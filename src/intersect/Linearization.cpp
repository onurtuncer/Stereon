// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import :Linearization;

/// The derivative products of a curve–surface intersection, through the
/// equations of its points (Linearization.cppm). The system of a point is
/// three by three, so it is solved by Cramer's rule on its columns C', -S_u
/// and -S_v, forward for the JVP and transposed for the VJP; the cross
/// products of the columns are the cofactors.
///
/// The curve enters through a side: its point and tangent at t, the
/// velocity of the point at fixed t for one binding, and the check of its
/// bindings. A line and a B-spline curve have one each; everything else is
/// written once over them.

namespace Stereon {

namespace {

/// The columns of F_z at a point, as C' = a, S_u and S_v, with the residual
/// and the determinant det F_z = C' . (S_u x S_v).
struct PointSystem {
    Vec3 Tangent;
    Vec3 Su;
    Vec3 Sv;
    Vec3 Residual;
    double Det;
    double Conditioning;
    /// One partial is zero within rounding against the other: a pole of
    /// the surface's parameters.
    bool AtPole;
};

/// A partial shorter than this fraction of the other is zero within the
/// rounding of its trigonometry: cos(pi/2) is 6e-17, not 0. Well below any
/// real surface's ratio of partials, well above that rounding.
constexpr double PoleRatio = 0x1p-50;

[[nodiscard]] Vec3 UnitAlong(const int k)
{
    return {k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0, k == 2 ? 1.0 : 0.0};
}

[[nodiscard]] Result<void> CheckParameter(const ParameterSet& set, const ParameterId parameter,
                                          const ParameterQuantity quantity,
                                          const std::string_view where)
{
    if (!set.Contains(parameter)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) + ": a binding's parameter is not of the set");
    }
    if (parameter.GetQuantity() != quantity) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) +
                        ": a coordinate is moved by a length and a weight by a dimensionless "
                        "parameter");
    }
    return {};
}

/// A line's side: dC/d(component) at t is a unit for the origin, t times a
/// unit for the direction.
struct LineSide {
    const Line& Geometry;
    std::span<const LineBinding> Bindings;

    [[nodiscard]] Point3 At(const double t) const { return Eval(Geometry, t); }
    [[nodiscard]] Vec3 Tangent(double /*t*/) const { return Geometry.Direction; }
    [[nodiscard]] static Vec3 Partial(const LineBinding& b, const double t)
    {
        const int k = static_cast<int>(b.Component);
        return k < 3 ? UnitAlong(k) : t * UnitAlong(k - 3);
    }
    [[nodiscard]] Result<void> Check(const ParameterSet& set, const std::string_view where) const
    {
        for (const LineBinding& b : Bindings) {
            STEREON_TRY_VOID(CheckParameter(set, b.Parameter, ParameterQuantity::Length, where));
        }
        return {};
    }
    [[nodiscard]] static bool Covers(double /*t*/) { return true; }
};

/// A B-spline curve's side: its derivative and its pole velocities, the
/// terms of EvalJvp (Stereon.Curve:PoleSensitivity).
struct SplineSide {
    const BSplineCurve& Geometry;
    std::span<const PoleBinding> Bindings;

    [[nodiscard]] Point3 At(const double t) const { return Eval(Geometry, t); }
    [[nodiscard]] Vec3 Tangent(const double t) const { return Deriv(Geometry, t, 1); }
    [[nodiscard]] Vec3 Partial(const PoleBinding& b, const double t) const
    {
        return PoleVelocity(Geometry, b, t);
    }
    [[nodiscard]] Result<void> Check(const ParameterSet& set, const std::string_view where) const
    {
        for (const PoleBinding& b : Bindings) {
            if (b.Pole >= Geometry.Poles.size()) {
                return Fail(ErrorCode::InvalidInput,
                            std::string(where) +
                                ": a binding names a pole the curve does not have");
            }
            const ParameterQuantity quantity = b.Component == PoleComponent::Weight
                                                   ? ParameterQuantity::Dimensionless
                                                   : ParameterQuantity::Length;
            STEREON_TRY_VOID(CheckParameter(set, b.Parameter, quantity, where));
        }
        return {};
    }
    [[nodiscard]] bool Covers(const double t) const
    {
        return std::isfinite(t) && Contains(Domain(Geometry), t);
    }
};

/// The surface's partials come from its own Deriv, so any surface fits. At
/// a pole of a sphere S_u vanishes, within rounding: the conditioning, the
/// sine of an angle, scales that out and stays near 1, so the pole is
/// detected on its own and the point is Indeterminate, a singularity of the
/// parameters, not of the geometry (ADR-0021, section 5).
template <typename Side>
[[nodiscard]] PointSystem SystemOf(const Side& curve, const Surface& surface,
                                   const CurveSurfacePoint& point)
{
    const Point2& uv = point.OnSurface;
    PointSystem s{curve.Tangent(point.Parameter),
                  Deriv(surface, uv, 1, 0),
                  Deriv(surface, uv, 0, 1),
                  curve.At(point.Parameter) - Eval(surface, uv),
                  0.0,
                  0.0,
                  false};
    const Vec3 normal = Cross(s.Su, s.Sv);
    s.Det = Dot(s.Tangent, normal);
    s.Conditioning = std::abs(s.Det) / (Length(s.Tangent) * Length(normal));
    const double shorter = std::min(Length(s.Su), Length(s.Sv));
    s.AtPole = !(shorter > PoleRatio * std::max(Length(s.Su), Length(s.Sv)));
    return s;
}

/// The surfaces whose defining values a SurfaceBinding can move.
[[nodiscard]] bool IsSupported(const Surface& surface)
{
    return surface.Target<Plane>() != nullptr || surface.Target<Sphere>() != nullptr ||
           surface.Target<Cylinder>() != nullptr || surface.Target<Cone>() != nullptr ||
           surface.Target<Torus>() != nullptr;
}

/// Whether the surface has the value a component moves: every surface its
/// origin, all but the plane a Radius, the torus a MinorRadius, the cone a
/// SemiAngle.
[[nodiscard]] bool Has(const Surface& surface, const SurfaceComponent component)
{
    switch (component) {
    case SurfaceComponent::Radius:
        return surface.Target<Plane>() == nullptr;
    case SurfaceComponent::MinorRadius:
        return surface.Target<Torus>() != nullptr;
    case SurfaceComponent::SemiAngle:
        return surface.Target<Cone>() != nullptr;
    default:
        return true;
    }
}

/// Whether a point has a derivative, and why not when it has none.
[[nodiscard]] std::pair<DerivativeStatus, std::string_view> Classify(const CurveSurfacePoint& point,
                                                                     const PointSystem& system)
{
    const double conditioning = system.Conditioning;
    if (point.Contact != ContactType::Crossing) {
        return {DerivativeStatus::Nondifferentiable,
                "a tangential contact: the points merge or vanish as the geometry moves"};
    }
    if (!point.Decided) {
        return {DerivativeStatus::Indeterminate,
                "a point the solver did not decide: its contact is not established"};
    }
    if (system.AtPole) {
        return {DerivativeStatus::Indeterminate,
                "a crossing at a pole of the surface's parameters: the angle there has no "
                "derivative"};
    }
    if (!(conditioning > 0.0) || !std::isfinite(conditioning)) { // STEREON-LINT-ALLOW(epsilon)
        return {DerivativeStatus::Indeterminate,
                "a decided crossing whose system is singular: at a pole of the surface's "
                "parameters, or tangent within rounding"};
    }
    if (conditioning < NearEventSine) {
        return {DerivativeStatus::NearEvent, "a crossing at a grazing angle"};
    }
    return {DerivativeStatus::Valid, {}};
}

/// Worse statuses rank higher: a result takes its worst point's.
[[nodiscard]] int Rank(const DerivativeStatus status)
{
    switch (status) {
    case DerivativeStatus::Valid:
        return 0;
    case DerivativeStatus::NearEvent:
        return 1;
    case DerivativeStatus::Indeterminate:
        return 2;
    case DerivativeStatus::Nondifferentiable:
        return 3;
    case DerivativeStatus::NotSupported:
        return 4;
    }
    return 4;
}

/// The systems of every point, the result's status and its diagnostics.
struct Analysis {
    std::vector<PointSystem> Systems;
    DerivativeStatus Status = DerivativeStatus::Valid;
    std::string Reason;
    double Residual = 0.0;
    double Conditioning = 1.0;
};

template <typename Side>
[[nodiscard]] Analysis Analyze(const Side& curve, const Surface& surface,
                               const CurveSurfaceIntersection& intersection)
{
    Analysis a;
    if (!intersection.Overlaps.empty()) {
        a.Status = DerivativeStatus::Nondifferentiable;
        a.Reason = "the curve lies in the surface: its points are not isolated";
        return a;
    }
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PointSystem s = SystemOf(curve, surface, intersection.Points[k]);
        const auto [status, why] = Classify(intersection.Points[k], s);
        if (Rank(status) > Rank(a.Status)) {
            a.Status = status;
            a.Reason = "point " + std::to_string(k) + ": " + std::string(why);
        }
        a.Residual = std::max(a.Residual, Length(s.Residual));
        a.Conditioning = std::min(a.Conditioning, s.Conditioning);
        a.Systems.push_back(s);
    }
    return a;
}

/// z' = (t', u', v') solving F_z z' = r, columns C', -S_u, -S_v.
[[nodiscard]] std::array<double, 3> Solve(const PointSystem& s, const Vec3& r)
{
    return {Dot(r, Cross(s.Su, s.Sv)) / s.Det, -Dot(s.Tangent, Cross(r, s.Sv)) / s.Det,
            -Dot(s.Tangent, Cross(s.Su, r)) / s.Det};
}

/// y solving F_z^T y = g.
[[nodiscard]] Vec3 SolveTransposed(const PointSystem& s, const std::array<double, 3>& g)
{
    return (g[0] * Cross(s.Su, s.Sv) + g[1] * Cross(s.Tangent, s.Sv) +
            g[2] * Cross(s.Su, s.Tangent)) /
           s.Det;
}

/// The unit radial direction cos u X + sin u Y of a surface of revolution.
template <typename S>
[[nodiscard]] Vec3 Radial(const S& surface, const double u)
{
    return std::cos(u) * surface.XAxis + std::sin(u) * surface.YAxis;
}

/// dS/d(component) at fixed (u, v) of a cone, O + (R + v sin a) radial(u) +
/// v cos a axis: the radius moves it along radial(u), the semi-angle by
/// v (cos a radial(u) - sin a axis).
[[nodiscard]] Vec3 ConePartial(const Cone& cone, const Point2& uv, const SurfaceComponent component)
{
    const Vec3 radial = Radial(cone, uv.X);
    if (component == SurfaceComponent::Radius) {
        return radial;
    }
    const double a = cone.SemiAngle;
    return uv.Y * (std::cos(a) * radial - std::sin(a) * Axis(cone));
}

/// dS/d(component) at fixed (u, v) of a torus, C + (R + r cos v) radial(u)
/// + r sin v axis: the major radius moves it along radial(u), the minor
/// one along the tube's normal cos v radial(u) + sin v axis.
[[nodiscard]] Vec3 TorusPartial(const Torus& torus, const Point2& uv,
                                const SurfaceComponent component)
{
    const Vec3 radial = Radial(torus, uv.X);
    if (component == SurfaceComponent::Radius) {
        return radial;
    }
    return std::cos(uv.Y) * radial + std::sin(uv.Y) * Axis(torus);
}

/// dS/d(component) at fixed (u, v), per unit of the component: the origin
/// moves every point alike; the radius moves a sphere's point along
/// (S - C) / r and a cylinder's along its radial direction; the cone's and
/// the torus's are ConePartial's and TorusPartial's.
[[nodiscard]] Vec3 SurfacePartial(const Surface& surface, const Point2& uv,
                                  const SurfaceComponent component)
{
    if (component == SurfaceComponent::OriginX || component == SurfaceComponent::OriginY ||
        component == SurfaceComponent::OriginZ) {
        return UnitAlong(static_cast<int>(component));
    }
    if (const Sphere* sphere = surface.Target<Sphere>()) {
        return (Eval(*sphere, uv) - sphere->Center) / sphere->Radius;
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        return (Eval(*cylinder, uv) - cylinder->Origin - uv.Y * Axis(*cylinder)) / cylinder->Radius;
    }
    if (const Cone* cone = surface.Target<Cone>()) {
        return ConePartial(*cone, uv, component);
    }
    if (const Torus* torus = surface.Target<Torus>()) {
        return TorusPartial(*torus, uv, component);
    }
    // A binding of a value the surface does not have was refused by
    // CheckSurfaceBindings; nothing moves here.
    return {};
}

[[nodiscard]] Result<void> CheckSurfaceBindings(const std::span<const SurfaceBinding> bindings,
                                                const Surface& surface, const ParameterSet& set,
                                                const std::string_view where)
{
    for (const SurfaceBinding& b : bindings) {
        const ParameterQuantity quantity = b.Component == SurfaceComponent::SemiAngle
                                               ? ParameterQuantity::Angle
                                               : ParameterQuantity::Length;
        STEREON_TRY_VOID(CheckParameter(set, b.Parameter, quantity, where));
        if (!Has(surface, b.Component)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) +
                            ": a binding of a value the surface does not have (a radius on a "
                            "plane, a minor radius off a torus, a semi-angle off a cone)");
        }
    }
    return {};
}

/// An intersection passed to the derivative API is an external primal
/// result. Check both its stored location and its current geometric residual
/// before using the inverse Jacobian; otherwise a stale result can produce a
/// plausible but incorrect sensitivity. The budget is for roundoff in the
/// geometry coordinates, not a geometric modelling tolerance.
[[nodiscard]] bool ConsistentContact(const Point3& stored, const Point3& onCurve,
                                     const Point3& onSurface)
{
    const auto size = [](const Point3& p) {
        return std::max({std::abs(p.X), std::abs(p.Y), std::abs(p.Z)});
    };
    const double scale = std::max({1.0, size(stored), size(onCurve), size(onSurface)});
    const double budget = 512.0 * std::numeric_limits<double>::epsilon() * scale;
    return std::isfinite(budget) && Distance(stored, onCurve) <= budget &&
           Distance(stored, onSurface) <= budget;
}

/// The curve's bindings, the surface's, and every point on the curve's
/// domain.
template <typename Side>
[[nodiscard]] Result<void>
CheckInput(const Side& curve, const std::span<const SurfaceBinding> surfaceBindings,
           const Surface& surface, const CurveSurfaceIntersection& intersection,
           const ParameterSet& set, const std::string_view where)
{
    STEREON_TRY_VOID(curve.Check(set, where));
    STEREON_TRY_VOID(CheckSurfaceBindings(surfaceBindings, surface, set, where));
    for (const CurveSurfacePoint& p : intersection.Points) {
        if (!std::isfinite(p.Parameter) || !std::isfinite(p.OnSurface.X) ||
            !std::isfinite(p.OnSurface.Y) || !curve.Covers(p.Parameter)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": point parameters are non-finite or out of domain");
        }
        if (!ConsistentContact(p.Position, curve.At(p.Parameter), Eval(surface, p.OnSurface))) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) +
                            ": stale intersection point does not lie on both current geometries");
        }
    }
    return {};
}

/// The answer for a surface no binding can move yet.
[[nodiscard]] Result<SensitivityResult> Unsupported(const std::string_view where)
{
    return SensitivityResult::Without(DerivativeStatus::NotSupported,
                                      std::string(where) +
                                          ": a curve against this surface type has no products "
                                          "yet; the elementary surfaces have theirs");
}

/// The regime: how many points, and how each meets the surface.
[[nodiscard]] std::uint64_t Signature(const CurveSurfaceIntersection& intersection)
{
    std::uint64_t hash = 14695981039346656037ull; // FNV-1a
    const auto mix = [&hash](const std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((value >> (8 * byte)) & 0xFFu)) * 1099511628211ull;
        }
    };
    mix(intersection.Points.size());
    mix(intersection.Overlaps.size());
    for (const CurveSurfacePoint& p : intersection.Points) {
        mix(static_cast<std::uint64_t>(p.Contact));
    }
    return hash;
}

/// The values as the analysis's status, unless the arithmetic overflowed,
/// which is Indeterminate rather than a derivative.
[[nodiscard]] Result<SensitivityResult> Published(const Analysis& a, std::vector<double> values,
                                                  const std::uint64_t signature,
                                                  const std::string_view where)
{
    for (std::size_t k = 0; k < values.size(); ++k) {
        if (!std::isfinite(values[k])) {
            return SensitivityResult::Without(DerivativeStatus::Indeterminate,
                                              std::string(where) +
                                                  ": the derivative overflowed: component " +
                                                  std::to_string(k) + " is not finite");
        }
    }
    STEREON_TRY(SensitivityResult result,
                SensitivityResult::Computed(a.Status, std::move(values), signature));
    result.SetDiagnostics(a.Residual, a.Conditioning);
    return result;
}

/// The point's velocity for the seed: C_p p' and S_p p'.
template <typename Side>
[[nodiscard]] std::pair<Vec3, Vec3>
Velocities(const Side& curve, const std::span<const SurfaceBinding> surfaceBindings,
           const Surface& surface, const SensitivityContext& sensitivity,
           const CurveSurfacePoint& point)
{
    const ParameterSet& set = sensitivity.GetParameters();
    const std::span<const double> seed = sensitivity.GetSeed();
    const auto rate = [&](const ParameterId p) { return seed[p.GetIndex()] * set.GetReference(p); };
    Vec3 onCurve{};
    Vec3 onSurface{};
    for (const auto& b : curve.Bindings) {
        onCurve += rate(b.Parameter) * curve.Partial(b, point.Parameter);
    }
    for (const SurfaceBinding& b : surfaceBindings) {
        onSurface += rate(b.Parameter) * SurfacePartial(surface, point.OnSurface, b.Component);
    }
    return {onCurve, onSurface};
}

template <typename Side>
[[nodiscard]] Result<CurveSurfaceLinearization>
LinearizeOn(const Side& curve, const Surface& surface, const CurveSurfacePoint& point)
{
    if (!std::isfinite(point.OnSurface.X) || !std::isfinite(point.OnSurface.Y) ||
        !std::isfinite(point.Parameter) || !curve.Covers(point.Parameter)) {
        return Fail(ErrorCode::InvalidInput,
                    "Linearize: the point's parameters must be finite and on the curve");
    }
    if (surface.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "Linearize: an empty surface");
    }
    const PointSystem s = SystemOf(curve, surface, point);
    return CurveSurfaceLinearization{s.Residual, s.Tangent,      s.Su,
                                     s.Sv,       s.Conditioning, Classify(point, s).first};
}

template <typename Side>
[[nodiscard]] Result<SensitivityResult> JvpOn(const Side& curve, const Surface& surface,
                                              const std::span<const SurfaceBinding> surfaceBindings,
                                              const CurveSurfaceIntersection& intersection,
                                              const SensitivityContext& sensitivity)
{
    if (surface.IsEmpty() || !IsSupported(surface)) {
        return Unsupported("IntersectJvp");
    }
    STEREON_TRY_VOID(CheckInput(curve, surfaceBindings, surface, intersection,
                                sensitivity.GetParameters(), "IntersectJvp"));
    if (sensitivity.GetSeed().size() != sensitivity.GetParameters().GetSize()) {
        return Fail(ErrorCode::InvalidInput, "IntersectJvp: the context has no seed");
    }
    const Analysis a = Analyze(curve, surface, intersection);
    if (!HasDerivative(a.Status)) {
        return SensitivityResult::Without(a.Status, "IntersectJvp: " + a.Reason);
    }
    std::vector<double> values;
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PointSystem& s = a.Systems[k];
        const auto [onCurve, onSurface] =
            Velocities(curve, surfaceBindings, surface, sensitivity, intersection.Points[k]);
        const std::array<double, 3> z = Solve(s, onSurface - onCurve);
        const Vec3 x = onCurve + z[0] * s.Tangent;
        values.insert(values.end(), {x.X, x.Y, x.Z, z[0], z[1], z[2]});
    }
    return Published(a, std::move(values), Signature(intersection), "IntersectJvp");
}

/// One point's share of the reverse product. With x' = C_p p' + C' t' and
/// F_z z' = S_p p' - C_p p', the pairing with (x_bar, z_bar) is
/// x_bar . C_p p' + y . (S_p p' - C_p p'), where
/// F_z^T y = (t_bar + C' . x_bar, u_bar, v_bar).
template <typename Side>
void AccumulatePoint(const Side& curve, const std::span<const SurfaceBinding> surfaceBindings,
                     const Surface& surface, const ParameterSet& set, const PointSystem& s,
                     const CurveSurfacePoint& point, const std::span<const double> bar,
                     std::vector<double>& values)
{
    const Vec3 xBar{bar[0], bar[1], bar[2]};
    const Vec3 y = SolveTransposed(s, {bar[3] + Dot(s.Tangent, xBar), bar[4], bar[5]});
    for (const auto& b : curve.Bindings) {
        values[b.Parameter.GetIndex()] +=
            set.GetReference(b.Parameter) * Dot(xBar - y, curve.Partial(b, point.Parameter));
    }
    for (const SurfaceBinding& b : surfaceBindings) {
        values[b.Parameter.GetIndex()] +=
            set.GetReference(b.Parameter) *
            Dot(y, SurfacePartial(surface, point.OnSurface, b.Component));
    }
}

template <typename Side>
[[nodiscard]] Result<SensitivityResult> VjpOn(const Side& curve, const Surface& surface,
                                              const std::span<const SurfaceBinding> surfaceBindings,
                                              const CurveSurfaceIntersection& intersection,
                                              const SensitivityContext& sensitivity)
{
    if (surface.IsEmpty() || !IsSupported(surface)) {
        return Unsupported("IntersectVjp");
    }
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY_VOID(
        CheckInput(curve, surfaceBindings, surface, intersection, set, "IntersectVjp"));
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != CurveSurfacePointValues * intersection.Points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "IntersectVjp: the adjoint must have six components per point");
    }
    const Analysis a = Analyze(curve, surface, intersection);
    if (!HasDerivative(a.Status)) {
        return SensitivityResult::Without(a.Status, "IntersectVjp: " + a.Reason);
    }
    std::vector<double> values(set.GetSize(), 0.0);
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        AccumulatePoint(curve, surfaceBindings, surface, set, a.Systems[k], intersection.Points[k],
                        adjoint.subspan(CurveSurfacePointValues * k, CurveSurfacePointValues),
                        values);
    }
    return Published(a, std::move(values), Signature(intersection), "IntersectVjp");
}

} // namespace

Result<CurveSurfaceLinearization> Linearize(const Line& line, const Surface& surface,
                                            const CurveSurfacePoint& point)
{
    return LinearizeOn(LineSide{line, {}}, surface, point);
}

Result<CurveSurfaceLinearization> Linearize(const BSplineCurve& curve, const Surface& surface,
                                            const CurveSurfacePoint& point)
{
    return LinearizeOn(SplineSide{curve, {}}, surface, point);
}

Result<SensitivityResult>
IntersectJvp(const Line& line, const std::span<const LineBinding> lineBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return JvpOn(LineSide{line, lineBindings}, surface, surfaceBindings, intersection, sensitivity);
}

Result<SensitivityResult>
IntersectJvp(const BSplineCurve& curve, const std::span<const PoleBinding> poleBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return JvpOn(SplineSide{curve, poleBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

Result<SensitivityResult>
IntersectVjp(const Line& line, const std::span<const LineBinding> lineBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return VjpOn(LineSide{line, lineBindings}, surface, surfaceBindings, intersection, sensitivity);
}

Result<SensitivityResult>
IntersectVjp(const BSplineCurve& curve, const std::span<const PoleBinding> poleBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return VjpOn(SplineSide{curve, poleBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

} // namespace Stereon
