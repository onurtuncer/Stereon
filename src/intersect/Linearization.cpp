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
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
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

/// Whether a conic has the radius a component moves: a circle its Radius,
/// an ellipse its MajorRadius and MinorRadius, every conic its centre.
template <typename Conic>
[[nodiscard]] bool HasRadius(const ConicComponent component)
{
    if (component == ConicComponent::Radius) {
        return std::is_same_v<Conic, Circle>;
    }
    if (component == ConicComponent::MajorRadius || component == ConicComponent::MinorRadius) {
        return std::is_same_v<Conic, Ellipse>;
    }
    return true;
}

/// A circle's or an ellipse's side: its derivative, and dC/d(component) at
/// the angle t at fixed frame: a unit for a centre coordinate, cos t X +
/// sin t Y for a circle's radius, cos t X and sin t Y for an ellipse's
/// major and minor radii.
template <typename Conic>
struct ConicSide {
    const Conic& Geometry;
    std::span<const ConicBinding> Bindings;

    [[nodiscard]] Point3 At(const double t) const { return Eval(Geometry, t); }
    [[nodiscard]] Vec3 Tangent(const double t) const { return Deriv(Geometry, t, 1); }
    [[nodiscard]] Vec3 Partial(const ConicBinding& b, const double t) const
    {
        switch (b.Component) {
        case ConicComponent::Radius:
            return std::cos(t) * Geometry.XAxis + std::sin(t) * Geometry.YAxis;
        case ConicComponent::MajorRadius:
            return std::cos(t) * Geometry.XAxis;
        case ConicComponent::MinorRadius:
            return std::sin(t) * Geometry.YAxis;
        default:
            return UnitAlong(static_cast<int>(b.Component));
        }
    }
    [[nodiscard]] Result<void> Check(const ParameterSet& set, const std::string_view where) const
    {
        for (const ConicBinding& b : Bindings) {
            STEREON_TRY_VOID(CheckParameter(set, b.Parameter, ParameterQuantity::Length, where));
            if (!HasRadius<Conic>(b.Component)) {
                return Fail(ErrorCode::InvalidInput,
                            std::string(where) +
                                ": a binding of a radius the conic does not have (a major or "
                                "minor radius on a circle, a radius on an ellipse)");
            }
        }
        return {};
    }
    [[nodiscard]] static bool Covers(const double t) { return std::isfinite(t); }
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
           surface.Target<Torus>() != nullptr || surface.Target<Ellipsoid>() != nullptr ||
           surface.Target<BSplineSurface>() != nullptr;
}

[[nodiscard]] bool IsPoleComponent(const SurfaceComponent component)
{
    return component == SurfaceComponent::PoleX || component == SurfaceComponent::PoleY ||
           component == SurfaceComponent::PoleZ || component == SurfaceComponent::PoleWeight;
}

/// Whether an analytic surface has the value a component moves: every one
/// its origin, all but the plane and the ellipsoid a Radius, the torus a
/// MinorRadius, the cone a SemiAngle, the ellipsoid its three semi-axes;
/// none a pole.
[[nodiscard]] bool HasAnalytic(const Surface& surface, const SurfaceComponent component)
{
    if (IsPoleComponent(component)) {
        return false;
    }
    switch (component) {
    case SurfaceComponent::Radius:
        return surface.Target<Plane>() == nullptr && surface.Target<Ellipsoid>() == nullptr;
    case SurfaceComponent::MinorRadius:
        return surface.Target<Torus>() != nullptr;
    case SurfaceComponent::SemiAngle:
        return surface.Target<Cone>() != nullptr;
    case SurfaceComponent::SemiAxisA:
    case SurfaceComponent::SemiAxisB:
    case SurfaceComponent::SemiAxisC:
        return surface.Target<Ellipsoid>() != nullptr;
    default:
        return true;
    }
}

/// Whether the surface has the value a binding moves: an analytic surface
/// as HasAnalytic says, a B-spline surface its poles and weights and
/// nothing else.
[[nodiscard]] bool Has(const Surface& surface, const SurfaceBinding& binding)
{
    if (const BSplineSurface* spline = surface.Target<BSplineSurface>()) {
        return IsPoleComponent(binding.Component) && binding.Pole < spline->Poles.size();
    }
    return HasAnalytic(surface, binding.Component);
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

/// dS/d(component) at fixed (u, v) of an ellipsoid, C + cos v (A cos u X +
/// B sin u Y) + C sin v axis: each semi-axis moves the point along its own
/// axis by the factor it multiplies.
[[nodiscard]] Vec3 EllipsoidPartial(const Ellipsoid& ellipsoid, const Point2& uv,
                                    const SurfaceComponent component)
{
    switch (component) {
    case SurfaceComponent::SemiAxisA:
        return (std::cos(uv.Y) * std::cos(uv.X)) * ellipsoid.XAxis;
    case SurfaceComponent::SemiAxisB:
        return (std::cos(uv.Y) * std::sin(uv.X)) * ellipsoid.YAxis;
    default:
        return std::sin(uv.Y) * Axis(ellipsoid);
    }
}

/// dS/d(component) at fixed (u, v), per unit of the component: the origin
/// moves every point alike; the radius moves a sphere's point along
/// its unit normal and a cylinder's along its radial direction. Form these
/// locally: subtracting world positions loses accuracy after translation.
/// The cone's and the torus's are ConePartial's and TorusPartial's.
[[nodiscard]] Vec3 SurfacePartial(const Surface& surface, const Point2& uv,
                                  const SurfaceBinding& binding)
{
    const SurfaceComponent component = binding.Component;
    if (const BSplineSurface* spline = surface.Target<BSplineSurface>()) {
        const auto pole = static_cast<PoleComponent>(static_cast<int>(component) -
                                                     static_cast<int>(SurfaceComponent::PoleX));
        return PoleVelocity(*spline, binding.Pole, pole, uv);
    }
    if (component == SurfaceComponent::OriginX || component == SurfaceComponent::OriginY ||
        component == SurfaceComponent::OriginZ) {
        return UnitAlong(static_cast<int>(component));
    }
    if (const Sphere* sphere = surface.Target<Sphere>()) {
        return std::cos(uv.Y) * Radial(*sphere, uv.X) + std::sin(uv.Y) * Axis(*sphere);
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        return Radial(*cylinder, uv.X);
    }
    if (const Cone* cone = surface.Target<Cone>()) {
        return ConePartial(*cone, uv, component);
    }
    if (const Torus* torus = surface.Target<Torus>()) {
        return TorusPartial(*torus, uv, component);
    }
    if (const Ellipsoid* ellipsoid = surface.Target<Ellipsoid>()) {
        return EllipsoidPartial(*ellipsoid, uv, component);
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
        ParameterQuantity quantity = ParameterQuantity::Length;
        if (b.Component == SurfaceComponent::SemiAngle) {
            quantity = ParameterQuantity::Angle;
        } else if (b.Component == SurfaceComponent::PoleWeight) {
            quantity = ParameterQuantity::Dimensionless;
        }
        STEREON_TRY_VOID(CheckParameter(set, b.Parameter, quantity, where));
        if (!Has(surface, b)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) +
                            ": a binding of a value the surface does not have (a radius on a "
                            "plane or an ellipsoid, a minor radius off a torus, a semi-angle "
                            "off a cone, a semi-axis off an ellipsoid, a pole off a B-spline "
                            "surface or past its last, an analytic value on one)");
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
        onSurface += rate(b.Parameter) * SurfacePartial(surface, point.OnSurface, b);
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
            set.GetReference(b.Parameter) * Dot(y, SurfacePartial(surface, point.OnSurface, b));
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

Result<CurveSurfaceLinearization> Linearize(const Circle& circle, const Surface& surface,
                                            const CurveSurfacePoint& point)
{
    return LinearizeOn(ConicSide<Circle>{circle, {}}, surface, point);
}

Result<CurveSurfaceLinearization> Linearize(const Ellipse& ellipse, const Surface& surface,
                                            const CurveSurfacePoint& point)
{
    return LinearizeOn(ConicSide<Ellipse>{ellipse, {}}, surface, point);
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

Result<SensitivityResult>
IntersectJvp(const Circle& circle, const std::span<const ConicBinding> conicBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return JvpOn(ConicSide<Circle>{circle, conicBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

Result<SensitivityResult>
IntersectJvp(const Ellipse& ellipse, const std::span<const ConicBinding> conicBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return JvpOn(ConicSide<Ellipse>{ellipse, conicBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

Result<SensitivityResult>
IntersectVjp(const Circle& circle, const std::span<const ConicBinding> conicBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return VjpOn(ConicSide<Circle>{circle, conicBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

Result<SensitivityResult>
IntersectVjp(const Ellipse& ellipse, const std::span<const ConicBinding> conicBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    return VjpOn(ConicSide<Ellipse>{ellipse, conicBindings}, surface, surfaceBindings, intersection,
                 sensitivity);
}

// The branch points of a surface–surface intersection (Linearization.cppm,
// BranchPointLinearization): four equations in (u_A, v_A, u_B, v_B), three
// from S_A = S_B and the normal-plane tracking constraint, solved by
// Gaussian elimination with partial pivoting, forward for the JVP and
// transposed for the VJP.

namespace {

using Row4 = std::array<double, 4>;
using Matrix4 = std::array<Row4, 4>;

/// Newton's method from a certified branch starts within its bound of the
/// point and converges quadratically: a few steps reach rounding. A cap,
/// not a tolerance.
constexpr int NewtonSteps = 8;

/// m and r reduced to upper triangular form, the pivot of each column the
/// largest left below the diagonal.
void Eliminate(Matrix4& m, Row4& r)
{
    for (std::size_t c = 0; c < 4; ++c) {
        std::size_t pivot = c;
        for (std::size_t i = c + 1; i < 4; ++i) {
            pivot = std::abs(m[i][c]) > std::abs(m[pivot][c]) ? i : pivot;
        }
        std::swap(m[c], m[pivot]);
        std::swap(r[c], r[pivot]);
        for (std::size_t i = c + 1; i < 4; ++i) {
            const double f = m[i][c] / m[c][c];
            for (std::size_t j = c; j < 4; ++j) {
                m[i][j] -= f * m[c][j];
            }
            r[i] -= f * r[c];
        }
    }
}

/// x solving m x = r; not finite when m is singular.
[[nodiscard]] Row4 Solve4(Matrix4 m, Row4 r)
{
    Eliminate(m, r);
    Row4 x{};
    for (std::size_t c = 4; c-- > 0;) {
        double sum = r[c];
        for (std::size_t j = c + 1; j < 4; ++j) {
            sum -= m[c][j] * x[j];
        }
        x[c] = sum / m[c][c];
    }
    return x;
}

[[nodiscard]] Matrix4 Transposed(const Matrix4& m)
{
    Matrix4 t{};
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            t[j][i] = m[i][j];
        }
    }
    return t;
}

/// The system of a branch point at its parameters: the position on A, the
/// residual S_A - S_B, the unit tangent T naming the tracking plane, the
/// partials, and the regularity.
struct BranchSystem {
    Point2 OnA;
    Point2 OnB;
    Point3 Position;
    Vec3 Residual;
    Vec3 Tangent;
    Vec3 SuA;
    Vec3 SvA;
    Vec3 SuB;
    Vec3 SvB;
    double Conditioning = 0.0;
    bool AtPole = false;
};

[[nodiscard]] bool AtPoleOf(const Vec3& su, const Vec3& sv)
{
    const double shorter = std::min(Length(su), Length(sv));
    return !(shorter > PoleRatio * std::max(Length(su), Length(sv)));
}

/// The system at (onA, onB), for the tangent already in s.
void Fill(const Surface& a, const Surface& b, const Point2& onA, const Point2& onB, BranchSystem& s)
{
    s.OnA = onA;
    s.OnB = onB;
    s.Position = Eval(a, onA);
    s.Residual = s.Position - Eval(b, onB);
    s.SuA = Deriv(a, onA, 1, 0);
    s.SvA = Deriv(a, onA, 0, 1);
    s.SuB = Deriv(b, onB, 1, 0);
    s.SvB = Deriv(b, onB, 0, 1);
    const Vec3 normalA = Cross(s.SuA, s.SvA);
    const Vec3 normalB = Cross(s.SuB, s.SvB);
    s.AtPole = AtPoleOf(s.SuA, s.SvA) || AtPoleOf(s.SuB, s.SvB);
    s.Conditioning =
        std::abs(Dot(s.Tangent, Cross(normalA, normalB))) / (Length(normalA) * Length(normalB));
}

/// F_z: rows S_Au, S_Av, -S_Bu, -S_Bv by coordinate, then the tracking row.
[[nodiscard]] Matrix4 MatrixOf(const BranchSystem& s)
{
    return {{{s.SuA.X, s.SvA.X, -s.SuB.X, -s.SvB.X},
             {s.SuA.Y, s.SvA.Y, -s.SuB.Y, -s.SvB.Y},
             {s.SuA.Z, s.SvA.Z, -s.SuB.Z, -s.SvB.Z},
             {Dot(s.Tangent, s.SuA), Dot(s.Tangent, s.SvA), 0.0, 0.0}}};
}

[[nodiscard]] bool AllFinite(const Row4& x)
{
    return std::ranges::all_of(x, [](const double v) { return std::isfinite(v); });
}

/// Newton's method on S_A - S_B = 0 and T . (S_A - x0) = 0 from the
/// parameters in s, which it leaves at the last iterate.
void Refine(const Surface& a, const Surface& b, const Point3& x0, BranchSystem& s)
{
    for (int step = 0; step < NewtonSteps; ++step) {
        const double tracking = Dot(s.Tangent, s.Position - x0);
        const Row4 dz =
            Solve4(MatrixOf(s), {-s.Residual.X, -s.Residual.Y, -s.Residual.Z, -tracking});
        if (!AllFinite(dz)) {
            return;
        }
        const Point2 onA{s.OnA.X + dz[0], s.OnA.Y + dz[1]};
        const Point2 onB{s.OnB.X + dz[2], s.OnB.Y + dz[3]};
        const bool settled = onA == s.OnA && onB == s.OnB;
        Fill(a, b, onA, onB, s);
        if (settled) {
            return;
        }
    }
}

/// Whether the refined point lies on both surfaces and in the tracking
/// plane to rounding, and near the branch's point: within its bound of
/// both surfaces, so within twice that over the sine of their angle of
/// the curve where they meet.
[[nodiscard]] bool OnBoth(const BranchSystem& s, const Point3& x0, const double bound)
{
    const auto size = [](const Point3& p) {
        return std::max({std::abs(p.X), std::abs(p.Y), std::abs(p.Z)});
    };
    const double scale = std::max({1.0, size(x0), size(s.Position)});
    const double budget = 512.0 * std::numeric_limits<double>::epsilon() * scale;
    const double reach = 2.0 * (bound + budget) / std::max(s.Conditioning, NearEventSine);
    return std::isfinite(budget) && Length(s.Residual) <= budget &&
           std::abs(Dot(s.Tangent, s.Position - x0)) <= budget && Distance(s.Position, x0) <= reach;
}

/// Whether a branch point has a derivative, and why not when it has none.
[[nodiscard]] std::pair<DerivativeStatus, std::string_view>
ClassifyBranch(const IntersectionBranch& branch, const BranchSystem& s)
{
    if (branch.Contact != ContactType::Crossing) {
        return {DerivativeStatus::Nondifferentiable,
                "a tangential or overlap branch: the surfaces touch along it, and it splits or "
                "vanishes as they move"};
    }
    if (!branch.IsCertified()) {
        return {DerivativeStatus::Indeterminate,
                "a branch the solver did not certify: its point is not established"};
    }
    if (s.AtPole) {
        return {DerivativeStatus::Indeterminate,
                "a point at a pole of a surface's parameters: the angle there has no derivative"};
    }
    if (!(s.Conditioning > 0.0) || !std::isfinite(s.Conditioning)) { // STEREON-LINT-ALLOW(epsilon)
        return {DerivativeStatus::Indeterminate,
                "a crossing whose system is singular: tangent within rounding, or the branch "
                "inside its tracking plane"};
    }
    if (s.Conditioning < NearEventSine) {
        return {DerivativeStatus::NearEvent, "a crossing at a grazing angle"};
    }
    return {DerivativeStatus::Valid, {}};
}

/// A branch point: its system and its status.
struct BranchPoint {
    BranchSystem System;
    DerivativeStatus Status = DerivativeStatus::Valid;
    std::string_view Reason;
};

[[nodiscard]] bool OnBranch(const IntersectionBranch& branch, const double t)
{
    return std::isfinite(t) && !branch.Curve3d.IsEmpty() && !branch.OnA.IsEmpty() &&
           !branch.OnB.IsEmpty() && Contains(branch.Range, t) &&
           Contains(Domain(branch.Curve3d), t) && Contains(Domain(branch.OnA), t) &&
           Contains(Domain(branch.OnB), t);
}

/// The branch's point at t, refined onto both surfaces when it has a
/// derivative.
[[nodiscard]] Result<BranchPoint> PointOn(const Surface& a, const Surface& b,
                                          const IntersectionBranch& branch, const double t,
                                          const std::string_view where)
{
    if (!OnBranch(branch, t)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) + ": the sample's parameter is not on its branch");
    }
    const Point3 x0 = Eval(branch.Curve3d, t);
    const Vec3 direction = Deriv(branch.Curve3d, t, 1);
    if (!(Length(direction) > 0.0) || !std::isfinite(Length(direction))) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) + ": the branch has no tangent at the sample");
    }
    BranchPoint p;
    p.System.Tangent = direction / Length(direction);
    Fill(a, b, Eval(branch.OnA, t), Eval(branch.OnB, t), p.System);
    std::tie(p.Status, p.Reason) = ClassifyBranch(branch, p.System);
    if (!HasDerivative(p.Status)) {
        return p;
    }
    Refine(a, b, x0, p.System);
    if (OnBoth(p.System, x0, branch.Bound)) {
        std::tie(p.Status, p.Reason) = ClassifyBranch(branch, p.System);
        return p;
    }
    if (p.System.Conditioning < NearEventSine) {
        p.Status = DerivativeStatus::Indeterminate;
        p.Reason = "a grazing crossing Newton's method could not bring onto both surfaces";
        return p;
    }
    return Fail(ErrorCode::InvalidInput,
                std::string(where) +
                    ": stale branch: its point does not lie on both current surfaces within "
                    "its bound");
}

/// The samples' systems, the result's status and its diagnostics.
struct BranchAnalysis {
    std::vector<BranchSystem> Systems;
    Analysis Summary;
};

[[nodiscard]] Result<BranchAnalysis> AnalyzeBranches(const Surface& a, const Surface& b,
                                                     const SurfaceSurfaceIntersection& intersection,
                                                     const std::span<const BranchSample> samples,
                                                     const std::string_view where)
{
    BranchAnalysis an;
    if (!intersection.Overlaps.empty()) {
        an.Summary.Status = DerivativeStatus::Nondifferentiable;
        an.Summary.Reason = "the surfaces coincide over a region: its boundary is not isolated";
        return an;
    }
    for (std::size_t k = 0; k < samples.size(); ++k) {
        if (samples[k].Branch >= intersection.Branches.size()) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a sample names a branch the result does not have");
        }
        STEREON_TRY(const BranchPoint p, PointOn(a, b, intersection.Branches[samples[k].Branch],
                                                 samples[k].Parameter, where));
        if (Rank(p.Status) > Rank(an.Summary.Status)) {
            an.Summary.Status = p.Status;
            an.Summary.Reason = "sample " + std::to_string(k) + ": " + std::string(p.Reason);
        }
        an.Summary.Residual = std::max(an.Summary.Residual, Length(p.System.Residual));
        an.Summary.Conditioning = std::min(an.Summary.Conditioning, p.System.Conditioning);
        an.Systems.push_back(p.System);
    }
    return an;
}

/// Both surfaces' bindings, or NotSupported as a status (nullopt when
/// both are supported).
[[nodiscard]] Result<std::optional<SensitivityResult>>
CheckSurfaces(const Surface& a, const std::span<const SurfaceBinding> bindingsA, const Surface& b,
              const std::span<const SurfaceBinding> bindingsB, const ParameterSet& set,
              const std::string_view where)
{
    if (a.IsEmpty() || b.IsEmpty() || !IsSupported(a) || !IsSupported(b)) {
        STEREON_TRY(SensitivityResult unsupported,
                    SensitivityResult::Without(
                        DerivativeStatus::NotSupported,
                        std::string(where) +
                            ": two surfaces have products only when both are planes, "
                            "cylinders, cones, spheres, tori or ellipsoids"));
        return std::optional<SensitivityResult>(std::move(unsupported));
    }
    STEREON_TRY_VOID(CheckSurfaceBindings(bindingsA, a, set, where));
    STEREON_TRY_VOID(CheckSurfaceBindings(bindingsB, b, set, where));
    return std::optional<SensitivityResult>{};
}

/// S_p p' at uv for the seed.
[[nodiscard]] Vec3 SurfaceVelocity(const Surface& surface,
                                   const std::span<const SurfaceBinding> bindings, const Point2& uv,
                                   const SensitivityContext& sensitivity)
{
    const ParameterSet& set = sensitivity.GetParameters();
    const std::span<const double> seed = sensitivity.GetSeed();
    Vec3 velocity{};
    for (const SurfaceBinding& b : bindings) {
        velocity += seed[b.Parameter.GetIndex()] * set.GetReference(b.Parameter) *
                    SurfacePartial(surface, uv, b);
    }
    return velocity;
}

/// The regime: the result's components, its branches' contacts and the
/// branches sampled.
[[nodiscard]] std::uint64_t BranchSignature(const SurfaceSurfaceIntersection& intersection,
                                            const std::span<const BranchSample> samples)
{
    std::uint64_t hash = 14695981039346656037ull; // FNV-1a
    const auto mix = [&hash](const std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((value >> (8 * byte)) & 0xFFu)) * 1099511628211ull;
        }
    };
    mix(intersection.Points.size());
    mix(intersection.Branches.size());
    mix(intersection.Overlaps.size());
    for (const IntersectionBranch& branch : intersection.Branches) {
        mix(static_cast<std::uint64_t>(branch.Contact));
    }
    for (const BranchSample& sample : samples) {
        mix(sample.Branch);
    }
    return hash;
}

/// One sample's forward values: F_z z' = (S_Bp - S_Ap, -T . S_Ap) p', and
/// x' = S_Ap p' + S_Au u_A' + S_Av v_A'.
void AppendBranchJvp(const BranchSystem& s, const Vec3& onA, const Vec3& onB,
                     std::vector<double>& values)
{
    const Vec3 r = onB - onA;
    const Row4 z = Solve4(MatrixOf(s), {r.X, r.Y, r.Z, -Dot(s.Tangent, onA)});
    const Vec3 x = onA + z[0] * s.SuA + z[1] * s.SvA;
    values.insert(values.end(), {x.X, x.Y, x.Z, z[0], z[1], z[2], z[3]});
}

/// One sample's share of the reverse product. The pairing of (x_bar,
/// z_bar) with the forward values is x_bar . S_Ap p' + y . r, where
/// F_z^T y = (S_Au . x_bar + z_bar_0, S_Av . x_bar + z_bar_1, z_bar_2,
/// z_bar_3) and r the forward right-hand side; S_Ap p' takes x_bar - y_xyz
/// - y_4 T, S_Bp p' takes y_xyz.
void AccumulateBranch(const Surface& a, const std::span<const SurfaceBinding> bindingsA,
                      const Surface& b, const std::span<const SurfaceBinding> bindingsB,
                      const ParameterSet& set, const BranchSystem& s,
                      const std::span<const double> bar, std::vector<double>& values)
{
    const Vec3 xBar{bar[0], bar[1], bar[2]};
    const Row4 y = Solve4(Transposed(MatrixOf(s)),
                          {Dot(s.SuA, xBar) + bar[3], Dot(s.SvA, xBar) + bar[4], bar[5], bar[6]});
    const Vec3 onB{y[0], y[1], y[2]};
    const Vec3 onA = xBar - onB - y[3] * s.Tangent;
    for (const SurfaceBinding& binding : bindingsA) {
        values[binding.Parameter.GetIndex()] +=
            set.GetReference(binding.Parameter) * Dot(onA, SurfacePartial(a, s.OnA, binding));
    }
    for (const SurfaceBinding& binding : bindingsB) {
        values[binding.Parameter.GetIndex()] +=
            set.GetReference(binding.Parameter) * Dot(onB, SurfacePartial(b, s.OnB, binding));
    }
}

} // namespace

Result<BranchPointLinearization> Linearize(const Surface& a, const Surface& b,
                                           const IntersectionBranch& branch, const double t)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "Linearize: an empty surface");
    }
    STEREON_TRY(const BranchPoint p, PointOn(a, b, branch, t, "Linearize"));
    const BranchSystem& s = p.System;
    return BranchPointLinearization{s.Position, s.OnA, s.OnB, s.Residual,     s.Tangent, s.SuA,
                                    s.SvA,      s.SuB, s.SvB, s.Conditioning, p.Status};
}

Result<SensitivityResult>
IntersectJvp(const Surface& a, const std::span<const SurfaceBinding> bindingsA, const Surface& b,
             const std::span<const SurfaceBinding> bindingsB,
             const SurfaceSurfaceIntersection& intersection,
             const std::span<const BranchSample> samples, const SensitivityContext& sensitivity)
{
    STEREON_TRY(
        std::optional<SensitivityResult> unsupported,
        CheckSurfaces(a, bindingsA, b, bindingsB, sensitivity.GetParameters(), "IntersectJvp"));
    if (unsupported) {
        return *std::move(unsupported);
    }
    if (sensitivity.GetSeed().size() != sensitivity.GetParameters().GetSize()) {
        return Fail(ErrorCode::InvalidInput, "IntersectJvp: the context has no seed");
    }
    STEREON_TRY(const BranchAnalysis an,
                AnalyzeBranches(a, b, intersection, samples, "IntersectJvp"));
    if (!HasDerivative(an.Summary.Status)) {
        return SensitivityResult::Without(an.Summary.Status, "IntersectJvp: " + an.Summary.Reason);
    }
    std::vector<double> values;
    for (const BranchSystem& s : an.Systems) {
        AppendBranchJvp(s, SurfaceVelocity(a, bindingsA, s.OnA, sensitivity),
                        SurfaceVelocity(b, bindingsB, s.OnB, sensitivity), values);
    }
    return Published(an.Summary, std::move(values), BranchSignature(intersection, samples),
                     "IntersectJvp");
}

Result<SensitivityResult>
IntersectVjp(const Surface& a, const std::span<const SurfaceBinding> bindingsA, const Surface& b,
             const std::span<const SurfaceBinding> bindingsB,
             const SurfaceSurfaceIntersection& intersection,
             const std::span<const BranchSample> samples, const SensitivityContext& sensitivity)
{
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY(std::optional<SensitivityResult> unsupported,
                CheckSurfaces(a, bindingsA, b, bindingsB, set, "IntersectVjp"));
    if (unsupported) {
        return *std::move(unsupported);
    }
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != BranchPointValues * samples.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "IntersectVjp: the adjoint must have seven components per sample");
    }
    STEREON_TRY(const BranchAnalysis an,
                AnalyzeBranches(a, b, intersection, samples, "IntersectVjp"));
    if (!HasDerivative(an.Summary.Status)) {
        return SensitivityResult::Without(an.Summary.Status, "IntersectVjp: " + an.Summary.Reason);
    }
    std::vector<double> values(set.GetSize(), 0.0);
    for (std::size_t k = 0; k < an.Systems.size(); ++k) {
        AccumulateBranch(a, bindingsA, b, bindingsB, set, an.Systems[k],
                         adjoint.subspan(BranchPointValues * k, BranchPointValues), values);
    }
    return Published(an.Summary, std::move(values), BranchSignature(intersection, samples),
                     "IntersectVjp");
}

} // namespace Stereon
