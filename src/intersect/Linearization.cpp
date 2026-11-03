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

/// The surface's partials come from its own Deriv, so any surface fits. At
/// a pole of a sphere S_u vanishes, within rounding: the conditioning, the
/// sine of an angle, scales that out and stays near 1, so the pole is
/// detected on its own and the point is Indeterminate, a singularity of the
/// parameters, not of the geometry (ADR-0021, section 5).
[[nodiscard]] PointSystem SystemOf(const Line& line, const Surface& surface,
                                   const CurveSurfacePoint& point)
{
    const Point2& uv = point.OnSurface;
    PointSystem s{line.Direction,
                  Deriv(surface, uv, 1, 0),
                  Deriv(surface, uv, 0, 1),
                  Eval(line, point.Parameter) - Eval(surface, uv),
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
           surface.Target<Cylinder>() != nullptr;
}

[[nodiscard]] bool HasRadius(const Surface& surface)
{
    return surface.Target<Sphere>() != nullptr || surface.Target<Cylinder>() != nullptr;
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

[[nodiscard]] Analysis Analyze(const Line& line, const Surface& surface,
                               const CurveSurfaceIntersection& intersection)
{
    Analysis a;
    if (!intersection.Overlaps.empty()) {
        a.Status = DerivativeStatus::Nondifferentiable;
        a.Reason = "the line lies in the surface: its points are not isolated";
        return a;
    }
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PointSystem s = SystemOf(line, surface, intersection.Points[k]);
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

[[nodiscard]] Vec3 UnitAlong(const int k)
{
    return {k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0, k == 2 ? 1.0 : 0.0};
}

/// dC/d(component) at t, per unit of the component: the origin moves every
/// point alike, the direction moves the point at t by t.
[[nodiscard]] Vec3 LinePartial(const LineComponent component, const double t)
{
    const int k = static_cast<int>(component);
    return k < 3 ? UnitAlong(k) : t * UnitAlong(k - 3);
}

/// dS/d(component) at fixed (u, v), per unit of the component: the origin
/// moves every point alike; the radius moves a sphere's point along
/// (S - C) / r and a cylinder's along its radial direction.
[[nodiscard]] Vec3 SurfacePartial(const Surface& surface, const Point2& uv,
                                  const SurfaceComponent component)
{
    if (component != SurfaceComponent::Radius) {
        return UnitAlong(static_cast<int>(component));
    }
    if (const Sphere* sphere = surface.Target<Sphere>()) {
        return (Eval(*sphere, uv) - sphere->Center) / sphere->Radius;
    }
    const Cylinder* cylinder = surface.Target<Cylinder>();
    return (Eval(*cylinder, uv) - cylinder->Origin - uv.Y * Axis(*cylinder)) / cylinder->Radius;
}

[[nodiscard]] Result<void> CheckBindings(const std::span<const LineBinding> lineBindings,
                                         const std::span<const SurfaceBinding> surfaceBindings,
                                         const Surface& surface, const ParameterSet& set,
                                         const std::string_view where)
{
    const auto check = [&](const ParameterId parameter) -> Result<void> {
        if (!set.Contains(parameter)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a binding's parameter is not of the set");
        }
        if (parameter.GetQuantity() != ParameterQuantity::Length) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a coordinate is moved by a length parameter");
        }
        return {};
    };
    for (const LineBinding& b : lineBindings) {
        STEREON_TRY_VOID(check(b.Parameter));
    }
    for (const SurfaceBinding& b : surfaceBindings) {
        STEREON_TRY_VOID(check(b.Parameter));
        if (b.Component == SurfaceComponent::Radius && !HasRadius(surface)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a radius binding on a surface without one");
        }
    }
    return {};
}

/// The answer for a surface no binding can move yet.
[[nodiscard]] Result<SensitivityResult> Unsupported(const std::string_view where)
{
    return SensitivityResult::Without(DerivativeStatus::NotSupported,
                                      std::string(where) +
                                          ": a line against this surface type has no products "
                                          "yet; planes, spheres and cylinders have theirs");
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
[[nodiscard]] std::pair<Vec3, Vec3>
Velocities(const std::span<const LineBinding> lineBindings,
           const std::span<const SurfaceBinding> surfaceBindings, const Surface& surface,
           const SensitivityContext& sensitivity, const CurveSurfacePoint& point)
{
    const ParameterSet& set = sensitivity.GetParameters();
    const std::span<const double> seed = sensitivity.GetSeed();
    const auto rate = [&](const ParameterId p) { return seed[p.GetIndex()] * set.GetReference(p); };
    Vec3 onCurve{};
    Vec3 onSurface{};
    for (const LineBinding& b : lineBindings) {
        onCurve += rate(b.Parameter) * LinePartial(b.Component, point.Parameter);
    }
    for (const SurfaceBinding& b : surfaceBindings) {
        onSurface += rate(b.Parameter) * SurfacePartial(surface, point.OnSurface, b.Component);
    }
    return {onCurve, onSurface};
}

} // namespace

Result<CurveSurfaceLinearization> Linearize(const Line& line, const Surface& surface,
                                            const CurveSurfacePoint& point)
{
    if (!std::isfinite(point.Parameter) || !std::isfinite(point.OnSurface.X) ||
        !std::isfinite(point.OnSurface.Y)) {
        return Fail(ErrorCode::InvalidInput, "Linearize: the point's parameters must be finite");
    }
    if (surface.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "Linearize: an empty surface");
    }
    const PointSystem s = SystemOf(line, surface, point);
    return CurveSurfaceLinearization{s.Residual, s.Tangent,      s.Su,
                                     s.Sv,       s.Conditioning, Classify(point, s).first};
}

Result<SensitivityResult>
IntersectJvp(const Line& line, const std::span<const LineBinding> lineBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    if (surface.IsEmpty() || !IsSupported(surface)) {
        return Unsupported("IntersectJvp");
    }
    STEREON_TRY_VOID(CheckBindings(lineBindings, surfaceBindings, surface,
                                   sensitivity.GetParameters(), "IntersectJvp"));
    if (sensitivity.GetSeed().size() != sensitivity.GetParameters().GetSize()) {
        return Fail(ErrorCode::InvalidInput, "IntersectJvp: the context has no seed");
    }
    const Analysis a = Analyze(line, surface, intersection);
    if (!HasDerivative(a.Status)) {
        return SensitivityResult::Without(a.Status, "IntersectJvp: " + a.Reason);
    }
    std::vector<double> values;
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PointSystem& s = a.Systems[k];
        const auto [onCurve, onSurface] =
            Velocities(lineBindings, surfaceBindings, surface, sensitivity, intersection.Points[k]);
        const std::array<double, 3> z = Solve(s, onSurface - onCurve);
        const Vec3 x = onCurve + z[0] * s.Tangent;
        values.insert(values.end(), {x.X, x.Y, x.Z, z[0], z[1], z[2]});
    }
    return Published(a, std::move(values), Signature(intersection), "IntersectJvp");
}

Result<SensitivityResult>
IntersectVjp(const Line& line, const std::span<const LineBinding> lineBindings,
             const Surface& surface, const std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity)
{
    if (surface.IsEmpty() || !IsSupported(surface)) {
        return Unsupported("IntersectVjp");
    }
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY_VOID(CheckBindings(lineBindings, surfaceBindings, surface, set, "IntersectVjp"));
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != CurveSurfacePointValues * intersection.Points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "IntersectVjp: the adjoint must have six components per point");
    }
    const Analysis a = Analyze(line, surface, intersection);
    if (!HasDerivative(a.Status)) {
        return SensitivityResult::Without(a.Status, "IntersectVjp: " + a.Reason);
    }
    // With x' = C_p p' + C' t' and F_z z' = S_p p' - C_p p', the pairing
    // with (x_bar, z_bar) is x_bar . C_p p' + y . (S_p p' - C_p p'), where
    // F_z^T y = (t_bar + C' . x_bar, u_bar, v_bar).
    std::vector<double> values(set.GetSize(), 0.0);
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PointSystem& s = a.Systems[k];
        const std::span<const double> bar =
            adjoint.subspan(CurveSurfacePointValues * k, CurveSurfacePointValues);
        const Vec3 xBar{bar[0], bar[1], bar[2]};
        const Vec3 y = SolveTransposed(s, {bar[3] + Dot(s.Tangent, xBar), bar[4], bar[5]});
        const CurveSurfacePoint& point = intersection.Points[k];
        for (const LineBinding& b : lineBindings) {
            values[b.Parameter.GetIndex()] +=
                set.GetReference(b.Parameter) *
                Dot(xBar - y, LinePartial(b.Component, point.Parameter));
        }
        for (const SurfaceBinding& b : surfaceBindings) {
            values[b.Parameter.GetIndex()] +=
                set.GetReference(b.Parameter) *
                Dot(y, SurfacePartial(surface, point.OnSurface, b.Component));
        }
    }
    return Published(a, std::move(values), Signature(intersection), "IntersectVjp");
}

} // namespace Stereon
