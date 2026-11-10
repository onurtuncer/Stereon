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
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import :Linearization;

/// The linearisation of points where two plane curves meet (ADR-0021,
/// section 9; Linearization.cppm, PlanarLinearization): F(s, t) = A(s) -
/// B(t) = 0, two equations in two unknowns, solved by Cramer's rule forward
/// and transposed for the reverse product. Each curve's side gives its point,
/// its tangent and the velocity of its point at fixed parameter for each
/// binding; the B-spline's is the curve library's PoleVelocity on the curve
/// lifted into the plane z = 0, where it is the same curve.

namespace Stereon {

namespace {

/// The 2D cross product, the determinant of the columns (a, b).
[[nodiscard]] double Cross2(const Vec2& a, const Vec2& b)
{
    return a.X * b.Y - a.Y * b.X;
}

[[nodiscard]] Vec2 Unit2(const int k)
{
    return {k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0};
}

/// A plane B-spline as the 3D curve in the plane z = 0.
[[nodiscard]] BSplineCurve Lifted(const BSplineCurve2d& curve)
{
    BSplineCurve lifted{curve.Degree, curve.Knots, {}, curve.Weights};
    lifted.Poles.reserve(curve.Poles.size());
    for (const Point2& pole : curve.Poles) {
        lifted.Poles.push_back({pole.X, pole.Y, 0.0});
    }
    return lifted;
}

[[nodiscard]] bool Among(const PlaneCurveComponent c, const PlaneCurveComponent first,
                         const PlaneCurveComponent last)
{
    return static_cast<int>(c) >= static_cast<int>(first) &&
           static_cast<int>(c) <= static_cast<int>(last);
}

/// The curve library's pole component for a plane B-spline's.
[[nodiscard]] PoleComponent PoleOf(const PlaneCurveBinding& b)
{
    switch (b.Component) {
    case PlaneCurveComponent::PoleY:
        return PoleComponent::Y;
    case PlaneCurveComponent::PoleWeight:
        return PoleComponent::Weight;
    default:
        return PoleComponent::X;
    }
}

/// One curve's side of the system.
struct PlaneSide {
    const Curve2d& Geometry;
    std::span<const PlaneCurveBinding> Bindings;
    std::optional<BSplineCurve> Spline; // the B-spline lifted, for PoleVelocity

    PlaneSide(const Curve2d& geometry, const std::span<const PlaneCurveBinding> bindings)
        : Geometry(geometry)
        , Bindings(bindings)
    {
        if (const BSplineCurve2d* spline = geometry.Target<BSplineCurve2d>()) {
            Spline = Lifted(*spline);
        }
    }

    /// The curves whose defining values a binding can move.
    [[nodiscard]] bool IsSupported() const
    {
        return Geometry.Target<Line2d>() != nullptr || Geometry.Target<Circle2d>() != nullptr ||
               Geometry.Target<Ellipse2d>() != nullptr || Spline.has_value();
    }

    [[nodiscard]] bool Has(const PlaneCurveBinding& b) const
    {
        using C = PlaneCurveComponent;
        if (Geometry.Target<Line2d>() != nullptr) {
            return Among(b.Component, C::OriginX, C::DirectionY);
        }
        if (Geometry.Target<Circle2d>() != nullptr) {
            return Among(b.Component, C::CenterX, C::Radius);
        }
        if (Geometry.Target<Ellipse2d>() != nullptr) {
            return Among(b.Component, C::CenterX, C::CenterY) ||
                   Among(b.Component, C::MajorRadius, C::MinorRadius);
        }
        return Spline && Among(b.Component, C::PoleX, C::PoleWeight) &&
               b.Pole < Spline->Poles.size();
    }

    [[nodiscard]] Point2 At(const double t) const { return Eval(Geometry, t); }
    [[nodiscard]] Vec2 Tangent(const double t) const { return Deriv(Geometry, t, 1); }
    [[nodiscard]] bool Covers(const double t) const
    {
        return std::isfinite(t) && (!Spline || Contains(Domain(Geometry), t));
    }

    /// dC/d(component) at fixed t, per unit of the component.
    [[nodiscard]] Vec2 Partial(const PlaneCurveBinding& b, const double t) const
    {
        using C = PlaneCurveComponent;
        if (Geometry.Target<Line2d>() != nullptr) {
            const int k = static_cast<int>(b.Component);
            return k < 2 ? Unit2(k) : t * Unit2(k - 2);
        }
        if (Spline) {
            const Vec3 v = PoleVelocity(*Spline, PoleBinding{b.Parameter, b.Pole, PoleOf(b)}, t);
            return {v.X, v.Y};
        }
        if (b.Component == C::CenterX || b.Component == C::CenterY) {
            return Unit2(static_cast<int>(b.Component) - static_cast<int>(C::CenterX));
        }
        if (const Circle2d* circle = Geometry.Target<Circle2d>()) {
            return std::cos(t) * circle->XAxis + std::sin(t) * circle->YAxis;
        }
        const Ellipse2d& ellipse = *Geometry.Target<Ellipse2d>();
        return b.Component == C::MajorRadius ? std::cos(t) * ellipse.XAxis
                                             : std::sin(t) * ellipse.YAxis;
    }

    [[nodiscard]] Result<void> Check(const ParameterSet& set, const std::string_view where) const
    {
        for (const PlaneCurveBinding& b : Bindings) {
            if (!set.Contains(b.Parameter)) {
                return Fail(ErrorCode::InvalidInput,
                            std::string(where) + ": a binding's parameter is not of the set");
            }
            const ParameterQuantity wanted = b.Component == PlaneCurveComponent::PoleWeight
                                                 ? ParameterQuantity::Dimensionless
                                                 : ParameterQuantity::Length;
            if (b.Parameter.GetQuantity() != wanted) {
                return Fail(ErrorCode::InvalidInput,
                            std::string(where) +
                                ": a coordinate or a radius is moved by a length and a weight by "
                                "a dimensionless parameter");
            }
            if (!Has(b)) {
                return Fail(ErrorCode::InvalidInput,
                            std::string(where) +
                                ": a binding of a value the curve does not have (a line's "
                                "origin and direction, a circle's centre and radius, an "
                                "ellipse's centre and radii, a B-spline's poles)");
            }
        }
        return {};
    }
};

/// One point's system: the tangents, the residual and the determinant of
/// F_z = [A', -B'].
struct PlanarSystem {
    Vec2 TangentA;
    Vec2 TangentB;
    Vec2 Residual;
    double Det = 0.0;
    double Conditioning = 0.0;
};

[[nodiscard]] PlanarSystem SystemOf(const PlaneSide& a, const PlaneSide& b, const PlanarPoint& p)
{
    PlanarSystem s;
    s.TangentA = a.Tangent(p.ParameterA);
    s.TangentB = b.Tangent(p.ParameterB);
    s.Residual = a.At(p.ParameterA) - b.At(p.ParameterB);
    s.Det = -Cross2(s.TangentA, s.TangentB);
    s.Conditioning = std::abs(s.Det) / (Length(s.TangentA) * Length(s.TangentB));
    return s;
}

/// Whether a point has a derivative, and why not when it has none.
[[nodiscard]] std::pair<DerivativeStatus, std::string_view> Classify(const PlanarPoint& point,
                                                                     const PlanarSystem& s)
{
    if (point.Contact != ContactType::Crossing) {
        return {DerivativeStatus::Nondifferentiable,
                "a tangential contact: the points merge or vanish as the geometry moves"};
    }
    if (!point.Decided) {
        return {DerivativeStatus::Indeterminate,
                "a point the solver did not decide: its contact is not established"};
    }
    if (!(s.Conditioning > 0.0) || !std::isfinite(s.Conditioning)) { // STEREON-LINT-ALLOW(epsilon)
        return {DerivativeStatus::Indeterminate,
                "a decided crossing whose system is singular: tangent within rounding"};
    }
    if (s.Conditioning < NearEventSine) {
        return {DerivativeStatus::NearEvent, "a crossing at a grazing angle"};
    }
    return {DerivativeStatus::Valid, {}};
}

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

struct PlanarAnalysis {
    std::vector<PlanarSystem> Systems;
    DerivativeStatus Status = DerivativeStatus::Valid;
    std::string Reason;
    double Residual = 0.0;
    double Conditioning = 1.0;
};

[[nodiscard]] PlanarAnalysis Analyze(const PlaneSide& a, const PlaneSide& b,
                                     const PlanarIntersection& intersection)
{
    PlanarAnalysis analysis;
    if (!intersection.Overlaps.empty()) {
        analysis.Status = DerivativeStatus::Nondifferentiable;
        analysis.Reason = "the curves overlap: their points are not isolated";
        return analysis;
    }
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PlanarSystem s = SystemOf(a, b, intersection.Points[k]);
        const auto [status, why] = Classify(intersection.Points[k], s);
        if (Rank(status) > Rank(analysis.Status)) {
            analysis.Status = status;
            analysis.Reason = "point " + std::to_string(k) + ": " + std::string(why);
        }
        analysis.Residual = std::max(analysis.Residual, Length(s.Residual));
        analysis.Conditioning = std::min(analysis.Conditioning, s.Conditioning);
        analysis.Systems.push_back(s);
    }
    return analysis;
}

/// The stored point lies on both current curves to rounding: a stale
/// result gives a plausible but wrong derivative, so it is refused.
[[nodiscard]] bool Consistent(const Point2& stored, const Point2& onA, const Point2& onB)
{
    const auto size = [](const Point2& p) { return std::max(std::abs(p.X), std::abs(p.Y)); };
    const double scale = std::max({1.0, size(stored), size(onA), size(onB)});
    const double budget = 512.0 * std::numeric_limits<double>::epsilon() * scale;
    return std::isfinite(budget) && Distance(stored, onA) <= budget &&
           Distance(stored, onB) <= budget;
}

[[nodiscard]] Result<void> CheckInput(const PlaneSide& a, const PlaneSide& b,
                                      const PlanarIntersection& intersection,
                                      const ParameterSet& set, const std::string_view where)
{
    STEREON_TRY_VOID(a.Check(set, where));
    STEREON_TRY_VOID(b.Check(set, where));
    for (const PlanarPoint& p : intersection.Points) {
        if (!a.Covers(p.ParameterA) || !b.Covers(p.ParameterB)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": point parameters are non-finite or out of domain");
        }
        if (!Consistent(p.Position, a.At(p.ParameterA), b.At(p.ParameterB))) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) +
                            ": stale intersection point does not lie on both current curves");
        }
    }
    return {};
}

/// The regime: how many points and overlaps, and how each point meets.
[[nodiscard]] std::uint64_t Signature(const PlanarIntersection& intersection)
{
    std::uint64_t hash = 14695981039346656037ull; // FNV-1a
    const auto mix = [&hash](const std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((value >> (8 * byte)) & 0xFFu)) * 1099511628211ull;
        }
    };
    mix(intersection.Points.size());
    mix(intersection.Overlaps.size());
    for (const PlanarPoint& p : intersection.Points) {
        mix(static_cast<std::uint64_t>(p.Contact));
    }
    return hash;
}

[[nodiscard]] Result<SensitivityResult> Published(const PlanarAnalysis& a,
                                                  std::vector<double> values,
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

[[nodiscard]] Result<SensitivityResult> Unsupported(const std::string_view where)
{
    return SensitivityResult::Without(DerivativeStatus::NotSupported,
                                      std::string(where) +
                                          ": a hyperbola or a parabola binds nothing yet; "
                                          "lines, circles, ellipses and B-splines do");
}

/// The velocity of a side's point at fixed parameter, for the seed.
[[nodiscard]] Vec2 Velocity(const PlaneSide& side, const SensitivityContext& sensitivity,
                            const double t)
{
    const ParameterSet& set = sensitivity.GetParameters();
    const std::span<const double> seed = sensitivity.GetSeed();
    Vec2 v{};
    for (const PlaneCurveBinding& b : side.Bindings) {
        v += seed[b.Parameter.GetIndex()] * set.GetReference(b.Parameter) * side.Partial(b, t);
    }
    return v;
}

/// z' = (s', t') solving [A', -B'] z' = r.
[[nodiscard]] std::array<double, 2> Solve(const PlanarSystem& s, const Vec2& r)
{
    const Vec2 minusB = -1.0 * s.TangentB;
    return {Cross2(r, minusB) / s.Det, Cross2(s.TangentA, r) / s.Det};
}

/// y solving [A', -B']^T y = g: the rows A' . y = g0 and -B' . y = g1.
[[nodiscard]] Vec2 SolveTransposed(const PlanarSystem& s, const std::array<double, 2>& g)
{
    const Vec2 a = s.TangentA;
    const Vec2 b = -1.0 * s.TangentB;
    // [a.X a.Y; b.X b.Y] y = g, whose determinant is Cross2(a, b) = s.Det.
    return {(g[0] * b.Y - a.Y * g[1]) / s.Det, (a.X * g[1] - g[0] * b.X) / s.Det};
}

} // namespace

Result<PlanarLinearization> Linearize(const Curve2d& a, const Curve2d& b, const PlanarPoint& point)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "Linearize: an empty curve");
    }
    const PlaneSide sideA(a, {});
    const PlaneSide sideB(b, {});
    if (!sideA.Covers(point.ParameterA) || !sideB.Covers(point.ParameterB)) {
        return Fail(ErrorCode::InvalidInput,
                    "Linearize: the point's parameters must be finite and on the curves");
    }
    const PlanarSystem s = SystemOf(sideA, sideB, point);
    return PlanarLinearization{s.Residual, s.TangentA, s.TangentB, s.Conditioning,
                               Classify(point, s).first};
}

Result<SensitivityResult>
IntersectJvp(const Curve2d& a, const std::span<const PlaneCurveBinding> bindingsA, const Curve2d& b,
             const std::span<const PlaneCurveBinding> bindingsB,
             const PlanarIntersection& intersection, const SensitivityContext& sensitivity)
{
    const PlaneSide sideA(a, bindingsA);
    const PlaneSide sideB(b, bindingsB);
    if (a.IsEmpty() || b.IsEmpty() || !sideA.IsSupported() || !sideB.IsSupported()) {
        return Unsupported("IntersectJvp");
    }
    STEREON_TRY_VOID(
        CheckInput(sideA, sideB, intersection, sensitivity.GetParameters(), "IntersectJvp"));
    if (sensitivity.GetSeed().size() != sensitivity.GetParameters().GetSize()) {
        return Fail(ErrorCode::InvalidInput, "IntersectJvp: the context has no seed");
    }
    const PlanarAnalysis analysis = Analyze(sideA, sideB, intersection);
    if (!HasDerivative(analysis.Status)) {
        return SensitivityResult::Without(analysis.Status, "IntersectJvp: " + analysis.Reason);
    }
    std::vector<double> values;
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PlanarPoint& p = intersection.Points[k];
        const PlanarSystem& s = analysis.Systems[k];
        const Vec2 onA = Velocity(sideA, sensitivity, p.ParameterA);
        const Vec2 onB = Velocity(sideB, sensitivity, p.ParameterB);
        const std::array<double, 2> z = Solve(s, onB - onA);
        const Vec2 x = onA + z[0] * s.TangentA;
        values.insert(values.end(), {x.X, x.Y, z[0], z[1]});
    }
    return Published(analysis, std::move(values), Signature(intersection), "IntersectJvp");
}

/// The reverse product. With x' = A_p p' + A' s' and F_z z' = B_p p' -
/// A_p p', the pairing with (x_bar, s_bar, t_bar) is x_bar . A_p p' +
/// y . (B_p p' - A_p p'), where F_z^T y = (s_bar + A' . x_bar, t_bar).
Result<SensitivityResult>
IntersectVjp(const Curve2d& a, const std::span<const PlaneCurveBinding> bindingsA, const Curve2d& b,
             const std::span<const PlaneCurveBinding> bindingsB,
             const PlanarIntersection& intersection, const SensitivityContext& sensitivity)
{
    const PlaneSide sideA(a, bindingsA);
    const PlaneSide sideB(b, bindingsB);
    if (a.IsEmpty() || b.IsEmpty() || !sideA.IsSupported() || !sideB.IsSupported()) {
        return Unsupported("IntersectVjp");
    }
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY_VOID(CheckInput(sideA, sideB, intersection, set, "IntersectVjp"));
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != PlanarPointValues * intersection.Points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "IntersectVjp: the adjoint must have four components per point");
    }
    const PlanarAnalysis analysis = Analyze(sideA, sideB, intersection);
    if (!HasDerivative(analysis.Status)) {
        return SensitivityResult::Without(analysis.Status, "IntersectVjp: " + analysis.Reason);
    }
    std::vector<double> values(set.GetSize(), 0.0);
    for (std::size_t k = 0; k < intersection.Points.size(); ++k) {
        const PlanarPoint& p = intersection.Points[k];
        const PlanarSystem& s = analysis.Systems[k];
        const std::span<const double> bar =
            adjoint.subspan(PlanarPointValues * k, PlanarPointValues);
        const Vec2 xBar{bar[0], bar[1]};
        const Vec2 y = SolveTransposed(s, {bar[2] + Dot(s.TangentA, xBar), bar[3]});
        for (const PlaneCurveBinding& binding : sideA.Bindings) {
            values[binding.Parameter.GetIndex()] +=
                set.GetReference(binding.Parameter) *
                Dot(xBar - y, sideA.Partial(binding, p.ParameterA));
        }
        for (const PlaneCurveBinding& binding : sideB.Bindings) {
            values[binding.Parameter.GetIndex()] +=
                set.GetReference(binding.Parameter) * Dot(y, sideB.Partial(binding, p.ParameterB));
        }
    }
    return Published(analysis, std::move(values), Signature(intersection), "IntersectVjp");
}

} // namespace Stereon
