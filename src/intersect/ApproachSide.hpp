// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// One curve of a closest approach (ClosestApproach.cppm), shared by the
// primal and the products: a line O + u D, a circle or an ellipse
// C + a cos u X + b sin u Y, or a polynomial B-spline curve, with its
// derivatives and bounds on them.
//
// Included in the purview of the closest approach's implementation units,
// after their imports, with <algorithm>, <cmath>, <cstdint>, <limits>,
// <string>, <string_view> and <vector> in their global module fragment. Everything is in an
// anonymous namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

struct ApproachSide {
    enum class Kind : std::uint8_t { Line, Circle, Ellipse, Spline };

    Kind Shape = Kind::Line;
    Point3 Origin; // a line's origin, a conic's centre, a B-spline's poles' box's centre
    Vec3 D;        // a line's direction
    Vec3 X;        // a conic's frame
    Vec3 Y;
    double Ra = 0.0; // a conic's radii along X and Y
    double Rb = 0.0;
    /// Bounds on |P'|, |P''| and |P'''| over the whole curve.
    double L1 = 0.0;
    double L2 = 0.0;
    double L3 = 0.0;
    /// A bound on |P - Origin| over the whole curve; infinite for a line.
    double Reach = std::numeric_limits<double>::infinity();
    /// A B-spline's curve, owned by the caller's Curve.
    const BSplineCurve* Spline = nullptr;
    /// A B-spline's interior knots where it is less than twice continuously
    /// differentiable (multiplicity at least degree - 1, every interior knot
    /// below degree 3): L1 and L2 bound phi's gradient's change, and L3 its
    /// Hessian's, only between them.
    std::vector<double> Cuts;

    [[nodiscard]] bool IsLine() const { return Shape == Kind::Line; }

    /// P(u) - Origin, formed from the frame.
    [[nodiscard]] Vec3 Offset(const double u) const
    {
        if (IsLine()) {
            return u * D;
        }
        if (Spline != nullptr) {
            return Eval(*Spline, u) - Origin;
        }
        return (Ra * std::cos(u)) * X + (Rb * std::sin(u)) * Y;
    }
    [[nodiscard]] Point3 At(const double u) const { return Origin + Offset(u); }
    [[nodiscard]] Vec3 D1(const double u) const
    {
        if (IsLine()) {
            return D;
        }
        if (Spline != nullptr) {
            return Deriv(*Spline, u, 1);
        }
        return (-Ra * std::sin(u)) * X + (Rb * std::cos(u)) * Y;
    }
    [[nodiscard]] Vec3 D2(const double u) const
    {
        if (IsLine()) {
            return Vec3{};
        }
        if (Spline != nullptr) {
            return Deriv(*Spline, u, 2);
        }
        return (-Ra * std::cos(u)) * X + (-Rb * std::sin(u)) * Y;
    }
};

[[nodiscard]] bool IsFiniteVec(const Vec3& v)
{
    return std::isfinite(v.X) && std::isfinite(v.Y) && std::isfinite(v.Z);
}

[[nodiscard]] ApproachSide ConicSideOf(const ApproachSide::Kind kind, const Point3& center,
                                       const double a, const double b, const Vec3& x, const Vec3& y)
{
    ApproachSide side;
    side.Shape = kind;
    side.Origin = center;
    side.X = x;
    side.Y = y;
    side.Ra = a;
    side.Rb = b;
    const double reach = std::max(a, b);
    side.L1 = side.L2 = side.L3 = side.Reach = reach;
    return side;
}

/// Bounds on the first three derivatives' lengths of a polynomial B-spline over its
/// domain: the largest pole of each derivative, itself a B-spline whose
/// poles p (P_i+1 - P_i) / (u_i+p+1 - u_i+1) bound it by the convex hull
/// property (a pole over an empty span multiplies a basis function that is
/// zero, and is dropped). Between the Cuts each derivative is continuous,
/// so the bounds are Lipschitz constants there. Widened by the rounding of
/// the differences and quotients.
[[nodiscard]] std::array<double, 3> DerivativeBounds(const BSplineCurve& curve)
{
    std::vector<Vec3> poles;
    poles.reserve(curve.Poles.size());
    for (const Point3& p : curve.Poles) {
        poles.push_back(p - curve.Poles.front());
    }
    std::vector<double> knots = curve.Knots;
    int degree = curve.Degree;
    std::array<double, 3> bounds{};
    for (double& bound : bounds) {
        if (degree == 0 || poles.size() < 2) {
            bound = 0.0;
            continue;
        }
        std::vector<Vec3> next;
        for (std::size_t i = 0; i + 1 < poles.size(); ++i) {
            const double span = knots[i + static_cast<std::size_t>(degree) + 1] - knots[i + 1];
            next.push_back(span > 0.0
                               ? (static_cast<double>(degree) / span) * (poles[i + 1] - poles[i])
                               : Vec3{});
            bound = std::max(bound, Length(next.back()));
        }
        bound *= 1.0 + 64.0 * static_cast<double>(degree) * std::numeric_limits<double>::epsilon();
        poles = std::move(next);
        knots = std::vector<double>(knots.begin() + 1, knots.end() - 1);
        --degree;
    }
    return bounds;
}

/// The interior knots of a B-spline where it is less than C2: multiplicity
/// at least degree - 1, and at least 1.
[[nodiscard]] std::vector<double> LowContinuityKnots(const BSplineCurve& curve)
{
    const Interval domain = Domain(curve);
    const auto least = static_cast<std::ptrdiff_t>(std::max(curve.Degree - 1, 1));
    std::vector<double> cuts;
    for (const double k : curve.Knots) {
        if (k > domain.GetLower() && k < domain.GetUpper() && (cuts.empty() || cuts.back() != k) &&
            std::ranges::count(curve.Knots, k) >= least) {
            cuts.push_back(k);
        }
    }
    return cuts;
}

/// The side of a polynomial B-spline: its poles' box's centre as origin, the
/// farthest pole as reach (the convex hull property), and its derivatives'
/// bounds. A rational one, unless every weight is the same, is not
/// supported: its derivatives are not bounded by its poles.
[[nodiscard]] Result<ApproachSide> SplineSideOf(const BSplineCurve& curve,
                                                const std::string_view where)
{
    const bool rational = std::ranges::any_of(
        curve.Weights, [&](const double w) { return w != curve.Weights.front(); });
    if (rational) {
        return Fail(ErrorCode::NotSupported,
                    std::string(where) +
                        ": a rational B-spline waits for bounds on its derivatives (ADR-0052, "
                        "addendum of 2026-11-12, decision 5); a polynomial one is supported");
    }
    ApproachSide side;
    side.Shape = ApproachSide::Kind::Spline;
    side.Spline = &curve;
    Point3 lo = curve.Poles.front();
    Point3 hi = lo;
    for (const Point3& p : curve.Poles) {
        if (!IsFiniteVec(p - Point3{})) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": a pole is not finite");
        }
        lo = Point3{std::min(lo.X, p.X), std::min(lo.Y, p.Y), std::min(lo.Z, p.Z)};
        hi = Point3{std::max(hi.X, p.X), std::max(hi.Y, p.Y), std::max(hi.Z, p.Z)};
    }
    side.Origin = Point3{0.5 * (lo.X + hi.X), 0.5 * (lo.Y + hi.Y), 0.5 * (lo.Z + hi.Z)};
    double reach = 0.0;
    for (const Point3& p : curve.Poles) {
        reach = std::max(reach, Distance(p, side.Origin));
    }
    side.Reach = reach * (1.0 + 64.0 * std::numeric_limits<double>::epsilon()) +
                 64.0 * std::numeric_limits<double>::epsilon() * Length(side.Origin - Point3{});
    const std::array<double, 3> bounds = DerivativeBounds(curve);
    side.L1 = bounds[0];
    side.L2 = bounds[1];
    side.L3 = bounds[2];
    side.Cuts = LowContinuityKnots(curve);
    return side;
}

/// The side of a line, a circle, an ellipse or a polynomial B-spline;
/// NotSupported for any other curve, InvalidInput for an empty one or a
/// degenerate line or conic.
[[nodiscard]] Result<ApproachSide> SideOf(const Curve& curve, const std::string_view where)
{
    if (curve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, std::string(where) + ": an empty curve");
    }
    if (const Line* line = curve.Target<Line>()) {
        const double length = Length(line->Direction);
        if (!(length > 0.0) || !std::isfinite(length) || !IsFiniteVec(line->Origin - Point3{})) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a line needs a finite origin and direction");
        }
        ApproachSide side;
        side.Origin = line->Origin;
        side.D = line->Direction;
        side.L1 = length;
        return side;
    }
    if (const Circle* circle = curve.Target<Circle>()) {
        return ConicSideOf(ApproachSide::Kind::Circle, circle->Center, circle->Radius,
                           circle->Radius, circle->XAxis, circle->YAxis);
    }
    if (const Ellipse* ellipse = curve.Target<Ellipse>()) {
        return ConicSideOf(ApproachSide::Kind::Ellipse, ellipse->Center, ellipse->MajorRadius,
                           ellipse->MinorRadius, ellipse->XAxis, ellipse->YAxis);
    }
    if (const BSplineCurve* spline = curve.Target<BSplineCurve>()) {
        return SplineSideOf(*spline, where);
    }
    return Fail(ErrorCode::NotSupported,
                std::string(where) +
                    ": lines, circles, ellipses and polynomial B-splines only (ADR-0052, "
                    "addendum of 2026-11-12)");
}

} // namespace
} // namespace Stereon::Detail
