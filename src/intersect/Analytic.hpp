// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the analytic handlers share once a configuration is decided: the
// nominal geometry of a branch, its range clipped to the requested domains
// with outward rounding, parameter-space boxes, the swap of a result's two
// sides, and numbers in messages.
//
// Included in the purview of the handlers' implementation units, after
// `import :Dispatch;`, with <array>, <charconv>, <cmath>, <cstddef>,
// <limits>, <string>, <string_view> and <utility> in their global module
// fragment, and <Stereon/Core/Try.hpp>. Everything is in an anonymous
// namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

inline constexpr double Infinity = std::numeric_limits<double>::infinity();

/// x in a message, shortest round-trip form.
[[nodiscard]] inline std::string Number(const double x)
{
    std::array<char, 32> buffer{};
    const std::to_chars_result r = std::to_chars(buffer.data(), buffer.data() + buffer.size(), x);
    return std::string(buffer.data(), r.ptr);
}

/// c - reach to c + reach, outward-rounded.
[[nodiscard]] inline Interval Around(const double c, const double reach)
{
    return Hull(Interval{c} - Interval{reach}, Interval{c} + Interval{reach});
}

/// The quadratic A t² + 2 B t + C = 0 that a line cuts out of a conic or a
/// quadric; its discriminant is B² - A C.
struct LineQuadratic {
    double A;
    double B;
    double C;
};

/// The roots, in order, for a non-negative discriminant: the one farther
/// from zero by the formula, the other from the product of the roots, so
/// that neither suffers the cancellation of the formula alone.
[[nodiscard]] inline std::array<double, 2> RootsOf(const LineQuadratic& q,
                                                   const double discriminant)
{
    const double root = std::sqrt(discriminant);
    const double head = -(q.B + (q.B < 0.0 ? -root : root));
    const double t1 = head / q.A;
    const double t2 = head == 0.0 ? t1 : q.C / head;
    return {std::min(t1, t2), std::max(t1, t2)};
}

/// The parameters t at which origin + t direction lies within range along
/// one axis, outward-rounded: Everything or Empty when the line runs along
/// the axis, inside or outside the range.
[[nodiscard]] inline Interval AxisRange(const double origin, const double direction,
                                        const Interval& range)
{
    if (direction == 0.0) {
        return Contains(range, origin) ? Interval::Everything() : Interval::Empty();
    }
    const auto at = [&](const double bound) {
        return (Interval{bound} - Interval{origin}) / Interval{direction};
    };
    const double lowerBound = direction > 0.0 ? range.GetLower() : range.GetUpper();
    const double upperBound = direction > 0.0 ? range.GetUpper() : range.GetLower();
    const double lower = std::isinf(lowerBound) ? -Infinity : at(lowerBound).GetLower();
    const double upper = std::isinf(upperBound) ? Infinity : at(upperBound).GetUpper();
    return Interval{lower, upper};
}

/// The parameters of the line inside the box, outward-rounded.
[[nodiscard]] inline Interval ClipLine(const Line2d& line, const Box2& box)
{
    return Intersection(AxisRange(line.Origin.X, line.Direction.X, box.X),
                        AxisRange(line.Origin.Y, line.Direction.Y, box.Y));
}

[[nodiscard]] inline bool IsBounded(const Box2& box)
{
    return box.X.IsBounded() && box.Y.IsBounded();
}

/// A box containing the image of box under map; Everything for an unbounded
/// box.
[[nodiscard]] inline Box2 ImageBox(const ParameterMap& map, const Box2& box)
{
    if (!IsBounded(box)) {
        return Box2::Everything();
    }
    return {Interval{map.Offset.X} + box.X * Interval{map.U.X} + box.Y * Interval{map.V.X},
            Interval{map.Offset.Y} + box.X * Interval{map.U.Y} + box.Y * Interval{map.V.Y}};
}

/// The map from B's parameters back to A's.
[[nodiscard]] inline ParameterMap Inverse(const ParameterMap& map)
{
    const double det = map.U.X * map.V.Y - map.U.Y * map.V.X;
    ParameterMap inverse;
    inverse.U = Vec2{map.V.Y, -map.U.Y} / det;
    inverse.V = Vec2{-map.V.X, map.U.X} / det;
    inverse.Offset = -(inverse.U * map.Offset.X + inverse.V * map.Offset.Y);
    return inverse;
}

/// A transverse branch: the surfaces cross, so the sides differ (ADR-0029
/// section 3). Uncertified until M4.
[[nodiscard]] inline IntersectionBranch Crossing(Curve curve, Curve2d onA, Curve2d onB,
                                                 const Interval& range)
{
    IntersectionBranch branch;
    branch.Curve3d = std::move(curve);
    branch.OnA = std::move(onA);
    branch.OnB = std::move(onB);
    branch.Range = range;
    branch.Contact = ContactType::Crossing;
    branch.Transition = SideTransition::OppositeSide;
    return branch;
}

/// A tangential branch between quadrics: a touch without a crossing, proved
/// from the closed form (ADR-0029 section 6), so the sides agree.
[[nodiscard]] inline IntersectionBranch Tangency(Curve curve, Curve2d onA, Curve2d onB,
                                                 const Interval& range)
{
    IntersectionBranch branch = Crossing(std::move(curve), std::move(onA), std::move(onB), range);
    branch.Contact = ContactType::Tangential;
    branch.Transition = SideTransition::SameSide;
    return branch;
}

/// The branch certified against both surfaces (ADR-0009 section 4): its
/// Bound is the larger of the two deviations, each proved within the
/// Context's resolution over the whole range. A branch over an unbounded
/// range stays uncertified: a face's loops bound it first (M6).
[[nodiscard]] inline Result<void> Certify(IntersectionBranch& branch, const Surface& a,
                                          const Surface& b, const Context& ctx)
{
    if (!branch.Range.IsBounded()) {
        return {};
    }
    const double target = ctx.GetResolution();
    STEREON_TRY(const DeviationCertificate onA,
                CertifyDeviation(branch.Curve3d, a, branch.OnA, branch.Range, target));
    STEREON_TRY(const DeviationCertificate onB,
                CertifyDeviation(branch.Curve3d, b, branch.OnB, branch.Range, target));
    branch.Bound = std::max(onA.Bound, onB.Bound);
    return {};
}

/// Certify then keep: the branch joins the result, or the failure is the
/// result's (ADR-0029: never a partial success).
[[nodiscard]] inline Result<void> Keep(SurfaceSurfaceIntersection& result,
                                       IntersectionBranch branch, const Surface& a,
                                       const Surface& b, const Context& ctx)
{
    STEREON_TRY_VOID(Certify(branch, a, b, ctx));
    result.Branches.push_back(std::move(branch));
    return {};
}

/// The same result seen from the other side: A and B exchanged.
[[nodiscard]] inline SurfaceSurfaceIntersection Swapped(SurfaceSurfaceIntersection result)
{
    for (IntersectionPoint& point : result.Points) {
        std::swap(point.OnA, point.OnB);
    }
    for (IntersectionBranch& branch : result.Branches) {
        std::swap(branch.OnA, branch.OnB);
    }
    for (OverlapRegion& overlap : result.Overlaps) {
        std::swap(overlap.OnA, overlap.OnB);
        overlap.AToB = Inverse(overlap.AToB);
    }
    return result;
}

} // namespace
} // namespace Stereon::Detail
