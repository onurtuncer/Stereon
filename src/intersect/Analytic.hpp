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
// fragment. Everything is in an anonymous namespace, so each unit gets its
// own copy.

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
