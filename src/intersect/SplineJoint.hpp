// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Where a B-spline's point is not a smooth function of the geometry: shared
// by the products of points on B-spline curves and surfaces
// (Linearization.cpp, PlanarLinearization.cpp; ADR-0021, section 9).
//
// Included in the purview of those implementation units, after their
// imports, with <algorithm>, <cmath> and <limits> in their global module
// fragment. Everything is in an anonymous namespace, so each unit gets its
// own copy.

#pragma once

namespace Stereon::Detail {
namespace {

/// One-sided first derivatives differing by more than this many units of
/// rounding, relative to the larger, are a corner: Split and Deriv each
/// round the tangents of a smooth joint by a few units.
constexpr double JointRounding = 64.0 * std::numeric_limits<double>::epsilon();

/// The multiplicity of the knot t in a knot vector: 0 when t is no knot.
[[nodiscard]] inline std::size_t MultiplicityOf(const std::vector<double>& knots, const double t)
{
    return static_cast<std::size_t>(std::ranges::count(knots, t));
}

/// Whether the parameter t of a B-spline curve is where a point of an
/// intersection cannot move smoothly with the geometry: one of the curve's
/// own ends, past which it has no point (a closed curve's too, where the
/// parameter jumps by the period), or a corner, an interior knot of
/// multiplicity at least the degree where the one-sided tangents differ,
/// so that the point changes arm as the geometry moves. A knot of lower
/// multiplicity, or a joint whose two tangents agree to rounding, is
/// smooth enough for a first derivative.
template <typename Spline>
[[nodiscard]] bool AtSplineJoint(const Spline& curve, const double t)
{
    const Interval domain = Domain(curve);
    if (!(t > domain.GetLower() && t < domain.GetUpper())) {
        return true;
    }
    if (MultiplicityOf(curve.Knots, t) < static_cast<std::size_t>(curve.Degree)) {
        return false;
    }
    const auto halves = Split(curve, t);
    if (!halves) {
        return true;
    }
    const auto left = Deriv(halves->first, t, 1);
    const auto right = Deriv(halves->second, t, 1);
    return Length(left - right) > JointRounding * std::max(Length(left), Length(right));
}

/// How a point of a B-spline surface at uv stands for a derivative: on the
/// surface's own boundary (the point leaves the surface as the geometry
/// moves, a closed surface's seam included), on an interior knot line of
/// multiplicity at least the degree in its direction (the partials across
/// it may jump, and whether they do is not established), or neither.
enum class SurfaceJoint : std::uint8_t { None, Boundary, KnotLine };

/// Where uv stands on a B-spline surface: its boundary, an interior knot line
/// of multiplicity at least the degree, or neither.
[[nodiscard]] inline SurfaceJoint SurfaceJointAt(const BSplineSurface& surface, const Point2& uv)
{
    const Box2 domain = Domain(surface);
    const auto inside = [](const Interval& range, const double x) {
        return x > range.GetLower() && x < range.GetUpper();
    };
    if (!inside(domain.X, uv.X) || !inside(domain.Y, uv.Y)) {
        return SurfaceJoint::Boundary;
    }
    if (MultiplicityOf(surface.KnotsU, uv.X) >= static_cast<std::size_t>(surface.DegreeU) ||
        MultiplicityOf(surface.KnotsV, uv.Y) >= static_cast<std::size_t>(surface.DegreeV)) {
        return SurfaceJoint::KnotLine;
    }
    return SurfaceJoint::None;
}

} // namespace
} // namespace Stereon::Detail
