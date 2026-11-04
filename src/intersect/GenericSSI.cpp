// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string_view>
#include <vector>

module Stereon.Intersect;

import :Dispatch;

/// Phase 2 generic SSI foundation (ADR-0010): a conservative broad phase.
/// This is NOT a NURBS intersection solver. The entire positive-weight
/// rational B-spline lies in the convex hull of its control net. Therefore
/// strictly separated control-net coordinate intervals prove the surfaces
/// disjoint, over *any* requested parameter subdomains. A hull overlap
/// proves nothing about actual intersection and must never return success.
///
/// The subsequent phases (Bézier patch decomposition, interval pruning,
/// complete branch discovery, marching and certification) are intentionally
/// not represented as success until implemented and qualified.

namespace Stereon::Detail {
namespace {

using Bounds = std::array<Interval, 3>;

/// Reject direct aggregate construction that violates a factory invariant.
/// Without those invariants the convex-hull proof does not apply.
[[nodiscard]] bool ValidKnotVector(const std::vector<double>& knots, const int degree)
{
    if (degree < 1 || knots.size() <= static_cast<std::size_t>(degree) + 1 ||
        !std::ranges::is_sorted(knots)) {
        return false;
    }
    for (const double knot : knots) {
        if (!std::isfinite(knot)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidPoleLayout(const BSplineSurface& s)
{
    const std::size_t nu = PoleCountU(s);
    const std::size_t nv = PoleCountV(s);
    return nv > 0 && nu > 0 && s.Poles.size() % nv == 0 && s.Poles.size() / nv == nu;
}

[[nodiscard]] bool ValidWeights(const BSplineSurface& s)
{
    if (!s.Weights.empty() && s.Weights.size() != s.Poles.size()) {
        return false;
    }
    for (const double weight : s.Weights) {
        if (!(std::isfinite(weight) && weight > 0.0)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidPoles(const BSplineSurface& s)
{
    for (const Point3& pole : s.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidNet(const BSplineSurface& s)
{
    const bool u = ValidKnotVector(s.KnotsU, s.DegreeU);
    const bool v = ValidKnotVector(s.KnotsV, s.DegreeV);
    if (!u || !v) {
        return false; // Do not compute pole counts from invalid knot vectors.
    }
    if (!ValidPoleLayout(s) || !ValidWeights(s) || !ValidPoles(s)) {
        return false;
    }
    // Apply the same degree, multiplicity, span and weight constraints
    // as the public B-spline factory before using the convex-hull theorem.
    const std::size_t nu = PoleCountU(s);
    const std::size_t nv = PoleCountV(s);
    constexpr std::string_view where = "NurbsHull";
    return CheckKnotVector(s.DegreeU, s.KnotsU, nu, where).has_value() &&
           CheckKnotVector(s.DegreeV, s.KnotsV, nv, where).has_value() &&
           CheckWeights(s.Weights, s.Poles.size(), where).has_value();
}

/// Bound each coordinate by the control points. Widen once by a floating
/// point successor on either end. The spline's nonnegative partition-of-
/// unity basis (and strictly positive rational weights) gives containment;
/// no sampled geometry, subtraction of close points or guessed tolerance
/// is involved in this proof.
[[nodiscard]] Bounds NetBounds(const BSplineSurface& s)
{
    const Point3& first = s.Poles.front();
    std::array<double, 3> low{first.X, first.Y, first.Z};
    std::array<double, 3> high = low;
    for (const Point3& pole : s.Poles) {
        const std::array<double, 3> xyz{pole.X, pole.Y, pole.Z};
        for (std::size_t i = 0; i < xyz.size(); ++i) {
            low[i] = std::min(low[i], xyz[i]);
            high[i] = std::max(high[i], xyz[i]);
        }
    }
    constexpr double inf = std::numeric_limits<double>::infinity();
    return {Interval{std::nextafter(low[0], -inf), std::nextafter(high[0], inf)},
            Interval{std::nextafter(low[1], -inf), std::nextafter(high[1], inf)},
            Interval{std::nextafter(low[2], -inf), std::nextafter(high[2], inf)}};
}

/// Use interval subtraction to round the gap outward; if its *lower* bound
/// exceeds the pair tolerance, neither nominal nor tolerantly coincident
/// surface points can meet. Touching or uncertain hulls go to the pending
/// generic solver, not to an empty successful result.
[[nodiscard]] bool ProvedApart(const Bounds& a, const Bounds& b, const double tolerance)
{
    for (std::size_t i = 0; i < a.size(); ++i) {
        if ((b[i] - a[i]).GetLower() > tolerance || (a[i] - b[i]).GetLower() > tolerance) {
            return true;
        }
    }
    return false;
}

} // namespace

using SSIResult = Result<SurfaceSurfaceIntersection>;

SSIResult NurbsHull(const Surface& a, const Box2&, const Surface& b, const Box2&,
                    const double tolerance, const Context&)
{
    const BSplineSurface* sa = a.Target<BSplineSurface>();
    const BSplineSurface* sb = b.Target<BSplineSurface>();
    if (sa == nullptr || sb == nullptr || !ValidNet(*sa) || !ValidNet(*sb)) {
        return Fail(ErrorCode::InvalidInput, "NurbsHull: expected valid B-spline surfaces");
    }
    if (ProvedApart(NetBounds(*sa), NetBounds(*sb), tolerance)) {
        return SurfaceSurfaceIntersection{}; // Proven empty, including trims.
    }
    return Fail(ErrorCode::NotImplemented,
                "BSpline x BSpline: control hulls not proved separate; "
                "generic SSI branch discovery and certification are not implemented "
                "(ADR-0010 Phase 2)");
}

} // namespace Stereon::Detail
