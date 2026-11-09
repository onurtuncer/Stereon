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
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Bernstein.hpp"

/// Phase 2 generic SSI (ADR-0010), its first stages: proofs of separation.
/// This is NOT yet a NURBS intersection solver. A positive-weight rational
/// B-spline lies in the convex hull of its control net, so strictly
/// separated control-net coordinate intervals prove two surfaces disjoint
/// over any requested subdomains: the whole nets first, then the Bezier
/// patches of the requested domains, split until each pair is proved apart.
/// The patches come from Boehm's knot insertion and de Casteljau's algorithm
/// in interval arithmetic, so every proof holds for the surfaces as given.
/// A pair that cannot be separated proves nothing and must never return
/// success: a pair within the tolerance is NotImplemented until branch
/// discovery, marching and certification exist, and a search past its budget
/// is IntersectionFailed.

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

// --- Enclosed control nets -----------------------------------------------

/// The pairs of patches one search may visit. Separating surfaces a gap g
/// apart costs about (extent / g)^2 pairs, so this settles gaps down to
/// about 1% of the patches' extent; closer ones fail rather than run on.
constexpr std::size_t PairBudget = std::size_t{1} << 14;

/// A grid of enclosed homogeneous poles, i along u: index i * CountV + j.
struct Net {
    std::size_t CountU = 0;
    std::size_t CountV = 0;
    std::vector<HPoint> Poles;
};

[[nodiscard]] std::size_t LineCount(const Net& net, const bool alongU)
{
    return alongU ? net.CountV : net.CountU;
}

[[nodiscard]] std::size_t LineLength(const Net& net, const bool alongU)
{
    return alongU ? net.CountU : net.CountV;
}

/// Where the pole `at` of line `line` is held.
[[nodiscard]] std::size_t PoleIndex(const Net& net, const bool alongU, const std::size_t line,
                                    const std::size_t at)
{
    return alongU ? at * net.CountV + line : line * net.CountV + at;
}

[[nodiscard]] std::vector<HPoint> LineOf(const Net& net, const bool alongU, const std::size_t line)
{
    std::vector<HPoint> poles(LineLength(net, alongU));
    for (std::size_t at = 0; at < poles.size(); ++at) {
        poles[at] = net.Poles[PoleIndex(net, alongU, line, at)];
    }
    return poles;
}

void SetLine(Net& net, const bool alongU, const std::size_t line, const std::vector<HPoint>& poles)
{
    for (std::size_t at = 0; at < poles.size(); ++at) {
        net.Poles[PoleIndex(net, alongU, line, at)] = poles[at];
    }
}

[[nodiscard]] Net NetOf(const BSplineSurface& s)
{
    Net net{PoleCountU(s), PoleCountV(s), {}};
    net.Poles.reserve(s.Poles.size());
    for (std::size_t i = 0; i < s.Poles.size(); ++i) {
        net.Poles.push_back(
            Homogeneous(s.Poles[i], Interval{s.Weights.empty() ? 1.0 : s.Weights[i]}));
    }
    return net;
}

/// A net with its lines along one direction replaced by `lines`.
[[nodiscard]] Net WithLines(const Net& net, const bool alongU,
                            const std::vector<std::vector<HPoint>>& lines)
{
    const std::size_t length = lines.front().size();
    Net refined{alongU ? length : net.CountU, alongU ? net.CountV : length, {}};
    refined.Poles.resize(refined.CountU * refined.CountV);
    for (std::size_t line = 0; line < lines.size(); ++line) {
        SetLine(refined, alongU, line, lines[line]);
    }
    return refined;
}

/// Every distinct knot of the domain raised to multiplicity p along one
/// direction, on all its lines, as Bernstein.hpp's SegmentsOf does for a curve.
[[nodiscard]] Result<void> RaiseKnots(Net& net, std::vector<double>& knots, const std::size_t p,
                                      const bool alongU)
{
    std::vector<double> values(knots.begin() + static_cast<std::ptrdiff_t>(p),
                               knots.end() - static_cast<std::ptrdiff_t>(p));
    const auto [first, last] = std::ranges::unique(values);
    values.erase(first, last);
    std::vector<std::vector<HPoint>> lines;
    for (std::size_t line = 0; line < LineCount(net, alongU); ++line) {
        lines.push_back(LineOf(net, alongU, line));
    }
    for (const double value : values) {
        while (std::ranges::count(knots, value) < static_cast<std::ptrdiff_t>(p)) {
            std::vector<double> refined;
            for (std::vector<HPoint>& line : lines) {
                refined = knots;
                STEREON_TRY_VOID(InsertKnot(refined, line, p, value));
            }
            knots = std::move(refined);
        }
    }
    net = WithLines(net, alongU, lines);
    return {};
}

/// The first pole index, past the degree, of each non-empty span of the domain.
[[nodiscard]] std::vector<std::size_t> SpansOf(const std::vector<double>& knots,
                                               const std::size_t p, const std::size_t poles,
                                               const Interval& domain)
{
    std::vector<std::size_t> spans;
    for (std::size_t j = p; j < poles; ++j) {
        if (knots[j] < knots[j + 1] && domain.GetLower() <= knots[j] &&
            knots[j + 1] <= domain.GetUpper()) {
            spans.push_back(j);
        }
    }
    return spans;
}

// --- Bezier patches ------------------------------------------------------

/// A Bezier patch: a stretch of the surface's parameters, its enclosed net,
/// and a box of 3D space holding it. NextAlongU alternates the splits.
struct Patch {
    Box2 Parameters;
    Net Poles;
    Bounds Extent;
    bool NextAlongU = true;
};

/// The convex hull's box: each pole's Cartesian coordinates, enclosed. A
/// weight not proved positive leaves nothing to prove with.
[[nodiscard]] Bounds ExtentOf(const Net& net)
{
    Bounds box;
    for (const HPoint& h : net.Poles) {
        if (!(h.W.GetLower() > 0.0)) {
            return {Interval::Everything(), Interval::Everything(), Interval::Everything()};
        }
        box = {Hull(box[0], h.X / h.W), Hull(box[1], h.Y / h.W), Hull(box[2], h.Z / h.W)};
    }
    return box;
}

[[nodiscard]] double LargestWidth(const Bounds& box)
{
    double largest = 0.0;
    for (const Interval& side : box) {
        largest = std::max(largest, Width(side));
    }
    return largest;
}

[[nodiscard]] Patch MakePatch(const Box2& parameters, Net net)
{
    Bounds extent = ExtentOf(net);
    return {parameters, std::move(net), extent};
}

/// The patch of spans ending at poles (ju, jv) of a net of Bezier form.
[[nodiscard]] Net SubNet(const Net& net, const std::size_t p, const std::size_t q,
                         const std::size_t ju, const std::size_t jv)
{
    Net patch{p + 1, q + 1, {}};
    for (std::size_t i = ju - p; i <= ju; ++i) {
        for (std::size_t j = jv - q; j <= jv; ++j) {
            patch.Poles.push_back(net.Poles[i * net.CountV + j]);
        }
    }
    return patch;
}

/// The surface's Bezier patches.
[[nodiscard]] Result<std::vector<Patch>> PatchesOf(const BSplineSurface& s)
{
    const auto p = static_cast<std::size_t>(s.DegreeU);
    const auto q = static_cast<std::size_t>(s.DegreeV);
    std::vector<double> knotsU = s.KnotsU;
    std::vector<double> knotsV = s.KnotsV;
    Net net = NetOf(s);
    STEREON_TRY_VOID(RaiseKnots(net, knotsU, p, true));
    STEREON_TRY_VOID(RaiseKnots(net, knotsV, q, false));
    const Box2 domain = Domain(s);
    std::vector<Patch> patches;
    for (const std::size_t ju : SpansOf(knotsU, p, net.CountU, domain.X)) {
        for (const std::size_t jv : SpansOf(knotsV, q, net.CountV, domain.Y)) {
            const Box2 span{Interval{knotsU[ju], knotsU[ju + 1]},
                            Interval{knotsV[jv], knotsV[jv + 1]}};
            patches.push_back(MakePatch(span, SubNet(net, p, q, ju, jv)));
        }
    }
    return patches;
}

/// de Casteljau's algorithm on one line at t, the two halves' poles.
void SplitLine(const std::vector<HPoint>& line, const Interval& t, std::vector<HPoint>& left,
               std::vector<HPoint>& right)
{
    std::vector<HPoint> work = line;
    const std::size_t n = line.size();
    left.assign(n, HPoint{});
    right.assign(n, HPoint{});
    left[0] = work[0];
    right[n - 1] = work[n - 1];
    for (std::size_t r = 1; r < n; ++r) {
        for (std::size_t i = 0; i + r < n; ++i) {
            work[i] = Blend(work[i], work[i + 1], t);
        }
        left[r] = work[0];
        right[n - 1 - r] = work[n - 1 - r];
    }
}

/// The patch cut at a parameter strictly inside its stretch; nothing when
/// the stretch is too narrow for its fraction to be enclosed.
[[nodiscard]] std::optional<std::pair<Patch, Patch>> Split(const Patch& patch, const bool alongU,
                                                           const double at)
{
    const Interval stretch = alongU ? patch.Parameters.X : patch.Parameters.Y;
    const Interval width = Interval{stretch.GetUpper()} - Interval{stretch.GetLower()};
    if (!(stretch.GetLower() < at && at < stretch.GetUpper() && width.GetLower() > 0.0)) {
        return std::nullopt;
    }
    const Interval t = (Interval{at} - Interval{stretch.GetLower()}) / width;
    std::vector<std::vector<HPoint>> lower(LineCount(patch.Poles, alongU));
    std::vector<std::vector<HPoint>> upper(lower.size());
    for (std::size_t line = 0; line < lower.size(); ++line) {
        SplitLine(LineOf(patch.Poles, alongU, line), t, lower[line], upper[line]);
    }
    Box2 below = patch.Parameters;
    Box2 above = patch.Parameters;
    (alongU ? below.X : below.Y) = Interval{stretch.GetLower(), at};
    (alongU ? above.X : above.Y) = Interval{at, stretch.GetUpper()};
    Patch first = MakePatch(below, WithLines(patch.Poles, alongU, lower));
    Patch second = MakePatch(above, WithLines(patch.Poles, alongU, upper));
    first.NextAlongU = second.NextAlongU = !patch.NextAlongU;
    return std::pair{std::move(first), std::move(second)};
}

/// The patch cut down to `wanted` along one direction; nothing when they
/// are disjoint. A cut too narrow to enclose is skipped: a larger patch
/// only makes a proof harder, never wrong.
[[nodiscard]] std::optional<Patch> Clipped(Patch patch, const bool alongU, const Interval& wanted)
{
    const Interval stretch = alongU ? patch.Parameters.X : patch.Parameters.Y;
    const Interval kept = Intersection(stretch, wanted);
    if (kept.IsEmpty()) {
        return std::nullopt;
    }
    if (auto cut = Split(patch, alongU, kept.GetLower())) {
        patch = std::move(cut->second);
    }
    if (auto cut = Split(patch, alongU, kept.GetUpper())) {
        patch = std::move(cut->first);
    }
    return patch;
}

/// The surface's patches over a requested domain.
[[nodiscard]] Result<std::vector<Patch>> PatchesOver(const BSplineSurface& s, const Box2& domain)
{
    STEREON_TRY(const std::vector<Patch> patches, PatchesOf(s));
    std::vector<Patch> kept;
    for (const Patch& patch : patches) {
        if (const auto u = Clipped(patch, true, domain.X)) {
            if (auto uv = Clipped(*u, false, domain.Y)) {
                kept.push_back(std::move(*uv));
            }
        }
    }
    return kept;
}

// --- The pairs of patches ------------------------------------------------

struct PatchPair {
    Patch A;
    Patch B;
};

/// The patch's two halves along its next direction, or its other one when
/// that is too narrow; nothing when neither can be cut.
[[nodiscard]] std::optional<std::pair<Patch, Patch>> Halves(const Patch& patch)
{
    for (const bool alongU : {patch.NextAlongU, !patch.NextAlongU}) {
        if (auto cut =
                Split(patch, alongU, Midpoint(alongU ? patch.Parameters.X : patch.Parameters.Y))) {
            return cut;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::string Describe(const PatchPair& pair)
{
    const auto text = [](const Box2& b) {
        return "[" + std::to_string(b.X.GetLower()) + ", " + std::to_string(b.X.GetUpper()) +
               "] x [" + std::to_string(b.Y.GetLower()) + ", " + std::to_string(b.Y.GetUpper()) +
               "]";
    };
    return text(pair.A.Parameters) + " against " + text(pair.B.Parameters);
}

[[nodiscard]] auto NotSeparated(const PatchPair& pair)
{
    return Fail(ErrorCode::NotImplemented,
                "BSpline x BSpline: patches " + Describe(pair) +
                    " could not be proved separate; generic SSI branch discovery and "
                    "certification are not implemented (ADR-0010 Phase 2)");
}

/// The larger patch's halves, each against the other patch, pushed; false
/// when neither patch can be cut.
[[nodiscard]] bool PushSplit(std::vector<PatchPair>& stack, const PatchPair& pair)
{
    const bool cutA = LargestWidth(pair.A.Extent) >= LargestWidth(pair.B.Extent);
    for (const bool first : {cutA, !cutA}) {
        if (auto halves = Halves(first ? pair.A : pair.B)) {
            const Patch& other = first ? pair.B : pair.A;
            stack.push_back(first ? PatchPair{std::move(halves->first), other}
                                  : PatchPair{other, std::move(halves->first)});
            stack.push_back(first ? PatchPair{std::move(halves->second), other}
                                  : PatchPair{other, std::move(halves->second)});
            return true;
        }
    }
    return false;
}

/// Every pair proved apart, or the reason one was not. Depth first, so that
/// a pair that cannot be separated is reached in about the depth of the
/// subdivision rather than after a whole level of pairs.
[[nodiscard]] Result<void> ProveAllApart(std::vector<PatchPair> stack, const double tolerance)
{
    std::size_t visited = 0;
    while (!stack.empty()) {
        const PatchPair pair = std::move(stack.back());
        stack.pop_back();
        if (ProvedApart(pair.A.Extent, pair.B.Extent, tolerance)) {
            continue;
        }
        if (++visited > PairBudget) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSpline x BSpline: the search for a separation ran out of its budget at "
                        "patches " +
                            Describe(pair));
        }
        const bool withinTolerance =
            std::max(LargestWidth(pair.A.Extent), LargestWidth(pair.B.Extent)) <= tolerance;
        if (withinTolerance || !PushSplit(stack, pair)) {
            return NotSeparated(pair);
        }
    }
    return {};
}

} // namespace

using SSIResult = Result<SurfaceSurfaceIntersection>;

SSIResult NurbsHull(const Surface& a, const Box2& domainA, const Surface& b, const Box2& domainB,
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
    STEREON_TRY(std::vector<Patch> patchesA, PatchesOver(*sa, domainA));
    STEREON_TRY(const std::vector<Patch> patchesB, PatchesOver(*sb, domainB));
    std::vector<PatchPair> pairs;
    for (const Patch& pa : patchesA) {
        for (const Patch& pb : patchesB) {
            pairs.push_back({pa, pb});
        }
    }
    STEREON_TRY_VOID(ProveAllApart(std::move(pairs), tolerance));
    return SurfaceSurfaceIntersection{};
}

} // namespace Stereon::Detail
