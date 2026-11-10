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
#include <memory>
#include <numbers>
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

/// The generic path of ADR-0010: two B-spline surfaces, rational or not.
///
/// A positive-weight B-spline lies in the convex hull of its control net, so
/// separated nets prove two surfaces disjoint over any domains. The nets
/// first, then the Bezier patches of the requested domains, enclosed by
/// Boehm's knot insertion and de Casteljau's algorithm in interval
/// arithmetic, so that every proof holds for the surfaces as given. Pairs of
/// patches are cut, off-centre, until each is settled:
///
/// - proved apart, the nets' boxes further than the tolerance;
/// - transverse: a direction L with L . (n_A x n_B) > 0 over the pair. Every
///   branch advances along L, so none closes up inside and each leaves through
///   the boundary. The roots of the pair's eight boundary curves against the
///   other patch, isolated by Krawczyk's operator, are therefore all its
///   branch ends: none is no branch, two are one arc, more need a smaller pair;
/// - an isolated touch, polynomial surfaces only: a unique solution of the
///   extremal system at which the relative second form is definite, within
///   the tolerance.
///
/// Arcs join at shared seeds into branches, every free end on a boundary of
/// the domains. Each arc is marched by Newton's method on S_A = S_B and
/// L . S_A = advance, halving chords until they hold, and the polylines in
/// space and on both surfaces are certified by CertifyDeviation, tightened
/// until they pass. Polylines, because a fitted cubic's pcurve can leave a
/// surface's domain.
///
/// What cannot be settled fails, never a partial or an unproved empty result:
/// a tangency along a curve, a branch point, an overlap, and any tangency of
/// rational surfaces, at the smallest pair, IntersectionFailed; a search past
/// its budget, the same.
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
/// apart costs about (extent / g)^2 pairs, so close ones fail rather than
/// run on.
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

struct PatchPair {
    Patch A;
    Patch B;
};

// --- Enclosures of a patch and its derivatives -----------------------------

/// Splits fall at this fraction of a stretch, not at its middle: symmetric
/// input would otherwise lay a branch along a split line, where a boundary
/// root is not isolated. sqrt(2) / 3.
constexpr double SplitFraction = std::numbers::sqrt2 / 3.0;

/// Newton's iterations for a point, and the depth of a chord's refinement.
constexpr int NewtonLimit = 16;
constexpr int RefineDepth = 26;
/// The narrowings of one Krawczyk search.
constexpr int KrawczykLimit = 48;
/// How often a branch's chords are tightened before certification fails.
constexpr int CertifyRetries = 5;

[[nodiscard]] Bounds Unbounded()
{
    return {Interval::Everything(), Interval::Everything(), Interval::Everything()};
}

[[nodiscard]] bool AllBounded(const Bounds& b)
{
    return b[0].IsBounded() && b[1].IsBounded() && b[2].IsBounded();
}

[[nodiscard]] Bounds CrossOf(const Bounds& a, const Bounds& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[nodiscard]] Interval DotOf(const Bounds& a, const Bounds& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[nodiscard]] Bounds MinusOf(const Bounds& a, const Bounds& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[nodiscard]] HPoint HullOf(const Net& net)
{
    HPoint hull;
    for (const HPoint& h : net.Poles) {
        hull = {Hull(hull.X, h.X), Hull(hull.Y, h.Y), Hull(hull.Z, h.Z), Hull(hull.W, h.W)};
    }
    return hull;
}

[[nodiscard]] Interval WidthOf(const Interval& stretch)
{
    return Interval{stretch.GetUpper()} - Interval{stretch.GetLower()};
}

/// The net of the partial derivative along one direction, its degree one
/// lower; a net of one zero pole across at degree 0. The stretch's width
/// must be proved positive.
[[nodiscard]] Net DerivativeNet(const Net& net, const bool alongU, const Interval& width)
{
    const std::size_t length = LineLength(net, alongU);
    if (length < 2) {
        Net zero{1, 1, {}};
        zero.Poles.push_back({Interval{0.0}, Interval{0.0}, Interval{0.0}, Interval{0.0}});
        return zero;
    }
    const Interval factor = Interval{static_cast<double>(length - 1)} / width;
    const auto minus = [&](const HPoint& high, const HPoint& low) {
        return HPoint{(high.X - low.X) * factor, (high.Y - low.Y) * factor,
                      (high.Z - low.Z) * factor, (high.W - low.W) * factor};
    };
    std::vector<std::vector<HPoint>> lines(LineCount(net, alongU));
    for (std::size_t line = 0; line < lines.size(); ++line) {
        const std::vector<HPoint> poles = LineOf(net, alongU, line);
        for (std::size_t at = 0; at + 1 < poles.size(); ++at) {
            lines[line].push_back(minus(poles[at + 1], poles[at]));
        }
    }
    return WithLines(net, alongU, lines);
}

/// A patch's hulls over all of it: its point and first partials.
struct Hulls {
    Bounds Point;
    Bounds Du;
    Bounds Dv;
};

[[nodiscard]] Bounds QuotientOf(const HPoint& d, const Bounds& point, const Interval& weight)
{
    return {(d.X - point[0] * d.W) / weight, (d.Y - point[1] * d.W) / weight,
            (d.Z - point[2] * d.W) / weight};
}

[[nodiscard]] Hulls HullsOf(const Patch& patch)
{
    const HPoint h = HullOf(patch.Poles);
    const Interval widthU = WidthOf(patch.Parameters.X);
    const Interval widthV = WidthOf(patch.Parameters.Y);
    if (!(h.W.GetLower() > 0.0 && widthU.GetLower() > 0.0 && widthV.GetLower() > 0.0)) {
        return {Unbounded(), Unbounded(), Unbounded()};
    }
    const HPoint du = HullOf(DerivativeNet(patch.Poles, true, widthU));
    const HPoint dv = HullOf(DerivativeNet(patch.Poles, false, widthV));
    return {patch.Extent, QuotientOf(du, patch.Extent, h.W), QuotientOf(dv, patch.Extent, h.W)};
}

/// The fraction of a stretch at which x lies, enclosed.
[[nodiscard]] Interval FractionOf(const Interval& stretch, const Interval& x)
{
    return (x - Interval{stretch.GetLower()}) / WidthOf(stretch);
}

[[nodiscard]] HPoint Reduced(std::vector<HPoint> line, const Interval& t)
{
    for (std::size_t r = line.size() - 1; r > 0; --r) {
        for (std::size_t i = 0; i < r; ++i) {
            line[i] = Blend(line[i], line[i + 1], t);
        }
    }
    return line.front();
}

/// The net's homogeneous value at enclosed parameters of its stretch.
[[nodiscard]] HPoint ValueOf(const Net& net, const Box2& stretch, const Interval& u,
                             const Interval& v)
{
    std::vector<HPoint> across;
    for (std::size_t i = 0; i < net.CountU; ++i) {
        across.push_back(Reduced(LineOf(net, false, i), FractionOf(stretch.Y, v)));
    }
    return Reduced(std::move(across), FractionOf(stretch.X, u));
}

[[nodiscard]] Bounds CartesianOf(const HPoint& h)
{
    return {h.X / h.W, h.Y / h.W, h.Z / h.W};
}

// --- Krawczyk's operator ---------------------------------------------------

template <std::size_t N>
using IntervalsN = std::array<Interval, N>;
template <std::size_t N>
using IntervalMatrix = std::array<IntervalsN<N>, N>;
template <std::size_t N>
using DoubleMatrix = std::array<std::array<double, N>, N>;

enum class Isolation : std::uint8_t { None, Unique, Unsettled };

template <std::size_t N>
struct Isolated {
    Isolation Outcome = Isolation::Unsettled;
    IntervalsN<N> Box{};
};

/// The row with the largest entry of a column, at or below it.
template <std::size_t N>
[[nodiscard]] std::size_t PivotRow(const DoubleMatrix<N>& a, const std::size_t column)
{
    std::size_t best = column;
    for (std::size_t r = column + 1; r < N; ++r) {
        if (std::abs(a[r][column]) > std::abs(a[best][column])) {
            best = r;
        }
    }
    return best;
}

/// Gauss-Jordan inversion in doubles: only a preconditioner, or a Newton step
/// checked afterwards, so its rounding proves nothing and harms nothing.
template <std::size_t N>
[[nodiscard]] std::optional<DoubleMatrix<N>> Inverted(DoubleMatrix<N> a)
{
    DoubleMatrix<N> inverse{};
    for (std::size_t i = 0; i < N; ++i) {
        inverse[i][i] = 1.0;
    }
    for (std::size_t c = 0; c < N; ++c) {
        const std::size_t pivot = PivotRow<N>(a, c);
        if (!(std::abs(a[pivot][c]) > 0.0) || !std::isfinite(a[pivot][c])) {
            return std::nullopt;
        }
        std::swap(a[c], a[pivot]);
        std::swap(inverse[c], inverse[pivot]);
        const double scale = 1.0 / a[c][c];
        for (std::size_t k = 0; k < N; ++k) {
            a[c][k] *= scale;
            inverse[c][k] *= scale;
        }
        for (std::size_t r = 0; r < N; ++r) {
            const double factor = r == c ? 0.0 : a[r][c];
            for (std::size_t k = 0; k < N; ++k) {
                a[r][k] -= factor * a[c][k];
                inverse[r][k] -= factor * inverse[c][k];
            }
        }
    }
    return inverse;
}

template <std::size_t N>
[[nodiscard]] bool BoundedMatrix(const IntervalMatrix<N>& j)
{
    return std::ranges::all_of(j, [](const IntervalsN<N>& row) {
        return std::ranges::all_of(row, [](const Interval& e) { return e.IsBounded(); });
    });
}

/// I - Y J, enclosed, and whether its largest row sum is below one, which
/// proves that the system has at most one root in the box.
template <std::size_t N>
struct Contraction {
    IntervalMatrix<N> Matrix;
    bool Contracts = true;
};

template <std::size_t N>
[[nodiscard]] Contraction<N> ContractionOf(const DoubleMatrix<N>& y, const IntervalMatrix<N>& j)
{
    Contraction<N> c{};
    for (std::size_t r = 0; r < N; ++r) {
        double rowSum = 0.0;
        for (std::size_t k = 0; k < N; ++k) {
            Interval e{r == k ? 1.0 : 0.0};
            for (std::size_t s = 0; s < N; ++s) {
                e = e - Interval{y[r][s]} * j[s][k];
            }
            c.Matrix[r][k] = e;
            rowSum += std::max(std::abs(e.GetLower()), std::abs(e.GetUpper()));
        }
        c.Contracts = c.Contracts && rowSum < 1.0;
    }
    return c;
}

template <std::size_t N>
[[nodiscard]] std::array<double, N> MidpointsOf(const IntervalsN<N>& x)
{
    std::array<double, N> m{};
    for (std::size_t i = 0; i < N; ++i) {
        m[i] = Midpoint(x[i]);
    }
    return m;
}

/// One application of the operator, m - Y F(m) + (I - Y J)(X - m), met with X.
template <std::size_t N>
[[nodiscard]] std::optional<IntervalsN<N>>
Narrowed(const IntervalsN<N>& x, const std::array<double, N>& m, const IntervalsN<N>& atMiddle,
         const DoubleMatrix<N>& y, const IntervalMatrix<N>& rest)
{
    IntervalsN<N> next{};
    for (std::size_t r = 0; r < N; ++r) {
        Interval k{m[r]};
        for (std::size_t s = 0; s < N; ++s) {
            k = k - Interval{y[r][s]} * atMiddle[s] + rest[r][s] * (x[s] - Interval{m[s]});
        }
        next[r] = Intersection(k, x[r]);
        if (next[r].IsEmpty()) {
            return std::nullopt;
        }
    }
    return next;
}

template <std::size_t N>
[[nodiscard]] double TotalWidth(const IntervalsN<N>& x)
{
    double total = 0.0;
    for (const Interval& e : x) {
        total += Width(e);
    }
    return total;
}

/// Whether a narrowed box is small enough for its middle to be the root.
template <std::size_t N>
[[nodiscard]] bool Tight(const IntervalsN<N>& x)
{
    double magnitude = 1.0;
    for (const Interval& e : x) {
        magnitude = std::max({magnitude, std::abs(e.GetLower()), std::abs(e.GetUpper())});
    }
    return TotalWidth<N>(x) <= 0x1p-34 * magnitude;
}

/// The roots of a system in the box, F enclosed at points and J over the
/// box. None: the operator left the box, so there is no root. Unique: J's
/// contraction proves at most one, and the box has narrowed around it; a
/// root within rounding of the box's boundary counts as in it. Unsettled:
/// the box is too large for either, or too large to narrow to the root.
template <std::size_t N, typename F>
[[nodiscard]] Isolated<N> Krawczyk(IntervalsN<N> x, const IntervalMatrix<N>& j, const F& f)
{
    if (!BoundedMatrix<N>(j)) {
        return {};
    }
    DoubleMatrix<N> middle{};
    for (std::size_t r = 0; r < N; ++r) {
        middle[r] = MidpointsOf<N>(j[r]);
    }
    const std::optional<DoubleMatrix<N>> y = Inverted<N>(middle);
    if (!y) {
        return {};
    }
    const Contraction<N> rest = ContractionOf<N>(*y, j);
    for (int k = 0; k < KrawczykLimit; ++k) {
        const std::array<double, N> m = MidpointsOf<N>(x);
        const std::optional<IntervalsN<N>> next = Narrowed<N>(x, m, f(m), *y, rest.Matrix);
        if (!next) {
            return {Isolation::None, x};
        }
        if (!rest.Contracts) {
            return {};
        }
        const bool stalled = !(TotalWidth<N>(*next) < 0.75 * TotalWidth<N>(x));
        x = *next;
        if (stalled) {
            break;
        }
    }
    if (!Tight<N>(x)) {
        return {};
    }
    return {Isolation::Unique, x};
}

// --- Normal cones ----------------------------------------------------------

/// A direction along which every intersection curve of the pair advances:
/// L . (n_A x n_B) > 0 over the whole pair, n the cross of each patch's
/// partials. The curves are then transverse, and each is monotone along L,
/// so none closes up inside the pair (Sederberg and Meyers' loop test, on
/// the cones of normals). Nothing when no such direction is proved.
[[nodiscard]] std::optional<Vec3> AdvanceDirection(const Hulls& a, const Hulls& b)
{
    const Bounds tangent = CrossOf(CrossOf(a.Du, a.Dv), CrossOf(b.Du, b.Dv));
    if (!AllBounded(tangent)) {
        return std::nullopt;
    }
    const Vec3 middle{Midpoint(tangent[0]), Midpoint(tangent[1]), Midpoint(tangent[2])};
    const double length = Length(middle);
    if (!(length > 0.0)) {
        return std::nullopt;
    }
    const Vec3 l = middle / length;
    const Bounds direction{Interval{l.X}, Interval{l.Y}, Interval{l.Z}};
    if (!(DotOf(direction, tangent).GetLower() > 0.0)) {
        return std::nullopt;
    }
    return l;
}

// --- Seeds on the patches' boundaries --------------------------------------

/// A point of both surfaces, with its parameters on each.
struct Seed {
    Point3 Position;
    Point2 OnA;
    Point2 OnB;
};

[[nodiscard]] Interval StretchOf(const Patch& patch, const bool alongU)
{
    return alongU ? patch.Parameters.X : patch.Parameters.Y;
}

/// One edge of a patch, its parameter running along a direction, the other
/// held at `Fixed`; its partial's hull is the enclosure along it.
struct Edge {
    const Patch& Carrier;
    bool AlongU;
    double Fixed;
};

[[nodiscard]] Bounds PointOnEdge(const Edge& e, const double t)
{
    const Interval along{t};
    const Interval fixed{e.Fixed};
    return CartesianOf(ValueOf(e.Carrier.Poles, e.Carrier.Parameters, e.AlongU ? along : fixed,
                               e.AlongU ? fixed : along));
}

[[nodiscard]] Point2 ParametersOnEdge(const Edge& e, const double t)
{
    return e.AlongU ? Point2{t, e.Fixed} : Point2{e.Fixed, t};
}

[[nodiscard]] Point3 Middle(const Bounds& b)
{
    return {Midpoint(b[0]), Midpoint(b[1]), Midpoint(b[2])};
}

/// A stretch reaching a few roundings past its ends. A root on a line that
/// two cells share would otherwise be refused by both, each by a rounding.
[[nodiscard]] Interval Reaching(const Interval& stretch)
{
    const double magnitude =
        std::max({1.0, std::abs(stretch.GetLower()), std::abs(stretch.GetUpper())});
    return Inflated(stretch, 0x1p-44 * magnitude);
}

/// A derivative's enclosure over a patch, widened to hold on the stretch
/// that reaches past it: the extra reach is some 2^-44 of the parameter, so
/// the derivative moves by far less than this.
[[nodiscard]] Interval Slackened(const Interval& d)
{
    if (!d.IsBounded()) {
        return d;
    }
    return Inflated(d, 0x1p-30 * std::max(std::abs(d.GetLower()), std::abs(d.GetUpper())));
}
/// Where the edge meets the other patch: C(t) - S(u, v) = 0 in (t, u, v).
[[nodiscard]] Isolated<3> EdgeRoot(const Edge& e, const Bounds& along, const Patch& other,
                                   const Hulls& hulls)
{
    IntervalMatrix<3> j{};
    for (std::size_t r = 0; r < 3; ++r) {
        j[r] = {Slackened(along[r]), Slackened(-hulls.Du[r]), Slackened(-hulls.Dv[r])};
    }
    const IntervalsN<3> box{Reaching(StretchOf(e.Carrier, e.AlongU)), Reaching(other.Parameters.X),
                            Reaching(other.Parameters.Y)};
    return Krawczyk<3>(box, j, [&](const std::array<double, 3>& m) {
        return MinusOf(PointOnEdge(e, m[0]), CartesianOf(ValueOf(other.Poles, other.Parameters,
                                                                 Interval{m[1]}, Interval{m[2]})));
    });
}

/// The seeds of the pair, without repeats within `join`; nothing when a root
/// could not be settled.
struct PairView {
    const PatchPair& Pair;
    const Hulls& OnA;
    const Hulls& OnB;
    double Join;
};

void AddSeed(std::vector<Seed>& seeds, const Seed& seed, const double join)
{
    const bool known = std::ranges::any_of(
        seeds, [&](const Seed& s) { return Distance(s.Position, seed.Position) <= join; });
    if (!known) {
        seeds.push_back(seed);
    }
}

[[nodiscard]] Point2 ClampedInto(const Point2& p, const Box2& box)
{
    return {std::clamp(p.X, box.X.GetLower(), box.X.GetUpper()),
            std::clamp(p.Y, box.Y.GetLower(), box.Y.GetUpper())};
}

/// The seed an isolated root stands for, its parameters kept in the cells.
[[nodiscard]] Seed SeedOfRoot(const Edge& edge, const Isolated<3>& root, const Patch& other,
                              const bool carrierIsA)
{
    const Interval stretch = StretchOf(edge.Carrier, edge.AlongU);
    const double along = std::clamp(Midpoint(root.Box[0]), stretch.GetLower(), stretch.GetUpper());
    const Point2 onCarrier = ParametersOnEdge(edge, along);
    const Point2 onOther =
        ClampedInto(Point2{Midpoint(root.Box[1]), Midpoint(root.Box[2])}, other.Parameters);
    return {Middle(PointOnEdge(edge, along)), carrierIsA ? onCarrier : onOther,
            carrierIsA ? onOther : onCarrier};
}

[[nodiscard]] bool EdgesOfCarrier(const PairView& v, const bool carrierIsA,
                                  std::vector<Seed>& seeds)
{
    const Patch& carrier = carrierIsA ? v.Pair.A : v.Pair.B;
    const Patch& other = carrierIsA ? v.Pair.B : v.Pair.A;
    const Hulls& own = carrierIsA ? v.OnA : v.OnB;
    const Hulls& theirs = carrierIsA ? v.OnB : v.OnA;
    for (const bool alongU : {true, false}) {
        const Interval across = StretchOf(carrier, !alongU);
        for (const double fixed : {across.GetLower(), across.GetUpper()}) {
            const Edge edge{carrier, alongU, fixed};
            const Isolated<3> root = EdgeRoot(edge, alongU ? own.Du : own.Dv, other, theirs);
            if (root.Outcome == Isolation::Unsettled) {
                return false;
            }
            if (root.Outcome == Isolation::Unique) {
                AddSeed(seeds, SeedOfRoot(edge, root, other, carrierIsA), v.Join);
            }
        }
    }
    return true;
}

[[nodiscard]] std::optional<std::vector<Seed>> SeedsOf(const PairView& v)
{
    std::vector<Seed> seeds;
    if (!EdgesOfCarrier(v, true, seeds) || !EdgesOfCarrier(v, false, seeds)) {
        return std::nullopt;
    }
    return seeds;
}

// --- Contact that is not transverse ----------------------------------------

/// A patch's derivative nets, for the distance's gradient and Hessian.
struct SecondNets {
    Net Du;
    Net Dv;
    Net Duu;
    Net Duv;
    Net Dvv;
};

[[nodiscard]] SecondNets SecondNetsOf(const Patch& p)
{
    const Interval wu = WidthOf(p.Parameters.X);
    const Interval wv = WidthOf(p.Parameters.Y);
    Net du = DerivativeNet(p.Poles, true, wu);
    Net dv = DerivativeNet(p.Poles, false, wv);
    Net duu = DerivativeNet(du, true, wu);
    Net duv = DerivativeNet(du, false, wv);
    Net dvv = DerivativeNet(dv, false, wv);
    return {std::move(du), std::move(dv), std::move(duu), std::move(duv), std::move(dvv)};
}

[[nodiscard]] Bounds XyzHull(const Net& net)
{
    const HPoint h = HullOf(net);
    return {h.X, h.Y, h.Z};
}

/// The stretch a derivative net is evaluated over: a net of one pole has
/// the same value everywhere, so any stretch serves.
[[nodiscard]] Bounds XyzAt(const Net& net, const Box2& stretch, const double u, const double v)
{
    if (net.Poles.size() == 1) {
        return {net.Poles[0].X, net.Poles[0].Y, net.Poles[0].Z};
    }
    const HPoint h = ValueOf(net, stretch, Interval{u}, Interval{v});
    return {h.X, h.Y, h.Z};
}

/// S x T for enclosed vectors, as a sum of the two product terms the
/// derivative of a normal has.
[[nodiscard]] Bounds PlusOf(const Bounds& a, const Bounds& b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/// B's normal S_u x S_v and its partials, enclosed.
struct NormalJet {
    Bounds N;
    Bounds Nu;
    Bounds Nv;
};

/// The Jacobian of the extremal system over the pair (see ExtremalPoint), in
/// (uA, vA, uB, vB, lambda), lambda enclosed in `reach`.
[[nodiscard]] IntervalMatrix<5> ExtremalJacobian(const SecondNets& a, const SecondNets& b,
                                                 const NormalJet& n, const Interval& reach)
{
    const Bounds au = XyzHull(a.Du);
    const Bounds av = XyzHull(a.Dv);
    const Bounds bu = XyzHull(b.Du);
    const Bounds bv = XyzHull(b.Dv);
    IntervalMatrix<5> j{};
    for (std::size_t r = 0; r < 3; ++r) {
        j[r] = {au[r], av[r], reach * (-n.Nu[r]) - bu[r], reach * (-n.Nv[r]) - bv[r], -n.N[r]};
    }
    const Interval zero{0.0};
    j[3] = {DotOf(XyzHull(a.Duu), n.N), DotOf(XyzHull(a.Duv), n.N), DotOf(au, n.Nu),
            DotOf(au, n.Nv), zero};
    j[4] = {DotOf(XyzHull(a.Duv), n.N), DotOf(XyzHull(a.Dvv), n.N), DotOf(av, n.Nu),
            DotOf(av, n.Nv), zero};
    return j;
}

[[nodiscard]] Bounds XyzValue(const Patch& p, const double u, const double v)
{
    const HPoint h = ValueOf(p.Poles, p.Parameters, Interval{u}, Interval{v});
    return {h.X, h.Y, h.Z};
}

[[nodiscard]] NormalJet NormalHull(const SecondNets& b)
{
    const Bounds bu = XyzHull(b.Du);
    const Bounds bv = XyzHull(b.Dv);
    return {CrossOf(bu, bv), PlusOf(CrossOf(XyzHull(b.Duu), bv), CrossOf(bu, XyzHull(b.Duv))),
            PlusOf(CrossOf(XyzHull(b.Duv), bv), CrossOf(bu, XyzHull(b.Dvv)))};
}

/// Where lambda may lie over the pair: R . n / n . n, enclosed.
[[nodiscard]] std::optional<Interval> LambdaReach(const PatchPair& pair, const Bounds& n)
{
    const Bounds r = MinusOf(XyzHull(pair.A.Poles), XyzHull(pair.B.Poles));
    const Interval square = DotOf(n, n);
    if (!(square.GetLower() > 0.0)) {
        return std::nullopt;
    }
    const Interval reach = DotOf(r, n) / square;
    return Inflated(reach, 0.5 * Width(reach));
}

/// The points where the two surfaces share a tangent plane, solved for
/// (uA, vA, uB, vB, lambda):
///     S_A - S_B - lambda n_B = 0,  A_u . n_B = 0,  A_v . n_B = 0.
/// Every tangential contact is one (with lambda 0), and at an isolated touch
/// the system is regular, as the distance's own gradient is not: its Hessian
/// is singular wherever the surfaces meet.
[[nodiscard]] Isolated<5> ExtremalPoint(const PatchPair& pair, const SecondNets& a,
                                        const SecondNets& b)
{
    const NormalJet n = NormalHull(b);
    const std::optional<Interval> reach = LambdaReach(pair, n.N);
    if (!reach) {
        return {};
    }
    const IntervalsN<5> box{pair.A.Parameters.X, pair.A.Parameters.Y, pair.B.Parameters.X,
                            pair.B.Parameters.Y, *reach};
    return Krawczyk<5>(box, ExtremalJacobian(a, b, n, *reach), [&](const std::array<double, 5>& m) {
        const Bounds normal = CrossOf(XyzAt(b.Du, pair.B.Parameters, m[2], m[3]),
                                      XyzAt(b.Dv, pair.B.Parameters, m[2], m[3]));
        const Bounds r = MinusOf(XyzValue(pair.A, m[0], m[1]), XyzValue(pair.B, m[2], m[3]));
        const Interval l{m[4]};
        return IntervalsN<5>{r[0] - l * normal[0], r[1] - l * normal[1], r[2] - l * normal[2],
                             DotOf(XyzAt(a.Du, pair.A.Parameters, m[0], m[1]), normal),
                             DotOf(XyzAt(a.Dv, pair.A.Parameters, m[0], m[1]), normal)};
    });
}

/// A symmetric 2 x 2 form: its entries (uu, uv, vv).
using Form2 = std::array<double, 3>;

[[nodiscard]] Vec3 MiddleVec(const Bounds& b)
{
    return {Midpoint(b[0]), Midpoint(b[1]), Midpoint(b[2])};
}

/// A surface's first derivatives and its second fundamental form at a point.
struct LocalForm {
    Vec3 Du;
    Vec3 Dv;
    Form2 Second;
};

[[nodiscard]] LocalForm LocalAt(const SecondNets& s, const Box2& box, const double u,
                                const double v, const Vec3& n)
{
    const auto at = [&](const Net& net) { return MiddleVec(XyzAt(net, box, u, v)); };
    return {at(s.Du), at(s.Dv), {Dot(n, at(s.Duu)), Dot(n, at(s.Duv)), Dot(n, at(s.Dvv))}};
}

/// Whether the surfaces separate on every side of the point of common
/// tangency: the relative second form, A's less B's pulled back through the
/// common tangent plane, is definite. Evaluated at the solution: Krawczyk's
/// regularity proves the form nonsingular there, so rounding cannot flip it
/// unless it is almost singular, which the margin refuses.
[[nodiscard]] bool RelativelyDefinite(const LocalForm& a, const LocalForm& b)
{
    const double e = Dot(b.Du, b.Du);
    const double f = Dot(b.Du, b.Dv);
    const double g = Dot(b.Dv, b.Dv);
    const double det = e * g - f * f;
    if (!(det > 0.0)) {
        return false;
    }
    const auto toB = [&](const Vec3& t) {
        const double x = Dot(t, b.Du);
        const double y = Dot(t, b.Dv);
        return std::array{(g * x - f * y) / det, (e * y - f * x) / det};
    };
    const std::array<double, 2> p = toB(a.Du);
    const std::array<double, 2> q = toB(a.Dv);
    const auto pulled = [&](const std::array<double, 2>& s, const std::array<double, 2>& t) {
        return b.Second[0] * s[0] * t[0] + b.Second[1] * (s[0] * t[1] + s[1] * t[0]) +
               b.Second[2] * s[1] * t[1];
    };
    const double m11 = a.Second[0] - pulled(p, p);
    const double m12 = a.Second[1] - pulled(p, q);
    const double m22 = a.Second[2] - pulled(q, q);
    const double scale = m11 * m11 + 2.0 * m12 * m12 + m22 * m22;
    return m11 * m22 - m12 * m12 > 0x1p-20 * scale;
}

// --- The search ------------------------------------------------------------

/// A stretch of one branch inside a pair: from one seed to the other,
/// advancing along Direction, within the patches' parameter boxes.
struct Arc {
    Seed From;
    Seed To;
    Vec3 Direction;
    Box2 BoxA;
    Box2 BoxB;
};

struct Findings {
    std::vector<Arc> Arcs;
    std::vector<IntersectionPoint> Touches;
};

struct SearchSettings {
    double Tolerance;
    double Join;
    bool Polynomial;
};

[[nodiscard]] std::optional<std::pair<Patch, Patch>> Halves(const Patch& patch)
{
    for (const bool alongU : {patch.NextAlongU, !patch.NextAlongU}) {
        const Interval s = StretchOf(patch, alongU);
        const double at = s.GetLower() + SplitFraction * (s.GetUpper() - s.GetLower());
        if (auto cut = Split(patch, alongU, at)) {
            return cut;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::string Describe(const Box2& b)
{
    return "[" + Number(b.X.GetLower()) + ", " + Number(b.X.GetUpper()) + "] x [" +
           Number(b.Y.GetLower()) + ", " + Number(b.Y.GetUpper()) + "]";
}

[[nodiscard]] std::string Describe(const PatchPair& pair)
{
    return Describe(pair.A.Parameters) + " against " + Describe(pair.B.Parameters);
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

/// What became of a pair: nothing more to do, or it must be cut.
enum class Settled : std::uint8_t { Done, Cut };

/// A pair with a proved advance direction: its branches enter and leave
/// through the boundary, so no seed means none, and one, a pair whose cells
/// meet only in that point, holds no arc longer than the join; two are one
/// arc; more need a smaller pair.
[[nodiscard]] Settled SettleTransverse(const PatchPair& pair, const Hulls& a, const Hulls& b,
                                       const Vec3& direction, const SearchSettings& s,
                                       Findings& found)
{
    const std::optional<std::vector<Seed>> seeds = SeedsOf({pair, a, b, s.Join});
    if (!seeds || seeds->size() > 2) {
        return Settled::Cut;
    }
    if (seeds->size() == 2) {
        found.Arcs.push_back(
            {(*seeds)[0], (*seeds)[1], direction, pair.A.Parameters, pair.B.Parameters});
    }
    return Settled::Done;
}

/// The largest distance a pair of enclosed points can be apart, and the
/// smallest.
[[nodiscard]] Interval DistanceBetween(const Bounds& r)
{
    const Interval squares = DotOf(r, r);
    const double upper = std::sqrt(squares.GetUpper());
    const double lower = std::sqrt(std::max(0.0, squares.GetLower()));
    return {std::nextafter(lower, 0.0), std::nextafter(upper, Infinity)};
}

/// A pair with no advance direction, polynomial. The one thing settled here
/// is an isolated touch: a unique point of common tangency, the surfaces
/// separating on every side of it, within the tolerance of each other. No
/// solution rules out only tangency, not a crossing, and a saddle has
/// branches through it; both are cut further, and a saddle ends in failure
/// at the smallest pair.
[[nodiscard]] Settled SettleTouch(const PatchPair& pair, const SearchSettings& s,
                                  [[maybe_unused]] Findings& found)
{
    if (!s.Polynomial) {
        return Settled::Cut;
    }
    const SecondNets a = SecondNetsOf(pair.A);
    const SecondNets b = SecondNetsOf(pair.B);
    const Isolated<5> extremal = ExtremalPoint(pair, a, b);
    if (extremal.Outcome != Isolation::Unique) {
        return Settled::Cut;
    }
    const std::array<double, 5> m = MidpointsOf<5>(extremal.Box);
    const Bounds onA = XyzValue(pair.A, m[0], m[1]);
    const Vec3 n = Normalized(MiddleVec(CrossOf(XyzAt(b.Du, pair.B.Parameters, m[2], m[3]),
                                                XyzAt(b.Dv, pair.B.Parameters, m[2], m[3]))));
    const Interval gap = DistanceBetween(MinusOf(onA, XyzValue(pair.B, m[2], m[3])));
    if (gap.GetUpper() > s.Tolerance ||
        !RelativelyDefinite(LocalAt(a, pair.A.Parameters, m[0], m[1], n),
                            LocalAt(b, pair.B.Parameters, m[2], m[3], n))) {
        return Settled::Cut;
    }
    // A nearby extremum does not exclude a surrounding crossing loop.
    // Until isolation of the zero set is proved, fail closed rather than
    // discard the remainder of this patch pair.
    return Settled::Cut;
}

[[nodiscard]] Settled Settle(const PatchPair& pair, const SearchSettings& s, Findings& found)
{
    const Hulls a = HullsOf(pair.A);
    const Hulls b = HullsOf(pair.B);
    if (const std::optional<Vec3> direction = AdvanceDirection(a, b)) {
        return SettleTransverse(pair, a, b, *direction, s, found);
    }
    return SettleTouch(pair, s, found);
}

/// The whole search, depth first, so that a pair that cannot be settled is
/// reached in about the depth of the subdivision.
[[nodiscard]] Result<Findings> SearchPairs(std::vector<PatchPair> stack, const SearchSettings& s)
{
    Findings found;
    std::size_t visited = 0;
    while (!stack.empty()) {
        const PatchPair pair = std::move(stack.back());
        stack.pop_back();
        if (ProvedApart(pair.A.Extent, pair.B.Extent, s.Tolerance)) {
            continue;
        }
        if (++visited > PairBudget) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSpline x BSpline: the search ran out of its budget at patches " +
                            Describe(pair));
        }
        if (Settle(pair, s, found) == Settled::Done) {
            continue;
        }
        const bool smallest =
            std::max(LargestWidth(pair.A.Extent), LargestWidth(pair.B.Extent)) <= s.Tolerance;
        if (smallest || !PushSplit(stack, pair)) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSpline x BSpline: a tangential contact, a branch point or an overlap "
                        "could not be resolved at patches " +
                            Describe(pair));
        }
    }
    return found;
}

// --- Marching ---------------------------------------------------------------

/// A point of a branch: its position and its parameters on each surface,
/// with how far it has advanced along the arc's direction.
struct Sample {
    Point3 Position;
    Point2 OnA;
    Point2 OnB;
    double Advance;
};

struct Marcher {
    const Surface& A;
    const Surface& B;
    double Chord;
};

/// Newton's system: S_A - S_B, and how far the point is from the plane
/// L . S_A = advance. The columns are A's and B's partials.
struct NewtonSystem {
    std::array<double, 4> F;
    DoubleMatrix<4> J;
};

[[nodiscard]] NewtonSystem SystemAt(const Marcher& m, const Arc& arc, const Sample& s)
{
    const Vec3 r = Eval(m.A, s.OnA) - Eval(m.B, s.OnB);
    const Vec3 au = Deriv(m.A, s.OnA, 1, 0);
    const Vec3 av = Deriv(m.A, s.OnA, 0, 1);
    const Vec3 bu = Deriv(m.B, s.OnB, 1, 0);
    const Vec3 bv = Deriv(m.B, s.OnB, 0, 1);
    const Vec3 l = arc.Direction;
    return {{r.X, r.Y, r.Z, Dot(l, ToVec(Eval(m.A, s.OnA))) - s.Advance},
            {{{au.X, av.X, -bu.X, -bv.X},
              {au.Y, av.Y, -bu.Y, -bv.Y},
              {au.Z, av.Z, -bu.Z, -bv.Z},
              {Dot(l, au), Dot(l, av), 0.0, 0.0}}}};
}

[[nodiscard]] double NormOf(const std::array<double, 4>& f)
{
    return std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2] + f[3] * f[3]);
}

/// The point of the arc at the sample's advance, by Newton's method from
/// the sample, kept within the arc's patches.
[[nodiscard]] std::optional<Sample> Solved(const Marcher& m, const Arc& arc, Sample s)
{
    const double converged = m.Chord * 0x1p-12;
    for (int k = 0; k < NewtonLimit; ++k) {
        const NewtonSystem system = SystemAt(m, arc, s);
        if (NormOf(system.F) <= converged) {
            s.Position = Midpoint(Eval(m.A, s.OnA), Eval(m.B, s.OnB));
            return s;
        }
        const std::optional<DoubleMatrix<4>> inverse = Inverted<4>(system.J);
        if (!inverse) {
            return std::nullopt;
        }
        std::array<double, 4> step{};
        for (std::size_t r = 0; r < 4; ++r) {
            for (std::size_t c = 0; c < 4; ++c) {
                step[r] -= (*inverse)[r][c] * system.F[c];
            }
        }
        s.OnA = ClampedInto(Point2{s.OnA.X + step[0], s.OnA.Y + step[1]}, arc.BoxA);
        s.OnB = ClampedInto(Point2{s.OnB.X + step[2], s.OnB.Y + step[3]}, arc.BoxB);
    }
    return std::nullopt;
}

/// Whether the chord between two samples stays within the marcher's bound
/// of both surfaces and of the arc, at its middle.
[[nodiscard]] bool ChordHolds(const Marcher& m, const Sample& a, const Sample& b, const Sample& mid)
{
    const Point3 chord = Midpoint(a.Position, b.Position);
    const Point2 onA{0.5 * (a.OnA.X + b.OnA.X), 0.5 * (a.OnA.Y + b.OnA.Y)};
    const Point2 onB{0.5 * (a.OnB.X + b.OnB.X), 0.5 * (a.OnB.Y + b.OnB.Y)};
    return Distance(chord, mid.Position) <= m.Chord && Distance(chord, Eval(m.A, onA)) <= m.Chord &&
           Distance(chord, Eval(m.B, onB)) <= m.Chord;
}

/// The samples strictly between a and b, in order: the chord is halved in
/// advance until it holds, the adaptive step of ADR-0010's marching.
[[nodiscard]] Result<void> Refine(const Marcher& m, const Arc& arc, const Sample& a,
                                  const Sample& b, const int depth, std::vector<Sample>& out)
{
    const Sample guess{Midpoint(a.Position, b.Position),
                       Point2{0.5 * (a.OnA.X + b.OnA.X), 0.5 * (a.OnA.Y + b.OnA.Y)},
                       Point2{0.5 * (a.OnB.X + b.OnB.X), 0.5 * (a.OnB.Y + b.OnB.Y)},
                       0.5 * (a.Advance + b.Advance)};
    const std::optional<Sample> mid = Solved(m, arc, guess);
    if (!mid || depth >= RefineDepth) {
        return Fail(ErrorCode::IntersectionFailed,
                    "BSpline x BSpline: marching lost a branch inside patches " +
                        Describe(arc.BoxA) + " against " + Describe(arc.BoxB));
    }
    if (ChordHolds(m, a, b, *mid)) {
        return {};
    }
    STEREON_TRY_VOID(Refine(m, arc, a, *mid, depth + 1, out));
    out.push_back(*mid);
    return Refine(m, arc, *mid, b, depth + 1, out);
}

[[nodiscard]] Sample SampleOf(const Seed& seed, const Vec3& direction)
{
    return {seed.Position, seed.OnA, seed.OnB, Dot(direction, ToVec(seed.Position))};
}

/// The arc's samples, its seeds at the ends, in the arc's own order.
[[nodiscard]] Result<std::vector<Sample>> Marched(const Marcher& m, const Arc& arc)
{
    const Sample first = SampleOf(arc.From, arc.Direction);
    const Sample last = SampleOf(arc.To, arc.Direction);
    std::vector<Sample> samples{first};
    STEREON_TRY_VOID(Refine(m, arc, first, last, 0, samples));
    samples.push_back(last);
    return samples;
}

// --- Branches from arcs -----------------------------------------------------

/// An arc's ends as nodes of the graph of arcs: arcs of neighbouring pairs
/// share a node where they share a seed.
struct ArcGraph {
    std::vector<Point3> Nodes;
    std::vector<std::array<std::size_t, 2>> Ends;
    std::vector<std::vector<std::size_t>> AtNode;
};

[[nodiscard]] std::size_t NodeOf(ArcGraph& g, const Point3& p, const double join)
{
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        if (Distance(g.Nodes[n], p) <= join) {
            return n;
        }
    }
    g.Nodes.push_back(p);
    g.AtNode.emplace_back();
    return g.Nodes.size() - 1;
}

[[nodiscard]] ArcGraph GraphOf(const std::vector<Arc>& arcs, const double join)
{
    ArcGraph g;
    for (std::size_t a = 0; a < arcs.size(); ++a) {
        const std::size_t from = NodeOf(g, arcs[a].From.Position, join);
        const std::size_t to = NodeOf(g, arcs[a].To.Position, join);
        g.Ends.push_back({from, to});
        g.AtNode[from].push_back(a);
        g.AtNode[to].push_back(a);
    }
    return g;
}

/// An arc used by a chain, and whether it is walked from its To end.
struct ChainStep {
    std::size_t Arc;
    bool Reversed;
};

using Chain = std::vector<ChainStep>;

/// The chain from a node along an unused arc, until it stops or closes.
[[nodiscard]] Chain WalkFrom(const ArcGraph& g, std::size_t node, std::vector<bool>& used)
{
    Chain chain;
    for (;;) {
        const auto next =
            std::ranges::find_if(g.AtNode[node], [&](const std::size_t a) { return !used[a]; });
        if (next == g.AtNode[node].end()) {
            return chain;
        }
        used[*next] = true;
        const bool reversed = g.Ends[*next][0] != node;
        chain.push_back({*next, reversed});
        node = g.Ends[*next][reversed ? 0 : 1];
    }
}

/// Every chain: those with free ends first, then the closed loops.
[[nodiscard]] std::vector<Chain> ChainsOf(const ArcGraph& g)
{
    std::vector<bool> used(g.Ends.size(), false);
    std::vector<Chain> chains;
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        if (g.AtNode[n].size() == 1 && !used[g.AtNode[n].front()]) {
            chains.push_back(WalkFrom(g, n, used));
        }
    }
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        Chain chain = WalkFrom(g, n, used);
        if (!chain.empty()) {
            chains.push_back(std::move(chain));
        }
    }
    return chains;
}

/// Whether the parameter lies on a side of its domain, within the join.
[[nodiscard]] bool OnSide(const double x, const Interval& domain)
{
    const double slack = 0x1p-30 * std::max(1.0, Width(domain));
    return x <= domain.GetLower() + slack || x >= domain.GetUpper() - slack;
}

[[nodiscard]] bool OnBoundary(const Seed& s, const Box2& a, const Box2& b)
{
    return OnSide(s.OnA.X, a.X) || OnSide(s.OnA.Y, a.Y) || OnSide(s.OnB.X, b.X) ||
           OnSide(s.OnB.Y, b.Y);
}

/// A node must join exactly two arcs, unless it is an end of the branch on
/// a boundary of the domains; anything else means the arcs do not fit.
[[nodiscard]] Result<void> CheckGraph(const ArcGraph& g, const std::vector<Arc>& arcs,
                                      const Box2& domainA, const Box2& domainB)
{
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        const std::size_t degree = g.AtNode[n].size();
        const Arc& arc = arcs[g.AtNode[n].front()];
        const Seed& seed = g.Ends[g.AtNode[n].front()][0] == n ? arc.From : arc.To;
        if (degree > 2 || (degree == 1 && !OnBoundary(seed, domainA, domainB))) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSpline x BSpline: arcs do not join into branches at (" +
                            Number(seed.Position.X) + ", " + Number(seed.Position.Y) + ", " +
                            Number(seed.Position.Z) + ")");
        }
    }
    return {};
}

/// The chain's samples in order, each arc's shared end used once.
[[nodiscard]] Result<std::vector<Sample>>
ChainSamples(const Marcher& m, const std::vector<Arc>& arcs, const Chain& chain)
{
    std::vector<Sample> samples;
    for (const ChainStep& step : chain) {
        STEREON_TRY(std::vector<Sample> part, Marched(m, arcs[step.Arc]));
        if (step.Reversed) {
            std::ranges::reverse(part);
        }
        samples.insert(samples.end(), part.begin() + (samples.empty() ? 0 : 1), part.end());
    }
    return samples;
}

/// A clamped degree-1 knot vector on the chord lengths of the points, and
/// the poles' parameters. Polylines stay in the domains, as a fitted cubic
/// need not, since Eval refuses a parameter outside them.
[[nodiscard]] std::vector<double> KnotsOf(const std::vector<Sample>& samples)
{
    std::vector<double> knots{0.0, 0.0};
    for (std::size_t k = 1; k < samples.size(); ++k) {
        const double step = Distance(samples[k - 1].Position, samples[k].Position);
        knots.push_back(knots.back() + step);
    }
    knots.push_back(knots.back());
    return knots;
}

/// The samples without repeats: knots must strictly increase.
[[nodiscard]] std::vector<Sample> Distinct(std::vector<Sample> samples)
{
    const auto [first, last] = std::ranges::unique(samples, [](const Sample& a, const Sample& b) {
        return !(Distance(a.Position, b.Position) > 0.0);
    });
    samples.erase(first, last);
    return samples;
}

[[nodiscard]] Result<IntersectionBranch> BranchOf(const std::vector<Sample>& samples)
{
    std::vector<Point3> points;
    std::vector<Point2> onA;
    std::vector<Point2> onB;
    for (const Sample& s : samples) {
        points.push_back(s.Position);
        onA.push_back(s.OnA);
        onB.push_back(s.OnB);
    }
    const std::vector<double> knots = KnotsOf(samples);
    STEREON_TRY(BSplineCurve curve, MakeBSplineCurve(1, knots, std::move(points)));
    STEREON_TRY(BSplineCurve2d curveA, MakeBSplineCurve(1, knots, std::move(onA)));
    STEREON_TRY(BSplineCurve2d curveB, MakeBSplineCurve(1, knots, std::move(onB)));
    return Crossing(Curve(std::move(curve)), Curve2d(std::move(curveA)), Curve2d(std::move(curveB)),
                    Interval{0.0, knots.back()});
}

/// The chain's branch, its chords tightened until it certifies.
[[nodiscard]] Result<void> AddBranch(SurfaceSurfaceIntersection& result, Marcher m,
                                     const std::vector<Arc>& arcs, const Chain& chain,
                                     const Context& ctx)
{
    for (int attempt = 0;; ++attempt) {
        STEREON_TRY(const std::vector<Sample> samples, ChainSamples(m, arcs, chain));
        STEREON_TRY(IntersectionBranch branch, BranchOf(Distinct(samples)));
        const Result<void> certified = Certify(branch, m.A, m.B, ctx);
        if (certified) {
            result.Branches.push_back(std::move(branch));
            return {};
        }
        if (certified.error().GetCode() != ErrorCode::CertificationFailed ||
            attempt >= CertifyRetries) {
            return certified;
        }
        m.Chord *= 0.25;
    }
}

/// The domains actually searched: the request met with the surface's own.
[[nodiscard]] Box2 SearchedDomain(const BSplineSurface& s, const Box2& requested)
{
    const Box2 domain = Domain(s);
    return {Intersection(domain.X, requested.X), Intersection(domain.Y, requested.Y)};
}

[[nodiscard]] Result<std::vector<PatchPair>> AllPairs(const std::vector<Patch>& a,
                                                      const std::vector<Patch>& b)
{
    std::vector<PatchPair> pairs;
    for (const Patch& pa : a) {
        for (const Patch& pb : b) {
            if (pairs.size() >= PairBudget) {
                return Fail(ErrorCode::IntersectionFailed, "Initial patch-pair budget exceeded");
            }
            pairs.push_back({pa, pb});
        }
    }
    return pairs;
}

/// Arcs and touches into a result.
[[nodiscard]] Result<SurfaceSurfaceIntersection>
Assembled(const Findings& found, const Surface& a, const Box2& domainA, const Surface& b,
          const Box2& domainB, const SearchSettings& s, const Context& ctx)
{
    SurfaceSurfaceIntersection result;
    result.Points = found.Touches;
    const ArcGraph graph = GraphOf(found.Arcs, s.Join);
    STEREON_TRY_VOID(CheckGraph(graph, found.Arcs, domainA, domainB));
    const Marcher marcher{a, b, 0.25 * ctx.GetResolution()};
    for (const Chain& chain : ChainsOf(graph)) {
        STEREON_TRY_VOID(AddBranch(result, marcher, found.Arcs, chain, ctx));
    }
    return result;
}

} // namespace

using SSIResult = Result<SurfaceSurfaceIntersection>;

SSIResult NurbsHull(const Surface& a, const Box2& domainA, const Surface& b, const Box2& domainB,
                    const double tolerance, const Context& ctx)
{
    const BSplineSurface* sa = a.Target<BSplineSurface>();
    const BSplineSurface* sb = b.Target<BSplineSurface>();
    if (sa == nullptr || sb == nullptr || !ValidNet(*sa) || !ValidNet(*sb)) {
        return Fail(ErrorCode::InvalidInput, "NurbsHull: expected valid B-spline surfaces");
    }
    if (ProvedApart(NetBounds(*sa), NetBounds(*sb), tolerance)) {
        return SurfaceSurfaceIntersection{}; // Proven empty, including trims.
    }
    STEREON_TRY(const std::vector<Patch> patchesA, PatchesOver(*sa, domainA));
    STEREON_TRY(const std::vector<Patch> patchesB, PatchesOver(*sb, domainB));
    const SearchSettings settings{tolerance, 0x1p-7 * ctx.GetResolution(),
                                  sa->Weights.empty() && sb->Weights.empty()};
    STEREON_TRY(std::vector<PatchPair> pairs, AllPairs(patchesA, patchesB));
    STEREON_TRY(const Findings found, SearchPairs(std::move(pairs), settings));
    return Assembled(found, a, SearchedDomain(*sa, domainA), b, SearchedDomain(*sb, domainB),
                     settings, ctx);
}

} // namespace Stereon::Detail
