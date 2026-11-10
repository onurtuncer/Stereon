// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Interval-enclosed Bezier patches of a B-spline surface and Krawczyk's
// operator on N-dimensional systems, shared by the generic surface-surface
// intersection (GenericSSI.cpp) and curves against B-spline surfaces
// (CurveSplineSurface.cpp). Moved from GenericSSI.cpp unchanged.
//
// Included in the purview of the handlers' implementation units, after
// Bernstein.hpp, with <algorithm>, <array>, <cmath>, <cstddef>, <cstdint>,
// <numbers>, <optional>, <utility> and <vector> in their global module
// fragment. Everything is in an anonymous namespace, so each unit gets its
// own copy.

#pragma once

namespace Stereon::Detail {
namespace {

using Bounds = std::array<Interval, 3>;

// --- Enclosed control nets -----------------------------------------------

/// A grid of enclosed homogeneous poles, i along u: index i * CountV + j.
struct Net {
    std::size_t CountU = 0;
    std::size_t CountV = 0;
    std::vector<HPoint> Poles;
};

[[maybe_unused]] [[nodiscard]] inline std::size_t LineCount(const Net& net, const bool alongU)
{
    return alongU ? net.CountV : net.CountU;
}

[[maybe_unused]] [[nodiscard]] inline std::size_t LineLength(const Net& net, const bool alongU)
{
    return alongU ? net.CountU : net.CountV;
}

/// Where the pole `at` of line `line` is held.
[[maybe_unused]] [[nodiscard]] inline std::size_t
PoleIndex(const Net& net, const bool alongU, const std::size_t line, const std::size_t at)
{
    return alongU ? at * net.CountV + line : line * net.CountV + at;
}

[[maybe_unused]] [[nodiscard]] inline std::vector<HPoint> LineOf(const Net& net, const bool alongU,
                                                                 const std::size_t line)
{
    std::vector<HPoint> poles(LineLength(net, alongU));
    for (std::size_t at = 0; at < poles.size(); ++at) {
        poles[at] = net.Poles[PoleIndex(net, alongU, line, at)];
    }
    return poles;
}

[[maybe_unused]] inline void SetLine(Net& net, const bool alongU, const std::size_t line,
                                     const std::vector<HPoint>& poles)
{
    for (std::size_t at = 0; at < poles.size(); ++at) {
        net.Poles[PoleIndex(net, alongU, line, at)] = poles[at];
    }
}

[[maybe_unused]] [[nodiscard]] inline Net NetOf(const BSplineSurface& s)
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
[[maybe_unused]] [[nodiscard]] inline Net WithLines(const Net& net, const bool alongU,
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
[[maybe_unused]] [[nodiscard]] inline Result<void>
RaiseKnots(Net& net, std::vector<double>& knots, const std::size_t p, const bool alongU)
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
[[maybe_unused]] [[nodiscard]] inline std::vector<std::size_t>
SpansOf(const std::vector<double>& knots, const std::size_t p, const std::size_t poles,
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
[[maybe_unused]] [[nodiscard]] inline Bounds ExtentOf(const Net& net)
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

[[maybe_unused]] [[nodiscard]] inline double LargestWidth(const Bounds& box)
{
    double largest = 0.0;
    for (const Interval& side : box) {
        largest = std::max(largest, Width(side));
    }
    return largest;
}

[[maybe_unused]] [[nodiscard]] inline Patch MakePatch(const Box2& parameters, Net net)
{
    Bounds extent = ExtentOf(net);
    return {parameters, std::move(net), extent};
}

/// The patch of spans ending at poles (ju, jv) of a net of Bezier form.
[[maybe_unused]] [[nodiscard]] inline Net SubNet(const Net& net, const std::size_t p,
                                                 const std::size_t q, const std::size_t ju,
                                                 const std::size_t jv)
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
[[maybe_unused]] [[nodiscard]] inline Result<std::vector<Patch>> PatchesOf(const BSplineSurface& s)
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
[[maybe_unused]] inline void SplitLine(const std::vector<HPoint>& line, const Interval& t,
                                       std::vector<HPoint>& left, std::vector<HPoint>& right)
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
[[maybe_unused]] [[nodiscard]] inline std::optional<std::pair<Patch, Patch>>
Split(const Patch& patch, const bool alongU, const double at)
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
[[maybe_unused]] [[nodiscard]] inline std::optional<Patch> Clipped(Patch patch, const bool alongU,
                                                                   const Interval& wanted)
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
[[maybe_unused]] [[nodiscard]] inline Result<std::vector<Patch>>
PatchesOver(const BSplineSurface& s, const Box2& domain)
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

// --- Enclosures of a patch and its derivatives -----------------------------

/// Splits fall at this fraction of a stretch, not at its middle: symmetric
/// input would otherwise lay a branch along a split line, where a boundary
/// root is not isolated. sqrt(2) / 3.
constexpr double SplitFraction = std::numbers::sqrt2 / 3.0;

/// The narrowings of one Krawczyk search.
constexpr int KrawczykLimit = 48;

[[maybe_unused]] [[nodiscard]] inline Bounds Unbounded()
{
    return {Interval::Everything(), Interval::Everything(), Interval::Everything()};
}

[[maybe_unused]] [[nodiscard]] inline bool AllBounded(const Bounds& b)
{
    return b[0].IsBounded() && b[1].IsBounded() && b[2].IsBounded();
}

[[maybe_unused]] [[nodiscard]] inline Bounds CrossOf(const Bounds& a, const Bounds& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[maybe_unused]] [[nodiscard]] inline Interval DotOf(const Bounds& a, const Bounds& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[maybe_unused]] [[nodiscard]] inline Bounds MinusOf(const Bounds& a, const Bounds& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[maybe_unused]] [[nodiscard]] inline HPoint HullOf(const Net& net)
{
    HPoint hull;
    for (const HPoint& h : net.Poles) {
        hull = {Hull(hull.X, h.X), Hull(hull.Y, h.Y), Hull(hull.Z, h.Z), Hull(hull.W, h.W)};
    }
    return hull;
}

[[maybe_unused]] [[nodiscard]] inline Interval WidthOf(const Interval& stretch)
{
    return Interval{stretch.GetUpper()} - Interval{stretch.GetLower()};
}

/// The net of the partial derivative along one direction, its degree one
/// lower; a net of one zero pole across at degree 0. The stretch's width
/// must be proved positive.
[[maybe_unused]] [[nodiscard]] inline Net DerivativeNet(const Net& net, const bool alongU,
                                                        const Interval& width)
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

[[maybe_unused]] [[nodiscard]] inline Bounds QuotientOf(const HPoint& d, const Bounds& point,
                                                        const Interval& weight)
{
    return {(d.X - point[0] * d.W) / weight, (d.Y - point[1] * d.W) / weight,
            (d.Z - point[2] * d.W) / weight};
}

[[maybe_unused]] [[nodiscard]] inline Hulls HullsOf(const Patch& patch)
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
[[maybe_unused]] [[nodiscard]] inline Interval FractionOf(const Interval& stretch,
                                                          const Interval& x)
{
    return (x - Interval{stretch.GetLower()}) / WidthOf(stretch);
}

[[maybe_unused]] [[nodiscard]] inline HPoint Reduced(std::vector<HPoint> line, const Interval& t)
{
    for (std::size_t r = line.size() - 1; r > 0; --r) {
        for (std::size_t i = 0; i < r; ++i) {
            line[i] = Blend(line[i], line[i + 1], t);
        }
    }
    return line.front();
}

/// The net's homogeneous value at enclosed parameters of its stretch.
[[maybe_unused]] [[nodiscard]] inline HPoint ValueOf(const Net& net, const Box2& stretch,
                                                     const Interval& u, const Interval& v)
{
    std::vector<HPoint> across;
    for (std::size_t i = 0; i < net.CountU; ++i) {
        across.push_back(Reduced(LineOf(net, false, i), FractionOf(stretch.Y, v)));
    }
    return Reduced(std::move(across), FractionOf(stretch.X, u));
}

[[maybe_unused]] [[nodiscard]] inline Bounds CartesianOf(const HPoint& h)
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
[[maybe_unused]] [[nodiscard]] inline std::size_t PivotRow(const DoubleMatrix<N>& a,
                                                           const std::size_t column)
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
[[maybe_unused]] [[nodiscard]] inline std::optional<DoubleMatrix<N>> Inverted(DoubleMatrix<N> a)
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
[[maybe_unused]] [[nodiscard]] inline bool BoundedMatrix(const IntervalMatrix<N>& j)
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
[[maybe_unused]] [[nodiscard]] inline Contraction<N> ContractionOf(const DoubleMatrix<N>& y,
                                                                   const IntervalMatrix<N>& j)
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
[[maybe_unused]] [[nodiscard]] inline std::array<double, N> MidpointsOf(const IntervalsN<N>& x)
{
    std::array<double, N> m{};
    for (std::size_t i = 0; i < N; ++i) {
        m[i] = Midpoint(x[i]);
    }
    return m;
}

/// One application of the operator, m - Y F(m) + (I - Y J)(X - m), met with X.
template <std::size_t N>
[[maybe_unused]] [[nodiscard]] inline std::optional<IntervalsN<N>>
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
[[maybe_unused]] [[nodiscard]] inline double TotalWidth(const IntervalsN<N>& x)
{
    double total = 0.0;
    for (const Interval& e : x) {
        total += Width(e);
    }
    return total;
}

/// Whether a narrowed box is small enough for its middle to be the root.
template <std::size_t N>
[[maybe_unused]] [[nodiscard]] inline bool Tight(const IntervalsN<N>& x)
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
[[maybe_unused]] [[nodiscard]] inline Isolated<N> Krawczyk(IntervalsN<N> x,
                                                           const IntervalMatrix<N>& j, const F& f)
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

// --- Patch helpers --------------------------------------------------------

[[maybe_unused]] [[nodiscard]] inline Interval StretchOf(const Patch& patch, const bool alongU)
{
    return alongU ? patch.Parameters.X : patch.Parameters.Y;
}

[[maybe_unused]] [[nodiscard]] inline Point3 Middle(const Bounds& b)
{
    return {Midpoint(b[0]), Midpoint(b[1]), Midpoint(b[2])};
}

/// A stretch reaching a few roundings past its ends. A root on a line that
/// two cells share would otherwise be refused by both, each by a rounding.
[[maybe_unused]] [[nodiscard]] inline Interval Reaching(const Interval& stretch)
{
    const double magnitude =
        std::max({1.0, std::abs(stretch.GetLower()), std::abs(stretch.GetUpper())});
    return Inflated(stretch, 0x1p-44 * magnitude);
}

/// A derivative's enclosure over a patch, widened to hold on the stretch
/// that reaches past it: the extra reach is some 2^-44 of the parameter, so
/// the derivative moves by far less than this.
[[maybe_unused]] [[nodiscard]] inline Interval Slackened(const Interval& d)
{
    if (!d.IsBounded()) {
        return d;
    }
    return Inflated(d, 0x1p-30 * std::max(std::abs(d.GetLower()), std::abs(d.GetUpper())));
}
[[maybe_unused]] [[nodiscard]] inline std::optional<std::pair<Patch, Patch>>
Halves(const Patch& patch)
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

} // namespace
} // namespace Stereon::Detail
