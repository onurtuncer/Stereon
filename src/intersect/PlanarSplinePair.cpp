// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
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
#include "Planar.hpp"

/// Two B-splines in the plane, polynomial or rational: the system
/// F(s, t) = a(s) - b(t) = 0 over pairs of pieces of their Bezier
/// segments, which Bernstein.hpp encloses (Boehm's insertion in
/// intervals). A pair of pieces is subdivided, both halved, until it is
///
/// - apart: the boxes of the pieces' poles are disjoint, or the hulls of
///   the poles' distances from either piece's chord are (with positive
///   weights each piece lies in the hull of its poles);
/// - empty or single by Krawczyk's operator over the pair's parameter box
///   X, inflated by an eighth of its width on every side:
///   K(X) = m - Y F(m) + (I - Y J(X)) (X - m), with m the box's middle,
///   J(X) the Jacobian [a', -b'] enclosed over X by interval de Casteljau
///   and Y the inverse of J's midpoint. K(X) disjoint from X proves no
///   root in X; K(X) inside X's interior proves exactly one, a transverse
///   crossing (J is then regular), and K is iterated on K(X) ∩ X to
///   bracket it to rounding.
///
/// The inflation lets a root on a piece's end be proved from either side,
/// and one beyond the end of a segment, a root of the segment's polynomial
/// continued, is dropped when its bracket is beyond the segment and
/// reported undecided when the bracket straddles a segment's end. Roots
/// proved twice, from overlapping boxes, are one. A pair that is none of
/// these when the pieces can no longer be halved, or a search past
/// PairBudget, fails with IntersectionFailed: a tangency of two
/// B-splines, an overlap, or a crossing too nearly tangential to separate,
/// is never dropped and never reported unproved (ADR-0029). Two identical
/// splines coincide: their overlap is the requested domains' common
/// stretch on the identity map (the spline's crossings with itself are not
/// reported then).

namespace Stereon::Detail {

namespace {

/// The pairs of pieces the search may visit.
constexpr std::size_t PairBudget = std::size_t{1} << 16;

/// The fraction of a box's width by which Krawczyk's box is inflated.
constexpr double Inflation = 0.125;

/// Krawczyk's refinement: at most this many steps, each shrinking the box.
constexpr int RefineSteps = 64;

// --- Pieces --------------------------------------------------------------

/// A piece [S0, S1] of segment Index's local parameter, dyadic, with its
/// homogeneous poles over it, enclosed, and its hodograph's.
struct Arc {
    std::size_t Index;
    double S0;
    double S1;
    std::vector<HPoint> Poles;
};

[[nodiscard]] HPoint Mid(const HPoint& a, const HPoint& b)
{
    const Interval half{0.5};
    return {(a.X + b.X) * half, (a.Y + b.Y) * half, Interval{0.0}, (a.W + b.W) * half};
}

/// de Casteljau at the middle, enclosed.
[[nodiscard]] std::pair<Arc, Arc> Halves(const Arc& arc)
{
    const std::size_t m = arc.Poles.size();
    std::vector<HPoint> work = arc.Poles;
    std::vector<HPoint> left(m);
    std::vector<HPoint> right(m);
    left[0] = work[0];
    right[m - 1] = work[m - 1];
    for (std::size_t level = 1; level < m; ++level) {
        for (std::size_t i = 0; i + level < m; ++i) {
            work[i] = Mid(work[i], work[i + 1]);
        }
        left[level] = work[0];
        right[m - 1 - level] = work[m - 1 - level];
    }
    const double mid = 0.5 * arc.S0 + 0.5 * arc.S1;
    return {Arc{arc.Index, arc.S0, mid, std::move(left)},
            Arc{arc.Index, mid, arc.S1, std::move(right)}};
}

/// b + s (c - b): de Casteljau's step in the form whose enclosure shrinks
/// with the poles' spread, for s anywhere, inside [0, 1] or not.
[[nodiscard]] Interval Step(const Interval& b, const Interval& c, const Interval& s)
{
    return b + s * (c - b);
}

/// The homogeneous point of the poles' Bezier polynomial at s, enclosed.
[[nodiscard]] HPoint EvalH(std::vector<HPoint> work, const Interval& s)
{
    for (std::size_t level = 1; level < work.size(); ++level) {
        for (std::size_t i = 0; i + level < work.size(); ++i) {
            work[i] = {Step(work[i].X, work[i + 1].X, s), Step(work[i].Y, work[i + 1].Y, s),
                       Interval{0.0}, Step(work[i].W, work[i + 1].W, s)};
        }
    }
    return work.front();
}

/// The hodograph's poles, p (H_i+1 - H_i), of degree p - 1; one zero pole
/// for a constant.
[[nodiscard]] std::vector<HPoint> Hodograph(const std::vector<HPoint>& poles)
{
    const Interval p{static_cast<double>(poles.size() - 1)};
    std::vector<HPoint> d;
    for (std::size_t i = 0; i + 1 < poles.size(); ++i) {
        d.push_back({p * (poles[i + 1].X - poles[i].X), p * (poles[i + 1].Y - poles[i].Y),
                     Interval{0.0}, p * (poles[i + 1].W - poles[i].W)});
    }
    if (d.empty()) {
        d.push_back({Interval{0.0}, Interval{0.0}, Interval{0.0}, Interval{0.0}});
    }
    return d;
}

/// The point x / w enclosed, or nothing when w is not proved positive.
[[nodiscard]] std::optional<std::array<Interval, 2>> PointOf(const HPoint& h)
{
    if (!(h.W.GetLower() > 0.0)) {
        return std::nullopt;
    }
    return std::array{h.X / h.W, h.Y / h.W};
}

/// The derivative (X' W - X W') / W² of the piece over s, enclosed.
[[nodiscard]] std::optional<std::array<Interval, 2>> SlopeOver(const Arc& arc, const Interval& s)
{
    const HPoint h = EvalH(arc.Poles, s);
    const HPoint d = EvalH(Hodograph(arc.Poles), s);
    if (!(h.W.GetLower() > 0.0)) {
        return std::nullopt;
    }
    const Interval ww = h.W * h.W;
    return std::array{(d.X * h.W - h.X * d.W) / ww, (d.Y * h.W - h.Y * d.W) / ww};
}

// --- Separation ----------------------------------------------------------

/// The box of the piece's Euclidean poles, which holds the piece.
[[nodiscard]] std::optional<std::array<Interval, 2>> BoxOf(const Arc& arc)
{
    std::array<Interval, 2> box{Interval::Empty(), Interval::Empty()};
    for (const HPoint& h : arc.Poles) {
        const std::optional<std::array<Interval, 2>> p = PointOf(h);
        if (!p) {
            return std::nullopt;
        }
        box[0] = Hull(box[0], (*p)[0]);
        box[1] = Hull(box[1], (*p)[1]);
    }
    return box;
}

/// The hull of n . (P_i - q) over the piece's poles: it holds the piece's
/// distances from the line through q across n, whatever n and q are.
[[nodiscard]] Interval DistancesOf(const Arc& arc, const Vec2& n, const Point2& q)
{
    Interval hull = Interval::Empty();
    for (const HPoint& h : arc.Poles) {
        if (const std::optional<std::array<Interval, 2>> p = PointOf(h)) {
            hull = Hull(hull, ((*p)[0] - Interval{q.X}) * Interval{n.X} +
                                  ((*p)[1] - Interval{q.Y}) * Interval{n.Y});
        } else {
            return Interval::Everything();
        }
    }
    return hull;
}

/// The middle of a pole's enclosure, as a point.
[[nodiscard]] Point2 Nominal(const HPoint& h)
{
    return {Midpoint(h.X) / Midpoint(h.W), Midpoint(h.Y) / Midpoint(h.W)};
}

/// Whether the hulls of both pieces' distances from `chord`'s chord, the
/// line through its end poles, are disjoint (a fat line, after Sederberg
/// and Nishita).
[[nodiscard]] bool ApartAcross(const Arc& chord, const Arc& other)
{
    const Point2 first = Nominal(chord.Poles.front());
    const Vec2 n = Perpendicular(Nominal(chord.Poles.back()) - first);
    if (!(Dot(n, n) > 0.0) || !std::isfinite(Dot(n, n))) {
        return false;
    }
    return !Intersects(DistancesOf(chord, n, first), DistancesOf(other, n, first));
}

[[nodiscard]] bool Apart(const Arc& a, const Arc& b)
{
    const std::optional<std::array<Interval, 2>> boxA = BoxOf(a);
    const std::optional<std::array<Interval, 2>> boxB = BoxOf(b);
    if (boxA && boxB &&
        (!Intersects((*boxA)[0], (*boxB)[0]) || !Intersects((*boxA)[1], (*boxB)[1]))) {
        return true;
    }
    return ApartAcross(a, b) || ApartAcross(b, a);
}

// --- Krawczyk's operator -------------------------------------------------

/// The outcome of Krawczyk's test on a box.
enum class Test : std::uint8_t { NoRoot, OneRoot, Unknown };

/// The parameter box of a pair: S on A's piece, T on B's, local to each
/// piece ([0, 1] is the piece itself).
struct ParameterBox {
    Interval S;
    Interval T;
};

/// K(X) for the pair, or nothing when F or J cannot be enclosed or J's
/// midpoint is singular.
[[nodiscard]] std::optional<ParameterBox> Krawczyk(const Arc& a, const Arc& b,
                                                   const ParameterBox& x)
{
    const double ms = Midpoint(x.S);
    const double mt = Midpoint(x.T);
    const auto pa = PointOf(EvalH(a.Poles, Interval{ms}));
    const auto pb = PointOf(EvalH(b.Poles, Interval{mt}));
    const auto da = SlopeOver(a, x.S);
    const auto db = SlopeOver(b, x.T);
    if (!pa || !pb || !da || !db) {
        return std::nullopt;
    }
    const std::array<Interval, 2> f{(*pa)[0] - (*pb)[0], (*pa)[1] - (*pb)[1]};
    // J = [[a'x, -b'x], [a'y, -b'y]].
    const std::array<Interval, 4> j{(*da)[0], -(*db)[0], (*da)[1], -(*db)[1]};
    const double c11 = Midpoint(j[0]);
    const double c12 = Midpoint(j[1]);
    const double c21 = Midpoint(j[2]);
    const double c22 = Midpoint(j[3]);
    const double det = c11 * c22 - c12 * c21;
    if (!(std::abs(det) > 0.0) || !std::isfinite(det)) {
        return std::nullopt;
    }
    const std::array<Interval, 4> y{Interval{c22 / det}, Interval{-c12 / det}, Interval{-c21 / det},
                                    Interval{c11 / det}};
    const Interval ds = x.S - Interval{ms};
    const Interval dt = x.T - Interval{mt};
    const auto row = [&](const std::size_t r, const double m, const double unitS,
                         const double unitT) {
        const Interval yf = y[2 * r] * f[0] + y[2 * r + 1] * f[1];
        const Interval rs = Interval{unitS} - (y[2 * r] * j[0] + y[2 * r + 1] * j[2]);
        const Interval rt = Interval{unitT} - (y[2 * r] * j[1] + y[2 * r + 1] * j[3]);
        return Interval{m} - yf + rs * ds + rt * dt;
    };
    return ParameterBox{row(0, ms, 1.0, 0.0), row(1, mt, 0.0, 1.0)};
}

[[nodiscard]] bool InsideInterior(const Interval& inner, const Interval& outer)
{
    return outer.GetLower() < inner.GetLower() && inner.GetUpper() < outer.GetUpper();
}

[[nodiscard]] Test Classify(const std::optional<ParameterBox>& k, const ParameterBox& x)
{
    if (!k) {
        return Test::Unknown;
    }
    if (!Intersects(k->S, x.S) || !Intersects(k->T, x.T)) {
        return Test::NoRoot;
    }
    return InsideInterior(k->S, x.S) && InsideInterior(k->T, x.T) ? Test::OneRoot : Test::Unknown;
}

/// The single root of a box Krawczyk's test proved to hold one, bracketed
/// by iterating X <- K(X) ∩ X, which keeps it, while the box shrinks.
[[nodiscard]] ParameterBox Refined(const Arc& a, const Arc& b, ParameterBox x)
{
    for (int step = 0; step < RefineSteps; ++step) {
        const std::optional<ParameterBox> k = Krawczyk(a, b, x);
        if (!k) {
            break;
        }
        const ParameterBox next{Intersection(k->S, x.S), Intersection(k->T, x.T)};
        if (next.S.IsEmpty() || next.T.IsEmpty() ||
            (Width(next.S) >= Width(x.S) && Width(next.T) >= Width(x.T))) {
            break;
        }
        x = next;
    }
    return x;
}

// --- The search ----------------------------------------------------------

/// A crossing found: its brackets in both curves' parameters, where it is
/// placed, and whether it is decided.
struct Meeting {
    Interval OnA;
    Interval OnB;
    bool Decided;
};

/// Both splines' segments, the search's roots and its count of pairs.
struct PairSearch {
    std::vector<Segment> A;
    std::vector<Segment> B;
    std::vector<Meeting> Roots;
    std::size_t Pairs = 0;
    std::string Handler;
};

/// A piece-local box in a segment's local parameter, then the curve's.
[[nodiscard]] Interval InSegment(const Arc& arc, const Interval& local)
{
    return Interval{arc.S0} + (Interval{arc.S1} - Interval{arc.S0}) * local;
}

[[nodiscard]] Interval OnCurve(const Segment& segment, const Interval& s)
{
    return Interval{segment.T0} + (Interval{segment.T1} - Interval{segment.T0}) * s;
}

/// Where a root's bracket in a segment's parameter lies against [0, 1]:
/// inside (a root of the curve), beyond (of the polynomial continued, not
/// of the segment), or across an end (undecided).
enum class Within : std::uint8_t { Inside, Beyond, Across };

[[nodiscard]] Within WithinSegment(const Interval& s)
{
    if (Contains(Interval{0.0, 1.0}, s)) {
        return Within::Inside;
    }
    return Intersects(s, Interval{0.0, 1.0}) ? Within::Across : Within::Beyond;
}

/// Records the proved root of the pair's inflated box.
void Record(PairSearch& search, const Arc& a, const Arc& b, const ParameterBox& x)
{
    const ParameterBox root = Refined(a, b, x);
    const Interval sa = InSegment(a, root.S);
    const Interval sb = InSegment(b, root.T);
    const Within wa = WithinSegment(sa);
    const Within wb = WithinSegment(sb);
    if (wa == Within::Beyond || wb == Within::Beyond) {
        return;
    }
    search.Roots.push_back({OnCurve(search.A[a.Index], sa), OnCurve(search.B[b.Index], sb),
                            wa == Within::Inside && wb == Within::Inside});
}

[[nodiscard]] KernelError Unresolved(const PairSearch& search, const Arc& a, const Arc& b,
                                     const std::string_view why)
{
    const double ta = Midpoint(OnCurve(search.A[a.Index], InSegment(a, Interval{0.5})));
    const double tb = Midpoint(OnCurve(search.B[b.Index], InSegment(b, Interval{0.5})));
    return Fail(ErrorCode::IntersectionFailed,
                search.Handler + ": " + std::string(why) + " near (" + Number(ta) + ", " +
                    Number(tb) +
                    ") (a tangency, an overlap or a nearly tangential crossing; ADR-0010: an "
                    "unresolved pair fails, it is never dropped)")
        .error();
}

/// One pair of pieces, apart, settled by Krawczyk's test, or halved.
[[nodiscard]] Result<void> VisitPair(PairSearch& search, const Arc& a, const Arc& b,
                                     const int depth)
{
    if (++search.Pairs > PairBudget) {
        return std::unexpected(Unresolved(search, a, b, "the subdivision budget is spent"));
    }
    if (Apart(a, b)) {
        return {};
    }
    const ParameterBox x{Interval{-Inflation, 1.0 + Inflation},
                         Interval{-Inflation, 1.0 + Inflation}};
    const Test test = Classify(Krawczyk(a, b, x), x);
    if (test == Test::NoRoot) {
        return {};
    }
    if (test == Test::OneRoot) {
        Record(search, a, b, x);
        return {};
    }
    if (depth == MaxDepth) {
        return std::unexpected(
            Unresolved(search, a, b, "a pair of pieces neither apart nor a single crossing"));
    }
    const std::pair<Arc, Arc> ha = Halves(a);
    const std::pair<Arc, Arc> hb = Halves(b);
    for (const Arc* pa : {&ha.first, &ha.second}) {
        for (const Arc* pb : {&hb.first, &hb.second}) {
            STEREON_TRY_VOID(VisitPair(search, *pa, *pb, depth + 1));
        }
    }
    return {};
}

/// Roots whose brackets meet on both curves are one, proved from boxes
/// that overlap; with an undecided one, undecided.
[[nodiscard]] std::vector<Meeting> Distinct(std::vector<Meeting> roots)
{
    std::ranges::stable_sort(roots, [](const Meeting& p, const Meeting& q) {
        return p.OnA.GetLower() < q.OnA.GetLower();
    });
    std::vector<Meeting> distinct;
    for (const Meeting& r : roots) {
        const auto same = std::ranges::find_if(distinct, [&](const Meeting& d) {
            return Intersects(d.OnA, r.OnA) && Intersects(d.OnB, r.OnB);
        });
        if (same == distinct.end()) {
            distinct.push_back(r);
            continue;
        }
        if (!(same->Decided && r.Decided)) {
            same->OnA = Hull(same->OnA, r.OnA);
            same->OnB = Hull(same->OnB, r.OnB);
            same->Decided = false;
        }
    }
    return distinct;
}

/// A root's parameter in the requested domain: its bracket's middle, or,
/// for a bracket across the domain's end, the end, undecided; nothing for
/// a bracket outside the domain.
[[nodiscard]] std::optional<double> Placed(const Interval& bracket, const Interval& domain,
                                           bool& decided)
{
    if (!Intersects(bracket, domain)) {
        return std::nullopt;
    }
    const double t = Midpoint(bracket);
    if (!Contains(domain, bracket)) {
        decided = false;
    }
    return std::clamp(t, domain.GetLower(), domain.GetUpper());
}

[[nodiscard]] std::vector<Segment> Reaching(std::vector<Segment> segments, const Interval& domain)
{
    std::erase_if(segments,
                  [&](const Segment& s) { return !Intersects(Interval{s.T0, s.T1}, domain); });
    return segments;
}

/// Every pair of segments, searched.
[[nodiscard]] Result<void> SearchAll(PairSearch& search)
{
    for (std::size_t i = 0; i < search.A.size(); ++i) {
        for (std::size_t j = 0; j < search.B.size(); ++j) {
            const Arc a{i, 0.0, 1.0, search.A[i].Poles};
            const Arc b{j, 0.0, 1.0, search.B[j].Poles};
            STEREON_TRY_VOID(VisitPair(search, a, b, 0));
        }
    }
    return {};
}

/// Identical splines coincide over their common requested stretch.
[[nodiscard]] PlanarIntersection Identical(const BSplineCurve2d& spline, const Interval& domainA,
                                           const Interval& domainB)
{
    PlanarIntersection result;
    const Interval common = Intersection(Intersection(domainA, domainB), Domain(spline));
    if (!common.IsEmpty()) {
        result.Overlaps.push_back({common, common, true});
    }
    return result;
}

[[nodiscard]] Result<void> CheckBoth(const BSplineCurve2d& a, const Interval& domainA,
                                     const BSplineCurve2d& b, const Interval& domainB,
                                     const std::string_view handler)
{
    STEREON_TRY_VOID(CheckPlaneSpline(a, handler));
    STEREON_TRY_VOID(CheckPlaneSpline(b, handler));
    STEREON_TRY_VOID(CheckPeriodicDomain(a, domainA, handler));
    return CheckPeriodicDomain(b, domainB, handler);
}

} // namespace

Result<PlanarIntersection> SplineSpline2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& /*ctx*/)
{
    const BSplineCurve2d* sa = a.Target<BSplineCurve2d>();
    const BSplineCurve2d* sb = b.Target<BSplineCurve2d>();
    if (sa == nullptr || sb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineSpline2d: dispatched for other curve types");
    }
    PairSearch search;
    search.Handler = "BSplineCurve2d x BSplineCurve2d";
    STEREON_TRY_VOID(CheckBoth(*sa, domainA, *sb, domainB, search.Handler));
    if (*sa == *sb) {
        return Identical(*sa, domainA, domainB);
    }
    STEREON_TRY(std::vector<Segment> segmentsA, SegmentsOf(*sa));
    STEREON_TRY(std::vector<Segment> segmentsB, SegmentsOf(*sb));
    search.A = Reaching(std::move(segmentsA), domainA);
    search.B = Reaching(std::move(segmentsB), domainB);
    STEREON_TRY_VOID(SearchAll(search));
    // A bracket across a spline's end is placed on the spline, in the domain.
    const Interval onA = Intersection(domainA, Domain(*sa));
    const Interval onB = Intersection(domainB, Domain(*sb));
    PlanarIntersection result;
    if (onA.IsEmpty() || onB.IsEmpty()) {
        return result;
    }
    for (const Meeting& root : Distinct(std::move(search.Roots))) {
        bool decided = root.Decided;
        const std::optional<double> ta = Placed(root.OnA, onA, decided);
        const std::optional<double> tb = Placed(root.OnB, onB, decided);
        if (ta && tb) {
            result.Points.push_back({*ta, *tb, Eval(*sa, *ta), ContactType::Crossing, decided});
        }
    }
    Order(result);
    return result;
}

} // namespace Stereon::Detail
