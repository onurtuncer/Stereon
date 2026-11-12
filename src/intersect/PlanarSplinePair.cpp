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
#include "BezierArc.hpp"
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
/// continued, is dropped when its bracket is beyond the segment. A bracket
/// that straddles a segment's end, a crossing at or within rounding of a
/// knot, is proved again over the same box in the curves' own parameters
/// with each spline taken whole, piecewise: F at the box's middle from the
/// segment holding it, J the hull of every segment's slope over its part
/// of the box (the end segments continued past the spline's ends). The
/// splines are continuous (an interior knot is repeated at most p times),
/// so F is continuous and piecewise polynomial, F(x) - F(m) lies in
/// J(X) (x - m) by the mean value theorem in integral form, and Krawczyk's
/// test is as sound across a knot, even a corner, as inside a segment.
/// Two brackets that meet on both curves are one root only when one lies
/// in the box the other was proved unique in, where the function proved is
/// the curves' own; otherwise their hull is reported undecided. Equal
/// parameters are never an argument. A pair that is none of
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

// --- Pieces (BezierArc.hpp holds the arcs themselves) ---------------------

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

/// K(X) for the pair, or nothing when F or J cannot be enclosed or J's
/// midpoint is singular.
[[nodiscard]] std::optional<ParameterBox> Krawczyk(const Arc& a, const Arc& b,
                                                   const ParameterBox& x)
{
    const auto pa = PointOf(EvalH(a.Poles, Interval{Midpoint(x.S)}));
    const auto pb = PointOf(EvalH(b.Poles, Interval{Midpoint(x.T)}));
    const auto da = SlopeOver(a, x.S);
    const auto db = SlopeOver(b, x.T);
    if (!pa || !pb || !da || !db) {
        return std::nullopt;
    }
    const std::array<Interval, 2> f{(*pa)[0] - (*pb)[0], (*pa)[1] - (*pb)[1]};
    // J = [[a'x, -b'x], [a'y, -b'y]].
    return KrawczykStep(x, f, {(*da)[0], -(*db)[0], (*da)[1], -(*db)[1]});
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
/// placed, and whether it is decided. A decided one carries a box, in both
/// curves' parameters and within both splines, in which it is proved the
/// only root of the curves' own F (empty when undecided).
struct Meeting {
    Interval OnA;
    Interval OnB;
    bool Decided;
    Interval ProofA = Interval::Empty();
    Interval ProofB = Interval::Empty();
};

/// Both splines' segments (those reaching the requested domains, and all of
/// them, in order), the search's roots and its count of pairs.
struct PairSearch {
    std::vector<Segment> A;
    std::vector<Segment> B;
    std::vector<Segment> AllA;
    std::vector<Segment> AllB;
    std::vector<Meeting> Roots;
    std::size_t Pairs = 0;
    std::string Handler;
};

// --- Across a knot: each spline whole, piecewise -------------------------

/// The local parameter of t in a segment, enclosed.
[[nodiscard]] Interval LocalOf(const Segment& segment, const Interval& t)
{
    return (t - Interval{segment.T0}) / (Interval{segment.T1} - Interval{segment.T0});
}

/// The spline's point at t, from the segment holding t (the first or the
/// last continued beyond the spline's ends); at a knot either neighbour
/// gives the same point, the spline being continuous.
[[nodiscard]] std::optional<std::array<Interval, 2>> PointAt(const std::vector<Segment>& segments,
                                                             const double t)
{
    std::size_t j = 0;
    while (j + 1 < segments.size() && segments[j].T1 < t) {
        ++j;
    }
    return PointOf(EvalH(segments[j].Poles, LocalOf(segments[j], Interval{t})));
}

/// The hull of the spline's slope dC/dt over t: each segment's over its
/// part of t, the end segments continued beyond the spline's ends.
[[nodiscard]] std::optional<std::array<Interval, 2>>
SlopeAcross(const std::vector<Segment>& segments, const Interval& t)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    std::array<Interval, 2> hull{Interval::Empty(), Interval::Empty()};
    for (std::size_t j = 0; j < segments.size(); ++j) {
        const Segment& segment = segments[j];
        const Interval reach{j == 0 ? -inf : segment.T0,
                             j + 1 == segments.size() ? inf : segment.T1};
        const Interval part = Intersection(t, reach);
        if (part.IsEmpty()) {
            continue;
        }
        const Arc whole{j, 0.0, 1.0, segment.Poles};
        const std::optional<std::array<Interval, 2>> d = SlopeOver(whole, LocalOf(segment, part));
        if (!d) {
            return std::nullopt;
        }
        const Interval width = Interval{segment.T1} - Interval{segment.T0};
        hull[0] = Hull(hull[0], (*d)[0] / width);
        hull[1] = Hull(hull[1], (*d)[1] / width);
    }
    return hull;
}

/// K(X) for the splines whole, X in both curves' parameters.
[[nodiscard]] std::optional<ParameterBox> KrawczykAcross(const PairSearch& search,
                                                         const ParameterBox& x)
{
    const auto pa = PointAt(search.AllA, Midpoint(x.S));
    const auto pb = PointAt(search.AllB, Midpoint(x.T));
    const auto da = SlopeAcross(search.AllA, x.S);
    const auto db = SlopeAcross(search.AllB, x.T);
    if (!pa || !pb || !da || !db) {
        return std::nullopt;
    }
    const std::array<Interval, 2> f{(*pa)[0] - (*pb)[0], (*pa)[1] - (*pb)[1]};
    return KrawczykStep(x, f, {(*da)[0], -(*db)[0], (*da)[1], -(*db)[1]});
}

/// The parameters the spline's segments cover, its domain.
[[nodiscard]] Interval Covered(const std::vector<Segment>& segments)
{
    return Interval{segments.front().T0, segments.back().T1};
}

/// The curve parameters a piece-local interval stands for, rounded inward:
/// each end's image enclosed and the side of that enclosure toward the
/// middle taken, so every point of the result is the image of a point of
/// `local` (empty when rounding leaves nothing).
[[nodiscard]] Interval InwardOnCurve(const Segment& segment, const Arc& arc, const Interval& local)
{
    const Interval lower = OnCurve(segment, InSegment(arc, Interval{local.GetLower()}));
    const Interval upper = OnCurve(segment, InSegment(arc, Interval{local.GetUpper()}));
    if (!(lower.GetUpper() <= upper.GetLower())) {
        return Interval::Empty();
    }
    return Interval{lower.GetUpper(), upper.GetLower()};
}

/// The box about c of half-widths h, in doubles (any box will do for the
/// test; it is the box the test proves things in).
[[nodiscard]] Interval About(const double c, const double h)
{
    return Interval{c - h, c + h};
}

/// A box in which Krawczyk's test, the splines whole, proves one root: the
/// pair's own box in the curves' parameters, else boxes about the
/// segment's bracket a factor of eight smaller each time, while they are
/// wider than the bracket; nothing when none is proved.
[[nodiscard]] std::optional<ParameterBox>
ProofBox(const PairSearch& search, const ParameterBox& pair, const ParameterBox& bracket)
{
    if (Classify(KrawczykAcross(search, pair), pair) == Test::OneRoot) {
        return pair;
    }
    const double cs = Midpoint(bracket.S);
    const double ct = Midpoint(bracket.T);
    double hs = Width(pair.S) / 2.0;
    double ht = Width(pair.T) / 2.0;
    for (int shrink = 0; shrink < RefineSteps; ++shrink) {
        hs /= 8.0;
        ht /= 8.0;
        if (!(hs > Width(bracket.S) && ht > Width(bracket.T))) {
            break;
        }
        const ParameterBox x{About(cs, hs), About(ct, ht)};
        if (Classify(KrawczykAcross(search, x), x) == Test::OneRoot) {
            return x;
        }
    }
    return std::nullopt;
}

/// A root whose segment-local bracket straddles a segment's end, proved
/// again with both splines whole; nothing when no box proves it, and the
/// root stays undecided.
[[nodiscard]] std::optional<Meeting> AcrossKnots(const PairSearch& search, const Arc& a,
                                                 const Arc& b, const ParameterBox& local,
                                                 const ParameterBox& bracket)
{
    const std::optional<ParameterBox> found =
        ProofBox(search,
                 {OnCurve(search.A[a.Index], InSegment(a, local.S)),
                  OnCurve(search.B[b.Index], InSegment(b, local.T))},
                 bracket);
    if (!found) {
        return std::nullopt;
    }
    const ParameterBox proof = *found;
    ParameterBox x = proof;
    for (int step = 0; step < RefineSteps; ++step) {
        const std::optional<ParameterBox> k = KrawczykAcross(search, x);
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
    // Past a spline's end its end segment continued was proved, not the
    // spline: a bracket across the end is undecided (and one beyond it is
    // dropped by the caller).
    const Interval onA = Covered(search.AllA);
    const Interval onB = Covered(search.AllB);
    if (!Contains(onA, x.S) || !Contains(onB, x.T)) {
        return Meeting{x.S, x.T, false};
    }
    return Meeting{x.S, x.T, true, Intersection(proof.S, onA), Intersection(proof.T, onB)};
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
    const Segment& segmentA = search.A[a.Index];
    const Segment& segmentB = search.B[b.Index];
    const ParameterBox bracket{OnCurve(segmentA, sa), OnCurve(segmentB, sb)};
    if (wa == Within::Inside && wb == Within::Inside) {
        // Unique in x for the segments' polynomials, which are the curves'
        // own over the segments.
        search.Roots.push_back(
            {bracket.S, bracket.T, true,
             Intersection(InwardOnCurve(segmentA, a, x.S), Interval{segmentA.T0, segmentA.T1}),
             Intersection(InwardOnCurve(segmentB, b, x.T), Interval{segmentB.T0, segmentB.T1})});
        return;
    }
    if (const std::optional<Meeting> across = AcrossKnots(search, a, b, x, bracket)) {
        if (Intersects(across->OnA, Covered(search.AllA)) &&
            Intersects(across->OnB, Covered(search.AllB))) {
            search.Roots.push_back(*across);
        }
        return;
    }
    search.Roots.push_back({bracket.S, bracket.T, false});
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

/// Whether `root`'s bracket lies in the box `proved` was proved the only
/// root of the curves' F in, so that any root of the curves `root` stands
/// for is `proved`'s.
[[nodiscard]] bool ProvedSame(const Meeting& proved, const Meeting& root)
{
    return proved.Decided && Contains(proved.ProofA, root.OnA) && Contains(proved.ProofB, root.OnB);
}

/// Roots whose brackets meet on both curves: one decided root when one's
/// bracket lies in the other's proof box (one crossing bracketed from
/// either side of a knot or from overlapping pieces, or a continued
/// segment's root that, if it is the curves' at all, is the proved one);
/// otherwise their hull, undecided. Brackets apart on either curve are
/// distinct roots.
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
        if (ProvedSame(*same, r)) {
            continue;
        }
        if (ProvedSame(r, *same)) {
            *same = r;
            continue;
        }
        *same = {Hull(same->OnA, r.OnA), Hull(same->OnB, r.OnB), false};
    }
    return distinct;
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
    search.AllA = segmentsA;
    search.AllB = segmentsB;
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
