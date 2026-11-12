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
/// intervals), both splines first moved by the middle of their poles when
/// every pole's move is exact (Centred), so the enclosures' rounding is
/// relative to the curves' size, not their distance from the origin. A
/// pair of pieces is subdivided, both halved, until it is
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
/// parameters are never an argument. A pair that is neither is tried as
/// two graphs over one direction before it is halved (Touches): a
/// tangency, where J is singular, is followed by the signs of the
/// derivatives of the graphs' height difference and reported bracketed,
/// Tangential and undecided. A pair that is none of these when the pieces
/// can no longer be halved, or a search past PairBudget, fails with
/// IntersectionFailed: an overlap, a touch no derivative of order up to
/// four settles, or a crossing too nearly tangential to separate, is never
/// dropped and never reported unproved (ADR-0029). Two identical splines,
/// or one curve in two knot vectors refined exactly (SameCurve), coincide:
/// their overlap is the requested domains' common stretch on the identity
/// map (the spline's crossings with itself are not reported then).

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
/// only root of the curves' own F (empty when undecided), and its contact:
/// a crossing, or a touch bracketed to rounding (Tangential, undecided).
struct Meeting {
    Interval OnA;
    Interval OnB;
    bool Decided;
    Interval ProofA = Interval::Empty();
    Interval ProofB = Interval::Empty();
    ContactType Contact = ContactType::Crossing;
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

// --- Touches: both pieces graphs over one direction ------------------------

/// The highest derivative a contact's order is looked for with.
constexpr std::size_t JetOrder = 4;

/// A function's derivatives of order 0 to JetOrder, each enclosed over one
/// parameter interval (all of it, so every pointwise identity between them
/// holds for the enclosures).
using Jet = std::array<Interval, JetOrder + 1>;

/// n / d as jets to the given order, by Leibniz's rule; nothing when d may
/// vanish.
[[nodiscard]] std::optional<Jet> Quotient(const Jet& n, const Jet& d, const std::size_t order)
{
    if (StrictSign(d[0]) == Strict::Unknown) {
        return std::nullopt;
    }
    Jet q;
    q.fill(Interval{0.0});
    for (std::size_t m = 0; m <= order; ++m) {
        Interval sum = n[m];
        for (std::size_t i = 1; i <= m; ++i) {
            sum = sum - Interval{Binomial(m, i)} * d[i] * q[m - i];
        }
        q[m] = sum / d[0];
    }
    return q;
}

/// The jet one order down: f' from f.
[[nodiscard]] Jet Shifted(const Jet& f)
{
    Jet d;
    d.fill(Interval{0.0});
    for (std::size_t i = 0; i < JetOrder; ++i) {
        d[i] = f[i + 1];
    }
    return d;
}

/// Coordinates along U and N from Origin: the frame both pieces are graphs
/// in, x along U and y along N.
struct Frame {
    Point2 Origin;
    Vec2 U;
    Vec2 N;
};

/// The piece's x and y as jets in its parameter over s; nothing when W is
/// not proved positive there.
[[nodiscard]] std::optional<std::pair<Jet, Jet>> CoordinatesOver(const Arc& arc, const Frame& f,
                                                                 const Interval& s)
{
    std::vector<HPoint> poles = arc.Poles;
    Jet x;
    Jet y;
    Jet w;
    for (std::size_t j = 0; j <= JetOrder; ++j) {
        // Derivatives of the homogeneous point, and of linear forms of it.
        const HPoint h = EvalH(poles, s);
        const Interval dx = h.X - h.W * Interval{f.Origin.X};
        const Interval dy = h.Y - h.W * Interval{f.Origin.Y};
        x[j] = dx * Interval{f.U.X} + dy * Interval{f.U.Y};
        y[j] = dx * Interval{f.N.X} + dy * Interval{f.N.Y};
        w[j] = h.W;
        poles = Hodograph(poles);
    }
    const std::optional<Jet> qx = Quotient(x, w, JetOrder);
    const std::optional<Jet> qy = Quotient(y, w, JetOrder);
    if (!qx || !qy) {
        return std::nullopt;
    }
    return std::pair{*qx, *qy};
}

/// The graph's d^j y / dx^j for j = 0 to JetOrder, enclosed over s, by
/// d/dx = (1 / x') d/ds; nothing when x' may vanish there.
[[nodiscard]] std::optional<Jet> GraphOver(const Arc& arc, const Frame& f, const Interval& s)
{
    const std::optional<std::pair<Jet, Jet>> c = CoordinatesOver(arc, f, s);
    if (!c) {
        return std::nullopt;
    }
    const Jet slope = Shifted(c->first);
    Jet graph;
    graph[0] = c->second[0];
    Jet g = c->second;
    for (std::size_t k = 1; k <= JetOrder; ++k) {
        const std::optional<Jet> next = Quotient(Shifted(g), slope, JetOrder - k);
        if (!next) {
            return std::nullopt;
        }
        g = *next;
        graph[k] = g[0];
    }
    return graph;
}

/// The piece's x at s, enclosed.
[[nodiscard]] Interval XAt(const Arc& arc, const Frame& f, const double s)
{
    const std::optional<std::array<Interval, 2>> p = PointOf(EvalH(arc.Poles, Interval{s}));
    if (!p) {
        return Interval::Everything();
    }
    return ((*p)[0] - Interval{f.Origin.X}) * Interval{f.U.X} +
           ((*p)[1] - Interval{f.Origin.Y}) * Interval{f.U.Y};
}

/// A piece whose x is strictly monotone, with that sense, +1 or -1.
struct Graph {
    const Arc* Of;
    double Sense;
};

/// x along the graph's sense, so that it increases with s.
[[nodiscard]] Interval Oriented(const Graph& g, const Interval& x)
{
    return g.Sense > 0.0 ? x : -x;
}

/// The piece's oriented x at s, enclosed.
[[nodiscard]] Interval RisingAt(const Graph& g, const Frame& f, const double s)
{
    return Oriented(g, XAt(*g.Of, f, s));
}

/// Bisection on [0, 1] for where `left` turns false, given it holds at 0
/// and not at 1: the last s found holding and the first found not.
template <typename Left>
[[nodiscard]] std::pair<double, double> Edge(const Left& left)
{
    double l = 0.0;
    double r = 1.0;
    for (double m = 0.5 * l + 0.5 * r; l < m && m < r; m = 0.5 * l + 0.5 * r) {
        (left(m) ? l : r) = m;
    }
    return {l, r};
}

/// The piece's parameters whose x may lie in `x`, rounded outward (every s
/// outside it has x proved outside `x`), by bisection on the monotone x;
/// empty when none.
[[nodiscard]] Interval Preimage(const Graph& g, const Frame& f, const Interval& x)
{
    const Interval target = Oriented(g, x);
    const auto below = [&](const double s) {
        return RisingAt(g, f, s).GetUpper() < target.GetLower();
    };
    const auto notAbove = [&](const double s) {
        return !(RisingAt(g, f, s).GetLower() > target.GetUpper());
    };
    if (below(1.0) || !notAbove(0.0)) {
        return Interval::Empty();
    }
    // The last s proved below, and the first proved above.
    const double lo = below(0.0) ? Edge(below).first : 0.0;
    const double hi = notAbove(1.0) ? 1.0 : Edge(notAbove).second;
    return lo <= hi ? Interval{lo, hi} : Interval::Empty();
}

/// The piece's parameters whose x is proved to lie in `x`, rounded inward;
/// empty when rounding leaves none.
[[nodiscard]] Interval InnerPreimage(const Graph& g, const Frame& f, const Interval& x)
{
    const Interval target = Oriented(g, x);
    const auto before = [&](const double s) {
        return RisingAt(g, f, s).GetLower() < target.GetLower();
    };
    const auto within = [&](const double s) {
        return RisingAt(g, f, s).GetUpper() <= target.GetUpper();
    };
    if (before(1.0) || !within(0.0)) {
        return Interval::Empty();
    }
    // The first s proved at or past the lower end (x being monotone, every
    // later one is), and the last proved at or before the upper end.
    const double lo = before(0.0) ? Edge(before).second : 0.0;
    const double hi = within(1.0) ? 1.0 : Edge(within).first;
    const Interval ends = Hull(RisingAt(g, f, lo), RisingAt(g, f, hi));
    return lo <= hi && Contains(target, ends) ? Interval{lo, hi} : Interval::Empty();
}

/// A pair of pieces, both graphs over the frame's x, whose meetings are the
/// roots of h(x) = y_a(x) - y_b(x) over the x both reach.
struct Touching {
    PairSearch& Search;
    const Arc& A;
    const Arc& B;
    Frame F;
    Graph GA;
    Graph GB;
    std::vector<Meeting> Found;
};

/// h's derivatives over the x-range x: each graph's over the parameters
/// that reach it, subtracted; nothing when x is beyond either piece or a
/// graph cannot be enclosed.
[[nodiscard]] std::optional<Jet> Difference(const Touching& t, const Interval& x)
{
    const Interval sa = Preimage(t.GA, t.F, x);
    const Interval sb = Preimage(t.GB, t.F, x);
    if (sa.IsEmpty() || sb.IsEmpty()) {
        return std::nullopt;
    }
    const std::optional<Jet> ga = GraphOver(t.A, t.F, sa);
    const std::optional<Jet> gb = GraphOver(t.B, t.F, sb);
    if (!ga || !gb) {
        return std::nullopt;
    }
    Jet h;
    for (std::size_t j = 0; j <= JetOrder; ++j) {
        h[j] = (*ga)[j] - (*gb)[j];
    }
    return h;
}

/// The sign of h^(j) at the point x, Unknown when it is not proved.
[[nodiscard]] Strict SignAt(const Touching& t, const double x, const std::size_t j)
{
    const std::optional<Jet> h = Difference(t, Interval{x});
    return h ? StrictSign((*h)[j]) : Strict::Unknown;
}

/// The meeting the x-range [x0, x1] stands for, in both curves'
/// parameters.
[[nodiscard]] Meeting MeetingOver(const Touching& t, const double x0, const double x1,
                                  const ContactType contact)
{
    const Interval x{x0, x1};
    const Interval sa = Preimage(t.GA, t.F, x);
    const Interval sb = Preimage(t.GB, t.F, x);
    return {OnCurve(t.Search.A[t.A.Index], InSegment(t.A, sa)),
            OnCurve(t.Search.B[t.B.Index], InSegment(t.B, sb)),
            false,
            Interval::Empty(),
            Interval::Empty(),
            contact};
}

/// h strictly monotone over [x0, x1] (h' of one sign): no root, one
/// crossing proved by opposite signs at the ends and bisected (decided, the
/// only meeting of the pair over the range), or a crossing within rounding
/// of an end (undecided).
[[nodiscard]] Result<void> Monotonic(Touching& t, const double x0, const double x1)
{
    const Strict start = SignAt(t, x0, 0);
    const Strict end = SignAt(t, x1, 0);
    if (start != Strict::Unknown && start == end) {
        return {};
    }
    if (start == Strict::Unknown || end == Strict::Unknown) {
        t.Found.push_back(MeetingOver(t, x0, x1, ContactType::Crossing));
        return {};
    }
    double l = x0;
    double r = x1;
    for (double m = 0.5 * l + 0.5 * r; l < m && m < r; m = 0.5 * l + 0.5 * r) {
        const Strict s = SignAt(t, m, 0);
        if (s == Strict::Unknown) {
            break;
        }
        (s == start ? l : r) = m;
    }
    Meeting root = MeetingOver(t, l, r, ContactType::Crossing);
    const Interval range{x0, x1};
    const Interval proofA = InnerPreimage(t.GA, t.F, range);
    if (!proofA.IsEmpty()) {
        const Segment& segmentA = t.Search.A[t.A.Index];
        const Segment& segmentB = t.Search.B[t.B.Index];
        root.Decided = true;
        root.ProofA =
            Intersection(InwardOnCurve(segmentA, t.A, proofA), Interval{segmentA.T0, segmentA.T1});
        root.ProofB = Intersection(InwardOnCurve(segmentB, t.B, Interval{0.0, 1.0}),
                                   Interval{segmentB.T0, segmentB.T1});
    }
    t.Found.push_back(root);
    return {};
}

[[nodiscard]] Result<void> Follow(Touching& t, double x0, double x1, std::size_t k, int depth);

/// The x-range [x0, x1] with h^(k), k >= 1, of one strict sign over it (on
/// it or on an enclosing range): Bernstein.hpp's Convex and the plane
/// handler's OfOrder in the graphs' common x. h proved of one sign: no
/// meeting. A lower derivative of one sign settles it at that order; h'
/// of one sign makes h monotone.
[[nodiscard]] Result<void> Chain(Touching& t, const double x0, const double x1, const std::size_t k,
                                 const int depth)
{
    if (++t.Search.Pairs > PairBudget) {
        return std::unexpected(Unresolved(t.Search, t.A, t.B, "the subdivision budget is spent"));
    }
    const std::optional<Jet> h = Difference(t, Interval{x0, x1});
    if (!h) {
        return {}; // Beyond either piece: no meeting of this pair.
    }
    if (StrictSign((*h)[0]) != Strict::Unknown) {
        return {};
    }
    for (std::size_t j = 1; j < k; ++j) {
        if (StrictSign((*h)[j]) != Strict::Unknown) {
            return j == 1 ? Monotonic(t, x0, x1) : Chain(t, x0, x1, j, depth);
        }
    }
    return k == 1 ? Monotonic(t, x0, x1) : Follow(t, x0, x1, k, depth);
}

/// h^(k-1), strictly monotone, at the range's ends: one strict sign there
/// is that sign throughout; else the range holding its single zero is
/// halved. Both ends within rounding of zero, and no lower derivative of
/// one sign, or MaxDepth: a contact of order up to k to rounding, h having
/// at most k roots there (Rolle), kept as one undecided tangential meeting
/// bracketed by the range, never dropped.
[[nodiscard]] Result<void> Follow(Touching& t, const double x0, const double x1,
                                  const std::size_t k, const int depth)
{
    const Strict start = SignAt(t, x0, k - 1);
    const Strict end = SignAt(t, x1, k - 1);
    if (start != Strict::Unknown && start == end) {
        return Chain(t, x0, x1, k - 1, depth);
    }
    const double m = 0.5 * x0 + 0.5 * x1;
    const bool flat = start == Strict::Unknown && end == Strict::Unknown;
    if (flat || depth == MaxDepth || !(x0 < m && m < x1)) {
        t.Found.push_back(MeetingOver(t, x0, x1, ContactType::Tangential));
        return {};
    }
    STEREON_TRY_VOID(Chain(t, x0, m, k, depth + 1));
    return Chain(t, m, x1, k, depth + 1);
}

/// The piece as a graph over the frame's x, if x' is proved of one sign
/// over it.
[[nodiscard]] std::optional<Graph> AsGraph(const Arc& arc, const Frame& f)
{
    const std::optional<std::pair<Jet, Jet>> c = CoordinatesOver(arc, f, Interval{0.0, 1.0});
    if (!c) {
        return std::nullopt;
    }
    const Strict sense = StrictSign(c->first[1]);
    if (sense == Strict::Unknown) {
        return std::nullopt;
    }
    return Graph{&arc, sense == Strict::Positive ? 1.0 : -1.0};
}

/// Marks every meeting whose brackets meet, on both curves, those of an
/// undecided one already joined to a touch; whether any was.
[[nodiscard]] bool JoinOnce(const std::vector<Meeting>& found, std::vector<bool>& joined)
{
    bool grown = false;
    for (std::size_t i = 0; i < found.size(); ++i) {
        for (std::size_t j = 0; j < found.size() && !joined[i]; ++j) {
            if (joined[j] && !found[j].Decided && Intersects(found[i].OnA, found[j].OnA) &&
                Intersects(found[i].OnB, found[j].OnB)) {
                joined[i] = true;
                grown = true;
            }
        }
    }
    return grown;
}

/// Whether the graphs' meetings stand: every undecided crossing (one within
/// rounding of a range's end, where the halving in x put it) is joined to a
/// touch by brackets that meet on both curves, so that Distinct merges it
/// into the touch's rounding. Otherwise the pair is halved as usual and
/// Krawczyk's test settles that crossing.
[[nodiscard]] bool Settled(const std::vector<Meeting>& found)
{
    std::vector<bool> joined(found.size());
    for (std::size_t i = 0; i < found.size(); ++i) {
        joined[i] = found[i].Decided || found[i].Contact == ContactType::Tangential;
    }
    while (JoinOnce(found, joined)) {
    }
    return std::ranges::all_of(joined, [](const bool b) { return b; });
}

/// Touches' pair, both pieces graphs: settled over the x both reach when h'
/// is not of one sign and some h^(k), k >= 2, is; nothing otherwise, or
/// when the graphs leave a crossing to Krawczyk's test.
[[nodiscard]] std::optional<Result<void>> AsGraphs(Touching t)
{
    const Interval xa = Hull(XAt(t.A, t.F, 0.0), XAt(t.A, t.F, 1.0));
    const Interval xb = Hull(XAt(t.B, t.F, 0.0), XAt(t.B, t.F, 1.0));
    const double x0 = std::max(xa.GetLower(), xb.GetLower());
    const double x1 = std::min(xa.GetUpper(), xb.GetUpper());
    if (!(x0 < x1)) {
        return std::nullopt;
    }
    const std::optional<Jet> h = Difference(t, Interval{x0, x1});
    // h' of one sign: no touch here, at most one transverse crossing, which
    // Krawczyk's test settles as the pieces are halved.
    if (!h || StrictSign((*h)[1]) != Strict::Unknown) {
        return std::nullopt;
    }
    for (std::size_t k = 2; k <= JetOrder; ++k) {
        if (StrictSign((*h)[k]) != Strict::Unknown) {
            STEREON_TRY_VOID(Chain(t, x0, x1, k, 0));
            if (!Settled(t.Found)) {
                return std::nullopt;
            }
            t.Search.Roots.insert(t.Search.Roots.end(), t.Found.begin(), t.Found.end());
            return Result<void>{};
        }
    }
    return std::nullopt;
}

/// A pair neither apart nor a single crossing, settled as graphs when both
/// pieces are graphs over the chord of the first and some derivative of h
/// of order two or more has one sign over their common x: a tangency of two
/// B-splines, of even order or odd, is then followed down one chain of
/// ranges and reported bracketed, Tangential and undecided, never dropped.
/// Nothing when it does not apply.
[[nodiscard]] std::optional<Result<void>> Touches(PairSearch& search, const Arc& a, const Arc& b)
{
    const Point2 origin = Nominal(a.Poles.front());
    const Vec2 u = Nominal(a.Poles.back()) - origin;
    if (!(Dot(u, u) > 0.0) || !std::isfinite(Dot(u, u))) {
        return std::nullopt;
    }
    const Frame f{origin, u, Perpendicular(u)};
    const std::optional<Graph> ga = AsGraph(a, f);
    const std::optional<Graph> gb = AsGraph(b, f);
    if (!ga || !gb) {
        return std::nullopt;
    }
    return AsGraphs(Touching{search, a, b, f, *ga, *gb, {}});
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
    if (std::optional<Result<void>> touch = Touches(search, a, b)) {
        return std::move(*touch);
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
        const bool touch =
            same->Contact == ContactType::Tangential || r.Contact == ContactType::Tangential;
        *same = {Hull(same->OnA, r.OnA),
                 Hull(same->OnB, r.OnB),
                 false,
                 Interval::Empty(),
                 Interval::Empty(),
                 touch ? ContactType::Tangential : ContactType::Crossing};
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

// --- One curve in two knot vectors ---------------------------------------

/// a + b, when it is exact in doubles (TwoSum's error term zero).
[[nodiscard]] std::optional<double> ExactSum(const double a, const double b)
{
    const double s = a + b;
    const double bb = s - a;
    if (!std::isfinite(s) || (a - (s - bb)) + (b - bb) != 0.0) {
        return std::nullopt;
    }
    return s;
}

/// a b, when it is exact in doubles (the fused residual zero).
[[nodiscard]] std::optional<double> ExactProduct(const double a, const double b)
{
    const double p = a * b;
    if (!std::isfinite(p) || (p != 0.0 && !std::isnormal(p)) || std::fma(a, b, -p) != 0.0) {
        return std::nullopt;
    }
    return p;
}

/// a / b, when it is exact in doubles.
[[nodiscard]] std::optional<double> ExactQuotient(const double a, const double b)
{
    if (b == 0.0) {
        return std::nullopt;
    }
    const double q = a / b;
    if (!std::isfinite(q) || (q != 0.0 && !std::isnormal(q)) || std::fma(q, b, -a) != 0.0) {
        return std::nullopt;
    }
    return q;
}

/// A homogeneous pole (w P, w) in doubles, held exactly.
using Exact3 = std::array<double, 3>;

/// (wa a + wb b) / w, when every step is exact: Boehm's blend with its ratio
/// never formed, so a pole that is exact is found even where the ratio
/// (2/3, say) is not a double.
[[nodiscard]] std::optional<Exact3> ExactBlend(const Exact3& a, const Exact3& b, const double wa,
                                               const double wb, const double w)
{
    Exact3 c{};
    for (std::size_t i = 0; i < 3; ++i) {
        const std::optional<double> l = ExactProduct(wa, a[i]);
        const std::optional<double> r = ExactProduct(wb, b[i]);
        const std::optional<double> s = l && r ? ExactSum(*l, *r) : std::nullopt;
        const std::optional<double> q = s ? ExactQuotient(*s, w) : std::nullopt;
        if (!q) {
            return std::nullopt;
        }
        c[i] = *q;
    }
    return c;
}

/// Boehm's insertion of u, once, when it is exact in doubles (Bernstein.hpp's
/// InsertKnot, without rounding).
[[nodiscard]] bool InsertExactly(std::vector<double>& knots, std::vector<Exact3>& poles,
                                 const std::size_t p, const double u)
{
    const std::size_t n = poles.size() - 1;
    const auto above = static_cast<std::size_t>(std::ranges::upper_bound(knots, u) - knots.begin());
    const std::size_t k = std::min(above - 1, n);
    std::vector<Exact3> refined(poles.begin(),
                                poles.begin() + static_cast<std::ptrdiff_t>(k - p + 1));
    for (std::size_t i = k - p + 1; i <= k; ++i) {
        // (t_i+p - u) P_i-1 + (u - t_i) P_i over t_i+p - t_i.
        const std::optional<double> right = ExactSum(u, -knots[i]);
        const std::optional<double> left = ExactSum(knots[i + p], -u);
        const std::optional<double> width = ExactSum(knots[i + p], -knots[i]);
        const std::optional<Exact3> pole =
            right && left && width ? ExactBlend(poles[i - 1], poles[i], *left, *right, *width)
                                   : std::nullopt;
        if (!pole) {
            return false;
        }
        refined.push_back(*pole);
    }
    refined.insert(refined.end(), poles.begin() + static_cast<std::ptrdiff_t>(k), poles.end());
    knots.insert(knots.begin() + static_cast<std::ptrdiff_t>(k + 1), u);
    poles = std::move(refined);
    return true;
}

/// The spline's homogeneous poles, when w P is exact for every pole.
[[nodiscard]] std::optional<std::vector<Exact3>> ExactPoles(const BSplineCurve2d& spline)
{
    std::vector<Exact3> poles;
    for (std::size_t i = 0; i < spline.Poles.size(); ++i) {
        const double w = spline.Weights.empty() ? 1.0 : spline.Weights[i];
        const std::optional<double> x = ExactProduct(w, spline.Poles[i].X);
        const std::optional<double> y = ExactProduct(w, spline.Poles[i].Y);
        if (!x || !y) {
            return std::nullopt;
        }
        poles.push_back({*x, *y, w});
    }
    return poles;
}

/// Raises every knot of `knots` to its multiplicity in `other`, exactly;
/// false when a knot lies outside the domain or an insertion would round.
[[nodiscard]] bool RaisedTo(std::vector<double>& knots, std::vector<Exact3>& poles,
                            const std::size_t p, const std::vector<double>& other,
                            const Interval& domain)
{
    for (const double u : other) {
        while (std::ranges::count(knots, u) < std::ranges::count(other, u)) {
            if (!Contains(domain, u) || !InsertExactly(knots, poles, p, u)) {
                return false;
            }
        }
    }
    return true;
}

/// Whether the splines are one curve on one parameter: of one degree and
/// domain, and, each refined to the other's knots by Boehm's insertion
/// carried out exactly in doubles (every sum, product and ratio checked
/// error-free), with equal knot vectors and equal homogeneous poles. Knot
/// insertion keeps the curve and its parameter, so the refined splines are
/// the inputs, and their overlap is the identity over the common stretch.
/// When any step would round it is not proved: false, and the pair is
/// searched (and, running along itself, fails closed).
[[nodiscard]] bool SameCurve(const BSplineCurve2d& a, const BSplineCurve2d& b)
{
    const Interval domain = Domain(a);
    if (a.Degree != b.Degree || domain.GetLower() != Domain(b).GetLower() ||
        domain.GetUpper() != Domain(b).GetUpper()) {
        return false;
    }
    const auto p = static_cast<std::size_t>(a.Degree);
    std::vector<double> ka = a.Knots;
    std::vector<double> kb = b.Knots;
    std::optional<std::vector<Exact3>> pa = ExactPoles(a);
    std::optional<std::vector<Exact3>> pb = ExactPoles(b);
    return pa && pb && RaisedTo(ka, *pa, p, b.Knots, domain) &&
           RaisedTo(kb, *pb, p, a.Knots, domain) && ka == kb && *pa == *pb;
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

/// Whether a - b is exact in doubles: Knuth's TwoSum of a and -b, whose
/// rounding error is itself exact, has none.
[[nodiscard]] bool ExactDifference(const double a, const double b)
{
    const double s = a - b;
    const double bb = s - a;
    const double error = (a - (s - bb)) + (-b - bb);
    return std::isfinite(s) && error == 0.0;
}

/// The spline with every pole moved by -c, if every coordinate's difference
/// is exact: then it is the same curve, translated, on the same parameter.
[[nodiscard]] std::optional<BSplineCurve2d> MovedExactly(const BSplineCurve2d& spline,
                                                         const Point2& c)
{
    BSplineCurve2d moved = spline;
    for (Point2& pole : moved.Poles) {
        if (!ExactDifference(pole.X, c.X) || !ExactDifference(pole.Y, c.Y)) {
            return std::nullopt;
        }
        pole = Point2{pole.X - c.X, pole.Y - c.Y};
    }
    return moved;
}

/// Both splines moved exactly by the middle of their poles' box, or nothing
/// when a move would round. Translating both by one vector keeps F and both
/// parameters, and the enclosures' rounding, relative to the coordinates'
/// size, then measures against the curves' size: a pair far from the
/// origin is searched as it is near it.
[[nodiscard]] std::optional<std::pair<BSplineCurve2d, BSplineCurve2d>>
Centred(const BSplineCurve2d& a, const BSplineCurve2d& b)
{
    double lowX = std::numeric_limits<double>::infinity();
    double lowY = lowX;
    double highX = -lowX;
    double highY = -lowX;
    for (const BSplineCurve2d* spline : {&a, &b}) {
        for (const Point2& pole : spline->Poles) {
            lowX = std::min(lowX, pole.X);
            lowY = std::min(lowY, pole.Y);
            highX = std::max(highX, pole.X);
            highY = std::max(highY, pole.Y);
        }
    }
    const Point2 c{0.5 * lowX + 0.5 * highX, 0.5 * lowY + 0.5 * highY};
    std::optional<BSplineCurve2d> movedA = MovedExactly(a, c);
    std::optional<BSplineCurve2d> movedB = MovedExactly(b, c);
    if (!movedA || !movedB) {
        return std::nullopt;
    }
    return std::pair{std::move(*movedA), std::move(*movedB)};
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

/// The splines' segments, and every pair of them searched. The search runs
/// on the pair moved to the middle of its poles when that move is exact:
/// the same F, enclosed at the size of the curves rather than of their
/// coordinates.
[[nodiscard]] Result<void> SearchPair(PairSearch& search, const BSplineCurve2d& a,
                                      const Interval& domainA, const BSplineCurve2d& b,
                                      const Interval& domainB)
{
    const std::optional<std::pair<BSplineCurve2d, BSplineCurve2d>> centred = Centred(a, b);
    STEREON_TRY(std::vector<Segment> segmentsA, SegmentsOf(centred ? centred->first : a));
    STEREON_TRY(std::vector<Segment> segmentsB, SegmentsOf(centred ? centred->second : b));
    search.AllA = segmentsA;
    search.AllB = segmentsB;
    search.A = Reaching(std::move(segmentsA), domainA);
    search.B = Reaching(std::move(segmentsB), domainB);
    return SearchAll(search);
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
    if (*sa == *sb || SameCurve(*sa, *sb)) {
        return Identical(*sa, domainA, domainB);
    }
    STEREON_TRY_VOID(SearchPair(search, *sa, domainA, *sb, domainB));
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
            result.Points.push_back({*ta, *tb, Eval(*sa, *ta), root.Contact, decided});
        }
    }
    Order(result);
    return result;
}

} // namespace Stereon::Detail
