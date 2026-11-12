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
#include <initializer_list>
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
#include "BezierArc.hpp"
#include "ConicArc.hpp"
#include "CurvePieces.hpp"

/// Two curves in space: B-splines, polynomial or rational, and lines and
/// conics as enclosed Bezier pieces of the same search, every pair but two
/// lines (LineLine.cpp). Two curves in space
/// meet only exceptionally: doubles placed to make them meet leave them a
/// rounding apart. So a meeting is defined by distance (ADR-0001: points
/// closer than the resolution cannot be told apart): a critical point of
/// the squared distance
///
///   G(s, t) = ( (a(s) - b(t)) . a'(s),  -(a(s) - b(t)) . b'(t) ) = 0,
///
/// the gradient of |a - b|² / 2, at which |a - b| is within the Context's
/// resolution. Every exact meeting is such a point, and so is the closest
/// approach of two curves that pass within the resolution.
///
/// The search runs over pairs of pieces of the curves' Bezier segments
/// (Bernstein.hpp encloses them; BezierArc.hpp halves them), each pair
///
/// - apart: the boxes of the pieces' poles, or their hulls along the common
///   normal of the pieces' chords (a fat slab, after Sederberg and
///   Nishita's fat line), are more than the resolution apart;
/// - without a critical point, or with exactly one, by Krawczyk's operator
///   on G over the pair's parameter box inflated by an eighth of its width,
///   the Jacobian
///
///     [ a'.a' + (a - b).a''      -a'.b'            ]
///     [ -a'.b'                    b'.b' - (a - b).b'' ]
///
///   enclosed by interval de Casteljau on the poles and their first and
///   second hodographs. Near a meeting it is the Gram matrix of the two
///   tangents, regular unless they are parallel;
/// - or halved.
///
/// A single critical point is bracketed to rounding, and its distance
/// enclosed over the bracket: within the resolution, a meeting; beyond it,
/// a near miss, not reported; across it, undecidable at the resolution, and
/// the intersection fails. Critical points outside the requested domains
/// are not meetings: there the distance's minimum over the domains lies on
/// their boundary, an end of one curve against the other's stretch, found
/// by the same search on the end's distance alone, or two ends. So two
/// edges sharing a vertex meet at it even when their ends are a rounding
/// apart, and the curves meet as two lines do (LineLine.cpp). Critical
/// points proved twice, from overlapping boxes, are one, and so are
/// meetings within the resolution of each other on both curves. A pair that is none of these
/// when the pieces can no longer be halved, or a search past the budget,
/// fails with IntersectionFailed: tangent curves, an overlap, or a meeting
/// too nearly tangential to separate is never dropped and never reported
/// unproved (ADR-0029). Two identical curves coincide: their overlap is
/// the requested domains' common stretch on the identity map.
///
/// A circle, an ellipse or a hyperbola is cut into rational quadratic arcs
/// (ConicArc.hpp); a meeting found on an arc's linear parameter is mapped
/// back to the conic's. An unbounded domain is clipped to where its curve
/// can come near the other, bounded curve; two unbounded domains are
/// NotSupported.

namespace Stereon::Detail {

namespace {

/// The pairs of pieces the search may visit.
constexpr std::size_t PairBudget = std::size_t{1} << 16;

/// The fraction of a box's width by which Krawczyk's box is inflated.
constexpr double Inflation = 0.125;

/// Krawczyk's refinement: at most this many steps, each shrinking the box.
constexpr int RefineSteps = 64;

using Vec3I = std::array<Interval, 3>;

[[nodiscard]] Vec3I Minus(const Vec3I& a, const Vec3I& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[nodiscard]] Interval DotI(const Vec3I& a, const Vec3I& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// A piece's point and first two derivatives over a parameter interval,
/// enclosed: with C = X / W, C' = (X' - C W') / W and C'' = (X'' - 2 C' W'
/// - C W'') / W. Nothing when W is not proved positive.
struct Jet3 {
    Vec3I Point;
    Vec3I First;
    Vec3I Second;
};

/// A piece's poles and their first and second hodographs.
struct Hodographs {
    std::vector<HPoint> Poles;
    std::vector<HPoint> First;
    std::vector<HPoint> Second;
};

[[nodiscard]] Hodographs HodographsOf(const Arc& arc)
{
    std::vector<HPoint> first = Hodograph(arc.Poles);
    std::vector<HPoint> second = Hodograph(first);
    return {arc.Poles, std::move(first), std::move(second)};
}

[[nodiscard]] std::optional<Jet3> JetOver(const Hodographs& h, const Interval& s)
{
    const HPoint p = EvalH(h.Poles, s);
    if (!(p.W.GetLower() > 0.0)) {
        return std::nullopt;
    }
    const HPoint d = EvalH(h.First, s);
    const HPoint dd = EvalH(h.Second, s);
    const Vec3I c{p.X / p.W, p.Y / p.W, p.Z / p.W};
    const Vec3I c1{(d.X - c[0] * d.W) / p.W, (d.Y - c[1] * d.W) / p.W, (d.Z - c[2] * d.W) / p.W};
    const Interval two{2.0};
    const auto second = [&](const Interval& x2, const std::size_t k) {
        return (x2 - two * c1[k] * d.W - c[k] * dd.W) / p.W;
    };
    return Jet3{c, c1, {second(dd.X, 0), second(dd.Y, 1), second(dd.Z, 2)}};
}

// --- Separation ----------------------------------------------------------

/// The box of the piece's Euclidean poles, which holds the piece.
[[nodiscard]] std::optional<Vec3I> BoxOf(const Arc& arc)
{
    Vec3I box{Interval::Empty(), Interval::Empty(), Interval::Empty()};
    for (const HPoint& h : arc.Poles) {
        if (!(h.W.GetLower() > 0.0)) {
            return std::nullopt;
        }
        box = {Hull(box[0], h.X / h.W), Hull(box[1], h.Y / h.W), Hull(box[2], h.Z / h.W)};
    }
    return box;
}

/// The hull of n . (P_i - q) over the piece's poles: it holds n . (C - q)
/// over the piece, whatever n and q are.
[[nodiscard]] Interval ProjectionsOf(const Arc& arc, const Vec3& n, const Point3& q)
{
    Interval hull = Interval::Empty();
    for (const HPoint& h : arc.Poles) {
        if (!(h.W.GetLower() > 0.0)) {
            return Interval::Everything();
        }
        hull = Hull(hull, (h.X / h.W - Interval{q.X}) * Interval{n.X} +
                              (h.Y / h.W - Interval{q.Y}) * Interval{n.Y} +
                              (h.Z / h.W - Interval{q.Z}) * Interval{n.Z});
    }
    return hull;
}

/// Whether two hulls are more than `gap` apart.
[[nodiscard]] bool FarApart(const Interval& a, const Interval& b, const double gap)
{
    return (Interval{a.GetUpper()} + Interval{gap}).GetUpper() < b.GetLower() ||
           (Interval{b.GetUpper()} + Interval{gap}).GetUpper() < a.GetLower();
}

/// Whether every point of one piece is more than `gap` from every point of
/// the other: along a coordinate axis, or along the common normal of their
/// chords (unit, so the hulls' separation is a distance).
[[nodiscard]] bool Apart(const Arc& a, const Arc& b, const double gap)
{
    const std::optional<Vec3I> boxA = BoxOf(a);
    const std::optional<Vec3I> boxB = BoxOf(b);
    if (boxA && boxB) {
        for (std::size_t k = 0; k < 3; ++k) {
            if (FarApart((*boxA)[k], (*boxB)[k], gap)) {
                return true;
            }
        }
    }
    const Point3 a0 = Nominal3(a.Poles.front());
    const Point3 b0 = Nominal3(b.Poles.front());
    const Vec3 n = Cross(Nominal3(a.Poles.back()) - a0, Nominal3(b.Poles.back()) - b0);
    const double length = Length(n);
    if (!(length > 0.0) || !std::isfinite(length)) {
        return false;
    }
    // A rounded unit normal: the hulls are of n . (P - q) with this n, so
    // they are distances only to within its rounding; the gap is widened by
    // that relative error, far below the resolution's own scale.
    const Vec3 unit = n / length;
    const double slack = 1.0 + 0x1p-40;
    return FarApart(ProjectionsOf(a, unit, a0), ProjectionsOf(b, unit, a0), gap * slack);
}

// --- Krawczyk's operator on the distance's gradient ----------------------

[[nodiscard]] std::optional<ParameterBox> Krawczyk(const Hodographs& a, const Hodographs& b,
                                                   const ParameterBox& x)
{
    const auto am = JetOver(a, Interval{Midpoint(x.S)});
    const auto bm = JetOver(b, Interval{Midpoint(x.T)});
    const auto ax = JetOver(a, x.S);
    const auto bx = JetOver(b, x.T);
    if (!am || !bm || !ax || !bx) {
        return std::nullopt;
    }
    const Vec3I gap = Minus(am->Point, bm->Point);
    const std::array<Interval, 2> g{DotI(gap, am->First), -DotI(gap, bm->First)};
    const Vec3I spread = Minus(ax->Point, bx->Point);
    const Interval cross = -DotI(ax->First, bx->First);
    const std::array<Interval, 4> j{DotI(ax->First, ax->First) + DotI(spread, ax->Second), cross,
                                    cross, DotI(bx->First, bx->First) - DotI(spread, bx->Second)};
    return KrawczykStep(x, g, j);
}

/// The single critical point of a box Krawczyk's test proved to hold one,
/// bracketed by iterating X <- K(X) ∩ X, which keeps it, while the box
/// shrinks.
[[nodiscard]] ParameterBox Refined(const Hodographs& a, const Hodographs& b, ParameterBox x)
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

/// |a(S) - b(T)| enclosed, or nothing when a weight is not proved positive.
[[nodiscard]] std::optional<Interval> DistanceOver(const Hodographs& a, const Interval& s,
                                                   const Hodographs& b, const Interval& t)
{
    const auto ja = JetOver(a, s);
    const auto jb = JetOver(b, t);
    if (!ja || !jb) {
        return std::nullopt;
    }
    const Vec3I d = Minus(ja->Point, jb->Point);
    const Interval squared = Squared(d[0]) + Squared(d[1]) + Squared(d[2]);
    return Sqrt(Interval{std::max(squared.GetLower(), 0.0), squared.GetUpper()});
}

// --- The search ----------------------------------------------------------

/// A meeting found: its brackets in both curves' parameters.
struct Meeting {
    Interval OnA;
    Interval OnB;
};

/// Both splines' segments, the requested domains, the resolution, the
/// meetings found and the count of pairs.
struct PairSearch {
    std::vector<Segment> A;
    std::vector<Segment> B;
    Interval DomainA;
    Interval DomainB;
    double Resolution = 0.0;
    std::vector<Meeting> Found;
    std::size_t Pairs = 0;
    bool ClosedA = false;
    bool ClosedB = false;
    std::string Handler;
};

[[nodiscard]] KernelError Unresolved(const PairSearch& search, const Arc& a, const Arc& b,
                                     const std::string_view why)
{
    const double ta = Midpoint(OnCurve(search.A[a.Index], InSegment(a, Interval{0.5})));
    const double tb = Midpoint(OnCurve(search.B[b.Index], InSegment(b, Interval{0.5})));
    return Fail(ErrorCode::IntersectionFailed,
                search.Handler + ": " + std::string(why) + " near (" + Number(ta) + ", " +
                    Number(tb) +
                    ") (tangent curves, an overlap or a nearly tangential meeting; ADR-0010: an "
                    "unresolved pair fails, it is never dropped)")
        .error();
}

/// A bracket in a curve's parameter, kept for a segment whose stretch in
/// the requested domain is `stretch`: itself when it reaches the stretch
/// (a bracket across a knot is proved from both segments and merged);
/// nothing when it lies past an end of the stretch: inside the domain the
/// neighbouring segment proves it, and past the domain's end the minimum
/// lies on the boundary, where SearchEnds finds it.
[[nodiscard]] std::optional<Interval> Kept(const Interval& bracket, const Interval& stretch,
                                           const Interval& domain)
{
    if (Intersects(bracket, stretch)) {
        return Intersection(bracket, domain);
    }
    return std::nullopt;
}

/// The segment-local parameter of a curve parameter.
[[nodiscard]] Interval LocalOf(const Segment& segment, const Interval& t)
{
    return (t - Interval{segment.T0}) / (Interval{segment.T1} - Interval{segment.T0});
}

/// Records the critical point of the pair's inflated box when it is a
/// meeting; fails when its distance straddles the resolution.
[[nodiscard]] Result<void> Record(PairSearch& search, const Arc& a, const Hodographs& ha,
                                  const Arc& b, const Hodographs& hb, const ParameterBox& x)
{
    const ParameterBox root = Refined(ha, hb, x);
    const Segment& segA = search.A[a.Index];
    const Segment& segB = search.B[b.Index];
    const std::optional<Interval> onA =
        Kept(OnCurve(segA, InSegment(a, root.S)),
             Intersection(Interval{segA.T0, segA.T1}, search.DomainA), search.DomainA);
    const std::optional<Interval> onB =
        Kept(OnCurve(segB, InSegment(b, root.T)),
             Intersection(Interval{segB.T0, segB.T1}, search.DomainB), search.DomainB);
    if (!onA || !onB) {
        return {};
    }
    // The distance at the kept brackets, over the whole segments' poles:
    // a bracket clipped to a domain's end may lie outside the piece.
    const std::optional<Interval> distance =
        DistanceOver(HodographsOf(Arc{a.Index, 0.0, 1.0, segA.Poles}), LocalOf(segA, *onA),
                     HodographsOf(Arc{b.Index, 0.0, 1.0, segB.Poles}), LocalOf(segB, *onB));
    if (!distance ||
        (distance->GetLower() <= search.Resolution && search.Resolution < distance->GetUpper())) {
        return std::unexpected(
            Unresolved(search, a, b, "a closest approach within rounding of the resolution"));
    }
    if (distance->GetUpper() <= search.Resolution) {
        search.Found.push_back({*onA, *onB});
    }
    return {};
}

/// A stationary distance requires both gradient components to contain zero.
/// Test the actual piece domains before inverting the near-singular Gram
/// matrix. Boundary minima are still handled independently by SearchEnds.
[[nodiscard]] bool NoStationaryPoint(const Hodographs& a, const Hodographs& b)
{
    const Interval domain{0.0, 1.0};
    const auto ja = JetOver(a, domain);
    const auto jb = JetOver(b, domain);
    if (!ja || !jb) {
        return false;
    }
    const Vec3I difference = Minus(ja->Point, jb->Point);
    return !Contains(DotI(difference, ja->First), 0.0) ||
           !Contains(DotI(difference, jb->First), 0.0);
}

/// One pair of pieces, apart, settled by Krawczyk's test, or halved.
[[nodiscard]] Result<void> VisitPair(PairSearch& search, const Arc& a, const Arc& b,
                                     const int depth)
{
    if (++search.Pairs > PairBudget) {
        return std::unexpected(Unresolved(search, a, b, "the subdivision budget is spent"));
    }
    if (Apart(a, b, search.Resolution)) {
        return {};
    }
    const Hodographs ha = HodographsOf(a);
    const Hodographs hb = HodographsOf(b);
    if (NoStationaryPoint(ha, hb)) {
        return {};
    }
    const ParameterBox x{Interval{-Inflation, 1.0 + Inflation},
                         Interval{-Inflation, 1.0 + Inflation}};
    const Test test = Classify(Krawczyk(ha, hb, x), x);
    if (test == Test::NoRoot) {
        return {};
    }
    if (test == Test::OneRoot) {
        return Record(search, a, ha, b, hb, x);
    }
    if (depth == MaxDepth) {
        return std::unexpected(
            Unresolved(search, a, b, "a pair of pieces neither apart nor a single meeting"));
    }
    const std::pair<Arc, Arc> halvesA = Halves(a);
    const std::pair<Arc, Arc> halvesB = Halves(b);
    for (const Arc* pa : {&halvesA.first, &halvesA.second}) {
        for (const Arc* pb : {&halvesB.first, &halvesB.second}) {
            STEREON_TRY_VOID(VisitPair(search, *pa, *pb, depth + 1));
        }
    }
    return {};
}

/// Meetings whose brackets meet on both curves are one, proved from boxes
/// that overlap.
[[nodiscard]] std::vector<Meeting> Distinct(std::vector<Meeting> found)
{
    std::ranges::stable_sort(found, [](const Meeting& p, const Meeting& q) {
        return p.OnA.GetLower() < q.OnA.GetLower();
    });
    std::vector<Meeting> distinct;
    for (const Meeting& m : found) {
        const auto same = std::ranges::find_if(distinct, [&](const Meeting& d) {
            return Intersects(d.OnA, m.OnA) && Intersects(d.OnB, m.OnB);
        });
        if (same == distinct.end()) {
            distinct.push_back(m);
        } else {
            same->OnA = Hull(same->OnA, m.OnA);
            same->OnB = Hull(same->OnB, m.OnB);
        }
    }
    return distinct;
}

// --- The domains' ends ---------------------------------------------------

/// An end of a requested domain: its parameter, the curve's point and
/// tangent there, enclosed, and its side: +1 for the upper end, -1 for the
/// lower. The tangent is the inside segment's, so a knot of reduced
/// continuity at the end does not matter.
struct DomainEnd {
    double T = 0.0;
    Vec3I Point;
    Vec3I Tangent;
    double Side = 1.0;
};

/// The ends of `domain` on a curve's segments, or nothing for an end whose
/// point a weight leaves unenclosed.
[[nodiscard]] std::vector<DomainEnd> EndsOf(const std::vector<Segment>& segments,
                                            const Interval& domain)
{
    std::vector<DomainEnd> ends;
    for (const double side : {-1.0, 1.0}) {
        const double t = side < 0.0 ? domain.GetLower() : domain.GetUpper();
        const auto inside = std::ranges::find_if(segments, [&](const Segment& s) {
            return side < 0.0 ? s.T0 <= t && t < s.T1 : s.T0 < t && t <= s.T1;
        });
        if (inside == segments.end()) {
            continue;
        }
        const Interval local =
            (Interval{t} - Interval{inside->T0}) / (Interval{inside->T1} - Interval{inside->T0});
        const auto jet = JetOver(HodographsOf(Arc{0, 0.0, 1.0, inside->Poles}), local);
        if (jet) {
            ends.push_back({t, jet->Point, jet->First, side});
        }
        if (domain.GetLower() == domain.GetUpper()) {
            break;
        }
    }
    return ends;
}

/// Whether the squared distance is proved to shrink from the end inward:
/// then the end is no minimum along its curve, and the meeting, if any,
/// lies inside, where the pair search finds it. Undecided keeps the end;
/// the merge below folds a meeting found both ways into one.
[[nodiscard]] bool ShrinksInward(const DomainEnd& end, const Vec3I& toOther)
{
    // d/dt |C - Q|² / 2 = (C - Q) . C'; the end is a minimum along C when
    // it decreases towards the end: (C - Q) . C' side < 0.
    const Interval slope = DotI(Minus(end.Point, toOther), end.Tangent) * Interval{end.Side};
    return slope.GetLower() >= 0.0;
}

/// One curve's end against the other's segments: the critical points of
/// the end's squared distance to the other curve, a root of
///
///   g(t) = (B(t) - P) . B'(t),   g'(t) = B' . B' + (B - P) . B'',
///
/// by Krawczyk's operator over halved pieces, as the pair search does.
struct EndSearch {
    const PairSearch* Pairs = nullptr;
    DomainEnd End;
    const std::vector<Segment>* Other = nullptr;
    Interval OtherDomain;
    bool Swapped = false;
    std::size_t Visited = 0;
    std::vector<Meeting> Found;
};

[[nodiscard]] KernelError EndUnresolved(const EndSearch& e, const std::string_view why)
{
    return Fail(ErrorCode::IntersectionFailed,
                e.Pairs->Handler + ": " + std::string(why) + " for the end at " + Number(e.End.T) +
                    " (a curve's end against a stretch of the other within the resolution; "
                    "ADR-0010: an unresolved pair fails, it is never dropped)")
        .error();
}

[[nodiscard]] std::optional<Interval> KrawczykEnd(const Hodographs& b, const Vec3I& p,
                                                  const Interval& x)
{
    const double m = Midpoint(x);
    const auto bm = JetOver(b, Interval{m});
    const auto bx = JetOver(b, x);
    if (!bm || !bx) {
        return std::nullopt;
    }
    const Interval g = DotI(Minus(bm->Point, p), bm->First);
    const Interval j = DotI(bx->First, bx->First) + DotI(Minus(bx->Point, p), bx->Second);
    const double c = Midpoint(j);
    if (!(std::abs(c) > 0.0) || !std::isfinite(c)) {
        return std::nullopt;
    }
    const Interval y{1.0 / c};
    return Interval{m} - y * g + (Interval{1.0} - y * j) * (x - Interval{m});
}

/// The single critical point of a piece's interval, bracketed as Refined does.
[[nodiscard]] Interval RefinedEnd(const Hodographs& h, const Vec3I& p, Interval x)
{
    for (int step = 0; step < RefineSteps; ++step) {
        const std::optional<Interval> k = KrawczykEnd(h, p, x);
        if (!k) {
            break;
        }
        const Interval next = Intersection(*k, x);
        if (next.IsEmpty() || Width(next) >= Width(x)) {
            break;
        }
        x = next;
    }
    return x;
}

/// Records the critical point of the piece's inflated interval when the
/// end lies within the resolution of the other curve there and is a
/// minimum along its own curve.
[[nodiscard]] Result<void> RecordEnd(EndSearch& e, const Arc& b, const Hodographs& h,
                                     const Interval& box)
{
    const Interval x = RefinedEnd(h, e.End.Point, box);
    const Segment& segment = (*e.Other)[b.Index];
    const Interval bracket = OnCurve(segment, InSegment(b, x));
    if (!Intersects(bracket, Intersection(Interval{segment.T0, segment.T1}, e.OtherDomain))) {
        return {}; // Past the domain's end, or the neighbouring segment's.
    }
    const Interval on = Intersection(bracket, e.OtherDomain);
    const auto jet =
        JetOver(HodographsOf(Arc{b.Index, 0.0, 1.0, segment.Poles}), LocalOf(segment, on));
    if (!jet) {
        return std::unexpected(EndUnresolved(e, "a weight not proved positive"));
    }
    const Vec3I d = Minus(e.End.Point, jet->Point);
    const Interval squared = Squared(d[0]) + Squared(d[1]) + Squared(d[2]);
    const double resolution = e.Pairs->Resolution;
    const Interval limit = Squared(Interval{resolution});
    if (squared.GetLower() > limit.GetUpper()) {
        return {};
    }
    if (!(squared.GetUpper() <= limit.GetLower())) {
        return std::unexpected(
            EndUnresolved(e, "a closest approach within rounding of the resolution"));
    }
    if (ShrinksInward(e.End, jet->Point)) {
        return {};
    }
    const Interval end{e.End.T};
    e.Found.push_back(e.Swapped ? Meeting{on, end} : Meeting{end, on});
    return {};
}

[[nodiscard]] bool ApartFromPoint(const Arc& b, const Vec3I& p, const double gap)
{
    const std::optional<Vec3I> box = BoxOf(b);
    if (!box) {
        return false;
    }
    for (std::size_t k = 0; k < 3; ++k) {
        if (FarApart((*box)[k], p[k], gap)) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] Result<void> VisitEnd(EndSearch& e, const Arc& b, const int depth)
{
    if (++e.Visited > PairBudget) {
        return std::unexpected(EndUnresolved(e, "the subdivision budget is spent"));
    }
    if (ApartFromPoint(b, e.End.Point, e.Pairs->Resolution)) {
        return {};
    }
    const Hodographs h = HodographsOf(b);
    const Interval x{-Inflation, 1.0 + Inflation};
    const std::optional<Interval> k = KrawczykEnd(h, e.End.Point, x);
    if (k && !Intersects(*k, x)) {
        return {};
    }
    if (k && InsideInterior(*k, x)) {
        return RecordEnd(e, b, h, x);
    }
    if (depth == MaxDepth) {
        return std::unexpected(EndUnresolved(e, "a piece neither apart nor a single minimum"));
    }
    const std::pair<Arc, Arc> halves = Halves(b);
    STEREON_TRY_VOID(VisitEnd(e, halves.first, depth + 1));
    return VisitEnd(e, halves.second, depth + 1);
}

/// Two ends within the resolution of each other, each a minimum along its
/// curve: a meeting at a corner of the parameter box.
[[nodiscard]] Result<void> Corner(PairSearch& search, const DomainEnd& a, const DomainEnd& b)
{
    const Vec3I d = Minus(a.Point, b.Point);
    const Interval squared = Squared(d[0]) + Squared(d[1]) + Squared(d[2]);
    const Interval limit = Squared(Interval{search.Resolution});
    if (squared.GetLower() > limit.GetUpper()) {
        return {};
    }
    if (!(squared.GetUpper() <= limit.GetLower())) {
        return Fail(ErrorCode::IntersectionFailed,
                    search.Handler + ": the ends at " + Number(a.T) + " and " + Number(b.T) +
                        " lie within rounding of the resolution apart");
    }
    if (!ShrinksInward(a, b.Point) && !ShrinksInward(b, a.Point)) {
        search.Found.push_back({Interval{a.T}, Interval{b.T}});
    }
    return {};
}

/// Each of one curve's ends against the other curve's stretch.
[[nodiscard]] Result<void> EndsAgainst(PairSearch& search, const std::vector<DomainEnd>& ends,
                                       const bool swapped)
{
    const std::vector<Segment>& other = swapped ? search.A : search.B;
    const Interval& domain = swapped ? search.DomainA : search.DomainB;
    for (const DomainEnd& end : ends) {
        EndSearch e{&search, end, &other, domain, swapped, 0, {}};
        for (std::size_t j = 0; j < other.size(); ++j) {
            STEREON_TRY_VOID(VisitEnd(e, Arc{j, 0.0, 1.0, other[j].Poles}, 0));
        }
        search.Found.insert(search.Found.end(), e.Found.begin(), e.Found.end());
    }
    return {};
}

/// The minima of the distance on the boundary of the parameter box: each
/// end of one curve against the other's stretch, and the corners. A
/// meeting inside the box is the pair search's; one on the boundary is a
/// critical point of the end's distance alone, which the pair search,
/// keeping only critical points in both domains, does not see: two edges
/// sharing a vertex a rounding apart, or nearly parallel segments touching
/// end to end, whose closest approach lies far past both.
[[nodiscard]] Result<void> SearchEnds(PairSearch& search)
{
    // A closed curve's whole domain has no ends: it runs on through its seam.
    const std::vector<DomainEnd> endsA =
        search.ClosedA ? std::vector<DomainEnd>{} : EndsOf(search.A, search.DomainA);
    const std::vector<DomainEnd> endsB =
        search.ClosedB ? std::vector<DomainEnd>{} : EndsOf(search.B, search.DomainB);
    STEREON_TRY_VOID(EndsAgainst(search, endsA, false));
    STEREON_TRY_VOID(EndsAgainst(search, endsB, true));
    for (const DomainEnd& a : endsA) {
        for (const DomainEnd& b : endsB) {
            STEREON_TRY_VOID(Corner(search, a, b));
        }
    }
    return {};
}

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

/// Both curves' pieces; an unbounded domain is first clipped to where its
/// curve can come near the other, bounded one.
struct BuiltPair {
    PieceCurve A;
    PieceCurve B;
};

[[nodiscard]] Result<PieceCurve> ClippedPieces(const Curve& c, const Interval& domain,
                                               const PieceCurve& other, const double resolution)
{
    STEREON_TRY(const Interval reach, Reach(c, BallOf(other.Segments, resolution)));
    const Interval clipped = Intersection(domain, reach);
    if (clipped.IsEmpty()) {
        return PieceCurve{{}, clipped, ParameterMap::Same, false, ""};
    }
    return PiecesOf(c, clipped);
}

[[nodiscard]] Result<BuiltPair> Build(const Curve& a, const Interval& da, const Curve& b,
                                      const Interval& db, const double resolution)
{
    if (!da.IsBounded() && !db.IsBounded()) {
        return Fail(ErrorCode::NotSupported,
                    "Intersect: two curves of unbounded domains; request a bounded one");
    }
    if (!da.IsBounded()) {
        STEREON_TRY(PieceCurve pb, PiecesOf(b, db));
        STEREON_TRY(PieceCurve pa, ClippedPieces(a, da, pb, resolution));
        return BuiltPair{std::move(pa), std::move(pb)};
    }
    STEREON_TRY(PieceCurve pa, PiecesOf(a, da));
    if (!db.IsBounded()) {
        STEREON_TRY(PieceCurve pb, ClippedPieces(b, db, pa, resolution));
        return BuiltPair{std::move(pa), std::move(pb)};
    }
    STEREON_TRY(PieceCurve pb, PiecesOf(b, db));
    return BuiltPair{std::move(pa), std::move(pb)};
}

/// Adds a meeting, or folds it into one already found within the
/// resolution of it on both curves (ADR-0001: points closer than that
/// cannot be told apart), keeping the nearer: a meeting at a domain's end
/// found both by the pair search and on the boundary, or at a closed
/// curve's seam from the pieces on both sides.
void Merge(std::vector<CurveCurvePoint>& points, const CurveCurvePoint& p, const Curve& a,
           const Curve& b, const double resolution)
{
    const auto gap = [&](const CurveCurvePoint& q) {
        return Distance(Eval(a, q.ParameterA), Eval(b, q.ParameterB));
    };
    for (CurveCurvePoint& q : points) {
        if (Distance(Eval(a, q.ParameterA), p.Position) <= resolution &&
            Distance(Eval(b, q.ParameterB), Eval(b, p.ParameterB)) <= resolution) {
            if (gap(p) < gap(q)) {
                q = p;
            }
            return;
        }
    }
    points.push_back(p);
}

template <typename T>
[[nodiscard]] bool SameCurve(const Curve& a, const Curve& b)
{
    const T* x = a.Target<T>();
    const T* y = b.Target<T>();
    return x != nullptr && y != nullptr && *x == *y;
}

/// Two identical curves coincide on the requested domains' common stretch.
/// A circle or an ellipse is compared in its base turn: outside it, the
/// same point has another parameter, so the search answers (it fails).
[[nodiscard]] bool Identical(const Curve& a, const Interval& da, const Curve& b, const Interval& db)
{
    const Interval turn{0.0, 2.0 * std::numbers::pi};
    const bool conic = SameCurve<Circle>(a, b) || SameCurve<Ellipse>(a, b);
    return SameCurve<BSplineCurve>(a, b) || SameCurve<Parabola>(a, b) ||
           SameCurve<Hyperbola>(a, b) || (conic && Contains(turn, da) && Contains(turn, db));
}

[[nodiscard]] CurveCurveIntersection Coincident(const Interval& da, const Interval& db)
{
    CurveCurveIntersection result;
    const Interval common = Intersection(da, db);
    if (!common.IsEmpty()) {
        result.Overlaps.push_back({common, common, true});
    }
    return result;
}

[[nodiscard]] Result<CurveCurveIntersection> Searched(const Curve& a, const PieceCurve& pa,
                                                      const Curve& b, const PieceCurve& pb,
                                                      const Context& ctx)
{
    PairSearch search;
    search.Handler = pa.Name + " x " + pb.Name;
    search.DomainA = pa.Domain;
    search.DomainB = pb.Domain;
    search.ClosedA = pa.Closed;
    search.ClosedB = pb.Closed;
    search.Resolution = ctx.GetResolution();
    search.A = pa.Segments;
    search.B = pb.Segments;
    STEREON_TRY_VOID(SearchAll(search));
    STEREON_TRY_VOID(SearchEnds(search));
    CurveCurveIntersection result;
    const auto placed = [](const PieceCurve& c, const Interval& on) {
        const double u = std::clamp(Midpoint(on), c.Domain.GetLower(), c.Domain.GetUpper());
        return std::clamp(CurveParameter(c, u), c.Domain.GetLower(), c.Domain.GetUpper());
    };
    for (const Meeting& m : Distinct(std::move(search.Found))) {
        const double ta = placed(pa, m.OnA);
        const double tb = placed(pb, m.OnB);
        Merge(result.Points, {ta, tb, Eval(a, ta)}, a, b, search.Resolution);
    }
    std::ranges::sort(result.Points, [](const CurveCurvePoint& p, const CurveCurvePoint& q) {
        return p.ParameterA < q.ParameterA;
    });
    return result;
}

} // namespace

Result<CurveCurveIntersection> CurvePair(const Curve& a, const Interval& domainA, const Curve& b,
                                         const Interval& domainB, const Context& ctx)
{
    STEREON_TRY(const Interval da, RequestedOf(a, domainA));
    STEREON_TRY(const Interval db, RequestedOf(b, domainB));
    if (da.IsEmpty() || db.IsEmpty()) {
        return CurveCurveIntersection{};
    }
    if (Identical(a, da, b, db)) {
        return Coincident(da, db);
    }
    STEREON_TRY(const BuiltPair built, Build(a, da, b, db, ctx.GetResolution()));
    if (built.A.Domain.IsEmpty() || built.B.Domain.IsEmpty()) {
        return CurveCurveIntersection{};
    }
    return Searched(a, built.A, b, built.B, ctx);
}

} // namespace Stereon::Detail
