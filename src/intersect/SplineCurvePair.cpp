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

/// Two B-spline curves in space, polynomial or rational. Two curves in space
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
/// the intersection fails. A critical point of the polynomials continued a
/// little past a requested domain's end is a meeting at that end when the
/// distance there is within the resolution: two edges sharing a vertex meet
/// at it even when their ends are a rounding apart. Critical points proved
/// twice, from overlapping boxes, are one. A pair that is none of these
/// when the pieces can no longer be halved, or a search past the budget,
/// fails with IntersectionFailed: tangent curves, an overlap, or a meeting
/// too nearly tangential to separate is never dropped and never reported
/// unproved (ADR-0029). Two identical splines coincide: their overlap is
/// the requested domains' common stretch on the identity map.

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
/// (a bracket across a knot is proved from both segments and merged); the
/// domain's end when it lies just past it (the polynomial continued: two
/// ends a rounding apart); nothing when it lies past an end inside the
/// domain, where the neighbouring segment proves it.
[[nodiscard]] std::optional<Interval> Kept(const Interval& bracket, const Interval& stretch,
                                           const Interval& domain)
{
    if (Intersects(bracket, stretch)) {
        return Intersection(bracket, domain);
    }
    if (bracket.GetUpper() < stretch.GetLower() && stretch.GetLower() == domain.GetLower()) {
        return Interval{domain.GetLower()};
    }
    if (bracket.GetLower() > stretch.GetUpper() && stretch.GetUpper() == domain.GetUpper()) {
        return Interval{domain.GetUpper()};
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
    // a bracket kept at a domain's end lies outside the piece.
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

[[nodiscard]] Result<void> CheckSpline(const BSplineCurve& curve, const Interval& domain,
                                       const std::string_view where)
{
    STEREON_TRY_VOID(CheckKnotVector(curve.Degree, curve.Knots, curve.Poles.size(), where));
    STEREON_TRY_VOID(CheckWeights(curve.Weights, curve.Poles.size(), where));
    for (const Point3& pole : curve.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return CheckPeriodicDomain(curve, domain, where);
}

} // namespace

Result<CurveCurveIntersection> SplineSpline(const Curve& a, const Interval& domainA, const Curve& b,
                                            const Interval& domainB, const Context& ctx)
{
    const BSplineCurve* sa = a.Target<BSplineCurve>();
    const BSplineCurve* sb = b.Target<BSplineCurve>();
    if (sa == nullptr || sb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineSpline: dispatched for other curve types");
    }
    PairSearch search;
    search.Handler = "BSplineCurve x BSplineCurve";
    STEREON_TRY_VOID(CheckSpline(*sa, domainA, search.Handler));
    STEREON_TRY_VOID(CheckSpline(*sb, domainB, search.Handler));
    search.DomainA = Intersection(domainA, Domain(*sa));
    search.DomainB = Intersection(domainB, Domain(*sb));
    CurveCurveIntersection result;
    if (search.DomainA.IsEmpty() || search.DomainB.IsEmpty()) {
        return result;
    }
    if (*sa == *sb) {
        result.Overlaps.push_back({Intersection(search.DomainA, search.DomainB),
                                   Intersection(search.DomainA, search.DomainB), true});
        if (result.Overlaps.front().OnA.IsEmpty()) {
            result.Overlaps.clear();
        }
        return result;
    }
    search.Resolution = ctx.GetResolution();
    STEREON_TRY(std::vector<Segment> segmentsA, SegmentsOf(*sa));
    STEREON_TRY(std::vector<Segment> segmentsB, SegmentsOf(*sb));
    search.A = Reaching(std::move(segmentsA), search.DomainA);
    search.B = Reaching(std::move(segmentsB), search.DomainB);
    STEREON_TRY_VOID(SearchAll(search));
    for (const Meeting& m : Distinct(std::move(search.Found))) {
        const double ta =
            std::clamp(Midpoint(m.OnA), search.DomainA.GetLower(), search.DomainA.GetUpper());
        const double tb =
            std::clamp(Midpoint(m.OnB), search.DomainB.GetLower(), search.DomainB.GetUpper());
        result.Points.push_back({ta, tb, Eval(*sa, ta)});
    }
    std::ranges::sort(result.Points, [](const CurveCurvePoint& p, const CurveCurvePoint& q) {
        return p.ParameterA < q.ParameterA;
    });
    return result;
}

} // namespace Stereon::Detail
