// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// An isolated touch, for the generic surface-surface intersection. Included
// once, in GenericSSI.cpp, after TangentBranch.hpp, whose helpers and those
// of "Branch points" it uses; it is a part of that unit kept in its own file.
//
// The surfaces share a tangent plane at a point and their relative second
// form there is definite: near it they part on every side. Such a point is
// settled on a pair (S, C), S's box R about it, in a frame (e1, e2, n) about
// C's normal there, by facts proved in interval arithmetic:
//
// 1. Coordinates and containment, as for a node (GenericSSI.cpp, "Branch
//    points", facts 1 and 2): the height h of S over C's graph is defined
//    on R, its zeros are the points the surfaces share there, and |h|
//    bounds how far a point of S is from C.
// 2. The critical point s* of h, isolated by Krawczyk's operator strictly
//    inside both boxes, the surfaces within the tolerance there (fact 3 of
//    a node).
// 3. Convexity. The Hessian of h, enclosed over R, is definite: h_uu
//    excludes zero and det H > 0 for every H in the enclosure. h is then
//    strictly convex (or, h_uu < 0, concave) on the convex box R, s* is its
//    only critical point there, and with lambda > 0 a lower bound on the
//    smaller eigenvalue's magnitude, sigma = sign h_uu,
//        sigma (h(s) - h(s*)) >= lambda |s - s*|^2 / 2.
//
// What follows, with h* = h(s*). If sigma h* > 0 (proved), h has no zero on
// R: the surfaces do not meet there, and are within the tolerance of each
// other at s*. If sigma h* may be zero or below, the zeros of h on R lie
// in the disc |s - s*| <= rho, rho = sqrt(2 max(0, -sigma h*) / lambda):
// every common point of the region is that close to s* in S's parameters,
// and so within |S_s| rho of S(s*) in space. The touch is reported as the
// point S(s*), its bound the larger of the gap at s* and that reach, only
// where that reach is within the resolution: a loop of the zero set any
// wider is the search's to find as a curve, never a point (a candidate for
// an event, ADR-0054; a gap within the tolerance is declared, ADR-0050).
//
// The search. A pair with no advance direction whose cells hold a candidate
// (Krawczyk's isolation of the critical point within the tolerance, the
// relative form definite in doubles) is centred on it: S's box reaches
// about s* as far as the surfaces need to part by TouchTolerances
// tolerances by the relative form's smaller eigenvalue, inside S's Bezier
// patch, and C's covers its projection inside C's; both are widened to a
// dyadic grid, and halved until the proof holds. The product of the two
// boxes is recorded as a settled region where nothing was found yet, and
// carved from the search as a node's is.
//
// Not settled, and failing closed: a touch on a knot line or a seam of
// either surface, or on the boundary of a requested domain (no single
// Bezier patch holds S's box about it); a touch at a singular point of a
// parameterization (a sphere's pole); a degenerate relative form; a loop
// of the zero set wider than the resolution that the search cannot resolve
// as a curve.

namespace Stereon::Detail {
namespace {

/// How many tolerances the surfaces part by at the sides of a touch's box.
constexpr double TouchTolerances = 8.0;
/// How often a touch's boxes are halved before the candidate is left to the
/// cutting.
constexpr int TouchNarrowings = 4;

/// The Hessian's enclosure definite: h_uu of one sign and det H > 0.
[[nodiscard]] bool Definite(const Matrix2& h)
{
    const Interval det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
    return ExcludesZero(h[0][0]) && det.GetLower() > 0.0;
}

/// A lower bound on the smaller eigenvalue's magnitude of every symmetric
/// matrix in a definite enclosure: det / |trace|, the product of the
/// eigenvalues over their sum.
[[nodiscard]] double SmallerEigenvalue(const Matrix2& h)
{
    const Interval det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
    const Interval trace = h[0][0] + h[1][1];
    const double largest = std::max(std::abs(trace.GetLower()), std::abs(trace.GetUpper()));
    return largest > 0.0 ? (Interval{det.GetLower()} / Interval{largest}).GetLower() : 0.0;
}

/// A proved touch in the pair's own order (S, C): the point, and the bound
/// on its distance from every common point of the region and from C.
struct TouchAt {
    Seed Point;
    double Bound;
};

/// How far from S(s*) in space a common point of R can lie: |S_s| rho, rho
/// the disc's radius in parameters, |S_s| bounded by the Frobenius norm of
/// S's partials over R. Zero where sigma h* is proved positive.
[[nodiscard]] double ReachOfLoop(const Interval& height, const Matrix2& h, const Jet& overS)
{
    const double sigma = h[0][0].GetLower() > 0.0 ? 1.0 : -1.0;
    const double below = -(sigma > 0.0 ? height.GetLower() : -height.GetUpper());
    if (!(below > 0.0)) {
        return 0.0;
    }
    const double lambda = SmallerEigenvalue(h);
    if (!(lambda > 0.0)) {
        return Infinity;
    }
    const Interval rho = Sqrt(Interval{2.0} * Interval{below} / Interval{lambda});
    const Interval speed = Sqrt(DotOf(overS.Pu, overS.Pu) + DotOf(overS.Pv, overS.Pv));
    return (rho * speed).GetUpper();
}

/// Facts 1 to 3 on a pair (small, cover) in a frame about `normal`, and the
/// touch's bound: nothing unless the reach of the zero set about s* is
/// within the resolution.
[[nodiscard]] std::optional<TouchAt> ProvedTouch(const PatchPair& pair, const Vec3& normal,
                                                 const SearchSettings& s)
{
    const std::optional<NodeFrame> frame = FrameAbout(normal);
    const std::optional<PairJets> j = frame ? JetsOfPair(pair) : std::nullopt;
    if (!j) {
        return std::nullopt;
    }
    const NodeSystem sys{pair, j->OnS, j->OnC, *frame, j->OverS, j->OverC};
    const std::optional<IntervalsN<4>> x = ProvedCritical(sys, s.Tolerance);
    const std::optional<Matrix2> inverse =
        x && Covered(sys) ? InverseProjection(sys) : std::nullopt;
    const std::optional<Matrix2> hessian =
        inverse ? std::optional<Matrix2>{HeightHessian(sys, *inverse)} : std::nullopt;
    if (!hessian || !Definite(*hessian)) {
        return std::nullopt;
    }
    const Bounds onS = XyzOver(pair.A, (*x)[0], (*x)[1]);
    const Interval height = DotOf(frame->N, MinusOf(onS, XyzOver(pair.B, (*x)[2], (*x)[3])));
    const double reach = ReachOfLoop(height, *hessian, j->OverS);
    // The resolution: the join is 2^-7 of it, exactly.
    if (!(reach <= 0x1p7 * s.Join)) {
        return std::nullopt;
    }
    const std::array<double, 4> m = MidpointsOf<4>(*x);
    const double gap = GapAt(pair, *x).GetUpper();
    return TouchAt{{Middle(onS), Point2{m[0], m[1]}, Point2{m[2], m[3]}}, std::max(gap, reach)};
}

/// Where a pair (small, cover) may hold a touch: a critical point isolated
/// within the tolerance, its relative form definite in doubles. The point's
/// parameters and the form, or nothing. A guess.
struct TouchGuess {
    Point2 OnS;
    Point2 OnC;
    Form2 Relative;
};

/// The candidate in the pair (small, cover), or nothing.
[[nodiscard]] std::optional<TouchGuess> TouchCandidate(const PatchPair& pair,
                                                       const SearchSettings& s)
{
    const std::optional<PairJets> j = JetsOfPair(pair);
    const std::optional<NodeFrame> frame =
        j ? FrameAbout(MiddleVec(CrossOf(j->OverC.Pu, j->OverC.Pv))) : std::nullopt;
    if (!frame) {
        return std::nullopt;
    }
    const Isolated<4> critical = CriticalPoint({pair, j->OnS, j->OnC, *frame, j->OverS, j->OverC});
    if (critical.Outcome != Isolation::Unique ||
        GapAt(pair, critical.Box).GetUpper() > s.Tolerance) {
        return std::nullopt;
    }
    const std::array<double, 4> m = MidpointsOf<4>(critical.Box);
    const Vec3 normal = NormalAt(pair.B, j->OnC, m[2], m[3]);
    if (!(Length(normal) > 0.0)) {
        return std::nullopt;
    }
    const Vec3 n = Normalized(normal);
    const std::optional<LocalForm> onS = FormAt(pair.A, j->OnS, m[0], m[1], n);
    const std::optional<LocalForm> onC = FormAt(pair.B, j->OnC, m[2], m[3], n);
    if (!onS || !onC || FormSign(*onS, *onC) <= 0) {
        return std::nullopt;
    }
    const std::optional<Form2> relative = RelativeForm(*onS, *onC);
    if (!relative) {
        return std::nullopt;
    }
    return TouchGuess{Point2{m[0], m[1]}, Point2{m[2], m[3]}, *relative};
}

/// The half-width of S's box about a touch, in its parameters: where the
/// surfaces part by TouchTolerances tolerances, by the relative form's
/// smaller eigenvalue.
[[nodiscard]] std::optional<double> TouchRadius(const Form2& m, const double tolerance)
{
    const double mean = 0.5 * (m[0] + m[2]);
    const double spread = std::hypot(0.5 * (m[0] - m[2]), m[1]);
    const double smaller = std::min(std::abs(mean - spread), std::abs(mean + spread));
    if (!(smaller > 0.0) || !std::isfinite(smaller)) {
        return std::nullopt;
    }
    return std::sqrt(2.0 * TouchTolerances * tolerance / smaller);
}

/// The boxes of a touch: S's about s*, `r` either way, and C's covering its
/// projection, by the projections' Jacobians at the point, each inside its
/// Bezier patch and widened to a dyadic grid; nothing when s* is not
/// strictly inside S's.
[[nodiscard]] std::optional<ValleyBoxes> TouchBoxes(const TracePatches& p, const TouchGuess& g,
                                                    const double r)
{
    const Box2& own = p.S->Parameters;
    const Box2 onS{Snapped(Intersection(Interval{g.OnS.X - r, g.OnS.X + r}, own.X), own.X),
                   Snapped(Intersection(Interval{g.OnS.Y - r, g.OnS.Y + r}, own.Y), own.Y)};
    const std::array<double, 2> margins = Margins(g.OnS, onS);
    const std::optional<Jet> s =
        JetAt(*p.S, SecondNetsOf(*p.S), Interval{g.OnS.X}, Interval{g.OnS.Y});
    const std::optional<Jet> c =
        JetAt(*p.C, SecondNetsOf(*p.C), Interval{g.OnC.X}, Interval{g.OnC.Y});
    if (!s || !c || !(std::min(margins[0], margins[1]) > 0.0)) {
        return std::nullopt;
    }
    // C's parameters per unit of S's: the least-squares map C_s = J^+ S_s,
    // in C's metric, at the point.
    const Vec3 cu = MiddleVec(c->Pu);
    const Vec3 cv = MiddleVec(c->Pv);
    const double e = Dot(cu, cu);
    const double f = Dot(cu, cv);
    const double gg = Dot(cv, cv);
    const double det = e * gg - f * f;
    if (!(det > 0.0)) {
        return std::nullopt;
    }
    const auto toC = [&](const Vec3& t, const double by) {
        const double x = Dot(t, cu);
        const double y = Dot(t, cv);
        return std::array{std::abs((gg * x - f * y) / det) * by,
                          std::abs((e * y - f * x) / det) * by};
    };
    const std::array<double, 2> du = toC(MiddleVec(s->Pu), Width(onS.X));
    const std::array<double, 2> dv = toC(MiddleVec(s->Pv), Width(onS.Y));
    const Box2& cover = p.C->Parameters;
    const Box2 onC{
        Snapped(Intersection(Inflated(Interval{g.OnC.X}, du[0] + dv[0]), cover.X), cover.X),
        Snapped(Intersection(Inflated(Interval{g.OnC.Y}, du[1] + dv[1]), cover.Y), cover.Y)};
    return ValleyBoxes{onS, onC};
}

/// A touch proved on the boxes about a candidate, halved until it holds;
/// the boxes and the touch, or nothing.
struct SettledTouch {
    ValleyBoxes Boxes;
    TouchAt Touch;
};

/// The touch about a candidate in the patches holding it, or nothing.
[[nodiscard]] std::optional<SettledTouch> TouchAround(const TracePatches& p, const TouchGuess& g,
                                                      const SearchSettings& s)
{
    const std::optional<double> radius = TouchRadius(g.Relative, s.Tolerance);
    if (!radius) {
        return std::nullopt;
    }
    const Vec3 normal = NormalAt(*p.C, SecondNetsOf(*p.C), g.OnC.X, g.OnC.Y);
    for (int k = 0; k < TouchNarrowings; ++k) {
        const std::optional<ValleyBoxes> boxes = TouchBoxes(p, g, std::ldexp(*radius, -k));
        std::optional<Patch> small = boxes ? ClippedTo(*p.S, boxes->OnS) : std::nullopt;
        std::optional<Patch> cover = boxes ? ClippedTo(*p.C, boxes->OnC) : std::nullopt;
        if (!small || !cover) {
            continue;
        }
        if (const std::optional<TouchAt> touch =
                ProvedTouch({std::move(*small), std::move(*cover)}, normal, s)) {
            return SettledTouch{*boxes, *touch};
        }
    }
    return std::nullopt;
}

/// The touch and its region into the findings, in the order (A, B).
void RecordTouch(const SettledTouch& t, const bool smallIsA, Findings& found)
{
    const Seed point = smallIsA ? t.Touch.Point : Swapped(t.Touch.Point);
    found.Touches.push_back({point.Position, point.OnA, point.OnB, t.Touch.Bound});
    const Box2& onA = smallIsA ? t.Boxes.OnS : t.Boxes.OnC;
    const Box2& onB = smallIsA ? t.Boxes.OnC : t.Boxes.OnS;
    found.Regions.push_back({{onA}, {onB}});
}

/// A touch of the pair, S the surface of A when smallIsA, or nothing.
[[nodiscard]] std::optional<SettledTouch> TouchOf(const PatchPair& pair, const SearchSettings& s,
                                                  const bool smallIsA)
{
    const PatchPair oriented = smallIsA ? pair : PatchPair{pair.B, pair.A};
    const std::optional<TouchGuess> guess = TouchCandidate(oriented, s);
    const std::optional<TracePatches> patches =
        guess ? PatchesHolding(pair, s, smallIsA) : std::nullopt;
    return patches ? TouchAround(*patches, *guess, s) : std::nullopt;
}

/// A pair with no advance direction that may hold an isolated touch, either
/// surface taken as S: the touch proved on boxes about it in the Bezier
/// patches holding the pair's cells, and recorded where nothing was found
/// yet. The pair goes back on the stack, to be carved by the region.
[[nodiscard]] bool SettleTouchPoint(const PatchPair& pair, const SearchSettings& s, Findings& found,
                                    std::vector<PatchPair>& stack)
{
    if (s.OnA == nullptr || s.OnB == nullptr) {
        return false;
    }
    for (const bool smallIsA : {true, false}) {
        const std::optional<SettledTouch> touch = TouchOf(pair, s, smallIsA);
        if (touch && NothingFoundIn(smallIsA ? touch->Boxes.OnS : touch->Boxes.OnC,
                                    smallIsA ? touch->Boxes.OnC : touch->Boxes.OnS, found)) {
            RecordTouch(*touch, smallIsA, found);
            stack.push_back(pair);
            return true;
        }
    }
    return false;
}

} // namespace
} // namespace Stereon::Detail
