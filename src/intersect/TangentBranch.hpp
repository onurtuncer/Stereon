// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Tangency along a curve, for the generic surface-surface intersection.
// Included once, in GenericSSI.cpp, after SeamNode.hpp, whose helpers and
// those of "Branch points" it uses; it is a part of that unit kept in its
// own file.
//
// Two surfaces tangent along a curve have no transverse pair of patches
// around it at any size, and no isolated point of common tangency. The
// branch is settled as a valley of one surface's height over the other.
//
// The setting. S is a Bezier patch of one surface and C one of the other,
// a the parameter of S across the curve and b the one along it. Each step
// along b is proved on two pairs of boxes in one frame (e1, e2, n): a wide
// one, R = [a0, a1] x [b0, b1] of S across the band the step settles with a
// box of C covering it, and a narrow one, R' = [a0', a1'] x [b0, b1] inside
// R with a box of C inside the wide one's, close about the valley, where a
// Hessian's enclosure is tight enough for facts 4 and 5.
//
// 1. Coordinates and containment, on both pairs, as for a node
//    (GenericSSI.cpp, "Branch points", facts 1 and 2): C's box projects one
//    to one onto the frame's plane, and the projection of S over the box
//    keeps clear of that of C's edges and meets it at the root below. The
//    height h(s) = n . S(s) - g(e . S(s)) of S over C's graph g is then
//    defined on R, the same on R' under the narrow box of C, its zeros are
//    exactly the points of S over R on C's box, and |h| bounds from above
//    how far the point of S is from C.
// 2. Convexity across. The Hessian of h is enclosed over R, and its entry
//    h_aa excludes zero: along every line b = const, h is strictly convex
//    or strictly concave on R, so h_a vanishes at most once on it.
// 3. A root. Krawczyk's operator isolates a point (a#, b0) of R''s edge
//    b = b0 with h_a = 0, the solution of
//        e . (S - C) = 0,  S_a . N_C = 0
//    in (a, u_C, v_C): h_a = nu . S_a for nu a normal of C at the foot.
// 4. The valley spans R'. By the mean value theorem on the convex box,
//        h_a(a, b) = h_aa (a - a#) + h_ab (b - b0)
//    for entries in the Hessian's enclosure over R'. This is proved of one
//    sign on the side a = a0' and of the other on a = a1', for every b.
//    With fact 2, each line b = const holds exactly one zero of h_a in R,
//    and it lies inside R': the valley, a curve a = phi(b) from the edge
//    b = b0 to the edge b = b1, phi' = -h_ab / h_aa.
// 5. Within the tolerance. Along the valley, by Taylor's theorem about the
//    root with grad h = (0, h_b) there and phi enclosed by fact 4, |h| is
//    proved at most the pair tolerance.
//
// What follows. On each line b = const, h has at most two zeros, one on
// either side of the valley's point, and between a zero and that point h
// lies between 0 and its value there. So every point the two surfaces
// share over R is joined to the valley, across it, by points of S within
// the tolerance of C: the valley is the tangential branch of ADR-0029
// section 6, and it accounts for every intersection in the region R x
// (C's box). Whether the surfaces touch, cross twice within the tolerance
// or stay apart within it is not told apart (ADR-0050); the branch's
// transition is left unresolved.
//
// The trace. The valley is followed across one whole Bezier patch of S,
// from its edge b = b0 to its edge b = b1, in boxes each proved as above;
// all of them, or none is kept. Each box reaches across at least as far as
// the band where the surfaces part by BandTolerances tolerances, so that the
// pairs left outside the region it settles can be proved apart. Each box
// becomes an arc, marched on the valley's own system (GenericSSI.cpp,
// "Marching a tangency along a curve"), and a settled region the search
// carves away. The ends lie on S's knot lines or on the boundary of the
// requested domain; the next patch's trace finds the same end from its
// side, since along the shared edge S_a is the edge's own tangent. A trace
// that fails over a pair of Bezier patches is not tried again from another
// pair of their cells.
//
// What is not settled fails closed as before: a curve of tangency that
// leaves S's patch across a (it turns along a parameter line, or closes
// up), that reaches an edge of C's patch (a knot line, a seam or the
// boundary of its domain), that runs along an edge of S, that ends inside
// both surfaces, or whose surfaces part by more than the tolerance; and
// surfaces tangent along curves of both kinds, across u and across v, or
// with either surface as S: the first trace fixes both for the search.

namespace Stereon::Detail {
namespace {

/// Newton's iterations for a point of the valley, in doubles.
constexpr int ValleyIterations = 16;
/// The boxes of one trace, and how often a step along b is halved.
constexpr std::size_t ValleyBoxLimit = std::size_t{1} << 11;
constexpr int ValleyHalvings = 14;
/// The steps of the walk back to the patch's edge, in doubles.
constexpr int ValleyWalk = 32;

/// The small surface's parameters from those across and along the valley.
[[nodiscard]] Point2 OnSmall(const bool acrossU, const double across, const double along)
{
    return acrossU ? Point2{across, along} : Point2{along, across};
}

/// A point of the valley in doubles: S's parameters and its foot on C.
struct ValleyPoint {
    double Across;
    double Along;
    Point2 OnC;
};

/// One trace's patches: a Bezier patch of each surface with its nets, which
/// surface is S, and which of its parameters runs across the valley.
struct ValleyTrace {
    const Patch& S;
    const SecondNets& OnS;
    const Patch& C;
    const SecondNets& OnC;
    bool SmallIsA;
    bool AcrossU;
    double Tolerance;
    double Join;
};

/// The plane (e1, e2) of a frame about C's normal at a foot, in doubles.
[[nodiscard]] std::optional<std::array<Vec3, 2>> PlaneAt(const ValleyTrace& t, const Point2& onC)
{
    const std::optional<NodeFrame> f = FrameAbout(NormalAt(t.C, t.OnC, onC.X, onC.Y));
    if (!f) {
        return std::nullopt;
    }
    return std::array{MiddleVec(f->E1), MiddleVec(f->E2)};
}

/// The two jets at a point of the valley's system, in doubles' enclosures.
struct ValleyJets {
    Jet S;
    Jet C;
};

/// S's and C's jets at a point of the valley's system.
[[nodiscard]] std::optional<ValleyJets> JetsAt(const ValleyTrace& t, const ValleyPoint& p)
{
    const Point2 uv = OnSmall(t.AcrossU, p.Across, p.Along);
    std::optional<Jet> s = JetAt(t.S, t.OnS, Interval{uv.X}, Interval{uv.Y});
    std::optional<Jet> c = JetAt(t.C, t.OnC, Interval{p.OnC.X}, Interval{p.OnC.Y});
    if (!s || !c) {
        return std::nullopt;
    }
    return ValleyJets{std::move(*s), std::move(*c)};
}

/// Newton's step of e . (S - C) = 0, S_a . N_C = 0 in (a, u_C, v_C), b held.
[[nodiscard]] std::optional<std::array<double, 3>>
ValleyStep(const ValleyTrace& t, const std::array<Vec3, 2>& e, const ValleyPoint& p)
{
    const std::optional<ValleyJets> at = JetsAt(t, p);
    if (!at) {
        return std::nullopt;
    }
    const auto m = [](const Bounds& b) { return MiddleVec(b); };
    const Jet& s = at->S;
    const Jet& c = at->C;
    const Vec3 sa = m(t.AcrossU ? s.Pu : s.Pv);
    const Vec3 saa = m(t.AcrossU ? s.Puu : s.Pvv);
    const Vec3 r = m(s.P) - m(c.P);
    const Vec3 n = Cross(m(c.Pu), m(c.Pv));
    const Vec3 nu = Cross(m(c.Puu), m(c.Pv)) + Cross(m(c.Pu), m(c.Puv));
    const Vec3 nv = Cross(m(c.Puv), m(c.Pv)) + Cross(m(c.Pu), m(c.Pvv));
    const DoubleMatrix<3> j{{{Dot(e[0], sa), -Dot(e[0], m(c.Pu)), -Dot(e[0], m(c.Pv))},
                             {Dot(e[1], sa), -Dot(e[1], m(c.Pu)), -Dot(e[1], m(c.Pv))},
                             {Dot(saa, n), Dot(sa, nu), Dot(sa, nv)}}};
    const std::array<double, 3> f{Dot(e[0], r), Dot(e[1], r), Dot(sa, n)};
    const std::optional<DoubleMatrix<3>> inverse = Inverted<3>(j);
    if (!inverse) {
        return std::nullopt;
    }
    std::array<double, 3> step{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t k = 0; k < 3; ++k) {
            step[row] -= (*inverse)[row][k] * f[k];
        }
    }
    return step;
}

/// The valley's point at b = along, from a guess, by Newton's method kept
/// in the patches; nothing unless the last step is a 2^-30th of them. A
/// guess: it proves nothing and is trusted with nothing.
[[nodiscard]] std::optional<ValleyPoint>
ValleyAt(const ValleyTrace& t, const std::array<Vec3, 2>& e, ValleyPoint p, const double along)
{
    const std::array<Interval, 3> box{StretchOf(t.S, t.AcrossU), t.C.Parameters.X,
                                      t.C.Parameters.Y};
    p.Along = along;
    for (int k = 0; k < ValleyIterations; ++k) {
        const std::optional<std::array<double, 3>> step = ValleyStep(t, e, p);
        if (!step) {
            return std::nullopt;
        }
        p.Across = std::clamp(p.Across + (*step)[0], box[0].GetLower(), box[0].GetUpper());
        p.OnC = {std::clamp(p.OnC.X + (*step)[1], box[1].GetLower(), box[1].GetUpper()),
                 std::clamp(p.OnC.Y + (*step)[2], box[2].GetLower(), box[2].GetUpper())};
        bool small = true;
        for (std::size_t r = 0; r < 3; ++r) {
            small = small && std::abs((*step)[r]) <= 0x1p-30 * Width(box[r]);
        }
        if (small) {
            return p;
        }
    }
    return std::nullopt;
}

/// The valley's point at b = along, in a frame about C's normal at the
/// guess's foot.
[[nodiscard]] std::optional<ValleyPoint> ValleyNear(const ValleyTrace& t, const ValleyPoint& guess,
                                                    const double along)
{
    const std::optional<std::array<Vec3, 2>> e = PlaneAt(t, guess.OnC);
    return e ? ValleyAt(t, *e, guess, along) : std::nullopt;
}

/// The root of fact 3 on the line b = along of the system's small patch:
/// Krawczyk's operator on e . (S - C) = 0, S_a . N_C = 0 in (a, u_C, v_C),
/// the Jacobian enclosed over the whole pair.
[[nodiscard]] Isolated<3> ValleyRoot(const NodeSystem& sys, const bool acrossU, const double along)
{
    const Patch& small = sys.Pair.A;
    const Patch& cover = sys.Pair.B;
    const NormalJet n = NormalOf(sys.OverC);
    const Jet& jc = sys.OverC;
    const NodeFrame& f = sys.Frame;
    const Bounds& sa = acrossU ? sys.OverS.Pu : sys.OverS.Pv;
    const Bounds& saa = acrossU ? sys.OverS.Puu : sys.OverS.Pvv;
    IntervalMatrix<3> j{};
    j[0] = {DotOf(f.E1, sa), -DotOf(f.E1, jc.Pu), -DotOf(f.E1, jc.Pv)};
    j[1] = {DotOf(f.E2, sa), -DotOf(f.E2, jc.Pu), -DotOf(f.E2, jc.Pv)};
    j[2] = {DotOf(saa, n.N), DotOf(sa, n.Nu), DotOf(sa, n.Nv)};
    const IntervalsN<3> box{StretchOf(small, acrossU), cover.Parameters.X, cover.Parameters.Y};
    return Krawczyk<3>(box, j, [&](const std::array<double, 3>& m) {
        const Point2 uv = OnSmall(acrossU, m[0], along);
        const std::optional<Jet> at = JetAt(small, sys.S, Interval{uv.X}, Interval{uv.Y});
        const std::optional<Jet> on = JetAt(cover, sys.C, Interval{m[1]}, Interval{m[2]});
        if (!at || !on) {
            const Interval all = Interval::Everything();
            return IntervalsN<3>{all, all, all};
        }
        const Bounds r = MinusOf(at->P, on->P);
        return IntervalsN<3>{DotOf(f.E1, r), DotOf(f.E2, r),
                             DotOf(acrossU ? at->Pu : at->Pv, CrossOf(on->Pu, on->Pv))};
    });
}

/// The root's box strictly inside the pair's: a across S's stretch, the
/// foot in C's box.
[[nodiscard]] bool RootInside(const NodeSystem& sys, const bool acrossU, const IntervalsN<3>& x)
{
    return Inside(x[0], StretchOf(sys.Pair.A, acrossU)) && Inside(x[1], sys.Pair.B.Parameters.X) &&
           Inside(x[2], sys.Pair.B.Parameters.Y);
}

/// The height and its slope along b at the root of fact 3: h = n . (S - C)
/// and h_b = nu . S_b, nu the co-normal of C at the foot.
struct ValleyHeight {
    Interval Height;
    Interval Slope;
};

/// h and h_b at the root, enclosed over its box.
[[nodiscard]] std::optional<ValleyHeight> HeightAtRoot(const NodeSystem& sys, const bool acrossU,
                                                       const IntervalsN<3>& x, const double along)
{
    const Interval b{along};
    const std::optional<Jet> at = JetAt(sys.Pair.A, sys.S, acrossU ? x[0] : b, acrossU ? b : x[0]);
    const std::optional<Jet> on = JetAt(sys.Pair.B, sys.C, x[1], x[2]);
    if (!at || !on) {
        return std::nullopt;
    }
    const std::optional<Matrix2> inverse = InverseOf(ProjectedJacobian(sys.Frame, on->Pu, on->Pv));
    if (!inverse) {
        return std::nullopt;
    }
    const Bounds nu = CoNormal(sys.Frame, *inverse, *on);
    return ValleyHeight{DotOf(sys.Frame.N, MinusOf(at->P, on->P)),
                        DotOf(nu, acrossU ? at->Pv : at->Pu)};
}

/// Facts 4 and 5 from the Hessian's enclosure over the box, the root's a
/// and the height there: the valley crosses neither side a = a0 nor
/// a = a1, and |h| along it is at most the tolerance.
[[nodiscard]] bool ValleySpans(const Matrix2& hessian, const bool acrossU, const Interval& across,
                               const Interval& along, const Interval& root, const ValleyHeight& at,
                               const double tolerance)
{
    const std::size_t a = acrossU ? 0 : 1;
    const std::size_t b = 1 - a;
    const Interval& haa = hessian[a][a];
    const Interval& hab = hessian[a][b];
    const Interval& hbb = hessian[b][b];
    if (!ExcludesZero(haa) || !hab.IsBounded() || !hbb.IsBounded()) {
        return false;
    }
    const Interval d{0.0, (Interval{along.GetUpper()} - Interval{along.GetLower()}).GetUpper()};
    const Interval low = haa * (Interval{across.GetLower()} - root) + hab * d;
    const Interval high = haa * (Interval{across.GetUpper()} - root) + hab * d;
    if (!ExcludesZero(low) || !ExcludesZero(high) || !((low * high).GetUpper() < 0.0)) {
        return false;
    }
    const Interval da = -(hab / haa) * d;
    const Interval height =
        at.Height + at.Slope * d +
        Interval{0.5} * (haa * Squared(da) + Interval{2.0} * hab * da * d + hbb * Squared(d));
    return height.GetUpper() <= tolerance && height.GetLower() >= -tolerance;
}

/// One proved box of a trace: its ends on the edges b = b0 and b = b1, the
/// boxes of S and C it settles, the plane of its frame and where the
/// valley leaves it, in doubles, for the next box.
struct ValleyBox {
    Seed From;
    Seed To;
    Box2 OnS;
    Box2 OnC;
    std::array<Vec3, 2> Plane;
    ValleyPoint End;
};

/// A seed at a root, in the pair's own order (S, C).
[[nodiscard]] Seed SeedAtRoot(const NodeSystem& sys, const bool acrossU, const IntervalsN<3>& x,
                              const double along)
{
    const std::array<double, 3> m = MidpointsOf<3>(x);
    const Point2 uv = OnSmall(acrossU, m[0], along);
    return {Middle(XyzOver(sys.Pair.A, Interval{uv.X}, Interval{uv.Y})), uv, Point2{m[1], m[2]}};
}

/// A patch clipped to a box of its parameters; nothing when a cut cannot be
/// enclosed or nothing is left.
[[nodiscard]] std::optional<Patch> ClippedTo(const Patch& p, const Box2& box)
{
    Result<std::optional<Patch>> u = Clipped(p, true, box.X);
    if (!u || !*u) {
        return std::nullopt;
    }
    Result<std::optional<Patch>> uv = Clipped(std::move(**u), false, box.Y);
    return uv && *uv ? std::move(*uv) : std::nullopt;
}

/// The valley's points at the ends and the middle of a step, in doubles.
struct ValleyGuesses {
    ValleyPoint First;
    ValleyPoint Middle;
    ValleyPoint Last;
};

/// The valley's points, in doubles, at the start, middle and end of a step.
[[nodiscard]] std::optional<ValleyGuesses> GuessesOver(const ValleyTrace& t,
                                                       const std::array<Vec3, 2>& e,
                                                       const ValleyPoint& start,
                                                       const Interval& along)
{
    const std::optional<ValleyPoint> first = ValleyAt(t, e, start, along.GetLower());
    const std::optional<ValleyPoint> middle =
        first ? ValleyAt(t, e, *first, Midpoint(along)) : std::nullopt;
    const std::optional<ValleyPoint> last =
        middle ? ValleyAt(t, e, *middle, along.GetUpper()) : std::nullopt;
    if (!last) {
        return std::nullopt;
    }
    return ValleyGuesses{*first, *middle, *last};
}

/// The interval widened outward to a dyadic grid of its stretch: each end
/// moved to the nearest multiple, from the stretch's lower end, of the
/// largest power of two no more than half the interval's width, and kept in
/// the stretch. The regions of consecutive steps then share most of their
/// sides, which the search cuts patches along: sides a rounding apart would
/// slice them into slivers. A larger box only asks more of the proof.
[[nodiscard]] Interval Snapped(const Interval& x, const Interval& stretch)
{
    const double width = Width(x);
    if (!(width > 0.0) || !std::isfinite(width)) {
        return x;
    }
    const double grid = std::ldexp(1.0, std::ilogb(width) - 1);
    const double origin = stretch.GetLower();
    const double lower = origin + std::floor((x.GetLower() - origin) / grid) * grid;
    const double upper = origin + std::ceil((x.GetUpper() - origin) / grid) * grid;
    return Intersection(Interval{std::min(lower, x.GetLower()), std::max(upper, x.GetUpper())},
                        stretch);
}

/// The hull of three numbers, inflated.
[[nodiscard]] Interval Around3(const double x, const double y, const double z, const double by)
{
    return Inflated(Interval{std::min({x, y, z}), std::max({x, y, z})}, by);
}

/// The boxes a step is proved on. S's: the step along b, and across a the
/// valley's guessed extent and, either side, as far again as the step is
/// long in space or as the band (BandAt), whichever is wider. C's: the
/// guessed feet and, either side, how far S's box reaches over C by the
/// projections' Jacobians at the middle. Guesses: the proof checks what it
/// needs of both boxes.
struct ValleyBoxes {
    Box2 OnS;
    Box2 OnC;
};

/// How many tolerances the surfaces part by at the band's sides.
constexpr double BandTolerances = 8.0;

/// The half-width across, in S's parameter a, of the band about the valley
/// past which the surfaces part by some BandTolerances tolerances: with h_aa
/// the relative normal curvature (II_S - II_C)(S_a, S_a) at the point,
/// h - h(valley) is about h_aa (a - phi)^2 / 2. The region a trace settles
/// spans the band, so that the pairs of patches left outside it are proved
/// apart; nearer, the surfaces are within the tolerance of each other, and
/// no pair there could be. A guess: the proof checks the box it sizes.
[[nodiscard]] std::optional<double> BandAt(const ValleyTrace& t, const ValleyJets& at)
{
    const auto m = [](const Bounds& b) { return MiddleVec(b); };
    const Vec3 cu = m(at.C.Pu);
    const Vec3 cv = m(at.C.Pv);
    const Vec3 normal = Cross(cu, cv);
    const double length = Length(normal);
    if (!(length > 0.0)) {
        return std::nullopt;
    }
    const Vec3 n = normal / length;
    const Vec3 sa = m(t.AcrossU ? at.S.Pu : at.S.Pv);
    const Vec3 saa = m(t.AcrossU ? at.S.Puu : at.S.Pvv);
    // S_a = alpha C_u + beta C_v in the tangent plane, by C's metric.
    const double e = Dot(cu, cu);
    const double f = Dot(cu, cv);
    const double g = Dot(cv, cv);
    const double det = e * g - f * f;
    if (!(det > 0.0)) {
        return std::nullopt;
    }
    const double x = Dot(sa, cu);
    const double y = Dot(sa, cv);
    const double alpha = (g * x - f * y) / det;
    const double beta = (e * y - f * x) / det;
    const double onC = alpha * alpha * Dot(m(at.C.Puu), n) +
                       2.0 * alpha * beta * Dot(m(at.C.Puv), n) + beta * beta * Dot(m(at.C.Pvv), n);
    const double curvature = std::abs(Dot(saa, n) - onC);
    if (!(curvature > 0.0) || !std::isfinite(curvature)) {
        return std::nullopt;
    }
    return std::sqrt(2.0 * BandTolerances * t.Tolerance / curvature);
}

/// S's box across: the guesses' extent and, either side, as far as the step
/// is long in space, or as the band where the box is wide.
[[nodiscard]] std::optional<Interval> AcrossOver(const ValleyTrace& t, const ValleyJets& at,
                                                 const ValleyGuesses& g, const Interval& along,
                                                 const bool wide)
{
    const Vec3 sa = MiddleVec(t.AcrossU ? at.S.Pu : at.S.Pv);
    const Vec3 sb = MiddleVec(t.AcrossU ? at.S.Pv : at.S.Pu);
    const double speed = Length(sa);
    const std::optional<double> band = speed > 0.0 ? BandAt(t, at) : std::nullopt;
    if (!band) {
        return std::nullopt;
    }
    const double step = Width(along) * Length(sb) / speed;
    const double reach = wide ? std::max(step, *band) : step;
    const Interval stretch = StretchOf(t.S, t.AcrossU);
    const Interval hull =
        Intersection(Around3(g.First.Across, g.Middle.Across, g.Last.Across, reach), stretch);
    if (hull.IsEmpty()) {
        return std::nullopt;
    }
    return wide ? Snapped(hull, stretch) : hull;
}

/// C's box: the guessed feet and, either side, how far S's box reaches over
/// C by the projections' Jacobians at the middle.
[[nodiscard]] std::optional<Box2> CoverOver(const ValleyTrace& t, const std::array<Vec3, 2>& e,
                                            const ValleyJets& at, const ValleyGuesses& g,
                                            const Box2& onS, const bool wide)
{
    const auto m = [](const Bounds& b) { return MiddleVec(b); };
    // dq = [e . C_u, e . C_v]; its inverse takes S's steps to C's.
    const double q00 = Dot(e[0], m(at.C.Pu));
    const double q01 = Dot(e[0], m(at.C.Pv));
    const double q10 = Dot(e[1], m(at.C.Pu));
    const double q11 = Dot(e[1], m(at.C.Pv));
    const double det = q00 * q11 - q01 * q10;
    if (!(std::abs(det) > 0.0)) {
        return std::nullopt;
    }
    const auto onC = [&](const Vec3& step, const double by) {
        const double x = Dot(e[0], step) * by;
        const double y = Dot(e[1], step) * by;
        return std::array{std::abs((q11 * x - q01 * y) / det), std::abs((q00 * y - q10 * x) / det)};
    };
    const std::array<double, 2> ru = onC(m(at.S.Pu), Width(onS.X));
    const std::array<double, 2> rv = onC(m(at.S.Pv), Width(onS.Y));
    const Interval cu = Intersection(
        Around3(g.First.OnC.X, g.Middle.OnC.X, g.Last.OnC.X, ru[0] + rv[0]), t.C.Parameters.X);
    const Interval cv = Intersection(
        Around3(g.First.OnC.Y, g.Middle.OnC.Y, g.Last.OnC.Y, ru[1] + rv[1]), t.C.Parameters.Y);
    if (cu.IsEmpty() || cv.IsEmpty()) {
        return std::nullopt;
    }
    return wide ? Box2{Snapped(cu, t.C.Parameters.X), Snapped(cv, t.C.Parameters.Y)} : Box2{cu, cv};
}

/// A step's boxes on S and C, narrow about the valley or wide across the band.
[[nodiscard]] std::optional<ValleyBoxes> BoxesOver(const ValleyTrace& t,
                                                   const std::array<Vec3, 2>& e,
                                                   const ValleyGuesses& g, const Interval& along,
                                                   const bool wide)
{
    const std::optional<ValleyJets> at = JetsAt(t, g.Middle);
    const std::optional<Interval> across = at ? AcrossOver(t, *at, g, along, wide) : std::nullopt;
    if (!across) {
        return std::nullopt;
    }
    const Box2 onS = t.AcrossU ? Box2{*across, along} : Box2{along, *across};
    const std::optional<Box2> onC = CoverOver(t, e, *at, g, onS, wide);
    if (!onC) {
        return std::nullopt;
    }
    return ValleyBoxes{onS, *onC};
}

/// A pair's derivative nets and its jets over the whole of each patch.
struct PairJets {
    SecondNets OnS;
    SecondNets OnC;
    Jet OverS;
    Jet OverC;
};

/// The pair's nets and jets; nothing when a jet cannot be enclosed.
[[nodiscard]] std::optional<PairJets> JetsOfPair(const PatchPair& pair)
{
    SecondNets small = SecondNetsOf(pair.A);
    SecondNets cover = SecondNetsOf(pair.B);
    std::optional<Jet> overS = JetOver(pair.A, small);
    std::optional<Jet> overC = JetOver(pair.B, cover);
    if (!overS || !overC) {
        return std::nullopt;
    }
    return PairJets{std::move(small), std::move(cover), std::move(*overS), std::move(*overC)};
}

/// Facts 1 and 2 on the wide pair, the band the step settles (S's box
/// across it, and C's covering it): C a graph, S over it, and h_aa of one
/// sign all over the band, so that each line across it holds at most one
/// zero of h_a. A Hessian's enclosure over the band is wide, so the valley
/// itself is proved on the narrow pair inside it (ProvedValley).
[[nodiscard]] bool ConvexAcross(const ValleyTrace& t, const PatchPair& pair, const NodeFrame& frame)
{
    const std::optional<PairJets> j = JetsOfPair(pair);
    if (!j) {
        return false;
    }
    const NodeSystem sys{pair, j->OnS, j->OnC, frame, j->OverS, j->OverC};
    const std::optional<Matrix2> inverse = Covered(sys) ? InverseProjection(sys) : std::nullopt;
    const std::size_t a = t.AcrossU ? 0 : 1;
    return inverse && ExcludesZero(HeightHessian(sys, *inverse)[a][a]);
}

/// Fact 3 at both ends of the step: each root isolated, strictly inside the
/// pair's boxes.
[[nodiscard]] std::optional<std::array<IntervalsN<3>, 2>>
RootsAtEnds(const NodeSystem& sys, const bool acrossU, const Interval& along)
{
    const Isolated<3> first = ValleyRoot(sys, acrossU, along.GetLower());
    const Isolated<3> last = ValleyRoot(sys, acrossU, along.GetUpper());
    const auto good = [&](const Isolated<3>& r) {
        return r.Outcome == Isolation::Unique && RootInside(sys, acrossU, r.Box);
    };
    if (!good(first) || !good(last)) {
        return std::nullopt;
    }
    return std::array{first.Box, last.Box};
}

/// Facts 1 to 5 on a clipped pair (S's box, C's box) in a frame, from b =
/// along's lower end to its upper: the two roots, or nothing. The pair is
/// the narrow one, inside the band ConvexAcross proved.
[[nodiscard]] std::optional<std::array<IntervalsN<3>, 2>>
ProvedValley(const ValleyTrace& t, const PatchPair& pair, const PairJets& j, const NodeFrame& frame,
             const Interval& along)
{
    const NodeSystem sys{pair, j.OnS, j.OnC, frame, j.OverS, j.OverC};
    const std::optional<Matrix2> inverse = Covered(sys) ? InverseProjection(sys) : std::nullopt;
    const std::optional<std::array<IntervalsN<3>, 2>> roots =
        inverse ? RootsAtEnds(sys, t.AcrossU, along) : std::nullopt;
    if (!roots) {
        return std::nullopt;
    }
    const std::optional<ValleyHeight> at =
        HeightAtRoot(sys, t.AcrossU, (*roots)[0], along.GetLower());
    if (!at || !ValleySpans(HeightHessian(sys, *inverse), t.AcrossU, StretchOf(pair.A, t.AcrossU),
                            along, (*roots)[0][0], *at, t.Tolerance)) {
        return std::nullopt;
    }
    return roots;
}

/// The narrow and the wide boxes of a step, the narrow inside the wide.
struct StepBoxes {
    ValleyBoxes Narrow;
    ValleyBoxes Wide;
};

/// Both pairs of boxes of a step, or nothing unless the narrow lie in the wide.
[[nodiscard]] std::optional<StepBoxes> StepBoxesOver(const ValleyTrace& t,
                                                     const std::array<Vec3, 2>& e,
                                                     const ValleyGuesses& g, const Interval& along)
{
    const std::optional<ValleyBoxes> narrow = BoxesOver(t, e, g, along, false);
    const std::optional<ValleyBoxes> wide = BoxesOver(t, e, g, along, true);
    if (!narrow || !wide || !Contains(wide->OnS, narrow->OnS) ||
        !Contains(wide->OnC, narrow->OnC)) {
        return std::nullopt;
    }
    return StepBoxes{*narrow, *wide};
}

/// The pair of patches clipped to a pair of boxes.
[[nodiscard]] std::optional<PatchPair> PairOver(const ValleyTrace& t, const ValleyBoxes& boxes)
{
    std::optional<Patch> small = ClippedTo(t.S, boxes.OnS);
    std::optional<Patch> cover = ClippedTo(t.C, boxes.OnC);
    if (!small || !cover) {
        return std::nullopt;
    }
    return PatchPair{std::move(*small), std::move(*cover)};
}

/// One step of a trace, from the valley's point `start` on b = along's
/// lower end to its upper end: the box proved, or nothing.
[[nodiscard]] std::optional<ValleyBox> ProvedStep(const ValleyTrace& t, const ValleyPoint& start,
                                                  const Interval& along)
{
    const std::optional<NodeFrame> frame =
        FrameAbout(NormalAt(t.C, t.OnC, start.OnC.X, start.OnC.Y));
    if (!frame) {
        return std::nullopt;
    }
    const std::array<Vec3, 2> e{MiddleVec(frame->E1), MiddleVec(frame->E2)};
    const std::optional<ValleyGuesses> guesses = GuessesOver(t, e, start, along);
    const std::optional<StepBoxes> boxes =
        guesses ? StepBoxesOver(t, e, *guesses, along) : std::nullopt;
    const std::optional<PatchPair> wide = boxes ? PairOver(t, boxes->Wide) : std::nullopt;
    const std::optional<PatchPair> pair =
        wide && ConvexAcross(t, *wide, *frame) ? PairOver(t, boxes->Narrow) : std::nullopt;
    const std::optional<PairJets> j = pair ? JetsOfPair(*pair) : std::nullopt;
    const std::optional<std::array<IntervalsN<3>, 2>> roots =
        j ? ProvedValley(t, *pair, *j, *frame, along) : std::nullopt;
    if (!roots) {
        return std::nullopt;
    }
    const NodeSystem sys{*pair, j->OnS, j->OnC, *frame, j->OverS, j->OverC};
    return ValleyBox{SeedAtRoot(sys, t.AcrossU, (*roots)[0], along.GetLower()),
                     SeedAtRoot(sys, t.AcrossU, (*roots)[1], along.GetUpper()),
                     boxes->Wide.OnS,
                     boxes->Wide.OnC,
                     e,
                     guesses->Last};
}

/// The valley's point on S's edge b = b0, walked back from a point of it in
/// doubles; nothing when the valley is lost on the way.
[[nodiscard]] std::optional<ValleyPoint> WalkedBack(const ValleyTrace& t, ValleyPoint p)
{
    const double target = StretchOf(t.S, !t.AcrossU).GetLower();
    const double from = p.Along;
    for (int k = 1; k <= ValleyWalk; ++k) {
        const double along = k == ValleyWalk ? target : from + (target - from) * k / ValleyWalk;
        const std::optional<ValleyPoint> next = ValleyNear(t, p, along);
        if (!next) {
            return std::nullopt;
        }
        p = *next;
    }
    return p;
}

/// The boxes of the valley across the whole of S's patch, each proved, the
/// step along b halved where one is not and doubled where one is; nothing
/// unless they reach from edge to edge.
[[nodiscard]] std::optional<std::vector<ValleyBox>> BoxesAcross(const ValleyTrace& t,
                                                                const ValleyPoint& start)
{
    const Interval range = StretchOf(t.S, !t.AcrossU);
    const double longest = 0.125 * Width(range);
    const double shortest = std::ldexp(Width(range), -ValleyHalvings);
    std::vector<ValleyBox> boxes;
    ValleyPoint p = start;
    double at = range.GetLower();
    double step = longest;
    while (at < range.GetUpper()) {
        if (boxes.size() >= ValleyBoxLimit || !(step >= shortest)) {
            return std::nullopt;
        }
        // No sliver at the end: a last step up to a quarter longer.
        const double to = at + 1.25 * step >= range.GetUpper() ? range.GetUpper() : at + step;
        std::optional<ValleyBox> box = ProvedStep(t, p, Interval{at, to});
        if (!box) {
            step *= 0.5;
            continue;
        }
        p = box->End;
        at = to;
        step = std::min(2.0 * step, longest);
        boxes.push_back(std::move(*box));
    }
    return boxes;
}

/// Whether nothing found so far lies in a box's region: no arc of a pair of
/// patches and no settled region shares area with it on both surfaces. A
/// branch already found there would otherwise be reported twice.
[[nodiscard]] bool NothingFoundIn(const Box2& onA, const Box2& onB, const Findings& found)
{
    const auto shares = [](const Box2& box, const std::vector<Box2>& boxes) {
        return std::ranges::any_of(boxes,
                                   [&](const Box2& other) { return HasAreaOverlap(box, other); });
    };
    const bool arc = std::ranges::any_of(found.Arcs, [&](const Arc& a) {
        return HasAreaOverlap(onA, a.BoxA) && HasAreaOverlap(onB, a.BoxB);
    });
    const bool region = std::ranges::any_of(found.Regions, [&](const NodeRegion& r) {
        return shares(onA, r.OnA) && shares(onB, r.OnB);
    });
    return !arc && !region;
}

/// The trace's boxes into arcs and settled regions. Each box's far end is
/// the next one's near end, found in the next one's frame, so that the arcs
/// join exactly.
void RecordValley(const ValleyTrace& t, std::vector<ValleyBox> boxes, Findings& found)
{
    for (std::size_t k = 0; k + 1 < boxes.size(); ++k) {
        boxes[k].To = boxes[k + 1].From;
    }
    const auto ordered = [&](const Seed& seed) { return t.SmallIsA ? seed : Swapped(seed); };
    for (const ValleyBox& box : boxes) {
        const Box2& onA = t.SmallIsA ? box.OnS : box.OnC;
        const Box2& onB = t.SmallIsA ? box.OnC : box.OnS;
        found.Arcs.push_back({ordered(box.From), ordered(box.To), Vec3{}, onA, onB,
                              TangentChart{t.SmallIsA, t.AcrossU, box.Plane[0], box.Plane[1]}});
        found.Regions.push_back({{onA}, {onB}});
    }
    found.Tangent = std::array{t.SmallIsA, t.AcrossU};
}

/// A trace from a point of the valley: walked back to S's edge, proved box
/// by box to the other edge, and recorded when nothing was found there.
[[nodiscard]] bool TraceValley(const ValleyTrace& t, const ValleyPoint& at, Findings& found)
{
    const std::optional<ValleyPoint> start = WalkedBack(t, at);
    std::optional<std::vector<ValleyBox>> boxes = start ? BoxesAcross(t, *start) : std::nullopt;
    if (!boxes || boxes->empty()) {
        return false;
    }
    for (const ValleyBox& box : *boxes) {
        if (!NothingFoundIn(t.SmallIsA ? box.OnS : box.OnC, t.SmallIsA ? box.OnC : box.OnS,
                            found)) {
            return false;
        }
    }
    RecordValley(t, std::move(*boxes), found);
    return true;
}

/// Where a pair's cells (S's, C's) may hold a valley within the tolerance:
/// its point on the middle line of S's cell, in the cells, the surfaces
/// there within the tolerance of each other and h's slope along b no more
/// than would part them by the tolerance over the patch. A guess, to spare
/// a trace where there is no tangency: a transverse pair has a valley too,
/// along which the surfaces part.
/// Whether the surfaces at a point of the valley are within the tolerance
/// of each other and h's slope along b would not part them by it over S's
/// patch.
[[nodiscard]] bool NearTangent(const ValleyTrace& t, const ValleyPoint& p)
{
    const std::optional<ValleyJets> at = JetsAt(t, p);
    if (!at) {
        return false;
    }
    const Vec3 normal = Cross(MiddleVec(at->C.Pu), MiddleVec(at->C.Pv));
    const double length = Length(normal);
    const double gap = Distance(Middle(at->S.P), Middle(at->C.P));
    if (!(length > 0.0) || !(gap <= t.Tolerance)) {
        return false;
    }
    const double slope = std::abs(Dot(MiddleVec(t.AcrossU ? at->S.Pv : at->S.Pu), normal)) / length;
    return slope * Width(StretchOf(t.S, !t.AcrossU)) <= t.Tolerance;
}

/// A candidate for the valley in a pair's cells (see NearTangent), or nothing.
[[nodiscard]] std::optional<ValleyPoint> ValleyCandidate(const ValleyTrace& t, const Box2& cellS,
                                                         const Box2& cellC)
{
    const Interval across = t.AcrossU ? cellS.X : cellS.Y;
    const Interval along = t.AcrossU ? cellS.Y : cellS.X;
    const ValleyPoint guess{Midpoint(across), Midpoint(along),
                            Point2{Midpoint(cellC.X), Midpoint(cellC.Y)}};
    const std::optional<ValleyPoint> p = ValleyNear(t, guess, guess.Along);
    const bool inside = p && Contains(across, p->Across) && Contains(cellC, p->OnC);
    return inside && NearTangent(t, *p) ? p : std::nullopt;
}

/// Whether a trace was tried over the same patches and failed.
[[nodiscard]] bool FailedBefore(const Findings& found, const FailedTrace& tried)
{
    return std::ranges::any_of(found.FailedTraces, [&](const FailedTrace& f) {
        return f.OnS == tried.OnS && f.OnC == tried.OnC && f.SmallIsA == tried.SmallIsA &&
               f.AcrossU == tried.AcrossU;
    });
}

/// The Bezier patches of S and C the pair's cells lie in, S the one of A
/// when smallIsA.
struct TracePatches {
    const Patch* S;
    const Patch* C;
};

/// The Bezier patches of the searched atlases that hold the pair's cells.
[[nodiscard]] std::optional<TracePatches>
PatchesHolding(const PatchPair& pair, const SearchSettings& s, const bool smallIsA)
{
    const Patch* small =
        Holding(smallIsA ? *s.OnA : *s.OnB, (smallIsA ? pair.A : pair.B).Parameters);
    const Patch* cover =
        Holding(smallIsA ? *s.OnB : *s.OnA, (smallIsA ? pair.B : pair.A).Parameters);
    if (small == nullptr || cover == nullptr) {
        return std::nullopt;
    }
    return TracePatches{small, cover};
}

/// One trace, S's parameter `acrossU` across, from a candidate in the
/// pair's cells; a failure is remembered.
[[nodiscard]] bool TriedTrace(const PatchPair& pair, const ValleyTrace& t, Findings& found)
{
    const FailedTrace tried{t.S.Parameters, t.C.Parameters, t.SmallIsA, t.AcrossU};
    if ((found.Tangent && *found.Tangent != std::array{t.SmallIsA, t.AcrossU}) ||
        FailedBefore(found, tried)) {
        return false;
    }
    const Box2& cellS = (t.SmallIsA ? pair.A : pair.B).Parameters;
    const Box2& cellC = (t.SmallIsA ? pair.B : pair.A).Parameters;
    const std::optional<ValleyPoint> at = ValleyCandidate(t, cellS, cellC);
    if (!at) {
        return false;
    }
    if (TraceValley(t, *at, found)) {
        return true;
    }
    found.FailedTraces.push_back(tried);
    return false;
}

/// A pair with no advance direction that may lie on a curve of tangency:
/// the valley traced across the Bezier patch of S the pair's cell lies in,
/// with either surface as S and either parameter across, as the search's
/// first trace fixed them. The pair goes back on the stack, to be carved by
/// the regions the trace settled.
[[nodiscard]] bool SettleTangent(const PatchPair& pair, const SearchSettings& s, Findings& found,
                                 std::vector<PatchPair>& stack)
{
    if (s.OnA == nullptr || s.OnB == nullptr) {
        return false;
    }
    for (const bool smallIsA : {true, false}) {
        const std::optional<TracePatches> patches = PatchesHolding(pair, s, smallIsA);
        if (!patches) {
            continue;
        }
        const SecondNets onS = SecondNetsOf(*patches->S);
        const SecondNets onC = SecondNetsOf(*patches->C);
        for (const bool acrossU : {true, false}) {
            const ValleyTrace t{*patches->S, onS,     *patches->C, onC,
                                smallIsA,    acrossU, s.Tolerance, s.Join};
            if (TriedTrace(pair, t, found)) {
                stack.push_back(pair);
                return true;
            }
        }
    }
    return false;
}

} // namespace
} // namespace Stereon::Detail
