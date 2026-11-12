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
#include "BezierPatch.hpp"

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
///   the tolerance (found, not yet settled: it fails closed);
/// - a transverse node, polynomial surfaces only: a critical point of one
///   surface's height over the other, isolated by Krawczyk's operator within
///   the tolerance, around which a Lyapunov function proved on the Hessian's
///   enclosure rules out every other piece of the intersection, so that the
///   four roots on the pair's boundary are the ends of the four
///   half-branches leaving it (see "Branch points" below).
///
/// Arcs join at shared seeds into branches, every free end on a boundary of
/// the domains or at a node, where four branches end. Each arc is marched by Newton's method on S_A
/// = S_B and L . S_A = advance, halving chords until they hold, and the polylines in space and on
/// both surfaces are certified by CertifyDeviation, tightened until they pass. Polylines, because a
/// fitted cubic's pcurve can leave a surface's domain.
///
/// What cannot be settled fails, never a partial or an unproved empty result:
/// a tangency along a curve, a branch point that is not a proved transverse
/// node (a higher-order singularity, a node on a line between Bezier patches,
/// a node of rational surfaces), an overlap, and any tangency of rational
/// surfaces, at the smallest pair, IntersectionFailed; a search past its
/// budget, the same.
namespace Stereon::Detail {
namespace {

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

/// The pairs of patches one search may visit. Separating surfaces a gap g
/// apart costs about (extent / g)^2 pairs, so close ones fail rather than
/// run on.
constexpr std::size_t PairBudget = std::size_t{1} << 14;

/// Newton's iterations for a point, and the depth of a chord's refinement.
constexpr int NewtonLimit = 16;
constexpr int RefineDepth = 26;
/// How often a branch's chords are tightened before certification fails.
constexpr int CertifyRetries = 5;

struct PatchPair {
    Patch A;
    Patch B;
};

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

/// The box of a line of poles, as ExtentOf boxes a net: the convex hull of
/// the edge it controls.
[[nodiscard]] Bounds ExtentOfLine(const std::vector<HPoint>& line)
{
    return ExtentOf(Net{line.size(), 1, line});
}

/// An edge's derivative along itself, enclosed from its own pole line as
/// HullsOf encloses a patch's, met with the patch's: the patch's alone can
/// hold zero, across a tangency inside it, where the edge's does not.
[[nodiscard]] Bounds EdgeDerivative(const std::vector<HPoint>& line, const Interval& width,
                                    const Bounds& patchWide)
{
    const Net net{line.size(), 1, line};
    const HPoint h = HullOf(net);
    if (!(h.W.GetLower() > 0.0 && width.GetLower() > 0.0)) {
        return patchWide;
    }
    const Bounds own = QuotientOf(HullOf(DerivativeNet(net, true, width)), ExtentOf(net), h.W);
    Bounds met;
    for (std::size_t i = 0; i < met.size(); ++i) {
        met[i] = Intersection(own[i], patchWide[i]);
        if (met[i].IsEmpty()) {
            return patchWide; // Both hold the derivative; only rounding parts them.
        }
    }
    return met;
}

/// The pole line controlling one edge of a patch, the lower or the upper.
[[nodiscard]] std::vector<HPoint> EdgeLine(const Patch& carrier, const bool alongU,
                                           const bool upper)
{
    return LineOf(carrier.Poles, alongU, upper ? LineCount(carrier.Poles, alongU) - 1 : 0);
}

/// The edge's root against the other patch: none at once when the edge's
/// hull misses the other patch's.
[[nodiscard]] Isolated<3> RootOnEdge(const Edge& edge, const std::vector<HPoint>& line,
                                     const Hulls& own, const Patch& other, const Hulls& theirs)
{
    if (ProvedApart(ExtentOfLine(line), other.Extent, 0.0)) {
        return {Isolation::None, {}};
    }
    const Bounds along = EdgeDerivative(line, WidthOf(StretchOf(edge.Carrier, edge.AlongU)),
                                        edge.AlongU ? own.Du : own.Dv);
    return EdgeRoot(edge, along, other, theirs);
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
        for (const bool upper : {false, true}) {
            const Edge edge{carrier, alongU, upper ? across.GetUpper() : across.GetLower()};
            const Isolated<3> root =
                RootOnEdge(edge, EdgeLine(carrier, alongU, upper), own, other, theirs);
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
[[nodiscard]] std::optional<Form2> RelativeForm(const LocalForm& a, const LocalForm& b)
{
    const double e = Dot(b.Du, b.Du);
    const double f = Dot(b.Du, b.Dv);
    const double g = Dot(b.Dv, b.Dv);
    const double det = e * g - f * f;
    if (!(det > 0.0)) {
        return std::nullopt;
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
    return Form2{a.Second[0] - pulled(p, p), a.Second[1] - pulled(p, q),
                 a.Second[2] - pulled(q, q)};
}

/// The relative form's determinant, against the margin rounding could
/// move it by: positive for a touch, negative for a crossing of two
/// branches.
[[nodiscard]] int FormSign(const LocalForm& a, const LocalForm& b)
{
    const std::optional<Form2> m = RelativeForm(a, b);
    if (!m) {
        return 0;
    }
    const auto [m11, m12, m22] = *m;
    const double scale = m11 * m11 + 2.0 * m12 * m12 + m22 * m22;
    const double det = m11 * m22 - m12 * m12;
    if (det > 0x1p-20 * scale) {
        return 1;
    }
    return det < -0x1p-20 * scale ? -1 : 0;
}

[[nodiscard]] bool RelativelyDefinite(const LocalForm& a, const LocalForm& b)
{
    return FormSign(a, b) > 0;
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
    /// Proved transverse nodes, where four arcs end.
    std::vector<Point3> Nodes;
};

struct SearchSettings {
    double Tolerance;
    double Join;
    bool Polynomial;
};

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

/// Whether the parameter lies on a side of its domain, within the join.
[[nodiscard]] bool OnSide(const double x, const Interval& domain)
{
    const double slack = 0x1p-30 * std::max(1.0, Width(domain));
    return x <= domain.GetLower() + slack || x >= domain.GetUpper() - slack;
}

// --- Branch points -----------------------------------------------------------
//
// A transverse node: the surfaces share a tangent plane at a point P and the
// relative second form there is indefinite and nondegenerate, so two smooth
// branches cross at P. It is settled on a pair of patches, a small one S and
// a covering one C, by four facts proved in interval arithmetic:
//
// 1. Coordinates. A frame (e1, e2, n), n near C's normal at P, projects C's
//    patch one to one onto the plane: the projection's Jacobian is regular
//    over the whole (convex) parameter box, so by the mean value theorem no
//    two parameters share a projection. C is then the graph z = g(p) over a
//    closed disc D, bounded by the projection of C's edges.
// 2. Containment. The projection of S's hull misses the projection of every
//    edge of C, and meets D at the critical point below, so it lies inside D.
//    The height h(s) = n . S(s) - g(e . S(s)) is then defined on S's whole
//    box, and its zeros are exactly the points of S on C (C's parameter
//    being determined).
// 3. The critical point. Krawczyk's operator isolates the unique solution in
//    the pair of  e . (S - C) = 0, S_u . N_C = S_v . N_C = 0,  a critical
//    point s* of h; the surfaces there are within the tolerance of each other.
// 4. No other structure. The Hessian of h is enclosed over S's box: with nu =
//    n - G1 e1 - G2 e2 for G the gradient of g,
//        Hess h = [nu . S_jk] - Dp^T Hess g Dp,  Hess g = Dq^-T [nu . C_jk] Dq^-1,
//    Dp and Dq the projections' Jacobians. Along the zero set, h's level
//    curve runs along J grad h, J a quarter turn, and grad h(s) = H (s - s*)
//    for H in that enclosure (mean value on the convex box, grad h(s*) = 0).
//    With M = sym(J H0), H0 the enclosure's middle, the symmetric part of
//    M J H is proved positive definite for every H in it, so
//    Q(s) = (s - s*)^T M (s - s*) strictly increases along every level curve
//    away from s*. (M J H0 has a positive definite symmetric part for any
//    indefinite nonsingular H0: J H0 is traceless with (J H0)^2 = -det(H0) I,
//    so sym(J H0) J H0 + (J H0)^T sym(J H0) = (J H0)^T (J H0) - det(H0) I.)
//    Hence h has no other critical point, the zero set has no closed loop,
//    and no arc leaves s* and returns to it; H is indefinite, so four
//    half-branches leave s* (Morse), each ending on S's boundary.
//
// Every arc of the zero set therefore ends on S's edges or at s*, and the
// edges' roots, each isolated by Krawczyk's operator, count them: four are
// the four half-branches and nothing else, s* joining them. Fewer or more,
// or anything unproved, and the node is not settled here.
//
// The proof places P within the tolerance of both surfaces, not on them: a
// gap within it is not told apart from a crossing (ADR-0050), and the node is
// reported as the shared end of the four branches, a candidate for an
// arrangement event, never as one by itself (ADR-0054).

/// A 2 x 2 matrix of intervals, by rows.
using Matrix2 = std::array<std::array<Interval, 2>, 2>;

[[nodiscard]] Matrix2 Times(const Matrix2& a, const Matrix2& b)
{
    Matrix2 c{};
    for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t k = 0; k < 2; ++k) {
            c[r][k] = a[r][0] * b[0][k] + a[r][1] * b[1][k];
        }
    }
    return c;
}

[[nodiscard]] Matrix2 Transposed(const Matrix2& a)
{
    return {{{a[0][0], a[1][0]}, {a[0][1], a[1][1]}}};
}

[[nodiscard]] bool ExcludesZero(const Interval& x)
{
    return x.GetLower() > 0.0 || x.GetUpper() < 0.0;
}

/// Linear coordinates for a node: across the plane (e1, e2), up it n. They
/// need not be orthonormal, only a basis, which rounding cannot spoil.
struct NodeFrame {
    Bounds E1;
    Bounds E2;
    Bounds N;
};

[[nodiscard]] Bounds Exactly(const Vec3& v)
{
    return {Interval{v.X}, Interval{v.Y}, Interval{v.Z}};
}

[[nodiscard]] std::optional<NodeFrame> FrameAbout(const Vec3& normal)
{
    const double length = Length(normal);
    if (!(length > 0.0) || !std::isfinite(length)) {
        return std::nullopt;
    }
    const Vec3 n = normal / length;
    const double ax = std::abs(n.X);
    const double ay = std::abs(n.Y);
    const double az = std::abs(n.Z);
    const Vec3 axis = ax <= ay && ax <= az ? Vec3{1.0, 0.0, 0.0}
                      : ay <= az           ? Vec3{0.0, 1.0, 0.0}
                                           : Vec3{0.0, 0.0, 1.0};
    const Vec3 e1 = Normalized(Cross(n, axis));
    const Vec3 e2 = Cross(n, e1);
    return NodeFrame{Exactly(e1), Exactly(e2), Exactly(n)};
}

/// The Jacobian of a patch's projection onto the frame's plane: rows e1
/// and e2, columns the partials along u and v.
[[nodiscard]] Matrix2 ProjectedJacobian(const NodeFrame& f, const Bounds& du, const Bounds& dv)
{
    return {{{DotOf(f.E1, du), DotOf(f.E1, dv)}, {DotOf(f.E2, du), DotOf(f.E2, dv)}}};
}

/// A small patch S and a patch C meant to cover it, as a pair (S, C).
struct NodeSystem {
    const PatchPair& Pair;
    const SecondNets& S;
    const SecondNets& C;
    const NodeFrame& Frame;
};

/// The critical point of S's height over C, solved for (u_S, v_S, u_C, v_C):
///     e1 . (S - C) = 0,  e2 . (S - C) = 0,  S_u . N_C = 0,  S_v . N_C = 0,
/// N_C = C_u x C_v. Polynomial patches only.
[[nodiscard]] Isolated<4> CriticalPoint(const NodeSystem& sys)
{
    const Patch& small = sys.Pair.A;
    const Patch& cover = sys.Pair.B;
    const NormalJet n = NormalHull(sys.C);
    const Bounds su = XyzHull(sys.S.Du);
    const Bounds sv = XyzHull(sys.S.Dv);
    const Bounds cu = XyzHull(sys.C.Du);
    const Bounds cv = XyzHull(sys.C.Dv);
    const NodeFrame& f = sys.Frame;
    IntervalMatrix<4> j{};
    j[0] = {DotOf(f.E1, su), DotOf(f.E1, sv), -DotOf(f.E1, cu), -DotOf(f.E1, cv)};
    j[1] = {DotOf(f.E2, su), DotOf(f.E2, sv), -DotOf(f.E2, cu), -DotOf(f.E2, cv)};
    j[2] = {DotOf(XyzHull(sys.S.Duu), n.N), DotOf(XyzHull(sys.S.Duv), n.N), DotOf(su, n.Nu),
            DotOf(su, n.Nv)};
    j[3] = {DotOf(XyzHull(sys.S.Duv), n.N), DotOf(XyzHull(sys.S.Dvv), n.N), DotOf(sv, n.Nu),
            DotOf(sv, n.Nv)};
    const IntervalsN<4> box{small.Parameters.X, small.Parameters.Y, cover.Parameters.X,
                            cover.Parameters.Y};
    return Krawczyk<4>(box, j, [&](const std::array<double, 4>& m) {
        const Bounds r = MinusOf(XyzValue(small, m[0], m[1]), XyzValue(cover, m[2], m[3]));
        const Bounds normal = CrossOf(XyzAt(sys.C.Du, cover.Parameters, m[2], m[3]),
                                      XyzAt(sys.C.Dv, cover.Parameters, m[2], m[3]));
        return IntervalsN<4>{DotOf(f.E1, r), DotOf(f.E2, r),
                             DotOf(XyzAt(sys.S.Du, small.Parameters, m[0], m[1]), normal),
                             DotOf(XyzAt(sys.S.Dv, small.Parameters, m[0], m[1]), normal)};
    });
}

/// The patch's points over enclosed parameters, polynomial.
[[nodiscard]] Bounds XyzOver(const Patch& p, const Interval& u, const Interval& v)
{
    const HPoint h = ValueOf(p.Poles, p.Parameters, u, v);
    return {h.X, h.Y, h.Z};
}

/// How far apart the two surfaces can be at the critical point.
[[nodiscard]] Interval GapAt(const PatchPair& pair, const IntervalsN<4>& x)
{
    return DistanceBetween(MinusOf(XyzOver(pair.A, x[0], x[1]), XyzOver(pair.B, x[2], x[3])));
}

/// e . P over a polynomial pole line or net, by the convex hull of its poles.
[[nodiscard]] Interval Across(const std::vector<HPoint>& poles, const Bounds& e)
{
    Interval range;
    for (const HPoint& h : poles) {
        range = Hull(range, DotOf(e, Bounds{h.X, h.Y, h.Z}));
    }
    return range;
}

/// Whether the projection of a pole line, cut in halves to a depth, misses
/// the box x by y of the plane.
[[nodiscard]] bool ProjectedApart(const std::vector<HPoint>& line, const NodeFrame& f,
                                  const Interval& x, const Interval& y, const int depth)
{
    if (Intersection(Across(line, f.E1), x).IsEmpty() ||
        Intersection(Across(line, f.E2), y).IsEmpty()) {
        return true;
    }
    if (depth == 0) {
        return false;
    }
    std::vector<HPoint> lower;
    std::vector<HPoint> upper;
    SplitLine(line, Interval{0.5}, lower, upper);
    return ProjectedApart(lower, f, x, y, depth - 1) && ProjectedApart(upper, f, x, y, depth - 1);
}

/// Facts 1 and 2: C projects one to one, and S's projection keeps clear of
/// the projection of C's edges.
[[nodiscard]] bool Covered(const NodeSystem& sys)
{
    const Patch& small = sys.Pair.A;
    const Patch& cover = sys.Pair.B;
    const Matrix2 dq = ProjectedJacobian(sys.Frame, XyzHull(sys.C.Du), XyzHull(sys.C.Dv));
    if (!ExcludesZero(dq[0][0] * dq[1][1] - dq[0][1] * dq[1][0])) {
        return false;
    }
    const Interval x = Across(small.Poles.Poles, sys.Frame.E1);
    const Interval y = Across(small.Poles.Poles, sys.Frame.E2);
    constexpr int depth = 6;
    for (const bool alongU : {true, false}) {
        const std::size_t last = LineCount(cover.Poles, alongU) - 1;
        for (const std::size_t line : {std::size_t{0}, last}) {
            if (!ProjectedApart(LineOf(cover.Poles, alongU, line), sys.Frame, x, y, depth)) {
                return false;
            }
        }
    }
    return true;
}

/// The inverse of C's projection Jacobian Dq, enclosed over C's box: each
/// entry an adjugate entry over the determinant, which excludes zero.
[[nodiscard]] std::optional<Matrix2> InverseProjection(const NodeSystem& sys)
{
    const Matrix2 dq = ProjectedJacobian(sys.Frame, XyzHull(sys.C.Du), XyzHull(sys.C.Dv));
    const Interval det = dq[0][0] * dq[1][1] - dq[0][1] * dq[1][0];
    if (!ExcludesZero(det)) {
        return std::nullopt;
    }
    return Matrix2{{{dq[1][1] / det, -dq[0][1] / det}, {-dq[1][0] / det, dq[0][0] / det}}};
}

/// Fact 4's enclosure of the Hessian of S's height over C, on S's box.
[[nodiscard]] Matrix2 HeightHessian(const NodeSystem& sys, const Matrix2& inverse)
{
    const NodeFrame& f = sys.Frame;
    // G = Dq^-T (n . C_u, n . C_v), the gradient of C's graph.
    const Interval w1 = DotOf(f.N, XyzHull(sys.C.Du));
    const Interval w2 = DotOf(f.N, XyzHull(sys.C.Dv));
    const Interval g1 = inverse[0][0] * w1 + inverse[1][0] * w2;
    const Interval g2 = inverse[0][1] * w1 + inverse[1][1] * w2;
    const Bounds nu{f.N[0] - g1 * f.E1[0] - g2 * f.E2[0], f.N[1] - g1 * f.E1[1] - g2 * f.E2[1],
                    f.N[2] - g1 * f.E1[2] - g2 * f.E2[2]};
    const auto second = [&](const SecondNets& s) {
        const Bounds uv = XyzHull(s.Duv);
        return Matrix2{{{DotOf(nu, XyzHull(s.Duu)), DotOf(nu, uv)},
                        {DotOf(nu, uv), DotOf(nu, XyzHull(s.Dvv))}}};
    };
    const Matrix2 graph = Times(Transposed(inverse), Times(second(sys.C), inverse));
    const Matrix2 dp = ProjectedJacobian(f, XyzHull(sys.S.Du), XyzHull(sys.S.Dv));
    const Matrix2 pulled = Times(Transposed(dp), Times(graph, dp));
    const Matrix2 own = second(sys.S);
    return Matrix2{{{own[0][0] - pulled[0][0], own[0][1] - pulled[0][1]},
                    {own[1][0] - pulled[1][0], own[1][1] - pulled[1][1]}}};
}

/// Fact 4's Lyapunov test: sym(M J H) positive definite over the enclosure.
[[nodiscard]] bool LevelCurvesLeave(const Matrix2& h)
{
    const double h00 = Midpoint(h[0][0]);
    const double h11 = Midpoint(h[1][1]);
    const double h01 = 0.5 * (Midpoint(h[0][1]) + Midpoint(h[1][0]));
    if (!std::isfinite(h00) || !std::isfinite(h11) || !std::isfinite(h01)) {
        return false;
    }
    // J H0 = [[-h01, -h11], [h00, h01]]; M is its symmetric part.
    const Interval off{0.5 * (h00 - h11)};
    const Matrix2 m{{{Interval{-h01}, off}, {off, Interval{h01}}}};
    const Matrix2 jh{{{-h[1][0], -h[1][1]}, {h[0][0], h[0][1]}}};
    const Matrix2 p = Times(m, jh);
    const Interval s01 = Interval{0.5} * (p[0][1] + p[1][0]);
    return p[0][0].GetLower() > 0.0 && (p[0][0] * p[1][1] - s01 * s01).GetLower() > 0.0;
}

/// Whether the enclosure lies strictly inside the stretch.
[[nodiscard]] bool Inside(const Interval& x, const Interval& stretch)
{
    return stretch.GetLower() < x.GetLower() && x.GetUpper() < stretch.GetUpper();
}

/// A settled node, in the pair's own order (small, cover): the node and the
/// four ends of its half-branches on the small patch's edges.
struct NodeEnds {
    Seed Node;
    std::vector<Seed> Ends;
};

/// What the proof knows of the critical point: its box (s*, c*), the height
/// h(s*) = n . (S - C) there, h's Hessian over S's box and Dq^-1 over C's.
struct NodeAt {
    IntervalsN<4> Box;
    Interval Height;
    Matrix2 Hessian;
    Matrix2 Inverse;
};

/// The deepest an edge of the small patch is halved while its roots are
/// counted.
constexpr int EdgeDepth = 12;

/// One edge of the small patch near a node, and the node's projection.
struct NodeEdge {
    const NodeSystem& Sys;
    const NodeAt& At;
    bool AlongU;
    bool Upper;
    double Fixed;
    Interval NodeX;
    Interval NodeY;
};

/// Whether Taylor's theorem keeps the height off zero on a piece of the
/// edge: h(s) in h(s*) + (s - s*)^T H (s - s*) / 2, grad h(s*) = 0.
[[nodiscard]] bool HeightApart(const NodeEdge& e, const Interval& piece)
{
    const Interval du = (e.AlongU ? piece : Interval{e.Fixed}) - e.At.Box[0];
    const Interval dv = (e.AlongU ? Interval{e.Fixed} : piece) - e.At.Box[1];
    const Matrix2& h = e.At.Hessian;
    return ExcludesZero(e.At.Height +
                        Interval{0.5} * (h[0][0] * Squared(du) + (h[0][1] + h[1][0]) * du * dv +
                                         h[1][1] * Squared(dv)));
}

/// The part of C that can hold a root on a piece of the edge, by the mean
/// value theorem: q(c) - q(c*) = M (c - c*) with M in Dq's enclosure, so
/// c lies in c* + Dq^-1 (p(piece) - p(s*)), as p(s*) = q(c*) exactly.
[[nodiscard]] std::optional<Patch> CoverReach(const NodeEdge& e, const std::vector<HPoint>& line)
{
    const Patch& cover = e.Sys.Pair.B;
    const Matrix2& inverse = e.At.Inverse;
    const Interval px = Across(line, e.Sys.Frame.E1) - e.NodeX;
    const Interval py = Across(line, e.Sys.Frame.E2) - e.NodeY;
    const Interval u =
        Intersection(e.At.Box[2] + inverse[0][0] * px + inverse[0][1] * py, cover.Parameters.X);
    const Interval v =
        Intersection(e.At.Box[3] + inverse[1][0] * px + inverse[1][1] * py, cover.Parameters.Y);
    if (u.IsEmpty() || v.IsEmpty()) {
        return std::nullopt;
    }
    std::optional<Patch> part = Clipped(cover, true, u);
    if (part) {
        part = Clipped(*part, false, v);
    }
    return part;
}

/// The root on one piece of the edge: none when the height keeps off zero
/// there or no part of C can hold one, else Krawczyk's operator's verdict.
[[nodiscard]] Isolation PieceRoot(const NodeEdge& e, const Interval& piece, const double join,
                                  std::vector<Seed>& seeds)
{
    if (HeightApart(e, piece)) {
        return Isolation::None;
    }
    const std::optional<Patch> carrier = Clipped(e.Sys.Pair.A, e.AlongU, piece);
    if (!carrier) {
        return Isolation::None;
    }
    const std::vector<HPoint> line = EdgeLine(*carrier, e.AlongU, e.Upper);
    const std::optional<Patch> part = CoverReach(e, line);
    if (!part) {
        return Isolation::None;
    }
    const Edge edge{*carrier, e.AlongU, e.Fixed};
    const Isolated<3> root = RootOnEdge(edge, line, HullsOf(*carrier), *part, HullsOf(*part));
    if (root.Outcome == Isolation::Unique) {
        AddSeed(seeds, SeedOfRoot(edge, root, *part, true), join);
    }
    return root.Outcome;
}

/// The roots of one edge of the small patch against C, each isolated by
/// Krawczyk's operator on a piece of the edge against the part of C that
/// can hold it, a piece halved while it is not settled; false when one is
/// still not settled at the deepest halving.
[[nodiscard]] bool NodeEdgeRoots(const NodeSystem& sys, const NodeAt& at, const bool alongU,
                                 const bool upper, const double join, std::vector<Seed>& seeds)
{
    const Interval across = StretchOf(sys.Pair.A, !alongU);
    const Bounds node = XyzOver(sys.Pair.A, at.Box[0], at.Box[1]);
    const NodeEdge e{sys,
                     at,
                     alongU,
                     upper,
                     upper ? across.GetUpper() : across.GetLower(),
                     DotOf(sys.Frame.E1, node),
                     DotOf(sys.Frame.E2, node)};
    std::vector<std::pair<Interval, int>> work{{StretchOf(sys.Pair.A, alongU), 0}};
    while (!work.empty()) {
        const auto [piece, depth] = work.back();
        work.pop_back();
        if (PieceRoot(e, piece, join, seeds) != Isolation::Unsettled) {
            continue;
        }
        if (depth >= EdgeDepth) {
            return false;
        }
        const double middle = Midpoint(piece);
        work.emplace_back(Interval{piece.GetLower(), middle}, depth + 1);
        work.emplace_back(Interval{middle, piece.GetUpper()}, depth + 1);
    }
    return true;
}

/// Fact 3: the critical point isolated, strictly inside both boxes, the
/// surfaces there within the tolerance.
[[nodiscard]] std::optional<IntervalsN<4>> ProvedCritical(const NodeSystem& sys,
                                                          const double tolerance)
{
    const Isolated<4> critical = CriticalPoint(sys);
    const Box2& small = sys.Pair.A.Parameters;
    const Box2& cover = sys.Pair.B.Parameters;
    const IntervalsN<4>& x = critical.Box;
    if (critical.Outcome != Isolation::Unique || !Inside(x[0], small.X) || !Inside(x[1], small.Y) ||
        !Inside(x[2], cover.X) || !Inside(x[3], cover.Y) ||
        GapAt(sys.Pair, x).GetUpper() > tolerance) {
        return std::nullopt;
    }
    return x;
}

/// The four ends on S's edges, or nothing when they are not exactly four or
/// one is not settled. By fact 2 no root lies on C's edges.
[[nodiscard]] std::optional<std::vector<Seed>> NodeEndsOf(const NodeSystem& sys, const NodeAt& at,
                                                          const double join)
{
    std::vector<Seed> ends;
    for (const bool alongU : {true, false}) {
        for (const bool upper : {false, true}) {
            if (!NodeEdgeRoots(sys, at, alongU, upper, join, ends)) {
                return std::nullopt;
            }
        }
    }
    if (ends.size() != 4) {
        return std::nullopt;
    }
    return ends;
}

/// Facts 1 to 4 on a pair (small, cover), with the tolerance the surfaces
/// may be apart at the node and the join of seeds.
[[nodiscard]] std::optional<NodeEnds> ProvedNode(const PatchPair& pair, const Vec3& normal,
                                                 const SearchSettings& s)
{
    const std::optional<NodeFrame> frame = FrameAbout(normal);
    if (!frame) {
        return std::nullopt;
    }
    const SecondNets small = SecondNetsOf(pair.A);
    const SecondNets cover = SecondNetsOf(pair.B);
    const NodeSystem sys{pair, small, cover, *frame};
    const std::optional<IntervalsN<4>> x = ProvedCritical(sys, s.Tolerance);
    const std::optional<Matrix2> inverse =
        x && Covered(sys) ? InverseProjection(sys) : std::nullopt;
    if (!inverse) {
        return std::nullopt;
    }
    const Matrix2 hessian = HeightHessian(sys, *inverse);
    if (!LevelCurvesLeave(hessian)) {
        return std::nullopt;
    }
    const Bounds onS = XyzOver(pair.A, (*x)[0], (*x)[1]);
    const Bounds gap = MinusOf(onS, XyzOver(pair.B, (*x)[2], (*x)[3]));
    std::optional<std::vector<Seed>> ends =
        NodeEndsOf(sys, NodeAt{*x, DotOf(frame->N, gap), hessian, *inverse}, s.Join);
    if (!ends) {
        return std::nullopt;
    }
    const std::array<double, 4> m = MidpointsOf<4>(*x);
    const Seed node{Middle(onS), Point2{m[0], m[1]}, Point2{m[2], m[3]}};
    return NodeEnds{node, std::move(*ends)};
}

/// A patch cut into a centre over `inner`, strictly inside it, and the four
/// pieces around it: two strips across u, and two pieces of the middle strip.
struct Centred {
    Patch Centre;
    std::vector<Patch> Rest;
};

[[nodiscard]] std::optional<Centred> CentredOn(const Patch& p, const Box2& inner)
{
    auto low = Split(p, true, inner.X.GetLower());
    if (!low) {
        return std::nullopt;
    }
    auto high = Split(low->second, true, inner.X.GetUpper());
    if (!high) {
        return std::nullopt;
    }
    auto below = Split(high->first, false, inner.Y.GetLower());
    if (!below) {
        return std::nullopt;
    }
    auto above = Split(below->second, false, inner.Y.GetUpper());
    if (!above) {
        return std::nullopt;
    }
    std::vector<Patch> rest{std::move(low->first), std::move(high->second), std::move(below->first),
                            std::move(above->second)};
    return Centred{std::move(above->first), std::move(rest)};
}

/// The box of half-widths r around a point.
/// The box reaching r beyond a point and 2 SplitFraction r, some 0.94 r,
/// short of it: off-centre for the reason the cuts are, as symmetric input
/// would otherwise end a branch exactly on a cut.
[[nodiscard]] Box2 Around(const Point2& c, const double ru, const double rv)
{
    constexpr double shorter = 2.0 * SplitFraction;
    return {Interval{c.X - shorter * ru, c.X + ru}, Interval{c.Y - shorter * rv, c.Y + rv}};
}

/// How far a point lies inside a box, along each direction.
[[nodiscard]] std::array<double, 2> Margins(const Point2& c, const Box2& box)
{
    return {std::min(c.X - box.X.GetLower(), box.X.GetUpper() - c.X),
            std::min(c.Y - box.Y.GetLower(), box.Y.GetUpper() - c.Y)};
}

[[nodiscard]] Seed Swapped(const Seed& seed)
{
    return {seed.Position, seed.OnB, seed.OnA};
}

[[nodiscard]] Vec3 NormalAt(const SecondNets& s, const Box2& box, const double u, const double v)
{
    return MiddleVec(CrossOf(XyzAt(s.Du, box, u, v), XyzAt(s.Dv, box, u, v)));
}

/// Where a pair (small, cover) may hold a node: a critical point isolated,
/// within the tolerance, its relative form indefinite. A touch or anything
/// unsettled is not one.
[[nodiscard]] std::optional<std::array<double, 4>> NodeCandidate(const PatchPair& pair,
                                                                 const SearchSettings& s)
{
    const SecondNets small = SecondNetsOf(pair.A);
    const SecondNets cover = SecondNetsOf(pair.B);
    const std::optional<NodeFrame> frame = FrameAbout(MiddleVec(NormalHull(cover).N));
    if (!frame) {
        return std::nullopt;
    }
    const Isolated<4> critical = CriticalPoint({pair, small, cover, *frame});
    if (critical.Outcome != Isolation::Unique ||
        GapAt(pair, critical.Box).GetUpper() > s.Tolerance) {
        return std::nullopt;
    }
    const std::array<double, 4> m = MidpointsOf<4>(critical.Box);
    const Vec3 normal = NormalAt(cover, pair.B.Parameters, m[2], m[3]);
    if (!(Length(normal) > 0.0)) {
        return std::nullopt;
    }
    const Vec3 n = Normalized(normal);
    if (FormSign(LocalAt(small, pair.A.Parameters, m[0], m[1], n),
                 LocalAt(cover, pair.B.Parameters, m[2], m[3], n)) >= 0) {
        return std::nullopt;
    }
    return m;
}

/// How often the covering and the small patch are narrowed around a node
/// candidate before it is left to the cutting.
constexpr int CoverNarrowings = 3;
constexpr int SmallNarrowings = 8;

/// A node settled on the centred pair: its four arcs, and the rest of the
/// pair, cut around the centred one, back on the stack.
void RecordNode(const NodeEnds& ends, const Centred& small, const Centred& cover,
                const PatchPair& oriented, const bool smallIsA, Findings& found,
                std::vector<PatchPair>& stack)
{
    const auto ordered = [&](const Seed& seed) { return smallIsA ? seed : Swapped(seed); };
    const Seed node = ordered(ends.Node);
    const Box2& boxA = smallIsA ? small.Centre.Parameters : cover.Centre.Parameters;
    const Box2& boxB = smallIsA ? cover.Centre.Parameters : small.Centre.Parameters;
    for (const Seed& e : ends.Ends) {
        const Seed end = ordered(e);
        const Vec3 d = ToVec(end.Position) - ToVec(node.Position);
        found.Arcs.push_back({node, end, d / Length(d), boxA, boxB});
    }
    found.Nodes.push_back(node.Position);
    for (const Patch& piece : small.Rest) {
        stack.push_back(smallIsA ? PatchPair{piece, oriented.B} : PatchPair{oriented.B, piece});
    }
    for (const Patch& piece : cover.Rest) {
        stack.push_back(smallIsA ? PatchPair{small.Centre, piece} : PatchPair{piece, small.Centre});
    }
}

/// A node candidate's place in an oriented pair (small, cover): its
/// parameters on each, and how far they lie inside the patches.
struct Placed {
    Point2 OnS;
    Point2 OnC;
    std::array<double, 2> MarginS;
    std::array<double, 2> MarginC;
};

/// The candidate's covering cell narrowed j times, and within it the small
/// cell narrowed until the node is proved, or nothing.
[[nodiscard]] bool SettleCentred(const PatchPair& oriented, const Placed& at, const int j,
                                 const bool smallIsA, const SearchSettings& s, Findings& found,
                                 std::vector<PatchPair>& stack)
{
    const double fc = std::ldexp(0.5, -j);
    const std::optional<Centred> cover =
        CentredOn(oriented.B, Around(at.OnC, fc * at.MarginC[0], fc * at.MarginC[1]));
    if (!cover) {
        return false;
    }
    const Vec3 normal =
        NormalAt(SecondNetsOf(cover->Centre), cover->Centre.Parameters, at.OnC.X, at.OnC.Y);
    for (int k = 0; k < SmallNarrowings; ++k) {
        const double fs = std::ldexp(0.5, -k);
        const std::optional<Centred> small =
            CentredOn(oriented.A, Around(at.OnS, fs * at.MarginS[0], fs * at.MarginS[1]));
        if (!small) {
            continue;
        }
        if (const std::optional<NodeEnds> ends =
                ProvedNode({small->Centre, cover->Centre}, normal, s)) {
            RecordNode(*ends, *small, *cover, oriented, smallIsA, found, stack);
            return true;
        }
    }
    return false;
}

/// A pair holding a node candidate, either patch taken as the small one:
/// the small patch's cell around the node narrowed until the covering
/// cell, also centred on it, is proved to cover it and facts 1 to 4 hold.
/// The pair is then the centred pair, settled, and the pieces around it.
[[nodiscard]] bool SettleNode(const PatchPair& pair, const SearchSettings& s, Findings& found,
                              std::vector<PatchPair>& stack)
{
    for (const bool smallIsA : {true, false}) {
        const PatchPair oriented = smallIsA ? pair : PatchPair{pair.B, pair.A};
        const std::optional<std::array<double, 4>> m = NodeCandidate(oriented, s);
        if (!m) {
            continue;
        }
        const Point2 onS{(*m)[0], (*m)[1]};
        const Point2 onC{(*m)[2], (*m)[3]};
        const Placed at{onS, onC, Margins(onS, oriented.A.Parameters),
                        Margins(onC, oriented.B.Parameters)};
        if (!(std::min(at.MarginS[0], at.MarginS[1]) > 0.0 &&
              std::min(at.MarginC[0], at.MarginC[1]) > 0.0)) {
            continue;
        }
        for (int j = 0; j < CoverNarrowings; ++j) {
            if (SettleCentred(oriented, at, j, smallIsA, s, found, stack)) {
                return true;
            }
        }
    }
    return false;
}

/// A pair with no advance direction, polynomial. Settled here: a proved
/// transverse node (SettleNode). An isolated touch, a unique point of common
/// tangency with the surfaces separating on every side of it, is found but
/// not yet settled. No solution rules out only tangency, not a crossing;
/// everything else is cut further, and fails at the smallest pair.
[[nodiscard]] Settled SettleTouch(const PatchPair& pair, const SearchSettings& s, Findings& found,
                                  std::vector<PatchPair>& stack)
{
    if (!s.Polynomial) {
        return Settled::Cut;
    }
    if (SettleNode(pair, s, found, stack)) {
        return Settled::Done;
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

[[nodiscard]] Settled Settle(const PatchPair& pair, const SearchSettings& s, Findings& found,
                             std::vector<PatchPair>& stack)
{
    const Hulls a = HullsOf(pair.A);
    const Hulls b = HullsOf(pair.B);
    if (const std::optional<Vec3> direction = AdvanceDirection(a, b)) {
        return SettleTransverse(pair, a, b, *direction, s, found);
    }
    return SettleTouch(pair, s, found, stack);
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
        if (Settle(pair, s, found, stack) == Settled::Done) {
            continue;
        }
        const bool smallest =
            std::max(LargestWidth(pair.A.Extent), LargestWidth(pair.B.Extent)) <= s.Tolerance;
        if (smallest || !PushSplit(stack, pair)) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSpline x BSpline: a tangential contact, a branch point not proved a "
                        "transverse node or an overlap could not be resolved at patches " +
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
    /// Whether a node of the graph is a proved branch point: chains end there.
    std::vector<bool> Branching;
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

[[nodiscard]] ArcGraph GraphOf(const Findings& found, const double join)
{
    ArcGraph g;
    const std::vector<Arc>& arcs = found.Arcs;
    for (std::size_t a = 0; a < arcs.size(); ++a) {
        const std::size_t from = NodeOf(g, arcs[a].From.Position, join);
        const std::size_t to = NodeOf(g, arcs[a].To.Position, join);
        g.Ends.push_back({from, to});
        g.AtNode[from].push_back(a);
        g.AtNode[to].push_back(a);
    }
    g.Branching.assign(g.Nodes.size(), false);
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        g.Branching[n] = std::ranges::any_of(
            found.Nodes, [&](const Point3& p) { return Distance(p, g.Nodes[n]) <= join; });
    }
    return g;
}

/// An arc used by a chain, and whether it is walked from its To end.
struct ChainStep {
    std::size_t Arc;
    bool Reversed;
};

using Chain = std::vector<ChainStep>;

/// The chain from a node along an unused arc, until it stops, closes or
/// reaches a branch point.
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
        if (g.Branching[node]) {
            return chain;
        }
    }
}

/// Every chain: those with free ends first, then those leaving branch
/// points, then the closed loops.
[[nodiscard]] std::vector<Chain> ChainsOf(const ArcGraph& g)
{
    std::vector<bool> used(g.Ends.size(), false);
    std::vector<Chain> chains;
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        if (!g.Branching[n] && g.AtNode[n].size() == 1 && !used[g.AtNode[n].front()]) {
            chains.push_back(WalkFrom(g, n, used));
        }
    }
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        while (g.Branching[n]) {
            Chain chain = WalkFrom(g, n, used);
            if (chain.empty()) {
                break;
            }
            chains.push_back(std::move(chain));
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

[[nodiscard]] bool OnBoundary(const Seed& s, const Box2& a, const Box2& b)
{
    return OnSide(s.OnA.X, a.X) || OnSide(s.OnA.Y, a.Y) || OnSide(s.OnB.X, b.X) ||
           OnSide(s.OnB.Y, b.Y);
}

/// A node must join exactly two arcs, unless it is an end of the branch on
/// a boundary of the domains, or a proved branch point joining its four;
/// anything else means the arcs do not fit.
[[nodiscard]] Result<void> CheckGraph(const ArcGraph& g, const std::vector<Arc>& arcs,
                                      const Box2& domainA, const Box2& domainB)
{
    for (std::size_t n = 0; n < g.Nodes.size(); ++n) {
        const std::size_t degree = g.AtNode[n].size();
        const Arc& arc = arcs[g.AtNode[n].front()];
        const Seed& seed = g.Ends[g.AtNode[n].front()][0] == n ? arc.From : arc.To;
        const bool fits = g.Branching[n]
                              ? degree == 4
                              : degree <= 2 && (degree == 2 || OnBoundary(seed, domainA, domainB));
        if (!fits) {
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
    const ArcGraph graph = GraphOf(found, s.Join);
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
