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
#include <expected>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

module Stereon.Intersect;

import :Dispatch;

#include "Analytic.hpp"
#include "Implicit.hpp"

/// The parametric-implicit tracer (ADR-0010's analytic table, the pairs
/// without a closed form; RS-0001 Tier B): one surface, the carrier, keeps
/// its parametrisation S(u, v); the other is its implicit equation F; and
/// the intersection is the zero set of f(u, v) = F(S(u, v)) on the
/// carrier's parameter plane (Hoffmann 1989). The zero set is found by
/// interval subdivision of the carrier's parameter box with certified
/// topology, in the manner of Snyder (1992) and Plantinga and Vegter (2004):
///
/// - a box whose enclosure of f misses zero holds no point of the
///   intersection, which proves separation there;
/// - a box over which f_v (or f_u) keeps one sign holds the zero set as the
///   graph of a function of u (of v), so each of its edges across that
///   direction meets it at most once, decided by the signs at its corners;
///   an edge along it is accepted when f keeps its sign there or is
///   monotone along it, and then the crossings, in order along the graph's
///   variable, pair into the arcs inside the box;
/// - any other box is split, within a budget; a box that cannot be
///   resolved (a tangency, a singular point, a branch through a pole of
///   the parametrisation) ends the search with IntersectionFailed naming
///   it, never with an unproved result.
///
/// Since f_u or f_v is nonzero over every accepted box, the gradient of F
/// is not normal to the carrier there: every traced branch is a proved
/// transverse crossing. Arcs of adjacent boxes join where their crossings
/// on a shared edge agree, into branches that end where the search box
/// ends: at the requested domain, or at the carrier chart's seam.
///
/// A branch's 3D curve and both pcurves are then cubic B-splines on one
/// parameter: a chord-length interpolant through the traced points serves
/// as a smooth guide, each guide point is moved onto the intersection along
/// a fixed direction in the carrier's parameters, and the points so found
/// are interpolated uniformly, doubling their count until CertifyDeviation
/// proves the curve within the Context's resolution of both surfaces
/// (ADR-0009, decision D3). A branch is cut where it crosses a seam of the
/// implicit side's chart.

namespace Stereon::Detail {

namespace {

constexpr double Pi = std::numbers::pi;

/// The boxes evaluated at most, and the deepest split: a work budget, not a
/// geometric tolerance (ADR-0010).
constexpr std::size_t BoxBudget = std::size_t{1} << 18;
constexpr int MaxDepth = 56;
/// A box is split a little off its middle, so that the corners of the
/// children avoid the symmetric points of symmetric configurations.
constexpr double SplitRatio = 0.5 + 0x1p-6;
/// Points sampled along each arc for the guide, samples per fitted piece
/// at first and at most, and the scan for the implicit side's seams.
constexpr std::size_t ArcSamples = 8;
constexpr std::size_t InitialSamples = 32;
constexpr std::size_t SampleBudget = std::size_t{1} << 14;
constexpr std::size_t SeamScan = 64;
/// Bisection steps when placing a root: to the last bit of a double.
constexpr int RootSteps = 80;

/// An interval's width and middle. Not Width and Midpoint: GCC 16.2 rejects
/// their contracts here as not constant (docs/toolchain/gcc16-feature-probe.md).
[[nodiscard]] double Extent(const Interval& i)
{
    return i.GetUpper() - i.GetLower();
}
[[nodiscard]] double Middle(const Interval& i)
{
    return 0.5 * (i.GetLower() + i.GetUpper());
}

// --- Roles and the search box ---------------------------------------------

/// The pair as the tracer takes it: the carrier, parametric, and the
/// implicit side, with their requested boxes and the carrier's search box.
struct Roles {
    const Surface* CarrierErased;
    const Surface* ImplicitErased;
    Elementary Carrier;
    Elementary ImplicitGeometry;
    Quadric Implicit;
    Box2 ImplicitBox;
    Box2 Search;
    bool Exchanged;
};

/// One period of an angle, in the caller's chart where it is narrower.
[[nodiscard]] Interval OneTurn(const Interval& angle)
{
    if (angle.IsBounded() && angle.GetUpper() - angle.GetLower() < CirclePeriod) {
        return angle;
    }
    const double start = std::isfinite(angle.GetLower())   ? angle.GetLower()
                         : std::isfinite(angle.GetUpper()) ? angle.GetUpper() - CirclePeriod
                                                           : 0.0;
    return {start, start + CirclePeriod};
}

/// Angles too large to count whole periods in exactly are refused.
[[nodiscard]] Result<void> CheckChart(const Box2& box)
{
    constexpr double MaxChart = 0x1p52;
    for (const Interval& side : {box.X, box.Y}) {
        const bool bounded = side.IsBounded();
        if (bounded && std::max(std::abs(side.GetLower()), std::abs(side.GetUpper())) >= MaxChart) {
            return Fail(ErrorCode::InvalidInput,
                        "Intersect: a requested parameter box is too far out to resolve whole "
                        "periods in");
        }
    }
    return {};
}

/// r widened for the rounding of the frames and the arithmetic: a relative
/// margin and an absolute one at the scale of the configuration.
[[nodiscard]] double Generous(const double r, const double scale)
{
    return r * (1.0 + 0x1p-20) + 0x1p-20 * scale;
}

/// How far the compact implicit side reaches from its centre.
[[nodiscard]] double ReachOf(const Quadric& q)
{
    if (q.Kind == QuadricKind::Ellipsoid) {
        return std::max({q.Radius, q.SemiB, q.SemiC});
    }
    return q.Kind == QuadricKind::Torus ? q.Major + q.Radius : q.Radius;
}

/// The carrier's unbounded parameters bounded by a ball that holds the
/// implicit side: what lies outside cannot meet it.
[[nodiscard]] std::optional<Box2> BallBound(const Elementary& carrier, const Quadric& ball)
{
    const double reach = ReachOf(ball);
    if (const Plane* p = std::get_if<Plane>(&carrier)) {
        const Vec3 w = ball.Origin - p->Origin;
        const double r = Generous(reach, Length(w) + reach);
        return Box2{Around(Dot(w, p->XAxis), r), Around(Dot(w, p->YAxis), r)};
    }
    if (const Cylinder* c = std::get_if<Cylinder>(&carrier)) {
        const Vec3 w = ball.Origin - c->Origin;
        const double r = Generous(reach + 0x1p-40 * c->Radius, Length(w) + reach + c->Radius);
        return Box2{Interval::Everything(), Around(Dot(w, Axis(*c)), r)};
    }
    if (const Cone* k = std::get_if<Cone>(&carrier)) {
        // The generators are unit-speed lines from the reference circle.
        const double far = Length(ball.Origin - k->Origin) + k->Radius + reach;
        return Box2{Interval::Everything(), Around(0.0, Generous(far, far))};
    }
    return std::nullopt;
}

/// The carrier's ruling at angle u, P0(u) + v w(u), against a quadric
/// implicit side: g(v) = F(P0 + v w) is a quadratic, read off at v = 0 and
/// +-1 over a range of u, and its roots lie within Cauchy's bound
/// 1 + max(|B|, |C|) / |A| wherever its leading coefficient A keeps away
/// from zero. Nothing when A may vanish: a ruling then runs along an
/// asymptotic direction of the implicit side, and the intersection may
/// reach infinity.
[[nodiscard]] std::optional<double> RulingBound(const Elementary& carrier, const Quadric& q,
                                                const Interval& u)
{
    const auto g = [&](const double v) {
        return CompositeOf(q, JetOver(carrier, IPoint2{I{u}, I{v}})).F.Value;
    };
    const Interval g0 = g(0.0);
    const Interval gPlus = g(1.0);
    const Interval gMinus = g(-1.0);
    const Interval lead = (gPlus + gMinus) * Interval{0.5} - g0;
    const Interval middle = (gPlus - gMinus) * Interval{0.5};
    if (Contains(lead, 0.0)) {
        return std::nullopt;
    }
    const double least = std::min(std::abs(lead.GetLower()), std::abs(lead.GetUpper()));
    const auto most = [](const Interval& x) {
        return std::max(std::abs(x.GetLower()), std::abs(x.GetUpper()));
    };
    return 1.0 + std::max(most(middle), most(g0)) / least;
}

/// Cauchy's bound over the whole turn, from the bound over each of its
/// pieces; nothing unless the carrier is ruled in v and the implicit side
/// a quadric.
[[nodiscard]] std::optional<Box2> RuledBound(const Elementary& carrier, const Quadric& q)
{
    const bool ruled =
        std::holds_alternative<Cylinder>(carrier) || std::holds_alternative<Cone>(carrier);
    const bool quadric = q.Kind == QuadricKind::Cylinder || q.Kind == QuadricKind::Cone;
    if (!ruled || !quadric) {
        return std::nullopt;
    }
    constexpr int Pieces = 256;
    double bound = 0.0;
    for (int k = 0; k < Pieces; ++k) {
        const Interval piece{CirclePeriod * k / Pieces, CirclePeriod * (k + 1) / Pieces};
        const std::optional<double> b = RulingBound(carrier, q, piece);
        if (!b) {
            return std::nullopt;
        }
        bound = std::max(bound, *b);
    }
    return Box2{Interval::Everything(), Around(0.0, Generous(bound, bound))};
}

/// The carrier's natural chart: one turn of each angle, the sphere's and the
/// ellipsoid's latitudes, everything along a length.
[[nodiscard]] Box2 ChartOf(const Elementary& carrier, const Box2& box)
{
    if (std::holds_alternative<Plane>(carrier)) {
        return box;
    }
    if (std::holds_alternative<Sphere>(carrier) || std::holds_alternative<Ellipsoid>(carrier)) {
        return {OneTurn(box.X),
                Intersection(box.Y, Interval{-SpherePoleLatitude, SpherePoleLatitude})};
    }
    if (std::holds_alternative<Torus>(carrier)) {
        return {OneTurn(box.X), OneTurn(box.Y)};
    }
    return {OneTurn(box.X), box.Y};
}

/// The search box: the chart met with the bound; nothing when it stays
/// unbounded.
[[nodiscard]] std::optional<Box2> SearchBox(const Elementary& carrier, const Box2& box,
                                            const Quadric& implicit)
{
    Box2 search = ChartOf(carrier, box);
    if (IsBounded(search)) {
        return search;
    }
    std::optional<Box2> bound =
        IsCompact(implicit.Kind) ? BallBound(carrier, implicit) : RuledBound(carrier, implicit);
    if (!bound) {
        return std::nullopt;
    }
    search = {Intersection(search.X, bound->X), Intersection(search.Y, bound->Y)};
    return IsBounded(search) ? std::optional<Box2>{search} : std::nullopt;
}

/// The preference for the carrier: no pole in its chart first, then the
/// simpler implicit equation for the other side.
[[nodiscard]] int CarrierRank(const Elementary& e)
{
    // Plane, Cylinder, Cone, Sphere, Torus, Ellipsoid: the ellipsoid last,
    // its chart having poles and its implicit equation being a quadric.
    constexpr std::array<int, 6> Rank{2, 1, 3, 4, 0, 5};
    return Rank[e.index()];
}

[[nodiscard]] std::optional<Roles> RolesOf(const Surface& carrier, const Box2& carrierBox,
                                           const Surface& other, const Box2& otherBox,
                                           const bool swapped)
{
    const std::optional<Elementary> c = ElementaryOf(carrier);
    const std::optional<Elementary> o = ElementaryOf(other);
    if (!c || !o) {
        return std::nullopt;
    }
    const Quadric implicit = QuadricOf(*o);
    const std::optional<Box2> search = SearchBox(*c, carrierBox, implicit);
    if (!search) {
        return std::nullopt;
    }
    return Roles{&carrier, &other, *c, *o, implicit, otherBox, *search, swapped};
}

/// The roles in order of preference, the first with a bounded search box.
[[nodiscard]] Result<Roles> ChooseRoles(const Surface& a, const Box2& boxA, const Surface& b,
                                        const Box2& boxB)
{
    const std::optional<Elementary> ea = ElementaryOf(a);
    const std::optional<Elementary> eb = ElementaryOf(b);
    if (!ea || !eb) {
        return Fail(ErrorCode::InvalidInput, "ImplicitTrace: dispatched for other surface types");
    }
    const bool bFirst = CarrierRank(*eb) < CarrierRank(*ea);
    std::optional<Roles> roles =
        bFirst ? RolesOf(b, boxB, a, boxA, true) : RolesOf(a, boxA, b, boxB, false);
    if (!roles) {
        roles = bFirst ? RolesOf(a, boxA, b, boxB, false) : RolesOf(b, boxB, a, boxA, true);
    }
    if (!roles) {
        return Fail(ErrorCode::NotImplemented,
                    "Intersect: the intersection may reach infinity along an asymptotic direction "
                    "of both surfaces; it needs a bounded parameter box on one of them (ADR-0010)");
    }
    return *roles;
}

// --- The implicit side's box ----------------------------------------------

/// Whether the angles of (x, y) provably miss the sector `range`: both or
/// either oriented boundary half-plane excluded, as the sector is at most a
/// half turn or more.
[[nodiscard]] bool SectorMissed(const I& x, const I& y, const Interval& range)
{
    if (!range.IsBounded() || range.GetUpper() - range.GetLower() >= CirclePeriod) {
        return false;
    }
    const I start{range.GetLower()};
    const I end{range.GetUpper()};
    const bool leftOut = (Cos(start) * y - Sin(start) * x).Value.GetUpper() < 0.0;
    const bool rightOut = (Sin(end) * x - Cos(end) * y).Value.GetUpper() < 0.0;
    return range.GetUpper() - range.GetLower() > Pi ? leftOut && rightOut : leftOut || rightOut;
}

/// The sphere's latitude of local height h provably outside the range.
[[nodiscard]] bool LatitudeMissed(const I& h, const double radius, const Interval& range)
{
    const Interval latitudes =
        Intersection(range, Interval{-SpherePoleLatitude, SpherePoleLatitude});
    if (latitudes.IsEmpty()) {
        return true;
    }
    const Interval sines =
        Hull(Sin(I{latitudes.GetLower()}).Value, Sin(I{latitudes.GetUpper()}).Value);
    return !Intersects(h.Value, sines * Interval{radius});
}

/// The cone's parameters of local coordinates: v from the height, and the
/// angle of (p, q) or of (-p, -q) on the far nappe; undecided where the
/// nappe is.
[[nodiscard]] bool ConeMissed(const Quadric& q, const IVec3& l, const Box2& box)
{
    const I v = l.Z / Cos(I{q.SemiAngle});
    if (!Intersects(v.Value, box.Y)) {
        return true;
    }
    const Interval radius = (I{q.Radius} + v * Sin(I{q.SemiAngle})).Value;
    if (radius.GetLower() > 0.0) {
        return SectorMissed(l.X, l.Y, box.X);
    }
    return radius.GetUpper() < 0.0 && SectorMissed(-l.X, -l.Y, box.X);
}

/// The torus's angles of local coordinates: about its axis, and about the
/// tube's circle.
[[nodiscard]] bool TorusMissed(const Quadric& q, const IVec3& l, const Box2& box)
{
    using std::sqrt;
    const I rho = sqrt(l.X * l.X + l.Y * l.Y);
    return SectorMissed(l.X, l.Y, box.X) || SectorMissed(rho - I{q.Major}, l.Z, box.Y);
}

/// The ellipsoid's parameters of local coordinates: in s = (x / A, y / B,
/// z / C) the ellipsoid is the unit sphere, its longitude and latitude the
/// sphere's.
[[nodiscard]] bool EllipsoidMissed(const Quadric& q, const IVec3& l, const Box2& box)
{
    const I sx = l.X / I{q.Radius};
    const I sy = l.Y / I{q.SemiB};
    const I sz = l.Z / I{q.SemiC};
    return LatitudeMissed(sz, 1.0, box.Y) || SectorMissed(sx, sy, box.X);
}

/// Whether the points enclosed by `x` provably lie outside the implicit
/// side's requested box: the clip is conservative, keeping whatever it
/// cannot exclude (see Intersect).
[[nodiscard]] bool OutsideImplicitBox(const Roles& r, const IPoint3& x)
{
    const Box2& box = r.ImplicitBox;
    const IVec3 l = LocalPoint(r.Implicit, x);
    switch (r.Implicit.Kind) {
    case QuadricKind::Plane:
        return !Intersects(l.X.Value, box.X) || !Intersects(l.Y.Value, box.Y);
    case QuadricKind::Cylinder:
        return !Intersects(l.Z.Value, box.Y) || SectorMissed(l.X, l.Y, box.X);
    case QuadricKind::Cone:
        return ConeMissed(r.Implicit, l, box);
    case QuadricKind::Sphere:
        return LatitudeMissed(l.Z, r.Implicit.Radius, box.Y) || SectorMissed(l.X, l.Y, box.X);
    case QuadricKind::Torus:
        return TorusMissed(r.Implicit, l, box);
    case QuadricKind::Ellipsoid:
        return EllipsoidMissed(r.Implicit, l, box);
    }
    return false;
}

// --- Certified boxes and their arcs -----------------------------------------

/// Where an arc meets an edge of its box: the edge's line, u = Fixed when
/// Vertical and v = Fixed otherwise, the edge's extent along the line, the
/// crossing's place along it, the side of the line the box is on, and the
/// crossing's parameters.
struct Endpoint {
    bool Vertical;
    double Fixed;
    Interval Segment;
    double Along;
    int Side;
    Point2 UV;
};

/// An arc inside a certified box: its ends, the box, and whether it is the
/// graph of a function of u (else of v).
struct Arc {
    std::size_t From;
    std::size_t To;
    Interval U;
    Interval V;
    bool OverU;
};

struct Trace {
    std::vector<Endpoint> Ends;
    std::vector<Arc> Arcs;
};

/// The certain signs of f at a box's corners, [i][j] at (u_i, v_j); a zero
/// where the enclosure holds zero.
using Corners = std::array<std::array<int, 2>, 2>;

[[nodiscard]] Corners CornerSigns(const Roles& r, const Interval& u, const Interval& v)
{
    Corners s{};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const Point2 uv{i == 0 ? u.GetLower() : u.GetUpper(),
                            j == 0 ? v.GetLower() : v.GetUpper()};
            s[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                CertainSign(ValueAt(r.Carrier, r.Implicit, uv));
        }
    }
    return s;
}

/// The parameters at a fraction s of the way along an edge from `from` to
/// `to`.
[[nodiscard]] Point2 Along(const Point2& from, const Point2& to, const double s)
{
    return {from.X + s * (to.X - from.X), from.Y + s * (to.Y - from.Y)};
}

/// The root of f on the edge from `from` to `to`, whose ends have the
/// certain signs given, which differ: bisection on f's rounded sign, to the
/// last bit.
[[nodiscard]] Point2 EdgeRoot(const Roles& r, const Point2& from, const Point2& to)
{
    const int fromSign = CertainSign(ValueAt(r.Carrier, r.Implicit, from));
    double lo = 0.0;
    double hi = 1.0;
    for (int step = 0; step < RootSteps; ++step) {
        const double m = 0.5 * (lo + hi);
        if (!(lo < m && m < hi)) {
            break;
        }
        const double value = NominalAt(r.Carrier, r.Implicit, Along(from, to, m)).F;
        ((value < 0.0) == (fromSign < 0) ? lo : hi) = m;
    }
    return Along(from, to, 0.5 * (lo + hi));
}

/// What an edge along the graph's direction holds: no crossing, one, or
/// undecided. `fixedU` is the coordinate held fixed on the edge.
enum class EdgeHolds : std::uint8_t { None, One, Undecided };

[[nodiscard]] EdgeHolds AlongEdge(const Roles& r, const Interval& u, const Interval& v,
                                  const bool alongU, const int startSign, const int endSign)
{
    const BoxValues e = ValuesOver(r.Carrier, r.Implicit, u, v);
    if (!Contains(e.F, 0.0)) {
        return EdgeHolds::None;
    }
    if (!Contains(alongU ? e.Fu : e.Fv, 0.0)) {
        return startSign != endSign ? EdgeHolds::One : EdgeHolds::None;
    }
    return EdgeHolds::Undecided;
}

/// A crossing of a certified box's edge, before it joins the trace.
struct EdgeCrossing {
    double Key; ///< Its place along the graph's variable.
    Endpoint End;
};

[[nodiscard]] EdgeCrossing CrossingOn(const Roles& r, const Point2& from, const Point2& to,
                                      const bool vertical, const int side, const bool overU)
{
    const Point2 uv = EdgeRoot(r, from, to);
    const double fixed = vertical ? from.X : from.Y;
    const Interval segment =
        vertical ? Hull(Interval{from.Y}, Interval{to.Y}) : Hull(Interval{from.X}, Interval{to.X});
    const double along = vertical ? uv.Y : uv.X;
    return {overU ? uv.X : uv.Y, Endpoint{vertical, fixed, segment, along, side, uv}};
}

/// The box with its axes exchanged, so that one routine certifies a graph
/// over u and over v: Swap reads (u, v) as (v, u).
struct Frame2 {
    Interval Graph;  ///< The graph's variable.
    Interval Across; ///< The other.
    bool OverU;

    [[nodiscard]] Point2 At(const double g, const double a) const
    {
        return OverU ? Point2{g, a} : Point2{a, g};
    }
    [[nodiscard]] Interval U(const Interval& g, const Interval& a) const { return OverU ? g : a; }
    [[nodiscard]] Interval V(const Interval& g, const Interval& a) const { return OverU ? a : g; }
};

/// The crossings of a box over which the zero set is a graph over the
/// frame's graph variable: the two edges across it hold one crossing each
/// where their corners' signs differ, and the two along it one where f is
/// monotone and their corners' signs differ. Nothing when an edge along it
/// is undecided or a corner's sign is not certain.
[[nodiscard]] std::optional<std::vector<EdgeCrossing>>
GraphCrossings(const Roles& r, const Frame2& f, const Corners& c)
{
    // c[g][a] in the frame's own order.
    const auto sign = [&](const int g, const int a) {
        return f.OverU ? c[static_cast<std::size_t>(g)][static_cast<std::size_t>(a)]
                       : c[static_cast<std::size_t>(a)][static_cast<std::size_t>(g)];
    };
    std::vector<EdgeCrossing> crossings;
    const std::array<double, 2> g{f.Graph.GetLower(), f.Graph.GetUpper()};
    const std::array<double, 2> a{f.Across.GetLower(), f.Across.GetUpper()};
    for (int i = 0; i < 2; ++i) {
        // The edge across the graph at g[i]: monotone, by the box's partial.
        if (sign(i, 0) != sign(i, 1)) {
            crossings.push_back(CrossingOn(r, f.At(g[static_cast<std::size_t>(i)], a[0]),
                                           f.At(g[static_cast<std::size_t>(i)], a[1]), f.OverU,
                                           i == 0 ? 1 : -1, f.OverU));
        }
        // The edge along the graph at a[i].
        const Interval fixed{a[static_cast<std::size_t>(i)]};
        const EdgeHolds holds =
            AlongEdge(r, f.U(f.Graph, fixed), f.V(f.Graph, fixed), f.OverU, sign(0, i), sign(1, i));
        if (holds == EdgeHolds::Undecided) {
            return std::nullopt;
        }
        if (holds == EdgeHolds::One) {
            crossings.push_back(CrossingOn(r, f.At(g[0], a[static_cast<std::size_t>(i)]),
                                           f.At(g[1], a[static_cast<std::size_t>(i)]), !f.OverU,
                                           i == 0 ? 1 : -1, f.OverU));
        }
    }
    return crossings;
}

/// The arcs of a certified box, joined to the trace: its crossings in order
/// along the graph's variable pair up, since the graph exists exactly
/// between consecutive ones.
[[nodiscard]] bool AddArcs(Trace& trace, std::vector<EdgeCrossing> crossings, const Interval& u,
                           const Interval& v, const bool overU)
{
    if (crossings.size() % 2 != 0) {
        return false;
    }
    std::ranges::sort(crossings, {}, &EdgeCrossing::Key);
    for (std::size_t k = 0; k < crossings.size(); k += 2) {
        const std::size_t from = trace.Ends.size();
        trace.Ends.push_back(crossings[k].End);
        trace.Ends.push_back(crossings[k + 1].End);
        trace.Arcs.push_back({from, from + 1, u, v, overU});
    }
    return true;
}

/// Tries to certify the box: a graph over u where f_v keeps its sign, over
/// v where f_u does. True when its arcs, if any, joined the trace.
[[nodiscard]] bool CertifyBox(const Roles& r, Trace& trace, const BoxValues& values,
                              const Interval& u, const Interval& v)
{
    const Corners corners = CornerSigns(r, u, v);
    const bool certain =
        std::ranges::all_of(corners, [](const auto& row) { return row[0] != 0 && row[1] != 0; });
    if (!certain) {
        return false;
    }
    for (const bool overU : {true, false}) {
        if (Contains(overU ? values.Fv : values.Fu, 0.0)) {
            continue;
        }
        const Frame2 frame{overU ? u : v, overU ? v : u, overU};
        std::optional<std::vector<EdgeCrossing>> crossings = GraphCrossings(r, frame, corners);
        if (crossings && AddArcs(trace, std::move(*crossings), u, v, overU)) {
            return true;
        }
    }
    return false;
}

// --- The subdivision ----------------------------------------------------------

struct Cell {
    Interval U;
    Interval V;
    int Depth;
};

/// The longer side in space, by the carrier's speeds at the middle.
[[nodiscard]] bool SplitAlongU(const Roles& r, const Cell& cell)
{
    const Point2 mid{Middle(cell.U), Middle(cell.V)};
    const Jet<double> jet = JetAt(r.Carrier, mid);
    const double alongU = Extent(cell.U) * std::max(Length(jet.Su), 0x1p-30);
    const double alongV = Extent(cell.V) * std::max(Length(jet.Sv), 0x1p-30);
    return alongU >= alongV;
}

/// The two halves of a cell, a little off its middle.
[[nodiscard]] std::optional<std::array<Cell, 2>> Halves(const Roles& r, const Cell& cell)
{
    const bool alongU = SplitAlongU(r, cell);
    const Interval& side = alongU ? cell.U : cell.V;
    const double cut = side.GetLower() + SplitRatio * (side.GetUpper() - side.GetLower());
    if (!(side.GetLower() < cut && cut < side.GetUpper()) || cell.Depth >= MaxDepth) {
        return std::nullopt;
    }
    const Interval low{side.GetLower(), cut};
    const Interval high{cut, side.GetUpper()};
    if (alongU) {
        return std::array<Cell, 2>{Cell{low, cell.V, cell.Depth + 1},
                                   Cell{high, cell.V, cell.Depth + 1}};
    }
    return std::array<Cell, 2>{Cell{cell.U, low, cell.Depth + 1},
                               Cell{cell.U, high, cell.Depth + 1}};
}

[[nodiscard]] std::unexpected<KernelError> Unresolved(const Cell& cell, const std::string_view why)
{
    return Fail(ErrorCode::IntersectionFailed,
                "Intersect: the parametric-implicit trace could not resolve the parameter box [" +
                    Number(cell.U.GetLower()) + ", " + Number(cell.U.GetUpper()) + "] x [" +
                    Number(cell.V.GetLower()) + ", " + Number(cell.V.GetUpper()) + "] (" +
                    std::string(why) +
                    "): a tangency, a singular point or a branch through a pole of a "
                    "parametrisation (ADR-0010)");
}

/// Every cell of the search box resolved: excluded, outside the implicit
/// side's box, or certified with its arcs; a stack, depth first.
[[nodiscard]] Result<Trace> Subdivide(const Roles& r)
{
    Trace trace;
    std::vector<Cell> work{Cell{r.Search.X, r.Search.Y, 0}};
    for (std::size_t evaluated = 0; !work.empty(); ++evaluated) {
        const Cell cell = work.back();
        work.pop_back();
        const BoxValues values = ValuesOver(r.Carrier, r.Implicit, cell.U, cell.V);
        if (!Contains(values.F, 0.0)) {
            continue;
        }
        const IPoint3 x = JetOver(r.Carrier, IPoint2{I{cell.U}, I{cell.V}}).Point;
        if (OutsideImplicitBox(r, x) || CertifyBox(r, trace, values, cell.U, cell.V)) {
            continue;
        }
        if (evaluated >= BoxBudget) {
            return Unresolved(cell, "box budget");
        }
        const std::optional<std::array<Cell, 2>> halves = Halves(r, cell);
        if (!halves) {
            return Unresolved(cell, "smallest box");
        }
        work.push_back((*halves)[1]);
        work.push_back((*halves)[0]);
    }
    return trace;
}

// --- Joining arcs into branches ---------------------------------------------

/// Whether two crossings on one line are the same: from opposite sides,
/// each within the other's edge.
[[nodiscard]] bool SameCrossing(const Endpoint& a, const Endpoint& b)
{
    const double slack = 0x1p-40 * std::max({1.0, std::abs(a.Along), Extent(a.Segment)});
    return a.Side != b.Side && Contains(Inflated(a.Segment, slack), b.Along) &&
           Contains(Inflated(b.Segment, slack), a.Along);
}

/// The crossing each crossing joins across its edge, if any: on the same
/// line, the nearest that is the same crossing.
[[nodiscard]] std::vector<std::size_t> Links(const Trace& t)
{
    constexpr std::size_t None = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> order(t.Ends.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    const auto key = [&](const std::size_t i) {
        const Endpoint& e = t.Ends[i];
        return std::make_tuple(e.Vertical, e.Fixed, e.Along);
    };
    std::ranges::sort(order, {}, key);
    std::vector<std::size_t> link(t.Ends.size(), None);
    for (std::size_t k = 0; k + 1 < order.size(); ++k) {
        const std::size_t i = order[k];
        const Endpoint& a = t.Ends[i];
        for (std::size_t m = k + 1; m < order.size() && m <= k + 3; ++m) {
            const std::size_t j = order[m];
            const Endpoint& b = t.Ends[j];
            if (b.Vertical != a.Vertical || b.Fixed != a.Fixed) {
                break;
            }
            if (link[i] == None && link[j] == None && SameCrossing(a, b)) {
                link[i] = j;
                link[j] = i;
            }
        }
    }
    return link;
}

/// One arc of a branch, in the direction the branch runs.
struct Step {
    std::size_t Arc;
    bool Forward;
};

using Chain = std::vector<Step>;

/// The chains: from every free end first, then the closed loops that are
/// left; each arc once.
[[nodiscard]] std::vector<Chain> Chains(const Trace& t)
{
    constexpr std::size_t None = std::numeric_limits<std::size_t>::max();
    const std::vector<std::size_t> link = Links(t);
    std::vector<std::size_t> arcOf(t.Ends.size());
    for (std::size_t a = 0; a < t.Arcs.size(); ++a) {
        arcOf[t.Arcs[a].From] = a;
        arcOf[t.Arcs[a].To] = a;
    }
    std::vector<bool> used(t.Arcs.size(), false);
    const auto walk = [&](std::size_t end) {
        Chain chain;
        while (end != None && !used[arcOf[end]]) {
            const Arc& arc = t.Arcs[arcOf[end]];
            used[arcOf[end]] = true;
            const bool forward = arc.From == end;
            chain.push_back({arcOf[end], forward});
            end = link[forward ? arc.To : arc.From];
        }
        return chain;
    };
    std::vector<Chain> chains;
    for (const bool freeOnly : {true, false}) {
        for (std::size_t e = 0; e < t.Ends.size(); ++e) {
            if (!used[arcOf[e]] && (!freeOnly || link[e] == None)) {
                chains.push_back(walk(e));
            }
        }
    }
    return chains;
}

/// The point of the arc at a value of its graph's variable: the other
/// coordinate by bisection on f's rounded sign across the box.
[[nodiscard]] Point2 OnArc(const Roles& r, const Arc& arc, const double g)
{
    const Interval& across = arc.OverU ? arc.V : arc.U;
    const auto at = [&](const double a) { return arc.OverU ? Point2{g, a} : Point2{a, g}; };
    const double low = NominalAt(r.Carrier, r.Implicit, at(across.GetLower())).F;
    const double high = NominalAt(r.Carrier, r.Implicit, at(across.GetUpper())).F;
    if ((low < 0.0) == (high < 0.0)) {
        // Rounding at an end of the arc: the end nearer zero.
        return at(std::abs(low) <= std::abs(high) ? across.GetLower() : across.GetUpper());
    }
    return EdgeRoot(r, at(across.GetLower()), at(across.GetUpper()));
}

/// The guide's parameters along a chain: each arc's ends and ArcSamples - 1
/// points between, the ends shared between arcs once.
[[nodiscard]] std::vector<Point2> GuideOf(const Roles& r, const Trace& t, const Chain& chain)
{
    std::vector<Point2> guide;
    for (const Step& step : chain) {
        const Arc& arc = t.Arcs[step.Arc];
        const Point2 from = t.Ends[step.Forward ? arc.From : arc.To].UV;
        const Point2 to = t.Ends[step.Forward ? arc.To : arc.From].UV;
        const double g0 = arc.OverU ? from.X : from.Y;
        const double g1 = arc.OverU ? to.X : to.Y;
        if (guide.empty()) {
            guide.push_back(from);
        }
        for (std::size_t k = 1; k < ArcSamples; ++k) {
            const double s = static_cast<double>(k) / static_cast<double>(ArcSamples);
            guide.push_back(OnArc(r, arc, g0 + s * (g1 - g0)));
        }
        guide.push_back(to);
    }
    return guide;
}

// --- The branches -----------------------------------------------------------

/// A traced branch: the chord-length interpolant through its points, the
/// smooth guide, and the parameters of its first point on the carrier.
struct Guide {
    BSplineCurve Curve;
    Point2 Start;
    double Length;
};

[[nodiscard]] Result<std::optional<Guide>> GuideCurve(const Roles& r, const std::vector<Point2>& uv)
{
    std::vector<Point3> points;
    double length = 0.0;
    for (const Point2& p : uv) {
        const Point3 x = PointAt(r.Carrier, p);
        if (!points.empty()) {
            const double step = Distance(points.back(), x);
            if (step <= 0x1p-40 * std::max(1.0, length)) {
                continue;
            }
            length += step;
        }
        points.push_back(x);
    }
    if (points.size() < 4) {
        return std::optional<Guide>{}; // Too short to fit: reported as a point (SliverPoint).
    }
    STEREON_TRY(CurveFit fit, Interpolate(points, 3, Parametrization::ChordLength));
    return std::optional<Guide>{Guide{std::move(fit.Curve), uv.front(), length}};
}

/// The angle a, shifted by whole periods to lie nearest `near`.
[[nodiscard]] double Unwrapped(const double a, const double near)
{
    return a - CirclePeriod * std::round((a - near) / CirclePeriod);
}

/// Which of the carrier's parameters are angles: (u, v).
[[nodiscard]] std::array<bool, 2> Periodic(const Elementary& e)
{
    if (std::holds_alternative<Plane>(e)) {
        return {false, false};
    }
    return {true, std::holds_alternative<Torus>(e)};
}

/// The parameters of p on the surface, each angle nearest `near`.
[[nodiscard]] Point2 ParametersNear(const Elementary& e, const Point3& p, const Point2& near)
{
    const std::array<bool, 2> periodic = Periodic(e);
    Point2 uv = ProjectOnto(e, p);
    if (periodic[0]) {
        uv.X = Unwrapped(uv.X, near.X);
    }
    if (periodic[1]) {
        uv.Y = Unwrapped(uv.Y, near.Y);
    }
    return uv;
}

/// The point of the intersection reached from `start` along the fixed
/// direction of f's gradient there: f(start + s g) = 0 by Newton's method
/// in s. A fixed direction makes the point a smooth function of the start.
[[nodiscard]] Result<Point2> OntoCurve(const Roles& r, const Point2& start)
{
    const Composite<double> c0 = NominalAt(r.Carrier, r.Implicit, start);
    const Vec2 g{c0.Fu, c0.Fv};
    const double g2 = Dot(g, g);
    if (!(g2 > 0.0)) {
        return Fail(ErrorCode::IntersectionFailed,
                    "Intersect: the traced curve meets a point where f's gradient vanishes");
    }
    double s = -c0.F / g2;
    for (int k = 0; k < 8; ++k) {
        const Composite<double> c = NominalAt(r.Carrier, r.Implicit, start + s * g);
        const double slope = c.Fu * g.X + c.Fv * g.Y;
        if (!(std::abs(slope) > 0.0)) {
            break;
        }
        const double step = c.F / slope;
        s -= step;
        if (std::abs(step) <= 0x1p-52 * std::abs(s)) {
            break;
        }
    }
    return start + s * g;
}

/// A point of a branch: in space, and on each surface.
struct BranchPoint {
    Point3 X;
    Point2 OnCarrier;
    Point2 OnImplicit;
};

/// The branch at the guide's parameter g, its angles nearest `near`.
[[nodiscard]] Result<BranchPoint> PointOf(const Roles& r, const Guide& guide, const double g,
                                          const BranchPoint& near)
{
    const Point3 y = Eval(guide.Curve, std::clamp(g, 0.0, 1.0));
    STEREON_TRY(const Point2 uv, OntoCurve(r, ParametersNear(r.Carrier, y, near.OnCarrier)));
    const Point3 x = PointAt(r.Carrier, uv);
    return BranchPoint{x, uv, ParametersNear(r.ImplicitGeometry, x, near.OnImplicit)};
}

/// The first point of the branch, its angles placed by the trace.
[[nodiscard]] Result<BranchPoint> FirstPoint(const Roles& r, const Guide& guide)
{
    const Point3 x0 = PointAt(r.Carrier, guide.Start);
    const BranchPoint seed{x0, guide.Start, ProjectOnto(r.ImplicitGeometry, x0)};
    return PointOf(r, guide, 0.0, seed);
}

/// Whether the implicit side's chart has a seam in its second parameter.
[[nodiscard]] bool SecondPeriodic(const Roles& r)
{
    return r.Implicit.Kind == QuadricKind::Torus;
}

/// The implicit side's angles at the guide's parameter g, in the chart.
[[nodiscard]] Result<Point2> ChartAnglesAt(const Roles& r, const Guide& guide, const double g)
{
    const Point3 y = Eval(guide.Curve, std::clamp(g, 0.0, 1.0));
    STEREON_TRY(const Point2 uv, OntoCurve(r, ProjectOnto(r.Carrier, y)));
    return ProjectOnto(r.ImplicitGeometry, PointAt(r.Carrier, uv));
}

/// An angle of the implicit side's chart followed continuously along the
/// guide: at G, unwrapped from the guide's start.
struct Followed {
    double G;
    double Angle;
};

/// The angle in coordinate `which` (0 for u, 1 for v) at g, unwrapped to
/// lie nearest `near`.
[[nodiscard]] Result<double> AngleNear(const Roles& r, const Guide& guide, const double g,
                                       const int which, const double near)
{
    STEREON_TRY(const Point2 at, ChartAnglesAt(r, guide, g));
    return Unwrapped(which == 0 ? at.X : at.Y, near);
}

/// The angle over (g0, g1], appended so that consecutive samples turn by
/// at most an eighth of a turn: halved where it turns faster, near a pole
/// of the chart. A branch through a pole turns by a half turn in no
/// distance at all, and is refused.
[[nodiscard]] Result<void> FollowAngle(const Roles& r, const Guide& guide, const int which,
                                       const Followed& from, const double g1, const int depth,
                                       std::vector<Followed>& out)
{
    STEREON_TRY(const double a1, AngleNear(r, guide, g1, which, from.Angle));
    if (std::abs(a1 - from.Angle) <= Pi / 4.0) {
        out.push_back({g1, a1});
        return {};
    }
    const double gm = 0.5 * (from.G + g1);
    if (depth >= 48 || !(from.G < gm && gm < g1)) {
        return Fail(ErrorCode::IntersectionFailed,
                    "Intersect: a traced branch passes through a pole of a parametrisation, "
                    "where its pcurve has no chart (ADR-0009)");
    }
    STEREON_TRY(const double am, AngleNear(r, guide, gm, which, from.Angle));
    STEREON_TRY_VOID(FollowAngle(r, guide, which, from, gm, depth + 1, out));
    STEREON_TRY_VOID(FollowAngle(r, guide, which, Followed{gm, am}, g1, depth + 1, out));
    return {};
}

/// The whole period an unwrapped angle lies in.
[[nodiscard]] double PeriodOf(const double angle)
{
    return std::floor(angle / CirclePeriod);
}

/// The guide's parameter between two followed samples where the angle
/// crosses into the next period: bisection on the period of the angle,
/// unwrapped from the first sample, which is within an eighth of a turn.
[[nodiscard]] Result<double> SeamBetween(const Roles& r, const Guide& guide, const int which,
                                         const Followed& a, const Followed& b)
{
    double lo = a.G;
    double hi = b.G;
    for (int step = 0; step < RootSteps; ++step) {
        const double m = 0.5 * (lo + hi);
        if (!(lo < m && m < hi)) {
            break;
        }
        STEREON_TRY(const double angle, AngleNear(r, guide, m, which, a.Angle));
        (PeriodOf(angle) == PeriodOf(a.Angle) ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
}

/// Where the angle in coordinate `which` crosses the chart's seam, a whole
/// number of periods, along the guide.
[[nodiscard]] Result<void> SeamCuts(const Roles& r, const Guide& guide, const int which,
                                    std::vector<double>& cuts)
{
    STEREON_TRY(const Point2 first, ChartAnglesAt(r, guide, 0.0));
    std::vector<Followed> followed{{0.0, which == 0 ? first.X : first.Y}};
    for (std::size_t i = 1; i <= SeamScan; ++i) {
        const double g = static_cast<double>(i) / static_cast<double>(SeamScan);
        const Followed from = followed.back();
        STEREON_TRY_VOID(FollowAngle(r, guide, which, from, g, 0, followed));
    }
    for (std::size_t i = 1; i < followed.size(); ++i) {
        if (PeriodOf(followed[i].Angle) != PeriodOf(followed[i - 1].Angle)) {
            STEREON_TRY(const double cut,
                        SeamBetween(r, guide, which, followed[i - 1], followed[i]));
            cuts.push_back(cut);
        }
    }
    return {};
}

/// The cuts less those that would leave a piece within the resolution of
/// its start: a seam found at a branch's end to rounding, or two found a
/// hair apart. The sliver joins its neighbour; the ends 0 and 1 stay.
[[nodiscard]] std::vector<double>
WithoutSlivers(const Guide& guide, const std::vector<double>& cuts, const double resolution)
{
    // The piece's reach from its start, over its end and three points between.
    const auto apart = [&](const double g0, const double g1) {
        const Point3 start = Eval(guide.Curve, g0);
        double reach = 0.0;
        for (const double f : {0.25, 0.5, 0.75, 1.0}) {
            reach = std::max(reach, Distance(start, Eval(guide.Curve, g0 + f * (g1 - g0))));
        }
        return reach > resolution;
    };
    std::vector<double> kept{0.0};
    for (std::size_t i = 1; i + 1 < cuts.size(); ++i) {
        if (cuts[i] > kept.back() && apart(kept.back(), cuts[i]) && apart(cuts[i], 1.0)) {
            kept.push_back(cuts[i]);
        }
    }
    kept.push_back(1.0);
    return kept;
}

/// The guide's parameters where the implicit side's chart wraps, with 0
/// and 1: the pieces each pcurve stays in one chart over.
[[nodiscard]] Result<std::vector<double>> PieceCuts(const Roles& r, const Guide& guide,
                                                    const double resolution)
{
    std::vector<double> cuts{0.0, 1.0};
    if (r.Implicit.Kind != QuadricKind::Plane) {
        STEREON_TRY_VOID(SeamCuts(r, guide, 0, cuts));
    }
    if (SecondPeriodic(r)) {
        STEREON_TRY_VOID(SeamCuts(r, guide, 1, cuts));
    }
    std::ranges::sort(cuts);
    return WithoutSlivers(guide, cuts, resolution);
}

struct Samples {
    std::vector<Point3> Points;
    std::vector<Point2> OnCarrier;
    std::vector<Point2> OnImplicit;
};

/// count + 1 points of the branch at even guide parameters over [g0, g1].
[[nodiscard]] Result<Samples> SamplePiece(const Roles& r, const Guide& guide, const double g0,
                                          const double g1, const std::size_t count,
                                          const BranchPoint& start)
{
    Samples s;
    BranchPoint near = start;
    for (std::size_t i = 0; i <= count; ++i) {
        const double g = g0 + (g1 - g0) * static_cast<double>(i) / static_cast<double>(count);
        STEREON_TRY(near, PointOf(r, guide, g, near));
        s.Points.push_back(near.X);
        s.OnCarrier.push_back(near.OnCarrier);
        s.OnImplicit.push_back(near.OnImplicit);
    }
    return s;
}

/// The samples' angles, continuous along the piece, into the chart: the
/// piece crosses no seam inside, so the shift by whole periods that brings
/// its middle sample into [0, 2 pi] brings every sample, and the ends,
/// which may lie past the seam to rounding, are clamped onto it.
void Chart(std::vector<Point2>& uv, const std::array<bool, 2>& periodic)
{
    for (std::size_t axis = 0; axis < 2; ++axis) {
        if (!periodic[axis]) {
            continue;
        }
        const auto get = [axis](Point2& p) -> double& { return axis == 0 ? p.X : p.Y; };
        const double shift = CirclePeriod * std::floor(get(uv[uv.size() / 2]) / CirclePeriod);
        for (Point2& p : uv) {
            get(p) = std::clamp(get(p) - shift, 0.0, CirclePeriod);
        }
    }
}

/// The implicit side's periodic parameters.
[[nodiscard]] std::array<bool, 2> ImplicitPeriodic(const Roles& r)
{
    return {r.Implicit.Kind != QuadricKind::Plane, SecondPeriodic(r)};
}

/// The fitted branch over a piece, on the parameter range given.
[[nodiscard]] Result<IntersectionBranch> FitPiece(const Roles& r, Samples s, const Interval& range)
{
    Chart(s.OnImplicit, ImplicitPeriodic(r));
    STEREON_TRY(const CurveFit fit, Interpolate(s.Points, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve curve, Reparametrized(fit.Curve, range));
    STEREON_TRY(const CurveFit2d fitC, Interpolate(s.OnCarrier, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d onC, Reparametrized(fitC.Curve, range));
    STEREON_TRY(const CurveFit2d fitI, Interpolate(s.OnImplicit, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d onI, Reparametrized(fitI.Curve, range));
    return Crossing(Curve(std::move(curve)), Curve2d(std::move(onC)), Curve2d(std::move(onI)),
                    range);
}

/// The largest deviation at the midpoints between the samples, of the fit
/// from the branch: a guide to the sample count; the certificate proves.
[[nodiscard]] Result<double> MidpointDeviation(const Roles& r, const Guide& guide, const double g0,
                                               const double g1, const IntersectionBranch& branch,
                                               const std::size_t count, const BranchPoint& start)
{
    double worst = 0.0;
    BranchPoint near = start;
    const double width = Extent(branch.Range);
    for (std::size_t i = 0; i < count; ++i) {
        const double s = (static_cast<double>(i) + 0.5) / static_cast<double>(count);
        STEREON_TRY(near, PointOf(r, guide, g0 + (g1 - g0) * s, near));
        const double t = branch.Range.GetLower() + width * s;
        const Point3 on = Eval(branch.Curve3d, t);
        worst = std::max({worst, Distance(on, near.X),
                          Distance(on, Eval(*r.CarrierErased, Eval(branch.OnA, t))),
                          Distance(on, Eval(*r.ImplicitErased, Eval(branch.OnB, t)))});
    }
    return worst;
}

[[nodiscard]] Result<void> CheckSizeBox(const Samples& s)
{
    const Box3 sizeBox = Context::GetSizeBox();
    for (const Point3& point : s.Points) {
        if (!Contains(sizeBox, point)) {
            return Fail(ErrorCode::CertificationFailed,
                        "Intersect: a traced branch leaves the size box, where no deviation "
                        "certifies within the resolution (ADR-0001)");
        }
    }
    return {};
}

/// The certified branch over the guide's parameters [g0, g1]: the samples
/// doubled until the midpoints agree with the branch and the certificate
/// holds, within the sample budget.
[[nodiscard]] Result<IntersectionBranch> FitCertified(const Roles& r, const Guide& guide,
                                                      const double g0, const double g1,
                                                      const BranchPoint& start, const Context& ctx)
{
    const double target = ctx.GetResolution();
    const Interval range{guide.Length * g0, guide.Length * g1};
    for (std::size_t count = InitialSamples;; count *= 2) {
        STEREON_TRY(const Samples samples, SamplePiece(r, guide, g0, g1, count, start));
        STEREON_TRY_VOID(CheckSizeBox(samples));
        STEREON_TRY(IntersectionBranch branch, FitPiece(r, samples, range));
        STEREON_TRY(const double deviation,
                    MidpointDeviation(r, guide, g0, g1, branch, count, start));
        if (deviation <= target / 4.0) {
            const Result<void> certified =
                Certify(branch, *r.CarrierErased, *r.ImplicitErased, ctx);
            if (certified) {
                return branch;
            }
            if (certified.error().GetCode() != ErrorCode::CertificationFailed) {
                return std::unexpected(certified.error());
            }
        }
        if (count >= SampleBudget) {
            return Fail(ErrorCode::CertificationFailed,
                        "Intersect: a traced branch exceeds the sample budget before it "
                        "certifies within the resolution " +
                            Number(target) + " mm");
        }
    }
}

/// The parameters that hold an arc, narrower than its certified box: over
/// its graph's variable g, the arc runs between its ends; across it, f_a
/// keeps its sign in the box, so the arc's slope is at most |f_g / f_a|
/// there and it stays within that slope times the run from its first end.
/// The ends are crossings found by bisection on an edge, padded by 2^-40.
[[nodiscard]] std::pair<Interval, Interval> ArcBox(const Roles& r, const Trace& t, const Arc& arc)
{
    const Point2 p = t.Ends[arc.From].UV;
    const Point2 q = t.Ends[arc.To].UV;
    const auto padded = [](const Interval& i) {
        const double pad =
            0x1p-40 * std::max({1.0, std::abs(i.GetLower()), std::abs(i.GetUpper())});
        return Interval{i.GetLower() - pad, i.GetUpper() + pad};
    };
    const double g0 = arc.OverU ? p.X : p.Y;
    const double g1 = arc.OverU ? q.X : q.Y;
    const double a0 = arc.OverU ? p.Y : p.X;
    const Interval graph = padded(Hull(Interval{g0}, Interval{g1}));
    const BoxValues e = ValuesOver(r.Carrier, r.Implicit, arc.U, arc.V);
    const Interval slope = arc.OverU ? e.Fu / e.Fv : e.Fv / e.Fu;
    const double steepest = std::max(std::abs(slope.GetLower()), std::abs(slope.GetUpper()));
    const double reach =
        steepest * (Interval{graph.GetUpper()} - Interval{graph.GetLower()}).GetUpper();
    const Interval across =
        Intersection(arc.OverU ? arc.V : arc.U, padded(Interval{a0 - reach, a0 + reach}));
    const Interval onGraph = Intersection(arc.OverU ? arc.U : arc.V, graph);
    return arc.OverU ? std::pair{onGraph, across} : std::pair{across, onGraph};
}

/// A branch too short to fit a guide through, its points within 2^-40 of
/// one another: a sliver of a branch at a corner of the search box, or a
/// loop as small. Dropped, it would leave an empty result claiming a
/// separation there (ADR-0029); it is reported as a point of contact on
/// the curve at its first point, bounded by the farthest corner of the
/// enclosures in space of the boxes that hold its arcs (ArcBox).
[[nodiscard]] Result<IntersectionPoint> SliverPoint(const Roles& r, const Trace& t,
                                                    const Chain& chain, const Point2& start)
{
    STEREON_TRY(const Point2 uv, OntoCurve(r, start));
    const Point3 x = PointAt(r.Carrier, uv);
    double bound = 0.0;
    for (const Step& step : chain) {
        const auto [u, v] = ArcBox(r, t, t.Arcs[step.Arc]);
        const IPoint3 over = JetOver(r.Carrier, IPoint2{I{u}, I{v}}).Point;
        const auto farther = [](const double c, const Interval& range) {
            return std::max(std::abs(c - range.GetLower()), std::abs(range.GetUpper() - c));
        };
        bound = std::max(bound, std::hypot(farther(x.X, over.X.Value), farther(x.Y, over.Y.Value),
                                           farther(x.Z, over.Z.Value)));
    }
    return IntersectionPoint{x, uv, ProjectOnto(r.ImplicitGeometry, x), bound};
}

/// Every piece of a chain's branch, certified, into the result; a branch
/// too short to fit, as a point (SliverPoint).
[[nodiscard]] Result<void> AddBranch(SurfaceSurfaceIntersection& result, const Roles& r,
                                     const Trace& t, const Chain& chain, const Context& ctx)
{
    const std::vector<Point2> samples = GuideOf(r, t, chain);
    STEREON_TRY(const std::optional<Guide> guide, GuideCurve(r, samples));
    if (!guide) {
        STEREON_TRY(IntersectionPoint point, SliverPoint(r, t, chain, samples.front()));
        result.Points.push_back(point);
        return {};
    }
    STEREON_TRY(const std::vector<double> cuts, PieceCuts(r, *guide, ctx.GetResolution()));
    STEREON_TRY(BranchPoint start, FirstPoint(r, *guide));
    for (std::size_t k = 1; k < cuts.size(); ++k) {
        if (k > 1) {
            STEREON_TRY(start, PointOf(r, *guide, cuts[k - 1], start));
            start.OnImplicit = ProjectOnto(r.ImplicitGeometry, start.X);
        }
        STEREON_TRY(IntersectionBranch branch,
                    FitCertified(r, *guide, cuts[k - 1], cuts[k], start, ctx));
        result.Branches.push_back(std::move(branch));
    }
    return {};
}

} // namespace

Result<SurfaceSurfaceIntersection> ImplicitTrace(const Surface& a, const Box2& boxA,
                                                 const Surface& b, const Box2& boxB,
                                                 const Context& ctx)
{
    STEREON_TRY_VOID(CheckChart(boxA));
    STEREON_TRY_VOID(CheckChart(boxB));
    STEREON_TRY(const Roles roles, ChooseRoles(a, boxA, b, boxB));
    STEREON_TRY(const Trace trace, Subdivide(roles));
    SurfaceSurfaceIntersection result;
    for (const Chain& chain : Chains(trace)) {
        STEREON_TRY_VOID(AddBranch(result, roles, trace, chain, ctx));
    }
    return roles.Exchanged ? Swapped(std::move(result)) : result;
}

} // namespace Stereon::Detail
