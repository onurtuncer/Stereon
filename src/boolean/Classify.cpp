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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Classify;
import Stereon.Robust;

namespace Stereon {

namespace {

/// How many times the seed probes' offset is halved: from a quarter of the
/// face's parameter extent down to 2^-40 of it, below any edge tolerance, so
/// a region narrower than that has no interior to seed.
constexpr int SeedHalvings = 40;

/// The k-th ray direction: a fixed turn about z by irrational fractions of
/// a turn, tilted by an irrational fraction too, so that no direction
/// aligns with the frames of boxes and cylinders.
[[nodiscard]] Vec3 Direction(const int k)
{
    const double turn = static_cast<double>(k) * std::numbers::phi * std::numbers::pi;
    const double tilt = 0.37 + 0.23 * std::numbers::sqrt2 * static_cast<double>(k);
    return Normalized(Vec3{Cos(turn), Sin(turn), Sin(tilt)});
}

// --- Regions ------------------------------------------------------------

/// The faces of the body grouped across every edge but the imprinted ones:
/// a union-find over the faces, joined by the edges they share.
[[nodiscard]] Result<std::vector<Region>> RegionsOf(const Shape& shape,
                                                    const std::set<std::uint32_t>& imprinted)
{
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(shape));
    std::map<std::uint32_t, std::uint32_t> parent;
    std::vector<FaceId> order;
    for (const FaceId face : faces) {
        parent[face.GetIndex()] = face.GetIndex();
        order.push_back(face);
    }
    const auto find = [&](std::uint32_t i) {
        while (parent[i] != i) {
            i = parent[i] = parent[parent[i]];
        }
        return i;
    };
    STEREON_TRY(const auto edges, AllEntities<EntityKind::Edge>(shape));
    for (const EdgeId edge : edges) {
        if (imprinted.contains(edge.GetIndex())) {
            continue;
        }
        STEREON_TRY(const std::vector<FaceId> adjacent, FacesOf(shape, edge));
        for (std::size_t i = 1; i < adjacent.size(); ++i) {
            parent[find(adjacent[i].GetIndex())] = find(adjacent.front().GetIndex());
        }
    }
    std::map<std::uint32_t, std::size_t> regionOf; // root -> region
    std::vector<Region> regions;
    for (const FaceId face : order) {
        const std::uint32_t root = find(face.GetIndex());
        if (!regionOf.contains(root)) {
            regionOf[root] = regions.size();
            regions.push_back({});
        }
        regions[regionOf[root]].Faces.push_back(face);
    }
    return regions;
}

// --- Seeds --------------------------------------------------------------

/// The midpoint of a coedge's pcurve and the unit normal to it there.
struct Midpoint2 {
    Point2 At;
    Vec2 Across;
};

/// Probes the box's centre, then offsets from the midpoints of the face's
/// pcurves, inward and outward, at a quarter of the face's parameter extent
/// halved SeedHalvings times: every edge at one distance before the next
/// smaller, so the widest seed is found first and a region however narrow
/// still gets one. The first probe FaceContains puts strictly inside.
[[nodiscard]] Result<std::optional<Point2>> FirstInside(const Shape& shape, const FaceId face,
                                                        const Box2& box,
                                                        const std::vector<Midpoint2>& middles,
                                                        const Context& ctx)
{
    const auto inside = [&](const Point2& probe) -> Result<bool> {
        STEREON_TRY(const std::optional<bool> in, FaceContains(shape, face, probe, ctx));
        return in && *in;
    };
    const Point2 centre{Midpoint(box.X), Midpoint(box.Y)};
    STEREON_TRY(const bool atCentre, inside(centre));
    if (atCentre) {
        return std::optional<Point2>{centre};
    }
    double step = 0.25 * std::max(Width(box.X), Width(box.Y));
    for (int halving = 0; halving < SeedHalvings; ++halving, step *= 0.5) {
        for (const Midpoint2& m : middles) {
            for (const Point2& probe : {m.At + step * m.Across, m.At - step * m.Across}) {
                STEREON_TRY(const bool in, inside(probe));
                if (in) {
                    return std::optional<Point2>{probe};
                }
            }
        }
    }
    return std::optional<Point2>{};
}

/// A point of the face's parameter domain inside its loops, off the
/// boundary (FirstInside). The box is the hull of the pcurves' ends and
/// midpoints, so a loop of one closed pcurve still spans it.
[[nodiscard]] Result<std::optional<Point2>> SeedIn(const Shape& shape, const FaceId face,
                                                   const Context& ctx)
{
    std::vector<Midpoint2> middles;
    Box2 box = Box2::Empty();
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const std::optional<Curve2d> pcurve = PCurveOf(shape, coedge);
            if (!pcurve) {
                continue;
            }
            const Interval range = RangeOf(shape, EdgeOf(shape, coedge));
            const double middle = Midpoint(range);
            const Point2 p = Eval(*pcurve, middle);
            middles.push_back({p, Perpendicular(Normalized(Deriv(*pcurve, middle, 1)))});
            box = Hull(
                Hull(Hull(box, Eval(*pcurve, range.GetLower())), Eval(*pcurve, range.GetUpper())),
                p);
        }
    }
    if (box.IsEmpty()) {
        return std::optional<Point2>{};
    }
    return FirstInside(shape, face, box, middles, ctx);
}

// --- Rays against the other body ----------------------------------------

/// A ray's meeting with a face of the other body: the ray parameter, or
/// nothing when the meeting is ambiguous (the surface grazed, the ray's
/// origin on it, or the hit on the face's boundary).
struct Hit {
    double Parameter;
};

/// A vector of outward-rounded intervals, for the ray's decisions.
using IVec3 = std::array<Interval, 3>;

[[nodiscard]] IVec3 Enclosed(const Vec3& v)
{
    return {Interval{v.X}, Interval{v.Y}, Interval{v.Z}};
}

/// p - q, enclosed.
[[nodiscard]] IVec3 Between(const Point3& p, const Point3& q)
{
    return {Interval{p.X} - Interval{q.X}, Interval{p.Y} - Interval{q.Y},
            Interval{p.Z} - Interval{q.Z}};
}

[[nodiscard]] Interval DotOf(const IVec3& a, const IVec3& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[nodiscard]] IVec3 CrossOf(const IVec3& a, const IVec3& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

/// Whether the interval is certainly positive, or certainly nonzero.
[[nodiscard]] bool Positive(const Interval& i)
{
    return SignOf(i) == Sign::Positive;
}

[[nodiscard]] bool NonZero(const Interval& i)
{
    const std::optional<Sign> sign = SignOf(i);
    return sign && *sign != Sign::Zero;
}

/// The parameters at which the line meets the plane, or nothing when the
/// line may run along it: d . (X x Y) is decided in Interval arithmetic on
/// the input doubles, and a sign the interval does not settle makes the ray
/// ambiguous, so another direction is tried; the parameter itself is a
/// construction, and FaceContains judges the hit within the tolerances.
[[nodiscard]] std::optional<std::vector<double>> MeetPlane(const Plane& plane, const Point3& origin,
                                                           const Vec3& d)
{
    const Interval along =
        DotOf(Enclosed(d), CrossOf(Enclosed(plane.XAxis), Enclosed(plane.YAxis)));
    if (!NonZero(along)) {
        return std::nullopt;
    }
    const Vec3 n = Cross(plane.XAxis, plane.YAxis);
    return std::vector<double>{Dot(plane.Origin - origin, n) / Dot(d, n)};
}

/// The parameters at which the line meets the cylinder, up to two, or
/// nothing when the line may run along the axis or graze the surface. With
/// A = X x Y, q = A . A and w = origin - O, the line o + s d meets the
/// cylinder where q |(w + s d)|^2 - ((w + s d) . A)^2 = q r^2, the quadratic
/// a s^2 + 2 b s + c = 0 of
///     a = q |d|^2 - (d . A)^2, b = q (w . d) - (w . A)(d . A),
///     c = q (|w|^2 - r^2) - (w . A)^2,
/// whose coefficients and discriminant b^2 - a c are enclosed in Interval
/// arithmetic straight from the input doubles, not from rounded values: a
/// certain negative discriminant is a miss, a certain positive one two
/// crossings, and anything the intervals do not settle is ambiguous.
[[nodiscard]] std::optional<std::vector<double>> MeetCylinder(const Cylinder& cylinder,
                                                              const Point3& origin, const Vec3& d)
{
    const IVec3 axis = CrossOf(Enclosed(cylinder.XAxis), Enclosed(cylinder.YAxis));
    const IVec3 dd = Enclosed(d);
    const IVec3 w = Between(origin, cylinder.Origin);
    const Interval q = DotOf(axis, axis);
    const Interval dA = DotOf(dd, axis);
    const Interval wA = DotOf(w, axis);
    const Interval radius{cylinder.Radius};
    const Interval a = q * DotOf(dd, dd) - Squared(dA);
    const Interval b = q * DotOf(w, dd) - wA * dA;
    const Interval c = q * (DotOf(w, w) - Squared(radius)) - Squared(wA);
    if (!Positive(a)) {
        return std::nullopt; // Along the axis, or not certainly off it.
    }
    const std::optional<Sign> discriminant = SignOf(Squared(b) - a * c);
    if (discriminant == Sign::Negative) {
        return std::vector<double>{}; // The line certainly misses the cylinder.
    }
    if (discriminant != Sign::Positive) {
        return std::nullopt; // A graze, or not certainly a crossing.
    }
    // The crossings' parameters, constructed in doubles; FaceContains then
    // judges each hit within the tolerances.
    const Vec3 unit = Normalized(Axis(cylinder));
    const Vec3 wv = origin - cylinder.Origin;
    const Vec3 dp = d - Dot(d, unit) * unit;
    const Vec3 wp = wv - Dot(wv, unit) * unit;
    const double qa = Dot(dp, dp);
    const double qb = Dot(dp, wp);
    const double qc = Dot(wp, wp) - cylinder.Radius * cylinder.Radius;
    const double root = std::sqrt(std::max(qb * qb - qa * qc, 0.0));
    return std::vector<double>{(-qb - root) / qa, (-qb + root) / qa};
}

/// The parameters at which the line meets the surface: a plane once, a
/// cylinder up to twice; nothing when a meeting is ambiguous, the line
/// running along the surface or touching it, or a sign the certified filter
/// does not settle, so that another direction is tried. A sign is never
/// taken from rounded doubles.
[[nodiscard]] std::optional<std::vector<double>> MeetSurface(const Surface& surface,
                                                             const Point3& origin, const Vec3& d)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        return MeetPlane(*plane, origin, d);
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        return MeetCylinder(*cylinder, origin, d);
    }
    return std::nullopt;
}

/// Whether the ray's meeting with the surface at parameter s crosses the
/// face: nothing when the meeting is ambiguous. The origin on the surface
/// is ambiguous only where the surface is the face; elsewhere the plane
/// merely passes through it.
[[nodiscard]] Result<std::optional<bool>> CrossesFace(const Shape& body, const FaceId face,
                                                      const Surface& surface, const Point3& origin,
                                                      const Vec3& d, const double s,
                                                      const Context& ctx)
{
    const bool atOrigin = std::abs(s) <= ctx.GetResolution();
    if (s < 0.0 && !atOrigin) { // STEREON-LINT-ALLOW(epsilon): a sign test, behind the origin
        return std::optional<bool>{false};
    }
    STEREON_TRY(const std::optional<bool> inside,
                FaceContains(body, face, Project(surface, origin + s * d), ctx));
    if (!inside || (atOrigin && *inside)) {
        return std::optional<bool>{}; // On the face's boundary, or on the face.
    }
    return std::optional<bool>{!atOrigin && *inside};
}

/// How many faces of the body the ray crosses, or nothing when a meeting is
/// ambiguous.
[[nodiscard]] Result<std::optional<int>> Crossings(const Shape& body, const Point3& origin,
                                                   const Vec3& d, const Context& ctx)
{
    int crossings = 0;
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(body));
    for (const FaceId face : faces) {
        const Surface surface = SurfaceOf(body, face);
        const std::optional<std::vector<double>> hits = MeetSurface(surface, origin, d);
        if (!hits) {
            return std::optional<int>{};
        }
        for (const double s : *hits) {
            STEREON_TRY(const std::optional<bool> crosses,
                        CrossesFace(body, face, surface, origin, d, s, ctx));
            if (!crosses) {
                return std::optional<int>{};
            }
            crossings += *crosses ? 1 : 0;
        }
    }
    return std::optional<int>{crossings};
}

/// The side by the parity of a ray's crossings of the other body, the ray
/// turned until no meeting is ambiguous; nothing when every direction of
/// the budget is.
[[nodiscard]] Result<std::optional<Side>> SeedByRays(Region& region, const Shape& other,
                                                     const int rays, const Context& ctx)
{
    for (int k = 0; k < rays; ++k) {
        ++region.Rays;
        STEREON_TRY(const std::optional<int> crossings,
                    Crossings(other, region.Seed, Direction(k), ctx));
        if (crossings) {
            return std::optional<Side>{*crossings % 2 == 1 ? Side::Inside : Side::Outside};
        }
    }
    return std::optional<Side>{};
}

/// The outward normal of a planar face, or nothing for another surface.
[[nodiscard]] std::optional<Vec3> PlaneNormal(const Shape& shape, const FaceId face)
{
    const Surface surface = SurfaceOf(shape, face); // Target points into it
    const Plane* plane = surface.Target<Plane>();
    if (plane == nullptr) {
        return std::nullopt;
    }
    const Vec3 n = Normalized(Cross(plane->XAxis, plane->YAxis));
    return IsReversed(shape, face) ? -n : n;
}

/// The largest tolerance of the face's edges, and the resolution: how far a
/// point may lie from the face's plane and still be on the face, the pair
/// tolerance of ADR-0001 as the face sees it.
[[nodiscard]] Result<double> FaceTolerance(const Shape& shape, const FaceId face,
                                           const Context& ctx)
{
    double tolerance = ctx.GetResolution();
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            tolerance = std::max(tolerance, ToleranceOf(shape, EdgeOf(shape, coedge)));
        }
    }
    return tolerance;
}

/// Whether the seed lies on the other body's planar face: within the
/// face's tolerance of its plane and strictly inside its loops.
[[nodiscard]] Result<bool> OnFace(const Shape& other, const FaceId face, const Point3& seed,
                                  const Context& ctx)
{
    const Surface surface = SurfaceOf(other, face); // Target points into it
    const Plane* plane = surface.Target<Plane>();
    STEREON_TRY(const double tolerance, FaceTolerance(other, face, ctx));
    const Vec3 n = Normalized(Cross(plane->XAxis, plane->YAxis));
    if (std::abs(Dot(seed - plane->Origin, n)) > tolerance) {
        return false;
    }
    STEREON_TRY(const std::optional<bool> inside,
                FaceContains(other, face, Project(*plane, seed), ctx));
    return inside && *inside;
}

/// OnSame or OnOpposite when the seed, on a planar face of its own body,
/// lies on a planar face of the other (a coincidence the imprint has bounded
/// with paired edges); nothing otherwise. A seed inside a region lies on a
/// face of the other body only if the two faces are coplanar: a transverse
/// face would meet the region along an imprinted edge that bounds it.
[[nodiscard]] Result<std::optional<Side>> CoincidentSide(const Shape& own, const FaceId face,
                                                         const Point3& seed, const Shape& other,
                                                         const Context& ctx)
{
    const std::optional<Vec3> normal = PlaneNormal(own, face);
    if (!normal) {
        return std::optional<Side>{};
    }
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(other));
    for (const FaceId candidate : faces) {
        const std::optional<Vec3> theirs = PlaneNormal(other, candidate);
        if (!theirs) {
            continue;
        }
        STEREON_TRY(const bool on, OnFace(other, candidate, seed, ctx));
        if (on) {
            // The normals are parallel: the sign of their product says which way.
            const bool same = Dot(*normal, *theirs) > 0.0; // STEREON-LINT-ALLOW(epsilon): a sign
            return std::optional<Side>{same ? Side::OnSame : Side::OnOpposite};
        }
    }
    return std::optional<Side>{};
}

/// The region's side by its seed: by a ray, else by the seed's bounded
/// winding number (ADR-0013, the decisions on acceptance); Unresolved when
/// neither settles it, naming the face, the rays tried and the last
/// interval (ADR-0029, section 6).
[[nodiscard]] Result<void> SeedRegion(Region& region, const Shape& own, const Shape& other,
                                      const ClassifyBudget& budget, const Context& ctx)
{
    for (const FaceId face : region.Faces) {
        STEREON_TRY(const std::optional<Point2> uv, SeedIn(own, face, ctx));
        if (!uv) {
            continue;
        }
        region.SeedFace = face;
        region.Seed = Eval(SurfaceOf(own, face), *uv);
        STEREON_TRY(const std::optional<Side> coincident,
                    CoincidentSide(own, face, region.Seed, other, ctx));
        if (coincident) {
            region.Side = *coincident;
            return {};
        }
        STEREON_TRY(const std::optional<Side> byRay, SeedByRays(region, other, budget.Rays, ctx));
        if (byRay) {
            region.Side = *byRay;
            return {};
        }
        STEREON_TRY(const WindingNumber winding,
                    WindingNumberOf(other, region.Seed, budget.Refinements, ctx));
        region.Winding = winding;
        if (winding.IsDecided()) {
            region.Side = winding.IsInside() ? Side::Inside : Side::Outside;
            return {};
        }
        return Fail(ErrorCode::Unresolved,
                    std::format("Classify: every ray from {} meets the other body ambiguously ({} "
                                "tried), and the winding number {:.6g} ± {:.3g} after {} "
                                "refinements does not settle it: the region is Unresolved "
                                "(decision D6)",
                                std::string(Describe(own, face)), region.Rays, winding.Value,
                                winding.Bound, winding.Refinements));
    }
    return Fail(ErrorCode::Unresolved,
                "Classify: no point inside any face of a region could be found for its seed");
}

// --- The cross-check ----------------------------------------------------

[[nodiscard]] std::size_t RegionOfFace(const std::vector<Region>& regions, const FaceId face)
{
    for (std::size_t i = 0; i < regions.size(); ++i) {
        if (std::ranges::find(regions[i].Faces, face) != regions[i].Faces.end()) {
            return i;
        }
    }
    return regions.size();
}

/// Every imprinted edge of the body: the regions on either side of it must
/// differ across an OppositeSide branch and agree across a SameSide one; an
/// Unresolved transition is not crossed, so it checks nothing.
[[nodiscard]] Result<void> CrossCheck(const Shape& shape, const std::vector<Region>& regions,
                                      const std::vector<std::pair<EdgeId, SideTransition>>& edges)
{
    for (const auto& [edge, transition] : edges) {
        if (transition == SideTransition::Unresolved) {
            continue;
        }
        STEREON_TRY(const std::vector<FaceId> adjacent, FacesOf(shape, edge));
        if (adjacent.size() != 2) {
            continue;
        }
        const Side first = regions[RegionOfFace(regions, adjacent[0])].Side;
        const Side second = regions[RegionOfFace(regions, adjacent[1])].Side;
        const bool differ = first != second;
        if (differ != (transition == SideTransition::OppositeSide)) {
            return Fail(ErrorCode::DegenerateResult,
                        "Classify: the regions on either side of " +
                            std::string(Describe(shape, edge)) +
                            " disagree with the branch's side transition, so the intersection, "
                            "the imprint or a seed went wrong (ADR-0013)");
        }
    }
    return {};
}

[[nodiscard]] Result<std::vector<Region>>
ClassifyBody(const Shape& own, const Shape& other,
             const std::vector<std::pair<EdgeId, SideTransition>>& edges,
             const ClassifyBudget& budget, const Context& ctx)
{
    std::set<std::uint32_t> imprinted;
    for (const auto& [edge, transition] : edges) {
        imprinted.insert(edge.GetIndex());
    }
    STEREON_TRY(std::vector<Region> regions, RegionsOf(own, imprinted));
    for (Region& region : regions) {
        STEREON_TRY_VOID(SeedRegion(region, own, other, budget, ctx));
    }
    STEREON_TRY_VOID(CrossCheck(own, regions, edges));
    return regions;
}

/// Whether the ray enters the surface's material side at a meeting: for a
/// plane, when d runs against X x Y, whose sign MeetPlane has settled; for
/// a cylinder, at the nearer of its two meetings, the one a line reaches
/// first from outside.
[[nodiscard]] bool Entering(const Surface& surface, const Vec3& d, const bool nearer)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        const double along = Dot(d, Cross(plane->XAxis, plane->YAxis));
        return along < 0.0; // STEREON-LINT-ALLOW(epsilon): a sign MeetPlane settled in Intervals
    }
    return nearer;
}

} // namespace

Result<std::optional<Point2>> Detail::SeedPointIn(const Shape& shape, const FaceId face,
                                                  const Context& ctx)
{
    return SeedIn(shape, face, ctx);
}

std::optional<std::vector<double>> Detail::RayMeets(const Surface& surface, const Point3& origin,
                                                    const Vec3& d)
{
    return MeetSurface(surface, origin, d);
}

Vec3 Detail::RayDirection(const int k)
{
    return Direction(k);
}

Result<std::optional<std::vector<Detail::FaceCrossing>>>
Detail::FaceCrossingsOf(const Shape& shape, const FaceId face, const Point3& origin, const Vec3& d,
                        const Context& ctx)
{
    const Surface surface = SurfaceOf(shape, face);
    const std::optional<std::vector<double>> hits = MeetSurface(surface, origin, d);
    if (!hits) {
        return std::optional<std::vector<FaceCrossing>>{};
    }
    std::vector<FaceCrossing> crossings;
    for (const double s : *hits) {
        STEREON_TRY(const std::optional<bool> crosses,
                    CrossesFace(shape, face, surface, origin, d, s, ctx));
        if (!crosses) {
            return std::optional<std::vector<FaceCrossing>>{};
        }
        if (*crosses) {
            crossings.push_back({s, Entering(surface, d, s <= *std::ranges::min_element(*hits))});
        }
    }
    return std::optional<std::vector<FaceCrossing>>{std::move(crossings)};
}

Result<Classification> Classify(const ImprintedBodies& bodies, const Context& ctx)
{
    return Classify(bodies, ctx, ClassifyBudget{});
}

Result<Classification> Classify(const ImprintedBodies& bodies, const Context& ctx,
                                const ClassifyBudget& budget)
{
    std::vector<std::pair<EdgeId, SideTransition>> inA;
    std::vector<std::pair<EdgeId, SideTransition>> inB;
    for (const SharedEdge& shared : bodies.Shared) {
        inA.emplace_back(shared.A, shared.Transition);
        inB.emplace_back(shared.B, shared.Transition);
    }
    Classification result;
    STEREON_TRY(result.A, ClassifyBody(bodies.A.Shape, bodies.B.Shape, inA, budget, ctx));
    STEREON_TRY(result.B, ClassifyBody(bodies.B.Shape, bodies.A.Shape, inB, budget, ctx));
    return result;
}

} // namespace Stereon
