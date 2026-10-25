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

namespace Stereon {

namespace {

/// How many points a face is probed at for a seed inside it.
constexpr int SeedProbes = 12;

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

/// A point of the face's parameter domain inside its loops, off the
/// boundary: probed at the box's centre and along deterministic offsets
/// from the midpoints of its edges' pcurves, inward and outward.
[[nodiscard]] Result<std::optional<Point2>> SeedIn(const Shape& shape, const FaceId face,
                                                   const Context& ctx)
{
    std::vector<Point2> probes;
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
            const double middle = 0.5 * (range.GetLower() + range.GetUpper());
            const Point2 p = Eval(*pcurve, middle);
            const Vec2 across = Perpendicular(Normalized(Deriv(*pcurve, middle, 1)));
            box = Hull(Hull(box, Eval(*pcurve, range.GetLower())), Eval(*pcurve, range.GetUpper()));
            for (int k = 1; k <= SeedProbes / 2; ++k) {
                const double step = 0.05 * static_cast<double>(k);
                probes.push_back(p + step * across);
                probes.push_back(p - step * across);
            }
        }
    }
    probes.insert(probes.begin(), Point2{0.5 * (box.X.GetLower() + box.X.GetUpper()),
                                         0.5 * (box.Y.GetLower() + box.Y.GetUpper())});
    for (const Point2& probe : probes) {
        STEREON_TRY(const std::optional<bool> inside, FaceContains(shape, face, probe, ctx));
        if (inside && *inside) {
            return std::optional<Point2>{probe};
        }
    }
    return std::optional<Point2>{};
}

// --- Rays against the other body ----------------------------------------

/// A ray's meeting with a face of the other body: the ray parameter, or
/// nothing when the meeting is ambiguous (the surface grazed, the ray's
/// origin on it, or the hit on the face's boundary).
struct Hit {
    double Parameter;
};

/// The parameters at which the line meets the surface: a plane once, a
/// cylinder up to twice; nothing when the line runs along the surface or
/// touches it (a vanishing discriminant, decided in doubles: a graze is
/// ambiguous and another direction is tried).
[[nodiscard]] std::optional<std::vector<double>> MeetSurface(const Surface& surface,
                                                             const Point3& origin, const Vec3& d)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        const Vec3 n = Cross(plane->XAxis, plane->YAxis);
        const double along = Dot(d, n);
        if (along == 0.0) { // STEREON-LINT-ALLOW(epsilon): a sign test
            return std::nullopt;
        }
        return std::vector<double>{Dot(plane->Origin - origin, n) / along};
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        const Vec3 axis = Normalized(Axis(*cylinder));
        const Vec3 w = origin - cylinder->Origin;
        const Vec3 dp = d - Dot(d, axis) * axis;
        const Vec3 wp = w - Dot(w, axis) * axis;
        const double a = Dot(dp, dp);
        const double b = Dot(dp, wp);
        const double c = Dot(wp, wp) - cylinder->Radius * cylinder->Radius;
        const double discriminant = b * b - a * c;
        if (discriminant < 0.0) {         // STEREON-LINT-ALLOW(epsilon): a sign test
            return std::vector<double>{}; // The line misses the cylinder.
        }
        if (a == 0.0 || discriminant == 0.0) { // STEREON-LINT-ALLOW(epsilon): a sign test
            return std::nullopt;               // Along the axis, or a graze.
        }
        const double root = std::sqrt(discriminant);
        return std::vector<double>{(-b - root) / a, (-b + root) / a};
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

} // namespace

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
