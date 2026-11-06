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
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Imprint;

namespace Stereon {

namespace {

/// Whether the branch's 3D curve is a closed form the edge can carry as
/// its own (Exact); a fitted B-spline is derived from the surface (Trim).
[[nodiscard]] bool IsExact(const Curve& curve)
{
    return curve.Target<Line>() != nullptr || curve.Target<Circle>() != nullptr ||
           curve.Target<Ellipse>() != nullptr;
}

/// The tolerance of what the imprint makes: the resolution, or the
/// branch's certified bound when that is larger.
[[nodiscard]] double ToleranceFor(const double bound, const Context& ctx)
{
    return std::max(ctx.GetResolution(),
                    bound < std::numeric_limits<double>::infinity() ? bound : 0.0);
}

/// The builder's vertex for each vertex of each network.
using MadeVertices = std::vector<std::vector<BuilderVertexId>>;

/// A network vertex on the boundary, by its network and index.
struct BoundaryVertex {
    std::size_t Network;
    std::size_t Vertex;
};

/// Splits every edge the networks' vertices cross at, once per distinct
/// crossing in order of the parameter along the edge: each split keeps the
/// lower stretch on the edge and makes the upper, so the next crossing on
/// the same edge splits the upper. Crossings of one edge within tolerance
/// of each other, from the two faces that share it, are one vertex. A
/// vertex at a loop vertex is that vertex.
[[nodiscard]] Result<void> SplitCrossedEdges(ShapeBuilder& builder, const Shape& shape,
                                             const std::span<const FaceNetwork> networks,
                                             MadeVertices& made, std::vector<Imprint>& imprints)
{
    std::map<std::uint32_t, std::vector<BoundaryVertex>> byEdge; // frozen edge index -> vertices
    for (std::size_t n = 0; n < networks.size(); ++n) {
        for (std::size_t i = 0; i < networks[n].Vertices.size(); ++i) {
            const std::optional<LoopCrossing>& c = networks[n].Vertices[i].Crossing;
            if (!c) {
                continue;
            }
            if (c->Vertex.IsValid()) {
                made[n][i] = builder.Handle(c->Vertex);
                continue;
            }
            byEdge[EdgeOf(shape, c->Coedge).GetIndex()].push_back({n, i});
        }
    }
    const auto at = [&](const BoundaryVertex& v) -> const NetworkVertex& {
        return networks[v.Network].Vertices[v.Vertex];
    };
    for (auto& [edgeIndex, vertices] : byEdge) {
        std::ranges::sort(vertices, [&](const BoundaryVertex& p, const BoundaryVertex& q) {
            return at(p).Crossing->OnEdge < at(q).Crossing->OnEdge;
        });
        const EdgeId frozen = EdgeOf(shape, at(vertices.front()).Crossing->Coedge);
        BuilderEdgeId current = builder.Handle(frozen);
        for (std::size_t k = 0; k < vertices.size();) {
            const NetworkVertex& first = at(vertices[k]);
            const double tolerance = std::max(ToleranceOf(shape, frozen), first.Tolerance);
            STEREON_TRY(const EdgeSplit split,
                        builder.SplitEdgeMakeVertex(current, first.Crossing->OnEdge, first.Position,
                                                    tolerance));
            // Every crossing of the edge within tolerance of this one meets here.
            for (; k < vertices.size() &&
                   Distance(at(vertices[k]).Position, first.Position) <= tolerance;
                 ++k) {
                made[vertices[k].Network][vertices[k].Vertex] = split.Vertex;
                imprints[vertices[k].Network].SplitVertices.push_back(split.Vertex);
            }
            current = split.Upper;
        }
    }
    return {};
}

/// The imprint of one network in progress: the builder's vertex for each
/// network vertex (invalid until made), the frozen shape the edit started
/// from (the builder reads no geometry beyond positions), the faces the
/// original has become, the pieces done so far, and the helper struts to
/// remove at the end.
struct Work {
    ShapeBuilder& Builder;
    const FaceNetwork& Network;
    const Context& Ctx;
    Imprint& Made;
    const Shape& Frozen;
    std::vector<BuilderFaceId> Faces;
    std::vector<bool> Done;
    std::vector<EdgeVertex> Helpers;   // the struts hung, with their uses
    std::vector<BuilderFaceId> HungOn; // the face each strut was hung on
};

/// The geometry of a piece's edge: Exact with the branch's curve, or Trim
/// derived from the face's surface through the pcurve.
[[nodiscard]] Result<EdgeGeometry> GeometryOf(Work& w, const NetworkEdge& edge)
{
    EdgeGeometry geometry;
    geometry.Range = edge.Range;
    geometry.Tolerance = ToleranceFor(edge.Bound, w.Ctx);
    STEREON_TRY(geometry.PCurve, w.Builder.AddPCurve(edge.PCurve));
    if (IsExact(edge.Curve3d)) {
        geometry.Kind = EdgeCurveKind::Exact;
        STEREON_TRY(geometry.Curve, w.Builder.AddCurve(edge.Curve3d));
    } else {
        geometry.Kind = EdgeCurveKind::Trim;
    }
    return geometry;
}

/// A Trim edge gets its cached curve; an Exact one carries its own.
[[nodiscard]] Result<void> CacheIfTrim(Work& w, const BuilderEdgeId edge,
                                       const EdgeGeometry& geometry)
{
    if (geometry.Kind == EdgeCurveKind::Trim) {
        STEREON_TRY_VOID(w.Builder.CacheEdgeCurve(edge));
    }
    return {};
}

/// The coedge ending at `vertex` among the loops of the faces, with its loop
/// and face.
struct UseAt {
    BuilderCoedgeId Coedge;
    BuilderLoopId Loop;
    BuilderFaceId Face;
};

/// Where a use ends in its face's parameters: its pcurve at the end of its
/// edge's range in the loop's direction; nothing without a pcurve or a
/// bounded range.
[[nodiscard]] std::optional<Point2> UseEnd(const Work& w, const BuilderCoedgeId coedge)
{
    const std::optional<Curve2d> pcurve = w.Builder.PCurveOf(coedge);
    const Interval range = w.Builder.RangeOf(w.Builder.EdgeOf(coedge));
    if (!pcurve || range.IsEmpty() || !range.IsBounded()) {
        return std::nullopt;
    }
    const bool along = w.Builder.SenseOf(coedge) == EdgeUseSense::Same;
    return Eval(*pcurve, along ? range.GetUpper() : range.GetLower());
}

/// A use ending at a vertex, with how far its end lies from the piece's end
/// in the face's parameters.
struct Candidate {
    UseAt Use;
    EdgeUseSense Sense;
    double Gap;
};

/// The uses of the faces' loops that end at the vertex, each with how far
/// its end lies from `uv`.
[[nodiscard]] std::vector<Candidate> UsesEndingAt(const Work& w, const BuilderVertexId vertex,
                                                  const Point2& uv)
{
    std::vector<Candidate> candidates;
    for (const BuilderFaceId face : w.Faces) {
        for (const BuilderLoopId loop : w.Builder.LoopsOf(face)) {
            for (const BuilderCoedgeId coedge : w.Builder.CoedgesOf(loop)) {
                const BuilderEdgeId edge = w.Builder.EdgeOf(coedge);
                const EdgeUseSense sense = w.Builder.SenseOf(coedge);
                const BuilderVertexId end =
                    sense == EdgeUseSense::Same ? w.Builder.EndOf(edge) : w.Builder.StartOf(edge);
                if (end != vertex) {
                    continue;
                }
                const std::optional<Point2> at = UseEnd(w, coedge);
                candidates.push_back(
                    {UseAt{coedge, loop, face}, sense,
                     at ? Distance(*at, uv) : std::numeric_limits<double>::infinity()});
            }
        }
    }
    return candidates;
}

/// Among several uses ending at a vertex, the one the piece continues: the
/// one whose pcurve ends where the piece's does, at `uv`, and of those, as
/// near as the face's resolution, the one of the preferred sense, the
/// crossed use's. On a periodic face the uses on the two sides of a seam end
/// a period apart in the parameters, though at one point in space; once a
/// piece is attached at a seam vertex its own reversed use ends there too,
/// with the sense of the seam's use on the far side, so the sense alone
/// can name the wrong side (the hole of a cylinder pierced across its seam,
/// docs/imprinting.rst).
[[nodiscard]] std::optional<UseAt> UseEndingAt(const Work& w, const BuilderVertexId vertex,
                                               const std::optional<EdgeUseSense> preferred,
                                               const Point2& uv)
{
    const std::vector<Candidate> candidates = UsesEndingAt(w, vertex, uv);
    double nearest = std::numeric_limits<double>::infinity();
    for (const Candidate& c : candidates) {
        nearest = std::min(nearest, c.Gap);
    }
    std::optional<UseAt> first;
    for (const Candidate& c : candidates) {
        if (std::isfinite(nearest) && c.Gap - nearest > w.Ctx.GetResolution()) {
            continue;
        }
        if (!preferred || c.Sense == *preferred) {
            return c.Use;
        }
        if (!first) {
            first = c.Use;
        }
    }
    return first;
}

/// Keeps the sense of the use `coedge` of the face if its pcurve, at the
/// crossing's parameter on the edge, comes nearer to `uv` than any before.
void NearerUse(const Work& w, const CoedgeId coedge, const LoopCrossing& crossing, const Point2& uv,
               double& nearest, std::optional<EdgeUseSense>& sense)
{
    if (EdgeOf(w.Frozen, coedge) != EdgeOf(w.Frozen, crossing.Coedge)) {
        return;
    }
    const std::optional<Curve2d> pcurve = PCurveOf(w.Frozen, coedge);
    const double distance = pcurve ? Distance(Eval(*pcurve, crossing.OnEdge), uv) : nearest;
    if (distance < nearest) {
        nearest = distance;
        sense = StoredEdgeUseSense(w.Frozen, coedge);
    }
}

/// The sense of the face's use of the crossed edge that the piece's own
/// pcurve reaches at `uv`, its end. Where a piece ends on a seam, both uses
/// of the edge end at the vertex, one from either side of the chart, and
/// the clip names one of them without telling which: the piece's pcurve
/// does. Nothing for a junction, which no use ends at.
[[nodiscard]] Result<std::optional<EdgeUseSense>> SenseAt(const Work& w, const std::size_t vertex,
                                                          const Point2& uv)
{
    const std::optional<LoopCrossing>& crossing = w.Network.Vertices[vertex].Crossing;
    if (!crossing) {
        return std::optional<EdgeUseSense>{};
    }
    std::optional<EdgeUseSense> sense{StoredEdgeUseSense(w.Frozen, crossing->Coedge)};
    double nearest = std::numeric_limits<double>::infinity();
    STEREON_TRY(const LoopRange loops, LoopsOf(w.Frozen, w.Network.Face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(w.Frozen, loop));
        for (const CoedgeId coedge : coedges) {
            NearerUse(w, coedge, *crossing, uv, nearest, sense);
        }
    }
    return sense;
}

/// A piece between two made vertices: a face split when they lie on one
/// loop, a bridge when on two loops of one face.
[[nodiscard]] Result<void> ImprintBetween(Work& w, const std::size_t index)
{
    const NetworkEdge& edge = w.Network.Edges[index];
    const BuilderVertexId start = w.Made.Vertices[edge.Start];
    const BuilderVertexId end = w.Made.Vertices[edge.End];
    STEREON_TRY(const std::optional<EdgeUseSense> startSense,
                SenseAt(w, edge.Start, Eval(edge.PCurve, edge.Range.GetLower())));
    STEREON_TRY(const std::optional<EdgeUseSense> endSense,
                SenseAt(w, edge.End, Eval(edge.PCurve, edge.Range.GetUpper())));
    const std::optional<UseAt> a =
        UseEndingAt(w, start, startSense, Eval(edge.PCurve, edge.Range.GetLower()));
    const std::optional<UseAt> b =
        UseEndingAt(w, end, endSense, Eval(edge.PCurve, edge.Range.GetUpper()));
    if (!a || !b || a->Face != b->Face) {
        return Fail(ErrorCode::IntersectionFailed,
                    "ImprintNetwork: a piece's ends are not both on the boundary of one face");
    }
    STEREON_TRY(const EdgeGeometry geometry, GeometryOf(w, edge));
    if (a->Loop == b->Loop) {
        STEREON_TRY(const EdgeFace made, w.Builder.MakeEdgeFace(a->Coedge, b->Coedge, geometry));
        STEREON_TRY_VOID(CacheIfTrim(w, made.Edge, geometry));
        w.Faces.push_back(made.Face);
        w.Made.Pieces[index] = {made.Edge, start, end, made.Face};
        return {};
    }
    STEREON_TRY(const EdgeRing made, w.Builder.MakeEdgeKillRing(a->Coedge, b->Coedge, geometry));
    STEREON_TRY_VOID(CacheIfTrim(w, made.Edge, geometry));
    w.Made.Pieces[index] = {made.Edge, start, end, BuilderFaceId{}};
    return {};
}

/// A piece from a made vertex to a junction not yet made: a strut.
[[nodiscard]] Result<void> ImprintStrut(Work& w, const std::size_t index)
{
    const NetworkEdge& edge = w.Network.Edges[index];
    const BuilderVertexId start = w.Made.Vertices[edge.Start];
    STEREON_TRY(const std::optional<EdgeUseSense> startSense,
                SenseAt(w, edge.Start, Eval(edge.PCurve, edge.Range.GetLower())));
    const std::optional<UseAt> at =
        UseEndingAt(w, start, startSense, Eval(edge.PCurve, edge.Range.GetLower()));
    if (!at) {
        return Fail(ErrorCode::IntersectionFailed,
                    "ImprintNetwork: a piece starts at a vertex on no loop of the face");
    }
    const NetworkVertex& to = w.Network.Vertices[edge.End];
    STEREON_TRY(const EdgeGeometry geometry, GeometryOf(w, edge));
    STEREON_TRY(const EdgeVertex made,
                w.Builder.MakeEdgeVertex(at->Coedge, to.Position, to.Tolerance, geometry));
    STEREON_TRY_VOID(CacheIfTrim(w, made.Edge, geometry));
    w.Made.Vertices[edge.End] = made.Vertex;
    w.Made.Pieces[index] = {made.Edge, start, made.Vertex, BuilderFaceId{}};
    return {};
}

/// A loop's polygon in its face's parameters: each use's pcurve sampled
/// over its edge's range in the loop's direction. Nothing when a use has no
/// pcurve.
[[nodiscard]] std::optional<std::vector<Point2>> PolygonOf(const ShapeBuilder& builder,
                                                           const BuilderLoopId loop)
{
    constexpr int samples = 16; // per use
    std::vector<Point2> polygon;
    for (const BuilderCoedgeId use : builder.CoedgesOf(loop)) {
        const std::optional<Curve2d> pcurve = builder.PCurveOf(use);
        const Interval range = builder.RangeOf(builder.EdgeOf(use));
        if (!pcurve || range.IsEmpty() || !range.IsBounded()) {
            return std::nullopt;
        }
        const bool along = builder.SenseOf(use) == EdgeUseSense::Same;
        for (int k = 0; k < samples; ++k) {
            const double s = static_cast<double>(along ? k : samples - k) / samples;
            polygon.push_back(Eval(*pcurve, range.GetLower() + Width(range) * s));
        }
    }
    return polygon;
}

/// Whether uv lies inside the polygon, by the parity of its crossings.
[[nodiscard]] bool Encloses(const std::vector<Point2>& polygon, const Point2& uv)
{
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const Point2& a = polygon[i];
        const Point2& b = polygon[j];
        if ((a.Y > uv.Y) != (b.Y > uv.Y) && uv.X < a.X + (uv.Y - a.Y) * (b.X - a.X) / (b.Y - a.Y)) {
            inside = !inside;
        }
    }
    return inside;
}

/// The absolute area a polygon encloses.
[[nodiscard]] double AreaOf(const std::vector<Point2>& polygon)
{
    double twice = 0.0;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        twice += polygon[j].X * polygon[i].Y - polygon[i].X * polygon[j].Y;
    }
    return std::abs(0.5 * twice);
}

/// The piece of a face, among the parts earlier networks split it into,
/// whose outer loop encloses uv, the innermost when one part's outer loop
/// also encloses another's; the first part when none does. A network
/// imprinted on the original face's handle after another had split the
/// face was hung wherever that handle had gone: the curve around one side
/// of a cylinder as a hole of the half disc on its other side
/// (docs/imprinting.rst).
[[nodiscard]] BuilderFaceId PartHolding(const ShapeBuilder& builder,
                                        const std::vector<BuilderFaceId>& parts, const Point2& uv)
{
    BuilderFaceId holding = parts.front();
    double smallest = std::numeric_limits<double>::infinity();
    for (const BuilderFaceId part : parts) {
        const std::vector<BuilderLoopId> loops = builder.LoopsOf(part);
        if (loops.empty()) {
            continue;
        }
        const std::optional<std::vector<Point2>> outer = PolygonOf(builder, loops.front());
        if (!outer || outer->size() < 3 || !Encloses(*outer, uv)) {
            continue;
        }
        const double area = AreaOf(*outer);
        if (area < smallest) {
            smallest = area;
            holding = part;
        }
    }
    return holding;
}

/// Where a network vertex lies in the face's parameters: the end of a piece
/// it bounds.
[[nodiscard]] Point2 PointAtVertex(const FaceNetwork& network, const std::size_t vertex)
{
    for (const NetworkEdge& edge : network.Edges) {
        if (edge.Start == vertex) {
            return Eval(edge.PCurve, edge.Range.GetLower());
        }
        if (edge.End == vertex) {
            return Eval(edge.PCurve, edge.Range.GetUpper());
        }
    }
    return {};
}

/// A helper strut from a loop vertex of the face to the network vertex, so
/// that a component touching no boundary can be hung on it; removed by
/// KillEdgeMakeRing once the component has closed. The strut hangs from the
/// part of the face that holds the component (PartHolding): pieces of this
/// network that touched the boundary may already have split the face, and
/// the part that kept the original face's handle need not be the one
/// around the component.
[[nodiscard]] Result<void> HangHelper(Work& w, const std::size_t vertex)
{
    const BuilderFaceId face = PartHolding(w.Builder, w.Faces, PointAtVertex(w.Network, vertex));
    const std::vector<BuilderLoopId> loops = w.Builder.LoopsOf(face);
    if (loops.empty() || w.Builder.CoedgesOf(loops.front()).empty()) {
        return Fail(ErrorCode::InvalidInput,
                    "ImprintNetwork: the face has no loop to hang a component on");
    }
    const BuilderCoedgeId at = w.Builder.CoedgesOf(loops.front()).front();
    const BuilderEdgeId atEdge = w.Builder.EdgeOf(at);
    const BuilderVertexId from = w.Builder.SenseOf(at) == EdgeUseSense::Same
                                     ? w.Builder.EndOf(atEdge)
                                     : w.Builder.StartOf(atEdge);
    const NetworkVertex& to = w.Network.Vertices[vertex];
    const Point3 origin = w.Builder.PositionOf(from);
    EdgeGeometry strut;
    strut.Kind = EdgeCurveKind::Spatial;
    strut.Range = Interval{0.0, 1.0};
    strut.Tolerance = to.Tolerance;
    STEREON_TRY(strut.Curve, w.Builder.AddCurve(Curve(Line{origin, to.Position - origin})));
    STEREON_TRY(const EdgeVertex made,
                w.Builder.MakeEdgeVertex(at, to.Position, to.Tolerance, strut));
    w.Made.Vertices[vertex] = made.Vertex;
    w.Helpers.push_back(made);
    w.HungOn.push_back(face);
    return {};
}

/// One round over the pieces not yet imprinted: those with both ends made,
/// then one with its start made. True when a piece was imprinted.
[[nodiscard]] Result<bool> ImprintRound(Work& w)
{
    bool progress = false;
    for (std::size_t i = 0; i < w.Network.Edges.size(); ++i) {
        const NetworkEdge& edge = w.Network.Edges[i];
        if (w.Done[i] || !w.Made.Vertices[edge.Start].IsValid() ||
            !w.Made.Vertices[edge.End].IsValid()) {
            continue;
        }
        STEREON_TRY_VOID(ImprintBetween(w, i));
        w.Done[i] = true;
        progress = true;
    }
    for (std::size_t i = 0; i < w.Network.Edges.size(); ++i) {
        const NetworkEdge& edge = w.Network.Edges[i];
        if (w.Done[i] || !w.Made.Vertices[edge.Start].IsValid() ||
            w.Made.Vertices[edge.End].IsValid()) {
            continue;
        }
        STEREON_TRY_VOID(ImprintStrut(w, i));
        w.Done[i] = true;
        return true; // The new junction may complete other pieces first.
    }
    return progress;
}

/// A vertex of a piece not yet imprinted whose ends are both unmade: the
/// start of a component touching no boundary; nothing when every piece is
/// done or waits on a made end.
[[nodiscard]] std::optional<std::size_t> UnmadeStart(const Work& w)
{
    for (std::size_t i = 0; i < w.Network.Edges.size(); ++i) {
        const NetworkEdge& edge = w.Network.Edges[i];
        if (!w.Done[i] && !w.Made.Vertices[edge.Start].IsValid() &&
            !w.Made.Vertices[edge.End].IsValid()) {
            return edge.Start;
        }
    }
    return std::nullopt;
}

[[nodiscard]] bool AllDone(const Work& w)
{
    return std::ranges::all_of(w.Done, [](const bool done) { return done; });
}

/// The signed area the ring's uses bound in the face's parameter space,
/// their pcurves sampled in the order and direction the loop runs:
/// positive when it turns counter-clockwise. Nothing when a use is not a
/// piece of this network, which a ring KEMR leaves never is.
[[nodiscard]] std::optional<double> TurnOf(const Work& w, const std::vector<BuilderCoedgeId>& uses)
{
    constexpr int samples = 64; // per piece
    std::vector<Point2> polygon;
    for (const BuilderCoedgeId use : uses) {
        const BuilderEdgeId edge = w.Builder.EdgeOf(use);
        const auto piece = std::ranges::find_if(
            w.Made.Pieces, [edge](const ImprintedPiece& p) { return p.Edge == edge; });
        if (piece == w.Made.Pieces.end()) {
            return std::nullopt;
        }
        const NetworkEdge& n =
            w.Network.Edges[static_cast<std::size_t>(piece - w.Made.Pieces.begin())];
        const bool along = w.Builder.SenseOf(use) == EdgeUseSense::Same;
        for (int k = 0; k < samples; ++k) {
            const double s = static_cast<double>(along ? k : samples - k) / samples;
            polygon.push_back(Eval(n.PCurve, n.Range.GetLower() + Width(n.Range) * s));
        }
    }
    double twice = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point2& p = polygon[i];
        const Point2& q = polygon[(i + 1) % polygon.size()];
        twice += p.X * q.Y - q.X * p.Y;
    }
    return 0.5 * twice;
}

/// A component's ring: the loop of its pieces, and its turn (TurnOf).
struct Ring {
    BuilderLoopId Loop;
    double Turn = 0.0;
};

/// The component's ring among the two loops KEMR leaves when the helper
/// goes, the one made of this network's pieces. KEMR keeps the helper's
/// anchor side in the original loop and returns the other, and which side
/// the component lay on depends on the way the strut was hung, so either
/// loop may be the ring; before this was settled the face's own boundary
/// was taken for the ring, and a ring on a reversed face kept the turn of
/// an outer loop.
[[nodiscard]] std::optional<Ring> RingAmong(const Work& w, const BuilderLoopId made,
                                            const BuilderLoopId kept)
{
    for (const BuilderLoopId loop : {made, kept}) {
        if (const std::optional<double> turn = TurnOf(w, w.Builder.CoedgesOf(loop))) {
            return Ring{loop, *turn};
        }
    }
    return std::nullopt;
}

/// The network's face whose loops include the loop.
[[nodiscard]] BuilderFaceId FaceHolding(const Work& w, const BuilderLoopId loop)
{
    for (const BuilderFaceId face : w.Faces) {
        if (std::ranges::contains(w.Builder.LoopsOf(face), loop)) {
            return face;
        }
    }
    return {};
}

/// A ring a closed component leaves in the original face must be a hole: it
/// runs clockwise about Su x Sv in a face oriented with its surface, the
/// other way in a reversed face, and the disc cut off runs against it. The
/// MakeEdgeFace that closed the component, from a coedge to itself for a
/// single closed piece or between two struts for several pieces, has no
/// loop order to take that from and may leave the disc's loop in the
/// original face; when the ring turns the wrong way the two loops change
/// faces, the disc's into the original face (KFMRH) and the ring's out into
/// a new face (MFKRH), which then is the face of the piece that closed the
/// component. `made` and `kept` are the loops KEMR left (RingAmong).
[[nodiscard]] Result<void> OrientRing(Work& w, const BuilderLoopId made, const BuilderLoopId kept,
                                      const BuilderFaceId hung)
{
    const std::optional<Ring> ring = RingAmong(w, made, kept);
    if (!ring) {
        return {};
    }
    const std::vector<BuilderCoedgeId> uses = w.Builder.CoedgesOf(ring->Loop);
    const auto closer = std::ranges::find_if(w.Made.Pieces, [&](const ImprintedPiece& p) {
        return p.Face.IsValid() && std::ranges::any_of(uses, [&](const BuilderCoedgeId use) {
                   return w.Builder.EdgeOf(use) == p.Edge;
               });
    });
    if (closer == w.Made.Pieces.end()) {
        return {};
    }
    const bool reversed = IsReversed(w.Frozen, w.Network.Face);
    if (reversed ? ring->Turn > 0.0 : ring->Turn < 0.0) { // STEREON-LINT-ALLOW(epsilon): a sense
        return {};
    }
    // The ring's face, the original one or the closer's, keeps the outer
    // boundary and takes the disc's loop as its hole; the disc's face goes,
    // and the ring's loop makes the new face of the piece that closed it.
    const BuilderFaceId ringFace = FaceHolding(w, ring->Loop);
    const BuilderFaceId discFace = ringFace == closer->Face ? hung : closer->Face;
    STEREON_TRY_VOID(w.Builder.KillFaceMakeRingHole(discFace, ringFace));
    STEREON_TRY(closer->Face, w.Builder.MakeFaceKillRingHole(ring->Loop));
    std::erase(w.Faces, discFace);
    std::erase(w.Faces, ringFace);
    w.Faces.insert(w.Faces.begin(), ringFace);
    w.Faces.push_back(closer->Face);
    return {};
}

/// One network's pieces, its boundary vertices already made, imprinted on
/// `part`, the piece of the original face that holds it (PartHolding). The
/// faces the network leaves, `part` among them, are appended to `parts`.
[[nodiscard]] Result<void> ImprintPieces(ShapeBuilder& builder, const Shape& shape,
                                         const FaceNetwork& network, Imprint& made,
                                         const Context& ctx, const BuilderFaceId part,
                                         std::vector<BuilderFaceId>& parts)
{
    if (!part.IsValid()) {
        return Fail(ErrorCode::InvalidInput, "ImprintNetwork: the face is not live in the builder");
    }
    Work w{builder, network, ctx, made, shape, {part}, {}, {}, {}};
    w.Done.assign(network.Edges.size(), false);
    while (!AllDone(w)) {
        STEREON_TRY(const bool progress, ImprintRound(w));
        if (progress) {
            continue;
        }
        const std::optional<std::size_t> start = UnmadeStart(w);
        if (!start) {
            return Fail(ErrorCode::IntersectionFailed,
                        "ImprintNetwork: a piece runs into the face from a junction no other piece "
                        "reaches; the network does not close");
        }
        STEREON_TRY_VOID(HangHelper(w, *start));
    }
    for (std::size_t h = 0; h < w.Helpers.size(); ++h) {
        const BuilderLoopId kept = builder.LoopOf(w.Helpers[h].Out);
        STEREON_TRY(const BuilderLoopId freed, builder.KillEdgeMakeRing(w.Helpers[h].Edge));
        STEREON_TRY_VOID(OrientRing(w, freed, kept, w.HungOn[h]));
    }
    for (const BuilderFaceId face : w.Faces) {
        if (!std::ranges::contains(parts, face)) {
            parts.push_back(face);
        }
    }
    return {};
}

/// A point of the network inside the face, off its boundary: the middle of
/// its first piece.
[[nodiscard]] Point2 InsidePointOf(const FaceNetwork& network)
{
    const NetworkEdge& first = network.Edges.front();
    return Eval(first.PCurve, Midpoint(first.Range));
}

} // namespace

FaceNetwork NetworkOf(const Shape& shape, const FaceId face, const IntersectionBranch& branch,
                      const BranchSide side, const FaceClip& clip, const Context& ctx)
{
    FaceNetwork network{face, {}, {}};
    const Curve2d& pcurve = side == BranchSide::A ? branch.OnA : branch.OnB;
    const double tolerance = ToleranceFor(branch.Bound, ctx);
    std::map<std::size_t, std::size_t> atCrossing; // crossing index -> vertex index
    const auto vertexAt = [&](const std::optional<std::size_t>& crossing, const double t) {
        if (crossing) {
            if (const auto found = atCrossing.find(*crossing); found != atCrossing.end()) {
                return found->second;
            }
            const LoopCrossing& c = clip.Crossings[*crossing];
            network.Vertices.push_back({Eval(branch.Curve3d, c.OnBranch), tolerance, c});
            atCrossing[*crossing] = network.Vertices.size() - 1;
            return network.Vertices.size() - 1;
        }
        network.Vertices.push_back({Eval(branch.Curve3d, t), tolerance, std::nullopt});
        return network.Vertices.size() - 1;
    };
    const auto senseAt = [&](const std::optional<std::size_t>& crossing) {
        return crossing ? std::optional{StoredEdgeUseSense(shape, clip.Crossings[*crossing].Coedge)}
                        : std::nullopt;
    };
    for (const FacePiece& piece : clip.Inside) {
        const std::size_t start = vertexAt(piece.Start, piece.Range.GetLower());
        // A loop only when the branch is a closed curve: a whole branch with
        // no crossing may be a segment lying inside the face.
        const bool closed = !piece.Start && !piece.End && piece.Range == branch.Range &&
                            Distance(Eval(branch.Curve3d, branch.Range.GetLower()),
                                     Eval(branch.Curve3d, branch.Range.GetUpper())) <= tolerance;
        const std::size_t end = closed ? start : vertexAt(piece.End, piece.Range.GetUpper());
        network.Edges.push_back({branch.Curve3d, pcurve, piece.Range, branch.Bound, start, end,
                                 senseAt(piece.Start), senseAt(piece.End)});
    }
    return network;
}

Result<std::vector<Imprint>> ImprintNetworks(ShapeBuilder& builder, const Shape& shape,
                                             const std::span<const FaceNetwork> networks,
                                             const Context& ctx)
{
    std::vector<Imprint> imprints(networks.size());
    MadeVertices made(networks.size());
    for (std::size_t n = 0; n < networks.size(); ++n) {
        made[n].assign(networks[n].Vertices.size(), BuilderVertexId{});
        imprints[n].Pieces.assign(networks[n].Edges.size(), ImprintedPiece{});
    }
    STEREON_TRY_VOID(SplitCrossedEdges(builder, shape, networks, made, imprints));
    // The pieces each face has been split into by the networks so far.
    std::map<std::uint32_t, std::vector<BuilderFaceId>> parts;
    for (std::size_t n = 0; n < networks.size(); ++n) {
        imprints[n].Vertices = std::move(made[n]);
        std::vector<BuilderFaceId>& faceParts = parts[networks[n].Face.GetIndex()];
        if (faceParts.empty()) {
            faceParts.push_back(builder.Handle(networks[n].Face));
        }
        const BuilderFaceId part =
            networks[n].Edges.empty() ? faceParts.front()
                                      : PartHolding(builder, faceParts, InsidePointOf(networks[n]));
        STEREON_TRY_VOID(
            ImprintPieces(builder, shape, networks[n], imprints[n], ctx, part, faceParts));
    }
    return imprints;
}

Result<Imprint> ImprintNetwork(ShapeBuilder& builder, const Shape& shape,
                               const FaceNetwork& network, const Context& ctx)
{
    STEREON_TRY(std::vector<Imprint> imprints,
                ImprintNetworks(builder, shape, std::span<const FaceNetwork>{&network, 1}, ctx));
    return std::move(imprints.front());
}

Result<Imprint> ImprintClip(ShapeBuilder& builder, const Shape& shape, const FaceId face,
                            const IntersectionBranch& branch, const BranchSide side,
                            const FaceClip& clip, const Context& ctx)
{
    return ImprintNetwork(builder, shape, NetworkOf(shape, face, branch, side, clip, ctx), ctx);
}

} // namespace Stereon
