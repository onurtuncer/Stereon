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
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Imprint;

namespace Stereon {

namespace {

/// The builder's vertex at each crossing a piece ends at: the loop vertex
/// the crossing is at, or the vertex its split made.
using CrossingVertices = std::map<std::size_t, BuilderVertexId>;

/// Whether the branch's 3D curve is a closed form the edge can carry as
/// its own (Exact); a fitted B-spline is derived from the surface (Trim).
[[nodiscard]] bool IsExact(const Curve& curve)
{
    return curve.Target<Line>() != nullptr || curve.Target<Circle>() != nullptr ||
           curve.Target<Ellipse>() != nullptr;
}

/// The tolerance of what the imprint makes: the resolution, or the
/// branch's certified bound when that is larger.
[[nodiscard]] double ToleranceFor(const IntersectionBranch& branch, const Context& ctx)
{
    return std::max(ctx.GetResolution(), branch.IsCertified() ? branch.Bound : 0.0);
}

/// Splits every edge the pieces' ends cross at, in order of the parameter
/// along the edge: each split keeps the lower stretch on the edge and makes
/// the upper one, so the next crossing on the same edge splits the upper.
[[nodiscard]] Result<void> SplitCrossedEdges(ShapeBuilder& builder, const Shape& shape,
                                             const IntersectionBranch& branch, const FaceClip& clip,
                                             CrossingVertices& vertices, Imprint& imprint,
                                             const Context& ctx)
{
    std::map<std::uint32_t, std::vector<std::size_t>> byEdge; // frozen edge index -> crossings
    std::vector<bool> wanted(clip.Crossings.size(), false);
    for (const FacePiece& piece : clip.Inside) {
        for (const std::optional<std::size_t>& end : {piece.Start, piece.End}) {
            if (end) {
                wanted[*end] = true;
            }
        }
    }
    for (std::size_t i = 0; i < clip.Crossings.size(); ++i) {
        const LoopCrossing& c = clip.Crossings[i];
        if (!wanted[i]) {
            continue;
        }
        if (c.Vertex.IsValid()) {
            vertices[i] = builder.Handle(c.Vertex);
            continue;
        }
        byEdge[EdgeOf(shape, c.Coedge).GetIndex()].push_back(i);
    }
    for (auto& [edgeIndex, crossings] : byEdge) {
        std::ranges::sort(crossings, [&](const std::size_t i, const std::size_t j) {
            return clip.Crossings[i].OnEdge < clip.Crossings[j].OnEdge;
        });
        const EdgeId frozen = EdgeOf(shape, clip.Crossings[crossings.front()].Coedge);
        BuilderEdgeId current = builder.Handle(frozen);
        const double tolerance = std::max(ToleranceOf(shape, frozen), ToleranceFor(branch, ctx));
        for (const std::size_t i : crossings) {
            const LoopCrossing& c = clip.Crossings[i];
            STEREON_TRY(const EdgeSplit split,
                        builder.SplitEdgeMakeVertex(current, c.OnEdge,
                                                    Eval(branch.Curve3d, c.OnBranch), tolerance));
            vertices[i] = split.Vertex;
            imprint.SplitVertices.push_back(split.Vertex);
            current = split.Upper;
        }
    }
    return {};
}

/// The geometry of a piece's edge: Exact with the branch's curve, or Trim
/// derived from the face's surface through the pcurve.
[[nodiscard]] Result<EdgeGeometry> GeometryOf(ShapeBuilder& builder,
                                              const IntersectionBranch& branch,
                                              const Curve2d& pcurve, const Interval& range,
                                              const Context& ctx)
{
    EdgeGeometry geometry;
    geometry.Range = range;
    geometry.Tolerance = ToleranceFor(branch, ctx);
    STEREON_TRY(geometry.PCurve, builder.AddPCurve(pcurve));
    if (IsExact(branch.Curve3d)) {
        geometry.Kind = EdgeCurveKind::Exact;
        STEREON_TRY(geometry.Curve, builder.AddCurve(branch.Curve3d));
    } else {
        geometry.Kind = EdgeCurveKind::Trim;
    }
    return geometry;
}

/// A Trim edge gets its cached curve; an Exact one carries its own.
[[nodiscard]] Result<void> CacheIfTrim(ShapeBuilder& builder, const BuilderEdgeId edge,
                                       const EdgeGeometry& geometry)
{
    if (geometry.Kind == EdgeCurveKind::Trim) {
        STEREON_TRY_VOID(builder.CacheEdgeCurve(edge));
    }
    return {};
}

/// The faces the original face has become: itself and every face a piece
/// cut off it.
struct Faces {
    std::vector<BuilderFaceId> All;
};

/// The coedge ending at `vertex` among the loops of the faces, with its loop.
struct UseAt {
    BuilderCoedgeId Coedge;
    BuilderLoopId Loop;
    BuilderFaceId Face;
};

[[nodiscard]] std::optional<UseAt> UseEndingAt(const ShapeBuilder& builder, const Faces& faces,
                                               const BuilderVertexId vertex)
{
    for (const BuilderFaceId face : faces.All) {
        for (const BuilderLoopId loop : builder.LoopsOf(face)) {
            for (const BuilderCoedgeId coedge : builder.CoedgesOf(loop)) {
                const BuilderEdgeId edge = builder.EdgeOf(coedge);
                const BuilderVertexId end = builder.SenseOf(coedge) == EdgeUseSense::Same
                                                ? builder.EndOf(edge)
                                                : builder.StartOf(edge);
                if (end == vertex) {
                    return UseAt{coedge, loop, face};
                }
            }
        }
    }
    return std::nullopt;
}

/// A piece between two vertices on the boundary: a face split when they
/// lie on one loop, a bridge when on two loops of one face.
[[nodiscard]] Result<ImprintedPiece> ImprintOpen(ShapeBuilder& builder, Faces& faces,
                                                 const BuilderVertexId start,
                                                 const BuilderVertexId end,
                                                 const EdgeGeometry& geometry)
{
    const std::optional<UseAt> a = UseEndingAt(builder, faces, start);
    const std::optional<UseAt> b = UseEndingAt(builder, faces, end);
    if (!a || !b || a->Face != b->Face) {
        return Fail(ErrorCode::IntersectionFailed,
                    "ImprintClip: a piece's ends are not both on the boundary of one face");
    }
    if (a->Loop == b->Loop) {
        STEREON_TRY(const EdgeFace made, builder.MakeEdgeFace(a->Coedge, b->Coedge, geometry));
        STEREON_TRY_VOID(CacheIfTrim(builder, made.Edge, geometry));
        faces.All.push_back(made.Face);
        return ImprintedPiece{made.Edge, start, end, made.Face};
    }
    STEREON_TRY(const EdgeRing made, builder.MakeEdgeKillRing(a->Coedge, b->Coedge, geometry));
    STEREON_TRY_VOID(CacheIfTrim(builder, made.Edge, geometry));
    return ImprintedPiece{made.Edge, start, end, BuilderFaceId{}};
}

/// A closed piece inside a face: a ring. A strut from a loop vertex to the
/// curve's start carries the closed edge in by MEF on the strut's own use,
/// and KEMR removes the strut, leaving the ring as a loop of the face and
/// the closed edge's other use as the new face's boundary.
[[nodiscard]] Result<ImprintedPiece> ImprintRing(ShapeBuilder& builder, Faces& faces,
                                                 const IntersectionBranch& branch,
                                                 const Interval& range,
                                                 const EdgeGeometry& geometry)
{
    const BuilderFaceId face = faces.All.front();
    const std::vector<BuilderLoopId> loops = builder.LoopsOf(face);
    if (loops.empty() || builder.CoedgesOf(loops.front()).empty()) {
        return Fail(ErrorCode::InvalidInput, "ImprintClip: the face has no loop to hang a ring on");
    }
    const BuilderCoedgeId at = builder.CoedgesOf(loops.front()).front();
    const BuilderEdgeId atEdge = builder.EdgeOf(at);
    const BuilderVertexId from =
        builder.SenseOf(at) == EdgeUseSense::Same ? builder.EndOf(atEdge) : builder.StartOf(atEdge);
    const Point3 origin = builder.PositionOf(from);
    const Point3 onCurve = Eval(branch.Curve3d, range.GetLower());
    EdgeGeometry strut;
    strut.Kind = EdgeCurveKind::Spatial;
    strut.Range = Interval{0.0, 1.0};
    strut.Tolerance = geometry.Tolerance;
    STEREON_TRY(strut.Curve, builder.AddCurve(Curve(Line{origin, onCurve - origin})));
    STEREON_TRY(const EdgeVertex made,
                builder.MakeEdgeVertex(at, onCurve, geometry.Tolerance, strut));
    STEREON_TRY(const EdgeFace closed, builder.MakeEdgeFace(made.Out, made.Out, geometry));
    STEREON_TRY_VOID(builder.KillEdgeMakeRing(made.Edge));
    STEREON_TRY_VOID(CacheIfTrim(builder, closed.Edge, geometry));
    faces.All.push_back(closed.Face);
    return ImprintedPiece{closed.Edge, made.Vertex, made.Vertex, closed.Face};
}

} // namespace

Result<Imprint> ImprintClip(ShapeBuilder& builder, const Shape& shape, const FaceId face,
                            const IntersectionBranch& branch, const BranchSide side,
                            const FaceClip& clip, const Context& ctx)
{
    const BuilderFaceId handle = builder.Handle(face);
    if (!handle.IsValid()) {
        return Fail(ErrorCode::InvalidInput, "ImprintClip: the face is not live in the builder");
    }
    Imprint imprint;
    CrossingVertices vertices;
    STEREON_TRY_VOID(SplitCrossedEdges(builder, shape, branch, clip, vertices, imprint, ctx));
    Faces faces{{handle}};
    const Curve2d& pcurve = side == BranchSide::A ? branch.OnA : branch.OnB;
    for (const FacePiece& piece : clip.Inside) {
        STEREON_TRY(const EdgeGeometry geometry,
                    GeometryOf(builder, branch, pcurve, piece.Range, ctx));
        if (!piece.Start && !piece.End) {
            STEREON_TRY(ImprintedPiece made,
                        ImprintRing(builder, faces, branch, piece.Range, geometry));
            imprint.Pieces.push_back(made);
            continue;
        }
        if (!piece.Start || !piece.End) {
            return Fail(ErrorCode::NotSupported,
                        "ImprintClip: a piece of the branch ends inside the face; imprinting "
                        "needs its ends on the boundary (the Boolean prototype's M8 joins "
                        "branches at nodes)");
        }
        STEREON_TRY(ImprintedPiece made, ImprintOpen(builder, faces, vertices.at(*piece.Start),
                                                     vertices.at(*piece.End), geometry));
        imprint.Pieces.push_back(made);
    }
    return imprint;
}

} // namespace Stereon
