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
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Bodies;

namespace Stereon {

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();

/// The box a face lies in: the hull of its vertices, each by its tolerance,
/// and of the bounding boxes of its edges' curves where those are bounded
/// (a circle's; a line's is everything, and its vertices bound it).
[[nodiscard]] Result<Box3> BoxOf(const Shape& shape, const FaceId face)
{
    Box3 box = Box3::Empty();
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const EdgeId edge = EdgeOf(shape, coedge);
            for (const VertexId v : {StartOf(shape, edge), EndOf(shape, edge)}) {
                const Point3 p = PositionOf(shape, v);
                const double tolerance = ToleranceOf(shape, v);
                box = Hull(box, Box3::FromCorners(p - Vec3{tolerance, tolerance, tolerance},
                                                  p + Vec3{tolerance, tolerance, tolerance}));
            }
            if (const std::optional<Curve> curve = CurveOf(shape, edge)) {
                const Box3 of = BoundingBox(*curve);
                if (of.X.IsBounded() && of.Y.IsBounded() && of.Z.IsBounded()) {
                    box = Hull(box, of);
                }
            }
        }
    }
    return box;
}

/// The face's parameter box: the hull of its pcurves over their ranges, a
/// line's by its ends, a circle's, an ellipse's or a B-spline's by its
/// bounding box, each widened by the resolution.
[[nodiscard]] Result<Box2> ParameterBoxOf(const Shape& shape, const FaceId face, const Context& ctx)
{
    Box2 box = Box2::Empty();
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const std::optional<Curve2d> pcurve = PCurveOf(shape, coedge);
            if (!pcurve) {
                return Fail(ErrorCode::InvalidInput, "ImprintBodies: a coedge has no pcurve: " +
                                                         std::string(Describe(shape, coedge)));
            }
            const Interval range = RangeOf(shape, EdgeOf(shape, coedge));
            if (pcurve->Target<Line2d>() != nullptr) {
                box = Hull(Hull(box, Eval(*pcurve, range.GetLower())),
                           Eval(*pcurve, range.GetUpper()));
            } else {
                box = Hull(box, BoundingBox(*pcurve));
            }
        }
    }
    const double margin = ctx.GetResolution();
    return Box2{Inflated(box.X, margin), Inflated(box.Y, margin)};
}

/// A face of one body with what the pairing needs of it.
struct FaceData {
    FaceId Face;
    Surface Surface;
    Box3 Box;
    Box2 Parameters;
};

[[nodiscard]] Result<std::vector<FaceData>> FacesOf(const Shape& shape, const Context& ctx)
{
    std::vector<FaceData> faces;
    STEREON_TRY(const auto all, AllEntities<EntityKind::Face>(shape));
    for (const FaceId face : all) {
        STEREON_TRY(const Box3 box, BoxOf(shape, face));
        STEREON_TRY(const Box2 parameters, ParameterBoxOf(shape, face, ctx));
        faces.push_back({face, SurfaceOf(shape, face), box, parameters});
    }
    return faces;
}

/// The networks a body's faces receive, built up pair by pair: a junction
/// within tolerance of one already in the face's network is that junction.
class Networks {
public:
    explicit Networks(const Shape& shape)
        : m_Shape(shape)
    {}

    /// The vertex of `face`'s network for the branch's point at t: the
    /// crossing's, or a junction merged by position.
    [[nodiscard]] std::size_t VertexOf(const FaceId face, const IntersectionBranch& branch,
                                       const std::optional<LoopCrossing>& crossing, const double t,
                                       const double tolerance)
    {
        FaceNetwork& network = Of(face);
        const Point3 position = Eval(branch.Curve3d, crossing ? crossing->OnBranch : t);
        for (std::size_t i = 0; i < network.Vertices.size(); ++i) {
            const NetworkVertex& v = network.Vertices[i];
            if (static_cast<bool>(v.Crossing) == static_cast<bool>(crossing) &&
                Distance(v.Position, position) <= std::max(v.Tolerance, tolerance)) {
                return i;
            }
        }
        network.Vertices.push_back({position, tolerance, crossing});
        return network.Vertices.size() - 1;
    }

    /// Adds a piece to `face`'s network; its index among the face's edges.
    std::size_t AddEdge(const FaceId face, NetworkEdge edge)
    {
        FaceNetwork& network = Of(face);
        network.Edges.push_back(std::move(edge));
        return network.Edges.size() - 1;
    }

    [[nodiscard]] std::vector<FaceNetwork>& All() { return m_Networks; }

private:
    FaceNetwork& Of(const FaceId face)
    {
        for (FaceNetwork& network : m_Networks) {
            if (network.Face == face) {
                return network;
            }
        }
        m_Networks.push_back({face, {}, {}});
        return m_Networks.back();
    }

    const Shape& m_Shape;
    std::vector<FaceNetwork> m_Networks;
};

/// A piece of the intersection: its edge in each body's networks.
struct SharedPiece {
    FaceId FaceA;
    std::size_t EdgeA;
    FaceId FaceB;
    std::size_t EdgeB;
};

/// The piece of a branch over `range`, inside both faces, added to both
/// networks: its ends are crossings where the range ends at a crossing of
/// that face, junctions otherwise.
void AddPiece(Networks& a, Networks& b, const Shape& sa, const Shape& sb, const FaceData& fa,
              const FaceData& fb, const IntersectionBranch& branch, const FaceClip& clipA,
              const FaceClip& clipB, const FacePiece& pa, const FacePiece& pb,
              const Interval& range, std::vector<SharedPiece>& shared, const Context& ctx)
{
    const double tolerance =
        std::max(ctx.GetResolution(), branch.IsCertified() ? branch.Bound : 0.0);
    const auto endOf = [&](const FaceClip& clip, const FacePiece& piece, const double t,
                           const bool lower) -> std::optional<LoopCrossing> {
        const std::optional<std::size_t>& at = lower ? piece.Start : piece.End;
        const double own = lower ? piece.Range.GetLower() : piece.Range.GetUpper();
        return at && t == own ? std::optional<LoopCrossing>{clip.Crossings[*at]} : std::nullopt;
    };
    const auto senseOf = [](const Shape& shape, const std::optional<LoopCrossing>& crossing) {
        return crossing ? std::optional{StoredEdgeUseSense(shape, crossing->Coedge)} : std::nullopt;
    };
    const auto edgeOn = [&](Networks& n, const Shape& shape, const FaceData& f,
                            const FaceClip& clip, const FacePiece& piece, const Curve2d& pcurve) {
        const std::optional<LoopCrossing> atStart = endOf(clip, piece, range.GetLower(), true);
        const std::optional<LoopCrossing> atEnd = endOf(clip, piece, range.GetUpper(), false);
        const std::size_t start = n.VertexOf(f.Face, branch, atStart, range.GetLower(), tolerance);
        const bool closed = range == branch.Range && !piece.Start && !piece.End;
        const std::size_t end =
            closed ? start : n.VertexOf(f.Face, branch, atEnd, range.GetUpper(), tolerance);
        return n.AddEdge(f.Face, {branch.Curve3d, pcurve, range, branch.Bound, start, end,
                                  senseOf(shape, atStart), senseOf(shape, atEnd)});
    };
    const std::size_t onA = edgeOn(a, sa, fa, clipA, pa, branch.OnA);
    const std::size_t onB = edgeOn(b, sb, fb, clipB, pb, branch.OnB);
    shared.push_back({fa.Face, onA, fb.Face, onB});
}

/// Every piece of the branch inside both faces: the common stretches of the
/// two clips' inside pieces.
[[nodiscard]] Result<void> AddBranch(Networks& a, Networks& b, const Shape& sa, const Shape& sb,
                                     const FaceData& fa, const FaceData& fb,
                                     const IntersectionBranch& branch,
                                     std::vector<SharedPiece>& shared, const Context& ctx)
{
    STEREON_TRY(const FaceClip clipA, ClipToFace(sa, fa.Face, branch, BranchSide::A, ctx));
    STEREON_TRY(const FaceClip clipB, ClipToFace(sb, fb.Face, branch, BranchSide::B, ctx));
    const double sliver =
        1e-9 * Width(branch.Range); // STEREON-LINT-ALLOW(epsilon): a parameter gap
    for (const FacePiece& pa : clipA.Inside) {
        for (const FacePiece& pb : clipB.Inside) {
            const Interval common = Intersection(pa.Range, pb.Range);
            if (common.IsEmpty() || Width(common) <= sliver) {
                continue;
            }
            AddPiece(a, b, sa, sb, fa, fb, branch, clipA, clipB, pa, pb, common, shared, ctx);
        }
    }
    return {};
}

/// Intersects every candidate face pair and gathers the networks.
[[nodiscard]] Result<void> Gather(Networks& a, Networks& b, const Shape& sa, const Shape& sb,
                                  std::vector<SharedPiece>& shared, const Context& ctx)
{
    STEREON_TRY(const std::vector<FaceData> facesA, FacesOf(sa, ctx));
    STEREON_TRY(const std::vector<FaceData> facesB, FacesOf(sb, ctx));
    for (const FaceData& fa : facesA) {
        for (const FaceData& fb : facesB) {
            if (!Intersects(fa.Box, fb.Box)) {
                continue;
            }
            STEREON_TRY(const SurfaceSurfaceIntersection meeting,
                        Intersect(fa.Surface, fa.Parameters, fb.Surface, fb.Parameters,
                                  ctx.GetResolution(), ctx));
            if (!meeting.Overlaps.empty()) {
                return Fail(ErrorCode::NotSupported,
                            "ImprintBodies: coincident faces (an overlap) wait for the assembly "
                            "of M8: " +
                                std::string(Describe(sa, fa.Face)) + " and " +
                                std::string(Describe(sb, fb.Face)));
            }
            for (const IntersectionBranch& branch : meeting.Branches) {
                STEREON_TRY_VOID(AddBranch(a, b, sa, sb, fa, fb, branch, shared, ctx));
            }
        }
    }
    return {};
}

/// The imprints of a body's networks, committed.
struct Imprinted {
    OperationResult Result;
    std::vector<std::pair<FaceId, Imprint>> ByFace;
};

[[nodiscard]] Result<Imprinted> ImprintBody(const Shape& shape, Networks& networks,
                                            const Context& ctx)
{
    STEREON_TRY(ShapeBuilder builder, ShapeBuilder::Open(shape));
    std::vector<std::pair<FaceId, Imprint>> byFace;
    STEREON_TRY_VOID(builder.Transaction([&](ShapeBuilder& candidate) -> Result<void> {
        byFace.clear();
        STEREON_TRY(std::vector<Imprint> made,
                    ImprintNetworks(candidate, shape, networks.All(), ctx));
        for (std::size_t n = 0; n < made.size(); ++n) {
            byFace.emplace_back(networks.All()[n].Face, std::move(made[n]));
        }
        return {};
    }));
    STEREON_TRY(OperationResult result, Commit(std::move(builder), ctx));
    return Imprinted{std::move(result), std::move(byFace)};
}

/// The frozen edge of a piece in an imprinted body.
[[nodiscard]] EdgeId EdgeOfPiece(const Imprinted& body, const FaceId face, const std::size_t edge)
{
    for (const auto& [imprintedFace, imprint] : body.ByFace) {
        if (imprintedFace == face) {
            return body.Result.Shape.Resolve(imprint.Pieces[edge].Edge);
        }
    }
    return {};
}

} // namespace

Result<ImprintedBodies> ImprintBodies(const Shape& a, const Shape& b, const Context& ctx)
{
    Networks networksA(a);
    Networks networksB(b);
    std::vector<SharedPiece> shared;
    STEREON_TRY_VOID(Gather(networksA, networksB, a, b, shared, ctx));
    STEREON_TRY(const Imprinted bodyA, ImprintBody(a, networksA, ctx));
    STEREON_TRY(const Imprinted bodyB, ImprintBody(b, networksB, ctx));
    ImprintedBodies result{bodyA.Result, bodyB.Result, {}};
    for (const SharedPiece& piece : shared) {
        result.Shared.emplace_back(EdgeOfPiece(bodyA, piece.FaceA, piece.EdgeA),
                                   EdgeOfPiece(bodyB, piece.FaceB, piece.EdgeB));
    }
    return result;
}

} // namespace Stereon
