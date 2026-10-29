// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Parallel.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <set>
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
    // Qualified: an unqualified type would change the meaning of the name
    // the member declares (-Wchanges-meaning).
    Stereon::Surface Surface;
    Box3 Box;
    Box2 Parameters;
    /// The face's tolerance: the largest of its boundary's edges and
    /// vertices, and at least the resolution. A face carries no tolerance of
    /// its own (ADR-0001); its boundary says how loosely it was made.
    double Tolerance = 0.0;
};

[[nodiscard]] Result<double> ToleranceOfFace(const Shape& shape, const FaceId face,
                                             const Context& ctx)
{
    double tolerance = ctx.GetResolution();
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const EdgeId edge = EdgeOf(shape, coedge);
            tolerance = std::max({tolerance, ToleranceOf(shape, edge),
                                  ToleranceOf(shape, StartOf(shape, edge)),
                                  ToleranceOf(shape, EndOf(shape, edge))});
        }
    }
    return tolerance;
}

[[nodiscard]] Result<std::vector<FaceData>> FacesOf(const Shape& shape, const Context& ctx)
{
    std::vector<FaceData> faces;
    STEREON_TRY(const auto all, AllEntities<EntityKind::Face>(shape));
    for (const FaceId face : all) {
        STEREON_TRY(const Box3 box, BoxOf(shape, face));
        STEREON_TRY(const Box2 parameters, ParameterBoxOf(shape, face, ctx));
        STEREON_TRY(const double tolerance, ToleranceOfFace(shape, face, ctx));
        faces.push_back({face, SurfaceOf(shape, face), box, parameters, tolerance});
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

/// Whether a stretch of a branch over `range` is a loop: all of the branch,
/// with no crossing at either end, and the branch a closed curve, its ends
/// within `tolerance` of each other. A coincident face's straight side lying
/// wholly inside the other face covers all of its branch with no crossing
/// too, and is no loop.
[[nodiscard]] bool IsLoop(const IntersectionBranch& branch, const Interval& range,
                          const FacePiece& piece, const double tolerance)
{
    return range == branch.Range && !piece.Start && !piece.End &&
           Distance(Eval(branch.Curve3d, branch.Range.GetLower()),
                    Eval(branch.Curve3d, branch.Range.GetUpper())) <= tolerance;
}

/// A stretch of a branch inside one face only, added to that face's network:
/// its ends are the clip's crossings where it ends at one, junctions
/// otherwise. A coincident face's boundary imprinted on the other body's
/// face is one: the other body already has the edge.
void AddOneSided(Networks& networks, const Shape& shape, const FaceData& face,
                 const IntersectionBranch& branch, const FaceClip& clip, const FacePiece& piece,
                 const Context& ctx)
{
    const double tolerance =
        std::max(ctx.GetResolution(), branch.IsCertified() ? branch.Bound : 0.0);
    const auto end = [&](const std::optional<std::size_t>& at) {
        return at ? std::optional<LoopCrossing>{clip.Crossings[*at]} : std::nullopt;
    };
    const auto senseOf = [&](const std::optional<LoopCrossing>& crossing) {
        return crossing ? std::optional{StoredEdgeUseSense(shape, crossing->Coedge)} : std::nullopt;
    };
    const std::optional<LoopCrossing> atStart = end(piece.Start);
    const std::optional<LoopCrossing> atEnd = end(piece.End);
    const std::size_t start =
        networks.VertexOf(face.Face, branch, atStart, piece.Range.GetLower(), tolerance);
    const bool closed = IsLoop(branch, piece.Range, piece, tolerance);
    const std::size_t stop =
        closed ? start
               : networks.VertexOf(face.Face, branch, atEnd, piece.Range.GetUpper(), tolerance);
    networks.AddEdge(face.Face, {branch.Curve3d, branch.OnA, piece.Range, branch.Bound, start, stop,
                                 senseOf(atStart), senseOf(atEnd)});
}

/// A piece of the intersection: its edge in each body's networks, and the
/// branch's side transition.
struct SharedPiece {
    FaceId FaceA;
    std::size_t EdgeA;
    FaceId FaceB;
    std::size_t EdgeB;
    SideTransition Transition;
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
        const bool closed = IsLoop(branch, range, piece, tolerance);
        const std::size_t end =
            closed ? start : n.VertexOf(f.Face, branch, atEnd, range.GetUpper(), tolerance);
        return n.AddEdge(f.Face, {branch.Curve3d, pcurve, range, branch.Bound, start, end,
                                  senseOf(shape, atStart), senseOf(shape, atEnd)});
    };
    const std::size_t onA = edgeOn(a, sa, fa, clipA, pa, branch.OnA);
    const std::size_t onB = edgeOn(b, sb, fb, clipB, pb, branch.OnB);
    shared.push_back({fa.Face, onA, fb.Face, onB, branch.Transition});
}

/// Every piece of the branch inside both faces: the common stretches of the
/// two clips' inside pieces.
void AddBranch(Networks& a, Networks& b, const Shape& sa, const Shape& sb, const FaceData& fa,
               const FaceData& fb, const IntersectionBranch& branch, const FaceClip& clipA,
               const FaceClip& clipB, std::vector<SharedPiece>& shared, const Context& ctx)
{
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
}

// --- Coincident faces ------------------------------------------------------------------

/// The parameters of a point of the plane: its offsets along the axes, the
/// plane's own Eval inverted (the axes orthonormal within FrameAllowance).
[[nodiscard]] Point2 OnPlane(const Plane& plane, const Point3& p)
{
    const Vec3 d = p - plane.Origin;
    return {Dot(d, plane.XAxis), Dot(d, plane.YAxis)};
}

[[nodiscard]] Vec2 AlongPlane(const Plane& plane, const Vec3& v)
{
    return {Dot(v, plane.XAxis), Dot(v, plane.YAxis)};
}

/// The pcurve on `plane` of a 3D line, circle or ellipse lying in it, on the
/// curve's own parameter: the curve's frame expressed in the plane's.
/// Nothing for any other curve.
[[nodiscard]] std::optional<Curve2d> PCurveIn(const Plane& plane, const Curve& curve)
{
    if (const Line* l = curve.Target<Line>()) {
        return Curve2d(Line2d{OnPlane(plane, l->Origin), AlongPlane(plane, l->Direction)});
    }
    if (const Circle* c = curve.Target<Circle>()) {
        return Curve2d(Circle2d{OnPlane(plane, c->Center), c->Radius, AlongPlane(plane, c->XAxis),
                                AlongPlane(plane, c->YAxis)});
    }
    if (const Ellipse* e = curve.Target<Ellipse>()) {
        return Curve2d(Ellipse2d{OnPlane(plane, e->Center), e->MajorRadius, e->MinorRadius,
                                 AlongPlane(plane, e->XAxis), AlongPlane(plane, e->YAxis)});
    }
    return std::nullopt;
}

/// The branch an edge of one face makes on a coincident face of the other
/// body: its own curve, its pcurve on the coincident plane, the edge's own
/// pcurve as the second, over the edge's range, its bound the edge's
/// tolerance and the planes' distance. NotSupported for an edge without a
/// line, circle or ellipse of its own.
[[nodiscard]] Result<IntersectionBranch> BoundaryBranch(const Shape& shape, const CoedgeId coedge,
                                                        const Plane& onto, const double distance)
{
    const EdgeId edge = EdgeOf(shape, coedge);
    const std::optional<Curve> curve = CurveOf(shape, edge);
    const std::optional<Curve2d> mapped = curve ? PCurveIn(onto, *curve) : std::nullopt;
    const std::optional<Curve2d> own = PCurveOf(shape, coedge);
    if (!mapped || !own) {
        return Fail(ErrorCode::NotSupported,
                    "ImprintBodies: a coincident face is bounded by " +
                        std::string(Describe(shape, edge)) +
                        ", which has no line, circle or ellipse of its own to imprint");
    }
    IntersectionBranch branch;
    branch.Curve3d = *curve;
    branch.OnA = *mapped;
    branch.OnB = *own;
    branch.Range = RangeOf(shape, edge);
    branch.Contact = ContactType::Overlap;
    branch.Transition = SideTransition::Unresolved;
    branch.Bound = ToleranceOf(shape, edge) + distance;
    return branch;
}

/// Imprints the boundary of `from`, a face of `source`, onto `onto`, the
/// coincident face of `target` (a plane): every edge of every loop as a
/// branch, clipped to the face; the stretches inside it join the face's
/// network, the stretches along its own edges are left to the pairing.
[[nodiscard]] Result<void> ImprintBoundary(Networks& networks, const Shape& target,
                                           const FaceData& onto, const Shape& source,
                                           const FaceId from, const double distance,
                                           const Context& ctx)
{
    const Plane* plane = onto.Surface.Target<Plane>();
    STEREON_TRY(const LoopRange loops, LoopsOf(source, from));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(source, loop));
        for (const CoedgeId coedge : coedges) {
            STEREON_TRY(const IntersectionBranch branch,
                        BoundaryBranch(source, coedge, *plane, distance));
            STEREON_TRY(const FaceClip clip,
                        ClipToFace(target, onto.Face, branch, BranchSide::A, ctx));
            for (const FacePiece& piece : clip.Inside) {
                AddOneSided(networks, target, onto, branch, clip, piece, ctx);
            }
        }
    }
    return {};
}

/// A coincidence of two faces: NotSupported unless both are planes (the
/// prototype imprints coincident planar faces; coincident cylinders wait),
/// otherwise each face's boundary imprinted onto the other.
[[nodiscard]] Result<void> AddCoincidence(Networks& a, Networks& b, const Shape& sa,
                                          const Shape& sb, const FaceData& fa, const FaceData& fb,
                                          const OverlapRegion& overlap, std::vector<Plane>& planes,
                                          const Context& ctx)
{
    const Plane* pa = fa.Surface.Target<Plane>();
    if (pa == nullptr || fb.Surface.Target<Plane>() == nullptr) {
        return Fail(ErrorCode::NotSupported,
                    "ImprintBodies: coincident curved faces are not handled yet: " +
                        std::string(Describe(sa, fa.Face)) + " and " +
                        std::string(Describe(sb, fb.Face)));
    }
    STEREON_TRY_VOID(ImprintBoundary(a, sa, fa, sb, fb.Face, overlap.Distance, ctx));
    STEREON_TRY_VOID(ImprintBoundary(b, sb, fb, sa, fa.Face, overlap.Distance, ctx));
    planes.push_back(*pa);
    return {};
}

// --- Gathering, in parallel ----------------------------------------------------------

/// A candidate face pair, by its faces' places in the two bodies' lists, in
/// input order (pipeline step 1).
struct Candidate {
    std::size_t A;
    std::size_t B;
};

/// A log record kept to be passed on later, in the pair's order.
struct KeptLog {
    LogLevel Level;
    std::string Message;
    SourceLocation Location;
};

/// What a pair's own work found, written only by the chunk of its index
/// (ADR-0007, G1): the surface intersection, both clips of each branch that
/// bounds no coincidence, and what the pair logged.
struct PairWork {
    Result<SurfaceSurfaceIntersection> Meeting;
    std::set<std::size_t> Bounding;
    std::vector<std::optional<Result<std::pair<FaceClip, FaceClip>>>> Clips;
    std::vector<KeptLog> Logs;
};

/// The candidate pairs: the faces whose boxes meet, A's faces in order and
/// B's within each. The broad phase finds them without testing every pair,
/// in the order testing every pair would list them.
[[nodiscard]] std::vector<Candidate> CandidatesOf(const std::vector<FaceData>& facesA,
                                                  const std::vector<FaceData>& facesB)
{
    const auto boxes = [](const std::vector<FaceData>& faces) {
        std::vector<Box3> out;
        out.reserve(faces.size());
        for (const FaceData& face : faces) {
            out.push_back(face.Box);
        }
        return out;
    };
    std::vector<Candidate> candidates;
    for (const BoxPair& pair : OverlappingPairs(boxes(facesA), boxes(facesB))) {
        candidates.push_back({pair.A, pair.B});
    }
    return candidates;
}

/// Both clips of a branch.
[[nodiscard]] Result<std::pair<FaceClip, FaceClip>> ClipsOf(const Shape& sa, const Shape& sb,
                                                            const FaceData& fa, const FaceData& fb,
                                                            const IntersectionBranch& branch,
                                                            const Context& ctx)
{
    STEREON_TRY(FaceClip clipA, ClipToFace(sa, fa.Face, branch, BranchSide::A, ctx));
    STEREON_TRY(FaceClip clipB, ClipToFace(sb, fb.Face, branch, BranchSide::B, ctx));
    return std::pair{std::move(clipA), std::move(clipB)};
}

/// One pair's work: the intersection at the pair tolerance, the larger of
/// the two faces' (ADR-0001, decision D2), within which exactly parallel
/// faces are declared coincident, and the clips of the branches that bound
/// no coincidence (a coincidence's own boundary branches bound the faces'
/// parameter boxes, not the faces: the faces' boundaries imprinted on each
/// other take their place, after, in order). It reads the two frozen shapes
/// only.
void WorkOnPair(PairWork& work, const Shape& sa, const Shape& sb, const FaceData& fa,
                const FaceData& fb, const Context& ctx)
{
    work.Meeting = Intersect(fa.Surface, fa.Parameters, fb.Surface, fb.Parameters,
                             std::max(fa.Tolerance, fb.Tolerance), ctx);
    if (!work.Meeting) {
        return;
    }
    for (const OverlapRegion& overlap : work.Meeting->Overlaps) {
        work.Bounding.insert(overlap.Boundary.begin(), overlap.Boundary.end());
    }
    for (std::size_t i = 0; i < work.Meeting->Branches.size(); ++i) {
        work.Clips.emplace_back();
        if (!work.Bounding.contains(i)) {
            work.Clips.back() = ClipsOf(sa, sb, fa, fb, work.Meeting->Branches[i], ctx);
        }
    }
}

/// What the parallel loop reads, and the slots it writes.
struct Gathering {
    const Shape* A;
    const Shape* B;
    const std::vector<FaceData>* FacesA;
    const std::vector<FaceData>* FacesB;
    const std::vector<Candidate>* Candidates;
    const Context* Ctx;
    std::vector<PairWork>* Work;
};

/// The parallel loop's body: the pairs of one chunk, each with a Context of
/// its own that runs serially, nothing nested in the pool, and keeps its log
/// records. It must not throw: anything thrown is the pair's failure.
void WorkOnChunk(void* state, const std::size_t begin, const std::size_t end)
{
    const Gathering& g = *static_cast<const Gathering*>(state);
    for (std::size_t k = begin; k < end; ++k) {
        PairWork& work = (*g.Work)[k];
        try {
            const Context own =
                g.Ctx->WithScheduler(Scheduler::Serial)
                    .WithLogSink([&work](const LogRecord& record) {
                        work.Logs.push_back(
                            {record.Level, std::string(record.Message), record.Location});
                    });
            const Candidate& c = (*g.Candidates)[k];
            WorkOnPair(work, *g.A, *g.B, (*g.FacesA)[c.A], (*g.FacesB)[c.B], own);
        } catch (const std::exception& e) {
            work.Meeting =
                Fail(ErrorCode::IntersectionFailed,
                     std::string("ImprintBodies: a face pair's work threw: ") + e.what());
        }
    }
}

/// Intersects every candidate face pair and gathers the networks: each
/// pair's intersection and clips on the Context's scheduler (ADR-0007), one
/// pair per chunk since their work is large and uneven; then, in pair
/// order, its log records passed on, its failure returned, or its
/// coincidences and branches added to the networks, so that nothing depends
/// on the thread count.
[[nodiscard]] Result<void> Gather(Networks& a, Networks& b, const Shape& sa, const Shape& sb,
                                  std::vector<SharedPiece>& shared, std::vector<Plane>& planes,
                                  const Context& ctx)
{
    STEREON_TRY(const std::vector<FaceData> facesA, FacesOf(sa, ctx));
    STEREON_TRY(const std::vector<FaceData> facesB, FacesOf(sb, ctx));
    const std::vector<Candidate> candidates = CandidatesOf(facesA, facesB);
    std::vector<PairWork> work(candidates.size());
    Gathering gathering{&sa, &sb, &facesA, &facesB, &candidates, &ctx, &work};
    Detail::ParallelFor(candidates.size(), 1, ctx.GetWorkerPool(), &gathering, &WorkOnChunk);
    for (std::size_t k = 0; k < candidates.size(); ++k) {
        for (const KeptLog& log : work[k].Logs) {
            ctx.Log(log.Level, log.Message, log.Location);
        }
        STEREON_TRY(const SurfaceSurfaceIntersection& meeting, std::move(work[k].Meeting));
        const FaceData& fa = facesA[candidates[k].A];
        const FaceData& fb = facesB[candidates[k].B];
        for (const OverlapRegion& overlap : meeting.Overlaps) {
            STEREON_TRY_VOID(AddCoincidence(a, b, sa, sb, fa, fb, overlap, planes, ctx));
        }
        for (std::size_t i = 0; i < meeting.Branches.size(); ++i) {
            if (work[k].Bounding.contains(i)) {
                continue;
            }
            STEREON_TRY(const auto& clips, std::move(*work[k].Clips[i]));
            AddBranch(a, b, sa, sb, fa, fb, meeting.Branches[i], clips.first, clips.second, shared,
                      ctx);
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

/// The point of the edge at parameter t: on its own curve, or for a kind
/// that evaluates through a surface, on its primary side.
[[nodiscard]] std::optional<Point3> EdgePoint(const Shape& shape, const EdgeId edge, const double t)
{
    if (const std::optional<Curve> curve = CurveOf(shape, edge)) {
        return Eval(*curve, t);
    }
    const CoedgeId primary = PrimaryCoedgeOf(shape, edge);
    const std::optional<Curve2d> pcurve =
        primary.IsValid() ? PCurveOf(shape, primary) : std::nullopt;
    if (!pcurve) {
        return std::nullopt;
    }
    return Eval(SurfaceOf(shape, FaceOf(shape, LoopOf(shape, primary))), Eval(*pcurve, t));
}

/// An edge as the pairing sees it: its ends and middle in space, and its
/// tolerance.
struct EdgeTrace {
    EdgeId Edge;
    Point3 Start;
    Point3 Middle;
    Point3 End;
    double Tolerance = 0.0;
};

[[nodiscard]] std::optional<EdgeTrace> TraceOf(const Shape& shape, const EdgeId edge)
{
    const std::optional<Point3> middle = EdgePoint(shape, edge, Midpoint(RangeOf(shape, edge)));
    if (!middle) {
        return std::nullopt;
    }
    return EdgeTrace{edge, PositionOf(shape, StartOf(shape, edge)), *middle,
                     PositionOf(shape, EndOf(shape, edge)), ToleranceOf(shape, edge)};
}

/// Whether the trace lies on the plane: its ends and middle within its
/// tolerance of it.
[[nodiscard]] bool OnPlane(const EdgeTrace& trace, const Plane& plane)
{
    const Vec3 n = Normalized(Cross(plane.XAxis, plane.YAxis));
    return std::ranges::all_of(
        std::array{trace.Start, trace.Middle, trace.End},
        [&](const Point3& p) { return std::abs(Dot(p - plane.Origin, n)) <= trace.Tolerance; });
}

/// The edges of the body that lie on one of the coincident planes, but for
/// those a piece of the intersection already pairs.
[[nodiscard]] Result<std::vector<EdgeTrace>> OnPlanes(const Shape& shape,
                                                      const std::vector<Plane>& planes,
                                                      const std::set<std::uint32_t>& paired)
{
    std::vector<EdgeTrace> traces;
    STEREON_TRY(const auto edges, AllEntities<EntityKind::Edge>(shape));
    for (const EdgeId edge : edges) {
        const std::optional<EdgeTrace> trace =
            paired.contains(edge.GetIndex()) ? std::nullopt : TraceOf(shape, edge);
        if (trace &&
            std::ranges::any_of(planes, [&](const Plane& p) { return OnPlane(*trace, p); })) {
            traces.push_back(*trace);
        }
    }
    return traces;
}

/// Whether two traces are one edge: ends together, either way round, and
/// middles together, within the larger tolerance.
[[nodiscard]] bool SameEdge(const EdgeTrace& a, const EdgeTrace& b)
{
    const double tolerance = std::max(a.Tolerance, b.Tolerance);
    const auto near = [&](const Point3& p, const Point3& q) { return Distance(p, q) <= tolerance; };
    const bool along = near(a.Start, b.Start) && near(a.End, b.End);
    const bool against = near(a.Start, b.End) && near(a.End, b.Start);
    return (along || against) && near(a.Middle, b.Middle);
}

/// Pairs the edges of the two imprinted bodies that coincide on the planes
/// of coincident faces: each such edge of A with the first edge of B that is
/// the same edge. Their transition is Unresolved: an overlap's boundary
/// separates regions without saying which side they are on.
[[nodiscard]] Result<void> PairCoincidentEdges(const Shape& a, const Shape& b,
                                               const std::vector<Plane>& planes,
                                               std::vector<SharedEdge>& shared)
{
    if (planes.empty()) {
        return {};
    }
    std::set<std::uint32_t> pairedA;
    std::set<std::uint32_t> pairedB;
    for (const SharedEdge& s : shared) {
        pairedA.insert(s.A.GetIndex());
        pairedB.insert(s.B.GetIndex());
    }
    STEREON_TRY(const std::vector<EdgeTrace> inA, OnPlanes(a, planes, pairedA));
    STEREON_TRY(const std::vector<EdgeTrace> inB, OnPlanes(b, planes, pairedB));
    for (const EdgeTrace& ea : inA) {
        for (const EdgeTrace& eb : inB) {
            if (!pairedB.contains(eb.Edge.GetIndex()) && SameEdge(ea, eb)) {
                shared.push_back({ea.Edge, eb.Edge, SideTransition::Unresolved});
                pairedB.insert(eb.Edge.GetIndex());
                break;
            }
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
    std::vector<Plane> planes; // the planes of coincident faces
    STEREON_TRY_VOID(Gather(networksA, networksB, a, b, shared, planes, ctx));
    STEREON_TRY(const Imprinted bodyA, ImprintBody(a, networksA, ctx));
    STEREON_TRY(const Imprinted bodyB, ImprintBody(b, networksB, ctx));
    ImprintedBodies result{bodyA.Result, bodyB.Result, {}};
    for (const SharedPiece& piece : shared) {
        result.Shared.push_back({EdgeOfPiece(bodyA, piece.FaceA, piece.EdgeA),
                                 EdgeOfPiece(bodyB, piece.FaceB, piece.EdgeB), piece.Transition});
    }
    STEREON_TRY_VOID(PairCoincidentEdges(result.A.Shape, result.B.Shape, planes, result.Shared));
    return result;
}

} // namespace Stereon
