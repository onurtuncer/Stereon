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
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Clip;

namespace Stereon {

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();

/// How many directions a ray from a point tries before a side is given up:
/// a ray meeting a loop vertex is ambiguous and turns by an irrational
/// fraction of a turn.
constexpr int RayDirections = 8;

/// A coedge of the face with what the clip needs of it: its pcurve on the
/// face's surface over the edge's range, and the vertices at the range's
/// ends (the edge's, whatever the use's sense).
struct LoopEdge {
    CoedgeId Coedge;
    Curve2d PCurve;
    Interval Range;
    VertexId AtLower;
    VertexId AtUpper;
};

[[nodiscard]] Result<std::vector<LoopEdge>> LoopEdgesOf(const Shape& shape, const FaceId face)
{
    std::vector<LoopEdge> edges;
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const EdgeId edge = EdgeOf(shape, coedge);
            std::optional<Curve2d> pcurve = PCurveOf(shape, coedge);
            if (!pcurve) {
                return Fail(ErrorCode::InvalidInput,
                            "ClipToFace: a coedge of the face has no pcurve: " +
                                std::string(Describe(shape, coedge)));
            }
            edges.push_back({coedge, std::move(*pcurve), RangeOf(shape, edge), StartOf(shape, edge),
                             EndOf(shape, edge)});
        }
    }
    return edges;
}

/// The branch's pcurve on the clipped face.
[[nodiscard]] const Curve2d& PCurveOf(const IntersectionBranch& branch, const BranchSide side)
{
    return side == BranchSide::A ? branch.OnA : branch.OnB;
}

/// The crossing at the vertex within whose tolerance the branch's point
/// lies, its edge parameter snapped to the range's end (ADR-0001); the
/// crossing as it is otherwise.
void SnapToVertex(const Shape& shape, const IntersectionBranch& branch, const LoopEdge& edge,
                  LoopCrossing& crossing)
{
    const Point3 point = Eval(branch.Curve3d, crossing.OnBranch);
    for (const auto& [vertex, end] : {std::pair{edge.AtLower, edge.Range.GetLower()},
                                      std::pair{edge.AtUpper, edge.Range.GetUpper()}}) {
        if (vertex.IsValid() &&
            Distance(point, PositionOf(shape, vertex)) <= ToleranceOf(shape, vertex)) {
            crossing.Vertex = vertex;
            crossing.OnEdge = end;
            return;
        }
    }
}

/// Every meeting of the branch's pcurve with the loops' pcurves: the points
/// as crossings and the overlaps as stretches along the boundary.
[[nodiscard]] Result<void> CollectCrossings(const Shape& shape, const IntersectionBranch& branch,
                                            const Curve2d& pcurve,
                                            const std::vector<LoopEdge>& edges, FaceClip& clip,
                                            const Context& ctx)
{
    for (const LoopEdge& edge : edges) {
        STEREON_TRY(const PlanarIntersection meeting,
                    Intersect(pcurve, branch.Range, edge.PCurve, edge.Range, ctx));
        for (const PlanarPoint& p : meeting.Points) {
            LoopCrossing crossing{edge.Coedge, VertexId{}, p.ParameterB, p.ParameterA,
                                  p.Position,  p.Contact,  p.Decided};
            SnapToVertex(shape, branch, edge, crossing);
            clip.Decided = clip.Decided && p.Decided;
            clip.Crossings.push_back(crossing);
        }
        for (const PlanarOverlap& overlap : meeting.Overlaps) {
            clip.OnBoundary.push_back(overlap.OnA);
        }
    }
    return {};
}

/// The parameter of the edge's pcurve nearest p, within the edge's range: a
/// periodic parameter is placed in the range when the range holds one of
/// its images, and a parameter beyond the range is its nearer end.
[[nodiscard]] double NearestOn(const LoopEdge& edge, const Point2& p)
{
    const double s = Project(edge.PCurve, p);
    if (edge.PCurve.Target<Line2d>() == nullptr) {
        for (const double shift : {0.0, CirclePeriod, -CirclePeriod}) {
            if (Contains(edge.Range, s + shift)) {
                return s + shift;
            }
        }
    }
    return std::clamp(s, edge.Range.GetLower(), edge.Range.GetUpper());
}

/// Whether a crossing of the coedge within rounding of t is already known.
[[nodiscard]] bool HasCrossingAt(const FaceClip& clip, const CoedgeId coedge, const double t,
                                 const double sliver)
{
    return std::ranges::any_of(clip.Crossings, [&](const LoopCrossing& c) {
        return c.Coedge == coedge && std::abs(c.OnBranch - t) <= sliver;
    });
}

/// A branch end lying within an edge's tolerance of that edge is on it
/// (ADR-0001), whether or not the plane curve intersection found a root
/// there: a B-spline pcurve's root at the very end of its domain may be
/// missed. The end becomes a crossing of the edge, bracketed, not decided.
void AddEndCrossings(const Shape& shape, const Surface& surface, const IntersectionBranch& branch,
                     const Curve2d& pcurve, const std::vector<LoopEdge>& edges, FaceClip& clip)
{
    const double sliver =
        1e-9 * Width(branch.Range); // STEREON-LINT-ALLOW(epsilon): a parameter gap
    for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
        const Point2 p = Eval(pcurve, t);
        const Point3 point = Eval(branch.Curve3d, t);
        for (const LoopEdge& edge : edges) {
            if (HasCrossingAt(clip, edge.Coedge, t, sliver)) {
                continue;
            }
            const double s = NearestOn(edge, p);
            const Point2 onEdge = Eval(edge.PCurve, s);
            const double tolerance = std::max(ToleranceOf(shape, EdgeOf(shape, edge.Coedge)),
                                              branch.IsCertified() ? branch.Bound : 0.0);
            if (Distance(point, Eval(surface, onEdge)) <= tolerance) {
                LoopCrossing crossing{edge.Coedge,           VertexId{}, s, t, onEdge,
                                      ContactType::Crossing, false};
                SnapToVertex(shape, branch, edge, crossing);
                clip.Decided = false;
                clip.Crossings.push_back(crossing);
            }
        }
    }
}

/// Crossings in order along the branch; the two edges meeting at a vertex
/// the branch passes through report it once each, which is one crossing.
void OrderCrossings(FaceClip& clip)
{
    std::ranges::sort(clip.Crossings, [](const LoopCrossing& p, const LoopCrossing& q) {
        return p.OnBranch < q.OnBranch;
    });
    std::vector<LoopCrossing> distinct;
    for (const LoopCrossing& c : clip.Crossings) {
        const bool sameVertex = !distinct.empty() && c.Vertex.IsValid() &&
                                c.Vertex == distinct.back().Vertex &&
                                c.OnBranch == distinct.back().OnBranch;
        if (!sameVertex) {
            distinct.push_back(c);
        }
    }
    clip.Crossings = std::move(distinct);
}

/// The ray from q in the direction of the given turn, as a plane curve over
/// [0, infinity).
[[nodiscard]] Curve2d RayFrom(const Point2& q, const int turn)
{
    // Irrational fractions of a turn, so that no two rays align with a
    // loop's edges in the same way.
    const double angle = static_cast<double>(turn) * std::numbers::phi * std::numbers::pi;
    return Curve2d{Line2d{q, Vec2{Cos(angle), Sin(angle)}}};
}

/// Whether the ray meets the edge's pcurve at one of its ends, where the
/// neighbouring edge reports the same point: ambiguous for a parity count.
[[nodiscard]] bool AtAnEnd(const PlanarPoint& p, const LoopEdge& edge)
{
    return p.ParameterB == edge.Range.GetLower() || p.ParameterB == edge.Range.GetUpper();
}

/// The parity of the ray's crossings of the loops, or nothing when the ray
/// meets a vertex; a tangential touch does not change the side. `decided`
/// is cleared when a crossing was bracketed to rounding only.
[[nodiscard]] Result<std::optional<bool>>
RayParity(const Curve2d& ray, const std::vector<LoopEdge>& edges, bool& decided, const Context& ctx)
{
    std::size_t crossings = 0;
    for (const LoopEdge& edge : edges) {
        STEREON_TRY(const PlanarIntersection meeting,
                    Intersect(ray, Interval{0.0, Infinity}, edge.PCurve, edge.Range, ctx));
        if (!meeting.Overlaps.empty()) {
            return std::optional<bool>{}; // Along an edge: another direction.
        }
        for (const PlanarPoint& p : meeting.Points) {
            // The ray's own origin on the edge: q is on the boundary.
            if (AtAnEnd(p, edge) || p.ParameterA == 0.0) { // STEREON-LINT-ALLOW(epsilon): exact
                return std::optional<bool>{};
            }
            if (p.Contact == ContactType::Crossing) {
                ++crossings;
                decided = decided && p.Decided;
            }
        }
    }
    return std::optional<bool>{crossings % 2 == 1};
}

/// Whether q lies inside the face's loops, by the parity of a ray's
/// crossings, the ray turned until it meets no vertex.
[[nodiscard]] Result<bool> InsideLoops(const Point2& q, const std::vector<LoopEdge>& edges,
                                       bool& decided, const Context& ctx)
{
    for (int turn = 0; turn < RayDirections; ++turn) {
        STEREON_TRY(const std::optional<bool> inside,
                    RayParity(RayFrom(q, turn), edges, decided, ctx));
        if (inside) {
            return *inside;
        }
    }
    return Fail(ErrorCode::IntersectionFailed,
                "ClipToFace: every ray from a point of the branch meets a vertex of the face's "
                "loops; the side cannot be settled");
}

/// Whether t lies within a stretch where the branch runs along the boundary.
[[nodiscard]] bool OnBoundary(const FaceClip& clip, const double t)
{
    return std::ranges::any_of(clip.OnBoundary,
                               [&](const Interval& stretch) { return Contains(stretch, t); });
}

/// The stretches between consecutive crossings (and the branch's ends) that
/// lie inside the face, each by the side of its middle point.
[[nodiscard]] Result<void> CollectPieces(const IntersectionBranch& branch, const Curve2d& pcurve,
                                         const std::vector<LoopEdge>& edges, FaceClip& clip,
                                         const Context& ctx)
{
    // A crossing within rounding of the one before, or of the branch's end,
    // cuts no piece of its own: a closed branch meets the seam it starts on
    // at both its ends, found to rounding.
    const double sliver =
        1e-9 * Width(branch.Range); // STEREON-LINT-ALLOW(epsilon): a parameter gap
    // A crossing at the branch's own end is that end: the piece starting or
    // ending there is on the boundary, not inside.
    std::vector<std::optional<std::size_t>> cuts{std::nullopt};
    std::vector<double> parameters{branch.Range.GetLower()};
    std::optional<std::size_t> last;
    for (std::size_t i = 0; i < clip.Crossings.size(); ++i) {
        const double t = clip.Crossings[i].OnBranch;
        if (clip.Crossings[i].Contact != ContactType::Crossing) {
            continue;
        }
        if (t - branch.Range.GetLower() <= sliver) {
            cuts.front() = i;
        } else if (branch.Range.GetUpper() - t <= sliver) {
            last = i;
        } else if (t - parameters.back() > sliver) {
            cuts.push_back(i);
            parameters.push_back(t);
        }
    }
    cuts.push_back(last);
    parameters.push_back(branch.Range.GetUpper());
    for (std::size_t i = 1; i < parameters.size(); ++i) {
        const Interval piece{parameters[i - 1], parameters[i]};
        const double middle = 0.5 * (piece.GetLower() + piece.GetUpper());
        if (piece.GetLower() == piece.GetUpper() || OnBoundary(clip, middle)) {
            continue;
        }
        STEREON_TRY(const bool inside, InsideLoops(Eval(pcurve, middle), edges, clip.Decided, ctx));
        if (inside) {
            clip.Inside.push_back({piece, cuts[i - 1], cuts[i]});
        }
    }
    return {};
}

} // namespace

Result<FaceClip> ClipToFace(const Shape& shape, const FaceId face, const IntersectionBranch& branch,
                            const BranchSide side, const Context& ctx)
{
    if (!branch.Range.IsBounded() || branch.Range.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "ClipToFace: the branch's range must be bounded");
    }
    STEREON_TRY(const std::vector<LoopEdge> edges, LoopEdgesOf(shape, face));
    const Curve2d& pcurve = PCurveOf(branch, side);
    FaceClip clip;
    STEREON_TRY_VOID(CollectCrossings(shape, branch, pcurve, edges, clip, ctx));
    AddEndCrossings(shape, SurfaceOf(shape, face), branch, pcurve, edges, clip);
    OrderCrossings(clip);
    STEREON_TRY_VOID(CollectPieces(branch, pcurve, edges, clip, ctx));
    return clip;
}

} // namespace Stereon
