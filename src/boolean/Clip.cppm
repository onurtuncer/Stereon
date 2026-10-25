// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

export module Stereon.Boolean:Clip;

import Stereon.Build;
import Stereon.Intersect;

/// Face clipping (the Boolean prototype's M6, pipeline step 3): an
/// intersection branch of two surfaces restricted to a face on one of them,
/// in that face's parameter space. The branch's pcurve on the face is met
/// with the pcurve of every coedge of the face's loops (Intersect on plane
/// curves), the meetings become crossings of the boundary, in order along
/// the branch, and the stretches between them are inside the face or
/// outside it by the parity of a ray from their middle. A crossing within
/// a loop vertex's tolerance of that vertex is at the vertex (ADR-0001);
/// everything else is decided on the nominal geometry, exactly where the
/// pcurves are lines, circles and ellipses (decision D1 of the plan).

export namespace Stereon {

/// Which of a branch's two pcurves lies on the face being clipped.
enum class BranchSide : std::uint8_t { A, B };

/// Where a branch meets the boundary of a face: the coedge whose edge it
/// meets, the loop vertex it meets at (invalid inside the edge), the
/// parameter on the edge and on the branch, the point in the face's
/// parameter space, how the curves meet, and whether the meeting was
/// decided exactly on the input doubles.
struct LoopCrossing {
    CoedgeId Coedge;
    VertexId Vertex;
    double OnEdge = 0.0;
    double OnBranch = 0.0;
    Point2 Position;
    ContactType Contact = ContactType::Crossing;
    bool Decided = true;
};

/// A stretch of a branch inside a face: its parameter range, and the
/// crossings at its ends as indices into the clip's Crossings; none at an
/// end of the branch itself, or where a closed branch lies wholly inside.
struct FacePiece {
    Interval Range;
    std::optional<std::size_t> Start;
    std::optional<std::size_t> End;
};

/// A branch restricted to a face: its crossings of the face's boundary in
/// order along the branch, the stretches inside the face, the stretches
/// along a loop edge (where a pcurve coincides with the edge's), and
/// whether every crossing and every side was decided exactly.
struct FaceClip {
    std::vector<LoopCrossing> Crossings;
    std::vector<FacePiece> Inside;
    std::vector<Interval> OnBoundary;
    bool Decided = true;
};

/// Whether the face contains the point at the parameters uv: inside its
/// loops, by the parity of a ray turned past the loops' vertices; nothing
/// when the point lies on the boundary, within an edge's tolerance of the
/// edge (ADR-0001). IntersectionFailed when every ray meets a vertex.
[[nodiscard]] Result<std::optional<bool>> FaceContains(const Shape& shape, FaceId face,
                                                       const Point2& uv, const Context& ctx);

/// The branch restricted to the face, whose surface is the one the pcurve
/// of `side` lies on. InvalidInput for a face the shape does not hold, or a
/// coedge without a pcurve; IntersectionFailed when a side cannot be
/// settled, a ray through the face's loops meeting a vertex in every
/// direction tried; otherwise what the plane curve intersections report.
[[nodiscard]] Result<FaceClip> ClipToFace(const Shape& shape, FaceId face,
                                          const IntersectionBranch& branch, BranchSide side,
                                          const Context& ctx);

} // namespace Stereon
