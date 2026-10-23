// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <vector>

export module Stereon.Boolean:Imprint;

import :Clip;

/// Imprinting (the Boolean prototype's M6, pipeline step 4): a clipped
/// branch written into a face through the controlled operators of
/// ADR-0027. The edges the branch crosses are split at the crossings
/// (SEMV), and each piece inside the face becomes an edge: splitting the
/// face in two when its ends lie on one loop (MEF), bridging two loops of
/// the face when they lie on different ones (MEKR), or a ring inside the
/// face when the piece is closed (a strut by MEV, the closed edge by MEF,
/// the strut removed by KEMR). The new edge is Exact where the branch's
/// 3D curve is a line, a circle or an ellipse, and otherwise a Trim edge
/// with its cached curve (ADR-0009, section 8), its tolerance the larger
/// of the resolution and the branch's certified bound.

export namespace Stereon {

/// What a piece became: its edge, the vertices at its ends (one vertex for
/// a ring), and the face it cut off the original face, invalid when it
/// bridged two loops instead.
struct ImprintedPiece {
    BuilderEdgeId Edge;
    BuilderVertexId Start;
    BuilderVertexId End;
    BuilderFaceId Face;
};

/// A clip imprinted: its pieces in the clip's order, and the vertices the
/// crossings made where the branch cut through an edge.
struct Imprint {
    std::vector<ImprintedPiece> Pieces;
    std::vector<BuilderVertexId> SplitVertices;
};

/// Imprints the clip of `branch` on `face` into `builder`, which was opened
/// on `shape` (ShapeBuilder::Open), inside the caller's transaction. Every
/// piece of the clip must end on the face's boundary or be closed:
/// NotSupported for a piece that ends inside the face. The faces the
/// pieces make are found through the builder, so several pieces of one
/// branch on one face imprint in order.
[[nodiscard]] Result<Imprint> ImprintClip(ShapeBuilder& builder, const Shape& shape, FaceId face,
                                          const IntersectionBranch& branch, BranchSide side,
                                          const FaceClip& clip, const Context& ctx);

} // namespace Stereon
