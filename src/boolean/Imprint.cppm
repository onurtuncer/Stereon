// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <vector>

export module Stereon.Boolean:Imprint;

import :Clip;

/// Imprinting (the Boolean prototype's M6, pipeline step 4): the curve
/// network a face receives from the intersection, written into the face
/// through the controlled operators of ADR-0027. The network's vertices are
/// crossings of the face's boundary, or junctions inside the face where the
/// intersection passes from one face of the other body to the next; its
/// edges are the pieces of branches between them. The edges the network
/// crosses are split at the crossings (SEMV); a piece reaching a junction
/// not yet made is added as a strut (MEV); a piece whose ends both exist
/// splits the face (MEF) or bridges two of its loops (MEKR); a component
/// that touches no boundary is hung on a helper strut from a loop vertex,
/// removed once the component has closed (KEMR). The new edges are Exact
/// where the branch's 3D curve is a line, a circle or an ellipse, and
/// otherwise Trim edges with their cached curves (ADR-0009, section 8),
/// their tolerance the larger of the resolution and the branch's certified
/// bound.

export namespace Stereon {

/// A vertex of a face's curve network: where it is, its tolerance, and the
/// crossing of the face's boundary it lies at; a junction inside the face
/// when there is none.
struct NetworkVertex {
    Point3 Position;
    double Tolerance = 0.0;
    std::optional<LoopCrossing> Crossing;
};

/// A piece of a branch inside the face, from the vertex Start to the
/// vertex End (one vertex, when closed): its 3D curve and its pcurve on the
/// face over the range, the branch's certified bound, and at an end on the
/// boundary the sense of the coedge the clip names as crossed there. Where
/// both uses of one edge end at the vertex (a seam) the clip cannot tell
/// them apart, so the imprint takes the use the pcurve's own end reaches,
/// and these senses are kept for the record.
struct NetworkEdge {
    Curve Curve3d;
    Curve2d PCurve;
    Interval Range;
    double Bound = std::numeric_limits<double>::infinity();
    std::size_t Start = 0;
    std::size_t End = 0;
    std::optional<EdgeUseSense> StartSense;
    std::optional<EdgeUseSense> EndSense;
};

struct FaceNetwork {
    FaceId Face;
    std::vector<NetworkVertex> Vertices;
    std::vector<NetworkEdge> Edges;
};

/// What a piece became: its edge, the vertices at its ends, and the face it
/// cut off the original face, invalid when it bridged two loops or reached
/// a new junction instead.
struct ImprintedPiece {
    BuilderEdgeId Edge;
    BuilderVertexId Start;
    BuilderVertexId End;
    BuilderFaceId Face;
};

/// A network imprinted: its pieces in the network's edge order, the
/// builder's vertex for each network vertex, and the vertices the splits
/// made where the network cut through an edge.
struct Imprint {
    std::vector<ImprintedPiece> Pieces;
    std::vector<BuilderVertexId> Vertices;
    std::vector<BuilderVertexId> SplitVertices;
};

/// One clip as a network: every inside piece with its ends at the clip's
/// crossings; a piece that ends inside the face, or a closed one, gets a
/// junction there.
[[nodiscard]] FaceNetwork NetworkOf(const Shape& shape, FaceId face,
                                    const IntersectionBranch& branch, BranchSide side,
                                    const FaceClip& clip, const Context& ctx);

/// Imprints the networks of a body's faces into `builder`, which was opened
/// on `shape` (ShapeBuilder::Open), inside the caller's transaction, one
/// result per network in order. An edge two faces share is split once for
/// the crossings of both networks, which meet at the same vertex.
/// IntersectionFailed when a piece's ends do not lie on one face's
/// boundary, or a component never closes.
[[nodiscard]] Result<std::vector<Imprint>> ImprintNetworks(ShapeBuilder& builder,
                                                           const Shape& shape,
                                                           std::span<const FaceNetwork> networks,
                                                           const Context& ctx);

/// One network, imprinted.
[[nodiscard]] Result<Imprint> ImprintNetwork(ShapeBuilder& builder, const Shape& shape,
                                             const FaceNetwork& network, const Context& ctx);

/// The clip of `branch` on `face` as a network, imprinted.
[[nodiscard]] Result<Imprint> ImprintClip(ShapeBuilder& builder, const Shape& shape, FaceId face,
                                          const IntersectionBranch& branch, BranchSide side,
                                          const FaceClip& clip, const Context& ctx);

} // namespace Stereon
