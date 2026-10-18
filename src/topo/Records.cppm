// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <compare>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <utility>

export module Stereon.Topo:Records;

import Stereon.Intersect;

/// The logical B-rep model of ADR-0004 (Golovanov's solid -> shell -> face ->
/// loop -> coedge -> edge -> vertex) as one record table per entity kind,
/// each a PagedArena so that a snapshot and its edit share every page the
/// edit does not touch (ADR-0005). Records reach one another by typed index;
/// geometry lives in separate pools and is referenced by index, so several
/// entities can share one curve or surface (ADR-0025).
///
/// The records are the storage of the topology module. Only ShapeBuilder's
/// controlled operations write them (ADR-0027), only the validator reads
/// them unchecked, and the public interface reads them through accessors
/// and ranges on a frozen Shape.

export namespace Stereon {

/// What defines an edge's geometry (ADR-0009, section 1). The number of
/// coedges on the edge is a separate, topological fact; the validator
/// enforces the combinations ADR-0009 allows.
enum class EdgeCurveKind : std::uint8_t {
    /// Two distinct surfaces and a pcurve on each.
    Intersection,
    /// A closed-form 3D curve, with exact pcurves on the adjacent surfaces.
    Exact,
    /// One pcurve on one surface, shared by both adjacent faces.
    Trim,
    /// One surface and one pcurve: a sheet boundary or a periodic seam.
    Boundary,
    /// An independent 3D curve; pcurves, where present, are projections.
    Spatial,
    /// A pole or apex: a pcurve with no 3D extent.
    Degenerate,
};

[[nodiscard]] constexpr std::string_view ToString(const EdgeCurveKind kind) noexcept
{
    switch (kind) {
    case EdgeCurveKind::Intersection:
        return "Intersection";
    case EdgeCurveKind::Exact:
        return "Exact";
    case EdgeCurveKind::Trim:
        return "Trim";
    case EdgeCurveKind::Boundary:
        return "Boundary";
    case EdgeCurveKind::Spatial:
        return "Spatial";
    case EdgeCurveKind::Degenerate:
        return "Degenerate";
    }
    return "EdgeCurveKind";
}

/// How a coedge uses its edge: along the edge's parameter or against it.
/// This is the stored edge-use sense of ADR-0037; a traversal that walks a
/// loop backwards reverses it once more, through the orientation helpers,
/// never by flipping this value.
enum class EdgeUseSense : std::uint8_t {
    Same,
    Opposite,
};

/// An index into one of a builder's or shape's geometry pools. Tag keeps
/// curves, plane curves and surfaces apart.
template <typename Tag>
class GeometryRef {
public:
    static constexpr std::uint32_t Invalid = 0xFFFFFFFFu;

    constexpr GeometryRef() = default;
    [[nodiscard]] static constexpr GeometryRef FromIndex(const std::uint32_t index) noexcept
    {
        GeometryRef ref;
        ref.m_Index = index;
        return ref;
    }

    [[nodiscard]] constexpr std::uint32_t GetIndex() const noexcept { return m_Index; }
    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Index != Invalid; }

    // Written out, as Id's are: MSVC cannot default <=> with a declared
    // return type in a module interface (C3774 on std::partial_ordering).
    constexpr bool operator==(const GeometryRef& other) const noexcept
    {
        return m_Index == other.m_Index;
    }
    constexpr std::strong_ordering operator<=>(const GeometryRef& other) const noexcept
    {
        return m_Index <=> other.m_Index;
    }

private:
    std::uint32_t m_Index = Invalid;
};

struct CurveTag;
struct PCurveTag;
struct SurfaceTag;

/// A 3D curve in the curve pool.
using CurveRef = GeometryRef<CurveTag>;
/// A curve in a surface's parameter plane, in the pcurve pool.
using PCurveRef = GeometryRef<PCurveTag>;
/// A surface in the surface pool.
using SurfaceRef = GeometryRef<SurfaceTag>;

} // namespace Stereon

export namespace Stereon::Detail {

/// Liveness and generation of a slot. Kept in the record rather than in a
/// separate bitset so that it shares pages with the record (ADR-0005): a
/// deletion copies one page, not a whole-table bitset.
struct SlotState {
    std::uint32_t Generation = 0;
    bool Alive = false;
};

struct VertexRecord {
    SlotState Slot;
    Point3 Position;
    double Tolerance = 0.0;
};

/// ADR-0004's edge record, following ADR-0009: kind, tolerance, the primary
/// coedge (for kinds that evaluate through a surface), the 3D curve where the
/// kind has one, and the shared parameter range of every pcurve. The radial
/// anchor starts the ring of coedges; it need not be the primary.
struct EdgeRecord {
    SlotState Slot;
    VertexId Start;
    VertexId End;
    CoedgeId RadialAnchor;
    CoedgeId Primary;
    CurveRef Curve;
    Interval Range;
    double Tolerance = 0.0;
    EdgeCurveKind Kind = EdgeCurveKind::Spatial;
};

/// A use of an edge in a loop: the loop cycle (Next, Prev), the radial ring
/// of uses of the same edge (Radial), and the pcurve of this use on the
/// face's surface, on the edge's shared parameter.
struct CoedgeRecord {
    SlotState Slot;
    CoedgeId Next;
    CoedgeId Prev;
    CoedgeId Radial;
    EdgeId Edge;
    LoopId Loop;
    PCurveRef PCurve;
    EdgeUseSense Sense = EdgeUseSense::Same;
};

struct LoopRecord {
    SlotState Slot;
    FaceId Face;
    /// The coedge a loop walk starts from.
    CoedgeId Anchor;
    LoopId NextInFace;
};

struct FaceRecord {
    SlotState Slot;
    ShellId Shell;
    LoopId FirstLoop;
    LoopId LastLoop;
    FaceId NextInShell;
    SurfaceRef Surface;
    /// The face's normal is the surface's reversed.
    bool Reversed = false;
};

struct ShellRecord {
    SlotState Slot;
    SolidId Solid;
    FaceId FirstFace;
    FaceId LastFace;
    ShellId NextInSolid;
};

struct SolidRecord {
    SlotState Slot;
    ShellId FirstShell;
    ShellId LastShell;
};

/// The record type of each entity kind.
template <EntityKind K>
struct RecordFor;
template <>
struct RecordFor<EntityKind::Vertex> {
    using Type = VertexRecord;
};
template <>
struct RecordFor<EntityKind::Edge> {
    using Type = EdgeRecord;
};
template <>
struct RecordFor<EntityKind::Coedge> {
    using Type = CoedgeRecord;
};
template <>
struct RecordFor<EntityKind::Loop> {
    using Type = LoopRecord;
};
template <>
struct RecordFor<EntityKind::Face> {
    using Type = FaceRecord;
};
template <>
struct RecordFor<EntityKind::Shell> {
    using Type = ShellRecord;
};
template <>
struct RecordFor<EntityKind::Solid> {
    using Type = SolidRecord;
};

template <EntityKind K>
using Record = typename RecordFor<K>::Type;

/// Every table of one topology store: the seven entity tables, the live
/// count of each, and the three geometry pools.
struct TopologyTables {
    PagedArena<VertexRecord> Vertices;
    PagedArena<EdgeRecord> Edges;
    PagedArena<CoedgeRecord> Coedges;
    PagedArena<LoopRecord> Loops;
    PagedArena<FaceRecord> Faces;
    PagedArena<ShellRecord> Shells;
    PagedArena<SolidRecord> Solids;
    std::uint32_t LiveCounts[EntityKindCount] = {};

    PagedArena<Curve> Curves;
    PagedArena<Curve2d> PCurves;
    PagedArena<Surface> Surfaces;

    template <EntityKind K>
    [[nodiscard]] const PagedArena<Record<K>>& Table() const noexcept
    {
        if constexpr (K == EntityKind::Vertex) {
            return Vertices;
        } else if constexpr (K == EntityKind::Edge) {
            return Edges;
        } else if constexpr (K == EntityKind::Coedge) {
            return Coedges;
        } else if constexpr (K == EntityKind::Loop) {
            return Loops;
        } else if constexpr (K == EntityKind::Face) {
            return Faces;
        } else if constexpr (K == EntityKind::Shell) {
            return Shells;
        } else {
            return Solids;
        }
    }

    template <EntityKind K>
    [[nodiscard]] PagedArena<Record<K>>& Table() noexcept
    {
        return const_cast<PagedArena<Record<K>>&>(std::as_const(*this).template Table<K>());
    }

    /// True when `index` names a live slot of the table of kind K.
    template <EntityKind K>
    [[nodiscard]] bool IsLive(const std::uint32_t index) const noexcept
    {
        const auto& table = Table<K>();
        return index < table.GetSize() && table[index].Slot.Alive;
    }

    template <EntityKind K>
    [[nodiscard]] std::uint32_t LiveCount() const noexcept
    {
        return LiveCounts[static_cast<std::size_t>(K)];
    }
};

} // namespace Stereon::Detail
