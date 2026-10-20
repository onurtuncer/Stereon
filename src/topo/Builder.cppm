// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

export module Stereon.Topo:Builder;

import Stereon.Intersect;
import :Records;
import :Report;
import :Shape;

/// ShapeBuilder: the mutable, single-threaded side of ADR-0004, and the
/// controlled construction path of ADR-0027. Its tables are private; every
/// connectivity change goes through an operation below that keeps loop
/// links, radial rings and parent lists consistent together, so no caller
/// writes a link field. Freeze() then runs the independent structural
/// validator and, only if it finds nothing, seals the tables into a
/// certified, immutable Shape without renumbering anything.
///
/// Operations check their handles: a handle that is not a live entity of
/// this builder, or a request the topology cannot take (a coedge of an
/// edge in a loop of another builder, a face added to a second shell), is
/// reported as InvalidInput and changes nothing. Every successful change
/// advances the builder's topology epoch (ADR-0038).
///
/// A builder opened on a frozen Shape (Open) starts from the shape's
/// tables and shares every page with it until it writes to one. Deletion
/// (the Remove operations) leaves tombstones: a removed slot is dead, its
/// generation bumped so that every handle to it is stale, and never reused
/// by this builder, so an index live in both the shape and its edit names
/// the same entity. Compact() is the separate, explicit renumbering.

export namespace Stereon {

using BuilderVertexId = BuilderId<EntityKind::Vertex>;
using BuilderEdgeId = BuilderId<EntityKind::Edge>;
using BuilderCoedgeId = BuilderId<EntityKind::Coedge>;
using BuilderLoopId = BuilderId<EntityKind::Loop>;
using BuilderFaceId = BuilderId<EntityKind::Face>;
using BuilderShellId = BuilderId<EntityKind::Shell>;
using BuilderSolidId = BuilderId<EntityKind::Solid>;

/// A new edge. Exact and Spatial edges carry their own 3D curve; every
/// kind has the shared parameter range [t0, t1] of its curve and pcurves
/// (ADR-0009). Start and End are the vertices at t0 and t1; a closed edge
/// has the same vertex at both.
struct EdgeSpec {
    EdgeCurveKind Kind = EdgeCurveKind::Spatial;
    BuilderVertexId Start{};
    BuilderVertexId End{};
    CurveRef Curve{};
    Interval Range{};
    double Tolerance = 0.0;
};

/// The stable identity of one builder (ADR-0038): it moves with the
/// builder, and a moved-from builder has none.
class BuilderToken {
public:
    constexpr BuilderToken() = default;
    constexpr explicit BuilderToken(const StoreId value) noexcept
        : m_Value(value)
    {}
    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Value != 0; }
    [[nodiscard]] constexpr StoreId GetValue() const noexcept { return m_Value; }
    constexpr bool operator==(const BuilderToken&) const = default;

private:
    StoreId m_Value = 0;
};

namespace Testing {
class TopologyFaults;
} // namespace Testing

class Compacted;
[[nodiscard]] Result<Compacted> Compact(const Shape& shape);

class ShapeBuilder {
public:
    ShapeBuilder();
    /// A builder on the shape's tables, sharing every page with it until it
    /// writes to one (ADR-0005). Its handles to the shape's entities come
    /// from Handle(); new entities take slots after the shape's last one.
    /// InvalidInput for an empty Shape.
    [[nodiscard]] static Result<ShapeBuilder> Open(const Shape& shape);
    ShapeBuilder(const ShapeBuilder&) = delete;
    ShapeBuilder& operator=(const ShapeBuilder&) = delete;
    /// The token, the epoch and the tables move; the source is left with no
    /// token, and every operation on it fails.
    ShapeBuilder(ShapeBuilder&& other) noexcept;
    ShapeBuilder& operator=(ShapeBuilder&& other) noexcept;
    ~ShapeBuilder();

    // --- Geometry pools ----------------------------------------------------
    // Values are shared by reference: several edges may use one curve, and
    // several faces one surface, each with its own bounds (ADR-0025).

    [[nodiscard]] CurveRef AddCurve(Curve curve);
    [[nodiscard]] PCurveRef AddPCurve(Curve2d pcurve);
    [[nodiscard]] SurfaceRef AddSurface(Surface surface);

    // --- Entities ----------------------------------------------------------

    [[nodiscard]] Result<BuilderVertexId> AddVertex(const Point3& position, double tolerance);

    /// An edge with no uses yet; coedges join its radial ring as loops take
    /// it. Its curve, where the kind has one, must already be in the pool.
    [[nodiscard]] Result<BuilderEdgeId> AddEdge(const EdgeSpec& spec);

    /// A face on a surface of the pool, in no shell and with no loops.
    [[nodiscard]] Result<BuilderFaceId> AddFace(SurfaceRef surface, bool reversed = false);

    /// An empty loop at the end of the face's loop list. A loop must have a
    /// coedge before Freeze().
    [[nodiscard]] Result<BuilderLoopId> AddLoop(BuilderFaceId face);

    /// A use of `edge` at the end of `loop`'s cycle, running along the edge
    /// or against it, with its pcurve on the face's surface. The coedge
    /// joins the edge's radial ring; the first use of an edge of a kind that
    /// evaluates through a surface (Intersection, Trim, Boundary,
    /// Degenerate) becomes the edge's primary coedge (ADR-0009), which
    /// SetPrimaryCoedge() can change explicitly.
    [[nodiscard]] Result<BuilderCoedgeId> AppendCoedge(BuilderLoopId loop, BuilderEdgeId edge,
                                                       EdgeUseSense sense, PCurveRef pcurve);

    /// Makes `coedge`, a use of `edge`, the edge's primary coedge.
    [[nodiscard]] Result<void> SetPrimaryCoedge(BuilderEdgeId edge, BuilderCoedgeId coedge);

    [[nodiscard]] Result<BuilderSolidId> AddSolid();
    /// A shell bounding `solid`, at the end of its shell list.
    [[nodiscard]] Result<BuilderShellId> AddShell(BuilderSolidId solid);
    /// A shell of no solid: a sheet or an open shell.
    [[nodiscard]] Result<BuilderShellId> AddShell();
    /// Adds a face of no shell at the end of the shell's face list.
    [[nodiscard]] Result<void> AddFaceToShell(BuilderShellId shell, BuilderFaceId face);

    // --- Editing -----------------------------------------------------------
    // Each operation unlinks what it removes from every list and ring it is
    // in and tombstones the slot; nothing else changes.

    /// Replaces the face's surface by another of the pool. The pcurves of its
    /// coedges are the caller's to keep consistent.
    [[nodiscard]] Result<void> SetSurface(BuilderFaceId face, SurfaceRef surface);

    /// Removes one use of an edge from its loop and its edge's radial ring.
    /// Refused when the coedge is its edge's primary and the edge has other
    /// uses: the caller reassigns the primary first (ADR-0027, section 4).
    [[nodiscard]] Result<void> RemoveCoedge(BuilderCoedgeId coedge);
    /// Removes a loop with its coedges from its face. A primary coedge among
    /// them passes the primary to the next use in the ring, since every use
    /// of this loop goes; an edge left with no use loses its primary.
    [[nodiscard]] Result<void> RemoveLoop(BuilderLoopId loop);
    /// Removes a face with its loops and coedges from its shell, as
    /// RemoveLoop does for each loop. Its edges and vertices stay.
    [[nodiscard]] Result<void> RemoveFace(BuilderFaceId face);
    /// Removes a shell from its solid; its faces stay, in no shell.
    [[nodiscard]] Result<void> RemoveShell(BuilderShellId shell);
    /// Removes a solid; its shells stay, of no solid.
    [[nodiscard]] Result<void> RemoveSolid(BuilderSolidId solid);
    /// Removes an edge with no uses from its vertices' lists; refused while
    /// a coedge uses it.
    [[nodiscard]] Result<void> RemoveEdge(BuilderEdgeId edge);
    /// Removes a vertex no edge starts or ends at; refused otherwise.
    [[nodiscard]] Result<void> RemoveVertex(BuilderVertexId vertex);

    // --- Inspection --------------------------------------------------------

    /// The builder's handle to an entity of the shape it was opened on: the
    /// same slot, at its current generation. Invalid when the handle names
    /// no live entity of this builder, in checked builds also when it comes
    /// from another shape.
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> Handle(const Id<K> h) const noexcept
    {
        if (!h.IsValid() || !m_Tables.template IsLive<K>(h.GetIndex())) {
            return {};
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            if (h.GetStore() != m_SourceStore) {
                return {};
            }
        }
        return BuilderId<K>::Make(h.GetIndex(),
                                  m_Tables.template Table<K>()[h.GetIndex()].Slot.Generation,
                                  m_Token.GetValue());
    }

    /// The structural check Freeze() will run, without freezing.
    [[nodiscard]] ValidityReport Validate() const;

    /// Live entities of kind K.
    template <EntityKind K>
    [[nodiscard]] std::uint32_t Count() const noexcept
    {
        return m_Tables.template LiveCount<K>();
    }

    [[nodiscard]] BuilderToken GetToken() const noexcept { return m_Token; }
    /// Advanced by every change to connectivity, membership or liveness.
    [[nodiscard]] std::uint64_t GetEpoch() const noexcept { return m_Epoch; }

    /// A builder on given tables, for Open() and Compact(): `source` is the
    /// store the tables came from, which Handle() checks.
    ShapeBuilder(Detail::TopologyTables tables, StoreId source, Detail::FreezeKey);

private:
    friend Result<Shape> Freeze(ShapeBuilder&& builder);
    friend class Testing::TopologyFaults;

    template <EntityKind K>
    [[nodiscard]] bool IsLive(BuilderId<K> h) const noexcept;
    template <EntityKind K>
    [[nodiscard]] Result<void> Require(BuilderId<K> h, const char* operation) const;
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> Allocate(Detail::Record<K> record);
    template <EntityKind K>
    void Release(std::uint32_t index);
    [[nodiscard]] Result<void> RequireToken(const char* operation) const;
    void Advance() noexcept;

    void LinkAtVertex(std::uint32_t edge, std::uint32_t vertex, bool atStart);
    void UnlinkAtVertex(std::uint32_t edge, std::uint32_t vertex);
    void DetachCoedge(std::uint32_t coedge, bool passPrimaryOn);
    void RemoveLoopAt(std::uint32_t loop);
    template <EntityKind P, EntityKind C>
    void UnlinkChild(std::uint32_t child, Id<P> Detail::Record<C>::* parent,
                     Id<C> Detail::Record<C>::* next, Id<C> Detail::Record<P>::* first,
                     Id<C> Detail::Record<P>::* last);

    Detail::TopologyTables m_Tables;
    BuilderToken m_Token;
    std::uint64_t m_Epoch = 0;
    StoreId m_SourceStore = 0;
};

/// The old-to-new slot map of one Compact() (ADR-0004), which feeds the
/// history of ADR-0011: for every kind, the new index of each old slot,
/// Invalid for a dead one, and the same for the three geometry pools, whose
/// unreferenced entries are dropped.
class RemapTable {
public:
    template <EntityKind K>
    [[nodiscard]] std::uint32_t NewIndex(const std::uint32_t old) const noexcept
    {
        const std::vector<std::uint32_t>& map = m_Entities[static_cast<std::size_t>(K)];
        return old < map.size() ? map[old] : Id<K>::Invalid;
    }
    [[nodiscard]] std::uint32_t NewCurve(const std::uint32_t old) const noexcept
    {
        return old < m_Curves.size() ? m_Curves[old] : CurveRef::Invalid;
    }
    [[nodiscard]] std::uint32_t NewPCurve(const std::uint32_t old) const noexcept
    {
        return old < m_PCurves.size() ? m_PCurves[old] : PCurveRef::Invalid;
    }
    [[nodiscard]] std::uint32_t NewSurface(const std::uint32_t old) const noexcept
    {
        return old < m_Surfaces.size() ? m_Surfaces[old] : SurfaceRef::Invalid;
    }
    /// Slots of kind K in the source, live or dead.
    template <EntityKind K>
    [[nodiscard]] std::size_t OldCount() const noexcept
    {
        return m_Entities[static_cast<std::size_t>(K)].size();
    }

private:
    friend Result<Compacted> Compact(const Shape& shape);

    std::array<std::vector<std::uint32_t>, EntityKindCount> m_Entities;
    std::vector<std::uint32_t> m_Curves;
    std::vector<std::uint32_t> m_PCurves;
    std::vector<std::uint32_t> m_Surfaces;
};

/// What Compact() gives: the compacted shape, certified like any other,
/// and the map from the source's slots to its own.
class Compacted {
public:
    Stereon::Shape Shape;
    RemapTable Remap;
    /// The store of the shape that was compacted, which Translate() checks
    /// in checked builds.
    StoreId SourceStore = 0;

    /// The compacted shape's handle to an entity of the source; invalid for a
    /// handle that named no live entity of the source.
    template <EntityKind K>
    [[nodiscard]] Id<K> Translate(const Id<K> old) const noexcept
    {
        if (!old.IsValid() || Shape.IsEmpty()) {
            return {};
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            if (old.GetStore() != SourceStore) {
                return {};
            }
        }
        const std::uint32_t index = Remap.template NewIndex<K>(old.GetIndex());
        return index == Id<K>::Invalid ? Id<K>{} : Shape.GetSnapshot()->template MakeId<K>(index);
    }
};

/// The explicit renumbering of ADR-0004: the shape's live entities with
/// consecutive indices in their old order, dead slots and unreferenced
/// geometry dropped, frozen and certified again, with the remap table. The
/// source shape is left as it is. InvalidInput for an empty Shape.
[[nodiscard]] Result<Compacted> Compact(const Shape& shape);

/// Seals the builder into a Shape: runs the structural validator over every
/// table (ADR-0038's Freeze invariants) and, if it finds nothing, moves the
/// tables into a certified snapshot, keeping every index. Otherwise reports
/// InvalidInput with the findings. Either way the builder is consumed.
[[nodiscard]] Result<Shape> Freeze(ShapeBuilder&& builder);

/// The structural check of a frozen shape, run again by the independent
/// validator: what Freeze() certified, re-derived (ADR-0036).
[[nodiscard]] ValidityReport ValidateTopology(const Shape& shape);

} // namespace Stereon

export namespace Stereon::Detail {

/// The independent structural validator (ADR-0018, ADR-0038): its own
/// bounded walks over the raw tables, checking every link before following
/// it, so it terminates and reads nothing out of range on any input.
[[nodiscard]] ValidityReport ValidateStructure(const TopologyTables& tables);

/// Describe() over raw tables, safe on corrupt ones: it checks every index
/// and bounds every walk, for the validator's findings.
[[nodiscard]] std::string DescribeEntity(const TopologyTables& tables, EntityKind kind,
                                         std::uint32_t index);

} // namespace Stereon::Detail

export namespace Stereon::Testing {

/// Seeded faults for the validator's tests (ADR-0027, section 8): direct
/// writes to link fields that the builder's operations never make. Not for
/// production code, which cannot otherwise reach the tables.
class TopologyFaults {
public:
    static void SetNext(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId next);
    static void SetPrev(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId prev);
    static void SetRadial(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId radial);
    static void SetCoedgeLoop(ShapeBuilder& b, BuilderCoedgeId c, BuilderLoopId loop);
    static void SetPrimary(ShapeBuilder& b, BuilderEdgeId e, BuilderCoedgeId primary);
    static void SetKind(ShapeBuilder& b, BuilderEdgeId e, EdgeCurveKind kind);
    static void SetEnd(ShapeBuilder& b, BuilderEdgeId e, BuilderVertexId end);
    static void SetPCurve(ShapeBuilder& b, BuilderCoedgeId c, PCurveRef pcurve);
    static void SetFaceShell(ShapeBuilder& b, BuilderFaceId f, BuilderShellId shell);
    static void SetLoopFace(ShapeBuilder& b, BuilderLoopId l, BuilderFaceId face);
    static void SetShellSolid(ShapeBuilder& b, BuilderShellId s, BuilderSolidId solid);
    static void SetFirstEdge(ShapeBuilder& b, BuilderVertexId v, BuilderEdgeId first);
    static void SetNextAtStart(ShapeBuilder& b, BuilderEdgeId e, BuilderEdgeId next);
    static void SetNextAtEnd(ShapeBuilder& b, BuilderEdgeId e, BuilderEdgeId next);
    /// Marks the slot dead without unlinking it.
    static void KillCoedge(ShapeBuilder& b, BuilderCoedgeId c);
    /// Points a link at a slot past the end of its table.
    static void SetNextOutOfRange(ShapeBuilder& b, BuilderCoedgeId c);
};

} // namespace Stereon::Testing
