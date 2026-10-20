// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Topo:Shape;

import Stereon.Intersect;
import :Records;
import :Report;

/// The frozen side of ADR-0004 and ADR-0038: an immutable TopologySnapshot,
/// certified by the structural validator when Freeze() made it, a Shape that
/// owns one, and the read-only ranges of ADR-0036 over it.
///
/// Every range holds one snapshot ownership token (a shared_ptr to the
/// snapshot), taken once at acquisition; its iterators borrow the range and
/// do no ownership work, allocate nothing and cost O(1) per step. A range
/// outlives the Shape it came from. Iterators must not outlive their range,
/// so the ranges are not borrowed ranges, and an iterator taken from a
/// temporary range dangles: use the temporary in a range-based for, whose
/// lifetime covers the loop, or keep the range.

export namespace Stereon {

class Shape;
class ShapeBuilder;
class Compacted;

[[nodiscard]] Result<Shape> Freeze(ShapeBuilder&& builder);
[[nodiscard]] Result<Compacted> Compact(const Shape& shape);

} // namespace Stereon

export namespace Stereon::Detail {

/// Only Freeze() makes snapshots, so a certificate cannot be forged; the
/// builder and Compact() hold the key to open a builder on a snapshot's
/// tables, which only Freeze() can turn into a snapshot again.
class FreezeKey {
    FreezeKey() = default;
    friend Result<Shape> Stereon::Freeze(ShapeBuilder&& builder);
    friend Result<Compacted> Stereon::Compact(const Shape& shape);
    friend class Stereon::ShapeBuilder;
};

/// Reports a handle that is not a live entity of the shape it is used with,
/// in checked builds, through the report of ADR-0004.
void ShapeHandleViolation(EntityKind kind, std::uint32_t index, std::uint32_t generation,
                          StoreId store, StoreId expectedStore, std::uint32_t size, bool alive);

/// The version recorded in every ValidationCertificate: incremented
/// whenever the structural rules change.
inline constexpr std::uint32_t StructuralValidatorVersion = 3;
/// The rule set the structural validator runs.
inline constexpr std::string_view StructuralInvariants = "Structural S01-S17";

} // namespace Stereon::Detail

export namespace Stereon {

/// What Freeze() proved about a snapshot (ADR-0038): which validator, which
/// invariants, which snapshot. A value the snapshot makes on request from
/// what it stores; never serialised, since a deserialised store is validated
/// again.
struct ValidationCertificate {
    /// Incremented whenever the structural rules change.
    std::uint32_t ValidatorVersion = 0;
    /// The rule set the validator ran.
    std::string_view Invariants;
    /// The snapshot it certifies.
    StoreId Snapshot = 0;
};

/// The immutable topology of one frozen Shape: its tables and its store ID,
/// which checked builds stamp on every handle it gives out. Shared, never
/// copied, by the Shape and every range over it.
///
/// Beside the tables it keeps two integers and builds its certificate from
/// them on request. On the Linux GCC 16.2 image, a certificate built with
/// designated initializers in a module implementation unit and passed by
/// value was read back with its fields scrambled by importing units
/// (docs/toolchain/gcc16-feature-probe.md); plain scalars survive the trip.
class TopologySnapshot {
public:
    TopologySnapshot(Detail::FreezeKey, Detail::TopologyTables tables, const StoreId store,
                     const StoreId sourceBuilder) noexcept
        : m_Tables(std::move(tables))
        , m_Store(store)
        , m_SourceBuilder(sourceBuilder)
    {}

    [[nodiscard]] const Detail::TopologyTables& GetTables() const noexcept { return m_Tables; }
    [[nodiscard]] ValidationCertificate GetCertificate() const noexcept
    {
        return {.ValidatorVersion = Detail::StructuralValidatorVersion,
                .Invariants = Detail::StructuralInvariants,
                .Snapshot = m_Store};
    }
    [[nodiscard]] StoreId GetStore() const noexcept { return m_Store; }
    /// The token of the builder Freeze() consumed: the store its handles
    /// carry in checked builds.
    [[nodiscard]] StoreId GetSourceBuilder() const noexcept { return m_SourceBuilder; }

    /// The handle to slot `index` of the table of kind K, stamped with this
    /// snapshot and the slot's generation in checked builds.
    template <EntityKind K>
    [[nodiscard]] Id<K> MakeId(const std::uint32_t index) const noexcept
    {
        if constexpr (STEREON_CHECKED_HANDLES) {
            return Id<K>::Make(index, GetStore(),
                               m_Tables.template Table<K>()[index].Slot.GetGeneration());
        } else {
            return Id<K>::FromIndex(index);
        }
    }

    /// True when h names a live entity of this snapshot; in checked builds
    /// it must also carry this snapshot's store and the slot's generation.
    template <EntityKind K>
    [[nodiscard]] bool Contains(const Id<K> h) const noexcept
    {
        if (!h.IsValid() || !m_Tables.template IsLive<K>(h.GetIndex())) {
            return false;
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            return h.GetStore() == GetStore() &&
                   h.GetGeneration() ==
                       m_Tables.template Table<K>()[h.GetIndex()].Slot.GetGeneration();
        }
        return true;
    }

private:
    Detail::TopologyTables m_Tables;
    StoreId m_Store = 0;
    StoreId m_SourceBuilder = 0;
};

/// A frozen, certified B-rep (ADR-0004): immutable, cheap to copy, safe to
/// read from any thread. A default-constructed Shape is empty: it has no
/// snapshot, and every query on it reports InvalidInput.
class Shape {
public:
    Shape() = default;
    Shape(Detail::FreezeKey, std::shared_ptr<const TopologySnapshot> snapshot) noexcept
        : m_Snapshot(std::move(snapshot))
    {}

    [[nodiscard]] bool IsEmpty() const noexcept { return m_Snapshot == nullptr; }

    /// The snapshot, which ranges share; null for an empty Shape.
    [[nodiscard]] const std::shared_ptr<const TopologySnapshot>& GetSnapshot() const noexcept
    {
        return m_Snapshot;
    }

    /// Live entities of kind K; 0 for an empty Shape.
    template <EntityKind K>
    [[nodiscard]] std::uint32_t Count() const noexcept
    {
        return m_Snapshot == nullptr ? 0 : m_Snapshot->GetTables().template LiveCount<K>();
    }

    /// The frozen handle of an entity the builder created: Freeze() keeps
    /// every index, so the slot is the same. Invalid when the builder's
    /// handle does not name a live entity of this shape, and in checked
    /// builds when it comes from another builder.
    template <EntityKind K>
    [[nodiscard]] Id<K> Resolve(const BuilderId<K> h) const noexcept
    {
        if (m_Snapshot == nullptr || !h.IsValid() ||
            !m_Snapshot->GetTables().template IsLive<K>(h.GetIndex())) {
            return {};
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            if (h.GetStore() != m_Snapshot->GetSourceBuilder()) {
                return {};
            }
        }
        const Detail::SlotState& slot =
            m_Snapshot->GetTables().template Table<K>()[h.GetIndex()].Slot;
        if (slot.GetGeneration() != h.GetGeneration()) {
            return {};
        }
        return m_Snapshot->template MakeId<K>(h.GetIndex());
    }

private:
    std::shared_ptr<const TopologySnapshot> m_Snapshot;
};

} // namespace Stereon

export namespace Stereon::Detail {

/// The record of a live entity of a non-empty shape. Handing an accessor a
/// handle that does not name one is a programming error (ADR-0006): checked
/// builds report it, release builds read whatever the slot holds.
template <EntityKind K>
[[nodiscard]] const Record<K>& RecordOf(const Shape& shape, const Id<K> h)
{
    const TopologySnapshot& snapshot = *shape.GetSnapshot();
    if constexpr (STEREON_CHECKED_HANDLES) {
        if (!snapshot.Contains(h)) {
            const auto& table = snapshot.GetTables().template Table<K>();
            const bool inRange = h.IsValid() && h.GetIndex() < table.GetSize();
            ShapeHandleViolation(K, h.GetIndex(), h.GetGeneration(), h.GetStore(),
                                 snapshot.GetStore(), static_cast<std::uint32_t>(table.GetSize()),
                                 inRange && table[h.GetIndex()].Slot.IsAlive());
        }
    }
    return snapshot.GetTables().template Table<K>()[h.GetIndex()];
}

/// The checks every range acquisition makes (ADR-0038): a certified
/// snapshot, and a handle that names a live entity of it.
template <EntityKind K>
[[nodiscard]] Result<void> CheckQuery(const Shape& shape, const Id<K> h,
                                      const std::string_view query)
{
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, std::string(query) + ": the Shape is empty");
    }
    if (!shape.GetSnapshot()->Contains(h)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(query) + ": the handle is not a live entity of this Shape");
    }
    return {};
}

} // namespace Stereon::Detail

export namespace Stereon {

/// Every live entity of kind K, once each, in ascending slot order
/// (ADR-0036, ADR-0037). Over a table of S slots it costs O(S): the scan
/// skips the dead slots, the first of them at construction.
template <EntityKind K>
class EntityRange {
public:
    class Iterator {
    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::input_iterator_tag;
        using value_type = Id<K>;
        using difference_type = std::ptrdiff_t;

        Iterator() = default;

        [[nodiscard]] Id<K> operator*() const { return m_Snapshot->template MakeId<K>(m_Index); }
        Iterator& operator++()
        {
            m_Index = NextLive(*m_Snapshot, m_Index + 1);
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator copy = *this;
            ++*this;
            return copy;
        }
        bool operator==(const Iterator& other) const noexcept { return m_Index == other.m_Index; }

    private:
        friend class EntityRange;
        Iterator(const TopologySnapshot* snapshot, const std::uint32_t index) noexcept
            : m_Snapshot(snapshot)
            , m_Index(index)
        {}

        const TopologySnapshot* m_Snapshot = nullptr;
        std::uint32_t m_Index = 0;
    };

    explicit EntityRange(std::shared_ptr<const TopologySnapshot> snapshot)
        : m_Snapshot(std::move(snapshot))
        , m_First(NextLive(*m_Snapshot, 0))
    {}

    [[nodiscard]] Iterator begin() const noexcept { return {m_Snapshot.get(), m_First}; }
    [[nodiscard]] Iterator end() const noexcept { return {m_Snapshot.get(), SizeOf(*m_Snapshot)}; }
    /// The live count, kept by the table: O(1).
    [[nodiscard]] std::size_t size() const noexcept
    {
        return m_Snapshot->GetTables().template LiveCount<K>();
    }
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }

private:
    [[nodiscard]] static std::uint32_t SizeOf(const TopologySnapshot& snapshot) noexcept
    {
        return static_cast<std::uint32_t>(snapshot.GetTables().template Table<K>().GetSize());
    }
    [[nodiscard]] static std::uint32_t NextLive(const TopologySnapshot& snapshot,
                                                std::uint32_t index) noexcept
    {
        const auto& table = snapshot.GetTables().template Table<K>();
        const std::uint32_t size = SizeOf(snapshot);
        while (index < size && !table[index].Slot.IsAlive()) {
            ++index;
        }
        return index;
    }

    std::shared_ptr<const TopologySnapshot> m_Snapshot;
    std::uint32_t m_First = 0;
};

namespace Detail {

/// A walk along a singly linked list or ring of entities of kind K: Step
/// gives the next index from the current record, and the walk ends at an
/// invalid index or, for a ring, back at the anchor. Each step reads one
/// record: O(1), no allocation, no ownership work.
template <EntityKind K, typename Step, bool Ring>
class LinkRange {
public:
    class Iterator {
    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::input_iterator_tag;
        using value_type = Id<K>;
        using difference_type = std::ptrdiff_t;

        Iterator() = default;

        [[nodiscard]] Id<K> operator*() const { return m_Snapshot->template MakeId<K>(m_Current); }
        Iterator& operator++()
        {
            const auto& record = m_Snapshot->GetTables().template Table<K>()[m_Current];
            const std::uint32_t next = Step{}(record);
            m_Current =
                (next == Id<K>::Invalid || (Ring && next == m_Anchor)) ? Id<K>::Invalid : next;
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator copy = *this;
            ++*this;
            return copy;
        }
        bool operator==(const Iterator& other) const noexcept
        {
            return m_Current == other.m_Current;
        }

    private:
        friend class LinkRange;
        Iterator(const TopologySnapshot* snapshot, const std::uint32_t current,
                 const std::uint32_t anchor) noexcept
            : m_Snapshot(snapshot)
            , m_Current(current)
            , m_Anchor(anchor)
        {}

        const TopologySnapshot* m_Snapshot = nullptr;
        std::uint32_t m_Current = Id<K>::Invalid;
        std::uint32_t m_Anchor = Id<K>::Invalid;
    };

    LinkRange(std::shared_ptr<const TopologySnapshot> snapshot, const std::uint32_t first)
        : m_Snapshot(std::move(snapshot))
        , m_First(first)
    {}

    [[nodiscard]] Iterator begin() const noexcept { return {m_Snapshot.get(), m_First, m_First}; }
    [[nodiscard]] Iterator end() const noexcept
    {
        return {m_Snapshot.get(), Id<K>::Invalid, m_First};
    }
    [[nodiscard]] bool empty() const noexcept { return m_First == Id<K>::Invalid; }

private:
    std::shared_ptr<const TopologySnapshot> m_Snapshot;
    std::uint32_t m_First = Id<K>::Invalid;
};

struct NextInFace {
    std::uint32_t operator()(const LoopRecord& r) const noexcept { return r.NextInFace.GetIndex(); }
};
struct NextInShell {
    std::uint32_t operator()(const FaceRecord& r) const noexcept
    {
        return r.NextInShell.GetIndex();
    }
};
struct NextInSolid {
    std::uint32_t operator()(const ShellRecord& r) const noexcept
    {
        return r.NextInSolid.GetIndex();
    }
};
struct NextInLoop {
    std::uint32_t operator()(const CoedgeRecord& r) const noexcept { return r.Next.GetIndex(); }
};
struct NextInRing {
    std::uint32_t operator()(const CoedgeRecord& r) const noexcept { return r.Radial.GetIndex(); }
};

} // namespace Detail

/// The loops of a face, in the order they were added; the first is the one
/// added first.
using LoopRange = Detail::LinkRange<EntityKind::Loop, Detail::NextInFace, false>;
/// The faces of a shell, in the order they were added.
using ShellFaceRange = Detail::LinkRange<EntityKind::Face, Detail::NextInShell, false>;
/// The shells of a solid, in the order they were added.
using SolidShellRange = Detail::LinkRange<EntityKind::Shell, Detail::NextInSolid, false>;
/// One walk around a loop along Next, from its anchor, each coedge once
/// (ADR-0037): k coedges cost O(k).
using CoedgeRange = Detail::LinkRange<EntityKind::Coedge, Detail::NextInLoop, true>;
/// Every use of an edge, once each, along its radial ring from the ring's
/// anchor (ADR-0037). The ring order is a connectivity order, not an angular
/// one, and the anchor need not be the primary coedge.
using RadialCoedgeRange = Detail::LinkRange<EntityKind::Coedge, Detail::NextInRing, true>;

/// The edges at a vertex, each once (ADR-0036: distinct, a closed edge
/// once), in ascending index order, read from the maintained
/// vertex-incidence index: each step reads one edge record, O(1), no
/// allocation, so a vertex of degree d costs O(d).
class VertexEdgeRange {
public:
    class Iterator {
    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::input_iterator_tag;
        using value_type = EdgeId;
        using difference_type = std::ptrdiff_t;

        Iterator() = default;

        [[nodiscard]] EdgeId operator*() const
        {
            return m_Snapshot->template MakeId<EntityKind::Edge>(m_Current);
        }
        Iterator& operator++()
        {
            const Detail::EdgeRecord& e = m_Snapshot->GetTables().Edges[m_Current];
            m_Current = (e.Start.GetIndex() == m_Vertex ? e.NextAtStart : e.NextAtEnd).GetIndex();
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator copy = *this;
            ++*this;
            return copy;
        }
        bool operator==(const Iterator& other) const noexcept
        {
            return m_Current == other.m_Current;
        }

    private:
        friend class VertexEdgeRange;
        Iterator(const TopologySnapshot* snapshot, const std::uint32_t vertex,
                 const std::uint32_t current) noexcept
            : m_Snapshot(snapshot)
            , m_Vertex(vertex)
            , m_Current(current)
        {}

        const TopologySnapshot* m_Snapshot = nullptr;
        std::uint32_t m_Vertex = VertexId::Invalid;
        std::uint32_t m_Current = EdgeId::Invalid;
    };

    VertexEdgeRange(std::shared_ptr<const TopologySnapshot> snapshot, const std::uint32_t vertex,
                    const std::uint32_t first)
        : m_Snapshot(std::move(snapshot))
        , m_Vertex(vertex)
        , m_First(first)
    {}

    [[nodiscard]] Iterator begin() const noexcept { return {m_Snapshot.get(), m_Vertex, m_First}; }
    [[nodiscard]] Iterator end() const noexcept
    {
        return {m_Snapshot.get(), m_Vertex, EdgeId::Invalid};
    }
    [[nodiscard]] bool empty() const noexcept { return m_First == EdgeId::Invalid; }

private:
    std::shared_ptr<const TopologySnapshot> m_Snapshot;
    std::uint32_t m_Vertex = VertexId::Invalid;
    std::uint32_t m_First = EdgeId::Invalid;
};

/// Every live entity of kind K in the shape: not the descendants of a root,
/// all of them (ADR-0036).
template <EntityKind K>
[[nodiscard]] Result<EntityRange<K>> AllEntities(const Shape& shape)
{
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AllEntities: the Shape is empty");
    }
    return EntityRange<K>(shape.GetSnapshot());
}

[[nodiscard]] Result<LoopRange> LoopsOf(const Shape& shape, FaceId face);
[[nodiscard]] Result<ShellFaceRange> FacesOf(const Shape& shape, ShellId shell);
[[nodiscard]] Result<SolidShellRange> ShellsOf(const Shape& shape, SolidId solid);
[[nodiscard]] Result<CoedgeRange> CoedgesOf(const Shape& shape, LoopId loop);
[[nodiscard]] Result<RadialCoedgeRange> RadialCoedgesOf(const Shape& shape, EdgeId edge);
/// The edges at a vertex, from the maintained index (ADR-0036).
[[nodiscard]] Result<VertexEdgeRange> EdgesOf(const Shape& shape, VertexId vertex);
/// The distinct faces of an edge's uses, once each, in ascending index order
/// (ADR-0036, ADR-0037): an eager result over the radial ring, k uses
/// costing O(k log k) and one allocation. A seam's face comes once; an edge
/// with no uses gives an empty vector.
[[nodiscard]] Result<std::vector<FaceId>> FacesOf(const Shape& shape, EdgeId edge);

// --- Accessors: O(1) reads of one record --------------------------------
//
// The handle must name a live entity of the shape (a programming error
// otherwise, reported in checked builds).

[[nodiscard]] inline Point3 PositionOf(const Shape& shape, const VertexId v)
{
    return Detail::RecordOf(shape, v).Position;
}
[[nodiscard]] inline double ToleranceOf(const Shape& shape, const VertexId v)
{
    return Detail::RecordOf(shape, v).Tolerance;
}

[[nodiscard]] inline EdgeCurveKind KindOf(const Shape& shape, const EdgeId e)
{
    return Detail::RecordOf(shape, e).Kind;
}
[[nodiscard]] inline double ToleranceOf(const Shape& shape, const EdgeId e)
{
    return Detail::RecordOf(shape, e).Tolerance;
}
/// The shared parameter range of the edge's curve and pcurves (ADR-0009).
[[nodiscard]] inline Interval RangeOf(const Shape& shape, const EdgeId e)
{
    return Detail::RecordOf(shape, e).Range;
}
[[nodiscard]] inline VertexId StartOf(const Shape& shape, const EdgeId e)
{
    return shape.GetSnapshot()->MakeId<EntityKind::Vertex>(
        Detail::RecordOf(shape, e).Start.GetIndex());
}
[[nodiscard]] inline VertexId EndOf(const Shape& shape, const EdgeId e)
{
    return shape.GetSnapshot()->MakeId<EntityKind::Vertex>(
        Detail::RecordOf(shape, e).End.GetIndex());
}
/// The coedge an edge of a surface-defined kind evaluates through
/// (ADR-0009); invalid for Exact and Spatial edges, which have their own
/// curve.
[[nodiscard]] inline CoedgeId PrimaryCoedgeOf(const Shape& shape, const EdgeId e)
{
    const CoedgeId primary = Detail::RecordOf(shape, e).Primary;
    return primary.IsValid() ? shape.GetSnapshot()->MakeId<EntityKind::Coedge>(primary.GetIndex())
                             : CoedgeId{};
}
// The geometry accessors return the geometry by value, wrapped: the pools
// keep the built-in types without a wrapper (GeometryPool), so there is no
// stored Curve to point to. A wrapper of a built-in type is a copy of a few
// dozen bytes; one of a B-spline shares its payload. Defined in Shape.cpp,
// as they instantiate the geometry modules' templates.

/// The edge's own 3D curve, or nothing for a kind without one.
[[nodiscard]] std::optional<Curve> CurveOf(const Shape& shape, EdgeId e);
/// The cached 3D curve of an Intersection or Trim edge (ADR-0009, section
/// 5): S(c(t)) of its primary side, approximated on the edge's parameter; a
/// cache, never the truth. Nothing when the edge keeps none.
[[nodiscard]] std::optional<Curve> CachedCurveOf(const Shape& shape, EdgeId e);
/// What is known of the cached curve's deviation (CacheCertificate).
[[nodiscard]] inline CacheCertificate CacheCertificateOf(const Shape& shape, const EdgeId e)
{
    return Detail::RecordOf(shape, e).CacheStatus;
}
/// The proven bound on the cached curve's deviation from S(c(t)) over the
/// whole range, for a Certified cache; infinite otherwise.
[[nodiscard]] inline double CacheBoundOf(const Shape& shape, const EdgeId e)
{
    return Detail::RecordOf(shape, e).CacheBound;
}

[[nodiscard]] inline EdgeId EdgeOf(const Shape& shape, const CoedgeId c)
{
    return shape.GetSnapshot()->MakeId<EntityKind::Edge>(
        Detail::RecordOf(shape, c).Edge.GetIndex());
}
[[nodiscard]] inline LoopId LoopOf(const Shape& shape, const CoedgeId c)
{
    return shape.GetSnapshot()->MakeId<EntityKind::Loop>(
        Detail::RecordOf(shape, c).Loop.GetIndex());
}
/// Whether the coedge runs along its edge's parameter or against it: the
/// stored edge-use sense of ADR-0037.
[[nodiscard]] inline EdgeUseSense StoredEdgeUseSense(const Shape& shape, const CoedgeId c)
{
    return Detail::RecordOf(shape, c).Sense;
}
/// The coedge's curve on its face's surface, on the edge's parameter;
/// nothing where the edge kind has none (a Spatial edge's optional
/// projection).
[[nodiscard]] std::optional<Curve2d> PCurveOf(const Shape& shape, CoedgeId c);

[[nodiscard]] inline FaceId FaceOf(const Shape& shape, const LoopId l)
{
    return shape.GetSnapshot()->MakeId<EntityKind::Face>(
        Detail::RecordOf(shape, l).Face.GetIndex());
}
/// The face's surface.
[[nodiscard]] Surface SurfaceOf(const Shape& shape, FaceId f);
/// True when the face's normal is its surface's reversed.
[[nodiscard]] inline bool IsReversed(const Shape& shape, const FaceId f)
{
    return Detail::RecordOf(shape, f).Reversed;
}
/// The shell the face belongs to; invalid for a face in no shell.
[[nodiscard]] inline ShellId ShellOf(const Shape& shape, const FaceId f)
{
    const ShellId shell = Detail::RecordOf(shape, f).Shell;
    return shell.IsValid() ? shape.GetSnapshot()->MakeId<EntityKind::Shell>(shell.GetIndex())
                           : ShellId{};
}
/// The solid the shell bounds; invalid for a shell in no solid.
[[nodiscard]] inline SolidId SolidOf(const Shape& shape, const ShellId s)
{
    const SolidId solid = Detail::RecordOf(shape, s).Solid;
    return solid.IsValid() ? shape.GetSnapshot()->MakeId<EntityKind::Solid>(solid.GetIndex())
                           : SolidId{};
}

/// The entity for logs, contract messages and checker reports (ADR-0004),
/// e.g. "Edge#4 [Exact Line, tol 1e-05, faces 0|2]"; "Edge#4 [not in this
/// shape]" for a handle that names no live entity of it.
[[nodiscard]] std::string Describe(const Shape& shape, VertexId h);
[[nodiscard]] std::string Describe(const Shape& shape, EdgeId h);
[[nodiscard]] std::string Describe(const Shape& shape, CoedgeId h);
[[nodiscard]] std::string Describe(const Shape& shape, LoopId h);
[[nodiscard]] std::string Describe(const Shape& shape, FaceId h);
[[nodiscard]] std::string Describe(const Shape& shape, ShellId h);
[[nodiscard]] std::string Describe(const Shape& shape, SolidId h);

} // namespace Stereon
