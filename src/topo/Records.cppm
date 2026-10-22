// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string_view>
#include <tuple>
#include <type_traits>
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

/// True for the kinds that keep a cached 3D curve of their primary side
/// (ADR-0009, sections 1 and 5): their 3D curve is S(c(t)), derived from a
/// pcurve and a surface rather than stored as the truth.
[[nodiscard]] constexpr bool HasCachedCurve(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Intersection || kind == EdgeCurveKind::Trim;
}

/// What is known of a cached edge curve's deviation from the curve it
/// approximates (ADR-0009, section 4).
enum class CacheCertificate : std::uint8_t {
    /// The edge keeps no cached curve: none was built, or a change to a
    /// pcurve, the primary coedge or a surface of the edge dropped it.
    None,
    /// A cached curve, built from the current pcurves and surfaces, whose
    /// deviation is not certified: the edge has no tolerance to certify
    /// within, or CertifyDeviation does not cover its types yet. Sampling
    /// the deviation is a regression check, never a certificate.
    Uncertified,
    /// A cached curve with a proven bound on its deviation over the whole
    /// range, from CertifyDeviation (ADR-0009 section 4).
    Certified,
};

[[nodiscard]] constexpr std::string_view ToString(const CacheCertificate certificate) noexcept
{
    switch (certificate) {
    case CacheCertificate::None:
        return "no cache";
    case CacheCertificate::Uncertified:
        return "cached, uncertified";
    case CacheCertificate::Certified:
        return "cached, certified";
    }
    return "CacheCertificate";
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
/// deletion copies one page, not a whole-table bitset. One word, the top bit
/// the liveness and the rest the generation, so that it costs every record
/// four bytes rather than eight with padding (B0.4).
class SlotState {
public:
    /// A live slot of generation 0.
    [[nodiscard]] static constexpr SlotState Live() noexcept { return SlotState(AliveBit); }

    [[nodiscard]] constexpr bool IsAlive() const noexcept { return (m_Word & AliveBit) != 0; }
    /// Bumped by every deletion: a handle made before it no longer matches.
    [[nodiscard]] constexpr std::uint32_t GetGeneration() const noexcept
    {
        return m_Word & ~AliveBit;
    }
    /// Deletes the slot: not alive, and a new generation.
    constexpr void Kill() noexcept { m_Word = (GetGeneration() + 1U) & ~AliveBit; }

    constexpr SlotState() = default;
    constexpr bool operator==(const SlotState&) const = default;

private:
    static constexpr std::uint32_t AliveBit = 0x8000'0000U;

    constexpr explicit SlotState(const std::uint32_t word) noexcept
        : m_Word(word)
    {}

    std::uint32_t m_Word = 0;
};

// The records' fields are ordered so that the four-byte slot and handles
// come first and fill the space before the first eight-byte field, and the
// one-byte fields share the tail: the sizes in a release build are checked
// in tests/topo/RecordLayoutTests.cpp (B0.4).

struct VertexRecord {
    SlotState Slot;
    /// The maintained vertex-incidence index (ADR-0004, ADR-0036): the first
    /// of the edges that start or end here, in ascending index order; each
    /// edge continues the list through NextAtStart or NextAtEnd, whichever
    /// of its ends this vertex is. A closed edge is listed once, through
    /// NextAtStart.
    EdgeId FirstEdge;
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
    /// The next edge in Start's incidence list, and in End's (ADR-0004).
    EdgeId NextAtStart;
    EdgeId NextAtEnd;
    CurveRef Curve;
    /// The cached 3D curve of an Intersection or Trim edge (ADR-0009,
    /// section 5), in the curve pool: an approximation of the primary side
    /// S(c(t)) on the edge's parameter. A cache, never the truth: dropped
    /// whenever a pcurve, the primary coedge or a surface of the edge
    /// changes (ADR-0004), and rebuilt from them deterministically.
    CurveRef Cached;
    EdgeCurveKind Kind = EdgeCurveKind::Spatial;
    /// What is known of the cache's deviation; None exactly when there is no
    /// cache. Plain scalars, as the snapshot's certificate is kept.
    CacheCertificate CacheStatus = CacheCertificate::None;
    Interval Range;
    double Tolerance = 0.0;
    /// The proven bound on the cache's deviation, for a Certified cache;
    /// infinite otherwise.
    double CacheBound = std::numeric_limits<double>::infinity();
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

/// A pool of geometry for a topology store (ADR-0004, B0.4): the common
/// built-in types stored by value, each in its own paged arena with no
/// wrapper around it, and everything else (B-splines, trimmed and swept
/// geometry, user types) as erased wrappers. A dense table of slots maps
/// each index the topology refers to onto its arena and its place there,
/// so that geometry references stay dense indices and every page, slots
/// included, is shared copy-on-write as the entity tables are (ADR-0005).
///
/// A line pcurve takes 32 bytes and its slot 4, against 80 for a Curve2d:
/// a shape keeps one curve per edge and one pcurve per coedge, so the pools
/// are most of its memory. Reading an entry gives the wrapper by value
/// (operator[]), or the stored value itself without one (Visit).
template <typename Wrapper, typename... Builtins>
class GeometryPool {
public:
    static constexpr std::size_t KindCount = sizeof...(Builtins) + 1;

    /// Appends the geometry and returns its index: in its type's arena when
    /// it is one of Builtins, among the wrappers otherwise.
    std::size_t PushBack(Wrapper value)
    {
        const std::uint32_t slot = Place<0>(std::move(value));
        return m_Slots.PushBack(slot);
    }

    /// Replaces the entry at index. The previous value stays in its arena,
    /// unreferenced, until Compact() rebuilds the pool.
    void Set(const std::size_t index, Wrapper value)
    {
        m_Slots.Set(index, Place<0>(std::move(value)));
    }

    /// Calls f with the stored value: the built-in type itself, or the
    /// wrapper. Every call of f must return the same type.
    template <typename F>
    decltype(auto) Visit(const std::size_t index, F&& f) const
    {
        const std::uint32_t slot = m_Slots[index];
        return VisitKind<0>(slot >> LocalBits, slot & LocalMask, f);
    }

    /// True when the entry holds geometry: a built-in value, or a wrapper
    /// that is not empty. After a failed ShapeBuilder::Transaction, later
    /// allocations pad reserved indices with empty wrappers so that no
    /// reference the aborted edit made names geometry added later.
    [[nodiscard]] bool HoldsGeometry(const std::size_t index) const
    {
        return Visit(index, []<typename T>(const T& value) {
            if constexpr (std::is_same_v<T, Wrapper>) {
                return !value.IsEmpty();
            } else {
                return true;
            }
        });
    }
    /// The entry, wrapped: a copy of the built-in value in a wrapper, or of
    /// the stored wrapper.
    [[nodiscard]] Wrapper operator[](const std::size_t index) const
    {
        return Visit(index, [](const auto& g) { return Wrapper(g); });
    }

    [[nodiscard]] std::size_t GetSize() const noexcept { return m_Slots.GetSize(); }

    /// The pages of the slots and of every arena.
    [[nodiscard]] std::size_t GetPageBytes() const noexcept
    {
        return m_Slots.GetPageBytes() + m_Other.GetPageBytes() +
               std::apply(
                   [](const auto&... a) { return (std::size_t{0} + ... + a.GetPageBytes()); },
                   m_Typed);
    }
    [[nodiscard]] std::size_t GetPageCount() const noexcept
    {
        return m_Slots.GetPageCount() + m_Other.GetPageCount() +
               std::apply(
                   [](const auto&... a) { return (std::size_t{0} + ... + a.GetPageCount()); },
                   m_Typed);
    }

    /// Pages two pools share, over the slots and every arena.
    friend std::size_t SharedPageCount(const GeometryPool& a, const GeometryPool& b) noexcept
    {
        return SharedPageCount(a.m_Slots, b.m_Slots) + SharedPageCount(a.m_Other, b.m_Other) +
               SharedTyped(a, b, std::index_sequence_for<Builtins...>{});
    }

private:
    /// A slot: the arena (Builtins in order, then the wrappers) in the top
    /// bits, the place in it below.
    static constexpr std::uint32_t LocalBits = 28;
    static constexpr std::uint32_t LocalMask = (1U << LocalBits) - 1U;
    static_assert(KindCount <= (1U << (32U - LocalBits)));

    template <std::size_t I>
    using BuiltinAt = std::tuple_element_t<I, std::tuple<Builtins...>>;

    [[nodiscard]] static std::uint32_t Encode(const std::size_t kind, const std::size_t local)
    {
        // An arena of 2^28 geometries of one type is past any body the
        // 32-bit entity indices can describe.
        STEREON_ASSERT(local <= LocalMask);
        return (static_cast<std::uint32_t>(kind) << LocalBits) | static_cast<std::uint32_t>(local);
    }

    template <std::size_t I>
    [[nodiscard]] std::uint32_t Place(Wrapper&& value)
    {
        if constexpr (I == sizeof...(Builtins)) {
            return Encode(I, m_Other.PushBack(std::move(value)));
        } else {
            if (const BuiltinAt<I>* builtin = value.template Target<BuiltinAt<I>>()) {
                return Encode(I, std::get<I>(m_Typed).PushBack(*builtin));
            }
            return Place<I + 1>(std::move(value));
        }
    }

    template <std::size_t I, typename F>
    decltype(auto) VisitKind(const std::uint32_t kind, const std::uint32_t local, F& f) const
    {
        if constexpr (I == sizeof...(Builtins)) {
            return std::invoke(f, m_Other[local]);
        } else {
            if (kind == I) {
                return std::invoke(f, std::get<I>(m_Typed)[local]);
            }
            return VisitKind<I + 1>(kind, local, f);
        }
    }

    template <std::size_t... I>
    [[nodiscard]] static std::size_t SharedTyped(const GeometryPool& a, const GeometryPool& b,
                                                 std::index_sequence<I...>) noexcept
    {
        return (std::size_t{0} + ... +
                SharedPageCount(std::get<I>(a.m_Typed), std::get<I>(b.m_Typed)));
    }

    PagedArena<std::uint32_t> m_Slots;
    std::tuple<PagedArena<Builtins>...> m_Typed;
    PagedArena<Wrapper> m_Other;
};

/// True when `ref` names geometry in `pool`: valid, in range and not an
/// empty entry.
template <typename Pool, typename Ref>
[[nodiscard]] bool InPool(const Pool& pool, const Ref ref)
{
    return ref.IsValid() && ref.GetIndex() < pool.GetSize() && pool.HoldsGeometry(ref.GetIndex());
}
/// The pools of a topology store. The built-in types stored by value are
/// the ones bodies are made of in bulk: the conics and analytic surfaces.
using CurvePool = GeometryPool<Curve, Line, Circle, Ellipse>;
using PCurvePool = GeometryPool<Curve2d, Line2d, Circle2d, Ellipse2d>;
using SurfacePool = GeometryPool<Surface, Plane, Cylinder, Cone, Sphere, Torus>;

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

    CurvePool Curves;
    PCurvePool PCurves;
    SurfacePool Surfaces;

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
        return index < table.GetSize() && table[index].Slot.IsAlive();
    }

    template <EntityKind K>
    [[nodiscard]] std::uint32_t LiveCount() const noexcept
    {
        return LiveCounts[static_cast<std::size_t>(K)];
    }
};

} // namespace Stereon::Detail
