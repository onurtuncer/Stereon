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

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

module Stereon.Topo;

/// Compact() of ADR-0004: the explicit renumbering that Freeze() never does.
/// It walks every table of a certified snapshot once, gives the live slots
/// of each kind consecutive indices in their old order, rewrites every link
/// through those maps, drops the geometry nothing references and the dead
/// slots, and freezes the result, so the compacted shape is certified like
/// any other. The maps come back as the RemapTable that feeds the history of
/// ADR-0011.

namespace Stereon {

namespace {

constexpr std::uint32_t None = 0xFFFFFFFFu;

/// The dense numbering of a table's live slots, in slot order.
template <typename T>
[[nodiscard]] std::vector<std::uint32_t> DenseNumbering(const PagedArena<T>& table)
{
    std::vector<std::uint32_t> map(table.GetSize(), None);
    std::uint32_t next = 0;
    for (std::uint32_t i = 0; i < table.GetSize(); ++i) {
        if (table[i].Slot.Alive) {
            map[i] = next++;
        }
    }
    return map;
}

/// The dense numbering of the pool entries that `used` marks.
[[nodiscard]] std::vector<std::uint32_t> DenseNumbering(const std::vector<bool>& used)
{
    std::vector<std::uint32_t> map(used.size(), None);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < used.size(); ++i) {
        if (used[i]) {
            map[i] = next++;
        }
    }
    return map;
}

template <EntityKind K>
[[nodiscard]] Id<K> Mapped(const Id<K> old, const std::vector<std::uint32_t>& map) noexcept
{
    return old.IsValid() && old.GetIndex() < map.size() ? Id<K>::FromIndex(map[old.GetIndex()])
                                                        : Id<K>{};
}

template <typename Tag>
[[nodiscard]] GeometryRef<Tag> Mapped(const GeometryRef<Tag> old,
                                      const std::vector<std::uint32_t>& map) noexcept
{
    return old.IsValid() && old.GetIndex() < map.size()
               ? GeometryRef<Tag>::FromIndex(map[old.GetIndex()])
               : GeometryRef<Tag>{};
}

/// Marks in `used` the pool entries that the live records' `member` names.
template <typename Record, typename Ref>
void MarkUsed(const PagedArena<Record>& records, Ref Record::* member, std::vector<bool>& used)
{
    for (const Record& record : records) {
        const Ref ref = record.*member;
        if (record.Slot.Alive && ref.IsValid() && ref.GetIndex() < used.size()) {
            used[ref.GetIndex()] = true;
        }
    }
}

/// Marks the pool entries the live records reference.
struct UsedGeometry {
    std::vector<bool> Curves;
    std::vector<bool> PCurves;
    std::vector<bool> Surfaces;

    explicit UsedGeometry(const Detail::TopologyTables& t)
        : Curves(t.Curves.GetSize(), false)
        , PCurves(t.PCurves.GetSize(), false)
        , Surfaces(t.Surfaces.GetSize(), false)
    {
        MarkUsed(t.Edges, &Detail::EdgeRecord::Curve, Curves);
        MarkUsed(t.Edges, &Detail::EdgeRecord::Cached, Curves);
        MarkUsed(t.Coedges, &Detail::CoedgeRecord::PCurve, PCurves);
        MarkUsed(t.Faces, &Detail::FaceRecord::Surface, Surfaces);
    }
};

/// Copies the live records of `from` into `to`, each rewritten by `rewrite`
/// and given a fresh slot state.
template <typename T, typename Rewrite>
void CopyLive(const PagedArena<T>& from, PagedArena<T>& to, Rewrite rewrite)
{
    for (const T& record : from) {
        if (record.Slot.Alive) {
            T copy = record;
            rewrite(copy);
            copy.Slot = {.Generation = 0, .Alive = true};
            static_cast<void>(to.PushBack(std::move(copy)));
        }
    }
}

template <typename T>
void CopyUsed(const PagedArena<T>& from, PagedArena<T>& to, const std::vector<bool>& used)
{
    for (std::size_t i = 0; i < from.GetSize(); ++i) {
        if (used[i]) {
            static_cast<void>(to.PushBack(from[i]));
        }
    }
}

/// The old-to-new maps Compact rewrites every link with.
struct LinkMaps {
    const std::vector<std::uint32_t>& Vertices;
    const std::vector<std::uint32_t>& Edges;
    const std::vector<std::uint32_t>& Coedges;
    const std::vector<std::uint32_t>& Loops;
    const std::vector<std::uint32_t>& Faces;
    const std::vector<std::uint32_t>& Shells;
    const std::vector<std::uint32_t>& Solids;
    const std::vector<std::uint32_t>& Curves;
    const std::vector<std::uint32_t>& PCurves;
    const std::vector<std::uint32_t>& Surfaces;
};

/// Copies the live entities of `old` into `fresh` with every link mapped.
void CopyEntities(const Detail::TopologyTables& old, Detail::TopologyTables& fresh,
                  const LinkMaps& m)
{
    CopyLive(old.Vertices, fresh.Vertices,
             [&](Detail::VertexRecord& v) { v.FirstEdge = Mapped(v.FirstEdge, m.Edges); });
    CopyLive(old.Edges, fresh.Edges, [&](Detail::EdgeRecord& e) {
        e.Start = Mapped(e.Start, m.Vertices);
        e.End = Mapped(e.End, m.Vertices);
        e.RadialAnchor = Mapped(e.RadialAnchor, m.Coedges);
        e.Primary = Mapped(e.Primary, m.Coedges);
        e.NextAtStart = Mapped(e.NextAtStart, m.Edges);
        e.NextAtEnd = Mapped(e.NextAtEnd, m.Edges);
        e.Curve = Mapped(e.Curve, m.Curves);
        e.Cached = Mapped(e.Cached, m.Curves);
    });
    CopyLive(old.Coedges, fresh.Coedges, [&](Detail::CoedgeRecord& c) {
        c.Next = Mapped(c.Next, m.Coedges);
        c.Prev = Mapped(c.Prev, m.Coedges);
        c.Radial = Mapped(c.Radial, m.Coedges);
        c.Edge = Mapped(c.Edge, m.Edges);
        c.Loop = Mapped(c.Loop, m.Loops);
        c.PCurve = Mapped(c.PCurve, m.PCurves);
    });
    CopyLive(old.Loops, fresh.Loops, [&](Detail::LoopRecord& l) {
        l.Face = Mapped(l.Face, m.Faces);
        l.Anchor = Mapped(l.Anchor, m.Coedges);
        l.NextInFace = Mapped(l.NextInFace, m.Loops);
    });
    CopyLive(old.Faces, fresh.Faces, [&](Detail::FaceRecord& f) {
        f.Shell = Mapped(f.Shell, m.Shells);
        f.FirstLoop = Mapped(f.FirstLoop, m.Loops);
        f.LastLoop = Mapped(f.LastLoop, m.Loops);
        f.NextInShell = Mapped(f.NextInShell, m.Faces);
        f.Surface = Mapped(f.Surface, m.Surfaces);
    });
    CopyLive(old.Shells, fresh.Shells, [&](Detail::ShellRecord& s) {
        s.Solid = Mapped(s.Solid, m.Solids);
        s.FirstFace = Mapped(s.FirstFace, m.Faces);
        s.LastFace = Mapped(s.LastFace, m.Faces);
        s.NextInSolid = Mapped(s.NextInSolid, m.Shells);
    });
    CopyLive(old.Solids, fresh.Solids, [&](Detail::SolidRecord& s) {
        s.FirstShell = Mapped(s.FirstShell, m.Shells);
        s.LastShell = Mapped(s.LastShell, m.Shells);
    });
}

} // namespace

Result<Compacted> Compact(const Shape& shape)
{
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "Compact: the Shape is empty");
    }
    const TopologySnapshot& snapshot = *shape.GetSnapshot();
    const Detail::TopologyTables& old = snapshot.GetTables();

    RemapTable remap;
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Vertex)] = DenseNumbering(old.Vertices);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Edge)] = DenseNumbering(old.Edges);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Coedge)] = DenseNumbering(old.Coedges);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Loop)] = DenseNumbering(old.Loops);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Face)] = DenseNumbering(old.Faces);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Shell)] = DenseNumbering(old.Shells);
    remap.m_Entities[static_cast<std::size_t>(EntityKind::Solid)] = DenseNumbering(old.Solids);
    const UsedGeometry used(old);
    remap.m_Curves = DenseNumbering(used.Curves);
    remap.m_PCurves = DenseNumbering(used.PCurves);
    remap.m_Surfaces = DenseNumbering(used.Surfaces);

    const auto entities = [&](const EntityKind k) -> const std::vector<std::uint32_t>& {
        return remap.m_Entities[static_cast<std::size_t>(k)];
    };
    Detail::TopologyTables fresh;
    CopyEntities(old, fresh,
                 LinkMaps{entities(EntityKind::Vertex), entities(EntityKind::Edge),
                          entities(EntityKind::Coedge), entities(EntityKind::Loop),
                          entities(EntityKind::Face), entities(EntityKind::Shell),
                          entities(EntityKind::Solid), remap.m_Curves, remap.m_PCurves,
                          remap.m_Surfaces});
    for (std::size_t k = 0; k < EntityKindCount; ++k) {
        fresh.LiveCounts[k] = old.LiveCounts[k];
    }
    CopyUsed(old.Curves, fresh.Curves, used.Curves);
    CopyUsed(old.PCurves, fresh.PCurves, used.PCurves);
    CopyUsed(old.Surfaces, fresh.Surfaces, used.Surfaces);

    // Freeze certifies the compacted tables as it certifies any builder's.
    STEREON_TRY(Shape compacted,
                Freeze(ShapeBuilder(std::move(fresh), snapshot.GetStore(), Detail::FreezeKey{})));
    // Built member by member rather than with designated initializers, which
    // the Linux GCC 16.2 image has scrambled in an implementation unit
    // (docs/toolchain/gcc16-feature-probe.md).
    Compacted result;
    result.Shape = std::move(compacted);
    result.Remap = std::move(remap);
    result.SourceStore = snapshot.GetStore();
    return result;
}

} // namespace Stereon
