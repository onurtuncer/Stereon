// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

module Stereon.Topo;

namespace Stereon {

namespace {

[[nodiscard]] bool Before(const HistoryEntry& a, const HistoryEntry& b) noexcept
{
    return std::tuple(static_cast<int>(a.Kind), a.Input, a.Source, a.Target,
                      static_cast<int>(a.Relation)) < std::tuple(static_cast<int>(b.Kind), b.Input,
                                                                 b.Source, b.Target,
                                                                 static_cast<int>(b.Relation));
}

} // namespace

History::History(std::vector<HistoryEntry> entries)
    : m_Entries(std::move(entries))
{
    std::ranges::sort(m_Entries, Before);
}

std::vector<HistoryEntry> History::OfSource(const EntityKind kind, const std::uint32_t source) const
{
    std::vector<HistoryEntry> found;
    for (const HistoryEntry& e : m_Entries) {
        if (e.Kind == kind && e.Input == 0 && e.Source == source && source != NoEntity) {
            found.push_back(e);
        }
    }
    return found;
}

std::vector<HistoryEntry> History::OfTarget(const EntityKind kind, const std::uint32_t target) const
{
    std::vector<HistoryEntry> found;
    for (const HistoryEntry& e : m_Entries) {
        if (e.Kind == kind && e.Target == target && target != NoEntity) {
            found.push_back(e);
        }
    }
    return found;
}

std::size_t History::Count(const HistoryRelation relation) const noexcept
{
    return static_cast<std::size_t>(std::ranges::count_if(
        m_Entries, [&](const HistoryEntry& e) { return e.Relation == relation; }));
}

std::size_t History::Count(const EntityKind kind, const HistoryRelation relation) const noexcept
{
    return static_cast<std::size_t>(std::ranges::count_if(m_Entries, [&](const HistoryEntry& e) {
        return e.Kind == kind && e.Relation == relation;
    }));
}

} // namespace Stereon

namespace Stereon::Detail {

namespace {

/// The members of a singly linked list, bounded by the table's size so that
/// a corrupt list cannot loop (the tables here are certified, but the walk
/// costs nothing to bound).
template <typename Table, typename Step>
[[nodiscard]] std::vector<std::uint32_t> Members(const Table& table, std::uint32_t first, Step step)
{
    std::vector<std::uint32_t> members;
    while (first < table.GetSize() && table[first].Slot.Alive &&
           members.size() <= table.GetSize()) {
        members.push_back(first);
        first = step(table[first]);
    }
    return members;
}

/// The coedges of a loop, from its anchor round its cycle.
[[nodiscard]] std::vector<std::uint32_t> CoedgesOfLoop(const TopologyTables& t,
                                                       const std::uint32_t loop)
{
    std::vector<std::uint32_t> coedges;
    const std::uint32_t anchor = t.Loops[loop].Anchor.GetIndex();
    for (std::uint32_t c = anchor;
         t.IsLive<EntityKind::Coedge>(c) && coedges.size() <= t.Coedges.GetSize();) {
        coedges.push_back(c);
        c = t.Coedges[c].Next.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    return coedges;
}

// Whether the entity at slot i is the same in both tables: its own fields
// and what it bounds, not the links that only chain it to its siblings in a
// parent's list (NextInShell and the like), which change when a neighbour
// goes without changing it. A cached edge curve is not the entity either:
// dropping one is reported as a warning, not as a modification.

[[nodiscard]] bool Same(const VertexRecord& a, const VertexRecord& b) noexcept
{
    return a.Position == b.Position && a.Tolerance == b.Tolerance;
}

[[nodiscard]] bool Same(const EdgeRecord& a, const EdgeRecord& b) noexcept
{
    return a.Start == b.Start && a.End == b.End && a.Primary == b.Primary && a.Curve == b.Curve &&
           a.Range == b.Range && a.Tolerance == b.Tolerance && a.Kind == b.Kind;
}

[[nodiscard]] bool Same(const CoedgeRecord& a, const CoedgeRecord& b) noexcept
{
    return a.Next == b.Next && a.Prev == b.Prev && a.Edge == b.Edge && a.Loop == b.Loop &&
           a.PCurve == b.PCurve && a.Sense == b.Sense;
}

template <EntityKind K>
[[nodiscard]] bool SameAt(const TopologyTables& x, const TopologyTables& y, const std::uint32_t i)
{
    if constexpr (K == EntityKind::Vertex || K == EntityKind::Edge || K == EntityKind::Coedge) {
        return Same(x.Table<K>()[i], y.Table<K>()[i]);
    } else if constexpr (K == EntityKind::Loop) {
        return x.Loops[i].Face == y.Loops[i].Face && CoedgesOfLoop(x, i) == CoedgesOfLoop(y, i);
    } else if constexpr (K == EntityKind::Face) {
        const FaceRecord& a = x.Faces[i];
        const FaceRecord& b = y.Faces[i];
        const auto next = [](const LoopRecord& l) { return l.NextInFace.GetIndex(); };
        return a.Surface == b.Surface && a.Reversed == b.Reversed && a.Shell == b.Shell &&
               Members(x.Loops, a.FirstLoop.GetIndex(), next) ==
                   Members(y.Loops, b.FirstLoop.GetIndex(), next);
    } else if constexpr (K == EntityKind::Shell) {
        const ShellRecord& a = x.Shells[i];
        const ShellRecord& b = y.Shells[i];
        const auto next = [](const FaceRecord& f) { return f.NextInShell.GetIndex(); };
        return a.Solid == b.Solid && Members(x.Faces, a.FirstFace.GetIndex(), next) ==
                                         Members(y.Faces, b.FirstFace.GetIndex(), next);
    } else {
        const auto next = [](const ShellRecord& s) { return s.NextInSolid.GetIndex(); };
        return Members(x.Shells, x.Solids[i].FirstShell.GetIndex(), next) ==
               Members(y.Shells, y.Solids[i].FirstShell.GetIndex(), next);
    }
}

[[nodiscard]] HistoryEntry Entry(const EntityKind kind, const HistoryRelation relation,
                                 const std::uint32_t source, const std::uint32_t target)
{
    HistoryEntry entry;
    entry.Kind = kind;
    entry.Relation = relation;
    entry.Source = source;
    entry.Target = target;
    return entry;
}

template <EntityKind K>
void AddKind(const TopologyTables* source, const TopologyTables& result,
             std::vector<HistoryEntry>& entries)
{
    const std::size_t sourceSize = source != nullptr ? source->Table<K>().GetSize() : 0;
    const std::size_t size = std::max(sourceSize, result.Table<K>().GetSize());
    for (std::uint32_t i = 0; i < size; ++i) {
        const bool before = source != nullptr && source->IsLive<K>(i);
        const bool after = result.IsLive<K>(i);
        if (before && after) {
            entries.push_back(Entry(K,
                                    SameAt<K>(*source, result, i) ? HistoryRelation::Unchanged
                                                                  : HistoryRelation::Modified,
                                    i, i));
        } else if (before) {
            entries.push_back(Entry(K, HistoryRelation::Deleted, i, NoEntity));
        } else if (after) {
            entries.push_back(Entry(K, HistoryRelation::Generated, NoEntity, i));
        }
    }
}

template <EntityKind K>
[[nodiscard]] bool KindIsComplete(const History& history, const TopologyTables* source,
                                  const TopologyTables& result)
{
    const auto& table = result.Table<K>();
    std::vector<std::uint8_t> targeted(table.GetSize(), 0);
    std::vector<std::uint8_t> sourced(source != nullptr ? source->Table<K>().GetSize() : 0, 0);
    for (const HistoryEntry& e : history.GetEntries()) {
        if (e.Kind != K) {
            continue;
        }
        const bool hasTarget = e.Relation != HistoryRelation::Deleted;
        const bool hasSource = e.Relation != HistoryRelation::Generated;
        if (hasTarget != (e.Target != NoEntity) || hasSource != (e.Source != NoEntity)) {
            return false;
        }
        if (hasTarget) {
            if (!result.IsLive<K>(e.Target)) {
                return false; // a published target must be a live entity
            }
            targeted[e.Target] = 1;
        }
        if (hasSource) {
            if (source == nullptr || !source->IsLive<K>(e.Source)) {
                return false;
            }
            sourced[e.Source] = 1;
        }
    }
    for (std::uint32_t i = 0; i < table.GetSize(); ++i) {
        if (table[i].Slot.Alive && targeted[i] == 0) {
            return false; // a live entity of the result without provenance
        }
    }
    for (std::uint32_t i = 0; i < sourced.size(); ++i) {
        if (source->IsLive<K>(i) && sourced[i] == 0) {
            return false; // a live input entity without a disposition
        }
    }
    return true;
}

} // namespace

History HistoryBetween(const TopologyTables* source, const TopologyTables& result)
{
    std::vector<HistoryEntry> entries;
    AddKind<EntityKind::Vertex>(source, result, entries);
    AddKind<EntityKind::Edge>(source, result, entries);
    AddKind<EntityKind::Coedge>(source, result, entries);
    AddKind<EntityKind::Loop>(source, result, entries);
    AddKind<EntityKind::Face>(source, result, entries);
    AddKind<EntityKind::Shell>(source, result, entries);
    AddKind<EntityKind::Solid>(source, result, entries);
    return History(std::move(entries));
}

bool HistoryIsComplete(const History& history, const TopologyTables* source,
                       const TopologyTables& result)
{
    return KindIsComplete<EntityKind::Vertex>(history, source, result) &&
           KindIsComplete<EntityKind::Edge>(history, source, result) &&
           KindIsComplete<EntityKind::Coedge>(history, source, result) &&
           KindIsComplete<EntityKind::Loop>(history, source, result) &&
           KindIsComplete<EntityKind::Face>(history, source, result) &&
           KindIsComplete<EntityKind::Shell>(history, source, result) &&
           KindIsComplete<EntityKind::Solid>(history, source, result);
}

} // namespace Stereon::Detail
