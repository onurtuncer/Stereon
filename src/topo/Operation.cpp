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
#include <format>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

module Stereon.Topo;

namespace Stereon {

namespace {

[[nodiscard]] bool Before(const HistoryEntry& a, const HistoryEntry& b) noexcept
{
    return std::tuple(static_cast<int>(a.Kind), a.Input, a.Source, a.Target,
                      static_cast<int>(a.Relation), a.Group, static_cast<int>(a.SourceKind)) <
           std::tuple(static_cast<int>(b.Kind), b.Input, b.Source, b.Target,
                      static_cast<int>(b.Relation), b.Group, static_cast<int>(b.SourceKind));
}

} // namespace

History::History(std::vector<HistoryEntry> entries)
    : m_Entries(std::move(entries))
{
    std::ranges::sort(m_Entries, Before);
}

std::vector<HistoryEntry> History::OfSource(const EntityKind kind, const std::uint32_t source,
                                            const std::uint8_t input) const
{
    std::vector<HistoryEntry> found;
    for (const HistoryEntry& e : m_Entries) {
        const bool generator = e.Relation == HistoryRelation::Generated;
        if (e.Kind == kind && !generator && e.Input == input && e.Source == source &&
            source != NoEntity) {
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
    while (first < table.GetSize() && table[first].Slot.IsAlive() &&
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

/// Marks the slot one side of an entry names, if the relation gives it
/// that side; false when the entry is malformed or names no live entity
/// (a published target must be live, ADR-0030).
template <EntityKind K>
[[nodiscard]] bool MarkSide(const bool present, const std::uint32_t index,
                            const TopologyTables* tables, std::vector<std::uint8_t>& marks)
{
    if (present != (index != NoEntity)) {
        return false;
    }
    if (!present) {
        return true;
    }
    if (tables == nullptr || !tables->IsLive<K>(index)) {
        return false;
    }
    marks[index] = 1;
    return true;
}

/// True when every live slot of the table of kind K is marked.
template <EntityKind K>
[[nodiscard]] bool EveryLiveMarked(const TopologyTables* tables,
                                   const std::vector<std::uint8_t>& marks)
{
    for (std::uint32_t i = 0; i < marks.size(); ++i) {
        if (tables->IsLive<K>(i) && marks[i] == 0) {
            return false;
        }
    }
    return true;
}

template <EntityKind K>
[[nodiscard]] bool KindIsComplete(const History& history, const TopologyTables* source,
                                  const TopologyTables& result)
{
    std::vector<std::uint8_t> targeted(result.Table<K>().GetSize(), 0);
    std::vector<std::uint8_t> sourced(source != nullptr ? source->Table<K>().GetSize() : 0, 0);
    for (const HistoryEntry& e : history.GetEntries()) {
        if (e.Kind == K &&
            (!MarkSide<K>(e.Relation != HistoryRelation::Deleted, e.Target, &result, targeted) ||
             !MarkSide<K>(e.Relation != HistoryRelation::Generated, e.Source, source, sourced))) {
            return false;
        }
    }
    // Every live entity of the result has provenance, and every live input
    // entity a disposition.
    return EveryLiveMarked<K>(&result, targeted) && EveryLiveMarked<K>(source, sourced);
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

namespace {

[[nodiscard]] std::string_view KindName(const EntityKind kind) noexcept
{
    constexpr std::string_view names[] = {"vertex", "edge",  "coedge", "loop",
                                          "face",   "shell", "solid"};
    const auto i = static_cast<std::size_t>(kind);
    return i < std::size(names) ? names[i] : "entity";
}

/// Whether slot i of the table of `kind` is live, the kind chosen at run time.
[[nodiscard]] bool LiveIn(const TopologyTables& t, const EntityKind kind, const std::uint32_t i)
{
    switch (kind) {
    case EntityKind::Vertex:
        return t.IsLive<EntityKind::Vertex>(i);
    case EntityKind::Edge:
        return t.IsLive<EntityKind::Edge>(i);
    case EntityKind::Coedge:
        return t.IsLive<EntityKind::Coedge>(i);
    case EntityKind::Loop:
        return t.IsLive<EntityKind::Loop>(i);
    case EntityKind::Face:
        return t.IsLive<EntityKind::Face>(i);
    case EntityKind::Shell:
        return t.IsLive<EntityKind::Shell>(i);
    case EntityKind::Solid:
        return t.IsLive<EntityKind::Solid>(i);
    }
    return false;
}

/// The slots of the table of `kind`, the kind chosen at run time.
[[nodiscard]] std::uint32_t SizeOf(const TopologyTables& t, const EntityKind kind)
{
    switch (kind) {
    case EntityKind::Vertex:
        return static_cast<std::uint32_t>(t.Vertices.GetSize());
    case EntityKind::Edge:
        return static_cast<std::uint32_t>(t.Edges.GetSize());
    case EntityKind::Coedge:
        return static_cast<std::uint32_t>(t.Coedges.GetSize());
    case EntityKind::Loop:
        return static_cast<std::uint32_t>(t.Loops.GetSize());
    case EntityKind::Face:
        return static_cast<std::uint32_t>(t.Faces.GetSize());
    case EntityKind::Shell:
        return static_cast<std::uint32_t>(t.Shells.GetSize());
    case EntityKind::Solid:
        return static_cast<std::uint32_t>(t.Solids.GetSize());
    }
    return 0;
}

constexpr EntityKind AllKinds[] = {EntityKind::Vertex, EntityKind::Edge, EntityKind::Coedge,
                                   EntityKind::Loop,   EntityKind::Face, EntityKind::Shell,
                                   EntityKind::Solid};

/// What the entries account for: the result's entities they target, and
/// each input's entities they dispose of.
struct Accounted {
    std::set<std::pair<int, std::uint32_t>> Targeted;
    std::set<std::tuple<int, int, std::uint32_t>> Disposed;
};

/// The entry's target: live in the result unless the entry is a deletion,
/// which has none.
[[nodiscard]] std::optional<std::string> TargetProblem(const HistoryEntry& e,
                                                       const TopologyTables& result)
{
    const bool deleted = e.Relation == HistoryRelation::Deleted;
    if (deleted != (e.Target == NoEntity)) {
        return std::format("a {} entry with{} a target", ToString(e.Relation),
                           deleted ? "" : "out");
    }
    if (!deleted && !LiveIn(result, e.Kind, e.Target)) {
        return std::format("an entry names {} {} of the result, which is not live",
                           KindName(e.Kind), e.Target);
    }
    return std::nullopt;
}

/// The entry's source: live in its input, of the generator's kind for a
/// Generated entry; only a Generated entry may have none, a creation root.
[[nodiscard]] std::optional<std::string> SourceProblem(const HistoryEntry& e,
                                                       const TopologyTables& input)
{
    const bool generated = e.Relation == HistoryRelation::Generated;
    if (!generated && e.Source == NoEntity) {
        return std::format("a {} entry without a source", ToString(e.Relation));
    }
    const EntityKind kind = generated ? e.SourceKind : e.Kind;
    if (e.Source != NoEntity && !LiveIn(input, kind, e.Source)) {
        return std::format("an entry names {} {} of input {}, which is not live", KindName(kind),
                           e.Source, e.Input);
    }
    return std::nullopt;
}

/// One entry checked (TargetProblem, SourceProblem), and what it accounts
/// for marked.
[[nodiscard]] std::optional<std::string> EntryProblem(const HistoryEntry& e,
                                                      std::span<const TopologyTables* const> inputs,
                                                      const TopologyTables& result, Accounted& seen)
{
    if (e.Input >= inputs.size() || inputs[e.Input] == nullptr) {
        return std::format("an entry names input {} of {}", e.Input, inputs.size());
    }
    if (std::optional<std::string> problem = TargetProblem(e, result)) {
        return problem;
    }
    if (std::optional<std::string> problem = SourceProblem(e, *inputs[e.Input])) {
        return problem;
    }
    if (e.Relation != HistoryRelation::Deleted) {
        seen.Targeted.insert({static_cast<int>(e.Kind), e.Target});
    }
    if (e.Relation != HistoryRelation::Generated) {
        seen.Disposed.insert({e.Input, static_cast<int>(e.Kind), e.Source});
    }
    return std::nullopt;
}

/// A Group's entries: one source continued by several targets (a split),
/// or several sources continued by one target (a merge), all Modified and
/// of one kind.
[[nodiscard]] std::optional<std::string> GroupProblem(const std::vector<HistoryEntry>& members)
{
    std::set<std::pair<int, std::uint32_t>> sources;
    std::set<std::uint32_t> targets;
    for (const HistoryEntry& e : members) {
        if (e.Relation != HistoryRelation::Modified || e.Kind != members.front().Kind) {
            return std::string("a split or merge group mixes kinds or has an entry that is not "
                               "Modified");
        }
        sources.insert({e.Input, e.Source});
        targets.insert(e.Target);
    }
    const bool split = sources.size() == 1 && targets.size() >= 2;
    const bool merge = targets.size() == 1 && sources.size() >= 2;
    if (!split && !merge) {
        return std::format("a group of {} sources and {} targets is neither a split nor a merge",
                           sources.size(), targets.size());
    }
    return std::nullopt;
}

/// Every live entity of the result targeted, every live entity of every
/// input disposed of.
[[nodiscard]] std::optional<std::string>
CoverageProblem(std::span<const TopologyTables* const> inputs, const TopologyTables& result,
                const Accounted& seen)
{
    for (const EntityKind kind : AllKinds) {
        for (std::uint32_t i = 0; i < SizeOf(result, kind); ++i) {
            if (LiveIn(result, kind, i) && !seen.Targeted.contains({static_cast<int>(kind), i})) {
                return std::format("{} {} of the result has no provenance", KindName(kind), i);
            }
        }
        for (std::size_t n = 0; n < inputs.size(); ++n) {
            for (std::uint32_t i = 0; i < SizeOf(*inputs[n], kind); ++i) {
                if (LiveIn(*inputs[n], kind, i) &&
                    !seen.Disposed.contains({static_cast<int>(n), static_cast<int>(kind), i})) {
                    return std::format("{} {} of input {} has no disposition", KindName(kind), i,
                                       n);
                }
            }
        }
    }
    return std::nullopt;
}

} // namespace

std::optional<std::string> HistoryProblem(const History& history,
                                          const std::span<const TopologyTables* const> inputs,
                                          const TopologyTables& result)
{
    Accounted seen;
    std::map<std::uint32_t, std::vector<HistoryEntry>> groups;
    for (const HistoryEntry& e : history.GetEntries()) {
        if (std::optional<std::string> problem = EntryProblem(e, inputs, result, seen)) {
            return problem;
        }
        if (e.Group != NoGroup) {
            groups[e.Group].push_back(e);
        }
    }
    for (const auto& [group, members] : groups) {
        if (std::optional<std::string> problem = GroupProblem(members)) {
            return problem;
        }
    }
    return CoverageProblem(inputs, result, seen);
}

} // namespace Stereon::Detail
