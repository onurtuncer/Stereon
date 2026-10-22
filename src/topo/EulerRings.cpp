// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
module Stereon.Topo;

namespace Stereon {
namespace {
constexpr std::uint32_t None = CoedgeId::Invalid;
[[nodiscard]] auto Refuse(const char* op, const char* why)
{
    return Fail(ErrorCode::InvalidInput, std::string(op) + ": " + why);
}
} // namespace

void ShapeBuilder::RecordEvent(const EntityKind kind, const TopologyEventKind event,
                               const std::uint32_t source, const std::uint32_t target)
{
    TopologyEvent e;
    e.Kind = kind;
    e.Event = event;
    e.Source = source;
    e.Target = target;
    m_Events.push_back(e);
}

std::vector<std::uint32_t> ShapeBuilder::UsesOf(const std::uint32_t edge) const
{
    std::vector<std::uint32_t> uses;
    const std::uint32_t anchor = m_Tables.Edges[edge].RadialAnchor.GetIndex();
    for (std::uint32_t c = anchor; c != None && uses.size() < m_Tables.Coedges.GetSize();) {
        uses.push_back(c);
        c = m_Tables.Coedges[c].Radial.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    return uses;
}

void ShapeBuilder::MoveLoop(const std::uint32_t loop, const std::uint32_t face)
{
    UnlinkChild<EntityKind::Face, EntityKind::Loop>(
        loop, &Detail::LoopRecord::Face, &Detail::LoopRecord::NextInFace,
        &Detail::FaceRecord::FirstLoop, &Detail::FaceRecord::LastLoop);
    const LoopId id = LoopId::FromIndex(loop);
    const std::uint32_t last = m_Tables.Faces[face].LastLoop.GetIndex();
    m_Tables.Loops.Modify(loop, [&](Detail::LoopRecord& l) {
        l.Face = FaceId::FromIndex(face);
        l.NextInFace = {};
    });
    if (last != LoopId::Invalid) {
        m_Tables.Loops.Modify(last, [&](Detail::LoopRecord& l) { l.NextInFace = id; });
    }
    m_Tables.Faces.Modify(face, [&](Detail::FaceRecord& f) {
        if (!f.FirstLoop.IsValid()) {
            f.FirstLoop = id;
        }
        f.LastLoop = id;
    });
    Advance();
    RecordEvent(EntityKind::Loop, TopologyEventKind::Modified, loop, loop);
}

Result<BuilderLoopId> ShapeBuilder::KillEdgeMakeRing(const BuilderEdgeId edge)
{
    constexpr const char* op = "KillEdgeMakeRing";
    STEREON_TRY_VOID(Require(edge, op));
    const std::uint32_t first = m_Tables.Edges[edge.GetIndex()].RadialAnchor.GetIndex();
    const std::uint32_t second = SecondOfTwoUses(edge.GetIndex());
    if (second == None) {
        return Refuse(op, "the bridge must have exactly two uses");
    }
    const auto a = m_Tables.Coedges[first];
    const auto b = m_Tables.Coedges[second];
    if (a.Loop != b.Loop || a.Sense == b.Sense || a.Next.GetIndex() == second ||
        b.Next.GetIndex() == first) {
        return Refuse(op, "opposite nonadjacent uses in one loop are required");
    }
    const std::uint32_t loop = a.Loop.GetIndex();
    const std::uint32_t face = m_Tables.Loops[loop].Face.GetIndex();
    const auto cycle = CycleFrom(first);
    const auto at = std::ranges::find(cycle, second);
    std::vector<std::uint32_t> left(cycle.begin() + 1, at);
    std::vector<std::uint32_t> right(at + 1, cycle.end());
    // Keep the old loop's anchor whenever it survives the removal.
    const auto anchor = m_Tables.Loops[loop].Anchor.GetIndex();
    if (std::ranges::find(right, anchor) != right.end()) {
        std::swap(left, right);
    }
    STEREON_TRY(const auto made,
                AddLoop(BuilderFaceId::Make(face, m_Tables.Faces[face].Slot.GetGeneration(),
                                            m_Token.GetValue())));
    RelinkCycle(loop, left);
    RelinkCycle(made.GetIndex(), right);
    KillEdgeAndUses(edge.GetIndex());
    RecordEvent(EntityKind::Edge, TopologyEventKind::Deleted, edge.GetIndex(), NoEntity);
    RecordEvent(EntityKind::Loop, TopologyEventKind::Split, loop, made.GetIndex());
    CheckNeighbourhood({loop, made.GetIndex()}, {});
    return made;
}

Result<EdgeRing> ShapeBuilder::MakeEdgeKillRing(const BuilderCoedgeId a, const BuilderCoedgeId b,
                                                const EdgeGeometry& geometry)
{
    constexpr const char* op = "MakeEdgeKillRing";
    STEREON_TRY_VOID(Require(a, op));
    STEREON_TRY_VOID(Require(b, op));
    STEREON_TRY_VOID(CheckEdgeGeometry(geometry, op));
    const std::uint32_t left = m_Tables.Coedges[a.GetIndex()].Loop.GetIndex();
    const std::uint32_t right = m_Tables.Coedges[b.GetIndex()].Loop.GetIndex();
    if (left == right || m_Tables.Loops[left].Face != m_Tables.Loops[right].Face) {
        return Refuse(op, "distinct loops of one face are required");
    }
    const auto from = EndVertexOf(a.GetIndex());
    const auto to = EndVertexOf(b.GetIndex());
    if (geometry.Kind == EdgeCurveKind::Degenerate && from != to) {
        return Refuse(op, "a degenerate bridge cannot join distinct vertices");
    }
    EdgeSpec spec;
    spec.Kind = geometry.Kind;
    spec.Curve = geometry.Curve;
    spec.Range = geometry.Range;
    spec.Tolerance = geometry.Tolerance;
    spec.Start = BuilderVertexId::Make(from, m_Tables.Vertices[from].Slot.GetGeneration(),
                                       m_Token.GetValue());
    spec.End =
        BuilderVertexId::Make(to, m_Tables.Vertices[to].Slot.GetGeneration(), m_Token.GetValue());
    const auto leftCycle = CycleFrom(a.GetIndex());
    const auto rightCycle = CycleFrom(b.GetIndex());
    EdgeRing made;
    STEREON_TRY(made.Edge, AddEdge(spec));
    made.Out = NewUse(made.Edge.GetIndex(), left, EdgeUseSense::Same, geometry.PCurve);
    made.Back = NewUse(made.Edge.GetIndex(), left, EdgeUseSense::Opposite, geometry.PCurve);
    std::vector<std::uint32_t> joined{a.GetIndex(), made.Out.GetIndex()};
    joined.insert(joined.end(), rightCycle.begin() + 1, rightCycle.end());
    joined.push_back(b.GetIndex());
    joined.push_back(made.Back.GetIndex());
    joined.insert(joined.end(), leftCycle.begin() + 1, leftCycle.end());
    RelinkCycle(left, joined);
    m_Tables.Loops.Modify(right, [](Detail::LoopRecord& l) { l.Anchor = {}; });
    RemoveLoopAt(right);
    RecordEvent(EntityKind::Edge, TopologyEventKind::Created, NoEntity, made.Edge.GetIndex());
    RecordEvent(EntityKind::Loop, TopologyEventKind::Merged, right, left);
    CheckNeighbourhood({left}, {made.Edge.GetIndex()});
    return made;
}

Result<void> ShapeBuilder::KillFaceMakeRingHole(const BuilderFaceId dead, const BuilderFaceId kept)
{
    constexpr const char* op = "KillFaceMakeRingHole";
    STEREON_TRY_VOID(Require(dead, op));
    STEREON_TRY_VOID(Require(kept, op));
    if (dead == kept) {
        return Refuse(op, "two distinct faces are required");
    }
    STEREON_TRY_VOID(CheckJoinable(dead.GetIndex(), kept.GetIndex(), op));
    const auto d = m_Tables.Faces[dead.GetIndex()];
    if (!d.FirstLoop.IsValid() || !m_Tables.Loops[d.FirstLoop.GetIndex()].Anchor.IsValid() ||
        !m_Tables.Faces[kept.GetIndex()].FirstLoop.IsValid()) {
        return Refuse(op,
                      "the removed face needs one nonempty loop and the kept face needs a loop");
    }
    MoveLoop(d.FirstLoop.GetIndex(), kept.GetIndex());
    if (d.Shell.IsValid()) {
        UnlinkChild<EntityKind::Shell, EntityKind::Face>(
            dead.GetIndex(), &Detail::FaceRecord::Shell, &Detail::FaceRecord::NextInShell,
            &Detail::ShellRecord::FirstFace, &Detail::ShellRecord::LastFace);
    }
    Release<EntityKind::Face>(dead.GetIndex());
    RecordEvent(EntityKind::Face, TopologyEventKind::Merged, dead.GetIndex(), kept.GetIndex());
    CheckNeighbourhood({d.FirstLoop.GetIndex()}, {});
    return {};
}

/// Splitting a face is shell-preserving only if the moved loop can reach
/// a remaining loop through radial neighbours in the same shell.
bool ShapeBuilder::CanSeparateLoop(const std::uint32_t moved) const
{
    const auto source = m_Tables.Loops[moved].Face;
    const auto shell = m_Tables.Faces[source.GetIndex()].Shell;
    if (!shell.IsValid()) {
        return true;
    }
    std::vector<bool> seen(m_Tables.Loops.GetSize(), false);
    std::vector<std::uint32_t> pending{moved};
    seen[moved] = true;
    while (!pending.empty()) {
        const auto loop = pending.back();
        pending.pop_back();
        for (const auto c : CycleFrom(m_Tables.Loops[loop].Anchor.GetIndex())) {
            for (const auto use : UsesOf(m_Tables.Coedges[c].Edge.GetIndex())) {
                const auto other = m_Tables.Coedges[use].Loop.GetIndex();
                if (other == moved) {
                    continue;
                }
                const auto face = m_Tables.Loops[other].Face;
                if (face == source) {
                    return true;
                }
                if (m_Tables.Faces[face.GetIndex()].Shell != shell) {
                    continue;
                }
                for (auto next = m_Tables.Faces[face.GetIndex()].FirstLoop; next.IsValid();
                     next = m_Tables.Loops[next.GetIndex()].NextInFace) {
                    if (!seen[next.GetIndex()]) {
                        seen[next.GetIndex()] = true;
                        pending.push_back(next.GetIndex());
                    }
                }
            }
        }
    }
    return false;
}

Result<BuilderFaceId> ShapeBuilder::MakeFaceKillRingHole(const BuilderLoopId loop)
{
    constexpr const char* op = "MakeFaceKillRingHole";
    STEREON_TRY_VOID(Require(loop, op));
    const auto source = m_Tables.Loops[loop.GetIndex()].Face.GetIndex();
    const auto face = m_Tables.Faces[source];
    if (face.FirstLoop == face.LastLoop || !m_Tables.Loops[loop.GetIndex()].Anchor.IsValid()) {
        return Refuse(op, "a nonempty loop of a face with several loops is required");
    }
    if (!CanSeparateLoop(loop.GetIndex())) {
        return Refuse(op, "separating the loop would disconnect its shell");
    }
    STEREON_TRY(const auto made, AddFace(face.Surface, face.Reversed));
    if (face.Shell.IsValid()) {
        const auto shell = face.Shell.GetIndex();
        STEREON_TRY_VOID(
            AddFaceToShell(BuilderShellId::Make(shell, m_Tables.Shells[shell].Slot.GetGeneration(),
                                                m_Token.GetValue()),
                           made));
    }
    MoveLoop(loop.GetIndex(), made.GetIndex());
    RecordEvent(EntityKind::Face, TopologyEventKind::Split, source, made.GetIndex());
    CheckNeighbourhood({loop.GetIndex()}, {});
    return made;
}
} // namespace Stereon
