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
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

module Stereon.Topo;

namespace Stereon {

namespace {

/// The kinds whose edge evaluates through a surface, and so needs a primary
/// coedge (ADR-0009, section 3).
[[nodiscard]] bool NeedsPrimary(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Intersection || kind == EdgeCurveKind::Trim ||
           kind == EdgeCurveKind::Boundary || kind == EdgeCurveKind::Degenerate;
}

/// The kinds with their own 3D curve.
[[nodiscard]] bool HasOwnCurve(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Exact || kind == EdgeCurveKind::Spatial;
}

template <EntityKind K>
[[nodiscard]] Id<K> Link(const BuilderId<K> h) noexcept
{
    return Id<K>::FromIndex(h.GetIndex());
}

} // namespace

ShapeBuilder::ShapeBuilder()
    : m_Token(Detail::NextStoreId())
{}

ShapeBuilder::ShapeBuilder(Detail::TopologyTables tables, const StoreId source, Detail::FreezeKey)
    : m_Tables(std::move(tables))
    , m_Token(Detail::NextStoreId())
    , m_SourceStore(source)
{}

Result<ShapeBuilder> ShapeBuilder::Open(const Shape& shape)
{
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "ShapeBuilder::Open: the Shape is empty");
    }
    const TopologySnapshot& snapshot = *shape.GetSnapshot();
    return ShapeBuilder(snapshot.GetTables(), snapshot.GetStore(), Detail::FreezeKey{});
}

ShapeBuilder::ShapeBuilder(ShapeBuilder&& other) noexcept
    : m_Tables(std::move(other.m_Tables))
    , m_Token(std::exchange(other.m_Token, BuilderToken{}))
    , m_Epoch(other.m_Epoch)
    , m_SourceStore(other.m_SourceStore)
{
    // Cursors on the source are invalid from now on (ADR-0038).
    other.Advance();
    other.m_Tables = {};
}

ShapeBuilder& ShapeBuilder::operator=(ShapeBuilder&& other) noexcept
{
    if (this != &other) {
        m_Tables = std::move(other.m_Tables);
        m_Token = std::exchange(other.m_Token, BuilderToken{});
        m_Epoch = other.m_Epoch + 1;
        m_SourceStore = other.m_SourceStore;
        other.Advance();
        other.m_Tables = {};
    }
    return *this;
}

ShapeBuilder::~ShapeBuilder() = default;

void ShapeBuilder::Advance() noexcept
{
    // ADR-0038: a wrapped epoch could match a stale cursor's. 2^64 changes
    // do not happen; the check makes sure of it.
    if (m_Epoch == std::numeric_limits<std::uint64_t>::max()) {
        std::fputs("Stereon: ShapeBuilder topology epoch exhausted\n", stderr);
        std::abort();
    }
    ++m_Epoch;
}

Result<void> ShapeBuilder::RequireToken(const char* operation) const
{
    if (!m_Token.IsValid()) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(operation) + ": the builder was moved from or frozen");
    }
    return {};
}

template <EntityKind K>
bool ShapeBuilder::IsLive(const BuilderId<K> h) const noexcept
{
    if (!h.IsValid() || !m_Tables.template IsLive<K>(h.GetIndex())) {
        return false;
    }
    if constexpr (STEREON_CHECKED_HANDLES) {
        if (h.GetStore() != m_Token.GetValue()) {
            return false;
        }
    }
    return m_Tables.template Table<K>()[h.GetIndex()].Slot.Generation == h.GetGeneration();
}

template <EntityKind K>
Result<void> ShapeBuilder::Require(const BuilderId<K> h, const char* operation) const
{
    STEREON_TRY_VOID(RequireToken(operation));
    if (!IsLive(h)) {
        return Fail(ErrorCode::InvalidInput, std::string(operation) + ": " +
                                                 std::string(ToString(K)) +
                                                 " handle is not a live entity of this builder");
    }
    return {};
}

template <EntityKind K>
BuilderId<K> ShapeBuilder::Allocate(Detail::Record<K> record)
{
    record.Slot = {.Generation = 0, .Alive = true};
    const std::size_t index = m_Tables.template Table<K>().PushBack(std::move(record));
    ++m_Tables.LiveCounts[static_cast<std::size_t>(K)];
    Advance();
    return BuilderId<K>::Make(static_cast<std::uint32_t>(index), 0, m_Token.GetValue());
}

/// Tombstones the slot: dead, its generation bumped so that every handle to
/// it is stale, never reused by this builder (ADR-0004).
template <EntityKind K>
void ShapeBuilder::Release(const std::uint32_t index)
{
    m_Tables.template Table<K>().Modify(index, [](Detail::Record<K>& r) {
        r.Slot.Alive = false;
        ++r.Slot.Generation;
    });
    --m_Tables.LiveCounts[static_cast<std::size_t>(K)];
    Advance();
}

// --- The vertex-incidence index -----------------------------------------------
//
// Each vertex lists the edges at it in ascending index order, an edge
// continuing the list through NextAtStart or NextAtEnd, whichever of its
// ends the vertex is; a closed edge is listed once, through NextAtStart.

namespace {

/// The link of edge e that continues vertex v's list.
[[nodiscard]] EdgeId& LinkAt(Detail::EdgeRecord& e, const std::uint32_t v) noexcept
{
    return e.Start.GetIndex() == v ? e.NextAtStart : e.NextAtEnd;
}

[[nodiscard]] std::uint32_t NextAt(const Detail::EdgeRecord& e, const std::uint32_t v) noexcept
{
    return (e.Start.GetIndex() == v ? e.NextAtStart : e.NextAtEnd).GetIndex();
}

} // namespace

void ShapeBuilder::LinkAtVertex(const std::uint32_t edge, const std::uint32_t vertex,
                                const bool atStart)
{
    // The list stays ascending: walk to the last edge below this one. A new
    // edge has the largest index, so this is the tail.
    std::uint32_t previous = EdgeId::Invalid;
    std::uint32_t current = m_Tables.Vertices[vertex].FirstEdge.GetIndex();
    for (std::size_t steps = 0;
         current != EdgeId::Invalid && current < edge && steps < m_Tables.Edges.GetSize();
         ++steps) {
        previous = current;
        current = NextAt(m_Tables.Edges[current], vertex);
    }
    const EdgeId after = EdgeId::FromIndex(current);
    m_Tables.Edges.Modify(
        edge, [&](Detail::EdgeRecord& e) { (atStart ? e.NextAtStart : e.NextAtEnd) = after; });
    if (previous == EdgeId::Invalid) {
        m_Tables.Vertices.Modify(
            vertex, [&](Detail::VertexRecord& v) { v.FirstEdge = EdgeId::FromIndex(edge); });
    } else {
        m_Tables.Edges.Modify(
            previous, [&](Detail::EdgeRecord& e) { LinkAt(e, vertex) = EdgeId::FromIndex(edge); });
    }
}

void ShapeBuilder::UnlinkAtVertex(const std::uint32_t edge, const std::uint32_t vertex)
{
    std::uint32_t previous = EdgeId::Invalid;
    std::uint32_t current = m_Tables.Vertices[vertex].FirstEdge.GetIndex();
    for (std::size_t steps = 0;
         current != EdgeId::Invalid && current != edge && steps < m_Tables.Edges.GetSize();
         ++steps) {
        previous = current;
        current = NextAt(m_Tables.Edges[current], vertex);
    }
    if (current != edge) {
        return; // not listed: nothing to unlink
    }
    const EdgeId after = EdgeId::FromIndex(NextAt(m_Tables.Edges[edge], vertex));
    if (previous == EdgeId::Invalid) {
        m_Tables.Vertices.Modify(vertex, [&](Detail::VertexRecord& v) { v.FirstEdge = after; });
    } else {
        m_Tables.Edges.Modify(previous, [&](Detail::EdgeRecord& e) { LinkAt(e, vertex) = after; });
    }
    m_Tables.Edges.Modify(edge, [&](Detail::EdgeRecord& e) { LinkAt(e, vertex) = EdgeId{}; });
}

Result<CurveRef> ShapeBuilder::AddCurve(Curve curve)
{
    STEREON_TRY_VOID(RequireToken("AddCurve"));
    if (curve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddCurve: the curve is empty");
    }
    return CurveRef::FromIndex(
        static_cast<std::uint32_t>(m_Tables.Curves.PushBack(std::move(curve))));
}

Result<PCurveRef> ShapeBuilder::AddPCurve(Curve2d pcurve)
{
    STEREON_TRY_VOID(RequireToken("AddPCurve"));
    if (pcurve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddPCurve: the pcurve is empty");
    }
    return PCurveRef::FromIndex(
        static_cast<std::uint32_t>(m_Tables.PCurves.PushBack(std::move(pcurve))));
}

Result<SurfaceRef> ShapeBuilder::AddSurface(Surface surface)
{
    STEREON_TRY_VOID(RequireToken("AddSurface"));
    if (surface.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddSurface: the surface is empty");
    }
    return SurfaceRef::FromIndex(
        static_cast<std::uint32_t>(m_Tables.Surfaces.PushBack(std::move(surface))));
}

Result<BuilderVertexId> ShapeBuilder::AddVertex(const Point3& position, const double tolerance)
{
    STEREON_TRY_VOID(RequireToken("AddVertex"));
    STEREON_TRY_VOID(Detail::CheckFinite(position, "the position", "AddVertex"));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the tolerance", "AddVertex"));
    Detail::VertexRecord record;
    record.Position = position;
    record.Tolerance = tolerance;
    return Allocate<EntityKind::Vertex>(record);
}

Result<BuilderEdgeId> ShapeBuilder::AddEdge(const EdgeSpec& spec)
{
    STEREON_TRY_VOID(Require(spec.Start, "AddEdge"));
    STEREON_TRY_VOID(Require(spec.End, "AddEdge"));
    STEREON_TRY_VOID(Detail::CheckNonNegative(spec.Tolerance, "the tolerance", "AddEdge"));
    if (spec.Range.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddEdge: an empty parameter range");
    }
    if (HasOwnCurve(spec.Kind) != spec.Curve.IsValid()) {
        return Fail(
            ErrorCode::InvalidInput,
            std::string("AddEdge: a ") + std::string(ToString(spec.Kind)) + " edge " +
                (HasOwnCurve(spec.Kind) ? "needs its 3D curve" : "has no 3D curve of its own"));
    }
    if (spec.Curve.IsValid() && spec.Curve.GetIndex() >= m_Tables.Curves.GetSize()) {
        return Fail(ErrorCode::InvalidInput, "AddEdge: the curve is not in this builder's pool");
    }
    if (spec.Kind == EdgeCurveKind::Degenerate && spec.Start != spec.End) {
        return Fail(ErrorCode::InvalidInput, "AddEdge: a degenerate edge has one vertex");
    }
    Detail::EdgeRecord record;
    record.Start = Link(spec.Start);
    record.End = Link(spec.End);
    record.Curve = spec.Curve;
    record.Range = spec.Range;
    record.Tolerance = spec.Tolerance;
    record.Kind = spec.Kind;
    const BuilderEdgeId edge = Allocate<EntityKind::Edge>(record);
    LinkAtVertex(edge.GetIndex(), spec.Start.GetIndex(), true);
    if (spec.End.GetIndex() != spec.Start.GetIndex()) {
        LinkAtVertex(edge.GetIndex(), spec.End.GetIndex(), false);
    }
    return edge;
}

Result<BuilderFaceId> ShapeBuilder::AddFace(const SurfaceRef surface, const bool reversed)
{
    STEREON_TRY_VOID(RequireToken("AddFace"));
    if (!surface.IsValid() || surface.GetIndex() >= m_Tables.Surfaces.GetSize()) {
        return Fail(ErrorCode::InvalidInput, "AddFace: the surface is not in this builder's pool");
    }
    Detail::FaceRecord record;
    record.Surface = surface;
    record.Reversed = reversed;
    return Allocate<EntityKind::Face>(record);
}

Result<BuilderLoopId> ShapeBuilder::AddLoop(const BuilderFaceId face)
{
    STEREON_TRY_VOID(Require(face, "AddLoop"));
    Detail::LoopRecord record;
    record.Face = Link(face);
    const BuilderLoopId loop = Allocate<EntityKind::Loop>(record);
    // Append to the face's list.
    const std::uint32_t last = m_Tables.Faces[face.GetIndex()].LastLoop.GetIndex();
    if (last == LoopId::Invalid) {
        m_Tables.Faces.Modify(face.GetIndex(), [&](Detail::FaceRecord& f) {
            f.FirstLoop = Link(loop);
            f.LastLoop = Link(loop);
        });
    } else {
        m_Tables.Loops.Modify(last, [&](Detail::LoopRecord& l) { l.NextInFace = Link(loop); });
        m_Tables.Faces.Modify(face.GetIndex(),
                              [&](Detail::FaceRecord& f) { f.LastLoop = Link(loop); });
    }
    return loop;
}

Result<BuilderCoedgeId> ShapeBuilder::AppendCoedge(const BuilderLoopId loop,
                                                   const BuilderEdgeId edge,
                                                   const EdgeUseSense sense, const PCurveRef pcurve)
{
    STEREON_TRY_VOID(Require(loop, "AppendCoedge"));
    STEREON_TRY_VOID(Require(edge, "AppendCoedge"));
    if (pcurve.IsValid() && pcurve.GetIndex() >= m_Tables.PCurves.GetSize()) {
        return Fail(ErrorCode::InvalidInput,
                    "AppendCoedge: the pcurve is not in this builder's pool");
    }
    const Detail::EdgeRecord& edgeRecord = m_Tables.Edges[edge.GetIndex()];
    if (!pcurve.IsValid() && edgeRecord.Kind != EdgeCurveKind::Spatial) {
        return Fail(ErrorCode::InvalidInput, std::string("AppendCoedge: a use of a ") +
                                                 std::string(ToString(edgeRecord.Kind)) +
                                                 " edge needs its pcurve");
    }

    Detail::CoedgeRecord record;
    record.Edge = Link(edge);
    record.Loop = Link(loop);
    record.PCurve = pcurve;
    record.Sense = sense;
    const BuilderCoedgeId coedge = Allocate<EntityKind::Coedge>(record);
    const CoedgeId self = Link(coedge);

    // Loop cycle: insert before the anchor, at the end of the walk.
    const CoedgeId anchor = m_Tables.Loops[loop.GetIndex()].Anchor;
    if (!anchor.IsValid()) {
        m_Tables.Coedges.Modify(coedge.GetIndex(), [&](Detail::CoedgeRecord& c) {
            c.Next = self;
            c.Prev = self;
        });
        m_Tables.Loops.Modify(loop.GetIndex(), [&](Detail::LoopRecord& l) { l.Anchor = self; });
    } else {
        const CoedgeId last = m_Tables.Coedges[anchor.GetIndex()].Prev;
        m_Tables.Coedges.Modify(coedge.GetIndex(), [&](Detail::CoedgeRecord& c) {
            c.Next = anchor;
            c.Prev = last;
        });
        m_Tables.Coedges.Modify(last.GetIndex(), [&](Detail::CoedgeRecord& c) { c.Next = self; });
        m_Tables.Coedges.Modify(anchor.GetIndex(), [&](Detail::CoedgeRecord& c) { c.Prev = self; });
    }

    // Radial ring: insert after the anchor.
    const CoedgeId radialAnchor = edgeRecord.RadialAnchor;
    if (!radialAnchor.IsValid()) {
        m_Tables.Coedges.Modify(coedge.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Radial = self; });
        m_Tables.Edges.Modify(edge.GetIndex(), [&](Detail::EdgeRecord& e) {
            e.RadialAnchor = self;
            if (NeedsPrimary(e.Kind)) {
                e.Primary = self;
            }
        });
    } else {
        const CoedgeId after = m_Tables.Coedges[radialAnchor.GetIndex()].Radial;
        m_Tables.Coedges.Modify(coedge.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Radial = after; });
        m_Tables.Coedges.Modify(radialAnchor.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Radial = self; });
    }
    return coedge;
}

Result<void> ShapeBuilder::SetPrimaryCoedge(const BuilderEdgeId edge, const BuilderCoedgeId coedge)
{
    STEREON_TRY_VOID(Require(edge, "SetPrimaryCoedge"));
    STEREON_TRY_VOID(Require(coedge, "SetPrimaryCoedge"));
    if (m_Tables.Coedges[coedge.GetIndex()].Edge.GetIndex() != edge.GetIndex()) {
        return Fail(ErrorCode::InvalidInput,
                    "SetPrimaryCoedge: the coedge is not a use of the edge");
    }
    m_Tables.Edges.Modify(edge.GetIndex(),
                          [&](Detail::EdgeRecord& e) { e.Primary = Link(coedge); });
    Advance();
    return {};
}

Result<BuilderSolidId> ShapeBuilder::AddSolid()
{
    STEREON_TRY_VOID(RequireToken("AddSolid"));
    return Allocate<EntityKind::Solid>({});
}

Result<BuilderShellId> ShapeBuilder::AddShell(const BuilderSolidId solid)
{
    STEREON_TRY_VOID(Require(solid, "AddShell"));
    Detail::ShellRecord record;
    record.Solid = Link(solid);
    const BuilderShellId shell = Allocate<EntityKind::Shell>(record);
    const std::uint32_t last = m_Tables.Solids[solid.GetIndex()].LastShell.GetIndex();
    if (last == ShellId::Invalid) {
        m_Tables.Solids.Modify(solid.GetIndex(), [&](Detail::SolidRecord& s) {
            s.FirstShell = Link(shell);
            s.LastShell = Link(shell);
        });
    } else {
        m_Tables.Shells.Modify(last, [&](Detail::ShellRecord& s) { s.NextInSolid = Link(shell); });
        m_Tables.Solids.Modify(solid.GetIndex(),
                               [&](Detail::SolidRecord& s) { s.LastShell = Link(shell); });
    }
    return shell;
}

Result<BuilderShellId> ShapeBuilder::AddShell()
{
    STEREON_TRY_VOID(RequireToken("AddShell"));
    return Allocate<EntityKind::Shell>({});
}

Result<void> ShapeBuilder::AddFaceToShell(const BuilderShellId shell, const BuilderFaceId face)
{
    STEREON_TRY_VOID(Require(shell, "AddFaceToShell"));
    STEREON_TRY_VOID(Require(face, "AddFaceToShell"));
    if (m_Tables.Faces[face.GetIndex()].Shell.IsValid()) {
        return Fail(ErrorCode::InvalidInput, "AddFaceToShell: the face is already in a shell");
    }
    m_Tables.Faces.Modify(face.GetIndex(), [&](Detail::FaceRecord& f) { f.Shell = Link(shell); });
    const std::uint32_t last = m_Tables.Shells[shell.GetIndex()].LastFace.GetIndex();
    if (last == FaceId::Invalid) {
        m_Tables.Shells.Modify(shell.GetIndex(), [&](Detail::ShellRecord& s) {
            s.FirstFace = Link(face);
            s.LastFace = Link(face);
        });
    } else {
        m_Tables.Faces.Modify(last, [&](Detail::FaceRecord& f) { f.NextInShell = Link(face); });
        m_Tables.Shells.Modify(shell.GetIndex(),
                               [&](Detail::ShellRecord& s) { s.LastFace = Link(face); });
    }
    Advance();
    return {};
}

// --- Editing --------------------------------------------------------------------

Result<void> ShapeBuilder::SetSurface(const BuilderFaceId face, const SurfaceRef surface)
{
    STEREON_TRY_VOID(Require(face, "SetSurface"));
    if (!surface.IsValid() || surface.GetIndex() >= m_Tables.Surfaces.GetSize()) {
        return Fail(ErrorCode::InvalidInput,
                    "SetSurface: the surface is not in this builder's pool");
    }
    m_Tables.Faces.Modify(face.GetIndex(), [&](Detail::FaceRecord& f) { f.Surface = surface; });
    Advance();
    return {};
}

/// Takes the coedge out of its loop's cycle and its edge's radial ring and
/// tombstones it. The edge's primary, when it is this coedge, passes to the
/// next use in the ring if `passPrimaryOn`, and otherwise is dropped, which
/// the caller allows only for an edge left with no use.
void ShapeBuilder::DetachCoedge(const std::uint32_t coedge, const bool passPrimaryOn)
{
    const Detail::CoedgeRecord r = m_Tables.Coedges[coedge];
    const CoedgeId self = CoedgeId::FromIndex(coedge);
    // The loop cycle.
    if (r.Next.GetIndex() == coedge) {
        m_Tables.Loops.Modify(r.Loop.GetIndex(), [](Detail::LoopRecord& l) { l.Anchor = {}; });
    } else {
        m_Tables.Coedges.Modify(r.Prev.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Next = r.Next; });
        m_Tables.Coedges.Modify(r.Next.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Prev = r.Prev; });
        if (m_Tables.Loops[r.Loop.GetIndex()].Anchor == self) {
            m_Tables.Loops.Modify(r.Loop.GetIndex(),
                                  [&](Detail::LoopRecord& l) { l.Anchor = r.Next; });
        }
    }
    // The radial ring, singly linked: find the use before this one.
    std::uint32_t before = coedge;
    for (std::size_t steps = 0;
         m_Tables.Coedges[before].Radial.GetIndex() != coedge && steps < m_Tables.Coedges.GetSize();
         ++steps) {
        before = m_Tables.Coedges[before].Radial.GetIndex();
    }
    const bool alone = before == coedge;
    if (!alone) {
        m_Tables.Coedges.Modify(before, [&](Detail::CoedgeRecord& c) { c.Radial = r.Radial; });
    }
    m_Tables.Edges.Modify(r.Edge.GetIndex(), [&](Detail::EdgeRecord& e) {
        if (e.RadialAnchor == self) {
            e.RadialAnchor = alone ? CoedgeId{} : r.Radial;
        }
        if (e.Primary == self) {
            e.Primary = !alone && passPrimaryOn ? r.Radial : CoedgeId{};
        }
    });
    Release<EntityKind::Coedge>(coedge);
}

Result<void> ShapeBuilder::RemoveCoedge(const BuilderCoedgeId coedge)
{
    STEREON_TRY_VOID(Require(coedge, "RemoveCoedge"));
    const Detail::CoedgeRecord& r = m_Tables.Coedges[coedge.GetIndex()];
    const Detail::EdgeRecord& e = m_Tables.Edges[r.Edge.GetIndex()];
    if (e.Primary.GetIndex() == coedge.GetIndex() && r.Radial.GetIndex() != coedge.GetIndex()) {
        return Fail(ErrorCode::InvalidInput,
                    "RemoveCoedge: the coedge is the primary of " +
                        Detail::DescribeEntity(m_Tables, EntityKind::Edge, r.Edge.GetIndex()) +
                        ", which has other uses; SetPrimaryCoedge first");
    }
    DetachCoedge(coedge.GetIndex(), false);
    return {};
}

/// Removes the loop's coedges and the loop itself from its face's list.
void ShapeBuilder::RemoveLoopAt(const std::uint32_t loop)
{
    const std::uint32_t anchor = m_Tables.Loops[loop].Anchor.GetIndex();
    if (anchor != CoedgeId::Invalid) {
        std::vector<std::uint32_t> coedges;
        std::uint32_t c = anchor;
        for (std::size_t steps = 0; steps < m_Tables.Coedges.GetSize(); ++steps) {
            coedges.push_back(c);
            c = m_Tables.Coedges[c].Next.GetIndex();
            if (c == anchor || c == CoedgeId::Invalid) {
                break;
            }
        }
        for (const std::uint32_t coedge : coedges) {
            DetachCoedge(coedge, true);
        }
    }
    UnlinkChild<EntityKind::Face, EntityKind::Loop>(
        loop, &Detail::LoopRecord::Face, &Detail::LoopRecord::NextInFace,
        &Detail::FaceRecord::FirstLoop, &Detail::FaceRecord::LastLoop);
    Release<EntityKind::Loop>(loop);
}

Result<void> ShapeBuilder::RemoveLoop(const BuilderLoopId loop)
{
    STEREON_TRY_VOID(Require(loop, "RemoveLoop"));
    RemoveLoopAt(loop.GetIndex());
    return {};
}

Result<void> ShapeBuilder::RemoveFace(const BuilderFaceId face)
{
    STEREON_TRY_VOID(Require(face, "RemoveFace"));
    std::vector<std::uint32_t> loops;
    for (std::uint32_t l = m_Tables.Faces[face.GetIndex()].FirstLoop.GetIndex();
         l != LoopId::Invalid && loops.size() < m_Tables.Loops.GetSize();
         l = m_Tables.Loops[l].NextInFace.GetIndex()) {
        loops.push_back(l);
    }
    for (const std::uint32_t loop : loops) {
        RemoveLoopAt(loop);
    }
    if (m_Tables.Faces[face.GetIndex()].Shell.IsValid()) {
        UnlinkChild<EntityKind::Shell, EntityKind::Face>(
            face.GetIndex(), &Detail::FaceRecord::Shell, &Detail::FaceRecord::NextInShell,
            &Detail::ShellRecord::FirstFace, &Detail::ShellRecord::LastFace);
    }
    Release<EntityKind::Face>(face.GetIndex());
    return {};
}

Result<void> ShapeBuilder::RemoveShell(const BuilderShellId shell)
{
    STEREON_TRY_VOID(Require(shell, "RemoveShell"));
    std::vector<std::uint32_t> faces;
    for (std::uint32_t f = m_Tables.Shells[shell.GetIndex()].FirstFace.GetIndex();
         f != FaceId::Invalid && faces.size() < m_Tables.Faces.GetSize();
         f = m_Tables.Faces[f].NextInShell.GetIndex()) {
        faces.push_back(f);
    }
    for (const std::uint32_t f : faces) {
        m_Tables.Faces.Modify(f, [](Detail::FaceRecord& r) {
            r.Shell = {};
            r.NextInShell = {};
        });
    }
    if (m_Tables.Shells[shell.GetIndex()].Solid.IsValid()) {
        UnlinkChild<EntityKind::Solid, EntityKind::Shell>(
            shell.GetIndex(), &Detail::ShellRecord::Solid, &Detail::ShellRecord::NextInSolid,
            &Detail::SolidRecord::FirstShell, &Detail::SolidRecord::LastShell);
    }
    Release<EntityKind::Shell>(shell.GetIndex());
    return {};
}

Result<void> ShapeBuilder::RemoveSolid(const BuilderSolidId solid)
{
    STEREON_TRY_VOID(Require(solid, "RemoveSolid"));
    std::vector<std::uint32_t> shells;
    for (std::uint32_t s = m_Tables.Solids[solid.GetIndex()].FirstShell.GetIndex();
         s != ShellId::Invalid && shells.size() < m_Tables.Shells.GetSize();
         s = m_Tables.Shells[s].NextInSolid.GetIndex()) {
        shells.push_back(s);
    }
    for (const std::uint32_t s : shells) {
        m_Tables.Shells.Modify(s, [](Detail::ShellRecord& r) {
            r.Solid = {};
            r.NextInSolid = {};
        });
    }
    Release<EntityKind::Solid>(solid.GetIndex());
    return {};
}

Result<void> ShapeBuilder::RemoveEdge(const BuilderEdgeId edge)
{
    STEREON_TRY_VOID(Require(edge, "RemoveEdge"));
    const Detail::EdgeRecord& e = m_Tables.Edges[edge.GetIndex()];
    if (e.RadialAnchor.IsValid()) {
        return Fail(
            ErrorCode::InvalidInput,
            "RemoveEdge: " + Detail::DescribeEntity(m_Tables, EntityKind::Edge, edge.GetIndex()) +
                " has uses; remove its coedges first");
    }
    const std::uint32_t start = e.Start.GetIndex();
    const std::uint32_t end = e.End.GetIndex();
    UnlinkAtVertex(edge.GetIndex(), start);
    if (end != start) {
        UnlinkAtVertex(edge.GetIndex(), end);
    }
    Release<EntityKind::Edge>(edge.GetIndex());
    return {};
}

Result<void> ShapeBuilder::RemoveVertex(const BuilderVertexId vertex)
{
    STEREON_TRY_VOID(Require(vertex, "RemoveVertex"));
    if (m_Tables.Vertices[vertex.GetIndex()].FirstEdge.IsValid()) {
        return Fail(ErrorCode::InvalidInput,
                    "RemoveVertex: " +
                        Detail::DescribeEntity(m_Tables, EntityKind::Vertex, vertex.GetIndex()) +
                        " has edges; remove them first");
    }
    Release<EntityKind::Vertex>(vertex.GetIndex());
    return {};
}

/// Takes `child` out of its parent's singly linked list (first, last, next)
/// and clears its parent and next links.
template <EntityKind P, EntityKind C>
void ShapeBuilder::UnlinkChild(const std::uint32_t child, Id<P> Detail::Record<C>::* parent,
                               Id<C> Detail::Record<C>::* next, Id<C> Detail::Record<P>::* first,
                               Id<C> Detail::Record<P>::* last)
{
    auto& parents = m_Tables.template Table<P>();
    auto& children = m_Tables.template Table<C>();
    const std::uint32_t p = (children[child].*parent).GetIndex();
    std::uint32_t previous = Id<C>::Invalid;
    for (std::uint32_t c = (parents[p].*first).GetIndex(), steps = 0;
         c != child && c != Id<C>::Invalid && steps < children.GetSize(); ++steps) {
        previous = c;
        c = (children[c].*next).GetIndex();
    }
    const Id<C> after = children[child].*next;
    if (previous == Id<C>::Invalid) {
        parents.Modify(p, [&](Detail::Record<P>& r) { r.*first = after; });
    } else {
        children.Modify(previous, [&](Detail::Record<C>& r) { r.*next = after; });
    }
    if ((parents[p].*last).GetIndex() == child) {
        parents.Modify(p, [&](Detail::Record<P>& r) { r.*last = Id<C>::FromIndex(previous); });
    }
    children.Modify(child, [&](Detail::Record<C>& r) {
        r.*parent = {};
        r.*next = {};
    });
}

ValidityReport ShapeBuilder::Validate() const
{
    return Detail::ValidateStructure(m_Tables);
}

Result<Shape> Freeze(ShapeBuilder&& builder)
{
    STEREON_TRY_VOID(builder.RequireToken("Freeze"));
    // The builder is consumed whatever happens (ADR-0027: a failed private
    // candidate is discarded, not reused).
    ShapeBuilder consumed(std::move(builder));
    const ValidityReport report = Detail::ValidateStructure(consumed.m_Tables);
    if (!report.IsValid()) {
        return Fail(ErrorCode::InvalidInput,
                    "Freeze: the topology is not structurally valid:\n" + report.ToText());
    }
    auto snapshot = std::make_shared<const TopologySnapshot>(
        Detail::FreezeKey{}, std::move(consumed.m_Tables), Detail::NextStoreId(),
        consumed.m_Token.GetValue());
    consumed.m_Token = BuilderToken{};
    return Shape(Detail::FreezeKey{}, std::move(snapshot));
}

ValidityReport ValidateTopology(const Shape& shape)
{
    if (shape.IsEmpty()) {
        return {};
    }
    return Detail::ValidateStructure(shape.GetSnapshot()->GetTables());
}

} // namespace Stereon

namespace Stereon::Testing {

void TopologyFaults::SetNext(ShapeBuilder& b, const BuilderCoedgeId c, const BuilderCoedgeId next)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(), [&](Detail::CoedgeRecord& r) { r.Next = Link(next); });
}

void TopologyFaults::SetPrev(ShapeBuilder& b, const BuilderCoedgeId c, const BuilderCoedgeId prev)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(), [&](Detail::CoedgeRecord& r) { r.Prev = Link(prev); });
}

void TopologyFaults::SetRadial(ShapeBuilder& b, const BuilderCoedgeId c,
                               const BuilderCoedgeId radial)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(),
                              [&](Detail::CoedgeRecord& r) { r.Radial = Link(radial); });
}

void TopologyFaults::SetCoedgeLoop(ShapeBuilder& b, const BuilderCoedgeId c,
                                   const BuilderLoopId loop)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(), [&](Detail::CoedgeRecord& r) { r.Loop = Link(loop); });
}

void TopologyFaults::SetPrimary(ShapeBuilder& b, const BuilderEdgeId e,
                                const BuilderCoedgeId primary)
{
    b.m_Tables.Edges.Modify(e.GetIndex(),
                            [&](Detail::EdgeRecord& r) { r.Primary = Link(primary); });
}

void TopologyFaults::SetKind(ShapeBuilder& b, const BuilderEdgeId e, const EdgeCurveKind kind)
{
    b.m_Tables.Edges.Modify(e.GetIndex(), [&](Detail::EdgeRecord& r) { r.Kind = kind; });
}

void TopologyFaults::SetEnd(ShapeBuilder& b, const BuilderEdgeId e, const BuilderVertexId end)
{
    b.m_Tables.Edges.Modify(e.GetIndex(), [&](Detail::EdgeRecord& r) { r.End = Link(end); });
}

void TopologyFaults::SetPCurve(ShapeBuilder& b, const BuilderCoedgeId c, const PCurveRef pcurve)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(), [&](Detail::CoedgeRecord& r) { r.PCurve = pcurve; });
}

void TopologyFaults::SetFaceShell(ShapeBuilder& b, const BuilderFaceId f,
                                  const BuilderShellId shell)
{
    b.m_Tables.Faces.Modify(f.GetIndex(), [&](Detail::FaceRecord& r) { r.Shell = Link(shell); });
}

void TopologyFaults::SetLoopFace(ShapeBuilder& b, const BuilderLoopId l, const BuilderFaceId face)
{
    b.m_Tables.Loops.Modify(l.GetIndex(), [&](Detail::LoopRecord& r) { r.Face = Link(face); });
}

void TopologyFaults::SetShellSolid(ShapeBuilder& b, const BuilderShellId s,
                                   const BuilderSolidId solid)
{
    b.m_Tables.Shells.Modify(s.GetIndex(), [&](Detail::ShellRecord& r) { r.Solid = Link(solid); });
}

void TopologyFaults::SetFirstEdge(ShapeBuilder& b, const BuilderVertexId v,
                                  const BuilderEdgeId first)
{
    b.m_Tables.Vertices.Modify(v.GetIndex(),
                               [&](Detail::VertexRecord& r) { r.FirstEdge = Link(first); });
}

void TopologyFaults::SetNextAtStart(ShapeBuilder& b, const BuilderEdgeId e,
                                    const BuilderEdgeId next)
{
    b.m_Tables.Edges.Modify(e.GetIndex(),
                            [&](Detail::EdgeRecord& r) { r.NextAtStart = Link(next); });
}

void TopologyFaults::SetNextAtEnd(ShapeBuilder& b, const BuilderEdgeId e, const BuilderEdgeId next)
{
    b.m_Tables.Edges.Modify(e.GetIndex(), [&](Detail::EdgeRecord& r) { r.NextAtEnd = Link(next); });
}

void TopologyFaults::KillCoedge(ShapeBuilder& b, const BuilderCoedgeId c)
{
    b.m_Tables.Coedges.Modify(c.GetIndex(), [](Detail::CoedgeRecord& r) {
        r.Slot.Alive = false;
        ++r.Slot.Generation;
    });
    --b.m_Tables.LiveCounts[static_cast<std::size_t>(EntityKind::Coedge)];
}

void TopologyFaults::SetNextOutOfRange(ShapeBuilder& b, const BuilderCoedgeId c)
{
    const auto past = static_cast<std::uint32_t>(b.m_Tables.Coedges.GetSize() + 7);
    b.m_Tables.Coedges.Modify(c.GetIndex(),
                              [&](Detail::CoedgeRecord& r) { r.Next = CoedgeId::FromIndex(past); });
}

void TopologyFaults::EmptyCurve(ShapeBuilder& b, const CurveRef curve)
{
    b.m_Tables.Curves.Set(curve.GetIndex(), Curve{});
}

void TopologyFaults::EmptyPCurve(ShapeBuilder& b, const PCurveRef pcurve)
{
    b.m_Tables.PCurves.Set(pcurve.GetIndex(), Curve2d{});
}

void TopologyFaults::EmptySurface(ShapeBuilder& b, const SurfaceRef surface)
{
    b.m_Tables.Surfaces.Set(surface.GetIndex(), Surface{});
}

} // namespace Stereon::Testing
