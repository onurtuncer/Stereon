// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include "Profiling.hpp"

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
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

/// Materialise reserved, dead indices before issuing another handle. A
/// failure here issues no handle and cannot release the reservation.
template <EntityKind K, typename Table>
void Tombstone(Table& table, const std::size_t size)
{
    while (table.GetSize() < size) {
        Detail::Record<K> record{};
        record.Slot.Kill();
        table.PushBack(std::move(record));
    }
}

/// Geometry references have no generation: reserve indices with empty
/// wrappers. This runs on allocation, never on failure cleanup.
template <typename Pool, typename Wrapper>
std::size_t AppendGeometry(Pool& pool, std::size_t& reserved, Wrapper value)
{
    while (pool.GetSize() < reserved) {
        static_cast<void>(pool.PushBack(Wrapper{}));
    }
    reserved = pool.GetSize() + 1;
    return pool.PushBack(std::move(value));
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
    // ADR-0004, safeguard 5: checked builds validate a public operation's
    // input again rather than trust the certificate it carries.
    if constexpr (STEREON_CHECKED_HANDLES) {
        const ValidityReport report = Detail::ValidateStructure(snapshot.GetTables());
        if (!report.IsValid()) {
            return Fail(ErrorCode::InvalidInput,
                        "ShapeBuilder::Open: the Shape is not structurally valid:\n" +
                            report.ToText());
        }
    }
    ShapeBuilder builder(snapshot.GetTables(), snapshot.GetStore(), Detail::FreezeKey{});
    builder.m_Source = shape.GetSnapshot();
    return builder;
}

ShapeBuilder::ShapeBuilder(ShapeBuilder&& other) noexcept
    : m_Tables(std::move(other.m_Tables))
    , m_Token(std::exchange(other.m_Token, BuilderToken{}))
    , m_ReservedEntities(other.m_ReservedEntities)
    , m_ReservedGeometry(other.m_ReservedGeometry)
    , m_Epoch(other.m_Epoch)
    , m_SourceStore(other.m_SourceStore)
    , m_Source(std::move(other.m_Source))
    , m_ToleranceChanges(std::move(other.m_ToleranceChanges))
    , m_Events(std::move(other.m_Events))
{
    // Cursors on the source are invalid from now on (ADR-0038).
    other.Advance();
    other.m_Tables = {};
    other.m_ToleranceChanges.clear();
    other.m_Events.clear();
}

ShapeBuilder& ShapeBuilder::operator=(ShapeBuilder&& other) noexcept
{
    if (this != &other) {
        m_Tables = std::move(other.m_Tables);
        m_Token = std::exchange(other.m_Token, BuilderToken{});
        m_ReservedEntities = other.m_ReservedEntities;
        m_ReservedGeometry = other.m_ReservedGeometry;
        m_Epoch = other.m_Epoch + 1;
        m_SourceStore = other.m_SourceStore;
        m_Source = std::move(other.m_Source);
        m_ToleranceChanges = std::move(other.m_ToleranceChanges);
        m_Events = std::move(other.m_Events);
        other.Advance();
        other.m_Tables = {};
        other.m_ToleranceChanges.clear();
        other.m_Events.clear();
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
    return m_Tables.template Table<K>()[h.GetIndex()].Slot.GetGeneration() == h.GetGeneration();
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
    auto& table = m_Tables.template Table<K>();
    std::size_t& reserved = m_ReservedEntities[static_cast<std::size_t>(K)];
    Tombstone<K>(table, reserved);
    reserved = table.GetSize() + 1;
    record.Slot = Detail::SlotState::Live();
    const std::size_t index = table.PushBack(std::move(record));
    ++m_Tables.LiveCounts[static_cast<std::size_t>(K)];
    Advance();
    return BuilderId<K>::Make(static_cast<std::uint32_t>(index), 0, m_Token.GetValue());
}

/// Tombstones the slot: dead, its generation bumped so that every handle to
/// it is stale, never reused by this builder (ADR-0004).
template <EntityKind K>
void ShapeBuilder::Release(const std::uint32_t index)
{
    m_Tables.template Table<K>().Modify(index, [](Detail::Record<K>& r) { r.Slot.Kill(); });
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
    return CurveRef::FromIndex(static_cast<std::uint32_t>(
        AppendGeometry(m_Tables.Curves, m_ReservedGeometry[0], std::move(curve))));
}

Result<PCurveRef> ShapeBuilder::AddPCurve(Curve2d pcurve)
{
    STEREON_TRY_VOID(RequireToken("AddPCurve"));
    if (pcurve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddPCurve: the pcurve is empty");
    }
    return PCurveRef::FromIndex(static_cast<std::uint32_t>(
        AppendGeometry(m_Tables.PCurves, m_ReservedGeometry[1], std::move(pcurve))));
}

Result<SurfaceRef> ShapeBuilder::AddSurface(Surface surface)
{
    STEREON_TRY_VOID(RequireToken("AddSurface"));
    if (surface.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "AddSurface: the surface is empty");
    }
    return SurfaceRef::FromIndex(static_cast<std::uint32_t>(
        AppendGeometry(m_Tables.Surfaces, m_ReservedGeometry[2], std::move(surface))));
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
    // An edge runs between vertices at its range's ends, so both are finite.
    if (!spec.Range.IsBounded()) {
        return Fail(ErrorCode::InvalidInput, "AddEdge: a parameter range that is not bounded");
    }
    if (HasOwnCurve(spec.Kind) != spec.Curve.IsValid()) {
        return Fail(
            ErrorCode::InvalidInput,
            std::string("AddEdge: a ") + std::string(ToString(spec.Kind)) + " edge " +
                (HasOwnCurve(spec.Kind) ? "needs its 3D curve" : "has no 3D curve of its own"));
    }
    if (spec.Curve.IsValid() && !Detail::InPool(m_Tables.Curves, spec.Curve)) {
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
    if (!Detail::InPool(m_Tables.Surfaces, surface)) {
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
    if (pcurve.IsValid() && !Detail::InPool(m_Tables.PCurves, pcurve)) {
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
    // A new use brings a new pcurve of the edge (ADR-0004).
    DropCache(edge.GetIndex());
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
    // The cache was the old primary side's S(c(t)).
    DropCache(edge.GetIndex());
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
    if (!Detail::InPool(m_Tables.Surfaces, surface)) {
        return Fail(ErrorCode::InvalidInput,
                    "SetSurface: the surface is not in this builder's pool");
    }
    m_Tables.Faces.Modify(face.GetIndex(), [&](Detail::FaceRecord& f) { f.Surface = surface; });
    // Every edge the face uses evaluated through the old surface on this
    // side; its cache was S(c(t)) on it.
    for (std::uint32_t l = m_Tables.Faces[face.GetIndex()].FirstLoop.GetIndex(), loops = 0;
         l != LoopId::Invalid && loops < m_Tables.Loops.GetSize();
         l = m_Tables.Loops[l].NextInFace.GetIndex(), ++loops) {
        const std::uint32_t anchor = m_Tables.Loops[l].Anchor.GetIndex();
        for (std::uint32_t c = anchor, steps = 0;
             c != CoedgeId::Invalid && steps < m_Tables.Coedges.GetSize(); ++steps) {
            DropCache(m_Tables.Coedges[c].Edge.GetIndex());
            c = m_Tables.Coedges[c].Next.GetIndex();
            if (c == anchor) {
                break;
            }
        }
    }
    Advance();
    return {};
}

void ShapeBuilder::DropCache(const std::uint32_t edge)
{
    const Detail::EdgeRecord& e = m_Tables.Edges[edge];
    if (!e.Cached.IsValid() && e.CacheStatus == CacheCertificate::None) {
        return; // nothing to drop: no page is copied
    }
    m_Tables.Edges.Modify(edge, [](Detail::EdgeRecord& r) {
        r.Cached = {};
        r.CacheStatus = CacheCertificate::None;
        r.CacheBound = std::numeric_limits<double>::infinity();
    });
}

template <EntityKind K>
Result<void> ShapeBuilder::SetToleranceOf(const BuilderId<K> h, const double tolerance,
                                          const std::string_view reason, const char* operation)
{
    STEREON_TRY_VOID(Require(h, operation));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the tolerance", operation));
    const double old = m_Tables.template Table<K>()[h.GetIndex()].Tolerance;
    m_Tables.template Table<K>().Modify(h.GetIndex(),
                                        [&](Detail::Record<K>& r) { r.Tolerance = tolerance; });
    ToleranceChange change;
    change.Kind = K;
    change.Index = h.GetIndex();
    change.Old = old;
    change.New = tolerance;
    change.Reason = std::string(reason);
    m_ToleranceChanges.push_back(std::move(change));
    Advance();
    return {};
}

Result<void> ShapeBuilder::SetTolerance(const BuilderVertexId vertex, const double tolerance,
                                        const std::string_view reason)
{
    return SetToleranceOf(vertex, tolerance, reason, "SetTolerance");
}

Result<void> ShapeBuilder::SetTolerance(const BuilderEdgeId edge, const double tolerance,
                                        const std::string_view reason)
{
    return SetToleranceOf(edge, tolerance, reason, "SetTolerance");
}

namespace {

/// True for a periodic direction: Period() and Periods() give exactly 0 for
/// one that is not, a marker rather than a tolerance.
[[nodiscard]] bool IsPeriodic(const double period)
{
    return period > 0.0; // STEREON-LINT-ALLOW(epsilon): 0 marks no period
}

/// True when the parameter lies in the domain in every direction that is
/// not periodic: where a pcurve may evaluate a surface.
[[nodiscard]] bool OnSurface(const Surface& surface, const Point2& uv)
{
    const Box2 domain = Domain(surface);
    const Vec2 periods = Periods(surface);
    return (IsPeriodic(periods.X) || Contains(domain.X, uv.X)) &&
           (IsPeriodic(periods.Y) || Contains(domain.Y, uv.Y));
}

} // namespace

Result<CurveRef> ShapeBuilder::CacheEdgeCurve(const BuilderEdgeId edge)
{
    STEREON_TRY_VOID(Require(edge, "CacheEdgeCurve"));
    const Detail::EdgeRecord& e = m_Tables.Edges[edge.GetIndex()];
    if (!HasCachedCurve(e.Kind)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string("CacheEdgeCurve: a ") + std::string(ToString(e.Kind)) +
                        " edge keeps no cached curve; only Intersection and Trim edges do");
    }
    if (!e.Primary.IsValid()) {
        return Fail(ErrorCode::InvalidInput,
                    "CacheEdgeCurve: the edge has no primary coedge to evaluate through");
    }
    const Detail::CoedgeRecord& use = m_Tables.Coedges[e.Primary.GetIndex()];
    const double low = e.Range.GetLower();
    const double high = e.Range.GetUpper();
    if (!use.PCurve.IsValid() || !std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        return Fail(
            ErrorCode::InvalidInput,
            "CacheEdgeCurve: the primary coedge needs a pcurve and the edge a finite range");
    }
    const Curve2d& pcurve = m_Tables.PCurves[use.PCurve.GetIndex()];
    const Surface& surface =
        m_Tables.Surfaces[m_Tables.Faces[m_Tables.Loops[use.Loop.GetIndex()].Face.GetIndex()]
                              .Surface.GetIndex()];
    if (!IsPeriodic(Period(pcurve)) && !Contains(Domain(pcurve), e.Range)) {
        return Fail(ErrorCode::InvalidInput,
                    "CacheEdgeCurve: the edge's range is not in the primary pcurve's domain");
    }
    // S(c(t)) at evenly spread parameters; the last one exactly the end.
    std::vector<Point3> points;
    points.reserve(CacheSamples);
    const double step = (high - low) / static_cast<double>(CacheSamples - 1);
    for (std::size_t k = 0; k < CacheSamples; ++k) {
        const double t = k + 1 == CacheSamples ? high : low + step * static_cast<double>(k);
        const Point2 uv = Eval(pcurve, t);
        if (!OnSurface(surface, uv)) {
            return Fail(ErrorCode::InvalidInput,
                        "CacheEdgeCurve: the primary pcurve leaves its surface's domain");
        }
        points.push_back(Eval(surface, uv));
    }
    // Uniform parametrisation puts sample k at k / (n - 1) of [0, 1], which
    // the reparametrisation maps to the edge's own parameter.
    STEREON_TRY(const CurveFit fit, Interpolate(points, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve spline, Reparametrized(fit.Curve, e.Range));
    Curve cache(std::move(spline));
    const CacheBound certificate = CertifyCache(cache, surface, pcurve, e.Range, e.Tolerance);
    STEREON_TRY(const CurveRef ref, AddCurve(std::move(cache)));
    m_Tables.Edges.Modify(edge.GetIndex(), [&](Detail::EdgeRecord& r) {
        r.Cached = ref;
        r.CacheStatus = certificate.Status;
        r.CacheBound = certificate.Bound;
    });
    return ref;
}

/// ADR-0009 section 4: the cache's deviation from S(c(t)) over the whole
/// range, proved within the edge's tolerance where interval evaluation
/// covers the types; otherwise the cache stays uncertified. An edge without
/// a tolerance has nothing to certify within.
ShapeBuilder::CacheBound ShapeBuilder::CertifyCache(const Curve& cache, const Surface& surface,
                                                    const Curve2d& pcurve, const Interval& range,
                                                    const double tolerance)
{
    CacheBound result{CacheCertificate::Uncertified, std::numeric_limits<double>::infinity()};
    if (tolerance > 0.0) { // STEREON-LINT-ALLOW(epsilon): a sign test
        if (const Result<DeviationCertificate> certificate =
                CertifyDeviation(cache, surface, pcurve, range, tolerance)) {
            result.Status = CacheCertificate::Certified;
            result.Bound = certificate->Bound;
        }
    }
    return result;
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
    // The use's pcurve leaves the edge, and its primary may have passed on.
    DropCache(r.Edge.GetIndex());
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

// The member templates the Euler operators (Euler.cpp) use, instantiated
// here where they are defined.
template Result<void> ShapeBuilder::Require<EntityKind::Vertex>(BuilderVertexId, const char*) const;
template Result<void> ShapeBuilder::Require<EntityKind::Edge>(BuilderEdgeId, const char*) const;
template Result<void> ShapeBuilder::Require<EntityKind::Coedge>(BuilderCoedgeId, const char*) const;
template Result<void> ShapeBuilder::Require<EntityKind::Loop>(BuilderLoopId, const char*) const;
template Result<void> ShapeBuilder::Require<EntityKind::Face>(BuilderFaceId, const char*) const;
template Result<void> ShapeBuilder::Require<EntityKind::Solid>(BuilderSolidId, const char*) const;
template BuilderCoedgeId
    ShapeBuilder::Allocate<EntityKind::Coedge>(Detail::Record<EntityKind::Coedge>);
template void ShapeBuilder::Release<EntityKind::Vertex>(std::uint32_t);
template void ShapeBuilder::Release<EntityKind::Edge>(std::uint32_t);
template void ShapeBuilder::Release<EntityKind::Coedge>(std::uint32_t);
template void ShapeBuilder::Release<EntityKind::Face>(std::uint32_t);
template void ShapeBuilder::Release<EntityKind::Shell>(std::uint32_t);
template void ShapeBuilder::UnlinkChild<EntityKind::Shell, EntityKind::Face>(
    std::uint32_t, Id<EntityKind::Shell> Detail::Record<EntityKind::Face>::*,
    Id<EntityKind::Face> Detail::Record<EntityKind::Face>::*,
    Id<EntityKind::Face> Detail::Record<EntityKind::Shell>::*,
    Id<EntityKind::Face> Detail::Record<EntityKind::Shell>::*);
template void ShapeBuilder::UnlinkChild<EntityKind::Solid, EntityKind::Shell>(
    std::uint32_t, Id<EntityKind::Solid> Detail::Record<EntityKind::Shell>::*,
    Id<EntityKind::Shell> Detail::Record<EntityKind::Shell>::*,
    Id<EntityKind::Shell> Detail::Record<EntityKind::Solid>::*,
    Id<EntityKind::Shell> Detail::Record<EntityKind::Solid>::*);

ShapeBuilder::ShapeBuilder(const ShapeBuilder& other, CandidateKey)
    : m_Tables(other.m_Tables)
    , m_Token(other.m_Token)
    , m_ReservedEntities(other.m_ReservedEntities)
    , m_ReservedGeometry(other.m_ReservedGeometry)
    , m_Epoch(other.m_Epoch)
    , m_SourceStore(other.m_SourceStore)
    , m_Source(other.m_Source)
    , m_ToleranceChanges(other.m_ToleranceChanges)
    , m_Events(other.m_Events)
{}

void ShapeBuilder::Abandon(const ShapeBuilder& candidate) noexcept
{
    for (std::size_t i = 0; i < m_ReservedEntities.size(); ++i) {
        m_ReservedEntities[i] = std::max(m_ReservedEntities[i], candidate.m_ReservedEntities[i]);
    }
    for (std::size_t i = 0; i < m_ReservedGeometry.size(); ++i) {
        m_ReservedGeometry[i] = std::max(m_ReservedGeometry[i], candidate.m_ReservedGeometry[i]);
    }
    Advance();
}

std::vector<BuilderCoedgeId> ShapeBuilder::CoedgesOf(const BuilderLoopId loop) const
{
    std::vector<BuilderCoedgeId> coedges;
    if (!IsLive(loop)) {
        return coedges;
    }
    const std::uint32_t anchor = m_Tables.Loops[loop.GetIndex()].Anchor.GetIndex();
    for (std::uint32_t c = anchor;
         c != CoedgeId::Invalid && coedges.size() < m_Tables.Coedges.GetSize();) {
        coedges.push_back(
            BuilderCoedgeId::Make(c, m_Tables.Coedges[c].Slot.GetGeneration(), m_Token.GetValue()));
        c = m_Tables.Coedges[c].Next.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    return coedges;
}

template <EntityKind K>
BuilderId<K> ShapeBuilder::HandleOfSlot(const std::uint32_t index) const noexcept
{
    if (index == Id<K>::Invalid || !m_Tables.template IsLive<K>(index)) {
        return {};
    }
    return BuilderId<K>::Make(index, m_Tables.template Table<K>()[index].Slot.GetGeneration(),
                              m_Token.GetValue());
}

BuilderEdgeId ShapeBuilder::EdgeOf(const BuilderCoedgeId coedge) const noexcept
{
    if (!IsLive(coedge)) {
        return {};
    }
    return HandleOfSlot<EntityKind::Edge>(m_Tables.Coedges[coedge.GetIndex()].Edge.GetIndex());
}

BuilderLoopId ShapeBuilder::LoopOf(const BuilderCoedgeId coedge) const noexcept
{
    if (!IsLive(coedge)) {
        return {};
    }
    return HandleOfSlot<EntityKind::Loop>(m_Tables.Coedges[coedge.GetIndex()].Loop.GetIndex());
}

EdgeUseSense ShapeBuilder::SenseOf(const BuilderCoedgeId coedge) const noexcept
{
    return IsLive(coedge) ? m_Tables.Coedges[coedge.GetIndex()].Sense : EdgeUseSense::Same;
}

BuilderVertexId ShapeBuilder::StartOf(const BuilderEdgeId edge) const noexcept
{
    if (!IsLive(edge)) {
        return {};
    }
    return HandleOfSlot<EntityKind::Vertex>(m_Tables.Edges[edge.GetIndex()].Start.GetIndex());
}

BuilderVertexId ShapeBuilder::EndOf(const BuilderEdgeId edge) const noexcept
{
    if (!IsLive(edge)) {
        return {};
    }
    return HandleOfSlot<EntityKind::Vertex>(m_Tables.Edges[edge.GetIndex()].End.GetIndex());
}

Point3 ShapeBuilder::PositionOf(const BuilderVertexId vertex) const noexcept
{
    return IsLive(vertex) ? m_Tables.Vertices[vertex.GetIndex()].Position : Point3{};
}

std::vector<BuilderLoopId> ShapeBuilder::LoopsOf(const BuilderFaceId face) const
{
    std::vector<BuilderLoopId> loops;
    if (!IsLive(face)) {
        return loops;
    }
    for (std::uint32_t l = m_Tables.Faces[face.GetIndex()].FirstLoop.GetIndex();
         l != LoopId::Invalid && loops.size() < m_Tables.Loops.GetSize();
         l = m_Tables.Loops[l].NextInFace.GetIndex()) {
        loops.push_back(HandleOfSlot<EntityKind::Loop>(l));
    }
    return loops;
}

ValidityReport ShapeBuilder::Validate() const
{
    return Detail::ValidateStructure(m_Tables);
}

Result<Shape> Freeze(ShapeBuilder&& builder)
{
    STEREON_PROFILE_SCOPE("topo.Freeze");
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

ValidityReport ValidateTopology(const Shape& shape, const ValidityLevel level)
{
    // An empty Shape, default-constructed or moved from, has no tables: it is
    // trivially valid at every level, as it is at the Structural one.
    ValidityReport report = ValidateTopology(shape);
    if (level == ValidityLevel::Structural || shape.IsEmpty() || !report.IsValid()) {
        return report;
    }
    return Detail::ValidateGeometry(shape.GetSnapshot()->GetTables());
}

// --- Shape-producing operations (ADR-0006) ------------------------------------

namespace {

[[nodiscard]] bool IsLiveOf(const Detail::TopologyTables& t, const EntityKind kind,
                            const std::uint32_t index)
{
    return kind == EntityKind::Vertex ? t.IsLive<EntityKind::Vertex>(index)
                                      : t.IsLive<EntityKind::Edge>(index);
}

/// The builder's SetTolerance() calls as one change per entity still live in
/// the result, ordered by kind and index: the first old value, the last new
/// one and every reason in order; none where the value came back to where it
/// started.
[[nodiscard]] std::vector<ToleranceChange> Coalesced(std::vector<ToleranceChange> calls,
                                                     const Detail::TopologyTables& result)
{
    std::vector<ToleranceChange> changes;
    for (ToleranceChange& call : calls) {
        const auto same = std::ranges::find_if(changes, [&](const ToleranceChange& c) {
            return c.Kind == call.Kind && c.Index == call.Index;
        });
        if (same == changes.end()) {
            changes.push_back(std::move(call));
        } else {
            same->New = call.New;
            same->Reason += "; " + call.Reason;
        }
    }
    std::erase_if(changes, [&](const ToleranceChange& c) {
        return c.New == c.Old || !IsLiveOf(result, c.Kind, c.Index);
    });
    std::ranges::sort(changes, [](const ToleranceChange& a, const ToleranceChange& b) {
        return std::pair(static_cast<int>(a.Kind), a.Index) <
               std::pair(static_cast<int>(b.Kind), b.Index);
    });
    return changes;
}

[[nodiscard]] std::string Named(const Shape& shape, const EntityKind kind,
                                const std::uint32_t index)
{
    return Detail::DescribeEntity(shape.GetSnapshot()->GetTables(), kind, index);
}

[[nodiscard]] Diagnostic MakeWarning(const DiagnosticCode code, std::string message,
                                     std::string entity)
{
    Diagnostic warning;
    warning.Code = code;
    warning.Message = std::move(message);
    warning.Entities.push_back(std::move(entity));
    return warning;
}

/// ToleranceExceeded for a change the Context does not allow (ADR-0001,
/// ADR-0006): above the ceiling, or any growth when growth fails.
[[nodiscard]] Result<void> CheckTolerances(const std::vector<ToleranceChange>& changes,
                                           const Context& context, const Shape& shape)
{
    for (const ToleranceChange& c : changes) {
        if (c.New > context.GetToleranceCeiling()) {
            return Fail(ErrorCode::ToleranceExceeded,
                        std::format("Commit: the tolerance of {} would be {}, above the "
                                    "Context's ceiling of {}",
                                    Named(shape, c.Kind, c.Index), c.New,
                                    context.GetToleranceCeiling()));
        }
        if (c.New > c.Old && context.GetToleranceGrowthPolicy() == ToleranceGrowth::Fail) {
            return Fail(ErrorCode::ToleranceExceeded,
                        std::format("Commit: the tolerance of {} grew from {} to {} ({}), and the "
                                    "Context fails any tolerance growth",
                                    Named(shape, c.Kind, c.Index), c.Old, c.New, c.Reason));
        }
    }
    return {};
}

/// ADR-0001's bounds on every live vertex and edge of a result, whatever its
/// change log says: an entity made with its tolerance, or kept from a shape
/// published under another Context, has no SetTolerance() entry to check.
/// ToleranceExceeded above the Context's ceiling; InvalidInput below its
/// resolution, which every local tolerance is at least.
[[nodiscard]] Result<void> CheckToleranceBounds(const Detail::TopologyTables& t,
                                                const Context& context, const Shape& shape,
                                                const char* operation)
{
    const auto check = [&](const EntityKind kind, const std::uint32_t index,
                           const double tolerance) -> Result<void> {
        if (tolerance > context.GetToleranceCeiling()) {
            return Fail(ErrorCode::ToleranceExceeded,
                        std::format("{}: the tolerance of {} is {}, above the Context's ceiling "
                                    "of {}",
                                    operation, Named(shape, kind, index), tolerance,
                                    context.GetToleranceCeiling()));
        }
        if (tolerance < context.GetResolution()) {
            return Fail(ErrorCode::InvalidInput,
                        std::format("{}: the tolerance of {} is {}, below the Context's "
                                    "resolution of {}",
                                    operation, Named(shape, kind, index), tolerance,
                                    context.GetResolution()));
        }
        return {};
    };
    for (std::uint32_t v = 0; v < t.Vertices.GetSize(); ++v) {
        if (t.IsLive<EntityKind::Vertex>(v)) {
            STEREON_TRY_VOID(check(EntityKind::Vertex, v, t.Vertices[v].Tolerance));
        }
    }
    for (std::uint32_t e = 0; e < t.Edges.GetSize(); ++e) {
        if (t.IsLive<EntityKind::Edge>(e)) {
            STEREON_TRY_VOID(check(EntityKind::Edge, e, t.Edges[e].Tolerance));
        }
    }
    return {};
}

/// Whether a publication checks the Geometric level: when the Context asks
/// for it, and always in checked builds, which check every level (ADR-0018).
[[nodiscard]] bool ChecksGeometry(const Context& context) noexcept
{
    return STEREON_CHECKED_HANDLES || context.GetValidityLevel() != ValidityLevel::Structural;
}

/// The warnings of a commit: tolerances that grew, cached edge curves the
/// edit dropped.
[[nodiscard]] std::vector<Diagnostic> WarningsOf(const std::vector<ToleranceChange>& changes,
                                                 const Detail::TopologyTables* before,
                                                 const Shape& shape)
{
    std::vector<Diagnostic> warnings;
    for (const ToleranceChange& c : changes) {
        if (c.New > c.Old) {
            warnings.push_back(
                MakeWarning(DiagnosticCode::ToleranceGrown,
                            std::format("tolerance grew from {} to {}: {}", c.Old, c.New, c.Reason),
                            Named(shape, c.Kind, c.Index)));
        }
    }
    if (before != nullptr) {
        const Detail::TopologyTables& after = shape.GetSnapshot()->GetTables();
        for (std::uint32_t e = 0; e < before->Edges.GetSize(); ++e) {
            if (before->IsLive<EntityKind::Edge>(e) && after.IsLive<EntityKind::Edge>(e) &&
                before->Edges[e].Cached.IsValid() && !after.Edges[e].Cached.IsValid()) {
                warnings.push_back(
                    MakeWarning(DiagnosticCode::CacheDropped,
                                "the cached 3D curve was dropped: a pcurve, the primary "
                                "coedge or a surface of the edge changed",
                                Named(shape, EntityKind::Edge, e)));
            }
        }
    }
    return warnings;
}

} // namespace

Result<OperationResult> Commit(ShapeBuilder&& builder, const Context& context)
{
    STEREON_PROFILE_SCOPE("topo.Commit");
    STEREON_TRY_VOID(builder.RequireToken("Commit"));
    // Consumed whatever happens, as Freeze() consumes it (ADR-0027).
    ShapeBuilder consumed(std::move(builder));
    const std::shared_ptr<const TopologySnapshot> source = consumed.m_Source;
    std::vector<ToleranceChange> calls = std::move(consumed.m_ToleranceChanges);
    std::vector<TopologyEvent> events = std::move(consumed.m_Events);
    if constexpr (STEREON_CHECKED_HANDLES) {
        if (source != nullptr) {
            const ValidityReport report = Detail::ValidateStructure(source->GetTables());
            if (!report.IsValid()) {
                return Fail(ErrorCode::InvalidInput,
                            "Commit: the shape the builder was opened on is not structurally "
                            "valid:\n" +
                                report.ToText());
            }
        }
    }
    Result<Shape> frozen = Freeze(std::move(consumed));
    if (!frozen) {
        return Fail(frozen.error().GetCode(),
                    "Commit: " + std::string(frozen.error().GetMessage()));
    }
    const Detail::TopologyTables& result = frozen->GetSnapshot()->GetTables();
    const Detail::TopologyTables* before = source != nullptr ? &source->GetTables() : nullptr;
    std::vector<ToleranceChange> changes = Coalesced(std::move(calls), result);
    STEREON_TRY_VOID(CheckTolerances(changes, context, *frozen));
    STEREON_TRY_VOID(CheckToleranceBounds(result, context, *frozen, "Commit"));
    // A result that fails the Geometric level is a kernel bug or a caller's
    // inconsistent geometry; either way it is not published (ADR-0006).
    const bool geometric = ChecksGeometry(context);
    if (geometric) {
        const ValidityReport report = Detail::ValidateGeometry(result);
        if (!report.IsValid()) {
            return Fail(ErrorCode::InvalidResult,
                        "Commit: the result is not geometrically valid:\n" + report.ToText());
        }
    }
    History history = Detail::HistoryBetween(before, result);
    // ADR-0030: missing history fails the operation before publication.
    if (!Detail::HistoryIsComplete(history, before, result)) {
        return Fail(ErrorCode::InvalidResult,
                    "Commit: the History does not account for every entity (a kernel bug)");
    }
    OperationResult published;
    published.TopologyEvents = std::move(events);
    published.Warnings = WarningsOf(changes, before, *frozen);
    published.ToleranceChanges = std::move(changes);
    published.History = std::move(history);
    published.Shape = *std::move(frozen);
    published.Checked = geometric ? ValidityLevel::Geometric : ValidityLevel::Structural;
    return published;
}

Result<OperationResult> CompactShape(const Shape& shape, const Context& context)
{
    // Compaction renumbers and changes no geometry or tolerance, but it
    // publishes under the Context's policy as Commit() does (ADR-0018): the
    // input carries no evidence of a Geometric check, since Freeze()
    // certifies structure only.
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "CompactShape: the Shape is empty");
    }
    if constexpr (STEREON_CHECKED_HANDLES) {
        const ValidityReport report = ValidateTopology(shape);
        if (!report.IsValid()) {
            return Fail(ErrorCode::InvalidInput,
                        "CompactShape: the Shape is not structurally valid:\n" + report.ToText());
        }
    }
    const Detail::TopologyTables& source = shape.GetSnapshot()->GetTables();
    STEREON_TRY_VOID(CheckToleranceBounds(source, context, shape, "CompactShape"));
    const bool geometric = ChecksGeometry(context);
    if (geometric) {
        const ValidityReport report = Detail::ValidateGeometry(source);
        if (!report.IsValid()) {
            return Fail(ErrorCode::InvalidInput,
                        "CompactShape: the Shape is not geometrically valid:\n" + report.ToText());
        }
    }
    STEREON_TRY(Compacted compacted, Compact(shape));
    std::vector<HistoryEntry> entries;
    const auto continuations = [&]<EntityKind K>() {
        for (std::uint32_t old = 0; old < compacted.Remap.template OldCount<K>(); ++old) {
            const std::uint32_t index = compacted.Remap.template NewIndex<K>(old);
            if (index != Id<K>::Invalid) {
                HistoryEntry entry;
                entry.Kind = K;
                entry.Relation = HistoryRelation::Unchanged;
                entry.Source = old;
                entry.Target = index;
                entries.push_back(entry);
            }
        }
    };
    continuations.template operator()<EntityKind::Vertex>();
    continuations.template operator()<EntityKind::Edge>();
    continuations.template operator()<EntityKind::Coedge>();
    continuations.template operator()<EntityKind::Loop>();
    continuations.template operator()<EntityKind::Face>();
    continuations.template operator()<EntityKind::Shell>();
    continuations.template operator()<EntityKind::Solid>();
    History history(std::move(entries));
    const Detail::TopologyTables& result = compacted.Shape.GetSnapshot()->GetTables();
    if (!Detail::HistoryIsComplete(history, &source, result)) {
        return Fail(ErrorCode::InvalidResult,
                    "CompactShape: the History does not account for every entity (a kernel bug)");
    }
    OperationResult published;
    published.History = std::move(history);
    published.Shape = std::move(compacted.Shape);
    published.Checked = geometric ? ValidityLevel::Geometric : ValidityLevel::Structural;
    return published;
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
    b.m_Tables.Coedges.Modify(c.GetIndex(), [](Detail::CoedgeRecord& r) { r.Slot.Kill(); });
    --b.m_Tables.LiveCounts[static_cast<std::size_t>(EntityKind::Coedge)];
}

void TopologyFaults::SetCache(ShapeBuilder& b, const BuilderEdgeId e, const CurveRef cached,
                              const CacheCertificate status, const double bound)
{
    b.m_Tables.Edges.Modify(e.GetIndex(), [&](Detail::EdgeRecord& r) {
        r.Cached = cached;
        r.CacheStatus = status;
        r.CacheBound = bound;
    });
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
