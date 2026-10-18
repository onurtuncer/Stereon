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

ShapeBuilder::ShapeBuilder(ShapeBuilder&& other) noexcept
    : m_Tables(std::move(other.m_Tables))
    , m_Token(std::exchange(other.m_Token, BuilderToken{}))
    , m_Epoch(other.m_Epoch)
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

CurveRef ShapeBuilder::AddCurve(Curve curve)
{
    return CurveRef::FromIndex(
        static_cast<std::uint32_t>(m_Tables.Curves.PushBack(std::move(curve))));
}

PCurveRef ShapeBuilder::AddPCurve(Curve2d pcurve)
{
    return PCurveRef::FromIndex(
        static_cast<std::uint32_t>(m_Tables.PCurves.PushBack(std::move(pcurve))));
}

SurfaceRef ShapeBuilder::AddSurface(Surface surface)
{
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
    return Allocate<EntityKind::Edge>(record);
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
    const ValidationCertificate certificate{
        .ValidatorVersion = Detail::StructuralValidatorVersion,
        .Invariants = "Structural S01-S15",
        .Snapshot = Detail::NextStoreId(),
    };
    auto snapshot =
        std::make_shared<const TopologySnapshot>(Detail::FreezeKey{}, std::move(consumed.m_Tables),
                                                 certificate, consumed.m_Token.GetValue());
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

} // namespace Stereon::Testing
