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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

module Stereon.Topo;

// The Euler operators of ADR-0027 (the Boolean prototype's M5): MVFS/KVFS,
// MEV/KEV and MEF/KEF. Each checks every precondition first and only then
// writes, so a refusal changes nothing. A loop's cycle is rebuilt whole from
// the list of its coedges (RelinkCycle), which keeps the splices of MEV, MEF
// and their inverses one operation rather than a case for each neighbour.

namespace Stereon {

namespace {

constexpr std::uint32_t None = CoedgeId::Invalid;

[[nodiscard]] bool NeedsPrimary(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Intersection || kind == EdgeCurveKind::Trim ||
           kind == EdgeCurveKind::Boundary || kind == EdgeCurveKind::Degenerate;
}

[[nodiscard]] bool HasOwnCurve(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Exact || kind == EdgeCurveKind::Spatial;
}

/// The refusal of an operator: InvalidInput, naming it. Converts to any
/// operator's Result.
[[nodiscard]] auto Refuse(const char* operation, const std::string& why)
{
    return Fail(ErrorCode::InvalidInput, std::string(operation) + ": " + why);
}

} // namespace

// --- Shared steps -------------------------------------------------------------

Result<void> ShapeBuilder::CheckEdgeGeometry(const EdgeGeometry& geometry,
                                             const char* operation) const
{
    STEREON_TRY_VOID(Detail::CheckNonNegative(geometry.Tolerance, "the tolerance", operation));
    if (geometry.Range.IsEmpty() || !geometry.Range.IsBounded()) {
        return Refuse(operation, "the edge's parameter range is empty or not bounded");
    }
    if (HasOwnCurve(geometry.Kind) != geometry.Curve.IsValid()) {
        return Refuse(
            operation,
            std::string("a ") + std::string(ToString(geometry.Kind)) + " edge " +
                (HasOwnCurve(geometry.Kind) ? "needs its 3D curve" : "has no 3D curve of its own"));
    }
    if (geometry.Curve.IsValid() && geometry.Curve.GetIndex() >= m_Tables.Curves.GetSize()) {
        return Refuse(operation, "the curve is not in this builder's pool");
    }
    if (!geometry.PCurve.IsValid() && geometry.Kind != EdgeCurveKind::Spatial) {
        return Refuse(operation, std::string("the uses of a ") +
                                     std::string(ToString(geometry.Kind)) +
                                     " edge need their pcurve");
    }
    if (geometry.PCurve.IsValid() && geometry.PCurve.GetIndex() >= m_Tables.PCurves.GetSize()) {
        return Refuse(operation, "the pcurve is not in this builder's pool");
    }
    return {};
}

/// The vertex at which a coedge ends: its edge's end along the edge, its
/// start against it.
std::uint32_t ShapeBuilder::EndVertexOf(const std::uint32_t coedge) const noexcept
{
    const Detail::CoedgeRecord& c = m_Tables.Coedges[coedge];
    const Detail::EdgeRecord& e = m_Tables.Edges[c.Edge.GetIndex()];
    return (c.Sense == EdgeUseSense::Same ? e.End : e.Start).GetIndex();
}

/// The coedges of a loop's cycle, starting at `coedge`.
std::vector<std::uint32_t> ShapeBuilder::CycleFrom(const std::uint32_t coedge) const
{
    std::vector<std::uint32_t> cycle;
    for (std::uint32_t c = coedge; c != None && cycle.size() < m_Tables.Coedges.GetSize();) {
        cycle.push_back(c);
        c = m_Tables.Coedges[c].Next.GetIndex();
        if (c == coedge) {
            break;
        }
    }
    return cycle;
}

/// Makes `cycle` the loop's cycle, in order: links, membership and anchor.
/// The anchor stays where it was if it is still in the cycle.
void ShapeBuilder::RelinkCycle(const std::uint32_t loop, const std::vector<std::uint32_t>& cycle)
{
    const std::size_t n = cycle.size();
    for (std::size_t i = 0; i < n; ++i) {
        const CoedgeId next = CoedgeId::FromIndex(cycle[(i + 1) % n]);
        const CoedgeId prev = CoedgeId::FromIndex(cycle[(i + n - 1) % n]);
        m_Tables.Coedges.Modify(cycle[i], [&](Detail::CoedgeRecord& c) {
            c.Next = next;
            c.Prev = prev;
            c.Loop = LoopId::FromIndex(loop);
        });
    }
    const std::uint32_t anchor = m_Tables.Loops[loop].Anchor.GetIndex();
    const bool kept = std::ranges::find(cycle, anchor) != cycle.end();
    if (!kept) {
        const CoedgeId first = n == 0 ? CoedgeId{} : CoedgeId::FromIndex(cycle.front());
        m_Tables.Loops.Modify(loop, [&](Detail::LoopRecord& l) { l.Anchor = first; });
    }
    Advance();
}

/// A new use of `edge` in `loop`, in the edge's radial ring after its
/// anchor, and its primary if the edge's kind needs one and has none. Its
/// cycle links are the caller's, through RelinkCycle.
BuilderCoedgeId ShapeBuilder::NewUse(const std::uint32_t edge, const std::uint32_t loop,
                                     const EdgeUseSense sense, const PCurveRef pcurve)
{
    Detail::CoedgeRecord record;
    record.Edge = EdgeId::FromIndex(edge);
    record.Loop = LoopId::FromIndex(loop);
    record.PCurve = pcurve;
    record.Sense = sense;
    const BuilderCoedgeId use = Allocate<EntityKind::Coedge>(record);
    const CoedgeId self = CoedgeId::FromIndex(use.GetIndex());
    const CoedgeId anchor = m_Tables.Edges[edge].RadialAnchor;
    if (!anchor.IsValid()) {
        m_Tables.Coedges.Modify(use.GetIndex(), [&](Detail::CoedgeRecord& c) { c.Radial = self; });
        m_Tables.Edges.Modify(edge, [&](Detail::EdgeRecord& e) { e.RadialAnchor = self; });
    } else {
        const CoedgeId after = m_Tables.Coedges[anchor.GetIndex()].Radial;
        m_Tables.Coedges.Modify(use.GetIndex(), [&](Detail::CoedgeRecord& c) { c.Radial = after; });
        m_Tables.Coedges.Modify(anchor.GetIndex(),
                                [&](Detail::CoedgeRecord& c) { c.Radial = self; });
    }
    m_Tables.Edges.Modify(edge, [&](Detail::EdgeRecord& e) {
        if (NeedsPrimary(e.Kind) && !e.Primary.IsValid()) {
            e.Primary = self;
        }
    });
    return use;
}

/// Tombstones an edge and every use in its ring, and takes it out of its
/// vertices' lists. The uses must already be out of their loops' cycles.
void ShapeBuilder::KillEdgeAndUses(const std::uint32_t edge)
{
    const std::uint32_t anchor = m_Tables.Edges[edge].RadialAnchor.GetIndex();
    std::vector<std::uint32_t> uses;
    for (std::uint32_t c = anchor; c != None && uses.size() < m_Tables.Coedges.GetSize();) {
        uses.push_back(c);
        c = m_Tables.Coedges[c].Radial.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    for (const std::uint32_t c : uses) {
        Release<EntityKind::Coedge>(c);
    }
    m_Tables.Edges.Modify(edge, [](Detail::EdgeRecord& e) {
        e.RadialAnchor = {};
        e.Primary = {};
    });
    const std::uint32_t start = m_Tables.Edges[edge].Start.GetIndex();
    const std::uint32_t end = m_Tables.Edges[edge].End.GetIndex();
    UnlinkAtVertex(edge, start);
    if (end != start) {
        UnlinkAtVertex(edge, end);
    }
    Release<EntityKind::Edge>(edge);
}

namespace {

/// True when edge e is in vertex v's list.
[[nodiscard]] bool Listed(const Detail::TopologyTables& t, const std::uint32_t e,
                          const std::uint32_t v)
{
    for (std::uint32_t x = t.Vertices[v].FirstEdge.GetIndex(), steps = 0;
         x != EdgeId::Invalid && steps < t.Edges.GetSize(); ++steps) {
        if (x == e) {
            return true;
        }
        const Detail::EdgeRecord& r = t.Edges[x];
        x = (r.Start.GetIndex() == v ? r.NextAtStart : r.NextAtEnd).GetIndex();
    }
    return false;
}

/// Loop l's cycle closes over live coedges of l with reciprocal links; adds
/// the edges it uses to `edges`.
void CheckLoop(const Detail::TopologyTables& t, const std::uint32_t l,
               std::vector<std::uint32_t>& edges)
{
    const std::uint32_t anchor = t.Loops[l].Anchor.GetIndex();
    std::size_t steps = 0;
    for (std::uint32_t c = anchor; c != None && steps <= t.Coedges.GetSize(); ++steps) {
        const Detail::CoedgeRecord& r = t.Coedges[c];
        STEREON_ASSERT(t.IsLive<EntityKind::Coedge>(c));
        STEREON_ASSERT(r.Loop.GetIndex() == l);
        STEREON_ASSERT(t.Coedges[r.Next.GetIndex()].Prev.GetIndex() == c);
        edges.push_back(r.Edge.GetIndex());
        c = r.Next.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    STEREON_ASSERT(steps <= t.Coedges.GetSize());
}

/// Edge e's ring closes over uses of e, and both its vertices list it.
void CheckEdge(const Detail::TopologyTables& t, const std::uint32_t e)
{
    const Detail::EdgeRecord& r = t.Edges[e];
    const std::uint32_t anchor = r.RadialAnchor.GetIndex();
    std::size_t uses = 0;
    for (std::uint32_t c = anchor; c != None && uses <= t.Coedges.GetSize(); ++uses) {
        STEREON_ASSERT(t.Coedges[c].Edge.GetIndex() == e);
        c = t.Coedges[c].Radial.GetIndex();
        if (c == anchor) {
            break;
        }
    }
    STEREON_ASSERT(uses <= t.Coedges.GetSize());
    STEREON_ASSERT(Listed(t, e, r.Start.GetIndex()));
    STEREON_ASSERT(Listed(t, e, r.End.GetIndex()));
}

} // namespace

/// ADR-0027, section 9: in checked builds, the loops and edges an operator
/// touched, and the radial rings and vertex lists of every edge those loops
/// use, are consistent. A failure is a bug in the operator, so a contract.
void ShapeBuilder::CheckNeighbourhood(const std::vector<std::uint32_t>& loops,
                                      const std::vector<std::uint32_t>& edges) const
{
    if constexpr (STEREON_CHECKED_HANDLES) {
        std::vector<std::uint32_t> touched = edges;
        for (const std::uint32_t l : loops) {
            CheckLoop(m_Tables, l, touched);
        }
        for (const std::uint32_t e : touched) {
            CheckEdge(m_Tables, e);
        }
    }
}

/// The edge's second use when it has exactly two, None otherwise.
std::uint32_t ShapeBuilder::SecondOfTwoUses(const std::uint32_t edge) const noexcept
{
    const std::uint32_t first = m_Tables.Edges[edge].RadialAnchor.GetIndex();
    if (first == None) {
        return None;
    }
    const std::uint32_t second = m_Tables.Coedges[first].Radial.GetIndex();
    const bool two = second != first && m_Tables.Coedges[second].Radial.GetIndex() == first;
    return two ? second : None;
}

/// True when `vertex` is an end of an edge that is not closed, and no other
/// edge starts or ends at it: the far end of a strut.
bool ShapeBuilder::IsFreeEnd(const std::uint32_t edge, const std::uint32_t vertex) const noexcept
{
    const Detail::EdgeRecord& e = m_Tables.Edges[edge];
    if (e.Start == e.End) {
        return false;
    }
    const bool atStart = e.Start.GetIndex() == vertex;
    if (!atStart && e.End.GetIndex() != vertex) {
        return false;
    }
    const EdgeId another = atStart ? e.NextAtStart : e.NextAtEnd;
    return m_Tables.Vertices[vertex].FirstEdge.GetIndex() == edge && !another.IsValid();
}

/// The face of a coedge's loop.
std::uint32_t ShapeBuilder::FaceOfUse(const std::uint32_t coedge) const noexcept
{
    return m_Tables.Loops[m_Tables.Coedges[coedge].Loop.GetIndex()].Face.GetIndex();
}

/// KEF's conditions on the face it removes and the face it keeps: one loop,
/// and the same surface, sense and shell.
Result<void> ShapeBuilder::CheckJoinable(const std::uint32_t dead, const std::uint32_t alive,
                                         const char* operation) const
{
    const Detail::FaceRecord& d = m_Tables.Faces[dead];
    const Detail::FaceRecord& a = m_Tables.Faces[alive];
    if (d.FirstLoop != d.LastLoop) {
        return Refuse(operation, "the face to remove must have one loop");
    }
    if (d.Surface != a.Surface || d.Reversed != a.Reversed || d.Shell != a.Shell) {
        return Refuse(operation,
                      "the two faces must lie on one surface, the same way round, in one shell");
    }
    return {};
}

// --- MVFS / KVFS --------------------------------------------------------------

Result<VertexFaceShell> ShapeBuilder::MakeVertexFaceShell(const Point3& position,
                                                          const double tolerance,
                                                          const SurfaceRef surface,
                                                          const bool reversed,
                                                          const BuilderSolidId solid)
{
    constexpr const char* op = "MakeVertexFaceShell";
    STEREON_TRY_VOID(RequireToken(op));
    STEREON_TRY_VOID(Detail::CheckFinite(position, "the position", op));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the tolerance", op));
    if (!surface.IsValid() || surface.GetIndex() >= m_Tables.Surfaces.GetSize()) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeVertexFaceShell: the surface is not in this builder's pool");
    }
    if (solid.IsValid()) {
        STEREON_TRY_VOID(Require(solid, op));
    }
    VertexFaceShell made;
    STEREON_TRY(made.Vertex, AddVertex(position, tolerance));
    STEREON_TRY(made.Face, AddFace(surface, reversed));
    STEREON_TRY(made.Loop, AddLoop(made.Face));
    STEREON_TRY(made.Shell, solid.IsValid() ? AddShell(solid) : AddShell());
    STEREON_TRY_VOID(AddFaceToShell(made.Shell, made.Face));
    return made;
}

Result<void> ShapeBuilder::KillVertexFaceShell(const BuilderVertexId vertex,
                                               const BuilderFaceId face)
{
    constexpr const char* op = "KillVertexFaceShell";
    STEREON_TRY_VOID(Require(vertex, op));
    STEREON_TRY_VOID(Require(face, op));
    if (m_Tables.Vertices[vertex.GetIndex()].FirstEdge.IsValid()) {
        return Refuse(op, "the vertex has edges");
    }
    const Detail::FaceRecord& f = m_Tables.Faces[face.GetIndex()];
    const std::uint32_t loop = f.FirstLoop.GetIndex();
    if (loop == LoopId::Invalid || f.LastLoop.GetIndex() != loop ||
        m_Tables.Loops[loop].Anchor.IsValid()) {
        return Refuse(op, "the face must have one loop, empty");
    }
    const std::uint32_t shell = f.Shell.GetIndex();
    if (shell == ShellId::Invalid ||
        m_Tables.Shells[shell].FirstFace != FaceId::FromIndex(face.GetIndex()) ||
        m_Tables.Shells[shell].LastFace != FaceId::FromIndex(face.GetIndex())) {
        return Refuse(op, "the face must be the only face of its shell");
    }
    RemoveLoopAt(loop);
    UnlinkChild<EntityKind::Shell, EntityKind::Face>(
        face.GetIndex(), &Detail::FaceRecord::Shell, &Detail::FaceRecord::NextInShell,
        &Detail::ShellRecord::FirstFace, &Detail::ShellRecord::LastFace);
    Release<EntityKind::Face>(face.GetIndex());
    if (m_Tables.Shells[shell].Solid.IsValid()) {
        UnlinkChild<EntityKind::Solid, EntityKind::Shell>(
            shell, &Detail::ShellRecord::Solid, &Detail::ShellRecord::NextInSolid,
            &Detail::SolidRecord::FirstShell, &Detail::SolidRecord::LastShell);
    }
    Release<EntityKind::Shell>(shell);
    Release<EntityKind::Vertex>(vertex.GetIndex());
    return {};
}

// --- MEV / KEV ----------------------------------------------------------------

/// The strut of MEV, every precondition checked: a vertex at `position`, an
/// edge to it from `from`, and its two uses after `after` in the loop (or
/// alone in it, for None).
Result<EdgeVertex> ShapeBuilder::GrowStrut(const std::uint32_t loop, const std::uint32_t after,
                                           const std::uint32_t from, const Point3& position,
                                           const double tolerance, const EdgeGeometry& geometry)
{
    EdgeVertex made;
    STEREON_TRY(made.Vertex, AddVertex(position, tolerance));
    EdgeSpec spec;
    spec.Kind = geometry.Kind;
    spec.Start = BuilderVertexId::Make(from, m_Tables.Vertices[from].Slot.GetGeneration(),
                                       m_Token.GetValue());
    spec.End = made.Vertex;
    spec.Curve = geometry.Curve;
    spec.Range = geometry.Range;
    spec.Tolerance = geometry.Tolerance;
    STEREON_TRY(made.Edge, AddEdge(spec));
    made.Out = NewUse(made.Edge.GetIndex(), loop, EdgeUseSense::Same, geometry.PCurve);
    made.Back = NewUse(made.Edge.GetIndex(), loop, EdgeUseSense::Opposite, geometry.PCurve);
    std::vector<std::uint32_t> cycle;
    if (after != None) {
        cycle = CycleFrom(after);
    }
    const auto at = cycle.empty() ? cycle.begin() : cycle.begin() + 1;
    cycle.insert(at, {made.Out.GetIndex(), made.Back.GetIndex()});
    RelinkCycle(loop, cycle);
    CheckNeighbourhood({loop}, {made.Edge.GetIndex()});
    return made;
}

Result<EdgeVertex> ShapeBuilder::MakeEdgeVertex(const BuilderCoedgeId at, const Point3& position,
                                                const double tolerance,
                                                const EdgeGeometry& geometry)
{
    constexpr const char* op = "MakeEdgeVertex";
    STEREON_TRY_VOID(Require(at, op));
    STEREON_TRY_VOID(Detail::CheckFinite(position, "the position", op));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the tolerance", op));
    STEREON_TRY_VOID(CheckEdgeGeometry(geometry, op));
    if (geometry.Kind == EdgeCurveKind::Degenerate) {
        return Refuse(op, "a degenerate edge has one vertex, and MEV makes a second");
    }
    const std::uint32_t loop = m_Tables.Coedges[at.GetIndex()].Loop.GetIndex();
    return GrowStrut(loop, at.GetIndex(), EndVertexOf(at.GetIndex()), position, tolerance,
                     geometry);
}

Result<EdgeVertex> ShapeBuilder::MakeEdgeVertex(const BuilderLoopId loop,
                                                const BuilderVertexId from, const Point3& position,
                                                const double tolerance,
                                                const EdgeGeometry& geometry)
{
    constexpr const char* op = "MakeEdgeVertex";
    STEREON_TRY_VOID(Require(loop, op));
    STEREON_TRY_VOID(Require(from, op));
    STEREON_TRY_VOID(Detail::CheckFinite(position, "the position", op));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the tolerance", op));
    STEREON_TRY_VOID(CheckEdgeGeometry(geometry, op));
    if (geometry.Kind == EdgeCurveKind::Degenerate) {
        return Refuse(op, "a degenerate edge has one vertex, and MEV makes a second");
    }
    if (m_Tables.Loops[loop.GetIndex()].Anchor.IsValid()) {
        return Refuse(op, "the loop is not empty; grow from one of its coedges");
    }
    if (m_Tables.Vertices[from.GetIndex()].FirstEdge.IsValid()) {
        return Refuse(op, "the vertex of an empty loop has no edges");
    }
    return GrowStrut(loop.GetIndex(), None, from.GetIndex(), position, tolerance, geometry);
}

Result<void> ShapeBuilder::KillEdgeVertex(const BuilderEdgeId edge, const BuilderVertexId vertex)
{
    constexpr const char* op = "KillEdgeVertex";
    STEREON_TRY_VOID(Require(edge, op));
    STEREON_TRY_VOID(Require(vertex, op));
    const std::uint32_t v = vertex.GetIndex();
    if (!IsFreeEnd(edge.GetIndex(), v)) {
        return Refuse(op, "the vertex must be an end of an edge that is not closed, with no "
                          "other edge at it");
    }
    const std::uint32_t first = m_Tables.Edges[edge.GetIndex()].RadialAnchor.GetIndex();
    const std::uint32_t second = SecondOfTwoUses(edge.GetIndex());
    if (second == None || m_Tables.Coedges[first].Loop != m_Tables.Coedges[second].Loop) {
        return Refuse(op, "the edge must have exactly two uses, in one loop");
    }
    // The use arriving at the vertex turns back along the other.
    const std::uint32_t in = EndVertexOf(first) == v ? first : second;
    const std::uint32_t out = in == first ? second : first;
    if (EndVertexOf(in) != v || m_Tables.Coedges[in].Next.GetIndex() != out) {
        return Refuse(op, "the edge's uses must meet at the vertex");
    }
    const std::uint32_t loop = m_Tables.Coedges[first].Loop.GetIndex();
    std::vector<std::uint32_t> cycle = CycleFrom(out);
    std::erase_if(cycle, [&](const std::uint32_t c) { return c == in || c == out; });
    RelinkCycle(loop, cycle);
    KillEdgeAndUses(edge.GetIndex());
    Release<EntityKind::Vertex>(v);
    CheckNeighbourhood({loop}, {});
    return {};
}

// --- MEF / KEF ----------------------------------------------------------------

Result<EdgeFace> ShapeBuilder::MakeEdgeFace(const BuilderCoedgeId a, const BuilderCoedgeId b,
                                            const EdgeGeometry& geometry)
{
    constexpr const char* op = "MakeEdgeFace";
    STEREON_TRY_VOID(Require(a, op));
    STEREON_TRY_VOID(Require(b, op));
    STEREON_TRY_VOID(CheckEdgeGeometry(geometry, op));
    const std::uint32_t loop = m_Tables.Coedges[a.GetIndex()].Loop.GetIndex();
    if (m_Tables.Coedges[b.GetIndex()].Loop.GetIndex() != loop) {
        return Refuse(op, "the two coedges must be in one loop");
    }
    const std::uint32_t from = EndVertexOf(a.GetIndex());
    const std::uint32_t to = EndVertexOf(b.GetIndex());
    if (geometry.Kind == EdgeCurveKind::Degenerate && from != to) {
        return Refuse(op, "a degenerate edge has one vertex");
    }
    const Detail::FaceRecord face = m_Tables.Faces[m_Tables.Loops[loop].Face.GetIndex()];
    const std::vector<std::uint32_t> cycle = CycleFrom(a.GetIndex());
    const auto bAt = std::ranges::find(cycle, b.GetIndex());

    EdgeFace made;
    EdgeSpec spec;
    spec.Kind = geometry.Kind;
    spec.Start = BuilderVertexId::Make(from, m_Tables.Vertices[from].Slot.GetGeneration(),
                                       m_Token.GetValue());
    spec.End =
        BuilderVertexId::Make(to, m_Tables.Vertices[to].Slot.GetGeneration(), m_Token.GetValue());
    spec.Curve = geometry.Curve;
    spec.Range = geometry.Range;
    spec.Tolerance = geometry.Tolerance;
    STEREON_TRY(made.Edge, AddEdge(spec));
    STEREON_TRY(made.Face, AddFace(face.Surface, face.Reversed));
    if (face.Shell.IsValid()) {
        const std::uint32_t s = face.Shell.GetIndex();
        STEREON_TRY_VOID(AddFaceToShell(
            BuilderShellId::Make(s, m_Tables.Shells[s].Slot.GetGeneration(), m_Token.GetValue()),
            made.Face));
    }
    STEREON_TRY(made.Loop, AddLoop(made.Face));
    made.Kept = NewUse(made.Edge.GetIndex(), loop, EdgeUseSense::Same, geometry.PCurve);
    made.New =
        NewUse(made.Edge.GetIndex(), made.Loop.GetIndex(), EdgeUseSense::Opposite, geometry.PCurve);

    // The original loop: a, the edge along, then what followed b. The new
    // loop: the edge against, then what followed a, up to b.
    std::vector<std::uint32_t> kept{a.GetIndex(), made.Kept.GetIndex()};
    kept.insert(kept.end(), bAt + 1, cycle.end());
    std::vector<std::uint32_t> split{made.New.GetIndex()};
    split.insert(split.end(), cycle.begin() + 1, bAt + 1);
    RelinkCycle(made.Loop.GetIndex(), split);
    RelinkCycle(loop, kept);
    CheckNeighbourhood({loop, made.Loop.GetIndex()}, {made.Edge.GetIndex()});
    return made;
}

Result<void> ShapeBuilder::KillEdgeFace(const BuilderEdgeId edge, const BuilderFaceId face)
{
    constexpr const char* op = "KillEdgeFace";
    STEREON_TRY_VOID(Require(edge, op));
    STEREON_TRY_VOID(Require(face, op));
    const std::uint32_t first = m_Tables.Edges[edge.GetIndex()].RadialAnchor.GetIndex();
    const std::uint32_t second = SecondOfTwoUses(edge.GetIndex());
    if (second == None) {
        return Refuse(op, "the edge must have exactly two uses");
    }
    const std::uint32_t gone = FaceOfUse(first) == face.GetIndex() ? first : second;
    const std::uint32_t kept = gone == first ? second : first;
    if (FaceOfUse(gone) != face.GetIndex() || FaceOfUse(kept) == face.GetIndex()) {
        return Refuse(op, "the edge must join the face to another face");
    }
    // A copy: the writes below may move the page this record is on.
    const Detail::FaceRecord dead = m_Tables.Faces[face.GetIndex()];
    STEREON_TRY_VOID(CheckJoinable(face.GetIndex(), FaceOfUse(kept), op));
    const std::uint32_t keptLoop = m_Tables.Coedges[kept].Loop.GetIndex();
    const std::uint32_t goneLoop = m_Tables.Coedges[gone].Loop.GetIndex();
    std::vector<std::uint32_t> joined = CycleFrom(kept);
    joined.erase(joined.begin());
    const std::vector<std::uint32_t> rest = CycleFrom(gone);
    joined.insert(joined.end(), rest.begin() + 1, rest.end());
    RelinkCycle(keptLoop, joined);
    m_Tables.Loops.Modify(goneLoop, [](Detail::LoopRecord& l) { l.Anchor = {}; });
    KillEdgeAndUses(edge.GetIndex());
    RemoveLoopAt(goneLoop);
    if (dead.Shell.IsValid()) {
        UnlinkChild<EntityKind::Shell, EntityKind::Face>(
            face.GetIndex(), &Detail::FaceRecord::Shell, &Detail::FaceRecord::NextInShell,
            &Detail::ShellRecord::FirstFace, &Detail::ShellRecord::LastFace);
    }
    Release<EntityKind::Face>(face.GetIndex());
    CheckNeighbourhood({keptLoop}, {});
    return {};
}

} // namespace Stereon
