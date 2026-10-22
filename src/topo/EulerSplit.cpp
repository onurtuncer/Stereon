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
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
module Stereon.Topo;
import Stereon.Robust;

namespace Stereon {
namespace {
[[nodiscard]] auto Refuse(const char* op, const char* why)
{
    return Fail(ErrorCode::InvalidInput, std::string(op) + ": " + why);
}
[[nodiscard]] bool IsPeriodic(const double period)
{
    return period > 0.0; // STEREON-LINT-ALLOW(epsilon): zero denotes no period
}
[[nodiscard]] Result<void> Within(const Point3& point, const Point3& center, const double tolerance)
{
    const auto comparison = CompareDistance(point, center, tolerance);
    if (!comparison.Supported || comparison.Value == Sign::Positive) {
        return Refuse("edge endpoint",
                      "the trace is outside the vertex tolerance or predicate domain");
    }
    return {};
}
} // namespace

/// Relink both endpoint lists, also when the old or new edge is closed.
void ShapeBuilder::RepointEdge(const std::uint32_t edge, const std::uint32_t start,
                               const std::uint32_t end)
{
    const auto old = m_Tables.Edges[edge];
    UnlinkAtVertex(edge, old.Start.GetIndex());
    if (old.End != old.Start) {
        UnlinkAtVertex(edge, old.End.GetIndex());
    }
    m_Tables.Edges.Modify(edge, [&](Detail::EdgeRecord& e) {
        e.Start = VertexId::FromIndex(start);
        e.End = VertexId::FromIndex(end);
        e.NextAtStart = {};
        e.NextAtEnd = {};
    });
    LinkAtVertex(edge, start, true);
    if (end != start) {
        LinkAtVertex(edge, end, false);
    }
    Advance();
}

Result<void> ShapeBuilder::CheckPointAt(const std::uint32_t edge, const double t,
                                        const Point3& position, const double tolerance) const
{
    const auto& e = m_Tables.Edges[edge];
    bool trace = false;
    if (e.Curve.IsValid()) {
        const Curve& curve = m_Tables.Curves[e.Curve.GetIndex()];
        if (!Contains(Domain(curve), t) && !IsPeriodic(Period(curve))) {
            return Refuse("edge endpoint", "the parameter leaves the curve domain");
        }
        STEREON_TRY_VOID(Within(Eval(curve, t), position, tolerance));
        trace = true;
    }
    for (const auto c : UsesOf(edge)) {
        const auto& use = m_Tables.Coedges[c];
        if (!use.PCurve.IsValid()) {
            continue;
        }
        STEREON_TRY_VOID(CheckUsePointAt(c, t, position, tolerance));
        trace = true;
    }
    if (!trace) {
        return Refuse("edge endpoint", "the edge has no evaluable trace");
    }
    return {};
}

Result<void> ShapeBuilder::CheckUsePointAt(const std::uint32_t c, const double t,
                                           const Point3& position, const double tolerance) const
{
    const auto& use = m_Tables.Coedges[c];
    const Curve2d& pcurve = m_Tables.PCurves[use.PCurve.GetIndex()];
    if (!Contains(Domain(pcurve), t) && !IsPeriodic(Period(pcurve))) {
        return Refuse("edge endpoint", "the parameter leaves a pcurve domain");
    }
    const auto uv = Eval(pcurve, t);
    const Surface& surface = m_Tables.Surfaces[m_Tables.Faces[FaceOfUse(c)].Surface.GetIndex()];
    const auto domain = Domain(surface);
    const auto periods = Periods(surface);
    if (!std::isfinite(uv.X) || !std::isfinite(uv.Y) ||
        (!IsPeriodic(periods.X) && !Contains(domain.X, uv.X)) ||
        (!IsPeriodic(periods.Y) && !Contains(domain.Y, uv.Y))) {
        return Refuse("edge endpoint", "the pcurve leaves the surface domain");
    }
    STEREON_TRY_VOID(Within(Eval(surface, uv), position, tolerance));
    return {};
}

/// Duplicate every use and splice next to its original, reversing the order
/// for opposite uses. The upper primary is the counterpart of the lower's.
void ShapeBuilder::SplitUses(const std::uint32_t lower, const std::uint32_t upper)
{
    const auto uses = UsesOf(lower);
    const auto primary = m_Tables.Edges[lower].Primary.GetIndex();
    for (const auto c : uses) {
        const auto use = m_Tables.Coedges[c];
        const auto made = NewUse(upper, use.Loop.GetIndex(), use.Sense, use.PCurve);
        auto cycle = CycleFrom(c);
        cycle.insert(use.Sense == EdgeUseSense::Same ? cycle.begin() + 1 : cycle.begin(),
                     made.GetIndex());
        RelinkCycle(use.Loop.GetIndex(), cycle);
        if (c == primary) {
            m_Tables.Edges.Modify(upper, [&](Detail::EdgeRecord& e) {
                e.Primary = CoedgeId::FromIndex(made.GetIndex());
            });
        }
        RecordEvent(EntityKind::Coedge, TopologyEventKind::Split, c, made.GetIndex());
    }
}

Result<EdgeSplit> ShapeBuilder::SplitEdgeMakeVertex(const BuilderEdgeId edge, const double t,
                                                    const Point3& position, const double tolerance)
{
    constexpr const char* op = "SplitEdgeMakeVertex";
    STEREON_TRY_VOID(Require(edge, op));
    STEREON_TRY_VOID(Detail::CheckFinite(position, "the split position", op));
    STEREON_TRY_VOID(Detail::CheckNonNegative(tolerance, "the vertex tolerance", op));
    const auto original = m_Tables.Edges[edge.GetIndex()];
    if (!std::isfinite(t) || !(original.Range.GetLower() < t && t < original.Range.GetUpper()) ||
        original.Kind == EdgeCurveKind::Degenerate) {
        return Refuse(op, "a nondegenerate edge and a strictly interior parameter are required");
    }
    STEREON_TRY_VOID(CheckPointAt(edge.GetIndex(), t, position, tolerance));
    EdgeSplit made;
    STEREON_TRY(made.Vertex, AddVertex(position, tolerance));
    EdgeSpec spec;
    spec.Kind = original.Kind;
    spec.Curve = original.Curve;
    spec.Start = made.Vertex;
    const auto end = original.End.GetIndex();
    spec.End =
        BuilderVertexId::Make(end, m_Tables.Vertices[end].Slot.GetGeneration(), m_Token.GetValue());
    spec.Range = Interval{t, original.Range.GetUpper()};
    spec.Tolerance = original.Tolerance;
    STEREON_TRY(made.Upper, AddEdge(spec));
    RepointEdge(edge.GetIndex(), original.Start.GetIndex(), made.Vertex.GetIndex());
    m_Tables.Edges.Modify(edge.GetIndex(), [&](Detail::EdgeRecord& e) {
        e.Range = Interval{original.Range.GetLower(), t};
    });
    SplitUses(edge.GetIndex(), made.Upper.GetIndex());
    m_Tables.Edges.Modify(made.Upper.GetIndex(), [&](Detail::EdgeRecord& e) {
        e.Cached = original.Cached;
        e.CacheStatus = original.CacheStatus;
        e.CacheBound = original.CacheBound;
    });
    RecordEvent(EntityKind::Vertex, TopologyEventKind::Created, NoEntity, made.Vertex.GetIndex());
    RecordEvent(EntityKind::Edge, TopologyEventKind::Split, edge.GetIndex(), made.Upper.GetIndex());
    std::vector<std::uint32_t> loops;
    for (const auto c : UsesOf(edge.GetIndex())) {
        loops.push_back(m_Tables.Coedges[c].Loop.GetIndex());
    }
    CheckNeighbourhood(loops, {edge.GetIndex(), made.Upper.GetIndex()});
    return made;
}

Result<void> ShapeBuilder::CheckEdgeJoin(const std::uint32_t lower, const std::uint32_t upper) const
{
    constexpr const char* op = "JoinEdgeKillVertex";
    const auto& a = m_Tables.Edges[lower];
    const auto& b = m_Tables.Edges[upper];
    if (lower == upper || a.End != b.Start || a.Start == a.End || b.Start == b.End ||
        a.Range.GetUpper() != b.Range.GetLower() || a.Kind != b.Kind || a.Curve != b.Curve ||
        a.Tolerance != b.Tolerance || a.Kind == EdgeCurveKind::Degenerate) {
        return Refuse(op, "contiguous edges of the same geometry and tolerance are required");
    }
    STEREON_TRY_VOID(CheckJoinVertex(lower, upper));
    return CheckJoinUses(lower, upper);
}

Result<void> ShapeBuilder::CheckJoinVertex(const std::uint32_t lower,
                                           const std::uint32_t upper) const
{
    constexpr const char* op = "JoinEdgeKillVertex";
    const auto& a = m_Tables.Edges[lower];
    const auto vertex = a.End.GetIndex();
    const auto first = m_Tables.Vertices[vertex].FirstEdge.GetIndex();
    const auto other = first == lower ? upper : lower;
    if ((first != lower && first != upper) ||
        (m_Tables.Edges[first].Start.GetIndex() == vertex ? m_Tables.Edges[first].NextAtStart
                                                          : m_Tables.Edges[first].NextAtEnd)
                .GetIndex() != other ||
        (m_Tables.Edges[other].Start.GetIndex() == vertex ? m_Tables.Edges[other].NextAtStart
                                                          : m_Tables.Edges[other].NextAtEnd)
            .IsValid()) {
        return Refuse(op, "the join vertex must have exactly these two incident edges");
    }
    return {};
}

Result<void> ShapeBuilder::CheckJoinUses(const std::uint32_t lower, const std::uint32_t upper) const
{
    constexpr const char* op = "JoinEdgeKillVertex";
    const auto& a = m_Tables.Edges[lower];
    const auto& b = m_Tables.Edges[upper];
    const auto uses = UsesOf(lower);
    if (uses.size() != UsesOf(upper).size()) {
        return Refuse(op, "radial valences differ");
    }
    for (const auto c : uses) {
        const auto& use = m_Tables.Coedges[c];
        const auto paired = (use.Sense == EdgeUseSense::Same ? use.Next : use.Prev).GetIndex();
        const auto& mate = m_Tables.Coedges[paired];
        if (mate.Edge.GetIndex() != upper || mate.Sense != use.Sense || mate.Loop != use.Loop ||
            mate.PCurve != use.PCurve) {
            return Refuse(op, "every use needs an adjacent counterpart with the same pcurve");
        }
        if (a.Primary.GetIndex() == c && b.Primary.GetIndex() != paired) {
            return Refuse(op, "the primary sides do not correspond");
        }
    }
    return {};
}

Result<void> ShapeBuilder::JoinEdgeKillVertex(const BuilderEdgeId lower, const BuilderEdgeId upper)
{
    constexpr const char* op = "JoinEdgeKillVertex";
    STEREON_TRY_VOID(Require(lower, op));
    STEREON_TRY_VOID(Require(upper, op));
    STEREON_TRY_VOID(CheckEdgeJoin(lower.GetIndex(), upper.GetIndex()));
    const auto a = m_Tables.Edges[lower.GetIndex()];
    const auto b = m_Tables.Edges[upper.GetIndex()];
    std::vector<std::uint32_t> loops;
    for (const auto c : UsesOf(upper.GetIndex())) {
        const auto loop = m_Tables.Coedges[c].Loop.GetIndex();
        auto cycle = CycleFrom(c);
        std::erase(cycle, c);
        RelinkCycle(loop, cycle);
        loops.push_back(loop);
    }
    KillEdgeAndUses(upper.GetIndex());
    RepointEdge(lower.GetIndex(), a.Start.GetIndex(), b.End.GetIndex());
    m_Tables.Edges.Modify(lower.GetIndex(), [&](Detail::EdgeRecord& e) {
        e.Range = Interval{a.Range.GetLower(), b.Range.GetUpper()};
        e.CacheBound = std::max(a.CacheBound, b.CacheBound);
    });
    if (a.Cached != b.Cached || a.CacheStatus != b.CacheStatus) {
        DropCache(lower.GetIndex());
    }
    Release<EntityKind::Vertex>(a.End.GetIndex());
    RecordEvent(EntityKind::Vertex, TopologyEventKind::Deleted, a.End.GetIndex(), NoEntity);
    RecordEvent(EntityKind::Edge, TopologyEventKind::Merged, upper.GetIndex(), lower.GetIndex());
    CheckNeighbourhood(loops, {lower.GetIndex()});
    return {};
}
Result<VertexMergeResult>
ShapeBuilder::MergeVertices(const std::span<const BuilderVertexId> vertices, const Context& context)
{
    for (const auto vertex : vertices) {
        STEREON_TRY_VOID(Require(vertex, "MergeVertices"));
    }
    return Transaction(
        [&](ShapeBuilder& candidate) { return candidate.MergeVerticesImpl(vertices, context); });
}

Result<VertexMergeResult>
ShapeBuilder::MergeVerticesImpl(const std::span<const BuilderVertexId> vertices,
                                const Context& context)
{
    std::vector<BuilderVertexId> sorted(vertices.begin(), vertices.end());
    std::ranges::sort(sorted, {}, &BuilderVertexId::GetIndex);
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::vector<std::vector<BuilderVertexId>> clusters;
    for (const auto vertex : sorted) {
        const auto& member = m_Tables.Vertices[vertex.GetIndex()];
        bool assigned = false;
        for (auto& cluster : clusters) {
            const auto& representative = m_Tables.Vertices[cluster.front().GetIndex()];
            const auto comparison =
                CompareDistance(member.Position, representative.Position,
                                std::max(member.Tolerance, representative.Tolerance));
            if (!comparison.Supported) {
                return Refuse("MergeVertices", "unsupported predicate domain");
            }
            if (comparison.Value != Sign::Positive) {
                cluster.push_back(vertex);
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            clusters.push_back({vertex});
        }
    }
    VertexMergeResult result;
    for (const auto& cluster : clusters) {
        STEREON_TRY_VOID(MergeVertexCluster(cluster, context));
        for (const auto source : cluster) {
            VertexMergeEntry entry;
            entry.Source = source;
            entry.Representative = cluster.front();
            result.Entries.push_back(entry);
        }
    }
    std::ranges::sort(result.Entries, [](const auto& a, const auto& b) {
        return a.Source.GetIndex() < b.Source.GetIndex();
    });
    return result;
}

Result<void> ShapeBuilder::MergeVertexCluster(const std::vector<BuilderVertexId>& cluster,
                                              const Context& context)
{
    if (cluster.size() < 2) {
        return {};
    }
    const auto kept = cluster.front();
    const auto original = m_Tables.Vertices[kept.GetIndex()];
    double tolerance = std::max(original.Tolerance, context.GetResolution());
    for (const auto vertex : cluster) {
        tolerance = std::max(tolerance, m_Tables.Vertices[vertex.GetIndex()].Tolerance);
    }
    if (tolerance > context.GetToleranceCeiling() ||
        (tolerance > original.Tolerance &&
         context.GetToleranceGrowthPolicy() == ToleranceGrowth::Fail)) {
        return Fail(ErrorCode::ToleranceExceeded, "MergeVertices: tolerance growth exceeds policy");
    }
    STEREON_TRY_VOID(CheckMergeEndpoints(cluster, original.Position, tolerance));
    if (tolerance != original.Tolerance) {
        STEREON_TRY_VOID(SetTolerance(kept, tolerance, "deterministic vertex merge"));
    }
    for (const auto vertex : cluster) {
        if (vertex == kept) {
            continue;
        }
        while (m_Tables.Vertices[vertex.GetIndex()].FirstEdge.IsValid()) {
            const auto edge = m_Tables.Vertices[vertex.GetIndex()].FirstEdge.GetIndex();
            const auto record = m_Tables.Edges[edge];
            RepointEdge(edge,
                        record.Start.GetIndex() == vertex.GetIndex() ? kept.GetIndex()
                                                                     : record.Start.GetIndex(),
                        record.End.GetIndex() == vertex.GetIndex() ? kept.GetIndex()
                                                                   : record.End.GetIndex());
        }
        Release<EntityKind::Vertex>(vertex.GetIndex());
        RecordEvent(EntityKind::Vertex, TopologyEventKind::Merged, vertex.GetIndex(),
                    kept.GetIndex());
    }
    return {};
}

Result<void> ShapeBuilder::CheckMergeEndpoints(const std::vector<BuilderVertexId>& cluster,
                                               const Point3& position, const double tolerance) const
{
    for (const auto vertex : cluster) {
        auto edge = m_Tables.Vertices[vertex.GetIndex()].FirstEdge;
        while (edge.IsValid()) {
            const auto index = edge.GetIndex();
            const auto record = m_Tables.Edges[index];
            const bool start = record.Start.GetIndex() == vertex.GetIndex();
            edge = start ? record.NextAtStart : record.NextAtEnd;
            if (record.Start.GetIndex() == vertex.GetIndex()) {
                STEREON_TRY_VOID(CheckPointAt(index, record.Range.GetLower(), position, tolerance));
            }
            if (record.End.GetIndex() == vertex.GetIndex()) {
                STEREON_TRY_VOID(CheckPointAt(index, record.Range.GetUpper(), position, tolerance));
            }
        }
    }
    return {};
}

} // namespace Stereon
