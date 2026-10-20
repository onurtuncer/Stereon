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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Topo;

/// The Geometric level of the validity checker (ADR-0018): ADR-0001's
/// containment invariant. A vertex's tolerance sphere holds the ends of its
/// edges (G01), and an edge's tolerance tube holds its 3D curve and every
/// pcurve mapped through its face's surface, all on the edge's shared
/// parameter (G02, ADR-0009).
///
/// The check is sampled, not certified: each trace of an edge is compared
/// with the edge's reference trace at evenly spread parameters, ends
/// included. A deviation between samples can pass, which certified
/// distance bounds will close when interval evaluation exists. It runs on
/// structurally valid tables only, and it never evaluates outside a
/// domain: a trace that does not cover the edge's range, or a pcurve that
/// leaves its surface's domain, is itself a finding.

namespace Stereon::Detail {

// Instantiated at namespace scope before any function uses it
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

/// Intervals each edge's range is sampled in; the samples are their ends.
constexpr int SampleIntervals = 16;

[[nodiscard]] bool IsPeriodic(const double period) noexcept
{
    return period > 0.0; // STEREON-LINT-ALLOW(epsilon): 0 marks no period
}

[[nodiscard]] bool Covers(const Interval domain, const double period, const Interval range)
{
    return IsPeriodic(period) ||
           (domain.GetLower() <= range.GetLower() && range.GetUpper() <= domain.GetUpper());
}

[[nodiscard]] bool OnSurface(const Surface& surface, const Point2& uv)
{
    const Box2 domain = Domain(surface);
    const Vec2 periods = Periods(surface);
    return (IsPeriodic(periods.X) || Contains(domain.X, uv.X)) &&
           (IsPeriodic(periods.Y) || Contains(domain.Y, uv.Y));
}

/// One way to evaluate an edge: its own 3D curve, or a coedge's pcurve
/// mapped through its face's surface.
struct Trace {
    std::optional<Curve> Curve3d;
    std::optional<Curve2d> PCurve;
    std::optional<Surface> On;
    std::uint32_t Coedge = 0;
};

class GeometryValidator {
public:
    explicit GeometryValidator(const TopologyTables& tables)
        : m_T(tables)
    {}

    [[nodiscard]] ValidityReport Run()
    {
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (m_T.IsLive<EntityKind::Edge>(e)) {
                CheckEdge(e);
            }
        }
        return std::move(m_Report);
    }

private:
    void Report(const TopologyRule rule, std::vector<std::string> entities, std::string message)
    {
        m_Report.Add(
            {.Rule = rule, .Entities = std::move(entities), .Message = std::move(message)});
    }

    [[nodiscard]] std::string Name(const EntityKind kind, const std::uint32_t index) const
    {
        return DescribeEntity(m_T, kind, index);
    }

    /// The edge's traces: its own curve first, then the primary coedge's
    /// pcurve, then the other uses' pcurves. Empty geometry, a structural
    /// matter, is left out.
    [[nodiscard]] std::vector<Trace> TracesOf(const std::uint32_t e) const
    {
        const EdgeRecord& edge = m_T.Edges[e];
        std::vector<Trace> traces;
        if (edge.Curve.IsValid()) {
            Curve curve = m_T.Curves[edge.Curve.GetIndex()];
            if (!curve.IsEmpty()) {
                Trace own;
                own.Curve3d = std::move(curve);
                traces.push_back(std::move(own));
            }
        }
        // The primary coedge's trace leads where the edge has no curve of
        // its own: the edge evaluates through it (ADR-0009).
        const auto primarySlot = static_cast<std::ptrdiff_t>(traces.size());
        const std::uint32_t anchor = edge.RadialAnchor.GetIndex();
        std::uint32_t c = anchor;
        for (std::size_t steps = 0; c != CoedgeId::Invalid && steps < m_T.Coedges.GetSize();
             ++steps) {
            if (std::optional<Trace> use = UseTrace(c)) {
                const bool primary = c == edge.Primary.GetIndex();
                traces.insert(primary ? traces.begin() + primarySlot : traces.end(),
                              *std::move(use));
            }
            c = m_T.Coedges[c].Radial.GetIndex();
            if (c == anchor) {
                break;
            }
        }
        return traces;
    }

    [[nodiscard]] std::optional<Trace> UseTrace(const std::uint32_t c) const
    {
        const CoedgeRecord& use = m_T.Coedges[c];
        if (!use.PCurve.IsValid()) {
            return std::nullopt;
        }
        const FaceRecord& face = m_T.Faces[m_T.Loops[use.Loop.GetIndex()].Face.GetIndex()];
        Trace trace;
        trace.PCurve = m_T.PCurves[use.PCurve.GetIndex()];
        trace.On = m_T.Surfaces[face.Surface.GetIndex()];
        trace.Coedge = c;
        if (trace.PCurve->IsEmpty() || trace.On->IsEmpty()) {
            return std::nullopt;
        }
        return trace;
    }

    /// The trace's points at every sample of the range, or nothing, after
    /// reporting why, when it cannot be evaluated there.
    [[nodiscard]] std::optional<std::vector<Point3>>
    PointsOf(const Trace& trace, const std::uint32_t e, const Interval range)
    {
        std::vector<Point3> points;
        if (trace.Curve3d) {
            if (!Covers(Domain(*trace.Curve3d), Period(*trace.Curve3d), range)) {
                Report(TopologyRule::EdgeContainment, {Name(EntityKind::Edge, e)},
                       "the edge's 3D curve does not cover its range");
                return std::nullopt;
            }
            for (int k = 0; k <= SampleIntervals; ++k) {
                points.push_back(Eval(*trace.Curve3d, SampleAt(range, k)));
            }
            return points;
        }
        if (!Covers(Domain(*trace.PCurve), Period(*trace.PCurve), range)) {
            Report(TopologyRule::EdgeContainment,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, trace.Coedge)},
                   "the pcurve does not cover the edge's range");
            return std::nullopt;
        }
        for (int k = 0; k <= SampleIntervals; ++k) {
            const Point2 uv = Eval(*trace.PCurve, SampleAt(range, k));
            if (!OnSurface(*trace.On, uv)) {
                Report(TopologyRule::EdgeContainment,
                       {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, trace.Coedge)},
                       "the pcurve leaves its surface's domain");
                return std::nullopt;
            }
            points.push_back(Eval(*trace.On, uv));
        }
        return points;
    }

    /// Sample k of the range: its ends exactly, evenly spread between.
    [[nodiscard]] static double SampleAt(const Interval range, const int k)
    {
        if (k == SampleIntervals) {
            return range.GetUpper();
        }
        const double s = static_cast<double>(k) / SampleIntervals;
        // Rounding must not step outside a domain that is not periodic.
        return std::clamp(range.GetLower() * (1.0 - s) + range.GetUpper() * s, range.GetLower(),
                          range.GetUpper());
    }

    void CheckEdge(const std::uint32_t e)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        const Interval range = edge.Range;
        if (!std::isfinite(range.GetLower()) || !std::isfinite(range.GetUpper()) ||
            range.GetLower() > range.GetUpper()) {
            return; // nothing to sample; the structural check owns the range
        }
        const std::vector<Trace> traces = TracesOf(e);
        std::optional<std::vector<Point3>> reference;
        for (const Trace& trace : traces) {
            std::optional<std::vector<Point3>> points = PointsOf(trace, e, range);
            if (!points) {
                continue;
            }
            if (!reference) {
                reference = std::move(points);
                CheckVertices(e, reference->front(), reference->back());
                continue;
            }
            CheckTube(e, trace, *reference, *points);
        }
    }

    void CheckTube(const std::uint32_t e, const Trace& trace, const std::vector<Point3>& reference,
                   const std::vector<Point3>& points)
    {
        const double tolerance = m_T.Edges[e].Tolerance;
        double worst = 0.0;
        std::size_t at = 0;
        for (std::size_t k = 0; k < points.size(); ++k) {
            const double d = Distance(points[k], reference[k]);
            if (!(d <= worst)) {
                worst = d;
                at = k;
            }
        }
        if (!(worst <= tolerance)) {
            Report(TopologyRule::EdgeContainment,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, trace.Coedge)},
                   std::format("the pcurve on its surface is {} from the edge at t = {}, outside "
                               "its tolerance {}",
                               worst, SampleAt(m_T.Edges[e].Range, static_cast<int>(at)),
                               tolerance));
        }
    }

    void CheckVertices(const std::uint32_t e, const Point3& start, const Point3& end)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        CheckVertex(e, edge.Start.GetIndex(), start, "start");
        CheckVertex(e, edge.End.GetIndex(), end, "end");
    }

    void CheckVertex(const std::uint32_t e, const std::uint32_t v, const Point3& point,
                     const char* which)
    {
        if (!m_T.IsLive<EntityKind::Vertex>(v)) {
            return; // a degenerate or dangling link; the structural check owns it
        }
        const VertexRecord& vertex = m_T.Vertices[v];
        const double d = Distance(point, vertex.Position);
        if (!(d <= vertex.Tolerance)) {
            Report(TopologyRule::VertexContainment,
                   {Name(EntityKind::Vertex, v), Name(EntityKind::Edge, e)},
                   std::format("the edge's {} is {} from the vertex, outside its tolerance {}",
                               which, d, vertex.Tolerance));
        }
    }

    const TopologyTables& m_T;
    ValidityReport m_Report;
};

} // namespace

ValidityReport ValidateGeometry(const TopologyTables& tables)
{
    return GeometryValidator(tables).Run();
}

} // namespace Stereon::Detail
