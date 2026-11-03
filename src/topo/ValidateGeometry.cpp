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
/// leaves its surface's domain, is itself a finding. A Degenerate edge's
/// traces must also stay at its vertex over the whole range (G02): a pole
/// has no 3D extent.

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
    if (!std::isfinite(uv.X) || !std::isfinite(uv.Y)) {
        return false;
    }
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
                CheckSurfaces(e);
            }
        }
        for (std::uint32_t l = 0; l < m_T.Loops.GetSize(); ++l) {
            if (m_T.IsLive<EntityKind::Loop>(l)) {
                CheckLoop(l);
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

    /// Every evaluated point must be finite, even on a lone Spatial trace.
    [[nodiscard]] bool AddPoint(std::vector<Point3>& points, const Point3& point,
                                const Trace& trace, const std::uint32_t e, const double t)
    {
        if (!std::isfinite(point.X) || !std::isfinite(point.Y) || !std::isfinite(point.Z)) {
            std::vector<std::string> entities{Name(EntityKind::Edge, e)};
            if (!trace.Curve3d) {
                entities.push_back(Name(EntityKind::Coedge, trace.Coedge));
            }
            Report(TopologyRule::EdgeContainment, std::move(entities),
                   std::format("the edge's trace is not finite at t = {}", t));
            return false;
        }
        points.push_back(point);
        return true;
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
                const double t = SampleAt(range, k);
                if (!AddPoint(points, Eval(*trace.Curve3d, t), trace, e, t)) {
                    return std::nullopt;
                }
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
            if (!AddPoint(points, Eval(*trace.On, uv), trace, e, SampleAt(range, k))) {
                return std::nullopt;
            }
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
        if (range.IsEmpty() || !range.IsBounded()) {
            // The structural check refuses such a range; on tables that
            // skipped it the containment cannot be sampled, which is a
            // finding, never a pass.
            Report(TopologyRule::EdgeContainment, {Name(EntityKind::Edge, e)},
                   "the edge's range is empty or not bounded, so its containment cannot be "
                   "checked");
            return;
        }
        const std::vector<Trace> traces = TracesOf(e);
        std::optional<std::vector<Point3>> reference;
        for (const Trace& trace : traces) {
            std::optional<std::vector<Point3>> points = PointsOf(trace, e, range);
            if (!points) {
                continue;
            }
            if (edge.Kind == EdgeCurveKind::Degenerate) {
                CheckCollapse(e, trace, *points);
            }
            if (!reference) {
                reference = std::move(points);
                CheckVertices(e, reference->front(), reference->back());
                continue;
            }
            CheckTube(e, trace, *reference, *points);
        }
        if (reference) {
            CheckCache(e, *reference);
        }
    }

    // --- G06: a cached curve stays in its tube and within its bound ---------

    /// The largest distance between the points and the reference, the first
    /// non-finite one stopping the search so that it is never erased.
    [[nodiscard]] static std::pair<double, std::size_t> Worst(const std::vector<Point3>& points,
                                                              const std::vector<Point3>& reference)
    {
        double worst = 0.0;
        std::size_t at = 0;
        for (std::size_t k = 0; k < points.size(); ++k) {
            const double d = Distance(points[k], reference[k]);
            if (!(d <= worst)) {
                worst = d;
                at = k;
            }
            if (!std::isfinite(d)) {
                break;
            }
        }
        return {worst, at};
    }

    /// An Intersection or Trim edge's cache approximates its primary side
    /// S(c(t)) (ADR-0009, section 5): it must lie in the edge's tube, as
    /// every trace of the edge does, and a Certified cache within its proven
    /// bound. A sampled deviation beyond the bound proves the certificate
    /// wrong or stale: the cache no longer matches the pcurve or surface it
    /// was built from.
    void CheckCache(const std::uint32_t e, const std::vector<Point3>& reference)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        if (!edge.Cached.IsValid()) {
            return;
        }
        const Curve cache = m_T.Curves[edge.Cached.GetIndex()];
        if (cache.IsEmpty()) {
            return; // the structural check owns an empty cache (S17)
        }
        if (!Covers(Domain(cache), Period(cache), edge.Range)) {
            Report(TopologyRule::CacheDeviation, {Name(EntityKind::Edge, e)},
                   "the cached curve does not cover the edge's range");
            return;
        }
        std::vector<Point3> points;
        for (int k = 0; k <= SampleIntervals; ++k) {
            points.push_back(Eval(cache, SampleAt(edge.Range, k)));
        }
        const auto [worst, at] = Worst(points, reference);
        const double t = SampleAt(edge.Range, static_cast<int>(at));
        if (edge.CacheStatus == CacheCertificate::Certified && !(worst <= edge.CacheBound)) {
            Report(TopologyRule::CacheDeviation, {Name(EntityKind::Edge, e)},
                   std::format("the certified cache is {} from the edge at t = {}, beyond its "
                               "proven bound {}: the certificate is stale",
                               worst, t, edge.CacheBound));
        } else if (!(worst <= edge.Tolerance)) {
            Report(TopologyRule::CacheDeviation, {Name(EntityKind::Edge, e)},
                   std::format("the cached curve is {} from the edge at t = {}, outside its "
                               "tolerance {}",
                               worst, t, edge.Tolerance));
        }
    }

    // --- G05: Trim and Intersection edges and their surfaces ----------------

    /// The surfaces of the faces that use the edge, in ring order.
    [[nodiscard]] std::vector<SurfaceRef> SurfacesOf(const std::uint32_t e) const
    {
        std::vector<SurfaceRef> surfaces;
        const std::uint32_t anchor = m_T.Edges[e].RadialAnchor.GetIndex();
        std::uint32_t c = anchor;
        for (std::size_t steps = 0; c != CoedgeId::Invalid && steps < m_T.Coedges.GetSize();
             ++steps) {
            const std::uint32_t face = m_T.Loops[m_T.Coedges[c].Loop.GetIndex()].Face.GetIndex();
            surfaces.push_back(m_T.Faces[face].Surface);
            c = m_T.Coedges[c].Radial.GetIndex();
            if (c == anchor) {
                break;
            }
        }
        return surfaces;
    }

    /// ADR-0009: a Trim edge's pcurve is shared by both its faces, so they
    /// lie on the same surface; an Intersection edge joins different ones.
    void CheckSurfaces(const std::uint32_t e)
    {
        const EdgeCurveKind kind = m_T.Edges[e].Kind;
        if (kind != EdgeCurveKind::Trim && kind != EdgeCurveKind::Intersection) {
            return;
        }
        const std::vector<SurfaceRef> surfaces = SurfacesOf(e);
        if (surfaces.empty()) {
            return;
        }
        const bool one = std::ranges::all_of(
            surfaces, [&](const SurfaceRef s) { return s == surfaces.front(); });
        if (kind == EdgeCurveKind::Trim && !one) {
            Report(TopologyRule::EdgeSurfaces, {Name(EntityKind::Edge, e)},
                   "the faces of a Trim edge lie on different surfaces");
        } else if (kind == EdgeCurveKind::Intersection && one) {
            Report(TopologyRule::EdgeSurfaces, {Name(EntityKind::Edge, e)},
                   "the faces of an Intersection edge all lie on one surface");
        }
    }

    // --- G03: loops closed in parameter space -------------------------------

    /// The point where the coedge's pcurve ends (`atEnd`) or starts in its
    /// loop's direction: the edge's range end for a use along the edge, its
    /// start for one against it. Nothing where the pcurve is missing or does
    /// not cover the range, which G02 reports.
    [[nodiscard]] std::optional<Point2> PcurvePoint(const std::uint32_t c, const bool atEnd) const
    {
        const CoedgeRecord& use = m_T.Coedges[c];
        if (!use.PCurve.IsValid()) {
            return std::nullopt;
        }
        const Curve2d pcurve = m_T.PCurves[use.PCurve.GetIndex()];
        const Interval range = m_T.Edges[use.Edge.GetIndex()].Range;
        if (pcurve.IsEmpty() || range.IsEmpty() || !range.IsBounded() ||
            !Covers(Domain(pcurve), Period(pcurve), range)) {
            return std::nullopt;
        }
        const bool upper = atEnd == (use.Sense == EdgeUseSense::Same);
        return Eval(pcurve, upper ? range.GetUpper() : range.GetLower());
    }

    /// How far the images of two pcurve ends meeting at a vertex may be
    /// apart: each lies within its edge's tolerance of the edge's end, which
    /// lies within the vertex's tolerance of the vertex (ADR-0001).
    [[nodiscard]] double JunctionLimit(const std::uint32_t c, const std::uint32_t n) const
    {
        const CoedgeRecord& use = m_T.Coedges[c];
        const EdgeRecord& edge = m_T.Edges[use.Edge.GetIndex()];
        const VertexId at = use.Sense == EdgeUseSense::Same ? edge.End : edge.Start;
        const double vertex = m_T.Vertices[at.GetIndex()].Tolerance;
        return 2.0 * vertex + edge.Tolerance + m_T.Edges[m_T.Coedges[n].Edge.GetIndex()].Tolerance;
    }

    /// Consecutive coedges of a loop meet in the surface's parameters, not
    /// only in space: the gap between one pcurve's end and the next one's
    /// start, carried onto the surface by its first derivatives, stays
    /// within the junction's tolerances. Not taken modulo a period: a loop
    /// crosses a seam through the seam's edge, whose two uses carry its two
    /// sides (ADR-0009), so a jump of a period is a finding, though its two
    /// ends are one point in space. A jump across a pole costs nothing, the
    /// derivative along it being zero; a Degenerate edge carries it there.
    void CheckJunction(const Surface& surface, const std::uint32_t l, const std::uint32_t c,
                       const std::uint32_t n)
    {
        const std::optional<Point2> end = PcurvePoint(c, true);
        const std::optional<Point2> start = PcurvePoint(n, false);
        if (!end || !start || !OnSurface(surface, *end)) {
            return;
        }
        const Vec2 gap = *start - *end;
        const double onSurface =
            Length(Deriv(surface, *end, 1, 0) * gap.X + Deriv(surface, *end, 0, 1) * gap.Y);
        const double limit = JunctionLimit(c, n);
        if (!(onSurface <= limit)) {
            Report(TopologyRule::LoopParameterGap,
                   {Name(EntityKind::Loop, l), Name(EntityKind::Coedge, c),
                    Name(EntityKind::Coedge, n)},
                   std::format("consecutive pcurves are ({}, {}) apart in the surface's "
                               "parameters, {} on the surface, beyond the junction's {}",
                               gap.X, gap.Y, onSurface, limit));
        }
    }

    void CheckLoop(const std::uint32_t l)
    {
        const LoopRecord& loop = m_T.Loops[l];
        const SurfaceRef ref = m_T.Faces[loop.Face.GetIndex()].Surface;
        if (!ref.IsValid()) {
            return;
        }
        const Surface surface = m_T.Surfaces[ref.GetIndex()];
        if (surface.IsEmpty()) {
            return;
        }
        const std::uint32_t anchor = loop.Anchor.GetIndex();
        std::uint32_t c = anchor;
        for (std::size_t steps = 0; c != CoedgeId::Invalid && steps < m_T.Coedges.GetSize();
             ++steps) {
            const std::uint32_t n = m_T.Coedges[c].Next.GetIndex();
            CheckJunction(surface, l, c, n);
            c = n;
            if (c == anchor) {
                break;
            }
        }
    }

    void CheckTube(const std::uint32_t e, const Trace& trace, const std::vector<Point3>& reference,
                   const std::vector<Point3>& points)
    {
        const double tolerance = m_T.Edges[e].Tolerance;
        const auto [worst, at] = Worst(points, reference);
        if (!(worst <= tolerance)) {
            Report(TopologyRule::EdgeContainment,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, trace.Coedge)},
                   std::format("the pcurve on its surface is {} from the edge at t = {}, outside "
                               "its tolerance {}",
                               worst, SampleAt(m_T.Edges[e].Range, static_cast<int>(at)),
                               tolerance));
        }
    }

    /// A Degenerate edge has no 3D extent (ADR-0009): each of its traces
    /// stays at its one vertex over the whole range, within the larger of the
    /// edge's and the vertex's tolerances. Coincident ends alone do not make
    /// a pole: a closed circle has them too.
    void CheckCollapse(const std::uint32_t e, const Trace& trace, const std::vector<Point3>& points)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        const std::uint32_t v = edge.Start.GetIndex();
        if (!m_T.IsLive<EntityKind::Vertex>(v)) {
            return; // a dangling link; the structural check owns it
        }
        const VertexRecord& vertex = m_T.Vertices[v];
        const double tolerance = std::max(edge.Tolerance, vertex.Tolerance);
        double worst = 0.0;
        std::size_t at = 0;
        for (std::size_t k = 0; k < points.size(); ++k) {
            const double d = Distance(points[k], vertex.Position);
            if (!(d <= worst)) {
                worst = d;
                at = k;
            }
            if (!std::isfinite(d)) {
                break; // a later finite sample must never erase this failure
            }
        }
        if (!(worst <= tolerance)) {
            std::vector<std::string> entities{Name(EntityKind::Edge, e)};
            if (!trace.Curve3d) {
                entities.push_back(Name(EntityKind::Coedge, trace.Coedge));
            }
            Report(TopologyRule::EdgeContainment, std::move(entities),
                   std::format("the degenerate edge's trace is {} from its vertex at t = {}, "
                               "outside its tolerance {}",
                               worst, SampleAt(edge.Range, static_cast<int>(at)), tolerance));
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
    STEREON_PROFILE_SCOPE("topo.ValidateGeometry");
    STEREON_PROFILE_COUNTER("topo.geometric.face_slots", tables.Faces.GetSize());
    STEREON_PROFILE_COUNTER("topo.geometric.edge_slots", tables.Edges.GetSize());
    return GeometryValidator(tables).Run();
}

namespace {

void CheckToleranceBound(ValidityReport& report, const TopologyTables& tables,
                         const EntityKind kind, const std::uint32_t index, const double tolerance,
                         const Context& context)
{
    if (tolerance >= context.GetResolution() && tolerance <= context.GetToleranceCeiling()) {
        return;
    }
    report.Add({.Rule = TopologyRule::ToleranceBounds,
                .Entities = {DescribeEntity(tables, kind, index)},
                .Message = std::format("the tolerance {} is outside the Context's [resolution, "
                                       "ceiling] of [{}, {}]",
                                       tolerance, context.GetResolution(),
                                       context.GetToleranceCeiling())});
}

} // namespace

void ValidateTolerances(const TopologyTables& tables, const Context& context,
                        ValidityReport& report)
{
    for (std::uint32_t v = 0; v < tables.Vertices.GetSize(); ++v) {
        if (tables.IsLive<EntityKind::Vertex>(v)) {
            CheckToleranceBound(report, tables, EntityKind::Vertex, v, tables.Vertices[v].Tolerance,
                                context);
        }
    }
    for (std::uint32_t e = 0; e < tables.Edges.GetSize(); ++e) {
        if (tables.IsLive<EntityKind::Edge>(e)) {
            CheckToleranceBound(report, tables, EntityKind::Edge, e, tables.Edges[e].Tolerance,
                                context);
        }
    }
}

void ValidateSolid(const TopologyTables& tables, ValidityReport& report)
{
    for (std::uint32_t s = 0; s < tables.Solids.GetSize(); ++s) {
        if (tables.IsLive<EntityKind::Solid>(s)) {
            report.Add({.Rule = TopologyRule::SolidUnchecked,
                        .Entities = {DescribeEntity(tables, EntityKind::Solid, s)},
                        .Message = "the Solid level's rules (ADR-0018) are not implemented: "
                                   "closed shells, orientation, signed volume, nesting and "
                                   "self-intersection are unchecked, so the solid is not "
                                   "reported valid at this level"});
        }
    }
}

} // namespace Stereon::Detail
