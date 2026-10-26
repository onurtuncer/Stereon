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
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Winding;

namespace Stereon {

namespace {

/// The first round's deviation target is the body's extent over this.
constexpr double FirstTargetDivisor = 64.0;
/// Segments an edge is sampled into at most, whatever the target asks; past
/// it the edge's chords stray further than the target, and the deviation
/// the sampling achieves, not the target, enters the bound.
constexpr double MaxSegments = 1 << 20;
/// Cells are split at most this deep, whatever the target asks.
constexpr int MaxDepth = 60;
/// Rounding allowed per triangle in the sum of solid angles over 4 pi. Each
/// atan2 is within a few units in the last place of pi, but a point near a
/// face sums thousands of triangles whose large angles cancel, and the sum
/// loses a few more digits; the allowance is generous, since the margin it
/// is held against is one half.
constexpr double RoundingPerTriangle = 0x1p-40;

/// A corner of a face's polygon: in the face's parameter space, in space,
/// and a bound on how far the position lies from the surface at the
/// parameters, beyond the boundary's own deviation. A boundary corner is on
/// the surface; a corner a cell split makes lies on the chord it cuts, so
/// that neighbouring cells share their sides exactly and the triangles close
/// up, and is off the surface by its chord's sagitta and its ends' offsets.
struct Corner {
    Point2 Uv;
    Point3 Position;
    double Off = 0.0;
};

/// A closed polygon of corners, one loop of a face.
using Polyline = std::vector<Corner>;

/// How a surface stretches its parameter space: the length in space of a
/// unit step in u and in v (so UScale + VScale bounds |dS| / |d(u, v)|), the
/// least stretch of any step (so a gap in space is at least MinScale times
/// the gap in parameters), and the radius a cylinder's chords stray by; 0
/// for a plane.
struct SurfaceScale {
    double UScale = 1.0;
    double VScale = 1.0;
    double MinScale = 1.0;
    double Radius = 0.0;

    [[nodiscard]] double Lipschitz() const noexcept { return UScale + VScale; }
};

/// Nothing for a surface outside the prototype's planes and cylinders.
[[nodiscard]] std::optional<SurfaceScale> ScaleOf(const Surface& surface)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        const double x = Length(plane->XAxis);
        const double y = Length(plane->YAxis);
        return SurfaceScale{x, y, std::min(x, y), 0.0};
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        const double r = cylinder->Radius;
        const double axis = Length(Axis(*cylinder));
        return SurfaceScale{r * (Length(cylinder->XAxis) + Length(cylinder->YAxis)), axis,
                            std::min(r, axis), r};
    }
    return std::nullopt;
}

/// How far a cylinder of radius r strays from the chord of an arc turning
/// du: r (1 - cos(du / 2)) up to a half turn, and never more than the
/// diameter.
[[nodiscard]] double Sagitta(const double r, const double du)
{
    return r * std::min(du * du / 8.0, 2.0);
}

// --- Bounds on the pcurves ---------------------------------------------------

/// The poles of a B-spline's derivative, as vectors, with the knots they are
/// on.
struct Hodograph {
    std::vector<Vec2> Poles;
    std::vector<double> Knots;
};

/// The derivative of a non-rational B-spline of this degree (The NURBS Book,
/// section 3.3): the pole differences scaled by the degree over the knot span
/// each straddles, on the knots without their first and last. Nothing when a
/// span vanishes, which a valid clamped curve has only at its ends.
[[nodiscard]] std::optional<Hodograph> Differentiated(const Hodograph& c, const int degree)
{
    Hodograph d;
    const auto p = static_cast<std::size_t>(degree);
    for (std::size_t i = 0; i + 1 < c.Poles.size(); ++i) {
        const double span = c.Knots[i + p + 1] - c.Knots[i + 1];
        if (!(span > 0.0)) { // STEREON-LINT-ALLOW(epsilon): a sign test on a knot span
            return std::nullopt;
        }
        d.Poles.push_back(static_cast<double>(degree) / span * (c.Poles[i + 1] - c.Poles[i]));
    }
    d.Knots.assign(c.Knots.begin() + 1, c.Knots.end() - 1);
    return d;
}

/// Bounds on a pcurve's derivatives over its whole range: on |u'|, the rate
/// it turns about a cylinder's axis, and on |p''|, its own curving in the
/// parameter plane. The edge in space is S(p(t)), whose second derivative is
/// within r u'^2 + (|Su| + |Sv|) |p''| on a cylinder of radius r (the only
/// second derivative of S is Suu, of length r), and within (|Su| + |Sv|)
/// |p''| on a plane.
struct PcurveRates {
    double Turn = 0.0;
    double Curving = 0.0;
};

/// The largest |X| component and the largest length among the poles.
[[nodiscard]] std::pair<double, double> PoleBounds(const Hodograph& h)
{
    double x = 0.0;
    double length = 0.0;
    for (const Vec2& pole : h.Poles) {
        x = std::max(x, std::abs(pole.X));
        length = std::max(length, Length(pole));
    }
    return {x, length};
}

/// The rates of a non-rational B-spline: each derivative curve lies in the
/// hull of its poles. Nothing for a rational one.
[[nodiscard]] std::optional<PcurveRates> BSplineRates(const BSplineCurve2d& b)
{
    if (!b.Weights.empty() || b.Degree < 1) {
        return std::nullopt;
    }
    Hodograph c{{}, b.Knots};
    for (const Point2& pole : b.Poles) {
        c.Poles.push_back(pole - Point2{});
    }
    const std::optional<Hodograph> first = Differentiated(c, b.Degree);
    if (!first) {
        return std::nullopt;
    }
    PcurveRates rates{PoleBounds(*first).first, 0.0};
    if (b.Degree < 2) {
        // A polyline curves nowhere within a knot span; its corners are at
        // the knots, where SamplingOf always places a point.
        return rates;
    }
    const std::optional<Hodograph> second = Differentiated(*first, b.Degree - 1);
    if (!second) {
        return std::nullopt;
    }
    rates.Curving = PoleBounds(*second).second;
    return rates;
}

/// Nothing for a rational B-spline or a type outside the prototype's.
[[nodiscard]] std::optional<PcurveRates> RatesOf(const Curve2d& pcurve)
{
    if (const Line2d* l = pcurve.Target<Line2d>()) {
        return PcurveRates{std::abs(l->Direction.X), 0.0};
    }
    if (const Circle2d* c = pcurve.Target<Circle2d>()) {
        return PcurveRates{c->Radius * (std::abs(c->XAxis.X) + std::abs(c->YAxis.X)),
                           c->Radius * (Length(c->XAxis) + Length(c->YAxis))};
    }
    if (const Ellipse2d* e = pcurve.Target<Ellipse2d>()) {
        return PcurveRates{e->MajorRadius * std::abs(e->XAxis.X) +
                               e->MinorRadius * std::abs(e->YAxis.X),
                           e->MajorRadius * Length(e->XAxis) + e->MinorRadius * Length(e->YAxis)};
    }
    if (const BSplineCurve2d* b = pcurve.Target<BSplineCurve2d>()) {
        return BSplineRates(*b);
    }
    return std::nullopt;
}

// --- The faces' boundary polygons ----------------------------------------------

/// Where an edge is sampled, in increasing order from one end of its range
/// to the other, and how far the chords between the samples may lie from
/// S(p(t)) on any face using it.
using Sampling = Detail::EdgeSampling;

/// The ends of the edge's range and every knot inside it of a B-spline
/// pcurve of any use: a B-spline is smooth only within a knot span, and at
/// a knot of a polyline it turns a corner the chords must keep.
[[nodiscard]] std::vector<double> BreaksOf(const std::vector<Curve2d>& pcurves,
                                           const Interval& range)
{
    std::vector<double> breaks{range.GetLower(), range.GetUpper()};
    for (const Curve2d& pcurve : pcurves) {
        if (const BSplineCurve2d* b = pcurve.Target<BSplineCurve2d>()) {
            for (const double knot : b->Knots) {
                if (knot > range.GetLower() && knot < range.GetUpper()) {
                    breaks.push_back(knot);
                }
            }
        }
    }
    std::ranges::sort(breaks);
    breaks.erase(std::ranges::unique(breaks).begin(), breaks.end());
    return breaks;
}

/// The sup over the edge's uses of |(S o p)''| (PcurveRates), and their
/// pcurves. NotSupported for a pcurve or surface outside the prototype's.
[[nodiscard]] Result<std::pair<double, std::vector<Curve2d>>> RateOf(const Shape& shape,
                                                                     const EdgeId edge)
{
    double rate = 0.0;
    std::vector<Curve2d> pcurves;
    STEREON_TRY(const RadialCoedgeRange uses, RadialCoedgesOf(shape, edge));
    for (const CoedgeId use : uses) {
        const std::optional<Curve2d> pcurve = PCurveOf(shape, use);
        const std::optional<PcurveRates> rates = pcurve ? RatesOf(*pcurve) : std::nullopt;
        const std::optional<SurfaceScale> scale =
            ScaleOf(SurfaceOf(shape, FaceOf(shape, LoopOf(shape, use))));
        if (!rates || !scale) {
            return Fail(ErrorCode::NotSupported,
                        "WindingNumber: a use of " + std::string(Describe(shape, edge)) +
                            " has a pcurve or a surface outside the prototype's lines, circles, "
                            "ellipses, non-rational B-splines, planes and cylinders");
        }
        rate = std::max(rate, scale->Radius * rates->Turn * rates->Turn +
                                  scale->Lipschitz() * rates->Curving);
        pcurves.push_back(*pcurve);
    }
    return std::pair{rate, std::move(pcurves)};
}

/// The edge's samples, the same for every use so the faces agree on the
/// points: every break (BreaksOf), and within each piece between breaks,
/// where every pcurve is smooth, enough uniform samples that a chord of
/// length h stays within h^2 rate / 8 <= target of S(p(t)). When the pieces
/// together would need more than MaxSegments, each takes its share of them
/// and the Deviation reports what that achieves, so the bound stays true.
[[nodiscard]] Result<Sampling> SamplingOf(const Shape& shape, const EdgeId edge,
                                          const double target)
{
    STEREON_TRY(const auto rated, RateOf(shape, edge));
    const auto& [rate, pcurves] = rated;
    const std::vector<double> breaks = BreaksOf(pcurves, RangeOf(shape, edge));
    std::vector<double> wanted; // segments each piece asks for
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double width = breaks[i + 1] - breaks[i];
        wanted.push_back(std::max(1.0, std::ceil(width * std::sqrt(rate / (8.0 * target)))));
        total += wanted.back();
    }
    const double share = std::min(1.0, MaxSegments / total);
    Sampling sampling;
    for (std::size_t i = 0; i < wanted.size(); ++i) {
        const double width = breaks[i + 1] - breaks[i];
        const double n = std::max(1.0, std::floor(wanted[i] * share));
        const double h = width / n;
        sampling.Deviation = std::max(sampling.Deviation, h * h * rate / 8.0);
        const auto count = static_cast<std::size_t>(n);
        for (std::size_t k = 0; k < count; ++k) {
            sampling.Parameters.push_back(breaks[i] + width * static_cast<double>(k) / n);
        }
    }
    sampling.Parameters.push_back(breaks.back());
    return sampling;
}

/// The use of the edge every face takes its points along the edge from: the
/// lowest coedge of the radial ring, so the choice is the same from every
/// face.
[[nodiscard]] Result<CoedgeId> AgreedUseOf(const Shape& shape, const EdgeId edge)
{
    STEREON_TRY(const RadialCoedgeRange uses, RadialCoedgesOf(shape, edge));
    std::optional<CoedgeId> agreed;
    for (const CoedgeId use : uses) {
        if (!agreed || use.GetIndex() < agreed->GetIndex()) {
            agreed = use;
        }
    }
    if (!agreed) {
        return Fail(ErrorCode::InvalidInput,
                    "WindingNumber: " + std::string(Describe(shape, edge)) + " has no use");
    }
    return *agreed;
}

/// Appends the coedge's corners along the loop's direction, all but its
/// last, which is the next coedge's first: the parameters are in the face's
/// own pcurve, the positions the edge's start vertex and then the agreed
/// use's surface at its pcurve, so that every face sharing the edge places
/// the same points and the triangles close up. Raises `tolerance` to the
/// edge's and the vertex's, which is how far those positions may lie from
/// this face's own boundary, and `deviation` to the chords' (SamplingOf).
/// Returns where the coedge ends in the face's parameter space, for the
/// loop's continuity check.
[[nodiscard]] Result<Point2> AppendCoedge(Polyline& polygon, double& tolerance, double& deviation,
                                          const Shape& shape, const CoedgeId coedge,
                                          const double target)
{
    const EdgeId edge = EdgeOf(shape, coedge);
    const std::optional<Curve2d> own = PCurveOf(shape, coedge);
    if (!own) {
        return Fail(ErrorCode::NotSupported,
                    "WindingNumber: " + std::string(Describe(shape, edge)) +
                        " has a use without a pcurve");
    }
    STEREON_TRY(const Sampling sampling, SamplingOf(shape, edge, target));
    STEREON_TRY(const CoedgeId agreed, AgreedUseOf(shape, edge));
    const Surface agreedSurface = SurfaceOf(shape, FaceOf(shape, LoopOf(shape, agreed)));
    const Curve2d agreedPcurve = *PCurveOf(shape, agreed); // SamplingOf saw every use has one
    const bool forward = StoredEdgeUseSense(shape, coedge) == EdgeUseSense::Same;
    const VertexId start = forward ? StartOf(shape, edge) : EndOf(shape, edge);
    const Interval range = RangeOf(shape, edge);
    tolerance = std::max({tolerance, ToleranceOf(shape, edge), ToleranceOf(shape, start)});
    deviation = std::max(deviation, sampling.Deviation);
    // The parameters run from the range's lower end to its upper end, both
    // exactly: a pcurve's domain may end where the edge's does.
    const std::vector<double>& ts = sampling.Parameters;
    const std::size_t segments = ts.size() - 1;
    for (std::size_t k = 0; k < segments; ++k) {
        const double t = ts[forward ? k : segments - k];
        const Point3 position =
            k == 0 ? PositionOf(shape, start) : Eval(agreedSurface, Eval(agreedPcurve, t));
        polygon.push_back({Eval(*own, t), position});
    }
    return Eval(*own, forward ? range.GetUpper() : range.GetLower());
}

/// A face's loops as polygons, their parameter box, how far the corners'
/// positions may lie from the face's own boundary, and how far the chords
/// between them may lie from its true edges.
struct FacePolygon {
    std::vector<Polyline> Loops;
    Box2 Box = Box2::Empty();
    double Tolerance = 0.0;
    double Deviation = 0.0;
};

/// Consecutive coedges of a loop must meet in the face's parameter space:
/// their pcurves end and start at one vertex, within the tolerances, so a
/// larger gap means the pcurves do not meet (a loop that jumps across a
/// periodic surface's seam, say) and the face has no polygon. InvalidResult
/// naming the two edges, so the shape's fault is not read as a side.
[[nodiscard]] Result<void> CheckMeeting(const Shape& shape, const FaceId face, const EdgeId from,
                                        const Point2& end, const EdgeId to, const Point2& start,
                                        const double tolerance, const SurfaceScale& scale)
{
    // Each pcurve reaches the vertex within the edge's and the vertex's
    // tolerance, so the two may miss each other in space by four of them.
    const double allowance = 4.0 * tolerance;
    if (scale.MinScale * Distance(end, start) <= allowance) {
        return {};
    }
    return Fail(ErrorCode::InvalidResult,
                std::format("WindingNumber: the loop of {} is not continuous in the face's "
                            "parameter space: {} ends at ({:.6g}, {:.6g}) and {} starts at "
                            "({:.6g}, {:.6g}), so the pcurves do not meet and the face cannot be "
                            "triangulated",
                            std::string(Describe(shape, face)), std::string(Describe(shape, from)),
                            end.X, end.Y, std::string(Describe(shape, to)), start.X, start.Y));
}

[[nodiscard]] Result<FacePolygon> PolygonOf(const Shape& shape, const FaceId face,
                                            const SurfaceScale& scale, const double target)
{
    FacePolygon polygon;
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        Polyline boundary;
        std::optional<std::pair<EdgeId, Point2>> previous; // the last coedge's edge and end
        STEREON_TRY(const CoedgeRange coedges, CoedgesOf(shape, loop));
        for (const CoedgeId coedge : coedges) {
            const std::size_t first = boundary.size();
            STEREON_TRY(const Point2 end, AppendCoedge(boundary, polygon.Tolerance,
                                                       polygon.Deviation, shape, coedge, target));
            if (previous) {
                STEREON_TRY_VOID(CheckMeeting(shape, face, previous->first, previous->second,
                                              EdgeOf(shape, coedge), boundary[first].Uv,
                                              polygon.Tolerance, scale));
            }
            previous = {EdgeOf(shape, coedge), end};
        }
        if (previous && !boundary.empty()) {
            STEREON_TRY_VOID(CheckMeeting(shape, face, previous->first, previous->second,
                                          EdgeOf(shape, *coedges.begin()), boundary.front().Uv,
                                          polygon.Tolerance, scale));
        }
        for (const Corner& corner : boundary) {
            polygon.Box = Hull(polygon.Box, corner.Uv);
        }
        polygon.Loops.push_back(std::move(boundary));
    }
    return polygon;
}

// --- Triangles -------------------------------------------------------------------

/// The solid angle the triangle abc subtends at p (Van Oosterom and
/// Strackee, 1983): positive when the triangle's right-hand normal points
/// away from p, so a body's outward faces add up to 4 pi at a point inside.
[[nodiscard]] double SolidAngle(const Point3& p, const Point3& a, const Point3& b, const Point3& c)
{
    const Vec3 ra = a - p;
    const Vec3 rb = b - p;
    const Vec3 rc = c - p;
    const double la = Length(ra);
    const double lb = Length(rb);
    const double lc = Length(rc);
    const double numerator = Dot(ra, Cross(rb, rc));
    const double denominator =
        la * lb * lc + Dot(ra, rb) * lc + Dot(ra, rc) * lb + Dot(rb, rc) * la;
    return 2.0 * std::atan2(numerator, denominator);
}

[[nodiscard]] double DistanceToSegment(const Point3& p, const Point3& a, const Point3& b)
{
    const Vec3 ab = b - a;
    const double ll = Dot(ab, ab);
    const double s = ll > 0.0 // STEREON-LINT-ALLOW(epsilon): a sign test; a point is its end
                         ? std::clamp(Dot(p - a, ab) / ll, 0.0, 1.0)
                         : 0.0;
    return Distance(p, a + s * ab);
}

/// The distance from p to the triangle abc: to its plane where the foot of
/// the perpendicular lies inside, else to the nearest edge.
[[nodiscard]] double DistanceToTriangle(const Point3& p, const Point3& a, const Point3& b,
                                        const Point3& c)
{
    const Vec3 n = Cross(b - a, c - a);
    const double nn = Dot(n, n);
    if (nn > 0.0) { // STEREON-LINT-ALLOW(epsilon): a sign test; a degenerate triangle is its edges
        // The foot lies inside when p is on the inner side of every edge's
        // plane through n; the component of p along n drops out of each.
        const double sides = std::min({Dot(Cross(b - a, p - a), n), Dot(Cross(c - b, p - b), n),
                                       Dot(Cross(a - c, p - c), n)});
        if (sides >= 0.0) { // STEREON-LINT-ALLOW(epsilon): a sign test
            return std::abs(Dot(p - a, n)) / std::sqrt(nn);
        }
    }
    return std::min(
        {DistanceToSegment(p, a, b), DistanceToSegment(p, b, c), DistanceToSegment(p, c, a)});
}

/// The fan of triangles from a polygon's first corner, which covers the
/// polygon's signed region once: its solid angle at p, the least distance
/// from p to any triangle, and their count.
struct Fan {
    double SolidAngle = 0.0;
    double Distance = std::numeric_limits<double>::infinity();
    std::size_t Triangles = 0;
};

[[nodiscard]] Fan FanOf(const Polyline& polygon, const Point3& p)
{
    Fan fan;
    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        const Point3& a = polygon.front().Position;
        const Point3& b = polygon[i].Position;
        const Point3& c = polygon[i + 1].Position;
        fan.SolidAngle += SolidAngle(p, a, b, c);
        fan.Distance = std::min(fan.Distance, DistanceToTriangle(p, a, b, c));
        ++fan.Triangles;
    }
    return fan;
}

// --- Cells of a face's parameter space ------------------------------------------

/// The polygon clipped to the half of the parameter plane where coordinate
/// `axis` (0 for u, 1 for v) lies above `at`, or below it
/// (Sutherland-Hodgman): the signed region it bounds is the polygon's region
/// within that half, joined by zero-width bridges where the cut parts it. A
/// new corner lies on the side it cuts, which the cell across that side
/// shares (a boundary chord with the face across the edge, an earlier cut
/// with the other half), so the triangles still close up into one
/// polyhedron; its offset from the surface is its ends' and the sagitta of
/// the side's turn.
[[nodiscard]] Polyline Clipped(const Polyline& polygon, const int axis, const double at,
                               const bool above, const double radius)
{
    const auto coordinate = [axis](const Corner& corner) {
        return axis == 0 ? corner.Uv.X : corner.Uv.Y;
    };
    const auto inside = [&](const Corner& corner) {
        return above ? coordinate(corner) >= at : coordinate(corner) <= at;
    };
    Polyline out;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Corner& a = polygon[i];
        const Corner& b = polygon[(i + 1) % polygon.size()];
        if (inside(a)) {
            out.push_back(a);
        }
        if (inside(a) != inside(b)) {
            const double s = (at - coordinate(a)) / (coordinate(b) - coordinate(a));
            const double off = std::max(a.Off, b.Off) + Sagitta(radius, std::abs(b.Uv.X - a.Uv.X));
            out.push_back(
                {a.Uv + s * (b.Uv - a.Uv), a.Position + s * (b.Position - a.Position), off});
        }
    }
    return out;
}

/// What the triangles of a face add up to, and how many of its cells could
/// not be separated from the point.
struct Tally {
    double SolidAngle = 0.0;
    std::size_t Triangles = 0;
    std::size_t Unseparated = 0;
};

/// What every cell of one face shares.
struct FaceFrame {
    Surface Geometry;
    SurfaceScale Scale;
    double Tolerance = 0.0;
    double Target = 0.0;
    double Chords = 0.0; // how far the boundary's chords lie from the edges
    Point3 Point;
};

/// How far the face may lie from the triangles of a cell turning du whose
/// corners lie up to `off` from the surface: a triangle strays from the
/// surface over its parameters by its corners' offsets and the sagitta of
/// its turn (taken three times, for the chord interpolation on a triangle
/// rather than a segment, with room to spare), the true boundary from the
/// chords by what the sampling achieved (SamplingOf: the target, or more
/// where MaxSegments capped it), and the positions from this face's own
/// boundary by the tolerances, at both ends of a chord and once more for a
/// cut.
[[nodiscard]] double DeviationOf(const FaceFrame& f, const double du, const double off)
{
    return off + 3.0 * Sagitta(f.Scale.Radius, du) + f.Chords + 4.0 * f.Tolerance;
}

/// The fans of all the loops together, and their corners' largest offset.
struct Cover {
    Fan Fans;
    double Off = 0.0;
};

[[nodiscard]] Cover CoverOf(const std::vector<Polyline>& loops, const Point3& p)
{
    Cover all;
    for (const Polyline& loop : loops) {
        const Fan fan = FanOf(loop, p);
        all.Fans.SolidAngle += fan.SolidAngle;
        all.Fans.Distance = std::min(all.Fans.Distance, fan.Distance);
        all.Fans.Triangles += fan.Triangles;
        for (const Corner& corner : loop) {
            all.Off = std::max(all.Off, corner.Off);
        }
    }
    return all;
}

/// The side to split a cell across: around a cylinder's axis first, while
/// the turn's sagitta is more than the target, so that cuts along the axis,
/// which turn as wide as their cell, stay short and their corners near the
/// surface; then across the longer side.
[[nodiscard]] int SplitAxisOf(const FaceFrame& f, const double du, const double dv)
{
    if (3.0 * Sagitta(f.Scale.Radius, du) > f.Target) {
        return 0;
    }
    return f.Scale.UScale * du >= f.Scale.VScale * dv ? 0 : 1;
}

/// One half of a cell, split at `at` across `axis`, with the loops clipped
/// to it; polygons left with fewer than three corners bound nothing.
struct Half {
    Box2 Cell;
    std::vector<Polyline> Loops;
};

[[nodiscard]] Half HalfOf(const Box2& cell, const std::vector<Polyline>& loops, const int axis,
                          const double at, const bool above, const double radius)
{
    Half half{cell, {}};
    const Interval& side = axis == 0 ? cell.X : cell.Y;
    (axis == 0 ? half.Cell.X : half.Cell.Y) =
        above ? Interval{at, side.GetUpper()} : Interval{side.GetLower(), at};
    for (const Polyline& loop : loops) {
        Polyline clipped = Clipped(loop, axis, at, above, radius);
        if (clipped.size() >= 3) {
            half.Loops.push_back(std::move(clipped));
        }
    }
    return half;
}

/// Sums the cell's triangles into the tally when the point is separated from
/// every one of them by more than the face can stray, or when the cell is
/// already as small as the target (then as unseparated); otherwise splits
/// the cell across its longer side and clips the loops to the halves.
void Accumulate(Tally& tally, const FaceFrame& f, const Box2& cell,
                const std::vector<Polyline>& loops, const int depth)
{
    const Cover all = CoverOf(loops, f.Point);
    if (all.Fans.Triangles == 0) {
        return;
    }
    const double du = Width(cell.X);
    const double dv = Width(cell.Y);
    const bool separated = all.Fans.Distance > DeviationOf(f, du, all.Off);
    const bool smallest = f.Scale.Lipschitz() * std::max(du, dv) <= f.Target || depth >= MaxDepth;
    if (separated || smallest) {
        tally.SolidAngle += all.Fans.SolidAngle;
        tally.Triangles += all.Fans.Triangles;
        tally.Unseparated += separated ? 0 : 1;
        return;
    }
    const int axis = SplitAxisOf(f, du, dv);
    const double at = Midpoint(axis == 0 ? cell.X : cell.Y);
    for (const bool above : {false, true}) {
        const Half half = HalfOf(cell, loops, axis, at, above, f.Scale.Radius);
        Accumulate(tally, f, half.Cell, half.Loops, depth + 1);
    }
}

[[nodiscard]] Result<void> AccumulateFace(Tally& tally, const Shape& body, const FaceId face,
                                          const Point3& point, const double target)
{
    const Surface surface = SurfaceOf(body, face);
    const std::optional<SurfaceScale> scale = ScaleOf(surface);
    if (!scale) {
        return Fail(ErrorCode::NotSupported,
                    "WindingNumber: " + std::string(Describe(body, face)) +
                        " lies on a surface outside the prototype's planes and cylinders");
    }
    STEREON_TRY(const FacePolygon polygon, PolygonOf(body, face, *scale, target));
    // A loop runs counter-clockwise about the face's own normal, reversed or
    // not, so its triangles, made from its corners in space, face the way
    // the face does: the orientation is in the corners' order, and no sign
    // for a face oriented against its surface may be applied on top of it.
    const FaceFrame frame{surface, *scale, polygon.Tolerance, target, polygon.Deviation, point};
    Accumulate(tally, frame, polygon.Box, polygon.Loops, 0);
    return {};
}

// --- The body -----------------------------------------------------------------------

/// The diagonal of the body's vertices' bounding box.
[[nodiscard]] Result<double> ExtentOf(const Shape& body)
{
    Box3 box = Box3::Empty();
    STEREON_TRY(const auto vertices, AllEntities<EntityKind::Vertex>(body));
    for (const VertexId vertex : vertices) {
        box = Hull(box, PositionOf(body, vertex));
    }
    if (box.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "WindingNumber: the body has no vertices");
    }
    return Length(Diagonal(box));
}

/// A decided number must be 0 or 1: anything else means the faces do not
/// close up, or one is oriented against its loops.
[[nodiscard]] Result<void> CheckWhole(const WindingNumber& winding)
{
    const double nearest = std::round(winding.Value);
    const bool whole = std::abs(winding.Value - nearest) <= winding.Bound;
    if (!whole || nearest < 0.0 || nearest > 1.0) { // STEREON-LINT-ALLOW(epsilon): whole numbers
        return Fail(ErrorCode::InvalidResult,
                    std::format("WindingNumber: the body's faces add up to {:.6g} ± {:.3g} at the "
                                "point, not 0 or 1: they do not close up, or a face is oriented "
                                "against its loops",
                                winding.Value, winding.Bound));
    }
    return {};
}

} // namespace

Result<Detail::EdgeSampling> Detail::WindingSamplingOf(const Shape& shape, const EdgeId edge,
                                                       const double target)
{
    return SamplingOf(shape, edge, target);
}

bool WindingNumber::IsDecided() const noexcept
{
    return std::abs(Value - 0.5) > Bound;
}

bool WindingNumber::IsInside() const noexcept
{
    return Value > 0.5; // STEREON-LINT-ALLOW(epsilon): the side of one half, not a tolerance
}

Result<WindingNumber> WindingNumberOf(const Shape& body, const Point3& point, const int refinements,
                                      const Context& ctx)
{
    if (refinements < 0) {
        return Fail(ErrorCode::InvalidInput, "WindingNumber: a negative refinement budget");
    }
    STEREON_TRY(const double extent, ExtentOf(body));
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(body));
    WindingNumber winding;
    for (int round = 0; round <= refinements; ++round) {
        const double target =
            std::max(extent / FirstTargetDivisor / std::pow(4.0, round), ctx.GetResolution());
        Tally tally;
        for (const FaceId face : faces) {
            STEREON_TRY_VOID(AccumulateFace(tally, body, face, point, target));
        }
        winding = {tally.SolidAngle / (4.0 * std::numbers::pi),
                   static_cast<double>(tally.Unseparated) +
                       static_cast<double>(tally.Triangles) * RoundingPerTriangle,
                   round, tally.Triangles};
        if (winding.IsDecided()) {
            STEREON_TRY_VOID(CheckWhole(winding));
            return winding;
        }
        if (target <= ctx.GetResolution()) {
            break; // Another round would be the same.
        }
    }
    return winding;
}

} // namespace Stereon
