// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Viewer.Scene;

namespace Stereon::Viewer {

// Instantiated at namespace scope before any function uses them
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));
static_assert(sizeof(Vec2) == 2 * sizeof(double));

namespace {

constexpr Rgba ControlColour{0.55F, 0.55F, 0.6F, 1.0F};
constexpr Rgba KnotColour{1.0F, 0.8F, 0.1F, 1.0F};
constexpr Rgba TickColour{0.85F, 0.85F, 0.85F, 1.0F};

/// The arrow's length as a share of the curve's drawn extent.
constexpr double ArrowShare = 0.06;
/// Where along its range a curve's arrow sits.
constexpr double ArrowAt = 0.6;

[[nodiscard]] Rgba FromDebug(const DebugColour& c, const Rgba fallback)
{
    if (!c.IsSet()) {
        return fallback;
    }
    constexpr float scale = 1.0F / 255.0F;
    return Rgba{static_cast<float>(c.R) * scale, static_cast<float>(c.G) * scale,
                static_cast<float>(c.B) * scale, static_cast<float>(c.A) * scale};
}

/// A unit vector perpendicular to the unit vector d.
[[nodiscard]] Vec3 Perpendicular(const Vec3& d)
{
    const Vec3 axis = std::abs(d.X) < std::abs(d.Y)
                          ? (std::abs(d.X) < std::abs(d.Z) ? Vec3{1, 0, 0} : Vec3{0, 0, 1})
                          : (std::abs(d.Y) < std::abs(d.Z) ? Vec3{0, 1, 0} : Vec3{0, 0, 1});
    const Vec3 n = Cross(d, axis);
    return (1.0 / Length(n)) * n;
}

/// `range` where it is bounded; around `centre`, `extent` (in model units
/// along a tangent of length `speed`) each way where it is not.
[[nodiscard]] Interval Around(const Interval range, const double centre, const double extent,
                              const double speed)
{
    if (range.IsBounded()) {
        return range;
    }
    const double reach = extent / std::max(speed, 1e-12);
    const double lower = std::max(range.GetLower(), centre - reach);
    const double upper = std::min(range.GetUpper(), centre + reach);
    return lower < upper ? Interval{lower, upper} : Interval{centre - reach, centre + reach};
}

[[nodiscard]] double DiagonalOf(const std::vector<Point3>& points)
{
    if (points.empty()) {
        return 0.0;
    }
    Point3 lo = points.front();
    Point3 hi = points.front();
    for (const Point3& p : points) {
        lo = Point3{std::min(lo.X, p.X), std::min(lo.Y, p.Y), std::min(lo.Z, p.Z)};
        hi = Point3{std::max(hi.X, p.X), std::max(hi.Y, p.Y), std::max(hi.Z, p.Z)};
    }
    return Length(hi - lo);
}

class Builder {
public:
    Builder(const SceneOptions& options, const Rgba colour,
            const std::optional<Placement> placement)
        : m_Options(options)
        , m_Colour(colour)
        , m_Placement(placement)
    {}

    void AddCurve(const Curve& curve)
    {
        const std::optional<Interval> range = CurveRange(curve);
        if (!range) {
            return;
        }
        const Result<Polyline> line = SampleCurve(curve, m_Options.Tolerance, *range);
        if (!line) {
            m_Out.Warnings.push_back("a curve: " + std::string(line.error().GetMessage()));
            return;
        }
        AddPolyline(line->Points, Layer::Geometry, m_Colour);
        AddTicks(curve, *range);
        AddArrow(curve, *range, ArrowShare * DiagonalOf(line->Points));
        if (const auto* bspline = curve.Target<BSplineCurve>()) {
            AddPolyline(bspline->Poles, Layer::ControlPolygon, ControlColour);
        }
        const std::vector<double> knots = Breakpoints(curve, *range);
        if (knots.size() > 2 || curve.Target<BSplineCurve>() != nullptr) {
            for (const double k : knots) {
                m_Out.Markers.push_back(Marker{Eval(curve, k), KnotColour, Layer::Knots});
            }
        }
    }

    void AddCurve2d(const Curve2d& curve)
    {
        const Interval domain = Domain(curve);
        const double centre = Project(curve, Point2{});
        const Interval range =
            Around(domain, centre, m_Options.DefaultExtent, Length(Deriv(curve, centre, 1)));
        const Result<Polyline2d> line = SampleCurve(curve, m_Options.Tolerance, range);
        if (!line) {
            m_Out.Warnings.push_back("a plane curve: " + std::string(line.error().GetMessage()));
            return;
        }
        for (std::size_t i = 0; i + 1 < line->Points.size(); ++i) {
            m_Out.PlaneSegments.push_back(
                PlaneSegment{line->Points[i], line->Points[i + 1], m_Colour, Layer::Geometry});
        }
        AddArrow2d(curve, range);
    }

    void AddSurface(const Surface& surface)
    {
        const std::optional<Box2> range = SurfaceRange(surface);
        if (!range) {
            return;
        }
        const Result<std::vector<IsoLine>> lines =
            SampleSurface(surface, m_Options.Tolerance, *range, m_Options.IsoLines);
        if (!lines) {
            m_Out.Warnings.push_back("a surface: " + std::string(lines.error().GetMessage()));
            return;
        }
        for (const IsoLine& iso : *lines) {
            AddPolyline(iso.Line.Points, Layer::Geometry, m_Colour);
        }
        if (const auto* bspline = surface.Target<BSplineSurface>()) {
            AddControlNet(*bspline);
        }
    }

    void AddDebug(const DebugChunk& debug)
    {
        for (const DebugPoint& p : debug.Points) {
            m_Out.Markers.push_back(
                Marker{p.Position, FromDebug(p.Colour, m_Colour), Layer::Debug});
        }
        for (const DebugPolyline& p : debug.Polylines) {
            AddPolyline(p.Points, Layer::Debug, FromDebug(p.Colour, m_Colour));
        }
        for (const DebugBox& b : debug.Boxes) {
            AddBox(b.Box, FromDebug(b.Colour, m_Colour));
        }
        for (const DebugTriangles& t : debug.Triangles) {
            const Rgba colour = FromDebug(t.Colour, m_Colour);
            for (std::size_t i = 0; i + 2 < t.Indices.size(); i += 3) {
                m_Out.Triangles.push_back(Triangle{t.Vertices[t.Indices[i]],
                                                   t.Vertices[t.Indices[i + 1]],
                                                   t.Vertices[t.Indices[i + 2]], colour});
            }
        }
        for (const DebugLabel& l : debug.Labels) {
            m_Out.Labels.push_back(Label{l.Position, l.Text, FromDebug(l.Colour, m_Colour)});
        }
    }

    [[nodiscard]] Primitives Take() { return std::move(m_Out); }

private:
    /// The curve's domain, or a range around the placement for an unbounded
    /// one; nothing (and HasUnbounded) without a placement.
    [[nodiscard]] std::optional<Interval> CurveRange(const Curve& curve)
    {
        const Interval domain = Domain(curve);
        if (domain.IsBounded()) {
            return domain;
        }
        m_Out.HasUnbounded = true;
        if (!m_Placement) {
            return std::nullopt;
        }
        const double centre = Project(curve, m_Placement->Centre);
        return Around(domain, centre, m_Placement->Extent, Length(Deriv(curve, centre, 1)));
    }

    [[nodiscard]] std::optional<Box2> SurfaceRange(const Surface& surface)
    {
        const Box2 domain = Domain(surface);
        if (domain.X.IsBounded() && domain.Y.IsBounded()) {
            return domain;
        }
        m_Out.HasUnbounded = true;
        if (!m_Placement) {
            return std::nullopt;
        }
        const Point2 centre = Project(surface, m_Placement->Centre);
        const double extent = m_Placement->Extent;
        return Box2{Around(domain.X, centre.X, extent, Length(Deriv(surface, centre, 1, 0))),
                    Around(domain.Y, centre.Y, extent, Length(Deriv(surface, centre, 0, 1)))};
    }

    void AddPolyline(const std::vector<Point3>& points, const Layer layer, const Rgba colour)
    {
        for (std::size_t i = 0; i + 1 < points.size(); ++i) {
            m_Out.Segments.push_back(Segment{points[i], points[i + 1], colour, layer});
        }
    }

    void AddTicks(const Curve& curve, const Interval range)
    {
        const std::size_t n = m_Options.Ticks;
        for (std::size_t i = 0; n > 1 && i < n; ++i) {
            const double t = range.GetLower() + (range.GetUpper() - range.GetLower()) *
                                                    static_cast<double>(i) /
                                                    static_cast<double>(n - 1);
            m_Out.Markers.push_back(Marker{Eval(curve, t), TickColour, Layer::ParameterTicks});
        }
    }

    void AddArrow(const Curve& curve, const Interval range, const double size)
    {
        const double t = range.GetLower() + ArrowAt * (range.GetUpper() - range.GetLower());
        const Vec3 tangent = Deriv(curve, t, 1);
        const double speed = Length(tangent);
        if (!(speed > 0.0) || !(size > 0.0)) {
            return;
        }
        const Vec3 d = (1.0 / speed) * tangent;
        const Vec3 n = Perpendicular(d);
        const Point3 tail = Eval(curve, t);
        const Point3 head = tail + size * d;
        const Point3 back = head - (0.4 * size) * d;
        for (const Point3& to : {tail, back + (0.25 * size) * n, back - (0.25 * size) * n}) {
            m_Out.Segments.push_back(Segment{head, to, m_Colour, Layer::Directions});
        }
    }

    void AddArrow2d(const Curve2d& curve, const Interval range)
    {
        const double t = range.GetLower() + ArrowAt * (range.GetUpper() - range.GetLower());
        const Vec2 tangent = Deriv(curve, t, 1);
        const double speed = Length(tangent);
        const double size =
            ArrowShare * Length(Eval(curve, range.GetUpper()) - Eval(curve, range.GetLower()));
        if (!(speed > 0.0) || !(size > 0.0)) {
            return;
        }
        const Vec2 d = (1.0 / speed) * tangent;
        const Vec2 n{-d.Y, d.X};
        const Point2 head = Eval(curve, t);
        const Point2 back = head - (0.4 * size) * d;
        for (const Point2& to : {back + (0.25 * size) * n, back - (0.25 * size) * n}) {
            m_Out.PlaneSegments.push_back(PlaneSegment{head, to, m_Colour, Layer::Directions});
        }
    }

    void AddControlNet(const BSplineSurface& s)
    {
        const std::size_t nu = PoleCountU(s);
        const std::size_t nv = PoleCountV(s);
        for (std::size_t i = 0; i < nu; ++i) {
            for (std::size_t j = 0; j < nv; ++j) {
                const Point3& p = s.Poles[i * nv + j];
                if (j + 1 < nv) {
                    m_Out.Segments.push_back(
                        Segment{p, s.Poles[i * nv + j + 1], ControlColour, Layer::ControlPolygon});
                }
                if (i + 1 < nu) {
                    m_Out.Segments.push_back(Segment{p, s.Poles[(i + 1) * nv + j], ControlColour,
                                                     Layer::ControlPolygon});
                }
            }
        }
    }

    void AddBox(const Box3& b, const Rgba colour)
    {
        const std::array<double, 2> xs{b.X.GetLower(), b.X.GetUpper()};
        const std::array<double, 2> ys{b.Y.GetLower(), b.Y.GetUpper()};
        const std::array<double, 2> zs{b.Z.GetLower(), b.Z.GetUpper()};
        const auto corner = [&](const std::size_t k) {
            return Point3{xs[k & 1U], ys[(k >> 1U) & 1U], zs[(k >> 2U) & 1U]};
        };
        // The twelve edges join corners that differ in one coordinate bit.
        for (std::size_t k = 0; k < 8; ++k) {
            for (const std::size_t bit : {std::size_t{1}, std::size_t{2}, std::size_t{4}}) {
                if ((k & bit) == 0) {
                    m_Out.Segments.push_back(
                        Segment{corner(k), corner(k | bit), colour, Layer::Debug});
                }
            }
        }
    }

    const SceneOptions& m_Options;
    Rgba m_Colour;
    std::optional<Placement> m_Placement;
    Primitives m_Out;
};

void Extend(Point3& lo, Point3& hi, const Point3& p)
{
    lo = Point3{std::min(lo.X, p.X), std::min(lo.Y, p.Y), std::min(lo.Z, p.Z)};
    hi = Point3{std::max(hi.X, p.X), std::max(hi.Y, p.Y), std::max(hi.Z, p.Z)};
}

} // namespace

Primitives BuildPrimitives(const GeometryArchive& archive, const SceneOptions& options,
                           const Rgba colour, const std::optional<Placement> placement)
{
    Builder builder(options, colour, placement);
    for (const Curve& curve : archive.Curves) {
        builder.AddCurve(curve);
    }
    for (const Curve2d& curve : archive.Curves2d) {
        builder.AddCurve2d(curve);
    }
    for (const Surface& surface : archive.Surfaces) {
        builder.AddSurface(surface);
    }
    if (archive.Debug) {
        builder.AddDebug(*archive.Debug);
    }
    return builder.Take();
}

Box3 BoundsOf(const Primitives& primitives)
{
    std::vector<Point3> points;
    for (const Segment& s : primitives.Segments) {
        points.push_back(s.A);
        points.push_back(s.B);
    }
    for (const Marker& m : primitives.Markers) {
        points.push_back(m.Position);
    }
    for (const Triangle& t : primitives.Triangles) {
        points.insert(points.end(), {t.A, t.B, t.C});
    }
    for (const Label& l : primitives.Labels) {
        points.push_back(l.Position);
    }
    if (points.empty()) {
        return Box3::Empty();
    }
    Point3 lo = points.front();
    Point3 hi = points.front();
    for (const Point3& p : points) {
        Extend(lo, hi, p);
    }
    return Box3{Interval{lo.X, hi.X}, Interval{lo.Y, hi.Y}, Interval{lo.Z, hi.Z}};
}

Rgba ColourFor(const std::string_view name)
{
    // FNV-1a, then the hash's low bits as a hue on a fixed saturation and
    // value, so that every name gets a clear colour and keeps it.
    std::uint32_t hash = 2166136261U;
    for (const char c : name) {
        hash = (hash ^ static_cast<std::uint8_t>(c)) * 16777619U;
    }
    const float hue = static_cast<float>(hash % 360U) / 60.0F;
    const auto sector = static_cast<int>(hue);
    const float f = hue - static_cast<float>(sector);
    constexpr float v = 0.95F;
    constexpr float p = v * 0.35F;
    const float q = v * (1.0F - 0.65F * f);
    const float t = v * (1.0F - 0.65F * (1.0F - f));
    const std::array<Rgba, 6> sectors{Rgba{v, t, p, 1.0F}, Rgba{q, v, p, 1.0F},
                                      Rgba{p, v, t, 1.0F}, Rgba{p, q, v, 1.0F},
                                      Rgba{t, p, v, 1.0F}, Rgba{v, p, q, 1.0F}};
    return sectors[static_cast<std::size_t>(sector % 6)];
}

} // namespace Stereon::Viewer
