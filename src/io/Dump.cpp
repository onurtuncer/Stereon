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
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

module Stereon.Io;

import :Debug;
import :Dump;
import :GeometryArchive;

namespace Stereon {

// Instantiated at namespace scope before any function uses it
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

using DumpResult = Result<std::optional<std::filesystem::path>>;

[[nodiscard]] bool IsNameCharacter(const char c) noexcept
{
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' ||
           c == '.' || c == '-';
}

[[nodiscard]] std::int64_t MillisecondsNow()
{
    using namespace std::chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

[[nodiscard]] DebugColour Rgb(const std::uint8_t r, const std::uint8_t g, const std::uint8_t b)
{
    DebugColour c;
    c.R = r;
    c.G = g;
    c.B = b;
    c.A = 255;
    return c;
}

[[nodiscard]] GeometryArchive ArchiveOf(const Context& context)
{
    GeometryArchive archive;
    archive.Policy = PolicyOf(context);
    return archive;
}

/// The edge drawn for display: its 3D curve over its range or, for a kind
/// without one (a sheet boundary, a seam, a pole), its primary coedge's
/// pcurve over the range, mapped onto the face's surface. Empty for an edge
/// with neither, or one sampling refuses.
[[nodiscard]] std::vector<Point3> EdgePoints(const Shape& shape, const EdgeId edge)
{
    const Interval range = RangeOf(shape, edge);
    if (const std::optional<Curve> curve = CurveOf(shape, edge)) {
        const Result<Polyline> line = SampleCurve(*curve, DisplayTolerance{}, range);
        return line ? line->Points : std::vector<Point3>{};
    }
    const CoedgeId coedge = PrimaryCoedgeOf(shape, edge);
    const std::optional<Curve2d> pcurve = coedge.IsValid() ? PCurveOf(shape, coedge) : std::nullopt;
    if (!pcurve) {
        return {};
    }
    const Result<Polyline2d> line = SampleCurve(*pcurve, DisplayTolerance{}, range);
    if (!line) {
        return {};
    }
    const Surface& surface = SurfaceOf(shape, FaceOf(shape, LoopOf(shape, coedge)));
    std::vector<Point3> points;
    points.reserve(line->Points.size());
    for (const Point2& uv : line->Points) {
        points.push_back(Eval(surface, uv));
    }
    return points;
}

[[nodiscard]] Result<void> AddVertices(const Shape& shape, DebugChunk& debug)
{
    STEREON_TRY(const auto vertices, AllEntities<EntityKind::Vertex>(shape));
    for (const VertexId v : vertices) {
        const Point3 p = PositionOf(shape, v);
        const double tolerance = ToleranceOf(shape, v);
        DebugPoint point;
        point.Position = p;
        debug.Points.push_back(point);
        DebugBox box;
        box.Box = Box3{Interval{p.X - tolerance, p.X + tolerance},
                       Interval{p.Y - tolerance, p.Y + tolerance},
                       Interval{p.Z - tolerance, p.Z + tolerance}};
        debug.Boxes.push_back(box);
        DebugLabel label;
        label.Position = p;
        label.Text = Describe(shape, v);
        debug.Labels.push_back(std::move(label));
    }
    return {};
}

/// An edge's curve trimmed to its range, or the curve itself when the trim
/// is refused (a range that is the whole of a periodic domain).
[[nodiscard]] Curve EdgeGeometry(const Curve& curve, const Interval range)
{
    Result<TrimmedCurve> trimmed = MakeTrimmedCurve(curve, range);
    return trimmed ? Curve(*std::move(trimmed)) : curve;
}

[[nodiscard]] Result<void> AddEdges(const Shape& shape, GeometryArchive& archive)
{
    STEREON_TRY(auto edges, AllEntities<EntityKind::Edge>(shape));
    for (const EdgeId e : edges) {
        if (const std::optional<Curve> curve = CurveOf(shape, e)) {
            archive.Curves.push_back(EdgeGeometry(*curve, RangeOf(shape, e)));
        }
        std::vector<Point3> points = EdgePoints(shape, e);
        DebugLabel label;
        label.Position =
            points.empty() ? PositionOf(shape, StartOf(shape, e)) : points[points.size() / 2];
        label.Text = Describe(shape, e);
        archive.Debug->Labels.push_back(std::move(label));
        if (!points.empty()) {
            DebugPolyline line;
            line.Points = std::move(points);
            line.Colour = EdgeCurveColour(KindOf(shape, e));
            archive.Debug->Polylines.push_back(std::move(line));
        }
    }
    return {};
}

/// Widens [lower, upper] to hold the points.
void Extend(Point2& lower, Point2& upper, const std::vector<Point2>& points)
{
    for (const Point2& p : points) {
        lower = Point2{std::min(lower.X, p.X), std::min(lower.Y, p.Y)};
        upper = Point2{std::max(upper.X, p.X), std::max(upper.Y, p.Y)};
    }
}

/// The parameter box the face's boundary covers: every pcurve of its loops,
/// sampled over its edge's range. Nothing for a face without pcurves or one
/// whose box is flat.
[[nodiscard]] std::optional<Box2> FaceBox(const Shape& shape, const FaceId face)
{
    const auto loops = LoopsOf(shape, face);
    if (!loops) {
        return std::nullopt;
    }
    constexpr double inf = std::numeric_limits<double>::infinity();
    Point2 lower{inf, inf};
    Point2 upper{-inf, -inf};
    for (const LoopId loop : *loops) {
        const auto coedges = CoedgesOf(shape, loop);
        if (!coedges) {
            continue;
        }
        for (const CoedgeId c : *coedges) {
            const std::optional<Curve2d> pcurve = PCurveOf(shape, c);
            if (!pcurve) {
                continue;
            }
            const Result<Polyline2d> line =
                SampleCurve(*pcurve, DisplayTolerance{}, RangeOf(shape, EdgeOf(shape, c)));
            if (line) {
                Extend(lower, upper, line->Points);
            }
        }
    }
    if (!(lower.X < upper.X && lower.Y < upper.Y)) {
        return std::nullopt;
    }
    return Box2{Interval{lower.X, upper.X}, Interval{lower.Y, upper.Y}};
}

/// The face's surface trimmed to the box its boundary covers, so that a
/// plane or a cylinder is drawn where the face is; the surface itself when
/// there is no such box or the trim is refused.
[[nodiscard]] Surface FaceGeometry(const Shape& shape, const FaceId face)
{
    const Surface& surface = SurfaceOf(shape, face);
    const std::optional<Box2> box = FaceBox(shape, face);
    if (!box) {
        return surface;
    }
    Result<TrimmedSurface> trimmed = MakeTrimmedSurface(surface, *box);
    return trimmed ? Surface(*std::move(trimmed)) : surface;
}

[[nodiscard]] Result<void> AddFacesAndCoedges(const Shape& shape, GeometryArchive& archive)
{
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(shape));
    for (const FaceId f : faces) {
        archive.Surfaces.push_back(FaceGeometry(shape, f));
    }
    STEREON_TRY(const auto coedges, AllEntities<EntityKind::Coedge>(shape));
    for (const CoedgeId c : coedges) {
        if (const std::optional<Curve2d> pcurve = PCurveOf(shape, c)) {
            archive.Curves2d.push_back(*pcurve);
        }
    }
    return {};
}

/// Bool for the debugger entry points.
[[nodiscard]] bool Written(const DumpResult& result)
{
    return result && result->has_value();
}

} // namespace

std::optional<std::filesystem::path> DumpDirectory()
{
    const char* value = std::getenv(std::string(DumpDirectoryVariable).c_str());
    if (value == nullptr || *value == '\0') {
        return std::nullopt;
    }
    return std::filesystem::path(value);
}

bool IsDumpName(const std::string_view name) noexcept
{
    if (name.empty() || name.front() == '.') {
        return false;
    }
    for (const char c : name) {
        if (!IsNameCharacter(c)) {
            return false;
        }
    }
    return true;
}

DumpResult DumpAt(const std::string_view name, GeometryArchive archive, const SourceLocation where)
{
    if (!IsDumpName(name)) {
        return Fail(ErrorCode::InvalidInput,
                    "Dump: \"" + std::string(name) +
                        "\" is not a dump name; use letters, digits, '_', '.' and '-', not "
                        "starting with '.'");
    }
    const std::optional<std::filesystem::path> directory = DumpDirectory();
    if (!directory) {
        return std::optional<std::filesystem::path>{};
    }
    if (!archive.Debug) {
        archive.Debug.emplace();
    }
    archive.Debug->Name = std::string(name);
    archive.Debug->File = where.FileName;
    archive.Debug->Line = static_cast<std::uint32_t>(where.Line);
    archive.Debug->Time = MillisecondsNow();
    std::error_code error;
    std::filesystem::create_directories(*directory, error);
    if (error) {
        return Fail(ErrorCode::InvalidInput,
                    "Dump: cannot create \"" + directory->string() + "\": " + error.message());
    }
    const std::filesystem::path path = *directory / (std::string(name) + ".stn");
    STEREON_TRY_VOID(SaveGeometry(path, archive));
    return std::optional<std::filesystem::path>{path};
}

DumpResult DumpAt(const std::string_view name, const Curve& curve, const SourceLocation where)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Curves.push_back(curve);
    return DumpAt(name, std::move(archive), where);
}

DumpResult DumpAt(const std::string_view name, const Curve2d& curve, const SourceLocation where)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Curves2d.push_back(curve);
    return DumpAt(name, std::move(archive), where);
}

DumpResult DumpAt(const std::string_view name, const Surface& surface, const SourceLocation where)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Surfaces.push_back(surface);
    return DumpAt(name, std::move(archive), where);
}

DumpResult DumpAt(const std::string_view name, const Shape& shape, const SourceLocation where)
{
    // Checked before the shape is walked, so that a dump switched off costs
    // the environment lookup only.
    if (IsDumpName(name) && !DumpDirectory()) {
        return std::optional<std::filesystem::path>{};
    }
    STEREON_TRY(GeometryArchive archive, ShapeArchive(shape));
    return DumpAt(name, std::move(archive), where);
}

DumpResult DumpAt(const std::string_view name, DebugChunk primitives, const SourceLocation where)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Debug = std::move(primitives);
    return DumpAt(name, std::move(archive), where);
}

Result<GeometryArchive> ShapeArchive(const Shape& shape)
{
    if (shape.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "ShapeArchive: the Shape is empty");
    }
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Debug.emplace();
    STEREON_TRY_VOID(AddFacesAndCoedges(shape, archive));
    STEREON_TRY_VOID(AddEdges(shape, archive));
    STEREON_TRY_VOID(AddVertices(shape, *archive.Debug));
    return archive;
}

DebugColour EdgeCurveColour(const EdgeCurveKind kind) noexcept
{
    switch (kind) {
    case EdgeCurveKind::Intersection:
        return Rgb(220, 60, 60);
    case EdgeCurveKind::Exact:
        return Rgb(60, 110, 220);
    case EdgeCurveKind::Trim:
        return Rgb(60, 170, 80);
    case EdgeCurveKind::Boundary:
        return Rgb(230, 150, 30);
    case EdgeCurveKind::Spatial:
        return Rgb(190, 70, 200);
    case EdgeCurveKind::Degenerate:
        return Rgb(140, 140, 140);
    }
    return Rgb(0, 0, 0);
}

bool Dump(const char* name, const Curve& curve)
{
    return name != nullptr && Written(DumpAt(name, curve, SourceLocation{}));
}

bool Dump(const char* name, const Curve2d& curve)
{
    return name != nullptr && Written(DumpAt(name, curve, SourceLocation{}));
}

bool Dump(const char* name, const Surface& surface)
{
    return name != nullptr && Written(DumpAt(name, surface, SourceLocation{}));
}

bool Dump(const char* name, const Shape& shape)
{
    return name != nullptr && Written(DumpAt(name, shape, SourceLocation{}));
}

} // namespace Stereon
