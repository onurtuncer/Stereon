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
#include <charconv>
#include <chrono>
#include <cstddef>
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

/// The face's loops as closed polygons in its surface's parameters: each
/// coedge's pcurve sampled over its edge's range, turned round for a use
/// against its edge, in the loop's order. Nothing when a use has no pcurve
/// or cannot be sampled.
[[nodiscard]] std::optional<std::vector<std::vector<Point2>>> FaceLoops(const Shape& shape,
                                                                        const FaceId face)
{
    const auto loops = LoopsOf(shape, face);
    if (!loops) {
        return std::nullopt;
    }
    std::vector<std::vector<Point2>> out;
    for (const LoopId loop : *loops) {
        const auto coedges = CoedgesOf(shape, loop);
        if (!coedges) {
            return std::nullopt;
        }
        std::vector<Point2> polygon;
        for (const CoedgeId c : *coedges) {
            const std::optional<Curve2d> pcurve = PCurveOf(shape, c);
            if (!pcurve) {
                return std::nullopt;
            }
            Result<Polyline2d> line =
                SampleCurve(*pcurve, DisplayTolerance{}, RangeOf(shape, EdgeOf(shape, c)));
            if (!line) {
                return std::nullopt;
            }
            if (StoredEdgeUseSense(shape, c) == EdgeUseSense::Opposite) {
                std::ranges::reverse(line->Points);
            }
            // Each use starts where the last ended: drop the repeated point.
            const std::size_t skip = polygon.empty() ? 0 : 1;
            polygon.insert(polygon.end(), line->Points.begin() + static_cast<std::ptrdiff_t>(skip),
                           line->Points.end());
        }
        if (polygon.size() >= 3) {
            out.push_back(std::move(polygon));
        }
    }
    return out;
}

/// The face drawn where it is: its surface's iso-lines inside its loops,
/// in FaceColour, as polylines of the debug chunk. False when the face's
/// loops or iso-lines cannot be sampled; the caller then draws the surface
/// over the box the boundary covers instead.
[[nodiscard]] bool AddFaceLines(const Shape& shape, const FaceId face, DebugChunk& debug,
                                const DebugColour colour = FaceColour())
{
    const std::optional<std::vector<std::vector<Point2>>> loops = FaceLoops(shape, face);
    const std::optional<Box2> box = FaceBox(shape, face);
    if (!loops || loops->empty() || !box) {
        return false;
    }
    const Result<std::vector<IsoLine>> lines =
        SampleTrimmedSurface(SurfaceOf(shape, face), DisplayTolerance{}, *loops, *box);
    if (!lines) {
        return false;
    }
    for (const IsoLine& iso : *lines) {
        DebugPolyline line;
        line.Points = iso.Line.Points;
        line.Colour = colour;
        debug.Polylines.push_back(std::move(line));
    }
    return true;
}

[[nodiscard]] Result<void> AddFacesAndCoedges(const Shape& shape, GeometryArchive& archive)
{
    STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(shape));
    for (const FaceId f : faces) {
        if (!AddFaceLines(shape, f, *archive.Debug)) {
            archive.Surfaces.push_back(FaceGeometry(shape, f));
        }
    }
    STEREON_TRY(const auto coedges, AllEntities<EntityKind::Coedge>(shape));
    for (const CoedgeId c : coedges) {
        if (const std::optional<Curve2d> pcurve = PCurveOf(shape, c)) {
            archive.Curves2d.push_back(*pcurve);
        }
    }
    return {};
}

/// The entity a finding names, found by its Describe() text: the
/// validator names entities as Describe does. Nothing for a name no live
/// vertex, edge, coedge or face has.
struct Named {
    std::optional<VertexId> Vertex;
    std::optional<EdgeId> Edge;
    std::optional<FaceId> Face;
};

template <EntityKind K, typename Id>
[[nodiscard]] std::optional<Id> Find(const Shape& shape, const std::string& name)
{
    const auto all = AllEntities<K>(shape);
    if (!all) {
        return std::nullopt;
    }
    for (const Id id : *all) {
        if (Describe(shape, id) == name) {
            return id;
        }
    }
    return std::nullopt;
}

[[nodiscard]] Named FindNamed(const Shape& shape, const std::string& name)
{
    Named named;
    named.Vertex = Find<EntityKind::Vertex, VertexId>(shape, name);
    named.Edge = Find<EntityKind::Edge, EdgeId>(shape, name);
    named.Face = Find<EntityKind::Face, FaceId>(shape, name);
    if (!named.Edge) {
        // A coedge stands for its edge.
        if (const auto coedge = Find<EntityKind::Coedge, CoedgeId>(shape, name)) {
            named.Edge = EdgeOf(shape, *coedge);
        }
    }
    return named;
}

/// "Rule (Edge#8, Shell#0): message": the entities by their short names,
/// the part of their Describe() text before the bracket.
[[nodiscard]] std::string FindingText(const ValidityFinding& finding)
{
    std::string text(RuleName(finding.Rule));
    std::string names;
    for (const std::string& entity : finding.Entities) {
        names += (names.empty() ? "" : ", ") + entity.substr(0, entity.find(" ["));
    }
    if (!names.empty()) {
        text += " (" + names + ")";
    }
    return text + ": " + finding.Message;
}

/// A finding drawn on the shape: every entity it names in FindingColour,
/// a vertex as a point, an edge as its polyline, a face as its lines, and
/// the rule and message as a label at the first of them.
void AddFinding(const Shape& shape, const ValidityFinding& finding, DebugChunk& debug)
{
    std::optional<Point3> anchor;
    for (const std::string& name : finding.Entities) {
        const Named named = FindNamed(shape, name);
        if (named.Vertex) {
            DebugPoint point;
            point.Position = PositionOf(shape, *named.Vertex);
            point.Colour = FindingColour();
            debug.Points.push_back(point);
            anchor = anchor.value_or(point.Position);
        }
        if (named.Edge) {
            std::vector<Point3> points = EdgePoints(shape, *named.Edge);
            if (!points.empty()) {
                anchor = anchor.value_or(points[points.size() / 2]);
                DebugPolyline line;
                line.Points = std::move(points);
                line.Colour = FindingColour();
                debug.Polylines.push_back(std::move(line));
            }
        }
        if (named.Face) {
            const std::size_t before = debug.Polylines.size();
            if (AddFaceLines(shape, *named.Face, debug, FindingColour()) &&
                debug.Polylines.size() > before && !debug.Polylines[before].Points.empty()) {
                const std::vector<Point3>& first = debug.Polylines[before].Points;
                anchor = anchor.value_or(first[first.size() / 2]);
            }
        }
    }
    DebugLabel label;
    label.Position = anchor.value_or(Point3{});
    label.Text = FindingText(finding);
    label.Colour = FindingColour();
    debug.Labels.push_back(std::move(label));
}

/// A double in the shortest form that reads back as itself.
[[nodiscard]] std::string Short(const double x)
{
    char buffer[32];
    const auto [end, error] = std::to_chars(buffer, buffer + sizeof buffer, x);
    return error == std::errc{} ? std::string(buffer, end) : std::string("?");
}

[[nodiscard]] std::string_view ContactName(const ContactType contact) noexcept
{
    switch (contact) {
    case ContactType::Crossing:
        return "Crossing";
    case ContactType::Tangential:
        return "Tangential";
    case ContactType::Overlap:
        return "Overlap";
    }
    return "Unknown";
}

[[nodiscard]] std::string_view TransitionName(const SideTransition transition) noexcept
{
    switch (transition) {
    case SideTransition::Unresolved:
        return "Unresolved";
    case SideTransition::SameSide:
        return "SameSide";
    case SideTransition::OppositeSide:
        return "OppositeSide";
    }
    return "Unknown";
}

/// The surface trimmed to the box when the box is bounded and the trim is
/// accepted; the surface itself otherwise, which the viewer draws around
/// the scene.
[[nodiscard]] Surface BoundedTo(const Surface& surface, const Box2& box)
{
    if (!box.X.IsBounded() || !box.Y.IsBounded()) {
        return surface;
    }
    Result<TrimmedSurface> trimmed = MakeTrimmedSurface(surface, box);
    return trimmed ? Surface(*std::move(trimmed)) : surface;
}

/// A branch: its 3D curve over its range as a polyline in its contact's
/// colour, its pcurves on A and B, and a label at its middle.
[[nodiscard]] Result<void> AddBranch(const IntersectionBranch& branch, const std::size_t k,
                                     GeometryArchive& archive)
{
    if (!branch.Range.IsBounded()) {
        return Fail(ErrorCode::NotSupported,
                    "IntersectionArchive: branch " + std::to_string(k) +
                        " has an unbounded range; bound the requested boxes to dump it");
    }
    STEREON_TRY(const Polyline line, SampleCurve(branch.Curve3d, DisplayTolerance{}, branch.Range));
    DebugPolyline polyline;
    polyline.Points = line.Points;
    polyline.Colour = ContactColour(branch.Contact);
    archive.Debug->Polylines.push_back(std::move(polyline));
    archive.Curves2d.push_back(branch.OnA);
    archive.Curves2d.push_back(branch.OnB);
    DebugLabel label;
    label.Position = line.Points[line.Points.size() / 2];
    label.Text = "branch " + std::to_string(k) + ": " + std::string(ContactName(branch.Contact)) +
                 ", " + std::string(TransitionName(branch.Transition)) + ", bound " +
                 Short(branch.Bound);
    label.Colour = ContactColour(branch.Contact);
    archive.Debug->Labels.push_back(std::move(label));
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

DumpResult DumpAt(const std::string_view name, const Shape& shape, const ValidityReport& report,
                  const SourceLocation where)
{
    if (IsDumpName(name) && !DumpDirectory()) {
        return std::optional<std::filesystem::path>{};
    }
    STEREON_TRY(GeometryArchive archive, ValidityArchive(shape, report));
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

DebugColour FaceColour() noexcept
{
    return Rgb(150, 150, 160);
}

DebugColour FindingColour() noexcept
{
    return Rgb(255, 40, 40);
}

DebugColour ContactColour(const ContactType contact) noexcept
{
    switch (contact) {
    case ContactType::Crossing:
        return Rgb(230, 50, 50);
    case ContactType::Tangential:
        return Rgb(240, 150, 30);
    case ContactType::Overlap:
        return Rgb(60, 190, 90);
    }
    return Rgb(0, 0, 0);
}

Result<GeometryArchive> IntersectionArchive(const Surface& a, const Box2& boxA, const Surface& b,
                                            const Box2& boxB,
                                            const SurfaceSurfaceIntersection& result)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Debug.emplace();
    archive.Surfaces.push_back(BoundedTo(a, boxA));
    archive.Surfaces.push_back(BoundedTo(b, boxB));
    for (std::size_t k = 0; k < result.Branches.size(); ++k) {
        STEREON_TRY_VOID(AddBranch(result.Branches[k], k, archive));
    }
    for (std::size_t k = 0; k < result.Points.size(); ++k) {
        const IntersectionPoint& p = result.Points[k];
        DebugPoint point;
        point.Position = p.Position;
        point.Colour = ContactColour(ContactType::Tangential);
        archive.Debug->Points.push_back(point);
        DebugLabel label;
        label.Position = p.Position;
        label.Text = "point " + std::to_string(k) + ": bound " + Short(p.Bound);
        label.Colour = point.Colour;
        archive.Debug->Labels.push_back(std::move(label));
    }
    for (std::size_t k = 0; k < result.Overlaps.size(); ++k) {
        const OverlapRegion& overlap = result.Overlaps[k];
        DebugLabel label;
        label.Position = Eval(a, Point2{overlap.OnA.X.GetLower(), overlap.OnA.Y.GetLower()});
        label.Text = "overlap " + std::to_string(k) +
                     (overlap.SameOrientation ? ": same orientation" : ": opposite orientation");
        label.Colour = ContactColour(ContactType::Overlap);
        archive.Debug->Labels.push_back(std::move(label));
    }
    return archive;
}

Result<GeometryArchive> IntersectionArchive(const Curve& curve, const Interval& domain,
                                            const Surface& surface, const Box2& box,
                                            const CurveSurfaceIntersection& result)
{
    GeometryArchive archive = ArchiveOf(Context{});
    archive.Debug.emplace();
    // The curve over its domain when bounded; the curve itself, drawn
    // around the scene, otherwise or when the trim is refused.
    archive.Curves.push_back(curve);
    if (domain.IsBounded()) {
        if (Result<TrimmedCurve> trimmed = MakeTrimmedCurve(curve, domain)) {
            archive.Curves.back() = Curve(*std::move(trimmed));
        }
    }
    archive.Surfaces.push_back(BoundedTo(surface, box));
    for (const CurveSurfacePoint& p : result.Points) {
        DebugPoint point;
        point.Position = p.Position;
        point.Colour = ContactColour(p.Contact);
        archive.Debug->Points.push_back(point);
        DebugLabel label;
        label.Position = p.Position;
        label.Text = "t = " + Short(p.Parameter) + ": " + std::string(ContactName(p.Contact)) +
                     (p.Decided ? "" : ", undecided");
        label.Colour = point.Colour;
        archive.Debug->Labels.push_back(std::move(label));
    }
    for (const CurveSurfaceOverlap& overlap : result.Overlaps) {
        STEREON_TRY(const Polyline line, SampleCurve(curve, DisplayTolerance{}, overlap.Range));
        DebugPolyline polyline;
        polyline.Points = line.Points;
        polyline.Colour = ContactColour(ContactType::Overlap);
        archive.Debug->Polylines.push_back(std::move(polyline));
        archive.Curves2d.push_back(overlap.OnSurface);
    }
    return archive;
}

Result<GeometryArchive> ValidityArchive(const Shape& shape, const ValidityReport& report)
{
    STEREON_TRY(GeometryArchive archive, ShapeArchive(shape));
    DebugChunk overlay;
    for (const ValidityFinding& finding : report.GetFindings()) {
        // Every Solid-level report has it: the check it names is not made.
        if (finding.Rule != TopologyRule::SelfIntersectionUnchecked) {
            AddFinding(shape, finding, overlay);
        }
    }
    if (overlay.Labels.empty()) {
        return archive;
    }
    // The findings' view: their labels alone, and their primitives first,
    // since a viewer draws the first of two lines at one depth on top.
    DebugChunk& debug = *archive.Debug;
    debug.Labels = std::move(overlay.Labels);
    overlay.Points.insert(overlay.Points.end(), debug.Points.begin(), debug.Points.end());
    debug.Points = std::move(overlay.Points);
    overlay.Polylines.insert(overlay.Polylines.end(), debug.Polylines.begin(),
                             debug.Polylines.end());
    debug.Polylines = std::move(overlay.Polylines);
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
