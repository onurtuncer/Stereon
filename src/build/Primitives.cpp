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

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Build;

import :Primitives;

namespace Stereon {

// Instantiated at namespace scope before any function uses them
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));
static_assert(sizeof(Vec2) == 2 * sizeof(double));

namespace {

/// A length larger than the resolution, as every dimension must be.
[[nodiscard]] bool IsDimension(const double length, const Context& context)
{
    return std::isfinite(length) && length > context.GetResolution();
}

[[nodiscard]] Result<void> CheckFrame(const Point3& origin, const Vec3& xAxis, const Vec3& yAxis,
                                      const std::string_view operation)
{
    if (!MakePlane(origin, xAxis, yAxis)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(operation) +
                        ": the origin must be finite and the axes orthonormal within "
                        "FrameAllowance");
    }
    return {};
}

/// The plane of a convex face through its first three corners, normal along
/// their counter-clockwise turn, so outward for corners listed
/// counter-clockwise from outside.
[[nodiscard]] Plane PlaneThrough(const Point3& a, const Point3& b, const Point3& c)
{
    const Vec3 x = Normalized(b - a);
    const Vec3 n = Normalized(Cross(b - a, c - a));
    Plane plane;
    plane.Origin = a;
    plane.XAxis = x;
    plane.YAxis = Cross(n, x);
    return plane;
}

[[nodiscard]] Point2 UvOn(const Plane& plane, const Point3& p)
{
    const Vec3 d = p - plane.Origin;
    return Point2{Dot(d, plane.XAxis), Dot(d, plane.YAxis)};
}

[[nodiscard]] EdgeSpec ExactEdge(const BuilderVertexId start, const BuilderVertexId end,
                                 const CurveRef curve, const Interval range, const double tolerance)
{
    // Member by member: a designated initializer of another module's
    // aggregate is not used across module units (gcc16-feature-probe.md).
    EdgeSpec spec;
    spec.Kind = EdgeCurveKind::Exact;
    spec.Start = start;
    spec.End = end;
    spec.Curve = curve;
    spec.Range = range;
    spec.Tolerance = tolerance;
    return spec;
}

/// The corners of the box, 0-3 counter-clockwise at the bottom and 4-7
/// above them, and its faces as corner lists counter-clockwise from
/// outside: bottom, top, then the four sides.
constexpr std::array<std::array<std::size_t, 4>, 6> BoxFaces{{
    {0, 3, 2, 1},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {2, 3, 7, 6},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
}};

/// The box's builder state: its vertices and the edges made so far, by
/// their corner pair.
class BoxBuilder {
public:
    BoxBuilder(std::array<Point3, 8> corners, const double tolerance)
        : m_Corners(corners)
        , m_Tolerance(tolerance)
    {}

    [[nodiscard]] Result<ShapeBuilder> Build() &&
    {
        for (std::size_t i = 0; i < m_Corners.size(); ++i) {
            STEREON_TRY(m_Vertices[i], m_B.AddVertex(m_Corners[i], m_Tolerance));
        }
        STEREON_TRY(const BuilderSolidId solid, m_B.AddSolid());
        STEREON_TRY(const BuilderShellId shell, m_B.AddShell(solid));
        for (const std::array<std::size_t, 4>& corners : BoxFaces) {
            STEREON_TRY_VOID(AddFace(shell, corners));
        }
        return std::move(m_B);
    }

private:
    [[nodiscard]] Result<void> AddFace(const BuilderShellId shell,
                                       const std::array<std::size_t, 4>& corners)
    {
        const Plane plane =
            PlaneThrough(m_Corners[corners[0]], m_Corners[corners[1]], m_Corners[corners[2]]);
        STEREON_TRY(const SurfaceRef surface, m_B.AddSurface(Surface(plane)));
        STEREON_TRY(const BuilderFaceId face, m_B.AddFace(surface, false));
        STEREON_TRY_VOID(m_B.AddFaceToShell(shell, face));
        STEREON_TRY(const BuilderLoopId loop, m_B.AddLoop(face));
        for (std::size_t i = 0; i < corners.size(); ++i) {
            const std::size_t a = corners[i];
            const std::size_t b = corners[(i + 1) % corners.size()];
            const std::size_t lo = a < b ? a : b;
            const std::size_t hi = a < b ? b : a;
            STEREON_TRY(const BuilderEdgeId edge, EdgeOf(lo, hi));
            const Point2 from = UvOn(plane, m_Corners[lo]);
            const Point2 to = UvOn(plane, m_Corners[hi]);
            STEREON_TRY(const PCurveRef pcurve, m_B.AddPCurve(Curve2d(Line2d{from, to - from})));
            STEREON_TRY_VOID(m_B.AppendCoedge(
                loop, edge, a == lo ? EdgeUseSense::Same : EdgeUseSense::Opposite, pcurve));
        }
        return {};
    }

    /// The line edge from corner lo to corner hi on [0, 1], made once.
    [[nodiscard]] Result<BuilderEdgeId> EdgeOf(const std::size_t lo, const std::size_t hi)
    {
        for (const auto& [key, edge] : m_Edges) {
            if (key.first == lo && key.second == hi) {
                return edge;
            }
        }
        const Point3 start = m_Corners[lo];
        const Point3 end = m_Corners[hi];
        STEREON_TRY(const CurveRef line, m_B.AddCurve(Curve(Line{start, end - start})));
        STEREON_TRY(const BuilderEdgeId edge,
                    m_B.AddEdge(ExactEdge(m_Vertices[lo], m_Vertices[hi], line, Interval{0.0, 1.0},
                                          m_Tolerance)));
        m_Edges.emplace_back(std::pair{lo, hi}, edge);
        return edge;
    }

    std::array<Point3, 8> m_Corners;
    double m_Tolerance;
    ShapeBuilder m_B;
    std::array<BuilderVertexId, 8> m_Vertices{};
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, BuilderEdgeId>> m_Edges;
};

/// The cylinder's three faces and their coedges, in the order the shape
/// lists them: the lateral face, then the bottom and top disks.
struct CylinderParts {
    BuilderEdgeId Bottom;
    BuilderEdgeId Top;
    BuilderEdgeId Seam;
};

[[nodiscard]] Result<void> AddLateralFace(ShapeBuilder& b, const BuilderShellId shell,
                                          const Cylinder& cylinder, const CylinderParts& parts,
                                          const double height)
{
    constexpr double tau = 2.0 * std::numbers::pi;
    STEREON_TRY(const SurfaceRef surface, b.AddSurface(Surface(cylinder)));
    STEREON_TRY(const BuilderFaceId face, b.AddFace(surface, false));
    STEREON_TRY_VOID(b.AddFaceToShell(shell, face));
    STEREON_TRY(const BuilderLoopId loop, b.AddLoop(face));
    // Around the lateral face: along the bottom circle at v = 0, up the
    // seam at u = 2 pi, back along the top circle at v = height, down the
    // seam at u = 0.
    const std::array<std::pair<BuilderEdgeId, std::pair<EdgeUseSense, Curve2d>>, 4> uses{{
        {parts.Bottom, {EdgeUseSense::Same, Curve2d(Line2d{Point2{0, 0}, Vec2{1, 0}})}},
        {parts.Seam, {EdgeUseSense::Same, Curve2d(Line2d{Point2{tau, 0}, Vec2{0, 1}})}},
        {parts.Top, {EdgeUseSense::Opposite, Curve2d(Line2d{Point2{0, height}, Vec2{1, 0}})}},
        {parts.Seam, {EdgeUseSense::Opposite, Curve2d(Line2d{Point2{0, 0}, Vec2{0, 1}})}},
    }};
    for (const auto& [edge, use] : uses) {
        STEREON_TRY(const PCurveRef pcurve, b.AddPCurve(use.second));
        STEREON_TRY_VOID(b.AppendCoedge(loop, edge, use.first, pcurve));
    }
    return {};
}

/// A disk at the end of the cylinder: the plane through `centre` with axes
/// `x` and `y` (outward normal x cross y), and the edge's circle in it.
[[nodiscard]] Result<void> AddDisk(ShapeBuilder& b, const BuilderShellId shell,
                                   const Point3& centre, const Vec3& x, const Vec3& y,
                                   const BuilderEdgeId edge, const Curve2d& circle,
                                   const EdgeUseSense sense)
{
    Plane plane;
    plane.Origin = centre;
    plane.XAxis = x;
    plane.YAxis = y;
    STEREON_TRY(const SurfaceRef surface, b.AddSurface(Surface(plane)));
    STEREON_TRY(const BuilderFaceId face, b.AddFace(surface, false));
    STEREON_TRY_VOID(b.AddFaceToShell(shell, face));
    STEREON_TRY(const BuilderLoopId loop, b.AddLoop(face));
    STEREON_TRY(const PCurveRef pcurve, b.AddPCurve(circle));
    STEREON_TRY_VOID(b.AppendCoedge(loop, edge, sense, pcurve));
    return {};
}

} // namespace

Result<OperationResult> MakeBoxSolid(const Point3& origin, const Vec3& xAxis, const Vec3& yAxis,
                                     const Vec3& size, const Context& context)
{
    STEREON_TRY_VOID(CheckFrame(origin, xAxis, yAxis, "MakeBoxSolid"));
    if (!IsDimension(size.X, context) || !IsDimension(size.Y, context) ||
        !IsDimension(size.Z, context)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeBoxSolid: every side must be finite and longer than the resolution");
    }
    const Vec3 x = size.X * xAxis;
    const Vec3 y = size.Y * yAxis;
    const Vec3 z = size.Z * Cross(xAxis, yAxis);
    const std::array<Point3, 8> corners{
        origin,     origin + x,     origin + x + y,     origin + y,
        origin + z, origin + x + z, origin + x + y + z, origin + y + z,
    };
    STEREON_TRY(ShapeBuilder builder, BoxBuilder(corners, context.GetResolution()).Build());
    return Commit(std::move(builder), context);
}

Result<OperationResult> MakeBoxSolid(const Box3& box, const Context& context)
{
    const Point3 lower{box.X.GetLower(), box.Y.GetLower(), box.Z.GetLower()};
    const Vec3 size{box.X.GetUpper() - box.X.GetLower(), box.Y.GetUpper() - box.Y.GetLower(),
                    box.Z.GetUpper() - box.Z.GetLower()};
    return MakeBoxSolid(lower, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, size, context);
}

Result<OperationResult> MakeCylinderSolid(const Point3& baseCentre, const Vec3& xAxis,
                                          const Vec3& yAxis, const double radius,
                                          const double height, const Context& context)
{
    STEREON_TRY_VOID(CheckFrame(baseCentre, xAxis, yAxis, "MakeCylinderSolid"));
    if (!IsDimension(radius, context) || !IsDimension(height, context)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeCylinderSolid: the radius and the height must be finite and larger than "
                    "the resolution");
    }
    constexpr double tau = 2.0 * std::numbers::pi;
    const double tolerance = context.GetResolution();
    const Vec3 z = Cross(xAxis, yAxis);
    const Point3 topCentre = baseCentre + height * z;
    ShapeBuilder b;
    STEREON_TRY(const BuilderVertexId v0, b.AddVertex(baseCentre + radius * xAxis, tolerance));
    STEREON_TRY(const BuilderVertexId v1, b.AddVertex(topCentre + radius * xAxis, tolerance));
    STEREON_TRY(const CurveRef bottomCircle,
                b.AddCurve(Curve(Circle{baseCentre, radius, xAxis, yAxis})));
    STEREON_TRY(const CurveRef topCircle,
                b.AddCurve(Curve(Circle{topCentre, radius, xAxis, yAxis})));
    CylinderParts parts;
    STEREON_TRY(parts.Bottom,
                b.AddEdge(ExactEdge(v0, v0, bottomCircle, Interval{0.0, tau}, tolerance)));
    STEREON_TRY(parts.Top, b.AddEdge(ExactEdge(v1, v1, topCircle, Interval{0.0, tau}, tolerance)));
    EdgeSpec seam;
    seam.Kind = EdgeCurveKind::Boundary;
    seam.Start = v0;
    seam.End = v1;
    seam.Range = Interval{0.0, height};
    seam.Tolerance = tolerance;
    STEREON_TRY(parts.Seam, b.AddEdge(seam));

    STEREON_TRY(const BuilderSolidId solid, b.AddSolid());
    STEREON_TRY(const BuilderShellId shell, b.AddShell(solid));
    Cylinder cylinder;
    cylinder.Origin = baseCentre;
    cylinder.Radius = radius;
    cylinder.XAxis = xAxis;
    cylinder.YAxis = yAxis;
    STEREON_TRY_VOID(AddLateralFace(b, shell, cylinder, parts, height));
    // The bottom disk faces down: axes x and -y, in which the circle
    // (r cos t, r sin t) reads (r cos t, -r sin t).
    STEREON_TRY_VOID(AddDisk(b, shell, baseCentre, xAxis, -yAxis, parts.Bottom,
                             Curve2d(Circle2d{Point2{}, radius, Vec2{1, 0}, Vec2{0, -1}}),
                             EdgeUseSense::Opposite));
    STEREON_TRY_VOID(AddDisk(b, shell, topCentre, xAxis, yAxis, parts.Top,
                             Curve2d(Circle2d{Point2{}, radius, Vec2{1, 0}, Vec2{0, 1}}),
                             EdgeUseSense::Same));
    return Commit(std::move(b), context);
}

} // namespace Stereon
