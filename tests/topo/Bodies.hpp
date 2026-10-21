// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Hand-built test bodies for the topology tests, built through ShapeBuilder's
// controlled operations only. Include after `import Stereon.Topo;`.
//
// - MakePolyhedron: planar faces given as outward counter-clockwise vertex
//   lists, Exact line edges with exact pcurves, one shell bounding one solid.
//   MakeBox and MakePrism use it.
// - MakeCylinder: a closed cylinder, lateral face with a Boundary seam edge
//   used twice, two Exact circles, two planar disks.
// - MakeSquareSheet, MakeHoleSheet: open sheets of Boundary edges.
// - MakeFan: k planar faces around one shared Exact edge, non-manifold for
//   k >= 3.
// - MakeWire: Spatial edges with no coedges.
// - MakeDegenerateFace: a loop that runs through a Degenerate edge.

#pragma once

// The standard and Catch2 headers this file needs come from BodiesStd.hpp,
// which every test includes before `import Stereon.Topo;`: GCC 16.2 rejects
// a standard header first included after an import of a module that uses it.

namespace TopoTest {

using Stereon::BuilderCoedgeId;
using Stereon::BuilderEdgeId;
using Stereon::BuilderFaceId;
using Stereon::BuilderLoopId;
using Stereon::BuilderShellId;
using Stereon::BuilderSolidId;
using Stereon::BuilderVertexId;
using Stereon::EdgeCurveKind;
using Stereon::EdgeUseSense;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;

/// The value of a successful Result; fails the test otherwise.
template <typename T>
T Must(Stereon::Result<T> result)
{
    if (!result) {
        FAIL("unexpected error: " << std::string(result.error().GetMessage()));
    }
    return std::move(*result);
}

inline void Must(const Stereon::Result<void>& result)
{
    if (!result) {
        FAIL("unexpected error: " << std::string(result.error().GetMessage()));
    }
}

/// The tolerance of every entity the bodies here make: the default Context's
/// resolution, since native construction is resolution-tight and no local
/// tolerance is below the resolution (ADR-0001), which Commit() checks.
inline constexpr double Tight = 1e-5;

/// A builder and the handles of everything it made, in creation order.
struct Body {
    Stereon::ShapeBuilder Builder;
    std::vector<BuilderVertexId> Vertices;
    std::vector<BuilderEdgeId> Edges;
    std::vector<BuilderFaceId> Faces;
    std::vector<BuilderLoopId> Loops;
    std::vector<BuilderCoedgeId> Coedges;
    std::vector<BuilderShellId> Shells;
    std::vector<BuilderSolidId> Solids;

    BuilderVertexId AddVertex(const Point3& p, const double tolerance = Tight)
    {
        return Vertices.emplace_back(Must(Builder.AddVertex(p, tolerance)));
    }
    BuilderEdgeId AddEdge(const Stereon::EdgeSpec& spec)
    {
        return Edges.emplace_back(Must(Builder.AddEdge(spec)));
    }
    BuilderFaceId AddFace(const Stereon::Surface& surface, const bool reversed = false)
    {
        return Faces.emplace_back(
            Must(Builder.AddFace(Must(Builder.AddSurface(surface)), reversed)));
    }
    BuilderLoopId AddLoop(const BuilderFaceId face)
    {
        return Loops.emplace_back(Must(Builder.AddLoop(face)));
    }
    BuilderCoedgeId AddCoedge(const BuilderLoopId loop, const BuilderEdgeId edge,
                              const EdgeUseSense sense, const Stereon::Curve2d& pcurve)
    {
        return Coedges.emplace_back(
            Must(Builder.AppendCoedge(loop, edge, sense, Must(Builder.AddPCurve(pcurve)))));
    }

    Stereon::Shape Freeze() && { return Must(Stereon::Freeze(std::move(Builder))); }
};

/// A pcurve for bodies whose geometry the test does not look at.
inline Stereon::Curve2d AnyPCurve()
{
    return Stereon::Curve2d(Stereon::Line2d{Point2{0.0, 0.0}, Vec2{1.0, 0.0}});
}

/// The plane through a convex face's first three corners, normal along
/// their counter-clockwise turn.
inline Stereon::Plane PlaneThrough(const Point3& a, const Point3& b, const Point3& c)
{
    const Vec3 x = Normalized(b - a);
    const Vec3 n = Normalized(Cross(b - a, c - a));
    return {a, x, Cross(n, x)};
}

inline Point2 UvOn(const Stereon::Plane& plane, const Point3& p)
{
    const Vec3 d = p - plane.Origin;
    return {Dot(d, plane.XAxis), Dot(d, plane.YAxis)};
}

/// A closed polyhedron: `faces` lists vertex indices counter-clockwise seen
/// from outside. Each pair of vertices joined by a face side becomes one
/// Exact line edge from the lower index to the higher, on t in [0, 1]; each
/// face side a coedge with the face's exact pcurve of that edge.
inline Body MakePolyhedron(std::span<const Point3> points,
                           std::span<const std::vector<std::size_t>> faces)
{
    Body body;
    for (const Point3& p : points) {
        body.AddVertex(p);
    }
    std::map<std::pair<std::size_t, std::size_t>, BuilderEdgeId> edges;
    const auto edgeOf = [&](const std::size_t a, const std::size_t b) {
        const auto key = std::minmax(a, b);
        if (const auto it = edges.find(key); it != edges.end()) {
            return it->second;
        }
        const Point3 start = points[key.first];
        const Point3 end = points[key.second];
        const BuilderEdgeId e = body.AddEdge({
            .Kind = EdgeCurveKind::Exact,
            .Start = body.Vertices[key.first],
            .End = body.Vertices[key.second],
            .Curve = Must(body.Builder.AddCurve(Stereon::Curve(Stereon::Line{start, end - start}))),
            .Range = Interval(0.0, 1.0),
            .Tolerance = Tight,
        });
        edges.emplace(key, e);
        return e;
    };

    const BuilderSolidId solid = body.Solids.emplace_back(Must(body.Builder.AddSolid()));
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell(solid)));
    for (const std::vector<std::size_t>& corners : faces) {
        const Stereon::Plane plane =
            PlaneThrough(points[corners[0]], points[corners[1]], points[corners[2]]);
        const BuilderFaceId face = body.AddFace(Stereon::Surface(plane));
        Must(body.Builder.AddFaceToShell(shell, face));
        const BuilderLoopId loop = body.AddLoop(face);
        for (std::size_t i = 0; i < corners.size(); ++i) {
            const std::size_t a = corners[i];
            const std::size_t b = corners[(i + 1) % corners.size()];
            const auto [lo, hi] = std::minmax(a, b);
            const Point2 from = UvOn(plane, points[lo]);
            const Point2 to = UvOn(plane, points[hi]);
            body.AddCoedge(loop, edgeOf(a, b),
                           a == lo ? EdgeUseSense::Same : EdgeUseSense::Opposite,
                           Stereon::Curve2d(Stereon::Line2d{from, to - from}));
        }
    }
    return body;
}

/// The box [0, a] x [0, b] x [0, c]: vertices 0-3 on z = 0, 4-7 above them.
inline Body MakeBox(const double a = 2.0, const double b = 3.0, const double c = 5.0)
{
    const std::vector<Point3> points{{0, 0, 0}, {a, 0, 0}, {a, b, 0}, {0, b, 0},
                                     {0, 0, c}, {a, 0, c}, {a, b, c}, {0, b, c}};
    const std::vector<std::vector<std::size_t>> faces{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                      {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    return MakePolyhedron(points, faces);
}

/// A prism over a regular n-gon of unit radius, height h: bottom vertices
/// 0..n-1, top n..2n-1.
inline Body MakePrism(const std::size_t n, const double h = 1.5)
{
    std::vector<Point3> points;
    for (std::size_t level = 0; level < 2; ++level) {
        for (std::size_t i = 0; i < n; ++i) {
            const double angle =
                2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
            points.push_back({std::cos(angle), std::sin(angle), static_cast<double>(level) * h});
        }
    }
    std::vector<std::vector<std::size_t>> faces;
    std::vector<std::size_t> bottom;
    std::vector<std::size_t> top;
    for (std::size_t i = 0; i < n; ++i) {
        bottom.push_back(n - 1 - i);
        top.push_back(n + i);
        faces.push_back({i, (i + 1) % n, n + (i + 1) % n, n + i});
    }
    faces.push_back(bottom);
    faces.push_back(top);
    return MakePolyhedron(points, faces);
}

/// The closed cylinder of radius r and height h about the z axis.
/// Vertices: 0 at (r, 0, 0), 1 at (r, 0, h). Edges: 0 the bottom circle,
/// 1 the top circle (Exact, t in [0, 2 pi]), 2 the seam (Boundary, the line
/// up from vertex 0, t in [0, h]). Faces: 0 the lateral face, 1 the bottom
/// disk, 2 the top disk. The lateral loop runs (u, v) = (0, 0) -> (2 pi, 0)
/// -> (2 pi, h) -> (0, h): bottom circle along, seam at u = 2 pi up, top
/// circle against, seam at u = 0 down.
inline Body MakeCylinder(const double r = 1.5, const double h = 4.0)
{
    constexpr double tau = 2.0 * std::numbers::pi;
    Body body;
    const BuilderVertexId v0 = body.AddVertex({r, 0, 0});
    const BuilderVertexId v1 = body.AddVertex({r, 0, h});
    const auto circle = [&](const double z) {
        return Must(body.Builder.AddCurve(
            Stereon::Curve(Stereon::Circle{{0, 0, z}, r, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
    };
    const BuilderEdgeId bottom = body.AddEdge({.Kind = EdgeCurveKind::Exact,
                                               .Start = v0,
                                               .End = v0,
                                               .Curve = circle(0.0),
                                               .Range = Interval(0.0, tau),
                                               .Tolerance = Tight});
    const BuilderEdgeId top = body.AddEdge({.Kind = EdgeCurveKind::Exact,
                                            .Start = v1,
                                            .End = v1,
                                            .Curve = circle(h),
                                            .Range = Interval(0.0, tau),
                                            .Tolerance = Tight});
    const BuilderEdgeId seam = body.AddEdge({.Kind = EdgeCurveKind::Boundary,
                                             .Start = v0,
                                             .End = v1,
                                             .Range = Interval(0.0, h),
                                             .Tolerance = Tight});

    const BuilderSolidId solid = body.Solids.emplace_back(Must(body.Builder.AddSolid()));
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell(solid)));

    const BuilderFaceId lateral = body.AddFace(
        Stereon::Surface(Stereon::Cylinder{{0, 0, 0}, r, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    const BuilderLoopId around = body.AddLoop(lateral);
    using Stereon::Curve2d;
    using Stereon::Line2d;
    body.AddCoedge(around, bottom, EdgeUseSense::Same, Curve2d(Line2d{{0, 0}, Vec2{1, 0}}));
    body.AddCoedge(around, seam, EdgeUseSense::Same, Curve2d(Line2d{{tau, 0}, Vec2{0, 1}}));
    body.AddCoedge(around, top, EdgeUseSense::Opposite, Curve2d(Line2d{{0, h}, Vec2{1, 0}}));
    body.AddCoedge(around, seam, EdgeUseSense::Opposite, Curve2d(Line2d{{0, 0}, Vec2{0, 1}}));

    // The bottom disk faces down: plane axes x and -y. In them the circle
    // (r cos t, r sin t) is (r cos t, -r sin t).
    const BuilderFaceId down =
        body.AddFace(Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, -1, 0}}));
    body.AddCoedge(body.AddLoop(down), bottom, EdgeUseSense::Opposite,
                   Curve2d(Stereon::Circle2d{{0, 0}, r, Vec2{1, 0}, Vec2{0, -1}}));
    const BuilderFaceId up =
        body.AddFace(Stereon::Surface(Stereon::Plane{{0, 0, h}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    body.AddCoedge(body.AddLoop(up), top, EdgeUseSense::Same,
                   Curve2d(Stereon::Circle2d{{0, 0}, r, Vec2{1, 0}, Vec2{0, 1}}));
    for (const BuilderFaceId f : {lateral, down, up}) {
        Must(body.Builder.AddFaceToShell(shell, f));
    }
    return body;
}

/// Adds a closed polygon of Boundary edges as a new loop of `face`, which
/// lies in the plane z = 0 with the identity frame.
inline void AddBoundaryPolygon(Body& body, const BuilderFaceId face,
                               std::span<const Point3> corners)
{
    std::vector<BuilderVertexId> vertices;
    for (const Point3& p : corners) {
        vertices.push_back(body.AddVertex(p));
    }
    const BuilderLoopId loop = body.AddLoop(face);
    for (std::size_t i = 0; i < corners.size(); ++i) {
        const std::size_t j = (i + 1) % corners.size();
        const BuilderEdgeId e = body.AddEdge({.Kind = EdgeCurveKind::Boundary,
                                              .Start = vertices[i],
                                              .End = vertices[j],
                                              .Range = Interval(0.0, 1.0),
                                              .Tolerance = Tight});
        const Point2 from{corners[i].X, corners[i].Y};
        const Point2 to{corners[j].X, corners[j].Y};
        body.AddCoedge(loop, e, EdgeUseSense::Same,
                       Stereon::Curve2d(Stereon::Line2d{from, to - from}));
    }
}

inline Stereon::Surface XyPlane()
{
    return Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}});
}

/// The unit square in z = 0 as an open sheet: one face, four Boundary
/// edges, a shell of no solid.
inline Body MakeSquareSheet()
{
    Body body;
    const BuilderFaceId face = body.AddFace(XyPlane());
    const std::vector<Point3> corners{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    AddBoundaryPolygon(body, face, corners);
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    Must(body.Builder.AddFaceToShell(shell, face));
    return body;
}

/// A square sheet [0, 3]^2 with a square hole [1, 2]^2: loop 0 the outer
/// boundary counter-clockwise, loop 1 the hole clockwise.
inline Body MakeHoleSheet()
{
    Body body;
    const BuilderFaceId face = body.AddFace(XyPlane());
    const std::vector<Point3> outer{{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 3, 0}};
    const std::vector<Point3> hole{{1, 1, 0}, {1, 2, 0}, {2, 2, 0}, {2, 1, 0}};
    AddBoundaryPolygon(body, face, outer);
    AddBoundaryPolygon(body, face, hole);
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    Must(body.Builder.AddFaceToShell(shell, face));
    return body;
}

/// k triangles around the shared Exact edge 0 from (0, 0, 0) to (0, 0, 1),
/// each in its own plane through the z axis, all in one open shell. The
/// shared edge has k uses: a non-manifold edge for k >= 3.
inline Body MakeFan(const std::size_t k)
{
    Body body;
    const BuilderVertexId v0 = body.AddVertex({0, 0, 0});
    const BuilderVertexId v1 = body.AddVertex({0, 0, 1});
    const BuilderEdgeId spine = body.AddEdge(
        {.Kind = EdgeCurveKind::Exact,
         .Start = v0,
         .End = v1,
         .Curve =
             Must(body.Builder.AddCurve(Stereon::Curve(Stereon::Line{{0, 0, 0}, Vec3{0, 0, 1}}))),
         .Range = Interval(0.0, 1.0),
         .Tolerance = Tight});
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    for (std::size_t i = 0; i < k; ++i) {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(k);
        const Point3 tip{std::cos(angle), std::sin(angle), 0.5};
        const Stereon::Plane plane = PlaneThrough({0, 0, 0}, {0, 0, 1}, tip);
        const BuilderVertexId w = body.AddVertex(tip);
        const BuilderFaceId face = body.AddFace(Stereon::Surface(plane));
        Must(body.Builder.AddFaceToShell(shell, face));
        const BuilderLoopId loop = body.AddLoop(face);
        const auto side = [&](const BuilderVertexId a, const BuilderVertexId b, const Point3& pa,
                              const Point3& pb) {
            const BuilderEdgeId e = body.AddEdge({.Kind = EdgeCurveKind::Boundary,
                                                  .Start = a,
                                                  .End = b,
                                                  .Range = Interval(0.0, 1.0),
                                                  .Tolerance = Tight});
            const Point2 from = UvOn(plane, pa);
            body.AddCoedge(loop, e, EdgeUseSense::Same,
                           Stereon::Curve2d(Stereon::Line2d{from, UvOn(plane, pb) - from}));
        };
        const Point2 from = UvOn(plane, {0, 0, 0});
        body.AddCoedge(loop, spine, EdgeUseSense::Same,
                       Stereon::Curve2d(Stereon::Line2d{from, UvOn(plane, {0, 0, 1}) - from}));
        side(v1, w, {0, 0, 1}, tip);
        side(w, v0, tip, {0, 0, 0});
    }
    return body;
}

/// A polyline of n Spatial edges, no faces: a wire body.
inline Body MakeWire(const std::size_t n)
{
    Body body;
    for (std::size_t i = 0; i <= n; ++i) {
        body.AddVertex({static_cast<double>(i), static_cast<double>(i % 2), 0});
    }
    for (std::size_t i = 0; i < n; ++i) {
        const Point3 a{static_cast<double>(i), static_cast<double>(i % 2), 0};
        const Point3 b{static_cast<double>(i + 1), static_cast<double>((i + 1) % 2), 0};
        body.AddEdge({.Kind = EdgeCurveKind::Spatial,
                      .Start = body.Vertices[i],
                      .End = body.Vertices[i + 1],
                      .Curve = Must(body.Builder.AddCurve(Stereon::Curve(Stereon::Line{a, b - a}))),
                      .Range = Interval(0.0, 1.0),
                      .Tolerance = Tight});
    }
    return body;
}

/// One face whose loop runs v0 -> v1 -> v0 along two Boundary edges, then
/// through a Degenerate edge at v0 (a pole), structurally: the pcurves are
/// placeholders, so the Geometric level refuses it. MakePoleCap is the
/// geometrically valid pole.
inline Body MakeDegenerateFace()
{
    Body body;
    const BuilderVertexId v0 = body.AddVertex({0, 0, 0});
    const BuilderVertexId v1 = body.AddVertex({1, 0, 0});
    const BuilderFaceId face = body.AddFace(XyPlane());
    const BuilderLoopId loop = body.AddLoop(face);
    const auto edge = [&](const EdgeCurveKind kind, const BuilderVertexId a,
                          const BuilderVertexId b) {
        return body.AddEdge(
            {.Kind = kind, .Start = a, .End = b, .Range = Interval(0.0, 1.0), .Tolerance = Tight});
    };
    const BuilderEdgeId out = edge(EdgeCurveKind::Boundary, v0, v1);
    const BuilderEdgeId back = edge(EdgeCurveKind::Boundary, v1, v0);
    const BuilderEdgeId pole = edge(EdgeCurveKind::Degenerate, v0, v0);
    body.AddCoedge(loop, out, EdgeUseSense::Same, AnyPCurve());
    body.AddCoedge(loop, back, EdgeUseSense::Same, AnyPCurve());
    body.AddCoedge(loop, pole, EdgeUseSense::Same, AnyPCurve());
    return body;
}

/// The northern hemisphere of the unit sphere about the origin, one face
/// whose loop runs in (u, v) as the cylinder's lateral loop does, with the
/// top side collapsed into the pole: (0, 0) -> (2 pi, 0) along the equator,
/// a Boundary edge closed at vertex 0 = (1, 0, 0); up the seam at u = 2 pi
/// to vertex 1 = (0, 0, 1); along the pole, a Degenerate edge at vertex 1,
/// from u = 2 pi back to u = 0; and down the seam at u = 0. Edges: 0 the
/// equator, 1 the seam (t in [0, pi/2]), 2 the pole (t in [0, 2 pi]).
inline Body MakePoleCap()
{
    using Stereon::Curve2d;
    using Stereon::Line2d;
    constexpr double tau = Stereon::CirclePeriod;
    constexpr double pole = Stereon::SpherePoleLatitude;
    Body body;
    const BuilderVertexId v0 = body.AddVertex({1, 0, 0});
    const BuilderVertexId v1 = body.AddVertex({0, 0, 1});
    const auto edge = [&](const EdgeCurveKind kind, const BuilderVertexId a,
                          const BuilderVertexId b, const double upper) {
        return body.AddEdge({.Kind = kind,
                             .Start = a,
                             .End = b,
                             .Range = Interval(0.0, upper),
                             .Tolerance = Tight});
    };
    const BuilderEdgeId equator = edge(EdgeCurveKind::Boundary, v0, v0, tau);
    const BuilderEdgeId seam = edge(EdgeCurveKind::Boundary, v0, v1, pole);
    const BuilderEdgeId apex = edge(EdgeCurveKind::Degenerate, v1, v1, tau);
    const BuilderFaceId face = body.AddFace(
        Stereon::Surface(Stereon::Sphere{{0, 0, 0}, 1.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    const BuilderLoopId loop = body.AddLoop(face);
    body.AddCoedge(loop, equator, EdgeUseSense::Same, Curve2d(Line2d{{0, 0}, Vec2{1, 0}}));
    body.AddCoedge(loop, seam, EdgeUseSense::Same, Curve2d(Line2d{{tau, 0}, Vec2{0, 1}}));
    body.AddCoedge(loop, apex, EdgeUseSense::Same, Curve2d(Line2d{{tau, pole}, Vec2{-1, 0}}));
    body.AddCoedge(loop, seam, EdgeUseSense::Opposite, Curve2d(Line2d{{0, 0}, Vec2{0, 1}}));
    return body;
}

/// A plane-curve segment from `from` to `to` on t in [0, 1].
inline Stereon::Curve2d Segment2(const Point2& from, const Point2& to)
{
    return Stereon::Curve2d(Stereon::Line2d{from, to - from});
}

/// Two unit squares meeting along the y axis, an open book: face 0 in the
/// plane z = 0 (x >= 0) and face 1 in the plane x = 0 (z >= 0), one open
/// shell. Edge 0 is the hinge, an Intersection edge from (0, 0, 0) to
/// (0, 1, 0) on t in [0, 1]; edges 1-3 and 4-6 are Boundary sides of faces 0
/// and 1. The hinge's first use, so its primary coedge, is face 0's, with
/// pcurve `hingeOnA` on `surfaceA`; by default face 0's plane, on which the
/// pcurve (0, t) gives (0, t, 0). `between` runs after face 0 is built and
/// before face 1's use of the hinge is appended.
template <typename Between>
Body MakeHinge(const Stereon::Surface& surfaceA, const Stereon::Curve2d& hingeOnA, Between between)
{
    Body body;
    const BuilderVertexId v0 = body.AddVertex({0, 0, 0});
    const BuilderVertexId v1 = body.AddVertex({0, 1, 0});
    const BuilderEdgeId hinge = body.AddEdge({.Kind = EdgeCurveKind::Intersection,
                                              .Start = v0,
                                              .End = v1,
                                              .Range = Interval(0.0, 1.0),
                                              .Tolerance = Tight});
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    const auto side = [&](const BuilderVertexId a, const BuilderVertexId b) {
        return body.AddEdge({.Kind = EdgeCurveKind::Boundary,
                             .Start = a,
                             .End = b,
                             .Range = Interval(0.0, 1.0),
                             .Tolerance = Tight});
    };
    // Face 0, in its plane's (x, y): v0 -> a1 -> a2 -> v1, then the hinge back.
    {
        const BuilderVertexId a1 = body.AddVertex({1, 0, 0});
        const BuilderVertexId a2 = body.AddVertex({1, 1, 0});
        const BuilderFaceId face = body.AddFace(surfaceA);
        Must(body.Builder.AddFaceToShell(shell, face));
        const BuilderLoopId loop = body.AddLoop(face);
        body.AddCoedge(loop, hinge, EdgeUseSense::Opposite, hingeOnA);
        body.AddCoedge(loop, side(v0, a1), EdgeUseSense::Same, Segment2({0, 0}, {1, 0}));
        body.AddCoedge(loop, side(a1, a2), EdgeUseSense::Same, Segment2({1, 0}, {1, 1}));
        body.AddCoedge(loop, side(a2, v1), EdgeUseSense::Same, Segment2({1, 1}, {0, 1}));
    }
    between(body);
    // Face 1, in its plane's (z, y): the hinge, then v1 -> b2 -> b1 -> v0.
    {
        const BuilderVertexId b2 = body.AddVertex({0, 1, 1});
        const BuilderVertexId b1 = body.AddVertex({0, 0, 1});
        const BuilderFaceId face =
            body.AddFace(Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{0, 0, 1}, Vec3{0, 1, 0}}));
        Must(body.Builder.AddFaceToShell(shell, face));
        const BuilderLoopId loop = body.AddLoop(face);
        body.AddCoedge(loop, hinge, EdgeUseSense::Same, Segment2({0, 0}, {0, 1}));
        body.AddCoedge(loop, side(v1, b2), EdgeUseSense::Same, Segment2({0, 1}, {1, 1}));
        body.AddCoedge(loop, side(b2, b1), EdgeUseSense::Same, Segment2({1, 1}, {1, 0}));
        body.AddCoedge(loop, side(b1, v0), EdgeUseSense::Same, Segment2({1, 0}, {0, 0}));
    }
    return body;
}

inline Body MakeHinge()
{
    return MakeHinge(XyPlane(), Segment2({0, 0}, {0, 1}), [](Body&) {});
}

/// The vertex a coedge starts from, along its loop: the edge's start for a
/// use in the same sense, its end otherwise.
inline Stereon::VertexId StartOfUse(const Stereon::Shape& shape, const Stereon::CoedgeId c)
{
    const Stereon::EdgeId e = EdgeOf(shape, c);
    return StoredEdgeUseSense(shape, c) == EdgeUseSense::Same ? StartOf(shape, e) : EndOf(shape, e);
}

inline Stereon::VertexId EndOfUse(const Stereon::Shape& shape, const Stereon::CoedgeId c)
{
    const Stereon::EdgeId e = EdgeOf(shape, c);
    return StoredEdgeUseSense(shape, c) == EdgeUseSense::Same ? EndOf(shape, e) : StartOf(shape, e);
}

/// The pcurve of a coedge, which every body here gives its uses; the test
/// fails on one without.
inline Stereon::Curve2d PCurveOfUse(const Stereon::Shape& shape, const Stereon::CoedgeId c)
{
    const std::optional<Stereon::Curve2d> pcurve = Stereon::PCurveOf(shape, c);
    if (!pcurve) {
        FAIL(Stereon::Describe(shape, c) << " has no pcurve");
        std::abort(); // FAIL throws; this keeps the compiler from seeing a null path
    }
    return *pcurve;
}

/// The object a pointer the test expects to be set points to; the test fails
/// on null. Dereferencing the call directly lets GCC 16.2 at -O3 see a null
/// path and warn (-Wnull-dereference).
template <typename T>
const T& Deref(const T* pointer)
{
    if (pointer == nullptr) {
        FAIL("unexpected null pointer");
        std::abort(); // FAIL throws; this keeps the compiler from seeing a null path
    }
    return *pointer;
}

/// The value of an optional the test expects to be set; the test fails on
/// an empty one. The geometry accessors return their geometry this way.
template <typename T>
T Deref(const std::optional<T>& value)
{
    if (!value.has_value()) {
        FAIL("unexpected empty value");
        std::abort(); // FAIL throws
    }
    return *value;
}

/// Collects a range into a vector.
template <typename Range>
auto Collect(const Range& range)
{
    std::vector<std::ranges::range_value_t<Range>> out;
    for (const auto h : range) {
        out.push_back(h);
    }
    return out;
}

} // namespace TopoTest
