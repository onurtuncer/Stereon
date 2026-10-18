// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The hand-built bodies of TODO.md section 7, a box (Exact edges) and a
// cylinder (Exact circles and a Boundary seam), read back through every
// accessor. Beyond the structure the validator certifies, the tests check
// the associations of ADR-0025: every pcurve agrees with its edge, loops
// close in parameter space, faces point outwards.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <map>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::FaceId;
using Stereon::LoopId;
using Stereon::Shape;
using Stereon::VertexId;

namespace {

/// The point of an edge at t: its own curve for Exact and Spatial edges,
/// the primary coedge's pcurve on its face's surface otherwise (ADR-0009).
Point3 EdgePoint(const Shape& shape, const EdgeId e, const double t)
{
    if (const Stereon::Curve* curve = Stereon::CurveOf(shape, e)) {
        return Eval(*curve, t);
    }
    const CoedgeId primary = Stereon::PrimaryCoedgeOf(shape, e);
    const FaceId face = Stereon::FaceOf(shape, Stereon::LoopOf(shape, primary));
    return Eval(Stereon::SurfaceOf(shape, face), Eval(PCurveOfUse(shape, primary), t));
}

/// Every use of every edge lies on the edge, at the same parameter, and its
/// ends are at the edge's vertices.
void RequireAssociationsAgree(const Shape& shape, const double tolerance)
{
    for (const CoedgeId c : Must(Stereon::AllEntities<EntityKind::Coedge>(shape))) {
        const EdgeId e = Stereon::EdgeOf(shape, c);
        const FaceId face = Stereon::FaceOf(shape, Stereon::LoopOf(shape, c));
        const Stereon::Surface& surface = Stereon::SurfaceOf(shape, face);
        const Stereon::Interval range = Stereon::RangeOf(shape, e);
        for (int i = 0; i <= 16; ++i) {
            const double t = range.GetLower() + (range.GetUpper() - range.GetLower()) * i / 16.0;
            const Point3 onFace = Eval(surface, Eval(PCurveOfUse(shape, c), t));
            INFO(Stereon::Describe(shape, c) << " at t = " << t);
            REQUIRE(Distance(onFace, EdgePoint(shape, e, t)) <= tolerance);
        }
    }
    for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
        const Stereon::Interval range = Stereon::RangeOf(shape, e);
        INFO(Stereon::Describe(shape, e));
        REQUIRE(Distance(EdgePoint(shape, e, range.GetLower()),
                         Stereon::PositionOf(shape, Stereon::StartOf(shape, e))) <= tolerance);
        REQUIRE(Distance(EdgePoint(shape, e, range.GetUpper()),
                         Stereon::PositionOf(shape, Stereon::EndOf(shape, e))) <= tolerance);
    }
}

/// Each coedge of a loop ends, in the face's parameter plane, where the
/// next begins: the loop closes in (u, v), seams included.
void RequireLoopsCloseInParameterSpace(const Shape& shape, const double tolerance)
{
    const auto uvAt = [&](const CoedgeId c, const bool atEnd) {
        const Stereon::Interval range = Stereon::RangeOf(shape, Stereon::EdgeOf(shape, c));
        const bool same = Stereon::StoredEdgeUseSense(shape, c) == EdgeUseSense::Same;
        return Eval(PCurveOfUse(shape, c), (atEnd == same) ? range.GetUpper() : range.GetLower());
    };
    for (const LoopId l : Must(Stereon::AllEntities<EntityKind::Loop>(shape))) {
        const auto loop = Collect(Must(Stereon::CoedgesOf(shape, l)));
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const CoedgeId next = loop[(i + 1) % loop.size()];
            INFO(Stereon::Describe(shape, loop[i]) << " -> " << Stereon::Describe(shape, next));
            REQUIRE(Distance(uvAt(loop[i], true), uvAt(next, false)) <= tolerance);
            REQUIRE(EndOfUse(shape, loop[i]) == StartOfUse(shape, next));
        }
    }
}

Vec3 FaceNormal(const Shape& shape, const FaceId f, const Point2& uv)
{
    const Stereon::Surface& s = Stereon::SurfaceOf(shape, f);
    const Vec3 n = Normalized(Cross(Deriv(s, uv, 1, 0), Deriv(s, uv, 0, 1)));
    return Stereon::IsReversed(shape, f) ? -n : n;
}

} // namespace

TEST_CASE("the hand-built box", "[topo][bodies][box]")
{
    Body body = MakeBox(2.0, 3.0, 5.0);
    const Shape shape = std::move(body).Freeze();

    SECTION("counts and the Euler characteristic of a sphere-like solid")
    {
        REQUIRE(shape.Count<EntityKind::Vertex>() == 8);
        REQUIRE(shape.Count<EntityKind::Edge>() == 12);
        REQUIRE(shape.Count<EntityKind::Coedge>() == 24);
        REQUIRE(shape.Count<EntityKind::Loop>() == 6);
        REQUIRE(shape.Count<EntityKind::Face>() == 6);
        REQUIRE(shape.Count<EntityKind::Shell>() == 1);
        REQUIRE(shape.Count<EntityKind::Solid>() == 1);
        REQUIRE(static_cast<int>(shape.Count<EntityKind::Vertex>()) -
                    static_cast<int>(shape.Count<EntityKind::Edge>()) +
                    static_cast<int>(shape.Count<EntityKind::Face>()) ==
                2);
    }
    SECTION("every edge is Exact, used twice in opposite senses by two faces")
    {
        for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
            REQUIRE(Stereon::KindOf(shape, e) == EdgeCurveKind::Exact);
            REQUIRE(Stereon::ToleranceOf(shape, e) == 1e-7);
            const auto ring = Collect(Must(Stereon::RadialCoedgesOf(shape, e)));
            REQUIRE(ring.size() == 2);
            REQUIRE(Stereon::StoredEdgeUseSense(shape, ring[0]) !=
                    Stereon::StoredEdgeUseSense(shape, ring[1]));
        }
    }
    SECTION("every face is a square-cornered quadrilateral in one shell")
    {
        const Stereon::ShellId shell = shape.Resolve(body.Shells[0]);
        REQUIRE(Stereon::SolidOf(shape, shell) == shape.Resolve(body.Solids[0]));
        const auto faces = Collect(Must(Stereon::FacesOf(shape, shell)));
        REQUIRE(faces.size() == 6);
        for (const FaceId f : faces) {
            REQUIRE(Stereon::ShellOf(shape, f) == shell);
            REQUIRE_FALSE(Stereon::IsReversed(shape, f));
            const auto loops = Collect(Must(Stereon::LoopsOf(shape, f)));
            REQUIRE(loops.size() == 1);
            REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, loops[0]))).size() == 4);
        }
    }
    SECTION("each vertex meets three edges")
    {
        std::map<VertexId, int> valence;
        for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
            ++valence[Stereon::StartOf(shape, e)];
            ++valence[Stereon::EndOf(shape, e)];
        }
        REQUIRE(valence.size() == 8);
        for (const auto& [v, n] : valence) {
            REQUIRE(n == 3);
        }
    }
    SECTION("the geometry agrees with the topology")
    {
        RequireAssociationsAgree(shape, 1e-12);
        RequireLoopsCloseInParameterSpace(shape, 1e-12);
    }
    SECTION("faces point away from the centre")
    {
        const Point3 centre{1.0, 1.5, 2.5};
        for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(shape))) {
            const Stereon::Surface& s = Stereon::SurfaceOf(shape, f);
            const Point3 p = Eval(s, Point2{0.1, 0.1});
            REQUIRE(Dot(FaceNormal(shape, f, {0.1, 0.1}), p - centre) > 0.0);
        }
    }
    SECTION("accessors read the vertices back")
    {
        REQUIRE(Stereon::PositionOf(shape, shape.Resolve(body.Vertices[6])) == Point3{2, 3, 5});
        REQUIRE(Stereon::ToleranceOf(shape, shape.Resolve(body.Vertices[0])) == 1e-7);
    }
    SECTION("Describe names the entity and its neighbourhood")
    {
        const EdgeId e = shape.Resolve(body.Edges[0]);
        const std::string edge = Stereon::Describe(shape, e);
        REQUIRE(edge.starts_with("Edge#0 [Exact Line, tol 1e-07, faces "));
        REQUIRE(edge.find('|') != std::string::npos);
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Vertices[1])) ==
                "Vertex#1 [(2, 0, 0), tol 1e-07]");
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Faces[0])) ==
                "Face#0 [Plane, Shell#0, 1 loop]");
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Loops[2])) ==
                "Loop#2 [Face#2, 4 coedges]");
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Shells[0])) ==
                "Shell#0 [Solid#0, 6 faces]");
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Solids[0])) == "Solid#0 [1 shell]");
        const std::string coedge = Stereon::Describe(shape, shape.Resolve(body.Coedges[0]));
        REQUIRE(coedge.starts_with("Coedge#0 [Edge#"));
        REQUIRE(coedge.find("Loop#0") != std::string::npos);
    }
}

TEST_CASE("the hand-built cylinder, seam included", "[topo][bodies][cylinder]")
{
    constexpr double r = 1.5;
    constexpr double h = 4.0;
    Body body = MakeCylinder(r, h);
    const Shape shape = std::move(body).Freeze();

    SECTION("counts")
    {
        REQUIRE(shape.Count<EntityKind::Vertex>() == 2);
        REQUIRE(shape.Count<EntityKind::Edge>() == 3);
        REQUIRE(shape.Count<EntityKind::Coedge>() == 6);
        REQUIRE(shape.Count<EntityKind::Loop>() == 3);
        REQUIRE(shape.Count<EntityKind::Face>() == 3);
        // V - E + F = 2 - 3 + 3: the seam keeps the count of a sphere.
        REQUIRE(static_cast<int>(shape.Count<EntityKind::Vertex>()) -
                    static_cast<int>(shape.Count<EntityKind::Edge>()) +
                    static_cast<int>(shape.Count<EntityKind::Face>()) ==
                2);
    }
    SECTION("edge kinds and uses")
    {
        const EdgeId bottom = shape.Resolve(body.Edges[0]);
        const EdgeId top = shape.Resolve(body.Edges[1]);
        const EdgeId seam = shape.Resolve(body.Edges[2]);
        REQUIRE(Stereon::KindOf(shape, bottom) == EdgeCurveKind::Exact);
        REQUIRE(Stereon::KindOf(shape, top) == EdgeCurveKind::Exact);
        REQUIRE(Stereon::KindOf(shape, seam) == EdgeCurveKind::Boundary);
        REQUIRE(Stereon::StartOf(shape, bottom) == Stereon::EndOf(shape, bottom));
        REQUIRE(Stereon::CurveOf(shape, seam) == nullptr);
        REQUIRE(Stereon::PrimaryCoedgeOf(shape, seam).IsValid());
        for (const EdgeId e : {bottom, top}) {
            const auto ring = Collect(Must(Stereon::RadialCoedgesOf(shape, e)));
            REQUIRE(ring.size() == 2);
            REQUIRE(Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[0])) !=
                    Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[1])));
        }
        REQUIRE(
            Stereon::Describe(shape, seam).starts_with("Edge#2 [Boundary, tol 1e-07, faces 0|0]"));
        REQUIRE(Stereon::Describe(shape, bottom).starts_with("Edge#0 [Exact Circle"));
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Faces[0])) ==
                "Face#0 [Cylinder, Shell#0, 1 loop]");
    }
    SECTION("the geometry agrees with the topology")
    {
        RequireAssociationsAgree(shape, 1e-12);
        RequireLoopsCloseInParameterSpace(shape, 1e-12);
    }
    SECTION("faces point outwards")
    {
        const FaceId lateral = shape.Resolve(body.Faces[0]);
        const FaceId down = shape.Resolve(body.Faces[1]);
        const FaceId up = shape.Resolve(body.Faces[2]);
        const Vec3 side = FaceNormal(shape, lateral, {0.3, 1.0});
        REQUIRE(Dot(side, Vec3{std::cos(0.3), std::sin(0.3), 0.0}) > 0.99);
        REQUIRE(FaceNormal(shape, down, {0.1, 0.1}).Z < -0.99);
        REQUIRE(FaceNormal(shape, up, {0.1, 0.1}).Z > 0.99);
    }
}

TEST_CASE("sheets, holes, fans, wires and poles are valid topology", "[topo][bodies]")
{
    SECTION("a square sheet")
    {
        Body body = MakeSquareSheet();
        const Shape shape = std::move(body).Freeze();
        RequireAssociationsAgree(shape, 1e-12);
        RequireLoopsCloseInParameterSpace(shape, 1e-12);
        REQUIRE_FALSE(Stereon::SolidOf(shape, shape.Resolve(body.Shells[0])).IsValid());
    }
    SECTION("a face with a hole")
    {
        Body body = MakeHoleSheet();
        const Shape shape = std::move(body).Freeze();
        REQUIRE(shape.Count<EntityKind::Loop>() == 2);
        REQUIRE(shape.Count<EntityKind::Edge>() == 8);
        RequireLoopsCloseInParameterSpace(shape, 1e-12);
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Faces[0])) ==
                "Face#0 [Plane, Shell#0, 2 loops]");
    }
    SECTION("a fan of seven faces around one edge")
    {
        Body body = MakeFan(7);
        const Shape shape = std::move(body).Freeze();
        RequireAssociationsAgree(shape, 1e-12);
        RequireLoopsCloseInParameterSpace(shape, 1e-12);
        REQUIRE(
            Stereon::Describe(shape, shape.Resolve(body.Edges[0])).find("faces 0|6|5|4|3|2|1") !=
            std::string::npos);
    }
    SECTION("a wire")
    {
        Body body = MakeWire(5);
        const Shape shape = std::move(body).Freeze();
        REQUIRE(shape.Count<EntityKind::Coedge>() == 0);
        RequireAssociationsAgree(shape, 1e-12);
        REQUIRE(Stereon::Describe(shape, shape.Resolve(body.Edges[0])) ==
                "Edge#0 [Spatial Line, tol 1e-07, faces none]");
    }
    SECTION("a loop through a pole")
    {
        Body body = MakeDegenerateFace();
        const Shape shape = std::move(body).Freeze();
        const EdgeId pole = shape.Resolve(body.Edges[2]);
        REQUIRE(Stereon::KindOf(shape, pole) == EdgeCurveKind::Degenerate);
        REQUIRE(Stereon::StartOf(shape, pole) == Stereon::EndOf(shape, pole));
        REQUIRE(Stereon::PrimaryCoedgeOf(shape, pole) == shape.Resolve(body.Coedges[2]));
        REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(body.Loops[0])))).size() == 3);
    }
}

TEST_CASE("two edges may share one curve over their own ranges", "[topo][bodies]")
{
    // ADR-0025, verification item 3: a shared support is not shared topology.
    Body body;
    const BuilderVertexId a = body.AddVertex({0, 0, 0});
    const BuilderVertexId b = body.AddVertex({1, 0, 0});
    const BuilderVertexId c = body.AddVertex({3, 0, 0});
    const Stereon::CurveRef line =
        body.Builder.AddCurve(Stereon::Curve(Stereon::Line{{0, 0, 0}, Vec3{1, 0, 0}}));
    const BuilderEdgeId first = body.AddEdge({.Kind = EdgeCurveKind::Spatial,
                                              .Start = a,
                                              .End = b,
                                              .Curve = line,
                                              .Range = Interval(0, 1)});
    const BuilderEdgeId second = body.AddEdge({.Kind = EdgeCurveKind::Spatial,
                                               .Start = b,
                                               .End = c,
                                               .Curve = line,
                                               .Range = Interval(1, 3)});
    // And two vertices at one point stay two vertices (item 4).
    const BuilderVertexId twin = body.AddVertex({1, 0, 0});
    const Shape shape = std::move(body).Freeze();
    REQUIRE(Stereon::CurveOf(shape, shape.Resolve(first)) ==
            Stereon::CurveOf(shape, shape.Resolve(second)));
    REQUIRE(Stereon::RangeOf(shape, shape.Resolve(first)).GetUpper() == 1.0);
    REQUIRE(Stereon::RangeOf(shape, shape.Resolve(second)).GetUpper() == 3.0);
    REQUIRE(shape.Resolve(twin) != shape.Resolve(b));
    REQUIRE(shape.Count<EntityKind::Vertex>() == 4);
    RequireAssociationsAgree(shape, 0.0);
}
