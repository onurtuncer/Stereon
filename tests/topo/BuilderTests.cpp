// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ShapeBuilder's controlled operations (ADR-0027): what each accepts, what
// it refuses and leaves unchanged, and the links it maintains.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::BuilderToken;
using Stereon::CoedgeId;
using Stereon::CurveRef;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::PCurveRef;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::SurfaceRef;

namespace {

template <typename T>
void RequireRefused(const Stereon::Result<T>& result, const std::string& fragment)
{
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    INFO(std::string(result.error().GetMessage()));
    REQUIRE(std::string(result.error().GetMessage()).find(fragment) != std::string::npos);
}

Stereon::EdgeSpec ExactLine(Body& body, const BuilderVertexId a, const BuilderVertexId b)
{
    return {.Kind = EdgeCurveKind::Exact,
            .Start = a,
            .End = b,
            .Curve = Must(
                body.Builder.AddCurve(Stereon::Curve(Stereon::Line{{0, 0, 0}, Vec3{1, 0, 0}}))),
            .Range = Interval(0.0, 1.0),
            .Tolerance = 1e-7};
}

} // namespace

TEST_CASE("AddVertex refuses positions and tolerances that are not finite", "[topo][builder]")
{
    ShapeBuilder b;
    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    RequireRefused(b.AddVertex({inf, 0, 0}, 0.0), "the position");
    RequireRefused(b.AddVertex({0, nan, 0}, 0.0), "the position");
    RequireRefused(b.AddVertex({0, 0, 0}, -1e-9), "tolerance");
    RequireRefused(b.AddVertex({0, 0, 0}, nan), "tolerance");
    RequireRefused(b.AddVertex({0, 0, 0}, inf), "tolerance");
    REQUIRE(b.Count<EntityKind::Vertex>() == 0);
    REQUIRE(b.AddVertex({1, 2, 3}, 0.0).has_value());
    REQUIRE(b.Count<EntityKind::Vertex>() == 1);
}

TEST_CASE("AddEdge checks its vertices, range, curve and kind", "[topo][builder]")
{
    Body body;
    const BuilderVertexId a = body.AddVertex({0, 0, 0});
    const BuilderVertexId c = body.AddVertex({1, 0, 0});

    SECTION("an invalid or foreign vertex handle")
    {
        Stereon::EdgeSpec spec = ExactLine(body, a, c);
        spec.Start = BuilderVertexId{};
        RequireRefused(body.Builder.AddEdge(spec), "not a live entity");
        spec.Start = BuilderVertexId::Make(57, 0, body.Builder.GetToken().GetValue());
        RequireRefused(body.Builder.AddEdge(spec), "not a live entity");
        spec.Start = BuilderVertexId::Make(a.GetIndex(), 3, body.Builder.GetToken().GetValue());
        RequireRefused(body.Builder.AddEdge(spec), "not a live entity");
        if constexpr (STEREON_CHECKED_HANDLES) {
            ShapeBuilder other;
            const BuilderVertexId foreign = Must(other.AddVertex({0, 0, 0}, 0.0));
            spec.Start = foreign;
            RequireRefused(body.Builder.AddEdge(spec), "not a live entity");
        }
    }
    SECTION("an empty range")
    {
        Stereon::EdgeSpec spec = ExactLine(body, a, c);
        spec.Range = Interval{};
        RequireRefused(body.Builder.AddEdge(spec), "empty parameter range");
    }
    SECTION("a tolerance that is negative or not finite")
    {
        Stereon::EdgeSpec spec = ExactLine(body, a, c);
        spec.Tolerance = -1.0;
        RequireRefused(body.Builder.AddEdge(spec), "tolerance");
        spec.Tolerance = std::numeric_limits<double>::quiet_NaN();
        RequireRefused(body.Builder.AddEdge(spec), "tolerance");
    }
    SECTION("Exact and Spatial edges need their curve; the others have none")
    {
        for (const EdgeCurveKind kind : {EdgeCurveKind::Exact, EdgeCurveKind::Spatial}) {
            Stereon::EdgeSpec spec = ExactLine(body, a, c);
            spec.Kind = kind;
            spec.Curve = CurveRef{};
            RequireRefused(body.Builder.AddEdge(spec), "needs its 3D curve");
        }
        for (const EdgeCurveKind kind :
             {EdgeCurveKind::Intersection, EdgeCurveKind::Trim, EdgeCurveKind::Boundary}) {
            Stereon::EdgeSpec spec = ExactLine(body, a, c);
            spec.Kind = kind;
            RequireRefused(body.Builder.AddEdge(spec), "no 3D curve of its own");
        }
    }
    SECTION("a curve that is not in the pool")
    {
        Stereon::EdgeSpec spec = ExactLine(body, a, c);
        spec.Curve = CurveRef::FromIndex(99);
        RequireRefused(body.Builder.AddEdge(spec), "not in this builder's pool");
    }
    SECTION("a degenerate edge has one vertex")
    {
        Stereon::EdgeSpec spec{
            .Kind = EdgeCurveKind::Degenerate, .Start = a, .End = c, .Range = Interval(0.0, 1.0)};
        RequireRefused(body.Builder.AddEdge(spec), "one vertex");
        spec.End = a;
        REQUIRE(body.Builder.AddEdge(spec).has_value());
    }
    SECTION("a closed edge has one vertex at both ends")
    {
        REQUIRE(body.Builder.AddEdge(ExactLine(body, a, a)).has_value());
    }
    REQUIRE(body.Builder.Count<EntityKind::Edge>() <= 1);
}

TEST_CASE("AddFace needs a surface of the pool", "[topo][builder]")
{
    ShapeBuilder b;
    RequireRefused(b.AddFace(SurfaceRef{}), "surface");
    RequireRefused(b.AddFace(SurfaceRef::FromIndex(0)), "surface");
    const SurfaceRef plane = Must(b.AddSurface(XyPlane()));
    const BuilderFaceId f = Must(b.AddFace(plane));
    const BuilderFaceId g = Must(b.AddFace(plane, true));
    REQUIRE(f != g);
    REQUIRE(b.Count<EntityKind::Face>() == 2);
}

TEST_CASE("the geometry pools refuse empty geometry and a spent builder", "[topo][builder][audit]")
{
    // Audit 2026-10-20, F3: an empty wrapper took a pool index that a face
    // or an edge then accepted, and the certified shape crashed Describe.
    ShapeBuilder b;
    RequireRefused(b.AddCurve(Stereon::Curve{}), "empty");
    RequireRefused(b.AddPCurve(Stereon::Curve2d{}), "empty");
    RequireRefused(b.AddSurface(Stereon::Surface{}), "empty");
    // Nothing was stored: the first geometry added takes index 0.
    REQUIRE(Must(b.AddSurface(XyPlane())).GetIndex() == 0);
    REQUIRE(Must(b.AddPCurve(AnyPCurve())).GetIndex() == 0);
    const Stereon::Curve line(Stereon::Line{{0, 0, 0}, Vec3{1, 0, 0}});
    REQUIRE(Must(b.AddCurve(line)).GetIndex() == 0);

    // Moved from or frozen, a builder refuses geometry as it refuses
    // everything else.
    ShapeBuilder c(std::move(b));
    RequireRefused(b.AddCurve(line), "moved from or frozen"); // NOLINT(bugprone-use-after-move)
    RequireRefused(b.AddPCurve(AnyPCurve()), "moved from or frozen");
    RequireRefused(b.AddSurface(XyPlane()), "moved from or frozen");
    REQUIRE(Stereon::Freeze(std::move(c)).has_value());
    RequireRefused(c.AddCurve(line), "moved from or frozen"); // NOLINT(bugprone-use-after-move)
    RequireRefused(c.AddPCurve(AnyPCurve()), "moved from or frozen");
    RequireRefused(c.AddSurface(XyPlane()), "moved from or frozen");
}

TEST_CASE("AppendCoedge needs a pcurve except for a Spatial edge", "[topo][builder]")
{
    Body body;
    const BuilderVertexId a = body.AddVertex({0, 0, 0});
    const BuilderFaceId face = body.AddFace(XyPlane());
    const BuilderLoopId loop = body.AddLoop(face);
    const BuilderEdgeId exact = body.AddEdge(ExactLine(body, a, a));
    RequireRefused(body.Builder.AppendCoedge(loop, exact, EdgeUseSense::Same, PCurveRef{}),
                   "needs its pcurve");
    RequireRefused(
        body.Builder.AppendCoedge(loop, exact, EdgeUseSense::Same, PCurveRef::FromIndex(4)),
        "not in this builder's pool");
    RequireRefused(body.Builder.AppendCoedge(BuilderLoopId{}, exact, EdgeUseSense::Same,
                                             Must(body.Builder.AddPCurve(AnyPCurve()))),
                   "not a live entity");
    REQUIRE(body.Builder.Count<EntityKind::Coedge>() == 0);

    Stereon::EdgeSpec spatial = ExactLine(body, a, a);
    spatial.Kind = EdgeCurveKind::Spatial;
    const BuilderEdgeId wire = body.AddEdge(spatial);
    REQUIRE(body.Builder.AppendCoedge(loop, wire, EdgeUseSense::Same, PCurveRef{}).has_value());
}

TEST_CASE("coedges join the loop cycle in order and the radial ring of their edge",
          "[topo][builder][ranges]")
{
    // A square sheet's loop, then the same face given a second loop over the
    // same edges in reverse: every edge then has two uses.
    Body body = MakeSquareSheet();
    const Shape shape = std::move(body).Freeze();
    const auto loops = Collect(Must(Stereon::LoopsOf(shape, shape.Resolve(body.Faces[0]))));
    REQUIRE(loops.size() == 1);
    const auto coedges = Collect(Must(Stereon::CoedgesOf(shape, loops[0])));
    REQUIRE(coedges.size() == 4);
    for (std::size_t i = 0; i < 4; ++i) {
        REQUIRE(coedges[i] == shape.Resolve(body.Coedges[i]));
        // Consecutive uses meet at a vertex.
        REQUIRE(EndOfUse(shape, coedges[i]) == StartOfUse(shape, coedges[(i + 1) % 4]));
    }
}

TEST_CASE("the first use of a surface-defined edge becomes its primary coedge",
          "[topo][builder][primary]")
{
    SECTION("Boundary edges")
    {
        Body body = MakeSquareSheet();
        const Shape shape = std::move(body).Freeze();
        for (std::size_t i = 0; i < 4; ++i) {
            REQUIRE(Stereon::PrimaryCoedgeOf(shape, shape.Resolve(body.Edges[i])) ==
                    shape.Resolve(body.Coedges[i]));
        }
    }
    SECTION("Exact edges have their own curve and no primary")
    {
        Body body = MakeBox();
        const Shape shape = std::move(body).Freeze();
        for (const BuilderEdgeId e : body.Edges) {
            REQUIRE_FALSE(Stereon::PrimaryCoedgeOf(shape, shape.Resolve(e)).IsValid());
            REQUIRE(Stereon::CurveOf(shape, shape.Resolve(e)).has_value());
        }
    }
    SECTION("the seam's first use, at u = 2 pi")
    {
        Body body = MakeCylinder();
        const Shape shape = std::move(body).Freeze();
        REQUIRE(Stereon::PrimaryCoedgeOf(shape, shape.Resolve(body.Edges[2])) ==
                shape.Resolve(body.Coedges[1]));
    }
}

TEST_CASE("SetPrimaryCoedge moves the primary without touching the ring",
          "[topo][builder][primary]")
{
    Body body = MakeCylinder();
    const BuilderEdgeId seam = body.Edges[2];
    RequireRefused(body.Builder.SetPrimaryCoedge(seam, body.Coedges[0]), "not a use of the edge");
    RequireRefused(body.Builder.SetPrimaryCoedge(seam, BuilderCoedgeId{}), "not a live entity");
    const std::uint64_t epoch = body.Builder.GetEpoch();
    Must(body.Builder.SetPrimaryCoedge(seam, body.Coedges[3]));
    REQUIRE(body.Builder.GetEpoch() > epoch);
    const Shape shape = std::move(body).Freeze();
    REQUIRE(Stereon::PrimaryCoedgeOf(shape, shape.Resolve(seam)) == shape.Resolve(body.Coedges[3]));
    const auto ring = Collect(Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(seam))));
    REQUIRE(ring ==
            std::vector<CoedgeId>{shape.Resolve(body.Coedges[1]), shape.Resolve(body.Coedges[3])});
}

TEST_CASE("shells and solids keep their members in order", "[topo][builder]")
{
    Body body;
    const BuilderSolidId solid = Must(body.Builder.AddSolid());
    const BuilderShellId outer = Must(body.Builder.AddShell(solid));
    const BuilderShellId inner = Must(body.Builder.AddShell(solid));
    const BuilderShellId free = Must(body.Builder.AddShell());
    RequireRefused(body.Builder.AddShell(BuilderSolidId{}), "not a live entity");

    const BuilderFaceId f = body.AddFace(XyPlane());
    Must(body.Builder.AddFaceToShell(free, f));
    RequireRefused(body.Builder.AddFaceToShell(outer, f), "already in a shell");
    RequireRefused(body.Builder.AddFaceToShell(BuilderShellId{}, f), "not a live entity");
    AddBoundaryPolygon(body, f, std::vector<Point3>{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});

    const Shape shape = std::move(body).Freeze();
    const auto shells = Collect(Must(Stereon::ShellsOf(shape, shape.Resolve(solid))));
    REQUIRE(shells == std::vector{shape.Resolve(outer), shape.Resolve(inner)});
    REQUIRE(Stereon::SolidOf(shape, shape.Resolve(outer)) == shape.Resolve(solid));
    REQUIRE_FALSE(Stereon::SolidOf(shape, shape.Resolve(free)).IsValid());
    REQUIRE(Collect(Must(Stereon::FacesOf(shape, shape.Resolve(free)))) ==
            std::vector{shape.Resolve(f)});
    REQUIRE(Collect(Must(Stereon::FacesOf(shape, shape.Resolve(outer)))).empty());
    REQUIRE(Stereon::ShellOf(shape, shape.Resolve(f)) == shape.Resolve(free));
}

TEST_CASE("every change advances the epoch; refusals do not", "[topo][builder][epoch]")
{
    ShapeBuilder b;
    std::uint64_t epoch = b.GetEpoch();
    const auto advanced = [&] {
        const bool moved = b.GetEpoch() > epoch;
        epoch = b.GetEpoch();
        return moved;
    };
    const BuilderVertexId v = Must(b.AddVertex({0, 0, 0}, 0.0));
    REQUIRE(advanced());
    const SurfaceRef s = Must(b.AddSurface(XyPlane()));
    REQUIRE_FALSE(advanced()); // geometry pools are not topology
    const BuilderFaceId f = Must(b.AddFace(s));
    REQUIRE(advanced());
    const BuilderLoopId l = Must(b.AddLoop(f));
    REQUIRE(advanced());
    const BuilderEdgeId e = Must(b.AddEdge(
        {.Kind = EdgeCurveKind::Boundary, .Start = v, .End = v, .Range = Interval(0.0, 1.0)}));
    REQUIRE(advanced());
    Must(b.AppendCoedge(l, e, EdgeUseSense::Same, Must(b.AddPCurve(AnyPCurve()))));
    REQUIRE(advanced());
    const BuilderShellId sh = Must(b.AddShell());
    REQUIRE(advanced());
    Must(b.AddFaceToShell(sh, f));
    REQUIRE(advanced());
    REQUIRE_FALSE(b.AddFaceToShell(sh, f).has_value());
    REQUIRE_FALSE(advanced());
    REQUIRE_FALSE(b.AddVertex({0, 0, 0}, -1.0).has_value());
    REQUIRE_FALSE(advanced());
}

TEST_CASE("the token moves with the builder and the source refuses everything",
          "[topo][builder][token]")
{
    ShapeBuilder a;
    const BuilderToken token = a.GetToken();
    REQUIRE(token.IsValid());
    const BuilderVertexId v = Must(a.AddVertex({0, 0, 0}, 0.0));
    const std::uint64_t epoch = a.GetEpoch();

    ShapeBuilder b(std::move(a));
    REQUIRE(b.GetToken() == token);
    REQUIRE(b.Count<EntityKind::Vertex>() == 1);
    REQUIRE_FALSE(a.GetToken().IsValid()); // NOLINT(bugprone-use-after-move): the point
    REQUIRE(a.GetEpoch() > epoch);         // NOLINT(bugprone-use-after-move)
    RequireRefused(a.AddVertex({0, 0, 0}, 0.0), "moved from");
    RequireRefused(a.AddSolid(), "moved from");
    RequireRefused(a.AddFace(SurfaceRef{}), "moved from");
    RequireRefused(Stereon::Freeze(std::move(a)), "moved from");

    ShapeBuilder c;
    const BuilderToken other = c.GetToken();
    REQUIRE(other != token);
    c = std::move(b);
    REQUIRE(c.GetToken() == token);
    REQUIRE(c.Count<EntityKind::Vertex>() == 1);
    REQUIRE(Must(c.AddEdge({.Kind = EdgeCurveKind::Boundary,
                            .Start = v,
                            .End = v,
                            .Range = Interval(0.0, 1.0)}))
                .IsValid());
}

TEST_CASE("each builder has its own token", "[topo][builder][token]")
{
    std::vector<ShapeBuilder> builders(16);
    for (std::size_t i = 0; i < builders.size(); ++i) {
        for (std::size_t j = i + 1; j < builders.size(); ++j) {
            REQUIRE(builders[i].GetToken() != builders[j].GetToken());
        }
    }
}

TEST_CASE("Count follows every kind of entity", "[topo][builder]")
{
    Body box = MakeBox();
    REQUIRE(box.Builder.Count<EntityKind::Vertex>() == 8);
    REQUIRE(box.Builder.Count<EntityKind::Edge>() == 12);
    REQUIRE(box.Builder.Count<EntityKind::Coedge>() == 24);
    REQUIRE(box.Builder.Count<EntityKind::Loop>() == 6);
    REQUIRE(box.Builder.Count<EntityKind::Face>() == 6);
    REQUIRE(box.Builder.Count<EntityKind::Shell>() == 1);
    REQUIRE(box.Builder.Count<EntityKind::Solid>() == 1);
    REQUIRE(box.Builder.Validate().IsValid());
}
