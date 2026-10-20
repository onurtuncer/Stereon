// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Editing: the Remove operations and SetSurface, what they refuse, the
// tombstones they leave, builders opened on a frozen Shape and the pages
// an edit shares with it (ADR-0004, ADR-0005, ADR-0027).

#include "BodiesStd.hpp"
#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::FaceId;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::TopologyRule;
using Stereon::VertexId;

namespace {

std::string MessageOf(const Stereon::Result<void>& result)
{
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    return std::string(result.error().GetMessage());
}

/// The uses of every edge, in slot order.
std::vector<std::size_t> UseCounts(const Shape& shape)
{
    std::vector<std::size_t> counts;
    for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
        counts.push_back(Collect(Must(Stereon::RadialCoedgesOf(shape, e))).size());
    }
    return counts;
}

/// Two planar triangles (a, b, c) and (b, a, d) sharing the Intersection
/// edge ab; the other sides are Boundary edges. Faces 0 and 1, edge 0 the
/// shared one, coedge 0 its first use (the primary), coedge 3 its second.
Body MakeTwoTriangles()
{
    Body body;
    const BuilderVertexId a = body.AddVertex({0, 0, 0});
    const BuilderVertexId b = body.AddVertex({1, 0, 0});
    const BuilderVertexId c = body.AddVertex({0, 1, 0});
    const BuilderVertexId d = body.AddVertex({0, -1, 0.5});
    const auto edge = [&](const EdgeCurveKind kind, const BuilderVertexId from,
                          const BuilderVertexId to) {
        return body.AddEdge({.Kind = kind,
                             .Start = from,
                             .End = to,
                             .Range = Interval(0.0, 1.0),
                             .Tolerance = 1e-7});
    };
    const BuilderEdgeId ab = edge(EdgeCurveKind::Intersection, a, b);
    const BuilderEdgeId bc = edge(EdgeCurveKind::Boundary, b, c);
    const BuilderEdgeId ca = edge(EdgeCurveKind::Boundary, c, a);
    const BuilderEdgeId ad = edge(EdgeCurveKind::Boundary, a, d);
    const BuilderEdgeId db = edge(EdgeCurveKind::Boundary, d, b);
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    const BuilderFaceId first = body.AddFace(XyPlane());
    const BuilderLoopId one = body.AddLoop(first);
    for (const BuilderEdgeId e : {ab, bc, ca}) {
        body.AddCoedge(one, e, EdgeUseSense::Same, AnyPCurve());
    }
    const BuilderFaceId second = body.AddFace(
        Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, -0.8, 0.6}}));
    const BuilderLoopId two = body.AddLoop(second);
    body.AddCoedge(two, ab, EdgeUseSense::Opposite, AnyPCurve());
    body.AddCoedge(two, ad, EdgeUseSense::Same, AnyPCurve());
    body.AddCoedge(two, db, EdgeUseSense::Same, AnyPCurve());
    Must(body.Builder.AddFaceToShell(shell, first));
    Must(body.Builder.AddFaceToShell(shell, second));
    return body;
}

} // namespace

TEST_CASE("removing a face takes its loops and coedges and leaves its edges", "[topo][edit]")
{
    Body box = MakeBox();
    ShapeBuilder& b = box.Builder;
    const std::uint64_t epoch = b.GetEpoch();
    Must(b.RemoveFace(box.Faces[0]));
    REQUIRE(b.GetEpoch() > epoch);
    REQUIRE(b.Count<EntityKind::Face>() == 5);
    REQUIRE(b.Count<EntityKind::Loop>() == 5);
    REQUIRE(b.Count<EntityKind::Coedge>() == 20);
    REQUIRE(b.Count<EntityKind::Edge>() == 12);
    REQUIRE(b.Count<EntityKind::Vertex>() == 8);
    REQUIRE(b.Validate().IsValid());
    // The removed face's handle is stale for every operation.
    REQUIRE(MessageOf(b.RemoveFace(box.Faces[0])).find("not a live entity") != std::string::npos);
    REQUIRE_FALSE(b.AddLoop(box.Faces[0]).has_value());
    REQUIRE_FALSE(b.RemoveCoedge(box.Coedges[0]).has_value());
    const Shape shape = std::move(box).Freeze();
    REQUIRE_FALSE(shape.Resolve(box.Faces[0]).IsValid());
    REQUIRE(shape.Resolve(box.Faces[1]).IsValid());
    // Four edges have one use now, the other eight still two.
    std::size_t single = 0;
    for (const std::size_t uses : UseCounts(shape)) {
        REQUIRE((uses == 1 || uses == 2));
        single += uses == 1 ? 1 : 0;
    }
    REQUIRE(single == 4);
    REQUIRE(Stereon::ValidateTopology(shape).IsValid());
}

TEST_CASE("removing a loop or a coedge keeps the cycle and the ring closed", "[topo][edit]")
{
    SECTION("the hole of a sheet, then its edges and vertices")
    {
        Body sheet = MakeHoleSheet();
        ShapeBuilder& b = sheet.Builder;
        Must(b.RemoveLoop(sheet.Loops[1]));
        REQUIRE(b.Count<EntityKind::Loop>() == 1);
        REQUIRE(b.Count<EntityKind::Coedge>() == 4);
        // The hole's Boundary edges have no use, which their kind forbids...
        Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::EdgeUseCount));
        REQUIRE(report.GetCount() == 4);
        // ...so they go, and then their vertices.
        for (std::size_t e = 4; e < 8; ++e) {
            Must(b.RemoveEdge(sheet.Edges[e]));
        }
        for (std::size_t v = 4; v < 8; ++v) {
            Must(b.RemoveVertex(sheet.Vertices[v]));
        }
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(sheet).Freeze();
        REQUIRE(shape.Count<EntityKind::Edge>() == 4);
        REQUIRE(shape.Count<EntityKind::Vertex>() == 4);
        REQUIRE(Collect(Must(Stereon::LoopsOf(shape, shape.Resolve(sheet.Faces[0])))).size() == 1);
    }
    SECTION("one coedge of a box face")
    {
        Body box = MakeBox();
        Must(box.Builder.RemoveCoedge(box.Coedges[0]));
        REQUIRE(box.Builder.Validate().IsValid());
        const Shape shape = std::move(box).Freeze();
        REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(box.Loops[0])))).size() == 3);
        REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(box.Loops[1])))).size() == 4);
        std::size_t single = 0;
        for (const std::size_t uses : UseCounts(shape)) {
            single += uses == 1 ? 1 : 0;
        }
        REQUIRE(single == 1);
    }
    SECTION("the anchor of a loop and of a ring")
    {
        // Coedge 0 anchors loop 0; removing it moves the anchor on, and the
        // walk still visits the rest once.
        Body box = MakeBox();
        Must(box.Builder.RemoveCoedge(box.Coedges[0]));
        Must(box.Builder.RemoveCoedge(box.Coedges[1]));
        const Shape shape = std::move(box).Freeze();
        REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(box.Loops[0])))).size() == 2);
        REQUIRE(Stereon::ValidateTopology(shape).IsValid());
    }
}

TEST_CASE("a primary coedge is not removed alone, and passes on with its loop",
          "[topo][edit][primary]")
{
    SECTION("the seam of a cylinder")
    {
        Body cylinder = MakeCylinder();
        ShapeBuilder& b = cylinder.Builder;
        const BuilderEdgeId seam = cylinder.Edges[2];
        const std::string refused = MessageOf(b.RemoveCoedge(cylinder.Coedges[1]));
        REQUIRE(refused.find("primary of Edge#2") != std::string::npos);
        REQUIRE(refused.find("SetPrimaryCoedge") != std::string::npos);
        REQUIRE(b.Count<EntityKind::Coedge>() == 6);
        Must(b.SetPrimaryCoedge(seam, cylinder.Coedges[3]));
        Must(b.RemoveCoedge(cylinder.Coedges[1]));
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(cylinder).Freeze();
        REQUIRE(PrimaryCoedgeOf(shape, shape.Resolve(seam)) == shape.Resolve(cylinder.Coedges[3]));
        REQUIRE(Collect(Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(seam)))).size() == 1);
    }
    SECTION("a loop takes its primaries with it")
    {
        Body triangles = MakeTwoTriangles();
        ShapeBuilder& b = triangles.Builder;
        REQUIRE(b.Validate().IsValid());
        Must(b.RemoveFace(triangles.Faces[0]));
        // The shared edge keeps a primary, the second face's use, and is
        // left with one use, which an Intersection edge may not have...
        const Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::EdgeUseCount));
        REQUIRE_FALSE(report.Has(TopologyRule::PrimaryCoedge));
        // ...until it is made a Boundary edge of the remaining face.
        Stereon::Testing::TopologyFaults::SetKind(b, triangles.Edges[0], EdgeCurveKind::Boundary);
        for (const std::size_t e : {1, 2}) {
            Must(b.RemoveEdge(triangles.Edges[e]));
        }
        Must(b.RemoveVertex(triangles.Vertices[2]));
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(triangles).Freeze();
        REQUIRE(PrimaryCoedgeOf(shape, shape.Resolve(triangles.Edges[0])) ==
                shape.Resolve(triangles.Coedges[3]));
    }
    SECTION("an edge left with no use loses its primary")
    {
        Body fan = MakeFan(3);
        ShapeBuilder& b = fan.Builder;
        Must(b.RemoveFace(fan.Faces[0]));
        // Its two Boundary sides have no use now: they must go.
        REQUIRE(b.Validate().Has(TopologyRule::EdgeUseCount));
        Must(b.RemoveEdge(fan.Edges[1]));
        Must(b.RemoveEdge(fan.Edges[2]));
        Must(b.RemoveVertex(fan.Vertices[2]));
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(fan).Freeze();
        REQUIRE(Must(Stereon::FacesOf(shape, shape.Resolve(fan.Edges[0]))).size() == 2);
    }
}

TEST_CASE("removals that would leave a dangling link are refused", "[topo][edit]")
{
    Body box = MakeBox();
    ShapeBuilder& b = box.Builder;
    const std::string edge = MessageOf(b.RemoveEdge(box.Edges[0]));
    REQUIRE(edge.find("RemoveEdge: Edge#0 [") != std::string::npos);
    REQUIRE(edge.find("has uses") != std::string::npos);
    const std::string vertex = MessageOf(b.RemoveVertex(box.Vertices[0]));
    REQUIRE(vertex.find("RemoveVertex: Vertex#0 [") != std::string::npos);
    REQUIRE(vertex.find("3 edges] has edges") != std::string::npos);
    REQUIRE(b.Count<EntityKind::Edge>() == 12);
    REQUIRE(b.Count<EntityKind::Vertex>() == 8);
    REQUIRE(b.Validate().IsValid());
    // A surface not in the pool.
    REQUIRE_FALSE(b.SetSurface(box.Faces[0], Stereon::SurfaceRef::FromIndex(99)).has_value());
}

TEST_CASE("removing a shell or a solid leaves the children where they are", "[topo][edit]")
{
    Body box = MakeBox();
    ShapeBuilder& b = box.Builder;
    SECTION("the shell")
    {
        Must(b.RemoveShell(box.Shells[0]));
        REQUIRE(b.Count<EntityKind::Shell>() == 0);
        REQUIRE(b.Count<EntityKind::Face>() == 6);
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(box).Freeze();
        for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(shape))) {
            REQUIRE_FALSE(ShellOf(shape, f).IsValid());
        }
        REQUIRE(Must(Stereon::ShellsOf(shape, shape.Resolve(box.Solids[0]))).empty());
    }
    SECTION("the solid")
    {
        Must(b.RemoveSolid(box.Solids[0]));
        REQUIRE(b.Validate().IsValid());
        const Shape shape = std::move(box).Freeze();
        REQUIRE_FALSE(SolidOf(shape, shape.Resolve(box.Shells[0])).IsValid());
        REQUIRE(Collect(Must(Stereon::FacesOf(shape, shape.Resolve(box.Shells[0])))).size() == 6);
    }
}

TEST_CASE("a builder opened on a shape shares every page the edit does not touch",
          "[topo][edit][sharing]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    ShapeBuilder b = Must(ShapeBuilder::Open(a));
    REQUIRE(b.GetToken().IsValid());
    REQUIRE(b.Count<EntityKind::Face>() == 6);
    const BuilderFaceId face = b.Handle(a.Resolve(box.Faces[2]));
    REQUIRE(face.IsValid());
    REQUIRE(face.GetIndex() == 2);
    REQUIRE_FALSE(b.Handle(FaceId{}).IsValid());
    REQUIRE_FALSE(b.Handle(FaceId::FromIndex(6)).IsValid());
    const BuilderVertexId kept = b.Handle(a.Resolve(box.Vertices[5]));
    Must(b.RemoveFace(face));
    // New entities take slots after the shape's last one.
    const BuilderVertexId added = Must(b.AddVertex({9, 9, 9}, 0.0));
    REQUIRE(added.GetIndex() == 8);
    REQUIRE_FALSE(b.Handle(a.Resolve(box.Faces[2])).IsValid());
    const Shape c = Must(Stereon::Freeze(std::move(b)));

    // The source is as it was, and the two share what the edit did not touch.
    REQUIRE(a.Count<EntityKind::Face>() == 6);
    REQUIRE(c.Count<EntityKind::Face>() == 5);
    REQUIRE(c.Count<EntityKind::Vertex>() == 9);
    const Stereon::Detail::TopologyTables& ta = a.GetSnapshot()->GetTables();
    const Stereon::Detail::TopologyTables& tc = c.GetSnapshot()->GetTables();
    REQUIRE(SharedPageCount(ta.Curves, tc.Curves) == ta.Curves.GetPageCount());
    REQUIRE(SharedPageCount(ta.PCurves, tc.PCurves) == ta.PCurves.GetPageCount());
    REQUIRE(SharedPageCount(ta.Surfaces, tc.Surfaces) == ta.Surfaces.GetPageCount());
    REQUIRE(SharedPageCount(ta.Solids, tc.Solids) == ta.Solids.GetPageCount());
    REQUIRE(SharedPageCount(ta.Faces, tc.Faces) == 0);
    REQUIRE(SharedPageCount(ta.Vertices, tc.Vertices) == 0); // the added vertex
    // An index live in both names the same entity.
    REQUIRE(PositionOf(c, c.Resolve(kept)) == PositionOf(a, a.Resolve(box.Vertices[5])));
    REQUIRE(PositionOf(c, c.Resolve(added)) == Point3{9, 9, 9});
    REQUIRE(Stereon::ValidateTopology(c).IsValid());

    REQUIRE(ShapeBuilder::Open(Shape{}).error().GetCode() == ErrorCode::InvalidInput);
    if constexpr (STEREON_CHECKED_HANDLES) {
        // A handle of another shape, even to a live slot.
        const Shape other = MakeBox().Freeze();
        ShapeBuilder d = Must(ShapeBuilder::Open(a));
        REQUIRE_FALSE(d.Handle(other.Resolve(box.Faces[0])).IsValid());
    }
}

TEST_CASE("SetSurface replaces a face's surface and copies one page", "[topo][edit][sharing]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    ShapeBuilder b = Must(ShapeBuilder::Open(a));
    const Stereon::SurfaceRef moved = Must(
        b.AddSurface(Stereon::Surface(Stereon::Plane{{0, 0, 7}, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
    const BuilderFaceId face = b.Handle(a.Resolve(box.Faces[0]));
    Must(b.SetSurface(face, moved));
    const Shape c = Must(Stereon::Freeze(std::move(b)));
    const auto* plane = SurfaceOf(c, c.Resolve(face)).Target<Stereon::Plane>();
    REQUIRE(plane != nullptr);
    REQUIRE(plane->Origin == Point3{0, 0, 7});
    REQUIRE(SurfaceOf(a, a.Resolve(box.Faces[0])).Target<Stereon::Plane>()->Origin ==
            Point3{0, 0, 0});
    const Stereon::Detail::TopologyTables& ta = a.GetSnapshot()->GetTables();
    const Stereon::Detail::TopologyTables& tc = c.GetSnapshot()->GetTables();
    REQUIRE(SharedPageCount(ta.Vertices, tc.Vertices) == ta.Vertices.GetPageCount());
    REQUIRE(SharedPageCount(ta.Edges, tc.Edges) == ta.Edges.GetPageCount());
    REQUIRE(SharedPageCount(ta.Coedges, tc.Coedges) == ta.Coedges.GetPageCount());
    REQUIRE(SharedPageCount(ta.Loops, tc.Loops) == ta.Loops.GetPageCount());
    REQUIRE(SharedPageCount(ta.Faces, tc.Faces) == 0);
}

TEST_CASE("random face removals keep the structure valid, and Compact keeps the rest",
          "[topo][edit][property]")
{
    rc::prop("on prisms of 3 to 7 sides", [] {
        const std::size_t n = *rc::gen::inRange<std::size_t>(3, 8);
        const std::size_t removals = *rc::gen::inRange<std::size_t>(1, n + 3);
        Body prism = MakePrism(n);
        ShapeBuilder& b = prism.Builder;
        std::vector<BuilderFaceId> live = prism.Faces;
        for (std::size_t k = 0; k < removals && !live.empty(); ++k) {
            const std::size_t pick = *rc::gen::inRange<std::size_t>(0, live.size());
            Must(b.RemoveFace(live[pick]));
            live.erase(live.begin() + static_cast<std::ptrdiff_t>(pick));
            // A face removal can leave the shell disconnected and an Exact
            // edge with no use, which its kind forbids; nothing else.
            for (const Stereon::ValidityFinding& f : b.Validate().GetFindings()) {
                RC_ASSERT((f.Rule == TopologyRule::ShellDisconnected ||
                           f.Rule == TopologyRule::EdgeUseCount));
            }
        }
        Must(b.RemoveShell(prism.Shells[0]));
        // The report names the edges left with no use: they go, and their
        // vertices may stay isolated.
        for (const Stereon::ValidityFinding& f : b.Validate().GetFindings()) {
            RC_ASSERT(f.Rule == TopologyRule::EdgeUseCount);
            RC_ASSERT(f.Entities[0].starts_with("Edge#"));
            const auto index = static_cast<std::size_t>(std::stoul(f.Entities[0].substr(5)));
            Must(b.RemoveEdge(prism.Edges[index]));
        }
        RC_ASSERT(b.Validate().IsValid());
        const Shape shape = Must(Stereon::Freeze(std::move(b)));
        RC_ASSERT(shape.Count<EntityKind::Face>() == live.size());
        const Stereon::Compacted compacted = Must(Stereon::Compact(shape));
        RC_ASSERT(compacted.Shape.Count<EntityKind::Face>() == live.size());
        RC_ASSERT(compacted.Shape.Count<EntityKind::Edge>() == shape.Count<EntityKind::Edge>());
        RC_ASSERT(Stereon::ValidateTopology(compacted.Shape).IsValid());
        for (const FaceId f : live | std::views::transform(
                                         [&](const BuilderFaceId h) { return shape.Resolve(h); })) {
            RC_ASSERT(compacted.Translate(f).IsValid());
        }
    });
}
