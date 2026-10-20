// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The maintained vertex-incidence index and the distinct-adjacency queries
// of ADR-0036, EdgesOf(vertex) and FacesOf(edge): against a referee that
// scans every record, on every test body; through deletions; and rule S16
// on seeded faults in the index.

#include "BodiesStd.hpp"

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
using Stereon::LoopId;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::TopologyRule;
using Stereon::VertexId;

namespace {

/// Every edge that starts or ends at v, in slot order: the referee for
/// EdgesOf, which must give the same without scanning.
std::vector<EdgeId> ScanEdgesOf(const Shape& shape, const VertexId v)
{
    std::vector<EdgeId> edges;
    for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
        if (StartOf(shape, e) == v || EndOf(shape, e) == v) {
            edges.push_back(e);
        }
    }
    return edges;
}

/// Every face with a coedge of e, in slot order, once: the referee for
/// FacesOf.
std::vector<FaceId> ScanFacesOf(const Shape& shape, const EdgeId e)
{
    std::vector<FaceId> faces;
    for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(shape))) {
        bool uses = false;
        for (const LoopId l : Must(Stereon::LoopsOf(shape, f))) {
            for (const CoedgeId c : Must(Stereon::CoedgesOf(shape, l))) {
                uses = uses || EdgeOf(shape, c) == e;
            }
        }
        if (uses) {
            faces.push_back(f);
        }
    }
    return faces;
}

void CheckAgainstReferee(const Shape& shape)
{
    for (const VertexId v : Must(Stereon::AllEntities<EntityKind::Vertex>(shape))) {
        REQUIRE(Collect(Must(Stereon::EdgesOf(shape, v))) == ScanEdgesOf(shape, v));
    }
    for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(shape))) {
        REQUIRE(Must(Stereon::FacesOf(shape, e)) == ScanFacesOf(shape, e));
    }
}

} // namespace

TEST_CASE("EdgesOf and FacesOf agree with a scan of every record on every body",
          "[topo][adjacency]")
{
    CheckAgainstReferee(MakeBox().Freeze());
    CheckAgainstReferee(MakePrism(5).Freeze());
    CheckAgainstReferee(MakeCylinder().Freeze());
    CheckAgainstReferee(MakeSquareSheet().Freeze());
    CheckAgainstReferee(MakeHoleSheet().Freeze());
    CheckAgainstReferee(MakeFan(4).Freeze());
    CheckAgainstReferee(MakeWire(3).Freeze());
    CheckAgainstReferee(MakeDegenerateFace().Freeze());
}

TEST_CASE("a box corner has three edges and a cylinder vertex its circle once and the seam",
          "[topo][adjacency]")
{
    SECTION("box")
    {
        const Shape box = MakeBox().Freeze();
        for (const VertexId v : Must(Stereon::AllEntities<EntityKind::Vertex>(box))) {
            const std::vector<EdgeId> edges = Collect(Must(Stereon::EdgesOf(box, v)));
            REQUIRE(edges.size() == 3);
            REQUIRE(edges[0] < edges[1]);
            REQUIRE(edges[1] < edges[2]);
        }
        for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(box))) {
            REQUIRE(Must(Stereon::FacesOf(box, e)).size() == 2);
        }
    }
    SECTION("cylinder")
    {
        Body body = MakeCylinder();
        const Shape shape = std::move(body).Freeze();
        const EdgeId bottom = shape.Resolve(body.Edges[0]);
        const EdgeId top = shape.Resolve(body.Edges[1]);
        const EdgeId seam = shape.Resolve(body.Edges[2]);
        // The closed bottom circle starts and ends at vertex 0: listed once.
        REQUIRE(Collect(Must(Stereon::EdgesOf(shape, shape.Resolve(body.Vertices[0])))) ==
                std::vector<EdgeId>{bottom, seam});
        REQUIRE(Collect(Must(Stereon::EdgesOf(shape, shape.Resolve(body.Vertices[1])))) ==
                std::vector<EdgeId>{top, seam});
        // The seam's two uses lie on one face.
        REQUIRE(Collect(Must(Stereon::RadialCoedgesOf(shape, seam))).size() == 2);
        REQUIRE(Must(Stereon::FacesOf(shape, seam)) ==
                std::vector<FaceId>{shape.Resolve(body.Faces[0])});
        REQUIRE(Must(Stereon::FacesOf(shape, bottom)) ==
                std::vector<FaceId>{shape.Resolve(body.Faces[0]), shape.Resolve(body.Faces[1])});
    }
    SECTION("fan and wire")
    {
        Body fan = MakeFan(5);
        const Shape shape = std::move(fan).Freeze();
        REQUIRE(Must(Stereon::FacesOf(shape, shape.Resolve(fan.Edges[0]))).size() == 5);
        REQUIRE(Collect(Must(Stereon::EdgesOf(shape, shape.Resolve(fan.Vertices[0])))).size() == 6);
        Body polyline = MakeWire(3);
        const Shape wire = std::move(polyline).Freeze();
        for (const EdgeId e : Must(Stereon::AllEntities<EntityKind::Edge>(wire))) {
            REQUIRE(Must(Stereon::FacesOf(wire, e)).empty());
        }
        REQUIRE(Stereon::Describe(wire, wire.Resolve(polyline.Vertices[1])).find("2 edges") !=
                std::string::npos);
    }
}

TEST_CASE("the index follows deletions", "[topo][adjacency][edit]")
{
    Body fan = MakeFan(4);
    const Shape before = std::move(fan).Freeze();
    ShapeBuilder b = Must(ShapeBuilder::Open(before));
    // The edited shape resolves the opened builder's handles, not the
    // original builder's.
    const BuilderEdgeId spine = b.Handle(before.Resolve(fan.Edges[0]));
    const BuilderVertexId origin = b.Handle(before.Resolve(fan.Vertices[0]));
    const BuilderVertexId tip = b.Handle(before.Resolve(fan.Vertices[3]));
    // Face 1 goes with its two Boundary side edges, which then have no use.
    Must(b.RemoveFace(b.Handle(before.Resolve(fan.Faces[1]))));
    Must(b.RemoveEdge(b.Handle(before.Resolve(fan.Edges[3]))));
    Must(b.RemoveEdge(b.Handle(before.Resolve(fan.Edges[4]))));
    const Shape after = Must(Stereon::Freeze(std::move(b)));
    CheckAgainstReferee(after);
    if constexpr (STEREON_CHECKED_HANDLES) {
        // Only the store tells the original builder's handle apart: the
        // spine's slot is live at the same generation in the edited shape.
        REQUIRE_FALSE(after.Resolve(fan.Edges[0]).IsValid());
    }
    REQUIRE(Must(Stereon::FacesOf(after, after.Resolve(spine))).size() == 3);
    REQUIRE(Collect(Must(Stereon::EdgesOf(after, after.Resolve(origin)))).size() == 4);
    // The tip of the removed face is now an isolated vertex.
    REQUIRE(Must(Stereon::EdgesOf(after, after.Resolve(tip))).empty());
    // The shape it was opened on is as it was.
    REQUIRE(Must(Stereon::FacesOf(before, before.Resolve(fan.Edges[0]))).size() == 4);
}

TEST_CASE("EdgesOf and FacesOf refuse the empty shape and handles not in the shape",
          "[topo][adjacency]")
{
    const Shape empty;
    REQUIRE(Stereon::EdgesOf(empty, VertexId::FromIndex(0)).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(Stereon::FacesOf(empty, EdgeId::FromIndex(0)).error().GetCode() ==
            ErrorCode::InvalidInput);
    const Shape box = MakeBox().Freeze();
    REQUIRE_FALSE(Stereon::EdgesOf(box, VertexId::FromIndex(8)).has_value());
    REQUIRE_FALSE(Stereon::FacesOf(box, EdgeId::FromIndex(12)).has_value());
    REQUIRE_FALSE(Stereon::EdgesOf(box, VertexId{}).has_value());
}

TEST_CASE("rule S16 catches a corrupted vertex-incidence index", "[topo][adjacency][validator]")
{
    using Stereon::Testing::TopologyFaults;
    Body wire = MakeWire(3); // v0 -e0- v1 -e1- v2 -e2- v3
    ShapeBuilder& b = wire.Builder;
    REQUIRE(b.Validate().IsValid());
    SECTION("an edge not at the vertex")
    {
        TopologyFaults::SetFirstEdge(b, wire.Vertices[0], wire.Edges[2]);
        const Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::VertexIncidence));
        REQUIRE(report.GetFindings()[0].Entities[0].starts_with("Vertex#0"));
        REQUIRE(report.GetFindings()[0].Message.find("does not start or end") != std::string::npos);
    }
    SECTION("a list that cycles")
    {
        TopologyFaults::SetNextAtStart(b, wire.Edges[0], wire.Edges[0]);
        const Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::VertexIncidence));
        REQUIRE(report.GetFindings()[0].Message.find("ascending") != std::string::npos);
    }
    SECTION("an edge missing from its vertex's list")
    {
        TopologyFaults::SetFirstEdge(b, wire.Vertices[0], BuilderEdgeId{});
        const Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::VertexIncidence));
        REQUIRE(report.GetFindings()[0].Entities[0].starts_with("Edge#0"));
        REQUIRE(report.GetFindings()[0].Message.find("start vertex's list") != std::string::npos);
    }
    SECTION("a link past the end of the table")
    {
        TopologyFaults::SetNextAtEnd(b, wire.Edges[0], BuilderEdgeId::Make(100, 0, 0));
        const Stereon::ValidityReport report = b.Validate();
        REQUIRE(report.Has(TopologyRule::DanglingLink));
        REQUIRE(report.GetFindings()[0].Message.find("NextAtEnd") != std::string::npos);
    }
    REQUIRE_FALSE(Stereon::Freeze(std::move(b)).has_value());
}
