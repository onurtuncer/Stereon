// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Geometric level of the validity checker (ADR-0018): ADR-0001's
// containment invariant. Every test body satisfies it; a vertex off its
// edge's end (G01) and a face surface moved away from its edges (G02) are
// found, named and refused by Commit; and a tolerance large enough to hold
// the gap makes the same geometry valid, which is what local tolerances are
// for.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::Context;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::TopologyRule;
using Stereon::ValidityLevel;
using Stereon::ValidityReport;

namespace {

ValidityReport Geometric(const Shape& shape)
{
    return Stereon::ValidateTopology(shape, ValidityLevel::Geometric);
}

/// The rules of a report's findings, in order.
std::vector<TopologyRule> RulesOf(const ValidityReport& report)
{
    std::vector<TopologyRule> rules;
    for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
        rules.push_back(finding.Rule);
    }
    return rules;
}

/// The box with face 2's surface replaced by its plane moved `offset` along
/// the plane's normal, so the face's pcurves no longer reach its edges.
ShapeBuilder BoxWithFaceMoved(const Shape& box, const Body& body, const double offset)
{
    const Stereon::FaceId face = box.Resolve(body.Faces[2]);
    Stereon::Plane plane = Deref(Stereon::SurfaceOf(box, face).Target<Stereon::Plane>());
    plane.Origin = plane.Origin + offset * Stereon::Cross(plane.XAxis, plane.YAxis);
    ShapeBuilder edit = Must(ShapeBuilder::Open(box));
    const Stereon::SurfaceRef moved = Must(edit.AddSurface(Stereon::Surface(plane)));
    Must(edit.SetSurface(edit.Handle(face), moved));
    return edit;
}

} // namespace

TEST_CASE("every test body holds its geometry within its tolerances", "[topo][validity]")
{
    std::vector<std::pair<std::string, Shape>> bodies;
    bodies.emplace_back("box", MakeBox().Freeze());
    bodies.emplace_back("prism", MakePrism(5).Freeze());
    bodies.emplace_back("cylinder", MakeCylinder().Freeze());
    bodies.emplace_back("square sheet", MakeSquareSheet().Freeze());
    bodies.emplace_back("sheet with a hole", MakeHoleSheet().Freeze());
    bodies.emplace_back("fan", MakeFan(3).Freeze());
    bodies.emplace_back("wire", MakeWire(3).Freeze());
    // Not MakeDegenerateFace: a structural fixture whose pcurves are
    // placeholders (Bodies.hpp), so its geometry is not meant to hold.
    bodies.emplace_back("hinge", MakeHinge().Freeze());
    for (const auto& [name, shape] : bodies) {
        INFO(name);
        const ValidityReport report = Geometric(shape);
        INFO(report.ToText());
        REQUIRE(report.IsValid());
    }
}

TEST_CASE("G01: a vertex off the end of its edge", "[topo][validity]")
{
    const auto wire = [](const double endTolerance) {
        Body body;
        const auto a = body.AddVertex({0, 0, 0});
        const auto b = body.AddVertex({1, 1e-3, 0}, endTolerance);
        const auto line = Must(
            body.Builder.AddCurve(Stereon::Curve(Stereon::Line{Point3{0, 0, 0}, Vec3{1, 0, 0}})));
        body.AddEdge({.Kind = Stereon::EdgeCurveKind::Spatial,
                      .Start = a,
                      .End = b,
                      .Curve = line,
                      .Range = Interval(0.0, 1.0),
                      .Tolerance = 1e-7});
        return std::move(body).Freeze();
    };
    // Structurally sound, so Freeze accepts it; the Geometric level does not.
    const ValidityReport report = Geometric(wire(1e-7));
    REQUIRE(RulesOf(report) == std::vector<TopologyRule>{TopologyRule::VertexContainment});
    REQUIRE(report.GetFindings()[0].Entities.at(0).starts_with("Vertex#1"));
    REQUIRE(report.GetFindings()[0].Message.find("end") != std::string::npos);
    // A vertex tolerance that holds the gap.
    REQUIRE(Geometric(wire(2e-3)).IsValid());
    // The Structural level alone does not look.
    REQUIRE(Stereon::ValidateTopology(wire(1e-7), ValidityLevel::Structural).IsValid());
}

TEST_CASE("G02: a face surface moved away from its edges", "[topo][validity]")
{
    Body body = MakeBox();
    const Shape box = std::move(body).Freeze();

    // Freeze certifies structure only, so the moved face publishes; the
    // Geometric level finds every use of the face's four edges.
    const Shape moved = Must(Stereon::Freeze(BoxWithFaceMoved(box, body, 1e-3)));
    const ValidityReport report = Geometric(moved);
    REQUIRE(report.GetCount() == 4);
    for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
        REQUIRE(finding.Rule == TopologyRule::EdgeContainment);
        REQUIRE(finding.Entities.size() == 2);
        REQUIRE(finding.Entities[1].starts_with("Coedge#"));
        REQUIRE(finding.Message.find("outside its tolerance") != std::string::npos);
    }

    // Commit refuses it with the rule IDs in the message.
    const auto refused = Stereon::Commit(BoxWithFaceMoved(box, body, 1e-3), Context::Default());
    REQUIRE_FALSE(refused);
    REQUIRE(refused.error().GetCode() == ErrorCode::InvalidResult);
    REQUIRE(std::string(refused.error().GetMessage()).find("G02") != std::string::npos);

    // Edge tolerances that hold the gap make the same geometry valid: the
    // local tolerances of ADR-0001.
    ShapeBuilder loosened = BoxWithFaceMoved(box, body, 1e-3);
    for (const auto e : body.Edges) {
        Must(loosened.SetTolerance(loosened.Handle(box.Resolve(e)), 2e-3, "a moved face"));
    }
    const auto accepted = Stereon::Commit(std::move(loosened), Context::Default());
    REQUIRE(accepted);
    REQUIRE(accepted->Checked == ValidityLevel::Geometric);
    REQUIRE(Geometric(accepted->Shape).IsValid());
}

TEST_CASE("Commit checks the Geometric level when asked, and always when handles are checked",
          "[topo][validity]")
{
    const Context structural = Context::Default().WithValidityLevel(ValidityLevel::Structural);
    const auto result = Stereon::Commit(std::move(MakeBox().Builder), structural);
    REQUIRE(result);
    REQUIRE(result->Checked ==
            (STEREON_CHECKED_HANDLES ? ValidityLevel::Geometric : ValidityLevel::Structural));
    const auto checked = Stereon::Commit(std::move(MakeBox().Builder), Context::Default());
    REQUIRE(checked);
    REQUIRE(checked->Checked == ValidityLevel::Geometric);
}

TEST_CASE("the Geometric rules have stable IDs", "[topo][validity]")
{
    REQUIRE(Stereon::RuleId(TopologyRule::VertexContainment) == "G01");
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeContainment) == "G02");
    REQUIRE(Stereon::RuleName(TopologyRule::VertexContainment) == "VertexContainment");
    REQUIRE(Stereon::RuleName(TopologyRule::EdgeContainment) == "EdgeContainment");
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeCache) == "S17");
}
