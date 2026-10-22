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

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"
#include "FaultGeometry.hpp"

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
    // MakePoleCap is the pole whose geometry does.
    bodies.emplace_back("pole cap", MakePoleCap().Freeze());
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
                      .Tolerance = Tight});
        return std::move(body).Freeze();
    };
    // Structurally sound, so Freeze accepts it; the Geometric level does not.
    const ValidityReport report = Geometric(wire(Tight));
    REQUIRE(RulesOf(report) == std::vector<TopologyRule>{TopologyRule::VertexContainment});
    REQUIRE(report.GetFindings()[0].Entities.at(0).starts_with("Vertex#1"));
    REQUIRE(report.GetFindings()[0].Message.find("end") != std::string::npos);
    // A vertex tolerance that holds the gap.
    REQUIRE(Geometric(wire(2e-3)).IsValid());
    // The Structural level alone does not look.
    REQUIRE(Stereon::ValidateTopology(wire(Tight), ValidityLevel::Structural).IsValid());
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

TEST_CASE("an empty Shape is valid at every level", "[topo][validity]")
{
    const Shape empty;
    Shape moved = MakeBox().Freeze();
    const Shape taken(std::move(moved));
    for (const Shape* shape : {&empty, static_cast<const Shape*>(&moved)}) {
        REQUIRE(shape->IsEmpty()); // NOLINT(bugprone-use-after-move)
        for (const ValidityLevel level :
             {ValidityLevel::Structural, ValidityLevel::Geometric, ValidityLevel::Solid}) {
            REQUIRE(Stereon::ValidateTopology(*shape, level).IsValid());
        }
    }
    REQUIRE_FALSE(taken.IsEmpty());
    REQUIRE(Geometric(taken).IsValid());
}

TEST_CASE("an edge's range is bounded", "[topo][validity]")
{
    const Interval infinite = Interval::Everything();
    const double inf = infinite.GetUpper();
    for (const Interval range : {infinite, Interval(0.0, inf), Interval(-inf, 1.0)}) {
        INFO(range.GetLower() << ", " << range.GetUpper());
        Body body;
        const auto a = body.AddVertex({0, 0, 0});
        const auto b = body.AddVertex({1, 0, 0});
        const auto line = Must(
            body.Builder.AddCurve(Stereon::Curve(Stereon::Line{Point3{0, 0, 0}, Vec3{1, 0, 0}})));
        const auto refused = body.Builder.AddEdge({.Kind = Stereon::EdgeCurveKind::Spatial,
                                                   .Start = a,
                                                   .End = b,
                                                   .Curve = line,
                                                   .Range = range,
                                                   .Tolerance = Tight});
        REQUIRE_FALSE(refused);
        REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
        REQUIRE(std::string(refused.error().GetMessage()).find("not bounded") != std::string::npos);
    }
    // Bounded ranges inside and outside the curve's domain: the first holds,
    // the second is a finding, never skipped.
    const auto arc = [](const Interval range) {
        Body body;
        const Stereon::Circle circle{{0, 0, 0}, 1.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
        const auto a = body.AddVertex(Eval(circle, range.GetLower()));
        const auto b = body.AddVertex(Eval(circle, range.GetUpper()));
        const auto curve = Must(body.Builder.AddCurve(Stereon::Curve(
            *Stereon::MakeTrimmedCurve(Stereon::Curve(circle), Interval(0.0, 1.0)))));
        body.AddEdge({.Kind = Stereon::EdgeCurveKind::Spatial,
                      .Start = a,
                      .End = b,
                      .Curve = curve,
                      .Range = range,
                      .Tolerance = Tight});
        return std::move(body).Freeze();
    };
    REQUIRE(Geometric(arc(Interval(0.25, 0.75))).IsValid());
    const ValidityReport outside = Geometric(arc(Interval(0.5, 2.0)));
    REQUIRE(RulesOf(outside) == std::vector<TopologyRule>{TopologyRule::EdgeContainment});
    REQUIRE(outside.GetFindings()[0].Message.find("does not cover") != std::string::npos);
}

TEST_CASE("G02: a degenerate edge's trace stays at its vertex", "[topo][validity]")
{
    // A pole of a sphere passes.
    REQUIRE(Geometric(MakePoleCap().Freeze()).IsValid());

    // A circle closed at its vertex, alleged a pole: its ends coincide with
    // the vertex, but its middle is 2 mm away. `uses` faces use it, on the
    // same plane with the same circle, so the uses agree with each other.
    const auto circle = [](const std::size_t uses) {
        Body body;
        const auto v = body.AddVertex({1, 0, 0});
        const auto pole = body.AddEdge({.Kind = Stereon::EdgeCurveKind::Degenerate,
                                        .Start = v,
                                        .End = v,
                                        .Range = Interval(0.0, Stereon::CirclePeriod),
                                        .Tolerance = Tight});
        for (std::size_t k = 0; k < uses; ++k) {
            const auto loop = body.AddLoop(body.AddFace(XyPlane()));
            body.AddCoedge(
                loop, pole, Stereon::EdgeUseSense::Same,
                Stereon::Curve2d(Stereon::Circle2d{{0, 0}, 1.0, Vec2{1, 0}, Vec2{0, 1}}));
        }
        return std::move(body).Freeze();
    };
    for (const std::size_t uses : {std::size_t{1}, std::size_t{2}}) {
        INFO(uses << " uses");
        const ValidityReport report = Geometric(circle(uses));
        REQUIRE(report.GetCount() == uses);
        for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
            REQUIRE(finding.Rule == TopologyRule::EdgeContainment);
            REQUIRE(finding.Message.find("degenerate edge's trace") != std::string::npos);
        }
    }
    // Commit refuses it.
    Body body;
    const auto v = body.AddVertex({1, 0, 0});
    const auto pole = body.AddEdge({.Kind = Stereon::EdgeCurveKind::Degenerate,
                                    .Start = v,
                                    .End = v,
                                    .Range = Interval(0.0, Stereon::CirclePeriod),
                                    .Tolerance = Tight});
    body.AddCoedge(body.AddLoop(body.AddFace(XyPlane())), pole, Stereon::EdgeUseSense::Same,
                   Stereon::Curve2d(Stereon::Circle2d{{0, 0}, 1.0, Vec2{1, 0}, Vec2{0, 1}}));
    REQUIRE(Stereon::Commit(std::move(body.Builder), Context::Default()).error().GetCode() ==
            ErrorCode::InvalidResult);
}

TEST_CASE("the Geometric rules have stable IDs", "[topo][validity]")
{
    REQUIRE(Stereon::RuleId(TopologyRule::VertexContainment) == "G01");
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeContainment) == "G02");
    REQUIRE(Stereon::RuleName(TopologyRule::VertexContainment) == "VertexContainment");
    REQUIRE(Stereon::RuleName(TopologyRule::EdgeContainment) == "EdgeContainment");
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeCache) == "S17");
}

TEST_CASE("G02 rejects every non-finite interior sample on a lone Spatial trace",
          "[topo][validity][property]")
{
    for (const double bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        for (int k = 1; k < 16; ++k) {
            Body body;
            const auto a = body.AddVertex({0, 0, 0});
            const auto z = body.AddVertex({1, 0, 0});
            FaultGeometry::Curve faulty;
            faulty.At = static_cast<double>(k) / 16;
            faulty.Value = bad;
            const auto curve = Must(body.Builder.AddCurve(Stereon::Curve(faulty)));
            body.AddEdge({.Kind = Stereon::EdgeCurveKind::Spatial,
                          .Start = a,
                          .End = z,
                          .Curve = curve,
                          .Range = Interval{0.0, 1.0},
                          .Tolerance = Tight});
            const auto shape = std::move(body).Freeze();
            const auto report = Geometric(shape);
            REQUIRE_FALSE(report.IsValid());
            REQUIRE(report.GetFindings().front().Rule == TopologyRule::EdgeContainment);
            auto edit = Must(ShapeBuilder::Open(shape));
            const auto result = Stereon::Commit(std::move(edit), Context::Default());
            REQUIRE_FALSE(result);
            REQUIRE(result.error().GetCode() == ErrorCode::InvalidResult);
        }
    }
}

TEST_CASE("G02 retains a non-finite surface sample before later finite samples", "[topo][validity]")
{
    for (const bool pole : {false, true}) {
        Body body = pole ? MakePoleCap() : MakeBox(1.0, 1.0, 1.0);
        const auto original = std::move(body).Freeze();
        const auto face = original.Resolve(body.Faces.front());
        FaultGeometry::Surface faulty;
        faulty.Base = Stereon::SurfaceOf(original, face);
        faulty.At = pole ? std::numbers::pi : 0.5;
        faulty.Value = std::numeric_limits<double>::quiet_NaN();
        auto edit = Must(ShapeBuilder::Open(original));
        const auto replacement = Must(edit.AddSurface(Stereon::Surface(faulty)));
        Must(edit.SetSurface(edit.Handle(face), replacement));
        const auto shape = Must(Stereon::Freeze(std::move(edit)));
        REQUIRE_FALSE(Geometric(shape).IsValid());
        REQUIRE(Geometric(original).IsValid());
    }
}
