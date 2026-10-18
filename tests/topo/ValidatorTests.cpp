// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The structural validator (ADR-0018, ADR-0038): a seeded fault for every
// rule, found by that rule; no finding on any valid body; termination and
// bounded reports on stores corrupted beyond repair.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::EntityKind;
using Stereon::PCurveRef;
using Stereon::TopologyRule;
using Stereon::ValidityReport;
using Faults = Stereon::Testing::TopologyFaults;

namespace {

/// The builder's report has a finding of `rule`, and Freeze() refuses it,
/// naming the rule.
void RequireFinds(Body& body, const TopologyRule rule)
{
    const ValidityReport report = body.Builder.Validate();
    INFO(report.ToText());
    REQUIRE(report.Has(rule));
    const auto frozen = Stereon::Freeze(std::move(body.Builder));
    REQUIRE_FALSE(frozen.has_value());
    REQUIRE(frozen.error().GetCode() == Stereon::ErrorCode::InvalidInput);
    const std::string text(frozen.error().GetMessage());
    REQUIRE(text.find(Stereon::RuleId(rule) + " " + std::string(Stereon::RuleName(rule))) !=
            std::string::npos);
}

} // namespace

TEST_CASE("rule IDs and names are stable", "[topo][validator]")
{
    REQUIRE(Stereon::RuleId(TopologyRule::DanglingLink) == "S01");
    REQUIRE(Stereon::RuleId(TopologyRule::EmptyLoop) == "S14");
    REQUIRE(Stereon::RuleId(TopologyRule::ShellDisconnected) == "S15");
    REQUIRE(Stereon::RuleName(TopologyRule::RadialRingOpen) == "RadialRingOpen");
    REQUIRE(Stereon::RuleName(TopologyRule::PrimaryCoedge) == "PrimaryCoedge");
}

TEST_CASE("no finding on any valid body", "[topo][validator]")
{
    std::vector<Body> bodies;
    bodies.push_back(MakeBox());
    bodies.push_back(MakeCylinder());
    bodies.push_back(MakePrism(3));
    bodies.push_back(MakePrism(17));
    bodies.push_back(MakeSquareSheet());
    bodies.push_back(MakeHoleSheet());
    bodies.push_back(MakeFan(1));
    bodies.push_back(MakeFan(5));
    bodies.push_back(MakeWire(4));
    bodies.push_back(MakeDegenerateFace());
    bodies.emplace_back();
    for (Body& body : bodies) {
        const ValidityReport report = body.Builder.Validate();
        INFO(report.ToText());
        REQUIRE(report.IsValid());
        REQUIRE(report.GetFindings().empty());
        const Stereon::Shape shape = std::move(body).Freeze();
        REQUIRE(Stereon::ValidateTopology(shape).IsValid());
    }
}

TEST_CASE("S01 DanglingLink: a link to a slot out of range or dead", "[topo][validator][S01]")
{
    SECTION("out of range")
    {
        Body body = MakeBox();
        Faults::SetNextOutOfRange(body.Builder, body.Coedges[5]);
        RequireFinds(body, TopologyRule::DanglingLink);
    }
    SECTION("a dead coedge still linked from its neighbours")
    {
        Body body = MakeBox();
        Faults::KillCoedge(body.Builder, body.Coedges[0]);
        RequireFinds(body, TopologyRule::DanglingLink);
    }
    SECTION("a primary that names nothing")
    {
        Body body = MakeSquareSheet();
        Faults::SetPrimary(body.Builder, body.Edges[0], Stereon::BuilderCoedgeId::Make(77, 0, 0));
        RequireFinds(body, TopologyRule::DanglingLink);
    }
}

TEST_CASE("S02 LoopLinksBroken: Next and Prev disagree", "[topo][validator][S02]")
{
    Body body = MakeBox();
    // Coedges 0..3 are the first face's loop.
    Faults::SetPrev(body.Builder, body.Coedges[1], body.Coedges[2]);
    RequireFinds(body, TopologyRule::LoopLinksBroken);
}

TEST_CASE("S03 LoopOpen: a walk that leaves its loop or circles short of the anchor",
          "[topo][validator][S03]")
{
    SECTION("into another loop")
    {
        Body body = MakeBox();
        Faults::SetNext(body.Builder, body.Coedges[1], body.Coedges[6]);
        RequireFinds(body, TopologyRule::LoopOpen);
    }
    SECTION("a cycle that does not pass the anchor")
    {
        Body body = MakeBox();
        Faults::SetNext(body.Builder, body.Coedges[2], body.Coedges[1]);
        RequireFinds(body, TopologyRule::LoopOpen);
    }
}

TEST_CASE("S04 CoedgeOutsideLoop: a coedge its loop never reaches", "[topo][validator][S04]")
{
    SECTION("a coedge that claims another loop")
    {
        Body body = MakeBox();
        // Coedge 0 claims the second face's loop, which never walks through it.
        Faults::SetCoedgeLoop(body.Builder, body.Coedges[0], body.Loops[1]);
        RequireFinds(body, TopologyRule::CoedgeOutsideLoop);
    }
    SECTION("an anchor that names itself next: the walk closes after one step")
    {
        Body body = MakePrism(5);
        Faults::SetNext(body.Builder, body.Coedges[0], body.Coedges[0]);
        const ValidityReport report = body.Builder.Validate();
        REQUIRE(report.Has(TopologyRule::LoopLinksBroken));
        REQUIRE_FALSE(report.Has(TopologyRule::LoopOpen));
        RequireFinds(body, TopologyRule::CoedgeOutsideLoop);
    }
}

TEST_CASE("S05 RadialRingOpen: a ring that reaches a use of another edge", "[topo][validator][S05]")
{
    Body body = MakeBox();
    Faults::SetRadial(body.Builder, body.Coedges[0], body.Coedges[1]);
    RequireFinds(body, TopologyRule::RadialRingOpen);
}

TEST_CASE("S06 CoedgeOutsideRing: a use the edge's ring skips", "[topo][validator][S06]")
{
    // The fan's spine has three uses: coedges 0, 3 and 6. Short-circuit the
    // ring at its anchor so that the other two drop out of it.
    Body body = MakeFan(3);
    Faults::SetRadial(body.Builder, body.Coedges[0], body.Coedges[0]);
    RequireFinds(body, TopologyRule::CoedgeOutsideRing);
}

TEST_CASE("S07 EdgeUseCount: coedge counts that do not fit the kind", "[topo][validator][S07]")
{
    SECTION("an Intersection edge with one use")
    {
        Body body = MakeSquareSheet();
        Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Intersection);
        RequireFinds(body, TopologyRule::EdgeUseCount);
    }
    SECTION("a Trim edge with three uses")
    {
        Body body = MakeFan(3);
        Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Trim);
        RequireFinds(body, TopologyRule::EdgeUseCount);
    }
    SECTION("an Exact edge with no use")
    {
        Body body = MakeWire(1);
        Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Exact);
        RequireFinds(body, TopologyRule::EdgeUseCount);
    }
    SECTION("a two-use Boundary edge on two faces is not a seam")
    {
        Body body = MakeBox();
        Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Boundary);
        RequireFinds(body, TopologyRule::EdgeUseCount);
    }
}

TEST_CASE("S08 PrimaryCoedge: missing, or not a use of the edge", "[topo][validator][S08]")
{
    SECTION("missing on a surface-defined edge")
    {
        Body body = MakeSquareSheet();
        Faults::SetPrimary(body.Builder, body.Edges[1], Stereon::BuilderCoedgeId{});
        RequireFinds(body, TopologyRule::PrimaryCoedge);
    }
    SECTION("a use of another edge")
    {
        Body body = MakeCylinder();
        Faults::SetPrimary(body.Builder, body.Edges[2], body.Coedges[0]);
        RequireFinds(body, TopologyRule::PrimaryCoedge);
    }
}

TEST_CASE("S09 LoopOutsideFace: a loop its face does not list", "[topo][validator][S09]")
{
    Body body = MakeBox();
    Faults::SetLoopFace(body.Builder, body.Loops[0], body.Faces[3]);
    RequireFinds(body, TopologyRule::LoopOutsideFace);
}

TEST_CASE("S10 FaceOutsideShell: a face its shell does not list", "[topo][validator][S10]")
{
    Body body = MakeBox();
    const BuilderShellId other = Must(body.Builder.AddShell());
    Faults::SetFaceShell(body.Builder, body.Faces[2], other);
    RequireFinds(body, TopologyRule::FaceOutsideShell);
}

TEST_CASE("S11 ShellOutsideSolid: a shell its solid does not list", "[topo][validator][S11]")
{
    Body body = MakeBox();
    const BuilderSolidId other = Must(body.Builder.AddSolid());
    Faults::SetShellSolid(body.Builder, body.Shells[0], other);
    RequireFinds(body, TopologyRule::ShellOutsideSolid);
}

TEST_CASE("S12 MissingGeometry: a use without its pcurve", "[topo][validator][S12]")
{
    SECTION("no pcurve")
    {
        Body body = MakeBox();
        Faults::SetPCurve(body.Builder, body.Coedges[3], PCurveRef{});
        RequireFinds(body, TopologyRule::MissingGeometry);
    }
    SECTION("a pcurve past the pool")
    {
        Body body = MakeBox();
        Faults::SetPCurve(body.Builder, body.Coedges[3], PCurveRef::FromIndex(4000));
        RequireFinds(body, TopologyRule::MissingGeometry);
    }
    SECTION("a Spatial edge's use may go without")
    {
        Body body = MakeBox();
        Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Spatial);
        for (const Stereon::BuilderCoedgeId c : body.Coedges) {
            Faults::SetPCurve(body.Builder, c, PCurveRef{});
        }
        const ValidityReport report = body.Builder.Validate();
        // Every other edge is Exact and still needs its pcurves.
        REQUIRE(report.Has(TopologyRule::MissingGeometry));
        REQUIRE(report.GetCount() == 24 - 2);
    }
}

TEST_CASE("S13 EdgeEndpoints: a degenerate edge with two vertices", "[topo][validator][S13]")
{
    Body body = MakeDegenerateFace();
    Faults::SetEnd(body.Builder, body.Edges[2], body.Vertices[1]);
    RequireFinds(body, TopologyRule::EdgeEndpoints);
}

TEST_CASE("S14 EmptyLoop: a loop without coedges", "[topo][validator][S14]")
{
    Body body = MakeSquareSheet();
    body.AddLoop(body.Faces[0]);
    RequireFinds(body, TopologyRule::EmptyLoop);
}

TEST_CASE("S15 ShellDisconnected: faces that share no edge chain", "[topo][validator][S15]")
{
    Body body;
    const BuilderShellId shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    for (int i = 0; i < 2; ++i) {
        const double x = 3.0 * i;
        const BuilderFaceId face = body.AddFace(XyPlane());
        AddBoundaryPolygon(body, face, std::vector<Point3>{{x, 0, 0}, {x + 1, 0, 0}, {x, 1, 0}});
        Must(body.Builder.AddFaceToShell(shell, face));
    }
    RequireFinds(body, TopologyRule::ShellDisconnected);
}

TEST_CASE("the validator terminates on a store corrupted throughout", "[topo][validator]")
{
    Body body = MakePrism(40);
    const std::size_t n = body.Coedges.size();
    for (std::size_t i = 0; i < n; ++i) {
        Faults::SetNext(body.Builder, body.Coedges[i], body.Coedges[(i * 7 + 3) % n]);
        Faults::SetRadial(body.Builder, body.Coedges[i], body.Coedges[(i * 13 + 1) % n]);
        if (i % 3 == 0) {
            Faults::SetCoedgeLoop(body.Builder, body.Coedges[i],
                                  body.Loops[(i * 5) % body.Loops.size()]);
        }
        if (i % 11 == 0) {
            Faults::KillCoedge(body.Builder, body.Coedges[i]);
        }
    }
    const ValidityReport report = body.Builder.Validate();
    REQUIRE_FALSE(report.IsValid());
    // The report keeps a bounded number of findings and counts the rest.
    REQUIRE(report.GetFindings().size() == ValidityReport::MaxFindings);
    REQUIRE(report.GetCount() > ValidityReport::MaxFindings);
    const std::string text = report.ToText();
    REQUIRE(text.find("more\n") != std::string::npos);
    REQUIRE_FALSE(Stereon::Freeze(std::move(body.Builder)).has_value());
}

TEST_CASE("findings name the entities with Describe", "[topo][validator]")
{
    Body body = MakeBox();
    Faults::SetPrev(body.Builder, body.Coedges[1], body.Coedges[2]);
    const ValidityReport report = body.Builder.Validate();
    REQUIRE_FALSE(report.GetFindings().empty());
    for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
        REQUIRE_FALSE(finding.Entities.empty());
        REQUIRE_FALSE(finding.Message.empty());
        for (const std::string& entity : finding.Entities) {
            REQUIRE(entity.find('#') != std::string::npos);
            REQUIRE(entity.find('[') != std::string::npos);
        }
    }
}
