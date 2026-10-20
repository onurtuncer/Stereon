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
#include <limits>
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
    bodies.push_back(MakeHinge());
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

TEST_CASE("S12 MissingGeometry: an empty curve, pcurve or surface in the pool",
          "[topo][validator][S12][audit]")
{
    // Audit 2026-10-20, F3: an index inside the pool is not a support unless
    // the entry holds geometry. The findings name the entities with
    // Describe, which must not visit the empty wrapper. In the bodies here
    // every face, edge and coedge added its own geometry, in order, so its
    // index in the pool is its creation index.
    const auto requireEmptyFound = [](Body& body, const std::string& named) {
        const ValidityReport report = body.Builder.Validate();
        const std::string text = report.ToText();
        INFO(text);
        REQUIRE(report.GetCount() == 1);
        REQUIRE(text.find("in the pool is empty") != std::string::npos);
        REQUIRE(text.find(named) != std::string::npos);
        RequireFinds(body, TopologyRule::MissingGeometry);
    };
    SECTION("an empty surface under a face")
    {
        Body body = MakeBox();
        Faults::EmptySurface(body.Builder, Stereon::SurfaceRef::FromIndex(2));
        requireEmptyFound(body, "Face#2 [empty geometry");
    }
    SECTION("an empty curve under an Exact edge")
    {
        Body body = MakeBox();
        Faults::EmptyCurve(body.Builder, Stereon::CurveRef::FromIndex(4));
        requireEmptyFound(body, "Edge#4 [Exact empty geometry");
    }
    SECTION("an empty curve under a Spatial edge")
    {
        Body body = MakeWire(3);
        Faults::EmptyCurve(body.Builder, Stereon::CurveRef::FromIndex(1));
        requireEmptyFound(body, "Edge#1 [Spatial empty geometry");
    }
    SECTION("an empty pcurve under a coedge")
    {
        Body body = MakeCylinder();
        Faults::EmptyPCurve(body.Builder, PCurveRef::FromIndex(5));
        requireEmptyFound(body, "Coedge#5");
    }
}

TEST_CASE("S17 EdgeCache: a cached curve and its certificate that disagree",
          "[topo][validator][S17][cache]")
{
    constexpr double Infinite = std::numeric_limits<double>::infinity();
    using Stereon::CacheCertificate;
    using Stereon::CurveRef;
    SECTION("a cache on an Exact edge, whose curve is its own")
    {
        Body body = MakeBox();
        Faults::SetCache(body.Builder, body.Edges[0], CurveRef::FromIndex(0),
                         CacheCertificate::Uncertified, Infinite);
        RequireFinds(body, TopologyRule::EdgeCache);
    }
    SECTION("a cache past the curve pool")
    {
        Body body = MakeHinge();
        Faults::SetCache(body.Builder, body.Edges[0], CurveRef::FromIndex(4000),
                         CacheCertificate::Uncertified, Infinite);
        RequireFinds(body, TopologyRule::EdgeCache);
    }
    SECTION("a status that claims a cache the edge does not have")
    {
        Body body = MakeHinge();
        Faults::SetCache(body.Builder, body.Edges[0], CurveRef{}, CacheCertificate::Uncertified,
                         Infinite);
        RequireFinds(body, TopologyRule::EdgeCache);
    }
    SECTION("an uncertified cache with a bound, a certified one without")
    {
        for (const auto& [status, bound] : {std::pair{CacheCertificate::Uncertified, 1e-6},
                                            std::pair{CacheCertificate::Certified, Infinite},
                                            std::pair{CacheCertificate::Certified, -1.0}}) {
            Body body = MakeHinge();
            const CurveRef cache = Must(body.Builder.CacheEdgeCurve(body.Edges[0]));
            Faults::SetCache(body.Builder, body.Edges[0], cache, status, bound);
            RequireFinds(body, TopologyRule::EdgeCache);
        }
    }
    SECTION("a cache as CacheEdgeCurve leaves it is valid")
    {
        Body body = MakeHinge();
        static_cast<void>(Must(body.Builder.CacheEdgeCurve(body.Edges[0])));
        REQUIRE(body.Builder.Validate().IsValid());
    }
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeCache) == "S17");
    REQUIRE(Stereon::RuleName(TopologyRule::EdgeCache) == "EdgeCache");
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

namespace {

/// The fan of k faces around one spine, its faces moved into two open
/// shells: face i into the second when bit i of `mask` is set. Face 0 was
/// added first, so its use of the spine is the edge's radial anchor.
Body FanInTwoShells(const std::size_t k, const unsigned mask)
{
    Body body = MakeFan(k);
    Must(body.Builder.RemoveShell(body.Shells[0]));
    const BuilderShellId first = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    const BuilderShellId second = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    for (std::size_t i = 0; i < k; ++i) {
        Must(body.Builder.AddFaceToShell(((mask >> i) & 1U) != 0 ? second : first, body.Faces[i]));
    }
    return body;
}

} // namespace

TEST_CASE("S15: an edge shared by two shells connects the faces of each, whichever holds its "
          "anchor",
          "[topo][validator][S15][audit]")
{
    // Audit 2026-10-20, F5: the faces around an edge were joined only when
    // they lay in the anchor face's shell, so a shell sharing a
    // non-manifold edge with another was reported as disconnected, or not,
    // depending on which shell's face came first. Every split of the fan's
    // faces between two shells keeps each shell connected through the
    // spine.
    for (const std::size_t k : {std::size_t{3}, std::size_t{4}}) {
        for (unsigned mask = 1; mask + 1 < (1U << k); ++mask) {
            Body body = FanInTwoShells(k, mask);
            const ValidityReport report = body.Builder.Validate();
            INFO("k = " << k << ", mask = " << mask << "\n" << report.ToText());
            REQUIRE_FALSE(report.Has(TopologyRule::ShellDisconnected));
            REQUIRE(report.IsValid());
            REQUIRE(Stereon::Freeze(std::move(body.Builder)).has_value());
        }
    }

    // A shell in two pieces is still found when the edge it shares lies in
    // another shell's anchor face: faces 1 and 2 share the spine with face
    // 0, but the square shares nothing with them.
    Body body = FanInTwoShells(3, 0b110U);
    const BuilderFaceId apart = body.AddFace(XyPlane());
    AddBoundaryPolygon(body, apart, std::vector<Point3>{{5, 5, 0}, {6, 5, 0}, {5, 6, 0}});
    Must(body.Builder.AddFaceToShell(body.Shells[2], apart));
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
