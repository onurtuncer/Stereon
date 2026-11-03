// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0018's verification: a seeded-fault corpus in which each fault is found
// by its intended rule at its intended level and is invisible to the levels
// below it, and a valid corpus with no finding at any implemented level. The
// Geometric rules G03 to G06 and the Solid level's L01 are tested here too:
// a jump of a period across a seam, a tolerance outside the Context's
// bounds, a Trim or Intersection edge on the wrong surfaces, a cached curve
// out of its tube or past its certified bound, and a solid that the Solid
// level cannot yet check, which is refused rather than passed.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::Context;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::TopologyRule;
using Stereon::ValidityLevel;
using Stereon::ValidityReport;
using Faults = Stereon::Testing::TopologyFaults;

namespace {

const Context ctx = Context::Default();

/// The cylinder with its seam's use at u = 2 pi moved to u = 0: the same
/// line in space, the wrong side of the seam in the lateral face's
/// parameters, so the loop jumps a period twice.
ShapeBuilder SeamOnWrongSide()
{
    Body body = MakeCylinder();
    const auto atZero = Must(body.Builder.AddPCurve(
        Stereon::Curve2d(Stereon::Line2d{Point2{0.0, 0.0}, Vec2{0.0, 1.0}})));
    Faults::SetPCurve(body.Builder, body.Coedges[1], atZero);
    return std::move(body.Builder);
}

/// A wire whose first vertex has the given tolerance.
ShapeBuilder WireWithVertexTolerance(const double tolerance)
{
    Body body;
    const auto a = body.AddVertex({0, 0, 0}, tolerance);
    const auto b = body.AddVertex({1, 0, 0});
    const auto line =
        Must(body.Builder.AddCurve(Stereon::Curve(Stereon::Line{Point3{0, 0, 0}, Vec3{1, 0, 0}})));
    body.AddEdge({.Kind = EdgeCurveKind::Spatial,
                  .Start = a,
                  .End = b,
                  .Curve = line,
                  .Range = Interval(0.0, 1.0),
                  .Tolerance = Tight});
    return std::move(body.Builder);
}

/// A wire whose end vertex is 1e-3 off its line's end.
ShapeBuilder VertexOffItsEdge()
{
    Body body;
    const auto a = body.AddVertex({0, 0, 0});
    const auto b = body.AddVertex({1, 1e-3, 0});
    const auto line =
        Must(body.Builder.AddCurve(Stereon::Curve(Stereon::Line{Point3{0, 0, 0}, Vec3{1, 0, 0}})));
    body.AddEdge({.Kind = EdgeCurveKind::Spatial,
                  .Start = a,
                  .End = b,
                  .Curve = line,
                  .Range = Interval(0.0, 1.0),
                  .Tolerance = Tight});
    return std::move(body.Builder);
}

/// The hinge with its Intersection edge relabelled Trim: a Trim edge
/// between faces on different planes.
ShapeBuilder TrimAcrossSurfaces()
{
    Body body = MakeHinge();
    Faults::SetKind(body.Builder, body.Edges[0], EdgeCurveKind::Trim);
    return std::move(body.Builder);
}

/// The hinge with both faces put on one plane: an Intersection edge on one
/// surface.
ShapeBuilder IntersectionOnOneSurface()
{
    Body body = MakeHinge();
    const auto plane = Must(body.Builder.AddSurface(
        Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
    Must(body.Builder.SetSurface(body.Faces[0], plane));
    Must(body.Builder.SetSurface(body.Faces[1], plane));
    return std::move(body.Builder);
}

/// The hinge with a cache of its Intersection edge `offset` above the edge,
/// marked as `status` with `bound`.
ShapeBuilder HingeWithCache(const double offset, const Stereon::CacheCertificate status,
                            const double bound)
{
    Body body = MakeHinge();
    const auto cache = Must(
        body.Builder.AddCurve(Stereon::Curve(Stereon::Line{Point3{0, 0, offset}, Vec3{0, 1, 0}})));
    Faults::SetCache(body.Builder, body.Edges[0], cache, status, bound);
    return std::move(body.Builder);
}

/// The box with face 2's plane moved 1e-3 off its edges, whose tolerance
/// does not hold the gap.
ShapeBuilder FaceOffItsEdges()
{
    Body body = MakeBox();
    const auto moved = Must(body.Builder.AddSurface(
        Stereon::Surface(Stereon::Plane{{0, 1e-3, 0}, Vec3{1, 0, 0}, Vec3{0, 0, 1}})));
    Must(body.Builder.SetSurface(body.Faces[2], moved));
    return std::move(body.Builder);
}

/// The hinge with its primary pcurve running out of a unit patch.
ShapeBuilder PcurveOutsideItsDomain()
{
    const Stereon::Surface patch(*Stereon::MakeTrimmedSurface(
        Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}),
        Stereon::Box2{Interval(0.0, 1.0), Interval(0.0, 1.0)}));
    Body body = MakeHinge(patch, Segment2({0, 0}, {0, 3}), [](Body&) {});
    return std::move(body.Builder);
}

/// The box with one face listed clockwise: its normal points in, which only
/// the Solid level's orientation rule can see.
ShapeBuilder FlippedFace()
{
    const std::vector<Point3> points{{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {0, 3, 0},
                                     {0, 0, 5}, {2, 0, 5}, {2, 3, 5}, {0, 3, 5}};
    const std::vector<std::vector<std::size_t>> faces{{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                      {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    return std::move(MakePolyhedron(points, faces).Builder);
}

/// The box without its top: a shell of a solid that does not close.
ShapeBuilder OpenShell()
{
    const std::vector<Point3> points{{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {0, 3, 0},
                                     {0, 0, 5}, {2, 0, 5}, {2, 3, 5}, {0, 3, 5}};
    const std::vector<std::vector<std::size_t>> faces{
        {0, 3, 2, 1}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    return std::move(MakePolyhedron(points, faces).Builder);
}

/// A seeded fault: what it is, how to make it, the rule that must find it
/// and the level at which it must be found.
struct SeededFault {
    std::string Name;
    std::function<ShapeBuilder()> Make;
    TopologyRule Rule;
    ValidityLevel Level;
};

std::vector<SeededFault> Corpus()
{
    const auto missingFromRing = [] {
        // The fan's spine has three uses, coedges 0, 3 and 6; the ring
        // short-circuits at its anchor and drops the other two.
        Body body = MakeFan(3);
        Faults::SetRadial(body.Builder, body.Coedges[0], body.Coedges[0]);
        return std::move(body.Builder);
    };
    const auto inTwoLoops = [] {
        Body body = MakeBox();
        Faults::SetCoedgeLoop(body.Builder, body.Coedges[0], body.Loops[1]);
        return std::move(body.Builder);
    };
    using L = ValidityLevel;
    using R = TopologyRule;
    using C = Stereon::CacheCertificate;
    return {
        {"a coedge missing from its edge's ring", missingFromRing, R::CoedgeOutsideRing,
         L::Structural},
        {"a coedge claimed by two loops", inTwoLoops, R::CoedgeOutsideLoop, L::Structural},
        {"a vertex off its edge's end", VertexOffItsEdge, R::VertexContainment, L::Geometric},
        {"an edge tolerance below the measured gap", FaceOffItsEdges, R::EdgeContainment,
         L::Geometric},
        {"a pcurve outside its surface's domain", PcurveOutsideItsDomain, R::EdgeContainment,
         L::Geometric},
        {"a seam use on the wrong side of the seam", SeamOnWrongSide, R::LoopParameterGap,
         L::Geometric},
        {"a tolerance above the Context's ceiling", [] { return WireWithVertexTolerance(2e-2); },
         R::ToleranceBounds, L::Geometric},
        {"a Trim edge between faces on different surfaces", TrimAcrossSurfaces, R::EdgeSurfaces,
         L::Geometric},
        {"an Intersection edge with both faces on one surface", IntersectionOnOneSurface,
         R::EdgeSurfaces, L::Geometric},
        {"an Intersection edge with a stale certificate",
         [] { return HingeWithCache(1e-3, C::Certified, 1e-9); }, R::CacheDeviation, L::Geometric},
        // The Solid level has no rules yet: its faults are refused by L01,
        // never passed, until orientation and closure rules replace it.
        {"a flipped face", FlippedFace, R::SolidUnchecked, L::Solid},
        {"an open shell", OpenShell, R::SolidUnchecked, L::Solid},
    };
}

/// Every level of the checker on a frozen shape, the Context's rules on.
ValidityReport At(const Shape& shape, const ValidityLevel level)
{
    return Stereon::ValidateTopology(shape, level, ctx);
}

} // namespace

TEST_CASE("each seeded fault is found by its rule, at its level and no lower",
          "[topo][validity][corpus]")
{
    for (const SeededFault& fault : Corpus()) {
        INFO(fault.Name);
        ShapeBuilder builder = fault.Make();
        const ValidityReport structural = builder.Validate();
        INFO(structural.ToText());
        if (fault.Level == ValidityLevel::Structural) {
            REQUIRE(structural.Has(fault.Rule));
            REQUIRE_FALSE(Stereon::Freeze(std::move(builder)).has_value());
            continue;
        }
        // Above Structural the fault is invisible to the structural check,
        // so Freeze certifies it, and the levels below its own pass it.
        REQUIRE(structural.IsValid());
        const Shape shape = Must(Stereon::Freeze(std::move(builder)));
        REQUIRE(At(shape, ValidityLevel::Structural).IsValid());
        if (fault.Level == ValidityLevel::Solid) {
            const ValidityReport geometric = At(shape, ValidityLevel::Geometric);
            INFO(geometric.ToText());
            REQUIRE(geometric.IsValid());
        }
        const ValidityReport found = At(shape, fault.Level);
        INFO(found.ToText());
        REQUIRE(found.Has(fault.Rule));
    }
}

TEST_CASE("no finding on the valid corpus at any implemented level", "[topo][validity][corpus]")
{
    std::vector<std::pair<std::string, Shape>> bodies;
    bodies.emplace_back("box", MakeBox().Freeze());
    bodies.emplace_back("prism", MakePrism(7).Freeze());
    bodies.emplace_back("cylinder", MakeCylinder().Freeze());
    bodies.emplace_back("square sheet", MakeSquareSheet().Freeze());
    bodies.emplace_back("sheet with a hole", MakeHoleSheet().Freeze());
    bodies.emplace_back("fan", MakeFan(3).Freeze());
    bodies.emplace_back("wire", MakeWire(3).Freeze());
    bodies.emplace_back("pole cap", MakePoleCap().Freeze());
    bodies.emplace_back("hinge", MakeHinge().Freeze());
    bodies.emplace_back(
        "hinge with a cache in its tube",
        Must(Stereon::Freeze(HingeWithCache(1e-6, Stereon::CacheCertificate::Certified, 2e-6))));
    for (const auto& [name, shape] : bodies) {
        INFO(name);
        const ValidityReport geometric = At(shape, ValidityLevel::Geometric);
        INFO(geometric.ToText());
        REQUIRE(geometric.IsValid());
        // At the Solid level a shape with no solid has nothing unchecked; one
        // with solids has exactly one L01 per solid and nothing else.
        const ValidityReport solid = At(shape, ValidityLevel::Solid);
        for (const Stereon::ValidityFinding& finding : solid.GetFindings()) {
            REQUIRE(finding.Rule == TopologyRule::SolidUnchecked);
        }
    }
}

TEST_CASE("G03: a loop jumps a period across a seam", "[topo][validity]")
{
    const Shape shape = Must(Stereon::Freeze(SeamOnWrongSide()));
    const ValidityReport report = At(shape, ValidityLevel::Geometric);
    INFO(report.ToText());
    // Into the moved use from the bottom circle at u = 2 pi, and out of it
    // into the top circle, which starts at u = 2 pi.
    REQUIRE(report.GetCount() == 2);
    for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
        REQUIRE(finding.Rule == TopologyRule::LoopParameterGap);
        REQUIRE(finding.Entities.size() == 3);
        REQUIRE(finding.Entities[0].starts_with("Loop#"));
    }
    // The seam is one line in space, so containment alone passes it.
    REQUIRE_FALSE(report.Has(TopologyRule::EdgeContainment));
}

TEST_CASE("G04: tolerances outside the Context's resolution and ceiling", "[topo][validity]")
{
    for (const double tolerance : {2e-2, 1e-7}) {
        INFO(tolerance);
        const Shape shape = Must(Stereon::Freeze(WireWithVertexTolerance(tolerance)));
        const ValidityReport report = At(shape, ValidityLevel::Geometric);
        REQUIRE(report.GetCount() == 1);
        REQUIRE(report.GetFindings()[0].Rule == TopologyRule::ToleranceBounds);
        REQUIRE(report.GetFindings()[0].Entities.at(0).starts_with("Vertex#0"));
        // Without the Context the bounds are not known, and not checked.
        REQUIRE(Stereon::ValidateTopology(shape, ValidityLevel::Geometric).IsValid());
    }
    // A Context whose ceiling holds the tolerance accepts it.
    const Shape loose = Must(Stereon::Freeze(WireWithVertexTolerance(2e-2)));
    REQUIRE(
        Stereon::ValidateTopology(loose, ValidityLevel::Geometric, ctx.WithToleranceCeiling(5e-2))
            .IsValid());
}

TEST_CASE("G05: Trim and Intersection edges and their surfaces", "[topo][validity]")
{
    const Shape trim = Must(Stereon::Freeze(TrimAcrossSurfaces()));
    const ValidityReport across = At(trim, ValidityLevel::Geometric);
    INFO(across.ToText());
    REQUIRE(across.Has(TopologyRule::EdgeSurfaces));
    REQUIRE(across.GetFindings()[0].Message.find("Trim") != std::string::npos);

    const Shape one = Must(Stereon::Freeze(IntersectionOnOneSurface()));
    const ValidityReport onOne = At(one, ValidityLevel::Geometric);
    INFO(onOne.ToText());
    REQUIRE(onOne.Has(TopologyRule::EdgeSurfaces));

    // A square split along its diagonal by a Trim edge, both halves on one
    // plane: the shared pcurve of ADR-0009, valid.
    Body body;
    const auto v0 = body.AddVertex({0, 0, 0});
    const auto v1 = body.AddVertex({1, 0, 0});
    const auto v2 = body.AddVertex({1, 1, 0});
    const auto v3 = body.AddVertex({0, 1, 0});
    const auto edge = [&](const EdgeCurveKind kind, const auto a, const auto b) {
        return body.AddEdge(
            {.Kind = kind, .Start = a, .End = b, .Range = Interval(0.0, 1.0), .Tolerance = Tight});
    };
    const auto diagonal = edge(EdgeCurveKind::Trim, v0, v2);
    const auto plane = Must(body.Builder.AddSurface(
        Stereon::Surface(Stereon::Plane{{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
    const auto shell = body.Shells.emplace_back(Must(body.Builder.AddShell()));
    const auto half = [&](const auto corner, const Point2 at, const bool first) {
        const auto face = body.Faces.emplace_back(Must(body.Builder.AddFace(plane)));
        Must(body.Builder.AddFaceToShell(shell, face));
        const auto loop = body.AddLoop(face);
        const Point2 origin{0, 0};
        const Point2 far{1, 1};
        if (first) {
            body.AddCoedge(loop, edge(EdgeCurveKind::Boundary, v0, corner), EdgeUseSense::Same,
                           Segment2(origin, at));
            body.AddCoedge(loop, edge(EdgeCurveKind::Boundary, corner, v2), EdgeUseSense::Same,
                           Segment2(at, far));
            body.AddCoedge(loop, diagonal, EdgeUseSense::Opposite, Segment2(origin, far));
        } else {
            body.AddCoedge(loop, diagonal, EdgeUseSense::Same, Segment2(origin, far));
            body.AddCoedge(loop, edge(EdgeCurveKind::Boundary, v2, corner), EdgeUseSense::Same,
                           Segment2(far, at));
            body.AddCoedge(loop, edge(EdgeCurveKind::Boundary, corner, v0), EdgeUseSense::Same,
                           Segment2(at, origin));
        }
    };
    half(v1, Point2{1, 0}, true);
    half(v3, Point2{0, 1}, false);
    const Shape split = std::move(body).Freeze();
    const ValidityReport shared = At(split, ValidityLevel::Geometric);
    INFO(shared.ToText());
    REQUIRE(shared.IsValid());
}

TEST_CASE("G06: a cached curve out of its tube or past its certified bound", "[topo][validity]")
{
    using Stereon::CacheCertificate;
    // Certified within 1e-9 yet 1e-3 away: a stale certificate.
    const Shape stale =
        Must(Stereon::Freeze(HingeWithCache(1e-3, CacheCertificate::Certified, 1e-9)));
    const ValidityReport bound = At(stale, ValidityLevel::Geometric);
    INFO(bound.ToText());
    REQUIRE(bound.GetCount() == 1);
    REQUIRE(bound.GetFindings()[0].Rule == TopologyRule::CacheDeviation);
    REQUIRE(bound.GetFindings()[0].Message.find("stale") != std::string::npos);

    // Uncertified and 1e-3 away, outside the edge's tolerance.
    const Shape loose = Must(Stereon::Freeze(HingeWithCache(
        1e-3, CacheCertificate::Uncertified, std::numeric_limits<double>::infinity())));
    const ValidityReport tube = At(loose, ValidityLevel::Geometric);
    INFO(tube.ToText());
    REQUIRE(tube.GetCount() == 1);
    REQUIRE(tube.GetFindings()[0].Message.find("outside its tolerance") != std::string::npos);

    // Within the tube and within the bound.
    REQUIRE(At(Must(Stereon::Freeze(HingeWithCache(1e-6, CacheCertificate::Certified, 2e-6))),
               ValidityLevel::Geometric)
                .IsValid());
}

TEST_CASE("L01: the Solid level refuses what it cannot check", "[topo][validity]")
{
    // A box has one solid: unchecked at the Solid level, never valid by
    // default; Geometric still passes it.
    const Shape box = MakeBox().Freeze();
    const ValidityReport solid = At(box, ValidityLevel::Solid);
    REQUIRE(solid.GetCount() == 1);
    REQUIRE(solid.GetFindings()[0].Rule == TopologyRule::SolidUnchecked);
    REQUIRE(solid.GetFindings()[0].Entities.at(0).starts_with("Solid#0"));
    REQUIRE(Stereon::ValidateTopology(box, ValidityLevel::Solid).Has(TopologyRule::SolidUnchecked));
    REQUIRE(At(box, ValidityLevel::Geometric).IsValid());
    // A sheet has no solid, so nothing is left unchecked.
    REQUIRE(At(MakeSquareSheet().Freeze(), ValidityLevel::Solid).IsValid());
}

TEST_CASE("the new rules have stable IDs and names", "[topo][validity]")
{
    REQUIRE(Stereon::RuleId(TopologyRule::LoopParameterGap) == "G03");
    REQUIRE(Stereon::RuleId(TopologyRule::ToleranceBounds) == "G04");
    REQUIRE(Stereon::RuleId(TopologyRule::EdgeSurfaces) == "G05");
    REQUIRE(Stereon::RuleId(TopologyRule::CacheDeviation) == "G06");
    REQUIRE(Stereon::RuleId(TopologyRule::SolidUnchecked) == "L01");
    REQUIRE(Stereon::RuleName(TopologyRule::LoopParameterGap) == "LoopParameterGap");
    REQUIRE(Stereon::RuleName(TopologyRule::CacheDeviation) == "CacheDeviation");
    REQUIRE(Stereon::RuleName(TopologyRule::SolidUnchecked) == "SolidUnchecked");
}
