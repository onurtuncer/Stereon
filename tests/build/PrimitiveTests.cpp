// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The primitive solids of the Boolean prototype (M1): a box and a cylinder,
// in any orthonormal frame, published by Commit with every entity
// Generated, valid at the Geometric level, faces oriented outward, and bad
// dimensions or frames refused.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

import Stereon.Build;

using Stereon::Context;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::HistoryRelation;
using Stereon::OperationResult;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ValidityLevel;
using Stereon::Vec3;

namespace {

/// A frame turned 30 degrees about z, then 45 about the new x: not aligned
/// with any axis, still orthonormal.
struct TiltedFrame {
    Vec3 X;
    Vec3 Y;
};

TiltedFrame Tilted()
{
    const double a = std::numbers::pi / 6.0;
    const double b = std::numbers::pi / 4.0;
    const Vec3 x{std::cos(a), std::sin(a), 0.0};
    const Vec3 y0{-std::sin(a), std::cos(a), 0.0};
    const Vec3 z0{0.0, 0.0, 1.0};
    return {x, std::cos(b) * y0 + std::sin(b) * z0};
}

void RequirePublished(const OperationResult& result, const std::size_t vertices,
                      const std::size_t edges, const std::size_t faces)
{
    const Shape& shape = result.Shape;
    REQUIRE(shape.Count<EntityKind::Vertex>() == vertices);
    REQUIRE(shape.Count<EntityKind::Edge>() == edges);
    REQUIRE(shape.Count<EntityKind::Face>() == faces);
    REQUIRE(shape.Count<EntityKind::Shell>() == 1);
    REQUIRE(shape.Count<EntityKind::Solid>() == 1);
    REQUIRE(result.Checked == ValidityLevel::Geometric);
    const Stereon::ValidityReport report =
        Stereon::ValidateTopology(shape, ValidityLevel::Geometric);
    INFO(report.ToText());
    REQUIRE(report.IsValid());
    // A new body: every entity Generated.
    REQUIRE(result.History.Count(HistoryRelation::Generated) == result.History.GetEntries().size());
}

/// The plane of a planar face; the test fails on another surface.
Stereon::Plane PlaneOf(const Shape& shape, const Stereon::FaceId face)
{
    const Stereon::Surface surface = Stereon::SurfaceOf(shape, face);
    const Stereon::Plane* plane = surface.Target<Stereon::Plane>();
    if (plane == nullptr) {
        FAIL("the face is not planar");
        return {};
    }
    return *plane;
}

/// The outward direction of a planar face: its plane's normal, flipped for
/// a reversed face.
Vec3 OutwardNormal(const Shape& shape, const Stereon::FaceId face)
{
    const Stereon::Plane plane = PlaneOf(shape, face);
    const Vec3 n = Stereon::Cross(plane.XAxis, plane.YAxis);
    return Stereon::IsReversed(shape, face) ? -n : n;
}

/// A closed, consistently oriented solid of positive volume: at the Solid
/// level the checker finds nothing but the self-intersection check, which
/// that level does not make.
void RequireClosedSolid(const Shape& shape)
{
    const Stereon::ValidityReport report =
        Stereon::ValidateTopology(shape, ValidityLevel::Solid, Context::Default());
    for (const Stereon::ValidityFinding& finding : report.GetFindings()) {
        INFO(finding.Message);
        CHECK(finding.Rule == Stereon::TopologyRule::SelfIntersectionUnchecked);
    }
}

} // namespace

TEST_CASE("a box: 8 vertices, 12 line edges, 6 faces facing out", "[build][primitives]")
{
    const auto box = Stereon::MakeBoxSolid(Stereon::Box3{Stereon::Interval{1.0, 3.0},
                                                         Stereon::Interval{-1.0, 2.0},
                                                         Stereon::Interval{0.0, 5.0}},
                                           Context::Default());
    REQUIRE(box);
    RequirePublished(*box, 8, 12, 6);
    const Shape& shape = box->Shape;
    const Point3 centre{2.0, 0.5, 2.5};
    for (const Stereon::FaceId f : *Stereon::AllEntities<EntityKind::Face>(shape)) {
        const Point3 onFace = PlaneOf(shape, f).Origin;
        REQUIRE(Stereon::Dot(onFace - centre, OutwardNormal(shape, f)) > 0.0);
    }
    for (const Stereon::EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
        REQUIRE(Stereon::KindOf(shape, e) == Stereon::EdgeCurveKind::Exact);
        REQUIRE(Stereon::ToleranceOf(shape, e) == Context::Default().GetResolution());
    }
}

TEST_CASE("a box in a tilted frame is as valid", "[build][primitives]")
{
    const TiltedFrame frame = Tilted();
    const auto box = Stereon::MakeBoxSolid(Point3{10.0, -4.0, 2.0}, frame.X, frame.Y,
                                           Vec3{2.0, 3.0, 0.5}, Context::Default());
    REQUIRE(box);
    RequirePublished(*box, 8, 12, 6);
}

TEST_CASE("a cylinder: two circles, a seam, a lateral face and two disks", "[build][primitives]")
{
    const auto cylinder = Stereon::MakeCylinderSolid(Point3{1.0, 2.0, 3.0}, Vec3{1, 0, 0},
                                                     Vec3{0, 1, 0}, 1.5, 4.0, Context::Default());
    REQUIRE(cylinder);
    RequirePublished(*cylinder, 2, 3, 3);
    const Shape& shape = cylinder->Shape;
    std::vector<Stereon::EdgeCurveKind> kinds;
    for (const Stereon::EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
        kinds.push_back(Stereon::KindOf(shape, e));
    }
    REQUIRE(kinds == std::vector<Stereon::EdgeCurveKind>{Stereon::EdgeCurveKind::Exact,
                                                         Stereon::EdgeCurveKind::Exact,
                                                         Stereon::EdgeCurveKind::Exact});
    // The disks face down and up.
    const auto faces = *Stereon::AllEntities<EntityKind::Face>(shape);
    std::vector<Stereon::FaceId> list(faces.begin(), faces.end());
    REQUIRE(OutwardNormal(shape, list[1]).Z < 0.0);
    REQUIRE(OutwardNormal(shape, list[2]).Z > 0.0);
}

TEST_CASE("a cylinder in a tilted frame is as valid", "[build][primitives]")
{
    const TiltedFrame frame = Tilted();
    const auto cylinder = Stereon::MakeCylinderSolid(Point3{-3.0, 0.0, 7.0}, frame.X, frame.Y, 0.25,
                                                     10.0, Context::Default());
    REQUIRE(cylinder);
    RequirePublished(*cylinder, 2, 3, 3);
}

TEST_CASE("primitives refuse bad dimensions and frames", "[build][primitives]")
{
    const Context context = Context::Default();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double tiny = context.GetResolution() / 2.0;
    const Vec3 x{1, 0, 0};
    const Vec3 y{0, 1, 0};
    for (const Vec3& size :
         {Vec3{0.0, 1.0, 1.0}, Vec3{1.0, -1.0, 1.0}, Vec3{1.0, 1.0, nan}, Vec3{1.0, tiny, 1.0}}) {
        const auto refused = Stereon::MakeBoxSolid(Point3{}, x, y, size, context);
        REQUIRE_FALSE(refused);
        REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
    }
    // Axes that are not orthonormal.
    REQUIRE_FALSE(Stereon::MakeBoxSolid(Point3{}, x, Vec3{1, 1, 0}, Vec3{1, 1, 1}, context));
    REQUIRE_FALSE(Stereon::MakeBoxSolid(Point3{}, Vec3{2, 0, 0}, y, Vec3{1, 1, 1}, context));
    REQUIRE_FALSE(Stereon::MakeCylinderSolid(Point3{}, x, Vec3{1, 1, 0}, 1.0, 1.0, context));
    for (const auto& [radius, height] :
         std::vector<std::pair<double, double>>{{0.0, 1.0}, {1.0, -2.0}, {nan, 1.0}, {1.0, tiny}}) {
        const auto refused = Stereon::MakeCylinderSolid(Point3{}, x, y, radius, height, context);
        REQUIRE_FALSE(refused);
        REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("a cylinder seam keeps an exact line and both periodic pcurves",
          "[build][primitives][regression]")
{
    const Context ctx = Context::Default();
    const Point3 base{1.0, -2.0, 3.0};
    for (const TiltedFrame& frame : {TiltedFrame{{1, 0, 0}, {0, 1, 0}}, Tilted()}) {
        const auto cylinder = Stereon::MakeCylinderSolid(base, frame.X, frame.Y, 1.5, 4.0, ctx);
        REQUIRE(cylinder);
        const Shape& shape = cylinder->Shape;
        std::size_t seams = 0;
        for (const Stereon::EdgeId edge : *AllEntities<EntityKind::Edge>(shape)) {
            if (StartOf(shape, edge) == EndOf(shape, edge)) {
                continue;
            }
            ++seams;
            REQUIRE(KindOf(shape, edge) == Stereon::EdgeCurveKind::Exact);
            const auto curve = CurveOf(shape, edge);
            REQUIRE(curve);
            REQUIRE(curve->Target<Stereon::Line>() != nullptr);
            const auto uses = RadialCoedgesOf(shape, edge);
            REQUIRE(uses);
            std::size_t count = 0;
            for (const Stereon::CoedgeId use : *uses) {
                ++count;
                const auto pcurve = PCurveOf(shape, use);
                REQUIRE(pcurve);
                const auto surface = SurfaceOf(shape, FaceOf(shape, LoopOf(shape, use)));
                for (const double t : {0.0, 1.0, 2.0, 4.0}) {
                    const Point3 expected = base + 1.5 * frame.X + t * Cross(frame.X, frame.Y);
                    CHECK(Distance(Eval(*curve, t), expected) <= ctx.GetResolution());
                    CHECK(Distance(Eval(surface, Eval(*pcurve, t)), expected) <=
                          ctx.GetResolution());
                }
            }
            CHECK(count == 2);
        }
        CHECK(seams == 1);
    }
}

TEST_CASE("a sphere: two poles, a seam, two pole edges and one face", "[build][primitives]")
{
    const auto sphere = Stereon::MakeSphereSolid(Point3{1.0, 2.0, 3.0}, Vec3{1, 0, 0},
                                                 Vec3{0, 1, 0}, 1.5, Context::Default());
    REQUIRE(sphere);
    RequirePublished(*sphere, 2, 3, 1);
    const Shape& shape = sphere->Shape;
    // The seam first, an exact meridian semicircle; then the two poles.
    std::vector<Stereon::EdgeCurveKind> kinds;
    for (const Stereon::EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
        kinds.push_back(Stereon::KindOf(shape, e));
    }
    REQUIRE(kinds == std::vector<Stereon::EdgeCurveKind>{Stereon::EdgeCurveKind::Exact,
                                                         Stereon::EdgeCurveKind::Degenerate,
                                                         Stereon::EdgeCurveKind::Degenerate});
    const auto edges = *Stereon::AllEntities<EntityKind::Edge>(shape);
    const Stereon::EdgeId seam = *edges.begin();
    REQUIRE(Stereon::CurveOf(shape, seam)->Target<Stereon::Circle>() != nullptr);
    // The poles lie on the frame's z axis, the seam runs from south to north.
    const auto vertices = *Stereon::AllEntities<EntityKind::Vertex>(shape);
    std::vector<Stereon::VertexId> poles(vertices.begin(), vertices.end());
    REQUIRE(Stereon::Distance(Stereon::PositionOf(shape, poles[0]), Point3{1.0, 2.0, 1.5}) <=
            1e-15);
    REQUIRE(Stereon::Distance(Stereon::PositionOf(shape, poles[1]), Point3{1.0, 2.0, 4.5}) <=
            1e-15);
    REQUIRE(Stereon::StartOf(shape, seam) == poles[0]);
    REQUIRE(Stereon::EndOf(shape, seam) == poles[1]);
    RequireClosedSolid(shape);
}

TEST_CASE("an ellipsoid: the sphere's topology with a half ellipse for a seam",
          "[build][primitives]")
{
    // Semi-axes 3, 2, 1: the seam is an ellipse of radii 3 along x and 1
    // along z. Then 1, 1.5, 3: the longer radius is along z, so the seam
    // is the ellipse of radii 3 along z and 1 along -x, and the latitude is
    // shifted on it by a quarter turn; the poles stay at +-c along z.
    for (const auto& [a, b, c] : {std::tuple{3.0, 2.0, 1.0}, std::tuple{1.0, 1.5, 3.0}}) {
        INFO("semi-axes " << a << ", " << b << ", " << c);
        const auto ellipsoid = Stereon::MakeEllipsoidSolid(Point3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0},
                                                           a, b, c, Context::Default());
        REQUIRE(ellipsoid);
        RequirePublished(*ellipsoid, 2, 3, 1);
        const Shape& shape = ellipsoid->Shape;
        const auto faces = *Stereon::AllEntities<EntityKind::Face>(shape);
        const Stereon::Surface surface = Stereon::SurfaceOf(shape, *faces.begin());
        REQUIRE(surface.Target<Stereon::Ellipsoid>() != nullptr);
        const auto edges = *Stereon::AllEntities<EntityKind::Edge>(shape);
        const Stereon::EdgeId seam = *edges.begin();
        const Stereon::Ellipse* meridian =
            Stereon::CurveOf(shape, seam)->Target<Stereon::Ellipse>();
        REQUIRE(meridian != nullptr);
        REQUIRE(meridian->MajorRadius == std::max(a, c));
        REQUIRE(meridian->MinorRadius == std::min(a, c));
        // The seam's ends are the poles, and its middle the vertex at +a x.
        const Stereon::Interval range = Stereon::RangeOf(shape, seam);
        const Stereon::Curve curve = *Stereon::CurveOf(shape, seam);
        REQUIRE(Stereon::Distance(Eval(curve, range.GetLower()), Point3{0.0, 0.0, -c}) <= 1e-15);
        REQUIRE(Stereon::Distance(Eval(curve, range.GetUpper()), Point3{0.0, 0.0, c}) <= 1e-15);
        REQUIRE(Stereon::Distance(Eval(curve, 0.5 * (range.GetLower() + range.GetUpper())),
                                  Point3{a, 0.0, 0.0}) <= 1e-15);
        RequireClosedSolid(shape);
    }
}

TEST_CASE("a sphere and an ellipsoid in a tilted frame are as valid", "[build][primitives]")
{
    const TiltedFrame frame = Tilted();
    const auto sphere =
        Stereon::MakeSphereSolid(Point3{-1.0, 0.5, 2.0}, frame.X, frame.Y, 2.0, Context::Default());
    REQUIRE(sphere);
    RequirePublished(*sphere, 2, 3, 1);
    RequireClosedSolid(sphere->Shape);
    const auto ellipsoid = Stereon::MakeEllipsoidSolid(Point3{-1.0, 0.5, 2.0}, frame.X, frame.Y,
                                                       2.0, 1.0, 0.5, Context::Default());
    REQUIRE(ellipsoid);
    RequirePublished(*ellipsoid, 2, 3, 1);
    RequireClosedSolid(ellipsoid->Shape);
}

TEST_CASE("spheres and ellipsoids refuse bad dimensions and frames", "[build][primitives]")
{
    const Context ctx = Context::Default();
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    REQUIRE(code(Stereon::MakeSphereSolid(Point3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, 0.0, ctx)) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeSphereSolid(Point3{}, Vec3{1, 0, 0}, Vec3{1, 1, 0}, 1.0, ctx)) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoidSolid(Point3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, 1.0, -2.0, 1.0,
                                             ctx)) == ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoidSolid(Point3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, 1.0, 2.0,
                                             std::numeric_limits<double>::quiet_NaN(), ctx)) ==
            ErrorCode::InvalidInput);
}
