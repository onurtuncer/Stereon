// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string_view>
#include <type_traits>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Transform;
using Stereon::Vec3;

namespace {

Cylinder SomeCylinder(const double radius)
{
    const Transform tilt = Stereon::MakeRotation(Vec3{0.3, 1.0, -0.4}, 0.8);
    return {Point3{1.0, 0.0, 2.0}, radius, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Plane SomePlane()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 0.0, 1.0}, 0.6);
    return {Point3{0.0, -1.0, 0.5}, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

/// The frame every other test surface uses, so that each component of each
/// result is non-trivial.
struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{-0.7, 0.2, 1.0}, 1.3);
    return {Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Cone SomeCone(const double radius)
{
    const Frame f = TiltedFrame();
    return {Point3{0.0, 1.0, -1.0}, radius, 0.4, f.X, f.Y};
}

Sphere SomeSphere(const double radius)
{
    const Frame f = TiltedFrame();
    return {Point3{2.0, -1.0, 0.5}, radius, f.X, f.Y};
}

Torus SomeTorus(const double minorRadius)
{
    const Frame f = TiltedFrame();
    return {Point3{-1.0, 0.0, 2.0}, minorRadius + 2.0, minorRadius, f.X, f.Y, false};
}

/// A bicubic B-spline patch over the parameter box the tests below use.
Stereon::BSplineSurface SomeBSplineSurface()
{
    std::vector<Point3> poles;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            poles.push_back({i * 1.0, j * 1.0, (i - 1.5) * (j - 1.5) * 0.25});
        }
    }
    return {3,
            3,
            {0.0, 0.0, 0.0, 0.0, 6.3, 6.3, 6.3, 6.3},
            {-1.6, -1.6, -1.6, -1.6, 1.6, 1.6, 1.6, 1.6},
            poles,
            {}};
}

} // namespace

TEST_CASE("VisitBuiltin hands built-in surfaces to the visitor as themselves", "[surface][visit]")
{
    const Surface plane(SomePlane());
    const Surface cylinder(SomeCylinder(1.5));
    REQUIRE(BuiltinIndexOf(plane) == 0);
    REQUIRE(BuiltinIndexOf(cylinder) == 1);
    REQUIRE(BuiltinIndexOf(Surface(SomeCone(1.0))) == 2);
    REQUIRE(BuiltinIndexOf(Surface(SomeSphere(1.0))) == 3);
    REQUIRE(BuiltinIndexOf(Surface(SomeTorus(1.0))) == 4);
    REQUIRE(BuiltinIndexOf(Surface(SomeBSplineSurface())) == 5);
    STATIC_REQUIRE(Stereon::BuiltinSurfaceTypes::Size == 6);
    STATIC_REQUIRE(Stereon::GeometryTypeName<Cone> == std::string_view{"Cone"});
    STATIC_REQUIRE(Stereon::GeometryTypeName<Sphere> == std::string_view{"Sphere"});
    STATIC_REQUIRE(Stereon::GeometryTypeName<Torus> == std::string_view{"Torus"});

    const bool sawCylinder =
        VisitBuiltin(cylinder, []<typename G>(const G&) { return std::is_same_v<G, Cylinder>; });
    REQUIRE(sawCylinder);
    STATIC_REQUIRE(Stereon::GeometryTypeName<Cylinder> == std::string_view{"Cylinder"});
}

TEST_CASE("the surface fast path and erased path give identical results",
          "[surface][visit][property]")
{
    rc::prop("for every built-in surface type and every operation", [] {
        // v inside the sphere's latitude range, the narrowest bounded domain.
        const Point2 uv{*rc::gen::inRange(0, 6283) / 1000.0,
                        *rc::gen::inRange(-1570, 1570) / 1000.0};
        const int du = *rc::gen::inRange(0, 4);
        const int dv = *rc::gen::inRange(du == 0 ? 1 : 0, 3);
        const Point3 p{*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                       *rc::gen::inRange(-50, 50) / 10.0};
        const double size = *rc::gen::inRange(1, 40) / 8.0;
        for (const Surface& s :
             {Surface(SomePlane()), Surface(SomeCylinder(size)), Surface(SomeCone(size)),
              Surface(SomeSphere(size)), Surface(SomeTorus(size)), Surface(SomeBSplineSurface())}) {
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Eval(g, uv); }) == Eval(s, uv));
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Deriv(g, uv, du, dv); }) ==
                      Deriv(s, uv, du, dv));
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Project(g, p); }) ==
                      Project(s, p));
        }
    });
}
