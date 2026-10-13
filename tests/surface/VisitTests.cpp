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

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Cylinder;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
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

} // namespace

TEST_CASE("VisitBuiltin hands built-in surfaces to the visitor as themselves", "[surface][visit]")
{
    const Surface plane(SomePlane());
    const Surface cylinder(SomeCylinder(1.5));
    REQUIRE(BuiltinIndexOf(plane) == 0);
    REQUIRE(BuiltinIndexOf(cylinder) == 1);
    STATIC_REQUIRE(Stereon::BuiltinSurfaceTypes::Size == 2);

    const bool sawCylinder =
        VisitBuiltin(cylinder, []<typename G>(const G&) { return std::is_same_v<G, Cylinder>; });
    REQUIRE(sawCylinder);
    STATIC_REQUIRE(Stereon::GeometryTypeName<Cylinder> == std::string_view{"Cylinder"});
}

TEST_CASE("the surface fast path and erased path give identical results",
          "[surface][visit][property]")
{
    rc::prop("for every built-in surface type and every operation", [] {
        const Point2 uv{*rc::gen::inRange(0, 6283) / 1000.0,
                        *rc::gen::inRange(-5000, 5000) / 1000.0};
        const int du = *rc::gen::inRange(0, 4);
        const int dv = *rc::gen::inRange(du == 0 ? 1 : 0, 3);
        const Point3 p{*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                       *rc::gen::inRange(-50, 50) / 10.0};
        for (const Surface& s :
             {Surface(SomePlane()), Surface(SomeCylinder(*rc::gen::inRange(1, 40) / 8.0))}) {
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Eval(g, uv); }) == Eval(s, uv));
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Deriv(g, uv, du, dv); }) ==
                      Deriv(s, uv, du, dv));
            RC_ASSERT(VisitBuiltin(s, [&](const auto& g) { return Project(g, p); }) ==
                      Project(s, p));
        }
    });
}
