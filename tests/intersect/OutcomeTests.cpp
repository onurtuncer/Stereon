// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The intersection outcomes of ADR-0029 (the Boolean prototype's M2): an
// empty successful result is a proved separation, a branch without a bound
// is uncertified, an overlap maps one surface's parameters onto the other's,
// and intersections are asked over requested domains, which must be valid.

#include <catch2/catch_test_macros.hpp>

#include <cmath>

import Stereon.Intersect;

using Stereon::Box2;
using Stereon::Context;
using Stereon::Curve;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

const Context ctx = Context::Default();

Curve SomeLine()
{
    return Curve(Stereon::Line{Point3{}, Vec3{1.0, 0.0, 0.0}});
}

Stereon::Surface SomePlane()
{
    return Stereon::Surface(Stereon::Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

} // namespace

TEST_CASE("an intersection with no component is a proved separation", "[intersect][outcome]")
{
    Stereon::SurfaceSurfaceIntersection result;
    REQUIRE(result.IsEmpty());
    result.Points.emplace_back();
    REQUIRE_FALSE(result.IsEmpty());

    Stereon::CurveCurveIntersection curves;
    REQUIRE(curves.IsEmpty());
    curves.Overlaps.emplace_back();
    REQUIRE_FALSE(curves.IsEmpty());
}

TEST_CASE("a branch is uncertified, and its side unresolved, until established",
          "[intersect][outcome]")
{
    Stereon::IntersectionBranch branch;
    REQUIRE_FALSE(branch.IsCertified());
    REQUIRE(branch.Transition == Stereon::SideTransition::Unresolved);
    branch.Bound = 1e-12;
    REQUIRE(branch.IsCertified());
    REQUIRE(std::isinf(Stereon::IntersectionPoint{}.Bound));
}

TEST_CASE("an overlap maps parameters from one surface to the other", "[intersect][outcome]")
{
    Stereon::ParameterMap identity;
    REQUIRE(identity(Point2{2.0, -3.0}) == Point2{2.0, -3.0});
    // A quarter turn and a shift, as two coplanar planes with turned frames
    // would have.
    const Stereon::ParameterMap turned{Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}, Vec2{5.0, 1.0}};
    REQUIRE(turned(Point2{2.0, 3.0}) == Point2{2.0, 3.0});
    REQUIRE(turned(Point2{1.0, 0.0}) == Point2{5.0, 2.0});

    Stereon::OverlapRegion overlap;
    REQUIRE(overlap.SameOrientation);
    REQUIRE(overlap.Evidence == Stereon::CoincidenceEvidence::Exact);
    REQUIRE(overlap.Distance == 0.0);
    REQUIRE(overlap.Boundary.empty());
}

TEST_CASE("intersections are asked over valid requested domains", "[intersect][outcome]")
{
    const Interval unit{0.0, 1.0};
    // An empty domain is refused before any handler runs. (A NaN-bounded
    // Interval cannot be made: its constructor's contract refuses it.)
    REQUIRE(Stereon::Intersect(SomeLine(), unit, SomeLine(), Interval{}, ctx).error().GetCode() ==
            ErrorCode::InvalidInput);
    const Box2 box{unit, unit};
    REQUIRE(Stereon::Intersect(SomePlane(), Box2{}, SomePlane(), box, ctx).error().GetCode() ==
            ErrorCode::InvalidInput);
    // Valid domains, including an unbounded one, reach the handler (a stub
    // until the analytic handlers land, M3).
    const auto stub = Stereon::Intersect(SomePlane(), Box2::Everything(), SomePlane(), box, ctx);
    REQUIRE(stub.error().GetCode() == ErrorCode::NotImplemented);
}
