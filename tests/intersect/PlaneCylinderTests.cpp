// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The plane–cylinder handler (the Boolean prototype's M3, ADR-0010,
// ADR-0029): a circle across the axis, two rulings, one tangent ruling or a
// proved separation along it, and an ellipse otherwise, every branch on the
// cylinder's angle with its pcurve on each surface; the ellipse's pcurve on
// the cylinder is a B-spline within the resolution of the sinusoid.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using IntersectTest::PcurveDeviation;
using Stereon::Box2;
using Stereon::Circle;
using Stereon::ContactType;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

/// Radius 5 about z.
const Surface aboutZ{Cylinder{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
/// Radius 5 about x, its angle from y towards z.
const Surface aboutX{Cylinder{Point3{}, 5.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true}};

Surface Horizontal(const double z)
{
    return Surface{Plane{Point3{0.0, 0.0, z}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
}

} // namespace

TEST_CASE("a plane across the axis cuts a circle on the cylinder's angle", "[intersect][cylinder]")
{
    const auto result = Intersect(Horizontal(3.0), aboutZ, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    const Circle* circle = branch.Curve3d.Target<Circle>();
    REQUIRE(circle != nullptr);
    REQUIRE(circle->Radius == 5.0);
    REQUIRE(Distance(circle->Center, Point3{0.0, 0.0, 3.0}) == 0.0);
    REQUIRE(branch.Contact == ContactType::Crossing);
    REQUIRE(branch.Transition == SideTransition::OppositeSide);
    REQUIRE(branch.Range == Interval{0.0, Stereon::CirclePeriod});
    REQUIRE(branch.OnA.Target<Stereon::Circle2d>() != nullptr);
    REQUIRE(branch.OnB.Target<Stereon::Line2d>() != nullptr);
    REQUIRE(PcurveDeviation(branch, Horizontal(3.0), aboutZ) <= 1e-14);
    REQUIRE(branch.IsCertified());
    REQUIRE(branch.Bound <= ctx.GetResolution());

    // Half the angular domain: half the circle. Outside the height: nothing.
    const Box2 half{Interval{0.0, pi}, Interval::Everything()};
    REQUIRE(
        Intersect(Horizontal(3.0), Box2::Everything(), aboutZ, half, ctx)->Branches.front().Range ==
        Interval{0.0, pi});
    const Box2 above{Interval{0.0, Stereon::CirclePeriod}, Interval{5.0, 9.0}};
    REQUIRE(Intersect(Horizontal(3.0), Box2::Everything(), aboutZ, above, ctx)->IsEmpty());
}

TEST_CASE("a left-handed cylinder's circle agrees on both surfaces", "[intersect][cylinder]")
{
    const Surface inward{Cylinder{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, false}};
    const auto result = Intersect(Horizontal(3.0), inward, ctx);
    REQUIRE(result->Branches.size() == 1);
    REQUIRE(PcurveDeviation(result->Branches.front(), Horizontal(3.0), inward) <= 1e-14);
    REQUIRE(Intersect(inward, Horizontal(3.0), ctx)->Branches.size() == 1);
}

TEST_CASE("a plane along the axis cuts two rulings, in order of angle", "[intersect][cylinder]")
{
    const auto result = Intersect(Horizontal(2.0), aboutX, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 2);
    double previous = -1.0;
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        const Stereon::Line* line = branch.Curve3d.Target<Stereon::Line>();
        REQUIRE(line != nullptr);
        REQUIRE_THAT(std::abs(line->Direction.X), WithinAbs(1.0, 1e-15));
        REQUIRE_THAT(line->Origin.Z, WithinAbs(2.0, 1e-14));
        REQUIRE_THAT(std::abs(line->Origin.Y), WithinAbs(std::sqrt(21.0), 1e-14));
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.Transition == SideTransition::OppositeSide);
        REQUIRE(PcurveDeviation(branch, Horizontal(2.0), aboutX) <= 1e-13);
        const double u = branch.OnB.Target<Stereon::Line2d>()->Origin.X;
        REQUIRE(u > previous);
        previous = u;
    }

    // Clipped by the cylinder's height and the plane's box.
    const Box2 planeBox{Interval{-1.0, 4.0}, Interval{-9.0, 9.0}};
    const Box2 cylinderBox{Interval{0.0, Stereon::CirclePeriod}, Interval{1.0, 7.0}};
    const auto clipped = Intersect(Horizontal(2.0), planeBox, aboutX, cylinderBox, ctx);
    REQUIRE(clipped->Branches.size() == 2);
    for (const Stereon::IntersectionBranch& branch : clipped->Branches) {
        REQUIRE_THAT(branch.Range.GetLower(), WithinAbs(1.0, 1e-14));
        REQUIRE_THAT(branch.Range.GetUpper(), WithinAbs(4.0, 1e-14));
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound < 1e-12);
    }
}

TEST_CASE("a plane touching the cylinder gives one tangent ruling, same side",
          "[intersect][cylinder]")
{
    const auto result = Intersect(Horizontal(5.0), aboutX, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    REQUIRE(branch.Contact == ContactType::Tangential);
    REQUIRE(branch.Transition == SideTransition::SameSide);
    REQUIRE(Distance(branch.Curve3d.Target<Stereon::Line>()->Origin, Point3{0.0, 0.0, 5.0}) <=
            1e-14);
    REQUIRE(PcurveDeviation(branch, Horizontal(5.0), aboutX) <= 1e-13);
    // From below as well.
    REQUIRE(Intersect(Horizontal(-5.0), aboutX, ctx)->Branches.front().Contact ==
            ContactType::Tangential);
}

TEST_CASE("a plane along the axis beyond the radius is proved separate", "[intersect][cylinder]")
{
    REQUIRE(Intersect(Horizontal(5.0004), aboutX, ctx)->IsEmpty());
    REQUIRE(Intersect(Horizontal(-7.0), aboutX, ctx)->IsEmpty());
    REQUIRE(Intersect(aboutX, Horizontal(5.0004), ctx)->IsEmpty());
}

TEST_CASE("an oblique plane cuts an ellipse, exact on the plane and fitted on the cylinder",
          "[intersect][cylinder]")
{
    const Surface oblique{IntersectTest::PlaneWithNormal(Point3{}, Vec3{0.0, 1.0, 1.0})};
    const auto result = Intersect(oblique, aboutZ, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    const Ellipse* ellipse = branch.Curve3d.Target<Ellipse>();
    REQUIRE(ellipse != nullptr);
    REQUIRE_THAT(ellipse->MajorRadius, WithinAbs(5.0 * std::sqrt(2.0), 1e-14));
    REQUIRE(ellipse->MinorRadius == 5.0);
    REQUIRE(Distance(ellipse->Center, Point3{}) <= 1e-15);
    REQUIRE(branch.OnA.Target<Stereon::Ellipse2d>() != nullptr);
    REQUIRE(branch.OnB.Target<Stereon::BSplineCurve2d>() != nullptr);
    REQUIRE_THAT(Width(branch.Range), WithinAbs(Stereon::CirclePeriod, 1e-14));
    // Certified over the whole turn, the fitted pcurve included.
    REQUIRE(branch.IsCertified());
    REQUIRE(branch.Bound <= ctx.GetResolution());
    // The exact pcurve on the plane, the fitted one on the cylinder.
    for (const double t : IntersectTest::SamplesOf(branch.Range, 500)) {
        const Point3 p = Eval(branch.Curve3d, t);
        REQUIRE(Distance(p, Eval(oblique, Eval(branch.OnA, t))) <= 1e-13);
        REQUIRE(Distance(p, Eval(aboutZ, Eval(branch.OnB, t))) <= ctx.GetResolution());
        const Stereon::Point2 uv = Eval(branch.OnB, t);
        REQUIRE(uv.X >= -1e-12);
        REQUIRE(uv.X <= Stereon::CirclePeriod + 1e-12);
    }

    // Half the angular domain.
    const Box2 half{Interval{0.0, pi}, Interval::Everything()};
    const auto part = Intersect(oblique, Box2::Everything(), aboutZ, half, ctx);
    REQUIRE(part->Branches.size() == 1);
    REQUIRE_THAT(Width(part->Branches.front().Range), WithinAbs(pi, 1e-14));
    // Heights the ellipse never reaches (it stays within |z| <= 5).
    const Box2 above{Interval{0.0, Stereon::CirclePeriod}, Interval{6.0, 9.0}};
    REQUIRE(Intersect(oblique, Box2::Everything(), aboutZ, above, ctx)->IsEmpty());
}

TEST_CASE("ellipse pcurve work is bounded before converting or allocating its samples",
          "[intersect][cylinder]")
{
    const Surface oblique{IntersectTest::PlaneWithNormal(Point3{}, Vec3{0.0, 1.0, 1.0})};
    // The first count is finite but beyond the work budget, the second beyond
    // size_t, and the last underflows the spacing calculation to zero.
    for (const double resolution : {1e-20, 1e-100, std::numeric_limits<double>::denorm_min()}) {
        const auto result = Intersect(oblique, aboutZ, ctx.WithResolution(resolution));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == Stereon::ErrorCode::CertificationFailed);
        REQUIRE(std::string(result.error().GetMessage()).find("sample budget") !=
                std::string::npos);
    }
}

TEST_CASE("the cylinder first gives the same branches with the sides swapped",
          "[intersect][cylinder]")
{
    const Surface oblique{
        IntersectTest::PlaneWithNormal(Point3{1.0, 2.0, 3.0}, Vec3{1.0, 0.0, 2.0})};
    const auto swapped = Intersect(aboutZ, oblique, ctx);
    INFO(IntersectTest::ErrorText(swapped));
    REQUIRE(swapped.has_value());
    REQUIRE(swapped->Branches.size() == 1);
    REQUIRE(swapped->Branches.front().OnA.Target<Stereon::BSplineCurve2d>() != nullptr);
    REQUIRE(swapped->Branches.front().OnB.Target<Stereon::Ellipse2d>() != nullptr);
    REQUIRE(PcurveDeviation(swapped->Branches.front(), aboutZ, oblique, 300) <=
            ctx.GetResolution());
}
