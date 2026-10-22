// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Cylinders with skew or intersecting axes (the Boolean prototype's M3b,
// ADR-0010): the configuration decided exactly from the axes' distance,
// the branch points from the discriminant's quartic, and each loop as
// certified B-spline pcurves on both cylinders with a cached 3D curve on
// the shared parameter.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using IntersectTest::PcurveDeviation;
using Stereon::Box2;
using Stereon::ContactType;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double TwoPi = 2.0 * std::numbers::pi;

/// Radius 3 about z through the origin.
const Surface aboutZ{Cylinder{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};

/// About x through `origin`.
Surface AboutX(const Point3& origin, const double radius)
{
    return Surface{Cylinder{origin, radius, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true}};
}

/// About y through `origin`.
Surface AboutY(const Point3& origin, const double radius)
{
    return Surface{Cylinder{origin, radius, Vec3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, true}};
}

/// The distance of p from the cylinder's axis.
double RadialDistance(const Surface& surface, const Point3& p)
{
    const Cylinder& c = *surface.Target<Cylinder>();
    const Vec3 offset = p - c.Origin;
    return Length(offset - Dot(offset, Axis(c)) * Axis(c));
}

/// Every branch certified, transverse, on both cylinders along its whole
/// turn, with its samples within the certified bound of both.
void RequireWellFormed(const Stereon::SurfaceSurfaceIntersection& result, const Surface& a,
                       const Surface& b)
{
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.Transition == SideTransition::OppositeSide);
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= ctx.GetResolution());
        REQUIRE(branch.Range == Interval{0.0, TwoPi});
        REQUIRE(branch.Curve3d.Target<Stereon::BSplineCurve>() != nullptr);
        REQUIRE(PcurveDeviation(branch, a, b, 400) <= branch.Bound);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 200)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE_THAT(RadialDistance(a, p), WithinAbs(a.Target<Cylinder>()->Radius, 2e-5));
            REQUIRE_THAT(RadialDistance(b, p), WithinAbs(b.Target<Cylinder>()->Radius, 2e-5));
        }
        // The cache closes: a loop.
        REQUIRE(Distance(Eval(branch.Curve3d, 0.0), Eval(branch.Curve3d, TwoPi)) <= 1e-9);
    }
}

} // namespace

TEST_CASE("skew cylinders further apart than their radii are proved separate", "[intersect][skew]")
{
    const auto result = Intersect(aboutZ, AboutX(Point3{0.0, 10.0, 0.0}, 2.0), ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->IsEmpty());
    // Just apart, as numbers: 5.0001 against 3 + 2.
    REQUIRE(Intersect(aboutZ, AboutX(Point3{0.0, 5.0001, 0.0}, 2.0), ctx)->IsEmpty());
}

TEST_CASE("skew cylinders touching from outside meet in one certified point", "[intersect][skew]")
{
    const Surface b = AboutX(Point3{0.0, 5.0, 0.0}, 2.0);
    const auto result = Intersect(aboutZ, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Points.size() == 1);
    const Stereon::IntersectionPoint& point = result->Points.front();
    REQUIRE(Distance(point.Position, Point3{0.0, 3.0, 0.0}) <= 1e-14);
    REQUIRE(point.Bound < 1e-12);
    REQUIRE(Distance(Eval(aboutZ, point.OnA), point.Position) <= 1e-14);
    REQUIRE(Distance(Eval(b, point.OnB), point.Position) <= 1e-14);
}

TEST_CASE("skew cylinders tangent along their curve are not supported yet", "[intersect][skew]")
{
    // The axes' distance equals the radii's difference: B inside A, touching.
    const auto inside = Intersect(aboutZ, AboutX(Point3{0.0, 1.0, 0.0}, 2.0), ctx);
    REQUIRE_FALSE(inside.has_value());
    REQUIRE(inside.error().GetCode() == ErrorCode::NotSupported);
    // Equal radii with meeting axes: the two ellipses, singular at two points.
    const auto crossed = Intersect(aboutZ, AboutX(Point3{}, 3.0), ctx);
    REQUIRE_FALSE(crossed.has_value());
    REQUIRE(crossed.error().GetCode() == ErrorCode::NotSupported);
}

TEST_CASE("a smaller cylinder through a larger one gives two loops either way round",
          "[intersect][skew]")
{
    const Surface through = AboutX(Point3{}, 2.0);
    // On the larger cylinder the loops lie between branch points: the
    // smaller one enters on one side and leaves on the other.
    const auto onLarger = Intersect(aboutZ, through, ctx);
    INFO(IntersectTest::ErrorText(onLarger));
    REQUIRE(onLarger.has_value());
    REQUIRE(onLarger->Branches.size() == 2);
    RequireWellFormed(*onLarger, aboutZ, through);
    // In order of the first branch point in [0, 2 pi): the loop around -x,
    // then the one around +x, which straddles the seam.
    const Point3 first = Eval(onLarger->Branches[0].Curve3d, 0.0);
    const Point3 second = Eval(onLarger->Branches[1].Curve3d, 0.0);
    REQUIRE(first.X < 0.0);
    REQUIRE(second.X > 0.0);

    // On the smaller cylinder every ruling crosses the larger twice: the
    // two branches v± over the whole turn.
    const auto onSmaller = Intersect(through, aboutZ, ctx);
    INFO(IntersectTest::ErrorText(onSmaller));
    REQUIRE(onSmaller.has_value());
    REQUIRE(onSmaller->Branches.size() == 2);
    RequireWellFormed(*onSmaller, through, aboutZ);
}

TEST_CASE("a cylinder entering another on one side gives one loop", "[intersect][skew]")
{
    const Surface b = AboutX(Point3{0.0, 4.0, 0.0}, 2.0);
    const auto result = Intersect(aboutZ, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    RequireWellFormed(*result, aboutZ, b);
    // The loop's angle on A stays around +y, between its branch points.
    for (const double t : IntersectTest::SamplesOf(Interval{0.0, TwoPi}, 100)) {
        const double u = Eval(result->Branches.front().OnA, t).X;
        REQUIRE(u > 0.5);
        REQUIRE(u < std::numbers::pi - 0.5);
    }
}

TEST_CASE("a loop across the first cylinder's seam keeps its pcurve continuous",
          "[intersect][skew]")
{
    // B's axis is nearest A's ruling at angle 0.
    const Surface b = AboutY(Point3{4.0, 0.0, 0.0}, 2.0);
    const auto result = Intersect(aboutZ, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    RequireWellFormed(*result, aboutZ, b);
    double previous = Eval(result->Branches.front().OnA, 0.0).X;
    bool beyond = false;
    for (const double t : IntersectTest::SamplesOf(Interval{0.0, TwoPi}, 200)) {
        const double u = Eval(result->Branches.front().OnA, t).X;
        REQUIRE(std::abs(u - previous) < 0.2);
        beyond = beyond || u < 0.0 || u > TwoPi;
        previous = u;
    }
    REQUIRE(beyond); // Past the seam, rather than wrapped.
}

TEST_CASE("nearly parallel cylinders are reported as not supported", "[intersect][skew]")
{
    // Axes a nanoradian apart, 0.5 mm apart: the loops run billions of
    // millimetres along them, out of the size box.
    const Surface nearly{
        IntersectTest::CylinderAbout(Point3{5.0, 0.5, 0.0}, Vec3{1e-9, 0.0, 1.0}, 2.0)};
    const auto result = Intersect(aboutZ, nearly, ctx);
    REQUIRE_FALSE(result.has_value());
    INFO(result.error().GetMessage());
    REQUIRE(result.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(std::string(result.error().GetMessage()).find("size box") != std::string::npos);
    // The corpus's pair, with equal radii and axes that meet far away, is
    // the singular case instead: tangent along its curve, decided exactly.
    const Surface meeting{
        IntersectTest::CylinderAbout(Point3{5.0, 0.0, 0.0}, Vec3{1e-9, 0.0, 1.0}, 3.0)};
    const auto singular = Intersect(aboutZ, meeting, ctx);
    REQUIRE(singular.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(std::string(singular.error().GetMessage()).find("tangent") != std::string::npos);
}

TEST_CASE("loops outside the requested domains are left out", "[intersect][skew]")
{
    const Surface b = AboutX(Point3{0.0, 4.0, 0.0}, 2.0);
    const Box2 whole = Box2::Everything();
    const Box2 below{Interval{0.0, TwoPi}, Interval{-10.0, -5.0}};
    const auto away = Intersect(aboutZ, below, b, whole, ctx);
    REQUIRE(away.has_value());
    REQUIRE(away->IsEmpty());
    const Box2 around{Interval{0.0, TwoPi}, Interval{0.0, 10.0}};
    REQUIRE(Intersect(aboutZ, around, b, whole, ctx)->Branches.size() == 1);
}
