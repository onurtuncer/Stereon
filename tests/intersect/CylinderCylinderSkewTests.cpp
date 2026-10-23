// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The cylinder–cylinder handler for skew and intersecting axes (the Boolean
// prototype's M3b, ADR-0010's decisions on acceptance, ADR-0029): the
// quartic as certified B-spline branches, one loop or two, through nodes
// where the cylinders are internally tangent, two ellipses for equal radii
// with meeting axes, one point for external tangency, and a proved
// separation; every configuration decided exactly on the input doubles,
// and the branches clipped to the requested boxes.

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
using IntersectTest::ConnectivityOf;
using IntersectTest::ctx;
using IntersectTest::DistanceFrom;
using Stereon::Box2;
using Stereon::ContactType;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

/// Radius 3 about z through the origin, its angle from x towards y.
const Surface a{Cylinder{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};

/// A cylinder along x through `through`, its angle from y towards z.
Surface AlongX(const Point3& through, const double radius)
{
    return Surface{Cylinder{through, radius, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true}};
}

/// Every branch is a certified transverse crossing whose samples lie on
/// both cylinders within its bound, with pcurves inside both charts.
void RequireBranchesOn(const Stereon::SurfaceSurfaceIntersection& result, const Surface& first,
                       const Surface& second)
{
    const Cylinder& ca = *first.Target<Cylinder>();
    const Cylinder& cb = *second.Target<Cylinder>();
    REQUIRE_FALSE(result.Branches.empty());
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.Transition == SideTransition::OppositeSide);
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= ctx.GetResolution());
        REQUIRE(branch.Range.IsBounded());
        REQUIRE(branch.Curve3d.Target<Stereon::BSplineCurve>() != nullptr);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 200)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(DistanceFrom(ca, p) <= branch.Bound + 1e-12);
            REQUIRE(DistanceFrom(cb, p) <= branch.Bound + 1e-12);
            REQUIRE(Distance(p, Eval(first, Eval(branch.OnA, t))) <= branch.Bound);
            REQUIRE(Distance(p, Eval(second, Eval(branch.OnB, t))) <= branch.Bound);
            for (const Point2 uv : {Eval(branch.OnA, t), Eval(branch.OnB, t)}) {
                REQUIRE(uv.X >= -1e-9);
                REQUIRE(uv.X <= Stereon::CirclePeriod + 1e-9);
            }
        }
    }
}

} // namespace

TEST_CASE("skew cylinders crossing partly meet in one loop", "[intersect][cylinder][skew]")
{
    // B's axis is 4 from A's: between |3 - 2| and 3 + 2.
    const Surface b = AlongX(Point3{0.0, 4.0, 0.0}, 2.0);
    const auto result = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Points.empty());
    REQUIRE(result->Overlaps.empty());
    RequireBranchesOn(*result, a, b);
    const IntersectTest::Connectivity shape = ConnectivityOf(*result);
    REQUIRE(shape.Curves == 1);
    REQUIRE(shape.Nodes == 0);
    // The loop stays where both cylinders are: y between 2 and 3.
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 50)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(p.Y >= 2.0 - 1e-6);
            REQUIRE(p.Y <= 3.0 + 1e-6);
        }
    }
    // The same geometry with the cylinders swapped.
    const auto swapped = Intersect(b, a, ctx);
    INFO(IntersectTest::ErrorText(swapped));
    REQUIRE(swapped.has_value());
    RequireBranchesOn(*swapped, b, a);
    REQUIRE(ConnectivityOf(*swapped).Curves == 1);
}

TEST_CASE("a cylinder passing through another gives two loops", "[intersect][cylinder][skew]")
{
    // B of radius 1 through A of radius 3, its axis 1 from A's: two arcs on A.
    const Surface b = AlongX(Point3{0.0, 1.0, 0.0}, 1.0);
    const auto larger = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(larger));
    REQUIRE(larger.has_value());
    RequireBranchesOn(*larger, a, b);
    REQUIRE(ConnectivityOf(*larger).Curves == 2);
    REQUIRE(ConnectivityOf(*larger).Nodes == 0);
    // From the smaller cylinder every ruling crosses the larger twice: two
    // sheets over the whole turn.
    const auto smaller = Intersect(b, a, ctx);
    INFO(IntersectTest::ErrorText(smaller));
    REQUIRE(smaller.has_value());
    RequireBranchesOn(*smaller, b, a);
    REQUIRE(ConnectivityOf(*smaller).Curves == 2);
    // Each loop is where B enters or leaves A: x of one sign each.
    std::size_t positive = 0;
    for (const Stereon::IntersectionBranch& branch : smaller->Branches) {
        const Point3 p = Eval(branch.Curve3d, Midpoint(branch.Range));
        REQUIRE(std::abs(p.X) > 2.0);
        if (p.X > 0.0) {
            ++positive;
        }
    }
    REQUIRE(positive > 0);
    REQUIRE(positive < smaller->Branches.size());

    // Axes that meet, with different radii, are the same two loops.
    const Surface through = AlongX(Point3{}, 2.0);
    const auto meeting = Intersect(a, through, ctx);
    INFO(IntersectTest::ErrorText(meeting));
    REQUIRE(meeting.has_value());
    RequireBranchesOn(*meeting, a, through);
    REQUIRE(ConnectivityOf(*meeting).Curves == 2);
    REQUIRE(ConnectivityOf(Intersect(through, a, ctx).value()).Curves == 2);
}

TEST_CASE("equal radii with meeting axes give two ellipses through two nodes",
          "[intersect][cylinder][skew]")
{
    const Surface b = AlongX(Point3{}, 3.0);
    std::vector<std::string> log;
    const auto result = Intersect(a, Box2::Everything(), b, Box2::Everything(), ctx.GetResolution(),
                                  IntersectTest::Logging(log));
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireBranchesOn(*result, a, b);
    REQUIRE(result->Branches.size() >= 4);
    const IntersectTest::Connectivity shape = ConnectivityOf(*result);
    REQUIRE(shape.Nodes == 2);
    REQUIRE(log.size() == 1);
    REQUIRE(log.front().find("two ellipses") != std::string::npos);
    // The ellipses lie in the planes z = x and z = -x.
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 50)) {
            const Point3 p = Eval(branch.Curve3d, t);
            const double toPlanes =
                std::min(std::abs(p.Z - p.X), std::abs(p.Z + p.X)) / std::sqrt(2.0);
            REQUIRE(toPlanes <= 2.0 * branch.Bound + 1e-9);
        }
    }
    // The nodes are where the axes' common normal meets the cylinders.
    bool atNode = false;
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        const Point3 start = Eval(branch.Curve3d, branch.Range.GetLower());
        atNode = atNode || Distance(start, Point3{0.0, 3.0, 0.0}) <= 1e-6 ||
                 Distance(start, Point3{0.0, -3.0, 0.0}) <= 1e-6;
    }
    REQUIRE(atNode);
}

TEST_CASE("externally tangent skew cylinders touch at one point", "[intersect][cylinder][skew]")
{
    // B's axis is 5 = 3 + 2 from A's, along y.
    const Surface b = AlongX(Point3{0.0, 5.0, 0.0}, 2.0);
    std::vector<std::string> log;
    const auto result = Intersect(a, Box2::Everything(), b, Box2::Everything(), ctx.GetResolution(),
                                  IntersectTest::Logging(log));
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Points.size() == 1);
    const Stereon::IntersectionPoint& point = result->Points.front();
    REQUIRE(Distance(point.Position, Point3{0.0, 3.0, 0.0}) <= 1e-12);
    REQUIRE_THAT(point.OnA.X, WithinAbs(pi / 2.0, 1e-15));
    REQUIRE_THAT(point.OnA.Y, WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(point.OnB.X, WithinAbs(pi, 1e-15));
    REQUIRE_THAT(point.OnB.Y, WithinAbs(0.0, 1e-15));
    // Certified: the point is on both cylinders to rounding.
    REQUIRE(point.Bound <= 1e-12);
    REQUIRE(log.size() == 1);
    REQUIRE(log.front().find("touch at one point") != std::string::npos);
    // Outside the requested heights, nothing.
    const Box2 above{Interval{0.0, Stereon::CirclePeriod}, Interval{1.0, 2.0}};
    REQUIRE(Intersect(a, above, b, Box2::Everything(), ctx)->IsEmpty());
}

TEST_CASE("internally tangent cylinders meet at a node where the loops touch",
          "[intersect][cylinder][skew]")
{
    // B's axis is 2 = 3 - 1 from A's: B inside A, touching at (0, 3, 0).
    const Surface b = AlongX(Point3{0.0, 2.0, 0.0}, 1.0);
    std::vector<std::string> log;
    const auto larger = Intersect(a, Box2::Everything(), b, Box2::Everything(), ctx.GetResolution(),
                                  IntersectTest::Logging(log));
    INFO(IntersectTest::ErrorText(larger));
    REQUIRE(larger.has_value());
    RequireBranchesOn(*larger, a, b);
    REQUIRE(ConnectivityOf(*larger).Nodes == 1);
    REQUIRE(log.size() == 1);
    REQUIRE(log.front().find("internally tangent") != std::string::npos);
    std::size_t endsAtNode = 0;
    for (const Stereon::IntersectionBranch& branch : larger->Branches) {
        for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
            if (Distance(Eval(branch.Curve3d, t), Point3{0.0, 3.0, 0.0}) <= 1e-7) {
                ++endsAtNode;
            }
        }
    }
    REQUIRE(endsAtNode == 4);
    // From the smaller cylinder: two sheets over the whole turn, touching at the node.
    const auto smaller = Intersect(b, a, ctx);
    INFO(IntersectTest::ErrorText(smaller));
    REQUIRE(smaller.has_value());
    RequireBranchesOn(*smaller, b, a);
    REQUIRE(ConnectivityOf(*smaller).Nodes == 1);
}

TEST_CASE("skew cylinders apart are proved separate", "[intersect][cylinder][skew]")
{
    REQUIRE(Intersect(a, AlongX(Point3{0.0, 6.0, 0.0}, 2.0), ctx)->IsEmpty());
    REQUIRE(Intersect(a, AlongX(Point3{0.0, -5.5, 0.0}, 2.0), ctx)->IsEmpty());
    // Just apart as numbers: separate, whatever the tolerance says.
    REQUIRE(Intersect(a, Box2::Everything(), AlongX(Point3{0.0, 5.0005, 0.0}, 2.0),
                      Box2::Everything(), 1e-3, ctx)
                ->IsEmpty());
    // Just overlapping: a small loop.
    const Surface near = AlongX(Point3{0.0, 4.9995, 0.0}, 2.0);
    const auto loop = Intersect(a, near, ctx);
    INFO(IntersectTest::ErrorText(loop));
    REQUIRE(loop.has_value());
    RequireBranchesOn(*loop, a, near);
    REQUIRE(ConnectivityOf(*loop).Curves == 1);
}

TEST_CASE("the threshold between contact and separation is decided exactly",
          "[intersect][cylinder][skew]")
{
    // B's axis direction (3, 4, 0) normalised gives the doubles 0.6 and 0.8,
    // and the point (-4, 3, 0) is 5 from A's axis along the common normal in
    // the reals, but 5 - 5e-33 in the input doubles: the exact stage decides
    // a loop too small to see, which no rounded evaluation could.
    const Surface tilted{
        IntersectTest::CylinderAbout(Point3{-4.0, 3.0, 0.0}, Vec3{3.0, 4.0, 0.0}, 2.0)};
    const auto loop = Intersect(a, tilted, ctx);
    INFO(IntersectTest::ErrorText(loop));
    REQUIRE(loop.has_value());
    REQUIRE(loop->Points.empty());
    REQUIRE_FALSE(loop->Branches.empty());
    for (const Stereon::IntersectionBranch& branch : loop->Branches) {
        REQUIRE(branch.IsCertified());
        for (const double t : IntersectTest::SamplesOf(branch.Range, 20)) {
            REQUIRE(Distance(Eval(branch.Curve3d, t), Point3{-2.4, 1.8, 0.0}) <= 1e-6);
        }
    }
    // A hair further out, the same decision says separate.
    const Surface apart{IntersectTest::CylinderAbout(
        Point3{-4.0 * (1.0 + 1e-15), 3.0 * (1.0 + 1e-15), 0.0}, Vec3{3.0, 4.0, 0.0}, 2.0)};
    const auto separate = Intersect(a, apart, ctx);
    INFO(IntersectTest::ErrorText(separate));
    REQUIRE(separate.has_value());
    REQUIRE(separate->IsEmpty());
}

TEST_CASE("skew branches are clipped to the requested boxes", "[intersect][cylinder][skew]")
{
    const Surface b = AlongX(Point3{0.0, 4.0, 0.0}, 2.0);
    const Box2 whole = Box2::Everything();
    // B's height is x: the half of the loop with x >= 0.
    const Box2 halfB{Interval{0.0, Stereon::CirclePeriod}, Interval{0.0, 10.0}};
    const auto half = Intersect(a, whole, b, halfB, ctx);
    INFO(IntersectTest::ErrorText(half));
    REQUIRE(half.has_value());
    RequireBranchesOn(*half, a, b);
    double furthest = 0.0;
    for (const Stereon::IntersectionBranch& branch : half->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 50)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(p.X >= -1e-3);
            furthest = std::max(furthest, p.X);
        }
    }
    REQUIRE(furthest > 2.0);
    // Heights the loop never reaches.
    REQUIRE(
        Intersect(a, whole, b, Box2{Interval{0.0, Stereon::CirclePeriod}, Interval{5.0, 10.0}}, ctx)
            ->IsEmpty());
    REQUIRE(
        Intersect(a, Box2{Interval{0.0, Stereon::CirclePeriod}, Interval{4.0, 9.0}}, b, whole, ctx)
            ->IsEmpty());
    // A's angles: the loop lies about the common normal, y, at angle pi/2,
    // within acos(2/3) of it.
    const Box2 awayA{Interval{0.0, 0.5}, Interval::Everything()};
    REQUIRE(Intersect(a, awayA, b, whole, ctx)->IsEmpty());
    const Box2 edgeA{Interval{0.0, 1.0}, Interval::Everything()};
    const auto edge = Intersect(a, edgeA, b, whole, ctx);
    INFO(IntersectTest::ErrorText(edge));
    REQUIRE(edge.has_value());
    REQUIRE_FALSE(edge->Branches.empty());
    for (const Stereon::IntersectionBranch& branch : edge->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 50)) {
            REQUIRE(Eval(branch.OnA, t).X <= 1.0 + 1e-6);
        }
    }
}

TEST_CASE("nearly parallel cylinders are the quartic within bounded heights and refused beyond "
          "the size box",
          "[intersect][cylinder][skew]")
{
    // Axes 1e-9 rad from parallel and 5 apart, meeting far away; equal
    // radii, so two ellipses of the order of 1e9 mm, which leave the size
    // box: refused over unbounded heights (ADR-0001), two nearly straight
    // pieces within bounded ones.
    const Surface nearly{
        IntersectTest::CylinderAbout(Point3{5.0, 0.0, 0.0}, Vec3{1e-9, 0.0, 1.0}, 3.0)};
    const auto unbounded = Intersect(a, nearly, ctx);
    REQUIRE_FALSE(unbounded.has_value());
    REQUIRE(unbounded.error().GetCode() == ErrorCode::CertificationFailed);
    REQUIRE(std::string(unbounded.error().GetMessage()).find("size box") != std::string::npos);

    const Box2 heights{Interval{0.0, Stereon::CirclePeriod}, Interval{-10.0, 10.0}};
    const auto bounded = Intersect(a, heights, nearly, heights, ctx);
    INFO(IntersectTest::ErrorText(bounded));
    REQUIRE(bounded.has_value());
    RequireBranchesOn(*bounded, a, nearly);
    // The clip follows the heights to within 2 pi / 2^40 of A's angle, which
    // the ellipses' slope of 3e9 along the axes turns into a few hundredths.
    const double y = std::sqrt(9.0 - 6.25);
    std::size_t upper = 0;
    for (const Stereon::IntersectionBranch& branch : bounded->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 20)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(std::abs(p.Z) <= 10.0 + 0.1);
            REQUIRE_THAT(p.X, WithinAbs(2.5, 1e-6));
            REQUIRE_THAT(std::abs(p.Y), WithinAbs(y, 1e-6));
        }
        if (Eval(branch.Curve3d, Midpoint(branch.Range)).Y > 0.0) {
            ++upper;
        }
    }
    REQUIRE(upper > 0);
    REQUIRE(upper < bounded->Branches.size());
}
