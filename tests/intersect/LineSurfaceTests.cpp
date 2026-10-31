// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Lines against planes, spheres and cylinders, v0.2's first curve–surface
// handlers. Configurations are decided exactly on the input doubles: a
// line along a plane a 2^-40 gap away is apart, a line at exactly the
// radius is a touch, and a crossing of a plane is kept or dropped by its
// exact rational parameter. Every point lies on both geometries at its
// parameters, in order of the line's parameter; a stretch where the line
// lies in the surface carries a pcurve that traces it; and reversing the
// line negates every parameter without moving a point.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Interval everything = Interval::Everything();

/// The plane z = 0 with the axes of space, the unit sphere and the unit
/// cylinder about the z axis: exact frames, so that what is decided exactly
/// can be set up exactly.
const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Sphere unitSphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Cylinder unitCylinder{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

template <typename S>
Result<CurveSurfaceIntersection> Meet(const Line& line, const S& surface,
                                      const Interval& domain = everything)
{
    return Intersect(Curve{line}, domain, Surface{surface}, Domain(surface), ctx);
}

template <typename S>
Result<CurveSurfaceIntersection> MeetIn(const Line& line, const S& surface, const Interval& domain,
                                        const Box2& box)
{
    return Intersect(Curve{line}, domain, Surface{surface}, box, ctx);
}

/// Every point lies on the line and on the surface at its parameters, the
/// points come in order of the line's parameter, and every overlap's
/// pcurve traces the line on the surface over its range.
template <typename S>
void RequireOnBoth(const CurveSurfaceIntersection& result, const Line& line, const S& surface,
                   const double within = 1e-12)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const CurveSurfacePoint& p : result.Points) {
        REQUIRE(Distance(Eval(line, p.Parameter), p.Position) <= within);
        REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= within);
        REQUIRE(p.Parameter >= previous);
        previous = p.Parameter;
    }
    for (const Stereon::CurveSurfaceOverlap& overlap : result.Overlaps) {
        for (const double t : IntersectTest::SamplesOf(overlap.Range, 16)) {
            REQUIRE(Distance(Eval(surface, Eval(overlap.OnSurface, t)), Eval(line, t)) <= within);
        }
    }
}

std::size_t Crossings(const CurveSurfaceIntersection& result)
{
    return static_cast<std::size_t>(std::ranges::count_if(result.Points, [](const auto& p) {
        return p.Contact == ContactType::Crossing && p.Decided;
    }));
}

/// A line along x at height y and depth z, from far on the negative side.
Line AlongX(const double y, const double z = 0.0)
{
    return Line{Point3{-5.0, y, z}, Vec3{1.0, 0.0, 0.0}};
}

} // namespace

TEST_CASE("a line crosses a plane at an exact rational parameter", "[intersect][line-surface]")
{
    // From (1, 2, -1) upwards at speed 3: the crossing is at t = 1/3,
    // whose double rounds down.
    const Line up{Point3{1.0, 2.0, -1.0}, Vec3{0.0, 0.0, 3.0}};
    const auto result = Meet(up, ground);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 1);
    REQUIRE(result->Overlaps.empty());
    const CurveSurfacePoint& p = result->Points.front();
    REQUIRE(p.Contact == ContactType::Crossing);
    REQUIRE(p.Decided);
    REQUIRE(Distance(p.Position, Point3{1.0, 2.0, 0.0}) <= 1e-15);
    REQUIRE(p.OnSurface == Point2{1.0, 2.0});
    RequireOnBoth(*result, up, ground);
    // The domain decides on the exact 1/3, not on the rounded parameter.
    const double below = 1.0 / 3.0;
    REQUIRE(Meet(up, ground, Interval{0.0, below})->IsEmpty());
    REQUIRE(Meet(up, ground, Interval{0.0, std::nextafter(below, 1.0)})->Points.size() == 1);
    REQUIRE(Meet(up, ground, Interval{below, 1.0})->Points.size() == 1);
    REQUIRE(Meet(up, ground, Interval{std::nextafter(below, 1.0), 1.0})->IsEmpty());
    // So does the box, on the exact parameters of the point on the plane.
    REQUIRE(MeetIn(up, ground, everything, Box2{Interval{0.0, 1.0}, Interval{0.0, 2.0}})
                ->Points.size() == 1);
    REQUIRE(MeetIn(up, ground, everything,
                   Box2{Interval{0.0, std::nextafter(1.0, 0.0)}, Interval{0.0, 2.0}})
                ->IsEmpty());
    REQUIRE(MeetIn(up, ground, everything,
                   Box2{Interval{0.0, 1.0}, Interval{std::nextafter(2.0, 3.0), 3.0}})
                ->IsEmpty());
}

TEST_CASE("a line along a plane is apart or lies in it, decided exactly",
          "[intersect][line-surface]")
{
    const double gap = 0x1p-40;
    for (const double z : {gap, -gap}) {
        const auto apart = Meet(Line{Point3{0.0, 0.0, z}, Vec3{1.0, 2.0, 0.0}}, ground);
        REQUIRE(apart.has_value());
        REQUIRE(apart->IsEmpty());
    }
    const Line inPlane{Point3{1.0, 1.0, 0.0}, Vec3{2.0, 0.0, 0.0}};
    const auto full = Meet(inPlane, ground);
    REQUIRE(full.has_value());
    REQUIRE(full->Points.empty());
    REQUIRE(full->Overlaps.size() == 1);
    REQUIRE(full->Overlaps[0].Range == everything);
    RequireOnBoth(*full, inPlane, ground);
    // Clipped by the domain and by the box, on the line's own parameter:
    // u = 1 + 2 t in [0, 3] gives t in [-1/2, 1], the domain ends at 1/2.
    const auto clipped =
        MeetIn(inPlane, ground, Interval{-2.0, 0.5}, Box2{Interval{0.0, 3.0}, Interval{0.0, 2.0}});
    REQUIRE(clipped.has_value());
    REQUIRE(clipped->Overlaps.size() == 1);
    REQUIRE(Contains(clipped->Overlaps[0].Range, Interval{-0.5, 0.5}));
    REQUIRE(Width(clipped->Overlaps[0].Range) <= 1.0 + 1e-12);
    RequireOnBoth(*clipped, inPlane, ground);
    // A box the line misses.
    REQUIRE(MeetIn(inPlane, ground, everything, Box2{Interval{0.0, 3.0}, Interval{2.0, 3.0}})
                ->IsEmpty());
}

TEST_CASE("a line meets a sphere in two crossings, a touch or not at all, decided exactly",
          "[intersect][line-surface]")
{
    const auto two = Meet(AlongX(0.5), unitSphere);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    REQUIRE(Crossings(*two) == 2);
    REQUIRE(two->Overlaps.empty());
    for (const CurveSurfacePoint& p : two->Points) {
        REQUIRE(std::abs(Length(Stereon::ToVec(p.Position)) - 1.0) <= 1e-15);
    }
    REQUIRE(std::abs(two->Points[0].Position.X + std::sqrt(0.75)) <= 1e-15);
    REQUIRE(std::abs(two->Points[1].Position.X - std::sqrt(0.75)) <= 1e-15);
    RequireOnBoth(*two, AlongX(0.5), unitSphere);
    // At the radius one touch, decided; a hair beyond, apart; a hair
    // inside, two crossings.
    const auto touch = Meet(AlongX(1.0), unitSphere);
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(touch->Points[0].Contact == ContactType::Tangential);
    REQUIRE(touch->Points[0].Decided);
    REQUIRE(Distance(touch->Points[0].Position, Point3{0.0, 1.0, 0.0}) <= 1e-15);
    REQUIRE(Meet(AlongX(1.0 + 0x1p-40), unitSphere)->IsEmpty());
    REQUIRE(Meet(AlongX(1.0 - 0x1p-40), unitSphere)->Points.size() == 2);
    // The domain keeps the second crossing alone.
    REQUIRE(Meet(AlongX(0.5), unitSphere, Interval{5.0, 10.0})->Points.size() == 1);
    // The box: the angle is placed in its range modulo the period, the
    // latitude taken as it is.
    const Box2 nearSide{Interval{0.0, pi}, Interval{-pi / 2.0, pi / 2.0}};
    REQUIRE(MeetIn(AlongX(-0.5), unitSphere, everything, nearSide)->IsEmpty());
    REQUIRE(MeetIn(AlongX(0.5), unitSphere, everything, nearSide)->Points.size() == 2);
    const Box2 nextTurn{Interval{2.0 * pi, 3.0 * pi}, Interval{-pi / 2.0, pi / 2.0}};
    const auto placed = MeetIn(AlongX(0.5), unitSphere, everything, nextTurn);
    REQUIRE(placed.has_value());
    REQUIRE(placed->Points.size() == 2);
    for (const CurveSurfacePoint& p : placed->Points) {
        REQUIRE(p.OnSurface.X >= 2.0 * pi);
        REQUIRE(p.OnSurface.X <= 3.0 * pi);
    }
    RequireOnBoth(*placed, AlongX(0.5), unitSphere);
    const Box2 north{Interval{0.0, 2.0 * pi}, Interval{0.0, pi / 2.0}};
    REQUIRE(MeetIn(AlongX(0.0, -0.5), unitSphere, everything, north)->IsEmpty());
    REQUIRE(MeetIn(AlongX(0.0, 0.5), unitSphere, everything, north)->Points.size() == 2);
}

TEST_CASE("a line meets a cylinder in two crossings, a touch, or lies along a ruling, decided "
          "exactly",
          "[intersect][line-surface]")
{
    const auto two = Meet(AlongX(0.5, 3.0), unitCylinder);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    REQUIRE(Crossings(*two) == 2);
    for (const CurveSurfacePoint& p : two->Points) {
        REQUIRE(std::abs(std::hypot(p.Position.X, p.Position.Y) - 1.0) <= 1e-15);
        REQUIRE(p.OnSurface.Y == 3.0);
    }
    RequireOnBoth(*two, AlongX(0.5, 3.0), unitCylinder);
    const auto touch = Meet(AlongX(1.0, -2.0), unitCylinder);
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(touch->Points[0].Contact == ContactType::Tangential);
    REQUIRE(touch->Points[0].Decided);
    REQUIRE(Distance(touch->Points[0].Position, Point3{0.0, 1.0, -2.0}) <= 1e-15);
    REQUIRE(Meet(AlongX(1.0 + 0x1p-40, 0.0), unitCylinder)->IsEmpty());
    // Oblique, through the axis at 45 degrees: the crossings at t = 2 and 4.
    const Line oblique{Point3{-3.0, 0.0, -3.0}, Vec3{1.0, 0.0, 1.0}};
    const auto slanted = Meet(oblique, unitCylinder);
    REQUIRE(slanted.has_value());
    REQUIRE(slanted->Points.size() == 2);
    REQUIRE(std::abs(slanted->Points[0].Parameter - 2.0) <= 1e-15);
    REQUIRE(std::abs(slanted->Points[1].Parameter - 4.0) <= 1e-15);
    RequireOnBoth(*slanted, oblique, unitCylinder);
    // Along the axis: inside, apart beyond the radius, or on a ruling.
    REQUIRE(Meet(Line{Point3{0.5, 0.0, 0.0}, Vec3{0.0, 0.0, 2.0}}, unitCylinder)->IsEmpty());
    REQUIRE(
        Meet(Line{Point3{1.0 + 0x1p-40, 0.0, 0.0}, Vec3{0.0, 0.0, 2.0}}, unitCylinder)->IsEmpty());
    const Line ruling{Point3{1.0, 0.0, 1.0}, Vec3{0.0, 0.0, 2.0}};
    const auto on = Meet(ruling, unitCylinder);
    REQUIRE(on.has_value());
    REQUIRE(on->Points.empty());
    REQUIRE(on->Overlaps.size() == 1);
    REQUIRE(on->Overlaps[0].Range == everything);
    RequireOnBoth(*on, ruling, unitCylinder);
    // Clipped in v by the box, on the line's parameter: v = 1 + 2 t in
    // [-1, 5] gives t in [-1, 2], and the domain starts at 0.
    const auto clipped = MeetIn(ruling, unitCylinder, Interval{0.0, 10.0},
                                Box2{Interval{0.0, 2.0 * pi}, Interval{-1.0, 5.0}});
    REQUIRE(clipped.has_value());
    REQUIRE(clipped->Overlaps.size() == 1);
    REQUIRE(Contains(clipped->Overlaps[0].Range, Interval{0.0, 2.0}));
    REQUIRE(Width(clipped->Overlaps[0].Range) <= 2.0 + 1e-12);
    RequireOnBoth(*clipped, ruling, unitCylinder);
    // The ruling's angle placed in the box's range, or outside it.
    const auto turned =
        MeetIn(ruling, unitCylinder, everything, Box2{Interval{2.0 * pi, 3.0 * pi}, everything});
    REQUIRE(turned.has_value());
    REQUIRE(turned->Overlaps.size() == 1);
    REQUIRE(Eval(turned->Overlaps[0].OnSurface, 0.0).X == 2.0 * pi);
    RequireOnBoth(*turned, ruling, unitCylinder);
    REQUIRE(
        MeetIn(ruling, unitCylinder, everything, Box2{Interval{1.0, 2.0}, everything})->IsEmpty());
}

TEST_CASE("reversing the line negates every parameter and moves no point",
          "[intersect][line-surface][property]")
{
    rc::prop("for random lines against the plane, the sphere and the cylinder", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Line line{Point3{coordinate(), coordinate(), coordinate()},
                        Vec3{coordinate(), coordinate(), coordinate()}};
        RC_PRE(line.Direction != Vec3{});
        const Line back{line.Origin, -line.Direction};
        const auto check = [&](const auto& surface) {
            const auto forward = Meet(line, surface);
            const auto reversed = Meet(back, surface);
            RC_ASSERT(forward.has_value());
            RC_ASSERT(reversed.has_value());
            RC_ASSERT(forward->Points.size() == reversed->Points.size());
            RC_ASSERT(forward->Overlaps.size() == reversed->Overlaps.size());
            const std::size_t n = forward->Points.size();
            for (std::size_t i = 0; i < n; ++i) {
                const CurveSurfacePoint& p = forward->Points[i];
                const CurveSurfacePoint& q = reversed->Points[n - 1 - i];
                RC_ASSERT(std::abs(q.Parameter + p.Parameter) <=
                          1e-12 * (1.0 + std::abs(p.Parameter)));
                RC_ASSERT(Distance(q.Position, p.Position) <=
                          1e-12 * (1.0 + Length(Stereon::ToVec(p.Position))));
                RC_ASSERT(q.Contact == p.Contact);
            }
        };
        check(ground);
        check(unitSphere);
        check(unitCylinder);
    });
}

TEST_CASE("a line through the inside of a quadric crosses it twice, on both geometries",
          "[intersect][line-surface][property]")
{
    rc::prop("for spheres and cylinders placed at random", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const double radius = *rc::gen::inRange(1, 33) / 8.0;
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(direction != Vec3{});
        // Through a point within half the radius of the centre, so the
        // chord is at least the radius and the roots are well conditioned.
        const auto offset = [&] { return *rc::gen::inRange(-4, 5) / 8.0 * radius; };
        const Vec3 inside{offset(), offset(), offset()};
        const Line line{center + inside - 3.0 * direction, direction};
        const double within = 1e-9 * (1.0 + Length(Stereon::ToVec(center)) + radius);
        const Sphere sphere{center, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const auto onSphere = Meet(line, sphere);
        RC_ASSERT(onSphere.has_value());
        RC_ASSERT(onSphere->Points.size() == 2);
        RC_ASSERT(onSphere->Overlaps.empty());
        for (const CurveSurfacePoint& p : onSphere->Points) {
            RC_ASSERT(p.Contact == ContactType::Crossing);
            RC_ASSERT(Distance(Eval(line, p.Parameter), p.Position) <= within);
            RC_ASSERT(Distance(Eval(sphere, p.OnSurface), p.Position) <= within);
            RC_ASSERT(std::abs(Distance(p.Position, center) - radius) <= within);
        }
        // A cylinder about a random axis through the centre, the line not
        // along it (exactly, on these dyadic coordinates).
        const Vec3 axis{coordinate(), coordinate(), coordinate()};
        RC_PRE(axis != Vec3{} && Cross(axis, direction) != Vec3{});
        const Cylinder cylinder = IntersectTest::CylinderAbout(center, axis, radius);
        const auto onCylinder = Meet(line, cylinder);
        RC_ASSERT(onCylinder.has_value());
        RC_ASSERT(onCylinder->Points.size() == 2);
        RC_ASSERT(onCylinder->Overlaps.empty());
        for (const CurveSurfacePoint& p : onCylinder->Points) {
            RC_ASSERT(p.Contact == ContactType::Crossing);
            RC_ASSERT(Distance(Eval(line, p.Parameter), p.Position) <= within);
            RC_ASSERT(Distance(Eval(cylinder, p.OnSurface), p.Position) <= within);
            RC_ASSERT(IntersectTest::DistanceFrom(cylinder, p.Position) <= within);
        }
    });
    rc::prop("for planes placed at random, one crossing on both", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Vec3 normal{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(normal != Vec3{} && direction != Vec3{});
        // Not within rounding of parallel: the frame's axes are rounded.
        RC_PRE(std::abs(Dot(normal, direction)) >= 0.125 * Length(normal) * Length(direction));
        const Point3 origin{coordinate(), coordinate(), coordinate()};
        const Plane plane = IntersectTest::PlaneWithNormal(origin, normal);
        const Line line{Point3{coordinate(), coordinate(), coordinate()}, direction};
        const auto result = Meet(line, plane);
        RC_ASSERT(result.has_value());
        RC_ASSERT(result->Points.size() == 1);
        RC_ASSERT(result->Overlaps.empty());
        const CurveSurfacePoint& p = result->Points.front();
        const double within = 1e-9 * (1.0 + Length(Stereon::ToVec(p.Position)));
        RC_ASSERT(Distance(Eval(line, p.Parameter), p.Position) <= within);
        RC_ASSERT(Distance(Eval(plane, p.OnSurface), p.Position) <= within);
    });
}
