// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Lines against planes, spheres, cylinders and ellipsoids: v0.2's first
// curve–surface handlers, and ADR-0048's first for the ellipsoid, checked
// against the sphere it is when its semi-axes are equal and against the
// sphere it is the image of when one axis is scaled.
// Configurations are decided exactly on the input doubles: a
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
#include <utility>
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
using Stereon::Ellipsoid;
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
/// The ellipsoid about the origin with the semi-axes 3, 2 and 1 along x, y
/// and z: an exact frame and dyadic semi-axes, so that a touch can be set
/// up exactly (ADR-0048).
const Ellipsoid stretched{Point3{}, 3.0, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

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

/// (x / A)² + (y / B)² + (z / C)² - 1 in the ellipsoid's frame.
double Level(const Ellipsoid& e, const Point3& p)
{
    const Vec3 w = p - e.Center;
    const double x = Dot(w, e.XAxis) / e.A;
    const double y = Dot(w, e.YAxis) / e.B;
    const double z = Dot(w, Axis(e)) / e.C;
    return x * x + y * y + z * z - 1.0;
}

/// An ellipsoid with its semi-axis c along the given direction.
Ellipsoid EllipsoidAbout(const Point3& center, const Vec3& axis, const double a, const double b,
                         const double c)
{
    const Cylinder frame = IntersectTest::CylinderAbout(center, axis, 1.0);
    return Ellipsoid{center, a, b, c, frame.XAxis, frame.YAxis};
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
    rc::prop("for random lines against the plane, the sphere, the cylinder and the ellipsoid", [] {
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
        check(stretched);
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

TEST_CASE("a crossing from a far line origin lies on the quadric to rounding",
          "[intersect][line-surface][regression]")
{
    // RapidCheck's counterexample of 2026-11-04 (the linearisation's
    // transpose property, seed 17734536904155436187): the line starts about
    // 24 units from a cylinder of radius 1/2, so the quadratic's
    // coefficients, formed from that offset, placed the crossings 1.4e-13
    // off the cylinder, about 600 units in the last place of the point, and
    // IntersectJvp rejected them as stale.
    const Point3 center{-11.0 / 8.0, -11.0 / 8.0, -12.0 / 8.0};
    constexpr double radius = 0.5;
    const Vec3 direction{-30.0 / 8.0, -40.0 / 8.0, -39.0 / 8.0};
    const Vec3 axis = Cross(direction, Vec3{-10.0 / 8.0, -20.0 / 8.0, -34.0 / 8.0});
    const Line line{center + Vec3{-0.125, -0.125, -0.1875} - 3.0 * direction, direction};
    const Sphere sphere{center, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Cylinder cylinder = IntersectTest::CylinderAbout(center, axis, radius);
    const auto onBoth = [&](const auto& surface) {
        const auto result = Intersect(Curve{line}, Surface{surface}, ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Points.size() == 2);
        for (const CurveSurfacePoint& p : result->Points) {
            const double ulps = 32.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, Length(Stereon::ToVec(p.Position)));
            REQUIRE(Distance(Eval(line, p.Parameter), p.Position) <= ulps);
            REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= ulps);
        }
    };
    onBoth(sphere);
    onBoth(cylinder);
}

TEST_CASE("a line's two crossings of a quadric stay apart however far its origin",
          "[intersect][line-surface][regression]")
{
    // The audit of 2026-11-05: from 10^8 away the quadratic's coefficients,
    // formed in doubles, cancelled the discriminant to zero, both crossings
    // were placed at the foot of the perpendicular from the centre, and the
    // polish, its step refused against a second root at no distance, left
    // them there. A line at distance 1/2 from the centre, its origin 13 L
    // back from the foot: the crossings sqrt(3) apart, each on the surface
    // to the rounding of the point, about L units in the last place.
    for (double length = 1e2; length <= 1e12; length *= 100.0) {
        const Vec3 direction{3.0, 4.0, 12.0};
        const Line line{Point3{-3.0 * length - 0.4, -4.0 * length + 0.3, -12.0 * length},
                        direction};
        const double within = 64.0 * std::numeric_limits<double>::epsilon() * 13.0 * length;
        const auto result = Intersect(Curve{line}, Surface{unitSphere}, ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Points.size() == 2);
        for (const CurveSurfacePoint& p : result->Points) {
            REQUIRE(p.Contact == ContactType::Crossing);
            REQUIRE(Distance(Eval(line, p.Parameter), p.Position) <= within);
            REQUIRE(Distance(Eval(unitSphere, p.OnSurface), p.Position) <= within);
        }
        REQUIRE(std::abs(Distance(result->Points[0].Position, result->Points[1].Position) -
                         std::sqrt(3.0)) <= 2.0 * within);
    }
    // The audit's configuration on the cylinder: a line along x, 1/2 from
    // the axis, from x = -L; the crossings at x = ±sqrt(3)/2.
    for (double length = 1e2; length <= 1e12; length *= 100.0) {
        const Line line{Point3{-length, 0.5, 0.25}, Vec3{1.0, 0.0, 0.0}};
        const auto result = Intersect(Curve{line}, Surface{unitCylinder}, ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Points.size() == 2);
        const double within = 64.0 * std::numeric_limits<double>::epsilon() * length;
        REQUIRE(std::abs(result->Points[0].Position.X + std::sqrt(0.75)) <= within);
        REQUIRE(std::abs(result->Points[1].Position.X - std::sqrt(0.75)) <= within);
        for (const CurveSurfacePoint& p : result->Points) {
            REQUIRE(Distance(Eval(unitCylinder, p.OnSurface), p.Position) <= within);
        }
    }
}

TEST_CASE("a line meets an ellipsoid in two crossings, a touch or not at all, decided exactly",
          "[intersect][line-surface][ellipsoid]")
{
    // Along x at y = 1/2: (x / 3)² + (1/4)² = 1 at x = ±3 sqrt(15) / 4.
    const double reach = 0.75 * std::sqrt(15.0);
    const auto two = Meet(AlongX(0.5), stretched);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    REQUIRE(Crossings(*two) == 2);
    REQUIRE(two->Overlaps.empty());
    REQUIRE(std::abs(two->Points[0].Position.X + reach) <= 4e-15);
    REQUIRE(std::abs(two->Points[1].Position.X - reach) <= 4e-15);
    for (const CurveSurfacePoint& p : two->Points) {
        REQUIRE(std::abs(Level(stretched, p.Position)) <= 4e-15);
    }
    RequireOnBoth(*two, AlongX(0.5), stretched);
    // At a semi-axis one touch, decided; a hair beyond, apart; a hair
    // inside, two crossings. Beside the equator, at y = 2:
    const auto touch = Meet(AlongX(2.0), stretched);
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(touch->Points[0].Contact == ContactType::Tangential);
    REQUIRE(touch->Points[0].Decided);
    REQUIRE(Distance(touch->Points[0].Position, Point3{0.0, 2.0, 0.0}) <= 1e-15);
    REQUIRE(Meet(AlongX(2.0 + 0x1p-40), stretched)->IsEmpty());
    REQUIRE(Meet(AlongX(2.0 - 0x1p-40), stretched)->Points.size() == 2);
    // Over the pole, at z = 1, where every longitude is the same point:
    const auto pole = Meet(AlongX(0.0, 1.0), stretched);
    REQUIRE(pole.has_value());
    REQUIRE(pole->Points.size() == 1);
    REQUIRE(pole->Points[0].Contact == ContactType::Tangential);
    REQUIRE(pole->Points[0].Decided);
    REQUIRE(Distance(pole->Points[0].Position, Point3{0.0, 0.0, 1.0}) <= 1e-15);
    REQUIRE(std::abs(pole->Points[0].OnSurface.Y - pi / 2.0) <= 1e-15);
    REQUIRE(Meet(AlongX(0.0, 1.0 + 0x1p-40), stretched)->IsEmpty());
    REQUIRE(Meet(AlongX(0.0, 1.0 - 0x1p-40), stretched)->Points.size() == 2);
    // And up the end of the longest semi-axis, at x = 3:
    const auto upAt = [](const double x) {
        return Line{Point3{x, 0.0, -5.0}, Vec3{0.0, 0.0, 1.0}};
    };
    const auto end = Meet(upAt(3.0), stretched);
    REQUIRE(end.has_value());
    REQUIRE(end->Points.size() == 1);
    REQUIRE(end->Points[0].Contact == ContactType::Tangential);
    REQUIRE(Distance(end->Points[0].Position, Point3{3.0, 0.0, 0.0}) <= 1e-15);
    REQUIRE(Meet(upAt(3.0 + 0x1p-40), stretched)->IsEmpty());
    REQUIRE(Meet(upAt(3.0 - 0x1p-40), stretched)->Points.size() == 2);
    // The domain keeps the second crossing alone.
    REQUIRE(Meet(AlongX(0.5), stretched, Interval{5.0, 10.0})->Points.size() == 1);
    // The box: the longitude is placed in its range modulo the period, the
    // latitude taken as it is.
    const Box2 nearSide{Interval{0.0, pi}, Interval{-pi / 2.0, pi / 2.0}};
    REQUIRE(MeetIn(AlongX(-0.5), stretched, everything, nearSide)->IsEmpty());
    REQUIRE(MeetIn(AlongX(0.5), stretched, everything, nearSide)->Points.size() == 2);
    const Box2 nextTurn{Interval{2.0 * pi, 3.0 * pi}, Interval{-pi / 2.0, pi / 2.0}};
    const auto placed = MeetIn(AlongX(0.5), stretched, everything, nextTurn);
    REQUIRE(placed.has_value());
    REQUIRE(placed->Points.size() == 2);
    for (const CurveSurfacePoint& p : placed->Points) {
        REQUIRE(p.OnSurface.X >= 2.0 * pi);
        REQUIRE(p.OnSurface.X <= 3.0 * pi);
    }
    RequireOnBoth(*placed, AlongX(0.5), stretched);
    const Box2 north{Interval{0.0, 2.0 * pi}, Interval{0.0, pi / 2.0}};
    REQUIRE(MeetIn(AlongX(0.0, -0.5), stretched, everything, north)->IsEmpty());
    REQUIRE(MeetIn(AlongX(0.0, 0.5), stretched, everything, north)->Points.size() == 2);
}

TEST_CASE("an ellipsoid with equal semi-axes gives the sphere's results",
          "[intersect][line-surface][ellipsoid][property]")
{
    // ADR-0048's oracle. On an exact frame the ellipsoid's quadratic is the
    // sphere's times r^4, so the count and the contact agree exactly, a
    // touch included; the points agree to rounding.
    const std::vector<std::pair<Vec3, Vec3>> frames{
        {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, {Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}},
        {Vec3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}}, {Vec3{0.0, 1.0, 0.0}, Vec3{1.0, 0.0, 0.0}},
        {Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 1.0, 0.0}}, {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};
    const auto same = [](const Line& line, const Sphere& sphere, const Ellipsoid& ellipsoid) {
        const auto onSphere = Meet(line, sphere);
        const auto onEllipsoid = Meet(line, ellipsoid);
        RC_ASSERT(onSphere.has_value());
        RC_ASSERT(onEllipsoid.has_value());
        RC_ASSERT(onSphere->Points.size() == onEllipsoid->Points.size());
        RC_ASSERT(onEllipsoid->Overlaps.empty());
        for (std::size_t i = 0; i < onSphere->Points.size(); ++i) {
            const CurveSurfacePoint& p = onSphere->Points[i];
            const CurveSurfacePoint& q = onEllipsoid->Points[i];
            const double within = 1e-12 * (1.0 + Length(Stereon::ToVec(p.Position)));
            RC_ASSERT(q.Contact == p.Contact);
            RC_ASSERT(q.Decided == p.Decided);
            RC_ASSERT(std::abs(q.Parameter - p.Parameter) <= 1e-12 * (1.0 + std::abs(p.Parameter)));
            RC_ASSERT(Distance(q.Position, p.Position) <= within);
            // The same point of the surface, whatever longitude a pole gets.
            RC_ASSERT(Distance(Eval(ellipsoid, q.OnSurface), Eval(sphere, p.OnSurface)) <= within);
            RC_ASSERT(std::abs(q.OnSurface.Y - p.OnSurface.Y) <= 1e-7);
        }
        return onSphere->Points.size();
    };
    rc::prop("for lines through and beside spheres on dyadic coordinates", [&] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const double radius = *rc::gen::inRange(1, 33) / 8.0;
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(direction != Vec3{});
        const auto& [x, y] = frames[static_cast<std::size_t>(*rc::gen::inRange(0, 6))];
        // Through a point within 1.5 radii of the centre along each axis:
        // some lines cross and some miss.
        const auto offset = [&] { return *rc::gen::inRange(-12, 13) / 8.0 * radius; };
        const Line line{center + Vec3{offset(), offset(), offset()} - 3.0 * direction, direction};
        const std::size_t points = same(line, Sphere{center, radius, x, y},
                                        Ellipsoid{center, radius, radius, radius, x, y});
        RC_CLASSIFY(points == 2, "crossing");
        RC_CLASSIFY(points == 0, "apart");
    });
    rc::prop("for lines touching the sphere at the end of an axis", [&] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const double radius = *rc::gen::inRange(1, 33) / 8.0;
        const auto& [x, y] = frames[static_cast<std::size_t>(*rc::gen::inRange(0, 6))];
        // Along x, a radius above the centre: tangent at the top.
        const Line line{center + Vec3{coordinate(), 0.0, radius}, Vec3{1.0, 0.0, 0.0}};
        const std::size_t points = same(line, Sphere{center, radius, x, y},
                                        Ellipsoid{center, radius, radius, radius, x, y});
        RC_ASSERT(points == 1U);
    });
}

TEST_CASE("scaling a sphere's configuration along an axis gives the ellipsoid's results",
          "[intersect][line-surface][ellipsoid][property]")
{
    // ADR-0048: with axis-aligned inputs, scaling by 2^k along x maps the
    // unit sphere onto an ellipsoid and a line onto a line, parameter for
    // parameter. Every coefficient of the quadratic scales by 4^k, exactly,
    // so the count and the contact are the sphere's; the parameters and the
    // points agree to rounding, the two handlers rounding their enclosures
    // through different operations.
    rc::prop("by 2^k, |k| <= 3, on dyadic coordinates", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(direction != Vec3{});
        const auto offset = [] { return *rc::gen::inRange(-12, 13) / 8.0; };
        const Line line{center + Vec3{offset(), offset(), offset()} - 3.0 * direction, direction};
        const double factor = std::ldexp(1.0, *rc::gen::inRange(-3, 4));
        const auto scaled = [&](const Point3& p) { return Point3{factor * p.X, p.Y, p.Z}; };
        const Sphere sphere{center, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const Ellipsoid ellipsoid{scaled(center),     factor, 1.0, 1.0, Vec3{1.0, 0.0, 0.0},
                                  Vec3{0.0, 1.0, 0.0}};
        const Line image{scaled(line.Origin), Vec3{factor * direction.X, direction.Y, direction.Z}};
        const auto onSphere = Meet(line, sphere);
        const auto onEllipsoid = Meet(image, ellipsoid);
        RC_ASSERT(onSphere.has_value());
        RC_ASSERT(onEllipsoid.has_value());
        RC_ASSERT(onSphere->Points.size() == onEllipsoid->Points.size());
        for (std::size_t i = 0; i < onSphere->Points.size(); ++i) {
            const CurveSurfacePoint& p = onSphere->Points[i];
            const CurveSurfacePoint& q = onEllipsoid->Points[i];
            RC_ASSERT(q.Contact == p.Contact);
            RC_ASSERT(q.Decided == p.Decided);
            RC_ASSERT(std::abs(q.Parameter - p.Parameter) <= 1e-13 * (1.0 + std::abs(p.Parameter)));
            RC_ASSERT(Distance(q.Position, scaled(p.Position)) <=
                      1e-13 * (1.0 + Length(Stereon::ToVec(q.Position))));
        }
    });
}

TEST_CASE("a line through the inside of an ellipsoid crosses it twice, on both geometries",
          "[intersect][line-surface][ellipsoid][property]")
{
    rc::prop("for ellipsoids placed and turned at random", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const auto semi = [] { return *rc::gen::inRange(2, 33) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const double a = semi();
        const double b = semi();
        const double c = semi();
        const Vec3 axis{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(axis != Vec3{} && direction != Vec3{});
        const Ellipsoid ellipsoid = EllipsoidAbout(center, axis, a, b, c);
        // Through a point within half the shortest semi-axis of the centre
        // along each axis of space: inside the sphere the ellipsoid holds.
        const double least = std::min({a, b, c});
        const auto offset = [&] { return *rc::gen::inRange(-4, 5) / 8.0 * least; };
        const Line line{center + Vec3{offset(), offset(), offset()} - 3.0 * direction, direction};
        const double within = 1e-9 * (1.0 + Length(Stereon::ToVec(center)) + std::max({a, b, c}));
        const auto result = Meet(line, ellipsoid);
        RC_ASSERT(result.has_value());
        RC_ASSERT(result->Points.size() == 2U);
        RC_ASSERT(result->Overlaps.empty());
        RC_ASSERT(result->Points[0].Parameter < result->Points[1].Parameter);
        for (const CurveSurfacePoint& p : result->Points) {
            RC_ASSERT(p.Contact == ContactType::Crossing);
            RC_ASSERT(p.Decided);
            RC_ASSERT(Distance(Eval(line, p.Parameter), p.Position) <= within);
            RC_ASSERT(Distance(Eval(ellipsoid, p.OnSurface), p.Position) <= within);
            RC_ASSERT(std::abs(Level(ellipsoid, p.Position)) <= 1e-9);
        }
    });
}

TEST_CASE("a line's two crossings of an ellipsoid stay apart however far its origin",
          "[intersect][line-surface][ellipsoid]")
{
    // The sphere's far-origin case (the audit of 2026-11-05) on the
    // ellipsoid: through (-0.4, 0.3, 0) inside it, the origin 13 L back.
    // The chord is the one the line has from that point itself, and each
    // crossing lies on both geometries to the rounding of the far origin.
    const Vec3 direction{3.0, 4.0, 12.0};
    const auto near = Meet(Line{Point3{-0.4, 0.3, 0.0}, direction}, stretched);
    REQUIRE(near.has_value());
    REQUIRE(near->Points.size() == 2);
    const double chord = Distance(near->Points[0].Position, near->Points[1].Position);
    for (double length = 1e2; length <= 1e12; length *= 100.0) {
        INFO("from " << 13.0 * length);
        const Line line{Point3{-3.0 * length - 0.4, -4.0 * length + 0.3, -12.0 * length},
                        direction};
        const double within = 64.0 * std::numeric_limits<double>::epsilon() * 13.0 * length;
        const auto result = Intersect(Curve{line}, Surface{stretched}, ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Points.size() == 2);
        for (const CurveSurfacePoint& p : result->Points) {
            REQUIRE(p.Contact == ContactType::Crossing);
            REQUIRE(Distance(Eval(line, p.Parameter), p.Position) <= within);
            REQUIRE(Distance(Eval(stretched, p.OnSurface), p.Position) <= within);
        }
        REQUIRE(std::abs(Distance(result->Points[0].Position, result->Points[1].Position) -
                         chord) <= 2.0 * within);
    }
}
