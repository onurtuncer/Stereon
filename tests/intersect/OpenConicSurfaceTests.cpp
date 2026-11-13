// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Hyperbolas and parabolas against the elementary surfaces (v0.2 ledger,
// curve-surface): the same proofs as circles', on the hyperbola's exact
// rational form in stretches of one unit of t and on the parabola's
// polynomial form, over a bounded requested domain; an unbounded one is
// NotSupported. The expected answers are closed forms. The hyperbola
// (cosh t, sinh t, 0) meets x = 2 at t = +-acosh 2, y = 0 at t = 0 (one of
// its cuts), and the sphere and the cylinder about z of radius 3 where
// cosh 2t = 9; the parabola (t^2, t, 0) meets x = 4 at t = +-2, y = 1 at
// t = 1 and the cylinder y^2 + z^2 = 1 about x at t = +-1. Metamorphic
// relations (ADR-0053): reversal (the y axis negated, t -> -t), translation
// and scaling by powers of two (a parabola's t is a length and scales).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Parabola;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// The points' parameters, required to succeed with decided crossings on
/// both geometries, in order.
std::vector<double> Crossings(const Curve& curve, const Interval& domain, const Surface& surface)
{
    const Result<CurveSurfaceIntersection> result =
        Intersect(curve, domain, surface, Domain(surface), ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    std::vector<double> t;
    for (const CurveSurfacePoint& p : result->Points) {
        REQUIRE(p.Decided);
        REQUIRE(p.Contact == ContactType::Crossing);
        const double scale = 1.0 + Distance(p.Position, Point3{});
        REQUIRE(Distance(Eval(curve, p.Parameter), p.Position) <= 1e-12 * scale);
        REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= 1e-12 * scale);
        t.push_back(p.Parameter);
    }
    return t;
}

void RequireAt(std::vector<double> t, std::vector<double> expected, const double within)
{
    std::ranges::sort(t);
    std::ranges::sort(expected);
    REQUIRE(t.size() == expected.size());
    for (std::size_t i = 0; i < t.size(); ++i) {
        REQUIRE_THAT(t[i], WithinAbs(expected[i], within));
    }
}

Point3 Moved(const Point3& at, const double k, const Vec3& v)
{
    return Point3{at.X + k * v.X, at.Y + k * v.Y, at.Z + k * v.Z};
}

} // namespace

TEST_CASE("a hyperbola crosses planes, a sphere and a cylinder where cosh and sinh say",
          "[intersect][conic-surface][open-conic]")
{
    const double t2 = std::acosh(2.0);
    const double t9 = 0.5 * std::acosh(9.0);
    for (const Point3 at : {Point3{}, Point3{1024.0, -512.0, 256.0}}) {
        for (const double k : {1.0, 0x1p10, 0x1p-10}) {
            for (const bool reversed : {false, true}) {
                DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ", " << at.Z << ") scaled by "
                                       << k << (reversed ? ", reversed" : ""))
                {
                    // MR-04: the y axis negated is t -> -t; the answers are
                    // symmetric but for y = 0's single one, at 0 either way.
                    const Curve h{Hyperbola{at, k, k, ex, reversed ? -ey : ey}};
                    const Interval domain{-2.0, 2.0};
                    const double within = 1e-12 * (1.0 + Distance(at, Point3{}) / k);
                    RequireAt(Crossings(h, domain, Surface{Plane{Moved(at, k, 2.0 * ex), ey, ez}}),
                              {-t2, t2}, within);
                    RequireAt(Crossings(h, domain, Surface{Plane{at, ez, ex}}), {0.0}, within);
                    RequireAt(Crossings(h, domain, Surface{Sphere{at, 3.0 * k, ex, ey}}), {-t9, t9},
                              within);
                    RequireAt(Crossings(h, domain,
                                        Surface{IntersectTest::CylinderAbout(at, ez, 3.0 * k)}),
                              {-t9, t9}, within);
                }
            }
        }
    }
}

TEST_CASE("a parabola crosses planes and a cylinder where its square says",
          "[intersect][conic-surface][open-conic]")
{
    for (const Point3 at : {Point3{}, Point3{1024.0, -512.0, 256.0}}) {
        for (const double k : {1.0, 0x1p10, 0x1p-10}) {
            for (const bool reversed : {false, true}) {
                DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ", " << at.Z << ") scaled by "
                                       << k << (reversed ? ", reversed" : ""))
                {
                    // (t^2 / k, t, 0) about `at`: f = k / 4, t a length.
                    const double sign = reversed ? -1.0 : 1.0;
                    const Curve p{Parabola{at, 0.25 * k, ex, reversed ? -ey : ey}};
                    const Interval domain{-4.0 * k, 4.0 * k};
                    const double within = 1e-12 * k * (1.0 + Distance(at, Point3{}) / k);
                    RequireAt(Crossings(p, domain, Surface{Plane{Moved(at, k, 4.0 * ex), ey, ez}}),
                              {-2.0 * k, 2.0 * k}, within);
                    RequireAt(Crossings(p, domain, Surface{Plane{Moved(at, k, ey), ez, ex}}),
                              {sign * k}, within);
                    RequireAt(
                        Crossings(p, domain, Surface{IntersectTest::CylinderAbout(at, ex, k)}),
                        {-k, k}, within);
                }
            }
        }
    }
}

TEST_CASE("a hyperbola or a parabola over an unbounded domain is NotSupported, not clipped",
          "[intersect][conic-surface][open-conic]")
{
    const Surface plane{Plane{Point3{2.0, 0.0, 0.0}, ey, ez}};
    for (const Curve& c :
         {Curve{Hyperbola{Point3{}, 1.0, 1.0, ex, ey}}, Curve{Parabola{Point3{}, 0.25, ex, ey}}}) {
        const auto whole = Intersect(c, Domain(c), plane, Domain(plane), ctx);
        REQUIRE_FALSE(whole.has_value());
        REQUIRE(whole.error().GetCode() == ErrorCode::NotSupported);
    }
}

TEST_CASE("a circle crossing a plane at an arc's middle, a dyadic point of its subdivision, is "
          "decided",
          "[intersect][conic-surface][open-conic]")
{
    // The circle (2 cos t, 0, 2 sin t) against the plane x = z, its normal
    // (1, 0, -1) / sqrt 2: t = pi/4 and 5pi/4, the middles of the arcs
    // [0, pi/2] and [pi, 3pi/2].
    const double r = std::sqrt(0.5);
    const Curve circle{Stereon::Circle{Point3{}, 2.0, ex, ez}};
    const Surface diagonal{Plane{Point3{}, ey, Vec3{r, 0.0, r}}};
    const double pi = std::acos(-1.0);
    RequireAt(Crossings(circle, Interval{0.0, 2.0 * pi}, diagonal), {pi / 4.0, 5.0 * pi / 4.0},
              1e-13);
}
