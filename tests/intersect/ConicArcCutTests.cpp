// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A circle's or an ellipse's crossing of a surface exactly on one of the
// cuts the kernel makes in its domain, every quarter turn from the
// domain's lower end (v0.2 ledger, curve-surface). Each arc's search
// brackets such a root across the arc's end; the handler searches it again
// on an arc that straddles the cut, where it is interior. The expected
// answers are closed forms on the exact geometry: the circle of radius 2
// in the xz-plane, (2 cos t, 0, 2 sin t), meets the plane x = 0 and the
// cylinder (x - 2)^2 + y^2 = 4 at t = pi/2 and 3pi/2, and the sphere about
// (2, 0, 2) of radius 2 where cos t + sin t = 1, t = 0 and pi/2; the
// ellipse (3 cos t, 0, sin t) meets x = 0 at pi/2 and 3pi/2. Metamorphic
// relations (ADR-0053): reversal (t -> 2pi - t), translation, scaling by
// powers of two.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::Circle;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// The points' parameters, required to succeed with decided crossings on
/// both geometries.
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

/// The angles, each to the rounding of the coordinates relative to the
/// conic's size: `far` is the distance from the origin over the radius.
void RequireAt(const std::vector<double>& t, const std::vector<double>& expected,
               const double far = 0.0)
{
    REQUIRE(t.size() == expected.size());
    for (std::size_t i = 0; i < t.size(); ++i) {
        REQUIRE_THAT(t[i], WithinAbs(expected[i], 1e-13 * (1.0 + far)));
    }
}

} // namespace

TEST_CASE("a circle crossing surfaces exactly on its quarter-turn cuts is decided",
          "[intersect][conic-surface][arc-cut]")
{
    for (const Point3 at : {Point3{}, Point3{1024.0, -512.0, 256.0}}) {
        for (const double k : {1.0, 0x1p10, 0x1p-10}) {
            DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ", " << at.Z << ") scaled by " << k)
            {
                const Circle circle{at, 2.0 * k, ex, ez};
                const Interval turn{0.0, 2.0 * pi};
                const double far = Distance(at, Point3{}) / k;
                // The plane x = 0 and the cylinder about z through x = 2.
                const Surface across{Plane{at, ey, ez}};
                const Surface tube{
                    IntersectTest::CylinderAbout(Point3{at.X + 2.0 * k, at.Y, at.Z}, ez, 2.0 * k)};
                for (const Surface& s : {across, tube}) {
                    RequireAt(Crossings(Curve{circle}, turn, s), {pi / 2.0, 3.0 * pi / 2.0}, far);
                    // MR-04: reversed, t -> 2pi - t, the same two angles.
                    RequireAt(Crossings(Curve{Reversed(circle)}, turn, s),
                              {pi / 2.0, 3.0 * pi / 2.0}, far);
                }
                // The sphere about (2, 0, 2), radius 2: t = 0 and pi/2, both
                // cuts of the domain [-pi/2, 3pi/2].
                const Surface ball{
                    Sphere{Point3{at.X + 2.0 * k, at.Y, at.Z + 2.0 * k}, 2.0 * k, ex, ey}};
                RequireAt(Crossings(Curve{circle}, Interval{-pi / 2.0, 3.0 * pi / 2.0}, ball),
                          {0.0, pi / 2.0}, far);
            }
        }
    }
}

TEST_CASE("an ellipse crossing a plane exactly on its quarter-turn cuts is decided",
          "[intersect][conic-surface][arc-cut]")
{
    for (const Point3 at : {Point3{}, Point3{1024.0, -512.0, 256.0}}) {
        DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ", " << at.Z << ")")
        {
            const Ellipse ellipse{at, 3.0, 1.0, ex, ez};
            const Surface across{Plane{at, ey, ez}};
            RequireAt(Crossings(Curve{ellipse}, Interval{0.0, 2.0 * pi}, across),
                      {pi / 2.0, 3.0 * pi / 2.0});
            RequireAt(Crossings(Curve{Reversed(ellipse)}, Interval{0.0, 2.0 * pi}, across),
                      {pi / 2.0, 3.0 * pi / 2.0});
        }
    }
}
