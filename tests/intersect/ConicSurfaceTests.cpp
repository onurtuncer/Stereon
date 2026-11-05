// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Circles and ellipses against the five elementary surfaces, and B-spline
// curves against the cone and the torus. As for the other curve-surface
// pairs, a successful empty result is a proof of separation: 2^-40 above a
// circle's top it is apart, 2^-40 below it crosses twice, and a conic on
// the surface fails with IntersectionFailed instead of coming back empty.
// Every point lies on both geometries at its parameters, the conic's
// parameter being its angle, and every sign change of the surface's
// implicit equation along the curve has a point in it.

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
using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

const Plane ground{Point3{}, ex, ey};
const Sphere unitSphere{Point3{}, 1.0, ex, ey};
/// The double cone x² + y² = z² about z, its apex at the origin.
const Cone rightCone{Point3{}, 0.0, pi / 4.0, ex, ey};
/// The torus about z with radii 2 and 1.
const Torus ring{Point3{}, 2.0, 1.0, ex, ey};

double ConeLevel(const Point3& p)
{
    return p.X * p.X + p.Y * p.Y - p.Z * p.Z;
}

double TorusLevel(const Point3& p)
{
    const double k = p.X * p.X + p.Y * p.Y + p.Z * p.Z + 3.0;
    return k * k - 16.0 * (p.X * p.X + p.Y * p.Y);
}

template <typename C, typename S>
Result<CurveSurfaceIntersection> Meet(const C& curve, const S& surface)
{
    return Intersect(Curve{curve}, Surface{surface}, ctx);
}

/// Every point on both geometries, in order of the curve's parameter, and
/// a point between any two of 512 samples where `level` changes sign.
template <typename C, typename S, typename Level>
void RequireComplete(const CurveSurfaceIntersection& result, const C& curve, const S& surface,
                     const Interval& domain, Level&& level)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const CurveSurfacePoint& p : result.Points) {
        REQUIRE(p.Decided);
        REQUIRE(Distance(Eval(curve, p.Parameter), p.Position) <= 1e-12);
        REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= 1e-12);
        REQUIRE(p.Parameter >= previous);
        previous = p.Parameter;
    }
    const std::vector<double> t = IntersectTest::SamplesOf(domain, 512);
    for (std::size_t i = 0; i + 1 < t.size(); ++i) {
        const double a = level(Eval(curve, t[i]));
        const double b = level(Eval(curve, t[i + 1]));
        if ((a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0)) {
            bool seen = false;
            for (const CurveSurfacePoint& p : result.Points) {
                seen = seen || (t[i] <= p.Parameter && p.Parameter <= t[i + 1]);
            }
            REQUIRE(seen);
        }
    }
}

} // namespace

TEST_CASE("a circle crosses a plane at its two angles", "[intersect][conic-surface]")
{
    // Radius 2 in the xz-plane against z = 1: sin t = 1/2.
    const Circle upright{Point3{}, 2.0, ex, ez};
    const Plane high{Point3{0.0, 0.0, 1.0}, ex, ey};
    const auto result = Meet(upright, high);
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    REQUIRE(result->Points.size() == 2);
    REQUIRE(std::abs(result->Points[0].Parameter - pi / 6.0) <= 1e-14);
    REQUIRE(std::abs(result->Points[1].Parameter - 5.0 * pi / 6.0) <= 1e-14);
    for (const CurveSurfacePoint& p : result->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
    }
    RequireComplete(*result, upright, high, Domain(upright),
                    [](const Point3& p) { return p.Z - 1.0; });

    // Over a quarter turn only the first; over two turns, each twice.
    const auto quarter =
        Intersect(Curve{upright}, Interval{0.0, pi / 2.0}, Surface{high}, Domain(high), ctx);
    REQUIRE(quarter.has_value());
    REQUIRE(quarter->Points.size() == 1);
    const auto twice =
        Intersect(Curve{upright}, Interval{0.0, 4.0 * pi}, Surface{high}, Domain(high), ctx);
    REQUIRE(twice.has_value());
    REQUIRE(twice->Points.size() == 4);
    REQUIRE(std::abs(twice->Points[3].Parameter - (2.0 * pi + 5.0 * pi / 6.0)) <= 1e-13);
}

TEST_CASE("an empty conic result is a proof: 2^-40 above the top is apart, below crosses twice",
          "[intersect][conic-surface][regression]")
{
    const Circle upright{Point3{}, 2.0, ex, ez};
    const double gap = std::ldexp(1.0, -40);
    const Plane above{Point3{0.0, 0.0, 2.0 + gap}, ex, ey};
    const auto apart = Meet(upright, above);
    REQUIRE(apart.has_value());
    REQUIRE(apart->Points.empty());
    REQUIRE(apart->Overlaps.empty());

    const Plane below{Point3{0.0, 0.0, 2.0 - gap}, ex, ey};
    const auto crossing = Meet(upright, below);
    REQUIRE(crossing.has_value());
    REQUIRE(crossing->Points.size() == 2);
    for (const CurveSurfacePoint& p : crossing->Points) {
        REQUIRE(p.Decided);
        REQUIRE(std::abs(p.Parameter - pi / 2.0) <= 1e-5);
    }

    // At the top exactly the circle touches: a decided touch or a failure,
    // never an empty success.
    const Plane top{Point3{0.0, 0.0, 2.0}, ex, ey};
    const auto touching = Meet(upright, top);
    if (touching.has_value()) {
        REQUIRE_FALSE(touching->Points.empty());
    } else {
        REQUIRE(touching.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}

TEST_CASE("a circle in the plane is an overlap whose pcurve traces it",
          "[intersect][conic-surface]")
{
    const Circle flat{Point3{1.0, 2.0, 0.0}, 1.5, ey, -ex};
    const auto result = Meet(flat, ground);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.empty());
    REQUIRE(result->Overlaps.size() == 1);
    const Stereon::CurveSurfaceOverlap& overlap = result->Overlaps.front();
    REQUIRE(overlap.Range == Domain(flat));
    for (const double t : IntersectTest::SamplesOf(overlap.Range, 16)) {
        REQUIRE(Distance(Eval(ground, Eval(overlap.OnSurface, t)), Eval(flat, t)) <= 1e-12);
    }

    // An ellipse likewise, in a plane whose parameters are not the world's.
    const Plane tilted{Point3{0.0, 0.0, 1.0}, ey, ez};
    const Ellipse upright{Point3{0.0, 0.5, 1.0}, 2.0, 0.5, ez, ey};
    const auto along = Meet(upright, tilted);
    REQUIRE(along.has_value());
    REQUIRE(along->Overlaps.size() == 1);
    for (const double t : IntersectTest::SamplesOf(along->Overlaps.front().Range, 16)) {
        REQUIRE(Distance(Eval(tilted, Eval(along->Overlaps.front().OnSurface, t)),
                         Eval(upright, t)) <= 1e-12);
    }

    // Leaving the requested box, it is not clipped yet.
    const auto clipped = Intersect(Curve{flat}, Domain(flat), Surface{ground},
                                   Stereon::Box2{Interval{0.0, 1.0}, Interval{0.0, 1.0}}, ctx);
    REQUIRE_FALSE(clipped.has_value());
    REQUIRE(clipped.error().GetCode() == ErrorCode::NotImplemented);
}

TEST_CASE("an ellipse crosses a plane, a sphere and a cylinder", "[intersect][conic-surface]")
{
    // Semi-axes 3 along x and 1/2 along y.
    const Ellipse flat{Point3{}, 3.0, 0.5, ex, ey};
    const Plane across = IntersectTest::PlaneWithNormal(Point3{1.5, 0.0, 0.0}, ex);
    const auto onPlane = Meet(flat, across);
    REQUIRE(onPlane.has_value());
    REQUIRE(onPlane->Points.size() == 2);
    REQUIRE(std::abs(onPlane->Points[0].Parameter - pi / 3.0) <= 1e-14);
    REQUIRE(std::abs(onPlane->Points[1].Parameter - 5.0 * pi / 3.0) <= 1e-14);
    RequireComplete(*onPlane, flat, across, Domain(flat),
                    [](const Point3& p) { return p.X - 1.5; });

    // 9 cos² t + sin² t / 4 = 1 four times a turn.
    const auto onSphere = Meet(flat, unitSphere);
    REQUIRE(onSphere.has_value());
    REQUIRE(onSphere->Points.size() == 4);
    const double c = std::sqrt(0.75 / 8.75);
    REQUIRE(std::abs(onSphere->Points[0].Parameter - std::acos(c)) <= 1e-13);
    RequireComplete(*onSphere, flat, unitSphere, Domain(flat),
                    [](const Point3& p) { return p.X * p.X + p.Y * p.Y + p.Z * p.Z - 1.0; });

    // About the y axis: x² + z² = 9 cos² t = 1 four times a turn.
    const Cylinder aboutY = IntersectTest::CylinderAbout(Point3{}, ey, 1.0);
    const auto onCylinder = Meet(flat, aboutY);
    REQUIRE(onCylinder.has_value());
    REQUIRE(onCylinder->Points.size() == 4);
    REQUIRE(std::abs(onCylinder->Points[0].Parameter - std::acos(1.0 / 3.0)) <= 1e-13);
    for (const CurveSurfacePoint& p : onCylinder->Points) {
        REQUIRE(IntersectTest::DistanceFrom(aboutY, p.Position) <= 1e-12);
    }
    RequireComplete(*onCylinder, flat, aboutY, Domain(flat),
                    [](const Point3& p) { return p.X * p.X + p.Z * p.Z - 1.0; });
}

TEST_CASE("a circle crosses a cone and a torus", "[intersect][conic-surface]")
{
    // In the xz-plane about (1/2, 0, 1), radius 1: on the cone where
    // 1/2 + cos t = 1 + sin t, cos t - sin t = 1/2, twice; the nappe below
    // would need cos t + sin t = -3/2, which no angle reaches.
    const Circle tilted{Point3{0.5, 0.0, 1.0}, 1.0, ex, ez};
    const auto onCone = Meet(tilted, rightCone);
    REQUIRE(onCone.has_value());
    REQUIRE(onCone->Points.size() == 2);
    for (const CurveSurfacePoint& p : onCone->Points) {
        const double t = p.Parameter;
        REQUIRE(std::abs(std::cos(t) - std::sin(t) - 0.5) <= 1e-14);
    }
    RequireComplete(*onCone, tilted, rightCone, Domain(tilted), ConeLevel);

    // A circle about the apex, through both nappes: (cos t, 0, sin t)
    // meets x = |z| at t = pi/4, 3pi/4, 5pi/4, 7pi/4.
    const Circle throughApex{Point3{}, 1.0, ex, ez};
    const auto bothNappes = Meet(throughApex, rightCone);
    REQUIRE(bothNappes.has_value());
    REQUIRE(bothNappes->Points.size() == 4);
    RequireComplete(*bothNappes, throughApex, rightCone, Domain(throughApex), ConeLevel);

    // In the xz-plane about (3, 0, 0), radius 1: on the tube where
    // (1 + cos t)² + sin² t = 1, cos t = -1/2.
    const Circle meridian{Point3{3.0, 0.0, 0.0}, 1.0, ex, ez};
    const auto onTorus = Meet(meridian, ring);
    REQUIRE(onTorus.has_value());
    REQUIRE(onTorus->Points.size() == 2);
    REQUIRE(std::abs(onTorus->Points[0].Parameter - 2.0 * pi / 3.0) <= 1e-13);
    REQUIRE(std::abs(onTorus->Points[1].Parameter - 4.0 * pi / 3.0) <= 1e-13);
    RequireComplete(*onTorus, meridian, ring, Domain(meridian), TorusLevel);

    // Horizontal at height 1/2, radius 2.5 about the axis: crossing no tube.
    const Circle inside{Point3{0.0, 0.0, 0.5}, 2.5, ex, ey};
    const auto apart = Meet(inside, ring);
    REQUIRE(apart.has_value());
    REQUIRE(apart->Points.empty());
}

TEST_CASE("a conic on the surface fails rather than coming back empty",
          "[intersect][conic-surface]")
{
    // The tube's own meridian and the cone's parallel at height 2.
    const Circle meridian{Point3{2.0, 0.0, 0.0}, 1.0, ex, ez};
    const auto onTorus = Meet(meridian, ring);
    REQUIRE_FALSE(onTorus.has_value());
    REQUIRE(onTorus.error().GetCode() == ErrorCode::IntersectionFailed);

    const Circle parallel{Point3{0.0, 0.0, 2.0}, 2.0, ex, ey};
    const auto onCone = Meet(parallel, rightCone);
    REQUIRE_FALSE(onCone.has_value());
    REQUIRE(onCone.error().GetCode() == ErrorCode::IntersectionFailed);

    const Circle equator{Point3{}, 1.0, ex, ey};
    const auto onSphere = Meet(equator, unitSphere);
    REQUIRE_FALSE(onSphere.has_value());
    REQUIRE(onSphere.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a B-spline crosses a cone and a torus", "[intersect][spline-surface]")
{
    // A cubic from inside the cone's upper nappe out through it, down
    // through the apex region and out through the lower one.
    auto made = Stereon::MakeBSplineCurve(
        3, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
        {{0.25, 0.0, 2.0}, {3.0, 0.5, 1.0}, {-1.0, -0.5, -1.0}, {0.25, 0.25, -2.0}});
    REQUIRE(made.has_value());
    const BSplineCurve& spline = *made;
    const auto onCone = Meet(spline, rightCone);
    REQUIRE(onCone.has_value());
    REQUIRE_FALSE(onCone->Points.empty());
    RequireComplete(*onCone, spline, rightCone, Domain(spline), ConeLevel);

    // Along the x axis from -4 to 4, bent upward: through the tube four
    // times.
    auto bent = Stereon::MakeBSplineCurve(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                                          {{-4.0, 0.0, 0.0}, {0.0, 0.25, 0.5}, {4.0, 0.0, 0.0}});
    REQUIRE(bent.has_value());
    const auto onTorus = Meet(*bent, ring);
    REQUIRE(onTorus.has_value());
    REQUIRE(onTorus->Points.size() == 4);
    RequireComplete(*onTorus, *bent, ring, Domain(*bent), TorusLevel);

    // Wholly inside the tube: apart, proved.
    auto within = Stereon::MakeBSplineCurve(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                                            {{2.0, 0.0, 0.0}, {2.0, 1.0, 0.25}, {1.5, 1.5, 0.0}});
    REQUIRE(within.has_value());
    const auto apart = Meet(*within, ring);
    REQUIRE(apart.has_value());
    REQUIRE(apart->Points.empty());
}

TEST_CASE("every point lies on the surface to the rounding of its coordinates",
          "[intersect][spline-surface][property][regression]")
{
    // The subdivision places a root to the rounding of the enclosed form,
    // which for a torus's quartic left a point 1.2e-13 off it, hundreds of
    // units in the last place of its coordinates, and IntersectJvp then
    // took Intersect's own point for a stale one (RapidCheck, 2026-11-04).
    rc::prop("for random cubic B-splines against a torus, a cone and a sphere", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-24, 25) / 8.0; };
        std::vector<Point3> poles;
        for (int i = 0; i < 5; ++i) {
            poles.push_back({coordinate(), coordinate(), coordinate()});
        }
        const auto made =
            Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0}, poles, {});
        RC_PRE(made.has_value());
        const auto check = [&](const auto& surface) {
            // Completeness and failing closed are tested elsewhere.
            const auto result = Intersect(Curve{*made}, Surface{surface}, ctx);
            if (!result) {
                return;
            }
            for (const CurveSurfacePoint& p : result->Points) {
                const double ulps = 64.0 * std::numeric_limits<double>::epsilon() *
                                    std::max(1.0, Length(Stereon::ToVec(p.Position)));
                RC_ASSERT(Distance(Eval(*made, p.Parameter), p.Position) <= ulps);
                RC_ASSERT(Distance(Eval(surface, p.OnSurface), p.Position) <= ulps);
            }
        };
        check(ring);
        check(Cone{Point3{}, 1.0, pi / 4.0, ex, ey});
        check(Sphere{Point3{0.25, 0.0, 0.0}, 1.5, ex, ey});
    });
}
