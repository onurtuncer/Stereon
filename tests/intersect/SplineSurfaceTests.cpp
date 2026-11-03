// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// B-spline curves against planes, spheres and cylinders. The safety
// property under test: a successful empty result is a proof of separation,
// so a curve 2^-40 from a plane is apart only when it is, a curve that
// cannot be resolved fails with IntersectionFailed instead of coming back
// empty, and every sign change of the surface's implicit equation along a
// random curve has a reported point in it. Every point lies on both
// geometries at its parameters, in order of the curve's parameter.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::BSplineCurve;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Sphere unitSphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Cylinder unitCylinder{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

BSplineCurve Spline(const int degree, std::vector<double> knots, std::vector<Point3> poles,
                    std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// A Bezier curve of degree poles - 1 on [0, 1].
BSplineCurve Bezier(std::vector<Point3> poles)
{
    const std::size_t n = poles.size();
    std::vector<double> knots(n, 0.0);
    knots.resize(2 * n, 1.0);
    return Spline(static_cast<int>(n) - 1, std::move(knots), std::move(poles));
}

template <typename S>
Result<CurveSurfaceIntersection> Meet(const BSplineCurve& spline, const S& surface)
{
    return Intersect(Curve{spline}, Surface{surface}, ctx);
}

template <typename S>
void RequireOnBoth(const CurveSurfaceIntersection& result, const BSplineCurve& spline,
                   const S& surface, const double within = 1e-9)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const CurveSurfacePoint& p : result.Points) {
        REQUIRE(Distance(Eval(spline, p.Parameter), p.Position) <= within);
        REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= within);
        REQUIRE(p.Parameter >= previous);
        previous = p.Parameter;
    }
}

/// The plane's signed distance and the quadrics' implicit equations at p.
double Height(const Point3& p)
{
    return p.Z;
}

double SphereLevel(const Point3& p)
{
    return Dot(Stereon::ToVec(p), Stereon::ToVec(p)) - 1.0;
}

} // namespace

TEST_CASE("a B-spline crosses a plane at each simple root", "[intersect][spline-surface]")
{
    // z(s) = -1 + 8 s - 8 s², zero at s = 1/2 -+ sqrt(2) / 4.
    const BSplineCurve arch = Bezier({{0.0, 0.0, -1.0}, {1.0, 0.0, 3.0}, {2.0, 0.0, -1.0}});
    const auto result = Meet(arch, ground);
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    REQUIRE(result->Points.size() == 2);
    const double offset = std::numbers::sqrt2 / 4.0;
    REQUIRE(std::abs(result->Points[0].Parameter - (0.5 - offset)) <= 1e-15);
    REQUIRE(std::abs(result->Points[1].Parameter - (0.5 + offset)) <= 1e-15);
    for (const CurveSurfacePoint& p : result->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
        REQUIRE(p.Decided);
    }
    RequireOnBoth(*result, arch, ground);

    // Over the first half only the first crossing is in the domain.
    const auto half =
        Intersect(Curve{arch}, Interval{0.0, 0.5}, Surface{ground}, Domain(ground), ctx);
    REQUIRE(half.has_value());
    REQUIRE(half->Points.size() == 1);
    REQUIRE(half->Points.front().Parameter < 0.5);
}

TEST_CASE("unresolved spline tangency and near miss never become successful contacts",
          "[intersect][spline-surface][regression]")
{
    const BSplineCurve touch = Bezier({{0.0, 0.0, 1.0}, {1.0, 0.0, -1.0}, {2.0, 0.0, 1.0}});
    const auto result = Meet(touch, ground);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    const double eps = std::numeric_limits<double>::epsilon();
    const BSplineCurve nearMiss =
        Bezier({{0.0, 0.0, 1.0 + eps}, {1.0, 0.0, -1.0 + eps}, {2.0, 0.0, 1.0 + eps}});
    const auto missed = Meet(nearMiss, ground);
    REQUIRE(!missed.has_value() || missed->IsEmpty());
    if (!missed) {
        REQUIRE(missed.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}

TEST_CASE("an empty result is a proof: 3 2^-40 above is apart, 3 2^-40 below crosses twice",
          "[intersect][spline-surface]")
{
    // (2 s - 1)² + c has Bernstein coefficients 1 + c, -1 + c, 1 + c. Not
    // 2^-40, whose roots 1/2 -+ 2^-21 below are dyadic: a subdivision lands
    // on them exactly, and a root at a piece's end is undecided.
    constexpr double c = 3.0 * 0x1p-40;
    const BSplineCurve above =
        Bezier({{0.0, 0.0, 1.0 + c}, {1.0, 0.0, -1.0 + c}, {2.0, 0.0, 1.0 + c}});
    const auto apart = Meet(above, ground);
    REQUIRE(apart.has_value());
    REQUIRE(apart->IsEmpty());

    const BSplineCurve below =
        Bezier({{0.0, 0.0, 1.0 - c}, {1.0, 0.0, -1.0 - c}, {2.0, 0.0, 1.0 - c}});
    const auto twice = Meet(below, ground);
    REQUIRE(twice.has_value());
    REQUIRE(twice->Points.size() == 2);
    for (const CurveSurfacePoint& p : twice->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
        REQUIRE(p.Decided);
    }
    // The roots are 1/2 -+ sqrt(3) 2^-21, apart by far more than rounding.
    REQUIRE(twice->Points[1].Parameter - twice->Points[0].Parameter > 0x1p-21);
}

TEST_CASE("an undecided spline root at a knot fails closed",
          "[intersect][spline-surface][regression]")
{
    const std::vector<Point3> poles{
        {0.0, 0.0, -1.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 1.0}};
    const BSplineCurve polyline = Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, poles);
    const auto result = Meet(polyline, ground);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a B-spline in the plane is an overlap whose pcurve traces it",
          "[intersect][spline-surface]")
{
    const BSplineCurve flat =
        Bezier({{0.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {3.0, -1.0, 0.0}, {4.0, 1.0, 0.0}});
    const auto result = Meet(flat, ground);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.empty());
    REQUIRE(result->Overlaps.size() == 1);
    const Stereon::CurveSurfaceOverlap& overlap = result->Overlaps.front();
    REQUIRE(overlap.Range == Interval{0.0, 1.0});
    for (const double t : IntersectTest::SamplesOf(overlap.Range, 16)) {
        REQUIRE(Distance(Eval(ground, Eval(overlap.OnSurface, t)), Eval(flat, t)) <= 1e-12);
    }
}

TEST_CASE("a B-spline that cannot be resolved fails rather than coming back empty",
          "[intersect][spline-surface]")
{
    // Along a ruling of the cylinder: the implicit equation is zero
    // throughout, so no piece is ever proved apart, monotone or a touch.
    const BSplineCurve ruling = Bezier({{1.0, 0.0, 0.0}, {1.0, 0.0, 5.0}});
    const auto result = Meet(ruling, unitCylinder);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a rational B-spline crosses where the circle it draws does",
          "[intersect][spline-surface]")
{
    // A quarter of the circle of radius 2 about the origin, against the
    // plane x = y: one crossing at (sqrt 2, sqrt 2, 0).
    const double w = std::numbers::sqrt2 / 2.0;
    const BSplineCurve arc =
        Spline(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
               {{2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0}}, {1.0, w, 1.0});
    const Plane diagonal = IntersectTest::PlaneWithNormal(Point3{}, Vec3{1.0, -1.0, 0.0});
    const auto result = Meet(arc, diagonal);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 1);
    const CurveSurfacePoint& p = result->Points.front();
    REQUIRE(p.Decided);
    REQUIRE(std::abs(p.Position.X - std::numbers::sqrt2) <= 1e-12);
    REQUIRE(std::abs(p.Position.Y - std::numbers::sqrt2) <= 1e-12);
    RequireOnBoth(*result, arc, diagonal);
}

TEST_CASE("a B-spline through a sphere and a cylinder crosses each twice",
          "[intersect][spline-surface]")
{
    const BSplineCurve through =
        Bezier({{-3.0, 0.25, 0.125}, {-1.0, -0.5, 0.25}, {1.0, 0.5, -0.25}, {3.0, -0.25, 0.0}});
    const auto onSphere = Meet(through, unitSphere);
    REQUIRE(onSphere.has_value());
    REQUIRE(onSphere->Points.size() == 2);
    for (const CurveSurfacePoint& p : onSphere->Points) {
        REQUIRE(p.Decided);
        REQUIRE(std::abs(Length(Stereon::ToVec(p.Position)) - 1.0) <= 1e-12);
    }
    RequireOnBoth(*onSphere, through, unitSphere);

    // About the y axis instead, so that the curve runs across it.
    const Cylinder aboutY = IntersectTest::CylinderAbout(Point3{}, Vec3{0.0, 1.0, 0.0}, 1.0);
    const auto onCylinder = Meet(through, aboutY);
    REQUIRE(onCylinder.has_value());
    REQUIRE(onCylinder->Points.size() == 2);
    for (const CurveSurfacePoint& p : onCylinder->Points) {
        REQUIRE(p.Decided);
        REQUIRE(IntersectTest::DistanceFrom(aboutY, p.Position) <= 1e-12);
    }
    RequireOnBoth(*onCylinder, through, aboutY);
}

TEST_CASE("every sign change along a random B-spline has a point in it",
          "[intersect][spline-surface][property]")
{
    // Completeness against sampling: wherever the implicit equation changes
    // sign between two samples, a reported point lies between them; and a
    // successful empty result never has a sign change at all.
    // Signs and parameters are sampled in doubles: a sign counts beyond
    // `noise`, and a point may sit `slack` past its samples.
    constexpr double noise = 1e-9;
    constexpr double slack = 1e-9;
    const auto check = [](const BSplineCurve& spline, const auto& surface, const auto& level) {
        const auto result = Meet(spline, surface);
        if (!result) {
            // A segment exactly in the plane cannot be resolved: it fails,
            // and is never an empty success.
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        const std::vector<double> t = IntersectTest::SamplesOf(Domain(spline), 256);
        for (std::size_t i = 0; i + 1 < t.size(); ++i) {
            const double a = level(Eval(spline, t[i]));
            const double b = level(Eval(spline, t[i + 1]));
            if ((a < -noise && b > noise) || (a > noise && b < -noise)) {
                RC_ASSERT(std::ranges::any_of(result->Points, [&](const CurveSurfacePoint& p) {
                    return t[i] - slack <= p.Parameter && p.Parameter <= t[i + 1] + slack;
                }));
            }
        }
        for (const CurveSurfacePoint& p : result->Points) {
            RC_ASSERT(Distance(Eval(spline, p.Parameter), p.Position) <= 1e-9);
            RC_ASSERT(Distance(Eval(surface, p.OnSurface), p.Position) <= 1e-9);
        }
    };
    rc::prop("for random cubic B-splines against the plane z = 0 and the unit sphere", [&] {
        const auto coordinate = [] { return *rc::gen::inRange(-24, 25) / 8.0; };
        std::vector<Point3> poles;
        for (int i = 0; i < 6; ++i) {
            poles.push_back(Point3{coordinate(), coordinate(), coordinate()});
        }
        const BSplineCurve spline = *Stereon::MakeBSplineCurve(
            3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0}, poles);
        check(spline, ground, Height);
        check(spline, unitSphere, SphereLevel);
    });
}
