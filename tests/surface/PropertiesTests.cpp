// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local differential properties of surfaces against the closed forms of the
// plane, cylinder, sphere and torus, and their invariance under the exact
// NURBS conversion, which changes the parametrisation but not the surface.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::FundamentalFormsResult;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::PrincipalCurvaturesResult;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-10;
constexpr double Pi = std::numbers::pi;

bool Near(const double a, const double b)
{
    return std::abs(a - b) <= Tolerance * (1.0 + std::abs(b));
}

bool Parallel(const Vec3& a, const Vec3& b)
{
    return Length(Cross(a, b)) <= Tolerance;
}

/// The directions are unit, tangent, orthogonal and right-handed about the
/// normal, and each is a principal direction: the normal curvature
/// II(d) / I(d) along it is its curvature.
void CheckFrame(const PrincipalCurvaturesResult& k, const Vec3& normal)
{
    CHECK(Near(Length(k.MinDirection), 1.0));
    CHECK(Near(Length(k.MaxDirection), 1.0));
    CHECK(Near(Dot(k.MinDirection, normal), 0.0));
    CHECK(Near(Dot(k.MaxDirection, normal), 0.0));
    CHECK(Parallel(Cross(k.MinDirection, k.MaxDirection), normal));
    CHECK(Near(Dot(Cross(k.MinDirection, k.MaxDirection), normal), 1.0));
    CHECK(k.Min <= k.Max);
}

} // namespace

TEST_CASE("a plane is flat everywhere", "[surface][properties]")
{
    const Plane plane{Point3{1.0, 2.0, 3.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const Point2 uv{0.3, -2.0};
    CHECK(Parallel(*UnitNormal(plane, uv), Vec3{1.0, 0.0, 0.0}));
    const FundamentalFormsResult f = *FundamentalForms(plane, uv);
    CHECK((f.E == 1.0 && f.F == 0.0 && f.G == 1.0));
    CHECK((f.L == 0.0 && f.M == 0.0 && f.N == 0.0));
    CHECK(*GaussianCurvature(plane, uv) == 0.0);
    CHECK(*MeanCurvature(plane, uv) == 0.0);
}

TEST_CASE("a sphere is umbilic with curvature -1/R", "[surface][properties]")
{
    const Sphere sphere{Point3{1.0, 0.0, -1.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    for (const Point2 uv : {Point2{0.0, 0.0}, Point2{1.0, 0.7}, Point2{4.0, -1.2}}) {
        const Vec3 normal = *UnitNormal(sphere, uv);
        CHECK(Parallel(normal, Eval(sphere, uv) - sphere.Center));
        CHECK(Dot(normal, Eval(sphere, uv) - sphere.Center) > 0.0);
        CHECK(Near(*GaussianCurvature(sphere, uv), 0.25));
        CHECK(Near(*MeanCurvature(sphere, uv), -0.5));
        const PrincipalCurvaturesResult k = *PrincipalCurvatures(sphere, uv);
        CHECK(Near(k.Min, -0.5));
        CHECK(Near(k.Max, -0.5));
        CheckFrame(k, normal);
    }
    // A left-handed frame turns the normal inward and the curvatures round.
    Sphere inward = sphere;
    inward.RightHanded = false;
    CHECK(Near(*MeanCurvature(inward, Point2{1.0, 0.7}), 0.5));
    CHECK(Near(*GaussianCurvature(inward, Point2{1.0, 0.7}), 0.25));
}

TEST_CASE("a cylinder bends round its axis only", "[surface][properties]")
{
    const Cylinder cylinder{Point3{0.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Point2 uv{0.8, 5.0};
    CHECK(Near(*GaussianCurvature(cylinder, uv), 0.0));
    CHECK(Near(*MeanCurvature(cylinder, uv), -1.0 / 6.0));
    const PrincipalCurvaturesResult k = *PrincipalCurvatures(cylinder, uv);
    CHECK(Near(k.Min, -1.0 / 3.0));
    CHECK(Near(k.Max, 0.0));
    CHECK(Parallel(k.MinDirection, Deriv(cylinder, uv, 1, 0)));
    CHECK(Parallel(k.MaxDirection, Vec3{0.0, 0.0, 1.0}));
    CheckFrame(k, *UnitNormal(cylinder, uv));
}

TEST_CASE("a torus is elliptic outside and hyperbolic inside", "[surface][properties]")
{
    const Torus torus{Point3{0.0, 0.0, 0.0}, 5.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    for (const double v : {0.0, 1.0, Pi / 2.0, 2.5, Pi}) {
        const Point2 uv{0.4, v};
        // Round the tube -1/r; along the parallels -cos v / (R + r cos v).
        const double tube = -1.0 / torus.MinorRadius;
        const double parallel =
            -std::cos(v) / (torus.MajorRadius + torus.MinorRadius * std::cos(v));
        CHECK(Near(*GaussianCurvature(torus, uv), tube * parallel));
        CHECK(Near(*MeanCurvature(torus, uv), 0.5 * (tube + parallel)));
        const PrincipalCurvaturesResult k = *PrincipalCurvatures(torus, uv);
        CHECK(Near(k.Min, std::min(tube, parallel)));
        CHECK(Near(k.Max, std::max(tube, parallel)));
        CHECK(Parallel(tube < parallel ? k.MinDirection : k.MaxDirection, Deriv(torus, uv, 0, 1)));
        CheckFrame(k, *UnitNormal(torus, uv));
    }
}

TEST_CASE("properties are undefined where the parametrisation is singular", "[surface][properties]")
{
    const Sphere sphere{Point3{0.0, 0.0, 0.0}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Point2 pole{0.5, Pi / 2.0};
    // cos(pi/2) is not exactly zero in floating point, so the sphere's pole
    // is only required not to produce a wrong normal silently; a cone of
    // radius 0 has its apex exactly at v = 0.
    const Cone cone{Point3{0.0, 0.0, 0.0}, 0.0, Pi / 4.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Point2 apex{0.5, 0.0};
    CHECK(UnitNormal(cone, apex).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(GaussianCurvature(cone, apex).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(PrincipalCurvatures(cone, apex).error().GetCode() == ErrorCode::DegenerateResult);
    const auto atPole = UnitNormal(sphere, pole);
    CHECK((!atPole || Parallel(*atPole, Vec3{0.0, 0.0, 1.0})));
    // A degenerate B-spline patch: one row of coincident poles.
    const BSplineSurface collapsed{1,
                                   1,
                                   {0.0, 0.0, 1.0, 1.0},
                                   {0.0, 0.0, 1.0, 1.0},
                                   {Point3{0.0, 0.0, 0.0}, Point3{0.0, 0.0, 0.0},
                                    Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}},
                                   {}};
    CHECK(UnitNormal(collapsed, Point2{0.0, 0.0}).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(MeanCurvature(collapsed, Point2{0.0, 0.0}).error().GetCode() ==
          ErrorCode::DegenerateResult);
}

TEST_CASE("curvatures survive the exact NURBS conversion", "[surface][properties][property]")
{
    const Torus torus{Point3{1.0, 2.0, 3.0}, 4.0, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Box2 range{Interval{0.2, 2.5}, Interval{-1.0, 3.0}};
    const Surface nurbs(*ToBSpline(torus, range));
    rc::prop("at every point of the converted patch", [&] {
        const Point2 uv{0.2 + 2.3 * *rc::gen::inRange(0, 1000) / 1000.0,
                        -1.0 + 4.0 * *rc::gen::inRange(0, 1000) / 1000.0};
        const Point3 p = Eval(torus, uv);
        const Point2 st = Project(nurbs, p);
        RC_ASSERT(Distance(Eval(nurbs, st), p) <= 1e-9);
        // Same surface, same orientation: equal curvatures, however
        // differently the two are parametrised.
        const PrincipalCurvaturesResult a = *PrincipalCurvatures(torus, uv);
        const PrincipalCurvaturesResult b = *PrincipalCurvatures(nurbs, st);
        RC_ASSERT(std::abs(a.Min - b.Min) <= 1e-7);
        RC_ASSERT(std::abs(a.Max - b.Max) <= 1e-7);
        RC_ASSERT(Dot(*UnitNormal(torus, uv), *UnitNormal(nurbs, st)) >= 1.0 - 1e-12);
    });
}
