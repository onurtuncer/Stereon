// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Exact NURBS conversion of the analytic surfaces: the converted surface
// covers exactly the original's patch over the range, agrees with it on the
// knot grid, and has the same normal orientation.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::CirclePeriod;
using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Transform;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-11;

struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{-1.0, 2.0, 1.0}, 0.7);
    return {Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

bool InRange(const double t, const Interval range, const double period)
{
    if (period == 0.0) {
        return t >= range.GetLower() - 1e-9 && t <= range.GetUpper() + 1e-9;
    }
    const double shifted = std::fmod(t - range.GetLower(), period);
    const double offset = shifted < 0.0 ? shifted + period : shifted;
    return offset <= range.GetUpper() - range.GetLower() + 1e-9 || offset >= period - 1e-9;
}

template <typename G>
Vec2 PeriodsOf(const G& g)
{
    if constexpr (Stereon::PeriodicSurface<G>) {
        return Periods(g);
    } else {
        static_cast<void>(g);
        return {};
    }
}

Vec3 Normal(const auto& s, const Point2& uv)
{
    return Cross(Deriv(s, uv, 1, 0), Deriv(s, uv, 0, 1));
}

template <typename G>
void CheckConversion(const G& g, const Box2& range)
{
    const auto converted = Stereon::ToBSpline(g, range);
    RC_ASSERT(converted.has_value());
    const BSplineSurface& s = *converted;
    RC_ASSERT(Domain(s) == range);
    for (const double u : s.KnotsU) {
        for (const double v : s.KnotsV) {
            RC_ASSERT(Distance(Eval(s, Point2{u, v}), Eval(g, Point2{u, v})) <= Tolerance);
        }
    }
    const auto at = [](const Interval i) {
        return i.GetLower() + (*rc::gen::inRange(1, 999) / 1000.0) * (i.GetUpper() - i.GetLower());
    };
    const Point2 st{at(range.X), at(range.Y)};
    const Point3 p = Eval(s, st);
    const Point2 uv = Project(g, p);
    const Vec2 periods = PeriodsOf(g);
    RC_ASSERT(Distance(Eval(g, uv), p) <= Tolerance);
    RC_ASSERT(InRange(uv.X, range.X, periods.X));
    RC_ASSERT(InRange(uv.Y, range.Y, periods.Y));
    RC_ASSERT(Dot(Normal(s, st), Normal(g, uv)) > 0.0);
    // Coverage: the original's point at st lies on the converted surface.
    const Point3 q = Eval(g, st);
    RC_ASSERT(Distance(Eval(s, Project(s, q)), q) <= 1e-9);
}

Interval RandomRange(const double start, const double maxLength)
{
    const double low = start * *rc::gen::inRange(-1000, 1000) / 1000.0;
    return {low, low + maxLength * *rc::gen::inRange(10, 1000) / 1000.0};
}

} // namespace

TEST_CASE("analytic surfaces convert exactly to NURBS", "[surface][conversion][property]")
{
    const Frame f = TiltedFrame();
    const Point3 o{0.5, -1.0, 2.0};
    rc::prop("planes, cylinders and cones", [&] {
        CheckConversion(Plane{o, f.X, f.Y}, Box2{RandomRange(3.0, 4.0), RandomRange(3.0, 4.0)});
        CheckConversion(Cylinder{o, 1.5, f.X, f.Y},
                        Box2{RandomRange(6.0, CirclePeriod), RandomRange(3.0, 4.0)});
        // v away from the apex at v = -Radius / sin(SemiAngle) = -3.
        CheckConversion(Cone{o, 1.5, 0.5, f.X, f.Y},
                        Box2{RandomRange(6.0, CirclePeriod),
                             Interval{-2.0, -2.0 + 4.0 * *rc::gen::inRange(1, 100) / 100.0}});
    });
    rc::prop("spheres and tori", [&] {
        const double band = *rc::gen::inRange(10, 280) / 100.0;
        CheckConversion(Sphere{o, 2.0, f.X, f.Y, false},
                        Box2{RandomRange(6.0, CirclePeriod), Interval{-1.4, -1.4 + band}});
        CheckConversion(Torus{o, 3.0, 1.0, f.X, f.Y},
                        Box2{RandomRange(6.0, CirclePeriod), RandomRange(6.0, CirclePeriod)});
    });
}

TEST_CASE("a plane converts with its own parametrisation", "[surface][conversion]")
{
    const Plane plane{Point3{1.0, 2.0, 3.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const BSplineSurface s =
        *Stereon::ToBSpline(plane, Box2{Interval{-1.0, 2.0}, Interval{0.0, 5.0}});
    REQUIRE(s.DegreeU == 1);
    REQUIRE(s.DegreeV == 1);
    REQUIRE(Distance(Eval(s, Point2{0.3, 4.1}), Eval(plane, Point2{0.3, 4.1})) <= 1e-14);
}

TEST_CASE("a B-spline surface converts to its restriction", "[surface][conversion]")
{
    const Sphere sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const BSplineSurface whole =
        *Stereon::ToBSpline(sphere, Box2{Interval{0.0, 3.0}, Interval{-1.0, 1.0}});
    const Box2 part{Interval{0.5, 2.0}, Interval{-0.25, 0.75}};
    const BSplineSurface piece = *Stereon::ToBSpline(whole, part);
    REQUIRE(Domain(piece) == part);
    for (int i = 0; i <= 6; ++i) {
        for (int j = 0; j <= 6; ++j) {
            const Point2 uv{0.5 + i * 0.25, -0.25 + j * (1.0 / 6.0)};
            REQUIRE(Distance(Eval(piece, uv), Eval(whole, uv)) <= 1e-13);
        }
    }
    REQUIRE(*Stereon::ToBSpline(whole, Domain(whole)) == whole);
}

TEST_CASE("surface conversion through the wrapper and its errors", "[surface][conversion]")
{
    const Torus torus{Point3{}, 3.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Box2 range{Interval{0.0, 1.0}, Interval{0.0, 2.0}};
    REQUIRE(*Stereon::ToBSpline(Surface(torus), range) == *Stereon::ToBSpline(torus, range));
    const auto invalid = ErrorCode::InvalidInput;
    const Sphere sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Stereon::ToBSpline(sphere, Box2{Interval{0.0, 1.0}, Interval{0.0, 2.0}})
                .error()
                .GetCode() == invalid); // beyond the north pole
    REQUIRE(
        Stereon::ToBSpline(torus, Box2{Interval{0.0, 7.0}, Interval{0.0, 1.0}}).error().GetCode() ==
        invalid); // more than a period
    REQUIRE(Stereon::ToBSpline(Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}},
                               Box2::Everything())
                .error()
                .GetCode() == invalid);
    // A whole sphere and a whole torus are fine.
    REQUIRE(Stereon::ToBSpline(sphere, Domain(sphere)).has_value());
    REQUIRE(Stereon::ToBSpline(torus, Domain(torus)).has_value());
}
