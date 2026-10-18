// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <numbers>

import Stereon.Surface;

#include "Formatting.hpp"

using namespace Stereon;

namespace {
constexpr double Tolerance = 1e-10;

BSplineCurve Profile()
{
    // A rational, non-planar cubic: neither a meridian nor a conic.
    return {
        3, {0, 0, 0, 0, 1, 1, 1, 1}, {{2, 0, -1}, {3, 1, 0}, {1, -1, 2}, {2, 0, 3}}, {1, 2, 3, 1}};
}

// Independent rational quadratic half-angle identity on a single arc.
double AngleAt(const double t, const double start, const double span)
{
    return start + 2.0 * std::atan2(t * std::sin(span / 2.0), 1.0 - t + t * std::cos(span / 2.0));
}
} // namespace

TEST_CASE("extrusion preserves rational profile parameters and mixed derivatives",
          "[surface][swept]")
{
    const auto profile = Profile();
    const Vec3 direction{1, -2, 3};
    const auto result = MakeExtrusionSurface(profile, direction, Interval{-2, 4});
    REQUIRE(result.has_value());
    const auto& surface = *result;
    REQUIRE(Domain(surface) == Box2{Interval{0, 1}, Interval{-2, 4}});
    for (const double u : {0.0, 0.17, 0.5, 0.83, 1.0}) {
        for (const double v : {-2.0, -0.3, 1.0, 4.0}) {
            REQUIRE(Distance(Eval(surface, Point2{u, v}), Eval(profile, u) + v * direction) <
                    Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 1, 0) - Deriv(profile, u, 1)) < Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 0, 1) - direction) < Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 1, 1)) < Tolerance);
        }
    }
    const Surface erased{surface};
    REQUIRE(Distance(Eval(erased, Point2{0.4, 2}), Eval(surface, Point2{0.4, 2})) == 0.0);
}

TEST_CASE("extrusion commutes with translation", "[surface][swept][property]")
{
    rc::prop("translation of the profile translates every generated point", [] {
        const int seed = *rc::gen::inRange(-100, 101);
        auto shifted = Profile();
        const Vec3 shift{static_cast<double>(seed), -2.0, 3.0};
        for (auto& pole : shifted.Poles) {
            pole = pole + shift;
        }
        const Vec3 direction{2, 3, -1};
        const auto a = MakeExtrusionSurface(Profile(), direction, Interval{-1, 2});
        const auto b = MakeExtrusionSurface(shifted, direction, Interval{-1, 2});
        RC_ASSERT(a.has_value());
        RC_ASSERT(b.has_value());
        const double u = static_cast<double>(*rc::gen::inRange(0, 101)) / 100.0;
        RC_ASSERT(Distance(Eval(*a, Point2{u, 0.7}) + shift, Eval(*b, Point2{u, 0.7})) < Tolerance);
    });
}

TEST_CASE("revolution matches Rodrigues rotation between angular knots", "[surface][swept]")
{
    const auto profile = Profile();
    const Point3 origin{1, -2, 0.5};
    const Vec3 axis{2, -1, 3};
    const double start = 5.5; // The sweep crosses the conventional 2 pi seam.
    const double span = 1.0;
    const auto result = MakeRevolutionSurface(profile, origin, axis, Interval{start, start + span});
    REQUIRE(result.has_value());
    for (const double t : {0.0, 0.13, 0.5, 0.79, 1.0}) {
        const double angle = AngleAt(t, start, span);
        const auto rotation = MakeRotation(origin, axis, angle);
        for (const double v : {0.0, 0.23, 0.6, 1.0}) {
            REQUIRE(Distance(Eval(*result, Point2{start + t * span, v}),
                             Transformed(Eval(profile, v), rotation)) < Tolerance);
            REQUIRE(Length(Deriv(*result, Point2{start + t * span, v}, 0, 1) -
                           Transformed(Deriv(profile, v, 1), rotation)) < Tolerance);
        }
    }
}

TEST_CASE("full revolution closes and admits poles on the axis", "[surface][swept]")
{
    const BSplineCurve profile{1, {0, 0, 1, 1}, {{0, 0, 0}, {2, 0, 3}}, {}};
    const auto result =
        MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 2}, Interval{0, CirclePeriod});
    REQUIRE(result.has_value());
    for (const double v : {0.0, 0.3, 1.0}) {
        REQUIRE(Distance(Eval(*result, Point2{0, v}), Eval(*result, Point2{CirclePeriod, v})) <
                Tolerance);
    }
    for (const double u : {0.0, 0.7, 2.0, 4.0, CirclePeriod}) {
        REQUIRE(Distance(Eval(*result, Point2{u, 0}), Point3{}) < Tolerance);
        const auto p = Eval(*result, Point2{u, 0.5});
        REQUIRE(std::abs(std::hypot(p.X, p.Y) - 1.0) < Tolerance);
        REQUIRE(std::abs(p.Z - 1.5) < Tolerance);
    }
}

TEST_CASE("swept factories reject malformed and overflowing inputs", "[surface][swept]")
{
    const auto profile = Profile();
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const Vec3 direction : {Vec3{}, Vec3{nan, 0, 1}, Vec3{0, infinity, 1}}) {
        REQUIRE_FALSE(MakeExtrusionSurface(profile, direction, Interval{0, 1}));
        REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{}, direction, Interval{0, 1}));
    }
    for (const Interval range : {Interval{1, 1}, Interval{0, infinity}}) {
        REQUIRE_FALSE(MakeExtrusionSurface(profile, Vec3{0, 0, 1}, range));
        REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 1}, range));
    }
    REQUIRE_FALSE(
        MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 1}, Interval{0, 2 * CirclePeriod}));
    REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{nan, 0, 0}, Vec3{0, 0, 1}, Interval{0, 1}));
    auto invalid = profile;
    invalid.Weights.pop_back();
    const auto rejected = MakeExtrusionSurface(invalid, Vec3{0, 0, 1}, Interval{0, 1});
    REQUIRE_FALSE(rejected);
    REQUIRE(rejected.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE_FALSE(MakeRevolutionSurface(invalid, Point3{}, Vec3{0, 0, 1}, Interval{0, 1}));
    const double huge = std::numeric_limits<double>::max();
    REQUIRE_FALSE(MakeExtrusionSurface(profile, Vec3{huge, 0, 0}, Interval{0, 2}));
}
