// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The fixtures of the adversarial audit of 2026-10-28 (A1-A5): valid
// inputs on which the B-spline rearrangements returned success with a wrong
// or invalid curve. Each operation must now return the curve it promises, or
// a structured error when no valid curve represents it; never an invalid
// success. And the properties that would have caught them: elevation keeps
// the curve at every degree and with strongly uneven knots, and every
// result passes the factory's checks.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

import Stereon.Curve;

using Stereon::BSplineCurve;
using Stereon::ErrorCode;
using Stereon::Point3;

namespace {

/// The largest distance between the two curves at n + 1 parameters spread
/// over the domain of `a`, mapped onto that of `b`.
double MaxDistance(const BSplineCurve& a, const BSplineCurve& b, const int n = 400)
{
    const Stereon::Interval da = Domain(a);
    const Stereon::Interval db = Domain(b);
    double worst = 0.0;
    for (int k = 0; k <= n; ++k) {
        const double s = static_cast<double>(k) / n;
        const double ta = da.GetLower() + s * (da.GetUpper() - da.GetLower());
        const double tb = db.GetLower() + s * (db.GetUpper() - db.GetLower());
        worst = std::max(worst, Distance(Eval(a, std::min(ta, da.GetUpper())),
                                         Eval(b, std::min(tb, db.GetUpper()))));
    }
    return worst;
}

double LargestPoleNorm(const BSplineCurve& c)
{
    double largest = 0.0;
    for (const Point3& p : c.Poles) {
        largest = std::max(largest, Length(Stereon::ToVec(p)));
    }
    return largest;
}

/// True when the curve passes the factory's checks of knots, poles and
/// weights.
bool Valid(const BSplineCurve& c)
{
    return Stereon::MakeBSplineCurve(c.Degree, c.Knots, c.Poles, c.Weights).has_value();
}

/// A1: degree 9, a first interior knot at 1e-8, poles in [0, 1] x [0, 2].
BSplineCurve AuditA1()
{
    std::vector<double> knots(10, 0.0);
    knots.insert(knots.end(), {1e-8, 0.5, 0.75});
    knots.insert(knots.end(), 10, 1.0);
    std::vector<Point3> poles;
    for (int i = 0; i < 13; ++i) {
        poles.push_back({static_cast<double>(i % 2), static_cast<double>(i % 3), 0.0});
    }
    return *Stereon::MakeBSplineCurve(9, knots, poles);
}

} // namespace

TEST_CASE("audit A1: elevation over a tiny knot interval keeps the curve",
          "[curve][bspline][audit]")
{
    const BSplineCurve c = AuditA1();
    const auto elevated = Stereon::ElevatedDegree(c, 10);
    REQUIRE(elevated.has_value());
    REQUIRE(elevated->Degree == 10);
    REQUIRE(Valid(*elevated));
    REQUIRE(Distance(Eval(*elevated, 0.155), Eval(c, 0.155)) <= 1e-12);
    REQUIRE(MaxDistance(c, *elevated) <= 1e-12);
    // The poles stay in the box the curve lives in, give or take rounding.
    REQUIRE(LargestPoleNorm(*elevated) <= LargestPoleNorm(c) * (1.0 + 1e-9));
}

TEST_CASE("audit A2: a curve far from the origin is periodic only if it is",
          "[curve][bspline][audit]")
{
    constexpr double origin = 1e16;
    const auto c = Stereon::MakeBSplineCurve(
        2, {origin, origin, origin, origin + 2.0, origin + 4.0, origin + 4.0, origin + 4.0},
        {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}, Point3{0.0, 0.0, 0.0},
         Point3{1.0, 0.0, 0.0}});
    REQUIRE(c.has_value());
    REQUIRE_FALSE(Stereon::IsPeriodic(*c));
    const auto periodic = Stereon::AsPeriodic(*c, 0.0);
    REQUIRE_FALSE(periodic.has_value());
    REQUIRE(periodic.error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("audit A3: concatenation refuses a piece whose span would vanish",
          "[curve][bspline][audit]")
{
    constexpr double origin = 1e16;
    const auto first = Stereon::MakeBSplineCurve(1, {origin, origin, origin + 4.0, origin + 4.0},
                                                 {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}});
    const auto second = Stereon::MakeBSplineCurve(1, {0.0, 0.0, 1.0, 1.0},
                                                  {Point3{1.0, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}});
    REQUIRE((first.has_value() && second.has_value()));
    const auto joined = Stereon::Concatenated(*first, *second, 0.0);
    REQUIRE_FALSE(joined.has_value());
    REQUIRE(joined.error().GetCode() == ErrorCode::NotSupported);
}

TEST_CASE("audit A4: concatenation keeps the weights in their range", "[curve][bspline][audit]")
{
    const double big = std::ldexp(1.0, 256);
    const double small = std::ldexp(1.0, -256);
    const auto first = Stereon::MakeBSplineCurve(
        1, {0.0, 0.0, 1.0, 1.0}, {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}}, {big, big});
    const auto second = Stereon::MakeBSplineCurve(
        1, {0.0, 0.0, 1.0, 1.0}, {Point3{1.0, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}}, {small, big});
    REQUIRE((first.has_value() && second.has_value()));
    const auto joined = Stereon::Concatenated(*first, *second, 0.0);
    REQUIRE(joined.has_value());
    REQUIRE(Valid(*joined));
    for (const double t : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        REQUIRE(Distance(Eval(*joined, t), Eval(*first, t)) <= 1e-15);
        REQUIRE(Distance(Eval(*joined, 1.0 + t), Eval(*second, t)) <= 1e-15);
    }
}

TEST_CASE("audit A5: the periodic form refuses knots it cannot extend", "[curve][bspline][audit]")
{
    const auto c = Stereon::MakeBSplineCurve(
        2, {0.0, 0.0, 0.0, 1e308, 1e308, 1e308},
        {Point3{0.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}, Point3{0.0, 0.0, 0.0}});
    REQUIRE(c.has_value());
    const auto periodic = Stereon::AsPeriodic(*c, 0.0);
    REQUIRE_FALSE(periodic.has_value());
    REQUIRE(periodic.error().GetCode() == ErrorCode::NotSupported);
}

TEST_CASE("elevation keeps the curve at every degree and with uneven knots",
          "[curve][bspline][audit][property]")
{
    // Degrees 1 to 24 raised by 1 to 25 - p, interior knots as close as
    // 1e-12 to each other or to an end, polynomial or rational: the result
    // is valid and within a rounding bound of the curve everywhere.
    rc::prop("for random curves and target degrees", [] {
        const int p = *rc::gen::inRange(1, Stereon::BSplineMaxDegree);
        const int q = *rc::gen::inRange(p + 1, Stereon::BSplineMaxDegree + 1);
        const bool rational = *rc::gen::arbitrary<bool>();
        const auto interior = static_cast<std::size_t>(*rc::gen::inRange(1, 5));
        std::vector<double> knots(static_cast<std::size_t>(p) + 1, 0.0);
        std::vector<double> inside;
        for (std::size_t k = 0; k < interior; ++k) {
            const int kind = *rc::gen::inRange(0, 3);
            const double tiny = std::pow(10.0, -*rc::gen::inRange(4, 13));
            inside.push_back(
                kind == 0 ? tiny : (kind == 1 ? 1.0 - tiny : *rc::gen::inRange(1, 1000) / 1000.0));
        }
        std::ranges::sort(inside);
        knots.insert(knots.end(), inside.begin(), inside.end());
        knots.insert(knots.end(), static_cast<std::size_t>(p) + 1, 1.0);
        const std::size_t count = knots.size() - static_cast<std::size_t>(p) - 1;
        std::vector<Point3> poles;
        std::vector<double> weights;
        for (std::size_t i = 0; i < count; ++i) {
            poles.push_back({*rc::gen::inRange(-50, 51) / 10.0, *rc::gen::inRange(-50, 51) / 10.0,
                             *rc::gen::inRange(-50, 51) / 10.0});
            if (rational) {
                weights.push_back(*rc::gen::inRange(2, 40) / 10.0);
            }
        }
        const auto made = Stereon::MakeBSplineCurve(p, knots, poles, weights);
        RC_PRE(made.has_value());
        const auto elevated = Stereon::ElevatedDegree(*made, q);
        RC_ASSERT(elevated.has_value());
        RC_ASSERT(elevated->Degree == q);
        RC_ASSERT(Valid(*elevated));
        const double scale = 1.0 + LargestPoleNorm(*made);
        RC_ASSERT(MaxDistance(*made, *elevated, 200) <= 1e-9 * scale);
    });
}
