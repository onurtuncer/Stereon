// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Two plane B-splines that touch (v0.2 ledger, 2D curve-curve): F's
// Jacobian is singular there, so Krawczyk's test never settles the pair.
// The handler takes both pieces as graphs over one direction and the
// difference of their heights h(x), whose derivative of some order is
// proved of one sign, and reports the touch bracketed, Tangential and
// undecided, never dropped; and two splines that are one curve in two knot
// vectors, refined exactly, overlap. The expected answers are closed
// forms, all poles exact: the cup y = x^2 against the cap y = -x^2 (order
// two at the origin), against y = -x^2 - g (apart) and y = g - x^2 (two
// crossings at x = +-sqrt(g / 2)), g = 3 2^-40; the quartic
// y = x^2 + (x / 3)^4 against the quadratic y = x^2 (order four); a
// spline and its copy with a knot inserted. Metamorphic relations
// (ADR-0053): exchange, reversal, translation, scaling by powers of two, a
// rational form.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::BSplineCurve2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Result;

namespace {

/// A Bezier curve on [0, 1] with poles moved by `at` and scaled by k.
BSplineCurve2d Bezier(std::vector<Point2> poles, const Point2& at = {}, const double k = 1.0,
                      std::vector<double> weights = {})
{
    const std::size_t n = poles.size();
    for (Point2& p : poles) {
        p = Point2{at.X + k * p.X, at.Y + k * p.Y};
    }
    std::vector<double> knots(n, 0.0);
    knots.resize(2 * n, 1.0);
    auto made = Stereon::MakeBSplineCurve(static_cast<int>(n) - 1, std::move(knots),
                                          std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// y = c + d x^2 over x in [-1, 1], x = 2s - 1: Bernstein heights c + d,
/// c - d, c + d.
BSplineCurve2d Parabola(const double c, const double d, const Point2& at = {}, const double k = 1.0,
                        std::vector<double> weights = {})
{
    return Bezier({{-1.0, c + d}, {0.0, c - d}, {1.0, c + d}}, at, k, std::move(weights));
}

/// The rounding a touch of the given order is placed to, in its parameter:
/// the order-th root of the rounding relative to the curve's size, which
/// grows with the distance from the origin.
double Reach(const Point2& at, const double k, const int order)
{
    const double rounding = 0x1p-52 * (1.0 + std::max(std::abs(at.X), std::abs(at.Y)) / k);
    return 16.0 * std::pow(rounding, 1.0 / order);
}

/// One touch, both ways round: Tangential, undecided, at s and t, `at`.
void RequireTouch(const Curve2d& a, const Curve2d& b, const double s, const double t,
                  const Point2& at, const double reach, const double k)
{
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "b x a" : "a x b"));
        const Result<PlanarIntersection> result =
            swapped ? Intersect(b, a, ctx) : Intersect(a, b, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        REQUIRE(result->Overlaps.empty());
        REQUIRE(result->Points.size() == 1);
        const PlanarPoint& p = result->Points.front();
        REQUIRE(p.Contact == ContactType::Tangential);
        REQUIRE_FALSE(p.Decided);
        REQUIRE_THAT(swapped ? p.ParameterB : p.ParameterA, WithinAbs(s, reach));
        REQUIRE_THAT(swapped ? p.ParameterA : p.ParameterB, WithinAbs(t, reach));
        REQUIRE_THAT(p.Position.X, WithinAbs(at.X, 8.0 * k * reach));
        REQUIRE_THAT(p.Position.Y, WithinAbs(at.Y, 8.0 * k * reach));
    }
}

} // namespace

TEST_CASE("two B-splines touching to order two are a bracketed touch, 3 2^-40 either side "
          "apart or two crossings",
          "[intersect][planar][spline][tangency]")
{
    for (const Point2 at : {Point2{0.0, 0.0}, Point2{1024.0, -512.0}}) {
        for (const double k : {1.0, 0x1p10, 0x1p-10}) {
            DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ") scaled by " << k)
            {
                const double reach = Reach(at, k, 2);
                const Curve2d cup{Parabola(0.0, 1.0, at, k)};
                const Curve2d cap{Parabola(0.0, -1.0, at, k)};
                RequireTouch(cup, cap, 0.5, 0.5, at, reach, k);
                // Reversed: s -> 1 - s keeps the touch at 1/2.
                RequireTouch(Curve2d{Reversed(Parabola(0.0, 1.0, at, k))}, cap, 0.5, 0.5, at, reach,
                             k);
                // The cap as a rational with equal weights: the same curve.
                RequireTouch(cup, Curve2d{Parabola(0.0, -1.0, at, k, {3.0, 3.0, 3.0})}, 0.5, 0.5,
                             at, reach, k);
            }
        }
    }
    SECTION("3 2^-40 apart and across")
    {
        const double g = 3.0 * 0x1p-40;
        const Curve2d cup{Parabola(0.0, 1.0)};
        const auto apart = Intersect(cup, Curve2d{Parabola(-g, -1.0)}, ctx);
        INFO(IntersectTest::ErrorText(apart));
        REQUIRE(apart.has_value());
        REQUIRE(apart->IsEmpty());
        const auto two = Intersect(cup, Curve2d{Parabola(g, -1.0)}, ctx);
        INFO(IntersectTest::ErrorText(two));
        REQUIRE(two.has_value());
        REQUIRE(two->Points.size() == 2);
        const double x = std::sqrt(g / 2.0);
        for (std::size_t i = 0; i < 2; ++i) {
            const PlanarPoint& p = two->Points[i];
            REQUIRE(p.Contact == ContactType::Crossing);
            REQUIRE_THAT(p.ParameterA, WithinAbs((1.0 + (i == 0 ? -x : x)) / 2.0, 1e-7));
        }
    }
}

TEST_CASE("two B-splines touching to order four are a bracketed touch",
          "[intersect][planar][spline][tangency]")
{
    for (const Point2 at : {Point2{0.0, 0.0}, Point2{1024.0, -512.0}}) {
        DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ")")
        {
            // y = x^2 + (x / 3)^4 and y = x^2 over x in [-3, 3], x = 3(2s - 1).
            const Curve2d quartic{
                Bezier({{-3.0, 10.0}, {-1.5, -1.0}, {0.0, -2.0}, {1.5, -1.0}, {3.0, 10.0}}, at)};
            const Curve2d quadratic{Bezier({{-3.0, 9.0}, {0.0, -9.0}, {3.0, 9.0}}, at)};
            RequireTouch(quartic, quadratic, 0.5, 0.5, at, Reach(at, 1.0, 4), 1.0);
        }
    }
}

TEST_CASE("a B-spline and its copy with knots inserted overlap; inexact copies fail closed",
          "[intersect][planar][spline][overlap]")
{
    const BSplineCurve2d arch = Parabola(0.0, 1.0);
    for (const double knot : {0.5, 0.25}) {
        const auto refined = InsertKnot(arch, knot, 1);
        REQUIRE(refined.has_value());
        for (const bool swapped : {false, true}) {
            const Curve2d a{swapped ? *refined : arch};
            const Curve2d b{swapped ? arch : *refined};
            const auto whole = Intersect(a, b, ctx);
            INFO(IntersectTest::ErrorText(whole));
            REQUIRE(whole.has_value());
            REQUIRE(whole->Points.empty());
            REQUIRE(whole->Overlaps.size() == 1);
            REQUIRE(whole->Overlaps.front().OnA == Interval{0.0, 1.0});
            REQUIRE(whole->Overlaps.front().OnB == Interval{0.0, 1.0});
            // The requested domains' common stretch.
            const auto part = Intersect(a, Interval{0.125, 0.75}, b, Interval{0.5, 1.0}, ctx);
            REQUIRE(part.has_value());
            REQUIRE(part->Overlaps.size() == 1);
            REQUIRE(part->Overlaps.front().OnA == Interval{0.5, 0.75});
        }
    }
    // A knot whose insertion rounds proves nothing: the pair runs along
    // itself and fails closed, never an empty or partial answer.
    const auto third = InsertKnot(arch, 1.0 / 3.0, 1);
    REQUIRE(third.has_value());
    const auto along = Intersect(Curve2d{arch}, Curve2d{*third}, ctx);
    REQUIRE_FALSE(along.has_value());
    REQUIRE(along.error().GetCode() == ErrorCode::IntersectionFailed);
}
