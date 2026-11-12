// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A B-spline's contact of order four with a line or a conic (v0.2 ledger,
// 2D curve-curve: touches reported bracketed, never dropped). The implicit
// equation along the spline vanishes there with its first three
// derivatives, so no piece about the touch is apart, monotone or convex in
// doubles; the handler isolates it by a derivative of higher order proved
// of one sign. The expected answers are closed forms: the quartic Bezier
// (3(2s - 1), 9(2s - 1)^2 + (2s - 1)^4) is y = x^2 + (x / 3)^4 with exact
// poles, touching the parabola y = x^2 to order four at s = 1/2, the
// origin; the quartic (2s - 1, (2s - 1)^4) is y = x^4, touching the axis
// to order four, missing y = -g and crossing y = g at
// x = +-g^(1/4), g = 3 2^-40. Metamorphic relations (ADR-0053): exchange,
// reversal, translation, scaling by powers of two, a rational form.

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
using Stereon::Line2d;
using Stereon::Parabola2d;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Result;
using Stereon::Vec2;

namespace {

/// A quartic Bezier on [0, 1] with poles moved by `at` and scaled by k.
BSplineCurve2d Quartic(std::vector<Point2> poles, const Point2& at, const double k,
                       std::vector<double> weights = {})
{
    for (Point2& p : poles) {
        p = Point2{at.X + k * p.X, at.Y + k * p.Y};
    }
    auto made = Stereon::MakeBSplineCurve(4, {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// y = x^2 + (x / 3)^4 over x in [-3, 3], x = 3(2s - 1).
BSplineCurve2d Cup(const Point2& at, const double k, std::vector<double> weights = {})
{
    return Quartic({{-3.0, 10.0}, {-1.5, -1.0}, {0.0, -2.0}, {1.5, -1.0}, {3.0, 10.0}}, at, k,
                   std::move(weights));
}

/// y = x^4 over x in [-1, 1], x = 2s - 1.
BSplineCurve2d Fourth(const Point2& at, const double k, std::vector<double> weights = {})
{
    return Quartic({{-1.0, 1.0}, {-0.5, -1.0}, {0.0, 1.0}, {0.5, -1.0}, {1.0, 1.0}}, at, k,
                   std::move(weights));
}

/// y = x^2 / k about `at`: points at + (-t, t^2 / (4 f)) with f = k / 4.
Curve2d ParabolaAt(const Point2& at, const double k)
{
    return Curve2d{Parabola2d{at, k / 4.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}}};
}

Curve2d HorizontalAt(const Point2& at)
{
    return Curve2d{Line2d{at, Vec2{1.0, 0.0}}};
}

/// The one point of the spline against the curve, both ways round: a touch
/// at the spline's s = 1/2, at `at`, tangential, never dropped. It is
/// Decided only if settled exactly, which this handler does not claim.
void RequireTouch(const Curve2d& spline, const Curve2d& other, const Point2& at, const double k)
{
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "other x spline" : "spline x other"));
        const Result<PlanarIntersection> result =
            swapped ? Intersect(other, spline, ctx) : Intersect(spline, other, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        REQUIRE(result->Overlaps.empty());
        REQUIRE(result->Points.size() == 1);
        const PlanarPoint& p = result->Points.front();
        REQUIRE(p.Contact == ContactType::Tangential);
        REQUIRE_FALSE(p.Decided);
        const double s = swapped ? p.ParameterB : p.ParameterA;
        // The touch is bracketed to where rounding hides the quartic: the
        // fourth root of the rounding relative to the curve's size, which
        // grows with the distance from the origin (2^-13 near it).
        const double rounding = 0x1p-52 * (1.0 + std::max(std::abs(at.X), std::abs(at.Y)) / k);
        const double reach = 8.0 * std::pow(rounding, 0.25);
        REQUIRE_THAT(s, WithinAbs(0.5, reach));
        REQUIRE_THAT(p.Position.X, WithinAbs(at.X, 8.0 * k * reach));
        REQUIRE_THAT(p.Position.Y, WithinAbs(at.Y, 8.0 * k * reach));
    }
}

} // namespace

TEST_CASE("a B-spline's contact of order four with a parabola is reported as a touch",
          "[intersect][planar][spline][osculation]")
{
    for (const Point2 at : {Point2{0.0, 0.0}, Point2{1024.0, -512.0}}) {
        for (const double k : {1.0, 0x1p10, 0x1p-10}) {
            DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ") scaled by " << k)
            {
                RequireTouch(Curve2d{Cup(at, k)}, ParabolaAt(at, k), at, k);
                // Reversed, s -> 1 - s keeps the touch at 1/2.
                RequireTouch(Curve2d{Reversed(Cup(at, k))}, ParabolaAt(at, k), at, k);
                // A rational form of the same curve, equal weights.
                RequireTouch(Curve2d{Cup(at, k, {3.0, 3.0, 3.0, 3.0, 3.0})}, ParabolaAt(at, k), at,
                             k);
            }
        }
    }
}

TEST_CASE("a B-spline's contact of order four with a line is reported as a touch; 3 2^-40 "
          "either side it is proved apart or two crossings",
          "[intersect][planar][spline][osculation]")
{
    for (const Point2 at : {Point2{0.0, 0.0}, Point2{1024.0, -512.0}}) {
        DYNAMIC_SECTION("at (" << at.X << ", " << at.Y << ")")
        {
            const Curve2d fourth{Fourth(at, 1.0)};
            RequireTouch(fourth, HorizontalAt(at), at, 1.0);
            RequireTouch(Curve2d{Reversed(Fourth(at, 1.0))}, HorizontalAt(at), at, 1.0);
            const double g = 3.0 * 0x1p-40;
            if (at.Y == 0.0) {
                // Below the vertex: proved apart, both ways round.
                const auto apart = Intersect(fourth, HorizontalAt(Point2{at.X, at.Y - g}), ctx);
                INFO(IntersectTest::ErrorText(apart));
                REQUIRE(apart.has_value());
                REQUIRE(apart->IsEmpty());
                // Above it: two proved crossings at x = +-g^(1/4).
                const auto two = Intersect(fourth, HorizontalAt(Point2{at.X, at.Y + g}), ctx);
                INFO(IntersectTest::ErrorText(two));
                REQUIRE(two.has_value());
                REQUIRE(two->Points.size() == 2);
                const double x = std::pow(g, 0.25);
                for (std::size_t i = 0; i < 2; ++i) {
                    const PlanarPoint& p = two->Points[i];
                    REQUIRE(p.Contact == ContactType::Crossing);
                    REQUIRE(p.Decided);
                    // Proved, but y' = 8 x^3 is 2^-26 there: the coefficients'
                    // rounding, 2^-52, leaves a bracket of a few 2^-26.
                    REQUIRE_THAT(p.ParameterA, WithinAbs((1.0 + (i == 0 ? -x : x)) / 2.0, 1e-7));
                }
            }
        }
    }
}
