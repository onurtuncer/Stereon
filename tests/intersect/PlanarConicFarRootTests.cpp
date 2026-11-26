// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The fixed replay of a counterexample of *plane crossings: second-order
// Taylor decay and the transpose*
// (tests/intersect/PlanarLinearizationTests.cpp, seed 10340478410048779848, CI
// on #344, 2026-11-26). Two parabolas share a frame turned by 34/64 of a turn,
// with focal length 1/4 and vertices
// (-1, -1) and (-1, -1/16). Translates of one shape in one frame meet
// exactly once: (s^2 - t^2) X + (s - t) Y = V_B - V_A has one solution.
// The conic handler takes B by its implicit equation y^2 = 4 f x in B's
// frame, which is orthonormal only to rounding (X . X != 1 in the
// doubles). Along A's chart that equation keeps a leading coefficient of
// -7.7e-19, and so a second root near u = 4.2e-19 of the far chart, at
// t ~ 2.4e18. The handler reported it as a decided crossing at a point
// ~5.6e36 that is no point of B. The products then refused the whole
// result as stale (InvalidInput), which was the property's failure.
// Every point the handler reports must lie on both parametrised curves
// within the resolution; a root that does not is beyond the doubles'
// reach, and the handler fails closed.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Curve2d;
using Stereon::DerivativeStatus;
using Stereon::ErrorCode;
using Stereon::Parabola2d;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::PlanarPoint;
using Stereon::PlaneCurveBinding;
using Stereon::PlaneCurveComponent;
using Stereon::Point2;
using Stereon::SensitivityContext;
using Stereon::Vec2;

TEST_CASE("a root of the implicit equation off the parametric conic fails closed",
          "[intersect][planar][conic]")
{
    const double angle = 34.0 * 2.0 * std::numbers::pi / 64.0;
    const Vec2 x{std::cos(angle), std::sin(angle)};
    const Vec2 y{-x.Y, x.X};
    const Curve2d lower{Parabola2d{Point2{-1.0, -1.0}, 0.25, x, y}};
    const Curve2d upper{Parabola2d{Point2{-1.0, -0.0625}, 0.25, x, y}};
    for (const bool exchanged : {false, true}) {
        INFO((exchanged ? "exchanged" : "as given"));
        const Curve2d& a = exchanged ? upper : lower;
        const Curve2d& b = exchanged ? lower : upper;
        const auto at = Stereon::Intersect(a, b, ctx);
        if (!at) {
            // Failing closed is an answer; a wrong point is not.
            CHECK(at.error().GetCode() == ErrorCode::IntersectionFailed);
            continue;
        }
        // Whatever is reported lies on both curves, and the products take
        // it: no point is refused as stale.
        for (const PlanarPoint& p : at->Points) {
            CHECK(Distance(Eval(a, p.ParameterA), p.Position) <= ctx.GetResolution());
            CHECK(Distance(Eval(b, p.ParameterB), p.Position) <= ctx.GetResolution());
        }
        ParameterSet set;
        const ParameterId id = *set.Add("p", ParameterQuantity::Length, 1.0);
        SensitivityContext sensitivity(std::move(set));
        REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
        const std::vector<PlaneCurveBinding> bindings{{id, PlaneCurveComponent::CenterX}};
        const auto jvp = Stereon::IntersectJvp(a, bindings, b, {}, *at, sensitivity);
        INFO(IntersectTest::ErrorText(jvp));
        CHECK(jvp.has_value());
    }
    // The one true crossing, s + t = (D . X) / (D . Y) and s - t = D . Y for
    // D = V_upper - V_lower, is found where nothing far out is in play: over
    // bounded domains.
    const Stereon::Interval near{-4.0, 4.0};
    const auto bounded = Stereon::Intersect(lower, near, upper, near, ctx);
    INFO(IntersectTest::ErrorText(bounded));
    REQUIRE(bounded.has_value());
    REQUIRE(bounded->Points.size() == 1);
    const Vec2 d{0.0, 0.9375};
    const double sum = Dot(d, x) / Dot(d, y);
    const double difference = Dot(d, y);
    CHECK(std::abs(bounded->Points.front().ParameterA - 0.5 * (sum + difference)) <= 1e-12);
    CHECK(std::abs(bounded->Points.front().ParameterB - 0.5 * (sum - difference)) <= 1e-12);
}
