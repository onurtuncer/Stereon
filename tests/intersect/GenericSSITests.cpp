// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Phase 2: a successful empty result requires a mathematical separation
// proof, never simply failure to find a seed or a small iteration budget.

#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <vector>

import Stereon.Intersect;

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Context;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();

/// A bilinear NURBS patch over a unit parameter square, with the optional
/// diagonal pole deliberately moved to stretch the convex hull.
BSplineSurface Patch(const double z, const bool rational = false)
{
    std::vector<Point3> poles;
    poles.push_back({0.0, 0.0, z});
    poles.push_back({0.0, 1.0, z});
    poles.push_back({1.0, 0.0, z});
    poles.push_back({1.0, 1.0, z});
    const std::vector<double> weights =
        rational ? std::vector<double>{1.0, 3.0, 2.0, 4.0} : std::vector<double>{};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

const Box2 middle{Interval{0.25, 0.75}, Interval{0.25, 0.75}};

} // namespace

TEST_CASE("separated polynomial and rational NURBS control nets prove empty SSI",
          "[intersect][generic][broadphase]")
{
    for (const bool rationalA : {false, true}) {
        for (const bool rationalB : {false, true}) {
            const Surface a{Patch(0.0, rationalA)};
            const Surface b{Patch(10.0, rationalB)};
            for (const Box2 box : {Box2::Everything(), middle}) {
                const auto forward = Intersect(a, box, b, box, 1e-6, ctx);
                const auto reverse = Intersect(b, box, a, box, 1e-6, ctx);
                REQUIRE(forward.has_value());
                REQUIRE(reverse.has_value());
                REQUIRE(forward->IsEmpty());
                REQUIRE(reverse->IsEmpty());
            }
        }
    }
}

TEST_CASE("touching or overlapping NURBS hulls never masquerade as proved separation",
          "[intersect][generic][broadphase]")
{
    const Surface base{Patch(0.0)};
    for (const double z : {0.0, 1e-7}) {
        const Surface other{Patch(z)};
        const auto result = Intersect(base, middle, other, middle, 1e-6, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::NotImplemented);
    }
    const auto rational = Intersect(base, middle, Surface{Patch(0.0, true)}, middle, 1e-6, ctx);
    REQUIRE_FALSE(rational.has_value());
    REQUIRE(rational.error().GetCode() == ErrorCode::NotImplemented);
}

TEST_CASE("tolerance prevents premature disjoint classification",
          "[intersect][generic][broadphase]")
{
    const Surface a{Patch(0.0)};
    const Surface b{Patch(2.0)};
    const auto separate = Intersect(a, middle, b, middle, 0.01, ctx);
    REQUIRE(separate.has_value());
    REQUIRE(separate->IsEmpty());
    const auto uncertain = Intersect(a, middle, b, middle, 3.0, ctx);
    REQUIRE_FALSE(uncertain.has_value());
    REQUIRE(uncertain.error().GetCode() == ErrorCode::NotImplemented);
}

TEST_CASE("malformed rational control nets fail rather than certify separation",
          "[intersect][generic][broadphase]")
{
    const Surface good{Patch(10.0)};
    BSplineSurface bad = Patch(0.0, true);
    bad.Weights[0] = -1.0;
    const auto negative = Intersect(Surface{bad}, good, ctx);
    REQUIRE_FALSE(negative.has_value());
    REQUIRE(negative.error().GetCode() == ErrorCode::InvalidInput);

    bad = Patch(0.0);
    bad.Poles[0].Z = std::numeric_limits<double>::quiet_NaN();
    // Use an explicit requested box: a malformed aggregate may not have
    // a valid computed domain, but must never reach the broad-phase proof.
    const auto nan = Intersect(Surface{bad}, middle, good, middle, 0.0, ctx);
    REQUIRE_FALSE(nan.has_value());
    REQUIRE(nan.error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("direct NURBS construction cannot bypass factory weight bounds",
          "[intersect][generic][broadphase][regression]")
{
    BSplineSurface invalid = Patch(0.0, true);
    invalid.Weights[0] = 0x1p257;
    const auto result =
        Intersect(Surface{invalid}, middle, Surface{Patch(100.0)}, middle, 0.0, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
}
