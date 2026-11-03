// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Phase 2: a successful empty result requires a mathematical separation
// proof, never simply failure to find a seed or a small iteration budget.

#include <catch2/catch_test_macros.hpp>

#include <array>
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
    const std::vector<Point3> poles{{0.0, 0.0, z}, {0.0, 1.0, z},
                                     {1.0, 0.0, z}, {1.0, 1.0, z}};
    const std::vector<double> weights =
        rational ? std::vector<double>{1.0, 3.0, 2.0, 4.0} : std::vector<double>{};
    const auto made = Stereon::MakeBSplineSurface(
        1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0}, poles, weights);
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
            for (const auto& [boxA, boxB] :
                 {std::pair{Box2::Everything(), Box2::Everything()},
                  std::pair{middle, middle}}) {
                const auto forward = Intersect(a, boxA, b, boxB, 1e-6, ctx);
                const auto reverse = Intersect(b, boxB, a, boxA, 1e-6, ctx);
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
    for (const Surface other : {Surface{Patch(0.0)}, Surface{Patch(1e-7)},
                                Surface{Patch(0.0, true)}}) {
        const auto result = Intersect(base, middle, other, middle, 1e-6, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::NotImplemented);
    }
}

TEST_CASE("tolerance prevents premature disjoint classification",
          "[intersect][generic][broadphase]")
{
    const Surface a{Patch(0.0)};
    const Surface b{Patch(2.0)};
    REQUIRE(Intersect(a, middle, b, middle, 0.01, ctx)->IsEmpty());
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
