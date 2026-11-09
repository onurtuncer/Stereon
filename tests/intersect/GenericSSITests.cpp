// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Phase 2: a successful empty result requires a mathematical separation
// proof, never simply failure to find a seed or a small iteration budget.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
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

namespace {

/// A biquadratic bowl over the unit square, curved along u: its height is
/// lift + 1 at u = 0 and 1 and lift + 1/2 - dip/2 at u = 1/2, while its net
/// reaches down to lift - dip. The optional weights vary along v only, so
/// they do not move the height.
BSplineSurface Bowl(const double lift, const double dip = 2.0, const bool rational = false)
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            poles.push_back({i / 2.0, j / 2.0, lift + (i == 1 ? -dip : 1.0)});
            weights.push_back(j == 1 ? 2.0 : 1.0);
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles,
                                                  rational ? weights : std::vector<double>{});
    REQUIRE(made.has_value());
    return *made;
}

/// A bilinear sheet over the unit square, rising with u from low to high.
BSplineSurface Slope(const double low, const double high)
{
    const std::vector<Point3> poles{
        {0.0, 0.0, low}, {0.0, 1.0, low}, {1.0, 0.0, high}, {1.0, 1.0, high}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// A bicubic with the given heights, x and y running over [0, 3].
BSplineSurface Bicubic(const std::vector<double>& z, const double lift)
{
    std::vector<Point3> poles;
    for (std::size_t k = 0; k < z.size(); ++k) {
        poles.push_back({static_cast<double>(k / 4), static_cast<double>(k % 4), z[k] + lift});
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(3, 3, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// The smallest distance between grids of points on the two surfaces.
double SampledDistance(const Surface& a, const Surface& b)
{
    constexpr int n = 12;
    double nearest = std::numeric_limits<double>::infinity();
    std::vector<Point3> onA;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            onA.push_back(Eval(a, Stereon::Point2{double(i) / n, double(j) / n}));
        }
    }
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            const Point3 q = Eval(b, Stereon::Point2{double(i) / n, double(j) / n});
            for (const Point3& p : onA) {
                nearest = std::min(nearest, Distance(p, q));
            }
        }
    }
    return nearest;
}

} // namespace

TEST_CASE("subdivision proves apart surfaces whose whole nets overlap",
          "[intersect][generic][subdivision]")
{
    for (const bool rational : {false, true}) {
        const Surface a{Bowl(0.0, 0.5, rational)};
        const Surface b{Bowl(0.5, 0.5, rational)};
        // The nets' boxes share z in [0, 1]; the surfaces stay 0.5 apart in
        // height everywhere.
        const auto forward = Intersect(a, b, ctx);
        const auto reverse = Intersect(b, a, ctx);
        INFO((forward ? std::string() : std::string(forward.error().GetMessage())));
        REQUIRE(forward.has_value());
        REQUIRE(reverse.has_value());
        REQUIRE(forward->IsEmpty());
        REQUIRE(reverse->IsEmpty());
    }
}

TEST_CASE("the requested domains decide what must be separated",
          "[intersect][generic][subdivision]")
{
    const Surface flat{Patch(0.0)};
    const Surface slope{Slope(-0.5, 0.5)}; // Through z = 0 at u = 1/2.
    const Box2 highU{Interval{0.75, 1.0}, Interval{0.0, 1.0}};
    const auto away = Intersect(flat, Box2::Everything(), slope, highU, 1e-6, ctx);
    REQUIRE(away.has_value());
    REQUIRE(away->IsEmpty());
    const auto swapped = Intersect(slope, highU, flat, Box2::Everything(), 1e-6, ctx);
    REQUIRE(swapped.has_value());
    REQUIRE(swapped->IsEmpty());

    const Box2 acrossU{Interval{0.25, 0.75}, Interval{0.0, 1.0}};
    const auto crossing = Intersect(flat, Box2::Everything(), slope, acrossU, 1e-6, ctx);
    REQUIRE_FALSE(crossing.has_value());
    REQUIRE(crossing.error().GetCode() == ErrorCode::NotImplemented);

    // Domains that miss the surface altogether leave nothing to intersect.
    const Box2 outside{Interval{2.0, 3.0}, Interval{0.0, 1.0}};
    const auto none = Intersect(flat, Box2::Everything(), slope, outside, 1e-6, ctx);
    REQUIRE(none.has_value());
    REQUIRE(none->IsEmpty());
}

TEST_CASE("a surface dipping through another is never proved apart",
          "[intersect][generic][subdivision]")
{
    // Lowest height 1/2 - 1.002/2 = -0.001: only a sliver crosses z = 0.
    const Surface bowl{Bowl(0.0, 1.002)};
    const Surface flat{Patch(0.0)};
    const auto result = Intersect(flat, bowl, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::NotImplemented);
}

TEST_CASE("a separation past the search's budget fails rather than runs on",
          "[intersect][generic][subdivision]")
{
    const Surface a{Bowl(0.0)};
    const Surface b{Bowl(1e-9)};
    const auto result = Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-12, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a proved empty result agrees with sampling",
          "[intersect][generic][subdivision][property]")
{
    rc::prop("for random bicubic pairs, lifted by a varying amount", [] {
        std::vector<double> z(16);
        std::vector<double> w(16);
        for (std::size_t k = 0; k < 16; ++k) {
            z[k] = *rc::gen::inRange(-16, 17) / 16.0;
            w[k] = *rc::gen::inRange(-16, 17) / 16.0;
        }
        const double lift = *rc::gen::inRange(0, 25) / 8.0;
        const Surface a{Bicubic(z, 0.0)};
        const Surface b{Bicubic(w, lift)};
        const auto result = Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-6, ctx);
        if (result) {
            RC_ASSERT(result->IsEmpty());
            RC_ASSERT(SampledDistance(a, b) > 1e-6);
        } else {
            RC_ASSERT(result.error().GetCode() == ErrorCode::NotImplemented ||
                      result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    });
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
