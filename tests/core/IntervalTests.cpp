// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <format>
#include <limits>
#include <numbers>

import Stereon.Core;

#include "Formatting.hpp"

using Stereon::Interval;

namespace {

constexpr double Inf = std::numeric_limits<double>::infinity();

// Multiples of 1/64 up to 16 in magnitude. Sums, differences and products of
// two of them are exact, so the real result of an operation is also its
// double result and enclosure can be tested exactly.
rc::Gen<double> SmallDyadic()
{
    return rc::gen::map(rc::gen::inRange(-1024, 1025),
                        [](const int value) { return static_cast<double>(value) / 64.0; });
}

rc::Gen<Interval> SmallInterval()
{
    return rc::gen::map(rc::gen::pair(SmallDyadic(), SmallDyadic()), [](const auto& bounds) {
        return Interval{std::min(bounds.first, bounds.second),
                        std::max(bounds.first, bounds.second)};
    });
}

// A dyadic point inside i, exactly.
double PointIn(const Interval& i)
{
    const double t = static_cast<double>(*rc::gen::inRange(0, 65)) / 64.0;
    return i.GetLower() + (i.GetUpper() - i.GetLower()) * t;
}

} // namespace

TEST_CASE("intervals are closed and may be empty", "[core][interval]")
{
    STATIC_REQUIRE(Interval{}.IsEmpty());
    STATIC_REQUIRE(Interval::Empty() == Interval{});
    STATIC_REQUIRE_FALSE(Interval{1, 2}.IsEmpty());
    STATIC_REQUIRE(Interval{1, 2}.GetLower() == 1.0);
    STATIC_REQUIRE(Interval{1, 2}.GetUpper() == 2.0);
    STATIC_REQUIRE(Interval{3.0} == Interval{3, 3});

    STATIC_REQUIRE(Interval{1, 2}.IsBounded());
    STATIC_REQUIRE(Interval{}.IsBounded());
    STATIC_REQUIRE_FALSE(Interval::Everything().IsBounded());
    STATIC_REQUIRE_FALSE(Interval{0, Inf}.IsBounded());

    STATIC_REQUIRE(Contains(Interval{1, 2}, 1.0));
    STATIC_REQUIRE(Contains(Interval{1, 2}, 2.0));
    STATIC_REQUIRE_FALSE(Contains(Interval{1, 2}, 2.5));
    STATIC_REQUIRE_FALSE(Contains(Interval{}, 0.0));
    STATIC_REQUIRE(Contains(Interval::Everything(), 1e300));

    STATIC_REQUIRE(Contains(Interval{0, 3}, Interval{1, 2}));
    STATIC_REQUIRE_FALSE(Contains(Interval{0, 3}, Interval{1, 4}));
    STATIC_REQUIRE(Contains(Interval{0, 3}, Interval{}));
    STATIC_REQUIRE(Contains(Interval{}, Interval{}));
}

TEST_CASE("hull, intersection and overlap", "[core][interval]")
{
    STATIC_REQUIRE(Hull(Interval{0, 1}, Interval{2, 3}) == Interval{0, 3});
    STATIC_REQUIRE(Hull(Interval{}, Interval{2, 3}) == Interval{2, 3});
    STATIC_REQUIRE(Hull(Interval{2, 3}, Interval{}) == Interval{2, 3});
    STATIC_REQUIRE(Hull(Interval{}, Interval{}).IsEmpty());
    STATIC_REQUIRE(Hull(Interval{1, 2}, 5.0) == Interval{1, 5});
    STATIC_REQUIRE(Hull(Interval{}, 5.0) == Interval{5.0});

    STATIC_REQUIRE(Intersection(Interval{0, 2}, Interval{1, 3}) == Interval{1, 2});
    STATIC_REQUIRE(Intersection(Interval{0, 1}, Interval{1, 2}) == Interval{1.0});
    // Disjoint intervals give the canonical empty interval.
    STATIC_REQUIRE(Intersection(Interval{0, 1}, Interval{2, 3}) == Interval{});
    STATIC_REQUIRE(Intersection(Interval{0, 1}, Interval{}) == Interval{});

    STATIC_REQUIRE(Intersects(Interval{0, 1}, Interval{1, 2}));
    STATIC_REQUIRE_FALSE(Intersects(Interval{0, 1}, Interval{2, 3}));
    STATIC_REQUIRE_FALSE(Intersects(Interval{}, Interval::Everything()));
}

TEST_CASE("the nearest image of a value in a periodic range", "[core][interval]")
{
    using Stereon::NearestInRange;
    const Interval range{0.2, 0.6};
    constexpr double tau = 2.0 * std::numbers::pi;
    // Not periodic: clamped.
    REQUIRE(NearestInRange(0.7, range, 0.0) == 0.6);
    REQUIRE(NearestInRange(-1.0, range, 0.0) == 0.2);
    REQUIRE(NearestInRange(0.3, range, 0.0) == 0.3);
    // Inside, periods away: brought in.
    REQUIRE(NearestInRange(0.3, range, tau) == 0.3);
    REQUIRE(NearestInRange(0.3 + tau, range, tau) == Catch::Approx(0.3).margin(1e-12));
    REQUIRE(NearestInRange(0.3 - 3.0 * tau, range, tau) == Catch::Approx(0.3).margin(1e-12));
    // A rounding error past an end stays at that end, whichever end and
    // however many periods away.
    REQUIRE(NearestInRange(0.6 + 1e-15, range, tau) == 0.6);
    REQUIRE(NearestInRange(0.2 - 1e-15, range, tau) == 0.2);
    REQUIRE(NearestInRange(0.6 + 1e-15 + tau, range, tau) == 0.6);
    REQUIRE(NearestInRange(0.2 - 1e-15 - tau, range, tau) == 0.2);
    // Well past an end: still the nearer end round the period.
    REQUIRE(NearestInRange(1.6, range, tau) == 0.6);
    REQUIRE(NearestInRange(-0.8, range, tau) == 0.2);
    REQUIRE(NearestInRange(-0.8 + 2.0 * tau, range, tau) == 0.2);
    // A range of a whole period: everything has an image inside.
    const Interval turn{0.0, tau};
    REQUIRE(NearestInRange(tau + 0.1, turn, tau) == Catch::Approx(0.1).margin(1e-12));
    REQUIRE(NearestInRange(-0.1, turn, tau) == Catch::Approx(tau - 0.1).margin(1e-12));
}

TEST_CASE("width, midpoint and inflation", "[core][interval]")
{
    STATIC_REQUIRE(Width(Interval{1, 4}) == 3.0);
    STATIC_REQUIRE(Width(Interval{0, Inf}) == Inf);
    STATIC_REQUIRE(Midpoint(Interval{1, 4}) == 2.5);
    STATIC_REQUIRE(Midpoint(Interval{-1e308, 1e308}) == 0.0);

    STATIC_REQUIRE(Inflated(Interval{1, 2}, 0.5) == Interval{0.5, 2.5});
    STATIC_REQUIRE(Inflated(Interval{1, 2}, -0.25) == Interval{1.25, 1.75});
    STATIC_REQUIRE(Inflated(Interval{1, 2}, -1.0) == Interval{});
    STATIC_REQUIRE(Inflated(Interval{}, 1.0) == Interval{});
}

TEST_CASE("arithmetic encloses the real result", "[core][interval][adr-0002]")
{
    // 0.1 + 0.2 is not 0.3 in double arithmetic, but the interval sum
    // contains it.
    const Interval sum = Interval{0.1} + Interval{0.2};
    REQUIRE(0.1 + 0.2 != 0.3);
    REQUIRE(Contains(sum, 0.3));
    REQUIRE(Contains(sum, 0.1 + 0.2));
    REQUIRE(Width(sum) < 1e-15);

    // Each bound moves exactly one representable value outward; this also
    // checks the bit-level stepping against the standard library.
    constexpr Interval widened = Interval{1, 2} + Interval{1, 1};
    STATIC_REQUIRE(widened.GetLower() < 2.0);
    STATIC_REQUIRE(widened.GetUpper() > 3.0);
    REQUIRE(widened.GetLower() == std::nextafter(2.0, -Inf));
    REQUIRE(widened.GetUpper() == std::nextafter(3.0, Inf));

    STATIC_REQUIRE(-Interval{1, 2} == Interval{-2, -1});
    constexpr Interval difference = Interval{1, 2} - Interval{1, 2};
    STATIC_REQUIRE(Contains(difference, -1.0));
    STATIC_REQUIRE(Contains(difference, 1.0));

    // Products take the extreme corner products: [-1, 2] x [3, 4] is [-4, 8].
    constexpr Interval product = Interval{-1, 2} * Interval{3, 4};
    STATIC_REQUIRE(Contains(product, -4.0));
    STATIC_REQUIRE(Contains(product, 8.0));
    STATIC_REQUIRE(product.GetLower() > -4.0 - 1e-14);
    STATIC_REQUIRE(product.GetUpper() < 8.0 + 1e-14);

    // 0 times anything, even an unbounded interval, is 0.
    constexpr Interval zero = Interval::Everything() * Interval{0.0};
    STATIC_REQUIRE(Contains(zero, 0.0));
    STATIC_REQUIRE(Width(zero) < 1e-300);

    // Division runs through the reciprocal; 1 / inf is 0.
    constexpr Interval quotient = Interval{1, 2} / Interval{2, 4};
    STATIC_REQUIRE(Contains(quotient, 0.25));
    STATIC_REQUIRE(Contains(quotient, 1.0));
    STATIC_REQUIRE(Width(quotient) < 0.75 + 1e-14);
    constexpr Interval tail = Interval{1, 2} / Interval{1, Inf};
    STATIC_REQUIRE(Contains(tail, 0.0));
    STATIC_REQUIRE(Contains(tail, 2.0));
    STATIC_REQUIRE(tail.IsBounded());
}

TEST_CASE("squaring and square roots", "[core][interval][adr-0002]")
{
    // Squared knows both factors are the same value; i * i does not.
    constexpr Interval squared = Squared(Interval{-1, 2});
    STATIC_REQUIRE(squared.GetLower() == 0.0);
    STATIC_REQUIRE(Contains(squared, 4.0));
    STATIC_REQUIRE(squared.GetUpper() < 4.0 + 1e-14);
    STATIC_REQUIRE((Interval{-1, 2} * Interval{-1, 2}).GetLower() < -1.9);

    STATIC_REQUIRE(Squared(Interval{-3, -2}).GetLower() > 3.9);
    STATIC_REQUIRE(Contains(Squared(Interval::Everything()), 1e300));

    const Interval root = Sqrt(Interval{4, 9});
    REQUIRE(Contains(root, 2.0));
    REQUIRE(Contains(root, 3.0));
    REQUIRE(Width(root) < 1.0 + 1e-14);
    REQUIRE(Sqrt(Interval{0, 1}).GetLower() == 0.0);
    REQUIRE(Sqrt(Interval{}).IsEmpty());
}

TEST_CASE("empty intervals absorb arithmetic", "[core][interval]")
{
    STATIC_REQUIRE((Interval{} + Interval{1, 2}).IsEmpty());
    STATIC_REQUIRE((Interval{1, 2} - Interval{}).IsEmpty());
    STATIC_REQUIRE((Interval{} * Interval::Everything()).IsEmpty());
    STATIC_REQUIRE((Interval{} / Interval{1, 2}).IsEmpty());
    STATIC_REQUIRE((-Interval{}).IsEmpty());
    STATIC_REQUIRE(Squared(Interval{}).IsEmpty());
}

TEST_CASE("intervals format as brackets", "[core][interval]")
{
    REQUIRE(std::format("{}", Interval{1, 2.5}) == "[1, 2.5]");
    REQUIRE(std::format("{}", Interval{}) == "[]");
    REQUIRE(std::format("{}", Interval::Everything()) == "[-inf, inf]");
}

TEST_CASE("interval arithmetic encloses every real result", "[core][interval][adr-0002][property]")
{
    rc::prop("sums, differences and products", [] {
        const Interval a = *SmallInterval();
        const Interval b = *SmallInterval();
        const double x = PointIn(a);
        const double y = PointIn(b);
        RC_ASSERT(Contains(a + b, x + y));
        RC_ASSERT(Contains(a - b, x - y));
        RC_ASSERT(Contains(a * b, x * y));
        RC_ASSERT(Contains(Squared(a), x * x));
        RC_ASSERT(Contains(-a, -x));
    });

    rc::prop("quotients and square roots", [] {
        const Interval a = *SmallInterval();
        const Interval b = *SmallInterval();
        const double x = PointIn(a);
        const double y = PointIn(b);
        RC_PRE(!Contains(b, 0.0));
        RC_ASSERT(Contains(a / b, x / y));
        if (a.GetLower() >= 0.0) {
            RC_ASSERT(Contains(Sqrt(a), std::sqrt(x)));
        }
    });

    rc::prop("hull and intersection are the set operations", [] {
        const Interval a = *SmallInterval();
        const Interval b = *SmallInterval();
        const double x = PointIn(a);
        const double z = *SmallDyadic();
        RC_ASSERT(Contains(Hull(a, b), x));
        RC_ASSERT(Contains(Hull(a, b), PointIn(b)));
        RC_ASSERT(Contains(Intersection(a, b), z) == (Contains(a, z) && Contains(b, z)));
        RC_ASSERT(Intersects(a, b) == !Intersection(a, b).IsEmpty());
    });
}
