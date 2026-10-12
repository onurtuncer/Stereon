// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <format>
#include <optional>

import Stereon.Robust;

#include "RobustFormatting.hpp"

using Stereon::FilteredSign;
using Stereon::Interval;
using Stereon::Sign;
using Stereon::SignOf;

TEST_CASE("signs negate, convert and format", "[robust][sign]")
{
    STATIC_REQUIRE(SignOf(2.5) == Sign::Positive);
    STATIC_REQUIRE(SignOf(-0.0) == Sign::Zero);
    STATIC_REQUIRE(SignOf(-1e-300) == Sign::Negative);
    STATIC_REQUIRE(-Sign::Positive == Sign::Negative);
    STATIC_REQUIRE(-Sign::Zero == Sign::Zero);
    STATIC_REQUIRE(Stereon::ToInt(Sign::Negative) == -1);
    STATIC_REQUIRE(Stereon::ToInt(Sign::Positive) == 1);
    if constexpr (Stereon::EnumNames) {
        REQUIRE(std::format("{}", Sign::Positive) == "Positive");
        REQUIRE(Stereon::ToString(Sign::Zero) == "Zero");
    }
}

TEST_CASE("an interval decides the sign when it excludes zero", "[robust][filter]")
{
    STATIC_REQUIRE(SignOf(Interval{1.0, 2.0}) == Sign::Positive);
    STATIC_REQUIRE(SignOf(Interval{-2.0, -1e-300}) == Sign::Negative);
    STATIC_REQUIRE(SignOf(Interval{0.0, 0.0}) == Sign::Zero);
    STATIC_REQUIRE(SignOf(Interval{-1.0, 1.0}) == std::nullopt);
    STATIC_REQUIRE(SignOf(Interval{0.0, 1.0}) == std::nullopt);
    STATIC_REQUIRE(SignOf(Interval{-1.0, 0.0}) == std::nullopt);
    STATIC_REQUIRE(SignOf(Interval{}) == std::nullopt);
    STATIC_REQUIRE(SignOf(Interval::Everything()) == std::nullopt);

    // Products of intervals that straddle zero fall through to the exact
    // computation; the filter never calls it otherwise.
    int exactCalls = 0;
    const auto exact = [&exactCalls] {
        ++exactCalls;
        return Sign::Negative;
    };
    REQUIRE(FilteredSign(Interval{3.0, 4.0} * Interval{1.0, 2.0}, exact) == Sign::Positive);
    REQUIRE(exactCalls == 0);
    REQUIRE(FilteredSign(Interval{-1.0, 1.0} * Interval{1.0, 2.0}, exact) == Sign::Negative);
    REQUIRE(exactCalls == 1);
    REQUIRE(FilteredSign(Interval{1.0, 1.0} - Interval{1.0, 1.0}, exact) == Sign::Negative);
    REQUIRE(exactCalls == 2);
}
