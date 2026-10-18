// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <string>

import Stereon.Core;

using Stereon::ErrorCode;

namespace {

constexpr double Inf = std::numeric_limits<double>::infinity();
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

void RequireInvalid(const Stereon::Result<void>& result, const std::string& fragment)
{
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(std::string(result.error().GetMessage()).find(fragment) != std::string::npos);
}

} // namespace

TEST_CASE("CheckNonNegative accepts zero and refuses negatives and non-finite values",
          "[core][checks]")
{
    using Stereon::Detail::CheckNonNegative;
    REQUIRE(CheckNonNegative(0.0, "the tolerance", "Here").has_value());
    REQUIRE(CheckNonNegative(-0.0, "the tolerance", "Here").has_value());
    REQUIRE(CheckNonNegative(1e-300, "the tolerance", "Here").has_value());
    RequireInvalid(CheckNonNegative(-1e-300, "the tolerance", "Here"), "the tolerance");
    RequireInvalid(CheckNonNegative(NaN, "the tolerance", "Here"), "not negative");
    RequireInvalid(CheckNonNegative(Inf, "the tolerance", "Here"), "finite");
    RequireInvalid(CheckNonNegative(-Inf, "the tolerance", "Here"), "Here");
}

TEST_CASE("CheckPositive refuses zero where CheckNonNegative accepts it", "[core][checks]")
{
    REQUIRE_FALSE(Stereon::Detail::CheckPositive(0.0, "the radius", "Here").has_value());
    REQUIRE(Stereon::Detail::CheckPositive(1.0, "the radius", "Here").has_value());
    REQUIRE(Stereon::Detail::CheckNonNegative(0.0, "the radius", "Here").has_value());
}

TEST_CASE("CheckFinite refuses a point with any non-finite coordinate", "[core][checks]")
{
    using Stereon::Point3;
    REQUIRE(Stereon::Detail::CheckFinite(Point3{1, 2, 3}, "the point", "Here").has_value());
    RequireInvalid(Stereon::Detail::CheckFinite(Point3{1, NaN, 3}, "the point", "Here"),
                   "must be finite");
    RequireInvalid(Stereon::Detail::CheckFinite(Point3{1, 2, Inf}, "the point", "Here"),
                   "the point");
}
