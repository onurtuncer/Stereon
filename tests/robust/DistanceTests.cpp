// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
import Stereon.Robust;

TEST_CASE("Distance comparison resolves boundaries and rescales extreme magnitudes",
          "[robust][distance]")
{
    using namespace Stereon;
    for (const double scale : {1.0, 0x1p-1020, 0x1p900}) {
        const auto exact = CompareDistance({3 * scale, 4 * scale, 0}, {0, 0, 0}, 5 * scale);
        REQUIRE(exact.Supported);
        REQUIRE(exact.Value == Sign::Zero);
        REQUIRE(
            CompareDistance({3 * scale, 4 * scale, 0}, {0, 0, 0}, std::nextafter(5 * scale, 0.0))
                .Value == Sign::Positive);
        REQUIRE(CompareDistance({3 * scale, 4 * scale, 0}, {0, 0, 0},
                                std::nextafter(5 * scale, std::numeric_limits<double>::infinity()))
                    .Value == Sign::Negative);
    }
    REQUIRE_FALSE(CompareDistance({0, 0, 0}, {0, 0, 0}, -1).Supported);
    REQUIRE_FALSE(
        CompareDistance({0, 0, 0}, {0, 0, 0}, std::numeric_limits<double>::infinity()).Supported);
    REQUIRE_FALSE(CompareDistance({0x1p900, 0x1p-900, 0}, {0, 0, 0}, 1).Supported);
}
