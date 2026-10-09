// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <Stereon/Core/Contracts.hpp>

#include <catch2/catch_test_macros.hpp>

import Stereon.Core;

namespace {

// Uses every contract macro. Parameters named in a postcondition must be const (P2900).
int Halve(const int value) STEREON_PRE(value % 2 == 0) STEREON_HOT_PRE(value >= 0)
    STEREON_POST(result : result * 2 == value)
{
    STEREON_ASSERT(value < 1'000'000);
    STEREON_HOT_ASSERT(value != 1);
    return value / 2;
}

} // namespace

TEST_CASE("contract macros compile and pass on valid input", "[core][contracts]")
{
    REQUIRE(Halve(8) == 4);
    REQUIRE(Halve(0) == 0);
}

TEST_CASE("hot-path contracts follow the build configuration", "[core][contracts]")
{
    STATIC_REQUIRE(Stereon::HotContracts == (STEREON_HOT_CONTRACTS != 0));
}
