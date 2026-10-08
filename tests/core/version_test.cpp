// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

import stereon.core;

TEST_CASE("version matches the CMake project version", "[core][version]")
{
    STATIC_REQUIRE(stereon::version == stereon::Version{STEREON_EXPECTED_VERSION_MAJOR,
                                                       STEREON_EXPECTED_VERSION_MINOR,
                                                       STEREON_EXPECTED_VERSION_PATCH});
}

TEST_CASE("checked handles follow the build configuration", "[core][handles]")
{
#if defined(STEREON_CHECKED_HANDLES)
    STATIC_REQUIRE(stereon::checked_handles);
#else
    STATIC_REQUIRE_FALSE(stereon::checked_handles);
#endif
}
