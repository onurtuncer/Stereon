# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Third-party dependencies, fetched at configure time. Every dependency must
# have an MPL-2.0-compatible licence (ADR-0008); record it next to its entry.
# Further dependencies (stdexec, xsimd, RapidCheck, Google Benchmark,
# nanobind) are added here when the first code that needs them lands.

include_guard(GLOBAL)
include(FetchContent)

if(STEREON_BUILD_TESTS)
    # Catch2 — BSL-1.0
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.1
        GIT_SHALLOW    TRUE
        SYSTEM
        FIND_PACKAGE_ARGS 3)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
