# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include(Catch)

# One test executable per library: tests/<library>/*.cpp -> stereon_<library>_tests.
# Test cases are tagged with their library, e.g. [core], plus topic tags.
function(stereon_add_tests library)
    file(GLOB sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${library}/*.cpp")
    add_executable(stereon_${library}_tests ${sources})
    target_link_libraries(stereon_${library}_tests
        PRIVATE stereon::${library} stereon_options Catch2::Catch2WithMain)
    target_compile_definitions(stereon_${library}_tests PRIVATE
        STEREON_EXPECTED_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        STEREON_EXPECTED_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        STEREON_EXPECTED_VERSION_PATCH=${PROJECT_VERSION_PATCH})
    catch_discover_tests(stereon_${library}_tests
        TEST_PREFIX "${library}: "
        DISCOVERY_MODE PRE_TEST)
endfunction()

stereon_add_tests(core)
stereon_add_tests(robust)
