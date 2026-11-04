# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# The worked examples of docs/examples.rst, which includes these files
# verbatim. They are built with the tests and run by ctest, so the manual's
# code always compiles and runs; each returns non-zero when its own check
# fails.
if(NOT STEREON_BUILD_TESTS)
    return()
endif()

foreach(example curves intersection boolean geometry_file)
    add_executable(stereon_example_${example} ${example}.cpp)
    target_link_libraries(stereon_example_${example} PRIVATE stereon::io stereon_options)
    stereon_anchor_contracts(stereon_example_${example})
    add_test(NAME "example: ${example}" COMMAND stereon_example_${example})
endforeach()
