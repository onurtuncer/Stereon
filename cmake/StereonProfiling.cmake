# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include_guard(GLOBAL)

# PRIVATE links to this target provide the adapter to implementation units.
# Static libraries propagate its link dependencies, not its compile interface.
add_library(stereon_profiling INTERFACE)
target_include_directories(stereon_profiling INTERFACE
    "${CMAKE_CURRENT_LIST_DIR}/../src/core/private")
target_compile_definitions(stereon_profiling INTERFACE
    STEREON_ENABLE_PROFILING=$<BOOL:${STEREON_ENABLE_PROFILING}>)
if(STEREON_ENABLE_PROFILING)
    target_link_libraries(stereon_profiling INTERFACE Tracy::TracyClient)
endif()
