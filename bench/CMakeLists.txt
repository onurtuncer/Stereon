# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Micro-benchmarks (Google Benchmark). The OCCT comparison harness from
# docs/benchmarks.md will live here too. Build with the bench preset.
file(GLOB sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable(stereon_bench ${sources})
target_link_libraries(stereon_bench PRIVATE stereon::core stereon_options benchmark::benchmark_main)
