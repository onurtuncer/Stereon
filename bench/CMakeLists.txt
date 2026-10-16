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
target_link_libraries(stereon_bench
    PRIVATE stereon::core stereon::robust stereon::surface stereon_options benchmark::benchmark_main)
stereon_anchor_contracts(stereon_bench)

# The Phase 1 comparison with OpenCASCADE (docs/benchmarks.md, B1.1-B1.5):
# built when OCCT is available, -DSTEREON_WITH_OCCT=ON with OCCT's install
# directory in CMAKE_PREFIX_PATH (tools/occt/build_occt.sh). OCCT is
# LGPL-2.1 and linked into this executable only (ADR-0008).
if(STEREON_WITH_OCCT)
    find_package(OpenCASCADE CONFIG REQUIRED)
    add_executable(stereon_occt_bench occt/OcctBench.cpp)
    target_link_libraries(stereon_occt_bench
        PRIVATE stereon::surface stereon_options benchmark::benchmark_main
                TKernel TKMath TKG3d TKGeomBase TKGeomAlgo)
    stereon_anchor_contracts(stereon_occt_bench)
endif()
