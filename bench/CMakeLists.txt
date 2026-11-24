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
    PRIVATE stereon::core stereon::robust stereon::surface stereon::topo stereon::boolean
            stereon::intersect stereon_options benchmark::benchmark_main
            $<$<PLATFORM_ID:Windows>:psapi>)
# The Boolean benchmark shares its fixture with the tests (BooleanCases.hpp);
# the intersection benchmarks (B2.1-B2.3, IntersectBench.cpp) read the corpus
# with the tests' JSON reader.
target_include_directories(stereon_bench PRIVATE "${PROJECT_SOURCE_DIR}/tests/support")
target_compile_definitions(stereon_bench PRIVATE STEREON_CORPUS_DIR="${PROJECT_SOURCE_DIR}/corpus")
stereon_anchor_contracts(stereon_bench)

# B1.3's convergence on a declared corpus (bench/convergence/ProjectionCorpus.cpp):
# a measurement with its own main, so not one of stereon_bench's benchmarks. The
# OCCT reference workflow runs it and hands its JSON to the release-gate check.
add_executable(stereon_b13_convergence convergence/ProjectionCorpus.cpp)
target_link_libraries(stereon_b13_convergence PRIVATE stereon::surface stereon_options)
stereon_anchor_contracts(stereon_b13_convergence)

# The comparison with OpenCASCADE (docs/benchmarks.md, B0.3-B0.5, B1.1-B1.6
# and B2.1-B2.3): built when OCCT is available, -DSTEREON_WITH_OCCT=ON with
# OCCT's install directory in CMAKE_PREFIX_PATH (tools/occt/build_occt.sh).
# The topology rows' Stereon twins are TopoBench.cpp and the intersection
# rows' IntersectBench.cpp, compiled in here too so that both libraries run
# interleaved in one process. OCCT is LGPL-2.1 and
# linked into this executable only (ADR-0008).
if(STEREON_WITH_OCCT)
    find_package(OpenCASCADE CONFIG REQUIRED)
    add_executable(stereon_occt_bench occt/OcctBench.cpp occt/OcctTopoBench.cpp TopoBench.cpp
                                      occt/OcctIntersectBench.cpp IntersectBench.cpp Memory.cpp)
    target_link_libraries(stereon_occt_bench
        PRIVATE stereon::surface stereon::topo stereon::intersect stereon_options
                benchmark::benchmark_main TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep
                TKTopAlgo $<$<PLATFORM_ID:Windows>:psapi>)
    target_include_directories(stereon_occt_bench PRIVATE "${PROJECT_SOURCE_DIR}/tests/support")
    target_compile_definitions(stereon_occt_bench
        PRIVATE STEREON_CORPUS_DIR="${PROJECT_SOURCE_DIR}/corpus")
    stereon_anchor_contracts(stereon_occt_bench)
    # GCC's -Wnull-dereference fires inside OCCT's NCollection_Array1 accessors
    # once they are inlined at -O3, where -isystem does not reach, as it does
    # in stdexec's queue (stereon_exec_sources()); with -Werror that stops the
    # reference run's build.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set_property(SOURCE occt/OcctBench.cpp occt/OcctTopoBench.cpp occt/OcctIntersectBench.cpp
                     APPEND PROPERTY
                     COMPILE_OPTIONS -Wno-null-dereference)
    endif()
endif()
