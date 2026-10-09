# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Third-party dependencies, fetched at configure time and pinned to a release
# tag (or a commit where the project has no releases). Every dependency must
# have an MPL-2.0-compatible licence (ADR-0008); record it next to its entry.
# Dependencies are SYSTEM, so their headers do not trip Stereon's warnings.

include_guard(GLOBAL)
include(FetchContent)

# --- Library polyfills (ADR-0016) -------------------------------------------
# Always fetched: the Stereon::Exec and Stereon::Simd wrappers pick the
# standard library or the polyfill with feature-test macros at compile time.

find_package(Threads REQUIRED)

# stdexec: Apache-2.0 WITH LLVM-exception. Header-only. Its own CMake build
# downloads rapids-cmake at configure time, so only the sources are fetched
# (SOURCE_SUBDIR points at a directory without a CMakeLists.txt) and the
# target is defined here. The parallel scheduler's default backend is
# compiled inline, so nothing needs linking besides threads.
FetchContent_Declare(stdexec
    GIT_REPOSITORY https://github.com/NVIDIA/stdexec.git
    GIT_TAG        nvhpc-26.05
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  stereon-headers-only)
FetchContent_MakeAvailable(stdexec)
add_library(stereon_stdexec INTERFACE)
add_library(stereon::stdexec ALIAS stereon_stdexec)
target_include_directories(stereon_stdexec SYSTEM INTERFACE "${stdexec_SOURCE_DIR}/include")
target_compile_definitions(stereon_stdexec INTERFACE STDEXEC_PARALLEL_SCHEDULER_HEADER_ONLY)
target_link_libraries(stereon_stdexec INTERFACE Threads::Threads)

# xsimd: BSD-3-Clause. Header-only.
FetchContent_Declare(xsimd
    GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
    GIT_TAG        14.3.0
    GIT_SHALLOW    TRUE
    SYSTEM)
FetchContent_MakeAvailable(xsimd)

# --- Tests -------------------------------------------------------------------

if(STEREON_BUILD_TESTS)
    # Catch2: BSL-1.0
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.1
        GIT_SHALLOW    TRUE
        SYSTEM
        FIND_PACKAGE_ARGS 3)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")

    # RapidCheck: BSD-2-Clause. No release tags upstream, so a commit is pinned.
    # Its RC_ENABLE_CATCH switch also builds a bundled Catch 2.x that clashes
    # with the Catch2 target above, so the header-only Catch adapter (rc::prop)
    # is defined here instead, and the bundled submodules are not cloned.
    set(RC_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(RC_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(rapidcheck
        GIT_REPOSITORY https://github.com/emil-e/rapidcheck.git
        GIT_TAG        2c3c4365aca21ef4e612768fdc95b1ce8b39a651
        GIT_SUBMODULES ""
        SYSTEM)
    FetchContent_MakeAvailable(rapidcheck)
    add_library(stereon_rapidcheck_catch INTERFACE)
    target_include_directories(stereon_rapidcheck_catch SYSTEM INTERFACE
        "${rapidcheck_SOURCE_DIR}/extras/catch/include")
    target_link_libraries(stereon_rapidcheck_catch INTERFACE rapidcheck Catch2::Catch2)
endif()

# --- Benchmarks --------------------------------------------------------------

if(STEREON_BUILD_BENCH)
    # Google Benchmark: Apache-2.0
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        v1.9.5
        GIT_SHALLOW    TRUE
        SYSTEM)
    FetchContent_MakeAvailable(benchmark)
endif()

# --- Python bindings ---------------------------------------------------------

if(STEREON_BUILD_PYTHON)
    find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
    # nanobind: BSD-3-Clause
    FetchContent_Declare(nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG        v3.1.0
        GIT_SHALLOW    TRUE
        SYSTEM)
    FetchContent_MakeAvailable(nanobind)
endif()
