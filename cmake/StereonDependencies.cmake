# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Third-party dependencies, fetched at configure time and pinned to a release
# tag (or a commit where the project has no releases). Every dependency must
# have an MPL-2.0-compatible licence (ADR-0008), recorded above its entry as
# "# Dependency licence: <SPDX expression>"; tools/check_dependencies.py
# checks the expression against the allow-list and, in CTest, against the
# licence file actually fetched.
# Dependencies are SYSTEM, so their headers do not trip Stereon's warnings.

include_guard(GLOBAL)
include(FetchContent)

# --- Library polyfills (ADR-0016) -------------------------------------------
# Always fetched: the Stereon::Exec and Stereon::Simd wrappers pick the
# standard library or the polyfill with feature-test macros at compile time.

find_package(Threads REQUIRED)

# stdexec. Header-only. Its own CMake build
# downloads rapids-cmake at configure time, so only the sources are fetched
# (SOURCE_SUBDIR points at a directory without a CMakeLists.txt) and the
# target is defined here. The parallel scheduler's default backend is
# compiled inline, so nothing needs linking besides threads.
# Dependency licence: Apache-2.0 WITH LLVM-exception
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

# xsimd. Header-only.
# Dependency licence: BSD-3-Clause
FetchContent_Declare(xsimd
    GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
    GIT_TAG        14.3.0
    GIT_SHALLOW    TRUE
    SYSTEM)
FetchContent_MakeAvailable(xsimd)

# --- Tests -------------------------------------------------------------------

if(STEREON_BUILD_TESTS)
    # Catch2.
    # Dependency licence: BSL-1.0
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.1
        GIT_SHALLOW    TRUE
        SYSTEM
        FIND_PACKAGE_ARGS 3)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")

    # RapidCheck. No release tags upstream, so a commit is pinned.
    # Its RC_ENABLE_CATCH switch also builds a bundled Catch 2.x that clashes
    # with the Catch2 target above, so the header-only Catch adapter (rc::prop)
    # is defined here instead, and the bundled submodules are not cloned.
    set(RC_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(RC_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    # Dependency licence: BSD-2-Clause
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
    # Google Benchmark.
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
    # Dependency licence: Apache-2.0
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
    # nanobind.
    # Dependency licence: BSD-3-Clause
    FetchContent_Declare(nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG        v3.1.0
        GIT_SHALLOW    TRUE
        SYSTEM)
    FetchContent_MakeAvailable(nanobind)
endif()

# --- Debug viewer --------------------------------------------------------------
# Only stn-view (viewer/app) may link these (ADR-0035, ADR-0034):
# tools/check_dependencies.py rejects a reference to them anywhere else.

if(STEREON_BUILD_VIEWER)
    # GLFW, for the window, input and the OpenGL context. X11 only on Linux,
    # so that the build needs the X11 development packages and not Wayland's.
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "" FORCE)
    # Dependency licence: Zlib
    FetchContent_Declare(glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG        3.4
        GIT_SHALLOW    TRUE
        SYSTEM)
    FetchContent_MakeAvailable(glfw)

    # Dear ImGui, for the panels, with its GLFW and OpenGL 3 backends. It has
    # no CMake build of its own, so the sources are fetched only and the
    # library defined here; its OpenGL 3 backend brings its own loader.
    # Dependency licence: MIT
    FetchContent_Declare(imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG        v1.91.8
        GIT_SHALLOW    TRUE
        SOURCE_SUBDIR  stereon-sources-only)
    FetchContent_MakeAvailable(imgui)
    add_library(stereon_imgui STATIC
        "${imgui_SOURCE_DIR}/imgui.cpp"
        "${imgui_SOURCE_DIR}/imgui_draw.cpp"
        "${imgui_SOURCE_DIR}/imgui_tables.cpp"
        "${imgui_SOURCE_DIR}/imgui_widgets.cpp"
        "${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp"
        "${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp")
    set_target_properties(stereon_imgui PROPERTIES CXX_SCAN_FOR_MODULES OFF)
    target_include_directories(stereon_imgui SYSTEM PUBLIC
        "${imgui_SOURCE_DIR}" "${imgui_SOURCE_DIR}/backends")
    target_link_libraries(stereon_imgui PUBLIC glfw ${CMAKE_DL_LIBS})
endif()