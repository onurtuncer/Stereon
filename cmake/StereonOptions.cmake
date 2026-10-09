# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Build options and the stereon_options interface target that every Stereon
# library and test links privately. Third-party code never sees these flags.

include_guard(GLOBAL)

option(STEREON_BUILD_TESTS "Build the unit and property tests" ON)
option(STEREON_BUILD_BENCH "Build the benchmark suite" OFF)
option(STEREON_BUILD_PYTHON "Build the Python bindings (nanobind)" OFF)
option(STEREON_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(STEREON_COVERAGE "Instrument Stereon libraries and tests for gcov coverage (GCC)" OFF)
option(STEREON_CHECKED_HANDLES
    "Handles carry store ID and generation; every access is checked (ADR-0004)" OFF)

option(STEREON_HOT_CONTRACTS
    "Compile hot-path contracts (STEREON_HOT_PRE/POST/ASSERT); release builds turn them off (ADR-0016)" ON)

set(STEREON_CONTRACT_SEMANTIC "enforce" CACHE STRING
    "Contract evaluation semantic (ADR-0016): ignore, observe, enforce or quick_enforce")
set_property(CACHE STEREON_CONTRACT_SEMANTIC
    PROPERTY STRINGS ignore observe enforce quick_enforce)

set(STEREON_SANITIZER "" CACHE STRING "Sanitizer to build with: empty, address or thread")
set_property(CACHE STEREON_SANITIZER PROPERTY STRINGS "" address thread)

set(_stereon_semantics ignore observe enforce quick_enforce)
if(NOT STEREON_CONTRACT_SEMANTIC IN_LIST _stereon_semantics)
    message(FATAL_ERROR "STEREON_CONTRACT_SEMANTIC must be one of: ${_stereon_semantics}")
endif()

add_library(stereon_options INTERFACE)

# GCC 16 is the reference compiler (ADR-0016). Clang and MSVC are built in
# non-blocking CI; they get C++26 language flags only when they accept them.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(stereon_options INTERFACE
        -freflection
        -fcontracts
        -fcontract-evaluation-semantic=${STEREON_CONTRACT_SEMANTIC})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-freflection STEREON_CLANG_HAS_REFLECTION)
    check_cxx_compiler_flag(-fcontracts STEREON_CLANG_HAS_CONTRACTS)
    target_compile_options(stereon_options INTERFACE
        $<$<BOOL:${STEREON_CLANG_HAS_REFLECTION}>:-freflection>
        $<$<BOOL:${STEREON_CLANG_HAS_CONTRACTS}>:-fcontracts>)
endif()

if(MSVC)
    # C4324 ("structure was padded due to alignment specifier") is
    # informational; stdexec triggers it in templates instantiated in Stereon
    # modules, where /external:W0 does not reach.
    target_compile_options(stereon_options INTERFACE
        /W4 /permissive- /utf-8 /Zc:__cplusplus /Zc:preprocessor /wd4324
        $<$<BOOL:${STEREON_WARNINGS_AS_ERRORS}>:/WX>)
else()
    target_compile_options(stereon_options INTERFACE
        -Wall -Wextra -Wpedantic
        -Wconversion -Wsign-conversion -Wshadow
        -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual
        -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough
        $<$<BOOL:${STEREON_WARNINGS_AS_ERRORS}>:-Werror>)
endif()

if(MSVC AND STEREON_SANITIZER STREQUAL "address")
    target_compile_options(stereon_options INTERFACE /fsanitize=address)
elseif(MSVC AND NOT STEREON_SANITIZER STREQUAL "")
    message(FATAL_ERROR "MSVC supports only STEREON_SANITIZER=address")
elseif(STEREON_SANITIZER STREQUAL "address")
    target_compile_options(stereon_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_link_options(stereon_options INTERFACE -fsanitize=address,undefined)
elseif(STEREON_SANITIZER STREQUAL "thread")
    target_compile_options(stereon_options INTERFACE -fsanitize=thread -fno-omit-frame-pointer)
    target_link_options(stereon_options INTERFACE -fsanitize=thread)
elseif(NOT STEREON_SANITIZER STREQUAL "")
    message(FATAL_ERROR "STEREON_SANITIZER must be empty, 'address' or 'thread'")
endif()

if(STEREON_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "STEREON_COVERAGE requires GCC (gcov)")
    endif()
    # Atomic counters keep coverage correct once operations run in parallel.
    target_compile_options(stereon_options INTERFACE --coverage -fprofile-update=atomic)
    target_link_options(stereon_options INTERFACE --coverage)
endif()

if(MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # std::print and std::println need libstdc++exp on MinGW
    # (docs/toolchain/gcc16-feature-probe.md).
    target_link_libraries(stereon_options INTERFACE stdc++exp)
endif()

# stereon_add_library(<name> MODULES <files...> [SOURCES <files...>] [DEPENDS <libs...>])
#
# Creates stereon_<name> with alias stereon::<name>.
#
# - MODULES are the module interface units (primary interface and partitions)
#   and go in a public CXX_MODULES file set.
# - SOURCES are module implementation units and other private sources.
# - DEPENDS names lower Stereon libraries, without the stereon:: prefix. They
#   are linked PUBLIC because a module interface may re-export them. The order
#   in src/CMakeLists.txt is the layering; tools/check_layering.py rejects any
#   upward dependency.
#
# Include visibility: src/<name>/include/ is the only public include
# directory (macro headers, since modules cannot export macros). Everything
# else under src/<name>/ is private to the library.
function(stereon_add_library name)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "MODULES;SOURCES;DEPENDS")
    if(arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "stereon_add_library(${name}): unknown arguments ${arg_UNPARSED_ARGUMENTS}")
    endif()
    if(NOT arg_MODULES)
        message(FATAL_ERROR "stereon_add_library(${name}): at least one module interface is required")
    endif()

    add_library(stereon_${name})
    add_library(stereon::${name} ALIAS stereon_${name})
    target_sources(stereon_${name}
        PUBLIC FILE_SET CXX_MODULES FILES ${arg_MODULES}
        PRIVATE ${arg_SOURCES})

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
        target_include_directories(stereon_${name}
            PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
    endif()
    target_include_directories(stereon_${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

    foreach(dep IN LISTS arg_DEPENDS)
        target_link_libraries(stereon_${name} PUBLIC stereon::${dep})
    endforeach()
    target_link_libraries(stereon_${name} PRIVATE stereon_options)
    target_compile_features(stereon_${name} PUBLIC cxx_std_${STEREON_CXX_STANDARD})
    set_target_properties(stereon_${name} PROPERTIES EXPORT_NAME ${name})
endfunction()

# stereon_exec_sources(<files...>)
#
# Marks source files that include <Stereon/Core/Exec.hpp>. They are excluded
# from module scanning, so CMake does not compile them with -fmodules: GCC 16.2
# crashes on stdexec under that flag. Such files cannot import modules. Call it
# in the directory that adds the files to their target.
#
# GCC's -Wnull-dereference is also turned off for them: it fires in stdexec's
# inlined intrusive queue after optimisation, where -isystem does not reach.
function(stereon_exec_sources)
    set_source_files_properties(${ARGN} PROPERTIES CXX_SCAN_FOR_MODULES OFF)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set_property(SOURCE ${ARGN} APPEND PROPERTY COMPILE_OPTIONS -Wno-null-dereference)
    endif()
endfunction()
