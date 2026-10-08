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
option(STEREON_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(STEREON_CHECKED_HANDLES
    "Handles carry store ID and generation; every access is checked (ADR-0004)" OFF)

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

target_compile_options(stereon_options INTERFACE
    -freflection
    -fcontracts
    -fcontract-evaluation-semantic=${STEREON_CONTRACT_SEMANTIC}
    -Wall -Wextra -Wpedantic
    -Wconversion -Wsign-conversion -Wshadow
    -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual
    -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough
    $<$<BOOL:${STEREON_WARNINGS_AS_ERRORS}>:-Werror>)

if(STEREON_SANITIZER STREQUAL "address")
    target_compile_options(stereon_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_link_options(stereon_options INTERFACE -fsanitize=address,undefined)
elseif(STEREON_SANITIZER STREQUAL "thread")
    target_compile_options(stereon_options INTERFACE -fsanitize=thread -fno-omit-frame-pointer)
    target_link_options(stereon_options INTERFACE -fsanitize=thread)
elseif(NOT STEREON_SANITIZER STREQUAL "")
    message(FATAL_ERROR "STEREON_SANITIZER must be empty, 'address' or 'thread'")
endif()

if(MINGW)
    # std::print and std::println need libstdc++exp on MinGW
    # (docs/toolchain/gcc16-feature-probe.md).
    target_link_libraries(stereon_options INTERFACE stdc++exp)
endif()

# stereon_add_library(<name> MODULES <files...>)
#
# Creates stereon_<name> with alias stereon::<name>. Module interface units
# go in a CXX_MODULES file set; private dependencies on other Stereon
# libraries are added by the caller with target_link_libraries.
function(stereon_add_library name)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "MODULES;SOURCES")
    add_library(stereon_${name})
    add_library(stereon::${name} ALIAS stereon_${name})
    target_sources(stereon_${name}
        PUBLIC FILE_SET CXX_MODULES FILES ${arg_MODULES}
        PRIVATE ${arg_SOURCES})
    target_link_libraries(stereon_${name} PRIVATE stereon_options)
    target_compile_features(stereon_${name} PUBLIC cxx_std_26)
    set_target_properties(stereon_${name} PROPERTIES EXPORT_NAME ${name})
endfunction()
