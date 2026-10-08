# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.30)

project(Stereon
    VERSION 0.0.1
    DESCRIPTION "A modern C++26 B-rep geometric modeling kernel"
    HOMEPAGE_URL "https://github.com/onurtuncer/Stereon"
    LANGUAGES CXX)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16)
    message(WARNING
        "Stereon's reference compiler is GCC 16 or newer (ADR-0016); "
        "found ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. "
        "Other compilers are not supported yet.")
endif()

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(StereonOptions)
include(StereonDependencies)

add_subdirectory(src)

if(STEREON_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
