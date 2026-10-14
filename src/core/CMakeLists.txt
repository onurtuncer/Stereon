# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(core
    MODULES
        Arena.cppm
        Box.cppm
        Context.cppm
        Core.cppm
        Enum.cppm
        EvaluatorChecks.cppm
        Error.cppm
        Handle.cppm
        Interval.cppm
        Journal.cppm
        Log.cppm
        Matrix.cppm
        Scalar.cppm
        Serialize.cppm
        Simd.cppm
        SmallBuffer.cppm
        Trigonometry.cppm
        Vector.cppm
        Version.cppm
    SOURCES
        Arena.cpp
        EvaluatorChecks.cpp
        Handle.cpp
        Serialize.cpp
        Trigonometry.cpp
        Log.cpp)

# stdexec is PUBLIC: <Stereon/Core/Exec.hpp> is a public header that includes
# it. xsimd is only included in the global module fragment of the Simd
# partition, so consumers never see its headers.
target_link_libraries(stereon_core PUBLIC stereon::stdexec PRIVATE xsimd)

target_compile_definitions(stereon_core PRIVATE
    STEREON_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    STEREON_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    STEREON_VERSION_PATCH=${PROJECT_VERSION_PATCH})

# These change declarations seen by every consumer (handle layout, which
# contracts exist), so they are PUBLIC and must agree across the build.
target_compile_definitions(stereon_core PUBLIC
    STEREON_CHECKED_HANDLES=$<BOOL:${STEREON_CHECKED_HANDLES}>
    STEREON_HOT_CONTRACTS=$<BOOL:${STEREON_HOT_CONTRACTS}>)
