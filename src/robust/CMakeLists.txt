# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(robust
    MODULES
        Expansion.cppm
        Filter.cppm
        Predicates.cppm
        Robust.cppm
        Sign.cppm
    DEPENDS
        core)

# The error-free transformations rely on every product and sum being rounded
# exactly once. A fused multiply-add skips a rounding and breaks them. The
# option is PUBLIC because the arithmetic is inline and its templates are
# instantiated in importing units, which must follow the same rule. GCC
# already turns contraction off for -std=c++26 (no GNU extensions); Clang
# contracts within expressions by default. MSVC contracts only under
# /fp:contract or /fp:fast, neither of which the build uses.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(stereon_robust PUBLIC -ffp-contract=off)
endif()
