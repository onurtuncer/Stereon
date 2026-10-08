# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(core
    MODULES
        Core.cppm
        Version.cppm)

target_compile_definitions(stereon_core PRIVATE
    STEREON_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    STEREON_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    STEREON_VERSION_PATCH=${PROJECT_VERSION_PATCH})

# Handle layout must agree between core and every consumer, so this is PUBLIC.
if(STEREON_CHECKED_HANDLES)
    target_compile_definitions(stereon_core PUBLIC STEREON_CHECKED_HANDLES=1)
endif()
