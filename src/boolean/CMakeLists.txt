# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(boolean
    MODULES
        Bodies.cppm
        Boolean.cppm
        Clip.cppm
        Imprint.cppm
    SOURCES
        Bodies.cpp
        Clip.cpp
        Imprint.cpp
    DEPENDS
        build
        intersect)
