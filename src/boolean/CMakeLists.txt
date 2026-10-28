# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(boolean
    MODULES
        Assemble.cppm
        Bodies.cppm
        Boolean.cppm
        Capabilities.cppm
        Classify.cppm
        Clip.cppm
        Imprint.cppm
        Winding.cppm
    SOURCES
        Assemble.cpp
        Bodies.cpp
        Classify.cpp
        Clip.cpp
        Imprint.cpp
        Winding.cpp
    DEPENDS
        build
        intersect)
