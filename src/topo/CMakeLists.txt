# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(topo
    MODULES
        Builder.cppm
        Capabilities.cppm
        Operation.cppm
        Records.cppm
        Report.cppm
        Shape.cppm
        Topo.cppm
    SOURCES
        Builder.cpp
        Compact.cpp
        Describe.cpp
        Operation.cpp
        Shape.cpp
        Validate.cpp
    DEPENDS
        intersect)
