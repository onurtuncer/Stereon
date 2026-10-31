# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(intersect
    MODULES
        Capabilities.cppm
        Dispatch.cppm
        Intersect.cppm
    SOURCES
        CylinderCylinder.cpp
        CylinderCylinderSkew.cpp
        Dispatch.cpp
        LineLine.cpp
        LineSurface.cpp
        PlanarConics.cpp
        PlanarLines.cpp
        PlanarSpline.cpp
        PlaneCylinder.cpp
        PlanePlane.cpp
    DEPENDS
        surface)

