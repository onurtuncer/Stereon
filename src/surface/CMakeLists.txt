# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(surface
    MODULES
        BSplineSurface.cppm
        Concepts.cppm
        Cone.cppm
        Cylinder.cppm
        Erased.cppm
        Plane.cppm
        Sphere.cppm
        Surface.cppm
        Torus.cppm
        Visit.cppm
    SOURCES
        BSplineSurface.cpp
        Factories.cpp
        Surface.cpp
    DEPENDS
        curve)
