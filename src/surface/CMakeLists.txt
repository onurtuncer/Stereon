# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(surface
    MODULES
        BSplineSurface.cppm
        Capabilities.cppm
        Certify.cppm
        Concepts.cppm
        Cone.cppm
        Conversion.cppm
        Cylinder.cppm
        Erased.cppm
        Fitting.cppm
        Plane.cppm
        Properties.cppm
        Sphere.cppm
        Swept.cppm
        Surface.cppm
        Torus.cppm
        Trimmed.cppm
        Visit.cppm
    SOURCES
        BSplineSurface.cpp
        Certify.cpp
        Conversion.cpp
        Factories.cpp
        Fitting.cpp
        Surface.cpp
        Swept.cpp
        Trimmed.cpp
    DEPENDS
        curve)
