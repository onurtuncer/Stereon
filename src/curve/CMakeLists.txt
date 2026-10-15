# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

stereon_add_library(curve
    MODULES
        BSplineBasis.cppm
        BSplineCurve.cppm
        Circle.cppm
        Concepts.cppm
        Conversion.cppm
        Curve.cppm
        Ellipse.cppm
        Erased.cppm
        Fitting.cppm
        Hyperbola.cppm
        Line.cppm
        Properties.cppm
        Visit.cppm
    SOURCES
        BSplineBasis.cpp
        BSplineCurve.cpp
        Conics.cpp
        Conversion.cpp
        Factories.cpp
        Fitting.cpp
        Properties.cpp
    DEPENDS
        robust)
