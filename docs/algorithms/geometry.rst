.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Analytic geometry, differential properties and projection
=========================================================

The native curves and surfaces evaluate in local frames. Derivatives give tangents, normals and curvature where the parameterisation is regular. Arc length integrates speed by adaptive Gauss–Legendre quadrature; inverse arc length refines the parameter using the accumulated integral. Projection must consider admissible extrema and parameter-domain boundaries. Trimming restricts a parent parameterisation, while a sweep introduces an additional parameter. The associated chapters give singularities, seam conventions and supported cases. See :doc:`/properties`, :doc:`/trimmed`, :doc:`/swept` and :doc:`/interchange`.

Method references: :cite:p:`docarmo1976differential,patrikalakis2002shape,davis1984numerical`.

Detailed treatment: :doc:`/analytic`.

Implementation notes
--------------------

.. algorithm-sources::
   src/curve/Capabilities.cppm
   src/curve/Circle.cppm
   src/curve/Concepts.cppm
   src/curve/Conics.cpp
   src/curve/Conversion.cpp
   src/curve/Conversion.cppm
   src/curve/Curve.cppm
   src/curve/Ellipse.cppm
   src/curve/Erased.cppm
   src/curve/Factories.cpp
   src/curve/Hyperbola.cppm
   src/curve/Interchange.cpp
   src/curve/Interchange.cppm
   src/curve/Line.cppm
   src/curve/Properties.cppm
   src/curve/Trimmed.cpp
   src/curve/Trimmed.cppm
   src/curve/Visit.cppm
   src/surface/Capabilities.cppm
   src/surface/Concepts.cppm
   src/surface/Cone.cppm
   src/surface/Conversion.cpp
   src/surface/Conversion.cppm
   src/surface/Cylinder.cppm
   src/surface/Ellipsoid.cpp
   src/surface/Ellipsoid.cppm
   src/surface/Erased.cppm
   src/surface/Factories.cpp
   src/surface/Plane.cppm
   src/surface/Properties.cppm
   src/surface/Sphere.cppm
   src/surface/Surface.cpp
   src/surface/Surface.cppm
   src/surface/Swept.cpp
   src/surface/Swept.cppm
   src/surface/Torus.cppm
   src/surface/Trimmed.cpp
   src/surface/Trimmed.cppm
   src/surface/Visit.cppm
