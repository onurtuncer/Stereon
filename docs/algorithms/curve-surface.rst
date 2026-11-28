.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Curve–surface intersections
===========================

For an analytic surface, substitution of the curve into an implicit equation reduces the dimension of the problem. For a spline surface, work instead with an enclosed curve piece and a surface patch. Disjoint hulls eliminate candidates; a root enclosure for C(t)−S(u,v) can certify an isolated meeting. Clip both parameter domains and distinguish tangencies and contained curves from regular crossings. The chapter gives the separate search and overlap policies for the supported families.

For ``LinePlane``, let n be the cross product of the plane axes, P its origin
and O + t d the line. A nonzero n·d gives t = n·(P−O)/(n·d).
Compare the rational parameter with the requested line interval and both
plane bounds before constructing a double-precision point. When n·d is
exactly zero, a nonzero numerator proves separation. If both vanish,
clip the projected line against the plane box and return the overlap.

Method references: :cite:p:`piegl1997nurbs,sederberg1990curve,krawczyk1969newton`.

Detailed treatment: :doc:`/intersection`.

Implementation notes
--------------------

.. algorithm-sources::
   src/intersect/BezierArc.hpp
   src/intersect/BoxClip.hpp
   src/intersect/ConicArc.hpp
   src/intersect/CurvePieces.hpp
   src/intersect/CurveSplineSurface.cpp
   src/intersect/LineSurface.cpp
   src/intersect/SplineSurface.cpp
