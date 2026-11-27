.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Intersections in the plane
==========================

Lines and conics admit specialised constructions. For two conics a rational chart of one curve turns the other's implicit equation into a univariate polynomial. A candidate root must map back to both requested domains. For spline pairs, enclose Bézier pieces, discard separated boxes and use a Krawczyk test to isolate a regular crossing. Failure of an inclusion test leaves a box undecided. Coincidence requires an overlap representation; a repeated polynomial root alone does not establish a transverse crossing.

Method references: :cite:p:`sederberg1990curve,farouki1987numerical,krawczyk1969newton`.

Detailed treatment: :doc:`/intersection`.

Implementation notes
--------------------

.. algorithm-sources::
   src/intersect/Planar.hpp
   src/intersect/PlanarConics.cpp
   src/intersect/PlanarLines.cpp
   src/intersect/PlanarSpline.cpp
   src/intersect/PlanarSplinePair.cpp
