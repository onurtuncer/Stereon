.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Sampling and display tessellation
=================================

Start curve sampling at knot spans and subdivide according to geometric departure from the chord. Initial subdivisions reduce the chance of missing an inflection whose midpoint lies on the chord. Face sampling adds vertices and triangles in the surface chart; the manual identifies the supported trim cases. The display tessellation's stopping criteria are separate from the interval certificate used when publishing an intersection.

Method references: :cite:p:`piegl1997nurbs,patrikalakis2002shape`.

Detailed treatment: :doc:`/meshing`.

Implementation notes
--------------------

.. algorithm-sources::
   src/mesh/Mesh.cppm
   src/mesh/Sampling.cpp
   src/mesh/Sampling.cppm
