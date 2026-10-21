.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Modeling operations
===================

Primitives, extrude, revolve, sweep, loft, sewing, Boolean operations, blends and operation history.

Primitive solids
----------------

``MakeBoxSolid`` and ``MakeCylinderSolid`` (``Stereon.Build``) make one solid
of one closed shell in any orthonormal frame: an origin and two axes, *x*
and *y*, with *z* their cross product, as for the analytic surfaces.

.. code-block:: cpp

   import Stereon.Build;

   const auto box = Stereon::MakeBoxSolid(
       Stereon::Box3{{0, 2}, {0, 3}, {0, 5}}, Stereon::Context::Default());
   const auto cylinder = Stereon::MakeCylinderSolid(
       Stereon::Point3{1, 1, 0}, {1, 0, 0}, {0, 1, 0}, /*radius*/ 0.5, /*height*/ 6,
       Stereon::Context::Default());

- A **box** has 8 vertices, 12 line edges and 6 planar faces. ``MakeBoxSolid``
  also takes an origin, two axes and the three side lengths.
- A **cylinder** has 2 vertices on the seam at angle 0, two circular edges, a
  seam, a lateral cylindrical face and two planar disks. The seam is a
  ``Boundary`` edge that the lateral face uses twice.

Every edge except the seam is ``Exact``, with exact pcurves on its faces
(ADR-0009). Faces are oriented outward. Vertex and edge tolerances equal the
``Context``'s resolution, because native construction is resolution-tight
(ADR-0001).

The result is published by ``Commit``, so it is certified at the Structural
level and checked at the Geometric level, and its History lists every entity
as Generated. A frame whose axes are not orthonormal, or a dimension that is
not finite or not larger than the resolution, is ``InvalidInput``.

These are the first operations of the Boolean prototype
(``docs/planning/boolean-prototype-plan.md``, milestone M1).

.. todo::

   Extrude, revolve, sweep, loft, sewing, Boolean operations, blends and
   operation history.

Related design decisions
------------------------

- `ADR-0011: Operation history and topological naming <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0011-history-and-naming.md>`_

