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

Face clipping and imprinting
----------------------------

The Boolean prototype (``docs/planning/boolean-prototype-plan.md``) builds
its pipeline in ``Stereon.Boolean`` milestone by milestone; M6 is the two
steps that take an intersection branch into a body's topology.

``ClipToFace(shape, face, branch, side, ctx)`` restricts a branch of two
surfaces to a face on one of them, in that face's parameter space. The
branch's pcurve on the face (``side`` says which of its two) is met with
the pcurve of every coedge of the face's loops through the plane curve
intersection of :doc:`intersection`; the meetings are the ``Crossings`` of
the ``FaceClip``, in order along the branch, each with the coedge, the
parameter on its edge and on the branch, and the loop vertex it is at when
it lies within that vertex's tolerance (ADR-0001). The stretches between
crossings are ``Inside`` the face or outside it by the parity of a ray from
their middle, the ray turned by irrational fractions of a turn until it
meets no loop vertex; stretches along a loop edge are ``OnBoundary``.
Everything else is decided on the nominal geometry, exactly where the
pcurves are lines, circles and ellipses (decision D1 of the plan), and
``Decided`` says whether every crossing and every side was.

``ImprintClip(builder, shape, face, branch, side, clip, ctx)`` writes the
clip into a ``ShapeBuilder`` opened on the shape, inside the caller's
transaction, through the Euler operators of :doc:`topology`: the edges the
branch crosses are split at the crossings (SEMV), and each inside piece
becomes an edge, splitting the face in two when its ends lie on one loop
(MEF), bridging two loops of the face when they lie on different ones
(MEKR), or a ring inside the face when the piece is closed (a strut by
MEV, the closed edge by MEF on the strut's own use, the strut removed by
KEMR, which leaves the ring as a loop of the face and the closed edge's
other use as the new face's boundary). The new edge is ``Exact`` with the
branch's curve where that is a line, a circle or an ellipse, and otherwise
a ``Trim`` edge with its cached curve (ADR-0009, section 8), its tolerance
the larger of the resolution and the branch's certified bound. The
``Imprint`` names each piece's edge, its end vertices and the face it cut
off, and the vertices the splits made. A piece that ends inside the face
is ``NotSupported`` here: branches meet at nodes, which the assembly of M8
joins.

.. code-block:: cpp

   const auto clip = Stereon::ClipToFace(shape, face, branch, Stereon::BranchSide::A, ctx);
   auto builder = *Stereon::ShapeBuilder::Open(shape);
   const auto imprint = builder.Transaction([&](Stereon::ShapeBuilder& candidate) {
       return Stereon::ImprintClip(candidate, shape, face, branch, Stereon::BranchSide::A, *clip, ctx);
   });
   const auto result = Stereon::Commit(std::move(builder), ctx);

Related design decisions
------------------------

- `ADR-0011: Operation history and topological naming <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0011-history-and-naming.md>`_

