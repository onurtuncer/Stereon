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

A face receives a *curve network* from the intersection: a ``FaceNetwork``
of vertices, each a crossing of the face's boundary or a junction inside
the face where the intersection passes from one face of the other body to
the next, and edges, the pieces of branches between them.
``ImprintNetworks(builder, shape, networks, ctx)`` writes a body's networks
into a ``ShapeBuilder`` opened on the shape, inside the caller's
transaction, through the Euler operators of :doc:`topology`. The edges the
networks cross are split at the crossings (SEMV), once per distinct
crossing, so that the two faces sharing an edge meet at the same vertex.
Each piece then becomes an edge: a strut when it reaches a junction not yet
made (MEV); a split of the face when both its ends exist on one loop (MEF),
or a bridge when on two loops of the face (MEKR); a component that touches
no boundary is hung on a helper strut from a loop vertex, removed once the
component has closed (KEMR), which leaves a ring as a loop of the face and
the ring's other use as the new face's boundary. The ring must run against
the face's outer loop, a hole; where the pcurve turns the other way the
ring and the disc change faces (KFMRH, then MFKRH), so the disc is the
piece's face. Where a piece starts and ends on one edge, a cylinder's
seam, both uses of that edge end at the vertex, one from either side of
the chart, and the imprint takes the use the piece's own pcurve reaches,
so the loops stay continuous in the face's parameter space. The new edge
is ``Exact`` with the
branch's curve where that is a line, a circle or an ellipse, and otherwise
a ``Trim`` edge with its cached curve (ADR-0009, section 8), its tolerance
the larger of the resolution and the branch's certified bound. The
``Imprint`` names each piece's edge, its end vertices and the face it cut
off, the builder's vertex for every network vertex, and the vertices the
splits made. ``NetworkOf`` turns one clip into a network, and
``ImprintClip`` imprints it.

``ImprintBodies(a, b, ctx)`` is the pipeline's first four steps for two
solids: candidate face pairs by bounding box in input order, the surface
intersection of each pair over the faces' parameter boxes, each branch
clipped to both faces, and both bodies imprinted with the networks their
faces receive, inside one transaction per body. The common stretch of a
branch inside both faces is a piece of the intersection; its end is a
crossing in the face whose boundary it meets and a junction in the other.
The result holds both ``OperationResult`` s and the pairs of edges, one in
each body, that are the same piece, for the assembly of M8 to identify;
coincident faces (an overlap) are ``NotSupported`` until then.

.. code-block:: cpp

   const auto bodies = Stereon::ImprintBodies(boxA, boxB, ctx);
   for (const auto& [edgeInA, edgeInB] : bodies->Shared) { /* one piece, twice */ }

   // One face, one branch, by hand.
   const auto clip = Stereon::ClipToFace(shape, face, branch, Stereon::BranchSide::A, ctx);
   auto builder = *Stereon::ShapeBuilder::Open(shape);
   const auto imprint = builder.Transaction([&](Stereon::ShapeBuilder& candidate) {
       return Stereon::ImprintClip(candidate, shape, face, branch, Stereon::BranchSide::A, *clip, ctx);
   });
   const auto result = Stereon::Commit(std::move(builder), ctx);

Classification
--------------

``Classify(bodies, ctx)`` is pipeline step 5 (M7; ADR-0013). The faces of
each imprinted body fall into *regions*, connected across every edge but
the imprinted ones (a union-find over the faces), and each region lies
``Inside`` the other solid or ``Outside`` it. One seed per region is
classified by deterministic rays: a point inside one of the region's faces
(probed at the parameter box's centre and then at offsets from the
midpoints of its edges' pcurves, inward and outward, starting at a quarter
of the face's parameter extent and halved forty times, every edge at one
offset before the next smaller, so the widest seed is found first and a
region however narrow still gets one; each probe settled by
``FaceContains``), then a ray
from it cast against every face of the other body, which it crosses an odd
number of times from inside. A plane is met once and a cylinder up to
twice, in closed form; a hit counts when ``FaceContains`` puts it inside
the face. A ray that runs along a surface, grazes a cylinder, starts on a
face or hits a face's boundary is ambiguous, and the next of sixteen
directions (irrational fractions of a turn, tilted) is tried; a seed whose
every direction is ambiguous falls back to its winding number, below.

When every direction is ambiguous, the seed's *bounded generalised winding
number* decides (``WindingNumberOf``, ADR-0013's decisions on acceptance).
Each face of the other body is triangulated in its parameter space, every
face sharing an edge taking the same points along it, so the triangles
close up into a polyhedron; the polyhedron's winding number at the seed is
the sum of the triangles' solid angles (Van Oosterom and Strackee) over
:math:`4\pi`. It equals the body's own winding number wherever the seed is
farther from every triangle than the face can stray from them: a surface
that agrees with the face inside and runs along the chords at its boundary
deforms onto the triangles within that distance without crossing the seed,
so the two numbers differ by a whole number, which must be zero. That
distance is summed from the cylinder's sagitta over a cell's turn, the
chords' deviation from the edges (:math:`h^2 / 8` of a bound on the second
derivative of :math:`S(p(t))`, :math:`r\,u'^2 + (\lvert S_u \rvert +
\lvert S_v \rvert)\,\lvert p'' \rvert`, with :math:`u'` and :math:`p''`
bounded in closed form for lines, circles and ellipses and by the poles of
the derivative curves for non-rational B-splines) and the tolerances of
the face's edges and vertices. An edge is sampled at every knot of a
B-spline pcurve of any of its uses, since a B-spline is smooth only within
a knot span and a polyline turns its corners at the knots, and uniformly
within each span. Its samples are capped at :math:`2^{20}`; where the cap
binds, the deviation the samples reach, not the target, enters the bound. A loop whose pcurves do not meet in the
face's parameter space, a face split across a seam say, is
``InvalidResult`` naming the two edges rather than a sum. A cell of a face's parameter space that
is not separated from the seed is split across its longer side and its
loops clipped to the halves (Sutherland–Hodgman), down to the round's
deviation target; a round that leaves a cell unseparated is followed by
one with the target quartered, from a sixty-fourth of the body's extent
down to the resolution, up to ``ClassifyBudget::Refinements`` rounds (eight
by default). The result carries its bound, the cells left unseparated plus
the rounding of the sum, and a side is read from it only when the number
lies further from one half than the bound (ADR-0029, section 4); a decided
number that is neither 0 nor 1 is ``InvalidResult``, since the faces then
do not close up or one is oriented against its loops. A seed undecided
after the budget leaves its region ``Unresolved`` and the operation fails
with ``ErrorCode::Unresolved`` naming the face, the rays tried and the
last interval (decision D6 of the plan), never guessing. The rays and the
rounds are a ``ClassifyBudget``, so a caller, or a test, can ask for no
rays and go straight to the winding number.

Every region is seeded on its own, and the seeds are then checked against
each other across the imprinted edges by the branches' side transitions
(ADR-0029, section 3): the regions on either side of an ``OppositeSide``
edge must differ, of a ``SameSide`` edge agree, and an ``Unresolved`` edge
checks nothing. A disagreement fails the operation with
``DegenerateResult``, since the intersection, the imprint or a seed went
wrong: this is the cross-check table of ADR-0013.

Assembly
--------

``Boolean(a, b, operation, ctx)`` is pipeline step 6 (M8; ADR-0013,
ADR-0030): the regularised ``Union``, ``Cut`` or ``Common`` of two
solids. The bodies are imprinted and their regions classified, and the
regions the operation keeps are copied into a new shape, A's before B's:
for a union each body's regions outside the other, for a common those
inside it, and for a cut A's regions outside B with B's regions inside A
turned inside out (the face reversed, its loops traversed backwards). A
piece of the intersection, an edge in each body carrying the same curve
over the same range, becomes one edge used by both bodies' faces, A's
where a kept face of A uses it; its end vertices are identified within the
pair tolerance, the larger of the two, A's position kept and any growth of
the tolerance logged (ADR-0001's representative policy). What no kept face
uses is never copied, so the result is regularised: a contact of lower
dimension leaves nothing behind, and two solids apart have nothing in
common, a successful empty result. The faces fall into shells by the
edges they share; a shell made only of B's reversed faces is a void of the
solid around it, and every other shell bounds a solid of its own.

The result is published through ``Commit`` with a History of both
operands, input 0 and input 1. Each entity of the result is traced through
the imprint, its History and its split events, to the operand entities it
continues: an operand entity continued by several is a split, several
continued by one a merge, each grouped; an intersection edge is
``Generated`` by the face of each operand that met there; and every
operand entity nothing continues is ``Deleted``. ``Commit`` checks the
History is complete before anything is published.

Coincident faces, ``OnSame`` and ``OnOpposite``, are still refused by the
imprint with ``NotSupported``; they are M8's next step.

Related design decisions
------------------------

- `ADR-0011: Operation history and topological naming <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0011-history-and-naming.md>`_

