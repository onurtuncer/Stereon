.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Topology — ``Stereon.Topo``
===========================

The boundary representation itself: solids, shells, faces, loops, coedges,
edges and vertices, and the geometry they carry. A B-rep lives in one of
two states:

- A **ShapeBuilder** is mutable and owned by one thread. Every change of
  connectivity goes through an operation that keeps the structure
  consistent, among them the classical Euler operators.
- A **Shape** is frozen. It is immutable, cheap to copy and safe to read
  from any thread. Only ``Freeze`` (or ``Commit``, for an operation's
  result) makes one, and only after the structural validator has passed
  it, so every ``Shape`` you hold is certified.

Editing a shape opens a builder on it. The builder shares every page of the
frozen tables until it writes to one, so an edit of one face copies almost
nothing (ADR-0004, ADR-0005).

.. code-block:: cpp

   import Stereon.Topo;

   using namespace Stereon;

   ShapeBuilder b = *ShapeBuilder::Open(shape);
   *b.RemoveFace(b.Handle(face));
   Result<Shape> edited = Freeze(std::move(b));   // shares every untouched page

   for (const FaceId f : *AllEntities<EntityKind::Face>(*edited)) {
       const Surface s = SurfaceOf(*edited, f);
   }

Start here
----------

``ShapeBuilder``
   Building and editing. Add geometry to its pools (``AddSurface``,
   ``AddCurve``, ``AddPCurve``), then topology. The low-level path is
   ``AddFace``, ``AddLoop``, ``AddEdge`` and ``AppendCoedge``; the Euler
   operators are ``MakeVertexFaceShell``, ``MakeEdgeVertex``,
   ``MakeEdgeFace`` and their inverses. ``Transaction`` runs a compound
   edit that either completes or leaves the builder as it was.

``Freeze``, ``Commit``
   Validate and publish. ``Commit`` also returns the ``History`` of what
   became of the input's entities, as an ``OperationResult``.

``Shape`` and its queries
   ``AllEntities``, ``FacesOf``, ``LoopsOf``, ``CoedgesOf``,
   ``RadialCoedgesOf``, ``EdgesOf``, ``SurfaceOf``, ``CurveOf``,
   ``PCurveOf``, ``ToleranceOf``, ``Describe``. Ranges come in a
   documented, deterministic order.

Handles
   ``VertexId`` ... ``SolidId`` address a frozen shape; ``BuilderVertexId``
   ... ``BuilderSolidId`` address a builder and carry a generation, so a
   stale one is caught. Their types are in :doc:`core`.

``ValidityReport``
   What the validator found: rule IDs and the entities involved, not a
   bare yes or no.

The guide chapter behind this library: :doc:`/topology`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: topo

How the entities fit together
-----------------------------

One table per entity kind, each record linking to others by handle. A
solid has shells and a shell has faces. A face is bounded by loops, and a
loop is a ring of coedges, each one use of an edge by the loop. An edge runs
from one vertex to another, and all the coedges that use it form its radial
ring. The ``Next...`` links thread the lists and rings that keep each
collection in order.

.. stereon-type-map:: topo
   :caption: The entity records of one store and the handles that link
             them, read from ``Records.cppm``: an arrow is a field holding
             a handle, labelled with the field's name.

Inside the module
-----------------

.. stereon-partition-graph:: topo

Reference
---------

.. stereon-api:: topo
