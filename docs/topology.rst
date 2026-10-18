.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

B-rep topology
==============

The topology store, typed handles, ShapeBuilder and immutable Shape, Euler operators, and the validity checker.

Handles
-------

Entities are addressed by typed 32-bit handles (ADR-0004). ``VertexId``,
``EdgeId``, ``CoedgeId``, ``LoopId``, ``FaceId``, ``ShellId`` and ``SolidId``
are the instances of ``Id<EntityKind>``, one type per kind, so passing an
edge handle where a vertex is expected, or a raw integer anywhere, is a
compile error. A default-constructed handle is invalid; ``IsValid()`` says so.
Handles compare and hash by index, so they key ``std::map`` and
``std::unordered_map`` side tables, and they format as ``Vertex#12``.

A ``ShapeBuilder`` hands out ``BuilderId<Kind>``, an index with a
generation. ``GenerationTable<Kind>`` owns the generations, the alive bitset
and the free list: releasing a slot bumps its generation, so every handle to
the released entity is stale and ``IsAlive()`` says no. Slots released in a
session are not reused in that session, which keeps an index that is live
before and after an edit on the same entity; ``BeginSession()`` makes them
reusable. Builder handles format as ``Vertex#12@3``.

``Oriented<Kind>`` carries a sense with the handle, as a loop uses a coedge
forwards or reversed. The sense lives in the top bit of the index and is
reached only through ``BaseOf(h)``, ``IsReversed(h)`` and ``Reversed(h)``;
nothing masks it by hand.

With ``STEREON_CHECKED_HANDLES`` (the debug and sanitizer presets) every
handle also records its store and generation, and ``GenerationTable::Check``
verifies store, bounds and generation. A mismatch prints a report naming the
handle, the entity kind and the expected store, then raises a contract
violation. Release builds compile ``Id`` down to its 4-byte index, which a
``static_assert`` guards. Handles are local to one store; identity across
operations comes from history maps (ADR-0011).

Storage pages
-------------

Entity tables and geometry pools keep their entries in ``PagedArena<T, N>``
(ADR-0005): fixed-size pages of ``N`` entries, 256 unless a table chooses
otherwise, behind ``std::shared_ptr``. Copying an arena copies the page
pointers, so a snapshot of a million-entry table costs a few thousand pointer
copies and shares every page with its source. The first write to a page
copies the page when another arena still holds it and writes in place when
the arena holds it alone; memory grows with what an edit touched, not with
the size of the table. ``PushBack`` appends, ``Resize`` grows with
value-initialised entries or drops whole pages, and ``GetPageSpan(page)``
hands hot loops one page at a time, so the page lookup and the reference
counting stay out of inner loops.

Writes are scoped. ``Set(index, value)`` replaces an entry, and
``Modify(index, fn)`` and ``ModifyPage(page, fn)`` call ``fn`` with a
reference to the entry or a span of the page, made unique first and valid
only for the duration of the call; whatever ``fn`` returns, by value, is
returned. No function returns a mutable reference. A reference kept past the
copy-on-write step would otherwise write into a page that a later copy
shares, and change a snapshot that is supposed to be immutable. While ``fn``
runs, the arena may be read but not copied, moved, assigned, cleared,
resized or edited again. Every build checks this with a plain branch and
reports a violation from a module implementation unit, so the check does
not depend on the compiler honouring contracts in imported templates
(``docs/toolchain/gcc16-feature-probe.md``).

.. code-block:: cpp

   PagedArena<FaceRecord> faces = frozen;      // shares every page
   faces.Set(17, record);                      // copies one page
   faces.Modify(42, [](FaceRecord& f) { f.Tolerance *= 2.0; });
   faces.ModifyPage(0, [](std::span<FaceRecord> page) {   // one copy-on-write step
       for (FaceRecord& f : page) { f.Tolerance = 1e-5; }
   });

A ``const`` arena can be read, and copied, from any thread; the mutating
operations belong to one thread, the builder's (ADR-0007).
``SharedPageCount(a, b)`` and ``IsPageShared(page)`` report what a snapshot
and its edit still share, which the tests and ``bench/ArenaBench.cpp`` use
to measure the copy cost. The page size is a tuning constant per table, fixed
by the Phase 0 measurement of ADR-0005 over 128, 256 and 1 024 entries;
nothing in the public API of the stores depends on it.

The store
---------

``Stereon.Topo`` keeps one record table per entity kind, each a
``PagedArena``: vertices (position, tolerance), edges, coedges, loops,
faces, shells and solids, following Golovanov's logical model (ADR-0004).
Records reach one another by index, and each record carries its own slot
state (alive, generation), so that a deletion copies one page rather than a
whole-table bitset (ADR-0005).

- An **edge** records its curve kind, its tolerance, its start and end
  vertices, the shared parameter range of its curve and pcurves, the anchor
  of its radial ring and, for the kinds that evaluate through a surface,
  its primary coedge (ADR-0009). ``Exact`` and ``Spatial`` edges carry
  their own 3D curve; ``Intersection``, ``Trim``, ``Boundary`` and
  ``Degenerate`` edges have none and are evaluated on the primary coedge's
  pcurve and surface.
- A **coedge** is one use of an edge in a loop: ``Next`` and ``Prev`` around
  the loop, ``Radial`` around the ring of uses of the same edge, its stored
  edge-use sense (``EdgeUseSense::Same`` along the edge's parameter or
  ``Opposite``) and its pcurve on the face's surface.
- A **loop** has its face and an anchor coedge; a **face** its surface, an
  orientation flag, its shell and its list of loops; a **shell** its solid
  and its faces; a **solid** its shells.

Geometry lives in three pools, of ``Curve``, ``Curve2d`` and ``Surface``
values, referenced by ``CurveRef``, ``PCurveRef`` and ``SurfaceRef``.
Several entities may share one value with bounds of their own: sharing a
support is not sharing topology (ADR-0025).

Building a shape
----------------

``ShapeBuilder`` is the mutable, single-threaded side, and the only writer
of the tables (ADR-0027). Its operations each keep every reciprocal link
consistent, so no caller writes a link field:

.. code-block:: cpp

   ShapeBuilder b;
   const SurfaceRef plane = b.AddSurface(Surface(Plane{...}));
   const BuilderFaceId face = *b.AddFace(plane);
   const BuilderLoopId loop = *b.AddLoop(face);
   const BuilderEdgeId edge = *b.AddEdge({.Kind = EdgeCurveKind::Exact,
                                          .Start = v0, .End = v1,
                                          .Curve = b.AddCurve(Curve(line)),
                                          .Range = Interval(0, 1)});
   *b.AppendCoedge(loop, edge, EdgeUseSense::Same, b.AddPCurve(pcurve));
   Result<Shape> shape = Freeze(std::move(b));

``AppendCoedge`` appends to the loop's cycle and inserts the coedge into its
edge's radial ring; the first use of a surface-defined edge becomes its
primary coedge, which only ``SetPrimaryCoedge`` changes. ``AddLoop``,
``AddShell`` and ``AddFaceToShell`` append to their parent's list, so
members come back in the order they were added. Every operation checks its
handles and its request and reports ``InvalidInput`` without changing
anything: a handle that is not a live entity of this builder, an empty
parameter range, a missing pcurve, an edge kind given a curve it does not
have, a face added to a second shell. Every change advances the builder's
topology epoch; its ``BuilderToken`` is stable across moves, and a
moved-from builder refuses every operation (ADR-0038).

Freeze and certification
------------------------

``Freeze(std::move(builder))`` is the certification gate of ADR-0038. It
runs the structural validator over every table and, only when it finds
nothing, moves the tables into an immutable ``TopologySnapshot`` with a
``ValidationCertificate`` (validator version, rule set, snapshot ID), and
returns a ``Shape`` that owns it. Indices are kept: a ``BuilderId`` names the
same slot as the frozen ``Id``, and ``shape.Resolve(builderId)`` gives it.
Otherwise ``Freeze`` reports ``InvalidInput`` with the findings, one line
each. The builder is consumed in both cases. ``builder.Validate()`` runs the
same check without freezing, and ``ValidateTopology(shape)`` runs it again
on a frozen shape.

Reading a shape
---------------

Queries follow ADR-0036 and ADR-0037. Each returns
``Result<Range>``: an empty shape, or a handle that is not a live entity
of the shape, gives ``InvalidInput``.

=====================================  ====================================================
Query                                  Result
=====================================  ====================================================
``AllEntities<K>(shape)``              every live entity of kind K, once, in slot order
``LoopsOf(shape, face)``               the face's loops, in the order they were added
``CoedgesOf(shape, loop)``             one walk around the loop from its anchor
``RadialCoedgesOf(shape, edge)``       every use of the edge, a seam's two included
``FacesOf(shape, shell)``              the shell's faces, in the order they were added
``ShellsOf(shape, solid)``             the solid's shells, in the order they were added
=====================================  ====================================================

The ranges are forward ranges of handles by value; they are not writable,
not bidirectional and not borrowed. A range holds one snapshot ownership
token, taken when it is made, and keeps the snapshot alive after the Shape
is gone; its iterators borrow it, so a step costs O(1) with no allocation
and no reference counting. Use a temporary range in a range-based for, or
keep the range; an iterator taken from a temporary range dangles. A
radial ring's order is a connectivity order, not an angular one, and its
anchor need not be the primary coedge.

O(1) accessors read one record: ``PositionOf``, ``ToleranceOf``,
``KindOf``, ``StartOf``, ``EndOf``, ``RangeOf``, ``PrimaryCoedgeOf``,
``CurveOf``, ``EdgeOf``, ``LoopOf``, ``StoredEdgeUseSense``, ``PCurveOf``,
``FaceOf``, ``SurfaceOf``, ``IsReversed``, ``ShellOf`` and ``SolidOf``.
Their handle must name a live entity of the shape; checked builds report
one that does not, and handles from a range carry their snapshot's store
there, so a handle of another shape is caught. ``Describe(shape, h)`` gives
the text logs and reports use, e.g. ``Edge#0 [Exact Line, tol 1e-07, faces
0|2]``.

The structural validator
------------------------

The validator (ADR-0018's ``Structural`` level) has its own walks over the
raw tables. It checks each link for bounds and liveness before following
it and bounds every walk by the size of its table, so it terminates on any
input; a finding stops the walk that made it and the others go on. The
report keeps the first ``ValidityReport::MaxFindings`` findings and counts
the rest.

====  ======================  ===================================================================
ID    Rule                    Finding
====  ======================  ===================================================================
S01   ``DanglingLink``        a link names a slot out of range or dead, or a required one is missing
S02   ``LoopLinksBroken``     ``Next`` and ``Prev`` are not mutual inverses
S03   ``LoopOpen``            a loop walk reaches another loop's coedge or revisits one early
S04   ``CoedgeOutsideLoop``   a live coedge is not reached by its loop's walk
S05   ``RadialRingOpen``      a radial walk reaches a use of another edge or revisits one early
S06   ``CoedgeOutsideRing``   a live coedge is not reached by its edge's radial walk
S07   ``EdgeUseCount``        the use count does not fit the kind (ADR-0009), or a two-use
                              ``Boundary`` edge is not a seam
S08   ``PrimaryCoedge``       a surface-defined edge lacks its primary, or it is another edge's use
S09   ``LoopOutsideFace``     a loop is not in its face's list exactly once
S10   ``FaceOutsideShell``    a face is not in its shell's list exactly once
S11   ``ShellOutsideSolid``   a shell is not in its solid's list exactly once
S12   ``MissingGeometry``     a face without a surface, an edge without the curve its kind needs,
                              a use without its pcurve, an empty parameter range
S13   ``EdgeEndpoints``       a degenerate edge with two different vertices
S14   ``EmptyLoop``           a loop without coedges
S15   ``ShellDisconnected``   the faces of a shell are not connected through shared edges
====  ======================  ===================================================================

``Testing::TopologyFaults`` writes link fields directly, for the
validator's seeded-fault tests (ADR-0027, section 8); production code has
no other way to reach the tables.

.. todo::

   Document the maintained vertex-incidence index (``EdgesOf``,
   ``FacesOf(edge)``), editing an existing shape, ``Compact()`` and the
   Euler operators as they land.

Related design decisions
------------------------

- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
- `ADR-0005: Immutability and structural sharing <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0005-immutability-and-sharing.md>`_

