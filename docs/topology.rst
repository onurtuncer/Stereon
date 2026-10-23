.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

B-rep topology
==============

The topology store, typed handles, ShapeBuilder and immutable Shape,
editing and compaction, and the validity checker.

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

- A **vertex** records its position and tolerance, and the head of its
  incidence list: the edges that start or end at it, in ascending index
  order, each continuing the list through its own ``NextAtStart`` or
  ``NextAtEnd``, whichever end the vertex is. This is the maintained
  vertex-incidence index of ADR-0004 and ADR-0036; the builder keeps it
  as edges come and go, and ``EdgesOf`` reads it without scanning.
- An **edge** records its curve kind, its tolerance, its start and end
  vertices, the shared parameter range of its curve and pcurves, the anchor
  of its radial ring, its two incidence links and, for the kinds that
  evaluate through a surface, its primary coedge (ADR-0009). ``Exact`` and
  ``Spatial`` edges carry their own 3D curve; ``Intersection``, ``Trim``,
  ``Boundary`` and ``Degenerate`` edges have none and are evaluated on the
  primary coedge's pcurve and surface.
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
   const SurfaceRef plane = *b.AddSurface(Surface(Plane{...}));
   const BuilderFaceId face = *b.AddFace(plane);
   const BuilderLoopId loop = *b.AddLoop(face);
   const BuilderEdgeId edge = *b.AddEdge({.Kind = EdgeCurveKind::Exact,
                                          .Start = v0, .End = v1,
                                          .Curve = *b.AddCurve(Curve(line)),
                                          .Range = Interval(0, 1)});
   *b.AppendCoedge(loop, edge, EdgeUseSense::Same, *b.AddPCurve(pcurve));
   Result<Shape> shape = Freeze(std::move(b));

``AppendCoedge`` appends to the loop's cycle and inserts the coedge into its
edge's radial ring; the first use of a surface-defined edge becomes its
primary coedge, which only ``SetPrimaryCoedge`` changes. ``AddLoop``,
``AddShell`` and ``AddFaceToShell`` append to their parent's list, so
members come back in the order they were added. Every operation checks its
handles and its request and reports ``InvalidInput`` without changing
anything: a handle that is not a live entity of this builder, an empty
curve, pcurve or surface offered to a geometry pool, an empty
parameter range, a missing pcurve, an edge kind given a curve it does not
have, a face added to a second shell. Every change advances the builder's
topology epoch; its ``BuilderToken`` is stable across moves, and a
moved-from builder refuses every operation (ADR-0038).

Editing a shape
---------------

``ShapeBuilder::Open(shape)`` gives a builder on a frozen shape's tables.
It shares every page with the shape until it writes to one (ADR-0005), so
an edit costs the pages it touches, not the size of the body; the shape
itself never changes. ``builder.Handle(id)`` is the builder's handle to
one of the shape's entities, the same slot; new entities take slots after
the shape's last one, and a removed slot is never reused by the builder,
so an index live in both the shape and its edit names the same entity,
edited or not (ADR-0004). ``SetSurface(face, surface)`` replaces a face's
surface by another of the pool.

Deletion is the family of ``Remove`` operations, each keeping every
reciprocal link consistent and leaving a tombstone: the slot is dead, its
generation bumped, so every handle to it is stale and every operation
refuses it.

- ``RemoveCoedge`` takes one use out of its loop's cycle and its edge's
  radial ring, moving the anchors on. It is refused when the coedge is its
  edge's primary and the edge has other uses: reassign the primary first
  with ``SetPrimaryCoedge`` (ADR-0027, section 4).
- ``RemoveLoop`` removes the loop's coedges and the loop from its face;
  ``RemoveFace`` removes the face's loops that way and the face from its
  shell. Every use of the loop goes at once, so a primary among them
  passes to the next use in the ring; an edge left with no use loses its
  primary. The edges and vertices stay. An edge whose kind needs a use
  (``Boundary``, ``Trim``, ``Intersection``) then fails the validator
  until it is removed or re-kinded, which keeps ADR-0009's rules in force.
- ``RemoveShell`` and ``RemoveSolid`` detach their children, which stay
  in no shell or of no solid.
- ``RemoveEdge`` is refused while a coedge uses the edge, and
  ``RemoveVertex`` while an edge starts or ends at the vertex: nothing is
  ever left pointing at a dead slot.

.. code-block:: cpp

   ShapeBuilder b = *ShapeBuilder::Open(shape);
   *b.RemoveFace(b.Handle(face));
   Result<Shape> edited = Freeze(std::move(b));   // shares every untouched page

Euler operators
---------------

The controlled edits of ADR-0027 are the classical Euler operators,
written for this store: radial rings with any number of uses, explicit
primary coedges, and pcurves on the edge's shared parameter. Each checks
all of its preconditions before it writes anything, so a refusal
(``InvalidInput``) leaves the builder as it was. Each has an inverse that
restores the connectivity up to the entities that survive: removed slots
are tombstones, so what an inverse brings back has new handles.

- ``MakeVertexFaceShell`` (MVFS) makes an isolated vertex, and a face with
  one empty loop in a new shell. ``KillVertexFaceShell`` (KVFS) removes
  them again.
- ``MakeEdgeVertex`` (MEV) grows a new vertex and an edge to it from the
  vertex where a coedge ends, as a strut of that coedge's loop: the edge's
  use along it and its use against it follow the coedge. On an empty loop
  it grows from the loop's isolated vertex. ``KillEdgeVertex`` (KEV)
  removes a strut and its far vertex.
- ``MakeEdgeFace`` (MEF) joins the vertices where two coedges of one loop
  end, splitting the loop: a new face on the same surface and in the same
  shell takes one side. ``KillEdgeFace`` (KEF) removes an edge between two
  such faces and one of the faces, joining its loop into the other.

The KEMR/MEKR, KFMRH/MFKRH and SEMV/JEKV pairs and ``MergeVertices``
are described in `M5 ring, split and merge operators`_ below.

.. code-block:: cpp

   ShapeBuilder b;
   auto [v0, face, loop, shell] = *b.MakeVertexFaceShell(origin, 0.0, plane);
   auto e1 = *b.MakeEdgeVertex(loop, v0, p1, 0.0, edge);    // v0 -> v1
   auto e2 = *b.MakeEdgeVertex(e1.Out, p2, 0.0, edge);      // v1 -> v2
   auto e3 = *b.MakeEdgeFace(e2.Out, e1.Back, edge);        // v2 -> v0: two faces

``Transaction(edit)`` runs a compound edit on a candidate copy of the
builder, which shares every page with it until it writes to one. The copy
replaces the builder only if the edit succeeds; on failure it is discarded,
and the builder is exactly as it was. Handles are valid in both. In
checked builds every operator also checks the loops, radial rings and
vertex lists it touched, a contract, while ``Freeze`` and ``Commit`` run
the independent validator in every build.

An edit composed of operators reads the connectivity it has made so far
back from the builder: ``EdgeOf``, ``LoopOf`` and ``SenseOf`` a coedge,
``StartOf`` and ``EndOf`` an edge, ``PositionOf`` a vertex, ``LoopsOf`` a
face and ``CoedgesOf`` a loop, each as handles at the current generation
(an invalid handle, or an empty list, for one that is not live). Geometry
beyond positions is read from the frozen ``Shape`` the edit started from;
the Booleans' imprinting (:doc:`modeling`) finds the two uses an MEF joins
this way after the splits before it.

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
``EdgesOf(shape, vertex)``             the edges at the vertex, once each, ascending
``FacesOf(shape, edge)``               the distinct faces of the edge's uses, ascending
``FacesOf(shape, shell)``              the shell's faces, in the order they were added
``ShellsOf(shape, solid)``             the solid's shells, in the order they were added
=====================================  ====================================================

``EdgesOf`` walks the maintained vertex-incidence index: a vertex of degree
d costs O(d), one record per step, no allocation, and a closed edge comes
once. ``FacesOf(shape, edge)`` is the eager result ADR-0036 allows for
distinct adjacency: it walks the radial ring, sorts the faces and drops
repeats, so k uses cost O(k log k) and one allocation, visible in its
``std::vector<FaceId>`` return type; a seam's face comes once.

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
``CurveOf``, ``CachedCurveOf``, ``CacheCertificateOf``, ``CacheBoundOf``,
``EdgeOf``,``LoopOf``, ``StoredEdgeUseSense``, ``PCurveOf``,
``FaceOf``, ``SurfaceOf``, ``IsReversed``, ``ShellOf`` and ``SolidOf``.
Their handle must name a live entity of the shape; checked builds report
one that does not, and handles from a range carry their snapshot's store
there, so a handle of another shape is caught. ``Describe(shape, h)`` gives
the text logs and reports use, e.g. ``Edge#0 [Exact Line, tol 1e-07, faces
0|2]`` or ``Vertex#3 [(0, 3, 0), tol 1e-07, 3 edges]``.

Compacting
----------

``Freeze`` keeps every index, so a long edit sequence leaves tombstones
behind; ``Compact(shape)`` is the explicit renumbering of ADR-0004. It
gives the live entities of each kind consecutive indices in their old
order, rewrites every link, drops the dead slots and the geometry nothing
references, and freezes the result, so the compacted shape is certified
like any other. The ``Compacted`` result carries the new ``Shape`` and the
``RemapTable``, the old-to-new index of every slot of every kind and of
the three pools (``Invalid`` for a dropped one), which the history of
ADR-0011 composes; ``compacted.Translate(id)`` is the compacted shape's
handle to an entity of the source, invalid for one that was dead. The
source shape is left as it is, and its handles do not reach the compacted
shape: checked builds tell the two stores apart. Traversal order is
preserved by the renumbering but not promised across it (ADR-0037).

Cached edge curves
------------------

An ``Intersection`` or ``Trim`` edge has no 3D curve of its own: it is
``S(c(t))``, its primary coedge's pcurve on that face's surface (ADR-0009).
ADR-0009 gives such an edge a cached B-spline approximation, for bounding
boxes, display and quick distance queries: a cache, never the truth.

``builder.CacheEdgeCurve(edge)`` builds it. It samples ``S(c(t))`` at
``ShapeBuilder::CacheSamples`` parameters spread evenly over the edge's
range and interpolates them with a cubic B-spline on the edge's own
parameter, so the cache passes through the primary side at every sample.
The same pcurve and surface give the same cache. ``CachedCurveOf(shape,
edge)`` reads the cache, ``CacheCertificateOf`` what is known of its
deviation, and ``CacheBoundOf`` the proven bound.

ADR-0009 asks for a *certificate*: a bound on the deviation over the whole
range, proven by interval evaluation of ``S(c(t))``. ``CacheEdgeCurve``
asks ``CertifyDeviation`` (see :doc:`intersection`) for one within the
edge's tolerance, and the cache is ``CacheCertificate::Certified`` with the
proven bound when it gets one. It stays ``CacheCertificate::Uncertified``,
with an infinite bound, when the edge has no tolerance to certify within,
when the certifier does not cover the surface or pcurve types yet, or when
the sampled cache deviates by more than the tolerance. Sampling the
deviation is a regression check in the tests, never a certificate.

Any change to what a cache was built from drops it, as ADR-0004 requires:

- a new use of the edge (``AppendCoedge``);
- a use removed (``RemoveCoedge``, ``RemoveLoop``, ``RemoveFace``);
- the primary coedge moved (``SetPrimaryCoedge``);
- a new surface on a face that uses the edge (``SetSurface``).

``Commit`` reports each dropped cache as a ``CacheDropped`` warning.
``Compact`` keeps caches and remaps them. Rule S17 checks that an edge's
cache fields agree with one another and with its kind.

Operation results
-----------------

A shape-producing operation returns an ``OperationResult`` (ADR-0006):

- ``Shape``, the result;
- ``History``, what became of the input's entities;
- ``ToleranceChanges``, the tolerances it changed;
- ``Warnings``, structured warnings that each carry a code, a message and
  the entities it names;
- ``Checked``, the validity level the result was checked at.

Everything is published together, and only after validation.
``Commit(std::move(builder), context)`` and ``CompactShape(shape,
context)`` return one. ``Freeze`` and ``Compact`` remain their building
blocks.

The History (ADR-0011, ADR-0030) covers one input. An entity of the input
continues ``Unchanged`` or ``Modified`` at one index of the result, or is
``Deleted``. An entity of the result with no source is ``Generated``: the
operation's creation root. Each entry carries the input's ordinal, the
kind and both indices. The entries are ordered by kind, input, source and
target, so two runs give the same History. ``History::OfSource`` and
``OfTarget`` look an entity up.

For ``Commit`` the History compares the result with the shape the builder
was opened on; a new builder's entities are all ``Generated``. An edit
keeps every index, so a slot live in both continues. It is ``Unchanged``
when the entity's own fields and what it bounds are as they were, and
``Modified`` otherwise. The links that only chain siblings in a parent's
list don't count, so a face is not modified because its neighbour was
deleted, but its shell is. ``CompactShape`` composes the remap: every
live entity continues ``Unchanged`` at its new index. Split and merge
events join with the Euler operators that make them.

``builder.SetTolerance(entity, value, reason)`` is the explicit operation
through which ADR-0001 lets a local tolerance change. ``Commit`` reports
one change per entity: the first old value, the last new one, and the
reasons in order. It warns with ``ToleranceGrown`` for a tolerance that
grew. It fails with ``ToleranceExceeded`` for a tolerance above the
``Context``'s ceiling, and for any growth when the ``Context`` sets
``ToleranceGrowth::Fail``.

The log is not the only check. ``Commit`` and ``CompactShape`` hold every
live vertex and edge of what they publish to ADR-0001's bounds under their
``Context``, including entities made with their tolerance and entities kept
unchanged from a shape published under another ``Context``. A tolerance
above the ceiling is ``ToleranceExceeded``, and one below the resolution is
``InvalidInput``. A model built at a finer resolution says so with
``Context::WithResolution``.

Validation happens where a shape is published: ``Freeze``, ``Commit``,
``Compact`` and ``CompactShape``. A builder's single steps are not
validated, because a valid edit passes through states that are not valid;
removing a face leaves its edges' other uses alone until the edit
continues. Publication always runs the structural validator. Checked
builds (the debug, ASan and fuzzing presets) also validate the input
again at ``ShapeBuilder::Open``, ``Commit`` and ``CompactShape``, rather
than trust its certificate (ADR-0004, safeguard 5). Before an operation
publishes, it checks that its History accounts for every live entity of
both sides and names only live ones; a gap is ``InvalidResult``, a kernel
bug.

``Commit`` and ``CompactShape`` also run the Geometric level (below) when
the ``Context``'s validity level asks for it, which is the default, and
always in checked builds. A ``Commit`` result that fails it is
``InvalidResult``. A shape ``CompactShape`` is given that fails it is
``InvalidInput``, because compaction changes no geometry. ``Checked`` says
which level ran. ``Freeze`` certifies structure only, so a builder can
publish geometry that is still being repaired. The Solid level does not
exist yet; asking for it checks the Geometric level.

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
                              a use without its pcurve, an empty or unbounded parameter range
S13   ``EdgeEndpoints``       a degenerate edge with two different vertices
S14   ``EmptyLoop``           a loop without coedges
S15   ``ShellDisconnected``   the faces of a shell are not connected through shared edges
S16   ``VertexIncidence``     a vertex's list names an edge not at it, twice or out of order,
                              or an edge is missing from one of its ends' lists
S17   ``EdgeCache``           a cached curve on a kind that keeps none or naming no pool entry,
                              a cache status without a cache, a bound its status does not allow
====  ======================  ===================================================================

The geometric validator
-----------------------

The Geometric level checks ADR-0001's containment invariant on a
structurally valid shape (``ValidateTopology(shape,
ValidityLevel::Geometric)``):

====  ======================  ===================================================================
ID    Rule                    Finding
====  ======================  ===================================================================
G01   ``VertexContainment``   an end of an edge lies outside its vertex's tolerance sphere
G02   ``EdgeContainment``     a pcurve on its surface, or the edge's 3D curve, leaves the edge's
                              tolerance tube; a trace does not cover the edge's range or a
                              pcurve leaves its surface's domain; or a ``Degenerate`` edge's
                              trace moves away from its vertex
====  ======================  ===================================================================

A ``Degenerate`` edge is a pole or an apex with no 3D extent. Each of its
traces must stay at its vertex over the whole range, within the larger of
the edge's and the vertex's tolerances. Ends that coincide with the vertex
are not enough, because a closed circle has them too. An empty ``Shape``,
default-constructed or moved from, is valid at every level.

Each edge has a reference trace: its own 3D curve, or for a kind without
one, its primary coedge's pcurve on that coedge's surface (ADR-0009). Every
other trace is compared with the reference at 17 evenly spread parameters,
ends included, on the edge's shared parameter. A gap between samples can
pass, so the check is sampled, not certified. A trace that cannot be
evaluated over the whole range is a finding, never a contract violation.
A local tolerance large enough to hold a gap makes the same geometry
valid, which is what ADR-0001's local tolerances are for.

``Testing::TopologyFaults`` writes link fields directly, for the
validator's seeded-fault tests (ADR-0027, section 8); production code has
no other way to reach the tables.

The seeded-fault suite of ADR-0004 (``tests/topo/FaultTests.cpp``) feeds
stale handles, handles from another builder or shape, off-by-one indices
and tombstones to every operation and query, mixes entity kinds at compile
time, and corrupts loop links, rings, slots and the incidence index: each
is refused through a ``Result`` naming the entity with ``Describe``, or by
``Freeze`` with the rule. The O(1) accessors report the same faults as a
contract violation in checked builds.

Related design decisions
------------------------

- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
- `ADR-0005: Immutability and structural sharing <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0005-immutability-and-sharing.md>`_
- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0009: Parameter-space representation <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0009-parameter-space-representation.md>`_
- `ADR-0030: Entity lineage in operation results <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0030-entity-lineage-in-operation-results.md>`_


M5 ring, split and merge operators
---------------------------------

The remaining classical operators are exposed on ``ShapeBuilder``. They
change incidence and membership; callers retain responsibility for geometric
nesting and face classification. None infers that a loop is a hole from its
orientation. Use ``Transaction`` for compound edits and ``Commit`` to publish
with structural/geometric validation and tolerance policy checks.

.. list-table:: Entity-count changes (inverse operations negate these)
   :header-rows: 1

   * - Operation
     - Vertices
     - Edges
     - Coedges
     - Loops
     - Faces
   * - ``KillEdgeMakeRing`` (KEMR)
     - 0
     - -1
     - -2
     - +1
     - 0
   * - ``MakeFaceKillRingHole`` (MFKRH)
     - 0
     - 0
     - 0
     - 0
     - +1
   * - ``SplitEdgeMakeVertex`` (SEMV), radial valence k
     - +1
     - +1
     - +k
     - 0
     - 0

KEMR requires two opposite, nonadjacent uses in one loop and produces two
nonempty loops. MEKR bridges the ends of uses in two loops of the same face.
KFMRH moves a face's sole nonempty loop to a face with identical surface,
orientation and shell; MFKRH moves one loop of a multi-loop face to a new face.
When the face belongs to a shell, MFKRH first checks for a radial-edge path
between the moved and remaining loops through faces of that shell. It refuses
to disconnect the shell, before changing any records. Faces outside a shell
have no such connectivity constraint.
These are incidence operations: genus and connected-component interpretations
require a manifold, geometrically classified input and are not inferred for
arbitrary radial topology.

SEMV requires a strictly interior parameter on a nondegenerate edge. Its new
vertex must contain every incident trace at that parameter. It restricts the
shared curve and pcurve parameter ranges, preserves certified cache bounds by
restriction, and duplicates every radial use, including seam sides and
non-manifold uses. JEKV requires adjacent uses with matching pcurves and primary
sides, contiguous ranges of identical geometry, and a degree-two join vertex.
It retains a common cache with the larger bound, or drops differing caches.

``MergeVertices`` sorts unique candidates by slot index and clusters against
the first candidate's original position and tolerance. It never averages or
chains through a previous cluster member. The representative's tolerance is
the maximum of member tolerances and the context resolution; the exact pair
predicate guarantees that this contains the cluster radius. Growth is logged,
the context ceiling and growth policy are enforced, and every incident endpoint
trace is checked before redirecting edges. Unsupported exact-predicate domains
are explicit failures. The operation uses a candidate transaction, so failed
merges leave topology and event history unchanged (rollback still reserves any
allocated identities and advances the epoch as specified by ADR-0040).

``GetTopologyEvents`` and ``OperationResult.TopologyEvents`` expose low-level
created, deleted, modified, split and merged relations from the Euler edits.
Indices may refer to intermediate entities and are not snapshot handles. These
events survive successful transactions and disappear on rollback; the enclosing
Boolean operation still supplies final multi-input semantic history in M8.

Regression tests cover inverse pairs, invalid preconditions, seam and closed
edges, arbitrary radial valence, fixed merge representatives and event rollback.
Native GCC 16 debug/sanitizer qualification and benchmark measurements remain
required before closing the M5 milestone.
