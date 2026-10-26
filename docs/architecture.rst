.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Architecture overview
=====================

The library layers under src/, the dependency rules between them, the Context object, and how the design decisions fit together.

The Context
-----------

Every operation takes a ``const Context&`` besides its inputs, and a result
is a function of the two alone: the kernel has no global or thread-local
state (ADR-0007). The ``Context`` carries the tolerance settings of ADR-0001
(``GetResolution``, ``GetAngularResolution``, ``GetToleranceCeiling``, the
±500 m ``GetSizeBox`` and the policy for tolerance growth), the scheduler
kind and thread count of ADR-0007, the validity level checked after every
operation (ADR-0018), the compaction threshold of ADR-0004, and the two
diagnostic sinks of ADR-0023: a ``LogSink`` that receives formatted log
records above a level threshold, and a ``JournalSink`` that receives the
entry an operation records before it runs and answers with the ``JournalId``
a ``KernelError`` carries for replay.

A ``Context`` is a value. Its settings are fixed at construction and changed
by making a modified copy, ``ctx.WithResolution(1e-4)``, so one context can be
shared by any number of operations on any number of threads. The settings are
the programmer's: a value outside its documented domain is a contract
violation, not a ``KernelError``. Without sinks a context is silent;
``StderrLogSink()`` and ``MemoryJournal`` are the ready-made sinks for tests
and tools, and an application attaches its own logging framework behind a
``LogSink``.

The libraries
-------------

The kernel is eleven libraries under ``src/``, one C++ module each
(``Stereon.Core``, ``Stereon.Robust``, ...), built as the CMake targets
``stereon::core``, ``stereon::robust``, and so on. They form a strict stack:
a library may depend only on the libraries below it, never on one above it or
beside it. The order is the order of the ``add_subdirectory()`` calls in
``src/CMakeLists.txt``, and ``tools/check_layering.py`` fails the lint
workflow on any upward ``import``, ``#include <Stereon/...>`` or CMake target
reference. Each module re-exports the one below it (``export import``), so an
application imports the highest layer it needs and sees the whole stack
beneath it.

.. graphviz::
   :align: center

   digraph layers {
       rankdir = BT;
       node [shape = box, fontname = "Helvetica", fontsize = 10, width = 3.6];
       edge [arrowsize = 0.6];
       core      [label = "core\nvectors, intervals, Result, Context,\narena and handles, SmallBuffer, journal"];
       robust    [label = "robust\nfloating-point filters, expansions,\nexact predicates"];
       curve     [label = "curve\nanalytic and B-spline curves,\nCurve, Curve2d, VisitBuiltin"];
       surface   [label = "surface\nanalytic, B-spline, swept surfaces,\nSurface, pcurve certification"];
       intersect [label = "intersect\ncurve-curve and surface-surface\nintersection, certified outcomes"];
       topo      [label = "topo\nShape, ShapeBuilder, Euler operators,\nCommit, History, validation"];
       build     [label = "build\nprimitive solids"];
       boolean   [label = "boolean\nimprint, classify, assemble"];
       blend     [label = "blend\n(reserved)"];
       mesh      [label = "mesh\ncurve and surface sampling"];
       io        [label = "io\ngeometry archive, dumps, debug JSON"];
       core -> robust -> curve -> surface -> intersect -> topo -> build
            -> boolean -> blend -> mesh -> io;
   }

From the bottom:

``core``
   Everything without geometry: scalars, vectors, matrices, frames, boxes and
   intervals; ``Result`` and ``KernelError`` (ADR-0006); the ``Context`` above;
   the log and journal sinks (ADR-0023); the paged arena and index handles the
   topology is stored in (ADR-0004, ADR-0005); ``SmallBuffer``, the type
   erasure that geometry wrappers are built on (ADR-0003); and serialisation
   primitives.

``robust``
   Adaptive-precision arithmetic: floating-point filters with certified error
   bounds, Shewchuk's expansion arithmetic and the exact orientation and
   in-circle predicates built from them (ADR-0002, :doc:`predicates`).

``curve`` and ``surface``
   The geometry: concrete types (``Line``, ``Circle``, ``BSplineCurve``,
   ``Plane``, ``Cylinder``, ``BSplineSurface``, ...) as plain values with free
   functions (``Eval``, ``Deriv``, ``Project``), the erased wrappers
   ``Curve``, ``Curve2d`` and ``Surface`` that topology stores, and
   ``VisitBuiltin``, which recovers the concrete type once so a loop runs on
   it without a virtual call per evaluation (ADR-0003). Fitting, properties,
   conversion and interchange formats live here too. ``surface`` also
   certifies pcurves: that a parameter-space curve lies on its 3D edge curve
   within the edge tolerance (ADR-0009).

``intersect``
   Curve-curve, planar and surface-surface intersection, dispatched on the
   pair of concrete types to analytic handlers where they exist (ADR-0010).
   Every result says whether it is certified and, if not, why (ADR-0029); see
   :doc:`intersection`.

``topo``
   The boundary representation (ADR-0025): solids, shells, faces, loops,
   coedges, edges and vertices in structure-of-arrays tables addressed by
   index handles (ADR-0004). A ``Shape`` is immutable and shares pages with
   the shapes it was derived from (ADR-0005). A ``ShapeBuilder`` is the only
   way to change one; it applies the Euler operators of ADR-0027, and
   ``Commit`` validates the result, freezes it and returns an
   ``OperationResult``: the new shape, its ``History`` (ADR-0011, ADR-0030),
   the tolerance changes and the warnings. See :doc:`topology`.

``build``
   Primitive solids, ``MakeBoxSolid`` and ``MakeCylinderSolid``, built through
   a ``ShapeBuilder`` and committed like any other operation.

``boolean``
   Union, cut and common (ADR-0013): imprinting both bodies with their
   intersection, classifying every region of the imprinted boundaries against
   the other body, and assembling the kept faces into the result. See
   :doc:`modeling`.

``blend``
   Reserved for fillets and chamfers; it re-exports ``boolean`` and has no
   content yet.

``mesh``
   Tolerance-driven sampling of curves and surfaces for display and export
   (``SampleCurve``, ``SampleSurface``); mesh generation proper is outside the
   kernel (ADR-0014). See :doc:`meshing`.

``io``
   The native geometry archive (ADR-0012), whole-shape archives for the debug
   viewer (ADR-0035), dump directories and debug JSON. See
   :doc:`data_exchange`.

Outside the stack and depending on it, never the other way round: the Python
bindings in ``python/`` (ADR-0020), the debug viewer in ``viewer/``, the
``stn-dump`` tool, the tests and the benchmarks. OCCT is linked only into the
differential tests and benchmarks (ADR-0008).

Why the order
~~~~~~~~~~~~~

The stack follows what each layer has to know. Predicates know numbers but
not geometry; curves need predicates but not surfaces; surfaces are built
from curves (swept surfaces, pcurves), so they come after. Intersection needs
both, and only geometry: it answers questions about curves and surfaces with
no idea of faces or edges. Topology stores curves and surfaces and certifies
the edge curves it caches against them (``CertifyDeviation`` in
``surface``); it sits above intersection rather than beside it so that its
capability catalogue (ADR-0021) can gather those of every geometric layer,
and so that the operations above it see one stack. Modelling operations (``build``, ``boolean``, ``blend``) sit on
topology and change shapes only through its builder. Sampling and file I/O
read finished shapes and come last.

Two consequences matter when adding code. A lower layer cannot call back into
a higher one, so a geometric routine that would like to know about faces
takes the data it needs as arguments instead. And a new kind of geometry
(another surface type, say) touches ``surface``, the dispatch tables in
``intersect`` and the capability tables of the layers above, which
``tools/check_capabilities.py`` keeps consistent.

Data flow
---------

An operation takes immutable shapes and a ``Context`` and returns a new shape
or an error; nothing is changed in place. A Boolean, from primitives to a
mesh, runs like this:

.. graphviz::
   :align: center

   digraph flow {
       rankdir = TB;
       node [fontname = "Helvetica", fontsize = 10];
       edge [arrowsize = 0.6, fontname = "Helvetica", fontsize = 9];
       ctx   [label = "Context\nresolution, tolerances, sinks", shape = note];
       prim  [label = "MakeBoxSolid, MakeCylinderSolid\n(build)", shape = box];
       a     [label = "Shape A", shape = ellipse];
       b     [label = "Shape B", shape = ellipse];
       imp   [label = "Imprint\nface pairs → surface-surface intersection (intersect)\n→ branches clipped to both faces → networks", shape = box];
       cls   [label = "Classify\nregions → inside, outside, OnSame, OnOpposite\nexact rays, winding-number fallback (robust)", shape = box];
       asm   [label = "Assemble\nkept faces → ShapeBuilder → Euler operators", shape = box];
       com   [label = "Commit (topo)\nvalidate, freeze, History", shape = box];
       res   [label = "OperationResult\nShape, History, tolerance changes, warnings", shape = ellipse];
       out   [label = "SampleSurface (mesh)\nShapeArchive (io)", shape = box];
       prim -> a; prim -> b;
       a -> imp; b -> imp;
       imp -> cls -> asm -> com -> res -> out;
       ctx -> imp [style = dashed]; ctx -> cls [style = dashed];
       ctx -> com [style = dashed];
   }

1. **Construction.** ``MakeBoxSolid`` opens a ``ShapeBuilder``, adds the
   surfaces, curves and pcurves, runs the Euler operators that make the
   faces, edges and vertices, and calls ``Commit``. The builder checks each
   operator's preconditions as it goes (ADR-0026, ADR-0039); ``Commit``
   validates the whole body at the context's validity level (ADR-0018) and
   freezes it into a ``Shape``.
2. **Imprinting.** ``Boolean(a, b, op, ctx)`` pairs the faces of the two
   bodies whose boxes overlap and intersects their surfaces. Each
   intersection branch is clipped to both faces' boundaries in parameter
   space and imprinted into both, so the two bodies' boundaries are split
   along the same curves with the same shared parameter. Coincident planar
   faces are imprinted with each other's boundaries instead
   (:doc:`modeling`).
3. **Classification.** Every region of each imprinted boundary is a
   connected set of faces between imprinted curves. A seed point in each is
   classified against the other body by exact ray casting, falling back to a
   certified winding number when every ray is ambiguous; an answer that
   neither settles is an error, not a guess (ADR-0013, ADR-0029).
4. **Assembly.** The operation picks the regions it keeps (union: both
   bodies' outside regions; common: the inside ones; cut: A's outside and
   B's inside, reversed) and copies their faces into a fresh builder, sharing
   the imprinted edges between the two sides.
5. **Commit.** ``Commit`` validates and freezes the result as in step 1 and
   records the ``History``: which input face or edge each output entity came
   from, split from or merged into, for both operands (ADR-0030).
6. **Consumers.** The result is a shape like any other: it can be the input
   of the next Boolean, sampled for display by ``mesh``, or written by
   ``io``.

The ``Context`` flows through every step unchanged. Each step returns a
``Result``; the first error stops the operation, carries the entities and
tolerances involved, and refers to the journal entry that replays it
(ADR-0006, ADR-0023). No step mutates its inputs, so a failed operation leaves
nothing to clean up, and the input shapes stay valid for another attempt with
a different context.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0007: Threading and determinism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0007-threading-and-determinism.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
- `ADR-0023: Logging and journal sinks <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0023-logging-and-journal-sinks.md>`_

