.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Introduction
============

Stereon is a geometric modelling kernel: a library that represents solid
shapes exactly enough to manufacture and analyse them, and that builds,
queries and modifies them. It is written from scratch in C++26 and
represents solids by their boundaries, the *boundary representation* or
B-rep that every mainstream CAD kernel uses. This chapter says what such a
kernel does, which guarantees Stereon makes and how it keeps them, what
exists today, and how the rest of the manual is organised.

The name is the ancient Greek *stereón*, "a solid", the word Euclid uses
in the eleventh book of the *Elements*:

   *A solid is that which has length, breadth and depth. And the boundary
   of a solid is a surface.* — Euclid, *Elements* XI, Definitions 1–2

The second definition is the idea this kernel is built on.

What a modelling kernel does
----------------------------

A B-rep describes a solid by the surfaces that bound it, stitched together
along curves that meet at points (:numref:`fig-intro-brep`). Two kinds of
information are kept apart:

- **Topology** says how the pieces connect: which *faces* bound which
  *shell*, which *edges* bound which face, as cycles of oriented edge uses
  called *coedges* grouped into *loops*, and which *vertices* end which
  edge. It is combinatorial and exact.
- **Geometry** says where the pieces are: the *surface* a face lies on,
  the *curve* an edge follows, the *point* at a vertex, and, for every use
  of an edge, its *pcurve*, the same curve drawn in the parameter plane of
  the adjacent surface. It is numerical, and so approximate.

.. _fig-intro-brep:

.. figure:: _static/tikz/intro-brep.png
   :width: 92%
   :alt: A box drawn as its boundary: a shaded top face with its loop of
         coedges, a highlighted edge and a vertex, each linked to the plane,
         line and point that carry it.

   A solid as its boundary. The topology (left) records that the top face
   is bounded by one loop of four coedges, each a use of an edge running
   one way or the other, and that edges end at vertices. The geometry
   (right) records where each entity is. A face refers to its surface, an
   edge to its curve and, through its coedges, to its pcurves, and a vertex
   to its point, which carries a tolerance.

Keeping the two apart is what makes a B-rep robust and economical. The
topology can be validated by counting and walking, with no arithmetic at
all. Geometry can be shared: a thousand faces may lie on one plane. And a
small numerical disagreement between a curve and the surface it should lie
on becomes a measured, bounded quantity, the entity's *tolerance*, instead
of a broken model. :doc:`topology` explains the separation in detail.

On this representation a kernel offers five families of operations:

- **evaluation and query**: points, derivatives, normals, curvature, arc
  length and the projection of a point onto a curve or a surface
  (:doc:`geometry`, :doc:`nurbs`, :doc:`properties`);
- **construction**: from analytic primitives, from fitted data, and by
  sweeping, lofting and sewing (:doc:`fitting`, :doc:`swept`, :doc:`modeling`);
- **intersection**: of curves with curves, curves with surfaces and surfaces
  with surfaces, the computation every other modelling operation rests on
  (:doc:`intersection`);
- **modification**: Boolean union, cut and common, and later blends and
  offsets, built from intersection, imprinting and classification
  (:doc:`imprinting`, :doc:`modeling`);
- **downstream services**: meshes for display and analysis, data exchange,
  and the derivatives of all of the above with respect to design parameters
  (:doc:`meshing`, :doc:`data_exchange`, :doc:`sensitivity`).

What Stereon guarantees
-----------------------

Mature kernels carry decades of design history: reference-counted object
graphs, global state, and tolerance behaviour that is hard to predict.
Stereon is narrower, and sharper about what it promises. Each promise below
is a recorded architecture decision (ADR), and each is enforced by code that
is described in this manual.

**A valid result or a structured error.**
   Every operation that can fail returns ``Result<T>``, a
   ``std::expected<T, KernelError>``. A shape leaves an operation only after
   the independent validity checker of ADR-0018 has passed it, and a failure
   names its cause and the entities involved. Nothing returns a silently
   broken shape, and nothing throws (:doc:`errors`, :doc:`topology`).

**An empty answer is a proof.**
   When an intersection reports that two pieces of geometry do not meet,
   that is a proved separation, not a failure to find a branch. Anything
   that cannot be resolved is reported as unresolved (ADR-0029,
   :doc:`intersection`).

**Exact decisions.**
   Every combinatorial decision, such as which side of a plane a point lies
   on, or whether two lines cross, is the sign of a polynomial in the input
   doubles, computed with floating-point filters and adaptive exact
   arithmetic, so that it is never wrong (ADR-0002, :doc:`predicates`).

**An explicit tolerance model.**
   A model resolution, the smallest distance the model distinguishes, plus
   local tolerances on edges and vertices that grow only through named,
   logged operations (ADR-0001, :doc:`tolerances`).

**Determinism.**
   The same input gives the same output, bit for bit, on any number of
   threads (ADR-0007).

**Value semantics and immutability.**
   Geometry types are plain structs with no base class, and a ``Shape`` is
   immutable once frozen. Edits go through a ``ShapeBuilder`` whose
   copy-on-write pages make snapshots and undo cheap, and shapes are safe to
   share between threads (ADR-0003, ADR-0004, ADR-0005, :doc:`architecture`).

**Declared differentiability.**
   Every operation declares which derivatives it provides and under which
   assumptions. Where it cannot differentiate, it says so with a status,
   never a silent zero (ADR-0021, :doc:`sensitivity`).

Where things stand
------------------

Stereon is pre-alpha software at version 0.1.0, the v0.1 "Geometry" release.
Its APIs change between minor versions until 1.0. The kernel is organised as
eleven libraries in a strict dependency order, each depending only on those
before it (:doc:`architecture`). Their state on ``main`` is:

.. list-table::
   :header-rows: 1
   :widths: 14 58 28

   * - Library
     - What it provides today
     - Chapters
   * - ``core``
     - Vectors, matrices, intervals and boxes; typed handles and paged
       arenas; ``Context``; errors and contracts; the parallel scheduler and
       SIMD wrappers; the sensitivity foundation.
     - :doc:`architecture`, :doc:`errors`
   * - ``robust``
     - Floating-point expansions, filtered exact predicates and interval
       arithmetic.
     - :doc:`predicates`
   * - ``curve``
     - Lines, circles, ellipses, hyperbolas and parabolas in 2D and 3D;
       B-spline and NURBS curves with their operations (knot insertion and
       removal, degree elevation, splitting, concatenation, periodic form);
       trimmed curves; fitting; differential properties; batch evaluation;
       derivatives with respect to poles and weights.
     - :doc:`analytic`, :doc:`nurbs`, :doc:`fitting`, :doc:`properties`
   * - ``surface``
     - Planes, cylinders, cones, spheres and tori; B-spline and NURBS
       surfaces; extrusions and revolutions; trimmed surfaces; projection;
       grid evaluation; certified deviation bounds.
     - :doc:`geometry`, :doc:`trimmed`, :doc:`swept`
   * - ``intersect``
     - The intersection services: analytic surface pairs (plane, cylinder,
       skew cylinders), lines and B-spline curves against planes, spheres
       and cylinders, plane curves for face clipping, and a proved-separation
       broad phase for B-spline surfaces.
     - :doc:`intersection`
   * - ``topo``
     - The B-rep store, the builder and its Euler operators, compaction,
       operation results with history, and the validity checker at the
       Structural, Geometric and Solid levels.
     - :doc:`topology`
   * - ``build``
     - Box and cylinder solids.
     - :doc:`modeling`
   * - ``boolean``
     - Union, cut and common for the solids of the Boolean prototype, whose
       faces are planes and cylinders: face clipping, imprinting and
       classification.
     - :doc:`imprinting`, :doc:`modeling`
   * - ``blend``
     - Nothing yet: fillets, chamfers and shelling are planned for v0.9.
     -
   * - ``mesh``
     - Display sampling of curves and surfaces.
     - :doc:`meshing`
   * - ``io``
     - The native ``.stn`` archive for geometry, and debug dumps.
     - :doc:`data_exchange`, :doc:`debugging`

A Python module, ``stereon``, builds with the kernel and reports its version
(:doc:`python`). A debug viewer displays any intermediate object
(:doc:`debugging`).

The v0.2 "Intersections" release is in progress. Its acceptance ledger in
``docs/planning/v0.2-acceptance-ledger.md`` tracks each intersection
service against its gate. The full release sequence, up to a stable kernel
at v1.0, is in the roadmap (``docs/planning/roadmap-and-versioning.md``).
Out of scope for v1.0 are parametric feature history, sheet metal, drafting
and variable-radius blends.

How this manual is organised
----------------------------

The manual is written to be read in order, though each chapter states what
it assumes.

- **Getting started** builds the kernel (:doc:`building`) and runs worked
  examples (:doc:`examples`).
- **Foundations** cover how the kernel is put together
  (:doc:`architecture`), how it treats numbers: tolerances
  (:doc:`tolerances`), exact predicates (:doc:`predicates`), and errors and
  contracts (:doc:`errors`).
- **Geometry** develops curves and surfaces from the analytic types
  (:doc:`geometry`, :doc:`analytic`) through B-splines and NURBS
  (:doc:`nurbs`) to trimmed and swept geometry (:doc:`trimmed`,
  :doc:`swept`). It continues with conversion between them
  (:doc:`interchange`), fitting (:doc:`fitting`), differential properties
  (:doc:`properties`) and design sensitivities (:doc:`sensitivity`), and
  ends with how all of it is made fast (:doc:`performance`).
- **Topology and operations** cover the B-rep store and its validity
  (:doc:`topology`), intersection (:doc:`intersection`), modelling
  (:doc:`modeling`) and imprinting (:doc:`imprinting`), then meshes and
  files (:doc:`meshing`, :doc:`data_exchange`).
- **Reference** has the Python bindings (:doc:`python`), the API by library
  (:doc:`api`), how the kernel is tested and debugged (:doc:`testing`,
  :doc:`debugging`), and the bibliography (:doc:`references`).

Each algorithm is named where it is used, with the literature it comes from
and the file that implements it.

Conventions
-----------

**Units.**
   Lengths are millimetres and angles radians (ADR-0001). A tolerance is a
   length, and the model resolution is :math:`10^{-5}` mm by default.

**Notation.**
   A curve is :math:`C(t)` on a parameter interval :math:`[a, b]`, and a
   surface :math:`S(u, v)` on a parameter box. A pcurve is
   :math:`c(t) = (u(t), v(t))`, and :math:`S(c(t))` is the curve it draws on
   the surface. Primes and subscripts denote derivatives: :math:`C'(t)`,
   :math:`S_u`, :math:`S_{uv}`. Vectors are upright (:math:`P`, :math:`d`)
   and scalars italic.

**Code.**
   Names in ``monospace`` are C++ identifiers in namespace ``Stereon``.
   Files are given from the repository root, as in
   ``src/curve/BSplineCurve.cppm``. Each library is a C++20 module named
   after it (``Stereon.Core``, ``Stereon.Curve``, ..., ``Stereon.Io``),
   which re-exports the libraries below it, so examples ``import`` the
   highest one they use.

**Decisions.**
   "ADR-*nnnn*" refers to the architecture decision record of that number in
   ``docs/adr/``, which states the decision, its alternatives and its
   verification. The manual explains what the code does; the ADRs explain
   why.

Licence and contributions
-------------------------

Stereon is free software under the Mozilla Public License 2.0 (ADR-0008).
Contributions follow ``CONTRIBUTING.md``. Every change comes with a test that
fails without it, and AI-assisted changes to the intersection, Boolean and
blend libraries are reviewed line by line (ADR-0015).

Related design decisions
------------------------

- `ADR-0001: Tolerance model <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0001-tolerance-model.md>`_
- `ADR-0002: Numeric robustness <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0002-numeric-robustness.md>`_
- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0007: Threading and determinism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0007-threading-and-determinism.md>`_
- `ADR-0008: Licence <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0008-licence.md>`_
- `ADR-0015: AI contribution policy <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0015-ai-contribution-policy.md>`_
- `ADR-0018: Validity checking <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0018-validity-checking.md>`_
- `ADR-0021: Shape derivatives <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
- `ADR-0029: Intersection and classification outcomes <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0029-intersection-and-classification-outcomes.md>`_
