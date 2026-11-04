.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Tolerances
==========

A B-rep kernel computes in floating-point arithmetic on geometry that is
itself approximate, yet it must make yes-or-no decisions about that geometry:
does this edge end at that vertex? Do these two faces meet? Is this curve
on that surface? This chapter sets out how Stereon reconciles the two. It
uses a small model resolution for everything the kernel constructs, local
tolerances where data demands more, an invariant that makes "approximate"
measurable, and a policy for identifying nearby entities that does not
depend on luck. The exact predicates that decide signs without any
tolerance have their own chapter, :doc:`predicates`.

Why a kernel needs tolerances
-----------------------------

Three sources of disagreement are unavoidable.

- **Representation.** A double carries 53 bits. The intersection of two
  planes, written back as a line through doubles, lies on neither plane
  exactly.
- **Construction.** Most constructions are not closed-form. The
  intersection curve of two B-spline surfaces is approximated by a spline,
  and a fitted curve meets its data only to a budget. The result lies near
  the exact one, not on it.
- **Exchange.** A model imported from another system was built to that
  system's tolerances. Gaps of :math:`10^{-3}` mm between the boundaries of
  adjacent faces are common, while Stereon's own constructions are accurate
  to :math:`10^{-7}` mm or better.

Exact arithmetic throughout would remove the first two sources. Some
kernels do this for polyhedra, but it is impractical for curved geometry,
whose intersections are algebraic curves of high degree. A single global
tolerance cannot represent imported data without healing all of it to that
value. Stereon follows the approach Parasolid made standard (ADR-0001): a
fixed **model resolution** governs every entity the kernel builds, and an
edge or a vertex may carry a larger **local tolerance** when its data
requires one.

Three quantities, kept apart
----------------------------

Much of the trouble kernels have with tolerances comes from one number
serving several purposes. Stereon keeps three apart (ADR-0001):

**Numerical error bound**
   How far a computed value may lie from the exact one. It is evidence
   attached to a construction (an intersection point, a fitted curve, a
   cached edge curve), never an adjustable setting.

**Modelling tolerance**
   How far apart two entities may be and still be *treated* as coincident.
   That is the model resolution, and the local tolerance of each vertex and
   edge.

**Healing budget**
   How large a local tolerance may become during one operation before the
   operation fails: the tolerance ceiling.

A construction's error bound must not exceed the modelling tolerance of the
entity it produces. If it would, the entity's tolerance grows, through the
explicit, reported path described below, or the operation fails.

ADR-0028 refines the picture with further, purpose-specific controls. Each
has a unit and an owner, and none may stand in for another:

.. list-table::
   :header-rows: 1
   :widths: 30 34 36

   * - Purpose
     - Units
     - Example in the code
   * - Modelling resolution and entity tolerance
     - mm
     - ``Context::GetResolution()``, ``ToleranceOf(shape, edge)``
   * - Angular resolution
     - rad
     - ``Context::GetAngularResolution()``
   * - Parametric allowance
     - the chart's own coordinates
     - parameter brackets of root isolation (:doc:`intersection`)
   * - Solver convergence target
     - per equation
     - the relative target of adaptive arc length,
       ``Detail::ArcLengthRelativeTarget`` (:doc:`properties`)
   * - Approximation budget
     - mm
     - the budget of ``Approximate`` (:doc:`fitting`)
   * - Display tolerance
     - mm and rad
     - ``DisplayTolerance{Chord, Angle}`` (:doc:`meshing`)
   * - Numerical error bound
     - the computed quantity's
     - an edge cache's ``CacheBound``, a certified deviation
       (:doc:`trimmed`)

Loosening a display or approximation setting never changes what counts as
coincident in the model. A parameter-space allowance is never converted to
millimetres without a derivative bound over the domain, because the
conversion :math:`\lVert S_u\,\delta u + S_v\,\delta v\rVert` varies across
a surface and vanishes at a pole.

The model resolution
--------------------

Lengths are in millimetres and angles in radians. All geometry must lie in
the **size box** :math:`[-5\times10^{5}, 5\times10^{5}]^3` mm, half a
kilometre either way of the origin (``Context::GetSizeBox()``). That keeps
the error bounds of the exact predicates, and the resolution itself,
meaningful. Placing geometry inside the box is the caller's obligation: the
kernel does not check it.

The **linear resolution** :math:`\varepsilon` is :math:`10^{-5}` mm by
default, and the **angular resolution** is :math:`10^{-11}` rad. Two points
closer than :math:`\varepsilon` cannot be told apart by a single comparison.
The resolution does **not** define an equivalence relation. "Within
:math:`\varepsilon` of" is not transitive, so it answers a question about
one pair, and whether entities are merged is decided by the merge policy
below.

These are values of the operation's ``Context``, an immutable object passed
to every operation, so two models with different resolutions can be
processed side by side:

.. code-block:: cpp

   import Stereon.Core;

   const Stereon::Context ctx = Stereon::Context::Default()
                                    .WithResolution(1e-6)          // mm
                                    .WithAngularResolution(1e-12)  // rad
                                    .WithToleranceCeiling(5e-2);   // mm

A model's resolution is saved with it and restored on reading, never
replaced by the reader's default (ADR-0028; :doc:`data_exchange`).

Local tolerances and the containment invariant
----------------------------------------------

Every vertex :math:`v` carries a tolerance :math:`\tau_v` and every edge
:math:`e` a tolerance :math:`\tau_e`, both at least :math:`\varepsilon`.
Entities the kernel constructs natively are *resolution-tight*:
:math:`\tau = \varepsilon`. A tolerance is a promise about the geometry,
the **containment invariant** (:numref:`fig-tolerance-containment`):

- **Edges (rule G02).** Every trace of an edge lies within :math:`\tau_e`
  of its reference trace, over the whole range :math:`[t_0, t_1]` of the
  edge's parameter. The reference is the edge's own 3D curve :math:`C(t)`
  where it has one. Otherwise it is its primary coedge's pcurve mapped
  through that face's surface. For every other use :math:`k` of the edge,
  with pcurve :math:`c_k` on surface :math:`S_k`,

  .. math::

     \max_{t \in [t_0, t_1]} \lVert S_k(c_k(t)) - C(t) \rVert \le \tau_e .

- **Vertices (rule G01).** A vertex's tolerance sphere holds the ends of
  every edge that meets it:
  :math:`\lVert C(t_0) - P_{v_0} \rVert \le \tau_{v_0}` and
  :math:`\lVert C(t_1) - P_{v_1} \rVert \le \tau_{v_1}`.

.. _fig-tolerance-containment:

.. figure:: _static/tikz/tolerance-containment.png
   :width: 85%
   :alt: An edge's tolerance tube around its reference curve, holding a
         second trace; the vertex tolerance spheres at its ends.

   The containment invariant. The tube of radius :math:`\tau_e` about the
   edge's reference trace (blue) holds every other trace of the edge, such
   as a pcurve on an adjacent surface (red). Each vertex's sphere holds the
   ends of the edges that meet there.

The invariant is what makes a B-rep with tolerances well defined. It turns
"the edge lies on both faces" from an unmeasurable wish into a checked
inequality. The validity checker evaluates it at the Geometric level
(:doc:`topology`), comparing every trace with the reference at 17 evenly
spread parameters, ends included. Being sampled, it is a regression check
rather than a proof: a deviation between samples can pass. Certified bounds
replace the samples where interval evaluation exists. Where two entities
meet, their **pair tolerance** is the larger of the two,
:math:`\tau_{AB} = \max(\tau_A, \tau_B)`.

Growing a tolerance
-------------------

A local tolerance changes only through a named operation that states why:

.. code-block:: cpp

   STEREON_TRY_VOID(builder.SetTolerance(edge, 2e-3, "import: face gap of 1.6 um"));

Publishing the builder with ``Commit(builder, ctx)`` then enforces the
policy:

- Every change appears in the ``OperationResult`` as a ``ToleranceChange``
  (entity, old value, new value, reason). A growth also adds a
  ``ToleranceGrown`` diagnostic, so a caller sees it even with logging off.
- A tolerance above ``Context::GetToleranceCeiling()`` (0.01 mm by default)
  fails the operation with ``ToleranceExceeded``. This is the healing
  budget: a model that needs looser tolerances must say so by raising the
  ceiling.
- A tolerance below the resolution fails with ``InvalidInput``.
- Under ``Context::WithToleranceGrowthPolicy(ToleranceGrowth::Fail)`` any
  growth fails the operation. Callers that must not heal use this.

The validity checker's rule G04 reports the same bounds on any shape
(:doc:`topology`).

Coincidence is not identity
---------------------------

Proximity is not transitive. If A lies within tolerance of B, and B within
tolerance of C, A may still be far from C. A kernel that merges whatever
lies close produces a topology that depends on the order in which it
discovered the coincidences, and whose clusters can grow without bound.
Wherever an operation may identify several vertices or edges as one, as
sewing, Boolean vertex coincidence and import healing do, Stereon follows
one policy (ADR-0001), illustrated in :numref:`fig-tolerance-merge`:

1. **Order.** Candidates are processed in a fixed order: by input entity
   index, then by an operation-defined secondary key, never by completion
   order. The result is therefore the same on any number of threads
   (ADR-0007).
2. **Representative.** The first candidate opens a cluster and becomes its
   representative. Its position is kept, never averaged, so repeating an
   operation does not drift.
3. **Joining.** A later candidate joins a cluster only if it lies within the
   pair tolerance of the cluster's *representative*, not of any member, and
   the cluster's radius after joining stays within the healing budget.
   Otherwise it opens a new cluster.
4. **Result.** The merged vertex takes the representative's position, and a
   tolerance equal to the largest of the members' tolerances and the
   cluster's radius, at least :math:`\varepsilon`. Any growth is logged.
5. **Edges** merge by the same procedure. The pair test is the containment
   of one edge's curve in the other's tube, at both ends and at the
   certified points of maximum deviation.

.. _fig-tolerance-merge:

.. figure:: _static/tikz/tolerance-merge.png
   :width: 92%
   :alt: Four vertices in a row, 0.8 tau apart, grouped into two clusters
         around the representatives A and C.

   The merge policy on a chain of vertices :math:`0.8\tau` apart. :math:`B`
   joins :math:`A`'s cluster. :math:`C` lies within :math:`\tau` of
   :math:`B` but :math:`1.6\tau` from the representative :math:`A`, so it
   opens a cluster of its own, which :math:`D` joins. Single-linkage
   clustering would have merged all four into one cluster
   :math:`2.4\tau` wide.

The policy is deliberately order-dependent in the rare non-transitive case,
but deterministically so, and every cluster stays within the healing
budget. An order-independent clustering cannot bound the cluster radius,
which is the property that matters. Single-linkage clustering ("merge
anything that touches anything already merged") is not used anywhere in the
kernel.

Tolerances in operations
------------------------

**Intersection.**
   ``Intersect(a, boxA, b, boxB, tolerance, ctx)`` takes the pair tolerance
   of the faces being intersected (:doc:`intersection`). It is used for one
   decision only: whether two *exactly parallel* supports are coincident.
   Two parallel planes whose separation is at most the pair tolerance are
   reported as an overlap, with ``CoincidenceEvidence::Tolerant`` and the
   separation in ``Distance``. Exactly coincident supports are reported
   with ``CoincidenceEvidence::Exact``. Every other decision, such as
   whether two surfaces cross, touch or miss, is made exactly on the input
   doubles by the predicates of :doc:`predicates`, without any tolerance.

**Booleans.**
   When a Boolean identifies a vertex of one body with a vertex of the
   other, the two must lie within their pair tolerance, or the operation
   fails with ``IntersectionFailed``. The identified vertex takes the pair
   tolerance, a growth that is logged with the reason "identified with the
   other body's vertex" (:doc:`imprinting`).

**No hidden epsilons.**
   Code in ``topo``, ``boolean`` and ``mesh`` never compares with a literal
   constant. It asks for ``ToleranceOf(shape, entity)`` or
   ``ctx.GetResolution()``, and ``tools/check_epsilon.py`` fails the build
   on a violation.

Implementation
--------------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Code
     - Role
   * - ``src/core/Context.cppm``
     - Resolution, angular resolution, ceiling, growth policy, size box.
   * - ``src/topo/Builder.cppm`` (``SetTolerance``)
     - The one way a local tolerance changes.
   * - ``src/topo/Builder.cpp`` (``Commit``)
     - Reports changes; enforces the ceiling, the resolution floor and the
       growth policy.
   * - ``src/topo/ValidateGeometry.cpp``
     - Rules G01, G02 and G04.
   * - ``src/boolean/Assemble.cpp``
     - Vertex identification at the pair tolerance.
   * - ``tools/check_epsilon.py``
     - The lint against literal tolerances.

Notes and further reading
-------------------------

Tolerant modelling with local tolerances on edges and vertices was
introduced in Parasolid and is now used, in some form, by every major
kernel. Mäntylä's text :cite:p:`mantyla1988solid` treats the
representation the tolerances are attached to. Patrikalakis and Maekawa
:cite:p:`patrikalakis2002shape` discuss the robustness of interrogation
algorithms under tolerances. The separation of exact combinatorial
decisions from approximate geometry follows the exact-geometric-computation
paradigm of Shewchuk :cite:p:`shewchuk1997adaptive` and others, which
:doc:`predicates` develops.

Related design decisions
------------------------

- `ADR-0001: Tolerance model <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0001-tolerance-model.md>`_
- `ADR-0002: Numeric robustness <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0002-numeric-robustness.md>`_
- `ADR-0018: Validity checking <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0018-validity-checking.md>`_
- `ADR-0028: Purpose-specific numerical tolerances <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0028-purpose-specific-numerical-tolerances.md>`_
- `ADR-0029: Intersection and classification outcomes <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0029-intersection-and-classification-outcomes.md>`_
