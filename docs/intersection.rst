.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Intersections
=============

Curve-curve, curve-surface and surface-surface intersection, including the marching and subdivision strategy.

Outcomes
--------

``Intersect`` returns a complete result or the reason there is none
(ADR-0029). Success means the completeness and accuracy obligations of
ADR-0010 are met over the requested domains. A result with no component is
therefore a **proved** separation, never "no seeds found". Anything less is a
``KernelError``: ``IntersectionFailed``, ``CertificationFailed`` or
``NotSupported``, never a partial success.

.. code-block:: cpp

   // Over the faces' parameter boxes, as a face-level operation asks.
   const auto result = Stereon::Intersect(surfaceA, boxA, surfaceB, boxB, ctx);
   // Over the whole domains.
   const auto whole = Stereon::Intersect(surfaceA, surfaceB, ctx);

A ``SurfaceSurfaceIntersection`` holds points and branches, which may
coexist:

- ``IntersectionPoint``: an isolated contact, with its parameters on both
  surfaces.
- ``IntersectionBranch``: a curve with its pcurves on both surfaces on one
  shared parameter (ADR-0009) and its parameter range inside the requested
  domains. It also carries how the surfaces meet (``Crossing``,
  ``Tangential``, or ``Overlap`` for an overlap's boundary), the
  ``SideTransition`` across it, and a certified deviation bound. An infinite
  bound means uncertified.

An ``OverlapRegion`` describes where the surfaces coincide: its domain on
each surface, the ``ParameterMap`` from one surface's parameters to the
other's (applied by ``Mapped``), whether the normals agree, and its
boundary branches as indices into the result's branches.

Contact and side transition are separate facts. A transverse crossing implies
``OppositeSide``, while a tangency alone implies neither side relation. For
quadrics a touch is proved ``SameSide`` from the closed form; otherwise the
transition is ``Unresolved``, and classification never propagates across it.

A coincidence is ``Exact`` when the supports coincide in the input doubles.
It is ``Tolerant`` when they are exactly parallel and within the pair
tolerance: declared, with its distance, and with no surface moved. Coincident
supports do not make overlapping faces, because the faces' trims decide that.

Exact analytic geometry stays exact: a plane–plane branch is a ``Line``, not
a B-spline.

Analytic handlers
-----------------

Plane–plane, plane–cylinder and cylinder–cylinder with parallel axes have
closed-form handlers (ADR-0010's analytic table, RS-0001 Tier A). Every
configuration decision is the sign of a polynomial in the input doubles
and the pair tolerance, estimated in interval arithmetic and settled by
expansion arithmetic from ``Stereon.Robust`` when the estimate straddles
zero (RS-0001 G5): whether two normals or axes are exactly parallel, a
height against a radius, the distance between parallel axes against the
sum and difference of the radii. Nothing depends on round-off, evaluation
order or thread count.

.. code-block:: cpp

   // With the pair's tolerance, the larger of the two faces' (ADR-0001).
   const auto result = Stereon::Intersect(plane, boxA, cylinder, boxB, tolerance, ctx);

The configurations:

- **Plane–plane:** a ``Line`` parametrised by arc length, with a ``Line2d``
  on each plane on the same parameter; exactly parallel planes are a proved
  separation beyond the pair tolerance and an overlap within it, ``Exact``
  when the height is exactly zero and ``Tolerant`` otherwise, with the map
  between the frames.
- **Plane–cylinder:** every branch is parametrised by the cylinder's angle.
  A plane across the axis cuts a ``Circle`` (a ``Line2d`` on the cylinder,
  a ``Circle2d`` on the plane); a plane along the axis cuts two rulings, one
  tangent ruling (``Tangential``, proved ``SameSide``) or nothing, decided by
  the axis's height over the plane against the radius; any other plane cuts
  an ``Ellipse``, exact as an ``Ellipse2d`` on the plane. Its pcurve on the
  cylinder is the sinusoid :math:`(u, v_c + B \cos t)`, which no exact 2D
  type represents: it is a clamped cubic B-spline through samples dense
  enough for the interpolation error bound to stay within half the
  resolution, and the certificate below proves what it achieved.
- **Cylinder–cylinder, parallel axes:** two rulings, an external or
  internal tangent ruling (``SameSide``), a proved separation apart or
  nested, or an overlap, exact or declared: the surfaces' Hausdorff
  distance, the axes' distance plus the radii's difference, within the
  pair tolerance. Skew and intersecting axes, the quartic, report
  ``NotImplemented`` until their handler (M3b).

A branch's ``Range`` is clipped to the requested domains with outward
rounding on the nominal curve: exactly for lines and for the angular
parameter of a cylinder's circles and ellipses, and otherwise only where
the branch's parameter box leaves the domain entirely. A face-level
operation restricts a branch to its loops itself; a box should include the
face's tolerance for the clip to keep every point within it. A declared
coincidence is logged at ``Info`` with its distance and the tolerance.

Certificates
------------

Every branch over a bounded range is *certified* (ADR-0009 section 4): its
``Bound`` is a proved bound on its deviation from both surfaces, through
their pcurves, over the whole range, within the Context's resolution. A
branch over an unbounded range stays uncertified until a face's loops
bound it. A branch that cannot be certified fails the intersection with
``CertificationFailed``, never a partial success (ADR-0029).

The proof is ``CertifyDeviation(curve, surface, pcurve, range, target)``
in ``Stereon.Surface``, which the topology's cached edge curves use too.
It rests on ``IntervalScalar``, a ``Scalar`` (ADR-0021) whose arithmetic
is core's outward-rounded ``Interval``: the geometry templates, given it,
evaluate over parameter intervals and return rigorous enclosures, with the
sine and cosine enclosing their extrema inside a span. The range is
bisected; on each span the deviation ``f(t) = c(t) - S(p(t))`` is enclosed
by the mean-value form ``f(m) + f'(span) (t - m)`` at the midpoint, and a
span is accepted once its enclosure is within the target. The largest
accepted enclosure is the bound, with the number of spans it took. A
B-spline's enclosures come from its poles: the hull of those active over
the span, narrowed within one knot span by the same mean-value form from
its derivative curves; ``Eval``'s rounded basis functions are covered by a
margin.

The search cannot be fooled and does not run away: a span whose midpoint
provably deviates by more than the target (the least magnitude of its
enclosure) ends it at once with that witness; a span that can no longer be
split, or a budget of 2^18 spans, ends it with the bound reached, which is
where a target below the rounding of the enclosures leads. Lines, circles,
ellipses and non-rational B-splines are covered as curves and pcurves, and
planes and cylinders as surfaces; anything else is ``NotSupported``, and a
range outside a B-spline's domain is ``InvalidInput``.

Every other pair reports ``NotImplemented`` until Phase 2.

.. todo::

   The generic path (subdivision, marching) when Phase 2 builds it.

Related design decisions
------------------------


