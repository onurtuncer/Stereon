.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Intersections
=============

Curve-curve, curve-surface and surface-surface intersection, including the marching and subdivision strategy; plane curve intersection for the face clipping of the Booleans.

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
  resolution, and the certificate below proves what it achieved. The fit
  is limited to 2^16 sample intervals; an excessive or non-finite count
  returns ``CertificationFailed`` before integer conversion or allocation.
- **Cylinder–cylinder, parallel axes:** two rulings, an external or
  internal tangent ruling (``SameSide``), a proved separation apart or
  nested, or an overlap, exact or declared: the surfaces' Hausdorff
  distance, the axes' distance plus the radii's difference, within the
  pair tolerance.
- **Cylinder–cylinder, skew or intersecting axes:** the quartic, as
  certified B-spline branches; see below.

Cylinders with skew axes
~~~~~~~~~~~~~~~~~~~~~~~~

Whether two axes are parallel is decided exactly, so nearly parallel
cylinders are the skew case (RS-0001 Tier B; ADR-0010's decisions on
acceptance). Along a ruling of the first cylinder at angle :math:`u`, the
second cylinder's equation is a quadratic in the height :math:`v`, so the
curve is two sheets :math:`v_\pm(u)` over the angles where the discriminant
is not negative. With :math:`\sigma` the angle from the axes' common normal
and :math:`\delta` the distance between the axes, the discriminant is
:math:`r_B^2 - (r_A \cos\sigma - \delta)^2`, and the configuration follows
from :math:`\delta` against :math:`r_A + r_B` and :math:`|r_A - r_B|`,
decided exactly on the input doubles (a polynomial of degree ten, on the
heap-allocated expansions of :doc:`predicates` when the interval filter
cannot settle it):

- :math:`\delta > r_A + r_B`: a proved separation; :math:`\delta = r_A +
  r_B`: one ``IntersectionPoint``, certified on both cylinders by interval
  evaluation.
- :math:`|r_A - r_B| < \delta < r_A + r_B`: one closed loop.
- :math:`\delta < |r_A - r_B|`: two loops, where the smaller cylinder
  enters and leaves the larger.
- :math:`\delta = |r_A - r_B|`: internally tangent at a *node*, a double
  root of the discriminant, where the loops meet; :math:`\delta = 0` with
  equal radii: two ellipses through two nodes. Both are logged at ``Info``.

Every branch is a ``Crossing`` with ``OppositeSide`` evidence; branches end
at nodes, since the surfaces are tangent there, and at the chart seams of
both cylinders, so that each pcurve stays within :math:`[0, 2\pi]`. A branch
is parametrised so that it is smooth: by the first cylinder's angle where
the discriminant stays positive, and otherwise so that the signed square
root of the discriminant is smooth across the arc's ends, the folds where
the curve's tangent runs along a ruling. Its 3D curve and both pcurves are
cubic B-splines through samples of the closed form, the samples doubled
from 32 until the midpoints agree with the closed form and the certificate
below holds, within a budget of 2^14; beyond it the intersection fails with
``CertificationFailed``, as it does at once for a branch that leaves the
size box of ADR-0001, where the rounding of a B-spline's poles alone
exceeds the resolution.

A branch is first clipped to the requested angular and height restrictions
on both cylinders by interval evaluation of the closed form over spans of
its parameter, bisected until a span is proved outside or inside. Half-infinite
height intervals also restrict the result. Angular sectors may cross a seam
or be shifted by whole turns. A bounded angular interval narrower than one
turn with either endpoint of magnitude at least 2^52 radians is refused with
``InvalidInput``: its whole-period shifts cannot be resolved safely. Two cylinders whose axes are nearly
parallel meet in loops of the order of :math:`r / \sin\gamma` along the
axes, which a face's bounded heights cut down to the two nearly straight
pieces inside them.

Coincident cylinders intersect the requested height domains through the
parameter map and intersect angular domains modulo the period. An overlap
crossing a chart seam is returned as separate regions. Disjoint requested
domains return no overlap, even when the supporting cylinders coincide.
``FullDomain`` describes coverage of both requested domains, including
bounded domains, rather than equality with the unbounded surface domains.

A branch's ``Range`` is clipped to the requested domains with outward
rounding on the nominal curve: exactly for lines and for the angular
parameter of a cylinder's circles and ellipses, and otherwise only where
the branch's parameter box leaves the domain entirely. A face-level
operation restricts a branch to its loops itself; a box should include the
face's tolerance for the clip to keep every point within it. A declared
coincidence is logged at ``Info`` with its distance and the tolerance.

Space curves
------------

The first v0.2 3D curve handler is ``Line`` against ``Line`` through the
existing ``Intersect(Curve, domain, Curve, domain, ctx)`` overload. Directions
need not be unit vectors; returned parameters retain each input's speed and
orientation. Other 3D curve pairs still return ``NotImplemented``.

Cross products and the scalar triple product decide parallelism, collinearity
and coplanarity exactly on the input doubles. Skew and parallel distinct
supports return an empty result even when closer than the model resolution.
There is no tolerant snapping. A nonparallel coplanar pair returns one point
only when its exact rational parameters belong to both requested closed
domains. Rounded parameters are never used to decide domain membership.

Collinear requested domains return an overlap, a single endpoint contact,
or proved separation. Endpoint-only contact is a ``CurveCurvePoint``, not a
zero-length overlap. Overlap intervals enclose the exact endpoint parameters
with outward rounding, restricted to the requested domains; ``SameDirection``
records orientation. Whole lines and half-infinite domains are supported.

The handler checks finite origins, finite nonzero directions and the numerical
domain of its degree-four expansion comparisons, including finite parameter
bounds. Inputs requiring rescaling or exceeding that domain return
``NotSupported`` in this initial implementation. Infinite-only domains are
``InvalidInput``. Point coordinates and rounded parameters must satisfy a
conservative interval error bound within the Context resolution; otherwise the
operation returns ``CertificationFailed``. This conservative check may reject
large-coordinate contacts even when a particular rounded point is exact.

Curves against surfaces
-----------------------

``Intersect(curve, domain, surface, box, ctx)`` meets a curve and a
surface over the requested domains (ADR-0029). A
``CurveSurfaceIntersection`` holds ``CurveSurfacePoint`` entries, in order
of the curve's parameter, each with its parameter on the curve, its
parameters on the surface, its position, the ``ContactType`` there
(``Crossing``, the curve passing through the surface, or ``Tangential``,
touching it and staying on one side) and whether the count and the contact
were *decided* exactly on the input doubles; and ``CurveSurfaceOverlap``
stretches where the curve lies in the surface, each the parameter interval
on the curve and the pcurve on the surface on the curve's own parameter
(ADR-0009). A point is kept when its parameters lie in the requested
domains, a periodic surface's angle taken modulo its period and placed in
the box. ``Intersect(curve, surface, ctx)`` takes the whole domains. How a
result moves with design parameters, and when it has no derivative, is in
:doc:`sensitivity` (*Products of intersections*): the ``Contact`` and
``Decided`` of a point decide its derivative's status.

The first handlers are a ``Line`` against a ``Plane``, a ``Sphere`` and a
``Cylinder``: ADR-0010's analytic table carried to the line–quadric pairs.
Every configuration decision is the sign of a polynomial in the input
doubles, written once over a scalar field, read as an enclosure and
settled exactly when the enclosure straddles zero (RS-0001 G5): whether
the line runs along the plane, :math:`(X \times Y) \cdot d = 0`, or along
the cylinder's axis, :math:`\lVert d \times (X \times Y) \rVert^2 = 0`;
whether it then lies in the surface or apart from it, however close; and
the discriminant of the quadratic the line cuts out of the sphere or the
cylinder, which settles whether there are two crossings, one touch or
none. The crossing of a plane has a rational parameter, and whether it
and the point's parameters on the plane lie in the requested domains is
decided exactly on that rational, as the line–line handler decides; the
roots of a quadric are irrational and are placed by their rounded values,
so a crossing within rounding of a domain's end may fall on either side
of it. A touch of a quadric is a ``Tangential`` point: the line stays
outside. A root that double precision cannot place, a line so nearly along
the axis that the quadratic's leading coefficient rounds to zero while it
is exactly positive, fails with ``IntersectionFailed`` rather than
dropping the point. ``tests/occt/IntersectionDifferentialTests.cpp``
compares the three handlers with OCCT's ``IntAna_IntConicQuad`` on
generic lines.

A ``BSplineCurve``, polynomial or rational, has handlers against the same
three surfaces. The surface's implicit equation composed with the curve is,
on each Bézier segment and in homogeneous coordinates, a polynomial
:math:`g` of degree :math:`p` against the plane and :math:`2p` against a
quadric, of the sign of the equation since every weight is positive. Its
Bernstein coefficients are enclosed in outward-rounded interval arithmetic,
from segments that are themselves enclosed: they come from Boehm's knot
insertion done in intervals, so that what is proved holds for the curve as
given, not for a rounded copy. The segment is subdivided at halves, and
every piece ends in one of four ways:

- **apart**, when :math:`g`'s coefficients all have one strict sign, a
  proof that there is no root (the Bernstein convex hull);
- **a crossing**, when the derivative's coefficients have one strict sign,
  so that there is at most one root, a simple one: ``Decided`` when
  :math:`g` has strictly opposite signs at the piece's ends, the root then
  bisected to rounding, undecided when an end value straddles zero;
- **a touch**, ``Tangential`` and undecided, when :math:`g` is proved
  strictly convex or concave and the piece holding its extremum is within
  rounding of it: a tangency, a near miss or two crossings within rounding.
  A touch is placed only to the square root of rounding, as any double
  root is;
- **unresolved**, a contact of higher order, a curve along a quadric, or a
  search past its budget: the intersection fails with
  ``IntersectionFailed``.

So an empty result is a proof that the curve misses the surface over the
requested domain, and nothing unresolved becomes one. A curve 3·2⁻⁴⁰ off a
plane is apart; a root that lands on a piece's end, a knot or a dyadic
parameter, is seen from both sides and reported once, undecided; a root
whose bracket reaches past the requested domain is kept, clamped into it,
rather than dropped on its rounded parameter. A curve whose every pole lies
exactly in the plane is an overlap, its pcurve the spline of the projected
poles. A periodic curve is intersected within its own domain only, a
curve in a plane that leaves the requested box is not clipped yet, and a
point's surface parameters are those of its rounded position, kept or
dropped against the box as the line handlers do. Every other pair reports
``NotImplemented``; the generic path for curves against B-spline surfaces
waits for Phase 2 with the surface–surface one.

Plane curves
------------

``Intersect(curve2d, domain, curve2d, domain, ctx)`` meets two plane
curves: the pcurves of a face's loop against the pcurve of a branch in the
face's parameter space, which is how a face-level operation restricts a
branch to its loops (ADR-0009; the Boolean prototype's M6). A
``PlanarIntersection`` holds ``PlanarPoint`` entries, in order of the first
curve's parameter, each with its parameter on both curves, its position,
the ``ContactType`` there (``Crossing`` or ``Tangential``) and whether the
count and the contact were *decided* exactly on the input doubles; and
``PlanarOverlap`` stretches where coincident curves share their domains,
with the parameter correspondence. A point is kept when both parameters
lie in the requested domains, a periodic curve's taken modulo its period
and placed there.

Every ordered pair of the six built-in plane curves, ``Line2d``,
``Circle2d``, ``Ellipse2d``, ``Hyperbola2d``, ``Parabola2d`` and
``BSplineCurve2d`` (polynomial or rational), has a handler; the table of
``src/intersect/Dispatch.cppm`` (``PlanarHandlerFor``) chooses it at
compile time, and a pair given the other way round is the same
computation with its result's sides exchanged (``SwappedPlanar``).

.. list-table:: Plane curve handlers
   :header-rows: 1
   :widths: 28 32 40

   * - Pair
     - Handler (file)
     - Method
   * - line–line, line–circle, line–ellipse, circle–circle
     - ``LineLine2d``, ``LineCircle2d``, ``LineEllipse2d``,
       ``CircleCircle2d`` (``PlanarLines.cpp``)
     - closed form; every decision exact
   * - line–hyperbola, line–parabola
     - ``LineHyperbola2d``, ``LineParabola2d`` (``PlanarLines.cpp``)
     - the conic's chart in the line's equation; every decision exact
   * - two conics, any of circle, ellipse, hyperbola and parabola, but two
       circles
     - ``ConicConic2d`` (``PlanarConics.cpp``)
     - Algorithm C, the exact chart quartic
   * - B-spline–line, B-spline–conic
     - ``SplineCurve2d``, ``CurveSpline2d`` (``PlanarSpline.cpp``)
     - Algorithm S, Bernstein subdivision of the implicit equation
   * - B-spline–B-spline
     - ``SplineSpline2d`` (``PlanarSplinePair.cpp``)
     - Algorithm K, subdivision with Krawczyk's test

What ``Decided`` means differs by row, and the contract is the same for
all: a successful result with no component is a proved separation, a
meeting that cannot be proved a single crossing or isolated as a touch is
never dropped, and anything left unresolved is ``IntersectionFailed``
(ADR-0029). ``Contact`` records whether the curves change sides there: a
tangency of odd order, which crosses, is a ``Crossing``; one of even
order, which does not, is ``Tangential``.

Lines against lines and conics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lines, circles and ellipses against each other decide every configuration
as the sign of a polynomial in the input doubles (RS-0001 G5): the
determinant of two directions, whether a parallel line passes through a
point, the discriminant of the quadratic a line cuts out of a circle or an
ellipse, and the distance between two centres against the sum and the
difference of the radii. The polynomials are written once over a scalar
field and read as an enclosure first and exactly, on heap expansions, only
when the enclosure straddles zero (``SignOfFormula`` in
``src/intersect/Planar.hpp``).

A hyperbola :math:`C(t) = c + a\cosh t\,X + b\sinh t\,Y` is taken in its
rational chart :math:`u = \tanh(t/2) \in (-1, 1)`, in which
:math:`(1-u^2)\,C = (1-u^2)\,c + a(1+u^2)\,X + 2bu\,Y` covers the branch
exactly once; the chart's ends :math:`u = \pm 1` are its points at
infinity. With :math:`n` across the line through :math:`O`,
:math:`k_0 = n\cdot(c-O)`, :math:`k_x = a\,n\cdot X` and
:math:`k_y = b\,n\cdot Y`, the line's equation along the chart is the
quadratic

.. math::

   f(u) = (k_x - k_0)\,u^2 + 2k_y\,u + (k_0 + k_x) = \alpha u^2 + 2\beta u + \gamma ,

whose roots in the open interval are the meetings. The discriminant
:math:`\beta^2 - \alpha\gamma` decides two crossings, one tangency or none;
where each root lies against :math:`\pm 1` is decided from the exact signs
of :math:`f(\pm 1)` against :math:`\alpha` (opposite between the roots,
alike outside) and of :math:`\alpha e + \beta` (which side of the vertex
:math:`e` is on, ``SideOf``). A line parallel to an asymptote makes
:math:`\alpha = 0` exactly: one crossing when
:math:`\gamma^2 < 4\beta^2`, none otherwise, which is how such a line
misses the branch or meets it once. A parabola
:math:`C(t) = V + t^2/(4f)\,X + t\,Y` gives the quadratic
:math:`(n\cdot X)\,t^2 + 2(2f\,n\cdot Y)\,t + 4f\,n\cdot(V-O)` over every
real :math:`t`, linear for a line along the axis. All of these decisions
are exact; the roots themselves are irrational and are placed by their
rounded values, as a line's roots on a circle are.

Algorithm C: two conics by the exact chart quartic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ConicConic2d`` (``src/intersect/PlanarConics.cpp``) takes one conic, A,
in a rational chart and the other, B, by its implicit equation in its frame
:cite:p:`patrikalakis2002shape`. A's point is the homogeneous
:math:`W(u)\,c_A + \alpha(u)\,X_A + \beta(u)\,Y_A` for :math:`u \in [-1, 1]`:

.. list-table:: Rational charts (``ChartOf``)
   :header-rows: 1
   :widths: 20 40 40

   * - A
     - chart
     - :math:`(W, \alpha, \beta)`
   * - ellipse, circle
     - :math:`u = \tan(s/2)` about :math:`s = 0`, and about :math:`s = \pi`
       (:math:`h = -1`)
     - :math:`(1+u^2,\; h a(1-u^2),\; 2hbu)`
   * - hyperbola
     - :math:`u = \tanh(t/2)`; :math:`u = \pm 1` at infinity
     - :math:`(1-u^2,\; a(1+u^2),\; 2bu)`
   * - parabola
     - :math:`t = Lu`, and :math:`t = L/u` (:math:`u = 0` at infinity),
       :math:`L` the power of two at or below :math:`4f`
     - :math:`(4f,\; L^2u^2,\; 4fLu)` and :math:`(4fu^2,\; L^2,\; 4fLu)`

With :math:`(\xi, \eta, W)` the homogeneous coordinates of A's point in
B's frame, :math:`\xi = (H - Wc_B)\cdot X_B` and
:math:`\eta = (H - Wc_B)\cdot Y_B`, B's equation along the chart is

.. math::

   Q(u) =
   \begin{cases}
   b^2\xi^2 + a^2\eta^2 - a^2b^2W^2 & \text{B an ellipse or a circle,}\\
   b^2\xi^2 - a^2\eta^2 - a^2b^2W^2 & \text{B a hyperbola,}\\
   \eta^2 - 4f\,\xi W & \text{B a parabola,}
   \end{cases}

a polynomial of degree four whose coefficients are polynomials in the input
doubles, formed exactly on heap expansions (``Composite``). Since a
hyperbola's equation holds both its branches, the open conic is the one
put in the chart when there is a choice (``ChartRank``: a hyperbola before
a parabola before an ellipse; between two of one rank, the first in a fixed
order of their data, so that the pair exchanged is the same computation,
``ChartedSecond``), and a root on B's other branch, whose
abscissa is at most :math:`-a`, :math:`2a` from the branch, is dropped by
the sign of its rounded abscissa. The steps:

1. **Coincidence.** If :math:`Q` vanishes identically, A lies on B's conic:
   two ellipses, hyperbolas or parabolas coincide (or A is the other branch
   of B's hyperbola, a proved separation). B's parameter is then
   :math:`\sigma t + \tau` for A's :math:`t`, :math:`\sigma = \pm 1` and
   :math:`\tau` read from A's point and tangent at 0, and the overlaps are
   the requested domains' common stretches under that map, modulo the
   period for ellipses (``AffineOverlap``, ``PeriodicOverlaps``).
2. **Points at infinity.** :math:`Q(\pm 1) = 0` on a hyperbola's chart, or
   :math:`Q(0) = 0` on a parabola's far chart, is a meeting at infinity (a
   curve along an asymptote's direction, or an axis's); it is divided out
   exactly by synthetic division by :math:`u - r`, which for a double
   :math:`r` needs only exact products and sums (``Deflated``).
3. **Isolation.** The chart is subdivided at halves
   :cite:p:`lane1981bounds`; on each piece :math:`[u_0, u_0 + w]` the
   Bernstein coefficients of :math:`Q` are formed from :math:`Q` itself:
   the Taylor coefficients :math:`c_k` at :math:`u_0` are exact, as are the
   shifted :math:`c_k w^k`, and only the change of basis
   :math:`b_i = \sum_{j\le i} \binom{i}{j}\binom{n}{j}^{-1} c_j w^j` is
   outward-rounded (``BernsteinOn``). The enclosures stay as tight as the
   rounding of their own size however small the piece, which the Bernstein
   basis's conditioning makes the right basis for
   :cite:p:`farouki1987numerical`; halving enclosed coefficients instead
   would leave a floor of the whole chart's rounding, below which a near
   miss or a contact of order four could not be resolved. Each piece is
   classified as Algorithm S classifies it, below, and a proved crossing's
   bracket, whose ends have strictly opposite signs, is narrowed by
   bisection on the doubles with :math:`Q`'s exact sign to two adjacent
   doubles (``Tighten``).
4. **Exact settlement.** A root left undecided, or not isolated at all, is
   settled when it lies at a double :math:`u^*` (its placement, its
   bracket's ends or middle, or 0): :math:`Q(u^*) = 0` is decided on
   expansions, the multiplicity :math:`m` is the index of the first
   non-zero exact Taylor coefficient at :math:`u^*`, and
   :math:`Q/(u-u^*)^m = \sum_{k\ge m} c_k (u-u^*)^{k-m}` is enclosed over
   the bracket by Horner's rule in intervals. When the enclosure excludes
   zero, :math:`u^*` is the bracket's only root, of multiplicity :math:`m`,
   and the point is ``Decided``: a ``Crossing`` for odd :math:`m`, a
   ``Tangential`` touch for even (``MultiplicityAt``, ``Settle``). The
   touch of a circle and an ellipse at a common vertex, and the contact of
   order four of an ellipse with its osculating circle there, are settled
   so.
5. **Charts to curve.** A root's parameter and bracket are mapped through
   the chart (:math:`s = 2\arctan u`, :math:`t = 2\,\mathrm{artanh}\,u`,
   :math:`t = Lu` or :math:`L/u`), and roots seen by two charts at their
   common end are merged (``MergeCharts``).

A root that is not settled is reported with ``Decided`` false: a touch, a
near miss or two crossings within rounding of each other, bracketed to
rounding; a root that is not isolated, a contact of higher order away from
any double, or a root within rounding of a point at infinity (beyond the
range of the parameter in doubles) fails with ``IntersectionFailed``.

Algorithm S: a B-spline against a line or a conic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SplineCurve2d`` (``src/intersect/PlanarSpline.cpp``) composes the other
curve's implicit equation with the spline. On each Bézier segment, in
homogeneous coordinates :math:`(X, Y, W) = (wP, w)`, it is a polynomial
:math:`g` of degree :math:`p` against a line, :math:`n\cdot(X - WO)`, and
of degree :math:`2p` against a conic, :math:`|X - Wc|^2 - r^2W^2` for a
circle and :math:`Q` above with :math:`(\xi, \eta, W)` the segment's
homogeneous coordinates in the conic's frame otherwise; :math:`g` has the
sign of the equation since every weight is positive, so rational splines
need nothing more. The segments are those of ``src/intersect/Bernstein.hpp``
(``SegmentsOf``), from Boehm's knot insertion done in interval arithmetic
:cite:p:`piegl1997nurbs`, so that every proof holds for the spline as given
rather than for a rounded copy, and :math:`g`'s Bernstein coefficients are
enclosed from them. Every piece of the subdivision ends in one of four
ways, the classification shared with the curve–surface handlers
(``Isolate``, ``Monotone``, ``Convex``):

.. graphviz::
   :caption: Classifying a piece of a polynomial in Bernstein form.

   digraph classify {
     rankdir=TB;
     node [shape=box, fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9];
     piece [label="piece [s0, s1]: coefficients b_i enclosed"];
     apart [label="apart: no root (convex hull)"];
     mono [label="monotone: at most one root"];
     crossing [label="ends of strictly opposite sign:\nproved crossing, bisected (Decided)"];
     edge_root [label="an end straddles zero:\ncrossing within rounding of it (undecided)"];
     convex [label="convex or concave: one extremum"];
     touch [label="extremum within rounding:\ntouch (Tangential, undecided)"];
     split [label="halve"];
     fail [label="at MaxDepth or past the budget:\nunresolved"];
     piece -> apart [label="all b_i strictly of one sign"];
     piece -> mono [label="all differences strictly of one sign"];
     mono -> crossing; mono -> edge_root;
     piece -> convex [label="all second differences\nstrictly of one sign"];
     convex -> touch; convex -> mono [label="sub-piece"];
     piece -> split [label="otherwise"];
     split -> piece;
     split -> fail;
   }

Unlike the curve–surface handlers, which fail on any root they cannot
decide, the plane curve handlers report an undecided root with ``Decided``
false, as the conic handlers do, and fail only on an unresolved one: a
contact of higher order, or the spline lying along the other curve (a
rational spline along the circle it draws fails so, never coming back
empty). A hyperbola's other branch is dropped as in Algorithm C. A spline
whose every pole lies exactly on a line lies along it: its overlap is
reported when the spline is clamped, is requested over its whole domain,
runs monotonically along the line
(exact signs of :math:`d\cdot(P_{i+1}-P_i)`; variation diminishing with
positive weights) and lies inside the line's requested domain, a proved
separation when it lies outside, and ``NotImplemented`` when the domains
cut it, which is not clipped yet. A periodic spline is intersected within
its own domain only (``NotSupported`` beyond it).

Algorithm K: two B-splines by subdivision with Krawczyk's test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SplineSpline2d`` (``src/intersect/PlanarSplinePair.cpp``) solves
:math:`F(s, t) = a(s) - b(t) = 0` over pairs of pieces of the two splines'
enclosed Bézier segments, halving both pieces of a pair until it is

- **apart**: the boxes of the pieces' poles are disjoint, or the hulls of
  both pieces' distances from either piece's chord are, a fat line in the
  sense of Sederberg and Nishita :cite:p:`sederberg1990curve`; with
  positive weights each piece lies in the hull of its poles, so either is a
  proof (``Apart``);
- **settled by Krawczyk's operator** :cite:p:`krawczyk1969newton` over the
  pair's parameter box :math:`X`:

  .. math::

     K(X) = m - Y F(m) + \bigl(I - Y J(X)\bigr)(X - m),

  with :math:`m` the box's middle, :math:`J(X) = [a'(S), -b'(T)]` the
  Jacobian enclosed over :math:`X` (the derivative
  :math:`(X'W - XW')/W^2` of each rational piece by interval de Casteljau
  on its poles and its hodograph's, in the form
  :math:`b_i + s(b_{i+1} - b_i)` whose enclosure shrinks with the poles'
  spread) and :math:`Y` the inverse of :math:`J`'s midpoint. If
  :math:`K(X) \cap X = \emptyset` there is no root in :math:`X`; if
  :math:`K(X)` lies in :math:`X`'s interior there is exactly one, and
  :math:`J` is regular there, so it is a transverse crossing
  :cite:p:`moore2009interval,neumaier1990interval`. The box is then
  shrunk to rounding by iterating :math:`X \leftarrow K(X) \cap X`, which
  keeps the root (``Krawczyk``, ``Refined``).

Krawczyk's box is the pair's own box inflated by an eighth of its width on
every side, so that a root on the end of a piece, where the halving put
it, is proved from either side. The inflation reaches past a segment's ends
into the segment's polynomial continued: a root whose bracket lies beyond
the segment is not the segment's and is dropped (the neighbouring segment
finds it, or it is no point of the curve); one whose bracket straddles a
segment's end is reported with ``Decided`` false; roots proved from
overlapping boxes are one. Every proved crossing is ``Decided``. A pair
that is neither apart nor settled when its pieces can no longer be halved,
or a search past :math:`2^{16}` pairs, fails with ``IntersectionFailed``:
a tangency of two B-splines (where :math:`J` is singular and no box ever
satisfies the test), an overlap, or a crossing too nearly tangential to
separate in doubles. Two identical splines coincide, and their overlap is
the requested domains' common stretch on the identity map; their crossings
with themselves are not reported then.

What is decided
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Pairs
     - ``Decided`` true when
   * - a line against a line or a conic; two circles
     - always: count and contact are exact decisions; irrational roots are
       placed to rounding
   * - two other conics
     - a root proved a single simple crossing, or settled exactly at a
       double with its multiplicity
   * - a B-spline against a line or a conic
     - a root proved a single simple crossing
   * - two B-splines
     - a crossing proved single by Krawczyk's test, its bracket inside both
       segments

A point whose parameter on a conic lies within rounding of a requested
domain's end may fall on either side of it; a B-spline root whose bracket
reaches past the requested domain is kept, clamped into it, undecided. The
tests are ``tests/intersect/PlanarTests.cpp`` and
``tests/intersect/PlanarPairTests.cpp``: exact configurations and their
neighbours a rounding or :math:`3\cdot 2^{-40}` either side, tangencies of
order two and four, coincidences, rational arcs, exchange and reversal of
the curves, and properties on random geometry (every sign change of a
conic's equation along a random cubic has a point; random conic pairs
agree with themselves exchanged).

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

Two B-spline surfaces
---------------------

The general intersection of two B-spline surfaces, by subdivision and
marching, is the generic path of ADR-0010. Its first stage exists. A
positive-weight B-spline surface lies in the convex hull of its control
net, so if the outward-rounded coordinate boxes of the two nets are further
apart than the pair tolerance along some axis,

.. math::

   \min_i P^{A}_{i,x} - \max_j P^{B}_{j,x} > \tau_{AB}
   \quad\text{(or the same with } A, B \text{ exchanged, or along } y \text{ or } z\text{)},

the surfaces are proved disjoint and the result is empty. That proof holds
over any requested domains. Otherwise the pair returns ``NotImplemented``,
never an empty result it has not proved (``src/intersect/GenericSSI.cpp``).

Every other pair reports ``NotImplemented``.

.. todo::

   The rest of the generic path (Bézier subdivision, interval and
   normal-cone pruning, marching) when Phase 2 builds it.

Related design decisions
------------------------



