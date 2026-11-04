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
- **Every other pair of planes, cylinders, cones, spheres and tori:** exact
  circles, lines and points where the surfaces share an axis or the plane
  passes through a cone's apex or contains a torus's axis, and certified
  B-spline branches traced otherwise; see *Elementary surfaces* below.

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

Elementary surfaces
~~~~~~~~~~~~~~~~~~~

Every other pair of planes, cylinders, cones, spheres and tori goes to
``ElementaryPair`` (``src/intersect/ElementaryPairs.cpp``). It tries three
closed forms in turn. The parametric–implicit trace takes every
configuration that has none, which covers the quartic curves of off-axis
pairs. The closed forms come first because they give exact geometry. They
also take the configurations the trace cannot certify: tangencies along a
curve, isolated points and coincidences.

**Coaxial surfaces of revolution.** Cylinders, cones and tori whose axes are
one line, a sphere centred on it, and a plane across it are all turns of a
*meridian profile* about the axis :cite:p:`hoffmann1989geometric`. In the
half-plane of the distance :math:`\rho` from the axis and the height
:math:`z` along it:

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Surface
     - Profile
     - Parameters
   * - Plane across the axis
     - the line :math:`z = h`
     - :math:`h` the plane's height
   * - Cylinder
     - the line :math:`\rho = r`
     -
   * - Cone
     - the lines :math:`\rho = \pm k (z - z_a)`, :math:`k = |\tan a|`
     - :math:`z_a` the apex's height; both nappes
   * - Sphere
     - the circle :math:`\rho^2 + (z - z_c)^2 = R^2`
     - :math:`z_c` the centre's height
   * - Torus
     - the circle :math:`(\rho - R)^2 + (z - z_c)^2 = r^2`
     - :math:`R \ge r`, so the circle stays in :math:`\rho \ge 0`

Turning a point :math:`(\rho, z)` of both profiles about the axis gives a
circle common to both surfaces. Where the profiles cross, the surfaces
cross along that circle. Where the profiles touch, the surfaces touch along
it: the circle is ``Tangential``, with ``SameSide`` evidence, which is the
configuration of a fillet's torus on the cylinder and the plane it blends.
A meeting on the axis, :math:`\rho = 0`, is an ``IntersectionPoint``, such
as a sphere through a cone's apex. Profiles that coincide give coincident
surfaces, which is possible only for two cones, two spheres or two tori
alike; that is an ``OverlapRegion``, with the affine map that turns one
angle into the other and carries one :math:`v` onto the other.

.. _fig-coaxial-profiles:

.. figure:: _static/tikz/coaxial-profiles.png
   :width: 90%
   :alt: The meridian half-plane with a torus's profile circle, a
         cylinder's vertical line touching it, a horizontal plane line
         crossing it twice and a cone's line through its apex crossing it
         twice.

   Coaxial surfaces in the meridian half-plane. The torus is a circle, the
   cylinder a line touching it (open marker), the plane a line crossing it
   twice and the cone a line from its apex crossing it twice (filled
   markers). Turned about the axis, the touch is the circle along which a
   fillet's torus meets its cylinder, and each crossing is a transverse
   circle.

The meetings are the elementary intersections of lines and circles
(``MeetProfiles``). Every circle is an exact ``Circle``, framed about the
axis. On a surface of revolution its pcurve is the line
:math:`(u_0 + s t, v_0)` of constant :math:`v`, with :math:`s = \pm 1` the
sense the surface's angle turns about the axis; on the plane it is an exact
``Circle2d``. A branch is cut where either pcurve crosses its chart's seam
and is clipped exactly to the requested angular and height ranges. A plane
keeps the whole circle when its box holds it and drops it when its box
misses it or lies inside the disk. Otherwise the plane would need a
certified clip, and the result is ``NotImplemented``, as for plane–cylinder.

Whether the axes coincide is decided exactly on the input doubles: the
cross products of the frames' axes, and of the offset between origins with
an axis, vanish or not. A sphere shares an axis with any line through its
centre, so two spheres always share the line of their centres, and a plane
and a sphere the normal through the centre. A sphere's circles are its
parallels of latitude only when its poles lie on the common axis. Otherwise
their pcurves are not lines, and such a pair is traced.

The meetings themselves cannot be decided exactly. A cone's profile comes
from the sine and cosine of its angle, which no double represents. So the
meetings are decided on the profiles' nominal values. A configuration
within the pair tolerance of a tangency is declared tangent and logged at
``Info``, as a coincidence within it is declared (ADR-0029 section 6).

**A sphere off a surface's axis** (``SphereTouch``). The points of a
cylinder, cone or torus nearest a point off its axis lie in that point's
meridian plane, where the surface is its profile. So the distance
:math:`\delta` from the sphere's centre to the surface is the distance to
the profile's nearest line or circle. When :math:`\delta` exceeds the radius
by more than the pair tolerance, the sphere is a proved separation. When it
equals the radius within the tolerance, the open ball misses the surface,
and the intersection is the nearest points alone: ``IntersectionPoint``\ s
of contact, logged as declared. This covers a ball against a wall, outside
or inside a cylinder, on a cone's generator or on a torus's tube. Otherwise
the sphere crosses the surface, and the pair is traced.

**A plane through a cone's apex** meets the cone in the generators whose
direction :math:`\sin a\, e(u) + \cos a\, \hat{A}` lies in the plane. With
:math:`n_x, n_y` the unit normal's components along the cone's axes, those
satisfy

.. math::

   n_x \cos u + n_y \sin u = -\cot a\, (n \cdot \hat{A}).

When :math:`\sqrt{n_x^2 + n_y^2}` exceeds the right side's magnitude there
are two crossing generators. At equality there is one tangent generator.
Below it the apex alone remains, an ``IntersectionPoint``. A generator is
an exact ``Line`` through the apex. Its pcurve on the cone is
:math:`(u, v_a + t)`, since the generators are unit-speed lines in
:math:`v`, and on the plane it is a ``Line2d``.

**A plane containing a torus's axis** meets the torus in two meridian
circles, of radius :math:`r` about :math:`C + R\, e(u)` at the two angles
:math:`u` where the plane cuts the axis's normal plane. This is decided
exactly. Their pcurves are :math:`(u, t)` on the torus and a ``Circle2d``
on the plane.

The parametric–implicit trace
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When no closed form applies, as for spheres, cones and tori off each
other's axes, planes cutting tori obliquely (spiric sections) or tori
crossing at an angle, the intersection is a curve of degree up to eight
with no closed form. ``ImplicitTrace`` (``src/intersect/ImplicitTrace.cpp``)
pairs one surface's parametrisation :math:`S(u, v)`, the *carrier*, with the
other's implicit equation :math:`F(x) = 0` :cite:p:`hoffmann1989geometric`.
The intersection is then the zero set of

.. math::

   f(u, v) = F(S(u, v))

in the carrier's parameter plane. With :math:`(p, q, h)` the coordinates of
:math:`x - O` along the implicit side's frame and
:math:`\rho^2 = p^2 + q^2`, the implicit equations are those of
``src/intersect/Implicit.hpp``:

.. math::

   F_{\text{plane}} = h, \qquad
   F_{\text{cylinder}} = \rho^2 - r^2, \qquad
   F_{\text{sphere}} = \rho^2 + h^2 - R^2,

   F_{\text{cone}} = \rho^2 \cos^2 a - (R \cos a + h \sin a)^2, \qquad
   F_{\text{torus}} = (\rho^2 + h^2 + R^2 - r^2)^2 - 4 R^2 \rho^2.

The local coordinates are a linear map :math:`L` of :math:`x`. So the
chain rule gives :math:`f_u = \nabla F \cdot L(S_u)` and
:math:`f_v = \nabla F \cdot L(S_v)` however the frame is rounded, and every
interval evaluation of :math:`f` is a rigorous enclosure for the nominal
:math:`F`.

*The search box.* The carrier is chosen so that the search is over a
bounded box. A torus or sphere is compact, so its chart is bounded. A plane,
cylinder or cone carrier against a compact implicit side is bounded by the
ball that holds that side. A cylinder or cone carrier against a cylinder or
cone is bounded by Cauchy's bound on the roots of
:math:`g(v) = F(P_0(u) + v\, w(u))`, a quadratic along each ruling, enclosed
over pieces of the turn. That bound fails where its leading coefficient may
vanish, along an asymptotic direction of the implicit side. There the
intersection may reach infinity, and only a bounded parameter box on one of
the surfaces makes the problem finite; without one the result is
``NotImplemented``.

*Certified subdivision.* The search box is subdivided in the manner of
Snyder :cite:p:`snyder1992interval` and Plantinga and Vegter
:cite:p:`plantinga2004isotopic`. Each box :math:`B` is enclosed with the
meet of its natural interval extension and its mean-value form about the
midpoint, and is then resolved in one of three ways:

- **Excluded.** When :math:`0 \notin f(B)` the box holds no point of the
  intersection. A box that provably lies outside the implicit side's
  requested box is also dropped.
- **Certified.** When :math:`f_v` keeps one sign over :math:`B`, the zero
  set is the graph :math:`v = g(u)` of a function over part of the box's
  :math:`u` range. Each edge across it, :math:`u` fixed, is monotone in
  :math:`v`, so it holds a crossing exactly when the certain signs of
  :math:`f` at its corners differ. An edge along it, :math:`v` fixed, is
  accepted when :math:`f` keeps its sign there, so it holds no crossing, or
  when :math:`f_u` does, so it holds one exactly when its corners' signs
  differ. The graph exists exactly between consecutive crossings in order
  of :math:`u`, so they pair into the box's arcs. The same holds with
  :math:`u` and :math:`v` exchanged when :math:`f_u` keeps its sign.
- **Split.** Otherwise the box is split a little off its middle, at
  :math:`33/64`, along its longer side in space, so that children's corners
  avoid the symmetric points of symmetric configurations.

A box still unresolved at depth 56, or after :math:`2^{18}` boxes, ends the
search with ``IntersectionFailed`` naming it. That happens at a tangency, at
a singular point of the curve such as the node of a lemniscate, or for a
branch through a pole of a parametrisation, and never yields an unproved
result.

Since :math:`f_u` or :math:`f_v` is nonzero on every certified box,
:math:`\nabla F` is not normal to the carrier there. Every traced branch is
therefore a *proved* transverse ``Crossing``, and its ``OppositeSide``
evidence is a theorem, not an assumption. An empty result is a proof of
separation over the search box.

.. _fig-implicit-trace:

.. figure:: _static/tikz/implicit-trace.png
   :width: 100%
   :alt: A flowchart: boxes from a stack are excluded, certified with their
         arcs, or split, and an unresolvable box fails; the arcs are then
         joined, guided, moved onto the curve, fitted and certified.

   The parametric–implicit trace. Left, the certified subdivision of the
   carrier's parameter box. Right, from arcs to certified branches.

*From arcs to branches.* Arcs join where their crossings on a shared edge
agree. The crossings come from opposite sides of the edge's line, each
within the other's edge, and on any one edge there is at most one. The
joined arcs are walked into branches, first from free ends, which lie on the
search box's boundary or at a dropped box, then round the closed loops left.
A branch's 3D curve and both pcurves must share one smooth parameter
(ADR-0009). A chord-length cubic interpolant through the traced points
serves as a *guide*. Each guide point is moved onto :math:`f = 0` by Newton's
method along the fixed direction of :math:`\nabla f` at its projection on
the carrier. Because the direction is fixed, the point found is a smooth
function of the guide's parameter. Where the implicit side's angle,
followed continuously with steps refined below an eighth of a turn, crosses
a whole number of periods, the branch is cut, so that each pcurve keeps to
one chart. A branch through a pole, where the angle turns a half turn in no
distance, is refused. The pieces are fitted uniformly on the guide's
parameter, scaled to arc length, with the samples doubled from 32 until the
midpoints agree with the branch and ``CertifyDeviation`` proves the curve
within the resolution of both surfaces. Beyond :math:`2^{14}` samples the
result is ``CertificationFailed``.

*Cost.* Pruning makes the number of boxes grow with the length of the
intersection over the box size at which arcs certify, not with the area of
the search box. In debug builds every configuration in
``tests/intersect/ImplicitTraceTests.cpp`` runs in under half a second,
fitting and certification included.

*Verification.* ``tests/intersect/ImplicitTraceTests.cpp`` checks the number
of curves where the topology is known: spiric sections of one oval, a
waisted curve and two ovals, a cylinder leaving a sphere through one loop
or two turns. It checks completeness against a dense-sampling oracle, which
brackets every sign change of the implicit function along a
:math:`121 \times 120` grid of the carrier's parameters and requires a
traced point within :math:`10^{-2}` mm of each root. It also checks that the
lemniscate's node fails closed. ``tests/intersect/ElementaryPairsTests.cpp``
covers the closed forms.

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
``PlanarOverlap`` stretches where coincident lines or circles share their
domains, with the parameter correspondence. A point is kept when both
parameters lie in the requested domains, a periodic curve's taken modulo
its period and placed there.

Lines, circles and ellipses against each other decide every configuration
as the sign of a polynomial in the input doubles (RS-0001 G5): the
determinant of two directions, whether a parallel line passes through a
point, the discriminant of the quadratic a line cuts out of a circle or an
ellipse, and the distance between two centres against the sum and the
difference of the radii. The polynomials are written once over a scalar
field and read as an enclosure first and exactly, on heap expansions, only
when the enclosure straddles zero. Two conics, a circle or an ellipse
against an ellipse, give a quartic in :math:`\tan(s/2)` over each half of
the turn; its simple roots are isolated by Bernstein subdivision and each
is certified a crossing by the exact signs of the quartic on either side,
while a root where the signs agree, a touch or a near miss, is reported as
a tangency bracketed to rounding and not decided. A B-spline against a
line or a circle is a polynomial per Bézier segment, whose roots are
bracketed to rounding; nothing there is decided, and a touch of even
multiplicity may be missed, which changes no side. Every other pair reports
``NotImplemented``; a rational B-spline is ``NotSupported``.

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
where a target below the rounding of the enclosures leads. Lines, conics
and non-rational B-splines are covered as curves, lines, circles, ellipses
and non-rational B-splines as pcurves, and planes, cylinders, cones,
spheres and tori as surfaces; anything else is ``NotSupported``, and a
range outside a B-spline's domain is ``InvalidInput``. A hyperbola's
enclosure takes ``sinh`` and ``cosh`` over a span from their monotonicity,
with the values at its ends widened by four units in the last place.

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

A B-spline surface against a plane, cylinder, cone, sphere or torus reports
``NotImplemented``.

.. todo::

   The rest of the generic path (Bézier subdivision, interval and
   normal-cone pruning, marching) when Phase 2 builds it.

Related design decisions
------------------------



