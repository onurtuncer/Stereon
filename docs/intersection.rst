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
serves as a *guide*. A branch whose traced points lie within
:math:`2^{-40}` of one another, too short for a guide, is a sliver of a
branch at a corner of the search box or a loop as small. It is reported as
a point of contact, not dropped, its ``Bound`` the farthest corner of its
arcs' certified boxes, which hold it. Each guide point is moved onto :math:`f = 0` by Newton's
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
of it. The quadratic's coefficients and discriminant are taken from their
enclosures when those are tight, and from their exact values otherwise:
formed in doubles, they cancel for a line whose origin lies far from the
surface, :math:`10^8` from a unit sphere say, and the two crossings
collapse onto one point. A touch of a quadric is a ``Tangential`` point: the line stays
outside. A root that double precision cannot place, a line so nearly along
the axis that the quadratic's leading coefficient rounds to zero while it
is exactly positive, fails with ``IntersectionFailed`` rather than
dropping the point. ``tests/occt/IntersectionDifferentialTests.cpp``
compares the three handlers with OCCT's ``IntAna_IntConicQuad`` on
generic lines.

A ``Line`` against a ``Cone`` or a ``Torus`` cannot be decided exactly on
the input doubles: the cone's equation holds the cosine and sine of its
semi-angle. These pairs use the proofs of the B-spline handlers below.

- **A cone over a bounded domain.** The line is a Bézier segment of degree
  1, its poles enclosed. The cone's equation along it is proved as a
  B-spline's is.
- **A cone over an unbounded domain.** The cone has no bound, so its
  quadratic is used, with coefficients enclosed in intervals.
  - A discriminant provably negative is a separation.
  - A discriminant provably positive, with a leading coefficient of strict
    sign, gives two simple crossings, placed in doubles and polished from
    the line re-centred at each.
  - Otherwise the call fails with ``IntersectionFailed``: a line along a
    generator, whose second root lies out near infinity, a line through
    the apex, or a line tangent within rounding. The same line over a
    bounded domain is resolved.
- **A torus.** Every point of the torus lies within :math:`R + r` of its
  centre. So the line is searched only over its stretch inside a sphere
  half a tube wider than that, an outward-rounded stretch whatever the
  domain, as a segment against the torus's quartic.

Either nappe of the cone is crossed, as its parametrisation runs through
the apex.

A ``BSplineCurve``, polynomial or rational, has handlers against all five
elementary surfaces. The surface's implicit equation composed with the
curve is, on each Bézier segment and in homogeneous coordinates, a
polynomial :math:`g` of degree :math:`p` against the plane, :math:`2p`
against a sphere, a cylinder or a cone and :math:`4p` against a torus, of
the sign of the equation since every weight is positive. The cone's
equation, :math:`\cos^2 a\,(x^2 + y^2) - (R \cos a + z \sin a)^2` in its
frame, vanishes on both nappes, as the cone's parametrisation runs through
its apex onto the second; the torus's is
:math:`(x^2 + y^2 + z^2 + R^2 - r^2)^2 - 4R^2 (x^2 + y^2)`. Both are read in
the surface's frame as given, the axis :math:`X \times Y` as the evaluator
rounds it, with :math:`\sin a` and :math:`\cos a` enclosed. Its
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
requested domain, and nothing unresolved becomes one.

Bisection to the rounding of :math:`g` leaves a crossing as far off the
surface as :math:`g`'s enclosed coefficients allow. For a torus's quartic,
that is hundreds of units in the last place of the point. So each decided
crossing then takes one Newton step on its signed distance from the
surface, :math:`(C(t) - S(u, v)) \cdot n / (C'(t) \cdot n)` at the
projection of its rounded point. The step is kept only inside the
crossing's proved bracket, and skipped at a pole of the surface's
parameters, where the normal vanishes. Every point then lies on both
geometries to the rounding of its coordinates.

A curve 3·2⁻⁴⁰ off a
plane is apart; a root that lands on a piece's end, a knot or a dyadic
parameter, is seen from both sides and reported once, undecided; a root
whose bracket reaches past the requested domain is kept, clamped into it,
rather than dropped on its rounded parameter. A curve whose every pole lies
exactly in the plane is an overlap, its pcurve the spline of the projected
poles, clipped to the requested box as described below. A periodic curve
is intersected within its own domain only, and a point's surface
parameters are those of its rounded position, kept or dropped against the
box as the line handlers do.

A ``Circle`` or an ``Ellipse`` against the five elementary surfaces goes
through the same proofs, on its exact rational form. The requested domain,
which must be bounded, is cut into equal arcs of at most a quarter turn,
and each arc :math:`[t_0, t_1]` is the rational quadratic Bézier segment
with end poles at :math:`t_0` and :math:`t_1`, weight 1, and middle pole
:math:`(\cos h\,C + A \cos m + B \sin m,\ \cos h)` in homogeneous
coordinates, :math:`h` the half-angle, :math:`m` the middle angle and
:math:`A`, :math:`B` the radii times the axes. Every one of these is
enclosed from the doubles :math:`t_0` and :math:`t_1`, so what is proved
holds for the conic as given. A root at Bézier parameter :math:`s` is the
angle :math:`m + 2 \arctan(\tan(h/2)(2s - 1))`, its bracket mapped through
that monotone map and widened by its rounding; points are reported on the
conic's angle, which may run over several turns. A conic whose centre and
axes lie exactly in a plane is an overlap there, its pcurve the same
conic in the plane's parameters. A conic on a curved surface, the equator
of a sphere or the meridian of a torus, fails with ``IntersectionFailed``.

A B-spline or a conic in a plane is kept whole when its control polygon,
or its bounding square, lies in the requested box. Otherwise it is
clipped to the box:

- Its range is cut wherever its pcurve meets the line of a finite edge,
  each found by the complete plane-curve intersection below
  (*Plane curves*). Between two cuts, the pcurve stays strictly on one
  side of every edge's line.
- The middle point of each piece decides it. If that point lies within
  rounding of an edge, the clip fails with ``IntersectionFailed`` rather
  than guess.
- A stretch along an edge's line is on the boundary. It is kept where the
  other coordinate keeps it in the box.
- A contact in the box with both neighbouring pieces outside, such as a
  touch from outside or a pass through a corner, is kept as an overlap of
  one parameter, as the line handlers keep one.
- Adjacent kept pieces are merged.
- A requested domain of one parameter has no piece: its point decides it,
  kept within rounding of an edge as a contact on the boundary is.

Every overlap shares the one pcurve and restricts its range.

Every other pair reports ``NotImplemented``; the generic path for curves
against B-spline surfaces waits for Phase 2 with the surface–surface one.

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
rounded values, as a line's roots on a circle are. They are solved from
the quadratic's coefficients and discriminant each to a few roundings
(``ValueOfFormula``): the middle of a tight enclosure, else the exact
value's estimate. In doubles, a line from :math:`(-10^8, 0)` along
:math:`x` gives the unit circle's :math:`C = 10^{16} - 1` as
:math:`10^{16}` and its discriminant, exactly 1, as zero, both crossings
at the centre.

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



