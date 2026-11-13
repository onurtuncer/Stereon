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

   // With the pair's tolerance, the sum of the two faces' (ADR-0050).
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
  Circular sections are clipped to a bounded plane box by the same complete
  planar edge-intersection clipper used for curves lying in a plane.
  Retained arcs keep the cylinder-angle parameter and are certified against
  both surfaces. Isolated box contacts are retained as singleton ranges,
  consistent with clipped rulings; uncertain piece classification fails
  with ``IntersectionFailed``. Partial ellipse clipping remains unsupported.
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

A declared configuration says so in the result:

- **A coincidence declared within the tolerance** is an overlap with
  ``Tolerant`` evidence. Its ``Distance`` is how far apart the profiles are:
  the offset of two parallel profile lines, or :math:`d + |R_1 - R_2|` for
  two profile circles whose centres are :math:`d` apart. Only profiles that
  are equal in the doubles give ``Exact`` evidence, as for parallel planes
  and cylinders.
- **A tangent circle declared within the tolerance** lies up to the gap
  away from one surface, by construction. It is certified against the pair
  tolerance, not the resolution, and its ``Bound`` records the deviation
  actually proved, with ``BoundCause`` ``TolerantCoincidence`` when that is
  past the resolution.

The generated corpus sweeps (``corpus/analytic-ssi/sweep-*``) found both
defects. Concentric spheres whose radii differed within the tolerance were
reported as an ``Exact`` overlap at distance 0. A sphere centred on a
cylinder's axis, the radii within the tolerance but not equal, failed the
whole intersection: the declared tangent circle was certified against the
resolution of 10\ :sup:`-5` mm and missed it. Regression tests in
``tests/intersect/ElementaryPairsTests.cpp`` cover both, and the sweeps
keep running in ``tests/intersect/CorpusTests.cpp``.

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
   F_{\text{torus}} = (\rho^2 + h^2 + R^2 - r^2)^2 - 4 R^2 \rho^2,

   F_{\text{ellipsoid}} = (BC\,p)^2 + (AC\,q)^2 + (AB\,h)^2 - (ABC)^2,

the last for an ``Ellipsoid`` (ADR-0048) with semi-axes :math:`A`,
:math:`B` and :math:`C`, multiplied through so that no division enters an
enclosure.

The local coordinates are a linear map :math:`L` of :math:`x`. So the
chain rule gives :math:`f_u = \nabla F \cdot L(S_u)` and
:math:`f_v = \nabla F \cdot L(S_v)` however the frame is rounded, and every
interval evaluation of :math:`f` is a rigorous enclosure for the nominal
:math:`F`.

*The search box.* The carrier is chosen so that the search is over a
bounded box. A torus, sphere or ellipsoid is compact, so its chart is
bounded. A plane, cylinder or cone carrier against a compact implicit side
is bounded by the ball that holds that side, the ellipsoid's of radius its
longest semi-axis. The ellipsoid's chart has poles, as the sphere's does,
and its implicit equation is a quadric, so it is preferred as the implicit
side and taken as the carrier only against another ellipsoid; on the
implicit side its requested box is tested in the coordinates
:math:`(x/A, y/B, z/C)`, where its longitude and latitude are the unit
sphere's. A cylinder or cone carrier against a cylinder or
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

A plane and an ellipsoid
~~~~~~~~~~~~~~~~~~~~~~~~

A ``Plane`` against an ``Ellipsoid`` (ADR-0048) has a closed form. In the
ellipsoid's frame :math:`(e_1, e_2, e_3)`, its axis as the evaluator rounds
it, the map :math:`s = (x/A, y/B, z/C)` takes the ellipsoid onto the unit
sphere and the plane :math:`n \cdot (x - O) = 0` onto the plane
:math:`m \cdot s = d`, with
:math:`m = (A\,n \cdot e_1,\ B\,n \cdot e_2,\ C\,n \cdot e_3)` and
:math:`d = n \cdot (O - c)`. The plane meets the unit sphere in a circle,
touches it at a point, or misses it, as :math:`d^2` is below, at or above
:math:`\lVert m \rVert^2`: a polynomial in the input doubles, decided
exactly. So a plane :math:`2^{-40}` beyond a vertex of the ellipsoid is a
proved separation, a plane through the vertex a touch at one
``IntersectionPoint``, and a plane :math:`2^{-40}` inside a section.

The circle maps back onto an ellipse, the affine image of a circle. Its
conjugate semi-diameters are turned to principal axes, and the ellipse is
the branch's 3D curve, exact, with the exact ``Ellipse2d`` as its pcurve
on the plane. Its pcurve on the ellipsoid is not a conic in
:math:`(u, v)`. It is the cubic interpolant through the section's
longitudes and latitudes, sampled more densely until ``CertifyDeviation``
proves it within half the resolution of the ellipsoid along the ellipse,
as the trace fits its pcurves (ADR-0009, D3). The section is cut where it
crosses the chart's seam, the half-plane :math:`y = 0, x > 0` of the
frame, found in closed form, so each branch's pcurve stays in the base
turn :math:`[0, 2\pi]`: a section round the axis is one branch over the
whole turn, one beside it zero or two. Each branch is certified on both
surfaces and clipped to both requested boxes.

Two cases fail closed:

- a plane through a pole of the chart, :math:`c \pm C e_3`, decided
  exactly, where the section's longitude has no value
  (``IntersectionFailed``);
- a section so near a pole that its pcurve cannot be certified within the
  sample budget (``IntersectionFailed`` too). A plane whose frame is
  rounded can miss the pole it was made through by a unit in the last
  place; it lands here.

A requested box on the ellipsoid outside its base turn in :math:`u` is
``NotImplemented``. With equal semi-axes the sections lie on the sphere,
as the sphere's own handler finds them.

Space curves
------------

The first v0.2 3D curve handler is ``Line`` against ``Line`` through the
existing ``Intersect(Curve, domain, Curve, domain, ctx)`` overload. Directions
need not be unit vectors; returned parameters retain each input's speed and
orientation. Every other pair of the built-in 3D curves has a handler too
(*Two B-splines in space* and *Lines and conics against other curves*,
below).

Two lines meet where they come within the Context's resolution of each
other, the rule two B-splines follow (*Two B-splines in space*, below), so
a ``Line`` and the same segment as a degree-1 ``BSplineCurve`` meet alike
(owner's decision, 2026-11-08). Over the requested domains the squared
distance of two segments is convex, so its minimum is one point, or a
stretch when the supports are parallel:

- **supports not parallel:** where the supports' closest approach lies in
  both domains it is the minimum, and the lines meet there when it is
  within the resolution. Otherwise the minimum lies on the domains'
  boundary: an end of one segment within the resolution of the other's
  interior, or of the other's end. The nearest such pair is the meeting,
  so two edges sharing a vertex meet at it even when their ends are a
  rounding apart, and a crossing a rounding past a domain's end meets at
  the end;
- **parallel supports** (ADR-0029, *tolerant coincidence*): exactly
  parallel supports within the resolution of each other cannot be told
  from collinear ones, and they overlap over the stretch both domains
  cover. Where that stretch is one parameter they meet at a point, and
  where the domains leave a gap their ends may meet across it. Supports
  that are not exactly parallel never overlap; they meet at a point, as
  above.

Every decision, whether the closest approach lies in a domain, whether a
distance is within the resolution and whether the supports are parallel,
is the sign of a polynomial of degree at most six in the input doubles and
the resolution. Each is read from an interval enclosure, then decided
exactly on expansions when the enclosure straddles zero. A gap of exactly
the resolution meets, and the next double above it does not. Rounded
parameters never decide domain membership. Only the placement is rounded.
A meeting is a ``CurveCurvePoint`` on the first line at its exact
parameter, rounded and clamped to the domain, so the second line's point
lies within the resolution of it. An endpoint-only contact is a point, not
a zero-length overlap. Overlap intervals enclose the exact endpoint
parameters with outward rounding, restricted to the requested domains;
``SameDirection`` records orientation. Whole lines and half-infinite
domains are supported.

OCCT's Boolean decides the same way. ``IntTools_EdgeEdge`` meets two edges,
lines included, by their distance against the edges' tolerances, and lines
that are parallel to within ``Precision::Angular()`` and within that
tolerance give a common part. Here the measure is the Context's
resolution, not an edge's tolerance, and only exactly parallel supports
overlap.

The handler checks finite origins, finite nonzero directions and the
numerical domain of its degree-six expansion comparisons, including finite
parameter bounds and the resolution. Inputs requiring rescaling or
exceeding that domain return ``NotSupported``. Infinite-only domains are
``InvalidInput``. The placed point and its rounded parameter must lie
within the resolution of the first line's exact point, by a conservative
interval bound; otherwise the operation returns ``CertificationFailed``.
This check may reject large-coordinate contacts even when a particular
rounded point is exact.

Two B-splines in space
~~~~~~~~~~~~~~~~~~~~~~

``CurvePair`` (``src/intersect/SplineCurvePair.cpp``) meets two
``BSplineCurve``\ s, polynomial or rational. Two curves in space meet only
exceptionally. Doubles placed to make them meet, such as two edges sharing
a vertex, leave them a rounding apart, so whether they meet exactly is
the wrong question. A meeting is defined by distance instead, at the
Context's resolution (ADR-0001: two points closer than that cannot be told
apart), as for two lines.

A meeting is a critical point of the squared distance, a root of its
gradient

.. math::

   G(s, t) = \begin{pmatrix} (a(s) - b(t)) \cdot a'(s) \\
             -(a(s) - b(t)) \cdot b'(t) \end{pmatrix} = 0,

at which :math:`|a(s) - b(t)|` is within the resolution. Every exact meeting
is such a point, and so is the closest approach of two curves that pass
within the resolution of each other.

The search is Algorithm K's, over pairs of pieces of the two splines'
enclosed Bézier segments. The machinery it shares with the plane handler
(``src/intersect/BezierArc.hpp``) is the pieces, their halving by interval
de Casteljau, their evaluation over a parameter interval and Krawczyk's
step for a two-by-two system. A pair is:

- **apart** when the boxes of the pieces' poles, or their hulls along the
  common normal of the two pieces' chords, are more than the resolution
  apart. This is a fat slab, the space analogue of the fat line;
- **settled by Krawczyk's operator** on :math:`G`, over the pair's box
  inflated by an eighth. The Jacobian

  .. math::

     \begin{pmatrix} a' \cdot a' + (a - b) \cdot a'' & -a' \cdot b' \\
                     -a' \cdot b' & b' \cdot b' - (a - b) \cdot b'' \end{pmatrix}

  is enclosed from the poles and their first and second hodographs, with
  :math:`C' = (X' - C W')/W` and :math:`C'' = (X'' - 2 C' W' - C W'')/W`.
  Near a meeting it is the Gram matrix of the two tangents, regular unless
  they are parallel;
- **halved** otherwise.

A single critical point is bracketed by iterating
:math:`X \leftarrow K(X) \cap X`. Its distance, enclosed over the
bracket, decides it:

- within the resolution, it is a meeting, a ``CurveCurvePoint`` placed at
  the bracket's middle;
- beyond the resolution, it is a near miss and is not reported;
- across the resolution, it cannot be decided at the resolution, and the
  intersection fails with ``IntersectionFailed``.

A bracket is placed only as precisely as the crossing's angle allows. The
gradient's rounding is amplified by the Jacobian's inverse, about
:math:`1/\sin^2\theta` for curves crossing at :math:`\theta`. Two cubics
crossing at 12.5° are placed to about :math:`10^{-11}`, six orders below
the resolution.

A critical point outside the requested domains is not a meeting. There
the distance's minimum over the domains lies on their boundary, and the
curves meet as two lines do (*Space curves*, above):

- **an end against the other curve's stretch.** The end's squared
  distance to the other curve has critical points at the roots of
  :math:`g(t) = (b(t) - P) \cdot b'(t)`, which Krawczyk's operator proves
  over the halved pieces of the other curve, with
  :math:`g' = b' \cdot b' + (b - P) \cdot b''`, and whose distance is
  decided against the resolution as above;
- **two ends** within the resolution of each other.

An end is dropped when the distance is proved to shrink from it inward,
the sign of :math:`(C - Q) \cdot C'` being the wrong one for that side:
the minimum then lies inside, where the pair search finds it. So two
edges sharing a vertex meet at it even when their ends are a rounding
apart. The same goes for a stem ending a rounding short of a bar, and for
nearly parallel segments touching end to end whose supports' closest
approach lies far past both. A critical point just past a segment's end
inside the domain belongs to the neighbouring segment, which proves it.
Critical points proved from overlapping boxes are one, and so are meetings
within the resolution of each other on both curves (ADR-0001), the nearer
kept: a meeting at an end the rounding leaves undecided is found both
inside and on the boundary.

Some pairs cannot be settled: a pair that is neither apart nor settled
when its pieces can no longer be halved, or a search past :math:`2^{16}`
pairs. These are tangent curves, an overlap, or a meeting too nearly
tangential to separate, and they fail with ``IntersectionFailed``, never
dropped and never reported unproved (ADR-0029). A parabola touching its
tangent line is the test case. Two identical splines coincide over the
requested domains' common stretch.

**Lines and splines meet alike.** A ``Line`` and the same segment as a
degree-1 B-spline give the same meetings, placed within the resolution
of each other on each curve. A random property compares them over
crossings, near misses around the resolution and vertices a rounding
apart, and a fixed case over nearly parallel segments touching end to end.
One difference remains, where the spline cannot settle a pair and never
answers differently: parallel supports within the resolution overlap as
lines, while the splines' search fails with ``IntersectionFailed``.

Lines and conics against other curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every pair of ``Line``, ``Circle``, ``Ellipse``, ``Hyperbola``, ``Parabola``
and ``BSplineCurve`` other than two lines goes through the same search and
the same rule. Each curve becomes a list of enclosed Bézier pieces,
computed from the doubles given, so a piece encloses the curve as given
and not a rounded NURBS form of it:

.. list-table::
   :header-rows: 1

   * - Curve
     - Pieces
     - Parameter returned
   * - ``Line``
     - one degree-1 piece
     - the line's own
   * - ``Parabola``
     - one quadratic, its middle pole the blossom of its end parameters
     - the parabola's own
   * - ``Circle``, ``Ellipse``
     - rational quadratic arcs of at most a quarter turn
     - the angle, by :math:`\tan((t-m)/2) = \tan(h/2)(2s-1)`
   * - ``Hyperbola``
     - rational quadratic arcs of width at most 1
     - by :math:`\tanh((t-m)/2) = \tanh(h/2)(2s-1)`

A whole turn of a circle or an ellipse, and the whole domain of a periodic
B-spline, has no ends: the curve runs on through its seam, and a meeting
there, reached from the pieces on both sides, is one. A circle or an ellipse
takes a bounded domain of at most one turn, and more is ``NotSupported``.

A curve with an unbounded domain is clipped first to the parameters that
can come within a ball round the other, bounded curve. The ball is wider
than the resolution, so a clipped end is never a meeting. Two unbounded
domains are ``NotSupported``, so a caller meeting a whole line with a whole
parabola or hyperbola requests a bounded domain on one of them.

Two identical curves of one type coincide over the requested domains'
common stretch; a circle or an ellipse is compared in its base turn. Anything
else the search cannot settle fails with ``IntersectionFailed``: a line
tangent to a circle, two circles equidistant all round, or two different
descriptions of one curve.

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

A ``Line`` against an ``Ellipsoid`` (ADR-0048) is decided the same way. In
the ellipsoid's frame, its axis as the evaluator rounds it, a point at
:math:`w` from the centre has the coordinates :math:`(x, y, z) = (w \cdot
X, w \cdot Y, w \cdot Z)` and lies on the surface where
:math:`(x/A)^2 + (y/B)^2 + (z/C)^2 = 1`. Multiplied through by
:math:`(ABC)^2`, that is a polynomial in the input doubles,

.. math::

   (BC\,x)^2 + (AC\,y)^2 + (AB\,z)^2 = (ABC)^2,

out of which the line cuts a quadratic whose leading coefficient is
positive for every direction. So a line never lies in an ellipsoid: it
meets it in two crossings, one touch or not at all, by the exact sign of
the discriminant, a touch at the end of a semi-axis or over a pole
included. The roots are placed as the sphere's are, from coefficients and
a discriminant good to a few roundings however far the line's origin, and
a point's longitude and latitude are ``Project``'s, those of the
ellipsoid's closest point to the placed one. With equal semi-axes on an
exact frame the quadratic is the sphere's times :math:`r^4`, and the
results are the sphere's: the same count and contact exactly, the same
points to rounding. Scaling a sphere's configuration by a power of two
along an axis of its frame gives the ellipsoid's, parameter for
parameter. B-spline curves, circles and ellipses against an ellipsoid are
below, and a plane against an ellipsoid under *A plane and an ellipsoid*.
An ellipsoid against a cylinder, a cone, a sphere, a torus or another
ellipsoid goes through the parametric–implicit trace; none of the closed
forms for surfaces of revolution applies to it, even when two of its
semi-axes are equal.

A ``Line`` against a ``Cone`` or a ``Torus`` cannot be decided exactly on
the input doubles: the cone's equation holds the cosine and sine of its
semi-angle. These pairs use the proofs of the B-spline handlers below.

- **A cone over a bounded domain.** The line is a Bézier segment of degree
  1, its poles enclosed. The cone's equation along it is proved as a
  B-spline's is.
- **A cone over an unbounded domain.** The cone has no bound, so its
  quadratic :math:`A t^2 + 2 B t + C` is used, with coefficients enclosed
  in intervals, in two steps.

  - From the line's origin: a discriminant :math:`B^2 - A C` provably
    negative is a separation, and a leading coefficient of no strict sign,
    a line along a generator within rounding, fails.
  - Then about the middle of the two roots, the line's point at
    :math:`-B / A` as rounded, its offset from the cone enclosed so that
    the proof holds for the line as given. The discriminant is the same
    number wherever the quadratic is centred. From an origin far from the
    crossings, though, :math:`B^2` and :math:`A C` each carry the square
    of that distance, and their enclosures cancel across zero however far
    apart the crossings are; about the middle of the roots they do not.
    So the decision does not depend on where along the line its origin
    lies. A discriminant provably negative is a separation. One provably
    positive gives two simple crossings, placed in doubles from the
    middle of that enclosure and polished from the line re-centred at
    each.
  - Otherwise the call fails with ``IntersectionFailed``: a line along a
    generator, whose second root lies out near infinity, a line through
    the apex, or a line tangent within rounding. The same line over a
    bounded domain is resolved.

  What remains of the origin is its own rounding: a line whose origin is
  so far that one rounding of its coordinates is comparable with the
  chord fails, as a tangent line does.
- **A torus.** Every point of the torus lies within :math:`R + r` of its
  centre. So the line is searched only over its stretch inside a sphere
  half a tube wider than that, an outward-rounded stretch whatever the
  domain, as a segment against the torus's quartic.

Either nappe of the cone is crossed, as its parametrisation runs through
the apex.

A ``BSplineCurve``, polynomial or rational, has handlers against all five
elementary surfaces and the ellipsoid. The surface's implicit equation
composed with the curve is, on each Bézier segment and in homogeneous
coordinates, a polynomial :math:`g` of degree :math:`p` against the plane,
:math:`2p` against a sphere, a cylinder, a cone or an ellipsoid and
:math:`4p` against a torus, of the sign of the equation since every weight
is positive. The ellipsoid's equation is the line handler's,
:math:`(BC\,x)^2 + (AC\,y)^2 + (AB\,z)^2 - (ABC)^2` in its frame
(ADR-0048). The cone's
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

A ``Circle`` or an ``Ellipse`` against the five elementary surfaces or the
ellipsoid goes through the same proofs, on its exact rational form. The
requested domain,
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

*Verification.* Besides each handler's own tests,
``tests/intersect/AdversarialTests.cpp`` runs properties over the whole
grid of lines, circles, ellipses and polynomial and rational B-splines
against the five surfaces. Each result is validated as geometry, not
counted: every point lies on both geometries to the rounding of its own
size, inside the requested domain and box and in order, and every
overlap's pcurve traces the curve; a handler that cannot decide fails
with ``IntersectionFailed`` and nothing else. The inputs are multiples of
:math:`1/8` in axis-aligned frames, so translating them by a multiple of
:math:`1/8` and scaling them by a power of two are exact, and every exact
decision must agree on the transformed problem while the points move with
it. Reversing the curve keeps the points; a line's origin moved up to
:math:`2^{30}` back along it shifts only its parameters; a periodic
surface's box a turn on gives the same points a turn on; a domain of one
parameter away from every point is a proved separation; a curve in a plane
is clipped to exactly the stretches in a random box; and a line
:math:`2^{-k}` either side of a tangency, :math:`20 \le k \le 50`, gets
the count its distance decides from a sphere or a cylinder. Against a
torus it may fail closed instead, which it does from :math:`k = 40` and
at an exact touch of the tube's top.

Curves against B-spline surfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CurveSplineSurface`` (``src/intersect/CurveSplineSurface.cpp``) meets any
built-in curve, ``Line``, ``Circle``, ``Ellipse``, ``Hyperbola``,
``Parabola`` or ``BSplineCurve``, with a ``BSplineSurface``, polynomial or
rational. The curve becomes the enclosed Bézier pieces of *Lines and conics
against other curves* (``CurvePieces.hpp``), and the surface the enclosed
Bézier patches of the generic SSI (``BezierPatch.hpp``, *Two B-spline
surfaces*). Every pair of a piece and a patch is settled as one of:

- **apart**, when the boxes of the piece's poles and of the patch's net are
  disjoint, so their convex hulls are;
- **one crossing or none**, by Krawczyk's operator on
  :math:`F(s, u, v) = C(s) - S(u, v)` over the pair's parameter box, with
  the Jacobian :math:`[C'(s), -S_u, -S_v]` enclosed over it. A unique root
  is a transverse crossing: the contraction proves the Jacobian regular, so
  the curve's tangent is not in the surface's tangent plane there;
- **cut**, the piece or the patch, whichever is larger in space.

A crossing is a ``CurveSurfacePoint`` with ``Crossing`` contact and
``Decided`` true, at the curve's own parameter (an arc's parameter mapped
back to a conic's) and the surface's parameters. A root found from two
pieces or two patches that share a boundary is one point. Unlike two
curves in space, a curve crosses a surface robustly, so meetings are exact
roots and not decided at the resolution.

What cannot be settled fails with ``IntersectionFailed``, never dropped:

- a pair that can no longer be cut, or a search past :math:`2^{16}` pairs.
  These are a tangent contact, a curve lying in the surface, or a crossing
  too nearly tangential to isolate;
- a crossing within rounding of the end of a requested domain, on the curve
  or on the surface, which cannot be placed on either side of it.

The exception is a seam. When the surface's first and last rows of poles
across a direction are the same, the surface closes on itself there, as a
NURBS cylinder does at angle 0. If the whole of that direction's domain is
requested, a crossing within rounding of its end is placed at the end. It
is reported once, though the patches on both sides find it.

A curve with an unbounded domain is clipped first to a ball around the
surface's patches, as against a bounded curve. The surface is searched over
the requested box met with its own domain; a periodic B-spline surface is
not wrapped.

The tests have three oracles:

- the saddle :math:`z = u^2 - v^2`, a degree-2 B-spline exactly, which a
  line through a known point crosses at that point's parameters;
- a NURBS cylinder from ``ToBSpline``, crossed where the analytic
  cylinder's own handler crosses it;
- a line tangent to the saddle, and one along it, which fail closed.

``tests/occt/IntersectionDifferentialTests.cpp`` compares lines through
known points of a wavy bicubic patch and lines across a NURBS cylinder with
OCCT's ``GeomAPI_IntCS``. Both kernels must find the same number of points,
each of Stereon's within :math:`10^{-7}` of one of OCCT's. The same 48
lines run without OCCT in the kernel's own tests.

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
enclosed Bézier segments, halving both pieces of a pair until it is apart
or settled, as below. Both splines are first moved by one vector, the middle
of their poles' box, when every pole's move is exact in doubles (checked by
Knuth's TwoSum, whose error term is zero exactly then; ``Centred``): the
translation keeps :math:`F` and both parameters, and the enclosures'
rounding, relative to the size of the coordinates, then measures against the
size of the curves. Without it a pair far from the origin, :math:`4\cdot
10^5` away and a few units across, is enclosed only to about
:math:`10^{-10}`, and two arcs :math:`3\cdot 2^{-34}` apart, or a crossing,
could not be settled. When some move would round, the splines are searched
where they are. Points are placed on the input splines. The two outcomes:

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
finds it, or it is no point of the curve).

A bracket that straddles a segment's end, a crossing at or within rounding
of an interior knot of either spline, is proved again with each spline
taken whole, in the curves' own parameters (``AcrossKnots``): :math:`F` at
the box's middle from the segment holding it, and :math:`J(X)` the hull of
every segment's slope over its part of the box, the end segments continued
past the spline's ends. An interior knot is repeated at most :math:`p`
times, so the splines are continuous and :math:`F` is continuous and
piecewise polynomial; then
:math:`F(x) - F(m) = \int_0^1 J\bigl(m + \tau(x - m)\bigr)\,d\tau\,(x - m)`
lies in :math:`J(X)(x - m)`, which is all Krawczyk's existence and
uniqueness argument needs, so the test is as sound across a knot, a corner
included, as inside a segment. The box tried first is the pair's own; if
that fails, boxes about the segment's bracket, an eighth as wide each time
while wider than the bracket. A bracket so proved and inside both splines
is ``Decided``; one across a spline's own end stays undecided (past it the
end segment continued was proved, not the spline), and a root no box
proves is reported with ``Decided`` false, never dropped.

Two brackets that meet on both curves are one root only *by proof*: when
one lies in the box the other was proved the only root of the curves' own
:math:`F` in (a segment's box restricted to the segment, rounded inward;
a whole-spline box restricted to the splines' domains). That merges the
brackets of one crossing found from either side of a knot, and absorbs an
unproved bracket of a continued segment, which if it is a root of the
curves at all is the proved one. Brackets that meet otherwise are reported
as their hull, undecided; equal parameters are never an argument. Every
proved crossing is ``Decided``. A pair
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
       segments, or, at or within rounding of an interior knot of either or
       both splines, proved single with the splines taken whole and its
       bracket inside both splines; a crossing at a spline's own end, or
       one no box proves, is reported undecided

A point whose parameter on a conic lies within rounding of a requested
domain's end may fall on either side of it; a B-spline root whose bracket
reaches past the requested domain is kept, clamped into it, undecided. The
tests are ``tests/intersect/PlanarTests.cpp``,
``tests/intersect/PlanarPairTests.cpp`` and
``tests/intersect/PlanarSplineKnotTests.cpp``: exact configurations and
their neighbours a rounding or :math:`3\cdot 2^{-40}` either side,
tangencies of order two and four, coincidences, rational arcs, exchange and
reversal of the curves, and properties on random geometry (every sign
change of a conic's equation along a random cubic has a point; random conic
pairs agree with themselves exchanged). The knot tests take two cubics
crossing at both their interior knots, polynomial and rational, a cubic
across a polyline's vertex at its own knot, a corner (a knot of full
multiplicity), and crossings a rounding and :math:`2^{-40}` either side of
a knot, with closed-form answers and an independent de Boor oracle in long
double, under exchange, reversal, knot insertion (the crossing's knot
raised to multiplicity three included), translation, a quarter turn and
scaling by :math:`2^{\pm 20}` (ADR-0053 MR-01 to MR-04).

Certificates
------------

Every branch over a bounded range is *certified* (ADR-0009 section 4): its
``Bound`` is a proved bound on its deviation from both surfaces, through
their pcurves, over the whole range, within the Context's resolution. A
branch over an unbounded range stays uncertified until a face's loops
bound it. A branch that cannot be certified fails the intersection with
``CertificationFailed``, never a partial success (ADR-0029).

The bound is the tolerance an edge made from the branch must have, and
``BoundCause`` says why it is that large (ADR-0050, step 0.2-T2):
``CertificateWidening`` for the error of a curve computed or fitted, within
the resolution, and ``TolerantCoincidence`` for a tangency declared within
the pair tolerance, whose bound is past the resolution. A tolerant overlap's
``Distance`` is the gap the edges on the other side grow to cover, cause
``TolerantCoincidence``. The corpus tests check this evidence on every case
(``IntersectTest::CheckToleranceEvidence``).

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
over any requested domains.

When the whole nets overlap, each surface is cut into its Bézier patches over
the requested domains, by Boehm's knot insertion and de Casteljau's algorithm
in outward-rounded interval arithmetic, so that the enclosed poles hold those
of the surface as given. Pairs of patches are examined depth first and cut,
off-centre so that a symmetric branch does not lie along a cut, until each
pair is settled in one of three ways.

*Apart.* The same test on the two patches' nets, as above.

*Transverse.* A direction :math:`L` with

.. math::

   L \cdot (\mathbf{n}_A \times \mathbf{n}_B) > 0

over the whole pair, the normals enclosed from the partial derivatives' nets.
Every branch inside advances along :math:`L`, so none closes up inside the
pair and each leaves through its boundary (the loop test of Sederberg and
Meyers on the cones of normals). The roots of the pair's eight boundary
curves against the other patch, isolated by Krawczyk's operator, are then all
the branch ends: none proves that the pair holds no branch, two are one arc,
and more need a smaller pair.

*An isolated touch*, for polynomial surfaces. Every tangential contact solves

.. math::

   S_A - S_B - \lambda\,\mathbf{n}_B = 0, \qquad
   \partial_u S_A \cdot \mathbf{n}_B = 0, \qquad
   \partial_v S_A \cdot \mathbf{n}_B = 0,

and at an isolated touch this system is regular, as the gradient of the
distance is not. A unique solution within the tolerance, where the relative
second form of the two surfaces is definite, is a point of the result.

Arcs join at shared seeds into branches, every free end on a boundary of the
domains. Each arc is marched by Newton's method on :math:`S_A = S_B` and
:math:`L \cdot S_A = s`, halving chords until they hold, and the polylines in
space and on both surfaces are certified by ``CertifyDeviation``, tightened
until they pass. They are polylines because a fitted cubic pcurve can leave a
surface's domain.

Whatever cannot be settled fails: a tangency along a curve, a branch point, an
overlap, any tangency of rational surfaces, and a search past its budget of
:math:`2^{14}` pairs return ``IntersectionFailed``, naming the patches, never
a partial or an unproved empty result (``src/intersect/GenericSSI.cpp``).

A B-spline surface against a plane, cylinder, cone, sphere, torus or
ellipsoid, in either order, goes to the parametric-implicit trace of
*Elementary surfaces*. The B-spline is the carrier and the analytic surface
supplies the implicit equation. The trace's enclosures of :math:`S`,
:math:`S_u` and :math:`S_v` over a box of parameters come from the
B-spline's Bézier patches, which are enclosed by Boehm's knot insertion in
interval arithmetic and cut down to the box by de Casteljau's algorithm.
Rational patches go through the quotient rule. The enclosure is the hull
over every patch the box meets, so a box across a knot holds both one-sided
partials, and the mean-value form stays valid at a knot where the surface
is only :math:`C^0`. At a point, or along an edge of no width, the patches'
nets are evaluated at the enclosed parameters.

The search box is the requested box met with the B-spline's domain, and a
B-spline chart has no seam. Everything else is the trace's: proved
separation, proved transverse arcs, branches fitted and certified by
``CertifyDeviation`` against the B-spline and the analytic surface, and
``IntersectionFailed`` at a box that cannot be resolved (a tangency, a
singular point) or past its budget. It never returns an unproved empty
result (``src/intersect/ImplicitTrace.cpp``,
``tests/intersect/SplineAnalyticTests.cpp``).

.. todo::

   The rest of the generic path (Bézier subdivision, interval and
   normal-cone pruning, marching) when Phase 2 builds it.

Related design decisions
------------------------




Generic SSI audit restriction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Isolated extrema currently fail closed: proximity and a definite relative
second form alone do not exclude a surrounding crossing loop. Isolated-touch
output is disabled pending a proof of the zero set over the patch pair.
Initial patch-pair allocation is capped by the same budget as subdivision.
