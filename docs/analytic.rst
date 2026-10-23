.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Analytic curves and surfaces
============================

The analytic types are the curves and surfaces with a closed-form
parametrisation: the line, circle, ellipse, hyperbola and parabola; the
plane, cylinder, cone, sphere and torus. They are most of what a mechanical part is made of, they
are what STEP files name explicitly, and their closed forms make evaluation,
derivatives and projection exact up to rounding rather than up to the
tolerance of an iteration. This chapter states each parametrisation, derives
its derivatives, normals and closest-point projection, and explains the
choices behind the conventions. :doc:`geometry` describes the concepts and
wrappers through which all of them are used.

Every type follows the conventions of STEP (ISO 10303-42) and of OCCT's
``Geom`` package, so that geometry exchanged with either keeps its
parameters, and ``tests/occt/GeometryDifferentialTests.cpp`` checks the
values, first and second derivatives and projections of every type against
OCCT 8.0.1.

Frames and handedness
---------------------

Each type is placed by a point and two orthonormal vectors
:math:`\mathbf{x}` and :math:`\mathbf{y}` (``XAxis`` and ``YAxis``). The
surfaces of revolution also need an axis :math:`\mathbf{z}`. Storing it would
cost 24 bytes and admit frames whose three vectors disagree, so instead each
surface stores a flag, ``RightHanded``, and computes

.. math::

   \mathbf{z} = \pm\, \mathbf{x} \times \mathbf{y}

with the sign of the flag (``FrameAxis``). The flag is what lets
``Reversed`` turn a surface's normal round without changing its point set:
negating :math:`\mathbf{y}` and flipping the flag together keep
:math:`\mathbf{z}` and run :math:`u` backwards,

.. math::

   S_{\text{reversed}}(u, v) = S(2\pi - u, v),

and since :math:`\partial_u` changes sign while :math:`\partial_v` does not,
the normal :math:`S_u \times S_v` changes sign. A reflection, whose linear
part has a negative determinant, flips the flag in ``Transformed`` so that
the new axis is the image of the old one.

All angle-parametrised types share one more device. Their :math:`n`-th
derivative in an angle :math:`t` involves :math:`\cos(t + n\pi/2)` and
:math:`\sin(t + n\pi/2)`. Evaluating those directly would round the multiple
of :math:`\pi/2` into the argument; instead ``QuarterTurned`` applies the
rotation :math:`(c, s) \mapsto (-s, c)` exactly :math:`n \bmod 4` times to
:math:`(\cos t, \sin t)`, so the derivatives of every order are exactly as
accurate as the point.

The types are templates over the ``Scalar`` concept (ADR-0021): ``Eval``,
``Deriv``, ``Reversed`` and ``Transformed`` run on the active scalar of an
automatic-differentiation backend as well as on ``double``, and the tests
check every derivative against forward-mode duals of ``Eval``.

Curves
------

The ellipse
~~~~~~~~~~~

An ``Ellipse`` with centre :math:`\mathbf{c}`, radii :math:`a \ge b > 0` and
frame :math:`(\mathbf{x}, \mathbf{y})` is

.. math::

   C(t) = \mathbf{c} + a\cos t\, \mathbf{x} + b\sin t\, \mathbf{y},
   \qquad t \in [0, 2\pi],

periodic with period :math:`2\pi`, and its derivatives are

.. math::

   C^{(n)}(t) = a\cos\!\left(t + \tfrac{n\pi}{2}\right) \mathbf{x}
              + b\sin\!\left(t + \tfrac{n\pi}{2}\right) \mathbf{y}.

The parameter :math:`t` is the *eccentric angle*: the polar angle of the
point on the circle of radius :math:`a` that projects onto :math:`C(t)`
along :math:`\mathbf{y}`. It equals the polar angle of :math:`C(t)` only at
the four vertices. With :math:`a = b` the ellipse is the ``Circle`` of the
same frame, and the tests require the two to agree in value, derivatives and
projection.

**Bounding box.** Along a coordinate axis :math:`i` the ellipse's offset
from its centre is :math:`a x_i \cos t + b y_i \sin t`, a sinusoid of
amplitude :math:`\sqrt{(a x_i)^2 + (b y_i)^2}`. ``BoundingBox`` computes this
amplitude in ``Interval`` arithmetic, which rounds outward, so the box
contains the exact ellipse and is tight to within a few units in the last
place.

**Projection.** The closest point of an ellipse in space to a point
:math:`\mathbf{p}` is the closest point to :math:`\mathbf{p}`'s projection
onto the ellipse's plane, because the squared distance splits into an
in-plane part and a constant normal part. In the plane, with
:math:`(y_0, y_1)` the coordinates of the offset in the frame, the problem
is symmetric in both axes, so it is solved for :math:`y_0, y_1 \ge 0` and
reflected back. Following Eberly :cite:p:`eberly2013ellipse`, the closest
point :math:`(x_0, x_1)` satisfies the Lagrange condition

.. math::

   (y_0 - x_0, y_1 - x_1) = s\, \left(\frac{x_0}{a^2}, \frac{x_1}{b^2}\right),

so :math:`x_0 = a^2 y_0 / (s + a^2)` and :math:`x_1 = b^2 y_1 / (s + b^2)`
for a multiplier :math:`s`. Substituting into the ellipse equation and
scaling by :math:`b^2` gives one equation in :math:`s`:

.. math::

   F(s) = \left(\frac{r_0 z_0}{s + r_0}\right)^2
        + \left(\frac{z_1}{s + 1}\right)^2 - 1 = 0,
   \qquad z_0 = \frac{y_0}{a},\; z_1 = \frac{y_1}{b},\; r_0 = \frac{a^2}{b^2}.

For :math:`y_0, y_1 > 0`, :math:`F` decreases strictly on
:math:`s > -1` from :math:`+\infty` to :math:`-1`, so it has exactly one root
there, and the root lies in :math:`[z_1 - 1,\ \sqrt{(r_0 z_0)^2 + z_1^2} - 1]`
when the point is outside the ellipse (:math:`F(0) > 0`) and in
:math:`[z_1 - 1, 0]` when it is inside. ``Detail::EllipseRoot`` bisects
this bracket until the midpoint equals an end point, which takes at most
1 074 steps (one per binary digit between the smallest subnormal and 2)
and needs no starting guess or convergence tolerance. Newton's method would
converge faster, but it can overshoot the pole at :math:`s = -1` or stall
near the evolute, and a projection that occasionally returns a wrong point
is worse than one that always takes a few microseconds.

The axes need their own cases. A point on the minor axis
(:math:`y_0 = 0`) projects to the minor vertex on its side. A point on the
major axis (:math:`y_1 = 0`) projects to the vertex when
:math:`a y_0 \ge a^2 - b^2`, that is, outside the *evolute*, the curve of the
ellipse's centres of curvature; inside it there are two closest points,
mirror images across the axis, at :math:`x_0 = a^2 y_0 / (a^2 - b^2)`, and
``Project`` returns the one with :math:`t \in (0, \pi)`. The centre is
equally close to both minor vertices and gets :math:`t = \pi/2`, except on a
circle, where every point is equally close and it gets :math:`0`, as
``Project(Circle)`` does. These choices make the result a function of the
input, which ADR-0007's determinism requires.

The parameter is recovered from the closest point as the polar angle of
:math:`(x_0 / a, x_1 / b)`, wrapped into :math:`[0, 2\pi)` by
``WrappedAngle``. A tiny negative angle plus :math:`2\pi` rounds to
:math:`2\pi` itself, the same direction as 0 but outside the half-open
period, and is mapped to 0.

The hyperbola
~~~~~~~~~~~~~

A ``Hyperbola`` with centre :math:`\mathbf{c}`, radii :math:`a, b > 0` (in
either order) and frame :math:`(\mathbf{x}, \mathbf{y})` is the branch

.. math::

   C(t) = \mathbf{c} + a\cosh t\, \mathbf{x} + b\sinh t\, \mathbf{y},
   \qquad t \in \R,

of :math:`(x/a)^2 - (y/b)^2 = 1` on the side of :math:`+\mathbf{x}`, with
asymptotes along :math:`a\mathbf{x} \pm b\mathbf{y}`. Differentiation swaps
:math:`\cosh` and :math:`\sinh`, so even orders are
:math:`a\cosh t\,\mathbf{x} + b\sinh t\,\mathbf{y}` and odd orders
:math:`a\sinh t\,\mathbf{x} + b\cosh t\,\mathbf{y}`. The domain and the
bounding box are everything.

**Projection.** With :math:`(x, y)` the in-plane coordinates of
:math:`\mathbf p`, the derivative of the squared distance vanishes where

.. math::

   (a^2 + b^2) \sinh t \cosh t = a x \sinh t + b y \cosh t.

Substituting :math:`s = e^t`, :math:`\cosh t = (s + s^{-1})/2`,
:math:`\sinh t = (s - s^{-1})/2`, and multiplying by :math:`4 s^2` turns it
into a quartic,

.. math::

   (a^2 + b^2)\, s^4 - 2 (a x + b y)\, s^3 + 2 (a x - b y)\, s - (a^2 + b^2) = 0,

whose positive roots are the stationary parameters. The closest point is no
farther from :math:`\mathbf p` than the vertex, at distance :math:`d_0`, so
its :math:`y` coordinate :math:`b\sinh t` lies within :math:`d_0` of
:math:`y`, and :math:`\abs{t} \le \operatorname{asinh}((\abs{y} + d_0)/b)`.
``Detail::PolynomialRoots`` converts the quartic to the Bernstein basis on
that range of :math:`s` and isolates its roots with the subdivision method
of :doc:`nurbs`; the nearest root wins. A point on the axis beyond the
centre of curvature of the vertex, :math:`x = a + b^2/a`, has two mirror
closest points.

The parabola
~~~~~~~~~~~~

A ``Parabola`` with vertex :math:`\mathbf{o}`, focal length :math:`f > 0`
and frame :math:`(\mathbf{x}, \mathbf{y})` is

.. math::

   C(t) = \mathbf{o} + \frac{t^2}{4f}\, \mathbf{x} + t\, \mathbf{y},
   \qquad t \in \R,

the parabola :math:`y^2 = 4 f x` with its focus at
:math:`\mathbf o + f\mathbf x`; :math:`t` is the coordinate along
:math:`\mathbf{y}`. Its derivatives are :math:`t/(2f)\,\mathbf x + \mathbf y`,
then :math:`\mathbf x / (2f)`, then zero.

**Projection.** The stationarity condition is the cubic

.. math::

   \frac{t^3}{8 f^2} + t \left(1 - \frac{x}{2f}\right) - y = 0,

and since the closest point is no farther than the vertex and its :math:`y`
coordinate is :math:`t` itself, its roots are searched on
:math:`[y - d_0, y + d_0]`. On the axis, the vertex is closest up to the
centre of curvature :math:`x = 2f`; beyond it two mirror points are.

Both conics follow STEP's ``hyperbola`` and ``parabola`` and OCCT's
``Geom_Hyperbola`` and ``Geom_Parabola``, with which the differential tests
compare their values, three derivatives and projections. ``Transformed``
takes similarities; for the parabola the parameter scales with the
similarity, because it is a length along :math:`\mathbf y`.

Surfaces
--------

The four surfaces of revolution share one structure. With the radial
direction

.. math::

   \mathbf{e}(u) = \cos u\, \mathbf{x} + \sin u\, \mathbf{y},
   \qquad \mathbf{e}^{(n)}(u) = \cos\!\left(u + \tfrac{n\pi}{2}\right)\mathbf{x}
                             + \sin\!\left(u + \tfrac{n\pi}{2}\right)\mathbf{y},

each surface is a *meridian* curve :math:`(\rho(v), h(v))` swept around the
axis:

.. math::

   S(u, v) = \mathbf{o} + \rho(v)\, \mathbf{e}(u) + h(v)\, \mathbf{z}.

The variables separate, so every mixed partial is a product:

.. math::

   \partial_u^{m} \partial_v^{n} S = \rho^{(n)}(v)\, \mathbf{e}^{(m)}(u)
   + [m = 0]\; h^{(n)}(v)\, \mathbf{z},

where :math:`[m = 0]` is 1 when :math:`m = 0` and 0 otherwise. This is how
``Deriv`` computes them, with ``RadialDirection`` for
:math:`\mathbf{e}^{(m)}`. The normal of a right-handed frame,

.. math::

   S_u \times S_v = \rho(v) \left( h'(v)\, \mathbf{e}(u) - \rho'(v)\, \mathbf{z} \right),

points away from the axis wherever the meridian rises (:math:`h' > 0`) at a
positive radius: outward on the cylinder, the sphere and the torus's outer
equator.

Projection reduces to the meridian in the same way. The surface is
symmetric under rotation about its axis, so the closest point to
:math:`\mathbf{p}` lies in the half-plane through the axis and
:math:`\mathbf{p}`, at the angle :math:`u_p` of :math:`\mathbf{p}` around the
axis; in that half-plane, with :math:`\mathbf{p}` at distance :math:`\rho_p`
from the axis and height :math:`h_p` along it, the problem is the 2D closest
point of the meridian to :math:`(\rho_p, h_p)`. A point on the axis is
equally close to every meridian and gets :math:`u = 0`.

.. list-table::
   :header-rows: 1

   * - Surface
     - Meridian :math:`(\rho(v), h(v))`
     - Domain
     - Periods
   * - ``Cylinder``
     - :math:`(r,\ v)`
     - :math:`[0, 2\pi] \times \R`
     - :math:`(2\pi, 0)`
   * - ``Cone``
     - :math:`(r + v\sin\alpha,\ v\cos\alpha)`
     - :math:`[0, 2\pi] \times \R`
     - :math:`(2\pi, 0)`
   * - ``Sphere``
     - :math:`(r\cos v,\ r\sin v)`
     - :math:`[0, 2\pi] \times [-\pi/2, \pi/2]`
     - :math:`(2\pi, 0)`
   * - ``Torus``
     - :math:`(R + r\cos v,\ r\sin v)`
     - :math:`[0, 2\pi]^2`
     - :math:`(2\pi, 2\pi)`

The plane, :math:`S(u, v) = \mathbf{o} + u\,\mathbf{x} + v\,\mathbf{y}`, is
not a surface of revolution; its derivatives are its axes and its
projection is the pair of dot products with them.

The sphere
~~~~~~~~~~

:math:`u` is the longitude and :math:`v` the latitude. The parametrisation
is singular at the poles :math:`v = \pm\pi/2`: every :math:`u` gives the
same point, :math:`S_u` vanishes, and the normal there is the limit
:math:`\pm\mathbf{z}` rather than :math:`S_u \times S_v`. The poles are the
``Degenerate`` edges of ADR-0009, and code that needs a normal must treat
them as such. The closest point of the sphere to :math:`\mathbf{p}` is on
the ray from the centre through :math:`\mathbf{p}`, so the projection is the
pair of angles

.. math::

   u = \operatorname{atan2}(p_y, p_x), \qquad
   v = \operatorname{atan2}\!\left(p_z, \sqrt{p_x^2 + p_y^2}\right)

of the offset in the frame; the centre, equally close to everything, gets
:math:`(0, 0)`. The bounding box is the centre plus and minus the radius,
tight.

The cone
~~~~~~~~

:math:`v` is the signed distance from the reference circle of radius
:math:`r` along a generator, so the generators are unit-speed lines and
every derivative of order 2 or more in :math:`v` vanishes. The semi-angle
:math:`\alpha` satisfies :math:`0 < \abs{\alpha} < \pi/2`; a negative angle
opens the cone towards :math:`-\mathbf{z}`. The apex is at
:math:`v = -r / \sin\alpha`, where the parametrisation is singular like the
sphere's poles, and the surface continues through it as the second
*nappe*, where :math:`\rho(v) < 0`: the point lies on the opposite side of
the axis from :math:`\mathbf{e}(u)`.

That second nappe is why the cone's projection needs care. The plane
through the axis and :math:`\mathbf{p}` cuts the cone in two full lines
through the apex: the generator at :math:`u_p`, with meridian points
:math:`(r + v\sin\alpha, v\cos\alpha)`, and the generator at
:math:`u_p + \pi`, which appears in the same signed coordinates as
:math:`(-(r + v\sin\alpha), v\cos\alpha)`. Both are unit-speed lines, so the
foot of the perpendicular on each is a dot product,

.. math::

   v_{\text{near}} = (\rho_p - r)\sin\alpha + h_p\cos\alpha, \qquad
   v_{\text{far}} = -(\rho_p + r)\sin\alpha + h_p\cos\alpha,

and ``Project`` keeps the nearer, with :math:`u = u_p + \pi` when it is the
second. Near the apex the second line can be the nearer one; for a point on
the axis the two are equally near and the first wins. The tests check the
distance against the independent formula
:math:`\min \abs{\rho\cos\alpha \mp (h - h_{\text{apex}})\sin\alpha}` for
the distance to two lines through the apex.

The cone is unbounded, so its bounding box is ``Box3::Everything()``.

The torus
~~~~~~~~~

``Torus`` supports ring tori (:math:`R > r > 0`) and horn tori
(:math:`R = r > 0`). ``MakeTorus`` rejects spindle tori (:math:`R < r`).
At the horn's inner point, :math:`v = \pi`, all values of :math:`u` map to
the centre. Normals and curvatures are undefined there and return
``DegenerateResult``. Nearby regular points use a cancellation-resistant
half-angle expression for the radial coordinate. Projection of the centre
chooses :math:`(u,v) = (0,\pi)`. Bounding boxes and rational B-spline
conversion also support horns. This is surface geometry support; it does
not introduce a toroidal solid primitive or Boolean support at the pinch.

The closest point lies in the half-plane at :math:`u_p`, and within it on
the tube's cross-section circle of centre :math:`(R, 0)` and radius
:math:`r`, so

.. math::

   v = \operatorname{atan2}(h_p,\ \rho_p - R).

The half-plane at :math:`u_p` holds the nearer of the two cross-sections the
plane through the axis cuts, because :math:`\rho_p \ge 0`. A point on the
tube's centre circle is equally close to the whole cross-section and gets
:math:`v = 0`.

**Bounding box.** Along coordinate axis :math:`i`, with
:math:`A = \sqrt{x_i^2 + y_i^2}`, the offset from the centre is
:math:`(R + r\cos v) A \cos(u - u_i) + r\sin v\, z_i` for some phase
:math:`u_i`. Since :math:`R + r\cos v \ge 0` on a ring or horn torus,

.. math::

   (R + r\cos v) A\cos(u - u_i) + r\sin v\, z_i
   \le R A + r (A\cos v + z_i \sin v)
   \le R A + r\sqrt{A^2 + z_i^2},

with equality reachable, so ``BoundingBox`` is tight. It evaluates the bound
in ``Interval`` arithmetic with the axis the evaluator computes, so the box
contains the exact torus even when the stored frame is orthonormal only up
to rounding.

Checked construction
--------------------

Each type is an aggregate, and kernel code that has established the
invariants builds it directly. Data from outside the kernel — a file, a user
interface, a script — goes through the type's checked factory instead
(ADR-0006, ADR-0016): ``MakeLine``, ``MakeCircle``, ``MakeEllipse``,
``MakePlane``, ``MakeCylinder``, ``MakeCone``, ``MakeSphere`` and
``MakeTorus`` return ``Result<T>``, with ``ErrorCode::InvalidInput`` and a
message naming the factory and the value when

- a coordinate, radius or angle is not finite;
- a radius is not positive (the cone's reference radius may be 0, putting
  the apex at the origin), the ellipse's radii are out of order, the cone's
  semi-angle is not strictly between 0 and :math:`\pi/2` in magnitude, or a
  torus is neither a ring nor a horn torus;
- the frame is not orthonormal within ``FrameAllowance``.

.. code-block:: cpp

   const Result<Sphere> ball = MakeSphere(center, radius, xAxis, yAxis);
   if (!ball) {
       std::println(stderr, "{}", ball.error().GetMessage());
       // "MakeSphere: Radius must be finite and positive"
   }

``FrameAllowance`` is :math:`2^{-40} \approx 9 \cdot 10^{-13}`, a
dimensionless relative allowance in the sense of ADR-0028:
:math:`\abs{\mathbf{x}\cdot\mathbf{x} - 1}`,
:math:`\abs{\mathbf{y}\cdot\mathbf{y} - 1}` and
:math:`\abs{\mathbf{x}\cdot\mathbf{y}}` must each be at most the allowance.
A unit vector computed in doubles carries a relative error of a few units of
:math:`2^{-53}` per operation, so frames that have been normalised, rotated
and transformed many times pass, while any frame that is visibly skewed or
scaled is rejected. The factories validate; they never repair, because a
silently orthonormalised frame would be a different surface from the one the
data described.

The evaluators themselves state their preconditions as hot-path checks
(:doc:`geometry`, *Evaluator preconditions*): a violation there is a bug in
the caller, not bad data.

Verification
------------

For every type, ``tests/curve`` and ``tests/surface`` check, over randomly
generated parameters and points (RapidCheck):

- each derivative against the forward-mode dual of the evaluator one order
  lower, in :math:`u`, in :math:`v` and mixed;
- that points satisfy the implicit equation of the surface
  (:math:`\abs{S - \mathbf{c}} = r` for the sphere,
  :math:`(\rho - R)^2 + h^2 = r^2` for the torus, :math:`\rho = \abs{r + h\tan\alpha}`
  for the cone);
- that ``Project`` inverts ``Eval`` on the domain, that its foot is a
  stationary point of the distance, and that its distance equals an
  independent closed form (or, for the ellipse, is no larger than the
  smallest distance over 20 000 samples);
- that ``Reversed`` runs the parameter backwards and turns the normal round,
  that ``Transformed`` commutes with ``Eval`` for similarities and
  reflections, and that the bounding box contains every evaluated point;
- that the factories accept valid data and reject each kind of invalid
  data, ``MakeHyperbola`` and ``MakeParabola`` included;
- that the fast path and the erased wrapper agree to the bit.

The OCCT differential tests compare the value and every first and second
partial at 41 parameters, and the projection distance of 343 points around
each geometry with OCCT's ``GeomAPI_ProjectPointOnCurve`` and
``GeomAPI_ProjectPointOnSurf``, to a relative :math:`10^{-9}`.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0009: Edge curves as surface-intersection curves <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0009-parameter-space-representation.md>`_ (degenerate edges at poles and apexes)
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
- `ADR-0028: Separate numerical tolerances by purpose <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0028-purpose-specific-numerical-tolerances.md>`_
