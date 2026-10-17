.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Trimmed curves
==============

Most curves in a model are only parts of curves: a line segment, a circular
arc, a stretch of a conic or of a B-spline between two parameters. The
exchange formats name these parts explicitly. STEP has ``trimmed_curve``,
IGES has the circular arc (entity 100), the line segment (110) and the
trimmed conic (104), and OCCT has ``Geom_TrimmedCurve``. ``TrimmedCurve``
(and ``TrimmedCurve2d`` in a surface's parameter plane) is Stereon's version.
It wraps any curve, built-in or not, in an erased ``Curve`` and restricts it
to a parameter range.

Restricting a curve looks trivial, and evaluation is. Projection is not. The
closest point of a piece is generally not the closest point of the whole
curve moved into the range. Most of this chapter is about getting projection
right, and doing it with the guarantees that :doc:`nurbs` gives for a whole
B-spline.

Definition
----------

Let :math:`C` be a curve with domain :math:`D` and :math:`[a, b]` a range
with :math:`a < b`, both ends finite. The trimmed curve is

.. math::

   C_T(t) = C(t), \qquad t \in [a, b],

and the reversed trimmed curve runs the same range backwards:

.. math::

   C_R(t) = C(a + b - t), \qquad t \in [a, b].

Both keep the basis's parametrisation. A parameter of the trimmed curve
names the same point as on the basis, so a parameter computed on one is
valid on the other with no conversion. This is STEP's convention. A
``trimmed_curve`` whose ``master_representation`` is ``parameter`` trims at
the basis's own parameter values, and ``sense_agreement = false`` is
``IsReversed()``. OCCT's ``Geom_TrimmedCurve`` keeps the parameter in the
same way for a forward trim.

Reversal is done by reflecting the range, not by reversing the basis. Each
type's ``Reversed`` chooses its own new parametrisation: a line negates its
direction, and a circle flips its frame. Reflecting the range keeps the
basis untouched, so the reflection :math:`t \mapsto a + b - t` is the only
map between the two curves. By the chain rule, the derivatives of the
reversed curve are

.. math::

   C_R^{(n)}(t) = (-1)^n \, C^{(n)}(a + b - t),

so ``Deriv`` negates the odd orders. A reversed curve's projection is the
reflection of the forward curve's projection,
:math:`t_R = a + b - t_T`, because both name the same point.

Valid ranges
------------

``MakeTrimmedCurve`` validates data that comes from outside the kernel
(ADR-0006). It returns ``InvalidInput`` unless all of these hold:

- the basis is not empty;
- the range is finite and non-empty;
- for a basis that is not periodic, the range lies inside its domain;
- for a periodic basis with period :math:`T` (the circle and the ellipse,
  :math:`T = 2\pi`), the range is at most one period long,
  :math:`b - a \le T`.

The periodic rule lets an arc cross the seam. A circle trimmed to
:math:`[5, 7]` is the arc from angle 5 through angle :math:`2\pi` to angle 7,
and the circle evaluates at 7 as at :math:`7 - 2\pi`. A range longer than
one period would cover part of the curve twice. Its points would then have
two parameters, projection would not be unique, and no exchange format can
express such a piece, so it is rejected. A range of exactly one period is
the whole closed curve.

The constructor of ``BasicTrimmedCurve`` performs no checks. Code inside
the kernel that already knows its range is valid uses it directly, as for
every geometry type.

The exact piece
---------------

Every built-in curve, restricted to a finite range, has an exact NURBS form
(:doc:`nurbs`, *Exact conversion*). For lines, parabolas and B-splines the
form keeps the parameter. A conic's form is a sequence of rational quadratic
arcs. Those arcs pass through the basis's point at every knot, but between
knots their parameter :math:`s` differs from the basis's :math:`t`. When it
is built, a trimmed curve whose basis is a built-in type stores the piece's
NURBS form :math:`N` over :math:`[a, b]`, together with a ``CurveProjector``
for it. A plane curve is first lifted to the plane :math:`z = 0` with
:math:`(x, y) \mapsto (x, y, 0)`. Distances are unchanged by that embedding,
so a closest point in the lifted plane is a closest point in the plane.

The piece is shared and immutable. Copies of a trimmed curve and its reversed
form share it, so copying a trimmed curve copies a pointer, not the
projector's prepared data.

Projection
----------

The problem
^^^^^^^^^^^

Projecting a point :math:`p` onto the trimmed curve means minimising

.. math::

   f(t) = \lVert C(t) - p \rVert^2, \qquad t \in [a, b].

:math:`f` is continuous and :math:`[a, b]` is compact, so a minimum exists.
If the minimum lies at an interior point, then :math:`f'(t) = 0` there:
:math:`C(t) - p` is orthogonal to :math:`C'(t)`. The minimiser is therefore
either one of the two ends or an interior *stationary* point. The first
method that comes to mind takes the basis's global projection :math:`t_g`
and clamps it to :math:`[a, b]`. That is correct only when :math:`t_g` falls
inside the range or :math:`f` has no local minimum inside it. Neither holds
in general.

A counterexample needs a curve whose :math:`f` has more than one local
minimum. The circle does not qualify. For a circle,
:math:`f(t) = r^2 + d^2 - 2rd\cos(t - \theta)` has exactly one minimum and
one maximum per period, so on an arc the minimum is at :math:`\theta` or at
an end. Take instead the ellipse with radii 3 and 1, the point
:math:`p = (0, 0.1)`, and the lower half :math:`[\pi, 2\pi]`. The global
closest point is the top vertex, at :math:`t_g = \pi/2`, which is outside
the range. Clamping gives an end, :math:`(\pm 3, 0)`, at distance about 3.
But the bottom vertex, at :math:`t = 3\pi/2`, is a local minimum of
:math:`f`, inside the range, at distance 1.1. Clamping returns a point
almost three times as far as the right answer. Points near the chord of an
arc that bends more than a circle, and points near any B-spline with
several bumps, fail the same way.

The method
^^^^^^^^^^

For a built-in basis, ``Project`` works on the exact piece:

1. Project :math:`p` onto the piece's NURBS form with its projector. The
   B-spline projection of :doc:`nurbs` is global: it finds a closest point
   of the whole piece, with no starting guess, by bounding each Bézier
   segment. This gives :math:`s^\ast` and the point
   :math:`q = N(s^\ast)`.

2. Map :math:`s^\ast` back to the basis's parameter :math:`t^\ast`:

   - for a line, a parabola or a B-spline, the piece keeps the parameter,
     so :math:`t^\ast = s^\ast`;
   - for a conic, :math:`q` lies on the basis up to rounding. The basis's
     own closed-form projection of a point on the basis returns that
     point's parameter, because a line or a conic does not cross itself on
     one period and a point is its own closest point. ``Project`` therefore
     sets :math:`t^\ast = \operatorname{Project}(C, q)`. For a periodic
     basis, the result is then moved by whole periods into :math:`[a, b]`:
     the basis answers in :math:`[0, 2\pi)`, while an arc across the seam
     lives in :math:`[5, 7]`.

   Rounding can leave :math:`t^\ast` a few units in the last place outside
   the range, so it is clamped. The error of the mapped parameter is of the
   order of the rounding of :math:`q` divided by the speed
   :math:`\lVert C'(t^\ast) \rVert`.

3. Compare :math:`t^\ast` with the two ends and return the closest of the
   three. Ties go to the earlier candidate, in the order
   :math:`a, t^\ast, b`, so the smaller parameter wins. This is the rule
   of the B-spline projection, and it keeps the answer deterministic
   (ADR-0007). Step 1
   already finds a global minimiser on the piece, so this comparison changes
   the answer only where rounding in steps 1 and 2 made an end slightly
   closer. It costs two evaluations.

On a reversed curve, the result is reflected: :math:`a + b - t^\ast`.

User curve types
^^^^^^^^^^^^^^^^

A basis whose type is not built in has no NURBS form the kernel knows, so
it has no piece, and ``GetPiece()`` is null. This includes a trimmed curve
wrapped in a ``Curve`` and trimmed again. ``Project`` then takes the
basis's global projection :math:`t_g`, moves it into the range (by periods,
then clamping), and compares it with the two ends. The answer is exact when
the global minimiser of :math:`f` lies in the range, or when :math:`f` has
no local minimum inside the range, because then the minimum is at an end.
It is approximate otherwise, as in the ellipse example. A user type that
needs exact projection on pieces should provide a NURBS form, or be
trimmed by its own code.

Bounding box
------------

A NURBS curve with positive weights lies in the convex hull of its poles,
because its rational basis functions are non-negative and sum to one. The
box of the piece's poles therefore contains the piece, and ``BoundingBox``
returns it. This box is finite even when the basis is unbounded: a line
segment's box is the segment's box, while the line's is
``Box3::Everything()``. It is also much tighter than the basis's box for a
short arc. For a user type, the basis's box is returned, which is valid but
loose.

Arc length
----------

The speed of the trimmed curve is the speed of the basis, and a reflection
does not change it, :math:`\lVert C_R'(t) \rVert = \lVert C'(a + b - t)
\rVert`. ``ArcLength`` integrates it with the adaptive Gauss–Legendre rule
of :doc:`properties`. That rule splits the range at the basis's knots,
where the speed may have a kink. For a reversed curve, the knots are
mapped through the reflection. A circle trimmed to :math:`[5, 7]` has
length :math:`2r`, and a piece and its reverse have the same length.

Transformation and conversion
-----------------------------

A rigid motion or similarity :math:`x` commutes with restriction:
:math:`x(C)` restricted to :math:`[a, b]` is :math:`x` applied to
:math:`C_T`, at the same parameters. ``Transformed`` therefore moves the
basis, keeps the range and the sense, and rebuilds the piece.

``ToBSpline`` returns the piece's NURBS form. For a reversed curve it
returns the form reversed with its knots reflected,
:math:`u \mapsto a + b - u`, so that its domain is still :math:`[a, b]`
and it runs as the trimmed curve runs. For a user type there is no form,
and the result is ``NotSupported``.

Correspondence with the exchange formats
----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Source
     - Stereon
   * - STEP ``trimmed_curve``, trimmed by parameter values
     - ``MakeTrimmedCurve(basis, {t1, t2})``, reversed when
       ``sense_agreement`` is false.
   * - STEP ``trimmed_curve``, trimmed by points
     - The reader projects each trimming point onto the basis to find its
       parameter, then trims as above.
   * - IGES 110, line
     - A ``Line`` from the start point to the end point, trimmed to
       :math:`[0, 1]`.
   * - IGES 100, circular arc
     - A ``Circle`` in the entity's plane, trimmed from the start angle to
       the end angle. The end angle is raised by :math:`2\pi` when it is not
       greater than the start angle, so that the arc runs counter-clockwise
       and crosses the seam when it must.
   * - IGES 104, conic arc
     - The conic from its implicit form, trimmed between the parameters of
       its start and end points.
   * - OCCT ``Geom_TrimmedCurve``
     - ``TrimmedCurve``. A forward trim has the same parameters. For a
       reversed trim, OCCT reverses the basis and maps the range through
       the basis's ``ReversedParameter``, while Stereon keeps the basis
       and reflects the range.
