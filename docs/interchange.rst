.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Interchange representations of curves
=====================================

Exchange formats store some curves in forms that the kernel does not keep as
types of their own. IGES has two that matter. The *parametric spline*
(entity 112) is a piecewise polynomial in the power basis. The *conic arc*
(entity 104) is a conic given by its implicit equation and two end points.
The ``:Interchange`` partition of ``Stereon.Curve`` converts both to kernel
types: the spline to a ``BSplineCurve``, exactly, and the conic to the
``Ellipse2d``, ``Hyperbola2d`` or ``Parabola2d`` that it describes. A
reader then trims the result with :doc:`trimmed`. This chapter derives the
conversions and states what they guarantee.

Piecewise polynomials
---------------------

The representation
^^^^^^^^^^^^^^^^^^

IGES 112 gives :math:`N` segments with breakpoints
:math:`t_0 < t_1 < \dots < t_N`. On segment :math:`i`, each coordinate is a
cubic in the *local* parameter :math:`s = t - t_i`:

.. math::

   x(t) = a_{i,0} + a_{i,1} s + a_{i,2} s^2 + a_{i,3} s^3,
   \qquad t \in [t_i, t_{i+1}],

and likewise for :math:`y` and :math:`z` (the entity's ``AX, BX, CX, DX``
and so on). ``PiecewisePolynomialToBSpline`` takes the breakpoints and, for
each piece, its vector coefficients :math:`\mathbf{a}_{i,j}`. Pieces may
have any degree. The curve takes the highest degree :math:`n` among them,
and at least 1. A piece of lower degree is raised by giving it zero leading
coefficients, which changes nothing in the power basis.

Power basis to Bernstein form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`h = t_{i+1} - t_i` and write :math:`s = h w` with
:math:`w \in [0, 1]`. The piece is then
:math:`\sum_j \mathbf{c}_j w^j` with :math:`\mathbf{c}_j = \mathbf{a}_j h^j`.
Each monomial expands in the Bernstein polynomials
:math:`B_{k,n}(w) = \binom{n}{k} w^k (1-w)^{n-k}` of degree :math:`n` as

.. math::

   w^j = \sum_{k=j}^{n} \frac{\binom{k}{j}}{\binom{n}{j}} B_{k,n}(w),

which follows from expanding :math:`w^j = w^j (w + (1 - w))^{n-j}`
:cite:p:`piegl1997nurbs`. Collecting terms gives the Bézier poles of the
piece:

.. math::

   \mathbf{b}_k = \sum_{j=0}^{k} \frac{\binom{k}{j}}{\binom{n}{j}}
                  \mathbf{a}_j h^j, \qquad k = 0, \dots, n.

The map is exact in exact arithmetic. Its coefficients are positive, and
:math:`\mathbf{b}_0 = \mathbf{a}_0` and
:math:`\mathbf{b}_n = \sum_j \mathbf{a}_j h^j`, the two ends of the piece.
In floating point, the error of each pole is a few units of rounding of
:math:`\sum_j \lvert \mathbf{a}_j \rvert h^j`, the size of the terms
summed. A spline whose pieces are long and whose high-order coefficients
are large loses precision in proportion. That is a property of the power
basis, not of the conversion, and the Bernstein form is the better
conditioned of the two from then on.

Joining the pieces
^^^^^^^^^^^^^^^^^^

Each piece becomes one Bézier segment, and the segments are joined by
interior knots of multiplicity :math:`n`:

.. math::

   U = (\underbrace{t_0, \dots, t_0}_{n+1},
        \underbrace{t_1, \dots, t_1}_{n}, \dots,
        \underbrace{t_{N-1}, \dots, t_{N-1}}_{n},
        \underbrace{t_N, \dots, t_N}_{n+1}).

A knot of multiplicity :math:`n` makes the B-spline interpolate the pole at
that knot. The basis there is only :math:`C^0`, so a segment's last pole is
the next segment's first, and the curve has :math:`nN + 1` poles. On each
span the B-spline is exactly the Bézier segment with the poles above, in
the parameter :math:`w = (t - t_i)/h`. The B-spline therefore evaluates as
the polynomial at every :math:`t`, and it keeps the parameter.

The data may be smoother than :math:`C^0`. IGES 112 usually stores a
:math:`C^2` cubic spline. The conversion keeps the full multiplicity
anyway, because removing a knot is only exact when the pieces really are
smooth there, and rounded data never are exactly. Knot removal with a
tolerance is a separate, approximating step.

Join gaps
^^^^^^^^^

Nothing in the format forces one piece to end where the next begins, and
exporters' data sometimes do not. Where the end
:math:`\mathbf{b}^{(i)}_n` of piece :math:`i` and the start
:math:`\mathbf{b}^{(i+1)}_0` of the next differ, the shared pole is their
mean. The largest distance between them is returned as ``MaxJoinGap``.
Moving a pole by :math:`\boldsymbol{\delta}` moves the curve by
:math:`B_{k,n}(w)\,\boldsymbol{\delta}`, and that is at most
:math:`\lVert\boldsymbol{\delta}\rVert` and zero at the far end of each
neighbouring segment. The curve is therefore within ``MaxJoinGap / 2`` of
the data near a join and exact elsewhere. The caller decides whether a gap
is acceptable: a reader compares it with its model tolerance (ADR-0028) and
reports it, rather than having the conversion fail or silently hide it.

Implicit conics
---------------

The equation
^^^^^^^^^^^^

IGES 104 defines a conic in its definition plane by

.. math::

   Q(x, y) = A x^2 + B x y + C y^2 + D x + E y + F = 0,

plus a start point and an end point on it, and a form number that says
whether it is an ellipse, a hyperbola or a parabola. In matrix form, with
:math:`\mathbf{x} = (x, y)`,

.. math::

   Q(\mathbf{x}) = \mathbf{x}^\mathsf{T} M \mathbf{x}
                 + \mathbf{g}^\mathsf{T} \mathbf{x} + F,
   \qquad
   M = \begin{pmatrix} A & B/2 \\ B/2 & C \end{pmatrix},
   \qquad \mathbf{g} = (D, E).

The equation is homogeneous: :math:`\lambda Q = 0` is the same conic for
every :math:`\lambda \neq 0`. ``ConicFromImplicit`` first divides the
coefficients by the largest magnitude among them. Its answer is then the
same for every scaling, and its allowances are relative to coefficients of
size 1.

Principal axes
^^^^^^^^^^^^^^

:math:`M` is symmetric, so a rotation diagonalises it. With

.. math::

   \theta = \tfrac{1}{2} \operatorname{atan2}(B, A - C),
   \qquad
   \mathbf{e}_1 = (\cos\theta, \sin\theta),
   \qquad
   \mathbf{e}_2 = (-\sin\theta, \cos\theta),

the eigenvalues are

.. math::

   \lambda_1 = A\cos^2\theta + B\sin\theta\cos\theta + C\sin^2\theta,
   \qquad
   \lambda_2 = A\sin^2\theta - B\sin\theta\cos\theta + C\cos^2\theta.

This closed form is stable: ``atan2`` is well defined for every
:math:`(B, A - C)`, including the circle :math:`B = 0, A = C`, where every
direction is principal. The eigenvalues are ordered so that
:math:`\lvert\lambda_1\rvert \ge \lvert\lambda_2\rvert`. In the eigenframe,
with coordinates :math:`x' = \mathbf{e}_1 \cdot \mathbf{x}` and
:math:`y' = \mathbf{e}_2 \cdot \mathbf{x}`, and with
:math:`g_k = \mathbf{g} \cdot \mathbf{e}_k`,

.. math::

   \lambda_1 x'^2 + \lambda_2 y'^2 + g_1 x' + g_2 y' + F = 0.

The type follows from the signs of the eigenvalues. Their product is
:math:`\det M = AC - B^2/4 = -\tfrac{1}{4}(B^2 - 4AC)`, the classical
discriminant up to its factor: positive for an ellipse, negative for a
hyperbola, zero for a parabola.

Central conics
^^^^^^^^^^^^^^

When neither eigenvalue is zero, completing both squares moves the origin to
the centre

.. math::

   \mathbf{c} = -\tfrac{1}{2} M^{-1} \mathbf{g}
              = -\frac{g_1}{2\lambda_1}\mathbf{e}_1
                -\frac{g_2}{2\lambda_2}\mathbf{e}_2,

and leaves :math:`\lambda_1 x''^2 + \lambda_2 y''^2 + F' = 0` with

.. math::

   F' = Q(\mathbf{c}) = F + \tfrac{1}{2}\,\mathbf{g} \cdot \mathbf{c}.

- **Ellipse**, where :math:`\lambda_1, \lambda_2` have the same sign. It is
  real only if :math:`F'` has the opposite sign. Its radii are
  :math:`r_k = \sqrt{-F'/\lambda_k}`. The larger radius belongs to the
  smaller eigenvalue, :math:`\lambda_2`, so ``XAxis`` is
  :math:`\mathbf{e}_2` and ``YAxis`` is a quarter turn counter-clockwise
  from it. The ellipse then runs counter-clockwise, as IGES 104's arcs run
  from the start point to the end point. Equal radii give an ellipse with
  :math:`a = b`. IGES writes a circle as entity 100, and a conic that turns
  out round stays the type its form number named.
- **Hyperbola**, where the eigenvalues have opposite signs. The transverse
  axis is the eigenvector whose eigenvalue :math:`\lambda_t` makes
  :math:`-F'/\lambda_t` positive. Then :math:`a = \sqrt{-F'/\lambda_t}`,
  and :math:`b = \sqrt{F'/\lambda_c}` for the conjugate eigenvalue
  :math:`\lambda_c`. A ``Hyperbola2d`` is one branch,
  :math:`\mathbf{c} + a\cosh t\,\mathbf{X} + b\sinh t\,\mathbf{Y}`, so
  the sign of :math:`\mathbf{X}` chooses the branch. Given a ``branch``
  point, such as the arc's start point, :math:`\mathbf{X}` is turned
  towards it: :math:`(\mathbf{p} - \mathbf{c}) \cdot \mathbf{X} > 0`.

Parabolas
^^^^^^^^^

When :math:`\lambda_2 = 0`, only one square can be completed:

.. math::

   \lambda_1 (x' - x_0)^2 + g_2 (y' - y_0) = 0,
   \qquad
   x_0 = -\frac{g_1}{2\lambda_1},
   \qquad
   y_0 = \frac{g_1^2 / (4\lambda_1) - F}{g_2}.

This is :math:`y' - y_0 = k (x' - x_0)^2` with
:math:`k = -\lambda_1 / g_2`. The vertex is
:math:`x_0\mathbf{e}_1 + y_0\mathbf{e}_2`, and the parabola opens along
:math:`\operatorname{sign}(k)\,\mathbf{e}_2`, which becomes ``XAxis``. A
``Parabola2d`` is :math:`\mathbf{o} + \frac{t^2}{4f}\mathbf{X} +
t\mathbf{Y}`, so matching :math:`\lvert k\rvert = 1/(4f)` gives the focal
length

.. math::

   f = \frac{\lvert g_2 \rvert}{4 \lvert \lambda_1 \rvert}.

The allowance
^^^^^^^^^^^^^

Rounded data almost never have :math:`\lambda_2 = 0` exactly. An exporter
writes a parabola's coefficients to about 16 significant digits, so its
discriminant is off zero by about :math:`10^{-16}` relative to the
coefficients. ``ConicAllowance`` (:math:`2^{-33} \approx 1.2 \times
10^{-10}`) is the named tolerance of the conversion (ADR-0028):

- a conic is read as a parabola when
  :math:`\lvert\lambda_2\rvert \le \texttt{ConicAllowance} \cdot
  \lvert\lambda_1\rvert`. The bound is far above rounding. It is also far
  below the eigenvalue ratio :math:`(b/a)^2` of any real ellipse or
  hyperbola, which would need an axis ratio below :math:`10^{-5}` to
  reach it;
- a central conic is degenerate when
  :math:`\lvert F' \rvert \le \texttt{ConicAllowance} \cdot \max(\lvert F
  \rvert, \tfrac{1}{2}\lvert g_1 c_1 \rvert, \tfrac{1}{2}\lvert g_2 c_2
  \rvert)`, that is, when :math:`F'` is lost in the cancellation of the
  terms that formed it;
- a parabola is degenerate when its axial linear term is
  :math:`\lvert g_2 \rvert \le \texttt{ConicAllowance}`, relative to the
  normalised coefficients.

Degenerate equations
^^^^^^^^^^^^^^^^^^^^

A quadratic equation can describe sets that are not conics. ``ConicFromImplicit``
returns ``InvalidInput`` for each of them rather than returning a curve that
does not satisfy the equation:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Case
     - Example
     - Detected by
   * - no real point
     - :math:`x^2 + y^2 + 1 = 0`
     - :math:`F'` has the eigenvalues' sign
   * - one point
     - :math:`x^2 + y^2 = 0`
     - :math:`F' \approx 0`, ellipse
   * - two crossing lines
     - :math:`x^2 - y^2 = 0`
     - :math:`F' \approx 0`, hyperbola
   * - parallel lines, one line, or none
     - :math:`x^2 - 1 = 0`
     - :math:`g_2 \approx 0`, parabola
   * - a line
     - :math:`x + 2y + 3 = 0`
     - :math:`M = 0`
   * - the whole plane
     - all coefficients zero
     - no largest coefficient

Each point of a converted conic satisfies the equation up to rounding.
``ImplicitResidual(conic, p)``, the value :math:`Q(p)` divided by the
largest coefficient, measures that independently of the scaling, and the
tests check it along every converted conic.

From the entity to a trimmed curve
----------------------------------

A complete IGES 104 reader does four things:

1. converts the coefficients with ``ConicFromImplicit``, passing the start
   point as ``branch``;
2. projects the start and end points onto the conic to find their
   parameters :math:`t_s, t_e`. For an ellipse, :math:`t_e` is raised by
   :math:`2\pi` when it is not above :math:`t_s`, because the arc runs
   counter-clockwise;
3. trims the conic to :math:`[t_s, t_e]` with ``MakeTrimmedCurve``;
4. lifts the plane curve into space by the entity's definition-plane
   offset :math:`z_T` and its transformation matrix.

For IGES 112, a reader collects the coefficients of each segment into
vectors, converts with ``PiecewisePolynomialToBSpline``, and checks
``MaxJoinGap`` against its tolerance.
