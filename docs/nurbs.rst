.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

B-spline and NURBS curves and surfaces
======================================

Free-form geometry in CAD is almost always a B-spline: a piecewise
polynomial, or with weights a piecewise rational, curve controlled by a
polygon of *poles* (control points). STEP stores every free-form curve and
surface this way, interpolation and approximation produce it, and the
rational form represents the conics exactly, so a kernel can carry circles,
ellipses and free-form shapes in one representation when it has to.

``BSplineCurve`` (and ``BSplineCurve2d`` in a surface's parameter plane)
and ``BSplineSurface`` implement the representation and the classical algorithms of Piegl and
Tiller :cite:p:`piegl1997nurbs`, which this chapter follows and cites by
algorithm number (A2.1 and so on). It adds a projection of its own that
finds the global closest point of a curve without a starting guess, and a
surface projection built on the same Bernstein-form reasoning.

Knot vectors and basis functions
--------------------------------

A B-spline of degree :math:`p` with poles :math:`P_0, \dots, P_n` has a
*knot vector*

.. math::

   U = (u_0, u_1, \dots, u_m), \qquad u_i \le u_{i+1}, \qquad m = n + p + 1,

and :math:`n + 1` *basis functions* defined by the Cox–de Boor recursion

.. math::

   N_{i,0}(t) &= \begin{cases} 1 & u_i \le t < u_{i+1} \\ 0 & \text{otherwise} \end{cases} \\
   N_{i,p}(t) &= \frac{t - u_i}{u_{i+p} - u_i} N_{i,p-1}(t)
              + \frac{u_{i+p+1} - t}{u_{i+p+1} - u_{i+1}} N_{i+1,p-1}(t),

with :math:`0/0 = 0`. Their properties are what everything else rests on:

- **Local support.** :math:`N_{i,p}` vanishes outside
  :math:`[u_i, u_{i+p+1})`, so on a knot span :math:`[u_j, u_{j+1})` only
  the :math:`p + 1` functions :math:`N_{j-p,p}, \dots, N_{j,p}` are non-zero.
  Moving one pole changes the curve on at most :math:`p + 1` spans.
- **Non-negativity and partition of unity.** On the domain
  :math:`[u_p, u_{n+1}]` every :math:`N_{i,p} \ge 0` and
  :math:`\sum_i N_{i,p}(t) = 1`.
- **Continuity.** At a knot of multiplicity :math:`k` the basis functions,
  and hence the curve, are :math:`C^{p-k}`. A knot repeated :math:`p` times
  inside the domain leaves the curve continuous with a possible corner; one
  repeated :math:`p + 1` times would disconnect it, which is why
  ``MakeBSplineCurve`` rejects it.

The **domain** is :math:`[u_p, u_{n+1}]`, where the partition of unity
holds. A *clamped* knot vector repeats each end knot :math:`p + 1` times;
the curve then starts at :math:`P_0` and ends at :math:`P_n`, tangent to
the first and last legs of the control polygon. Unclamped knot vectors are
also valid, and are what a periodic B-spline becomes when it is written out
without its period. A *Bezier curve* is the clamped B-spline with a single
span, :math:`U = (0^{p+1}, 1^{p+1})`, whose basis functions are the
Bernstein polynomials
:math:`B_{i,p}(s) = \binom{p}{i} s^i (1 - s)^{p-i}`.

``FindSpan`` (A2.1) locates the span by binary search. The end of the
domain, which no half-open span contains, is given the last non-empty span,
so that :math:`C(u_{n+1})` is the limit from the left. ``BasisFunctions``
(A2.2) computes the :math:`p + 1` non-zero values on span :math:`s` by
building the triangle of :math:`N_{\cdot, r}` for :math:`r = 0, \dots, p`
level by level. With :math:`\ell_j = t - u_{s+1-j}` and
:math:`\rho_j = u_{s+j} - t`, both recursion terms of level :math:`r` share
the denominator :math:`\rho_{k+1} + \ell_{r-k}`, so each level costs
:math:`r` divisions and the whole triangle :math:`O(p^2)` operations; for a
valid knot vector no denominator of a non-zero term vanishes.
``BasisDerivatives`` (A2.3) computes all their derivatives up to a given
order from the same triangle, using

.. math::

   N^{(k)}_{i,p} = \frac{p!}{(p-k)!} \sum_{j=0}^{k} a_{k,j} N_{i+j,p-k},

with coefficients :math:`a_{k,j}` built from knot differences. Derivatives
above the degree vanish. The tests compare both against the textbook
recursions above, implemented independently.

The degree is limited to ``BSplineMaxDegree`` = 25, OCCT's limit, so that
geometry exchanged with OCCT always fits and the evaluators keep their work
arrays on the stack.

Polynomial and rational curves
------------------------------

A ``BSplineCurve`` is

.. math::

   C(t) = \frac{\sum_{i=0}^{n} N_{i,p}(t)\, w_i P_i}{\sum_{i=0}^{n} N_{i,p}(t)\, w_i},
   \qquad t \in [u_p, u_{n+1}],

with positive weights :math:`w_i`. Without weights (``Weights`` empty, every
:math:`w_i = 1`) the denominator is 1 and :math:`C` is a polynomial
B-spline; with them it is a *NURBS* curve. The rational form is a polynomial
B-spline in homogeneous coordinates: the curve

.. math::

   C^w(t) = \sum_i N_{i,p}(t)\, (w_i P_i, w_i)

in one more dimension, projected back by dividing by its last coordinate.
That is how the algorithms below treat both kinds alike.

The weights earn their keep with the conics. The quarter of the unit circle
from :math:`(1, 0)` to :math:`(0, 1)` is the rational quadratic Bezier curve
with poles :math:`(1, 0), (1, 1), (0, 1)` and weights
:math:`1, \tfrac{\sqrt 2}{2}, 1`, exactly; no polynomial represents it.

Every point of the curve is a convex combination of :math:`p + 1` poles,
because the coefficients :math:`N_{i,p} w_i / \sum_j N_{j,p} w_j` are
non-negative and sum to 1. Two consequences:

- **Convex hull.** The curve lies in the convex hull of its poles, and on
  each span in the hull of the :math:`p + 1` poles active there.
  ``BoundingBox`` is the box of the poles, widened by a bound on the
  rounding of ``Eval`` (``Detail::PoleBoxMargin``,
  :math:`8 (p + 1)^2 \cdot 2^{-53}` times the largest coordinate): the
  computed basis functions sum to 1 only up to a few units in the last place
  per recurrence step, so a computed point can fall that far outside the
  hull of the exact one.
- **Affine invariance.** An affine map :math:`x \mapsto Ax + b` commutes
  with convex combinations, so ``Transformed`` maps the poles and keeps the
  knots and weights, and accepts any affine map, including non-uniform
  scaling and shear, which the analytic types cannot.

Derivatives
-----------

For a polynomial curve the derivatives are those of the basis functions
(A3.2):

.. math::

   C^{(k)}(t) = \sum_{j=0}^{p} N^{(k)}_{s-p+j,p}(t)\, P_{s-p+j},

on the span :math:`s` containing :math:`t`, and zero for :math:`k > p`. A
rational curve :math:`C = A / w` with :math:`A = \sum N_i w_i P_i` and
:math:`w = \sum N_i w_i` has :math:`A = w C`; differentiating with Leibniz's
rule and solving for the highest term gives (A4.2)

.. math::

   C^{(k)} = \frac{1}{w} \left( A^{(k)} - \sum_{i=1}^{k} \binom{k}{i} w^{(i)} C^{(k-i)} \right).

:math:`A^{(i)}` and :math:`w^{(i)}` vanish above the degree, but
:math:`C^{(k)}` does not: a rational curve has non-zero derivatives of every
order, and ``Deriv`` returns them.

Parameters are ``double``; poles and weights may be any ``Scalar``
(ADR-0021). The knot span is a discrete choice made on the parameter, while
the point depends smoothly on the poles, so evaluating a curve whose poles
are active scalars gives :math:`\partial C / \partial P_i`, which for a
polynomial curve is :math:`N_{i,p}(t)` per coordinate. The tests check this
with forward-mode duals.

Knot insertion and what it builds
---------------------------------

Inserting a knot :math:`\bar t \in [u_k, u_{k+1})` once replaces the poles
:math:`P_{k-p+1}, \dots, P_k` by

.. math::

   Q_i = \alpha_i P_i + (1 - \alpha_i) P_{i-1}, \qquad
   \alpha_i = \frac{\bar t - u_i}{u_{i+p} - u_i},

and adds one pole, without changing the curve (Boehm's algorithm; A5.1
inserts a knot several times in one pass). In homogeneous coordinates the
same formula applies to rational curves. Each insertion is a convex
combination, so it is numerically benign. ``InsertKnot`` accepts a knot in
the domain whose multiplicity stays at most the degree.

Repeated insertion builds the other operations:

- **Bezier extraction** (A5.6). Once every knot in the domain has
  multiplicity :math:`p`, the :math:`p + 1` poles of each span are the
  Bernstein coefficients of that span: ``ToBezierSegments`` returns one
  clamped single-span curve per non-empty span.
- **Splitting.** With :math:`\bar t` raised to multiplicity :math:`p`, the
  poles before and after it form two curves meeting at
  :math:`C(\bar t)`; ``Split`` returns them with the original
  parametrisation.

``Reversed`` reflects the knots, :math:`u_i \mapsto (a + b) - u_{m-i}` on
the domain :math:`[a, b]`, and reverses the poles and weights.
``Reparametrized`` maps the knots affinely onto a new domain, keeping the
domain's end knots exact so that a clamped curve stays clamped.

Global projection
-----------------

``Project(c, p)`` must return a parameter of a closest point of the whole
curve to :math:`\mathbf p`. Most kernels sample the curve and refine the
best samples with Newton's method, which can miss a narrow closest point
between samples, and a starting point beyond an inflection can converge to
the wrong local minimum. Stereon instead reduces the problem on each Bezier
segment to the roots of one polynomial, which it can isolate completely.

**The stationarity polynomial.** On a segment, with :math:`s \in [0, 1]`,
write the homogeneous curve in the Bernstein basis,
:math:`A(s) = \sum_i B_{i,p}(s)\, w_i P_i` and
:math:`w(s) = \sum_i B_{i,p}(s)\, w_i`, and let
:math:`g(s) = \norm{A(s) - w(s)\,\mathbf p}^2`. The squared distance is
:math:`f = g / w^2`, and

.. math::

   f'(s) = \frac{g'(s)\, w(s) - 2\, g(s)\, w'(s)}{w(s)^3}.

Since :math:`w > 0`, the stationary parameters are exactly the roots of
:math:`h = g' w - 2 g w'`, a polynomial of degree :math:`3p - 1`
(:math:`2p - 1` for a polynomial curve, where :math:`w' = 0`). Its Bernstein
coefficients follow from two exact identities: the product of Bernstein
polynomials of degrees :math:`m` and :math:`n`,

.. math::

   \Bigl(\sum_i a_i B_{i,m}\Bigr)\Bigl(\sum_j b_j B_{j,n}\Bigr)
   = \sum_{k=0}^{m+n} \Biggl( \sum_{i+j=k} \frac{\binom{m}{i}\binom{n}{j}}{\binom{m+n}{k}} a_i b_j \Biggr) B_{k,m+n},

and the derivative
:math:`\bigl(\sum_i a_i B_{i,m}\bigr)' = m \sum_{i} (a_{i+1} - a_i) B_{i,m-1}`.

**Root isolation.** The Bernstein form has the *variation-diminishing*
property: the number of roots of :math:`h` in :math:`(0, 1)` is at most the
number of sign changes :math:`V` in its coefficient sequence, and has the
same parity (Descartes' rule of signs in the Bernstein basis). So

- :math:`V = 0`: no root in the interval;
- :math:`V = 1` with end values of opposite sign: exactly one root, which
  bisection finds to adjacent doubles;
- otherwise: subdivide at :math:`s = 1/2` with de Casteljau's algorithm,
  whose triangle yields the Bernstein coefficients of both halves, and
  repeat. A root exactly at the midpoint, which neither half would see as
  an interior sign change, is recorded there.

Subdivision shrinks the coefficient sequence towards the polynomial's
values, so every simple root is isolated after finitely many halvings;
multiple roots (where the distance function has an inflection) stop at a
depth of 53 halvings, the precision of a double, and contribute the
midpoint as a candidate.

**Assembly.** ``Project`` extracts the Bezier segments, skips any whose
pole box is farther from :math:`\mathbf p` than the best candidate so far
(the convex hull property makes that a valid lower bound), and evaluates
the actual distance at every isolated root and at the segment ends. The
smallest wins; an exact tie goes to the smaller parameter, which keeps the
result deterministic (ADR-0007).

Every candidate is a true stationary point or end point, and every
stationary point is a candidate, so the result is the global minimum up to
the rounding of the coefficients and of the evaluation. The tests require
the projected distance to be no larger than the smallest distance over 4 000
samples, and the OCCT differential tests require it to match
``GeomAPI_ProjectPointOnCurve`` on 1 372 probe points across clamped,
unclamped, rational and degree-5 curves.

Surfaces
--------

A ``BSplineSurface`` is the tensor product of two B-spline bases,

.. math::

   S(u, v) = \frac{\sum_{i=0}^{n} \sum_{j=0}^{m} N_{i,p}(u)\, M_{j,q}(v)\, w_{ij} P_{ij}}
                  {\sum_{i=0}^{n} \sum_{j=0}^{m} N_{i,p}(u)\, M_{j,q}(v)\, w_{ij}},

over :math:`[u_p, u_{n+1}] \times [v_q, v_{m+1}]`, with degree :math:`p`
and knots ``KnotsU`` along :math:`u`, degree :math:`q` and ``KnotsV`` along
:math:`v`, and the poles stored row by row, pole :math:`(i, j)` at index
:math:`i (m + 1) + j`. Everything said about curves carries over per
direction: local support on a pair of spans, the partition of unity, the
convex hull (``BoundingBox`` is the widened pole box, with a margin of
:math:`8 (p + q + 2)^2 \cdot 2^{-53}` times the largest coordinate), affine
invariance, and the homogeneous form :math:`(w_{ij} P_{ij}, w_{ij})` in which
a NURBS surface is a polynomial one. Rational surfaces represent the
quadrics exactly; ``tests/surface/BSplineSurfaceTests.cpp`` builds a quarter
cylinder from the rational quarter circle.

**Derivatives** (A3.6 and A4.4). With
:math:`A = \sum N_i M_j w_{ij} P_{ij}` and
:math:`w = \sum N_i M_j w_{ij}`, the homogeneous derivatives are sums of
products of basis derivatives,
:math:`A^{(k,l)} = \sum_{ij} N^{(k)}_i M^{(l)}_j w_{ij} P_{ij}`. Leibniz's
rule applied to :math:`A = w S` in both parameters gives

.. math::

   S^{(k,l)} = \frac{1}{w} \Biggl( A^{(k,l)}
     - \sum_{\substack{i \le k,\ j \le l \\ (i,j) \ne (0,0)}}
       \binom{k}{i}\binom{l}{j}\, w^{(i,j)} S^{(k-i,l-j)} \Biggr),

which ``Deriv`` evaluates for the requested order from a table of the lower
ones. Small tables stay on the stack (``Detail::SurfaceTable``).

**Knot insertion** in :math:`u` is curve knot insertion applied to every
column of poles (each column is a curve in :math:`u` with the same knots),
and in :math:`v` to every row (A5.3). ``InsertKnot`` takes the direction,
``ToBezierPatches`` raises every knot of both domains to its degree and
returns the patches, and ``IsoCurve`` returns the curve
:math:`v \mapsto S(u_0, v)` or :math:`u \mapsto S(u, v_0)`, whose
homogeneous poles are :math:`\sum_i N_i(u_0) (w_{ij} P_{ij}, w_{ij})`.

``Reversed`` reverses :math:`u`, which turns the normal
:math:`S_u \times S_v` round; ``Transformed`` maps the poles.

Surface projection
~~~~~~~~~~~~~~~~~~

The closest point of a surface to :math:`\mathbf p` is either on the
boundary of the domain or an interior stationary point of the squared
distance. ``Project(BSplineSurface, Point3)`` finds both kinds patch by
patch, keeping the best candidate so far:

1. **Pruning.** A Bezier patch whose pole box is farther from
   :math:`\mathbf p` than the best candidate cannot hold a closer point
   (convex hull) and is skipped.
2. **Boundaries.** Each edge of a patch is a rational Bezier curve, a row or
   column of its poles; the curve projection above finds its closest point
   globally. Edges between patches are included, so a closest point on a
   patch boundary is found directly.
3. **Interior.** On a patch, with the homogeneous form :math:`A / w`,

   .. math::

      F = (A - w\mathbf p) \cdot (A_u w - A w_u), \qquad
      G = (A - w\mathbf p) \cdot (A_v w - A w_v)

   are, up to the positive factor :math:`2 / w^3`, the partial derivatives
   of :math:`\norm{A/w - \mathbf p}^2`, so the interior stationary points are
   the common zeros of :math:`F` and :math:`G`: polynomials of degree
   :math:`(3p - 1, 3q)` and :math:`(3p, 3q - 1)` whose tensor-product
   Bernstein coefficients follow from the bivariate product formula. The
   search subdivides the patch at its midpoints, carrying the four
   homogeneous coordinate polynomials and :math:`F` and :math:`G` along by
   de Casteljau's algorithm, and discards any sub-patch on which :math:`F`
   or :math:`G` has coefficients of one strict sign (then it has no zero
   there) or whose pole box is too far. This is the exclusion test of the
   Bernstein subdivision solvers of Patrikalakis and Maekawa
   :cite:p:`patrikalakis2002shape`.
4. **Refinement.** Sub-patches that survive six levels of subdivision
   (``LeafDepth``, 1/64 of the patch per direction) start Newton's method
   on :math:`(F, G) = 0` in the form
   :math:`((S - \mathbf p)\cdot S_u, (S - \mathbf p)\cdot S_v) = 0`, whose
   Jacobian is the Hessian of half the squared distance, from their centre,
   kept inside the sub-patch.

The boundary part is certified like the curve projection. The interior part
is certified in what it excludes: a sub-patch is discarded only when it
provably contains no stationary point or no closer point. It is not certified
inside the surviving leaves, where Newton's method converges to the
stationary point when the leaf lies in its basin of attraction; a leaf is a
small fraction of a Bezier patch, but no bound on the basin is computed. The
tests require the projected point to be no farther than the closest of
6 561 grid samples on random surfaces, and the OCCT differential tests
require it to be no farther than ``GeomAPI_ProjectPointOnSurf``'s on 1 029
probes; neither has found a counterexample.

One case is expensive rather than wrong: when many points of a patch are at
almost the same distance, as from the centre of a spherical patch, the
coefficients of :math:`F` and :math:`G` are rounding noise, nothing is
excluded, and all :math:`4^6` leaves start Newton's method. Every candidate
is then equally good, and the result is still correct.

``MakeBSplineSurface`` validates each knot vector as ``MakeBSplineCurve``
does, requires a full grid of finite poles, and the weights to be absent or
one positive finite weight per pole.

Checked construction
--------------------

``MakeBSplineCurve(degree, knots, poles, weights)`` validates data from
outside the kernel (ADR-0006, ADR-0016) and returns ``InvalidInput`` with a
message unless

- :math:`1 \le p \le 25` and there are at least :math:`p + 1` poles;
- there are exactly :math:`n + p + 2` knots, finite and non-decreasing, with
  a non-empty domain;
- no knot value is repeated more than :math:`p + 1` times, and none inside
  the domain more than :math:`p` times;
- the poles are finite, and the weights are absent or one positive, finite
  weight per pole.

``ExpandKnots(values, multiplicities)`` converts the knot form STEP and OCCT
store, distinct values with multiplicities, into the full knot vector.
``InsertKnot``, ``Split`` and ``Reparametrized`` validate their parameters
the same way. The evaluators assume a valid curve and check only their
hot-path preconditions (supported degree, consistent sizes, a finite
parameter in the domain) in debug builds.

A ``BSplineCurve`` holds vectors, so copying it allocates. The type-erased
``Curve`` therefore stores it as a shared immutable payload, whatever its
size: ``FitsInline`` admits only trivially copyable types.

Current limits
--------------

- Periodic B-splines have no periodic flag; they are represented by their
  unclamped equivalent and are not ``PeriodicCurve``\ s or
  ``PeriodicSurface``\ s.
- Surface projection is not certified inside its final sub-patches (see
  *Surface projection*).
- Degree elevation, knot removal and curve fitting are not yet implemented;
  fitting is planned for v0.1.
- Evaluation of a cubic costs about 105 ns (GCC 16.2, -O3, MinGW), measured
  by ``bench/GeometryBench.cpp``; the comparison with OCCT that benchmark
  B1.1 asks for is not yet made.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
