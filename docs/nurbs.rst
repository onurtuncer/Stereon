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

For a polynomial curve the derivatives are those of the basis functions,

.. math::

   C^{(k)}(t) = \sum_{j=0}^{p} N^{(k)}_{s-p+j,p}(t)\, P_{s-p+j},

on the span :math:`s` containing :math:`t`, and zero for :math:`k > p`.
``Deriv`` does not differentiate the basis functions, though (A3.2); it uses
the fact that the derivative of a B-spline is again a B-spline, one degree
lower, on the same knots (eq. 3.8):

.. math::

   C'(t) = \sum_{i=0}^{n-1} N_{i+1,p-1}(t)\, Q_i, \qquad
   Q_i = \frac{p}{u_{i+p+1} - u_{i+1}}\, (P_{i+1} - P_i).

Applied :math:`k` times to the :math:`p + 1` poles active on the span, this
gives the poles of :math:`C^{(k)}` there (A3.3). Their basis functions, of
degree :math:`p - k`, are a level of the triangle that the recurrence for
the point builds anyway, so ``BasisFunctionsOfAllDegrees`` keeps the whole
triangle and every derivative costs only :math:`O(p^2)` pole differences
after it. On a non-empty span every denominator is positive. The combined
call ``Derivatives(curve, t, n, out)`` returns the point and all
derivatives up to order :math:`n` from one evaluation, the counterpart of
OCCT's ``D1``, ``D2`` and ``D3``; ``Deriv`` is that call for one order. The
work arrays are sized for degrees below 8, which covers the curves met in
practice, because vectors zero themselves and arrays sized for degree 25
would cost more to clear than the arithmetic.

A rational curve :math:`C = A / w` with :math:`A = \sum N_i w_i P_i` and
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
  clamped single-span curve per non-empty span. For a curve clamped at both
  ends it walks the knots once, raising each to multiplicity :math:`p` by
  updating only the current segment's poles and copying the shared ones
  forward, rather than inserting knot by knot into the whole pole array.
- **Splitting.** With :math:`\bar t` raised to multiplicity :math:`p`, the
  poles before and after it form two curves meeting at
  :math:`C(\bar t)`; ``Split`` returns them with the original
  parametrisation.

``Reversed`` reflects the knots, :math:`u_i \mapsto (a + b) - u_{m-i}` on
the domain :math:`[a, b]`, and reverses the poles and weights.
``Reparametrized`` maps the knots affinely onto a new domain, keeping the
domain's end knots exact so that a clamped curve stays clamped.

Degree elevation, knot removal and the other rearrangements
-----------------------------------------------------------

The remaining chapter-5 operations live in ``BSplineOperations.cppm``. All of
them work in homogeneous coordinates, so the rational case costs nothing
extra, and all but knot removal leave the point set unchanged up to
rounding.

- **Clamping.** ``Clamped`` inserts the two domain ends up to multiplicity
  :math:`p` and drops the knots and poles beyond them: on the first
  non-empty span the basis functions depend only on the knots from the
  second copy of :math:`a` on, so the knot before it may be set to
  :math:`a` and everything earlier discarded, and likewise at :math:`b`. The
  result starts at its first pole and ends at its last, on the same domain,
  without any new arithmetic beyond the insertions. The periodic fits of
  :doc:`fitting` produce unclamped curves; this is how to hand them to code
  that expects end poles on the curve.
- **The periodic form.** Stereon has no periodic flag; a closed curve is
  periodic by its data. ``IsPeriodic`` recognises the form: with
  :math:`n + p` poles, the last :math:`p` repeat the first :math:`p`, and the
  knots repeat with the width :math:`T` of the domain,
  :math:`u_{i+n} = u_i + T` up to a few units of the rounding of doubles
  of the knots' size. That rounding must stay below :math:`1/256` of the
  shortest knot interval, or the knots cannot resolve the period: a curve
  whose knots sit so far from the origin (about :math:`10^{13}` for unit
  intervals) is not recognised, whatever its values. ``AsPeriodic``
  builds it from any closed curve: the curve is clamped, its end poles,
  which must lie within the given tolerance, are merged at their midpoint,
  the seam knot takes multiplicity :math:`p` and the poles are rotated by
  :math:`p - 1` so that the seam falls where the clamped curve began. The
  curve on the domain is unchanged except for the merge, which moves it by
  at most half the gap. The extended knot vector is checked as a factory
  checks one: when a shifted copy overflows or the width is not finite,
  ``AsPeriodic`` returns ``NotSupported`` rather than an invalid curve.
  A clamped closed curve is :math:`C^0` at its seam
  and so is its periodic form; the fits produce seams of multiplicity 1,
  :math:`C^{p-1}`, which ``RemoveKnot`` at the seam would recover where the
  curve allows. A rational curve needs equal weights at its two ends;
  scaling one side's weights by a constant changes no point.

Periodic curves
~~~~~~~~~~~~~~~

A curve in the periodic form is periodic, with no flag to say so:
``Period`` returns the width of its domain, and 0 for every other B-spline
curve, so ``BSplineCurve`` and ``BSplineCurve2d`` are ``PeriodicCurve``\ s
whose period their data decide. Everything that reads a period then treats
the curve as it treats a circle:

- **Evaluation** (``Eval``, ``Deriv``, ``Derivatives``, ``EvalAll``,
  ``EvalAllFast``, and the span polynomials of ``ToSpanPolynomials``, which
  record the period) takes every finite parameter. One outside the domain is
  carried into it by whole periods (``NearestInRange``); the periodicity is
  checked only then, so evaluation inside the domain costs what it did.
- **Arc length** (``ArcLength``, ``ParameterAtArcLength``) over any finite
  range, counting the turns it covers and splitting the quadrature at the
  images of the knots in every period, not at the stored knots of one
  (:doc:`properties`).
- **Knot insertion** keeps the form. A parameter stands for all its images,
  and the knot goes in at every one of them: the form is written out over a
  few periods on each side, every image inserted by Boehm's algorithm, and
  one period of knots and poles taken back from the middle, its last
  :math:`p` poles copies of its first :math:`p`. The result has
  :math:`n + r` distinct poles for :math:`r` insertions, the same domain
  and the same period, and is the same curve up to rounding. The
  multiplicity limit counts one period, so the seam of ``AsPeriodic``,
  already at multiplicity :math:`p`, takes no more.
- **Pieces** (``ToBSpline`` over a range of at most one period) may cross
  the seam and lie in any period: the range is moved by whole periods to
  start in the domain, cut from the form written out over one more period
  on each side, and its knots moved back. A piece is a clamped curve, so it
  is not periodic itself. Trimmed curves on a periodic B-spline cross its
  seam the same way, as they cross a circle's (:doc:`trimmed`).
- **Opening and splitting.** ``Opened(curve, t)`` opens the curve at
  :math:`t`, in any period: the clamped curve over :math:`[t, t + T]` that
  traces it once round, its last pole made its first, so that it closes
  exactly and ``AsPeriodic`` turns it back into the form. ``Split`` takes
  any finite :math:`t` and cuts the domain at its image, as it cuts any
  curve; at the seam there is nothing on one side, and it refuses.
- **Degree elevation** (A5.9). ``ElevatedDegree`` splits the clamped curve
  into its Bezier segments, raises each by the convex combination
  :math:`Q_i = \tfrac{i}{q+1} P_{i-1} + (1 - \tfrac{i}{q+1}) P_i` one degree
  at a time, joins the segments with every breakpoint at the new degree's
  multiplicity and then removes each interior knot :math:`p - s` times,
  where :math:`s` was its multiplicity: the elevated curve is as smooth
  there as the original, so these removals are exact in real arithmetic.
  In doubles they are not always: a removal inverts knot insertion, and
  dividing by the ratio of a tiny knot interval to a span magnifies the
  rounding of the poles by that inverse. A degree-9 curve with its first
  interior knot at :math:`10^{-8}` and poles inside a :math:`1 \times 2` box
  moved by :math:`3.8 \cdot 10^7` under the unchecked removals (the audit of
  2026-10-28, A1). Each removal is therefore kept only while the running
  bound on the curve's movement, the one ``RemoveKnot`` reports, stays
  within :math:`2^{-40} (1 + \max |P|)` and the weights stay valid; where it
  would not, the knot keeps its remaining copies. The result is as small as
  the arithmetic allows and within that bound of the curve everywhere; with
  ordinary knot spacings the removals succeed and it is minimal. Every
  interior multiplicity grows by the degree increase; ``Continuity`` is
  unchanged.
- **Knot removal** (A5.8). Inserting a knot is a convex combination of the
  poles; removing one inverts it. ``RemoveKnot`` runs the inversion from
  both ends of the :math:`p - s + 1` poles the knot touches towards the
  middle, where the two sides must agree; how far they disagree is the
  distance by which one pole of the re-inserted curve would miss the
  original, and the curve moves by at most that much (for a rational
  curve, scaled by :math:`(1 + \max |P|) / \min w`, since
  :math:`\mathbf C = \mathbf A / w`). The removal is kept while the running
  sum of those bounds stays within the caller's tolerance; the result
  reports how many copies went and the bound ``Deviation``, so that a
  caller can log the change as :doc:`tolerances` asks. A tolerance of 0
  removes only what rounding allows, which is how the degree elevation
  above uses it internally.
- **Concatenation.** ``Concatenated`` clamps both curves, elevates the one
  of lower degree, shifts the second's knots so that its domain starts
  where the first's ends, scales the second's weights so that the two
  junction poles carry the same weight, merges those poles at their
  midpoint (the two ends must lie within the tolerance) and joins the knot
  vectors with the junction at multiplicity :math:`p`, so the parametrisation
  of each piece is kept and the join is :math:`C^0`. Where the pieces are
  in fact smoother, ``RemoveKnot`` at the junction with a tolerance lets
  the extra copies go; splitting a curve and concatenating the halves,
  then removing the junction knot, gives the original back. Two joins have
  no valid curve, and return ``NotSupported`` instead of one (the audit of
  2026-10-28, A3 and A4): a second curve whose knot intervals round away
  when shifted onto the first's knot values, since its parametrisation,
  not a stretched one, is the contract; and junction weights whose ratio
  takes the scaled weights outside :math:`[2^{-256}, 2^{256}]`. Before that
  check the joined weights are rescaled by the common power of two that
  centres their exponents in the range, which changes no point, so a join
  like :math:`2^{256}` against :math:`2^{-256}` succeeds.

Every rearrangement that returns a curve publishes it through the checks
``MakeBSplineCurve`` makes of data from outside the kernel: knots finite,
non-decreasing and of finite width with the multiplicities the degree
allows, weights in range, poles finite. An operation on a valid curve can
still produce knots or weights no valid curve has, as the cases above show,
and the check turns that into ``NotSupported`` rather than a successful
invalid curve. ``tests/curve/BSplineOperationsAuditTests.cpp`` holds the
five audit fixtures and a property over every degree with knots as close
as :math:`10^{-12}`.

Two queries read the structure rather than change it. ``Continuity`` is
the order of continuity the knots guarantee, :math:`p - s` for the largest
interior multiplicity :math:`s` (the degree for a single span), at one
parameter or over the whole curve; OCCT's ``IsCN(n)`` is
``Continuity(c) >= n``. ``IsClosed`` compares the points at the two ends
within a tolerance.

Batch and grid evaluation
-------------------------

``EvalAll(curve, parameters, points)`` evaluates a curve at many
parameters, and ``EvalGrid(surface, u, v, points)`` a surface at every pair
of a list of :math:`u` and a list of :math:`v`, row by row along
:math:`u`. Both return, bit for bit, the points ``Eval`` returns one at a
time. That is the whole contract: no tolerance, nothing that depends on the
SIMD width of the machine, the order of the parameters, or a thread count,
so a batch can replace a loop anywhere without changing a result
(ADR-0007).

They get there by doing the same arithmetic in fewer and wider steps.

- **Lanes.** A vector of :math:`W` doubles (``Simd::Width<double>``: 2 on
  the baseline x86-64 and arm64 builds, 4 with AVX2) holds :math:`W`
  consecutive parameters. When they share a knot span, the knots of the
  basis recurrence (A2.2) and the poles are the same for every lane and are
  broadcast, while each lane runs the recurrence and the weighted sum of
  ``Eval`` on its own parameter, in ``Eval``'s order. IEEE arithmetic gives
  the same result in a lane as in a scalar, and the build does not contract
  a multiply and an add into one rounding (``-ffp-contract=off``), so every
  lane reproduces ``Eval``. A run that straddles a span, and the last
  :math:`n \bmod W` parameters, go through ``Eval`` itself. Sorted
  parameters share spans and so take the vector path almost always.
- **Shared bases.** On a grid the span and the :math:`q + 1` basis functions
  of each :math:`v` depend on :math:`v` alone, and those of each :math:`u`
  on :math:`u` alone. ``EvalGrid`` computes them once per column and once
  per row, :math:`O(p^2 n_u + q^2 n_v)` in all instead of
  :math:`O((p^2 + q^2)\, n_u n_v)`, and leaves only the
  :math:`(p + 1)(q + 1)` weighted sum per point, run in lanes along
  :math:`v` with each column's basis functions stored contiguously. The
  sum keeps ``Eval``'s order, the product :math:`N_r(u) N_s(v)` formed
  first and then weighted, which is what makes the points identical.

Faster schemes exist that give up the identity; ``EvalAllFast`` below is
one. ``tests/curve/BatchEvaluationTests.cpp`` and
``tests/surface/GridEvaluationTests.cpp`` check the identity with ``==`` on
random curves and surfaces of degrees 1 to 7, polynomial and rational,
clamped and not, over sorted, shuffled and degenerate parameter lists.

Fast batch evaluation
---------------------

``EvalAllFast(curve, parameters, points)`` trades the identity for speed,
as OCCT's grid evaluators do. ``ToSpanPolynomials(curve)`` converts every
non-empty knot span :math:`[a, b]` once into a polynomial in the local
parameter :math:`u = (t - m) / h`, with :math:`m` the span's midpoint and
:math:`h = b - a`, so :math:`|u| \le 1/2`. Its coefficients are the
homogeneous curve's derivatives at :math:`m`,

.. math::

   (w\mathbf C, w)(t) = \sum_{k=0}^{p} \mathbf c_k\, u^k, \qquad
   \mathbf c_k = \frac{h^k}{k!}\, (w\mathbf C, w)^{(k)}(m),

computed by de Boor differencing (A3.3) in :math:`u` itself: the active
poles are differenced level by level with the factor
:math:`\frac{p - k}{k + 1}\, \frac{h}{\Delta}`, :math:`\Delta` the knot
difference of the level, which contains the span and so is at least
:math:`h`. Every factor is therefore at most :math:`p`, and the
coefficients of a span narrow against its poles come out finite, where
the derivatives in :math:`t` would overflow while :math:`h^k / k!`
underflows and their product is NaN (a defect of the release audit of
2026-10-30; ``tests/curve/FastBatchEvaluationTests.cpp`` keeps the span
:math:`h = 10^{-200}` that showed it, and the bound below on random
curves whose knots are scaled by :math:`2^{e}`, :math:`|e| \le 1000`,
or moved as far as :math:`2^{40}` from the origin). A point is then
:math:`p` multiply-adds per coordinate by Horner's rule, and one division
for a rational curve, against the :math:`O(p^2)` basis recurrence of
``Eval``.
Parameters that share a span run in the lanes of one vector, by the same
steps as the scalar fallback, so the result still does not depend on the
SIMD width, the parameters' order or how they fall into spans. Degrees 1
to 7 each get a kernel of their own, for polynomial and rational curves, so
that the compiler unrolls Horner's loop. An expansion statement
(``template for``, C++26) generates the kernels, and higher degrees share
one with a run-time degree. All kernels perform the same operations, so the
choice changes no result. A caller that
evaluates one curve many times keeps the ``SpanPolynomials`` and calls
``EvalAllFast(polynomials, parameters, points)``. The curve overload instead
converts a span when a parameter first lands in it, into a buffer on the
stack, as OCCT's evaluator rebuilds its span cache. Sorted parameters convert
each span once. Shuffled ones may convert a span again, which costs time
but not bits, since both overloads convert with the same function.

**How close it is to Eval.** With :math:`\varepsilon = 2^{-52}` and the
scale :math:`S = \max_i \lVert P_i \rVert \cdot \max_i w_i / \min_i w_i`
(the weight ratio is 1 for a polynomial curve), the tests hold

.. math::

   \lVert \texttt{EvalAllFast} - \texttt{Eval} \rVert \le
   \bigl(8(p + 1) + 2^{p-1}\bigr)\, \varepsilon S

for every degree :math:`p` from 1 to 25. This bound is tested, not proved.
``tests/curve/FastBatchEvaluationTests.cpp`` checks it on random clamped and
unclamped, polynomial and rational curves with repeated knots, at sorted
parameters including every knot and both ends. Its ``[measure]`` case,
2000 curves per degree up to 12 and 300 above, each at about 200
parameters, gives these maxima in units of :math:`\varepsilon S` (measured
2026-10-31, with the coefficients computed in the span's coordinate; the
earlier computation gave the same orders of magnitude):

=========== ============ ========= ==========
Degree      Polynomial   Rational  Bound
=========== ============ ========= ==========
1           1.3          1.8       17
2           2.1          1.7       26
3           2.7          2.3       36
5           5.4          3.0       64
7           13           9.8       128
9           48           39        336
12          331          187       2 152
16          3 396        800       32 904
20          36 978       18 238    524 456
25          1.1 × 10⁶    3.1 × 10⁵ 1.7 × 10⁷
=========== ============ ========= ==========

Up to degree 5 the deviation grows like :math:`p + 1`, as Horner's rule
alone would give. Beyond that it roughly doubles with each degree. That
growth is the conditioning of the power basis on a span: the coefficients
:math:`\mathbf c_k` may be up to about :math:`2^p` times the poles and
cancel in the sum. For a cubic with poles within 100 mm, the deviation is
at most :math:`3.1\, \varepsilon S \approx 7 \cdot 10^{-14}` mm, far below
any modelling resolution. At degree 25 it reaches about
:math:`3 \cdot 10^{-10}` relative to :math:`S`. Use ``EvalAll`` where the
points must equal ``Eval``'s, such as results compared bit for bit, or
results that feed a robust predicate. ``EvalAllFast`` suits display,
sampling and machining-style batches at ordinary degrees.

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

**Preparation.** The point enters the stationarity polynomial only through
:math:`g`, and

.. math::

   g = \sum_d A_d^2 - 2 \sum_d p_d\, (A_d w) + \norm{\mathbf p}^2 w^2

is linear in the four polynomials :math:`\sum_d A_d^2`, :math:`A_d w` and
:math:`w^2` with the coordinates of :math:`\mathbf p` as weights. The map
:math:`L(X) = X' w - 2 X w'` that turns :math:`g` into :math:`h` is linear
too, so :math:`h` is the same combination of :math:`L` of those four,
which do not depend on the point. A ``CurveProjector`` extracts the Bezier
segments and computes them once; a query only forms the combination, a few
multiply-adds per coefficient. Build one to project many points onto the
same curve; ``Project(curve, p)`` builds one for its single point.
``ProjectAll(projector, points, context)`` projects a batch on the
``Context``'s scheduler and thread count, each point's answer in its own
slot, so it returns the same bits serially and on any number of threads
(ADR-0007; :doc:`testing`, *Determinism*); ``ProjectAll`` on a
``SurfaceProjector`` does the same for surfaces.

**Assembly.** ``Project`` extracts the Bezier segments, skips any whose
pole box is farther from :math:`\mathbf p` than the best candidate so far
(the convex hull property makes that a valid lower bound), and evaluates
the actual distance at every isolated root and at the segment ends. The
segments are visited in increasing order of their bound, so the closest
ones are searched first and the search stops at the first segment whose
bound exceeds the best distance found. The smallest distance wins; an exact
tie goes to the smaller parameter, which keeps the result deterministic
(ADR-0007).

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
   Jacobian is the Hessian of half the squared distance, from their centre.
   The iteration is safeguarded: where the Hessian is not positive definite
   it uses the first fundamental form instead (Gauss–Newton), every step is
   halved until it does not increase the distance, and a coordinate on the
   patch boundary that the step would push outward stays there while the
   other moves. It is kept inside the Bezier patch, not the leaf: a step
   clamped to a leaf edge can settle there, pinned at a point that is no
   stationary point, which is how a property test once found a surface point
   projecting 3.2e-3 mm away from itself. A run that stalls, stopping
   strictly inside the patch at no stationary point, subdivides its leaf up
   to two more levels (``ExtraLeafDepth``) and restarts from the quarters'
   centres; a run that ends on the patch boundary does not, since the
   boundary's minima come from the certified edge projections.

The boundary part is certified like the curve projection. The interior part
is certified in what it excludes: a sub-patch is discarded only when it
provably contains no stationary point or no closer point. It is not certified
inside the surviving leaves, where Newton's method converges to the
stationary point when the leaf lies in its basin of attraction; a leaf is a
small fraction of a Bezier patch, but no bound on the basin is computed. The
tests require the projected point to be no farther than the closest of
6 561 grid samples on random surfaces, and every surface point to project
to itself within 1e-8 mm; the OCCT differential tests require it to be no
farther than ``GeomAPI_ProjectPointOnSurf``'s on 1 029 probes. The one
counterexample found so far, the pinned Newton run above, is a regression
test.

**Preparation and order.** As for curves, the point enters :math:`F` and
:math:`G` linearly,
:math:`F = \sum_d A_d H_d - \sum_d p_d (w H_d)` with
:math:`H = A_u w - A w_u`, so a ``SurfaceProjector`` computes the constant
and the three linear parts of both polynomials once per patch, with
projectors for its four edges. Patches, and within a patch its four
quarters at every level, are searched in increasing order of their distance
bound, and the search stops at the first that cannot hold a closer point.
``Project(surface, p)`` prepares only the patches its search reaches. Each
Newton step takes the point and all five partials it needs from one
``TotalDerivatives(surface, uv, 2, out)`` call, which computes only the
partials of total order at most 2.

``bench/reports/phase1-2026-10-15.md`` measures the evaluators, projectors
and knot operations against OCCT 8.0.1 (benchmarks B1.1–B1.5): evaluation is
2.5–3 times faster, evaluation with derivatives 1.1–1.6 times, curve
projection about 4 times, and surface projection, interpolation and Bezier
decomposition are at parity.

One case is expensive rather than wrong: when many points of a patch are at
almost the same distance, as from the centre of a spherical patch, the
coefficients of :math:`F` and :math:`G` are rounding noise, nothing is
excluded, and all :math:`4^6` leaves start Newton's method. Every candidate
is then equally good, and the result is still correct.

``MakeBSplineSurface`` validates each knot vector as ``MakeBSplineCurve``
does, requires a full grid of finite poles, and the weights to be absent or
one positive finite weight per pole.

Exact conversion of the analytic types
--------------------------------------

Every built-in curve and surface, restricted to a finite parameter range,
*is* a NURBS curve or surface: ``ToBSpline(geometry, range)`` returns it
exactly, not approximately. This is the tier-3 capability of ADR-0003 for
the built-in types: an operation that only knows NURBS (a generic
intersection, a tessellator, a STEP writer) can take any built-in geometry.
The result's domain is the range itself, it runs the same way, and it agrees
with the original at every knot; whether it agrees between knots depends on
the type.

**Lines and parabolas** are polynomial in their parameter, so their
B-splines keep the parametrisation exactly. A line over :math:`[a, b]` is
the degree-1 curve through :math:`C(a)` and :math:`C(b)`; a parabola is the
quadratic Bezier curve with poles :math:`C(a)`,
:math:`C(a) + \tfrac{b-a}{2} C'(a)` and :math:`C(b)`, the middle pole being
where the end tangents meet.

**Circular arcs** (algorithm A7.1). The rational quadratic Bezier curve
with poles :math:`(\cos a, \sin a)`,
:math:`(\cos m, \sin m) / \cos\delta` and :math:`(\cos b, \sin b)`, middle
weight :math:`\cos\delta`, :math:`m = (a + b)/2` and
:math:`\delta = (b - a)/2`, is exactly the unit circle's arc from angle
:math:`a` to :math:`b`: its middle pole is where the end tangents meet, and
at its parameter :math:`1/2` it passes through

.. math::

   \frac{\tfrac14 (\cos a, \sin a) + \tfrac12 \cos\delta \frac{(\cos m, \sin m)}{\cos\delta}
         + \tfrac14 (\cos b, \sin b)}{\tfrac12 + \tfrac12 \cos\delta}
   = (\cos m, \sin m),

using :math:`\cos a + \cos b = 2 \cos m \cos\delta` and the same for
:math:`\sin`. A conic is fixed by two points, their tangents and one more
point, so the curve is the circle. An arc wider than a quarter turn is split
into equal pieces of at most a quarter turn (``Detail::MakeUnitArc``), which
keeps every middle weight at least :math:`\cos(\pi/4)` and every middle pole
near the circle; a full turn is four pieces, nine poles, and its last pole
is set to its first so that it closes exactly. The knots are the piece end
angles, so the arc passes through the circle's point of angle :math:`t` at
every knot :math:`t`. Between knots its parameter is not the angle: no
rational polynomial is proportional to arc length on a circle.

**Ellipses** are the image of the unit circle under the affine map
:math:`(x, y) \mapsto \mathbf c + a x \mathbf x + b y \mathbf y`, and NURBS
are affinely invariant, so an elliptic arc is the circular arc with mapped
poles and the same weights.

**Hyperbolic arcs** follow the hyperbolic analogue. With
:math:`m = (t_0 + t_1)/2` and :math:`h = (t_1 - t_0)/2`, the rational
quadratic Bezier curve with poles :math:`H(t_0)`,
:math:`\mathbf c + (a\cosh m\, \mathbf x + b\sinh m\, \mathbf y) / \cosh h`
and :math:`H(t_1)`, middle weight :math:`\cosh h`, passes through
:math:`H(m)` at its parameter :math:`1/2`, by
:math:`\cosh t_0 + \cosh t_1 = 2 \cosh m \cosh h` and
:math:`\sinh t_0 + \sinh t_1 = 2 \sinh m \cosh h`. The range is split into
pieces of parameter length at most 1, so that :math:`\cosh h` stays near 1.

**Surfaces of revolution** (algorithm A8.1). The cylinder, cone, sphere and
torus are all :math:`\mathbf o + \rho(v)\, \mathbf e(u) + h(v)\, \mathbf z`
(:doc:`analytic`), a meridian :math:`(\rho, h)` swept around the axis. If
the meridian is a NURBS curve with poles :math:`(\rho_j, h_j)` and weights
:math:`w_j`, and the unit arc in :math:`u` has poles :math:`U_i` and weights
:math:`\omega_i`, the tensor-product surface with poles
:math:`\mathbf o + \rho_j U_i + h_j \mathbf z` and weights
:math:`\omega_i w_j` is

.. math::

   \mathbf o
   + \frac{\sum_j M_j w_j \rho_j}{\sum_j M_j w_j}\,
     \frac{\sum_i N_i \omega_i U_i}{\sum_i N_i \omega_i}
   + \frac{\sum_j M_j w_j h_j}{\sum_j M_j w_j}\, \mathbf z
   = \mathbf o + \rho(v)\, \mathbf e(u) + h(v)\, \mathbf z,

because the sums separate. The meridian is a line for the cylinder and the
cone, and a circular arc for the sphere (about the axis) and the torus
(about the centre of the tube). It holds for :math:`\rho < 0` as well, so a
cone range may cross its apex.

**Planes** over a box are the bilinear patch through the four corners, with
their own parametrisation.

**B-splines** convert to their restriction to the range: a curve is split at
the range ends, a surface cut in each direction by raising the cut values
to the degree and dropping the poles beyond them. A range equal to the
domain returns the geometry unchanged.

``ToBSpline`` validates the range as outside data: ``InvalidInput`` unless
it is finite, non-empty and inside the domain, or in a periodic parameter at
most one period long anywhere. Through the wrappers ``Curve`` and
``Surface`` it dispatches with ``VisitBuiltin``; a user type gets
``NotSupported``, since the kernel knows no exact form for it (ADR-0003
expects such types to supply their own ``ToBSpline`` with a certified
deviation).

The tests check every type over random ranges. Every point of the
conversion must lie on the original geometry inside the range, the two must
agree on every knot (and knot pair), tangents and normals must point the same
way, and every point of the original's range must project onto the
conversion at distance zero, which also exercises the curve and surface
projections on rational input.

Checked construction
--------------------

``MakeBSplineCurve(degree, knots, poles, weights)`` validates data from
outside the kernel (ADR-0006, ADR-0016) and returns ``InvalidInput`` with a
message unless

- :math:`1 \le p \le 25` and there are at least :math:`p + 1` poles;
- there are exactly :math:`n + p + 2` knots, finite and non-decreasing, with
  a non-empty domain, and the last minus the first is finite: the basis
  functions divide by differences of knots, and :math:`10^{308} -
  (-10^{308})` is not a double;
- no knot value is repeated more than :math:`p + 1` times, and none inside
  the domain more than :math:`p` times;
- the poles are finite, and the weights are absent or one positive weight
  per pole in :math:`[2^{-256}, 2^{256}]` (``Detail::MinWeight`` and
  ``Detail::MaxWeight``).

The weight range is the numerical domain evaluation can serve. Weights of
:math:`10^{308}` are finite, but the weighted poles :math:`w_i P_i`
overflow even for a curve well inside the model; within the range neither
they nor the denominator overflow or lose precision to underflow. Inside it
a rational curve is invariant under a common scale of its weights, and a
scale by a power of two changes no computed point or derivative, bit for
bit. ``MakeBSplineSurface`` applies the same rules to its weights.

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

- B-spline *surfaces* have no periodic form: no fit produces one and nothing
  recognises one, so a surface is evaluated inside its domain only. Curves
  in the periodic form are periodic (*Periodic curves*, above).
- Surface projection is not certified inside its final sub-patches (see
  *Surface projection*).
- Degree elevation, knot removal, clamping and concatenation exist for
  curves only; surfaces have knot insertion alone.
- Batch and grid evaluation give points only; derivatives are evaluated a
  point at a time.
- ``EvalAllFast`` is for curves only; grids keep ``EvalGrid``, which is
  bit-identical. Its accuracy bound is tested, not proved.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
