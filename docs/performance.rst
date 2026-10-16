.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Performance against OpenCASCADE
===============================

Stereon's geometry is measured against OpenCASCADE (OCCT), operation by
operation, on the same geometry, the same machine and the same compiler.
This chapter explains how the comparison is made, what each operation costs
in theory and why, and what the measurements show. The catalogue of
benchmarks and their targets is ``docs/benchmarks.md``; the dated reports,
with the exact versions and machine, are in ``bench/reports/``.

A speed result only counts when both libraries compute the same thing
correctly. The differential tests (``tests/occt/GeometryDifferentialTests.cpp``)
check that values and derivatives agree with OCCT to a relative
:math:`10^{-13}`, and projections are never farther than OCCT's; the arc
lengths below agree to :math:`2 \cdot 10^{-14}` relative.

Method
------

The harness is ``bench/occt/OcctBench.cpp`` (target ``stereon_occt_bench``,
built when ``STEREON_WITH_OCCT`` is on), written with Google Benchmark.
Every Stereon benchmark has an OCCT twin with the same name after the last
slash. Both libraries are release builds (``-O3 -DNDEBUG``) of GCC 16.2;
Stereon's contracts are compiled out, as in an application's release
build; OCCT 8.0.1 is built by ``tools/occt/build_occt.sh`` from a pinned
commit.

Each benchmark processes a batch of :math:`n` items per iteration: 4096
parameters for evaluation and local properties, 256 query points for
projection, 64 lengths for the point-at-length search, one fit, one
decomposition. Google Benchmark repeats the batch until the timing is
stable, giving a CPU time :math:`T_r` per batch in repetition
:math:`r = 1, \dots, R` with :math:`R = 5`. The reported time per item is
the median over repetitions,

.. math::

   t = \operatorname{median}_{r}\, \frac{T_r}{n},

and the comparison is the speed ratio

.. math::

   \rho = \frac{t_{\mathrm{OCCT}}}{t_{\mathrm{Stereon}}},

so :math:`\rho > 1` means Stereon is faster and :math:`\rho = 2` means it
takes half the time. The median rather than the mean keeps a stray slow
repetition, from the operating system or the CPU's frequency scaling, out
of the result. On the laptop used here repeated runs move a ratio by up
to about 10 %, so a ratio between 0.9 and 1.1 is reported as parity.

Each library gets its best usage where it has one. Projections reuse an
initialised projector, ``GeomAPI_ProjectPointOnCurve::Perform`` against a
reused ``CurveProjector``. Evaluation with derivatives is one ``D2`` call
against one ``Derivatives`` call. Curve properties are one
``GeomLProp_CLProps::SetParameter`` against one ``FrenetFrame``. Bezier
decomposition asks OCCT for every arc, because OCCT builds them lazily.

Cost of each operation
----------------------

The costs below count floating-point operations per item for a curve of
degree :math:`p` with :math:`n + 1` poles in :math:`d` dimensions, or a
surface of degrees :math:`p, q`. They explain the measured ratios: where
both libraries do the same arithmetic the ratio comes from memory traffic
and overheads, and where one does less arithmetic, it shows.

Evaluation (B1.1)
~~~~~~~~~~~~~~~~~

The point of a B-spline on its span :math:`s`,
:math:`\mathbf C(t) = \sum_{j=0}^{p} N_{s-p+j,p}(t)\, \mathbf P_{s-p+j}`,
needs the :math:`p + 1` non-zero basis functions, which the triangular
Cox–de Boor recurrence produces in about :math:`p^2` multiply–adds (The NURBS
Book, A2.2), and a combination of :math:`p + 1` poles,
:math:`(p + 1)\, d` more. Finding the span is a binary search,
:math:`O(\log n)`. A rational curve combines homogeneous poles
:math:`(w_i \mathbf P_i, w_i)` and divides once. A surface needs two basis
triangles and a :math:`(p + 1) \times (q + 1)` combination:

.. math::

   \mathrm{cost}_{D0} \approx p^2 + (p + 1)\, d
   \quad\text{(curve)}, \qquad
   p^2 + q^2 + (p + 1)(q + 1)\, d
   \quad\text{(surface)}.

The arithmetic is the same in both libraries; Stereon's lead of
:math:`\rho \approx 2.5` comes from what surrounds it: poles in one
contiguous array of plain structs, work arrays on the stack sized to the
degree, and no handle or adaptor indirection.

Derivatives use the fact that the derivative of a B-spline is a B-spline
of degree :math:`p - 1` on the same knots, with differenced poles
:math:`\mathbf Q_i = p\, (\mathbf P_{i+1} - \mathbf P_i) / (u_{i+p+1} -
u_{i+1})`. The basis triangle already holds the functions of every lower
degree, so all derivatives up to order :math:`k \le p` cost

.. math::

   \mathrm{cost}_{Dk} \approx p^2 + k\,(p + 1)\, d + (k + 1)(p + 1)\, d

after one span search: one triangle, :math:`k` rounds of pole
differences, and one combination per order (:doc:`nurbs`). A rational
curve then applies the quotient rule (A4.2), :math:`O(k^2 d)` more. The
lead over OCCT is smaller here, :math:`\rho \approx 1.4`–:math:`1.6`, as
the arithmetic both libraries must do is a larger share of the work.

Projection (B1.3)
~~~~~~~~~~~~~~~~~

The closest point of a curve to :math:`\mathbf x` is a root of the
stationarity function :math:`(\mathbf C(t) - \mathbf x) \cdot \mathbf
C'(t)`. On each Bezier segment of a rational curve, with homogeneous
numerator :math:`\mathbf A` and weight :math:`w`, it becomes the
polynomial

.. math::

   h(t) = \sum_d \bigl(A_d' w - A_d w'\bigr)\bigl(A_d - x_d w\bigr),

whose roots are isolated by Bernstein subdivision with Descartes' rule of
signs. The point enters linearly, so the products that do not involve it
are computed once per curve by a ``CurveProjector``; a query forms
:math:`h` from :math:`d + 1` precomputed polynomials in
:math:`O(d\, m)` operations for a segment polynomial of degree :math:`m`,
then isolates its roots. Segments are visited best-first by a lower bound
on their distance, and a segment whose bound exceeds the best distance so
far is never expanded. OCCT's ``Extrema`` package instead samples the
curve and refines candidates with Newton's method, evaluating the curve
many times per query; the measured ratio is :math:`\rho \approx 4`.

On surfaces both libraries are dominated by the subdivision search, and
they are at parity.

Fitting (B1.4)
~~~~~~~~~~~~~~

Interpolation solves the collocation system :math:`N \mathbf P = \mathbf
Q`, whose matrix :math:`N_{ki} = N_{i,p}(t_k)` is banded with bandwidth
:math:`p` and totally positive, so Gaussian elimination without pivoting
is stable and costs :math:`O(n\, p^2)`. Both libraries do this, at parity.
Least squares forms the normal equations :math:`(N^\mathsf{T} N)\mathbf P =
N^\mathsf{T}\mathbf Q`, banded too. ``Approximate`` searches for the
fewest poles that meet an ``ApproximationBudget``; its ratio of about 40 is
**not** like for like, since OCCT's ``GeomAPI_PointsToBSpline`` also
optimises parameters and continuity (:doc:`fitting`).

Knot operations (B1.5)
~~~~~~~~~~~~~~~~~~~~~~

Inserting a knot once changes :math:`p` poles, :math:`O(p\, d)` work (A5.1),
plus moving the array tail. Bezier decomposition raises every interior
knot to multiplicity :math:`p`; Stereon walks the knots once and updates
only the current segment's poles (A5.6), :math:`O(n\, p^2 d)` in total,
rather than inserting knot by knot into the whole array, which would move
:math:`O(n)` poles each time.

Local properties (B1.6)
~~~~~~~~~~~~~~~~~~~~~~~

With :math:`\mathbf C'` and :math:`\mathbf C''` from one evaluation, the
Frenet frame and the curvature (:doc:`properties`)

.. math::

   \mathbf T = \frac{\mathbf C'}{\norm{\mathbf C'}}, \quad
   \mathbf B = \frac{\mathbf C' \times \mathbf C''}{\norm{\mathbf C' \times \mathbf C''}}, \quad
   \mathbf N = \mathbf B \times \mathbf T, \quad
   \kappa = \frac{\norm{\mathbf C' \times \mathbf C''}}{\norm{\mathbf C'}^3}

cost two cross products, two square roots and a few divisions on top of
:math:`\mathrm{cost}_{D2}`. ``FrenetFrame`` is that one evaluation, as is
one ``SetParameter`` of OCCT's ``GeomLProp_CLProps``; the ratio,
:math:`\rho \approx 1.8`, is close to the one for D0–D2 evaluation, which
dominates both. Every property function evaluates the curve once, through
the type's combined ``Derivatives`` where it has one; a first version
called ``Deriv`` once per order and was twice as slow as OCCT.

On surfaces the principal curvatures need
:math:`\mathbf S_u, \mathbf S_v, \mathbf S_{uu}, \mathbf S_{uv}, \mathbf
S_{vv}`, the fundamental forms (six dot products), and the roots
:math:`H \pm \sqrt{H^2 - K}`. Stereon obtains the partials from one call of
``Derivatives(surface, uv, 2, 2, out)``, which also computes
:math:`\mathbf S_{uuv}, \mathbf S_{uvv}, \mathbf S_{uuvv}`, three partials
the curvature does not use. That surplus is why this is the one row where
OCCT is ahead (:math:`\rho \approx 0.8`); a triangular evaluator returning
only orders :math:`k + l \le 2` would remove it.

Arc length (B1.6)
~~~~~~~~~~~~~~~~~

``ArcLength`` integrates the speed with 16-point Gauss–Legendre rules on
pieces that are halved until the one-piece value :math:`Q_1` and the
two-piece value :math:`Q_2` of each agree,

.. math::

   \abs{Q_2 - Q_1} \le 2^{-40}\, s_{\text{so far}},

after splitting at the knots, where the speed has kinks. Its cost is
:math:`16 \times (\text{pieces evaluated}) \times \mathrm{cost}_{D1}`; on
this curve the result agrees with OCCT's ``GCPnts_AbscissaPoint::Length``
at tolerance :math:`10^{-12}` to :math:`2 \cdot 10^{-14}` relative, slightly
faster (:math:`\rho \approx 1.1`).

``ParameterAtArcLength`` solves :math:`s(t) = \ell` by Newton's method,
:math:`t_{k+1} = t_k - (s(t_k) - \ell) / \norm{\mathbf C'(t_k)}`. The
length to each iterate is accumulated, :math:`s(t_{k+1}) = s(t_k) +
s(t_k, t_{k+1})`, so each step integrates only the short stretch between
two iterates, whose length shrinks quadratically as Newton converges:

.. math::

   \mathrm{cost} \approx \mathrm{cost}\bigl(s(t_0, t_1)\bigr)
   + \sum_{k \ge 1} \mathrm{cost}\bigl(s(t_k, t_{k+1})\bigr)
   \approx \mathrm{cost}(s(\text{start}, t_1)) + O(\text{iterations}).

Integrating from the start at every step, as a first version did, cost
more than twice as much and was slower than OCCT; accumulated, it is
:math:`\rho \approx 1.5`.

Results
-------

Median CPU time per item, GCC 16.2 ``-O3``, OCCT 8.0.1, Intel Core
i5-10210U (Windows 11 laptop), single thread, 2026-10-16.

.. list-table::
   :header-rows: 1
   :widths: 8 40 14 14 10

   * - ID
     - Benchmark
     - Stereon
     - OCCT
     - :math:`\rho`
   * - B1.1
     - Curve D0, cubic, 24 poles
     - 110 ns
     - 291 ns
     - 2.6
   * - B1.1
     - Curve D0, rational cubic
     - 126 ns
     - 360 ns
     - 2.9
   * - B1.1
     - Curve D0, degree 7
     - 238 ns
     - 627 ns
     - 2.6
   * - B1.1
     - Curve D0–D2, cubic
     - 266 ns
     - 383 ns
     - 1.4
   * - B1.1
     - Curve D0–D2, rational cubic
     - 327 ns
     - 523 ns
     - 1.6
   * - B1.1
     - Surface D0, bicubic 12 × 12
     - 286 ns
     - 665 ns
     - 2.3
   * - B1.1
     - Surface D0, rational bicubic
     - 293 ns
     - 734 ns
     - 2.5
   * - B1.1
     - Surface D0–D1, bicubic
     - 959 ns
     - 991 ns
     - parity
   * - B1.3
     - Curve projection, rational cubic
     - 18.8 µs
     - 76.3 µs
     - 4.1
   * - B1.3
     - Surface projection, bicubic
     - 537 µs
     - 519 µs
     - parity
   * - B1.4
     - Interpolation, 1000 points, cubic
     - 0.49 ms
     - 0.48 ms
     - parity
   * - B1.4
     - Approximation, 1000 points, 1e-3 mm
     - 6.0 ms
     - 234 ms
     - 39 (not like for like)
   * - B1.5
     - 100 knot insertions
     - 6.0 µs
     - 11.2 µs
     - 1.9
   * - B1.5
     - Bezier decomposition, 21 segments
     - 36.8 µs
     - 44.3 µs
     - 1.2
   * - B1.6
     - Frenet frame and curvature, rational cubic
     - 375 ns
     - 692 ns
     - 1.8
   * - B1.6
     - Principal curvatures and directions, bicubic
     - 1.88 µs
     - 1.49 µs
     - 0.79
   * - B1.6
     - Arc length of the whole curve, rational cubic
     - 356 µs
     - 399 µs
     - 1.1
   * - B1.6
     - Parameter at a given length
     - 617 µs
     - 922 µs
     - 1.5

Not yet measured: the SIMD batch evaluation of B1.1, grid evaluation
(B1.2), B1.3's convergence rate over a corpus, and degree elevation (B1.5),
which Stereon does not have yet. All numbers come from one laptop; the
release report for v0.1 needs the CI reference machine.

Reproducing the numbers
-----------------------

With OCCT installed by ``tools/occt/build_occt.sh``:

.. code-block:: console

   cmake -S . -B build/bench-occt -G Ninja -DCMAKE_BUILD_TYPE=Release \
         -DSTEREON_BUILD_BENCH=ON -DSTEREON_WITH_OCCT=ON \
         -DSTEREON_CONTRACT_SEMANTIC=ignore -DSTEREON_HOT_CONTRACTS=OFF \
         -DCMAKE_PREFIX_PATH=<occt install prefix>
   cmake --build build/bench-occt --target stereon_occt_bench
   build/bench-occt/bench/stereon_occt_bench --benchmark_repetitions=5 \
         --benchmark_report_aggregates_only=true

``--benchmark_filter=B1.6`` runs one group. The ``length_mm`` counter of
the arc-length rows prints each library's result, so the accuracy check
travels with the timing.
