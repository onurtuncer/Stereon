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

Two kinds of run are told apart. A *development run* is made on whatever
machine is to hand, as the table below was, and may be quoted in a pull
request as long as its machine is named. The *reference run* is made by CI,
``.github/workflows/occt-reference.yml``: the same ``gcc:16.2`` container
and pinned OCCT build as the differential tests, on a GitHub-hosted
``ubuntu-24.04`` runner, with 15 repetitions of the two libraries
interleaved. Each run uploads the raw Google Benchmark JSON, a record of the
machine (CPU model, cores, memory, kernel, compiler, OCCT tag and commit,
Stereon commit) and the report that ``tools/bench_report.py`` renders from
them, with the median per row and the largest deviation of any repetition
from it. The report a release gate cites (``docs/benchmarks.md``,
*Reporting*) is a reference run committed to ``bench/reports/``; a hosted
runner's CPU differs from run to run and is shared with other jobs, so a
row whose spread is large is run again rather than quoted.

A third kind, the *trend run*, follows Stereon alone from commit to commit.
Every push to main runs ``stereon_bench`` (all of it but the
:math:`10^6`-face topology rows and the :math:`10^8`-evaluation dispatch
rows) with five interleaved repetitions in the same container, and
``tools/bench_trend.py`` records each benchmark's median per item on the
``bench-data`` branch, one file per run. The job's summary lists every
benchmark more than 10 % slower than in the run before, as
``docs/benchmarks.md`` asks a pull request to explain, and the
`trend page <benchmarks/trend/index.html>`_ plots the last 60 runs. It warns
and never fails: the runner's CPU, recorded with every run, changes between
runs, and a step in a chart is read against it.

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
:math:`\rho \approx 2`–:math:`3` comes from what surrounds it: poles in one
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
lead over OCCT is smaller here, :math:`\rho \approx 1.2`–:math:`1.6`, as
the arithmetic both libraries must do is a larger share of the work.

A surface's first partials have a kernel of their own
(``TotalDerivatives(s, uv, 1, out)``, and ``Deriv`` of first order): the
basis functions and their first derivatives of both directions (A2.3), then
one pass over the :math:`(p + 1)(q + 1)` active poles with three
multiply–adds each, for :math:`\mathbf S`, :math:`\mathbf S_u` and
:math:`\mathbf S_v` together, relative to the pole with the largest basis
product for the rounding reason of the general path. That is what OCCT's
``D1`` returns, and what the *Surface D0–D1* row compares; the mixed
partial that ``Derivatives(s, uv, 1, 1, out)`` also gives is not part of
either.

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
many times per query; the measured ratio is :math:`\rho \approx 3.5`–:math:`4`.

On surfaces both libraries are dominated by the search. Two of its costs
are avoided where they prove nothing: the pole-box bound of a sub-patch is
taken without a division per coefficient when the surface is polynomial,
since :math:`w` is then exactly 1 in every coefficient and stays so under
subdivision; and a patch's edges, whose four certified curve projections
were the dearest step of its search, are projected only when the edge's own
bounding box is nearer than the best point so far. Neither changes an
answer: only candidates provably farther than the best go unconsidered. The
leaf Newton iteration stays safeguarded (see the results below).

Fitting (B1.4)
~~~~~~~~~~~~~~

Interpolation solves the collocation system :math:`N \mathbf P = \mathbf
Q`, whose matrix :math:`N_{ki} = N_{i,p}(t_k)` is banded with bandwidth
:math:`p` and totally positive, so Gaussian elimination without pivoting
is stable and costs :math:`O(n\, p^2)`. Both libraries do this. Stereon
builds the collocation rows once: the span of each parameter is found by
scanning forward from the previous row's, since the parameters increase,
the basis functions at it fill the band matrix, and after the solve the same
rows give the fit's deviation at its points, which ``MaxDeviationOf`` would
otherwise re-evaluate the curve for, point by point, with a span search and
a basis triangle each.
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
:math:`\rho \approx 1.4`–:math:`1.8` across runs, is close to the one for
D0–D2 evaluation, which dominates both. Every property function evaluates the curve once, through
the type's combined ``Derivatives`` where it has one; a first version
called ``Deriv`` once per order and was twice as slow as OCCT.

On surfaces the principal curvatures need
:math:`\mathbf S_u, \mathbf S_v, \mathbf S_{uu}, \mathbf S_{uv}, \mathbf
S_{vv}`, the fundamental forms (six dot products), and the roots
:math:`H \pm \sqrt{H^2 - K}`. Stereon obtains the partials from one call of
``TotalDerivatives(surface, uv, 2, out)``, which computes only the partials
of total order :math:`k + l \le 2`, through the surface evaluator below.

A first version called ``Derivatives(surface, uv, 2, 2, out)``, which also
computes :math:`\mathbf S_{uuv}, \mathbf S_{uvv}, \mathbf S_{uuvv}`, and
formed each of the nine partials as a full sum over the
:math:`(p + 1)(q + 1)` window: :math:`\rho = 0.79`. Measuring it stage by
stage (``bench/SurfaceDerivativeBench.cpp``) showed the basis derivatives
of algorithm A2.3 taking over half the time, most of it in divisions that
depend on one another.

**Surface derivatives without basis derivatives.** A tensor-product surface
:math:`\mathbf S(u, v) = \sum_r \sum_s N_r(u) M_s(v)\, \mathbf P_{rs}` is
differentiated one direction at a time. For each of the :math:`p + 1`
active rows :math:`r`, the row is a curve in :math:`v`; its :math:`l`-th
derivative is a B-spline of degree :math:`q - l` on the same knots whose
poles are the row's poles differenced :math:`l` times (eq. 3.8),

.. math::

   \mathbf R_{r,l}(v) = \sum_{s} M_{s,q-l}(v)\, \Delta_v^{l} \mathbf P_{rs},
   \qquad
   \Delta_v \mathbf Q_s = \frac{q - k}{v_{s+q+1} - v_{s+k+1}}
                          (\mathbf Q_{s+1} - \mathbf Q_s)

at difference level :math:`k`, and its basis functions are a level of the
triangle that Cox–de Boor builds for the point anyway. The values
:math:`\mathbf R_{r,l}` are in turn the poles in :math:`u` of
:math:`\partial^l \mathbf S / \partial v^l`, so the same step along
:math:`u` gives

.. math::

   \frac{\partial^{k+l} \mathbf S}{\partial u^k\, \partial v^l}
   = \sum_r N_{r,p-k}(u)\, \Delta_u^{k} \mathbf R_{r,l}(v).

A rational surface does the same on homogeneous poles
:math:`(w \mathbf P, w)` and applies the quotient rule (A4.4). The cost is
two basis triangles, about :math:`p^2 + q^2`, a few divisions for the
difference factors, which depend only on the knots and are computed once
per point, and a contraction of about
:math:`(p + 1)(q + 1)\, k_v + k_u k_v (p + 1)` vector multiply–adds instead
of :math:`k_u k_v (p + 1)(q + 1)`. With it, principal curvatures went from
:math:`\rho = 0.79` to :math:`\rho \approx 1.1`, and every other surface
derivative (``Deriv``, ``Derivatives``, the projection's Newton steps)
became faster too.

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
:math:`\rho \approx 1.5`–:math:`1.7`.

Results
-------

The reference run of 2026-10-19 (``bench/reports/phase1-reference-2026-10-19.md``, run
39236132501): median CPU time per item over 15 repetitions with the
two libraries' runs interleaved, GCC 16.2 ``-O3``, OCCT 8.0.1, the
``gcc:16.2`` container on a GitHub-hosted ``ubuntu-24.04`` runner, single
thread. The report lists each row's spread across the repetitions.

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
     - 51.4 ns
     - 84.9 ns
     - 1.65
   * - B1.1
     - Curve D0, rational cubic
     - 57.3 ns
     - 88.3 ns
     - 1.54
   * - B1.1
     - Curve D0, degree 7
     - 94.6 ns
     - 153 ns
     - 1.61
   * - B1.1
     - Curve D0–D2, cubic
     - 94.7 ns
     - 98.7 ns
     - 1.04
   * - B1.1
     - Curve D0–D2, rational cubic
     - 104 ns
     - 139 ns
     - 1.33
   * - B1.1
     - Surface D0, bicubic 12 × 12
     - 109 ns
     - 172 ns
     - 1.58
   * - B1.1
     - Surface D0, rational bicubic
     - 115 ns
     - 198 ns
     - 1.73
   * - B1.1
     - Surface D0–D1, bicubic
     - 191 ns
     - 255 ns
     - 1.33
   * - B1.3
     - Curve projection, rational cubic
     - 5.5 µs
     - 20.5 µs
     - 3.7
   * - B1.3
     - Surface projection, bicubic
     - 104 µs
     - 128 µs
     - 1.24
   * - B1.4
     - Approximation, 1000 points, 1e-3 mm
     - 1.95 ms
     - 36.7 ms
     - 19
   * - B1.4
     - Interpolation, 1000 points, cubic
     - 84.1 µs
     - 113 µs
     - 1.34
   * - B1.5
     - Bezier decomposition, 21 segments
     - 4.97 µs
     - 7.3 µs
     - 1.47
   * - B1.5
     - 100 knot insertions
     - 88 µs
     - 227 µs
     - 2.6
   * - B1.6
     - Arc length of the whole curve, rational cubic
     - 106 µs
     - 112 µs
     - 1.06
   * - B1.6
     - Frenet frame and curvature, rational cubic
     - 135 ns
     - 189 ns
     - 1.41
   * - B1.6
     - Parameter at a given length
     - 182 µs
     - 293 µs
     - 1.61
   * - B1.6
     - Principal curvatures and directions, bicubic
     - 442 ns
     - 409 ns
     - 0.93

No row is below parity. Surface projection, the one row where OCCT had been
ahead (0.89 on the first reference run), reads 1.24 since the search
skips the divisions and edge projections that prove nothing (*Projection*
above); the safeguards of the leaf Newton iteration (:doc:`nurbs`, *Surface
projection*), which cost about 10 % when they were added, stay, since being
farther from the right answer is not on offer at any speed. Surface first
derivatives (1.33) and interpolation (1.34) were below parity on that run
too and have the kernels described under *Evaluation* and *Fitting*.

Not yet measured: the SIMD batch evaluation of B1.1, grid evaluation
(B1.2), B1.3's convergence rate over a corpus, and degree elevation (B1.5),
which Stereon does not have yet. The reference host varies between runs by
up to about 1.5× in absolute time on both sides, so a ratio is read within
one run, and one within a few percent of 1 is parity.

Where OCCT gets it wrong
------------------------

Speed only counts with the right answer, so the comparison also hunts for
inputs on which OCCT is wrong. ``tests/occt/RobustnessHunt.cpp`` generates
adversarial geometry with a fixed seed, runs both kernels, and judges both
against a reference that shares no code with either:

- lengths of polylines (degree-1 B-splines) by exact arithmetic on their
  points;
- lengths of smooth curves by an inscribed polyline of :math:`2^k` chords,
  doubled until two values agree to :math:`10^{-12}` and extrapolated
  once, since its error falls as :math:`1/n^2`;
- distances by dense sampling, which bounds the true distance from above:
  a projection farther than the sampled bound is wrong.

OCCT is given a tolerance of :math:`10^{-10}` throughout. The hunt runs
on request (``stereon_occt_tests "[.hunt]"``) and reports what it finds;
the findings worth keeping become ordinary tests in
``tests/occt/OcctRobustnessTests.cpp``, which CI runs against the pinned
OCCT 8.0.1.

.. list-table:: Results of the hunt, OCCT 8.0.1
   :header-rows: 1
   :widths: 30 34 18 18

   * - Job
     - Adversarial input
     - OCCT wrong
     - Stereon wrong
   * - Point at a given length (``GCPnts_AbscissaPoint``)
     - Polylines with stalls (repeated points) and knot spans from
       :math:`10^{-6}` to :math:`10^{2}`
     - 10 of 400, worst 4.3 % of the length
     - 0
   * - Arc length (``GCPnts_AbscissaPoint::Length``)
     - Rational cubics with weights from :math:`10^{-3}` to :math:`10^{3}`
     - 5 of 60, worst 1.4 %
     - 0
   * - Point at a given length
     - The same rational cubics
     - 4 of 60, worst 4.5 %
     - 0
   * - Surface projection (``GeomAPI_ProjectPointOnSurf``)
     - Bicubic sheets with one narrow bump, points beside it
     - 10 of 120, worst 0.048 mm too far
     - 0
   * - Curve projection (``GeomAPI_ProjectPointOnCurve``)
     - Narrow weighted spikes; hairpins with gaps down to :math:`10^{-4}`
     - 0 of 800
     - 0
   * - Arc length
     - The polylines above
     - 0 of 400
     - 0
   * - Principal curvatures (``GeomLProp_SLProps``)
     - The exact NURBS sphere, :math:`10^{-14}` to :math:`10^{-2}` of the
       parameter range from an edge where a row of poles collapses to the
       pole
     - 114 of 300 wrong (by up to 50 % already at :math:`10^{-8}`), 116 no
       answer
     - 107 of 300 beyond :math:`10^{-6}`, all closer than :math:`10^{-10}`;
       worst 1.1 %
   * - Surface projection
     - Points within :math:`10^{-12}` to :math:`10^{-1}` of the axis of the
       NURBS sphere
     - 8 of 300, worst :math:`3.1 \cdot 10^{-7}` mm too far
     - 0
   * - Surface projection
     - The same points on the analytic sphere
     - 0 of 300
     - 0
   * - Curve and surface projection
     - Points within :math:`10^{-15}` to :math:`10^{-3}` rad of the seam of a
       circle, a closed NURBS circle, a cylinder and a torus
     - 0 of 1 200
     - 0
   * - Arc length and point at a length
     - A circle, starting just before its seam, over up to five turns
     - 0 of 600
     - 0

The tests record what OCCT returns, not why. The pattern fits the two
designs. Stereon's arc length halves every piece until the halves agree,
after splitting at every knot, so a speed that varies by six orders of
magnitude only costs it more pieces; OCCT's error control evidently does not
see such variation, since it returns a length off by 1.4 % at a requested
tolerance of 1e-10. OCCT's surface extrema begin from a sampling of the
surface and refine with Newton's method, so a bump between samples can hide
the closest point, and the refinement settles at a farther local minimum;
Stereon's Bernstein subdivision discards a region only once it has shown
that the region holds no stationary point or no closer point.

The curvature row needs a word. Near an edge where the poles collapse,
:math:`\mathbf S_u` shrinks with the distance :math:`d` to the edge while
the second fundamental form keeps its size, so a curvature computed from
this parametrisation in double precision carries an error of order
:math:`\varepsilon / d` at best. Stereon reaches that bound: within
:math:`10^{-6}` down to :math:`10^{-10}` of the range from the edge, 1.1 %
at :math:`10^{-14}`. Getting there took two fixes the hunt forced, both
recorded as regression tests:

- ``ToBSpline(sphere)`` did not collapse the pole rows exactly. Its meridian
  ended at :math:`\cos(\pi/2) = 6.1 \cdot 10^{-17}` rather than 0, so the
  converted sphere had a ring of radius :math:`1.5 \cdot 10^{-16} R` at each
  pole, whose real curvature grows like :math:`10^{-16} / d^2`: errors up to
  :math:`2 \cdot 10^{11}` at :math:`d = 10^{-14}`. Unit-circle points are
  now reduced by quarter turns before their cosine and sine are taken, so
  arcs that end on a quarter turn end exactly on it.
- The surface evaluator worked in absolute coordinates, where a difference
  of poles near the point carries an error of order :math:`\varepsilon
  \lvert \mathbf P \rvert` however small the difference is. It now works
  relative to the pole with the largest basis product at the point, so near
  a collapsed edge every term shrinks with :math:`d`.

The hunt also found a weakness in Stereon: ``ParameterAtArcLength`` used to
stop when its Newton step was small relative to the whole parameter range,
which on a curve whose knot spans range from :math:`10^{-6}` to
:math:`10^{2}` left 155 of 400 answers off by up to :math:`10^{-5}` of the
length. It now stops on the length residual, and all 400 are within
:math:`10^{-9}`. That is the point of hunting both ways.

Curve dispatch (B0.2)
---------------------

B0.2 times :math:`10^8` evaluations per measurement of a line, a circle and
ADR-0003's cubic B-spline, along three paths, each against its OCCT
counterpart in the same process:

- *Direct*: ``Eval`` on the concrete type, against ``EvalD0`` through
  ``occ::handle<Geom_Line>`` (or ``Geom_Circle``, ``Geom_BSplineCurve``),
  where ``EvalD0`` is ``final`` and the call is not virtual.
- *Dispatch*: ``VisitBuiltin``, one type test and then the concrete loop,
  against ``EvalD0`` through ``occ::handle<Geom_Curve>``, a virtual call per
  evaluation. B0.2 asks Stereon to be faster here.
- *Virtual*: the erased ``Curve``, a virtual call per evaluation, against
  the same OCCT call, like for like.

B0.2's second target, ``VisitBuiltin`` within :math:`1.1\times` of the direct
call, compares two Stereon rows; the report renders it as a table of its
own. OCCT parametrises a line by arc length, so its parameters are scaled
by :math:`\lVert d \rVert`, and the ``sum_x_mm`` counter, the
:math:`x`-coordinates summed over one pass, shows that both libraries
evaluate the same points. The rows run with B1 in the reference run;
results are quoted here from the first reference run that includes them.

Topology (B0.3–B0.5)
--------------------

The topology benchmarks build a sheet of :math:`n \times n` planar quads
with each library: Stereon's tables and TopoDS faces sharing their edges and
vertices. At 10⁶ faces this is the body that the gates of ADR-0004 and
ADR-0005 name. The figures below come from the reference run of 2026-10-20
(``bench/reports/phase0-topology-reference-2026-10-20.md``, run
39399395722, commit ``2c24d79``): the same container and runner as above,
5 interleaved repetitions for the topology rows, and both libraries in one
process.

.. list-table::
   :header-rows: 1
   :widths: 10 46 14 14 16

   * - ID
     - Measure, 10⁶ faces
     - Stereon
     - OCCT
     - Gate
   * - B0.3
     - Faces of every edge (OCCT builds its ancestor map)
     - 123 ms
     - 2516 ms
     - ≥ 5×: **20×, met**
   * - B0.3
     - Edges of every vertex
     - 53 ms
     - 3158 ms
     - ≥ 5×: **60×, met**
   * - B0.4
     - Live heap per face
     - 1255 B
     - 1601 B
     - < 50%: **38%, met** (78% at first; see below)
   * - B0.5
     - Pages copied by one surface replacement
     - 0.003 %
     - —
     - < 1%: **met**
   * - B0.5
     - Pages copied by one face deletion
     - 0.014 %
     - —
     - < 1%: **met**
   * - ADR-0005
     - Body and its copy after 100 edits, over the body
     - 1.008×
     - —
     - < 1.2×: **met**

The adjacency queries read Stereon's maintained incidence: ``EdgesOf``
walks the vertex's list, and ``FacesOf`` walks the radial ring. OCCT has to
build ``TopExp::MapShapesAndAncestors`` before it can answer the first
query, and the comparison counts that build, as a query session must pay
for it. A local edit opened on the frozen body copies a few pages of the
tables and nothing else. The whole edit adds about 1 MB of heap against
OCCT's 48 MB for ``BRepTools_ReShape``.

**Memory per face is not met.** The measure is the live heap the sheet
holds with the shape alive, read through glibc's ``mallinfo2`` the same way
for both libraries (``bench/Memory.hpp``), geometry included on both sides.
Stereon's live heap per face is almost exactly its table pages per face
(1253 B), so the cost is in the tables and geometry pools, not in anything
else.

The first run estimated that the geometry pools dominated: each of a face's
four coedges has its own pcurve and each of its two edges its own curve, and
every one was an erased ``Curve`` or ``Curve2d`` of 144 bytes (a 128-byte
buffer and a pointer, aligned to 16). Two changes followed: smaller records
(liveness folded into the generation word, fields ordered against padding)
and one buffer size per wrapper kind, so that every built-in type still fits
inline (:doc:`geometry`, *The buffer size, revisited*). B0.4 now reports its
bytes per face by table. A run of only these rows on the reference image
(run 39433824401, three repetitions) measured, at 10⁶ faces:

.. list-table::
   :header-rows: 1

   * - Table
     - Bytes per face
   * - Vertices
     - 40
   * - Edges
     - 144
   * - Coedges
     - 128
   * - Loops and faces
     - 44
   * - Curve pool (2 lines)
     - 208
   * - Pcurve pool (4 plane lines)
     - 320
   * - Live heap, all of it
     - 887, against OCCT's 1601: **55%, not met**

The geometry pools were still 60% of the total, and their wrappers most of
it: a ``Line2d`` is 32 bytes, its ``Curve2d`` 80. The pools now keep the
built-in types by value, each in its own arena, with the erased wrapper kept
for B-splines, composite and user geometry (``GeometryPool``,
``src/topo/Records.cppm``). Run 39437110770 measured, at 10⁶ faces:

.. list-table::
   :header-rows: 1

   * - Table
     - Bytes per face
   * - Records (vertices, edges, coedges, loops and faces)
     - 356
   * - Curve pool (2 lines and their slots)
     - 104
   * - Pcurve pool (4 plane lines and their slots)
     - 144
   * - Live heap, all of it
     - 608, against OCCT's 1601: **38%, met**
Splitting a face or an edge is part of ADR-0005's edit sequence and is not
measured. It needs the Euler operators of ADR-0027, which do not exist yet.

The run of 2026-10-20 predates the cached-curve fields of the edge record;
the measurement above includes them.

Booleans on the pool (M10)
--------------------------

A Boolean spends nearly all its time on its candidate face pairs: each
pair's surface intersection and the clip of every branch to both faces.
That work reads the two frozen operands only, so it runs on the Context's
scheduler, one pair per chunk since pairs cost unevenly, each writing its
own slot with a serial Context of its own that keeps its log records; the
networks are then built from the slots in pair order, the records passed on
and the first failure returned, so that the result does not depend on the
thread count (ADR-0007). The seeds of the classification's regions run the
same way, one region per chunk. ``Detail::ParallelFor`` takes the grain for
this: the chunks are cut by the count and the grain alone.

The fixture is the plan's: a slab of 10 x 10 cells per hole, 10 mm thick,
cut by one operand made of a cylinder of radius 2 through every cell
(``tests/boolean/SpeedupTests.cpp``, hidden behind ``[speedup]``). On a
6-core, 12-thread development machine, the debug build, the serial
scheduler against 8 threads:

.. list-table::
   :header-rows: 1

   * - Holes
     - Serial
     - 8 threads
     - Speed-up
   * - 25
     - 7.58 s
     - 1.54 s
     - 4.9x
   * - 49
     - 16.07 s
     - 2.93 s
     - 5.5x
   * - 100
     - 32.53 s
     - 6.16 s
     - 5.3x

These are a development machine's figures, not the measurement of record:
ADR-0007 asks for at least 3x on the 4 cores of the OCCT reference
workflow's runner, on a body of more than 1 000 faces, measured as that
workflow measures.

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

``--benchmark_filter=B1.6`` runs one group, ``--benchmark_filter=B0.2``
the dispatch rows, about a second per row. The ``length_mm`` counter of
the arc-length rows prints each library's result, so the accuracy check
travels with the timing.
