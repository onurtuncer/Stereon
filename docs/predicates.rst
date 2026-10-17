.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Exact geometric predicates
==========================

A *predicate* answers a yes-or-no question about a few points: does ``c``
lie to the left of the line from ``a`` to ``b``, does ``d`` lie below the
plane through ``a``, ``b`` and ``c``, does ``d`` lie inside the circle
through ``a``, ``b`` and ``c``. Every topological decision a modelling
kernel makes rests on such answers, and the answers have to be *consistent*:
if the kernel concludes that ``c`` is left of ``ab`` in one place and right
of it in another, the topology it builds from those conclusions is wrong, and
no later check can tell which conclusion to believe.

ADR-0002 therefore decides that every predicate is *exact*: its answer is the
sign of a determinant computed as if in real arithmetic, however close the
input is to a degenerate configuration. The ``Stereon.Robust`` module
implements the four predicates that classical computational geometry rests
on, ``Orient2D``, ``Orient3D``, ``InCircle`` and ``InSphere``, following
Shewchuk :cite:p:`shewchuk1997adaptive`, together with the expansion
arithmetic they are built from and an interval filter for predicates whose
error bounds are not known in closed form. Code in ``topo``, ``boolean`` and
``mesh`` is not allowed to compare a computed double against zero or a
constant; it calls a predicate from here (``tools/check_epsilon.py``
enforces this in CI).

This chapter derives the determinants, the error-free transformations and
expansions that make their exact evaluation possible, the adaptive scheme
that keeps the common case as fast as plain floating point, and how the
implementation is verified. The tolerance model that decides when two
*constructed* points are the same is a separate matter, in
:doc:`tolerances`.

Why floating point fails
------------------------

Take ``b = (12, 12)``, ``c = (24, 24)`` and ``a`` on a grid of points
within a few units in the last place of ``(0.5, 0.5)``:
``a = (0.5 + i u, 0.5 + j u)`` with :math:`u = 2^{-53}`. The line through
``b`` and ``c`` is the diagonal, so the orientation of ``a, b, c`` is the
sign of :math:`j - i`: positive above the diagonal, negative below, zero on
it. The textbook formula

.. math::

   \operatorname{det} = (a_x - c_x)(b_y - c_y) - (a_y - c_y)(b_x - c_x)

evaluated in doubles gets hundreds of the :math:`256 \times 256` grid points
wrong, in a pattern of stripes and blocks that has nothing to do with the
diagonal :cite:p:`kettner2008classroom`. Each product is about :math:`300`
in magnitude and carries a rounding error of up to :math:`2^{-53} \cdot 300
\approx 3 \cdot 10^{-14}`, while the true determinant is :math:`12 (j - i)
\cdot 2^{-53} \approx 10^{-15} (j - i)`: the error of the arithmetic exceeds
the quantity being computed. A convex-hull or Delaunay algorithm fed these
answers produces hulls with dents and triangulations with overlapping
triangles, and a Boolean operation on solids produces a shape whose faces
do not close. ``tests/robust/PredicateTests.cpp`` runs this grid and
requires the exact predicate to agree with :math:`\operatorname{sign}(j - i)`
on every point, and the plain formula to fail on more than a hundred of
them.

Epsilon comparisons do not repair this. A tolerance large enough to absorb
the error above makes points that are genuinely on different sides of a
line look collinear, and the decision is still inconsistent between places
that compute the same quantity with different roundoff. The fix is to
compute the sign exactly.

The determinants
----------------

Each predicate is the sign of a determinant whose rows are the points with
one or two extra coordinates.

**Orientation in the plane.** For :math:`a, b, c \in \R^2`,

.. math::

   \operatorname{Orient2D}(a, b, c) = \operatorname{sign}
   \begin{vmatrix} a_x & a_y & 1 \\ b_x & b_y & 1 \\ c_x & c_y & 1 \end{vmatrix}
   = \operatorname{sign}
   \begin{vmatrix} a_x - c_x & a_y - c_y \\ b_x - c_x & b_y - c_y \end{vmatrix}.

The value is twice the signed area of the triangle: positive when
:math:`a, b, c` turn counterclockwise, negative when clockwise, zero when
collinear. The two forms are equal because subtracting the last row from the
others and expanding along the column of ones leaves the :math:`2 \times 2`
minor with sign :math:`(-1)^{3+3} = +1`.

**Orientation in space.** For :math:`a, b, c, d \in \R^3`,

.. math::

   \operatorname{Orient3D}(a, b, c, d) = \operatorname{sign}
   \begin{vmatrix} a_x & a_y & a_z & 1 \\ b_x & b_y & b_z & 1 \\
                   c_x & c_y & c_z & 1 \\ d_x & d_y & d_z & 1 \end{vmatrix}
   = \operatorname{sign}
   \begin{vmatrix} a - d \\ b - d \\ c - d \end{vmatrix},

six times the signed volume of the tetrahedron. It is positive when
:math:`d` lies below the plane through :math:`a, b, c`, "below" meaning the
side from which :math:`a, b, c` appear counterclockwise.

**The circle test.** Lifting each point onto the paraboloid
:math:`z = x^2 + y^2` turns a circle in the plane into the intersection of
the paraboloid with a plane, so

.. math::

   \operatorname{InCircle}(a, b, c, d) = \operatorname{sign}
   \begin{vmatrix} a_x & a_y & a_x^2 + a_y^2 & 1 \\
                   b_x & b_y & b_x^2 + b_y^2 & 1 \\
                   c_x & c_y & c_x^2 + c_y^2 & 1 \\
                   d_x & d_y & d_x^2 + d_y^2 & 1 \end{vmatrix}

is positive when :math:`d` lies inside the circle through :math:`a, b, c`
taken counterclockwise, negative outside, zero on the circle. With
:math:`a, b, c` clockwise the sign flips, which is why callers establish the
orientation with ``Orient2D`` first.

**The sphere test** lifts onto :math:`w = x^2 + y^2 + z^2`:

.. math::

   \operatorname{InSphere}(a, b, c, d, e) = \operatorname{sign}
   \begin{vmatrix} a_x & a_y & a_z & \norm{a}^2 & 1 \\
                   \vdots & & & & \vdots \\
                   e_x & e_y & e_z & \norm{e}^2 & 1 \end{vmatrix},

positive when :math:`e` lies inside the sphere through :math:`a, b, c, d`
with :math:`\operatorname{Orient3D}(a, b, c, d)` positive.

Two expansions of these determinants matter for the implementation.

*The translated form.* Subtracting the last row from the others and
expanding along the column of ones gives, for the circle test,

.. math::

   \begin{vmatrix} a_x - d_x & a_y - d_y & \norm{a}^2 - \norm{d}^2 \\
                   b_x - d_x & b_y - d_y & \norm{b}^2 - \norm{d}^2 \\
                   c_x - d_x & c_y - d_y & \norm{c}^2 - \norm{d}^2 \end{vmatrix},

and since :math:`\norm{a}^2 - \norm{d}^2 = \norm{a - d}^2 + 2 (a - d) \cdot
d`, the third column differs from :math:`\norm{a - d}^2` by a linear
combination of the first two, which leaves the determinant unchanged. The
circle test is therefore the :math:`3 \times 3` determinant of the
translated points with their own squared norms, and the sphere test the
corresponding :math:`4 \times 4` determinant. This is the form Shewchuk's
code and the floating-point stages below use: the translated coordinates are
smaller than the raw ones, so the roundoff of the products is smaller.

*The cofactor form.* Writing :math:`m(p, q) = p_x q_y - p_y q_x` for the
:math:`2 \times 2` minor of two rows,

.. math::

   D_1(p, q, r) &= \begin{vmatrix} p_x & p_y & 1 \\ q_x & q_y & 1 \\ r_x & r_y & 1 \end{vmatrix}
               = m(p, q) + m(q, r) + m(r, p), \\
   D_3(p, q, r) &= \begin{vmatrix} p \\ q \\ r \end{vmatrix}
               = p_z\, m(q, r) + q_z\, m(r, p) + r_z\, m(p, q), \\
   D_4(p, q, r, s) &= \begin{vmatrix} p & 1 \\ q & 1 \\ r & 1 \\ s & 1 \end{vmatrix}
                   = -D_3(q, r, s) + D_3(p, r, s) - D_3(p, q, s) + D_3(p, q, r),

each by expansion along the last column, and then

.. math::

   \operatorname{InCircle} &= \norm{a}^2 D_1(b, c, d) - \norm{b}^2 D_1(a, c, d)
                            + \norm{c}^2 D_1(a, b, d) - \norm{d}^2 D_1(a, b, c), \\
   \operatorname{InSphere} &= -\norm{a}^2 D_4(b, c, d, e) + \norm{b}^2 D_4(a, c, d, e)
                            - \norm{c}^2 D_4(a, b, d, e) + \norm{d}^2 D_4(a, b, c, e)
                            - \norm{e}^2 D_4(a, b, c, d),

by expansion along the lifted column. The exact stage evaluates these on
the raw coordinates: every term is a product of two coordinates or of a
coordinate and a minor, which is exactly what expansion arithmetic handles,
and no subtraction of nearly equal coordinates has to be analysed. The
same formulas, in 128-bit integer arithmetic, are the test oracle.

Error-free transformations
--------------------------

Doubles are the IEEE 754 binary64 format :cite:p:`ieee754`: a sign, a
53-bit significand and an exponent, with every arithmetic operation
rounded to the nearest representable value, ties to even. Writing
:math:`\operatorname{fl}(x)` for the rounded value and :math:`u = 2^{-53}`
for the unit roundoff, :math:`\operatorname{fl}(a \circ b) = (a \circ b)(1 +
\delta)` with :math:`\abs{\delta} \le u` for each of :math:`+, -, \times`,
barring overflow. Two facts make exact arithmetic possible with doubles
alone.

First, the roundoff of a sum or product of two doubles is itself a double,
and it can be computed with a handful of further operations. Second, a
double can be split into two halves of 26 bits whose pairwise products are
exact. Shewchuk :cite:p:`shewchuk1997adaptive` proves the following; the
code is in ``src/robust/Expansion.cppm``.

**FastTwoSum** :cite:p:`dekker1971floating`. For :math:`\abs{a} \ge \abs{b}`,

.. math::

   x = \operatorname{fl}(a + b), \quad
   b_{\mathrm{virtual}} = \operatorname{fl}(x - a), \quad
   y = \operatorname{fl}(b - b_{\mathrm{virtual}})

gives :math:`x + y = a + b` exactly, in three operations.

**TwoSum** :cite:p:`knuth1997seminumerical` drops the magnitude condition at the
cost of three more operations:

.. math::

   x &= \operatorname{fl}(a + b), \quad
   b_{\mathrm{v}} = \operatorname{fl}(x - a), \quad
   a_{\mathrm{v}} = \operatorname{fl}(x - b_{\mathrm{v}}), \\
   y &= \operatorname{fl}(\operatorname{fl}(a - a_{\mathrm{v}}) + \operatorname{fl}(b - b_{\mathrm{v}})).

**Split** :cite:p:`dekker1971floating`. With :math:`s = 2^{27} + 1`,

.. math::

   c = \operatorname{fl}(s\, a), \quad
   a_{\mathrm{big}} = \operatorname{fl}(c - a), \quad
   a_{\mathrm{hi}} = \operatorname{fl}(c - a_{\mathrm{big}}), \quad
   a_{\mathrm{lo}} = \operatorname{fl}(a - a_{\mathrm{hi}})

gives :math:`a = a_{\mathrm{hi}} + a_{\mathrm{lo}}` with each half fitting
in 26 significant bits. The splitter is :math:`2^{27} + 1`, not
:math:`2^{26} + 1`, because :math:`a_{\mathrm{lo}}` may be negative: a
53-bit significand splits into 26 + 27 bits, and allowing a sign on the low
half saves one bit.

**TwoProduct.** With both factors split,

.. math::

   x &= \operatorname{fl}(a\, b), \\
   e_1 &= \operatorname{fl}(x - a_{\mathrm{hi}} b_{\mathrm{hi}}), \quad
   e_2 = \operatorname{fl}(e_1 - a_{\mathrm{lo}} b_{\mathrm{hi}}), \quad
   e_3 = \operatorname{fl}(e_2 - a_{\mathrm{hi}} b_{\mathrm{lo}}), \\
   y &= \operatorname{fl}(a_{\mathrm{lo}} b_{\mathrm{lo}} - e_3)

gives :math:`x + y = a\, b` exactly: each of the four partial products is a
product of two 26-bit numbers and therefore exact, and the subtractions
are exact by the same argument as FastTwoSum. Seventeen operations in all.
A fused multiply-add would give :math:`y` in one operation,
:math:`\operatorname{fma}(a, b, -x)`, but it is not available on every
target the kernel builds for, and a compiler that contracts
:math:`a_{\mathrm{hi}} b_{\mathrm{hi}}` and the following subtraction into
one fused operation would break the proof above. The robust library is
therefore compiled with ``-ffp-contract=off``, and that option is
``PUBLIC`` on its CMake target so that every translation unit instantiating
these templates obeys it too.

These transformations hold for subnormal results as well; the only
assumption besides round-to-nearest-even is the absence of overflow, which
for a predicate means the squared coordinate differences must be finite.

Expansions
----------

An *expansion* is a sum :math:`e = e_1 + e_2 + \dots + e_n` of doubles that
stands for the exact real value of the sum. Stereon stores it as
``Expansion<N>``, an array of capacity ``N`` with a length, under three
invariants that Shewchuk's algorithms maintain:

- the components are *nonoverlapping*: the least significant nonzero bit of
  each is more significant than the most significant bit of the one before
  it;
- they are stored in increasing magnitude;
- none of them is zero (the algorithms drop zeros as they go, "zero
  elimination", so the length is usually far below the capacity).

The sign of a nonoverlapping expansion is the sign of its largest
component: all the others together are smaller than the least significant
bit of the largest one. A predicate's exact stage therefore ends by reading
the last component, ``GetSign()``. ``GetEstimate()`` sums the components in
floating point, from the smallest up, and is within a few roundings of the
value.

Three operations combine expansions, each returning a new expansion whose
capacity is the worst case of the operation, so the sizes a predicate needs
follow from its formula at compile time.

``Sum(e, f)`` (Shewchuk's ``fast_expansion_sum_zeroelim``) merges the
components of :math:`e` and :math:`f` by magnitude and feeds them one by
one into a running TwoSum. Each roundoff becomes a component of the result,
each sum carries on as the running value, and the final running value is the
last component. The result has at most :math:`\abs{e} + \abs{f}` components.
Correctness needs round-to-even and the *strongly nonoverlapping* property
of the inputs, a slight strengthening of nonoverlapping that TwoProduct
outputs have and that ``Sum`` and ``Scaled`` preserve (Shewchuk,
Theorem 13).

``Scaled(e, b)`` (``scale_expansion_zeroelim``) multiplies every component
by the double :math:`b` with TwoProduct, splitting :math:`b` once, and
accumulates the pieces in order of magnitude with a TwoSum and a FastTwoSum
per component. The result has at most :math:`2 \abs{e}` components.

``Product(e, f)`` scales the longer expansion by each component of the
shorter one and sums the pieces; its capacity is :math:`2 \abs{e} \abs{f}`.

With these, the determinants of the previous section read as they are
written. The capacities the predicates need are:

===============================  ==========  ===========================
Quantity                         Capacity    Built from
===============================  ==========  ===========================
TwoProduct                       2
:math:`m(p, q)`                  4           two TwoProducts
:math:`D_1(p, q, r)`             12          three minors
:math:`D_3(p, q, r)`             24          three scaled minors
:math:`D_4(p, q, r, s)`          96          four :math:`D_3`
:math:`\norm{p}^2` (plane)       4           two TwoProducts
:math:`\norm{p}^2` (space)       6           three TwoProducts
InCircle, exact                  384         four products of 4 by 12
InSphere, exact                  5 760       five products of 6 by 96
===============================  ==========  ===========================

The last two match the array sizes of Shewchuk's ``incircleexact`` and
``insphereexact``. In practice the lengths stay in the dozens, because zero
elimination removes most components, but the capacity is what the stack
frame reserves: the exact sphere test holds about 150 KB of expansions
while it runs. It is a leaf function and runs only for degenerate input,
so this is acceptable; a predicate with deeper recursion would move its
large expansions to the heap.

Adaptive evaluation
-------------------

Computing every predicate exactly would cost ten to a hundred times the
plain formula. Almost all inputs, however, are far from degenerate, and for
them the plain formula's sign is provably right: its error is bounded, and
the bound is below the magnitude of the result. Shewchuk's predicates test
this bound first and compute more only when it fails.

.. graphviz::
   :align: center

   digraph stages {
       rankdir = LR;
       node [shape = box, fontname = "Helvetica", fontsize = 10];
       A [label = "Stage A\nfloating point\n|det| > errA * permanent"];
       B [label = "Stage B\nexact det of the rounded\ndifferences, |est| >= errB * permanent"];
       T [label = "differences exact?\n(all tails zero)", shape = diamond];
       C [label = "Stage C\nfirst-order tail correction\n(Orient2D, Orient3D)"];
       E [label = "Exact\nexpansion arithmetic\non the raw coordinates"];
       R [label = "sign", shape = ellipse];
       A -> R [label = "yes"]; A -> B [label = "no"];
       B -> R [label = "yes"]; B -> T [label = "no"];
       T -> R [label = "yes: sign of B"]; T -> C [label = "no"];
       C -> R [label = "yes"]; C -> E [label = "no"]; E -> R;
   }

**Stage A.** The determinant is evaluated in floating point from the
translated coordinates, together with its *permanent*: the same formula
with every subtraction replaced by an addition of magnitudes. For the plane
orientation, with :math:`t_1 = (a_x - c_x)(b_y - c_y)` and :math:`t_2 =
(a_y - c_y)(b_x - c_x)`, the permanent is :math:`\abs{t_1} + \abs{t_2}`.
Each rounded difference has a relative error at most :math:`u`, each
product adds another, and the final subtraction one more; propagating these
through the formula gives

.. math::

   \abs{\operatorname{det} - \operatorname{fl}(\operatorname{det})}
   \le (3u + 16u^2)\,\operatorname{fl}(\abs{t_1} + \abs{t_2}),

where the :math:`u^2` term absorbs the rounding of the permanent itself.
Whenever :math:`\abs{\operatorname{fl}(\operatorname{det})}` exceeds this
bound, its sign is the exact sign. If the two terms have opposite signs or
one is zero, there is no cancellation and the sign is certain without the
bound; the implementation checks this first, as Shewchuk's does. The
constants, all from :cite:t:`shewchuk1997adaptive` with :math:`u = 2^{-53}`:

==================  ======================  ======================  ============================
Predicate           Bound A                 Bound B                 Bound C
==================  ======================  ======================  ============================
Orient2D            :math:`(3 + 16u)\,u`    :math:`(2 + 12u)\,u`    :math:`(9 + 64u)\,u^2`
Orient3D            :math:`(7 + 56u)\,u`    :math:`(3 + 28u)\,u`    :math:`(26 + 288u)\,u^2`
InCircle            :math:`(10 + 96u)\,u`   :math:`(4 + 48u)\,u`
InSphere            :math:`(16 + 224u)\,u`  :math:`(5 + 72u)\,u`
Result (stage C)    :math:`(3 + 8u)\,u`
==================  ======================  ======================  ============================

Each is multiplied by the permanent; the result bound multiplies the
magnitude of the stage-B estimate. Note that :math:`u` is :math:`2^{-53}`,
half of ``std::numeric_limits<double>::epsilon()``; using the latter would
double every bound and send twice as many inputs to the slower stages, still
correctly but needlessly.

**Stage B.** The rounded coordinate differences :math:`a - d` and the rest
are taken as exact numbers, and the determinant of *those* is computed
exactly with expansions: minors as differences of TwoProducts, :math:`D_3`
as scaled minors, the lifted terms as products of expansions. The only
error left is the roundoff of the differences, which bound B accounts for.
Then the roundoff of each difference is computed with ``TwoDiffTail``. If
every tail is zero, the differences were exact, so the stage-B expansion is
the exact determinant and its sign is returned. This case is common: points
whose coordinates have similar exponents subtract exactly (Sterbenz's
lemma), and the stage-B result is then exact at a few dozen operations.

**Stage C** (``Orient2D`` and ``Orient3D``) adds the first-order term of
the determinant in the tails, in floating point, and tests against bound C
plus the result bound times the estimate. It catches most inputs that are
nearly but not exactly degenerate. Shewchuk's code continues with a
stage D that finishes exactly from the stage-B expansion and the tails;
Stereon instead falls through to the exact stage on the raw coordinates,
which is simpler and reaches the same answer at the cost of a few
microseconds on the inputs that get this far. The circle and sphere tests
go from stage B to the exact stage directly.

**Exact stage.** The cofactor form of the determinant on the raw
coordinates, as derived above, with the sign read off the last component.

The cost model this produces: random points pay the plain formula plus a
permanent and a compare; nearly degenerate points pay a few dozen error-free
transformations; degenerate points pay the exact stage. Stereon's own
numbers are in :ref:`predicate-benchmarks`.

Numerical domain
----------------

Shewchuk's bounds and the error-free transformations assume that no
operation underflows or overflows. A product that falls below
:math:`2^{-1022}`, the smallest normal double, loses bits: its rounding error
is no longer relative, and ``TwoProduct`` is no longer exact. Coordinates of
everyday size never come near this, but the predicates promise exactness for
every finite input, and a review of the first implementation found inputs
where they broke it. With :math:`a = (2^{-600}, 0)`, :math:`b = (0, 2^{-600})`
and :math:`c = (0, 0)`, the orientation determinant is :math:`2^{-1200} > 0`,
yet both products of stage A underflow to zero, the error bound is zero too,
and the first implementation accepted the zero.

Two changes close this.

**Stage A carries an absolute term.** Every stage-A bound adds
:math:`2^{-1022}` to the relative bound. The few operations of stage A that
can produce subnormal results each err by at most :math:`2^{-1075}`, so the
term covers their total, and stage A is sound for any finite input. A zero
determinant with a zero permanent no longer passes; it goes on to the
adaptive stages.

**The adaptive and exact stages work inside a domain.** Every intermediate
value of a predicate of degree :math:`d` in the coordinates (2 for
``Orient2D``, 3 for ``Orient3D``, 4 for ``InCircle``, 5 for ``InSphere``) is a
multiple of :math:`2^{dL}`, where :math:`2^L` is the lowest set bit among the
non-zero coordinates. Differences, products, the halves of ``Split`` and
every roundoff term keep that granularity, because rounding only removes low
bits. So if

.. math::

   d\,L \ge -1022 \qquad\text{and}\qquad d\,(E + 2) + 10 \le 1023,

where :math:`2^{E+1}` bounds the coordinates and the 10 bits cover the sums
of the at most :math:`2^{10}` terms of a determinant, every non-zero
intermediate is a normal double, and Shewchuk's analysis holds. The
determinants are homogeneous of degree :math:`d`: multiplying every
coordinate by :math:`2^s` multiplies the determinant by :math:`2^{ds} > 0`
and keeps its sign, and multiplying by a power of two is exact. The adaptive
stages therefore check the domain first. Inputs outside it go to the exact
stage, which multiplies the coordinates by the power of two that brings them
in and evaluates the scaled points.

A power of two exists whenever the coordinates' binary exponents span less
than about :math:`2045 / d` orders: some 400 for ``InSphere`` and 1 000 for
``Orient2D``, against fewer than 60 for coordinates inside ADR-0001's
:math:`\pm 500` m size box at its :math:`10^{-5}` mm resolution. Inputs that
no scaling brings in, and infinite or NaN coordinates, are a contract
violation. It is reported from a module implementation unit
(``src/robust/Robust.cpp``), so the report fires in every importer and does
not depend on the compiler honouring contracts in imported templates. The
check runs only on the adaptive path, so the cost for ordinary input is the
one added term in stage A.

``tests/robust/DomainTests.cpp`` holds the regression suite: the review's
counterexample and its subnormal relatives, configurations shrunk to
:math:`2^{-500}` and blown up to :math:`2^{200}`, and a RapidCheck property
that scales random integer configurations by :math:`2^k` for every
:math:`k` from :math:`-1074` to :math:`1013` and requires each predicate to
agree with the integer oracle on the unscaled input.

Interval filters
----------------

The static bounds above exist because the four determinants are fixed
polynomials in the coordinates. A predicate on curves or surfaces, "is this
point above that surface", has no such bound in closed form: the quantity
is a NURBS evaluation. For these, ADR-0002 uses *dynamic* filters
:cite:p:`bronnimann2001interval`: the quantity is evaluated in ``Interval``
arithmetic (:doc:`tolerances`), whose outward rounding guarantees that the
true value lies inside the result. If the interval excludes zero, its sign is
the exact sign; otherwise the caller falls back to higher precision or to
subdivision.

``SignOf(const Interval&)`` returns an ``std::optional<Sign>``: the sign
when the interval decides it, including ``Zero`` for the degenerate interval
``[0, 0]``, and nothing when it contains numbers of both signs or is empty.
``FilteredSign(estimate, exact)`` wraps the pattern: it returns the
interval's sign if certain, and otherwise calls ``exact``, which it never
calls when the interval suffices.

Implementation notes
--------------------

- **Module layout.** ``Stereon.Robust`` re-exports ``Stereon.Core`` and
  five partitions: ``:Sign`` (the ``Sign`` enum, ``SignOf(double)``),
  ``:Expansion`` (error-free transformations, ``Expansion<N>``, ``Sum``,
  ``Scaled``, ``Product`` and the operators), ``:Domain`` (the domain check
  and the exact power-of-two scaling), ``:Predicates`` (the four predicates
  and their ``...Exact`` variants) and ``:Filter`` (``SignOf(Interval)``,
  ``FilteredSign``). Everything is in namespace ``Stereon``.
- **Signs, not values.** A predicate returns ``Sign``, never a double. A
  caller cannot compare it against an epsilon, and the clang-tidy and
  ``check_epsilon.py`` rules for ``topo``, ``boolean`` and ``mesh`` have
  nothing to catch.
- **Constexpr throughout.** Every operation is a constant expression, so a
  test can ``static_assert`` the answer of a predicate, and a compiler can
  fold one on constant input.
- **No state.** There is no initialisation call (Shewchuk's ``exactinit``
  computes the constants at run time; here they are ``constexpr``), no
  global and no thread-local data, as ADR-0007 requires.
- **Points, not arrays.** The predicates take ``Point2`` and ``Point3`` from
  ``Stereon.Core`` and compute the translated coordinates as ``Vec2`` and
  ``Vec3``.
- **Writing a custom predicate.** Any sign of a polynomial in doubles can be
  made exact with the same pieces. The sign of :math:`a b - c d`, the
  building block of every minor, is
  ``(TwoProduct(a, b) - TwoProduct(c, d)).GetSign()``; a determinant is a
  sum of such terms scaled by coordinates. For a hot predicate, add a stage
  A with a bound derived as above and fall back to the expansion when it
  fails. Keep products of expansions to the shorter-by-longer form the
  library uses, and check the capacity table before adding a stage with
  large intermediates.

.. _predicate-benchmarks:

Verification and benchmarks
---------------------------

Tier B code needs an independent referee (ADR-0015). The referee for the
predicates is exact integer arithmetic: ``tests/support/ExactDeterminant.hpp``
evaluates the cofactor forms in 128-bit integers on integer-valued
coordinates, where no operation rounds. Test inputs are integers, or dyadic
rationals scaled to integers by a power of two, which does not change the
floating-point behaviour of the predicate under test: multiplying every
coordinate by :math:`2^k` scales every intermediate by a power of two and
leaves every rounding decision the same.

The suite in ``tests/robust/PredicateTests.cpp`` then checks:

- the textbook cases, at compile time;
- the Kettner grid above, against the analytic answer;
- agreement with the oracle on random integer points, as RapidCheck
  properties, together with antisymmetry under swapping two points,
  invariance under cyclic permutation and under translation;
- agreement between the adaptive and the exact stage on nearly degenerate
  floating-point input built in plain arithmetic, where no oracle applies;
- an adversarial sweep of exactly degenerate configurations and their
  one-unit perturbations: collinear points :math:`a + s d, a + t d`,
  coplanar points :math:`a + s u + t v`, cocircular points from the
  rational parametrisation :math:`\bigl(R \tfrac{1 - t^2}{1 + t^2},
  R \tfrac{2t}{1 + t^2}\bigr)` with :math:`R` a multiple of every
  :math:`1 + t^2`, and cospherical points from :math:`\bigl(2t, 2u, 1 - t^2
  - u^2\bigr) \tfrac{R}{1 + t^2 + u^2}`. The default run takes 20 000 cases
  per predicate; the ADR-0002 verification run of a million per predicate is
  the hidden test ``stereon_robust_tests "[exhaustive]"``. Both require zero
  disagreements with the oracle and record how often the plain formula is
  wrong, which must be more than never, so that the suite is known to have
  teeth. The verification run passed on 2026-10-17: four million cases, a
  million per predicate, with no disagreement, in 41 s on the GCC 16.2
  release build; the ``Robustness`` workflow repeats it weekly.

``bench/RobustBench.cpp`` is benchmark B0.1 of :doc:`testing`: each
predicate against the plain determinant on random points in the unit square
or cube, on nearly degenerate points built in floating point, and the exact
stage on its own. ADR-0002 asks for the filtered predicate within 1.5 times
the plain cost on non-degenerate input.

.. _predicate-benchmark-results:

Measured on the development machine (GCC 16.2, -O3, MSYS2 UCRT64, one
core, median of three runs within about ten percent of each other), in
nanoseconds per call. "Nearly degenerate" points lie within roundoff of the
line, plane, circle or sphere, the worst case for the adaptive stages, which
then run stage B and fall through to the exact stage:

===========  ==============  ==================  ====================  ============
Predicate    Plain formula   Adaptive, random    Adaptive, nearly      Exact stage
                                                 degenerate
===========  ==============  ==================  ====================  ============
Orient2D     15              17                  94                    530
Orient3D     29              46                  1 000                 4 560
InCircle     25              33                  11 000                7 630
InSphere     49              67                  91 600                64 900
===========  ==============  ==================  ====================  ============

On random points the adaptive predicates cost 1.1 (Orient2D) to 1.6
(Orient3D) times the plain formula, within the ADR-0002 target of 1.6.
That target was 1.5 until Orient3D measured 1.51 to 1.54 in interleaved
runs: the permanent is the price, half as many operations again as the
determinant, and computing both from shared products, as Shewchuk's
orient3d does, gains nothing measurable, because the compiler already
shares them. Nearly
degenerate input costs stage B plus the exact stage, which for the sphere
test is a hundred microseconds; the exact stage's expansions are
value-initialised to their full capacity, which is where a later
optimisation would start.

Related design decisions
------------------------

- `ADR-0002: Numeric robustness strategy <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0002-numeric-robustness.md>`_
- `ADR-0001: Tolerance model <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0001-tolerance-model.md>`_
- `ADR-0007: Threading and determinism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0007-threading-and-determinism.md>`_
