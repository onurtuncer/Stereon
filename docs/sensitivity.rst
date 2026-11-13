.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Shape sensitivities
===================

The types derivatives of shapes are expressed in (ADR-0021, stage 1), and
the first products (stage 2): forward and reverse products of a B-spline
curve's points in its poles and weights. Every other operation declares what
a derivative would depend on and answers ``NotSupported``, never a silent
zero.

Parameters
----------

A ``ParameterSet`` lists a model's design parameters in order, and that order
is the convention of every vector over them: component *k* of a seed, a
gradient or a Jacobian row belongs to the parameter whose ``ParameterId`` has
index *k*. Each parameter has a quantity (``Length`` in millimetres,
``Angle`` in radians, or ``Dimensionless``, per ADR-0001) and a reference
scale in that unit. Seeds and derivatives are per reference, so a step of 1
moves a parameter by its reference and millimetre and radian parameters stay
comparable.

A ``ParameterId`` is bound to the set that made it, not just to an index.
Two sets whose first parameters are both lengths give different ids, and
neither set accepts the other's. For a foreign id, ``GetName`` and
``GetUnit`` return an empty string and ``GetReference`` returns NaN, never
another parameter's values. A set can be moved but not copied, so two sets
that grew apart never share an id. A moved-from set is empty and has an
identity of its own.

.. code-block:: cpp

   Stereon::ParameterSet parameters;
   const auto radius = parameters.Add("radius", Stereon::ParameterQuantity::Length, 10.0);
   const auto angle = parameters.Add("angle", Stereon::ParameterQuantity::Angle);

Statuses and results
--------------------

A completed evaluation returns ``Result<SensitivityResult>``: a
``KernelError`` for invalid input (a seed of the wrong length, a value that is
not finite), otherwise a result with one of five statuses.

================== ==============================================================
Status             Meaning
================== ==============================================================
Valid              computed under the declared assumptions and checks
NearEvent          usable but qualified: near a branch event or poorly conditioned
Nondifferentiable  no ordinary derivative under the selected semantics
NotSupported       the operation or product has no implementation
Indeterminate      correspondence, convergence or regularity not established
================== ==============================================================

Only ``Valid`` and ``NearEvent`` carry values; for the others ``GetValues()``
is empty, so an absent derivative cannot be read as zero. A result also
carries the branch and topology signature, the entities it concerns, the
tolerance policy held fixed and the solver residual and conditioning, each
NaN until known.

Capability declarations
-----------------------

Every public geometric operation declares, before any evaluation, an
``OperationCapability``: what its result depends on (defining values, input
geometry, input points, the operation's own values, tolerances, topology),
the conventions its samples may be tracked by, the assumptions under which a
derivative would be valid, and the products it implements (JVP, VJP, dense
Jacobian; so far the JVP and VJP of B-spline curve evaluation). Each library keeps a catalogue in its ``Capabilities``
partition: ``CurveCapabilities()``, ``SurfaceCapabilities()``,
``IntersectCapabilities()``, ``TopoCapabilities()``, ``BuildCapabilities()``
and ``BooleanCapabilities()``. ``KernelCapabilities()``, in
``Stereon.Boolean``, the library above all the others, returns them together.

.. code-block:: cpp

   const auto* project = Stereon::FindCapability(Stereon::CurveCapabilities(), "Project");
   // project->Supports(Stereon::DerivativeProduct::Jvp) is false: ask, and the
   // answer is SensitivityResult::NotSupported(*project, DerivativeProduct::Jvp).

``tools/check_capabilities.py`` enforces the rule. It reads every function
exported from ``src/curve``, ``src/surface``, ``src/intersect``,
``src/topo``, ``src/build`` and ``src/boolean`` and fails when one is
neither declared in its library's catalogue nor listed in the script's
``NOT_OPERATIONS`` with a reason (an evaluator, a query, an accessor, a
validator, a classification); it also fails on declarations
of functions that no longer exist. It runs in the lint workflow and in CTest.

Correspondence
--------------

A ``TrackedSample`` is a point an analysis follows as parameters move: on a
face at (*u*, *v*) or on an edge at *t*, with its orientation and its
tracking convention (fixed surface coordinates, fixed edge parameter, moving
trim coordinates or unique projection). It names its entity by the store of
the snapshot it belongs to and the slot there, never by the slot alone, so a
compacted or edited snapshot, whose slots mean something else, never matches
it. Carrying a sample across an operation is the History's job (ADR-0011).

The sensitivity context
-----------------------

A ``SensitivityContext`` holds what one evaluation owns: its parameter set,
the seed of a forward product, the adjoint of a reverse one and solver
workspace. The caller owns it; it is separate from ``Context`` and never
stored in a ``Shape``. It cannot be copied, so two contexts share no mutable
state and evaluations in separate contexts can run concurrently.

Products of B-spline curve evaluation
-------------------------------------

``EvalJvp`` and ``EvalVjp`` differentiate the points ``Eval`` gives a
``BSplineCurve`` at fixed parameters with respect to design parameters bound
to its poles and weights. A ``PoleBinding`` binds a parameter to one
coordinate of one pole (a ``Length``) or to its weight (``Dimensionless``).
A parameter may move several components, as a translation moves every pole,
and a component may be moved by several parameters; contributions add. The
knots, the degree and the pole count are structural and stay fixed, so the
map is smooth wherever the weights stay positive and every evaluation is
``Valid``. A weight bound on a polynomial curve is a weight of 1 that may
move.

With basis functions :math:`N_i`, :math:`W = \sum_j N_j w_j` and
:math:`R_i = N_i w_i / W`, the point at :math:`t` moves with a coordinate of
pole :math:`i` by :math:`R_i(t)` and with its weight by
:math:`N_i(t) (P_i - C(t)) / W(t)`, each times the parameter's reference.

- ``EvalJvp`` takes the context's seed and returns three values per point
  (x, y, z), in the order of the parameters given.
- ``EvalVjp`` takes the context's adjoint, three components per point under
  the Euclidean pairing, and returns one value per parameter of the set, in
  its order, accumulated point by point in order (ADR-0007).
- Both return the same signature, a digest of the degree, pole count,
  rationality and knots, so two results with different signatures were not
  computed in one regime.

.. code-block:: cpp

   Stereon::ParameterSet parameters;
   const auto lift = *parameters.Add("lift", Stereon::ParameterQuantity::Length);
   const std::vector<Stereon::PoleBinding> bindings{{lift, 2, Stereon::PoleComponent::Z},
                                                    {lift, 3, Stereon::PoleComponent::Z}};
   Stereon::SensitivityContext sensitivity(std::move(parameters));
   (void)sensitivity.SetSeed(std::vector<double>{1.0});
   const auto velocity = Stereon::EvalJvp(curve, bindings, samples, sensitivity);

The tests compare ``EvalJvp`` with forward-mode dual numbers pushed through
the scalar-generic ``Eval``, an independent derivative, and check that
``EvalVjp`` is its transpose. Products of other curve types, of surfaces and
of batch and grid evaluation are not implemented yet; they answer
``NotSupported``.

Products of intersections
-------------------------

An intersection is differentiated through the equations that define its
result, never through the solver that found it (ADR-0021, section 6; the
v0.2 gate of section 9's stage 4). A point where a curve meets a surface
solves :math:`F(z; p) = C(t; p) - S(u, v; p) = 0` in :math:`z = (t, u, v)`:
three equations in three unknowns, so a point needs no tracking constraint.
With :math:`F_z = [C', -S_u, -S_v]`, a design change :math:`\dot p` moves it
by :math:`F_z \dot z = S_p \dot p - C_p \dot p` and its position by
:math:`\dot x = C_p \dot p + C' \dot t`. The determinant
:math:`C' \cdot (S_u \times S_v)` vanishes exactly at a tangency, and
:math:`|\det F_z| / (|C'| |S_u \times S_v|)`, the sine of the angle between
curve and surface, is the conditioning.

``Linearize(line, surface, point)`` exposes this for one point of
``Intersect``: the residual :math:`F`, the columns of :math:`F_z` (``Tangent``,
``Su``, ``Sv``), the conditioning and the status. Whether a derivative exists
is read from what the primal solver proved:

=====================================  ================================================
Point                                  Status
=====================================  ================================================
a crossing the solver decided          ``Valid``, or ``NearEvent`` below ``NearEventSine``
                                       (2\ :sup:`-20`), qualified
a ``Tangential`` point, or an overlap  ``Nondifferentiable``: the points merge, vanish
                                       or stop being isolated as the geometry moves
a point the solver left undecided      ``Indeterminate``
a crossing at a pole of the surface's  ``Indeterminate``: the angle there has no
parameters (a sphere's)                derivative, though the geometry is regular
=====================================  ================================================

A result takes its worst point's status, and carries values only for
``Valid`` and ``NearEvent``.

``IntersectJvp`` and ``IntersectVjp`` implement the products for a line or
a B-spline curve against a plane, a cylinder, a cone, a sphere, a torus or
an ellipsoid.
``LineBinding`` binds a parameter to a coordinate of the line's origin or
direction. ``SurfaceBinding`` binds one to a value of the surface, each a
``Length`` but the semi-angle, an ``Angle``:

- ``OriginX``, ``OriginY``, ``OriginZ``: a coordinate of its origin (a
  sphere's or a torus's centre);
- ``Radius``: a sphere's or a cylinder's radius, a cone's reference radius,
  a torus's major radius;
- ``MinorRadius``: a torus's tube;
- ``SemiAngle``: a cone's opening;
- ``SemiAxisA``, ``SemiAxisB``, ``SemiAxisC``: an ellipsoid's semi-axes,
  along its ``XAxis``, its ``YAxis`` and its axis.

A binding of a value the surface does not have is ``InvalidInput``. The
surface's partials :math:`S_u`, :math:`S_v` come from its own ``Deriv``.
At fixed :math:`(u, v)`, with :math:`\hat\rho = \cos u\,X + \sin u\,Y` and
:math:`A` the axis:

- the origin moves every point alike;
- the radius moves a sphere's point along :math:`(S - C) / r`, and a
  cylinder's, a cone's and a torus's along :math:`\hat\rho`;
- a cone's semi-angle moves it by
  :math:`v (\cos a\,\hat\rho - \sin a\,A)`;
- a torus's minor radius moves it along the tube's normal,
  :math:`\cos v\,\hat\rho + \sin v\,A`;
- an ellipsoid's semi-axes move it by the factor each multiplies,
  :math:`\cos v \cos u\,X`, :math:`\cos v \sin u\,Y` and
  :math:`\sin v\,A`; its poles are detected as the sphere's are.

``Intersect`` has no line handler for the cone and the torus yet, so their
points come from B-spline curves, or from a caller's own closed form, which
``IntersectJvp`` accepts once it lies on both geometries. The cone's apex is
a pole, detected as the sphere's are. A pole is detected on its own: the
conditioning, the sine of an angle, scales out the vanishing :math:`S_u`, so
a partial within rounding of zero against the other (below 2\ :sup:`-50` of
it) makes the point ``Indeterminate``. Each point gives six values,
``CurveSurfacePointValues``: :math:`\dot x` (three), then :math:`\dot t`,
:math:`\dot u` and :math:`\dot v`. The VJP takes an adjoint of six components
per point and returns one value per parameter, accumulated point by point in
order. The 3 × 3 system is solved by Cramer's rule, transposed for the VJP.
Both report the largest residual and the smallest conditioning over the
points, and a signature of the point count and contacts.

The tests check the forward product against the crossing's closed form,
the second-order decay of the Taylor residual
:math:`\lVert y(p + h v) - y(p) - h J v \rVert` on random lines and planes,
that the VJP is the JVP's transpose, and every status: an overlap, a touch,
an undecided point and a grazing crossing never carry an unqualified
derivative. The sphere and the cylinder are checked against the closed form
of a line crossing a circle, :math:`x' = r / x` for a growing radius, and
the Taylor and transpose properties run on random lines through both.

A ``BSplineCurve`` has the same products, its poles and weights moved by
``PoleBinding`` exactly as ``EvalJvp`` moves them (a coordinate by a
``Length``, a weight by a ``Dimensionless`` parameter, at fixed knots). The
curve's side of the system is its derivative :math:`C'(t)` and
``PoleVelocity(curve, binding, t)``, the velocity of the point at fixed
:math:`t` that ``EvalJvp`` sums; the rest is written once for every curve.
A degree-1 B-spline from :math:`P_0` to :math:`P_1` is checked against the
line from :math:`P_0` along :math:`P_1 - P_0`, moving :math:`P_1` against
moving the line's direction, through every surface; the Taylor and
transpose properties run on random cubic B-splines, rational half the
time with weight bindings, through spheres, cylinders and planes; and a
B-spline touching a plane, the undecided touch the B-spline handler
reports, is ``Nondifferentiable``.

A ``Circle`` and an ``Ellipse`` have the same products, the point's
parameter their angle, their defining values moved by ``ConicBinding``: a
coordinate of the centre, a circle's ``Radius``, an ellipse's
``MajorRadius`` and ``MinorRadius``, each a ``Length``, the frame fixed as a
surface's is. At the angle :math:`t`, :math:`C = O + r(\cos t X + \sin t
Y)` moves by a unit with a centre coordinate and by :math:`\cos t X + \sin
t Y` with its radius; an ellipse :math:`O + a \cos t X + b \sin t Y` by
:math:`\cos t X` with :math:`a` and :math:`\sin t Y` with :math:`b`. The
tests check a circle of radius 2 and an ellipse of semi-axes 3 and 2, both
cut by the plane :math:`z = 1`, against their closed forms, and the Taylor
and transpose properties on random circles and ellipses through spheres,
cylinders and planes. A conic meets a quadric up to four times, so those
properties discard a step that changes the number of crossings (a near
miss elsewhere becoming two crossings), pair crossings by position (one
near :math:`t = 0` may pass the seam), and keep away from a sphere's poles,
where :math:`u` moves too fast for the steps.

A ``BSplineSurface`` binds its poles and weights the same way: a
``SurfaceBinding`` with ``SurfaceComponent::PoleX``, ``PoleY`` or ``PoleZ``
(a ``Length``) or ``PoleWeight`` (``Dimensionless``) and the ``Pole`` it
moves, the index into ``Poles``, row by row. The surface's side of the
system is ``PoleVelocity(surface, pole, component, uv)``, the velocity of
:math:`S(u, v)` at fixed :math:`(u, v)`: :math:`R_{ij} e_c` for a coordinate
and :math:`N_i N_j (P_{ij} - S) / W` for a weight. It is checked against the
derivative of evaluating the surface over dual numbers, every component of
random polynomial and rational surfaces. The products are checked against a
closed form, a vertical line through a flat bilinear patch rising with a
corner pole by :math:`R_{11}(u, v)`, and the Taylor and transpose
properties run on random lines through random bicubic patches, rational half
the time, with three poles bound. Other surfaces (extrusions, revolutions)
answer ``NotSupported``; each new curve–surface handler binds its
geometry's parameters to the same system.

Points of two plane curves
~~~~~~~~~~~~~~~~~~~~~~~~~~

Two plane curves meet where :math:`F(s, t) = A(s) - B(t) = 0`: two
equations in two unknowns, with :math:`F_z = [A'(s), -B'(t)]`, singular
exactly where the curves are tangent. Its conditioning is the sine of the
angle between them, :math:`|A' \times B'| / (|A'| |B'|)`. A design change
moves the parameters by :math:`F_z z' = B_p p' - A_p p'` and the point by
:math:`x' = A_p p' + A' s'`; the reverse product solves the transposed
system once per point. ``Linearize(a, b, point)``, ``IntersectJvp`` and
``IntersectVjp`` take the two ``Curve2d`` and a ``PlaneCurveBinding`` per
parameter: a line's origin and direction, a circle's centre and radius, an
ellipse's or a hyperbola's centre and radii, a parabola's vertex and focal
length, or a B-spline's poles and weights, whose velocity is the curve
library's ``PoleVelocity`` on the curve lifted into the plane
:math:`z = 0`. Every frame is held fixed. At fixed parameter the
hyperbola :math:`C + a \cosh t\, X + b \sinh t\, Y` moves by
:math:`\cosh t\, X` per unit of ``MajorRadius`` and :math:`\sinh t\, Y`
per unit of ``MinorRadius``; the parabola
:math:`V + t^2 / (4 f)\, X + t\, Y` by :math:`-t^2 / (4 f^2)\, X` per unit
of ``FocalLength``, its vertex bound through ``CenterX`` and ``CenterY``.
The statuses are a curve and a surface's: a decided crossing ``Valid``, or
``NearEvent`` below ``NearEventSine``; a tangential contact or an overlap
``Nondifferentiable``; an undecided one ``Indeterminate``. Every built-in
plane curve binds; ``NotSupported`` is kept for a kind no binding reaches.

The tests check the x-axis against a circle, an ellipse and a segment
against their closed forms, and the Taylor and transpose properties on
random pairs of lines, circles, ellipses, hyperbolas, parabolas and cubic
B-splines, rational half the time, pairing points by position across a
circle's seam and discarding a step that changes their number or that the
primal solver does not answer. Four closed forms cover the open conics,
each as either curve of its pair. The line :math:`y = 4` cuts
:math:`y = x^2` (:math:`f = 1/4`) at :math:`(\mp 2, 4)`, :math:`t = \pm 2`:
raising the vertex by :math:`v` gives :math:`t = \pm 2 \sqrt{f (4 - v)}`,
so :math:`x' = \pm 1/4`, and a longer focal length
:math:`dt/df = \pm \sqrt{4 / f} = \pm 4`, so :math:`x' = \mp 4`. The line
:math:`x = 2` meets the same parabola at :math:`(2, 4)`, :math:`t = -2`,
where the focal length moves the height by :math:`-t^2 / (4 f^2) = -16` and
the vertex's :math:`x` by :math:`-2 (x - u) = -4`. The line :math:`x = 5`
cuts the rectangular hyperbola of semi-axes 3 at :math:`(5, \pm 4)`,
:math:`\cosh t = 5/3`: :math:`a \cosh t = 5` gives
:math:`t' = -5 / (a^2 \sinh t) = \mp 5/12` and
:math:`y' = b \cosh t\, t' = \mp 25/12` per unit of :math:`a`,
:math:`y' = \sinh t = \pm 4/3` per unit of :math:`b`. The line :math:`y = 4` meets it once, at :math:`(5, 4)`,
:math:`\sinh t = 4/3`: :math:`x' = \cosh t = 5/3` per unit of :math:`a`
and :math:`t' = -4 / (b^2 \cosh t) = -4/15`, :math:`x' = -16/15` per unit of
:math:`b`. Each new partial's sign, flipped, fails these. A parabola's and a
hyperbola's tangencies with a line (their exact discriminants vanish) and
each against itself are ``Nondifferentiable``; a line a hair inside each
vertex, ``NearEvent``.

Two space curves meet only where their closest approach is within the
resolution, which a generic change of either curve undoes: a meeting is
itself an event, of codimension one. ``IntersectJvp(a, b, intersection,
context)`` and ``IntersectVjp`` for two ``Curve`` therefore answer
``Nondifferentiable`` for a result with a point or an overlap, saying why,
and ``Valid`` with no values (zeros, in reverse) for a proved separation,
which a small enough change keeps (ADR-0052, the owner's decision of
2026-11-10).
Tracking the closest approach smoothly, :math:`(A - B) \cdot A' = 0` and
:math:`(A - B) \cdot B' = 0` in :math:`(s, t)`, is a separate product,
``ClosestApproach``, still to come, with its own tracking convention; it is
not hidden inside ``Intersect``'s derivatives.

Branch points of two surfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Where two surfaces meet in a curve, a point has four unknowns,
:math:`z = (u_A, v_A, u_B, v_B)`, and only three equations,
:math:`S_A(u_A, v_A) - S_B(u_B, v_B) = 0`: it can slide along the branch.
The fourth equation is a tracking constraint, and Stereon uses the normal
plane:

.. math::

   T \cdot \bigl(S_A(u_A, v_A) - X_0\bigr) = 0,

where :math:`X_0` and :math:`T` are the branch's point and unit tangent at
the sample's parameter :math:`t`, the shared edge parameter of ADR-0009. As
the geometry moves, the point stays in the plane through :math:`X_0` normal
to the branch, so :math:`T \cdot \dot x = 0`. The constraint is geometric:
it does not depend on either surface's parametrisation, unlike holding one
pcurve coordinate fixed, and it is well conditioned wherever the branch is
transverse. With the rows
:math:`[S_{A,u}, S_{A,v}, -S_{B,u}, -S_{B,v}]` and
:math:`[T \cdot S_{A,u}, T \cdot S_{A,v}, 0, 0]`,

.. math::

   F_z \dot z = \begin{pmatrix} S_{B,p} \dot p - S_{A,p} \dot p \\
                -T \cdot S_{A,p} \dot p \end{pmatrix}, \qquad
   \dot x = S_{A,p} \dot p + S_{A,u} \dot u_A + S_{A,v} \dot v_A .

:math:`F_z` is singular where the surfaces are tangent or where the branch
lies in its own tracking plane. The conditioning is
:math:`|T \cdot (\hat n_A \times \hat n_B)|`, the sine of the angle between
the surfaces when :math:`T` follows the branch.

A branch's curves are fitted, within its certified ``Bound``, so a sample
is first refined onto both surfaces. Newton's method on the four equations
starts from the branch's pcurves at :math:`t`. The refined point must lie on
both surfaces and in the tracking plane to rounding. It must also lie
within :math:`2(\text{Bound} + \text{rounding}) / \sin\theta` of
:math:`X_0`; otherwise the result came from other geometry and is
``InvalidInput``.

``Linearize(a, b, branch, t)`` exposes the refined point, its parameters,
the residual, :math:`T`, the four partials, the conditioning and the
status. ``IntersectJvp`` and ``IntersectVjp`` take the result, a
``SurfaceBinding`` list for each surface and a list of ``BranchSample``
values (a branch's index and a parameter in its range). One parameter may
move both surfaces. Each sample gives seven values, ``BranchPointValues``:
:math:`\dot x` (three), then :math:`\dot u_A`, :math:`\dot v_A`,
:math:`\dot u_B` and :math:`\dot v_B`. The 4 × 4 system is solved by
Gaussian elimination with partial pivoting, transposed for the VJP.

=========================================  ============================================
Sample                                     Status
=========================================  ============================================
on a certified ``Crossing`` branch         ``Valid``, or ``NearEvent`` below
                                           ``NearEventSine``
on a ``Tangential`` or ``Overlap`` branch  ``Nondifferentiable``: the branch splits or
                                           vanishes as the surfaces move
any sample of a result with an overlap     ``Nondifferentiable``
on an uncertified branch                   ``Indeterminate``
at a pole of either surface's parameters   ``Indeterminate``
=========================================  ============================================

Both surfaces must be planes, cylinders, cones, spheres, tori,
ellipsoids or B-spline surfaces, whose poles and weights bind as above, so
the branches of the generic SSI have products too; other pairs answer
``NotSupported``. Two crossing bilinear patches, one flat, are checked
against the closed form of raising the flat one's corner pole. The tests check the forward product against
a closed form. A plane at height :math:`h` cuts a sphere of radius
:math:`R` in a circle of radius :math:`\rho = \sqrt{R^2 - h^2}`, whose
point moves radially by :math:`(R \dot R - h \dot h) / \rho` and up by
:math:`\dot h`. They also check the second-order Taylor decay on random
planes through spheres, planes across cylinders and cylinders through
spheres, the last one traced and certified, not exact. The moved point is
found in the same tracking plane and held there, and on both moved
surfaces, to rounding. Further tests check the transpose identity, the
statuses of a tangential branch, an overlap and an uncertified branch, and
the refusal of stale results.

The AD backend: CppAD
---------------------

CppAD is Stereon's production AD backend (ADR-0045). ``CppAD::AD<double>``
satisfies ``Stereon::Scalar``, so the scalar-generic evaluators take it
unchanged: poles, weights and transform parameters become active, while knots
and parameters stay ``double``. Their branches depend on knots and parameters
only, so a tape recorded once holds for every value of the active inputs.

.. code-block:: cpp

   std::vector<CppAD::AD<double>> x = /* poles, then weights */;
   CppAD::Independent(x);
   Stereon::BasicBSplineCurve<CppAD::AD<double>, 3> active{curve.Degree, curve.Knots, {}, {}};
   // ... fill active.Poles and active.Weights from x ...
   std::vector<CppAD::AD<double>> y = /* x, y, z of Eval(active, t) for each t */;
   CppAD::ADFun<double> f(x, y);
   const std::vector<double> jacobian = f.Jacobian(x0);

CppAD is EPL-2.0, so under ADR-0008 it is an optional adapter: no kernel
library links or includes it, and ``tools/check_dependencies.py`` enforces
that. Configuring with ``-DSTEREON_WITH_CPPAD=ON`` fetches the pinned release
and builds ``stereon_cppad_tests`` (``tests/ad/CppADTests.cpp``), which checks
CppAD's Jacobians of curve and surface points and of a composed rigid motion
against closed forms and against forward-mode duals, and checks that
``EvalJvp`` and ``EvalVjp`` equal CppAD's forward and reverse modes. An application that
tapes Stereon's evaluators links CppAD itself, under CppAD's licence.
