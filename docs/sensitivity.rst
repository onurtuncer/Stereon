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
Jacobian; none yet). Each library keeps a catalogue in its ``Capabilities``
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
