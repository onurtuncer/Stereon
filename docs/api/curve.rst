.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Curves — ``Stereon.Curve``
==========================

Curves in space and in a surface's parameter plane. There are the
analytic curves (``Line``, ``Circle``, ``Ellipse``, ``Hyperbola`` and
``Parabola``), B-spline and NURBS curves (``BSplineCurve``), trimmed curves,
and the type-erased ``Curve``, which holds any of them. Every curve type
also comes in a 2D form (``Line2d``, ``Circle2d``, ...) for the pcurves that
describe an edge in the parameter plane of a face's surface.

Around the types are the operations a modeller expects:

- evaluation and derivatives of any order, also in batches;
- closest-point projection;
- exact conversion of any analytic curve to NURBS;
- the B-spline toolkit: knot insertion and removal, degree elevation,
  splitting, Bezier decomposition, the periodic form;
- fitting curves to points;
- local properties: tangent, curvature, torsion, the Frenet frame, and arc
  length.

.. code-block:: cpp

   import Stereon.Curve;

   using namespace Stereon;

   const Circle rim{Point3{0, 0, 0}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
   const Point3 top = Eval(rim, std::numbers::pi / 2);   // (0, 5, 0)
   const double t = Project(rim, Point3{0, -9, 3});      // 3 pi / 2

   const Curve edge(rim);                  // any curve, behind one value type
   if (const Circle* c = edge.Target<Circle>()) { /* the circle itself */ }

   const Result<ArcLengthResult> half = ArcLength(edge, Interval{0.0, std::numbers::pi});
   // half->Length is 5 pi; half->ErrorEstimate is the quadrature's own estimate

Start here
----------

Pick the type for your curve
   ``Line``, ``Circle``, ``Ellipse``, ``Hyperbola``, ``Parabola``,
   ``BSplineCurve`` or ``TrimmedCurve``. They are plain aggregates: a
   ``Circle`` is a centre, a radius and two axes. For data from outside the
   kernel use the checked factory (``MakeCircle``, ``MakeBSplineCurve``,
   ...), which returns ``InvalidInput`` with a reason instead of a broken
   curve.

The same operations work on every curve
   ``Eval``, ``Deriv``, ``Domain``, ``BoundingBox``, ``Project``,
   ``Reversed``, ``Transformed`` and ``ToBSpline``. The *Operations at a
   glance* table in the reference lists them with the types each is
   declared for.

``Curve`` and ``Curve2d``: any curve in one variable
   Copyable values that hold any curve type, built in or your own. A call
   through them costs one virtual call; ``VisitBuiltin`` runs a hot loop on
   the concrete type instead.

B-splines
   ``InsertKnot``, ``RemoveKnot``, ``ElevatedDegree``, ``Split``,
   ``ToBezierSegments``, ``Clamped``, ``AsPeriodic``, ``Concatenated``. For
   many evaluations at once: ``EvalAll``, which matches ``Eval`` to the bit,
   and ``EvalAllFast``, which is faster within a tested bound.

Fitting
   ``Interpolate``, ``FitLeastSquares`` and ``Approximate``, each also
   periodic. Every fit returns the curve together with the deviation it
   achieved.

Properties
   ``UnitTangent``, ``Curvature``, ``Torsion``, ``FrenetFrame``,
   ``ArcLength``, ``ParameterAtArcLength``, for any curve type.

The guide chapters behind this library: :doc:`/geometry` (the type-erasure
design), :doc:`/analytic`, :doc:`/nurbs`, :doc:`/trimmed`,
:doc:`/interchange`, :doc:`/fitting`, :doc:`/properties` and
:doc:`/sensitivity`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: curve

How the types fit together
--------------------------

There is no class hierarchy. A curve type is any type with the free
functions the ``CurveGeometry`` concept asks for (``Curve2dGeometry`` in the
plane), and ``Curve`` holds any such type by value (ADR-0003). The six
built-in types are a closed list that ``VisitBuiltin`` and the intersection
tables dispatch over without a virtual call; any other type goes through
the virtual path. Each built-in name is an alias of a class template written
once for both dimensions.

.. stereon-type-map:: curve
   :caption: The wrappers (blue), the built-in types they dispatch over
             (solid arrows), the class templates behind them (dotted) and
             the concepts any other type must model (dashed).

Inside the module
-----------------

.. stereon-partition-graph:: curve

Reference
---------

.. stereon-api:: curve
