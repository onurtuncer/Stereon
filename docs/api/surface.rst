.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Surfaces — ``Stereon.Surface``
==============================

Surfaces, organised the way :doc:`curve` organises curves. There are the
analytic surfaces (``Plane``, ``Cylinder``, ``Cone``, ``Sphere`` and
``Torus``), B-spline and NURBS surfaces (``BSplineSurface``), extrusions and
revolutions of any profile curve, rectangular trimmed surfaces, and the
type-erased ``Surface``, which holds any of them.

A surface maps a parameter point ``(u, v)`` to a point in space.
``Deriv(s, uv, i, j)`` is the partial derivative of order ``i`` in u and
``j`` in v. On top of evaluation come:

- projection;
- exact conversion to NURBS;
- normals, fundamental forms and principal, Gaussian and mean curvatures;
- grid evaluation and fitting of B-spline surfaces;
- interval certificates that bound how far a curve strays from a surface.

.. code-block:: cpp

   import Stereon.Surface;

   using namespace Stereon;

   const Cylinder wall{Point3{0, 0, 0}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
   const Vec3 outward = Cross(Deriv(wall, {0.0, 2.0}, 1, 0),
                              Deriv(wall, {0.0, 2.0}, 0, 1));

   const Surface face(wall);               // any surface, behind one value type

Start here
----------

Pick the type for your surface
   ``Plane``, ``Cylinder``, ``Cone``, ``Sphere``, ``Torus``,
   ``BSplineSurface``, ``ExtrusionSurface``, ``RevolutionSurface`` or
   ``TrimmedSurface``. Data from outside the kernel goes through the
   checked factory (``MakeCylinder``, ``MakeBSplineSurface``,
   ``MakeExtrusion``, ...).

The same operations work on every surface
   ``Eval``, ``Deriv``, ``Domain``, ``Periods``, ``BoundingBox``,
   ``Project``, ``Reversed``, ``Transformed`` and ``ToBSpline``.

``Surface``: any surface in one variable
   As ``Curve`` is for curves: a copyable value with a fast path through
   ``VisitBuiltin``.

Properties
   ``UnitNormal``, ``FundamentalForms``, ``PrincipalCurvatures``,
   ``GaussianCurvature``, ``MeanCurvature``. At a singular point (a cone's
   apex, a sphere's pole) they return ``DegenerateResult`` rather than a
   meaningless value.

B-spline surfaces
   ``InsertKnot``, ``ToBezierPatches``, ``IsoCurve``, ``EvalGrid``, and the
   fits ``InterpolateSurface``, ``FitSurfaceLeastSquares`` and
   ``ApproximateSurface``.

Certificates
   ``CertifyDeviation`` bounds the distance from a curve to a surface over
   a whole parameter range, not just at samples.

The guide chapters behind this library: :doc:`/geometry`,
:doc:`/analytic`, :doc:`/nurbs`, :doc:`/swept`, :doc:`/trimmed`,
:doc:`/fitting` and :doc:`/properties`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: surface

How the types fit together
--------------------------

As for curves, there is no class hierarchy: a surface type is any type with
the free functions of the ``SurfaceGeometry`` concept, and ``Surface``
holds any such type by value. The six built-in types are dispatched
without a virtual call. Extrusions, revolutions and trimmed surfaces wrap
other curves and surfaces and travel the general path.

.. stereon-type-map:: surface
   :caption: The wrapper (blue), the built-in types it dispatches over
             (solid arrows), the class templates behind them (dotted) and
             the concept any other type must model (dashed).

Inside the module
-----------------

.. stereon-partition-graph:: surface

Reference
---------

.. stereon-api:: surface
