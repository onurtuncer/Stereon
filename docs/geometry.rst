.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Curves and surfaces
===================

Analytic, NURBS and procedural curves and surfaces, the geometry concepts, and the type-erased Curve and Surface wrappers.

Vectors, points and transforms
------------------------------

``Stereon.Core`` provides the fixed-size types that every geometry type is
built from: ``Vec2``, ``Vec3``, ``Point2``, ``Point3``, the square matrices
``Mat2``, ``Mat3`` and ``Mat4``, and the affine ``Transform``. All of their
operations are ``constexpr``.

Points are positions and vectors are displacements, so ``Point - Point`` is a
``Vec``, ``Point + Vec`` is a ``Point``, and ``Point + Point`` does not
compile. ``ToVec`` and ``ToPoint`` convert explicitly.

.. code-block:: cpp

   import Stereon.Core;

   using namespace Stereon;

   constexpr Point3 corner{10, 0, 0};
   constexpr Vec3 up{0, 0, 1};

   // (a * b) applies b first: rotate about Z, then lift by 5 mm.
   const Transform lift = MakeTranslation(5.0 * up) * MakeRotation(up, 0.5);
   const Point3 moved = Transformed(corner, lift);
   const Vec3 direction = Transformed(up, lift); // vectors ignore the translation

Each type is the ``double`` instance of a scalar-generic template
(``BasicVec3<T>``, ``BasicMat<T, N>``, ``BasicTransform<T>``) constrained by
the ``Scalar`` concept, so the same code runs on the active scalar type of an
automatic-differentiation backend (ADR-0021). Generic code keeps the scalar
type: constants are written ``T{2}``, and mathematical functions are called
through ADL, ``using std::sqrt; sqrt(x);``, so that an active type can supply
its own. Lengths are millimetres and angles radians.

``Determinant`` and ``Inverse`` are closed-form cofactor expansions, so their
results do not depend on a pivoting strategy (ADR-0007). A singular matrix or
a zero vector passed to ``Normalized`` is a contract violation, not a
``KernelError``: code that may meet one checks first.

Every type formats with ``std::format``: ``"{}"`` gives ``(1, 2, 3)`` for a
vector or point and ``[[1, 0], [0, 1]]`` for a matrix.

Curves and surfaces
-------------------

.. todo::

   Write this section. It is a placeholder for the documentation effort.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_

