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
``KernelError``: code that may meet one checks first. These primitives run in
inner loops, so their contracts are hot-path contracts (ADR-0016), compiled in
debug and sanitizer builds and left out of release builds.

Every type formats with ``std::format``: ``"{}"`` gives ``(1, 2, 3)`` for a
vector or point and ``[[1, 0], [0, 1]]`` for a matrix.

Intervals and bounding boxes
----------------------------

``Interval`` is a closed interval of doubles with outward rounding: every
operation returns an interval that contains every real result of the
operation on real numbers from its operands. It is the enclosure arithmetic
behind the predicate filters of ADR-0002. Rounding steps each bound one
representable value outward on the bit pattern, so the arithmetic is
``constexpr``, carries no rounding-mode state, and gives the same bits on
every platform (ADR-0007); the price is that an exact result is also widened
by one step.

.. code-block:: cpp

   const Interval sum = Interval{0.1} + Interval{0.2};
   assert(0.1 + 0.2 != 0.3);      // double arithmetic rounds
   assert(Contains(sum, 0.3));    // the enclosure does not lose it

``Interval{}`` is the empty interval, the identity of ``Hull``; ``Squared``
is tighter than ``i * i`` when the interval contains zero; and division by an
interval containing zero is a contract violation rather than a result.

``Box2`` and ``Box3`` are products of one ``Interval`` per axis (``box.X``,
``box.Y``, ``box.Z``). ``Contains``, ``Intersects``, ``Hull`` and
``Intersection`` are the per-axis interval operations, so a box is empty as
soon as one axis is, and ``Box3{}`` is the empty box. Boxes are closed:
touching faces intersect. ``Box3::FromCorners`` accepts two opposite corners
in either order. Both intervals and boxes are ``double`` only; they bound and
decide, which stays discrete (ADR-0021).

Curves and surfaces
-------------------

Geometry types are plain aggregates with no base class (ADR-0003). What a
type can do is given by free functions found by argument-dependent lookup,
and three concepts name the set a type must provide:

- ``CurveGeometry``, a curve in space: ``Eval(c, t)``, ``Deriv(c, t, n)``,
  ``Domain(c)``, ``BoundingBox(c)``, ``Project(c, p)``, ``Reversed(c)`` and
  ``Transformed(c, x)``;
- ``Curve2dGeometry``, a curve in a surface's parameter plane (a pcurve,
  ADR-0009): the same without ``Transformed``, since there is no 2D transform
  type yet;
- ``SurfaceGeometry``: ``Eval(s, uv)``, ``Deriv(s, uv, du, dv)``, ``Domain``
  (a ``Box2``), ``BoundingBox``, ``Project`` (returning the parameter point),
  ``Reversed`` and ``Transformed``.

These are tier 1 of ADR-0003. The bounding operations of tier 2 and the
intersection capability of tier 3 arrive with the intersection work.

The built-in types so far are ``Line`` and ``Circle`` (and ``Line2d`` and
``Circle2d`` in the plane), ``Plane`` and ``Cylinder``. A line is
``Origin + t * Direction``; a circle is ``Center + Radius * (cos t * XAxis +
sin t * YAxis)`` on ``[0, 2π]``; a plane is ``Origin + u * XAxis + v *
YAxis``; a cylinder is the circle swept along its axis. Like the core math
types they are instances of templates over the ``Scalar`` concept, so
``Eval``, ``Deriv``, ``Reversed`` and ``Transformed`` also run on the active
scalar type of an automatic-differentiation backend (ADR-0021).
``Domain``, ``BoundingBox`` and ``Project`` bound and decide, and stay
``double`` like ``Interval``.

.. code-block:: cpp

   import Stereon.Surface;

   using namespace Stereon;

   const Circle rim{Point3{0, 0, 0}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
   const Point3 top = Eval(rim, std::numbers::pi / 2);   // (0, 5, 0)
   const double t = Project(rim, Point3{0, -9, 3});      // 3 pi / 2

   const Cylinder wall{Point3{0, 0, 0}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
   const Vec3 normal = Cross(Deriv(wall, {0.0, 2.0}, 1, 0),
                             Deriv(wall, {0.0, 2.0}, 0, 1)); // points outward

Some conventions worth knowing:

- ``Deriv`` of order n of a circle turns the radius vector by n quarter
  turns, selected by ``n mod 4``, so that no rounded multiple of π/2 enters
  the argument.
- ``Project`` on a circle or cylinder returns an angle in ``[0, 2π)``; a point
  on the axis, which is equally close to every point of the circle, gets 0.
- ``Reversed`` runs a curve the other way and turns a surface's normal round:
  for a cylinder it runs u backwards, which is why ``Cylinder`` carries a
  ``RightHanded`` flag and its frame may have either handedness.
- ``Transformed`` accepts similarities: rigid motions, reflections and
  uniform scaling. A non-uniform scaling would turn a circle into an
  ellipse, which the type cannot hold.
- Circle bounding boxes are computed in ``Interval`` arithmetic, so they
  enclose the exact circle. Lines, planes and cylinders are unbounded and
  return ``Everything()``.

Type-erased wrappers
~~~~~~~~~~~~~~~~~~~~

``Curve``, ``Curve2d`` and ``Surface`` hold any type satisfying the matching
concept behind one copyable value type, and satisfy the concept themselves, so
code written against the concept works on either. A wrapped value smaller
than the 96-byte buffer, which every built-in type is, lives inside the
wrapper and allocates nothing; a larger one, such as a NURBS curve with many
poles, is held immutable behind a shared pointer and shared between copies.
``Target<T>()`` returns the wrapped value when it is a ``T``.

.. code-block:: cpp

   const Curve edge(rim);                       // explicit: wrapping is a choice
   const Point3 p = Eval(edge, 1.0);            // one virtual call
   if (const Circle* c = edge.Target<Circle>()) { /* the circle itself */ }
   const Curve back = Reversed(edge);           // still a Circle inside

The wrappers are ``double`` only and are ``std::semiregular``, so the paged
stores of ADR-0004 can hold them.

The fast path
~~~~~~~~~~~~~

A call through ``Curve`` costs one virtual call, which is a lot next to
evaluating a line. Hot loops therefore dispatch once with
``VisitBuiltin(curve, f)``: when the curve wraps one of the built-in types
(``BuiltinCurveTypes``, ``BuiltinCurve2dTypes``, ``BuiltinSurfaceTypes``),
``f`` receives that type and its loop makes direct calls; any other type
reaches ``f`` as the wrapper, whose calls stay virtual. The type test is
generated with an expansion statement and pack indexing over the type list
(ADR-0016) and runs once per ``VisitBuiltin``.

.. code-block:: cpp

   const double length = VisitBuiltin(edge, [&](const auto& g) {
       double total = 0.0;
       for (std::size_t i = 1; i < samples.size(); ++i) {   // runs on Circle itself
           total += Length(Eval(g, samples[i]) - Eval(g, samples[i - 1]));
       }
       return total;
   });

Every call of ``f`` must return the same type. The two paths run the same
code, and property tests require them to agree to the bit on every built-in
type. ``BuiltinIndexOf`` gives a wrapped value's position in its type list,
which the intersection tables of ``Stereon.Intersect`` use: ``Intersect(a,
b, ctx)`` dispatches built-in pairs through a table generated at compile
time, one handler per ordered pair. Every handler is a stub returning
``NotImplemented`` until the analytic handlers of ADR-0010 arrive, and types
that are not built in return ``NotSupported``.

Measured in ``bench/GeometryBench.cpp`` (GCC 16.2, -O3), in nanoseconds per
evaluation:

===================================  ======  ======
Path                                 Line    Circle
===================================  ======  ======
Direct call on the concrete type     2.2     49
``VisitBuiltin``, loop inside        2.0     46
``VisitBuiltin`` around every call   6.4     50
``Curve``, one virtual call each     5.6     50
===================================  ======  ======

``VisitBuiltin`` meets ADR-0003's target of 10% over direct calls. The
virtual path misses its target of 2x for a line, where the evaluation costs
about as much as the call; for a circle the trigonometry dominates. The copy
benchmark settles the buffer size: at 64 bytes a circle no longer fits and a
copy costs 42 ns instead of 12, while 128 bytes gains nothing, so
``GeometryBufferSize`` stays at 96, which also holds ``Cylinder`` (88 bytes).

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_

