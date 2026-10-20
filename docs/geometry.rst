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

Every transform that fixes a point is a linear map conjugated by a
translation, :math:`x \mapsto A(x - p) + p`, which ``AboutPoint(A, p)``
builds. The helpers that OCCT's ``gp_Trsf`` offers are built on it:
``MakeRotation(point, axis, angle)`` turns about an arbitrary line,
``MakeScaling(point, factor)`` scales about a point, and the three mirrors
reflect through a point (:math:`A = -I`), in a line
(:math:`A = 2\mathbf k\mathbf k^\mathsf{T} - I`, the half-turn about it) and
in a plane (:math:`A = I - 2\mathbf n\mathbf n^\mathsf{T}`):
``MakePointMirror(point)``, ``MakeAxisMirror(point, axis)`` and
``MakePlaneMirror(point, normal)``. The axis and the normal need not be unit
vectors. The point and plane mirrors reverse orientation (determinant
:math:`-1`); ``Transformed`` takes all of them.

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

The built-in analytic types are ``Line``, ``Circle``, ``Ellipse``,
``Hyperbola`` and ``Parabola`` (and their ``2d`` forms in the plane), and
``Plane``,
``Cylinder``, ``Cone``, ``Sphere`` and ``Torus``. A line is
``Origin + t * Direction``; a circle is ``Center + Radius * (cos t * XAxis +
sin t * YAxis)`` on ``[0, 2π]``; a plane is ``Origin + u * XAxis + v *
YAxis``; a cylinder is the circle swept along its axis. :doc:`analytic`
gives every parametrisation with its derivatives, singular points and
projection. Free-form geometry is ``BSplineCurve``, ``BSplineCurve2d`` and
``BSplineSurface``, polynomial or rational (NURBS), described in
:doc:`nurbs`. Like the core math
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
- Bounding boxes of the bounded types (circle, ellipse, sphere, torus) are
  computed in ``Interval`` arithmetic, so they enclose the exact geometry.
  Lines, planes, cylinders and cones are unbounded and return
  ``Everything()``.
- Data from outside the kernel goes through a checked factory,
  ``MakeCircle``, ``MakeCylinder`` and so on, which returns ``InvalidInput``
  instead of building an invalid value (:doc:`analytic`, *Checked
  construction*).

Type-erased wrappers
~~~~~~~~~~~~~~~~~~~~

``Curve``, ``Curve2d`` and ``Surface`` hold any type satisfying the matching
concept behind one copyable value type, and satisfy the concept themselves, so
code written against the concept works on either. A wrapped value smaller
than the 128-byte buffer, which every built-in analytic type is, lives inside the
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
copy costs 42 ns instead of 12, while 128 bytes gains nothing for these
types, so ``GeometryBufferSize`` was set to 96, which also holds
``Cylinder`` (88 bytes), until ``Cone`` and ``Torus`` arrived (see *The
buffer size, revisited*).

Evaluator preconditions
~~~~~~~~~~~~~~~~~~~~~~~

``Eval``, ``Deriv`` and ``Project`` require finite parameters and points and
derivative orders that exist. The erased wrappers also require a parameter
inside ``Domain`` unless the geometry is periodic in it: a type that is
periodic provides ``Period(c)`` for curves or ``Periods(s)`` for surfaces
(``PeriodicCurve``, ``PeriodicSurface``), and every finite parameter is then
valid for it. ``Circle`` and ``Ellipse`` have period 2π; ``Cylinder``, ``Cone`` and
``Sphere`` periods (2π, 0); ``Torus`` periods (2π, 2π).

The evaluators are templates, and GCC 16.2 silently drops contracts in
templates imported from a module, so these checks are not ``pre``
contracts. They are plain branches that report through a function defined
in a core implementation unit, ``Stereon::Detail::EvaluatorViolation``, which
prints the function and the problem and raises a contract violation in every
importer. They are hot-path checks (ADR-0016): compiled in the debug and
sanitizer presets, compiled out in release builds, where evaluators run in
inner loops. The wrappers check at the boundary to user geometry, so a user
type need not check anything itself.

Serialisation
~~~~~~~~~~~~~

``Serialized(value)``, ``Deserialized<T>(bytes)`` and ``ToJson(value)`` are
generated by static reflection: they walk a struct's members in declaration
order with ``std::meta::nonstatic_data_members_of``, recursively, down to
arithmetic values and enumerations, with no hand-written code per type
(ADR-0012, ADR-0016). The binary form is little-endian, doubles as their
IEEE 754 bit pattern, bools as one byte and a ``std::vector`` as its element
count (8 bytes) followed by its elements, so a ``Cylinder`` is 81 bytes, a
B-spline carries its knots, poles and weights, and both read back bit for
bit. Truncated or overlong input, a bool byte other than 0 or 1 and an
element count larger than the remaining bytes could hold are a
``ParseError``, because file contents are data (ADR-0006); the count is
checked before anything is allocated.
The JSON view names every member, for debugging and for the ``stn-dump``
tool of ADR-0012:

.. code-block:: json

   {"Origin":{"X":1,"Y":2,"Z":3},"Direction":{"X":0.5,"Y":0,"Z":-1}}

The geometry records of the native ``.stn`` format are these bytes; the
chunks, checksums, versioning and the stored numerical policy around them
are the ``io`` library's (:doc:`data_exchange`, *Saving and loading
geometry*).
Compilers without reflection get the same functions reporting
``NotSupported``.

The buffer size, revisited
~~~~~~~~~~~~~~~~~~~~~~~~~~

``Cone`` and ``Torus`` take 96 bytes, 104 with the model's virtual table
pointer, so at the original 96-byte buffer they were shared rather than
stored inline. The copy benchmark settled it (GCC 16.2, -O3, nanoseconds per
wrapper copy):

========  =======  ========
Type      96 B     128 B
========  =======  ========
Line      10.5     12.6
Circle    12.4     15.0
Cylinder  13.1     13.5
Cone      45.1     17.3
Torus     42.1     15.0
========  =======  ========

A shared payload costs an atomic reference count per copy, almost three
times an inline copy of the same type, against about 2.5 ns more for the
types that already fit. The buffer was therefore 128 bytes, and every
built-in analytic type was stored inline.

Topology later made the size matter far more than a copy does: a shape keeps
one erased curve per edge and one erased pcurve per coedge, and at 144 bytes
each (128 plus the model pointer, aligned to 16) they were about two thirds
of its memory (B0.4, :doc:`performance`). Each kind of wrapper now has its
own buffer, sized to its largest built-in type plus the table pointer, at the
alignment of a double:

.. list-table::
   :header-rows: 1

   * - Wrapper
     - Largest built-in type
     - Buffer
     - Wrapper size
   * - ``Curve``
     - ``Ellipse``, ``Hyperbola`` (88 bytes)
     - 96
     - 104
   * - ``Curve2d``
     - ``Ellipse2d``, ``Hyperbola2d`` (64 bytes)
     - 72
     - 80
   * - ``Surface``
     - ``Cone``, ``Torus`` (96 bytes)
     - 104
     - 112

Every built-in analytic type is still stored inline, so no copy got slower;
``static_assert``\ s next to each type and in
``tests/topo/RecordLayoutTests.cpp`` keep it so.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_

