.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. _api:

API reference
=============

Stereon is eleven C++ modules, one per library, stacked so that each uses
only the ones below it. This reference has a page per library. Each page
starts with what the library is for and where to begin, then shows how its
types fit together and how its module is divided, and ends with the full
reference, generated from the source comments.

The chapters of the *Kernel guide* explain the mathematics and the design
decisions; this part answers "what is it called, what does it take, and
what does it return". Each library page links to the chapters behind it.

.. stereon-library-graph::
   :caption: The libraries. An arrow points to a library used directly;
             what that one uses comes with it. Click a library to open its
             page.

.. toctree::
   :maxdepth: 1
   :caption: Libraries, from the bottom up

   api/core
   api/robust
   api/curve
   api/surface
   api/intersect
   api/topo
   api/build
   api/boolean
   api/blend
   api/mesh
   api/io

Reading the reference
---------------------

A few conventions run through the whole API. Knowing them makes the pages
much shorter to read.

**One namespace, one import per library.** Everything public is in
``Stereon``. Importing a library gives you it and everything below it, so
``import Stereon.Boolean;`` is enough to build and combine solids:

.. code-block:: cpp

   import Stereon.Surface;          // also brings in Core, Robust and Curve

   using namespace Stereon;

   const Circle rim{Point3{0, 0, 0}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
   const Point3 top = Eval(rim, std::numbers::pi / 2);   // (0, 5, 0)
   const double t = Project(rim, Point3{0, -9, 3});      // 3 pi / 2

**Operations are free functions, overloaded per type.** ``Eval``,
``Deriv``, ``Project``, ``BoundingBox``, ``Reversed``, ``Transformed`` and
``ToBSpline`` work the same way on every curve and surface: geometry types
are plain aggregates with no base class, and a call finds the right overload
by argument-dependent lookup (ADR-0003). Each library page opens its
reference with an *Operations at a glance* table listing every function and
the types it is declared for.

**Templates and their everyday names.** The geometry is written once for
any scalar type and dimension, as class templates named ``Basic...``;
the names you use are aliases of them. ``Line`` is ``BasicLine<double, 3>``
and ``Line2d`` is ``BasicLine<double, 2>``, the same line in a surface's
parameter plane. The type maps on the :doc:`api/curve` and
:doc:`api/surface` pages show every pair.

**Failure is a value.** An operation that can fail returns
``Result<T>``, a ``std::expected<T, KernelError>``. The ``KernelError``
says what went wrong (an ``ErrorCode``), why, in words, and which
tolerances were involved. The kernel does not throw; only
``std::bad_alloc`` may still propagate (see :doc:`errors`).

**Checked factories for outside data.** Geometry from a file or a user
goes through a ``Make...`` function, such as ``MakeCircle``,
``MakeBSplineCurve`` or ``MakeTrimmedSurface``. These validate it and
return ``InvalidInput`` with a reason when it is malformed. Aggregate
initialisation is for values the kernel itself produced (ADR-0006).

**The Context carries everything else.** Tolerances, the scheduler, the
checking level and the log and journal sinks travel in a ``Context``
argument instead of global state (ADR-0001, ADR-0007).
``Context::Default()`` is a sensible start.

**Lengths are millimetres** unless a type says otherwise, and every
numerical allowance has its own type (``ApproximationBudget``,
``DisplayTolerance``, ...), so one cannot be passed for another
(ADR-0028, :doc:`tolerances`).

Where to find things
--------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - To ...
     - Look in
   * - use points, vectors, transforms, intervals and boxes
     - :doc:`api/core`, partitions *Vector*, *Matrix*, *Interval*, *Box*
   * - decide orientation or in-circle questions exactly
     - :doc:`api/robust`, ``Orient2D``, ``Orient3D``, ``InCircle``, ``InSphere``
   * - make, evaluate or project a curve
     - :doc:`api/curve`: the type for your curve, then ``Eval``, ``Deriv``, ``Project``
   * - work with B-splines: insert knots, elevate the degree, split, fit
     - :doc:`api/curve`, *BSplineCurve*, *BSplineOperations*, *Fitting*
   * - measure a curve: curvature, Frenet frame, arc length
     - :doc:`api/curve`, *Properties*
   * - make, evaluate or project a surface; normals and curvatures
     - :doc:`api/surface`, then *Properties*
   * - hold any curve or surface in one variable
     - ``Curve``, ``Curve2d`` and ``Surface``, in :doc:`api/curve` and :doc:`api/surface`, *Erased*
   * - intersect curves and surfaces
     - :doc:`api/intersect`, ``Intersect``
   * - build or read a B-rep
     - :doc:`api/topo`: ``ShapeBuilder`` to build, ``Freeze`` to publish, ``Shape`` to read
   * - make a box or a cylinder
     - :doc:`api/build`, ``MakeBoxSolid``, ``MakeCylinderSolid``
   * - unite, cut or intersect solids
     - :doc:`api/boolean`, ``Boolean``
   * - sample geometry for display
     - :doc:`api/mesh`, ``SampleCurve``, ``SampleSurface``
   * - save and load geometry
     - :doc:`api/io`, ``SaveGeometry``, ``LoadGeometry``

What is left out
----------------

The reference covers the API a program uses. Implementation details are
left out: the ``Stereon::Detail`` and ``Stereon::Testing`` namespaces,
private members, and the module partitions a library keeps to itself (such
as ``Stereon.Curve:Homogeneous``). The source links on each page lead to
them. Stereon is pre-1.0, and the API changes from release to release; the
roadmap's versioning policy says how (``docs/planning/roadmap-and-versioning.md``).
