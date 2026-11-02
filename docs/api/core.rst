.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Core — ``Stereon.Core``
=======================

The foundation every other library stands on. Core has no geometry of its
own; it has the vocabulary the geometry is written in:

- **math types**: points, vectors, matrices and affine transforms;
- **enclosures**: closed intervals with outward rounding, and the boxes made
  of them;
- **failure**: the ``Result<T>`` every fallible operation returns, and the
  ``KernelError`` inside it;
- **the Context**: tolerances, the scheduler, the checking level and the
  log and journal sinks, passed to operations instead of living in globals;
- **storage for topology**: typed handles and the copy-on-write paged arena
  behind the B-rep tables of :doc:`topo`;
- **groundwork**: numerical allowance types, enum names and serialisation by
  reflection, and the declarations of what each operation can differentiate.

Parallel execution is a header rather than a module,
``<Stereon/Core/Exec.hpp>``, because the compilers in use crash on the
library behind it inside modules.

.. code-block:: cpp

   import Stereon.Core;

   using namespace Stereon;

   constexpr Point3 corner{10, 0, 0};
   constexpr Vec3 up{0, 0, 1};

   // (a * b) applies b first: rotate about Z, then lift by 5 mm.
   const Transform lift = MakeTranslation(5.0 * up) * MakeRotation(up, 0.5);
   const Point3 moved = Transformed(corner, lift);
   const Vec3 direction = Transformed(up, lift);   // vectors ignore the translation

   const Interval sum = Interval{0.1} + Interval{0.2};
   assert(Contains(sum, 0.3));                      // the enclosure does not lose it

Start here
----------

``Point3``, ``Vec3``, ``Point2``, ``Vec2``
   Positions and displacements are different types: a point minus a point
   is a vector, and a transform moves points but only turns vectors.
   ``Dot``, ``Cross``, ``Length``, ``Normalized`` and ``Distance`` are in
   *Vector*.

``Transform``
   An affine map. ``MakeTranslation``, ``MakeRotation``, ``MakeScaling`` and
   the mirrors build one, ``AboutPoint`` moves its fixed point, and
   ``Transformed(x, t)`` applies it to anything that has an overload.

``Interval``, ``Box2``, ``Box3``
   Closed intervals whose arithmetic rounds outward, so a computed interval
   always contains the exact result; boxes are one interval per axis.
   ``Interval{}`` is empty, the identity of ``Hull``.

``Result<T>``, ``KernelError``, ``ErrorCode``, ``Fail``
   How every operation reports failure (:doc:`/errors`).

``Context``
   What an operation may consult besides its inputs. ``Context::Default()``
   is the place to start.

``Id<Kind>`` and ``VertexId`` ... ``SolidId``
   The handles a frozen shape gives out (:doc:`topo`).

``ApproximationBudget`` and the other allowances
   Each purpose gets its own type, so that a fitting budget cannot be
   passed where a modelling tolerance belongs (:doc:`/tolerances`).

The guide chapters behind this library: :doc:`/geometry` (math types,
intervals, boxes), :doc:`/errors`, :doc:`/tolerances` and
:doc:`/sensitivity`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: core

Inside the module
-----------------

.. stereon-partition-graph:: core

Reference
---------

.. stereon-api:: core
