.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Examples
========

Four complete programs, from evaluating a curve to drilling a solid. Each
is a file in ``examples/``, included here verbatim. They are built with the
tests and run by ``ctest`` on every change (``example: curves`` and so
on), and each checks its own result, so the code on this page always
compiles and does what it says. The output shown is what they print.

Build and run them from any build with tests (:doc:`building`):

.. code-block:: bash

   cmake --build --preset gcc16-debug
   build/gcc16-debug/examples/stereon_example_curves

A B-spline curve
----------------

A cubic B-spline is built through the checked factory, which validates the
degree, the knots and the poles and returns ``InvalidInput`` for anything
it cannot represent (:doc:`errors`). It is then evaluated, differentiated
and projected onto. At an interior closest point the residual
:math:`q - C(t)` is perpendicular to the tangent :math:`C'(t)`, the
first-order condition of the minimisation, and the program checks it
(:doc:`nurbs`).

.. literalinclude:: ../examples/curves.cpp
   :language: cpp
   :lines: 11-

.. code-block:: text

   C(0.5)   = (1.468750, 1.687500, 0.281250)
   C'(0.5)  = (2.812500, 1.125000, 0.937500)
   C''(0.5) = (-0.750000, -7.500000, 0.750000)
   closest parameter t = 0.683134609377, distance 1.224028657
   (q - C(t)) . C'(t)  = 1.804e-16

Curves against surfaces
-----------------------

``Intersect`` returns every point where a curve meets a surface, each with
its parameter on the curve, its parameters on the surface, its position,
how the two meet there and whether that was decided exactly; and every
stretch where the curve lies in the surface. A result with no points is a
proof that the two do not meet (:doc:`intersection`).

.. literalinclude:: ../examples/intersection.cpp
   :language: cpp
   :lines: 13-

.. code-block:: text

   line x sphere: 2 point(s), 0 overlap(s)
     t = 1.200000000  (u, v) = (2.498092, 0.000000)  at (-0.800000, 0.600000, 0.000000)  crossing
     t = 2.800000000  (u, v) = (0.643501, 0.000000)  at (0.800000, 0.600000, 0.000000)  crossing
   arch x plane: 2 point(s), 0 overlap(s)
     t = 0.146446609  (u, v) = (0.292893, 0.000000)  at (0.292893, 0.000000, -0.000000)  crossing
     t = 0.853553391  (u, v) = (1.707107, 0.000000)  at (1.707107, 0.000000, -0.000000)  crossing

The line crosses the sphere where :math:`x^2 + 0.6^2 = 1`, at
:math:`x = \mp 0.8`. The arch is the quadratic Bézier
:math:`z(s) = -1 + 8s - 8s^2`, which vanishes at
:math:`s = \tfrac12 \mp \tfrac{\sqrt 2}{4}`, that is 0.146447 and 0.853553.

A drilled box
-------------

Two primitive solids and a Boolean cut. ``MakeBoxSolid`` and
``MakeCylinderSolid`` return an ``OperationResult``: the shape, with the
history and diagnostics of the operation that made it. ``Boolean``
imprints each body on the other, classifies the pieces and assembles the
result (:doc:`imprinting`). Each step either succeeds with a shape that has
passed the validity checker or returns a ``KernelError``.

.. literalinclude:: ../examples/boolean.cpp
   :language: cpp
   :lines: 14-

.. code-block:: text

   box: 6 faces, 12 edges, 8 vertices
   box - cylinder: 7 faces, 15 edges, 10 vertices
   valid at the Geometric level: true

The counts check against the Euler–Poincaré formula for a solid of genus
:math:`g` with :math:`F` faces, :math:`L` loops and one shell,
:math:`V - E + 2F - L = 2(1 - g)` (:doc:`topology`). The drilled box has
the box's 6 faces and the hole's cylindrical face; 12 box edges, the two
circles where the hole leaves the box, and the hole's seam; the box's 8
vertices and one on each circle. Its top and bottom faces each have an
inner loop around the hole, so :math:`L = 9`, and
:math:`10 - 15 + 14 - 9 = 0 = 2(1 - 1)`: one hole, genus one.

The native file
---------------

Geometry is saved to a ``.stn`` archive with the model's numerical policy
(its units, resolution and tolerance ceiling) and read back. Every double
is written as its IEEE 754 bit pattern, so the round trip loses nothing,
and the program checks it bit for bit (:doc:`data_exchange`).

.. literalinclude:: ../examples/geometry_file.cpp
   :language: cpp
   :lines: 13-

.. code-block:: text

   format version 3, written by 0.1.0, resolution 1e-05 mm
   curve and surface evaluate to the same bits after the round trip: true
