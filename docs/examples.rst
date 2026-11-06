.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Examples
========

Five complete programs, from evaluating a curve to drilling a solid and
intersecting three cylinders. Each is a file in ``examples/``, included
here verbatim. They are built with the tests and run by ``ctest`` on every
change (``example: curves`` and so on), and each checks its own result, so
the code on this page always compiles and does what it says. The output
shown is what they print.

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

The Steinmetz solids
--------------------

Two cylinders of one radius :math:`r` whose axes cross at right angles
have the *Steinmetz solid*, or bicylinder, in common; three have the
*tricylinder*. They are hard cases for a Boolean and good checks of one.
They are hard because the cylinders are tangent to each other at points,
where their curves of intersection cross, so that four pieces of curve
meet at one vertex inside a face. They are good checks because their
volumes and areas are known in closed form.

The program unites two such cylinders, intersects them and cuts each from
the other, then intersects their common with a third. Each result is
counted, checked closed by the validity checker at the Solid level, and
asked by its own winding number whether it holds eight points, one for
each way of being inside some of the three cylinders and outside the
others.

.. literalinclude:: ../examples/steinmetz.cpp
   :language: cpp
   :lines: 16-

.. code-block:: text

   z u x        solids 1  faces  8  edges 14  vertices  8  closed true  members true
   z n x        solids 1  faces  5  edges  7  vertices  4  closed true  members true
   z - x        solids 2  faces  6  edges 10  vertices  6  closed true  members true
   x - z        solids 2  faces  7  edges 11  vertices  6  closed true  members true
   z n x n y    solids 1  faces 12  edges 24  vertices 14  closed true  members true

The cylinders are named by their axes; ``u`` is a union, ``n`` a common
and ``-`` a cut.

The curves
~~~~~~~~~~

Take the cylinder about :math:`z`, :math:`x^2 + y^2 = r^2`, and the one
about :math:`x`, :math:`y^2 + z^2 = r^2`. Subtracting gives
:math:`x^2 = z^2`: the cylinders meet in the planes :math:`x = z` and
:math:`x = -z`, in two ellipses of semi-axes :math:`r` and
:math:`r\sqrt{2}`. The ellipses cross each other at the *nodes*
:math:`(0, \pm r, 0)`, where both cylinders have the normal :math:`\pm y`
and are tangent. The cylinder–cylinder handler returns the two ellipses as
branches that end at the nodes and at the seams of both cylinders
(:doc:`intersection`, *Cylinders with skew axes*): six pieces, which the
imprint writes into both lateral faces, four of them meeting at each node
(:doc:`imprinting`, *Two cylinders of one radius whose axes cross*).

Two cylinders
~~~~~~~~~~~~~

The two ellipses cut each lateral face into four regions about the nodes.
On the cylinder about :math:`z` they are the parts above and below both
ellipses and two *lenses* between them, about :math:`+x` and :math:`-x`;
the lens about :math:`+x` holds that cylinder's seam and is two faces. On
the cylinder about :math:`x` they are the parts towards either end and the
lenses above and below; its seam passes through a node and lies in the end
parts. Each result keeps some of these faces and some of the four caps:

.. list-table::
   :header-rows: 1
   :widths: 14 46 8 8 8 8

   * - Result
     - Its faces
     - Solids
     - Faces
     - Edges
     - Vertices
   * - :math:`z \cup x`
     - The four caps, and of each lateral face the two parts outside the
       other cylinder
     - 1
     - 8
     - 14
     - 8
   * - :math:`z \cap x`
     - The four lenses, the one about :math:`+x` in its two halves: the
       Steinmetz solid
     - 1
     - 5
     - 7
     - 4
   * - :math:`z - x`
     - Above and below: a cap, the part of the lateral face beyond the
       ellipses and a lens of the other cylinder
     - 2
     - 6
     - 10
     - 6
   * - :math:`x - z`
     - At either end: a cap, the end part of the lateral face and a lens of
       the other cylinder, the one about :math:`+x` in its two halves
     - 2
     - 7
     - 11
     - 6

The six pieces of the ellipses are ``Intersection`` edges in every result:
each joins a face of one cylinder to a face of the other. A cut leaves two
lobes, two solids, which touch at the two nodes and share a vertex at
each: the cylinders are tangent there, so what is left of one pinches to a
point. Every face has one loop, so the counts check against the
Euler–Poincaré formula (:doc:`topology`) as :math:`V - E + F`. The common
has :math:`4 - 7 + 5 = 2`, one closed surface of genus zero, and so has
the union, :math:`8 - 14 + 8`. The cuts have :math:`6 - 10 + 6 = 2` and
:math:`6 - 11 + 7 = 2`, not the 4 of two separate spheres, by the two
vertices their lobes share.

Three cylinders
~~~~~~~~~~~~~~~

A third cylinder about :math:`y`, :math:`x^2 + z^2 = r^2`, is tangent to
the one about :math:`z` at :math:`(\pm r, 0, 0)` and to the one about
:math:`x` at :math:`(0, 0, \pm r)`, and meets each in two more ellipses.
The six ellipses meet three at a time at the eight points
:math:`(\pm r, \pm r, \pm r)/\sqrt{2}`. The common of the three is the
Steinmetz solid of two cut down by the third. It has the faces, edges and
vertices of a rhombic dodecahedron: 12 faces, four on each cylinder; 24
edges, every one an ``Intersection`` edge; and 14 vertices, the six nodes,
where four faces meet, and the eight points, where three do. Again
:math:`14 - 24 + 12 = 2`.

The closed forms
~~~~~~~~~~~~~~~~

A plane :math:`y = y_0` cuts the cylinder about :math:`z` in the strip
:math:`\lvert x \rvert \le \sqrt{r^2 - y_0^2}` and the cylinder about
:math:`x` in the strip :math:`\lvert z \rvert \le \sqrt{r^2 - y_0^2}`, so
it cuts their common in a square of side :math:`2\sqrt{r^2 - y_0^2}`:

.. math::

   V_2 = \int_{-r}^{r} 4\,(r^2 - y^2)\,dy = \frac{16}{3}\,r^3 .

On the cylinder about :math:`z`, at the angle :math:`u` about its axis,
the common keeps the heights
:math:`\lvert z \rvert \le r\lvert\cos u\rvert`, an area of
:math:`\int_0^{2\pi} 2 r \lvert\cos u\rvert\, r\,du = 8 r^2`. The other
cylinder gives as much:

.. math::

   S_2 = 16\,r^2 .

For the tricylinder,

.. math::

   V_3 = 8\,(2 - \sqrt{2})\,r^3, \qquad S_3 = 24\,(2 - \sqrt{2})\,r^2 ,

a third of that area on each cylinder. Every face of both solids lies at
the distance :math:`r` from the origin along its normal, so the divergence
theorem gives :math:`V = r S / 3` for each, which these forms satisfy.

The other results follow from these. With :math:`V_z, S_z` the volume and
area of the cylinder about :math:`z`, and likewise for the others:

- the union of two has the volume :math:`V_z + V_x - V_2` and the area
  :math:`S_z + S_x - S_2`;
- a cut has the volume :math:`V_z - V_2` and the area :math:`S_z` exactly:
  it gives up the :math:`8 r^2` of its own lateral face that lay inside
  the other cylinder and gains the :math:`8 r^2` of the other's that lay
  inside it;
- the union of three has, by inclusion and exclusion, the volume
  :math:`V_z + V_x + V_y - 3 V_2 + V_3` and the area
  :math:`S_z + S_x + S_y - 3 S_2 + S_3`.

For :math:`r = 3`, the cylinder about :math:`z` 10 long and the other two
12 long, as in the program:

.. list-table::
   :header-rows: 1
   :widths: 20 30 30

   * - Result
     - Volume
     - Area
   * - :math:`z \cup x`
     - :math:`198\pi - 144 \approx 478.035`
     - :math:`168\pi - 144 \approx 383.788`
   * - :math:`z \cap x`
     - :math:`144`
     - :math:`144`
   * - :math:`z - x`
     - :math:`90\pi - 144 \approx 138.743`
     - :math:`78\pi \approx 245.044`
   * - :math:`x - z`
     - :math:`108\pi - 144 \approx 195.292`
     - :math:`90\pi \approx 282.743`
   * - :math:`z \cap x \cap y`
     - :math:`216\,(2 - \sqrt{2}) \approx 126.530`
     - :math:`216\,(2 - \sqrt{2}) \approx 126.530`
   * - :math:`z \cup x \cup y`
     - :math:`306\pi - 432 + 216\,(2 - \sqrt{2}) \approx 655.857`
     - :math:`258\pi - 432 + 216\,(2 - \sqrt{2}) \approx 505.061`

What holds the cases
~~~~~~~~~~~~~~~~~~~~

``tests/boolean/MeasureTests.cpp`` measures the volume and the area of
every result in the table, by integrating over its faces, and compares
them with these values. The ellipses are traced and fitted curves, known
to the result's tolerances, so the test allows the result's area times its
largest tolerance on the volume, and its edges' length times that
tolerance on the area. What is measured agrees with the closed forms
within :math:`2 \times 10^{-6}`: eight digits.

Beside the program and the measures, ``tests/boolean/BodiesTests.cpp``
checks the imprint of two cylinders, the four pieces at each node and the
six shared pairs. ``tests/boolean/AssembleTests.cpp`` checks the counts,
the closure, the ``Intersection`` edges and the membership of the union,
the common, a cut and the tricylinder.
``tests/boolean/IntegrationTests.cpp`` holds the pair of cylinders, in all
four operations, to what any correct Boolean does: a closed solid, no
duplicate edge or vertex, the volume identities
:math:`\operatorname{vol}(A \cup B) + \operatorname{vol}(A \cap B) =
\operatorname{vol}(A) + \operatorname{vol}(B)` and
:math:`\operatorname{vol}(A - B) = \operatorname{vol}(A) -
\operatorname{vol}(A \cap B)`, and membership over a grid of points.

What is not covered: cylinders of one radius whose axes cross at another
angle meet in two ellipses in the same way and go through the same code,
but no test holds them yet; and cylinders on one axis and radius, whose
lateral faces coincide, are refused with ``NotSupported``
(:doc:`modeling`).

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
