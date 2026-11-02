.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Construction — ``Stereon.Build``
================================

Solids made from a description rather than edited face by face. The library
is meant to hold primitives, extrusion, revolution, sweeping, lofting and
sewing; today it has the two primitives the Boolean prototype needs, a box
and a cylinder. Each is built through the controlled operations of
``ShapeBuilder`` and published with ``Commit``, so the result is a
certified ``Shape`` with its history, like any other operation's.

.. code-block:: cpp

   import Stereon.Build;

   using namespace Stereon;

   const Context ctx = Context::Default();
   const Result<OperationResult> box = MakeBoxSolid(Box3{{0, 2}, {0, 3}, {0, 5}}, ctx);
   const Result<OperationResult> cylinder =
       MakeCylinderSolid(Point3{1, 1, 0}, {1, 0, 0}, {0, 1, 0}, /*radius*/ 0.5, /*height*/ 6, ctx);

Start here
----------

``MakeBoxSolid``
   A box from an axis-aligned ``Box3``, or from an origin, two axes and
   three side lengths: 8 vertices, 12 line edges, 6 planar faces.

``MakeCylinderSolid``
   A cylinder from its base centre, two axes, a radius and a height: a
   lateral face that meets itself along a seam at angle 0, and two disks.

The guide chapter behind this library: :doc:`/modeling`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: build

Inside the module
-----------------

.. stereon-partition-graph:: build

Reference
---------

.. stereon-api:: build
