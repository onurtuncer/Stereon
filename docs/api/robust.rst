.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Robust — ``Stereon.Robust``
===========================

Exact answers to the yes-or-no questions geometry asks: which side of a
plane a point is on, whether a point lies inside a circle. Every topological
decision in :doc:`topo`, :doc:`boolean` and :doc:`mesh` goes through a
predicate from here (ADR-0002). In floating point, a decision made from a
rounded determinant can contradict another made a moment earlier, and a
B-rep built on contradictions falls apart. These predicates are never
wrong, and on most inputs they cost little more than the rounded
computation, because they only compute exactly when a fast filter cannot
decide.

.. code-block:: cpp

   import Stereon.Robust;

   using namespace Stereon;

   const Sign side = Orient3D(a, b, c, d);
   if (side == Sign::Zero) {
       // d lies exactly in the plane of a, b and c: no tolerance involved
   }

Start here
----------

``Orient2D``, ``Orient3D``, ``InCircle``, ``InSphere``
   The four classical predicates, after Shewchuk (1997). Each returns the
   exact ``Sign`` of a determinant, never the determinant itself.

``CompareDistance``
   Whether two points are closer than a radius, exactly.

``Sign``
   ``Negative``, ``Zero`` or ``Positive``: the only answer a predicate
   gives.

``SignOf``, ``FilteredSign``
   For predicates on curves and surfaces: evaluate in interval arithmetic,
   and fall back to an exact evaluation only when the interval contains
   zero.

``Expansion``, ``DynamicExpansion``
   The exact arithmetic underneath: a real number as a sum of
   non-overlapping doubles. For writing a new exact predicate.

The guide chapter behind this library: :doc:`/predicates`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: robust

Inside the module
-----------------

.. stereon-partition-graph:: robust

Reference
---------

.. stereon-api:: robust
