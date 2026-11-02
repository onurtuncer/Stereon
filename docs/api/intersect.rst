.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Intersection — ``Stereon.Intersect``
====================================

Where curves and surfaces meet. One function, ``Intersect``, is
overloaded for every pairing:

- curve and curve;
- two curves in a parameter plane;
- curve and surface;
- surface and surface.

Each overload comes with or without explicit parameter domains and
tolerance. The result types say everything a modeller needs to trust the
answer: isolated points, branches with their pcurves on both surfaces, and
overlapping regions, each with a bound on its error (ADR-0010, ADR-0029).

A successful result is complete. Within the requested domains, nothing is
missing, so a result with no component proves the inputs are apart. A pair
that is not handled yet says so rather than guessing:

- a built-in pair without a handler returns ``NotImplemented``, naming
  both types;
- a type that is not built in gets ``NotSupported``.

.. code-block:: cpp

   import Stereon.Intersect;

   using namespace Stereon;

   const Context ctx = Context::Default();
   const Result<SurfaceSurfaceIntersection> meet = Intersect(Surface(plane), Surface(cylinder), ctx);
   if (meet) {
       for (const IntersectionBranch& branch : meet->Branches) { /* a curve on both */ }
   }

Start here
----------

``Intersect``
   The entry point, overloaded by what you pass. Without domains it
   intersects over the whole of each input, at the Context's resolution.

``SurfaceSurfaceIntersection``, ``IntersectionBranch``, ``IntersectionPoint``, ``OverlapRegion``
   What two surfaces have in common.

``CurveSurfaceIntersection``, ``CurveCurveIntersection``, ``PlanarIntersection``
   The same for the other pairings: points of contact and coincident
   stretches.

``ContactType``, ``SideTransition``
   How the inputs meet along a branch (``Crossing``, ``Tangential`` or
   ``Overlap``), and whether the inside and outside swap across it. The
   Booleans read these.

The guide chapter behind this library: :doc:`/intersection`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: intersect

Inside the module
-----------------

.. stereon-partition-graph:: intersect

Reference
---------

.. stereon-api:: intersect
