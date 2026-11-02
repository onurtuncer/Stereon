.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Meshing — ``Stereon.Mesh``
==========================

Tessellation, for display and, later, for analysis (ADR-0014). Today it is
display sampling. It turns curves into polylines and surfaces into
iso-lines, for the debug viewer and any other front end. A
``DisplayTolerance`` steers it: a chord height and an angle. That is a
purpose of its own, not a modelling tolerance (ADR-0028).

.. code-block:: cpp

   import Stereon.Mesh;

   using namespace Stereon;

   const Result<Polyline> outline = SampleCurve(Curve(rim), DisplayTolerance{.Chord = 0.01, .Angle = 0.1});
   if (outline) {
       // outline->Points and outline->Parameters, one per sample
   }

Start here
----------

``SampleCurve``
   A curve, or a range of it, as a ``Polyline``: points with the parameter
   of each, within the tolerance.

``SampleSurface``
   A surface as ``IsoLine``\ s: evenly spaced lines in both parameter
   directions, plus one at every interior knot of a B-spline surface.

``Breakpoints``
   The parameters a sampled curve must pass through whatever the
   tolerance: the ends of the range and every distinct knot of a B-spline.

The guide chapters behind this library: :doc:`/meshing` and
:doc:`/debugging` (the debug viewer).

Where it sits
-------------

.. stereon-library-graph::
   :highlight: mesh

Reference
---------

.. stereon-api:: mesh
