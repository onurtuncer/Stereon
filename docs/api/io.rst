.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Data exchange — ``Stereon.Io``
==============================

Reading and writing models (ADR-0012). The library is meant to cover STEP,
STL, OBJ, glTF and Stereon's native ``.stn`` format. Today it has the
geometry subset of ``.stn``:

- curves, plane curves and surfaces of the built-in types;
- the model's numerical policy (unit, resolutions, tolerance ceiling),
  stored with them, because the geometry means nothing without it;
- everything read back bit for bit.

It also writes dumps for the debug viewer: one ``.stn`` file per call, with
a debug chunk naming where and when it was made.

.. code-block:: cpp

   import Stereon.Io;

   using namespace Stereon;

   GeometryArchive archive;
   archive.Policy = PolicyOf(context);           // the model's numerical policy
   archive.Curves.push_back(Curve(circle));
   archive.Surfaces.push_back(Surface(nurbs));

   if (const Result<void> saved = SaveGeometry("part.stn", archive); !saved) {
       // InvalidInput (cannot write, bad policy) or NotSupported (a user type);
       // an existing part.stn is left as it was
   }
   const Result<LoadedGeometry> loaded = LoadGeometry("part.stn");

Start here
----------

``GeometryArchive``, ``ModelPolicy``, ``PolicyOf``
   What a file holds: the geometry, and the policy it was built under.

``SaveGeometry``, ``LoadGeometry``
   To and from a file. ``WriteGeometry`` and ``ReadGeometry`` do the same
   with bytes in memory. A damaged chunk or an unknown required one is
   refused; an unknown optional chunk is skipped, and
   ``LoadedGeometry::Warnings`` says so.

``GeometryToJson``
   A readable dump of an archive, for inspection.

``Dump``, ``DumpAt``
   Write geometry and shapes for the debug viewer, under
   ``STEREON_DUMP_DIR``.

The guide chapters behind this library: :doc:`/data_exchange` and
:doc:`/debugging`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: io

Inside the module
-----------------

.. stereon-partition-graph:: io

Reference
---------

.. stereon-api:: io
