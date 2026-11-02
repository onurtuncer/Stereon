.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

.. cpp:namespace:: Stereon

Booleans — ``Stereon.Boolean``
==============================

Union, cut and common of solids (ADR-0013). ``Boolean`` is the whole
operation; the rest of the library is its pipeline, public so that each
step can be tested, inspected and reused on its own:

#. find the pairs of faces whose boxes meet;
#. intersect their surfaces;
#. clip each intersection branch to the faces (*Clip*);
#. imprint the pieces into both bodies (*Imprint*, *Bodies*);
#. classify the regions of each body as inside or outside the other
   (*Classify*, with the winding number of *Winding* as the fallback);
#. assemble the regions the operation keeps (*Assemble*).

The result is a certified ``Shape`` with the history of what became of
each input entity.

.. code-block:: cpp

   import Stereon.Boolean;

   using namespace Stereon;

   const Context ctx = Context::Default();
   const auto box = MakeBoxSolid(Box3{{0, 2}, {0, 3}, {0, 5}}, ctx);
   const auto pin = MakeCylinderSolid(Point3{1, 1, 0}, {1, 0, 0}, {0, 1, 0}, 0.5, 6, ctx);

   const Result<OperationResult> drilled = Boolean(box->Shape, pin->Shape, BooleanOperation::Cut, ctx);
   if (!drilled) {
       // e.g. Unresolved: a region could not be classified within its budget
   }

Start here
----------

``Boolean``, ``BooleanOperation``
   The operation: ``Union``, ``Cut`` or ``Common`` of two solids. It
   returns an ``OperationResult``: the shape, its ``History``, the tolerance
   changes it made and its warnings.

``ImprintBodies``
   Steps 1 to 4 for two solids: both bodies with the intersection imprinted,
   and the edges each piece became in each.

``Classify``, ``ClassifyBudget``
   Step 5: the regions of both bodies and their sides. A region that
   cannot be decided within the budget is ``Unresolved`` rather than
   guessed.

``ClipToFace``, ``ImprintClip``, ``WindingNumberOf``
   The single steps, for one face and one branch.

``KernelCapabilities``
   What every operation of every library can differentiate. It lives here
   because this library sits above all the others.

The guide chapters behind this library: :doc:`/modeling` and
:doc:`/imprinting`.

Where it sits
-------------

.. stereon-library-graph::
   :highlight: boolean

Inside the module
-----------------

.. stereon-partition-graph:: boolean

Reference
---------

.. stereon-api:: boolean
