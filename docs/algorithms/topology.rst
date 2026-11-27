.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Boundary representation, Euler operators and validation
=======================================================

Euler edits update incidence among vertices, edges, coedges, loops and faces. Each operator has topological preconditions; preserving an Euler count alone cannot establish that a loop closes or that its coedges agree geometrically. The builder checks connectivity before publishing an immutable shape. Compaction assigns dense indices to live slots, rewrites links through the resulting maps, then validates the new snapshot. The remap is essential for clients following entity lineage. Primitive construction supplies explicit geometry and oriented boundary uses; see :doc:`/modeling`.

Method references: :cite:p:`mantyla1988solid,hoffmann1989geometric`.

Detailed treatment: :doc:`/topology`.

Implementation notes
--------------------

.. algorithm-sources::
   src/blend/Blend.cppm
   src/build/Build.cppm
   src/build/Capabilities.cppm
   src/build/Primitives.cpp
   src/build/Primitives.cppm
   src/topo/Builder.cpp
   src/topo/Builder.cppm
   src/topo/Capabilities.cppm
   src/topo/Compact.cpp
   src/topo/Describe.cpp
   src/topo/Euler.cpp
   src/topo/EulerRings.cpp
   src/topo/EulerSplit.cpp
   src/topo/Operation.cpp
   src/topo/Operation.cppm
   src/topo/Records.cppm
   src/topo/Report.cppm
   src/topo/Shape.cpp
   src/topo/Shape.cppm
   src/topo/Topo.cppm
   src/topo/Validate.cpp
   src/topo/ValidateGeometry.cpp
   src/topo/ValidateSolid.cpp
