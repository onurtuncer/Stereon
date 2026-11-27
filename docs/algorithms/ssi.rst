.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Surface–surface intersections and branch assembly
=================================================

Plane and elementary-pair handlers exploit closed forms or profile intersections. Skew cylinders are parameterised by angle on one cylinder; the quadratic in axial height has two sheets, which merge where the discriminant vanishes. A parametric–implicit pair traces the zero set of F(S(u,v)). Search the requested charts, including their seams, before joining arcs.

For two spline surfaces, bound pairs of Bézier patches and prove separation or a transverse direction. If a branch advances monotonically along that direction it must reach the patch boundary. Isolate the boundary roots, connect compatible endpoints, march in three dimensions and certify the resulting space curve and pcurves. Nodes require their own local proof. An undecided patch pair can prevent publication of the whole result, because publishing only the discovered branches would lose completeness. The procedures below include the edge, node, seam and chart-specific arguments used in the implementation.

Method references: :cite:p:`hoffmann1989geometric,patrikalakis2002shape,krawczyk1969newton,snyder1992interval`.

Detailed treatment: :doc:`/intersection`.

Implementation notes
--------------------

.. algorithm-sources::
   src/intersect/Analytic.hpp
   src/intersect/Bernstein.hpp
   src/intersect/BezierPatch.hpp
   src/intersect/Capabilities.cppm
   src/intersect/CylinderCylinder.cpp
   src/intersect/CylinderCylinderSkew.cpp
   src/intersect/Dispatch.cpp
   src/intersect/Dispatch.cppm
   src/intersect/ElementaryPairs.cpp
   src/intersect/Exact.hpp
   src/intersect/GenericSSI.cpp
   src/intersect/Implicit.hpp
   src/intersect/ImplicitTrace.cpp
   src/intersect/Intersect.cppm
   src/intersect/PlaneCylinder.cpp
   src/intersect/PlaneEllipsoid.cpp
   src/intersect/PlanePlane.cpp
   src/intersect/SeamNode.hpp
   src/intersect/TangentBranch.hpp
