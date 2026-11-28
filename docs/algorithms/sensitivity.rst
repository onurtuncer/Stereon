.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Derivatives of evaluation, intersections and closest approach
=============================================================

At a regular point, differentiate the defining equations F(z,p)=0. A forward product solves F_z dz = −F_p dp. A reverse product uses the transposed system to pair output sensitivities with input perturbations. A point sampled on a branch needs a tracking equation to fix motion along the branch. Closest approach instead differentiates stationarity in the free variables and retains active endpoint constraints. The primal regularity status determines whether a derivative can be reported. The cited AD survey gives the forward/reverse distinction; the Stereon chapter specifies each residual and the event policy.

Method references: :cite:p:`baydin2018ad,patrikalakis2002shape`.

Detailed treatment: :doc:`/sensitivity`.

Implementation notes
--------------------

.. algorithm-sources::
   src/core/Sensitivity.cpp
   src/core/Sensitivity.cppm
   src/curve/PoleSensitivity.cpp
   src/curve/PoleSensitivity.cppm
   src/intersect/ClosestApproachLinearization.cpp
   src/intersect/Linearization.cpp
   src/intersect/Linearization.cppm
   src/intersect/PlanarLinearization.cpp
   src/intersect/SplineJoint.hpp
