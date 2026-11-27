.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Spline evaluation, subdivision and control-net editing
======================================================

Locate the knot span and evaluate the locally nonzero basis functions. For rational geometry accumulate weighted coordinates and weight before division. A tensor-product surface repeats the basis construction in the two parameter directions. Repeated knots affect which one-sided derivatives exist. Knot insertion changes the control net while preserving the represented curve; splitting at inserted knots gives Bézier pieces suitable for enclosure. Degree elevation, knot removal and periodic rearrangement have different acceptance conditions. The chapter identifies the operations that can change the geometry and the bounds they check.

Method references: :cite:p:`piegl1997nurbs,deboor2001practical`.

Detailed treatment: :doc:`/nurbs`.

Implementation notes
--------------------

.. algorithm-sources::
   src/curve/BSplineBasis.cpp
   src/curve/BSplineBasis.cppm
   src/curve/BSplineCurve.cpp
   src/curve/BSplineCurve.cppm
   src/curve/BSplineOperations.cpp
   src/curve/BSplineOperations.cppm
   src/curve/BatchEvaluation.cpp
   src/curve/BatchEvaluation.cppm
   src/curve/FastBatchEvaluation.cpp
   src/curve/Homogeneous.cppm
   src/surface/BSplineSurface.cpp
   src/surface/BSplineSurface.cppm
   src/surface/GridEvaluation.cpp
   src/surface/GridEvaluation.cppm
