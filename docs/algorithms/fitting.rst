.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Interpolation, approximation and periodic fitting
=================================================

Assign sample parameters, choose knots, and assemble only nonzero basis coefficients into the collocation rows. Interpolation solves for poles reproducing the samples. Least squares solves a different objective; residuals at the samples do not bound unsampled geometry. Prescribed end derivatives add constraints. Periodic fits couple the first and last unknowns. Surface fitting processes rows and columns with shared parameters in each direction. Degenerate parameter sequences and singular solves need explicit outcomes.

Method references: :cite:p:`piegl1997nurbs,lee1989choosing,schoenberg1953polya`.

Detailed treatment: :doc:`/fitting`.

Implementation notes
--------------------

.. algorithm-sources::
   src/curve/Fitting.cpp
   src/curve/Fitting.cppm
   src/curve/PeriodicFitting.cpp
   src/surface/Fitting.cpp
   src/surface/Fitting.cppm
