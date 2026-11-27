.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Deviation certification
=======================

An intersection branch carries a space curve and parameter curves on its two surfaces. Certification bounds the deviation between each surface composed with its parameter curve and the space curve over a parameter span. Knots delimit intervals of smoothness. Derivative enclosures and Taylor remainders cover the interior between samples; subdivision tightens a bound. A certificate justifies the reported curve's accuracy. Establishing that every branch was found requires a separate argument in the primal intersection solver.

In ``CertifyDeviation``, require a finite positive target, a bounded
parameter interval and spline domains covering the requested range.
Try a whole-span closed-form bound for the applicable plane and cylinder
cases. If that bound exceeds the target, interval bisection may tighten it
or produce a witness against the claim. A missing interval evaluator gives
``NotSupported``. This distinction matters when interpreting a failed
certification run.

Method references: :cite:p:`moore2009interval,snyder1992interval,piegl1997nurbs`.

Detailed treatment: :doc:`/intersection`.

Implementation notes
--------------------

.. algorithm-sources::
   src/surface/Certify.cpp
   src/surface/Certify.cppm
