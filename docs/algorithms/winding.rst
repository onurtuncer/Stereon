.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Bounded winding-number classification
=====================================

For a triangle with vertex offsets a, b, c from the query point, SolidAngle evaluates twice the oriented atan2 of the scalar triple product and the symmetric length/dot-product denominator. Summing signed solid angles and dividing by 4π yields the winding value. Reversing orientation changes the sign.

A curved face introduces an approximation error beyond rounding the triangle formula. Stereon samples boundary curves, clips loops while subdividing parameter cells, and tracks the deviation of the resulting fans. Cuts retain matching corners on both sides so internal edges cancel in space. A classified point must be separated sufficiently for the accumulated bound to decide the result. The triangle formula, generalised winding theory and polygon clipping each have a distinct reference below. The curved-face deviation budget is Stereon's construction.

For offset vectors of lengths A, B and C the oriented contribution is

.. math::

   \Omega=2\operatorname{atan2}\!\left(
     a\cdot(b\times c),\;
     ABC+C(a\cdot b)+A(b\cdot c)+B(c\cdot a)\right).

Using both arguments of atan2 preserves the quadrant when the denominator
is negative. In ``WindingNumberOf``, refinement also accounts for cells
whose separation from the query has not been proved. ``CheckWhole``
checks a decided value against the admissible whole value and its error
bound.

Method references: :cite:p:`vanoosterom1983solid,jacobson2013winding,sutherland1974clip`.

Detailed treatment: :doc:`/modeling`.

Implementation notes
--------------------

.. algorithm-sources::
   src/boolean/Winding.cpp
   src/boolean/Winding.cppm
