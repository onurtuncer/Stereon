.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Exact signs, filters and enclosures
===================================

Orientation and in-circle or in-sphere decisions depend on the sign of a determinant. A fast sign is usable only outside its error bound. An interval containing zero leaves the decision open. Expansion arithmetic retains the low parts lost by rounded sums and products so the fallback can determine the sign on the supported input domain. Zero has geometric meaning (for example, collinearity); substituting a fixed epsilon changes that meaning. The chapter derives the determinants, filter stages and numerical-domain restrictions.

In ``Orient3DAdaptive``, translate three input points by the fourth.
Compute the determinant of the rounded differences as an expansion and
compare its estimate with the error bound. If all subtraction tails vanish,
the expansion's sign suffices. Otherwise add the first-order tail correction,
test the tighter bound, and invoke the full exact determinant if necessary.
The original coordinates enter the exact fallback; recomputing on the rounded
differences alone would lose the data that triggered it.

Method references: :cite:p:`shewchuk1997adaptive,dekker1971floating,priest1991algorithms,bronnimann2001interval`.

Detailed treatment: :doc:`/predicates`.

Implementation notes
--------------------

.. algorithm-sources::
   src/core/Interval.cppm
   src/robust/Domain.cppm
   src/robust/Expansion.cppm
   src/robust/Filter.cppm
   src/robust/Predicates.cppm
   src/robust/Robust.cpp
   src/robust/Robust.cppm
   src/robust/Sign.cppm
