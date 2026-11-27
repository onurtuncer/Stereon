.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Space-curve intersection and minimum distance
=============================================

Two lines admit a closed-form configuration test. For general supported pairs, bounding Bézier pieces separates distant portions before testing possible close meetings. Closest approach minimises half the squared separation over the two parameter domains. Bounds discard boxes that cannot beat the incumbent, and an active-set refinement handles a minimum on an endpoint. A local stationary point alone does not prove a global minimum; the solver also has to exclude competitors or report ambiguity.

The routine ``Searched`` clusters surviving parameter boxes before choosing
the result status. A single cluster still needs a proof of uniqueness.
With several clusters it refines each and reports the best candidate with
undecided uniqueness. Exhausting the search budget takes a separate path;
the numerical candidate survives, but the exhausted search cannot certify
that it excluded all competing minima.

Method references: :cite:p:`patrikalakis2002shape,neumaier1990interval`.

Detailed treatment: :doc:`/intersection`.

Implementation notes
--------------------

.. algorithm-sources::
   src/intersect/ApproachSide.hpp
   src/intersect/ClosestApproach.cpp
   src/intersect/ClosestApproach.cppm
   src/intersect/LineLine.cpp
   src/intersect/SplineCurvePair.cpp
