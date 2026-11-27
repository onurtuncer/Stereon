.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Boolean clipping, imprinting, classification and assembly
=========================================================

Pair candidate faces using their conservative boxes, intersect their supports, and clip the resulting branches to each face's loops. Split crossed edges and imprint the resulting networks on both operands. Ray classification needs an unambiguous crossing: a hit at a loop vertex can require another ray or a refusal. The Boolean operation selects oriented fragments by their classifications. Assembly pairs compatible boundary uses and tracks the original entities through splits. The worked examples in the imprinting chapter explain why coincident faces and seams require special handling.

In ``NetworkOf``, a clipped section piece contributes vertices and an edge
for the appropriate face. A whole branch without a trim crossing is closed
only if its two evaluated endpoints agree to the applicable tolerance.
A segment wholly inside a face can also have no crossings; marking every
such segment closed would join distinct endpoints during imprinting.

Method references: :cite:p:`mantyla1988solid,hoffmann1989geometric`.

Detailed treatment: :doc:`/imprinting`.

Implementation notes
--------------------

.. algorithm-sources::
   src/boolean/Assemble.cpp
   src/boolean/Assemble.cppm
   src/boolean/Bodies.cpp
   src/boolean/Bodies.cppm
   src/boolean/Boolean.cppm
   src/boolean/Capabilities.cppm
   src/boolean/Classify.cpp
   src/boolean/Classify.cppm
   src/boolean/Clip.cpp
   src/boolean/Clip.cppm
   src/boolean/Imprint.cpp
   src/boolean/Imprint.cppm
