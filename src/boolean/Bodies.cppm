// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <utility>
#include <vector>

export module Stereon.Boolean:Bodies;

import :Imprint;

/// The Boolean pipeline's first four steps for two solids (the plan,
/// section 2): candidate face pairs by bounding box in input order, the
/// surface intersection of each pair over the faces' parameter boxes, each
/// branch clipped to both faces, and both bodies imprinted with the curve
/// networks their faces receive, inside one transaction per body. Where a
/// branch leaves one body's face it is a junction in the other body's
/// face; the two bodies' imprints of one piece of the intersection are
/// paired for the assembly (M8) to identify.

export namespace Stereon {

/// One piece of the intersection as an edge in each body, with the branch's
/// side transition across it (ADR-0029, section 3).
struct SharedEdge {
    EdgeId A;
    EdgeId B;
    SideTransition Transition = SideTransition::Unresolved;
};

/// Both bodies imprinted, and the pieces of the intersection, each an edge
/// in each result.
struct ImprintedBodies {
    OperationResult A;
    OperationResult B;
    std::vector<SharedEdge> Shared;
};

/// Intersects every face pair of the two solids and imprints both.
/// NotSupported for coincident faces (an overlap), which the assembly of
/// M8 handles; otherwise what the intersection, the clip and the imprint
/// report.
[[nodiscard]] Result<ImprintedBodies> ImprintBodies(const Shape& a, const Shape& b,
                                                    const Context& ctx);

} // namespace Stereon
