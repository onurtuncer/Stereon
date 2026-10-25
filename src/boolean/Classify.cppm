// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

export module Stereon.Boolean:Classify;

import :Bodies;
import :Winding;

/// Classification (the Boolean prototype's M7, pipeline step 5; ADR-0013).
/// The faces of an imprinted body fall into regions, connected across every
/// edge but the imprinted ones; each region lies inside the other solid or
/// outside it. One seed point per region is classified by deterministic
/// rays: a ray from a point inside one of the region's faces, cast against
/// every face of the other body, crosses it an odd number of times from
/// inside. A ray that meets a face's boundary, grazes its surface or starts
/// on it is ambiguous and the next of up to sixteen directions is tried;
/// when every direction is ambiguous the seed's bounded winding number
/// decides (:Winding), and when that is undecided too the region is
/// Unresolved and the operation fails naming it, the rays tried and the
/// last interval (decision D6 of the plan), never guessing. Every region is
/// seeded on its own, and the seeds are then checked against each other
/// across the imprinted edges by the branches' side transitions (ADR-0029,
/// section 3): regions on either side of an OppositeSide edge must differ,
/// of a SameSide edge agree, and an Unresolved edge is never crossed. A
/// disagreement fails the operation: the cross-check table of ADR-0013.

export namespace Stereon {

/// Where a region lies relative to the other solid.
enum class Side : std::uint8_t {
    Inside,
    Outside,
    /// Neither a ray nor the winding number could settle it (decision D6:
    /// the operation fails).
    Unresolved,
};

/// How much work a seed may take before its region is Unresolved
/// (ADR-0013, the decisions on acceptance): the ray directions tried, then
/// the rounds the winding number's tessellation is refined.
struct ClassifyBudget {
    int Rays = 16;
    int Refinements = 8;
};

/// A connected set of a body's faces bounded by imprinted edges, with its
/// side, the face and point its seed was cast from, how many rays it took,
/// and the winding number when the rays could not settle it.
struct Region {
    std::vector<FaceId> Faces;
    // Qualified: an unqualified type would change the meaning of the name
    // the member declares (-Wchanges-meaning).
    Stereon::Side Side = Stereon::Side::Unresolved;
    FaceId SeedFace;
    Point3 Seed;
    int Rays = 0;
    /// The winding number that decided the region once every ray was
    /// ambiguous; nothing when a ray decided it.
    std::optional<WindingNumber> Winding;
};

/// The regions of both bodies.
struct Classification {
    std::vector<Region> A;
    std::vector<Region> B;
};

/// Classifies the regions of both imprinted bodies against each other with
/// the default budget. Unresolved naming the region when neither its rays
/// nor its winding number settle it; DegenerateResult when two regions
/// disagree across an imprinted edge with the branch's side transition,
/// which means the intersection, the imprint or a seed went wrong
/// (ADR-0013).
[[nodiscard]] Result<Classification> Classify(const ImprintedBodies& bodies, const Context& ctx);

/// As above, within `budget`; a budget without rays goes straight to the
/// winding number.
[[nodiscard]] Result<Classification> Classify(const ImprintedBodies& bodies, const Context& ctx,
                                              const ClassifyBudget& budget);

} // namespace Stereon
