// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <vector>

export module Stereon.Boolean:Classify;

import :Bodies;

/// Classification (the Boolean prototype's M7, pipeline step 5; ADR-0013).
/// The faces of an imprinted body fall into regions, connected across every
/// edge but the imprinted ones; each region lies inside the other solid or
/// outside it. One seed point per region is classified by deterministic
/// rays: a ray from a point inside one of the region's faces, cast against
/// every face of the other body, crosses it an odd number of times from
/// inside. A ray that meets a face's boundary, grazes its surface or starts
/// on it is ambiguous and the next of up to sixteen directions is tried;
/// when every direction is ambiguous the region is Unresolved, and the
/// operation fails naming it (decision D6 of the plan), never guessing.
/// Every region is seeded on its own, and the seeds are then checked
/// against each other across the imprinted edges by the branches' side
/// transitions (ADR-0029, section 3): regions on either side of an
/// OppositeSide edge must differ, of a SameSide edge agree, and an
/// Unresolved edge is never crossed. A disagreement fails the operation:
/// the cross-check table of ADR-0013.

export namespace Stereon {

/// Where a region lies relative to the other solid.
enum class Side : std::uint8_t {
    Inside,
    Outside,
    /// No ray could settle it (decision D6: the operation fails).
    Unresolved,
};

/// A connected set of a body's faces bounded by imprinted edges, with its
/// side, the face and point its seed was cast from, and how many rays it
/// took.
struct Region {
    std::vector<FaceId> Faces;
    Side Side = Side::Unresolved;
    FaceId SeedFace;
    Point3 Seed;
    int Rays = 0;
};

/// The regions of both bodies.
struct Classification {
    std::vector<Region> A;
    std::vector<Region> B;
};

/// Classifies the regions of both imprinted bodies against each other.
/// IntersectionFailed naming the region when one is Unresolved, or when
/// two regions disagree across an imprinted edge with the branch's side
/// transition.
[[nodiscard]] Result<Classification> Classify(const ImprintedBodies& bodies, const Context& ctx);

} // namespace Stereon
