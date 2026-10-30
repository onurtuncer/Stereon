// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The broad phase: which boxes of one list meet which of another, without
// testing every pair. Two lists of ten thousand boxes are a hundred million
// pairs; a box hierarchy over one list finds the k pairs that meet in about
// (m + n) log n + k tests.

module;

#include <compare>
#include <cstddef>
#include <span>
#include <vector>

export module Stereon.Core:BoxPairs;

import :Box;

export namespace Stereon {

/// The places of two boxes that meet, A's in its list and B's in its own.
struct BoxPair {
    std::size_t A;
    std::size_t B;

    constexpr auto operator<=>(const BoxPair&) const = default;
};

/// Every pair (i, j) for which Intersects(a[i], b[j]), ordered by i and then
/// by j: exactly what testing every pair in that order would list, so a
/// caller that depends on the order of its pairs (ADR-0007) can take this in
/// its place. Boxes meet as Intersects says: closed, so that boxes touching
/// at a face, edge or corner meet, and an empty box meets nothing.
///
/// A hierarchy of B's boxes, split at the median of their centres along the
/// longest axis, is queried with each of A's; boxes with an unbounded side,
/// which have no centre, are tested against every query. The hierarchy is
/// built in O(n log n) and each query visits the nodes whose boxes meet it,
/// so the work is about (m + n) log n + k for k pairs, whatever the boxes'
/// layout, where testing every pair is m n.
[[nodiscard]] std::vector<BoxPair> OverlappingPairs(std::span<const Box3> a,
                                                    std::span<const Box3> b);

} // namespace Stereon
