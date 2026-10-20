// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Intersect:Capabilities;

import Stereon.Surface;

/// The derivative capabilities of the intersection library's operations
/// (ADR-0021, section 5), as Stereon.Curve:Capabilities declares the curve
/// library's: none implements a product yet.

export namespace Stereon {

inline constexpr OperationCapability IntersectCatalogue[] = {
    Declare("Intersect", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::MovingTrimCoordinates, TrackingConvention::UniqueProjection},
            "transversal contact only: Tangential and Overlap contacts (ADR-0010) are "
            "nondifferentiable; the branch count and the tolerance stay fixed"),
};

/// The intersection library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> IntersectCapabilities() noexcept
{
    return IntersectCatalogue;
}

} // namespace Stereon
