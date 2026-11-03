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
/// library's. Intersect implements the forward and reverse products for a
/// line against a plane (IntersectJvp, IntersectVjp, through Linearize);
/// every other pair answers NotSupported.

export namespace Stereon {

inline constexpr OperationCapability IntersectCatalogue[] = {
    Declare("Intersect", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::MovingTrimCoordinates, TrackingConvention::UniqueProjection},
            "a line against a plane, through IntersectJvp and IntersectVjp, in the line's origin "
            "and direction and the plane's origin: decided transversal crossings only, NearEvent "
            "below NearEventSine; Tangential and Overlap contacts (ADR-0010) are "
            "nondifferentiable, undecided points indeterminate; the point count and the "
            "tolerance stay fixed; other pairs answer NotSupported",
            {DerivativeProduct::Jvp, DerivativeProduct::Vjp}),
};

/// The intersection library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> IntersectCapabilities() noexcept
{
    return IntersectCatalogue;
}

} // namespace Stereon
