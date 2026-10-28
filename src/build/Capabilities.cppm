// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Build:Capabilities;

import Stereon.Topo;

/// The derivative capabilities of the construction library's operations
/// (ADR-0021, section 5), as Stereon.Curve:Capabilities declares the curve
/// library's: none implements a product yet.

export namespace Stereon {

inline constexpr OperationCapability BuildCatalogue[] = {
    Declare("MakeBoxSolid", {Dependency::DefiningValues, Dependency::Topology},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter},
            "smooth in the origin, the frame and the sides; the topology is fixed (8 vertices, 12 "
            "edges, 6 faces) and the frame stays orthonormal"),
    Declare("MakeCylinderSolid", {Dependency::DefiningValues, Dependency::Topology},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter},
            "smooth in the base centre, the frame, the radius and the height; the seam stays at "
            "angle 0 and the radius and height stay positive"),
};

/// The construction library's declarations. KernelCapabilities(), in
/// Stereon.Boolean, collects them with every other library's.
[[nodiscard]] inline std::span<const OperationCapability> BuildCapabilities() noexcept
{
    return BuildCatalogue;
}

} // namespace Stereon
