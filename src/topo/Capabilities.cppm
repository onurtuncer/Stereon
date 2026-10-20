// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>
#include <string_view>
#include <vector>

export module Stereon.Topo:Capabilities;

import Stereon.Intersect;

/// The derivative capabilities of the topology library's operations
/// (ADR-0021, section 5), and the catalogue of every library together: the
/// query an optimiser makes before it asks for a product. None implements a
/// product yet.

export namespace Stereon {

inline constexpr OperationCapability TopoCatalogue[] = {
    Declare("Freeze", {Dependency::InputGeometry, Dependency::Topology},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter},
            "the topology is fixed: a derivative moves geometry, never connectivity, and the "
            "frozen shape stays valid along the step"),
    Declare("Commit", {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter,
             TrackingConvention::MovingTrimCoordinates},
            "as Freeze; a tolerance change between two evaluations is a branch event (ADR-0021, "
            "section 8)"),
    Declare("Compact", {Dependency::Topology}, {},
            "renumbering only: no geometry moves, and samples are carried by the History, not by "
            "indices (ADR-0021, section 4)"),
    Declare("CompactShape", {Dependency::Topology}, {}, "as Compact"),
};

/// The topology library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> TopoCapabilities() noexcept
{
    return TopoCatalogue;
}

/// Every library's declarations: curve, surface, intersection, topology.
[[nodiscard]] inline std::vector<OperationCapability> KernelCapabilities()
{
    std::vector<OperationCapability> all;
    for (const std::span<const OperationCapability> catalogue :
         {CurveCapabilities(), SurfaceCapabilities(), IntersectCapabilities(),
          TopoCapabilities()}) {
        all.insert(all.end(), catalogue.begin(), catalogue.end());
    }
    return all;
}

} // namespace Stereon
