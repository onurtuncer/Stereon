// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>
#include <vector>

export module Stereon.Boolean:Capabilities;

import Stereon.Build;
import Stereon.Intersect;

/// The derivative capabilities of the Boolean library's operations
/// (ADR-0021, section 5), and the catalogue of every library together: the
/// Boolean library sits above all the others, so only it can collect them.
/// None implements a product yet. A Boolean derivative, once one exists, is
/// that of one regular local construction (ADR-0021, section 8): a change in
/// the regions kept, a tangency or a coincidence is a branch event.

export namespace Stereon {

inline constexpr OperationCapability BooleanCatalogue[] = {
    Declare("Boolean", {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter,
             TrackingConvention::MovingTrimCoordinates},
            "every face pair meets transversally, the regions kept and their classification stay "
            "the same, and no tolerance grows: the intersection edges move with the operands, and "
            "a region appearing or vanishing, a tangency or a coincidence is a branch event"),
    Declare("ImprintBodies",
            {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::FixedEdgeParameter,
             TrackingConvention::MovingTrimCoordinates},
            "as Boolean, before classification: the set of intersecting face pairs and the "
            "number of branches on each stay the same"),
    Declare("ClipToFace", {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter, TrackingConvention::MovingTrimCoordinates},
            "the branch crosses the face's loops transversally and away from their vertices; a "
            "crossing that appears, vanishes or reaches a vertex is a branch event"),
    Declare("NetworkOf", {Dependency::InputGeometry, Dependency::Topology},
            {TrackingConvention::MovingTrimCoordinates},
            "as ClipToFace: the pieces and their junctions stay the same"),
    Declare("ImprintClip",
            {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::MovingTrimCoordinates},
            "as ClipToFace, and the edges split by the clip are split at moving parameters"),
    Declare("ImprintNetwork",
            {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::MovingTrimCoordinates},
            "as ImprintClip, for one network"),
    Declare("ImprintNetworks",
            {Dependency::InputGeometry, Dependency::Topology, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates, TrackingConvention::MovingTrimCoordinates},
            "as ImprintNetwork for each network; an edge two networks split stays split at the "
            "same vertex"),
};

/// The Boolean library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> BooleanCapabilities() noexcept
{
    return BooleanCatalogue;
}

/// Every library's declarations: curve, surface, intersection, topology,
/// construction and Boolean.
[[nodiscard]] inline std::vector<OperationCapability> KernelCapabilities()
{
    std::vector<OperationCapability> all;
    for (const std::span<const OperationCapability> catalogue :
         {CurveCapabilities(), SurfaceCapabilities(), IntersectCapabilities(), TopoCapabilities(),
          BuildCapabilities(), BooleanCapabilities()}) {
        all.insert(all.end(), catalogue.begin(), catalogue.end());
    }
    return all;
}

} // namespace Stereon
