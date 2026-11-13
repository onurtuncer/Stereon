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
/// line, a circle, an ellipse or a B-spline curve against a plane, a cylinder, a cone, a sphere,
/// a torus, an ellipsoid or a B-spline surface, and for branch points of
/// two of them, and for points where two plane curves meet (IntersectJvp,
/// IntersectVjp, through Linearize); every other pair answers NotSupported.

export namespace Stereon {

inline constexpr OperationCapability IntersectCatalogue[] = {
    Declare("Intersect", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::MovingTrimCoordinates, TrackingConvention::UniqueProjection},
            "a line, a circle, an ellipse or a B-spline curve against a plane, a cylinder, a "
            "cone, a sphere, a torus, "
            "an ellipsoid or a B-spline surface, "
            "through IntersectJvp and IntersectVjp, in the line's origin and direction, a "
            "conic's centre and radii at a fixed frame, the "
            "curve's poles and weights at fixed knots, and the surface's origin, its radii, a "
            "cone's semi-angle, an ellipsoid's semi-axes and a B-spline surface's poles and "
            "weights at fixed knots: decided transversal "
            "crossings only, NearEvent below NearEventSine; "
            "Tangential and Overlap contacts (ADR-0010) are nondifferentiable, undecided points "
            "and poles of the surface's parameters indeterminate; the point count and the "
            "tolerance stay fixed; branch points of two such surfaces under the normal-plane "
            "tracking constraint (BranchSample), certified crossing branches only; "
            "points where two plane curves (lines, circles, ellipses, B-splines) cross, in "
            "their defining values (PlaneCurveBinding); a meeting of two space curves is an "
            "event of codimension one, nondifferentiable, a separation valid with no values; "
            "other pairs answer NotSupported",
            {DerivativeProduct::Jvp, DerivativeProduct::Vjp}),
    Declare("ClosestApproach", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::UniqueProjection},
            "two lines, circles or ellipses in space (ADR-0052, addendum of 2026-11-12): the "
            "pair of closest points, proved the unique global minimiser over the requested "
            "domains, tracked through ClosestApproachJvp and ClosestApproachVjp in a line's "
            "origin and direction and a conic's centre and radii at a fixed frame; a parameter "
            "at an end of its domain is held there; a set of minimisers, a meeting and an end "
            "the minimum may leave are nondifferentiable, an undecided minimiser or a singular "
            "Hessian indeterminate, NearEvent below NearEventSine on the normalised Hessian; "
            "B-spline curves answer NotSupported",
            {DerivativeProduct::Jvp, DerivativeProduct::Vjp}),
};

/// The intersection library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> IntersectCapabilities() noexcept
{
    return IntersectCatalogue;
}

} // namespace Stereon
