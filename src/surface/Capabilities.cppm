// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Surface:Capabilities;

import Stereon.Curve;

/// The derivative capabilities of the surface library's operations
/// (ADR-0021, section 5), as Stereon.Curve:Capabilities declares the
/// curve library's: none implements a product yet.

export namespace Stereon {

inline constexpr OperationCapability SurfaceCatalogue[] = {
    // Evaluation, which the curve library differentiates for B-splines; the
    // surface products are not implemented yet.
    Declare("Eval", {Dependency::DefiningValues}, {TrackingConvention::FixedSurfaceCoordinates},
            "at fixed (u, v), knots, degrees and pole counts, smooth in the defining values"),
    // Construction from defining values.
    Declare("MakePlane", {Dependency::DefiningValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the origin and frame; the frame stays orthonormal"),
    Declare("MakeCylinder", {Dependency::DefiningValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the origin, radius and frame; the radius stays positive"),
    Declare("MakeCone", {Dependency::DefiningValues}, {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the origin, radius, semi-angle and frame; samples stay off the apex, a "
            "coordinate singularity (ADR-0009)"),
    Declare("MakeSphere", {Dependency::DefiningValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the centre, radius and frame; samples stay off the poles, coordinate "
            "singularities"),
    Declare("MakeTorus", {Dependency::DefiningValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the centre, radii and frame; the minor radius stays below the major"),
    Declare("MakeBSplineSurface", {Dependency::DefiningValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "smooth in the poles and weights at fixed knot vectors and degrees"),
    Declare("MakeTrimmedSurface", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "the basis's regime, with the box held fixed or moved as a parameter"),
    Declare("MakeExtrusion", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "the profile's regime, linear in the direction"),
    Declare("MakeRevolution", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "the profile's regime; samples stay off the axis, where the surface degenerates"),
    Declare("MakeExtrusionSurface", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "linear in the profile's poles and the direction over a fixed range"),
    Declare("MakeRevolutionSurface", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "at a fixed arc count over a fixed angular range"),
    // Conversion, refinement and rearrangement of a surface.
    Declare("ToBSpline", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "at a fixed patch layout; parameter correspondence holds at the knots for rational "
            "arcs"),
    Declare("ToBezierPatches", {Dependency::InputGeometry},
            {TrackingConvention::FixedSurfaceCoordinates},
            "linear in the poles at fixed knot vectors"),
    Declare("InsertKnot", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates},
            "linear in the poles; the inserted knot stays strictly inside a span"),
    Declare("IsoCurve", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles; the iso value stays inside its direction's domain"),
    Declare("Reversed", {Dependency::InputGeometry}, {TrackingConvention::FixedSurfaceCoordinates},
            "a fixed reflection of u"),
    Declare("Transformed", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedSurfaceCoordinates}, "linear in the transform's entries"),
    // Fitting.
    Declare("InterpolateSurface", {Dependency::InputPoints},
            {TrackingConvention::FixedSurfaceCoordinates}, "at fixed parameters and knots"),
    Declare("FitSurfaceLeastSquares", {Dependency::InputPoints},
            {TrackingConvention::FixedSurfaceCoordinates},
            "at fixed pole counts, parameters and knots; well conditioned normal equations"),
    Declare("ApproximateSurface", {Dependency::InputPoints, Dependency::Tolerances},
            {TrackingConvention::FixedSurfaceCoordinates},
            "at the pole counts the budget chose; a change of them is a branch event"),
    // Implicit constructions (ADR-0021, section 6).
    Declare("Project", {Dependency::InputGeometry, Dependency::InputPoints},
            {TrackingConvention::UniqueProjection},
            "a unique closest point, interior to the domain, with a regular parametrisation and "
            "a non-singular Hessian of the squared distance"),
    Declare("ProjectAll", {Dependency::InputGeometry, Dependency::InputPoints},
            {TrackingConvention::UniqueProjection}, "as Project, point by point"),
};

/// The surface library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> SurfaceCapabilities() noexcept
{
    return SurfaceCatalogue;
}

} // namespace Stereon
