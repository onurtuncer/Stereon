// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Curve:Capabilities;

import Stereon.Robust;

/// The derivative capabilities of the curve library's operations (ADR-0021,
/// section 5): one declaration per public operation, by function name, its
/// overloads together. Eval of a B-spline curve implements the forward and
/// reverse products in its poles and weights (EvalJvp, EvalVjp); every other
/// declaration has an empty product set, and an evaluation asked of it is
/// NotSupported. tools/check_capabilities.py checks that every public
/// function of the library is declared here or listed there as no
/// operation (an evaluator, a query, a utility).

export namespace Stereon {

inline constexpr OperationCapability CurveCatalogue[] = {
    // Evaluation, differentiated in the defining values (ADR-0021, stage 2).
    Declare("Eval", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "B-spline curves only, through EvalJvp and EvalVjp: at fixed parameters, knots, "
            "degree and pole count, smooth in the poles and in weights that stay positive; "
            "other curve types answer NotSupported",
            {DerivativeProduct::Jvp, DerivativeProduct::Vjp}),
    // Construction from defining values.
    Declare("MakeLine", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "smooth in the origin and direction; the direction stays non-zero"),
    Declare("MakeCircle", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "smooth in the centre, radius and frame; the radius stays positive and the frame "
            "orthonormal"),
    Declare("MakeEllipse", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "smooth in the centre, radii and frame; the radii stay positive and the frame "
            "orthonormal"),
    Declare("MakeHyperbola", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "smooth in the centre, radii and frame; the radii stay positive and the frame "
            "orthonormal"),
    Declare("MakeParabola", {Dependency::DefiningValues}, {TrackingConvention::FixedEdgeParameter},
            "smooth in the vertex, focal length and frame; the focal length stays positive"),
    Declare("MakeBSplineCurve", {Dependency::DefiningValues},
            {TrackingConvention::FixedEdgeParameter},
            "smooth in the poles and weights at a fixed knot vector and degree; knots and "
            "multiplicities are structural (ADR-0021, section 4)"),
    Declare("MakeTrimmedCurve", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "the basis's regime, with the trim range held fixed or moved as a parameter"),
    Declare("ConicFromImplicit", {Dependency::DefiningValues},
            {TrackingConvention::FixedEdgeParameter},
            "within one conic type: a change of type (ellipse, parabola, hyperbola) is a branch "
            "event"),
    Declare("PiecewisePolynomialToBSpline", {Dependency::DefiningValues},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the coefficients at fixed breakpoints"),
    // Conversion, refinement and rearrangement of a curve.
    Declare("ToBSpline", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "at a fixed segment count; parameter correspondence holds at the knots for conics"),
    Declare("ToBezierSegments", {Dependency::InputGeometry},
            {TrackingConvention::FixedEdgeParameter}, "linear in the poles at a fixed knot vector"),
    Declare("InsertKnot", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles; the inserted knot stays strictly inside a span"),
    Declare("Split", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "the split parameter stays strictly inside the domain"),
    Declare("Reparametrized", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter}, "an affine map of the parameter"),
    Declare("Clamped", {Dependency::InputGeometry}, {TrackingConvention::FixedEdgeParameter},
            "linear in the poles at a fixed knot vector: knot insertion at the domain ends"),
    Declare("AsPeriodic", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles at a fixed knot vector; the end poles stay within the tolerance"),
    Declare("ElevatedDegree", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles at a fixed knot vector and target degree"),
    Declare("RemoveKnot",
            {Dependency::InputGeometry, Dependency::OperationValues, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles at a fixed knot vector; how many copies the tolerance lets go "
            "is a branch event"),
    Declare("Concatenated", {Dependency::InputGeometry, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter},
            "linear in the poles of both curves at fixed knot vectors; the junction stays within "
            "the tolerance"),
    Declare("Reversed", {Dependency::InputGeometry}, {TrackingConvention::FixedEdgeParameter},
            "a fixed reflection of the parameter"),
    Declare("Transformed", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter}, "linear in the transform's entries"),
    // Fitting: solves whose discrete choices are branch events.
    Declare("Interpolate", {Dependency::InputPoints}, {TrackingConvention::FixedEdgeParameter},
            "at fixed parameters and knots; parametrisation by chord length moves them with the "
            "points"),
    Declare("InterpolatePeriodic", {Dependency::InputPoints},
            {TrackingConvention::FixedEdgeParameter},
            "at fixed parameters and knots, as Interpolate"),
    Declare("FitLeastSquares", {Dependency::InputPoints}, {TrackingConvention::FixedEdgeParameter},
            "at a fixed pole count, parameters and knots; the normal equations stay well "
            "conditioned"),
    Declare("FitLeastSquaresPeriodic", {Dependency::InputPoints},
            {TrackingConvention::FixedEdgeParameter}, "as FitLeastSquares"),
    Declare("Approximate", {Dependency::InputPoints, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter},
            "at the pole count the budget chose; a change of that count is a branch event"),
    Declare("ApproximatePeriodic", {Dependency::InputPoints, Dependency::Tolerances},
            {TrackingConvention::FixedEdgeParameter}, "as Approximate"),
    // Implicit constructions (ADR-0021, section 6).
    Declare("Project", {Dependency::InputGeometry, Dependency::InputPoints},
            {TrackingConvention::UniqueProjection},
            "a unique closest point, interior to the domain, where the curve is not tangent to a "
            "circle about the point (no ambiguous or degenerate projection)"),
    Declare("ProjectAll", {Dependency::InputGeometry, Dependency::InputPoints},
            {TrackingConvention::UniqueProjection}, "as Project, point by point"),
    Declare("ParameterAtArcLength", {Dependency::InputGeometry, Dependency::OperationValues},
            {TrackingConvention::FixedEdgeParameter},
            "a regular curve: the speed stays away from zero along the solved length"),
    // Functionals.
    Declare("ArcLength", {Dependency::InputGeometry, Dependency::OperationValues}, {},
            "a regular curve over a fixed range; quadrature tolerance held fixed"),
};

/// The curve library's declarations.
[[nodiscard]] inline std::span<const OperationCapability> CurveCapabilities() noexcept
{
    return CurveCatalogue;
}

} // namespace Stereon
