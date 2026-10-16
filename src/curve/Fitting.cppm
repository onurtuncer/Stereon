// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

export module Stereon.Curve:Fitting;

import Stereon.Robust;
import :BSplineCurve;

/// Fitting B-spline curves to points (The NURBS Book, chapter 9): global
/// interpolation, least-squares approximation with a given number of poles,
/// and approximation to a budget. Every fit reports the parameters it gave
/// the points and the deviation it achieved at them, separately from what
/// was asked (ADR-0028).

export namespace Stereon {

/// How the points are given parameters in [0, 1] before fitting
/// (The NURBS Book, section 9.2.1). Chord length suits evenly sampled data;
/// centripetal, which uses the square roots of the chord lengths, follows
/// sharp turns more faithfully; uniform ignores the spacing.
enum class Parametrization : std::uint8_t {
    Uniform,
    ChordLength,
    Centripetal,
};

/// A fitted curve with what it achieved.
template <std::size_t D>
struct BasicCurveFit {
    BasicBSplineCurve<double, D> Curve;
    /// The parameter given to each point, in [0, 1].
    std::vector<double> Parameters;
    /// The largest distance |C(t_k) - Q_k| over the points, in millimetres.
    /// Between the points the data say nothing, so this is the deviation
    /// the fit can certify, not a bound on the curve's distance from
    /// whatever the points were sampled from.
    double MaxDeviation = 0.0;
};

using CurveFit = BasicCurveFit<3>;
using CurveFit2d = BasicCurveFit<2>;

/// The curve of the given degree through all the points, in order (algorithm
/// A9.1): parameters by `parametrization`, knots by averaging the parameters
/// (eq. 9.8), which keeps the collocation matrix non-singular, and poles
/// from a banded linear system. The domain is [0, 1] and the curve is
/// clamped, so it starts at the first point and ends at the last.
/// InvalidInput unless 1 <= degree <= BSplineMaxDegree, there are at least
/// degree + 1 finite points, and, for chord-length and centripetal
/// parameters, no two consecutive points coincide.
[[nodiscard]] Result<CurveFit>
Interpolate(std::span<const Point3> points, int degree,
            Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
Interpolate(std::span<const Point2> points, int degree,
            Parametrization parametrization = Parametrization::ChordLength);

/// The curve of the given degree and pole count closest to the points in
/// the least-squares sense, with its ends fixed at the first and last point
/// (algorithm A9.7): the interior poles solve the normal equations
/// N^T N P = R, a banded symmetric positive-definite system. Knots are
/// placed so that every span holds parameters (eq. 9.69). InvalidInput
/// unless the points suit Interpolate and degree + 1 <= poleCount <=
/// points.size().
[[nodiscard]] Result<CurveFit>
FitLeastSquares(std::span<const Point3> points, int degree, std::size_t poleCount,
                Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
FitLeastSquares(std::span<const Point2> points, int degree, std::size_t poleCount,
                Parametrization parametrization = Parametrization::ChordLength);

/// The least-squares fit with the fewest poles whose deviation at the points
/// is within the budget: pole counts are tried by doubling from degree + 1
/// and then by bisection, so the count is the smallest that passes among
/// those tried, and with as many poles as points the fit interpolates.
/// InvalidInput for invalid points or a budget that is not finite and
/// positive; ToleranceExceeded, with the budget and the best deviation
/// reached, if rounding keeps even the interpolant above the budget.
[[nodiscard]] Result<CurveFit>
Approximate(std::span<const Point3> points, int degree, ApproximationBudget budget,
            Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
Approximate(std::span<const Point2> points, int degree, ApproximationBudget budget,
            Parametrization parametrization = Parametrization::ChordLength);

} // namespace Stereon

export namespace Stereon::Detail {

/// The parameters of the points in [0, 1] (section 9.2.1), or empty when a
/// chord is zero under chord-length or centripetal parametrisation.
[[nodiscard]] std::vector<double> ParametersOf(std::span<const Point3> points,
                                               Parametrization parametrization);

/// Knots of degree p for interpolation at these parameters, by averaging
/// (eq. 9.8): clamped on [0, 1], interior knot j + p the mean of
/// parameters j, ..., j + p - 1.
[[nodiscard]] std::vector<double> AveragedKnots(std::span<const double> parameters, int degree);

/// The poles of the clamped B-spline with these knots and degree that
/// passes through each point at its parameter: the collocation system
/// sum_i N_i(t_k) P_i = Q_k, banded with half-bandwidth p, solved by
/// Gaussian elimination without pivoting, which is stable for B-spline
/// collocation matrices because they are totally positive (de Boor).
/// Surfaces interpolate their rows and columns with it.
[[nodiscard]] std::vector<Point3> InterpolatingPoles(std::span<const Point3> points,
                                                     std::span<const double> parameters,
                                                     std::span<const double> knots, int degree);

} // namespace Stereon::Detail
