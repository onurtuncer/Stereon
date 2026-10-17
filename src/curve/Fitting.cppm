// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <string_view>
#include <vector>

export module Stereon.Curve:Fitting;

import Stereon.Robust;
import :BSplineCurve;

/// Fitting B-spline curves to points (The NURBS Book, chapter 9): global
/// interpolation, least-squares approximation with a given number of poles,
/// and approximation to a budget, each with or without prescribed end
/// derivatives. Every fit reports the parameters it gave the points and the
/// deviation it achieved at them, separately from what was asked (ADR-0028).

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

/// Derivatives prescribed at the two ends of a fit, with respect to the
/// curve's parameter on [0, 1] (The NURBS Book, section 9.2.2). Their
/// lengths matter as much as their directions: a chord-length parametrised
/// curve of total chord L runs at a speed of about L, so a unit tangent
/// direction T is prescribed as L T; a shorter vector slows the curve into
/// its end and a longer one makes it overshoot.
template <std::size_t D>
struct BasicEndDerivatives {
    BasicVecN<double, D> Start;
    BasicVecN<double, D> End;
};

using EndDerivatives = BasicEndDerivatives<3>;
using EndDerivatives2d = BasicEndDerivatives<2>;

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

/// The interpolant with its end derivatives prescribed as well (algorithm
/// A9.2, for any degree): two more poles than points, the second and the
/// second to last fixed by the derivatives, P_1 = Q_0 + u_{p+1} D_0 / p and
/// P_{n-1} = Q_m - (1 - u_n) D_m / p (eq. 9.10), with knots that average p
/// parameters shifted back by one (eq. 9.11), which keeps the remaining
/// collocation system non-singular. InvalidInput as for Interpolate, for a
/// non-finite derivative, and for degree 1, whose shifted knots would repeat
/// the ends.
[[nodiscard]] Result<CurveFit>
Interpolate(std::span<const Point3> points, int degree, EndDerivatives ends,
            Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
Interpolate(std::span<const Point2> points, int degree, EndDerivatives2d ends,
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

/// The least-squares fit with its end derivatives prescribed as well: the
/// second and the second to last pole are fixed by the derivatives as in
/// the constrained interpolant, and the poles between them solve the normal
/// equations. InvalidInput as for FitLeastSquares and the constrained
/// interpolant, and unless poleCount >= 4, so that the four fixed poles are
/// distinct.
[[nodiscard]] Result<CurveFit>
FitLeastSquares(std::span<const Point3> points, int degree, std::size_t poleCount,
                EndDerivatives ends,
                Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
FitLeastSquares(std::span<const Point2> points, int degree, std::size_t poleCount,
                EndDerivatives2d ends,
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

/// Approximation with the end derivatives prescribed as well: the same
/// search over the constrained least-squares fits, from 4 poles up to as
/// many as points, with the constrained interpolant (two more poles than
/// points) as the last resort.
[[nodiscard]] Result<CurveFit>
Approximate(std::span<const Point3> points, int degree, ApproximationBudget budget,
            EndDerivatives ends, Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
Approximate(std::span<const Point2> points, int degree, ApproximationBudget budget,
            EndDerivatives2d ends, Parametrization parametrization = Parametrization::ChordLength);

/// The closed curve of the given degree through all the points, in order,
/// returning to the first point after the last with the smoothness of the
/// spline, C^(p-1), at the seam. The points are given once: the closing
/// chord from the last point back to the first takes part in the
/// parametrisation like any other. The curve is a periodic B-spline of
/// period 1 written out as its unclamped equivalent (see BSplineCurve): the
/// domain is [0, 1], the knots repeat with period 1 and the last `degree`
/// poles repeat the first. Knots are placed at the parameters for an odd
/// degree and halfway between them for an even degree, the rule under which
/// periodic spline interpolation is uniquely solvable; the parameters of
/// the points increase from 0 (odd degree) or from half the closing chord's
/// share (even degree) and stay below 1. The poles solve a cyclic banded
/// system. InvalidInput unless 1 <= degree <= BSplineMaxDegree, there are
/// at least degree + 1 finite points and, for chord-length and centripetal
/// parameters, no two consecutive points coincide, the last and the first
/// included.
[[nodiscard]] Result<CurveFit>
InterpolatePeriodic(std::span<const Point3> points, int degree,
                    Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
InterpolatePeriodic(std::span<const Point2> points, int degree,
                    Parametrization parametrization = Parametrization::ChordLength);

/// The closed curve of the given degree and number of distinct poles
/// closest to the points in the least-squares sense: every pole is free,
/// and they solve the periodic normal equations, a cyclic banded symmetric
/// positive-definite system. Knots are spread over the period so that
/// every span holds parameters, moved back by half a step for an even
/// degree as in InterpolatePeriodic; with as many poles as points the fit
/// interpolates. The curve has poleCount + degree poles, the last `degree`
/// repeating the first. InvalidInput unless the points suit
/// InterpolatePeriodic and degree + 1 <= poleCount <= points.size().
[[nodiscard]] Result<CurveFit>
FitLeastSquaresPeriodic(std::span<const Point3> points, int degree, std::size_t poleCount,
                        Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
FitLeastSquaresPeriodic(std::span<const Point2> points, int degree, std::size_t poleCount,
                        Parametrization parametrization = Parametrization::ChordLength);

/// Approximate for a closed curve: the periodic least-squares fit with the
/// fewest distinct poles, among those tried by doubling and bisection as in
/// Approximate, whose deviation at the points is within the budget. The
/// same errors as Approximate, with InterpolatePeriodic's conditions on the
/// points.
[[nodiscard]] Result<CurveFit>
ApproximatePeriodic(std::span<const Point3> points, int degree, ApproximationBudget budget,
                    Parametrization parametrization = Parametrization::ChordLength);
[[nodiscard]] Result<CurveFit2d>
ApproximatePeriodic(std::span<const Point2> points, int degree, ApproximationBudget budget,
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

/// Knots of degree p for a least-squares fit with poleCount poles at these
/// parameters (eq. 9.69): clamped on [0, 1], with d = (m + 1) / (n - p + 1)
/// interior knot p + j blends the parameters around j d, so that every knot
/// span holds at least one parameter and the normal equations are positive
/// definite. Needs degree + 1 <= poleCount <= parameters.size().
[[nodiscard]] std::vector<double> SpreadKnots(std::span<const double> parameters, int degree,
                                              std::size_t poleCount);

/// The poles of the least-squares fit with these knots and degree: the
/// first and last are the first and last point, the rest solve the normal
/// equations over the interior points (algorithm A9.7). Surfaces fit their
/// rows and columns with it.
[[nodiscard]] std::vector<Point3> LeastSquaresPoles(std::span<const Point3> points,
                                                    std::span<const double> parameters,
                                                    std::span<const double> knots, int degree,
                                                    std::size_t poleCount);

/// The smallest pole count among those tried whose deviation is within the
/// budget: counts are tried by doubling from `fewest` until one passes or
/// `most` is reached, then by bisection between the last failure and the
/// first success. `deviationOf(count)` fits with that many poles and
/// returns its deviation, or the error that ends the search.
/// ToleranceExceeded, naming `where` and carrying the budget and the best
/// deviation reached, when even `most` poles fail.
[[nodiscard]] Result<std::size_t>
FewestPolesWithin(std::size_t fewest, std::size_t most, double budget,
                  const std::function<Result<double>(std::size_t)>& deviationOf,
                  std::string_view where);

/// The knot vector of the unclamped equivalent of a periodic B-spline of
/// degree p whose distinct knots over one period are v_0 < ... < v_n, with
/// v_n = v_0 + period: the n + 2p + 1 knots v_{n-p} - period, ...,
/// v_{n-1} - period, v_0, ..., v_n, v_1 + period, ..., v_p + period, whose
/// domain is [v_0, v_n]. `periodKnots` holds v_0, ..., v_n with n >= p + 1.
[[nodiscard]] std::vector<double> PeriodicKnotVector(std::span<const double> periodKnots,
                                                     int degree);

/// The poles of that equivalent: the n distinct poles followed by the first
/// p of them again, so that the curve on [v_0, v_n] closes with C^(p-1)
/// continuity at the seam. `poles` holds n >= p + 1 poles.
[[nodiscard]] std::vector<Point3> WrappedPoles(std::span<const Point3> poles, int degree);

} // namespace Stereon::Detail
