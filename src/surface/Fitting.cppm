// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <span>
#include <vector>

export module Stereon.Surface:Fitting;

import Stereon.Curve;
import :BSplineSurface;

/// Fitting B-spline surfaces to a grid of points (The NURBS Book, sections
/// 9.2.5 and 9.4.3): interpolation, least squares with given pole counts,
/// and approximation to a budget. As for curves, every fit reports the
/// parameters it used and the deviation it achieved (ADR-0028).

export namespace Stereon {

/// A fitted surface with what it achieved.
struct SurfaceFit {
    BSplineSurface Surface;
    /// The parameter of each grid column along u and each row along v.
    std::vector<double> ParametersU;
    std::vector<double> ParametersV;
    /// The largest distance |S(u_k, v_l) - Q_kl| over the grid, in
    /// millimetres; as for curves, a deviation at the data, not between.
    double MaxDeviation = 0.0;
};

/// The surface of degrees (degreeU, degreeV) through every point of a grid
/// of countU x countV points, stored row by row like surface poles (point
/// (k, l) at grid[k * countV + l], k along u). Parameters are the means of
/// the lines' parameters (algorithm A9.3), skipping lines whose points all
/// coincide, such as the collapsed row at a sphere's pole; knots are
/// averaged as for curves; the poles come from interpolating every line
/// along u and then every line along v (A9.4). The domain is [0, 1]^2.
/// InvalidInput unless the degrees are valid, countU * countV finite points
/// form at least (degreeU + 1) x (degreeV + 1), and each direction has a
/// usable line.
[[nodiscard]] Result<SurfaceFit>
InterpolateSurface(std::span<const Point3> grid, std::size_t countU, std::size_t countV,
                   int degreeU, int degreeV,
                   Parametrization parametrization = Parametrization::ChordLength);

/// The surface of degrees (degreeU, degreeV) with poleCountU x poleCountV
/// poles closest to the grid in the least-squares sense (section 9.4.3):
/// parameters as for InterpolateSurface, knots spread as for
/// FitLeastSquares, and the poles from a least-squares curve fit of every
/// grid line along u followed by one of every resulting line along v, each
/// with its ends fixed. Because the grid shares its parameters along each
/// direction, the two stages give the exact least-squares minimiser over the
/// interior poles once the boundary poles are fixed, and the boundary poles
/// are those of the least-squares fits of the four boundary lines, so the
/// corners are exact. A direction with as many poles as grid lines is
/// interpolated instead, as Approximate does for curves. InvalidInput unless
/// the grid suits InterpolateSurface and each pole count is between its
/// degree + 1 and its grid count.
[[nodiscard]] Result<SurfaceFit>
FitSurfaceLeastSquares(std::span<const Point3> grid, std::size_t countU, std::size_t countV,
                       int degreeU, int degreeV, std::size_t poleCountU, std::size_t poleCountV,
                       Parametrization parametrization = Parametrization::ChordLength);

/// The least-squares surface with the fewest poles among those tried whose
/// deviation at the grid is within the budget. The pole count along u is
/// searched first, by doubling from degreeU + 1 and then bisection as
/// Approximate does for curves, with every line along v interpolated so that
/// only u limits the fit; then the count along v with that u count. The
/// search stops at the first passing count of each direction among those
/// tried, so a surface that is simple along one direction gets few poles
/// there whatever the other needs. InvalidInput for an invalid grid or a
/// budget that is not finite and positive; ToleranceExceeded, with the
/// budget and the best deviation reached, if even the interpolant misses it.
[[nodiscard]] Result<SurfaceFit>
ApproximateSurface(std::span<const Point3> grid, std::size_t countU, std::size_t countV,
                   int degreeU, int degreeV, ApproximationBudget budget,
                   Parametrization parametrization = Parametrization::ChordLength);

} // namespace Stereon
