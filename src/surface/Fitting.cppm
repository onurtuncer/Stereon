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

/// Fitting B-spline surfaces to a grid of points (The NURBS Book, section
/// 9.2.5).

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

} // namespace Stereon
