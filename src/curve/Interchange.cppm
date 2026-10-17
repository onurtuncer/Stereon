// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

export module Stereon.Curve:Interchange;

import Stereon.Robust;
import :BSplineCurve;
import :Erased;

/// Conversions from the curve representations of the exchange formats that
/// the kernel does not keep as types of their own (docs/interchange.rst):
///
/// - piecewise polynomials in the power basis, IGES's parametric spline
///   (entity 112), to B-splines, exactly;
/// - conics given by their implicit equation, IGES's conic arc (entity
///   104), to the Ellipse2d, Hyperbola2d or Parabola2d they describe.

export namespace Stereon {

/// A piecewise polynomial converted to a B-spline, with how well its pieces
/// met.
template <std::size_t D>
struct BasicPiecewiseConversion {
    /// Piece i on [t_i, t_{i+1}] as one Bezier segment, joined to the next
    /// by a knot of multiplicity equal to the degree.
    BasicBSplineCurve<double, D> Curve;
    /// The largest distance, in millimetres, between the end of one piece
    /// and the start of the next, over the joins. A join that does not
    /// meet gets the mean of the two points as its pole; the curve is then
    /// within MaxJoinGap / 2 of the data near the join and exact elsewhere.
    double MaxJoinGap = 0.0;
};

using PiecewiseConversion = BasicPiecewiseConversion<3>;
using PiecewiseConversion2d = BasicPiecewiseConversion<2>;

/// The curve that is, on [breakpoints[i], breakpoints[i + 1]],
///
///   x(t) = sum_j pieces[i][j] s^j,   s = t - breakpoints[i],
///
/// IGES 112's form: the coefficients are in the power basis in the local
/// parameter s, and pieces[i][0] is the position of the piece's start as a
/// vector from the origin. Pieces of lower degree are raised to the
/// highest; the degree is at least 1. The conversion to Bernstein form is
/// exact up to rounding, and it keeps the parameter.
///
/// InvalidInput unless there is at least one piece, every piece has at
/// least one coefficient, there is one more breakpoint than pieces, the
/// breakpoints are finite and strictly increasing, and every coefficient
/// is finite.
[[nodiscard]] Result<PiecewiseConversion>
PiecewisePolynomialToBSpline(std::span<const double> breakpoints,
                             std::span<const std::vector<Vec3>> pieces);
[[nodiscard]] Result<PiecewiseConversion2d>
PiecewisePolynomialToBSpline(std::span<const double> breakpoints,
                             std::span<const std::vector<Vec2>> pieces);

/// The coefficients of A x^2 + B xy + C y^2 + D x + E y + F = 0, in IGES
/// 104's order.
struct ImplicitConic {
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    double E = 0.0;
    double F = 0.0;
};

/// How close to zero, relative to the largest, the smaller eigenvalue of
/// the quadratic part may be for the conic to be read as a parabola, and
/// how close to zero, relative to the terms that cancel in it, a quantity
/// may be before the conic is degenerate (ADR-0028). An exporter writes a
/// parabola's coefficients rounded, so its discriminant B^2 - 4AC is
/// rarely exactly zero; the allowance is far above that rounding (about
/// 1e-16) and far below the eigenvalue ratio of any ellipse or hyperbola a
/// CAD model holds, whose axis ratio would be 1e-5 at this bound.
inline constexpr double ConicAllowance = 0x1p-33;

/// The conic the equation describes, in the plane: an Ellipse2d, a
/// Hyperbola2d or a Parabola2d, whose points satisfy it up to rounding. An
/// ellipse runs counter-clockwise, as IGES 104's arcs do, and has its
/// major axis as XAxis; an ellipse with equal radii stays an ellipse. A
/// hyperbola is one branch: the one on the side of `branch` when it is
/// given (IGES 104's start point), otherwise the one along the positive
/// direction of its transverse axis's eigenvector as computed. A parabola
/// opens along XAxis.
///
/// InvalidInput when a coefficient is not finite, and when the equation
/// describes no curve the kernel can hold: no point (x^2 + y^2 + 1 = 0), a
/// single point, a line or a pair of lines, or all of the plane.
[[nodiscard]] Result<Curve2d> ConicFromImplicit(const ImplicitConic& conic,
                                                std::optional<Point2> branch = std::nullopt);

/// The equation's value at p, scaled by its largest coefficient: zero on
/// the conic, and comparable between equations of different scale.
[[nodiscard]] double ImplicitResidual(const ImplicitConic& conic, const Point2& p);

} // namespace Stereon
