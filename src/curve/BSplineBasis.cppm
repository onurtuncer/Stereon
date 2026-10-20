// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

export module Stereon.Curve:BSplineBasis;

import Stereon.Robust;

/// Knot vectors and B-spline basis functions, shared by B-spline curves and
/// surfaces (Piegl and Tiller, The NURBS Book, chapter 2).
///
/// A B-spline of degree p with n + 1 poles has a non-decreasing knot vector
/// U = (u_0, ..., u_m) with m = n + p + 1. Its parameter domain is
/// [u_p, u_{n+1}]: on every non-empty knot span [u_i, u_{i+1}) inside it,
/// exactly the p + 1 basis functions N_{i-p,p}, ..., N_{i,p} are non-zero.
/// The functions here are double only: knots and parameters decide which
/// span is active, which is a discrete choice (ADR-0021), while the poles a
/// curve combines them with may be any Scalar.
///
/// They are defined in BSpline.cpp: non-inline functions keep the module
/// interface free of instantiations of core's templates (see
/// docs/toolchain/gcc16-feature-probe.md).

export namespace Stereon {

/// The highest degree a B-spline may have. 25 is OCCT's
/// Geom_BSplineCurve::MaxDegree(), so geometry exchanged with OCCT always
/// fits; the bound lets the evaluators keep their work arrays on the stack.
inline constexpr int BSplineMaxDegree = 25;

/// The index i of the knot span [u_i, u_{i+1}) that contains t, with
/// degree <= i <= poleCount - 1 and u_i < u_{i+1} (algorithm A2.1). A t
/// at or beyond the end of the domain gets the last non-empty span, and one
/// at or before its start the first, so the end points of the domain belong
/// to the spans inside it. The knot vector must be valid for the degree and
/// pole count (CheckKnotVector).
[[nodiscard]] std::size_t FindSpan(std::span<const double> knots, int degree, std::size_t poleCount,
                                   double t) noexcept;

/// The degree + 1 basis functions N_{span-degree,degree}(t), ...,
/// N_{span,degree}(t) that are non-zero on the span, into `values`
/// (algorithm A2.2). They are non-negative and sum to 1.
void BasisFunctions(std::span<const double> knots, int degree, std::size_t span, double t,
                    std::span<double> values) noexcept;

/// The basis functions of every degree j = 0, ..., degree that are non-zero
/// on the span: table[j * (degree + 1) + r] = N_{span-j+r,j}(t) for r <= j.
/// They are the intermediate levels of BasisFunctions' triangle, kept: the
/// k-th derivative of a B-spline of degree p is a B-spline of degree p - k
/// on the same knots, so the levels below p are what derivatives need
/// (see Deriv(BSplineCurve)). `table` holds (degree + 1)^2 entries.
void BasisFunctionsOfAllDegrees(std::span<const double> knots, int degree, std::size_t span,
                                double t, std::span<double> table) noexcept;

/// The basis functions of the span and their derivatives up to `order`:
/// values[k * (degree + 1) + j] is the k-th derivative of
/// N_{span-degree+j,degree} at t (algorithm A2.3). Derivatives above the
/// degree are zero. `values` holds (order + 1) * (degree + 1) entries.
void BasisDerivatives(std::span<const double> knots, int degree, std::size_t span, double t,
                      int order, std::span<double> values) noexcept;

/// How many times `value` occurs in the knot vector.
[[nodiscard]] std::size_t KnotMultiplicity(std::span<const double> knots, double value) noexcept;

/// The full knot vector from distinct knot values and their multiplicities,
/// the form STEP and OCCT store: ExpandKnots({0, 1, 2}, {3, 1, 3}) is
/// (0, 0, 0, 1, 2, 2, 2). InvalidInput unless the two have the same length,
/// the values are finite and strictly increasing, and every multiplicity is
/// positive.
[[nodiscard]] Result<std::vector<double>> ExpandKnots(std::span<const double> values,
                                                      std::span<const int> multiplicities);

} // namespace Stereon

export namespace Stereon::Detail {

/// InvalidInput, naming `where`, unless the degree and knots describe a
/// B-spline with poleCount poles: 1 <= degree <= BSplineMaxDegree,
/// poleCount >= degree + 1, poleCount + degree + 1 knots, all finite and
/// non-decreasing, a non-empty domain [u_p, u_n+1], no knot value repeated
/// more than degree + 1 times, and no value strictly inside the domain
/// repeated more than degree times, so that the curve is continuous.
[[nodiscard]] Result<void> CheckKnotVector(int degree, std::span<const double> knots,
                                           std::size_t poleCount, std::string_view where);

/// The real roots in `range` of the polynomial sum_k coefficients[k] t^k,
/// in increasing order. The polynomial is converted to the Bernstein basis
/// on the range and its roots are isolated by de Casteljau subdivision with
/// Descartes' rule of signs, then refined by bisection, as in
/// Project(BSplineCurve): simple roots are all found, to adjacent doubles. A
/// root of even multiplicity, where the polynomial touches zero without
/// changing sign, may be reported as a nearby parameter or missed; callers
/// that minimise a distance evaluate every candidate anyway. The range must
/// be finite and non-empty.
[[nodiscard]] std::vector<double> PolynomialRoots(std::span<const double> coefficients,
                                                  Interval range);

/// The range of weights a rational B-spline may have. Evaluation forms the
/// weighted poles w P and the sums of w N over a span; a weight of 1e308
/// overflows them although the curve itself is well inside the model, and
/// one of 1e-320 loses all precision. Within [2^-256, 2^256] neither
/// happens for any pole inside the size box, and scaling every weight by a
/// common power of two in the range changes no computed point, bit for bit.
inline constexpr double MinWeight = 0x1p-256;
inline constexpr double MaxWeight = 0x1p256;

/// InvalidInput unless there are no weights (a polynomial B-spline) or one
/// finite, positive weight per pole, each in [MinWeight, MaxWeight].
/// Positive weights keep a rational B-spline inside the convex hull of its
/// poles and free of poles of its own.
[[nodiscard]] Result<void> CheckWeights(std::span<const double> weights, std::size_t poleCount,
                                        std::string_view where);

} // namespace Stereon::Detail
