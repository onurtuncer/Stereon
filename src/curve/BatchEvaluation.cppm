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

export module Stereon.Curve:BatchEvaluation;

import Stereon.Robust;
import :BSplineCurve;

/// Evaluation of a B-spline curve at many parameters at once, data-parallel
/// across the parameters (ADR-0007: std::simd, through Stereon::Simd). See
/// docs/nurbs.rst, *Batch and grid evaluation*.

export namespace Stereon {

/// points[k] = Eval(curve, parameters[k]) for every k, bit for bit: a run of
/// Simd::Width<double> consecutive parameters in one knot span is evaluated
/// in the lanes of one vector, by the same operations in the same order as
/// Eval performs on each, and any other parameter by Eval itself. The result
/// therefore does not depend on the SIMD width, the order of the
/// parameters, or how they fall into spans; sorted parameters, which share
/// spans, are the fast case. Preconditions as Eval's for every parameter,
/// and as many points as parameters.
void EvalAll(const BSplineCurve& curve, std::span<const double> parameters,
             std::span<Point3> points);
void EvalAll(const BSplineCurve2d& curve, std::span<const double> parameters,
             std::span<Point2> points);

/// Every non-empty knot span of a B-spline curve as a polynomial, for
/// EvalAllFast: on the span [a, b], with midpoint m and width h, the
/// homogeneous curve (w P, w) is sum_k c_k u^k in u = (t - m) / h, |u| <= 1/2,
/// with c_k its k-th derivative at m times h^k / k!. Built once, it serves any
/// number of evaluations. The fields are the evaluator's.
template <std::size_t D>
struct BasicSpanPolynomials {
    int Degree = 1;
    bool Rational = false;
    std::vector<double> Knots;
    std::size_t PoleCount = 0;
    /// Per knot index s, the offset of span [Knots[s], Knots[s + 1]] in
    /// Midpoints, InverseWidths and Coefficients / ((Degree + 1) (D + 1)),
    /// or PoleCount + 1 for an empty span or one outside the domain.
    std::vector<std::size_t> Slots;
    std::vector<double> Midpoints;
    std::vector<double> InverseWidths;
    /// Per span, (Degree + 1) rows of D + 1 values: c_k's coordinates and
    /// weight. The weight is 1 for k = 0 and 0 otherwise on a polynomial curve.
    std::vector<double> Coefficients;
};

using SpanPolynomials = BasicSpanPolynomials<3>;
using SpanPolynomials2d = BasicSpanPolynomials<2>;

/// The polynomial form of every non-empty span of the curve.
[[nodiscard]] SpanPolynomials ToSpanPolynomials(const BSplineCurve& curve);
[[nodiscard]] SpanPolynomials2d ToSpanPolynomials(const BSplineCurve2d& curve);

/// points[k] = the curve at parameters[k], from its span polynomials by
/// Horner's rule, the parameters of a span in the lanes of one vector. Not
/// bit for bit Eval: the polynomial form rounds differently from de Boor's,
/// and the difference is bounded in docs/nurbs.rst, *Fast batch
/// evaluation*. The result does not depend on the SIMD width, the order of
/// the parameters or how they fall into spans. Preconditions as Eval's for
/// every parameter, and as many points as parameters.
void EvalAllFast(const SpanPolynomials& polynomials, std::span<const double> parameters,
                 std::span<Point3> points);
void EvalAllFast(const SpanPolynomials2d& polynomials, std::span<const double> parameters,
                 std::span<Point2> points);

/// As above, converting the curve first; for a curve evaluated once.
void EvalAllFast(const BSplineCurve& curve, std::span<const double> parameters,
                 std::span<Point3> points);
void EvalAllFast(const BSplineCurve2d& curve, std::span<const double> parameters,
                 std::span<Point2> points);

} // namespace Stereon
