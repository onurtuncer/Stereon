// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Curve:BSplineCurve;

import Stereon.Robust;
import :BSplineBasis;
import :Concepts;

/// B-spline and NURBS curves, in space and in the parameter plane.

export namespace Stereon {

/// The curve C(t) = sum_i N_{i,p}(t) w_i P_i / sum_i N_{i,p}(t) w_i over the
/// domain [Knots[p], Knots[n+1]], with p = Degree, poles P_0, ..., P_n and
/// the basis functions N_{i,p} of the knot vector (see BSplineBasis). With
/// no weights it is the polynomial B-spline sum_i N_{i,p}(t) P_i; with
/// weights, all positive, it is a NURBS curve, which can also represent the
/// conics exactly.
///
/// The knot vector may be clamped (the end knots repeated p + 1 times, so
/// the curve starts at P_0 and ends at P_n) or not; a Bezier curve is the
/// clamped B-spline with a single span. Periodic B-splines are represented
/// by their non-periodic equivalent. The invariants are those
/// MakeBSplineCurve checks; kernel code that builds the aggregate directly
/// must establish them.
///
/// Poles and weights may be any Scalar (ADR-0021), so the derivative of a
/// point with respect to its poles is available. Knots and parameters stay
/// double: they select the active span, a discrete choice.
template <Scalar T, std::size_t D>
struct BasicBSplineCurve {
    int Degree = 1;
    std::vector<double> Knots;
    std::vector<BasicPointN<T, D>> Poles;
    std::vector<T> Weights;

    bool operator==(const BasicBSplineCurve&) const = default;
};

using BSplineCurve = BasicBSplineCurve<double, 3>;
using BSplineCurve2d = BasicBSplineCurve<double, 2>;

/// The curve after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless the knot vector suits the degree and pole count
/// (Detail::CheckKnotVector), the poles are finite and the weights are
/// either absent or one positive finite weight per pole.
[[nodiscard]] Result<BSplineCurve> MakeBSplineCurve(int degree, std::vector<double> knots,
                                                    std::vector<Point3> poles,
                                                    std::vector<double> weights = {});
[[nodiscard]] Result<BSplineCurve2d> MakeBSplineCurve(int degree, std::vector<double> knots,
                                                      std::vector<Point2> poles,
                                                      std::vector<double> weights = {});

/// True for a NURBS curve, one with weights.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr bool IsRational(const BasicBSplineCurve<T, D>& curve) noexcept
{
    return !curve.Weights.empty();
}

/// [Knots[p], Knots[n+1]].
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain(const BasicBSplineCurve<T, D>& curve)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    return Interval{curve.Knots[p], curve.Knots[curve.Poles.size()]};
}

} // namespace Stereon

namespace Stereon::Detail {

/// The work arrays of one evaluation, on the stack.
inline constexpr std::size_t BSplineBasisCount = BSplineMaxDegree + 1;

/// The evaluator preconditions of a B-spline curve (hot-path checks,
/// ADR-0016): a supported degree, consistent sizes, and a finite parameter
/// inside the domain. The full knot-vector check is MakeBSplineCurve's.
template <Scalar T, std::size_t D>
constexpr void CheckEvaluation(const BasicBSplineCurve<T, D>& curve, const double t,
                               const std::string_view function)
{
    if constexpr (EvaluatorChecks) {
        const std::size_t poles = curve.Poles.size();
        Require(curve.Degree >= 1 && curve.Degree <= BSplineMaxDegree, function,
                "an unsupported degree");
        Require(poles > static_cast<std::size_t>(curve.Degree) &&
                    curve.Knots.size() == poles + static_cast<std::size_t>(curve.Degree) + 1 &&
                    (curve.Weights.empty() || curve.Weights.size() == poles),
                function, "inconsistent knot, pole or weight counts");
        RequireFinite(t, function, "a non-finite parameter");
        Require(Contains(Domain(curve), t), function, "a parameter outside the domain");
    }
}

/// The weight of pole i: 1 for a polynomial curve.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr T WeightOf(const BasicBSplineCurve<T, D>& curve, const std::size_t i)
{
    return curve.Weights.empty() ? T{1.0} : curve.Weights[i];
}

/// The derivative of order k of a rational curve from the derivatives of its
/// homogeneous numerator A = sum N w P and denominator w = sum N w (The
/// NURBS Book, eq. 4.8): C^(k) = (A^(k) - sum_{i=1..k} binom(k, i) w^(i)
/// C^(k-i)) / w. aw and w hold the derivatives up to the degree, which are
/// all that are non-zero; c receives C^(0), ..., C^(n).
template <Scalar T, std::size_t D>
[[nodiscard]] BasicVecN<T, D> RationalDerivative(std::span<const BasicVecN<T, D>> aw,
                                                 std::span<const T> w, const int n,
                                                 std::span<BasicVecN<T, D>> c)
{
    const std::size_t known = aw.size();
    for (std::size_t k = 0; k <= static_cast<std::size_t>(n); ++k) {
        BasicVecN<T, D> v = k < known ? aw[k] : BasicVecN<T, D>{};
        double binomial = 1.0;
        for (std::size_t i = 1; i <= k && i < known; ++i) {
            binomial = binomial * static_cast<double>(k - i + 1) / static_cast<double>(i);
            v -= (T{binomial} * w[i]) * c[k - i];
        }
        c[k] = v / w[0];
    }
    return c[static_cast<std::size_t>(n)];
}

} // namespace Stereon::Detail

export namespace Stereon {

template <Scalar T, std::size_t D>
[[nodiscard]] BasicPointN<T, D> Eval(const BasicBSplineCurve<T, D>& curve, const double t)
{
    Detail::CheckEvaluation(curve, t, "Eval(BSplineCurve)");
    const std::size_t span = FindSpan(curve.Knots, curve.Degree, curve.Poles.size(), t);
    const auto p = static_cast<std::size_t>(curve.Degree);
    std::array<double, Detail::BSplineBasisCount> basis; // filled by BasisFunctions
    BasisFunctions(curve.Knots, curve.Degree, span, t, std::span{basis}.first(p + 1));
    BasicVecN<T, D> numerator{};
    T denominator{};
    for (std::size_t j = 0; j <= p; ++j) {
        const std::size_t i = span - p + j;
        const T weighted = T{basis[j]} * Detail::WeightOf(curve, i);
        numerator += weighted * ToVec(curve.Poles[i]);
        denominator = denominator + weighted;
    }
    return IsRational(curve) ? ToPoint(numerator / denominator) : ToPoint(numerator);
}

/// The n-th derivative, n >= 1, from the derivatives of the basis
/// functions (algorithm A3.2), and for a rational curve from those of its
/// homogeneous form (algorithm A4.2). A polynomial curve's derivatives
/// above its degree are zero; a rational curve's are not.
template <Scalar T, std::size_t D>
[[nodiscard]] BasicVecN<T, D> Deriv(const BasicBSplineCurve<T, D>& curve, const double t,
                                    const int n)
{
    Detail::CheckEvaluation(curve, t, "Deriv(BSplineCurve)");
    Detail::Require(n >= 1, "Deriv(BSplineCurve)", "a derivative order below 1");
    if (!IsRational(curve) && n > curve.Degree) {
        return {};
    }
    const std::size_t span = FindSpan(curve.Knots, curve.Degree, curve.Poles.size(), t);
    const auto p = static_cast<std::size_t>(curve.Degree);
    const std::size_t known = std::min(static_cast<std::size_t>(n), p) + 1;
    std::array<double, Detail::BSplineBasisCount * Detail::BSplineBasisCount> ders; // filled below
    BasisDerivatives(curve.Knots, curve.Degree, span, t, static_cast<int>(known) - 1,
                     std::span{ders}.first(known * (p + 1)));
    std::array<BasicVecN<T, D>, Detail::BSplineBasisCount> aw;
    std::array<T, Detail::BSplineBasisCount> w;
    for (std::size_t k = 0; k < known; ++k) {
        aw[k] = {};
        w[k] = T{};
        for (std::size_t j = 0; j <= p; ++j) {
            const std::size_t i = span - p + j;
            const T weighted = T{ders[k * (p + 1) + j]} * Detail::WeightOf(curve, i);
            aw[k] += weighted * ToVec(curve.Poles[i]);
            w[k] = w[k] + weighted;
        }
    }
    if (!IsRational(curve)) {
        return aw[static_cast<std::size_t>(n)];
    }
    const auto order = static_cast<std::size_t>(n);
    if (order < Detail::BSplineBasisCount) {
        std::array<BasicVecN<T, D>, Detail::BSplineBasisCount> c;
        return Detail::RationalDerivative<T, D>(std::span{aw}.first(known),
                                                std::span{w}.first(known), n,
                                                std::span{c}.first(order + 1));
    }
    std::vector<BasicVecN<T, D>> c(order + 1); // orders beyond any degree: rare
    return Detail::RationalDerivative<T, D>(std::span{aw}.first(known), std::span{w}.first(known),
                                            n, c);
}

} // namespace Stereon

namespace Stereon::Detail {

/// A bound on how far Eval's rounding can move a computed point beyond the
/// convex hull of the poles, relative to the largest pole coordinate. The
/// basis functions are computed with relative errors of a few units of
/// 2^-53 per recurrence step, p steps for degree p, so they sum to 1 only
/// up to about p^2 such units; the weighted sums add p + 1 more. 8 (p + 1)^2
/// units covers both with room to spare: 5e-13 relative at degree 25.
[[nodiscard]] constexpr double PoleBoxMargin(const int degree, const double magnitude)
{
    const auto terms = static_cast<double>(degree + 1);
    return 8.0 * terms * terms * 0x1p-53 * magnitude;
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The box of the poles, widened by Detail::PoleBoxMargin. Every point of
/// the curve is a convex combination of p + 1 poles (the basis functions are
/// non-negative and sum to 1, and the weights are positive), so the exact
/// curve lies in their convex hull and hence in their box; the margin also
/// admits the rounding of Eval, so the box contains every evaluated point
/// as it does for the other built-in types.
template <std::size_t D>
[[nodiscard]] BoxN<D> BoundingBox(const BasicBSplineCurve<double, D>& curve)
{
    BoxN<D> box;
    double magnitude = 0.0;
    for (const BasicPointN<double, D>& pole : curve.Poles) {
        box = Hull(box, pole);
        for (std::size_t d = 0; d < D; ++d) {
            magnitude = std::max(magnitude, pole[d] < 0.0 ? -pole[d] : pole[d]);
        }
    }
    return Inflated(box, Detail::PoleBoxMargin(curve.Degree, magnitude));
}

/// The parameter of a point of the curve closest to p, over the whole
/// domain. The curve is split into its Bezier segments; a segment whose
/// pole box is farther from p than the best point so far is skipped, and on
/// every other the parameters where the squared distance is stationary are
/// the roots of one polynomial in Bernstein form, isolated by subdivision
/// with Descartes' rule of signs and refined by bisection. The closest of
/// those roots and the segment ends wins, the smaller parameter on a tie.
/// No starting point or convergence tolerance is involved, so the result is
/// the global closest point up to rounding (docs/nurbs.rst).
[[nodiscard]] double Project(const BSplineCurve& curve, const Point3& p);
[[nodiscard]] double Project(const BSplineCurve2d& curve, const Point2& p);

/// The same curve traversed the other way:
/// Eval(Reversed(c), a + b - t) == Eval(c, t) on the domain [a, b], up to
/// the rounding of the reflected knots.
template <Scalar T, std::size_t D>
[[nodiscard]] BasicBSplineCurve<T, D> Reversed(const BasicBSplineCurve<T, D>& curve)
{
    const Interval domain = Domain(curve);
    const double sum = domain.GetLower() + domain.GetUpper();
    BasicBSplineCurve<T, D> reversed{curve.Degree, {}, curve.Poles, curve.Weights};
    reversed.Knots.reserve(curve.Knots.size());
    for (auto it = curve.Knots.rbegin(); it != curve.Knots.rend(); ++it) {
        reversed.Knots.push_back(sum - *it);
    }
    std::ranges::reverse(reversed.Poles);
    std::ranges::reverse(reversed.Weights);
    return reversed;
}

/// The image of the curve under any affine map x: B-splines are affinely
/// invariant, because their basis functions sum to 1, so transforming the
/// poles transforms every point. The parameter and the weights are kept.
template <Scalar T>
[[nodiscard]] BasicBSplineCurve<T, 3> Transformed(const BasicBSplineCurve<T, 3>& curve,
                                                  const BasicTransform<T>& x)
{
    BasicBSplineCurve<T, 3> moved = curve;
    for (BasicPoint3<T>& pole : moved.Poles) {
        pole = Transformed(pole, x);
    }
    return moved;
}

/// The curve with the knot t inserted `times` times (Boehm's algorithm,
/// A5.1): the same point set and parametrisation with more poles.
/// InvalidInput unless t lies in the domain, times >= 1 and the multiplicity
/// of t stays at most the degree.
[[nodiscard]] Result<BSplineCurve> InsertKnot(const BSplineCurve& curve, double t, int times = 1);
[[nodiscard]] Result<BSplineCurve2d> InsertKnot(const BSplineCurve2d& curve, double t,
                                                int times = 1);

/// The curve's Bezier segments, in order: one clamped single-span curve of
/// the same degree per non-empty knot span of the domain, which together
/// trace the curve (A5.6, by knot insertion).
[[nodiscard]] std::vector<BSplineCurve> ToBezierSegments(const BSplineCurve& curve);
[[nodiscard]] std::vector<BSplineCurve2d> ToBezierSegments(const BSplineCurve2d& curve);

/// The two pieces of the curve on either side of t, each with the same
/// parametrisation as the original. InvalidInput unless t lies strictly
/// inside the domain.
[[nodiscard]] Result<std::pair<BSplineCurve, BSplineCurve>> Split(const BSplineCurve& curve,
                                                                  double t);
[[nodiscard]] Result<std::pair<BSplineCurve2d, BSplineCurve2d>> Split(const BSplineCurve2d& curve,
                                                                      double t);

/// The same curve with its domain mapped affinely onto `domain`:
/// Eval(Reparametrized(c, [a', b']), a' + (t - a)(b' - a')/(b - a)) ==
/// Eval(c, t). InvalidInput unless the new domain is finite and non-empty.
[[nodiscard]] Result<BSplineCurve> Reparametrized(const BSplineCurve& curve, Interval domain);
[[nodiscard]] Result<BSplineCurve2d> Reparametrized(const BSplineCurve2d& curve, Interval domain);

static_assert(CurveGeometry<BSplineCurve>);
static_assert(Curve2dGeometry<BSplineCurve2d>);

} // namespace Stereon
