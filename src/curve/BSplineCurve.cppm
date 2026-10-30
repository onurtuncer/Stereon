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
#include <cmath>
#include <concepts>
#include <cstddef>
#include <memory>
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

/// The period of a curve in the periodic form (IsPeriodic), the width of its
/// domain, and 0, "not periodic", for any other: a B-spline curve is a
/// PeriodicCurve whose period its data decide, with no flag to disagree with
/// them. A curve in the form continues across either end of its domain, so
/// it evaluates at every finite parameter, carried into the domain by whole
/// periods, and its knot insertion and conversion keep the form. O(n): the
/// evaluators ask for it only for a parameter outside the domain.
[[nodiscard]] double Period(const BSplineCurve& curve);
[[nodiscard]] double Period(const BSplineCurve2d& curve);

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

/// t as the evaluators use it: a parameter outside the domain of a curve in
/// the periodic form is carried into the domain by whole periods
/// (NearestInRange); any other is returned as it is, for CheckEvaluation to
/// judge. Only a curve of doubles has a period; the others, interval and
/// dual curves among them, are evaluated inside their domains.
template <Scalar T, std::size_t D>
[[nodiscard]] double IntoDomain(const BasicBSplineCurve<T, D>& curve, const double t)
{
    if constexpr (std::same_as<T, double>) {
        const Interval domain = Domain(curve);
        if (!Contains(domain, t) && std::isfinite(t)) {
            if (const double period = Period(curve); period > 0.0) {
                return NearestInRange(t, domain, period);
            }
        }
    }
    return t;
}

/// The parameters as the batch evaluators use them: when some lie outside
/// the domain of a curve in the periodic form, all carried into the domain,
/// in `storage`; otherwise the parameters themselves, with nothing copied.
template <std::size_t D>
[[nodiscard]] std::span<const double> IntoDomain(const BasicBSplineCurve<double, D>& curve,
                                                 const std::span<const double> parameters,
                                                 std::vector<double>& storage)
{
    const Interval domain = Domain(curve);
    const bool outside =
        std::ranges::any_of(parameters, [&](const double t) { return !Contains(domain, t); });
    if (!outside || !(Period(curve) > 0.0)) {
        return parameters;
    }
    storage.reserve(parameters.size());
    for (const double t : parameters) {
        storage.push_back(IntoDomain(curve, t));
    }
    return storage;
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

/// Work arrays sized for low degrees. Vectors zero themselves on
/// construction, so arrays sized for the largest degree would clear about
/// 2 KB per evaluation; degrees up to SmallDegreeCapacity - 1, which covers
/// the curves met in practice, use arrays of this size instead.
inline constexpr std::size_t SmallDegreeCapacity = 8;

/// The derivatives A^(k) = sum N^(k) w P and w^(k) = sum N^(k) w of the
/// homogeneous curve for k < known (known <= p + 1 <= Capacity), by
/// differencing the active poles rather than differentiating the basis (The
/// NURBS Book, eq. 3.8 and A3.3): the k-th derivative of a degree-p
/// B-spline is the degree-(p - k) B-spline on the same knots whose poles are
/// Q^(k)_i = (p - k + 1) (Q^(k-1)_(i+1) - Q^(k-1)_i) / (u_(i+p+1) - u_(i+k)),
/// and its basis functions on the span are a level of the triangle that
/// BasisFunctionsOfAllDegrees computes for the point anyway. On a non-empty
/// span every denominator is positive.
template <std::size_t Capacity, Scalar T, std::size_t D>
void HomogeneousDerivatives(const BasicBSplineCurve<T, D>& curve, const double t,
                            const std::size_t known, std::array<BasicVecN<T, D>, Capacity>& aw,
                            std::array<T, Capacity>& w)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    const std::size_t span = FindSpan(curve.Knots, curve.Degree, curve.Poles.size(), t);
    std::array<double, Capacity * Capacity> table; // filled below
    BasisFunctionsOfAllDegrees(curve.Knots, curve.Degree, span, t,
                               std::span{table}.first((p + 1) * (p + 1)));
    // The active window of homogeneous poles, differenced in place.
    std::array<BasicVecN<T, D>, Capacity> qv;
    std::array<T, Capacity> qw;
    for (std::size_t i = 0; i <= p; ++i) {
        const T weight = WeightOf(curve, span - p + i);
        qv[i] = weight * ToVec(curve.Poles[span - p + i]);
        qw[i] = weight;
    }
    for (std::size_t k = 0; k < known; ++k) {
        const std::size_t degree = p - k;
        aw[k] = {};
        w[k] = T{};
        for (std::size_t i = 0; i <= degree; ++i) {
            const T basis{table[degree * (p + 1) + i]};
            aw[k] += basis * qv[i];
            w[k] = w[k] + basis * qw[i];
        }
        for (std::size_t i = 0; k + 1 < known && i < degree; ++i) {
            const double scale = static_cast<double>(degree) /
                                 (curve.Knots[span + i + 1] - curve.Knots[span - p + i + k + 1]);
            qv[i] = T{scale} * (qv[i + 1] - qv[i]);
            qw[i] = T{scale} * (qw[i + 1] - qw[i]);
        }
    }
}

/// Derivatives with work arrays of the given capacity, Capacity > degree.
template <std::size_t Capacity, Scalar T, std::size_t D>
void DerivativesWithCapacity(const BasicBSplineCurve<T, D>& curve, const double t,
                             const std::size_t order, std::span<BasicVecN<T, D>> out)
{
    const std::size_t known = std::min(order, static_cast<std::size_t>(curve.Degree)) + 1;
    std::array<BasicVecN<T, D>, Capacity> aw;
    std::array<T, Capacity> w;
    HomogeneousDerivatives<Capacity>(curve, t, known, aw, w);
    if (!IsRational(curve)) {
        for (std::size_t k = 0; k <= order; ++k) {
            out[k] = k < known ? aw[k] : BasicVecN<T, D>{};
        }
        return;
    }
    static_cast<void>(RationalDerivative<T, D>(std::span{aw}.first(known),
                                               std::span{w}.first(known), static_cast<int>(order),
                                               out.first(order + 1)));
}

} // namespace Stereon::Detail

export namespace Stereon {

template <Scalar T, std::size_t D>
[[nodiscard]] BasicPointN<T, D> Eval(const BasicBSplineCurve<T, D>& curve, const double parameter)
{
    const double t = Detail::IntoDomain(curve, parameter);
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

/// The point and its derivatives up to order n, into out[0], ..., out[n]:
/// out[0] is the position vector of C(t), out[k] = C^(k)(t). Cheaper than
/// n + 1 separate calls, because every order comes from one table of basis
/// functions (see Detail::HomogeneousDerivatives). `out` holds n + 1
/// entries.
template <Scalar T, std::size_t D>
void Derivatives(const BasicBSplineCurve<T, D>& curve, const double parameter, const int n,
                 std::span<BasicVecN<T, D>> out)
{
    const double t = Detail::IntoDomain(curve, parameter);
    Detail::CheckEvaluation(curve, t, "Derivatives(BSplineCurve)");
    Detail::Require(n >= 0 && out.size() > static_cast<std::size_t>(n), "Derivatives(BSplineCurve)",
                    "an order below 0 or too small an output");
    const auto order = static_cast<std::size_t>(n);
    if (static_cast<std::size_t>(curve.Degree) < Detail::SmallDegreeCapacity) {
        Detail::DerivativesWithCapacity<Detail::SmallDegreeCapacity>(curve, t, order, out);
    } else {
        Detail::DerivativesWithCapacity<Detail::BSplineBasisCount>(curve, t, order, out);
    }
}

/// The n-th derivative, n >= 1. A polynomial curve's derivatives above its
/// degree are zero; a rational curve's are not (algorithm A4.2).
template <Scalar T, std::size_t D>
[[nodiscard]] BasicVecN<T, D> Deriv(const BasicBSplineCurve<T, D>& curve, const double parameter,
                                    const int n)
{
    const double t = Detail::IntoDomain(curve, parameter);
    Detail::CheckEvaluation(curve, t, "Deriv(BSplineCurve)");
    Detail::Require(n >= 1, "Deriv(BSplineCurve)", "a derivative order below 1");
    if (!IsRational(curve) && n > curve.Degree) {
        return {};
    }
    const auto order = static_cast<std::size_t>(n);
    if (order < Detail::SmallDegreeCapacity) {
        std::array<BasicVecN<T, D>, Detail::SmallDegreeCapacity> all;
        Derivatives(curve, t, n, std::span{all}.first(order + 1));
        return all[order];
    }
    std::vector<BasicVecN<T, D>> all(order + 1); // high orders: rare
    Derivatives(curve, t, n, std::span{all});
    return all[order];
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

} // namespace Stereon

export namespace Stereon::Detail {

struct CurveProjectorData;

} // namespace Stereon::Detail

export namespace Stereon {

/// Projects many points onto one B-spline curve. Construction does the work
/// that does not depend on the point: the Bezier segments, their pole boxes,
/// and the point-independent parts of each segment's stationarity
/// polynomial, of which the polynomial for a point is a linear combination.
/// A query then visits the segments in order of their distance bound and
/// stops at the first that cannot hold a closer point. The answers are those
/// of Project(curve, p), which builds a projector for one point. A planar
/// curve is projected in the plane z = 0. Copies share the prepared data.
class CurveProjector {
public:
    explicit CurveProjector(const BSplineCurve& curve);
    explicit CurveProjector(const BSplineCurve2d& curve);

    /// The parameter of a point of the curve closest to p.
    [[nodiscard]] double Project(const Point3& p) const;
    [[nodiscard]] double Project(const Point2& p) const;

    /// ProjectAll reuses one query's work arrays across the points of a
    /// chunk, which takes the prepared data.
    friend std::vector<double> ProjectAll(const CurveProjector& projector,
                                          std::span<const Point3> points, const Context& context);

private:
    std::shared_ptr<const Detail::CurveProjectorData> m_Data;
};

/// The parameter of a closest point of the curve to each point, in the
/// points' order: what Project gives one by one, computed on the Context's
/// scheduler with its thread count (ADR-0007). Each point's answer goes to
/// its own slot, so the result is the same, bit for bit, for the serial
/// scheduler and for every thread count.
[[nodiscard]] std::vector<double>
ProjectAll(const CurveProjector& projector, std::span<const Point3> points, const Context& context);

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
/// of t stays at most the degree. A curve in the periodic form takes any
/// finite t, which stands for all its images a period apart: the knot is
/// inserted at every one of them, so the result is in the form again, with
/// `times` more distinct poles and the same period.
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
/// A curve in the periodic form is split as any other, on its domain: to
/// open it at t instead, take ToBSpline over [t, t + Period(curve)].
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
