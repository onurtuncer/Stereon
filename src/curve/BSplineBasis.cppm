// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <string_view>
#include <utility>
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
/// FindSpan and the basis functions are constexpr, so that small fixtures
/// can be checked by static_assert (ADR-0042, section 5); their kernels are
/// the ones in Stereon::Detail below, which use the standard library only.
/// Everything else is defined in BSplineBasis.cpp: non-inline functions keep
/// the module interface free of instantiations of core's templates (see
/// docs/toolchain/gcc16-feature-probe.md), which is also why the per-degree
/// dispatch of BasisFunctions stays there.

export namespace Stereon {

/// The highest degree a B-spline may have. 25 is OCCT's
/// Geom_BSplineCurve::MaxDegree(), so geometry exchanged with OCCT always
/// fits; the bound lets the evaluators keep their work arrays on the stack.
inline constexpr int BSplineMaxDegree = 25;

} // namespace Stereon

namespace Stereon::Detail {

/// The most basis functions that are non-zero on a span.
inline constexpr std::size_t MaxBasisCount = static_cast<std::size_t>(BSplineMaxDegree) + 1;

/// Work-array capacity for degrees below 8, which covers the curves and
/// surfaces met in practice; those degrees also get basis kernels of their
/// own (BasisFunctions).
inline constexpr std::size_t SmallBasisCapacity = 8;

/// The span index for a t at or outside an end of the domain: the first or
/// the last non-empty span.
[[nodiscard]] constexpr std::size_t EndSpan(const std::span<const double> knots,
                                            const std::size_t p, const std::size_t n,
                                            const bool atStart) noexcept
{
    if (atStart) {
        std::size_t i = p;
        while (knots[i] == knots[i + 1]) {
            ++i;
        }
        return i;
    }
    std::size_t i = n;
    while (knots[i] == knots[i + 1]) {
        --i;
    }
    return i;
}

/// A2.2 at degree Fixed, when it is not 0, so that the compiler unrolls
/// both loops and keeps the work arrays in registers; otherwise at the
/// run-time degree. The operations are the same in either case, so the
/// results are too.
template <std::size_t Fixed>
constexpr void BasisFunctionsAt(const std::span<const double> knots, const std::size_t degree,
                                const std::size_t span, const double t, double* values) noexcept
{
    const std::size_t p = Fixed != 0 ? Fixed : degree;
    constexpr std::size_t Size = Fixed != 0 ? Fixed + 1 : MaxBasisCount;
    // Written before read, so not zeroed: these run in every evaluation.
    std::array<double, Size> left;
    std::array<double, Size> right;
    values[0] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            const double temp = values[r] / (right[r + 1] + left[j - r]);
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
    }
}

/// BasisFunctions at run time: a kernel per degree below SmallBasisCapacity,
/// the run-time one above. Defined in BSplineBasis.cpp.
void BasisFunctionsDispatched(std::span<const double> knots, std::size_t degree, std::size_t span,
                              double t, double* values) noexcept;

/// The table of algorithm A2.3 for one evaluation: entry (j, r) holds the
/// basis function N_{span-j+r,j} on and above the diagonal and a knot
/// difference below it. Stored compactly with stride p + 1, so a low-degree
/// evaluation touches a few cache lines, and not zeroed: every entry read
/// is written first.
template <std::size_t Capacity>
class BasisTable {
public:
    constexpr explicit BasisTable(const std::size_t stride) noexcept
        : m_Stride(stride)
    {}
    [[nodiscard]] constexpr double& operator()(const std::size_t i, const std::size_t j) noexcept
    {
        return m_Entries[i * m_Stride + j];
    }
    [[nodiscard]] constexpr double operator()(const std::size_t i,
                                              const std::size_t j) const noexcept
    {
        return m_Entries[i * m_Stride + j];
    }

private:
    std::size_t m_Stride;
    std::array<double, Capacity * Capacity> m_Entries; // NOLINT: written before read
};

template <std::size_t Capacity>
constexpr void FillBasisTable(const std::span<const double> knots, const std::size_t p,
                              const std::size_t span, const double t,
                              BasisTable<Capacity>& ndu) noexcept
{
    // Written before read, so not zeroed: these run in every evaluation.
    std::array<double, Capacity> left;
    std::array<double, Capacity> right;
    ndu(0, 0) = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            ndu(j, r) = right[r + 1] + left[j - r]; // lower triangle: knot differences
            const double temp = ndu(r, j - 1) / ndu(j, r);
            ndu(r, j) = saved + right[r + 1] * temp; // upper triangle: basis functions
            saved = left[j - r] * temp;
        }
        ndu(j, j) = saved;
    }
}

/// The derivatives of orders 1 to n (n <= p) of basis function r of the span,
/// unscaled, into values[k * (p + 1) + r]: the inner loops of A2.3, which
/// build the coefficients a_{k,j} of each order from those of the order
/// below, alternating between two rows.
template <std::size_t Capacity>
constexpr void BasisFunctionDerivatives(const BasisTable<Capacity>& ndu, const std::size_t p,
                                        const std::size_t r, const std::size_t n,
                                        const std::span<double> values) noexcept
{
    std::array<std::array<double, Capacity>, 2> a; // written before read
    std::size_t s1 = 0;
    std::size_t s2 = 1;
    a[0][0] = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        double d = 0.0;
        const auto rk = static_cast<long long>(r) - static_cast<long long>(k);
        const std::size_t pk = p - k;
        if (r >= k) {
            a[s2][0] = a[s1][0] / ndu(pk + 1, static_cast<std::size_t>(rk));
            d = a[s2][0] * ndu(static_cast<std::size_t>(rk), pk);
        }
        const std::size_t j1 = rk >= -1 ? 1 : static_cast<std::size_t>(-rk);
        const std::size_t j2 = r <= pk + 1 ? k - 1 : p - r;
        for (std::size_t j = j1; j <= j2; ++j) {
            const auto row = static_cast<std::size_t>(rk + static_cast<long long>(j));
            a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu(pk + 1, row);
            d += a[s2][j] * ndu(row, pk);
        }
        if (r <= pk) {
            a[s2][k] = -a[s1][k - 1] / ndu(pk + 1, r);
            d += a[s2][k] * ndu(r, pk);
        }
        values[k * (p + 1) + r] = d;
        std::swap(s1, s2);
    }
}

/// A2.3 with a table of the given capacity: the values and the unscaled
/// derivatives up to order n (n <= p) of every basis function of the span.
template <std::size_t Capacity>
constexpr void UnscaledBasisDerivatives(const std::span<const double> knots, const std::size_t p,
                                        const std::size_t span, const double t, const std::size_t n,
                                        const std::span<double> values) noexcept
{
    BasisTable<Capacity> ndu(p + 1);
    FillBasisTable(knots, p, span, t, ndu);
    for (std::size_t r = 0; r <= p; ++r) {
        values[r] = ndu(r, p);
        BasisFunctionDerivatives(ndu, p, r, n, values);
    }
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The index i of the knot span [u_i, u_{i+1}) that contains t, with
/// degree <= i <= poleCount - 1 and u_i < u_{i+1} (algorithm A2.1). A t
/// at or beyond the end of the domain gets the last non-empty span, and one
/// at or before its start the first, so the end points of the domain belong
/// to the spans inside it. The knot vector must be valid for the degree and
/// pole count (CheckKnotVector).
[[nodiscard]] constexpr std::size_t FindSpan(const std::span<const double> knots, const int degree,
                                             const std::size_t poleCount, const double t) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    if (t >= knots[n + 1]) {
        return Detail::EndSpan(knots, p, n, false);
    }
    if (t <= knots[p]) {
        return Detail::EndSpan(knots, p, n, true);
    }
    // Binary search keeping knots[low] <= t < knots[high] (A2.1).
    std::size_t low = p;
    std::size_t high = n + 1;
    std::size_t mid = (low + high) / 2;
    while (t < knots[mid] || t >= knots[mid + 1]) {
        if (t < knots[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return mid;
}

/// The degree + 1 basis functions N_{span-degree,degree}(t), ...,
/// N_{span,degree}(t) that are non-zero on the span, into `values`
/// (algorithm A2.2). They are non-negative and sum to 1.
///
/// At run time a kernel per degree below 8 does the work, and in a constant
/// expression the run-time-degree kernel does: the same operations in the
/// same order, so the same values.
constexpr void BasisFunctions(const std::span<const double> knots, const int degree,
                              const std::size_t span, const double t,
                              const std::span<double> values) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    if consteval {
        Detail::BasisFunctionsAt<0>(knots, p, span, t, values.data());
    } else {
        Detail::BasisFunctionsDispatched(knots, p, span, t, values.data());
    }
}

/// The basis functions of every degree j = 0, ..., degree that are non-zero
/// on the span: table[j * (degree + 1) + r] = N_{span-j+r,j}(t) for r <= j.
/// They are the intermediate levels of BasisFunctions' triangle, kept: the
/// k-th derivative of a B-spline of degree p is a B-spline of degree p - k
/// on the same knots, so the levels below p are what derivatives need
/// (see Deriv(BSplineCurve)). `table` holds (degree + 1)^2 entries.
constexpr void BasisFunctionsOfAllDegrees(const std::span<const double> knots, const int degree,
                                          const std::size_t span, const double t,
                                          const std::span<double> table) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t stride = p + 1;
    // Written before read, so not zeroed: these run in every evaluation.
    std::array<double, Detail::MaxBasisCount> left;
    std::array<double, Detail::MaxBasisCount> right;
    table[0] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            const double temp = table[(j - 1) * stride + r] / (right[r + 1] + left[j - r]);
            table[j * stride + r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        table[j * stride + j] = saved;
    }
}

/// The basis functions of the span and their derivatives up to `order`:
/// values[k * (degree + 1) + j] is the k-th derivative of
/// N_{span-degree+j,degree} at t (algorithm A2.3). Derivatives above the
/// degree are zero. `values` holds (order + 1) * (degree + 1) entries.
constexpr void BasisDerivatives(const std::span<const double> knots, const int degree,
                                const std::size_t span, const double t, const int order,
                                const std::span<double> values) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const auto n = static_cast<std::size_t>(order);
    const std::size_t computed = std::min(n, p);
    // Small tables for the degrees met in practice: a table sized for degree
    // 25 is 5 KB of stack, which costs more than the arithmetic.
    if (p < Detail::SmallBasisCapacity) {
        Detail::UnscaledBasisDerivatives<Detail::SmallBasisCapacity>(knots, p, span, t, computed,
                                                                     values);
    } else {
        Detail::UnscaledBasisDerivatives<Detail::MaxBasisCount>(knots, p, span, t, computed,
                                                                values);
    }
    // Scale order k by p! / (p - k)!; orders above the degree vanish.
    double factor = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        factor *= static_cast<double>(p + 1 - std::min(k, p + 1));
        for (std::size_t r = 0; r <= p; ++r) {
            values[k * (p + 1) + r] = k > p ? 0.0 : values[k * (p + 1) + r] * factor;
        }
    }
}

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
/// non-decreasing, every non-empty knot interval having a finite reciprocal,
/// a non-empty domain [u_p, u_n+1], no knot value repeated
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
