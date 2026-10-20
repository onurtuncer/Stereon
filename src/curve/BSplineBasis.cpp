// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md): every module
// implementation unit of a geometry library includes it, because one that
// calls core's inline functions with contracts and parses none itself
// leaves __tu_has_violation undefined, which the macOS linker rejects.
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;

namespace Stereon {

namespace {

constexpr std::size_t MaxBasis = static_cast<std::size_t>(BSplineMaxDegree) + 1;

/// Work-array capacity for degrees below 8, which covers the curves and
/// surfaces met in practice.
constexpr std::size_t SmallCapacity = 8;

/// "<where>: <problem>", concatenated (see EvaluatorChecks.cpp).
[[nodiscard]] std::unexpected<KernelError> Invalid(const std::string_view where,
                                                   const std::string_view problem)
{
    std::string message{where};
    message += ": ";
    message += problem;
    return Fail(ErrorCode::InvalidInput, std::move(message));
}

/// The span index for a t at or outside an end of the domain: the first or
/// the last non-empty span.
[[nodiscard]] std::size_t EndSpan(const std::span<const double> knots, const std::size_t p,
                                  const std::size_t n, const bool atStart) noexcept
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

/// The table of algorithm A2.3 for one evaluation: entry (j, r) holds the
/// basis function N_{span-j+r,j} on and above the diagonal and a knot
/// difference below it. Stored compactly with stride p + 1, so a low-degree
/// evaluation touches a few cache lines, and not zeroed: every entry read
/// is written first.
template <std::size_t Capacity>
class BasisTable {
public:
    explicit BasisTable(const std::size_t stride) noexcept
        : m_Stride(stride)
    {}
    [[nodiscard]] double& operator()(const std::size_t i, const std::size_t j) noexcept
    {
        return m_Entries[i * m_Stride + j];
    }
    [[nodiscard]] double operator()(const std::size_t i, const std::size_t j) const noexcept
    {
        return m_Entries[i * m_Stride + j];
    }

private:
    std::size_t m_Stride;
    std::array<double, Capacity * Capacity> m_Entries; // NOLINT: written before read
};

template <std::size_t Capacity>
void FillBasisTable(const std::span<const double> knots, const std::size_t p,
                    const std::size_t span, const double t, BasisTable<Capacity>& ndu) noexcept
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
void BasisFunctionDerivatives(const BasisTable<Capacity>& ndu, const std::size_t p,
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

} // namespace

std::size_t FindSpan(const std::span<const double> knots, const int degree,
                     const std::size_t poleCount, const double t) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    if (t >= knots[n + 1]) {
        return EndSpan(knots, p, n, false);
    }
    if (t <= knots[p]) {
        return EndSpan(knots, p, n, true);
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

void BasisFunctions(const std::span<const double> knots, const int degree, const std::size_t span,
                    const double t, const std::span<double> values) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    // Written before read, so not zeroed: these run in every evaluation.
    std::array<double, MaxBasis> left;
    std::array<double, MaxBasis> right;
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

void BasisFunctionsOfAllDegrees(const std::span<const double> knots, const int degree,
                                const std::size_t span, const double t,
                                const std::span<double> table) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t stride = p + 1;
    // Written before read, so not zeroed: these run in every evaluation.
    std::array<double, MaxBasis> left;
    std::array<double, MaxBasis> right;
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

void BasisDerivatives(const std::span<const double> knots, const int degree, const std::size_t span,
                      const double t, const int order, const std::span<double> values) noexcept
{
    const auto p = static_cast<std::size_t>(degree);
    const auto n = static_cast<std::size_t>(order);
    const std::size_t computed = std::min(n, p);
    // Small tables for the degrees met in practice: a table sized for degree
    // 25 is 5 KB of stack, which costs more than the arithmetic.
    if (p < SmallCapacity) {
        BasisTable<SmallCapacity> ndu(p + 1);
        FillBasisTable(knots, p, span, t, ndu);
        for (std::size_t r = 0; r <= p; ++r) {
            values[r] = ndu(r, p);
            BasisFunctionDerivatives(ndu, p, r, computed, values);
        }
    } else {
        BasisTable<MaxBasis> ndu(p + 1);
        FillBasisTable(knots, p, span, t, ndu);
        for (std::size_t r = 0; r <= p; ++r) {
            values[r] = ndu(r, p);
            BasisFunctionDerivatives(ndu, p, r, computed, values);
        }
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

std::size_t KnotMultiplicity(const std::span<const double> knots, const double value) noexcept
{
    return static_cast<std::size_t>(std::ranges::count(knots, value));
}

Result<std::vector<double>> ExpandKnots(const std::span<const double> values,
                                        const std::span<const int> multiplicities)
{
    constexpr std::string_view where = "ExpandKnots";
    if (values.size() != multiplicities.size()) {
        return Invalid(where, "values and multiplicities must have the same length");
    }
    std::vector<double> knots;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i]) || (i > 0 && !(values[i] > values[i - 1]))) {
            return Invalid(where, "knot values must be finite and strictly increasing");
        }
        if (multiplicities[i] < 1) {
            return Invalid(where, "multiplicities must be positive");
        }
        knots.insert(knots.end(), static_cast<std::size_t>(multiplicities[i]), values[i]);
    }
    return knots;
}

} // namespace Stereon

namespace Stereon::Detail {

namespace {

[[nodiscard]] bool IsNonDecreasingAndFinite(const std::span<const double> knots)
{
    for (std::size_t i = 0; i < knots.size(); ++i) {
        if (!std::isfinite(knots[i]) || (i > 0 && knots[i] < knots[i - 1])) {
            return false;
        }
    }
    return true;
}

/// The largest multiplicity of a value strictly inside [low, high], and of
/// any value at all.
struct Multiplicities {
    std::size_t Inside = 0;
    std::size_t Any = 0;
};

[[nodiscard]] Multiplicities MultiplicitiesOf(const std::span<const double> knots, const double low,
                                              const double high)
{
    Multiplicities result;
    std::size_t run = 0;
    for (std::size_t i = 0; i < knots.size(); ++i) {
        run = (i > 0 && knots[i] == knots[i - 1]) ? run + 1 : 1;
        result.Any = std::max(result.Any, run);
        if (knots[i] > low && knots[i] < high) {
            result.Inside = std::max(result.Inside, run);
        }
    }
    return result;
}

} // namespace

Result<void> CheckKnotVector(const int degree, const std::span<const double> knots,
                             const std::size_t poleCount, const std::string_view where)
{
    if (degree < 1 || degree > BSplineMaxDegree) {
        return Invalid(where, "the degree must be between 1 and BSplineMaxDegree (25)");
    }
    const auto p = static_cast<std::size_t>(degree);
    if (poleCount < p + 1) {
        return Invalid(where, "a B-spline needs at least degree + 1 poles");
    }
    if (knots.size() != poleCount + p + 1) {
        return Invalid(where, "the knot vector must have poles + degree + 1 knots");
    }
    if (!IsNonDecreasingAndFinite(knots)) {
        return Invalid(where, "knots must be finite and non-decreasing");
    }
    // The basis functions divide by differences of knots, which must not
    // overflow: [-1e308, 1e308] is finite at both ends but 2e308 wide.
    if (!std::isfinite(knots.back() - knots.front())) {
        return Invalid(where, "the knot vector must span a finite width");
    }
    const double low = knots[p];
    const double high = knots[poleCount];
    if (!(low < high)) {
        return Invalid(where, "the parameter domain [u_p, u_n+1] must not be empty");
    }
    const Multiplicities m = MultiplicitiesOf(knots, low, high);
    if (m.Any > p + 1) {
        return Invalid(where, "no knot may be repeated more than degree + 1 times");
    }
    if (m.Inside > p) {
        return Invalid(where, "a knot inside the domain may be repeated at most degree times");
    }
    return {};
}

Result<void> CheckWeights(const std::span<const double> weights, const std::size_t poleCount,
                          const std::string_view where)
{
    if (weights.empty()) {
        return {};
    }
    if (weights.size() != poleCount) {
        return Invalid(where, "a rational B-spline needs one weight per pole");
    }
    for (const double w : weights) {
        if (!std::isfinite(w) || !(w > 0.0)) {
            return Invalid(where, "weights must be finite and positive");
        }
        if (w < MinWeight || w > MaxWeight) {
            return Invalid(where, "weights must lie in [MinWeight, MaxWeight], [2^-256, 2^256]");
        }
    }
    return {};
}

} // namespace Stereon::Detail
