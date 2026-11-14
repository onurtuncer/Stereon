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

/// "<where>: <problem>", concatenated (see EvaluatorChecks.cpp).
[[nodiscard]] std::unexpected<KernelError> Invalid(const std::string_view where,
                                                   const std::string_view problem)
{
    std::string message{where};
    message += ": ";
    message += problem;
    return Fail(ErrorCode::InvalidInput, std::move(message));
}

} // namespace

namespace Detail {

void BasisFunctionsDispatched(const std::span<const double> knots, const std::size_t degree,
                              const std::size_t span, const double t, double* const values) noexcept
{
    // A kernel per degree below SmallBasisCapacity; higher degrees share the
    // run-time one.
    DispatchFixed<SmallBasisCapacity>(degree, [&]<std::size_t Fixed>() {
        BasisFunctionsAt<Fixed>(knots, degree, span, t, values);
    });
}

void BasisDerivativesAtRunTime(const std::span<const double> knots, const int degree,
                               const std::size_t span, const double t, const int order,
                               const std::span<double> values) noexcept
{
    BasisDerivativesKernel(knots, degree, span, t, order, values);
}

} // namespace Detail

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
    // Evaluation divides by knot differences, and the fast evaluator caches
    // inverse span widths. A positive subnormal interval can overflow both
    // even though its knots pass the finite/order checks above.
    for (std::size_t i = 1; i < knots.size(); ++i) {
        const double width = knots[i] - knots[i - 1];
        if (width > 0.0 && !std::isfinite(1.0 / width)) {
            return Invalid(where, "non-empty knot intervals must have finite reciprocals");
        }
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
