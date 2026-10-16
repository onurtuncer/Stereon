// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>

module Stereon.Curve;

import :Properties;

namespace Stereon::Detail {

namespace {

constexpr std::size_t GaussPoints = 16;

/// Newton steps on the Legendre polynomial; from the Chebyshev-like initial
/// guess the iteration converges in four or five, the bound is generous.
constexpr int LegendreNewtonSteps = 100;

/// P_n(x) and P_n'(x) by the three-term recurrence.
std::pair<double, double> Legendre(const double x)
{
    double previous = 1.0;
    double current = x;
    for (std::size_t k = 2; k <= GaussPoints; ++k) {
        const double kd = static_cast<double>(k);
        const double next = ((2.0 * kd - 1.0) * x * current - (kd - 1.0) * previous) / kd;
        previous = current;
        current = next;
    }
    const double n = static_cast<double>(GaussPoints);
    return {current, n * (x * current - previous) / (x * x - 1.0)};
}

std::array<std::pair<double, double>, GaussPoints> ComputeGaussLegendre()
{
    std::array<std::pair<double, double>, GaussPoints> rule{};
    const double n = static_cast<double>(GaussPoints);
    for (std::size_t i = 0; i < GaussPoints / 2; ++i) {
        // The i-th largest root lies close to cos(pi (i + 3/4) / (n + 1/2)).
        double x = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) / (n + 0.5));
        for (int step = 0; step < LegendreNewtonSteps; ++step) {
            const auto [p, dp] = Legendre(x);
            const double next = x - p / dp;
            if (next == x) {
                break;
            }
            x = next;
        }
        const double dp = Legendre(x).second;
        const double weight = 2.0 / ((1.0 - x * x) * dp * dp);
        rule[i] = {-x, weight};
        rule[GaussPoints - 1 - i] = {x, weight};
    }
    return rule;
}

} // namespace

const std::array<std::pair<double, double>, 16>& GaussLegendre16() noexcept
{
    static const std::array<std::pair<double, double>, GaussPoints> rule = ComputeGaussLegendre();
    return rule;
}

} // namespace Stereon::Detail
