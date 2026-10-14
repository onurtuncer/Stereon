// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// A forward-mode dual number: Value is the primal value, Deriv its derivative
/// with respect to one parameter. It stands in for an automatic-differentiation
/// scalar in tests of the scalar-generic math (ADR-0021); the kernel itself
/// selects no AD backend.

#pragma once

#include <cmath>

namespace StereonTest {

struct Dual {
    double Value{};
    double Deriv{};

    constexpr Dual() = default;
    // Implicit, so that constants written T{2} in generic code are duals with
    // zero derivative.
    constexpr Dual(const double value)
        : Value(value)
    {} // NOLINT(google-explicit-constructor)
    constexpr Dual(const double value, const double deriv)
        : Value(value)
        , Deriv(deriv)
    {}

    constexpr bool operator==(const Dual&) const = default;

    // Ordering compares primal values, as a discrete decision on an active
    // scalar would.
    friend constexpr bool operator<(const Dual& a, const Dual& b) { return a.Value < b.Value; }
    friend constexpr bool operator>(const Dual& a, const Dual& b) { return a.Value > b.Value; }

    friend constexpr Dual operator+(const Dual& a, const Dual& b)
    {
        return {a.Value + b.Value, a.Deriv + b.Deriv};
    }
    friend constexpr Dual operator-(const Dual& a, const Dual& b)
    {
        return {a.Value - b.Value, a.Deriv - b.Deriv};
    }
    friend constexpr Dual operator-(const Dual& a) { return {-a.Value, -a.Deriv}; }
    friend constexpr Dual operator*(const Dual& a, const Dual& b)
    {
        return {a.Value * b.Value, a.Deriv * b.Value + a.Value * b.Deriv};
    }
    friend constexpr Dual operator/(const Dual& a, const Dual& b)
    {
        return {a.Value / b.Value, (a.Deriv * b.Value - a.Value * b.Deriv) / (b.Value * b.Value)};
    }

    // Found by ADL from generic code that writes `using std::sqrt; sqrt(x);`.
    friend Dual sqrt(const Dual& a)
    {
        const double root = std::sqrt(a.Value);
        return {root, a.Deriv / (2.0 * root)};
    }
    friend Dual sin(const Dual& a) { return {std::sin(a.Value), std::cos(a.Value) * a.Deriv}; }
    friend Dual cos(const Dual& a) { return {std::cos(a.Value), -std::sin(a.Value) * a.Deriv}; }
    friend Dual sinh(const Dual& a) { return {std::sinh(a.Value), std::cosh(a.Value) * a.Deriv}; }
    friend Dual cosh(const Dual& a) { return {std::cosh(a.Value), std::sinh(a.Value) * a.Deriv}; }
};

/// The parameter being differentiated with respect to: value p, derivative 1.
[[nodiscard]] constexpr Dual Variable(const double p)
{
    return {p, 1.0};
}

} // namespace StereonTest
