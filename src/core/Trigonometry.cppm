// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cmath>

export module Stereon.Core:Trigonometry;

import :Scalar;

/// Sin, Cos, Atan2, Sqrt and Hypot for geometry modules above core, and
/// QuarterTurned for the derivatives of every angle-parametrised type.
///
/// A geometry module that includes <cmath> in its own global module fragment
/// while importing another module that also did can produce a module file
/// that GCC 16.2 cannot read back on glibc ("failed to load binding
/// '::sin'"; docs/toolchain/gcc16-feature-probe.md). Modules above core
/// therefore call these instead of std::sin and friends. The double
/// overloads are defined in Trigonometry.cpp, so no interface outside core
/// refers to the C library's functions; the templates serve active scalar
/// types (ADR-0021), whose functions are found by argument-dependent lookup.

export namespace Stereon {

[[nodiscard]] double Sin(double x) noexcept;
[[nodiscard]] double Cos(double x) noexcept;
/// The angle of (x, y) in (-pi, pi].
[[nodiscard]] double Atan2(double y, double x) noexcept;
/// The square root; x must not be negative.
[[nodiscard]] double Sqrt(double x) noexcept;
/// sqrt(x^2 + y^2) without overflow or underflow in the squares.
[[nodiscard]] double Hypot(double x, double y) noexcept;

template <Scalar T>
[[nodiscard]] T Sin(const T& x)
{
    using std::sin;
    return sin(x);
}

template <Scalar T>
[[nodiscard]] T Cos(const T& x)
{
    using std::cos;
    return cos(x);
}

/// A cosine and sine pair, (cos a, sin a) for some angle a.
template <Scalar T>
struct BasicCosSin {
    T Cos;
    T Sin;

    constexpr bool operator==(const BasicCosSin&) const = default;
};

/// (cos(t + n pi/2), sin(t + n pi/2)) from c = cos t and s = sin t, for
/// n >= 0. Each quarter turn maps (c, s) to (-s, c), so the result is
/// selected by n mod 4 with no rounded multiple of pi/2 in an argument.
/// The n-th derivative of cos t is cos(t + n pi/2), and of sin t is
/// sin(t + n pi/2), which is what the angle-parametrised curves and surfaces
/// use it for.
template <Scalar T>
[[nodiscard]] constexpr BasicCosSin<T> QuarterTurned(const T& c, const T& s, const int n)
{
    switch (n % 4) {
    case 1:
        return {-s, c};
    case 2:
        return {-c, -s};
    case 3:
        return {s, -c};
    default:
        return {c, s};
    }
}

} // namespace Stereon
