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

/// Sin, Cos and Atan2 for geometry modules above core.
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

} // namespace Stereon
