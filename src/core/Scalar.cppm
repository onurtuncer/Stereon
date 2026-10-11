// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <concepts>
#include <format>

export module Stereon.Core:Scalar;

export namespace Stereon {

/// A real-number-like type that the math types and the built-in geometry
/// evaluators are written for (ADR-0021): double, float and the active scalar
/// types of an automatic-differentiation backend. Integers are excluded; the
/// kernel's arithmetic is on lengths, angles and parameters.
///
/// Generic code keeps the scalar type: constants are written T{2}, never
/// converted through double, and mathematical functions are called through
/// ADL (using std::sqrt; sqrt(x);) so that an active type can supply its own.
template <typename T>
concept Scalar = std::regular<T> && !std::integral<T> && std::constructible_from<T, double> &&
                 requires(const T a, const T b) {
                     { a + b } -> std::convertible_to<T>;
                     { a - b } -> std::convertible_to<T>;
                     { a * b } -> std::convertible_to<T>;
                     { a / b } -> std::convertible_to<T>;
                     { -a } -> std::convertible_to<T>;
                     { a < b } -> std::convertible_to<bool>;
                 };

} // namespace Stereon

namespace Stereon::Detail {

/// Shared by the std::formatter specialisations of the math types, none of
/// which take a format specification of their own.
struct EmptySpecParser {
    template <typename ParseContext>
    static constexpr auto Parse(ParseContext& ctx)
    {
        const auto it = ctx.begin();
        if (it != ctx.end() && *it != '}') {
            throw std::format_error("Stereon math types take no format specification");
        }
        return it;
    }
};

} // namespace Stereon::Detail
