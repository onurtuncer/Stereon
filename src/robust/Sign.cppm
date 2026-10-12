// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdint>
#include <format>
#include <string_view>

export module Stereon.Robust:Sign;

import Stereon.Core;

/// The answer of a predicate (ADR-0002). A predicate returns the sign of a
/// determinant, never the determinant: the sign is exact, a floating-point
/// value of it would not be, and a caller that compared such a value against
/// zero would be making the decision the predicate exists to make.

export namespace Stereon {

enum class Sign : std::int8_t {
    Negative = -1,
    Zero = 0,
    Positive = 1,
};

/// The sign of a floating-point number. NaN is outside the domain of every
/// predicate; here it compares as Zero.
[[nodiscard]] constexpr Sign SignOf(const double x) noexcept
{
    if (x > 0.0) {
        return Sign::Positive;
    }
    return x < 0.0 ? Sign::Negative : Sign::Zero;
}

/// The opposite sign; Zero stays Zero.
[[nodiscard]] constexpr Sign operator-(const Sign s) noexcept
{
    return static_cast<Sign>(-static_cast<int>(s));
}

/// -1, 0 or 1.
[[nodiscard]] constexpr int ToInt(const Sign s) noexcept
{
    return static_cast<int>(s);
}

/// "Positive", "Zero" or "Negative"; "Sign" where the compiler has no
/// reflection.
[[nodiscard]] constexpr std::string_view ToString(const Sign s) noexcept
{
    const std::string_view name = EnumName(s);
    return name.empty() ? std::string_view{"Sign"} : name;
}

} // namespace Stereon

template <>
struct std::formatter<Stereon::Sign> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        const auto it = ctx.begin();
        if (it != ctx.end() && *it != '}') {
            throw std::format_error("Sign takes no format specification");
        }
        return it;
    }
    template <typename FormatContext>
    auto format(const Stereon::Sign s, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "{}", Stereon::ToString(s));
    }
};
