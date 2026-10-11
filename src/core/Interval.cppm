// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>

export module Stereon.Core:Interval;

import :Scalar;

/// Closed intervals of doubles with outward rounding: the enclosure
/// arithmetic behind the predicate filters of ADR-0002. Every operation
/// returns an interval containing every real result of the operation on real
/// numbers taken from its operands.
///
/// Rounding steps each bound one representable value outward, on the bit
/// representation, instead of changing the processor rounding mode. The
/// arithmetic is therefore constexpr, has no global state and gives the same
/// bits on every platform (ADR-0007). The price is that an exact result is
/// also widened by one step.
///
/// Interval is double only. It bounds and decides, which stays discrete
/// (ADR-0021).

namespace Stereon::Detail {

inline constexpr double Infinity = std::numeric_limits<double>::infinity();

/// The next representable double above x. +inf and NaN come back unchanged;
/// both zeros step to the smallest positive subnormal.
[[nodiscard]] constexpr double NextUp(const double x) noexcept
{
    if (x != x || x == Infinity) {
        return x;
    }
    if (x == 0.0) {
        return std::numeric_limits<double>::denorm_min();
    }
    const auto bits = std::bit_cast<std::uint64_t>(x);
    return std::bit_cast<double>(x > 0.0 ? bits + 1 : bits - 1);
}

[[nodiscard]] constexpr double NextDown(const double x) noexcept
{
    return -NextUp(-x);
}

/// Multiplication with 0 * inf = 0, the interval-arithmetic convention: a zero
/// factor is an exact zero however large the other factor.
[[nodiscard]] constexpr double Multiply(const double a, const double b) noexcept
{
    return (a == 0.0 || b == 0.0) ? 0.0 : a * b;
}

} // namespace Stereon::Detail

export namespace Stereon {

class Interval {
public:
    /// The empty interval, the identity of Hull.
    constexpr Interval() = default;

    /// The closed interval [lower, upper]. Bounds may be infinite; NaN is not
    /// a bound.
    constexpr Interval(const double lower, const double upper) STEREON_PRE(lower <= upper)
        : m_Lower(lower)
        , m_Upper(upper)
    {}

    /// The degenerate interval [x, x].
    explicit constexpr Interval(const double x)
        : m_Lower(x)
        , m_Upper(x)
    {}

    [[nodiscard]] static constexpr Interval Empty() { return {}; }
    [[nodiscard]] static constexpr Interval Everything()
    {
        return {-Detail::Infinity, Detail::Infinity};
    }

    [[nodiscard]] constexpr double GetLower() const noexcept { return m_Lower; }
    [[nodiscard]] constexpr double GetUpper() const noexcept { return m_Upper; }
    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return m_Lower > m_Upper; }
    /// True when neither bound is infinite; the empty interval is bounded.
    [[nodiscard]] constexpr bool IsBounded() const noexcept
    {
        return -Detail::Infinity < m_Lower && m_Upper < Detail::Infinity;
    }

    constexpr bool operator==(const Interval&) const = default;

    friend constexpr Interval operator+(const Interval& a, const Interval& b)
    {
        if (a.IsEmpty() || b.IsEmpty()) {
            return {};
        }
        return {Detail::NextDown(a.m_Lower + b.m_Lower), Detail::NextUp(a.m_Upper + b.m_Upper)};
    }

    friend constexpr Interval operator-(const Interval& a, const Interval& b)
    {
        if (a.IsEmpty() || b.IsEmpty()) {
            return {};
        }
        return {Detail::NextDown(a.m_Lower - b.m_Upper), Detail::NextUp(a.m_Upper - b.m_Lower)};
    }

    /// Exact: negation never rounds.
    friend constexpr Interval operator-(const Interval& a)
    {
        if (a.IsEmpty()) {
            return {};
        }
        return {-a.m_Upper, -a.m_Lower};
    }

    friend constexpr Interval operator*(const Interval& a, const Interval& b)
    {
        if (a.IsEmpty() || b.IsEmpty()) {
            return {};
        }
        const double ll = Detail::Multiply(a.m_Lower, b.m_Lower);
        const double lu = Detail::Multiply(a.m_Lower, b.m_Upper);
        const double ul = Detail::Multiply(a.m_Upper, b.m_Lower);
        const double uu = Detail::Multiply(a.m_Upper, b.m_Upper);
        return {Detail::NextDown(std::min({ll, lu, ul, uu})),
                Detail::NextUp(std::max({ll, lu, ul, uu}))};
    }

    /// Division by an interval containing zero has no finite enclosure, so it
    /// is a precondition, not a result.
    friend constexpr Interval operator/(const Interval& a, const Interval& b)
        STEREON_PRE(b.m_Lower > 0.0 || b.m_Upper < 0.0)
    {
        if (a.IsEmpty() || b.IsEmpty()) {
            return {};
        }
        // 1 / [l, u] = [1 / u, 1 / l] when 0 is outside; 1 / inf is 0.
        const Interval reciprocal{Detail::NextDown(1.0 / b.m_Upper),
                                  Detail::NextUp(1.0 / b.m_Lower)};
        return a * reciprocal;
    }

private:
    double m_Lower = Detail::Infinity;
    double m_Upper = -Detail::Infinity;
};

[[nodiscard]] constexpr bool Contains(const Interval& i, const double x) noexcept
{
    return i.GetLower() <= x && x <= i.GetUpper();
}

/// Every interval contains the empty one.
[[nodiscard]] constexpr bool Contains(const Interval& outer, const Interval& inner) noexcept
{
    return outer.GetLower() <= inner.GetLower() && inner.GetUpper() <= outer.GetUpper();
}

/// Closed intervals that share only an endpoint intersect; the empty interval
/// intersects nothing. (The emptiness checks are needed: the empty interval's
/// bounds, +inf and -inf, pass the comparisons against unbounded intervals.)
[[nodiscard]] constexpr bool Intersects(const Interval& a, const Interval& b) noexcept
{
    return !a.IsEmpty() && !b.IsEmpty() && a.GetLower() <= b.GetUpper() &&
           b.GetLower() <= a.GetUpper();
}

/// The smallest interval containing both.
[[nodiscard]] constexpr Interval Hull(const Interval& a, const Interval& b)
{
    if (a.IsEmpty()) {
        return b;
    }
    if (b.IsEmpty()) {
        return a;
    }
    return {std::min(a.GetLower(), b.GetLower()), std::max(a.GetUpper(), b.GetUpper())};
}

[[nodiscard]] constexpr Interval Hull(const Interval& i, const double x)
{
    return Hull(i, Interval{x});
}

/// Empty when the intervals do not intersect.
[[nodiscard]] constexpr Interval Intersection(const Interval& a, const Interval& b)
{
    const double lower = std::max(a.GetLower(), b.GetLower());
    const double upper = std::min(a.GetUpper(), b.GetUpper());
    if (lower > upper) {
        return {};
    }
    return {lower, upper};
}

[[nodiscard]] constexpr double Width(const Interval& i) STEREON_PRE(!i.IsEmpty())
{
    return i.GetUpper() - i.GetLower();
}

[[nodiscard]] constexpr double Midpoint(const Interval& i)
    STEREON_PRE(!i.IsEmpty() && i.IsBounded())
{
    // Never overflows, unlike (lower + upper) / 2.
    return 0.5 * i.GetLower() + 0.5 * i.GetUpper();
}

/// Each bound moved outward by margin; a negative margin shrinks the interval
/// and may empty it. The empty interval stays empty.
[[nodiscard]] constexpr Interval Inflated(const Interval& i, const double margin)
{
    if (i.IsEmpty()) {
        return {};
    }
    const double lower = i.GetLower() - margin;
    const double upper = i.GetUpper() + margin;
    if (lower > upper) {
        return {};
    }
    return {lower, upper};
}

/// The square of the interval as a set, which is tighter than i * i when the
/// interval contains zero: [-1, 2] squared is [0, 4], not [-2, 4].
[[nodiscard]] constexpr Interval Squared(const Interval& i)
{
    if (i.IsEmpty()) {
        return {};
    }
    const double ll = Detail::Multiply(i.GetLower(), i.GetLower());
    const double uu = Detail::Multiply(i.GetUpper(), i.GetUpper());
    const double upper = Detail::NextUp(std::max(ll, uu));
    if (Contains(i, 0.0)) {
        return {0.0, upper};
    }
    return {Detail::NextDown(std::min(ll, uu)), upper};
}

[[nodiscard]] constexpr Interval Sqrt(const Interval& i)
    STEREON_PRE(i.IsEmpty() || i.GetLower() >= 0.0)
{
    if (i.IsEmpty()) {
        return {};
    }
    return {std::max(0.0, Detail::NextDown(std::sqrt(i.GetLower()))),
            Detail::NextUp(std::sqrt(i.GetUpper()))};
}

} // namespace Stereon

// std::format support: "[lower, upper]", or "[]" when empty.

template <>
struct std::formatter<Stereon::Interval> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::Interval& i, FormatContext& ctx) const
    {
        if (i.IsEmpty()) {
            return std::format_to(ctx.out(), "[]");
        }
        return std::format_to(ctx.out(), "[{}, {}]", i.GetLower(), i.GetUpper());
    }
};
