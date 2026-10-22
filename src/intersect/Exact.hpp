// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The exact decisions of the analytic handlers (RS-0001 G5, ADR-0002): a
// configuration decision is the sign of a polynomial in the input doubles
// and the tolerance. Each quantity is carried two ways, as an outward-rounded
// Interval for the fast path and as an Expansion for the exact one, and a
// decision reads the interval first and builds the exact polynomial only
// when the interval straddles zero (FilteredSign). Nothing here depends on
// round-off, evaluation order or thread count.
//
// Included in the purview of the handlers' implementation units, after
// `import Stereon.Robust;`, with <array> and <cstddef> in their global
// module fragment. Everything is in an anonymous namespace, so each unit
// gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

/// A scalar given exactly, with a rounded enclosure.
template <std::size_t N>
struct ExactScalar {
    Interval Hull;
    Expansion<N> Value;
};

/// A vector given exactly, with a rounded enclosure.
template <std::size_t N>
struct ExactVector {
    std::array<Interval, 3> Hull;
    std::array<Expansion<N>, 3> Value;
};

[[nodiscard]] inline ExactScalar<1> Exactly(const double x)
{
    return {Interval{x}, Expansion<1>(x)};
}

[[nodiscard]] inline ExactVector<1> Exactly(const Vec3& v)
{
    return {{Interval{v.X}, Interval{v.Y}, Interval{v.Z}},
            {Expansion<1>(v.X), Expansion<1>(v.Y), Expansion<1>(v.Z)}};
}

/// p - q.
[[nodiscard]] inline ExactVector<2> ExactOffset(const Point3& p, const Point3& q)
{
    return {{Interval{p.X} - Interval{q.X}, Interval{p.Y} - Interval{q.Y},
             Interval{p.Z} - Interval{q.Z}},
            {TwoDiff(p.X, q.X), TwoDiff(p.Y, q.Y), TwoDiff(p.Z, q.Z)}};
}

/// a - b.
[[nodiscard]] inline ExactScalar<2> ExactDifference(const double a, const double b)
{
    return {Interval{a} - Interval{b}, TwoDiff(a, b)};
}

/// a + b.
[[nodiscard]] inline ExactScalar<2> ExactSum(const double a, const double b)
{
    return {Interval{a} + Interval{b}, TwoSum(a, b)};
}

template <std::size_t N>
[[nodiscard]] inline ExactVector<N> ExactNegated(const ExactVector<N>& v)
{
    return {{-v.Hull[0], -v.Hull[1], -v.Hull[2]},
            {Negated(v.Value[0]), Negated(v.Value[1]), Negated(v.Value[2])}};
}

template <std::size_t A, std::size_t B>
[[nodiscard]] inline ExactVector<4 * A * B> ExactCross(const ExactVector<A>& a,
                                                       const ExactVector<B>& b)
{
    return {{a.Hull[1] * b.Hull[2] - a.Hull[2] * b.Hull[1],
             a.Hull[2] * b.Hull[0] - a.Hull[0] * b.Hull[2],
             a.Hull[0] * b.Hull[1] - a.Hull[1] * b.Hull[0]},
            {Difference(Product(a.Value[1], b.Value[2]), Product(a.Value[2], b.Value[1])),
             Difference(Product(a.Value[2], b.Value[0]), Product(a.Value[0], b.Value[2])),
             Difference(Product(a.Value[0], b.Value[1]), Product(a.Value[1], b.Value[0]))}};
}

template <std::size_t A, std::size_t B>
[[nodiscard]] inline ExactScalar<6 * A * B> ExactDot(const ExactVector<A>& a,
                                                     const ExactVector<B>& b)
{
    return {a.Hull[0] * b.Hull[0] + a.Hull[1] * b.Hull[1] + a.Hull[2] * b.Hull[2],
            Sum(Sum(Product(a.Value[0], b.Value[0]), Product(a.Value[1], b.Value[1])),
                Product(a.Value[2], b.Value[2]))};
}

/// The sign of s.
template <std::size_t N>
[[nodiscard]] inline Sign SignOf(const ExactScalar<N>& s)
{
    return FilteredSign(s.Hull, [&] { return s.Value.GetSign(); });
}

/// True when v is exactly the zero vector.
template <std::size_t N>
[[nodiscard]] inline bool IsExactlyZero(const ExactVector<N>& v)
{
    for (std::size_t i = 0; i < 3; ++i) {
        if (const std::optional<Sign> sign = SignOf(v.Hull[i]); sign && *sign != Sign::Zero) {
            return false;
        }
    }
    return v.Value[0].IsZero() && v.Value[1].IsZero() && v.Value[2].IsZero();
}

/// The sign of value - factor² scale: how a squared length compares with a
/// radius or tolerance squared, scaled by a squared normal's length.
template <std::size_t A, std::size_t B, std::size_t C>
[[nodiscard]] inline Sign Compare(const ExactScalar<A>& value, const ExactScalar<B>& factor,
                                  const ExactScalar<C>& scale)
{
    const Interval estimate = value.Hull - Squared(factor.Hull) * scale.Hull;
    return FilteredSign(estimate, [&] {
        return Difference(value.Value, Product(Product(factor.Value, factor.Value), scale.Value))
            .GetSign();
    });
}

/// The sign of value² - factor² scale: Compare for a value that is a signed
/// length rather than a squared one.
template <std::size_t A, std::size_t B, std::size_t C>
[[nodiscard]] inline Sign CompareSquared(const ExactScalar<A>& value, const ExactScalar<B>& factor,
                                         const ExactScalar<C>& scale)
{
    const Interval estimate = Squared(value.Hull) - Squared(factor.Hull) * scale.Hull;
    return FilteredSign(estimate, [&] {
        return Difference(Product(value.Value, value.Value),
                          Product(Product(factor.Value, factor.Value), scale.Value))
            .GetSign();
    });
}

} // namespace
} // namespace Stereon::Detail
