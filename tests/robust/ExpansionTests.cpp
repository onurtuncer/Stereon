// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

import Stereon.Robust;

#include "ExactDeterminant.hpp"
#include "RobustFormatting.hpp"

using Stereon::Expansion;
using Stereon::Sign;
using Stereon::TwoDiff;
using Stereon::TwoProduct;
using Stereon::TwoSum;

namespace {

/// An integer-valued double with up to `bits` significant bits, odd when
/// possible so that products and sums of two of them need rounding.
double RandomInteger(const int bits)
{
    const auto magnitude = *rc::gen::inRange<std::int64_t>(0, std::int64_t{1} << bits);
    const bool negative = *rc::gen::arbitrary<bool>();
    const std::int64_t value = (magnitude | 1) * (negative ? -1 : 1);
    return static_cast<double>(value);
}

#if STEREON_HAS_EXACT_ORACLE
/// The exact value of an expansion whose components are integers.
template <std::size_t N>
Stereon::Tests::Wide WideValue(const Expansion<N>& e)
{
    Stereon::Tests::Wide total = 0;
    for (const double c : e.GetComponents()) {
        RC_ASSERT(Stereon::Tests::IsIntegerValued(c));
        total += Stereon::Tests::ToWide(c);
    }
    return total;
}
#endif

} // namespace

TEST_CASE("expansion capacities follow the formula", "[robust][expansion]")
{
    STATIC_REQUIRE(std::is_same_v<decltype(TwoSum(1.0, 2.0)), Expansion<2>>);
    STATIC_REQUIRE(
        std::is_same_v<decltype(TwoProduct(1.0, 2.0) - TwoProduct(3.0, 4.0)), Expansion<4>>);
    STATIC_REQUIRE(std::is_same_v<decltype(Expansion<4>{} * 2.0), Expansion<8>>);
    STATIC_REQUIRE(std::is_same_v<decltype(Expansion<4>{} * Expansion<6>{}), Expansion<48>>);
    STATIC_REQUIRE(
        std::is_same_v<decltype(Expansion<8>{} + Expansion<8>{} + Expansion<8>{}), Expansion<24>>);
    STATIC_REQUIRE(Expansion<96>::Capacity == 96);
    STATIC_REQUIRE(Stereon::Splitter == 0x1p27 + 1.0);
    // STEREON-LINT-ALLOW(literal): the decimal Shewchuk's paper gives
    STATIC_REQUIRE(Stereon::Splitter == 134217729.0);
}

TEST_CASE("expansions are built at compile time and read by their last component",
          "[robust][expansion]")
{
    constexpr Expansion<2> six = TwoProduct(2.0, 3.0);
    STATIC_REQUIRE(six.GetLength() == 1);
    STATIC_REQUIRE(six.GetMostSignificant() == 6.0);
    STATIC_REQUIRE(six.GetEstimate() == 6.0);
    STATIC_REQUIRE(six.GetSign() == Sign::Positive);

    STATIC_REQUIRE(TwoProduct(0.0, 3.0).IsZero());
    STATIC_REQUIRE(TwoProduct(0.0, 3.0).GetSign() == Sign::Zero);
    STATIC_REQUIRE(Expansion<4>{}.IsZero());
    STATIC_REQUIRE(Expansion<4>{-2.5}.GetSign() == Sign::Negative);
    STATIC_REQUIRE(Expansion<4>{0.0}.IsZero());

    // 1 + 2^-60 is not a double, so the sum keeps its roundoff.
    constexpr Expansion<2> sum = TwoSum(1.0, 0x1p-60);
    STATIC_REQUIRE(sum.GetLength() == 2);
    STATIC_REQUIRE(sum[0] == 0x1p-60);
    STATIC_REQUIRE(sum[1] == 1.0);
    STATIC_REQUIRE(TwoDiff(1.0, 0x1p-60)[0] == -0x1p-60);
    STATIC_REQUIRE(Stereon::TwoDiffTail(1.0, 0x1p-60, 1.0 - 0x1p-60) == -0x1p-60);

    // (2^27 + 1)^2 = 2^54 + 2^28 + 1 needs 55 bits: the product carries a
    // roundoff of 1, and the exact difference from 2^27 * (2^27 + 2) is 1
    // where floating point sees 0.
    constexpr double a = 0x1p27 + 1.0;
    constexpr double c = 0x1p27;
    constexpr double d = 0x1p27 + 2.0;
    STATIC_REQUIRE(a * a - c * d == 0.0);
    STATIC_REQUIRE(TwoProduct(a, a).GetLength() == 2);
    STATIC_REQUIRE(TwoProduct(a, a)[0] == 1.0);
    STATIC_REQUIRE((TwoProduct(a, a) - TwoProduct(c, d)).GetSign() == Sign::Positive);
    STATIC_REQUIRE((TwoProduct(a, a) - TwoProduct(c, d)).GetMostSignificant() == 1.0);
    STATIC_REQUIRE((TwoProduct(c, d) - TwoProduct(a, a)).GetSign() == Sign::Negative);
    STATIC_REQUIRE((TwoProduct(a, a) - TwoProduct(a, a)).IsZero());

    // Widening keeps the components; narrowing is the identity on a short
    // expansion.
    constexpr Expansion<8> wide = sum;
    STATIC_REQUIRE(wide.GetLength() == 2);
    STATIC_REQUIRE(wide[1] == 1.0);
    STATIC_REQUIRE(Expansion<2>::Narrowed(wide)[0] == 0x1p-60);
}

TEST_CASE("splitting a double gives two exact halves", "[robust][expansion][property]")
{
    rc::prop("hi + lo == a and the halves have at most 26 bits", [] {
        const double a = RandomInteger(52);
        const Stereon::Halves halves = Stereon::Split(a);
        RC_ASSERT(halves.Hi + halves.Lo == a);
        // A 26-bit integer times 2^27 is still exact in a double.
        RC_ASSERT(std::fabs(halves.Hi) < 0x1p53);
        RC_ASSERT(std::fabs(halves.Lo) <= 0x1p26);
        RC_ASSERT(std::fabs(halves.Lo) <= std::fabs(halves.Hi) || halves.Hi == 0.0);
    });
}

#if STEREON_HAS_EXACT_ORACLE

TEST_CASE("error-free transformations are exact", "[robust][expansion][property]")
{
    using Stereon::Tests::ToWide;
    rc::prop("TwoSum, TwoDiff and TwoProduct keep every bit", [] {
        const double a = RandomInteger(52);
        const double b = RandomInteger(52) * 0x1p20; // different scale: the sum rounds
        RC_ASSERT(WideValue(TwoSum(a, b)) == ToWide(a) + ToWide(b));
        RC_ASSERT(WideValue(TwoDiff(a, b)) == ToWide(a) - ToWide(b));
        RC_ASSERT(WideValue(TwoProduct(a, b)) == ToWide(a) * ToWide(b));
        RC_ASSERT(Stereon::TwoDiffTail(a, b, a - b) ==
                  static_cast<double>(ToWide(a) - ToWide(b) - ToWide(a - b)));
    });
}

TEST_CASE("expansion arithmetic is exact", "[robust][expansion][property]")
{
    using Stereon::Tests::ToWide;
    using Stereon::Tests::Wide;

    rc::prop("sums and differences of products", [] {
        const Expansion<2> p = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Expansion<2> q = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Expansion<2> r = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Wide exact = WideValue(p) - WideValue(q) + WideValue(r);
        const Expansion<6> sum = p - q + r;
        RC_ASSERT(WideValue(sum) == exact);
        RC_ASSERT(sum.GetSign() == Stereon::Tests::SignOfWide(exact));
        for (std::size_t i = 1; i < sum.GetLength(); ++i) {
            RC_ASSERT(std::fabs(sum[i - 1]) < std::fabs(sum[i])); // increasing, nonoverlapping
            RC_ASSERT(sum[i - 1] != 0.0);
        }
    });

    rc::prop("scaling and products", [] {
        const Expansion<4> e = TwoProduct(RandomInteger(40), RandomInteger(40)) -
                               TwoProduct(RandomInteger(40), RandomInteger(40));
        const double b = RandomInteger(40);
        RC_ASSERT(WideValue(e * b) == WideValue(e) * ToWide(b));
        RC_ASSERT(WideValue(b * e) == WideValue(e) * ToWide(b));

        const Expansion<2> f = TwoProduct(RandomInteger(20), RandomInteger(20));
        RC_ASSERT(WideValue(e * f) == WideValue(e) * WideValue(f));
        RC_ASSERT(WideValue(f * e) == WideValue(e) * WideValue(f));
        RC_ASSERT(WideValue(-e) == -WideValue(e));
        RC_ASSERT((e - e).IsZero());
        RC_ASSERT((e * 0.0).IsZero());
    });

    rc::prop("the heap expansion agrees with the stack one, component for component", [] {
        const Expansion<4> e = TwoProduct(RandomInteger(40), RandomInteger(40)) -
                               TwoProduct(RandomInteger(40), RandomInteger(40));
        const Expansion<2> f = TwoProduct(RandomInteger(20), RandomInteger(20));
        const double b = RandomInteger(40);
        const Stereon::DynamicExpansion de{e};
        const Stereon::DynamicExpansion df{f};
        const auto same = [](const auto& stack, const Stereon::DynamicExpansion& heap) {
            RC_ASSERT(stack.GetLength() == heap.GetLength());
            for (std::size_t i = 0; i < heap.GetLength(); ++i) {
                RC_ASSERT(stack[i] == heap[i]);
            }
            RC_ASSERT(stack.GetSign() == heap.GetSign());
            RC_ASSERT(stack.GetEstimate() == heap.GetEstimate());
        };
        same(e + f, Sum(de, df));
        same(e - f, Difference(de, df));
        same(e * b, Scaled(de, b));
        same(e * f, Product(de, df));
        same(-e, Negated(de));
        RC_ASSERT(Difference(de, de).IsZero());
        RC_ASSERT(Scaled(de, 0.0).IsZero());
    });

    rc::prop("a heap expansion is exact at a degree beyond any stack capacity", [] {
        // Ten factors of twelve bits: the exact product has 120 bits, within
        // the oracle, and its stack capacity would be 2^10 components per
        // monomial.
        Stereon::DynamicExpansion product{1.0};
        Wide exact = 1;
        for (int i = 0; i < 10; ++i) {
            const double factor = RandomInteger(12);
            product = Scaled(product, factor);
            exact *= ToWide(factor);
        }
        const Expansion<2> other = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Stereon::DynamicExpansion difference =
            Difference(product, Stereon::DynamicExpansion{other});
        Wide total = 0;
        for (const double c : difference.GetComponents()) {
            RC_ASSERT(Stereon::Tests::IsIntegerValued(c));
            total += ToWide(c);
        }
        RC_ASSERT(total == exact - WideValue(other));
        RC_ASSERT(difference.GetSign() == Stereon::Tests::SignOfWide(total));
        for (std::size_t i = 1; i < difference.GetLength(); ++i) {
            RC_ASSERT(std::fabs(difference[i - 1]) < std::fabs(difference[i]));
        }
    });

    rc::prop("the estimate is the value to within a few roundings", [] {
        const Expansion<2> p = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Expansion<2> q = TwoProduct(RandomInteger(52), RandomInteger(52));
        const Expansion<4> d = p - q;
        const double exact = static_cast<double>(WideValue(d));
        const double estimate = d.GetEstimate();
        RC_ASSERT(std::fabs(estimate - exact) <= 4.0 * 0x1p-53 * std::fabs(exact));
        // The last component alone is only within a factor of two.
        RC_ASSERT(d.IsZero() || std::fabs(d.GetMostSignificant() - exact) <
                                    2.0 * std::fabs(d.GetMostSignificant()));
    });
}

#endif // STEREON_HAS_EXACT_ORACLE
