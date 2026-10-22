// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Exact arithmetic without a capacity: DynamicExpansion agrees with
// Expansion<N> component for component and is exact where doubles are not;
// ExactPolynomial evaluates signs exactly; SturmSequence counts the distinct
// real roots in any interval, whatever the leading coefficient's sign.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

import Stereon.Robust;

using Stereon::DynamicExpansion;
using Stereon::ExactPolynomial;
using Stereon::Expansion;
using Stereon::Sign;
using Stereon::SturmSequence;

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();

bool SameComponents(const DynamicExpansion& d, const std::span<const double> fixed)
{
    const std::span<const double> dc = d.GetComponents();
    if (dc.size() != fixed.size()) {
        return false;
    }
    for (std::size_t i = 0; i < dc.size(); ++i) {
        if (dc[i] != fixed[i]) {
            return false;
        }
    }
    return true;
}

/// (x - r1)(x - r2)...: the monic polynomial with the given roots.
ExactPolynomial WithRoots(const std::vector<double>& roots)
{
    ExactPolynomial p(std::vector<DynamicExpansion>{DynamicExpansion{1.0}});
    for (const double r : roots) {
        const ExactPolynomial factor(
            std::vector<DynamicExpansion>{DynamicExpansion{-r}, DynamicExpansion{1.0}});
        p = Product(p, factor);
    }
    return p;
}

} // namespace

TEST_CASE("dynamic expansions agree with fixed ones and are exact", "[robust][dynamic]")
{
    rc::prop("sum, product and scaling match Expansion<N>", [] {
        const double a = *rc::gen::arbitrary<double>();
        const double b = *rc::gen::arbitrary<double>();
        const double c = *rc::gen::arbitrary<double>();
        RC_PRE(std::isfinite(a) && std::isfinite(b) && std::isfinite(c));
        RC_PRE(std::abs(a) < 1e100 && std::abs(b) < 1e100 && std::abs(c) < 1e100);
        RC_PRE(std::abs(a) > 1e-100 && std::abs(b) > 1e-100 && std::abs(c) > 1e-100);
        const Expansion<2> ab = Stereon::TwoProduct(a, b);
        const DynamicExpansion dab = Product(DynamicExpansion{a}, DynamicExpansion{b});
        RC_ASSERT(SameComponents(dab, ab.GetComponents()));
        const Expansion<4> sum = Sum(ab, Stereon::TwoSum(a, c));
        RC_ASSERT(SameComponents(Sum(dab, Sum(DynamicExpansion{a}, DynamicExpansion{c})),
                                 sum.GetComponents()));
        RC_ASSERT(SameComponents(Scaled(dab, c), Scaled(ab, c).GetComponents()));
        RC_ASSERT(SameComponents(DynamicExpansion{sum}, sum.GetComponents()));
    });
    // Exact where doubles are not: (a + b) - a - b and the difference of
    // squares identity vanish exactly.
    const DynamicExpansion a{0.1};
    const DynamicExpansion b{0.2};
    REQUIRE(((a + b) - a - b).IsZero());
    REQUIRE(((a + b) * (a - b) - (a * a - b * b)).IsZero());
    REQUIRE((a + b).GetSign() == Sign::Positive);
    REQUIRE((-(a + b)).GetSign() == Sign::Negative);
    REQUIRE(DynamicExpansion{}.GetSign() == Sign::Zero);
    REQUIRE(std::abs((a + b).GetEstimate() - 0.3) < 1e-16);
}

TEST_CASE("exact polynomials evaluate their signs exactly", "[robust][dynamic]")
{
    // (x - 1)(x - 2)(x - 3)(x - 4)
    const ExactPolynomial p = WithRoots({1.0, 2.0, 3.0, 4.0});
    REQUIRE(p.GetDegree() == 4);
    REQUIRE(p.SignAt(0.0) == Sign::Positive);
    REQUIRE(p.SignAt(1.0) == Sign::Zero);
    REQUIRE(p.SignAt(1.5) == Sign::Negative);
    REQUIRE(p.SignAt(2.5) == Sign::Positive);
    REQUIRE(p.SignAt(3.5) == Sign::Negative);
    REQUIRE(p.SignAt(5.0) == Sign::Positive);
    REQUIRE(p.SignAtInfinity() == Sign::Positive);
    REQUIRE(p.SignAtMinusInfinity() == Sign::Positive);
    const ExactPolynomial d = Derivative(p);
    REQUIRE(d.GetDegree() == 3);
    REQUIRE(d.SignAtMinusInfinity() == Sign::Negative);
    // 4x^3 - 30x^2 + 70x - 50 at 1: -6.
    REQUIRE(d.SignAt(1.0) == Sign::Negative);
    // Leading zeros are dropped; the zero polynomial has degree -1.
    const ExactPolynomial trimmed(
        std::vector<DynamicExpansion>{DynamicExpansion{2.0}, DynamicExpansion{}});
    REQUIRE(trimmed.GetDegree() == 0);
    REQUIRE(ExactPolynomial{}.GetDegree() == -1);
    REQUIRE(Difference(p, p).IsZero());
}

TEST_CASE("Sturm sequences count distinct real roots in any interval", "[robust][dynamic]")
{
    const SturmSequence sturm(WithRoots({1.0, 2.0, 3.0, 4.0}));
    REQUIRE(sturm.CountRoots(-Infinity, Infinity) == 4);
    REQUIRE(sturm.CountRoots(1.5, 3.5) == 2);
    REQUIRE(sturm.CountRoots(0.0, 1.0) == 1); // (0, 1] holds the root 1.
    REQUIRE(sturm.CountRoots(1.0, 2.0) == 1);
    REQUIRE(sturm.CountRoots(4.0, Infinity) == 0);
    REQUIRE(sturm.CountRoots(-Infinity, 0.5) == 0);

    // A negative leading coefficient and odd reductions: the signed
    // pseudo-remainder keeps the count right.
    const SturmSequence negated(Scaled(WithRoots({-1.0, 0.5, 7.0}), DynamicExpansion{-3.0}));
    REQUIRE(negated.CountRoots(-Infinity, Infinity) == 3);
    REQUIRE(negated.CountRoots(0.0, 1.0) == 1);

    // No real roots, and a multiple root counted once.
    const ExactPolynomial circle(std::vector<DynamicExpansion>{
        DynamicExpansion{1.0}, DynamicExpansion{}, DynamicExpansion{1.0}});
    REQUIRE(SturmSequence(circle).CountRoots(-Infinity, Infinity) == 0);
    REQUIRE(SturmSequence(WithRoots({1.0, 1.0, -2.0})).CountRoots(-Infinity, Infinity) == 2);

    // Degree 0 and the zero polynomial have no roots.
    REQUIRE(SturmSequence(ExactPolynomial(std::vector<DynamicExpansion>{DynamicExpansion{5.0}}))
                .CountRoots(-Infinity, Infinity) == 0);
    REQUIRE(SturmSequence(ExactPolynomial{}).CountRoots(-Infinity, Infinity) == 0);

    rc::prop("the count of a product of random linear factors is the number of distinct roots", [] {
        const auto count = *rc::gen::inRange<std::size_t>(1, 6);
        std::vector<double> roots;
        for (std::size_t i = 0; i < count; ++i) {
            roots.push_back(static_cast<double>(*rc::gen::inRange<int>(-20, 20)) / 4.0);
        }
        std::vector<double> distinct = roots;
        std::ranges::sort(distinct);
        const auto [first, last] = std::ranges::unique(distinct);
        distinct.erase(first, last);
        const SturmSequence s(WithRoots(roots));
        RC_ASSERT(s.CountRoots(-Infinity, Infinity) == static_cast<int>(distinct.size()));
        const double a = static_cast<double>(*rc::gen::inRange<int>(-25, 25)) / 4.0;
        const double b = a + static_cast<double>(*rc::gen::inRange<int>(1, 20)) / 4.0;
        int within = 0;
        for (const double r : distinct) {
            within += (a < r && r <= b) ? 1 : 0;
        }
        RC_ASSERT(s.CountRoots(a, b) == within);
    });
}
