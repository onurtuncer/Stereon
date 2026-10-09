// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

import Stereon.Core;

namespace Simd = Stereon::Simd;

namespace {

constexpr std::size_t Lanes = Simd::Width<double>;

// Multiples of 1/64 below 2^14 in magnitude: sums, products and square roots
// of squares are exact in double, so lanes can be compared with ==.
rc::Gen<double> DyadicDouble()
{
    return rc::gen::map(rc::gen::inRange(-1'000'000, 1'000'001),
                        [](int value) { return static_cast<double>(value) / 64.0; });
}

std::vector<double> RandomLanes()
{
    return *rc::gen::container<std::vector<double>>(Lanes, DyadicDouble());
}

std::vector<double> Stored(const Simd::Vec<double>& value)
{
    std::vector<double> lanes(Lanes);
    Simd::Store(value, lanes.data());
    return lanes;
}

} // namespace

TEST_CASE("Simd vectors have at least one lane", "[core][simd]")
{
    STATIC_REQUIRE(Simd::Width<double> >= 1);
    STATIC_REQUIRE(Simd::Width<float> >= Simd::Width<double>);
}

TEST_CASE("Simd operations match scalar code lane by lane", "[core][simd][property]")
{
    rc::prop("load then store round-trips", [] {
        const auto lanes = RandomLanes();
        RC_ASSERT(Stored(Simd::Load(lanes.data())) == lanes);
    });

    rc::prop("broadcast fills every lane", [] {
        const double value = *DyadicDouble();
        RC_ASSERT(Stored(Simd::Broadcast(value)) == std::vector<double>(Lanes, value));
    });

    rc::prop("arithmetic operators act lane-wise", [] {
        const auto a = RandomLanes();
        const auto b = RandomLanes();
        const auto sum = Stored(Simd::Load(a.data()) + Simd::Load(b.data()));
        const auto product = Stored(Simd::Load(a.data()) * Simd::Load(b.data()));
        for (std::size_t i = 0; i < Lanes; ++i) {
            RC_ASSERT(sum[i] == a[i] + b[i]);
            RC_ASSERT(product[i] == a[i] * b[i]);
        }
    });

    rc::prop("ReduceAdd equals the scalar sum", [] {
        const auto lanes = RandomLanes();
        double expected = 0.0;
        for (const double value : lanes) {
            expected += value;
        }
        RC_ASSERT(Simd::ReduceAdd(Simd::Load(lanes.data())) == expected);
    });

    rc::prop("Min, Max, Abs, Sqrt and Fma act lane-wise", [] {
        const auto a = RandomLanes();
        const auto b = RandomLanes();
        const auto c = RandomLanes();
        const auto va = Simd::Load(a.data());
        const auto vb = Simd::Load(b.data());
        const auto vc = Simd::Load(c.data());
        const auto low = Stored(Simd::Min(va, vb));
        const auto high = Stored(Simd::Max(va, vb));
        const auto magnitude = Stored(Simd::Abs(va));
        const auto root = Stored(Simd::Sqrt(va * va));
        const auto fused = Stored(Simd::Fma(va, vb, vc));
        for (std::size_t i = 0; i < Lanes; ++i) {
            RC_ASSERT(low[i] == std::min(a[i], b[i]));
            RC_ASSERT(high[i] == std::max(a[i], b[i]));
            RC_ASSERT(magnitude[i] == std::abs(a[i]));
            RC_ASSERT(root[i] == std::abs(a[i]));
            RC_ASSERT(fused[i] == std::fma(a[i], b[i], c[i]));
        }
    });
}
