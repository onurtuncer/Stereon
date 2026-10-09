// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <vector>

import Stereon.Core;

namespace {

namespace Simd = Stereon::Simd;

std::vector<double> MakeData(std::size_t count)
{
    std::vector<double> data(count);
    for (std::size_t i = 0; i < count; ++i) {
        data[i] = static_cast<double>(i % 1024) * 0.25;
    }
    return data;
}

void ScalarSum(benchmark::State& state)
{
    const auto data = MakeData(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        double total = 0.0;
        for (const double value : data) {
            total += value;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * state.range(0));
}

void SimdSum(benchmark::State& state)
{
    const auto data = MakeData(static_cast<std::size_t>(state.range(0)));
    constexpr std::size_t lanes = Simd::Width<double>;
    for (auto _ : state) {
        auto partial = Simd::Broadcast(0.0);
        std::size_t i = 0;
        for (; i + lanes <= data.size(); i += lanes) {
            partial = partial + Simd::Load(data.data() + i);
        }
        double total = Simd::ReduceAdd(partial);
        for (; i < data.size(); ++i) {
            total += data[i];
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * state.range(0));
}

} // namespace

BENCHMARK(ScalarSum)->Arg(1 << 16);
BENCHMARK(SimdSum)->Arg(1 << 16);
