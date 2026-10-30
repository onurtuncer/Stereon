// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The broad phase against testing every pair, on the boxes of two 100 x 100
// grids of unit cubes two apart, the second shifted by half a unit: 10 000
// boxes a side, 100 million pairs, 10 000 that meet. And the opposite case,
// where testing every pair is the answer: 1 000 boxes a side that all meet,
// a million pairs, every one listed.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <vector>

import Stereon.Core;

namespace {

using Stereon::Box3;
using Stereon::BoxPair;
using Stereon::Point3;

std::vector<Box3> Grid(const int n, const double shift)
{
    std::vector<Box3> boxes;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const Point3 lo{2.0 * i + shift, 2.0 * j, 0.0};
            boxes.push_back(Box3::FromCorners(lo, Point3{lo.X + 1.0, lo.Y + 1.0, 1.0}));
        }
    }
    return boxes;
}

std::vector<Box3> Heap(const int n)
{
    std::vector<Box3> boxes;
    for (int k = 0; k < n; ++k) {
        const double d = 1e-3 * k;
        boxes.push_back(Box3::FromCorners(Point3{d, d, d}, Point3{10.0 + d, 10.0 + d, 10.0 + d}));
    }
    return boxes;
}

/// The loop the broad phase replaced.
std::vector<BoxPair> AllPairs(const std::vector<Box3>& a, const std::vector<Box3>& b)
{
    std::vector<BoxPair> pairs;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            if (Intersects(a[i], b[j])) {
                pairs.push_back({i, j});
            }
        }
    }
    return pairs;
}

struct Lists {
    std::vector<Box3> A;
    std::vector<Box3> B;
};

const Lists& Disjoint()
{
    static const Lists lists{Grid(100, 0.0), Grid(100, 0.5)};
    return lists;
}

const Lists& Overlapping()
{
    static const Lists lists{Heap(1000), Heap(1000)};
    return lists;
}

void BM_BroadPhase(benchmark::State& state, const Lists& lists)
{
    for (auto _ : state) {
        benchmark::DoNotOptimize(OverlappingPairs(lists.A, lists.B));
    }
}

void BM_AllPairs(benchmark::State& state, const Lists& lists)
{
    for (auto _ : state) {
        benchmark::DoNotOptimize(AllPairs(lists.A, lists.B));
    }
}

} // namespace

BENCHMARK_CAPTURE(BM_BroadPhase, disjoint, Disjoint())
    ->Name("Broad/Pairs/10000x10000-disjoint/hierarchy")
    ->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BM_AllPairs, disjoint, Disjoint())
    ->Name("Broad/Pairs/10000x10000-disjoint/all-pairs")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(BM_BroadPhase, overlapping, Overlapping())
    ->Name("Broad/Pairs/1000x1000-overlapping/hierarchy")
    ->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BM_AllPairs, overlapping, Overlapping())
    ->Name("Broad/Pairs/1000x1000-overlapping/all-pairs")
    ->Unit(benchmark::kMillisecond);
