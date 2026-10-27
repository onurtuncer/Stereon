// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The speed-up of ADR-0007's verification (the Boolean prototype's M10): a
// Boolean on a body of more than 1 000 faces, the slab drilled with a grid
// of 20 x 20 holes (one operand of 400 cylinders, 1 200 faces), cut on the
// serial scheduler and on the pool of 4 threads, the hosted runner's cores.
// The OCCT reference workflow runs both with interleaved repetitions and
// reports the ratio of their medians; ADR-0007 asks for at least 3.
// The fixture is tests/support/BooleanCases.hpp's, made once.

#include <benchmark/benchmark.h>

#include <bit>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"

namespace {

/// The holes across the slab: 20 x 20, 1 200 faces in the drills.
constexpr int Holes = 20;

struct Fixture {
    Stereon::Shape Slab;
    Stereon::Shape Drills;
};

const Fixture& DrilledSlab()
{
    static const Fixture fixture{BooleanCases::Slab(Holes), BooleanCases::Drills(Holes)};
    return fixture;
}

/// The cut under `context`; the thread count 0 is the serial scheduler.
void Cut(benchmark::State& state, const Stereon::Context& context)
{
    const Fixture& f = DrilledSlab();
    for (auto _ : state) {
        auto result = Boolean(f.Slab, f.Drills, Stereon::BooleanOperation::Cut, context);
        if (!result) {
            state.SkipWithError(std::string(result.error().GetMessage()));
            return;
        }
        benchmark::DoNotOptimize(result);
    }
    state.counters["faces"] =
        static_cast<double>(BooleanCases::CountOf<Stereon::EntityKind::Face>(f.Drills) +
                            BooleanCases::CountOf<Stereon::EntityKind::Face>(f.Slab));
}

void BM_DrilledSlabCutSerial(benchmark::State& state)
{
    Cut(state, Stereon::Context::Default().WithScheduler(Stereon::Scheduler::Serial));
}

void BM_DrilledSlabCutParallel(benchmark::State& state)
{
    // One Context for every repetition: its pool is started once (ADR-0007).
    static const Stereon::Context parallel =
        Stereon::Context::Default().WithScheduler(Stereon::Scheduler::Parallel).WithThreadCount(4);
    Cut(state, parallel);
}

} // namespace

BENCHMARK(BM_DrilledSlabCutSerial)
    ->Name("M10/DrilledSlabCut/serial")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Iterations(1);
BENCHMARK(BM_DrilledSlabCutParallel)
    ->Name("M10/DrilledSlabCut/threads:4")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Iterations(1);
