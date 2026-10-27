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
//
// The control: batch projection of 200 000 points onto the determinism
// checks' rational bicubic surface (tests/support/DeterminismCases.hpp), the
// same way, serially and on 4 threads. It is parallel through and through,
// no serial part and no merge, so its ratio is what the machine's 4 threads
// can give any work: the ceiling the Boolean's ratio is read against.

#include <benchmark/benchmark.h>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "DeterminismCases.hpp"

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

/// The control's points, made once.
constexpr std::size_t ControlPoints = 200'000;

void Project(benchmark::State& state, const Stereon::Context& context)
{
    static const std::vector<Stereon::Point3> points = DeterminismTest::Points(ControlPoints);
    static const Stereon::SurfaceProjector projector(DeterminismTest::Surface());
    for (auto _ : state) {
        benchmark::DoNotOptimize(ProjectAll(projector, points, context));
    }
}

void BM_ControlSerial(benchmark::State& state)
{
    Project(state, Stereon::Context::Default().WithScheduler(Stereon::Scheduler::Serial));
}

void BM_ControlParallel(benchmark::State& state)
{
    static const Stereon::Context parallel =
        Stereon::Context::Default().WithScheduler(Stereon::Scheduler::Parallel).WithThreadCount(4);
    Project(state, parallel);
}

} // namespace

BENCHMARK(BM_ControlSerial)
    ->Name("M10/Control/serial")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Iterations(1);
BENCHMARK(BM_ControlParallel)
    ->Name("M10/Control/threads:4")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Iterations(1);

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
