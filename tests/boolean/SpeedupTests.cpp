// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The many-face fixture of the Boolean prototype's M10 (ADR-0007): a slab
// with a grid of cylindrical holes, cut by one operand that is every
// cylinder of the grid (tests/support/BooleanCases.hpp). The result must
// not depend on the thread count, at 1, 2, 8 and 32 threads and over
// repeated runs; and, hidden from the default run since it measures the
// machine, the cut is timed serially and on the pool, stage by stage.
// bench/BooleanBench.cpp is the measurement of record.

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "Fixtures.hpp"

using BooleanCases::Drills;
using BooleanCases::Slab;
using BooleanTest::CountOf;
using Stereon::BooleanOperation;
using Stereon::Context;
using Stereon::EntityKind;
using Stereon::OperationResult;
using Stereon::Scheduler;
using Stereon::Shape;

namespace {

const Context serial = Context::Default().WithScheduler(Scheduler::Serial);

OperationResult Cut(const Shape& slab, const Shape& drills, const Context& context)
{
    auto result = Boolean(slab, drills, BooleanOperation::Cut, context);
    const std::string error = result ? std::string() : std::string(result.error().GetMessage());
    INFO(error);
    REQUIRE(result.has_value());
    return std::move(*result);
}

double Seconds(const std::chrono::steady_clock::duration d)
{
    return std::chrono::duration<double>(d).count();
}

} // namespace

TEST_CASE("a slab drilled with a grid of holes is the same at every thread count",
          "[boolean][determinism]")
{
    const int n = 3;
    const Shape slab = Slab(n);
    const Shape drills = Drills(n);
    REQUIRE(CountOf<EntityKind::Solid>(drills) == static_cast<std::size_t>(n * n));
    const OperationResult reference = Cut(slab, drills, serial);
    // The slab's six faces, and one hole's face per drill.
    REQUIRE(CountOf<EntityKind::Face>(reference.Shape) == static_cast<std::size_t>(6 + n * n));
    REQUIRE(CountOf<EntityKind::Solid>(reference.Shape) == 1);
    const std::vector<std::uint64_t> bits = BooleanCases::Bits(reference);
    for (const unsigned threads : {1U, 2U, 8U, 32U, 0U}) {
        CAPTURE(threads);
        const Context parallel =
            Context::Default().WithScheduler(Scheduler::Parallel).WithThreadCount(threads);
        for (int run = 0; run < 2; ++run) {
            REQUIRE(BooleanCases::Bits(Cut(slab, drills, parallel)) == bits);
        }
    }
}

TEST_CASE("the drilled slab's speed-up, stage by stage", "[.][boolean][speedup]")
{
    for (const int n : {7, 10}) {
        using Clock = std::chrono::steady_clock;
        const auto built = Clock::now();
        const Shape slab = Slab(n);
        const Shape drills = Drills(n);
        std::printf("holes %4d: the drills built in %7.3f s\n", n * n,
                    Seconds(Clock::now() - built));
        for (const unsigned threads : {1U, 8U}) {
            const Context context = threads == 1U ? serial
                                                  : Context::Default()
                                                        .WithScheduler(Scheduler::Parallel)
                                                        .WithThreadCount(threads);
            const auto t0 = Clock::now();
            const auto imprinted = ImprintBodies(slab, drills, context);
            REQUIRE(imprinted.has_value());
            const auto t1 = Clock::now();
            const auto classified = Classify(*imprinted, context);
            REQUIRE(classified.has_value());
            const auto t2 = Clock::now();
            const OperationResult cut = Cut(slab, drills, context);
            const auto t3 = Clock::now();
            std::printf("holes %4d threads %2u: imprint %7.3f s, classify %7.3f s, whole cut "
                        "%7.3f s (%zu faces)\n",
                        n * n, threads, Seconds(t1 - t0), Seconds(t2 - t1), Seconds(t3 - t2),
                        CountOf<EntityKind::Face>(cut.Shape));
        }
    }
}
