// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detail::ParallelFor (Stereon/Core/Parallel.hpp), the loop the parallel
// operations run on (ADR-0007): every index once, in chunks cut the same way
// at every thread count, serially on the calling thread when asked.

#include <Stereon/Core/Parallel.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <thread>
#include <vector>

namespace {

using Stereon::Detail::ParallelFor;
using Stereon::Detail::ParallelGrain;

/// What one run records: how often each index ran, and the chunk each was in.
struct Visits {
    std::vector<int> Count;
    std::vector<std::size_t> ChunkStart;
};

Visits Run(const std::size_t count, const unsigned threads, const bool serial)
{
    Visits visits{std::vector<int>(count, 0), std::vector<std::size_t>(count, 0)};
    ParallelFor(count, threads, serial, &visits,
                [](void* state, const std::size_t begin, const std::size_t end) {
                    Visits& v = *static_cast<Visits*>(state);
                    for (std::size_t i = begin; i < end; ++i) {
                        ++v.Count[i];
                        v.ChunkStart[i] = begin;
                    }
                });
    return visits;
}

} // namespace

TEST_CASE("ParallelFor runs every index once, in chunks independent of the thread count",
          "[core][parallel][determinism]")
{
    for (const std::size_t count : {std::size_t{0}, std::size_t{1}, ParallelGrain - 1,
                                    ParallelGrain, ParallelGrain + 1, std::size_t{5000}}) {
        INFO("count " << count);
        const Visits serial = Run(count, 1, true);
        REQUIRE(serial.Count == std::vector<int>(count, 1));
        for (std::size_t i = 0; i < count; ++i) {
            REQUIRE(serial.ChunkStart[i] == i - i % ParallelGrain);
        }
        for (const unsigned threads : {1U, 2U, 8U, 32U, 0U}) {
            INFO(threads << " threads");
            const Visits parallel = Run(count, threads, false);
            REQUIRE(parallel.Count == serial.Count);
            REQUIRE(parallel.ChunkStart == serial.ChunkStart);
        }
    }
}

TEST_CASE("ParallelFor's serial mode runs on the calling thread", "[core][parallel]")
{
    struct State {
        std::thread::id Caller;
        bool Elsewhere = false;
    } state{std::this_thread::get_id()};
    ParallelFor(3 * ParallelGrain, 32, true, &state, [](void* raw, std::size_t, std::size_t) {
        State& s = *static_cast<State*>(raw);
        s.Elsewhere = s.Elsewhere || std::this_thread::get_id() != s.Caller;
    });
    REQUIRE_FALSE(state.Elsewhere);
}
