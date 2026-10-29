// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detail::ParallelFor (Stereon/Core/Parallel.hpp), the loop the parallel
// operations run on (ADR-0007): every index once, in chunks cut the same way
// at every thread count, serially on the calling thread without a pool; and
// the pool a Context owns, shared by its copies, started on first use, one
// per thread count.

#include <Stereon/Core/Parallel.hpp>

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <memory>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>

import Stereon.Core;

namespace {

using Stereon::Context;
using Stereon::Scheduler;
using Stereon::Detail::MakePoolHandle;
using Stereon::Detail::ParallelFor;
using Stereon::Detail::ParallelGrain;
using Stereon::Detail::WorkerPool;

/// What one run records: how often each index ran, and the chunk each was in.
struct Visits {
    std::vector<int> Count;
    std::vector<std::size_t> ChunkStart;
};

Visits Run(const std::size_t count, const std::size_t grain, WorkerPool* const pool)
{
    Visits visits{std::vector<int>(count, 0), std::vector<std::size_t>(count, 0)};
    ParallelFor(count, grain, pool, &visits,
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
    for (const std::size_t grain : {std::size_t{1}, std::size_t{7}, ParallelGrain}) {
        for (const std::size_t count : {std::size_t{0}, std::size_t{1}, ParallelGrain - 1,
                                        ParallelGrain, ParallelGrain + 1, std::size_t{5000}}) {
            INFO("count " << count << ", grain " << grain);
            const Visits serial = Run(count, grain, nullptr);
            REQUIRE(serial.Count == std::vector<int>(count, 1));
            for (std::size_t i = 0; i < count; ++i) {
                REQUIRE(serial.ChunkStart[i] == i - i % grain);
            }
            for (const unsigned threads : {1U, 2U, 8U, 32U, 0U}) {
                INFO(threads << " threads");
                const auto handle = MakePoolHandle(threads);
                const Visits parallel = Run(count, grain, &handle->Get());
                REQUIRE(parallel.Count == serial.Count);
                REQUIRE(parallel.ChunkStart == serial.ChunkStart);
            }
        }
    }
}

TEST_CASE("ParallelFor joins every started chunk before an exception leaves it", "[core][parallel]")
{
    // ADR-0046, section 5 (AC-0046-05): the loop completes or joins all
    // started work before unwinding past the caller's state, which here is
    // on the stack. One chunk throws; every other chunk that started
    // finishes; nothing runs once the exception has propagated.
    struct State {
        std::atomic<int> Started{0};
        std::atomic<int> Finished{0};
        std::atomic<bool> Returned{false};
        std::atomic<bool> RanAfterReturn{false};
    };
    for (const unsigned threads : {1U, 4U}) {
        INFO(threads << " threads");
        State state;
        const auto handle = MakePoolHandle(threads);
        bool thrown = false;
        try {
            ParallelFor(64, 1, &handle->Get(), &state,
                        [](void* raw, const std::size_t begin, std::size_t) {
                            State& s = *static_cast<State*>(raw);
                            if (s.Returned.load()) {
                                s.RanAfterReturn = true;
                            }
                            ++s.Started;
                            std::this_thread::sleep_for(std::chrono::microseconds(50));
                            if (begin == 7) {
                                throw std::runtime_error("chunk 7");
                            }
                            ++s.Finished;
                        });
        } catch (const std::runtime_error&) {
            thrown = true;
        }
        state.Returned = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        REQUIRE(thrown);
        REQUIRE(state.Started.load() == state.Finished.load() + 1);
        REQUIRE_FALSE(state.RanAfterReturn.load());
    }
}

TEST_CASE("ParallelFor without a pool runs on the calling thread", "[core][parallel]")
{
    struct State {
        std::thread::id Caller;
        bool Elsewhere = false;
    } state{std::this_thread::get_id()};
    ParallelFor(3 * ParallelGrain, nullptr, &state, [](void* raw, std::size_t, std::size_t) {
        State& s = *static_cast<State*>(raw);
        s.Elsewhere = s.Elsewhere || std::this_thread::get_id() != s.Caller;
    });
    REQUIRE_FALSE(state.Elsewhere);
    // So does a serial Context's work: it has no pool to offer.
    REQUIRE(Context::Default().WithScheduler(Scheduler::Serial).GetWorkerPool() == nullptr);
}

TEST_CASE("a Context owns one pool, shared by its copies, one per thread count", "[core][parallel]")
{
    const Context four = Context::Default().WithThreadCount(4);
    WorkerPool* const pool = four.GetWorkerPool();
    REQUIRE(pool != nullptr);
    // Copies share it, whatever else they change; asking again starts nothing new.
    REQUIRE(four.GetWorkerPool() == pool);
    REQUIRE(four.WithResolution(1e-4).GetWorkerPool() == pool);
    REQUIRE(four.WithThreadCount(4).GetWorkerPool() == pool);
    // The serial scheduler keeps the pool and does not use it; back on the
    // parallel scheduler, it is the same pool.
    const Context serial = four.WithScheduler(Scheduler::Serial);
    REQUIRE(serial.GetWorkerPool() == nullptr);
    REQUIRE(serial.WithScheduler(Scheduler::Parallel).GetWorkerPool() == pool);
    // Another thread count is another pool.
    const Context eight = four.WithThreadCount(8);
    REQUIRE(eight.GetWorkerPool() != pool);
    // And it runs no more workers than its count: at most four distinct threads.
    std::set<std::thread::id> seen;
    std::atomic<bool> busy{false};
    struct State {
        std::set<std::thread::id>* Seen;
        std::atomic<bool>* Busy;
    } state{&seen, &busy};
    ParallelFor(400, 1, pool, &state, [](void* raw, std::size_t, std::size_t) {
        State& s = *static_cast<State*>(raw);
        bool expected = false;
        while (!s.Busy->compare_exchange_weak(expected, true)) {
            expected = false;
        }
        s.Seen->insert(std::this_thread::get_id());
        s.Busy->store(false);
    });
    REQUIRE(seen.size() >= 1);
    REQUIRE(seen.size() <= 4);
}

TEST_CASE("a Context's pool is started on first use from any thread, once", "[core][parallel]")
{
    const Context shared = Context::Default().WithThreadCount(3);
    std::vector<WorkerPool*> pools(8, nullptr);
    {
        std::vector<std::jthread> threads;
        for (std::size_t i = 0; i < pools.size(); ++i) {
            threads.emplace_back([&shared, &pools, i] { pools[i] = shared.GetWorkerPool(); });
        }
    }
    for (WorkerPool* const pool : pools) {
        REQUIRE(pool != nullptr);
        REQUIRE(pool == pools.front());
    }
}
