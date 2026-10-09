// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Includes Stereon/Core/Exec.hpp, so this file is not scanned for modules and
// does not import Stereon.Core (stereon_exec_sources in tests/CMakeLists.txt).
#include <Stereon/Core/Exec.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

namespace Exec = Stereon::Exec;

TEST_CASE("then composes on just", "[core][exec]")
{
    auto work = Exec::just(20) | Exec::then([](int value) { return value + 22; });
    auto [result] = Exec::sync_wait(std::move(work)).value();
    REQUIRE(result == 42);
}

TEST_CASE("when_all joins values in argument order", "[core][exec]")
{
    auto [count, length] = Exec::sync_wait(Exec::when_all(Exec::just(3), Exec::just(2.5))).value();
    REQUIRE(count == 3);
    REQUIRE(length == 2.5);
}

TEST_CASE("run_loop runs scheduled work on the thread that drives it", "[core][exec]")
{
    // The serial context ADR-0007 requires for debugging: work started on a
    // run_loop runs on whichever thread calls run().
    Exec::run_loop loop;
    std::thread::id worker;
    auto work = Exec::just() | Exec::then([&] { worker = std::this_thread::get_id(); });
    std::thread waiter([&] {
        Exec::sync_wait(Exec::starts_on(loop.get_scheduler(), std::move(work)));
        loop.finish();
    });
    loop.run();
    waiter.join();
    REQUIRE(worker == std::this_thread::get_id());
}

TEST_CASE("bulk on the parallel scheduler fills per-index slots", "[core][exec][determinism]")
{
    // ADR-0007 pattern: parallel stages write per-index buffers, then reduce in index order.
    constexpr std::size_t count = 10'000;
    std::vector<double> slots(count, 0.0);
    auto work = Exec::schedule(Exec::get_parallel_scheduler()) |
                Exec::bulk(Exec::par, count, [&slots](std::size_t index) {
                    slots[index] = 0.5 * static_cast<double>(index);
                });
    Exec::sync_wait(std::move(work));

    double total = 0.0;
    for (const double value : slots) {
        total += value;
    }
    REQUIRE(total == 0.5 * static_cast<double>(count * (count - 1) / 2));
}
