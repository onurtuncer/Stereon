// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Includes Stereon/Core/Exec.hpp, so this file is not scanned for modules and
// imports none (stereon_exec_sources in src/core/CMakeLists.txt).
#include <Stereon/Core/Exec.hpp>
#include <Stereon/Core/Parallel.hpp>

#include <algorithm>
#include <cstddef>
#include <thread>

namespace Stereon::Detail {

void ParallelFor(const std::size_t count, const unsigned threads, const bool serial,
                 void* const context, const ParallelBody body)
{
    ParallelFor(count, ParallelGrain, threads, serial, context, body);
}

void ParallelFor(const std::size_t count, const std::size_t grain, const unsigned threads,
                 const bool serial, void* const context, const ParallelBody body)
{
    const std::size_t step = std::max<std::size_t>(1, grain);
    const std::size_t chunks = (count + step - 1) / step;
    const auto run = [&](const std::size_t chunk) {
        const std::size_t begin = chunk * step;
        body(context, begin, std::min(count, begin + step));
    };
    if (chunks == 0) {
        return;
    }
    if (serial) {
        for (std::size_t chunk = 0; chunk < chunks; ++chunk) {
            run(chunk);
        }
        return;
    }
    // A pool of exactly the Context's thread count, one worker included, so
    // that a run at one thread goes through the same sender graph as one at
    // 32 and differs from the serial reference only in where it runs.
    const unsigned workers =
        threads != 0 ? threads : std::max(1U, std::thread::hardware_concurrency());
    Exec::static_thread_pool pool(workers);
    Exec::sync_wait(Exec::schedule(pool.get_scheduler()) | Exec::bulk(Exec::par, chunks, run));
}

} // namespace Stereon::Detail
