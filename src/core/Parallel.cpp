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
#include <memory>
#include <mutex>
#include <thread>

namespace Stereon::Detail {

class WorkerPool {
public:
    explicit WorkerPool(const unsigned workers)
        : Pool(workers)
    {}

    Exec::static_thread_pool Pool;
};

namespace {

[[nodiscard]] unsigned WorkersFor(const unsigned threads) noexcept
{
    return threads != 0 ? threads : std::max(1U, std::thread::hardware_concurrency());
}

} // namespace

PoolHandle::PoolHandle(const unsigned threads) noexcept
    : m_Threads(threads)
{}

PoolHandle::~PoolHandle() = default;

WorkerPool& PoolHandle::Get() &
{
    std::call_once(m_Started,
                   [this] { m_Pool = std::make_unique<WorkerPool>(WorkersFor(m_Threads)); });
    return *m_Pool;
}

unsigned PoolHandle::GetWorkerCount() const noexcept
{
    return WorkersFor(m_Threads);
}

std::shared_ptr<PoolHandle> MakePoolHandle(const unsigned threads)
{
    return std::make_shared<PoolHandle>(threads);
}

void ParallelFor(const std::size_t count, const std::size_t grain, WorkerPool* const pool,
                 void* const context, const ParallelBody body)
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
    if (pool == nullptr) {
        for (std::size_t chunk = 0; chunk < chunks; ++chunk) {
            run(chunk);
        }
        return;
    }
    Exec::sync_wait(Exec::schedule(pool->Pool.get_scheduler()) |
                    Exec::bulk(Exec::par, chunks, run));
}

void ParallelFor(const std::size_t count, WorkerPool* const pool, void* const context,
                 const ParallelBody body)
{
    ParallelFor(count, ParallelGrain, pool, context, body);
}

} // namespace Stereon::Detail
