// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// The parallel loop operations run their per-element work on, and the pool
/// of worker threads it runs on (ADR-0007).
///
/// Sender graphs cannot live in module units (Stereon/Core/Exec.hpp), so the
/// pool and the loop are defined in src/core/Parallel.cpp, a plain source
/// file, and module code reaches them through this header, included in the
/// global module fragment of a unit: the declarations belong to the global
/// module, as the definitions do. The pool is opaque here, so that this
/// header includes no stdexec.
///
/// One pool per Context (ADR-0007): a Context holds a PoolHandle through a
/// shared pointer, made with its thread count; copies of the Context share
/// it, a different thread count makes a different one, and the pool's
/// threads start the first time an operation asks for them. No operation
/// creates threads of its own.
///
/// Determinism (ADR-0007, G1): [0, count) is cut into chunks of `grain`
/// indices whatever the pool, and a body writes only the slots of its own
/// indices, so an operation that does no more than that gives the same bytes
/// for every scheduler and thread count. A merge across indices is the
/// caller's, in index order, after the loop.

#pragma once

#include <cstddef>
#include <memory>
#include <mutex>

namespace Stereon::Detail {

/// A pool of worker threads; defined in Parallel.cpp.
class WorkerPool;

/// The pool a Context owns: made cheaply with a thread count, 0 for the
/// hardware's, its threads started on the first Get() from any thread and
/// joined when the last Context sharing it is destroyed. That destruction
/// must not happen on one of its own threads, and a body run on the pool must
/// not wait on the same pool.
class PoolHandle {
public:
    explicit PoolHandle(unsigned threads) noexcept;
    ~PoolHandle();
    PoolHandle(const PoolHandle&) = delete;
    PoolHandle& operator=(const PoolHandle&) = delete;

    /// The pool, its threads started on the first call. Thread-safe. A
    /// borrow of this handle (ADR-0046): asked of a named handle only.
    [[nodiscard]] WorkerPool& Get() &;
    WorkerPool& Get() && = delete;
    /// The worker count the pool has or will have: the thread count, or the
    /// hardware's for 0.
    [[nodiscard]] unsigned GetWorkerCount() const noexcept;

private:
    unsigned m_Threads;
    std::once_flag m_Started;
    std::unique_ptr<WorkerPool> m_Pool;
};

/// A handle with that thread count, its pool not yet started.
[[nodiscard]] std::shared_ptr<PoolHandle> MakePoolHandle(unsigned threads);

/// The work of one chunk, the indices [begin, end), on the state `context`
/// points to. It must not throw, and must write nothing another chunk
/// writes or reads.
using ParallelBody = void (*)(void* context, std::size_t begin, std::size_t end);

/// The indices of one chunk for cheap, even elements, fixed so that the cut
/// does not depend on the thread count.
inline constexpr std::size_t ParallelGrain = 64;

/// Runs `body` over every chunk of [0, count), chunks of `grain` indices (at
/// least one), exactly once each. With no pool, on the calling thread in
/// index order (the serial scheduler of ADR-0007, for debugging and as the
/// determinism tests' reference); otherwise on the pool's workers, through a
/// bulk sender, which a run at one thread goes through as one at 32 does.
/// Returns when every chunk has run. A coarse grain suits elements whose
/// work is large and uneven, a face pair's intersection say.
void ParallelFor(std::size_t count, std::size_t grain, WorkerPool* pool, void* context,
                 ParallelBody body);

/// In chunks of ParallelGrain.
void ParallelFor(std::size_t count, WorkerPool* pool, void* context, ParallelBody body);

} // namespace Stereon::Detail
