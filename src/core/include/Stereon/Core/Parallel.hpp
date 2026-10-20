// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// The parallel loop operations run their per-element work on (ADR-0007).
///
/// Sender graphs cannot live in module units (Stereon/Core/Exec.hpp), so the
/// loop is defined in src/core/Parallel.cpp, a plain source file, and module
/// code reaches it through this header, included in the global module
/// fragment of an implementation unit: the declaration belongs to the global
/// module, as the definition does. It takes the scheduler kind and thread
/// count of a Context as plain values, since this header imports no module.
///
/// Determinism (ADR-0007, G1): [0, count) is cut into chunks of ParallelGrain
/// indices whatever the thread count, and a body writes only the slots of its
/// own indices, so an operation that does no more than that gives the same
/// bytes for every scheduler and thread count. A merge across indices is the
/// caller's, in index order, after the loop.

#pragma once

#include <cstddef>

namespace Stereon::Detail {

/// The work of one chunk, the indices [begin, end), on the state `context`
/// points to. It must not throw, and must write nothing another chunk
/// writes or reads.
using ParallelBody = void (*)(void* context, std::size_t begin, std::size_t end);

/// The indices of one chunk, fixed so that the cut does not depend on the
/// thread count.
inline constexpr std::size_t ParallelGrain = 64;

/// Runs `body` over every chunk of [0, count) exactly once. With `serial`,
/// on the calling thread in index order (the serial scheduler of ADR-0007,
/// for debugging and as the determinism tests' reference); otherwise on a
/// pool of `threads` workers, the hardware's count for 0, through a bulk
/// sender. Returns when every chunk has run.
void ParallelFor(std::size_t count, unsigned threads, bool serial, void* context,
                 ParallelBody body);

} // namespace Stereon::Detail
