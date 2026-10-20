// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

// Defined in Memory.cpp, which imports no module: <windows.h> brings in the
// x86 intrinsics, which on MinGW conflict with the ones a Stereon module
// already declares (GCC 16.2).

namespace StereonBench {

/// The process's resident set, in bytes; 0 where the platform is not known.
/// It does not grow when an allocation reuses pages a previous one freed, so
/// a benchmark that builds and frees a body in every iteration sees 0 after
/// the first; HeapBytes() is the measure for memory per entity.
[[nodiscard]] std::size_t ResidentBytes();

/// The bytes the heap holds in live allocations, from glibc's mallinfo2():
/// the arenas' in-use blocks and the mmapped ones. The difference across a
/// build, with the result still alive, is what the result occupies, freed
/// pages reused or not, and it counts every allocation through malloc, so
/// OCCT's (built without TBB, whose allocator would bypass it) and Stereon's
/// alike. Falls back to ResidentBytes() without glibc.
[[nodiscard]] std::size_t HeapBytes();

/// True when HeapBytes() counts live allocations rather than the resident
/// set: a report says which measure its memory rows used.
[[nodiscard]] bool HeapBytesAreExact();

} // namespace StereonBench
