// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ResidentMemory.hpp"

#if defined(_WIN32)
#include <windows.h>
// After windows.h.
#include <psapi.h>
#elif defined(__linux__)
#include <unistd.h>

#include <fstream>
#endif

namespace StereonBench {

std::size_t ResidentBytes()
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS counters{};
    counters.cb = sizeof counters;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, counters.cb) != 0) {
        return counters.WorkingSetSize;
    }
    return 0;
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    std::size_t pages = 0;
    std::size_t resident = 0;
    if (statm >> pages >> resident) {
        return resident * static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    }
    return 0;
#else
    return 0;
#endif
}

} // namespace StereonBench
