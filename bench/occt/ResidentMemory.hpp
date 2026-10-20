// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

namespace StereonBench {

/// The process's resident set, in bytes; 0 where the platform is not known.
/// Defined in ResidentMemory.cpp, which imports no module: <windows.h>
/// brings in the x86 intrinsics, which on MinGW conflict with the ones a
/// Stereon module already declares (GCC 16.2).
[[nodiscard]] std::size_t ResidentBytes();

} // namespace StereonBench
