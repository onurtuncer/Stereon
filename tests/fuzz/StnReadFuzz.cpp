// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Fuzz target: the .stn reader (ADR-0012), the kernel's first entry point for
// bytes from outside. Any input must give a LoadedGeometry or a KernelError,
// never a crash or a contract violation (ADR-0006, ADR-0016). An archive
// that loads must write again and read back: a file the reader accepts is
// one the writer can store.

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <span>
#include <vector>

import Stereon.Io;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, const std::size_t size)
{
    const std::span<const std::byte> bytes(reinterpret_cast<const std::byte*>(data), // NOLINT
                                           size);
    const Stereon::Result<Stereon::LoadedGeometry> loaded = Stereon::ReadGeometry(bytes);
    static_cast<void>(Stereon::GeometryToJson(bytes));
    if (!loaded) {
        return 0;
    }
    const Stereon::Result<std::vector<std::byte>> written = Stereon::WriteGeometry(loaded->Archive);
    if (!written || !Stereon::ReadGeometry(*written)) {
        std::abort(); // the reader accepted what the writer cannot store
    }
    return 0;
}
