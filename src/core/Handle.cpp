// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <string>

module Stereon.Core;

import :Handle;

namespace Stereon::Detail {

void HandleCheckFailed(const HandleReport& report)
{
    // The report of ADR-0004: the handle, the entity kind and the expected
    // store, then whatever the slot itself says.
    std::string text = std::format(
        "Stereon: invalid handle {}#{}@{} from store {}: the table is store {} with {} slots",
        ToString(report.Kind), report.Index, report.Generation, report.Store, report.ExpectedStore,
        report.Size);
    if (report.Store != report.ExpectedStore) {
        text += " (handle from another store)";
    } else if (report.Index >= report.Size) {
        text += " (index out of range)";
    } else if (!report.SlotAlive) {
        text += std::format(" (slot released, now at generation {})", report.SlotGeneration);
    } else {
        text += std::format(" (slot reused, now at generation {})", report.SlotGeneration);
    }
    text += '\n';
    std::fputs(text.c_str(), stderr);

    const bool handleIsValid = false;
    STEREON_ASSERT(handleIsValid);
#if !STEREON_HAS_CONTRACTS
    // Without contracts the build has no way to report a programming error
    // but this; an unchecked stale handle would corrupt a shape silently.
    std::abort();
#endif
}

StoreId NextStoreId() noexcept
{
    // 0 is "no store".
    static std::atomic<StoreId> s_Next{1};
    return s_Next.fetch_add(1, std::memory_order_relaxed);
}

} // namespace Stereon::Detail
