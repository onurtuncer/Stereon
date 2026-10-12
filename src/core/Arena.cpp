// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstdio>
#include <cstdlib>
#include <format>
#include <string>
#include <string_view>

module Stereon.Core;

import :Arena;

namespace Stereon::Detail {

void ArenaEditViolation(const std::string_view operation) noexcept
{
    // The arena's borrowing rule (Arena.cppm): while Modify() or ModifyPage()
    // runs its callback, the reference it handed out must stay the arena's
    // own, so nothing may copy, move, resize or edit the arena again.
    const std::string text = std::format(
        "Stereon: PagedArena {} during an edit callback: the callback's reference would be "
        "shared or invalidated\n",
        operation);
    std::fputs(text.c_str(), stderr);

    const bool arenaIsNotBeingEdited = false;
    STEREON_ASSERT(arenaIsNotBeingEdited);
    // Under the observe semantic, or without contracts, execution reaches
    // here; continuing would corrupt a snapshot silently.
    std::abort();
}

} // namespace Stereon::Detail
