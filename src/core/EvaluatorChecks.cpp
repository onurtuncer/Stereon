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
#include <string>
#include <string_view>

module Stereon.Core;

import :EvaluatorChecks;

namespace Stereon::Detail {

void EvaluatorViolation(const std::string_view function, const std::string_view problem) noexcept
{
    // Concatenated, not std::format: formatting string_views here and in an
    // importer makes GCC 16.2 on MinGW define a libstdc++ static twice
    // (docs/toolchain/gcc16-feature-probe.md).
    std::string text = "Stereon: ";
    text += function;
    text += " called with ";
    text += problem;
    text += '\n';
    std::fputs(text.c_str(), stderr);

    const bool evaluatorPreconditionHolds = false;
    STEREON_ASSERT(evaluatorPreconditionHolds);
    // Under observe, or without contracts: evaluating anyway would return a
    // meaningless point to code that trusts it.
    std::abort();
}

} // namespace Stereon::Detail
