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

module Stereon.Robust;

import :Domain;

namespace Stereon::Detail {

void PredicateDomainViolation(const char* predicate, const DomainStatus status) noexcept
{
    // The domain of the exact predicates (Domain.cppm): finite coordinates
    // whose binary exponents span few enough orders that one power of two
    // keeps every intermediate value a normal double.
    const char* reason =
        status == DomainStatus::NotFinite
            ? "a coordinate is infinite or NaN"
            : "the coordinates' exponents span too many orders to evaluate exactly";
    const std::string text =
        std::format("Stereon: {} outside the predicate domain: {}\n", predicate, reason);
    std::fputs(text.c_str(), stderr);

    const bool coordinatesAreInDomain = false;
    STEREON_ASSERT(coordinatesAreInDomain);
    // Under observe, or without contracts: an inexact sign would make
    // topology inconsistent silently.
    std::abort();
}

} // namespace Stereon::Detail
