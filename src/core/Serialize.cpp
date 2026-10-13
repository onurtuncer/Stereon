// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <format>
#include <string>

module Stereon.Core;

import :Serialize;

namespace Stereon::Detail {

std::string TruncatedInputMessage(const std::size_t needed, const std::size_t offset,
                                  const std::size_t available)
{
    return std::format("truncated input: {} bytes needed at offset {}, {} available", needed,
                       offset, available);
}

std::string InvalidBoolMessage(const std::uint64_t value, const std::size_t offset)
{
    return std::format("invalid bool byte {} at offset {}", value, offset);
}

std::string LeftoverBytesMessage(const std::size_t count)
{
    return std::format("{} bytes left over after the value", count);
}

void AppendJsonDouble(const double x, std::string& out)
{
    if (x != x) {
        out += "\"nan\"";
    } else if (x - x != 0.0) {
        out += x > 0.0 ? "\"inf\"" : "\"-inf\"";
    } else {
        out += std::format("{}", x); // the shortest form that reads back exactly
    }
}

void AppendJsonInteger(const long long x, std::string& out)
{
    out += std::format("{}", x);
}

void AppendJsonUnsigned(const unsigned long long x, std::string& out)
{
    out += std::format("{}", x);
}

} // namespace Stereon::Detail
