// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Catch2 output for the robust library's types. Include it after
/// `import Stereon.Robust;`.

#pragma once

#include <catch2/catch_tostring.hpp>

#include <cstddef>
#include <format>
#include <string>

namespace Catch {

template <>
struct StringMaker<Stereon::Sign> {
    static std::string convert(const Stereon::Sign s) { return std::format("{}", s); }
};

template <std::size_t N>
struct StringMaker<Stereon::Expansion<N>> {
    static std::string convert(const Stereon::Expansion<N>& e)
    {
        std::string text = "[";
        for (const double c : e.GetComponents()) {
            text += std::format("{}{}", text.size() > 1 ? ", " : "", c);
        }
        return text + "]";
    }
};

} // namespace Catch
