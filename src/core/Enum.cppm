// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <version>

#if defined(__cpp_impl_reflection) && defined(__cpp_lib_reflection)
#include <meta>
#define STEREON_ENUM_NAMES 1
#else
#define STEREON_ENUM_NAMES 0
#endif

export module Stereon.Core:Enum;

/// Enumerator names by static reflection (ADR-0016), for diagnostics and
/// formatting. A new enumerator is named automatically; nothing is written
/// twice. Clang and MSVC have no reflection yet, so there EnumNames is false
/// and EnumName() returns an empty view, and callers fall back to numbers.

export namespace Stereon {

/// True when EnumName() reports names.
inline constexpr bool EnumNames = STEREON_ENUM_NAMES != 0;

} // namespace Stereon

#if STEREON_ENUM_NAMES
namespace Stereon::Detail {

template <typename E>
struct EnumEntry {
    E Value{};
    std::string_view Name;
};

// Built with an ordinary loop: `template for` raises a false -Wshadow
// warning in GCC 16.2 (docs/toolchain/gcc16-feature-probe.md).
template <typename E>
consteval auto MakeEnumNameTable()
{
    constexpr auto enumerators = std::define_static_array(std::meta::enumerators_of(^^E));
    std::array<EnumEntry<E>, enumerators.size()> table{};
    for (std::size_t i = 0; i < enumerators.size(); ++i) {
        table[i] = {std::meta::extract<E>(enumerators[i]),
                    std::define_static_string(std::meta::identifier_of(enumerators[i]))};
    }
    return table;
}

template <typename E>
inline constexpr auto EnumNameTable = MakeEnumNameTable<E>();

} // namespace Stereon::Detail
#endif

export namespace Stereon {

/// The enumerator's name, or an empty view when the value is not an
/// enumerator or the compiler has no reflection.
template <typename E>
    requires std::is_enum_v<E>
[[nodiscard]] constexpr std::string_view EnumName([[maybe_unused]] const E value) noexcept
{
#if STEREON_ENUM_NAMES
    for (const auto& entry : Detail::EnumNameTable<E>) {
        if (entry.Value == value) {
            return entry.Name;
        }
    }
#endif
    return {};
}

} // namespace Stereon
