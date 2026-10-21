// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Private implementation header: include in the global module fragment of
// implementation units only, never in an exported interface or public header.
#ifndef STEREON_ENABLE_PROFILING
#define STEREON_ENABLE_PROFILING 0
#endif

#if STEREON_ENABLE_PROFILING

#include <tracy/TracyC.h>

#include <cstdint>
#include <type_traits>

namespace Stereon::Detail::Profiling {

// The C API's literal context lets this RAII guard exist during constant
// evaluation. Tracy's runtime C++ ScopedZone is not a literal type.
class Scope final {
public:
    constexpr explicit Scope(const ___tracy_source_location_data* location) noexcept
    {
        if (!std::is_constant_evaluated()) {
            m_Context = ___tracy_emit_zone_begin(location, 1);
        }
    }

    constexpr ~Scope()
    {
        if (!std::is_constant_evaluated()) {
            ___tracy_emit_zone_end(m_Context);
        }
    }

    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&&) = delete;
    Scope& operator=(Scope&&) = delete;

private:
    TracyCZoneCtx m_Context{};
};

} // namespace Stereon::Detail::Profiling

#define STEREON_PROFILE_JOIN_IMPL(a, b) a##b
#define STEREON_PROFILE_JOIN(a, b) STEREON_PROFILE_JOIN_IMPL(a, b)

// Names must be string literals. Static descriptors outlive asynchronous
// capture. One scope per source line; braces delimit nested scopes.
// Do not span a coroutine suspension or move a scope between threads.
#define STEREON_PROFILE_SCOPE(name)                                                                \
    static constexpr ___tracy_source_location_data STEREON_PROFILE_JOIN(                           \
        stereonLocation, __LINE__){name, __func__, __FILE__, static_cast<std::uint32_t>(__LINE__), \
                                   0};                                                             \
    [[maybe_unused]] const ::Stereon::Detail::Profiling::Scope STEREON_PROFILE_JOIN(               \
        stereonScope, __LINE__)(&STEREON_PROFILE_JOIN(stereonLocation, __LINE__))

// A sample, not an increment of a shared counter. Arguments must be
// observational and the value representable as int64_t. No evaluation in
// constant evaluation or a disabled build.
#define STEREON_PROFILE_COUNTER(name, value)                                                       \
    do {                                                                                           \
        if (!std::is_constant_evaluated()) {                                                       \
            ___tracy_emit_plot_int(name, static_cast<std::int64_t>(value));                        \
        }                                                                                          \
    } while (false)

#else

#define STEREON_PROFILE_SCOPE(name) ((void)0)
#define STEREON_PROFILE_COUNTER(name, value) ((void)0)

#endif
