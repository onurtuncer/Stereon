// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdint>
#include <format>
#include <functional>
#include <string_view>

export module Stereon.Core:Log;

import :Enum;
import :Scalar;

/// Logging (ADR-0007, ADR-0023). The kernel formats its messages with
/// std::format and hands them to a LogSink owned by the Context; it has no
/// logging library of its own and never writes anywhere by itself. The
/// application supplies the sink, which may wrap spdlog or anything else.
/// Log output is diagnostic only: no result depends on whether a sink is
/// attached or what it does.

export namespace Stereon {

/// Severity, in increasing order. Off is a threshold that drops everything.
enum class LogLevel : std::uint8_t {
    Trace = 0,
    Debug,
    Info,
    Warning,
    Error,
    Off,
};

/// "Info", "Warning", ...; "LogLevel" where the compiler has no reflection.
[[nodiscard]] constexpr std::string_view ToString(const LogLevel level) noexcept
{
    const std::string_view name = EnumName(level);
    return name.empty() ? std::string_view{"LogLevel"} : name;
}

/// Where a log call was made. The kernel's own stand-in for
/// std::source_location, whose type GCC 16.2 treats as TU-local inside a
/// module interface (docs/toolchain/gcc16-feature-probe.md); Current() uses
/// the same compiler builtins, evaluated at the call site when it appears
/// as a default argument.
struct SourceLocation {
    const char* FileName = "";
    const char* FunctionName = "";
    std::uint_least32_t Line = 0;

    [[nodiscard]] static constexpr SourceLocation
    Current(const char* file = __builtin_FILE(), const char* function = __builtin_FUNCTION(),
            const std::uint_least32_t line = __builtin_LINE()) noexcept
    {
        return {file, function, line};
    }
};

/// One message on its way to the sink. The views are valid for the duration
/// of the sink call only; a sink that keeps the message copies it.
struct LogRecord {
    LogLevel Level = LogLevel::Info;
    std::string_view Message;
    SourceLocation Location;
};

/// Where log records go. Called on the thread that logs, which will be a
/// worker thread once operations run in parallel, and possibly from several
/// threads at once; a sink is thread-safe or serialises itself.
using LogSink = std::function<void(const LogRecord&)>;

/// A sink that writes "level: message (file:line)" lines to stderr, for
/// tests and tools. Defined in Log.cpp.
[[nodiscard]] LogSink StderrLogSink();

} // namespace Stereon

template <>
struct std::formatter<Stereon::LogLevel> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::LogLevel level, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "{}", Stereon::ToString(level));
    }
};
