// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cmath>
#include <cstdint>
#include <string_view>
#include <utility>

export module Stereon.Core:Context;

import :Box;
import :Enum;
import :Journal;
import :Log;
import :Vector;

/// The Context (ADR-0001, ADR-0007): everything an operation may consult
/// besides its inputs. Tolerances, the scheduler, the checking policy and
/// the diagnostic sinks live here and nowhere else; the kernel has no global
/// or thread-local state, so a result is a function of the inputs and the
/// Context alone.
///
/// A Context is a value. Its settings are fixed at construction and changed
/// by making a modified copy (`ctx.WithResolution(1e-4)`), so a `const
/// Context&` may be shared by any number of operations on any number of
/// threads. The sinks are the only things it calls out to, and they are
/// documented as thread-safe.
///
/// Settings are the programmer's, not the user's: a value outside its
/// documented domain is a contract violation (ADR-0006, ADR-0016), not a
/// KernelError. Data that arrives from outside the kernel is validated by the
/// operation that receives it, with GetSizeBox() and GetResolution() as the
/// yardsticks.

export namespace Stereon {

/// Which execution context operations run their sender graphs on (ADR-0007).
/// The serial scheduler is always available, for debugging and for the
/// determinism tests; results do not depend on the choice.
enum class Scheduler : std::uint8_t {
    Serial,
    Parallel,
};

/// The validity level checked after every operation in release builds
/// (ADR-0018). Structural is always checked; debug and sanitizer presets check
/// every level regardless of this setting.
enum class ValidityLevel : std::uint8_t {
    Structural,
    Geometric,
    Solid,
};

/// Whether an operation that had to grow a tolerance succeeds with a warning
/// or fails (ADR-0006: "any tolerance growth fails").
enum class ToleranceGrowth : std::uint8_t {
    Allow,
    Fail,
};

class Context {
public:
    /// ADR-0001: lengths in millimetres, angles in radians.
    static constexpr double DefaultResolution = 1e-5;
    static constexpr double DefaultAngularResolution = 1e-11;
    static constexpr double DefaultToleranceCeiling = 1e-2;
    /// Half the side of the size box: all geometry lies within ±500 m.
    static constexpr double SizeBoxHalfWidth = 500'000.0;
    /// ADR-0004: Compact() is due when fewer than this share of the slots of a
    /// table are alive.
    static constexpr double DefaultCompactAliveRatio = 0.5;

    /// The defaults of the ADRs, no sinks, parallel scheduler.
    Context() = default;
    [[nodiscard]] static Context Default() { return {}; }

    // --- Tolerances (ADR-0001) -----------------------------------------------

    /// The linear resolution in millimetres: closer than this, two points
    /// cannot be told apart by one comparison.
    [[nodiscard]] double GetResolution() const noexcept { return m_Resolution; }
    [[nodiscard]] double GetAngularResolution() const noexcept { return m_AngularResolution; }
    /// A local tolerance above this fails the operation with ToleranceExceeded.
    [[nodiscard]] double GetToleranceCeiling() const noexcept { return m_ToleranceCeiling; }
    [[nodiscard]] ToleranceGrowth GetToleranceGrowthPolicy() const noexcept
    {
        return m_ToleranceGrowth;
    }
    /// The box all geometry must lie in.
    [[nodiscard]] static constexpr Box3 GetSizeBox() noexcept
    {
        return Box3::FromCorners({-SizeBoxHalfWidth, -SizeBoxHalfWidth, -SizeBoxHalfWidth},
                                 {SizeBoxHalfWidth, SizeBoxHalfWidth, SizeBoxHalfWidth});
    }

    /// Raise the ceiling before the resolution when both grow.
    [[nodiscard]] Context WithResolution(const double resolution) const
        STEREON_PRE(std::isfinite(resolution) && resolution > 0.0 &&
                    resolution <= GetToleranceCeiling())
    {
        Context c = *this;
        c.m_Resolution = resolution;
        return c;
    }
    [[nodiscard]] Context WithAngularResolution(const double resolution) const
        STEREON_PRE(std::isfinite(resolution) && resolution > 0.0)
    {
        Context c = *this;
        c.m_AngularResolution = resolution;
        return c;
    }
    [[nodiscard]] Context WithToleranceCeiling(const double ceiling) const
        STEREON_PRE(std::isfinite(ceiling) && ceiling >= GetResolution())
    {
        Context c = *this;
        c.m_ToleranceCeiling = ceiling;
        return c;
    }
    [[nodiscard]] Context WithToleranceGrowthPolicy(const ToleranceGrowth policy) const
    {
        Context c = *this;
        c.m_ToleranceGrowth = policy;
        return c;
    }

    // --- Threading (ADR-0007) ------------------------------------------------

    [[nodiscard]] Scheduler GetScheduler() const noexcept { return m_Scheduler; }
    /// Worker threads for the parallel scheduler; 0 means the hardware's
    /// count. The determinism tests run the same input at 1, 2, 8 and 32.
    [[nodiscard]] unsigned GetThreadCount() const noexcept { return m_ThreadCount; }

    [[nodiscard]] Context WithScheduler(const Scheduler scheduler) const
    {
        Context c = *this;
        c.m_Scheduler = scheduler;
        return c;
    }
    [[nodiscard]] Context WithThreadCount(const unsigned count) const
    {
        Context c = *this;
        c.m_ThreadCount = count;
        return c;
    }

    // --- Checking (ADR-0018, ADR-0004) ---------------------------------------

    [[nodiscard]] ValidityLevel GetValidityLevel() const noexcept { return m_ValidityLevel; }
    [[nodiscard]] double GetCompactAliveRatio() const noexcept { return m_CompactAliveRatio; }

    [[nodiscard]] Context WithValidityLevel(const ValidityLevel level) const
    {
        Context c = *this;
        c.m_ValidityLevel = level;
        return c;
    }
    [[nodiscard]] Context WithCompactAliveRatio(const double ratio) const
        STEREON_PRE(ratio > 0.0 && ratio <= 1.0)
    {
        Context c = *this;
        c.m_CompactAliveRatio = ratio;
        return c;
    }

    // --- Logging (ADR-0023) --------------------------------------------------

    /// Records below this level are dropped before they are formatted.
    [[nodiscard]] LogLevel GetLogLevel() const noexcept { return m_LogLevel; }
    [[nodiscard]] bool HasLogSink() const noexcept { return static_cast<bool>(m_LogSink); }
    /// True when a record at `level` would reach the sink; test it before
    /// formatting an expensive message.
    [[nodiscard]] bool IsLogging(const LogLevel level) const noexcept
    {
        return m_LogSink && level >= m_LogLevel && level != LogLevel::Off;
    }

    /// Hands the message to the sink when IsLogging(level).
    void Log(const LogLevel level, const std::string_view message,
             const SourceLocation location = SourceLocation::Current()) const
    {
        if (IsLogging(level)) {
            m_LogSink(LogRecord{level, message, location});
        }
    }

    [[nodiscard]] Context WithLogLevel(const LogLevel level) const
    {
        Context c = *this;
        c.m_LogLevel = level;
        return c;
    }
    [[nodiscard]] Context WithLogSink(LogSink sink) const
    {
        Context c = *this;
        c.m_LogSink = std::move(sink);
        return c;
    }

    // --- Journal (ADR-0004, ADR-0006) ----------------------------------------

    [[nodiscard]] bool HasJournal() const noexcept { return static_cast<bool>(m_JournalSink); }

    /// Records the entry and returns its id, or the invalid id when there is
    /// no journal. Operations call this before they run.
    [[nodiscard]] JournalId Record(const JournalEntry& entry) const
    {
        return m_JournalSink ? m_JournalSink(entry) : JournalId{};
    }

    [[nodiscard]] Context WithJournalSink(JournalSink sink) const
    {
        Context c = *this;
        c.m_JournalSink = std::move(sink);
        return c;
    }

private:
    double m_Resolution = DefaultResolution;
    double m_AngularResolution = DefaultAngularResolution;
    double m_ToleranceCeiling = DefaultToleranceCeiling;
    double m_CompactAliveRatio = DefaultCompactAliveRatio;
    ToleranceGrowth m_ToleranceGrowth = ToleranceGrowth::Allow;
    Scheduler m_Scheduler = Scheduler::Parallel;
    unsigned m_ThreadCount = 0;
    ValidityLevel m_ValidityLevel = ValidityLevel::Geometric;
    LogLevel m_LogLevel = LogLevel::Info;
    LogSink m_LogSink;
    JournalSink m_JournalSink;
};

} // namespace Stereon
