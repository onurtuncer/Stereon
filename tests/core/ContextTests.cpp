// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

import Stereon.Core;

#include "Formatting.hpp"

using Stereon::Context;
using Stereon::JournalEntry;
using Stereon::JournalId;
using Stereon::LogLevel;
using Stereon::LogRecord;
using Stereon::MemoryJournal;
using Stereon::Scheduler;
using Stereon::ToleranceGrowth;
using Stereon::ValidityLevel;

TEST_CASE("the default context carries the ADR-0001 and ADR-0007 defaults", "[core][context]")
{
    const Context ctx = Context::Default();
    REQUIRE(ctx.GetResolution() == 1e-5);
    REQUIRE(ctx.GetAngularResolution() == 1e-11);
    REQUIRE(ctx.GetToleranceCeiling() == 1e-2);
    REQUIRE(ctx.GetToleranceGrowthPolicy() == ToleranceGrowth::Allow);
    REQUIRE(ctx.GetScheduler() == Scheduler::Parallel);
    REQUIRE(ctx.GetThreadCount() == 0);
    REQUIRE(ctx.GetValidityLevel() == ValidityLevel::Geometric);
    REQUIRE(ctx.GetCompactAliveRatio() == 0.5);
    REQUIRE(ctx.GetLogLevel() == LogLevel::Info);
    REQUIRE_FALSE(ctx.HasLogSink());
    REQUIRE_FALSE(ctx.HasJournal());

    // The size box is ±500 m, closed.
    STATIC_REQUIRE(Context::SizeBoxHalfWidth == 500'000.0);
    constexpr Stereon::Box3 box = Context::GetSizeBox();
    STATIC_REQUIRE(Contains(box, Stereon::Point3{500'000.0, -500'000.0, 0.0}));
    STATIC_REQUIRE_FALSE(Contains(box, Stereon::Point3{500'000.5, 0.0, 0.0}));
    STATIC_REQUIRE(Contains(box, Stereon::Point3{}));

    STATIC_REQUIRE(std::is_copy_constructible_v<Context>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<Context>);
}

TEST_CASE("settings change by making a modified copy", "[core][context]")
{
    const Context base = Context::Default();
    const Context loose = base.WithToleranceCeiling(0.1)
                              .WithResolution(1e-4)
                              .WithAngularResolution(1e-9)
                              .WithToleranceGrowthPolicy(ToleranceGrowth::Fail)
                              .WithScheduler(Scheduler::Serial)
                              .WithThreadCount(8)
                              .WithValidityLevel(ValidityLevel::Solid)
                              .WithCompactAliveRatio(0.25)
                              .WithLogLevel(LogLevel::Debug);

    REQUIRE(loose.GetToleranceCeiling() == 0.1);
    REQUIRE(loose.GetResolution() == 1e-4);
    REQUIRE(loose.GetAngularResolution() == 1e-9);
    REQUIRE(loose.GetToleranceGrowthPolicy() == ToleranceGrowth::Fail);
    REQUIRE(loose.GetScheduler() == Scheduler::Serial);
    REQUIRE(loose.GetThreadCount() == 8);
    REQUIRE(loose.GetValidityLevel() == ValidityLevel::Solid);
    REQUIRE(loose.GetCompactAliveRatio() == 0.25);
    REQUIRE(loose.GetLogLevel() == LogLevel::Debug);

    // The original is a value and did not change.
    REQUIRE(base.GetToleranceCeiling() == 1e-2);
    REQUIRE(base.GetResolution() == 1e-5);
    REQUIRE(base.GetScheduler() == Scheduler::Parallel);
    REQUIRE(base.GetLogLevel() == LogLevel::Info);
}

TEST_CASE("log records reach the sink at or above the threshold with their location",
          "[core][context][log]")
{
    std::vector<std::pair<LogLevel, std::string>> received;
    std::uint_least32_t line = 0;
    const Context ctx = Context::Default().WithLogSink([&](const LogRecord& record) {
        received.emplace_back(record.Level, std::string(record.Message));
        line = record.Location.Line;
    });

    REQUIRE(ctx.HasLogSink());
    REQUIRE(ctx.IsLogging(LogLevel::Info));
    REQUIRE(ctx.IsLogging(LogLevel::Error));
    REQUIRE_FALSE(ctx.IsLogging(LogLevel::Debug));
    REQUIRE_FALSE(ctx.IsLogging(LogLevel::Trace));
    REQUIRE_FALSE(ctx.IsLogging(LogLevel::Off));

    ctx.Log(LogLevel::Debug, "dropped");
    ctx.Log(LogLevel::Info, std::format("resolution {}", ctx.GetResolution()));
    const auto expectedLine = static_cast<std::uint_least32_t>(__LINE__ + 1);
    ctx.Log(LogLevel::Warning, "kept");
    REQUIRE(received.size() == 2);
    REQUIRE(received[0].first == LogLevel::Info);
    REQUIRE(received[0].second == "resolution 1e-05");
    REQUIRE(received[1].first == LogLevel::Warning);
    REQUIRE(received[1].second == "kept");
    REQUIRE(line == expectedLine);

    // The threshold is part of the context; Off drops everything.
    const Context verbose = ctx.WithLogLevel(LogLevel::Trace);
    REQUIRE(verbose.IsLogging(LogLevel::Trace));
    verbose.Log(LogLevel::Trace, "now kept");
    REQUIRE(received.size() == 3);
    const Context silent = ctx.WithLogLevel(LogLevel::Off);
    REQUIRE_FALSE(silent.IsLogging(LogLevel::Error));
    silent.Log(LogLevel::Error, "never");
    REQUIRE(received.size() == 3);

    // Without a sink nothing is logged and nothing is formatted.
    const Context none = Context::Default();
    REQUIRE_FALSE(none.IsLogging(LogLevel::Error));
    none.Log(LogLevel::Error, "nowhere");

    if constexpr (Stereon::EnumNames) {
        REQUIRE(std::format("{}", LogLevel::Warning) == "Warning");
        REQUIRE(Stereon::ToString(LogLevel::Off) == "Off");
    }
    REQUIRE(static_cast<bool>(Stereon::StderrLogSink()));
}

TEST_CASE("journal entries are recorded in order and carry ids into errors",
          "[core][context][journal]")
{
    const Context none = Context::Default();
    REQUIRE_FALSE(none.Record(JournalEntry{"Test.Noop", ""}).IsValid());

    const auto journal = std::make_shared<MemoryJournal>();
    const Context ctx = Context::Default().WithJournalSink(Stereon::SinkFor(journal));
    REQUIRE(ctx.HasJournal());

    const JournalId first = ctx.Record(JournalEntry{"Boolean.Union", "a=box b=cylinder"});
    const JournalId second = ctx.Record(JournalEntry{"Boolean.Cut", "a=box b=sphere"});
    REQUIRE(first == JournalId{1});
    REQUIRE(second == JournalId{2});
    REQUIRE(first.IsValid());
    REQUIRE_FALSE(JournalId{}.IsValid());
    REQUIRE(journal->GetSize() == 2);

    const std::vector<MemoryJournal::Entry> entries = journal->GetEntries();
    REQUIRE(entries.size() == 2);
    REQUIRE(entries[0].Id == first);
    REQUIRE(entries[0].Operation == "Boolean.Union");
    REQUIRE(entries[0].Parameters == "a=box b=cylinder");
    REQUIRE(entries[1].Operation == "Boolean.Cut");
    REQUIRE(journal->ToText() ==
            "#1 Boolean.Union a=box b=cylinder\n#2 Boolean.Cut a=box b=sphere\n");
    REQUIRE(std::format("{}", second) == "J#2");
    REQUIRE(std::format("{}", JournalId{}) == "J#-");

    // The id travels with the error of the operation it belongs to.
    const Stereon::KernelError error =
        Stereon::KernelError(Stereon::ErrorCode::IntersectionFailed, "no curve")
            .WithJournalId(second);
    REQUIRE(error.GetJournalId() == second);
    REQUIRE_FALSE(
        Stereon::KernelError(Stereon::ErrorCode::InvalidInput, "x").GetJournalId().IsValid());

    // The journal outlives the context that recorded into it.
    {
        const Context scoped = Context::Default().WithJournalSink(Stereon::SinkFor(journal));
        REQUIRE(scoped.Record(JournalEntry{"Test.Scoped", "-"}) == JournalId{3});
    }
    REQUIRE(journal->GetSize() == 3);
}

TEST_CASE("several threads record into one journal without losing entries",
          "[core][context][journal][threads]")
{
    const auto journal = std::make_shared<MemoryJournal>();
    const Context ctx = Context::Default().WithJournalSink(Stereon::SinkFor(journal));
    constexpr std::size_t threadCount = 8;
    constexpr std::size_t perThread = 200;

    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < threadCount; ++t) {
        threads.emplace_back([&ctx, t] {
            const std::string name = std::format("Thread{}", t);
            for (std::size_t i = 0; i < perThread; ++i) {
                const JournalId id = ctx.Record(JournalEntry{name, "-"});
                if (!id.IsValid()) {
                    return;
                }
            }
        });
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    REQUIRE(journal->GetSize() == threadCount * perThread);
    std::set<std::uint64_t> ids;
    for (const MemoryJournal::Entry& entry : journal->GetEntries()) {
        REQUIRE(ids.insert(entry.Id.Value).second);
    }
    REQUIRE(*ids.begin() == 1);
    REQUIRE(*ids.rbegin() == threadCount * perThread);
}

TEST_CASE("the log threshold filters exactly by level", "[core][context][log][property]")
{
    rc::prop("IsLogging is level >= threshold, except that Off drops everything", [] {
        const auto threshold = static_cast<LogLevel>(*rc::gen::inRange(0, 6));
        const auto level = static_cast<LogLevel>(*rc::gen::inRange(0, 6));
        int delivered = 0;
        const Context ctx = Context::Default().WithLogLevel(threshold).WithLogSink(
            [&](const LogRecord&) { ++delivered; });
        const bool expected = level >= threshold && level != LogLevel::Off;
        RC_ASSERT(ctx.IsLogging(level) == expected);
        ctx.Log(level, "message");
        RC_ASSERT(delivered == (expected ? 1 : 0));
    });
}
