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
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Core:Journal;

import :Scalar;

/// The operation journal (ADR-0004 safeguard 6, ADR-0006, ADR-0023). Before
/// an operation runs it records its name and its inputs and parameters, as
/// text, with the Context's JournalSink; the sink answers with a JournalId
/// that a KernelError carries, so that a failure reported from the field can
/// be replayed in isolation. The kernel owns the record format and the sink
/// interface; where the entries go is the application's choice.
///
/// Replay itself arrives with the first operations: the text of an entry is
/// the operation's own serialisation of its inputs, which the native format
/// of ADR-0012 will define for shapes.

export namespace Stereon {

/// The number a sink gave an entry. 0 means the operation was not
/// journaled, because the Context had no sink.
struct JournalId {
    std::uint64_t Value = 0;

    [[nodiscard]] constexpr bool IsValid() const noexcept { return Value != 0; }
    constexpr bool operator==(const JournalId&) const = default;
};

/// What an operation records. The views are valid for the duration of the
/// sink call only; a sink that keeps an entry copies it.
struct JournalEntry {
    /// The operation's name, e.g. "Boolean.Union".
    std::string_view Operation;
    /// Inputs and parameters in the operation's text form.
    std::string_view Parameters;
};

/// Where journal entries go; returns the id the entry got. Thread safety as
/// for LogSink: a sink may be called from several threads at once.
using JournalSink = std::function<JournalId(const JournalEntry&)>;

/// A journal kept in memory, entries in order of recording and ids counting
/// from 1. For tests, and for writing the replay files of CONTRIBUTING.md.
/// Safe to record into from several threads.
class MemoryJournal {
public:
    struct Entry {
        JournalId Id;
        std::string Operation;
        std::string Parameters;
    };

    MemoryJournal() = default;
    MemoryJournal(const MemoryJournal&) = delete;
    MemoryJournal& operator=(const MemoryJournal&) = delete;

    JournalId Record(const JournalEntry& entry)
    {
        const std::scoped_lock lock(m_Mutex);
        const JournalId id{static_cast<std::uint64_t>(m_Entries.size()) + 1};
        m_Entries.push_back({id, std::string(entry.Operation), std::string(entry.Parameters)});
        return id;
    }

    [[nodiscard]] std::size_t GetSize() const
    {
        const std::scoped_lock lock(m_Mutex);
        return m_Entries.size();
    }

    /// A copy of the entries recorded so far.
    [[nodiscard]] std::vector<Entry> GetEntries() const
    {
        const std::scoped_lock lock(m_Mutex);
        return m_Entries;
    }

    /// One line per entry: "#<id> <operation> <parameters>".
    [[nodiscard]] std::string ToText() const
    {
        std::string text;
        for (const Entry& entry : GetEntries()) {
            text += std::format("#{} {} {}\n", entry.Id.Value, entry.Operation, entry.Parameters);
        }
        return text;
    }

private:
    mutable std::mutex m_Mutex;
    std::vector<Entry> m_Entries;
};

/// A sink that records into `journal`, which it keeps alive.
[[nodiscard]] inline JournalSink SinkFor(std::shared_ptr<MemoryJournal> journal)
{
    return [journal = std::move(journal)](const JournalEntry& entry) {
        return journal->Record(entry);
    };
}

} // namespace Stereon

// std::format support: "J#12", or "J#-" for the invalid id.

template <>
struct std::formatter<Stereon::JournalId> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::JournalId id, FormatContext& ctx) const
    {
        if (!id.IsValid()) {
            return std::format_to(ctx.out(), "J#-");
        }
        return std::format_to(ctx.out(), "J#{}", id.Value);
    }
};
