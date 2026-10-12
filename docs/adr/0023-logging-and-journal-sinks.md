# ADR-0023: Logging and journal sinks

- **Status:** Proposed
- **Date:** 2026-10-12
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** ADR-0004 (replayable journals), ADR-0006 (journal ID in `KernelError`), ADR-0007 (logging lives in the `Context`), ADR-0008 (licences), ADR-0019 (ABI), ADR-0020 (bindings)

## Context

ADR-0007 puts logging in the `Context`, with no global state, and ADR-0004 asks every operation to be able to record its inputs and parameters so that a failure can be replayed. Both need a destination: a file, a console, an application's own logging framework, a Python logger. The question is whether the kernel should adopt a logging library and, if so, which.

Candidates in 2026: spdlog (MIT; mature, synchronous by default, a `std::format` mode that drops its bundled fmt), Quill and fmtlog (MIT; asynchronous, formatting on a backend thread, a few nanoseconds on the logging thread), BqLog (Apache-2.0; asynchronous), Boost.Log (BSL). None offers a C++20 module interface; all would enter the kernel through the global module fragment, as xsimd does. No standard logging facility exists in C++26.

Three constraints from existing decisions shape the answer. ADR-0019 and ADR-0020 keep third-party types out of the public API and the C ABI, so a logging library's types could not appear in `Context` without becoming part of the kernel's interface. ADR-0007 forbids global mutable state, which rules out the default-logger registries these libraries provide. And a Python binding or a host application usually has a logging framework already; the kernel's job is to deliver records to it, not to choose it.

## Options

1. **spdlog inside the kernel.** Familiar; brings fmt or `std::format`, sinks and formatting patterns. Its types in `Context` would leak through the public API, and its registry is global state the kernel would have to avoid by convention.
2. **An asynchronous library (Quill, fmtlog) inside the kernel.** The lowest per-call cost; a backend thread the kernel does not control, output reordered relative to execution, and the same interface leak as option 1.
3. **Kernel-owned sinks.** The kernel formats with `std::format` and hands records to `std::function` sinks owned by the `Context`. Applications attach whatever library they use. The kernel depends on nothing.
4. **No logging in the kernel.** Diagnostics only through `OperationResult` warnings and `KernelError`. Loses the trace of what an operation did on the way to a failure.

## Decision

Adopt option 3.

- **Sinks.** `Stereon.Core:Log` defines `LogLevel`, `LogRecord` (level, message, source location) and `LogSink = std::function<void(const LogRecord&)>`; `Stereon.Core:Journal` defines `JournalId`, `JournalEntry` (operation name and parameters as text) and `JournalSink = std::function<JournalId(const JournalEntry&)>`. Both live in the `Context`, which has no sinks by default and therefore writes nothing by itself.
- **Formatting.** Kernel code formats with `std::format`, and only after `ctx.IsLogging(level)` says the record will be delivered, so a context without a sink pays a compare per log site.
- **Threading.** Sinks are called on the thread that logs, which will be a worker thread once operations parallelise, and possibly from several threads at once. A sink is thread-safe or serialises itself; the kernel does not serialise for it. Log output is diagnostic only and never influences a result (ADR-0007).
- **Journal.** Every shape-producing operation records a `JournalEntry` before it runs and attaches the returned `JournalId` to its `KernelError` on failure (ADR-0006). The text of an entry is the operation's own serialisation of its inputs; the native format of ADR-0012 defines it for shapes. `MemoryJournal` collects entries in memory for tests and for writing the replay files that CONTRIBUTING.md asks for. Replay itself arrives with the first operations.
- **Applications.** The viewer, the tests and the Python package may use any library behind a sink. spdlog with `SPDLOG_USE_STD_FORMAT` is the recommended default where a library is wanted: it is MIT (ADR-0008), synchronous, and reuses the kernel's formatting. An asynchronous library is appropriate for a host whose profiling shows logging on a hot path; the kernel does not care. Python attaches a sink that forwards to the `logging` module (ADR-0020).
- **No third-party logging dependency in `src/`.** `tools/check_layering.py` already forbids upward dependencies; a review item forbids `spdlog`, `quill` and similar includes anywhere under `src/`.

## Consequences

- The kernel's public API and C ABI stay free of logging types; a C sink is a function pointer plus a user-data pointer (ADR-0020).
- Applications keep one logging configuration for their own code and the kernel.
- The kernel does not get spdlog's pattern formatting, rotation or filtering; the sink does, and the `Context` offers only a level threshold.
- `std::function` costs an indirect call per delivered record, after the level test. Acceptable: records are rare compared to arithmetic, and the test is a compare.
- A kernel without a sink is silent. Tests and tools attach `StderrLogSink()` when they want output.

## Verification

- `tests/core/ContextTests.cpp`: a context without a sink logs nothing and `IsLogging` is false at every level; with a sink, records below the threshold are dropped and records at or above it arrive with their message and source location; several threads recording into a `MemoryJournal` get distinct consecutive ids and lose no entry (the TSan preset runs this).
- The first shape-producing operation (the ADR-0001 prototype) records a journal entry and its failure path carries the `JournalId`; the replay of a recorded entry reproduces the failure.
- No file under `src/` includes a logging library header (review, and the layering check once it lists forbidden headers).
