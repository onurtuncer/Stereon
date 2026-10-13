// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The paged arena against std::vector on a 10^6-record table, the size of
// the ADR-0005 body, for page sizes 128, 256 and 1 024 (ADR-0005 verification,
// storage level; the topology-level run needs the stores of ADR-0004):
//
// - building the table, and reading it by index, by page span, by iterator
//   and at random positions;
// - taking a snapshot: pointer copies against a deep copy;
// - 100 local edits on a snapshot, reporting the pages copied, the fraction
//   of the body copied per edit (ADR-0005 asks < 1%) and the memory of body
//   plus edited copy relative to the body alone (ADR-0005 asks < 1.2x);
// - an undo history of 100 snapshots, one edit each, and its memory.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

import Stereon.Core;

namespace {

using Stereon::PagedArena;

/// A face-like record, about what the face table of ADR-0004 will hold.
struct FaceRecord {
    std::uint32_t Surface = 0;
    std::uint32_t FirstLoop = 0;
    std::uint32_t Shell = 0;
    std::uint32_t Flags = 0;
    double Tolerance = 0.0;
    double Area = 0.0;
};
static_assert(sizeof(FaceRecord) == 32);

template <std::size_t N>
using FaceArena = PagedArena<FaceRecord, N>;

/// The 10^6-face body of ADR-0005.
constexpr std::size_t BodyFaces = 1'000'000;
/// Local edits in the ADR-0005 verification sequence.
constexpr std::size_t EditCount = 100;
/// Random reads per iteration.
constexpr std::size_t RandomReads = 1 << 16;

FaceRecord MakeRecord(const std::size_t i)
{
    return {.Surface = static_cast<std::uint32_t>(i / 4),
            .FirstLoop = static_cast<std::uint32_t>(i),
            .Shell = 0,
            .Flags = static_cast<std::uint32_t>(i % 3),
            .Tolerance = 1e-6,
            .Area = static_cast<double>(i % 1024) * 0.25};
}

template <std::size_t N>
FaceArena<N> MakeBody()
{
    FaceArena<N> body;
    body.Reserve(BodyFaces);
    for (std::size_t i = 0; i < BodyFaces; ++i) {
        body.PushBack(MakeRecord(i));
    }
    return body;
}

std::vector<FaceRecord> MakeVectorBody()
{
    std::vector<FaceRecord> body;
    body.reserve(BodyFaces);
    for (std::size_t i = 0; i < BodyFaces; ++i) {
        body.push_back(MakeRecord(i));
    }
    return body;
}

void SetItems(benchmark::State& state, const std::size_t perIteration)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(perIteration));
}

/// Fixed-seed indices into the body, the same for every page size.
std::vector<std::size_t> RandomIndices(const std::size_t count)
{
    std::mt19937_64 rng{20261012};
    std::uniform_int_distribution<std::size_t> pick(0, BodyFaces - 1);
    std::vector<std::size_t> indices(count);
    for (std::size_t& index : indices) {
        index = pick(rng);
    }
    return indices;
}

/// Distinct pages over a set of snapshots: the memory they hold together.
template <std::size_t N>
std::size_t UniquePageCount(const std::vector<FaceArena<N>>& snapshots)
{
    std::unordered_set<const typename FaceArena<N>::Page*> pages;
    for (const FaceArena<N>& snapshot : snapshots) {
        for (std::size_t page = 0; page < snapshot.GetPageCount(); ++page) {
            pages.insert(snapshot.GetPage(page).get());
        }
    }
    return pages.size();
}

// --- Building ----------------------------------------------------------------

template <std::size_t N>
void ArenaPushBack(benchmark::State& state)
{
    for (auto _ : state) {
        FaceArena<N> body = MakeBody<N>();
        benchmark::DoNotOptimize(body);
    }
    SetItems(state, BodyFaces);
    state.counters["bytes/page"] = static_cast<double>(FaceArena<N>::BytesPerPage);
}

void VectorPushBack(benchmark::State& state)
{
    for (auto _ : state) {
        std::vector<FaceRecord> body = MakeVectorBody();
        benchmark::DoNotOptimize(body);
    }
    SetItems(state, BodyFaces);
}

// --- Reading -----------------------------------------------------------------

template <std::size_t N>
void ArenaIndexRead(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    for (auto _ : state) {
        double total = 0.0;
        for (std::size_t i = 0; i < body.GetSize(); ++i) {
            total += body[i].Area;
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, BodyFaces);
}

/// The access pattern ADR-0005 asks of hot loops: one span per page.
template <std::size_t N>
void ArenaPageSpanRead(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    for (auto _ : state) {
        double total = 0.0;
        for (std::size_t page = 0; page < body.GetPageCount(); ++page) {
            for (const FaceRecord& face : body.GetPageSpan(page)) {
                total += face.Area;
            }
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, BodyFaces);
}

template <std::size_t N>
void ArenaIteratorRead(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    for (auto _ : state) {
        double total = 0.0;
        for (const FaceRecord& face : body) {
            total += face.Area;
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, BodyFaces);
}

void VectorRead(benchmark::State& state)
{
    const std::vector<FaceRecord> body = MakeVectorBody();
    for (auto _ : state) {
        double total = 0.0;
        for (const FaceRecord& face : body) {
            total += face.Area;
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, BodyFaces);
}

template <std::size_t N>
void ArenaRandomRead(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    const std::vector<std::size_t> indices = RandomIndices(RandomReads);
    for (auto _ : state) {
        double total = 0.0;
        for (const std::size_t index : indices) {
            total += body[index].Area;
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, RandomReads);
}

void VectorRandomRead(benchmark::State& state)
{
    const std::vector<FaceRecord> body = MakeVectorBody();
    const std::vector<std::size_t> indices = RandomIndices(RandomReads);
    for (auto _ : state) {
        double total = 0.0;
        for (const std::size_t index : indices) {
            total += body[index].Area;
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state, RandomReads);
}

// --- Snapshots ---------------------------------------------------------------

template <std::size_t N>
void ArenaSnapshot(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    for (auto _ : state) {
        FaceArena<N> snapshot = body;
        benchmark::DoNotOptimize(snapshot);
    }
    state.counters["pages"] = static_cast<double>(body.GetPageCount());
}

void VectorSnapshot(benchmark::State& state)
{
    const std::vector<FaceRecord> body = MakeVectorBody();
    for (auto _ : state) {
        std::vector<FaceRecord> snapshot = body;
        benchmark::DoNotOptimize(snapshot);
    }
}

// --- Local edits (ADR-0005 verification, storage level) ----------------------

/// A snapshot of the body with EditCount single-record edits. Counters: the
/// pages the edits copied, the fraction of the body's storage that is per
/// edit, and the memory of body plus edited snapshot over the body alone.
template <std::size_t N>
void ArenaLocalEdits(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    const std::vector<std::size_t> edits = RandomIndices(EditCount);
    std::size_t copiedPages = 0;
    for (auto _ : state) {
        FaceArena<N> builder = body;
        for (const std::size_t index : edits) {
            builder.Modify(index, [](FaceRecord& face) { face.Tolerance *= 2.0; });
        }
        benchmark::DoNotOptimize(builder);
        copiedPages = builder.GetPageCount() - SharedPageCount(body, builder);
    }
    const double bodyBytes = static_cast<double>(body.GetPageBytes());
    const double copiedBytes = static_cast<double>(copiedPages * FaceArena<N>::BytesPerPage);
    state.counters["pages copied"] = static_cast<double>(copiedPages);
    state.counters["copied/edit"] = copiedBytes / bodyBytes / static_cast<double>(EditCount);
    state.counters["memory/body"] = (bodyBytes + copiedBytes) / bodyBytes;
    SetItems(state, EditCount);
}

/// The same edits on a deep copy, the cost ADR-0005's option 2 would pay.
void VectorLocalEdits(benchmark::State& state)
{
    const std::vector<FaceRecord> body = MakeVectorBody();
    const std::vector<std::size_t> edits = RandomIndices(EditCount);
    for (auto _ : state) {
        std::vector<FaceRecord> builder = body;
        for (const std::size_t index : edits) {
            builder[index].Tolerance *= 2.0;
        }
        benchmark::DoNotOptimize(builder);
    }
    state.counters["copied/edit"] = 1.0 / static_cast<double>(EditCount);
    state.counters["memory/body"] = 2.0;
    SetItems(state, EditCount);
}

/// Undo as ADR-0005 defines it: every edit keeps the previous snapshot, so
/// EditCount edits leave EditCount + 1 live arenas. The memory counter is the
/// distinct pages they hold together, relative to the body alone.
template <std::size_t N>
void ArenaEditHistory(benchmark::State& state)
{
    const FaceArena<N> body = MakeBody<N>();
    const std::vector<std::size_t> edits = RandomIndices(EditCount);
    const auto edit = [&] {
        std::vector<FaceArena<N>> history;
        history.reserve(EditCount + 1);
        history.push_back(body);
        for (const std::size_t index : edits) {
            FaceArena<N> next = history.back();
            next.Modify(index, [](FaceRecord& face) { face.Tolerance *= 2.0; });
            history.push_back(std::move(next));
        }
        return history;
    };
    for (auto _ : state) {
        std::vector<FaceArena<N>> history = edit();
        benchmark::DoNotOptimize(history);
    }
    const std::vector<FaceArena<N>> history = edit();
    state.counters["memory/body"] =
        static_cast<double>(UniquePageCount<N>(history)) / static_cast<double>(body.GetPageCount());
    SetItems(state, EditCount);
}

} // namespace

BENCHMARK(ArenaPushBack<128>);
BENCHMARK(ArenaPushBack<256>);
BENCHMARK(ArenaPushBack<1024>);
BENCHMARK(VectorPushBack);

BENCHMARK(ArenaIndexRead<128>);
BENCHMARK(ArenaIndexRead<256>);
BENCHMARK(ArenaIndexRead<1024>);
BENCHMARK(ArenaPageSpanRead<128>);
BENCHMARK(ArenaPageSpanRead<256>);
BENCHMARK(ArenaPageSpanRead<1024>);
BENCHMARK(ArenaIteratorRead<128>);
BENCHMARK(ArenaIteratorRead<256>);
BENCHMARK(ArenaIteratorRead<1024>);
BENCHMARK(VectorRead);
BENCHMARK(ArenaRandomRead<128>);
BENCHMARK(ArenaRandomRead<256>);
BENCHMARK(ArenaRandomRead<1024>);
BENCHMARK(VectorRandomRead);

BENCHMARK(ArenaSnapshot<128>);
BENCHMARK(ArenaSnapshot<256>);
BENCHMARK(ArenaSnapshot<1024>);
BENCHMARK(VectorSnapshot);

BENCHMARK(ArenaLocalEdits<128>);
BENCHMARK(ArenaLocalEdits<256>);
BENCHMARK(ArenaLocalEdits<1024>);
BENCHMARK(VectorLocalEdits);
BENCHMARK(ArenaEditHistory<128>);
BENCHMARK(ArenaEditHistory<256>);
BENCHMARK(ArenaEditHistory<1024>);
