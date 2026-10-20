// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <ranges>
#include <set>
#include <span>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Core;

using Stereon::DefaultPageSize;
using Stereon::PagedArena;

namespace {

/// Pages of four, so that a handful of entries already spans several pages.
using SmallArena = PagedArena<int, 4>;

template <std::size_t N>
concept ValidPageSize = requires { typename PagedArena<int, N>; };

std::vector<int> ToVector(const SmallArena& arena)
{
    return {arena.begin(), arena.end()};
}

/// 0, 1, ..., count - 1.
SmallArena Filled(const int count)
{
    SmallArena arena;
    for (int i = 0; i < count; ++i) {
        arena.PushBack(i);
    }
    return arena;
}

} // namespace

TEST_CASE("paged arenas have the layout ADR-0005 specifies", "[core][arena]")
{
    STATIC_REQUIRE(DefaultPageSize == 256);
    STATIC_REQUIRE(PagedArena<double>::PageSize == 256);
    STATIC_REQUIRE(PagedArena<double>::BytesPerPage == 256 * sizeof(double));
    STATIC_REQUIRE(SmallArena::PageSize == 4);
    STATIC_REQUIRE(SmallArena::BytesPerPage == 4 * sizeof(int));

    // The ADR-0005 sweep sizes are valid; anything but a power of two is not.
    STATIC_REQUIRE(ValidPageSize<128>);
    STATIC_REQUIRE(ValidPageSize<256>);
    STATIC_REQUIRE(ValidPageSize<1024>);
    STATIC_REQUIRE_FALSE(ValidPageSize<100>);
    STATIC_REQUIRE_FALSE(ValidPageSize<0>);

    // A read-only random-access range: writes go through Set(), Modify() and
    // ModifyPage(), which do the copy-on-write step.
    STATIC_REQUIRE(std::random_access_iterator<SmallArena::ConstIterator>);
    STATIC_REQUIRE(std::ranges::random_access_range<SmallArena>);
    STATIC_REQUIRE(std::ranges::sized_range<SmallArena>);
    STATIC_REQUIRE(std::is_same_v<std::ranges::range_value_t<SmallArena>, int>);
    STATIC_REQUIRE(std::is_same_v<std::ranges::range_reference_t<SmallArena>, const int&>);
}

TEST_CASE("entries are appended page by page and read back", "[core][arena]")
{
    SmallArena arena;
    REQUIRE(arena.IsEmpty());
    REQUIRE(arena.GetSize() == 0);
    REQUIRE(arena.GetPageCount() == 0);
    REQUIRE(arena.GetPageBytes() == 0);
    REQUIRE(arena.begin() == arena.end());

    for (int i = 0; i < 10; ++i) {
        REQUIRE(arena.PushBack(i * 10) == static_cast<std::size_t>(i));
    }
    REQUIRE_FALSE(arena.IsEmpty());
    REQUIRE(arena.GetSize() == 10);
    REQUIRE(arena.GetPageCount() == 3);
    REQUIRE(arena.GetPageBytes() == 3 * SmallArena::BytesPerPage);
    REQUIRE(arena[0] == 0);
    REQUIRE(arena[5] == 50);
    REQUIRE(arena[9] == 90);

    // Page spans cover the live entries only.
    REQUIRE(arena.GetPageSpan(0).size() == 4);
    REQUIRE(arena.GetPageSpan(1).size() == 4);
    REQUIRE(arena.GetPageSpan(2).size() == 2);
    REQUIRE(arena.GetPageSpan(1)[1] == 50);
    REQUIRE(arena.GetPageSpan(2)[1] == 90);

    // Iteration agrees with indexing, in every form the iterator offers.
    REQUIRE(ToVector(arena) == std::vector<int>{0, 10, 20, 30, 40, 50, 60, 70, 80, 90});
    REQUIRE(std::ranges::size(arena) == 10);
    REQUIRE(arena.end() - arena.begin() == 10);
    REQUIRE(*(arena.begin() + 7) == 70);
    REQUIRE(*(arena.end() - 1) == 90);
    REQUIRE(arena.begin()[3] == 30);
    REQUIRE(arena.begin() < arena.end());
    auto it = arena.begin();
    it += 4;
    REQUIRE(*it == 40);
    REQUIRE(*it++ == 40);
    REQUIRE(*it == 50);
    REQUIRE(*--it == 40);

    // Nothing shares a page with a freshly built arena.
    for (std::size_t page = 0; page < arena.GetPageCount(); ++page) {
        REQUIRE_FALSE(arena.IsPageShared(page));
    }

    arena.Clear();
    REQUIRE(arena.IsEmpty());
    REQUIRE(arena.GetPageCount() == 0);
    REQUIRE(arena.begin() == arena.end());
}

TEST_CASE("a copy shares every page and the first write to a page copies it", "[core][arena]")
{
    const SmallArena frozen = Filled(10); // pages [0..3] [4..7] [8, 9]
    SmallArena builder = frozen;

    REQUIRE(builder == frozen);
    REQUIRE(SharedPageCount(frozen, builder) == 3);
    REQUIRE(frozen.IsPageShared(1));
    REQUIRE(builder.IsPageShared(1));
    REQUIRE(builder.GetPage(1) == frozen.GetPage(1));

    // Writing through the builder copies the page it writes to, and only that
    // one; the snapshot does not change.
    builder.Set(5, 500);
    REQUIRE(builder[5] == 500);
    REQUIRE(frozen[5] == 5);
    REQUIRE(builder != frozen);
    REQUIRE(SharedPageCount(frozen, builder) == 2);
    REQUIRE_FALSE(frozen.IsPageShared(1));
    REQUIRE_FALSE(builder.IsPageShared(1));
    REQUIRE(frozen.IsPageShared(0));
    REQUIRE(builder.GetPage(0) == frozen.GetPage(0));
    REQUIRE(builder.GetPage(1) != frozen.GetPage(1));

    // A second write to the same page finds it unique and writes in place.
    const SmallArena::Page* const copied = builder.GetPage(1).get();
    builder.Set(6, 600);
    REQUIRE(builder.GetPage(1).get() == copied);
    REQUIRE(builder[6] == 600);
    REQUIRE(frozen[6] == 6);

    // Appending onto the shared last page copies that page only.
    REQUIRE(builder.PushBack(10) == 10);
    REQUIRE(builder.GetSize() == 11);
    REQUIRE(frozen.GetSize() == 10);
    REQUIRE(SharedPageCount(frozen, builder) == 1);
    REQUIRE(frozen.GetPageSpan(2).size() == 2);
    REQUIRE(builder.GetPageSpan(2).size() == 3);

    // Appending into a new page touches nothing the snapshot holds.
    builder.PushBack(11);
    builder.PushBack(12);
    REQUIRE(builder.GetPageCount() == 4);
    REQUIRE(frozen.GetPageCount() == 3);
    REQUIRE(SharedPageCount(frozen, builder) == 1);
    REQUIRE(ToVector(frozen) == std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    REQUIRE(ToVector(builder) == std::vector<int>{0, 1, 2, 3, 4, 500, 600, 7, 8, 9, 10, 11, 12});
}

TEST_CASE("a held page keeps what it read while the arena writes", "[core][arena]")
{
    SmallArena arena = Filled(4);
    const std::shared_ptr<const SmallArena::Page> held = arena.GetPage(0);
    REQUIRE(arena.IsPageShared(0));

    arena.Set(0, 42);
    REQUIRE((*held)[0] == 0);
    REQUIRE(arena[0] == 42);
    REQUIRE_FALSE(arena.IsPageShared(0));
    REQUIRE(arena.GetPage(0) != held);
}

TEST_CASE("page spans edit a page in one copy-on-write step", "[core][arena]")
{
    const SmallArena frozen = Filled(6);
    SmallArena builder = frozen;

    const std::size_t length = builder.ModifyPage(1, [](const std::span<int> page) {
        for (int& value : page) {
            value *= 100;
        }
        return page.size();
    });
    REQUIRE(length == 2);
    REQUIRE(ToVector(builder) == std::vector<int>{0, 1, 2, 3, 400, 500});
    REQUIRE(ToVector(frozen) == std::vector<int>{0, 1, 2, 3, 4, 5});
    REQUIRE(SharedPageCount(frozen, builder) == 1);

    // The whole page is the span, so a page the snapshot does not hold is
    // edited in place.
    SmallArena own = Filled(8);
    const int* const before = own.GetPageSpan(1).data();
    own.ModifyPage(1, [](const std::span<int> page) { page[0] = -4; });
    REQUIRE(own.GetPageSpan(1).data() == before);
    REQUIRE(own[4] == -4);
}

TEST_CASE("edits are scoped to a callback and never hand out a reference", "[core][arena]")
{
    const SmallArena frozen = Filled(6);
    SmallArena builder = frozen;

    // Modify returns what the callback returns, by value.
    const int old = builder.Modify(5, [](int& entry) {
        const int previous = entry;
        entry = 50;
        return previous;
    });
    REQUIRE(old == 5);
    REQUIRE(builder[5] == 50);
    REQUIRE(frozen[5] == 5);
    REQUIRE(SharedPageCount(frozen, builder) == 1);

    // The arena knows while a callback runs, and may be read meanwhile.
    REQUIRE_FALSE(builder.IsEditing());
    builder.Modify(0, [&builder](int& entry) {
        REQUIRE(builder.IsEditing());
        entry = builder[1] + builder[5]; // reads are allowed
    });
    REQUIRE_FALSE(builder.IsEditing());
    REQUIRE(builder[0] == 51);

    // A page edit returns its result the same way.
    const int sum = builder.ModifyPage(0, [](const std::span<int> page) {
        int total = 0;
        for (int& value : page) {
            value += 1;
            total += value;
        }
        return total;
    });
    REQUIRE(sum == 52 + 2 + 3 + 4);
    REQUIRE(ToVector(frozen) == std::vector<int>{0, 1, 2, 3, 4, 5});

    // Copying, moving, resizing or taking a page with GetPage() inside a
    // callback is reported and stops the program; tests/probes/ContractProbe.cpp
    // runs those cases.
}

TEST_CASE("resizing grows with value-initialised entries and drops whole pages", "[core][arena]")
{
    SmallArena arena = Filled(6);
    arena.Resize(9);
    REQUIRE(arena.GetSize() == 9);
    REQUIRE(arena.GetPageCount() == 3);
    REQUIRE(ToVector(arena) == std::vector<int>{0, 1, 2, 3, 4, 5, 0, 0, 0});

    arena.Resize(5);
    REQUIRE(arena.GetSize() == 5);
    REQUIRE(arena.GetPageCount() == 2);
    REQUIRE(arena.GetPageSpan(1).size() == 1);
    REQUIRE(ToVector(arena) == std::vector<int>{0, 1, 2, 3, 4});

    // What a shrink dropped does not come back.
    arena.Resize(7);
    REQUIRE(ToVector(arena) == std::vector<int>{0, 1, 2, 3, 4, 0, 0});
    arena.Resize(7);
    REQUIRE(ToVector(arena) == std::vector<int>{0, 1, 2, 3, 4, 0, 0});

    arena.Resize(0);
    REQUIRE(arena.IsEmpty());
    REQUIRE(arena.GetPageCount() == 0);

    // Shrinking copies nothing, even when the page is shared; growing back
    // into the shared page copies it and resets the entries.
    const SmallArena frozen = Filled(6);
    SmallArena builder = frozen;
    builder.Resize(5);
    REQUIRE(SharedPageCount(frozen, builder) == 2);
    REQUIRE(ToVector(builder) == std::vector<int>{0, 1, 2, 3, 4});
    builder.Resize(7);
    REQUIRE(ToVector(builder) == std::vector<int>{0, 1, 2, 3, 4, 0, 0});
    REQUIRE(ToVector(frozen) == std::vector<int>{0, 1, 2, 3, 4, 5});
    REQUIRE(SharedPageCount(frozen, builder) == 1);

    SmallArena reserved;
    reserved.Reserve(100);
    REQUIRE(reserved.IsEmpty());
    REQUIRE(reserved.GetPageCount() == 0);
}

TEST_CASE("copying a page shares the payloads its entries point to", "[core][arena]")
{
    using Payload = std::shared_ptr<const std::string>;
    PagedArena<Payload, 2> frozen;
    const Payload surface = std::make_shared<const std::string>("cylinder");
    frozen.PushBack(surface);
    frozen.PushBack(nullptr);
    REQUIRE(surface.use_count() == 2);

    // The copy shares the page, so the payload is not referenced again.
    PagedArena<Payload, 2> builder = frozen;
    REQUIRE(surface.use_count() == 2);

    // Copying the page copies the pointers, not what they point to.
    builder.Set(1, std::make_shared<const std::string>("plane"));
    REQUIRE(surface.use_count() == 3);
    REQUIRE(builder[0].get() == surface.get());
    REQUIRE(*builder[1] == "plane");
    REQUIRE(frozen[1] == nullptr);

    // Dropped entries release their payload when the arena holds the page
    // alone.
    builder.Clear();
    REQUIRE(surface.use_count() == 2);
    PagedArena<Payload, 2> own;
    own.PushBack(surface);
    REQUIRE(surface.use_count() == 3);
    own.Resize(0);
    REQUIRE(surface.use_count() == 2);
}

TEST_CASE("arenas compare by their entries", "[core][arena]")
{
    const SmallArena a = Filled(6);
    SmallArena b = Filled(6);
    REQUIRE(SharedPageCount(a, b) == 0);
    REQUIRE(a == b);

    b.Set(5, 55);
    REQUIRE(a != b);

    SmallArena shorter = a;
    shorter.Resize(5);
    REQUIRE(a != shorter);

    REQUIRE(SmallArena{} == SmallArena{});
    REQUIRE(SmallArena{} != a);
}

TEST_CASE("a frozen arena is read from several threads while a copy is edited",
          "[core][arena][threads]")
{
    using Arena = PagedArena<std::uint64_t, 64>;
    Arena frozen;
    for (std::uint64_t i = 0; i < 1000; ++i) {
        frozen.PushBack(i);
    }
    constexpr std::uint64_t expected = 999 * 1000 / 2;

    constexpr std::size_t readerCount = 4;
    std::vector<std::uint64_t> sums(readerCount, 0);
    std::vector<std::thread> threads;
    for (std::size_t reader = 0; reader < readerCount; ++reader) {
        threads.emplace_back([&frozen, &sums, reader] {
            for (int round = 0; round < 50; ++round) {
                const Arena snapshot = frozen; // copying is a read
                std::uint64_t sum = 0;
                for (std::size_t page = 0; page < snapshot.GetPageCount(); ++page) {
                    for (const std::uint64_t value : snapshot.GetPageSpan(page)) {
                        sum += value;
                    }
                }
                sums[reader] = sum;
            }
        });
    }
    threads.emplace_back([&frozen] {
        Arena builder = frozen;
        for (std::size_t i = 0; i < builder.GetSize(); ++i) {
            builder.Set(i, 0);
        }
        builder.PushBack(7);
        builder.Resize(10);
    });
    for (std::thread& thread : threads) {
        thread.join();
    }

    for (const std::uint64_t sum : sums) {
        REQUIRE(sum == expected);
    }
    REQUIRE(frozen.GetSize() == 1000);
    REQUIRE(frozen[999] == 999);
}

namespace {

/// A SmallArena with the std::vector it must agree with.
struct ArenaModel {
    SmallArena Arena;
    std::vector<int> Model;
};

void PushRandom(ArenaModel& m)
{
    const int value = *rc::gen::arbitrary<int>();
    m.Arena.PushBack(value);
    m.Model.push_back(value);
}

void ModifyRandom(ArenaModel& m)
{
    if (m.Model.empty()) {
        return;
    }
    const auto index = *rc::gen::inRange<std::size_t>(0, m.Model.size());
    const int value = *rc::gen::arbitrary<int>();
    m.Arena.Set(index, value);
    m.Model[index] = value;
}

void BumpRandomPage(ArenaModel& m)
{
    if (m.Model.empty()) {
        return;
    }
    const auto page = *rc::gen::inRange<std::size_t>(0, m.Arena.GetPageCount());
    m.Arena.ModifyPage(page, [](const std::span<int> entries) {
        for (int& value : entries) {
            ++value;
        }
    });
    const std::size_t first = page * SmallArena::PageSize;
    const std::size_t last = std::min(first + SmallArena::PageSize, m.Model.size());
    for (std::size_t i = first; i < last; ++i) {
        ++m.Model[i];
    }
}

void ResizeRandom(ArenaModel& m)
{
    const auto size = *rc::gen::inRange<std::size_t>(0, m.Model.size() + 6);
    m.Arena.Resize(size);
    m.Model.resize(size);
}

void RandomStep(ArenaModel& m)
{
    switch (*rc::gen::inRange(0, 5)) {
    case 0:
    case 1:
        PushRandom(m);
        break;
    case 2:
        ModifyRandom(m);
        break;
    case 3:
        BumpRandomPage(m);
        break;
    default:
        ResizeRandom(m);
        break;
    }
}

} // namespace

TEST_CASE("paged arenas behave as vectors with immutable snapshots", "[core][arena][property]")
{
    rc::prop("every operation agrees with std::vector and never changes a snapshot", [] {
        ArenaModel m;
        std::vector<ArenaModel> snapshots;
        const int steps = *rc::gen::inRange(1, 80);
        for (int step = 0; step < steps; ++step) {
            if (*rc::gen::inRange(0, 4) == 0) {
                snapshots.push_back(m);
            }
            RandomStep(m);

            RC_ASSERT(m.Arena.GetSize() == m.Model.size());
            RC_ASSERT(m.Arena.GetPageCount() == (m.Model.size() + 3) / 4);
            RC_ASSERT(ToVector(m.Arena) == m.Model);
            for (const ArenaModel& snapshot : snapshots) {
                RC_ASSERT(ToVector(snapshot.Arena) == snapshot.Model);
            }
        }
    });

    rc::prop("an edit copies only the pages it touches", [] {
        const int count = *rc::gen::inRange(1, 200);
        const SmallArena frozen = Filled(count);
        SmallArena builder = frozen;
        const int edits = *rc::gen::inRange(0, 20);
        std::set<std::size_t> touched;
        for (int edit = 0; edit < edits; ++edit) {
            const auto index = *rc::gen::inRange<std::size_t>(0, static_cast<std::size_t>(count));
            builder.Set(index, -1);
            touched.insert(index / SmallArena::PageSize);
        }

        RC_ASSERT(SharedPageCount(frozen, builder) == frozen.GetPageCount() - touched.size());
        for (std::size_t page = 0; page < frozen.GetPageCount(); ++page) {
            RC_ASSERT(frozen.IsPageShared(page) == !touched.contains(page));
        }
        for (int i = 0; i < count; ++i) {
            RC_ASSERT(frozen[static_cast<std::size_t>(i)] == i);
        }
    });

    rc::prop("a snapshot outlives the arena it was taken from", [] {
        const int count = *rc::gen::inRange(0, 50);
        SmallArena snapshot;
        {
            SmallArena source = Filled(count);
            snapshot = source;
            source.Clear();
        }
        RC_ASSERT(snapshot.GetSize() == static_cast<std::size_t>(count));
        for (int i = 0; i < count; ++i) {
            RC_ASSERT(snapshot[static_cast<std::size_t>(i)] == i);
        }
    });
}
