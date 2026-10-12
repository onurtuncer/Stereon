// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

export module Stereon.Core:Arena;

/// Paged storage with copy-on-write pages (ADR-0005): the container behind
/// the entity tables and geometry pools of ADR-0004.
///
/// A PagedArena keeps its entries in fixed-size pages behind shared_ptr.
/// Copying an arena copies the page pointers, so a snapshot of a 10^6-entry
/// table costs a few thousand pointer copies and shares every page with its
/// source. The first write to a page makes it unique: a page that another
/// arena still holds is copied, a page this arena holds alone is written in
/// place. Memory therefore grows with what an edit touched, not with the size
/// of the table; a frozen Shape and the builder editing a copy of it share
/// all but the pages the edit wrote to.
///
/// A const arena can be read, and copied, from any thread. The mutating
/// operations belong to one thread at a time, the builder's (ADR-0004,
/// ADR-0007).
///
/// Writes are scoped. Set(index, value) replaces an entry; Modify(index, fn)
/// and ModifyPage(page, fn) hand fn a reference to the entry or a span of
/// the page, made unique first, valid only for the duration of the call.
/// No mutable reference is ever returned, because a reference that outlived
/// the copy-on-write step would write into a page that a later copy shares,
/// changing a snapshot that is meant to be immutable. During fn the arena
/// must not be copied, moved, assigned, cleared, resized or edited again;
/// each of those checks this on every build and reports a violation from a
/// module implementation unit (Arena.cpp), so the check does not depend on
/// GCC 16.2 honouring contracts in imported templates.
///
/// The page size is a tuning constant per table, fixed by the Phase 0
/// measurement of ADR-0005 (bench/ArenaBench.cpp); nothing in the public
/// API of the stores depends on it. Hot loops take the span of a page once
/// with GetPageSpan() instead of indexing entry by entry, which keeps the
/// page lookup and the reference counting out of the inner loop.

namespace Stereon::Detail {

/// Reports an arena operation attempted while an edit callback runs, then
/// raises a contract violation (and aborts without contracts). Defined in
/// Arena.cpp, a module implementation unit, so it fires in every importer.
[[noreturn]] void ArenaEditViolation(std::string_view operation) noexcept;

} // namespace Stereon::Detail

export namespace Stereon {

/// Entries per page unless a table chooses otherwise (ADR-0005).
inline constexpr std::size_t DefaultPageSize = 256;

/// A sequence of T in pages of N entries, shared between copies until written.
/// N is a power of two so that an index splits into page and offset by shift
/// and mask.
template <std::semiregular T, std::size_t N = DefaultPageSize>
    requires(N > 0 && std::has_single_bit(N))
class PagedArena {
public:
    using Page = std::array<T, N>;
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using const_reference = const T&;

    static constexpr std::size_t PageSize = N;
    static constexpr std::size_t BytesPerPage = sizeof(Page);

    /// Reads entries by index; there is no mutable iteration, since every
    /// write has to go through the copy-on-write step of Modify().
    class ConstIterator {
    public:
        using iterator_concept = std::random_access_iterator_tag;
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using reference = const T&;
        using pointer = const T*;

        constexpr ConstIterator() = default;

        [[nodiscard]] reference operator*() const
            STEREON_HOT_PRE(m_Arena != nullptr && m_Index < m_Arena->GetSize())
        {
            return (*m_Arena)[m_Index];
        }
        [[nodiscard]] pointer operator->() const { return &**this; }
        [[nodiscard]] reference operator[](const difference_type n) const { return *(*this + n); }

        ConstIterator& operator++() noexcept
        {
            ++m_Index;
            return *this;
        }
        ConstIterator operator++(int) noexcept
        {
            const ConstIterator copy = *this;
            ++m_Index;
            return copy;
        }
        ConstIterator& operator--() noexcept
        {
            --m_Index;
            return *this;
        }
        ConstIterator operator--(int) noexcept
        {
            const ConstIterator copy = *this;
            --m_Index;
            return copy;
        }
        ConstIterator& operator+=(const difference_type n) noexcept
        {
            m_Index = static_cast<std::size_t>(static_cast<difference_type>(m_Index) + n);
            return *this;
        }
        ConstIterator& operator-=(const difference_type n) noexcept { return *this += -n; }

        friend ConstIterator operator+(ConstIterator it, const difference_type n) noexcept
        {
            return it += n;
        }
        friend ConstIterator operator+(const difference_type n, ConstIterator it) noexcept
        {
            return it += n;
        }
        friend ConstIterator operator-(ConstIterator it, const difference_type n) noexcept
        {
            return it -= n;
        }
        friend difference_type operator-(const ConstIterator& a, const ConstIterator& b) noexcept
        {
            return static_cast<difference_type>(a.m_Index) -
                   static_cast<difference_type>(b.m_Index);
        }

        /// Iterators of one arena order by index.
        constexpr bool operator==(const ConstIterator&) const = default;
        constexpr std::strong_ordering operator<=>(const ConstIterator& other) const noexcept
        {
            return m_Index <=> other.m_Index;
        }

    private:
        friend class PagedArena;

        constexpr ConstIterator(const PagedArena* arena, const std::size_t index) noexcept
            : m_Arena(arena)
            , m_Index(index)
        {}

        const PagedArena* m_Arena = nullptr;
        std::size_t m_Index = 0;
    };

    using const_iterator = ConstIterator;
    using iterator = ConstIterator;

    /// Copies share every page with their source; the copy is cheap and
    /// neither arena can observe a write made through the other.
    PagedArena() = default;
    PagedArena(const PagedArena& other)
        : m_Pages((other.CheckNotEditing("copy"), other.m_Pages))
        , m_Size(other.m_Size)
    {}
    PagedArena(PagedArena&& other) noexcept
        : m_Pages((other.CheckNotEditing("move"), std::move(other.m_Pages)))
        , m_Size(std::exchange(other.m_Size, 0))
    {}
    ~PagedArena() = default;

    /// Copy-and-swap, so that an assignment copies page pointers into a
    /// fresh vector rather than over the old one. The old pages are released
    /// at the end. GCC 16.2 at -O3 also reports a false -Wnull-dereference
    /// inside std::vector's copy assignment of shared_ptr elements, which
    /// this path never takes (docs/toolchain/gcc16-feature-probe.md).
    PagedArena& operator=(PagedArena other) noexcept
    {
        CheckNotEditing("assignment");
        m_Pages.swap(other.m_Pages);
        std::swap(m_Size, other.m_Size);
        return *this;
    }

    // --- Reading: on a const arena, from any thread ----------------------

    [[nodiscard]] std::size_t GetSize() const noexcept { return m_Size; }
    [[nodiscard]] bool IsEmpty() const noexcept { return m_Size == 0; }
    [[nodiscard]] std::size_t GetPageCount() const noexcept { return m_Pages.size(); }
    /// Bytes held in this arena's pages, shared or not.
    [[nodiscard]] std::size_t GetPageBytes() const noexcept
    {
        return m_Pages.size() * BytesPerPage;
    }

    [[nodiscard]] const T& operator[](const std::size_t index) const
        STEREON_HOT_PRE(index < GetSize())
    {
        return (*m_Pages[PageOf(index)])[OffsetOf(index)];
    }

    /// The live entries of a page: PageSize of them for every page but the
    /// last. Hot loops read a page through this, once per page.
    [[nodiscard]] std::span<const T> GetPageSpan(const std::size_t page) const
        STEREON_PRE(page < GetPageCount())
    {
        return std::span<const T>{m_Pages[page]->data(), PageLength(page)};
    }

    /// The page itself, for code that holds pages (a serialiser, a diff).
    /// While the pointer lives the page counts as shared, so a write to it
    /// through this arena copies it and the holder keeps what it read.
    [[nodiscard]] std::shared_ptr<const Page> GetPage(const std::size_t page) const
        STEREON_PRE(page < GetPageCount())
    {
        return m_Pages[page];
    }

    /// True when another arena, or a GetPage() holder, holds this page too,
    /// so that the next write to it copies.
    [[nodiscard]] bool IsPageShared(const std::size_t page) const STEREON_PRE(page < GetPageCount())
    {
        return m_Pages[page].use_count() > 1;
    }

    [[nodiscard]] ConstIterator begin() const noexcept { return {this, 0}; }
    [[nodiscard]] ConstIterator end() const noexcept { return {this, m_Size}; }

    /// Equal when they hold the same entries. Pages that are the same object
    /// compare in one step, so a snapshot and its edit compare in time
    /// proportional to what the edit touched.
    [[nodiscard]] bool operator==(const PagedArena& other) const
        requires std::equality_comparable<T>
    {
        if (m_Size != other.m_Size) {
            return false;
        }
        for (std::size_t page = 0; page < m_Pages.size(); ++page) {
            if (m_Pages[page] != other.m_Pages[page] &&
                !std::ranges::equal(GetPageSpan(page), other.GetPageSpan(page))) {
                return false;
            }
        }
        return true;
    }

    /// Pages at the same index that a and b hold as the same object: what a
    /// snapshot and its edit still share.
    friend std::size_t SharedPageCount(const PagedArena& a, const PagedArena& b) noexcept
    {
        const std::size_t pages = std::min(a.m_Pages.size(), b.m_Pages.size());
        std::size_t shared = 0;
        for (std::size_t page = 0; page < pages; ++page) {
            if (a.m_Pages[page] == b.m_Pages[page]) {
                ++shared;
            }
        }
        return shared;
    }

    // --- Writing: one thread, the builder's -----------------------------

    /// Appends value and returns its index. Opens a new page when the last
    /// one is full and otherwise copies the last page first if it is shared.
    std::size_t PushBack(T value)
    {
        CheckNotEditing("PushBack");
        const std::size_t index = m_Size;
        if (OffsetOf(index) == 0) {
            m_Pages.push_back(std::make_shared<Page>());
        }
        Unshare(PageOf(index))[OffsetOf(index)] = std::move(value);
        ++m_Size;
        return index;
    }

    /// Replaces the entry, copying its page first when it is shared.
    void Set(const std::size_t index, T value) STEREON_PRE(index < GetSize())
    {
        CheckNotEditing("Set");
        Unshare(PageOf(index))[OffsetOf(index)] = std::move(value);
    }

    /// Calls fn(entry) with the entry made writable: its page is copied first
    /// when it is shared, so no other arena sees the write. The reference is
    /// valid only during the call; fn returns a value, never a reference,
    /// and must not keep the address. For many entries of one page,
    /// ModifyPage() does the copy-on-write step once.
    template <typename F>
        requires std::invocable<F&, T&>
    decltype(auto) Modify(const std::size_t index, F&& fn) STEREON_PRE(index < GetSize())
    {
        static_assert(!std::is_reference_v<std::invoke_result_t<F&, T&>>,
                      "an arena edit must not return a reference into the arena");
        CheckNotEditing("Modify");
        const EditScope scope(*this);
        return std::invoke(fn, Unshare(PageOf(index))[OffsetOf(index)]);
    }

    /// Calls fn(entries) with the live entries of a page made writable, the
    /// page copied first when it is shared; the span is valid only during
    /// the call, as for Modify().
    template <typename F>
        requires std::invocable<F&, std::span<T>>
    decltype(auto) ModifyPage(const std::size_t page, F&& fn) STEREON_PRE(page < GetPageCount())
    {
        static_assert(!std::is_reference_v<std::invoke_result_t<F&, std::span<T>>>,
                      "an arena edit must not return a reference into the arena");
        CheckNotEditing("ModifyPage");
        const EditScope scope(*this);
        return std::invoke(fn, std::span<T>{Unshare(page).data(), PageLength(page)});
    }

    /// True while a Modify() or ModifyPage() callback runs.
    [[nodiscard]] bool IsEditing() const noexcept { return m_Editing; }

    /// Room for the page pointers of size entries; pages themselves are
    /// allocated as entries arrive.
    void Reserve(const std::size_t size) { m_Pages.reserve(PagesFor(size)); }

    /// Grows with value-initialised entries, or drops the entries past size
    /// together with the pages that held only those. Shrinking never copies
    /// a page: dropped entries in a page another arena still holds stay in
    /// that page, unobservable, until the arena grows into them again.
    void Resize(const std::size_t size)
    {
        CheckNotEditing("Resize");
        if (size < m_Size) {
            Shrink(size);
        } else {
            Grow(size);
        }
    }

    /// Releases every page.
    void Clear() noexcept
    {
        CheckNotEditing("Clear");
        m_Pages.clear();
        m_Size = 0;
    }

private:
    using PagePtr = std::shared_ptr<Page>;

    /// Marks the arena as being edited for the lifetime of the scope, also
    /// when the callback exits with an exception (std::bad_alloc).
    class EditScope {
    public:
        explicit EditScope(PagedArena& arena) noexcept
            : m_Arena(arena)
        {
            m_Arena.m_Editing = true;
        }
        EditScope(const EditScope&) = delete;
        EditScope& operator=(const EditScope&) = delete;
        ~EditScope() { m_Arena.m_Editing = false; }

    private:
        PagedArena& m_Arena;
    };

    /// Every build: an operation during an edit callback would invalidate or
    /// share the reference the callback holds.
    void CheckNotEditing(const std::string_view operation) const noexcept
    {
        if (m_Editing) {
            Detail::ArenaEditViolation(operation);
        }
    }

    static constexpr std::size_t Shift = static_cast<std::size_t>(std::countr_zero(N));
    static constexpr std::size_t Mask = N - 1;

    [[nodiscard]] static constexpr std::size_t PageOf(const std::size_t index) noexcept
    {
        return index >> Shift;
    }
    [[nodiscard]] static constexpr std::size_t OffsetOf(const std::size_t index) noexcept
    {
        return index & Mask;
    }
    /// Pages that hold size entries; m_Pages.size() is always PagesFor(m_Size).
    [[nodiscard]] static constexpr std::size_t PagesFor(const std::size_t size) noexcept
    {
        return (size + Mask) >> Shift;
    }
    /// Live entries in a page.
    [[nodiscard]] std::size_t PageLength(const std::size_t page) const noexcept
    {
        return std::min(N, m_Size - page * N);
    }

    /// The page, made this arena's alone: copied when another arena or a
    /// GetPage() holder has it, otherwise returned as it is.
    [[nodiscard]] Page& Unshare(const std::size_t page)
    {
        PagePtr& ptr = m_Pages[page];
        if (ptr.use_count() > 1) {
            ptr = std::make_shared<Page>(*ptr);
        }
        return *ptr;
    }

    void Grow(const std::size_t size)
    {
        const std::size_t oldSize = m_Size;
        if (OffsetOf(oldSize) != 0) {
            // The last page may hold entries a Shrink() left behind.
            const std::size_t end = std::min(N, size - PageOf(oldSize) * N);
            std::ranges::fill(std::span<T>{Unshare(PageOf(oldSize))}.subspan(
                                  OffsetOf(oldSize), end - OffsetOf(oldSize)),
                              T{});
        }
        m_Pages.reserve(PagesFor(size));
        while (m_Pages.size() < PagesFor(size)) {
            m_Pages.push_back(std::make_shared<Page>());
        }
        m_Size = size;
    }

    void Shrink(const std::size_t size)
    {
        const std::size_t oldSize = m_Size;
        m_Pages.resize(PagesFor(size));
        m_Size = size;
        if (OffsetOf(size) != 0 && m_Pages.back().use_count() == 1) {
            // The page is this arena's alone: release what the dropped
            // entries held.
            const std::size_t end = std::min(N, oldSize - PageOf(size) * N);
            std::ranges::fill(
                std::span<T>{*m_Pages.back()}.subspan(OffsetOf(size), end - OffsetOf(size)), T{});
        }
    }

    std::vector<PagePtr> m_Pages;
    std::size_t m_Size = 0;
    bool m_Editing = false;
};

static_assert(std::random_access_iterator<PagedArena<int>::ConstIterator>);
static_assert(std::ranges::random_access_range<PagedArena<int>>);
static_assert(std::ranges::sized_range<PagedArena<int>>);

} // namespace Stereon
