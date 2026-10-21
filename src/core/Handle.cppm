// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <compare>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <string_view>
#include <vector>

export module Stereon.Core:Handle;

import :Enum;
import :Scalar;

/// Typed handles into the topology stores (ADR-0004).
///
/// A frozen Shape addresses its entities with Id<Kind>: a 32-bit index wrapped
/// in one type per entity kind, so that mixing kinds or passing a raw integer
/// is a compile error. ShapeBuilder uses BuilderId<Kind>, an index with a
/// generation, so a handle to a removed entity is stale and detected.
/// Oriented<Kind> packs a sense bit into a handle, for coedges, reachable
/// only through BaseOf() and IsReversed().
///
/// With STEREON_CHECKED_HANDLES every handle also records the store it came
/// from and its generation, and GenerationTable::Check() verifies store,
/// bounds and generation; a mismatch prints a report naming the handle, the
/// entity kind and the expected store, then raises a contract violation.
/// Release builds compile Id down to its 4-byte index.
///
/// Handles are local to one store. Identity across operations comes from
/// history maps (ADR-0011), never from handles.

export namespace Stereon {

/// The entity kinds of the logical model (ADR-0004). The C ABI exposes the
/// numbers (ADR-0020), so values are never changed, only appended.
enum class EntityKind : std::uint8_t {
    Vertex = 0,
    Edge,
    Coedge,
    Loop,
    Face,
    Shell,
    Solid,
};

inline constexpr std::size_t EntityKindCount = 7;

/// "Vertex", "Edge", ...; "Entity" where the compiler has no reflection.
[[nodiscard]] constexpr std::string_view ToString(const EntityKind kind) noexcept
{
    const std::string_view name = EnumName(kind);
    return name.empty() ? std::string_view{"Entity"} : name;
}

/// Identifies one store (a Shape or a ShapeBuilder) in checked builds. 0 is
/// "no store": the value of handles made without one, and of every handle in
/// release builds.
using StoreId = std::uint32_t;

/// Handle to an entity of a frozen Shape: its slot index in the store of its
/// kind. Stores make handles; other code only receives and compares them.
template <EntityKind K>
class Id {
public:
    static constexpr EntityKind Kind = K;
    /// The index of the default-constructed, invalid handle.
    static constexpr std::uint32_t Invalid = 0xFFFFFFFFu;
    /// The largest index Oriented<K> can carry, since it keeps the top bit for
    /// the sense and the all-ones pattern for its own invalid value.
    static constexpr std::uint32_t MaxOrientable = 0x7FFFFFFEu;

    constexpr Id() = default;

    /// The handle to slot `index` of `store`, at `generation`. Release builds
    /// keep only the index.
    [[nodiscard]] static constexpr Id Make(const std::uint32_t index,
                                           [[maybe_unused]] const StoreId store,
                                           [[maybe_unused]] const std::uint32_t generation) noexcept
    {
        Id id;
        id.m_Index = index;
#if STEREON_CHECKED_HANDLES
        id.m_Store = store;
        id.m_Generation = generation;
#endif
        return id;
    }

    /// A handle with no store, for tests and for readers of serialised data.
    [[nodiscard]] static constexpr Id FromIndex(const std::uint32_t index) noexcept
    {
        return Make(index, 0, 0);
    }

    [[nodiscard]] constexpr std::uint32_t GetIndex() const noexcept { return m_Index; }
    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Index != Invalid; }

    /// 0 in release builds.
    [[nodiscard]] constexpr StoreId GetStore() const noexcept
    {
#if STEREON_CHECKED_HANDLES
        return m_Store;
#else
        return 0;
#endif
    }
    /// 0 in release builds.
    [[nodiscard]] constexpr std::uint32_t GetGeneration() const noexcept
    {
#if STEREON_CHECKED_HANDLES
        return m_Generation;
#else
        return 0;
#endif
    }

    /// Handles are equal when they name the same slot. The checked-build
    /// payload is diagnostic and never changes a comparison, so debug and
    /// release builds take the same decisions (ADR-0007).
    constexpr bool operator==(const Id& other) const noexcept { return m_Index == other.m_Index; }
    constexpr std::strong_ordering operator<=>(const Id& other) const noexcept
    {
        return m_Index <=> other.m_Index;
    }

private:
    std::uint32_t m_Index = Invalid;
#if STEREON_CHECKED_HANDLES
    StoreId m_Store = 0;
    std::uint32_t m_Generation = 0;
#endif
};

using VertexId = Id<EntityKind::Vertex>;
using EdgeId = Id<EntityKind::Edge>;
using CoedgeId = Id<EntityKind::Coedge>;
using LoopId = Id<EntityKind::Loop>;
using FaceId = Id<EntityKind::Face>;
using ShellId = Id<EntityKind::Shell>;
using SolidId = Id<EntityKind::Solid>;

static_assert(sizeof(VertexId) == (STEREON_CHECKED_HANDLES ? 12 : 4),
              "release-build handles are exactly 4 bytes (ADR-0004)");

/// Handle to an entity of a ShapeBuilder: a slot index with the generation
/// the slot had when the handle was made. Releasing the slot bumps its
/// generation, so every handle to the released entity is stale.
template <EntityKind K>
class BuilderId {
public:
    static constexpr EntityKind Kind = K;
    static constexpr std::uint32_t Invalid = Id<K>::Invalid;

    constexpr BuilderId() = default;

    [[nodiscard]] static constexpr BuilderId Make(const std::uint32_t index,
                                                  const std::uint32_t generation,
                                                  [[maybe_unused]] const StoreId store) noexcept
    {
        BuilderId id;
        id.m_Index = index;
        id.m_Generation = generation;
#if STEREON_CHECKED_HANDLES
        id.m_Store = store;
#endif
        return id;
    }

    [[nodiscard]] constexpr std::uint32_t GetIndex() const noexcept { return m_Index; }
    [[nodiscard]] constexpr std::uint32_t GetGeneration() const noexcept { return m_Generation; }
    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Index != Invalid; }
    /// 0 in release builds.
    [[nodiscard]] constexpr StoreId GetStore() const noexcept
    {
#if STEREON_CHECKED_HANDLES
        return m_Store;
#else
        return 0;
#endif
    }

    /// The generation is part of the identity: a stale handle differs from
    /// the live handle to the same slot.
    constexpr bool operator==(const BuilderId& other) const noexcept
    {
        return m_Index == other.m_Index && m_Generation == other.m_Generation;
    }
    constexpr std::strong_ordering operator<=>(const BuilderId& other) const noexcept
    {
        if (const auto byIndex = m_Index <=> other.m_Index; byIndex != 0) {
            return byIndex;
        }
        return m_Generation <=> other.m_Generation;
    }

private:
    std::uint32_t m_Index = Invalid;
    std::uint32_t m_Generation = 0;
#if STEREON_CHECKED_HANDLES
    StoreId m_Store = 0;
#endif
};

static_assert(sizeof(BuilderId<EntityKind::Vertex>) == (STEREON_CHECKED_HANDLES ? 12 : 8));

/// A handle with a sense: the entity used forwards or reversed, as a loop
/// uses a coedge. The sense lives in the top bit of the index and is reached
/// only through the functions below; nothing else masks it (ADR-0004).
template <EntityKind K>
class Oriented {
public:
    static constexpr EntityKind Kind = K;

    constexpr Oriented() = default;

    /// `base` must be orientable: its index at most Id<K>::MaxOrientable, or
    /// invalid, which gives the invalid oriented handle.
    [[nodiscard]] static constexpr Oriented Make(const Id<K> base, const bool reversed) noexcept
        STEREON_HOT_PRE(!base.IsValid() || base.GetIndex() <= Id<K>::MaxOrientable)
    {
        Oriented h;
        if (base.IsValid()) {
            h.m_Bits = base.GetIndex() | (reversed ? SenseBit : 0u);
        }
#if STEREON_CHECKED_HANDLES
        h.m_Store = base.GetStore();
        h.m_Generation = base.GetGeneration();
#endif
        return h;
    }

    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Bits != Id<K>::Invalid; }

    /// The handle without its sense.
    friend constexpr Id<K> BaseOf(const Oriented h) noexcept
    {
        if (!h.IsValid()) {
            return {};
        }
        return Id<K>::Make(h.m_Bits & ~SenseBit, h.GetStore(), h.GetGeneration());
    }

    friend constexpr bool IsReversed(const Oriented h) noexcept
    {
        return h.IsValid() && (h.m_Bits & SenseBit) != 0;
    }

    /// The same entity with the opposite sense.
    friend constexpr Oriented Reversed(Oriented h) noexcept
    {
        if (h.IsValid()) {
            h.m_Bits ^= SenseBit;
        }
        return h;
    }

    /// 0 in release builds.
    [[nodiscard]] constexpr StoreId GetStore() const noexcept
    {
#if STEREON_CHECKED_HANDLES
        return m_Store;
#else
        return 0;
#endif
    }
    /// 0 in release builds.
    [[nodiscard]] constexpr std::uint32_t GetGeneration() const noexcept
    {
#if STEREON_CHECKED_HANDLES
        return m_Generation;
#else
        return 0;
#endif
    }

    /// Equal when they name the same slot with the same sense.
    constexpr bool operator==(const Oriented& other) const noexcept
    {
        return m_Bits == other.m_Bits;
    }
    constexpr std::strong_ordering operator<=>(const Oriented& other) const noexcept
    {
        return m_Bits <=> other.m_Bits;
    }

private:
    static constexpr std::uint32_t SenseBit = 0x80000000u;

    std::uint32_t m_Bits = Id<K>::Invalid;
#if STEREON_CHECKED_HANDLES
    StoreId m_Store = 0;
    std::uint32_t m_Generation = 0;
#endif
};

using OrientedCoedgeId = Oriented<EntityKind::Coedge>;

static_assert(sizeof(OrientedCoedgeId) == sizeof(CoedgeId));

} // namespace Stereon

// Exported for the topology store, whose frozen tables check handles the way
// GenerationTable does and draw their store IDs from the same counter.
export namespace Stereon::Detail {

/// What a failed handle check knows; HandleCheckFailed() turns it into the
/// report of ADR-0004.
struct HandleReport {
    EntityKind Kind{};
    std::uint32_t Index = 0;
    std::uint32_t Generation = 0;
    StoreId Store = 0;
    StoreId ExpectedStore = 0;
    std::uint32_t Size = 0;
    std::uint32_t SlotGeneration = 0;
    bool SlotAlive = false;
};

/// Prints the report and raises a contract violation. Defined in Handle.cpp,
/// a module implementation unit, so that the contract is compiled in the
/// module's own unit and fires in every importer at every optimisation level
/// (GCC 16.2; see Contracts.hpp).
void HandleCheckFailed(const HandleReport& report);

/// A fresh identity: of a store in checked builds, of a builder, of a
/// ParameterSet. The counter is the only global state in the kernel; it
/// tells identities apart and never influences a result's value (ADR-0007).
[[nodiscard]] StoreId NextStoreId() noexcept;

} // namespace Stereon::Detail

export namespace Stereon {

/// Slot bookkeeping for the entities of one kind in a ShapeBuilder
/// (ADR-0004): a generation per slot, the alive bitset and the free list.
/// Releasing a slot bumps its generation, so every BuilderId to the released
/// entity is stale. Slots released in a session are not reused in that
/// session, so an index live before and after an edit names the same entity;
/// BeginSession() makes them reusable, which a builder opened on an existing
/// Shape calls first.
template <EntityKind K>
class GenerationTable {
public:
    using Handle = BuilderId<K>;

    GenerationTable()
    {
#if STEREON_CHECKED_HANDLES
        m_Store = Detail::NextStoreId();
#endif
    }

    /// A live handle to a reusable slot, or to a new one.
    [[nodiscard]] Handle Allocate()
    {
        std::uint32_t index;
        if (m_Free.empty()) {
            index = static_cast<std::uint32_t>(m_Generations.size());
            m_Generations.push_back(0);
            m_Alive.push_back(true);
        } else {
            index = m_Free.back();
            m_Free.pop_back();
            m_Alive[index] = true;
        }
        ++m_AliveCount;
        return Handle::Make(index, m_Generations[index], m_Store);
    }

    /// Marks the slot dead and makes every handle to it stale. A handle that
    /// is not live is reported (checked builds) and otherwise ignored.
    void Release(const Handle h)
    {
        if (!IsAlive(h)) {
            Check(h);
            return;
        }
        const std::uint32_t index = h.GetIndex();
        m_Alive[index] = false;
        ++m_Generations[index];
        m_Pending.push_back(index);
        --m_AliveCount;
    }

    /// Slots released before this call become reusable.
    void BeginSession()
    {
        m_Free.insert(m_Free.end(), m_Pending.begin(), m_Pending.end());
        m_Pending.clear();
    }

    /// True when h names a live slot of this table. Asks; never reports.
    [[nodiscard]] bool IsAlive(const Handle h) const noexcept
    {
#if STEREON_CHECKED_HANDLES
        if (h.GetStore() != m_Store) {
            return false;
        }
#endif
        const std::uint32_t index = h.GetIndex();
        return index < m_Generations.size() && m_Alive[index] &&
               m_Generations[index] == h.GetGeneration();
    }

    /// Verifies store, bounds and generation in checked builds and reports a
    /// mismatch as a contract violation; does nothing in release builds.
    void Check([[maybe_unused]] const Handle h) const
    {
#if STEREON_CHECKED_HANDLES
        if (IsAlive(h)) {
            return;
        }
        Detail::HandleReport report{
            .Kind = K,
            .Index = h.GetIndex(),
            .Generation = h.GetGeneration(),
            .Store = h.GetStore(),
            .ExpectedStore = m_Store,
            .Size = static_cast<std::uint32_t>(m_Generations.size()),
        };
        if (h.GetIndex() < m_Generations.size()) {
            report.SlotGeneration = m_Generations[h.GetIndex()];
            report.SlotAlive = m_Alive[h.GetIndex()];
        }
        Detail::HandleCheckFailed(report);
#endif
    }

    /// Slots, live or dead.
    [[nodiscard]] std::uint32_t GetSize() const noexcept
    {
        return static_cast<std::uint32_t>(m_Generations.size());
    }
    [[nodiscard]] std::uint32_t GetAliveCount() const noexcept { return m_AliveCount; }
    [[nodiscard]] bool IsAliveAt(const std::uint32_t index) const noexcept
    {
        return index < m_Alive.size() && m_Alive[index];
    }
    [[nodiscard]] std::uint32_t GetGeneration(const std::uint32_t index) const
        STEREON_HOT_PRE(index < GetSize())
    {
        return m_Generations[index];
    }
    /// 0 in release builds.
    [[nodiscard]] StoreId GetStore() const noexcept { return m_Store; }

private:
    std::vector<std::uint32_t> m_Generations;
    std::vector<bool> m_Alive;
    std::vector<std::uint32_t> m_Free;
    std::vector<std::uint32_t> m_Pending;
    std::uint32_t m_AliveCount = 0;
    StoreId m_Store = 0;
};

} // namespace Stereon

// std::format support: "Vertex#12" (invalid: "Vertex#-"), "Vertex#12@3" for a
// builder handle at generation 3, "Coedge#12 (reversed)" for an oriented one.

namespace Stereon::Detail {

template <typename FormatContext>
auto FormatHandle(const EntityKind kind, const std::uint32_t index, const bool valid,
                  FormatContext& ctx)
{
    if (!valid) {
        return std::format_to(ctx.out(), "{}#-", ToString(kind));
    }
    return std::format_to(ctx.out(), "{}#{}", ToString(kind), index);
}

} // namespace Stereon::Detail

template <Stereon::EntityKind K>
struct std::formatter<Stereon::Id<K>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::Id<K>& h, FormatContext& ctx) const
    {
        return Stereon::Detail::FormatHandle(K, h.GetIndex(), h.IsValid(), ctx);
    }
};

template <Stereon::EntityKind K>
struct std::formatter<Stereon::BuilderId<K>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BuilderId<K>& h, FormatContext& ctx) const
    {
        auto out = Stereon::Detail::FormatHandle(K, h.GetIndex(), h.IsValid(), ctx);
        if (h.IsValid()) {
            out = std::format_to(out, "@{}", h.GetGeneration());
        }
        return out;
    }
};

template <Stereon::EntityKind K>
struct std::formatter<Stereon::Oriented<K>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::Oriented<K>& h, FormatContext& ctx) const
    {
        const Stereon::Id<K> base = BaseOf(h);
        auto out = Stereon::Detail::FormatHandle(K, base.GetIndex(), base.IsValid(), ctx);
        if (IsReversed(h)) {
            out = std::format_to(out, " (reversed)");
        }
        return out;
    }
};

// std::hash, for the side tables keyed by handle (ADR-0004). Consistent with
// operator==: Id hashes its index only.

template <Stereon::EntityKind K>
struct std::hash<Stereon::Id<K>> {
    std::size_t operator()(const Stereon::Id<K>& h) const noexcept
    {
        return std::hash<std::uint32_t>{}(h.GetIndex());
    }
};

template <Stereon::EntityKind K>
struct std::hash<Stereon::BuilderId<K>> {
    std::size_t operator()(const Stereon::BuilderId<K>& h) const noexcept
    {
        const std::uint64_t bits = (std::uint64_t{h.GetIndex()} << 32) | h.GetGeneration();
        return std::hash<std::uint64_t>{}(bits);
    }
};

template <Stereon::EntityKind K>
struct std::hash<Stereon::Oriented<K>> {
    std::size_t operator()(const Stereon::Oriented<K>& h) const noexcept
    {
        const Stereon::Id<K> base = BaseOf(h);
        return std::hash<std::uint32_t>{}(base.GetIndex() | (IsReversed(h) ? 0x80000000u : 0u));
    }
};
