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
#include <set>
#include <type_traits>
#include <unordered_set>
#include <vector>

import Stereon.Core;

#include "Formatting.hpp"

using Stereon::BuilderId;
using Stereon::CheckedHandles;
using Stereon::CoedgeId;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::GenerationTable;
using Stereon::OrientedCoedgeId;
using Stereon::VertexId;

using VertexBuilderId = BuilderId<EntityKind::Vertex>;

TEST_CASE("handles have the sizes ADR-0004 specifies", "[core][handles]")
{
    STATIC_REQUIRE(sizeof(VertexId) == (CheckedHandles ? 12u : 4u));
    STATIC_REQUIRE(sizeof(OrientedCoedgeId) == sizeof(CoedgeId));
    STATIC_REQUIRE(sizeof(VertexBuilderId) == (CheckedHandles ? 12u : 8u));
    STATIC_REQUIRE(std::is_trivially_copyable_v<VertexId>);
    STATIC_REQUIRE(std::is_trivially_copyable_v<VertexBuilderId>);
    STATIC_REQUIRE(std::is_trivially_copyable_v<OrientedCoedgeId>);
}

TEST_CASE("entity kinds are distinct types with no integer conversions", "[core][handles]")
{
    STATIC_REQUIRE_FALSE((std::is_convertible_v<VertexId, EdgeId>));
    STATIC_REQUIRE_FALSE((std::is_constructible_v<VertexId, EdgeId>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<std::uint32_t, VertexId>));
    STATIC_REQUIRE_FALSE((std::is_constructible_v<VertexId, int>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<VertexId, std::uint32_t>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<VertexId, bool>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<VertexBuilderId, VertexId>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<OrientedCoedgeId, CoedgeId>));

    STATIC_REQUIRE(VertexId::Kind == EntityKind::Vertex);
    STATIC_REQUIRE(EdgeId::Kind == EntityKind::Edge);
    // The C ABI exposes these numbers (ADR-0020).
    STATIC_REQUIRE(static_cast<int>(EntityKind::Vertex) == 0);
    STATIC_REQUIRE(static_cast<int>(EntityKind::Solid) == 6);
    STATIC_REQUIRE(Stereon::EntityKindCount == 7);
}

TEST_CASE("handles compare by slot and format by kind", "[core][handles]")
{
    constexpr VertexId a = VertexId::FromIndex(12);
    constexpr VertexId b = VertexId::FromIndex(12);
    constexpr VertexId c = VertexId::FromIndex(13);
    STATIC_REQUIRE(a == b);
    STATIC_REQUIRE(a != c);
    STATIC_REQUIRE(a < c);
    STATIC_REQUIRE(a.GetIndex() == 12);
    STATIC_REQUIRE(a.IsValid());
    STATIC_REQUIRE_FALSE(VertexId{}.IsValid());
    STATIC_REQUIRE(VertexId{}.GetIndex() == VertexId::Invalid);

    // The checked-build payload never takes part in a comparison.
    constexpr VertexId checked = VertexId::Make(12, 7, 3);
    REQUIRE(checked == a);
    if constexpr (CheckedHandles) {
        REQUIRE(checked.GetStore() == 7);
        REQUIRE(checked.GetGeneration() == 3);
    } else {
        REQUIRE(checked.GetStore() == 0);
        REQUIRE(checked.GetGeneration() == 0);
    }

    const std::set<VertexId> ordered{c, a, b};
    REQUIRE(ordered.size() == 2);
    REQUIRE(*ordered.begin() == a);
    const std::unordered_set<VertexId> hashed{a, b, c, checked};
    REQUIRE(hashed.size() == 2);

    if constexpr (Stereon::EnumNames) {
        REQUIRE(std::format("{}", a) == "Vertex#12");
        REQUIRE(std::format("{}", EdgeId::FromIndex(4)) == "Edge#4");
        REQUIRE(std::format("{}", VertexId{}) == "Vertex#-");
        REQUIRE(Stereon::ToString(EntityKind::Coedge) == "Coedge");
    } else {
        REQUIRE(std::format("{}", a) == "Entity#12");
    }
}

TEST_CASE("oriented handles expose the sense only through named functions", "[core][handles]")
{
    constexpr CoedgeId base = CoedgeId::FromIndex(42);
    constexpr OrientedCoedgeId forward = OrientedCoedgeId::Make(base, false);
    constexpr OrientedCoedgeId reversed = OrientedCoedgeId::Make(base, true);

    STATIC_REQUIRE(BaseOf(forward) == base);
    STATIC_REQUIRE(BaseOf(reversed) == base);
    STATIC_REQUIRE_FALSE(IsReversed(forward));
    STATIC_REQUIRE(IsReversed(reversed));
    STATIC_REQUIRE(forward != reversed);
    STATIC_REQUIRE(Reversed(forward) == reversed);
    STATIC_REQUIRE(Reversed(Reversed(forward)) == forward);
    STATIC_REQUIRE(forward.IsValid());

    STATIC_REQUIRE_FALSE(OrientedCoedgeId{}.IsValid());
    STATIC_REQUIRE_FALSE(BaseOf(OrientedCoedgeId{}).IsValid());
    STATIC_REQUIRE_FALSE(IsReversed(OrientedCoedgeId{}));
    STATIC_REQUIRE(Reversed(OrientedCoedgeId{}) == OrientedCoedgeId{});
    STATIC_REQUIRE_FALSE(OrientedCoedgeId::Make(CoedgeId{}, true).IsValid());

    // The largest orientable index survives the round trip with either sense.
    constexpr CoedgeId largest = CoedgeId::FromIndex(CoedgeId::MaxOrientable);
    STATIC_REQUIRE(BaseOf(OrientedCoedgeId::Make(largest, true)) == largest);
    STATIC_REQUIRE(OrientedCoedgeId::Make(largest, true).IsValid());
    STATIC_REQUIRE(IsReversed(OrientedCoedgeId::Make(largest, true)));

    const std::unordered_set<OrientedCoedgeId> hashed{forward, reversed, Reversed(reversed)};
    REQUIRE(hashed.size() == 2);

    if constexpr (Stereon::EnumNames) {
        REQUIRE(std::format("{}", forward) == "Coedge#42");
        REQUIRE(std::format("{}", reversed) == "Coedge#42 (reversed)");
        REQUIRE(std::format("{}", OrientedCoedgeId{}) == "Coedge#-");
    }
}

TEST_CASE("generation tables detect stale handles", "[core][handles]")
{
    GenerationTable<EntityKind::Vertex> table;
    REQUIRE(table.GetSize() == 0);

    const VertexBuilderId a = table.Allocate();
    const VertexBuilderId b = table.Allocate();
    REQUIRE(a != b);
    REQUIRE(a.GetIndex() == 0);
    REQUIRE(b.GetIndex() == 1);
    REQUIRE(a.GetGeneration() == 0);
    REQUIRE(table.GetSize() == 2);
    REQUIRE(table.GetAliveCount() == 2);
    REQUIRE(table.IsAlive(a));
    REQUIRE(table.IsAlive(b));

    table.Release(a);
    REQUIRE_FALSE(table.IsAlive(a));
    REQUIRE(table.IsAlive(b));
    REQUIRE(table.GetAliveCount() == 1);
    REQUIRE(table.GetGeneration(0) == 1);
    REQUIRE_FALSE(table.IsAliveAt(0));

    // Within the session the released slot is not reused.
    const VertexBuilderId c = table.Allocate();
    REQUIRE(c.GetIndex() == 2);
    REQUIRE(table.GetSize() == 3);

    // After BeginSession() it is, at the next generation, and the old handle
    // stays stale.
    table.BeginSession();
    const VertexBuilderId d = table.Allocate();
    REQUIRE(d.GetIndex() == 0);
    REQUIRE(d.GetGeneration() == 1);
    REQUIRE(d != a);
    REQUIRE(table.IsAlive(d));
    REQUIRE_FALSE(table.IsAlive(a));
    REQUIRE(table.GetSize() == 3);
    REQUIRE(table.GetAliveCount() == 3);

    // Releasing a stale handle changes nothing (and is reported in checked
    // builds; see the contract probe).
    REQUIRE_FALSE(table.IsAlive(VertexBuilderId{}));
    REQUIRE_FALSE(table.IsAlive(VertexBuilderId::Make(99, 0, table.GetStore())));

    GenerationTable<EntityKind::Vertex> other;
    const VertexBuilderId foreign = other.Allocate();
    REQUIRE(other.IsAlive(foreign));
    if constexpr (CheckedHandles) {
        REQUIRE(other.GetStore() != table.GetStore());
        REQUIRE(other.GetStore() != 0);
        REQUIRE_FALSE(table.IsAlive(foreign));
    }

    if constexpr (Stereon::EnumNames) {
        REQUIRE(std::format("{}", d) == "Vertex#0@1");
        REQUIRE(std::format("{}", VertexBuilderId{}) == "Vertex#-");
    }
    const std::unordered_set<VertexBuilderId> hashed{a, b, c, d};
    REQUIRE(hashed.size() == 4);
}

TEST_CASE("generation tables keep their invariants under random use", "[core][handles][property]")
{
    rc::prop("live handles stay live, released handles stay stale", [] {
        GenerationTable<EntityKind::Edge> table;
        std::vector<BuilderId<EntityKind::Edge>> live;
        std::vector<BuilderId<EntityKind::Edge>> stale;
        const auto steps = *rc::gen::inRange(1, 60);

        for (int step = 0; step < steps; ++step) {
            const int action = *rc::gen::inRange(0, 10);
            if (action < 5 || live.empty()) {
                live.push_back(table.Allocate());
            } else if (action < 9) {
                const auto which = *rc::gen::inRange<std::size_t>(0, live.size());
                stale.push_back(live[which]);
                table.Release(live[which]);
                live.erase(live.begin() + static_cast<std::ptrdiff_t>(which));
            } else {
                table.BeginSession();
            }

            RC_ASSERT(table.GetAliveCount() == live.size());
            for (const auto& h : live) {
                RC_ASSERT(table.IsAlive(h));
            }
            for (const auto& h : stale) {
                RC_ASSERT(!table.IsAlive(h));
            }
        }

        // Every live handle names a distinct slot, and no slot exceeds the
        // table.
        std::set<std::uint32_t> indices;
        for (const auto& h : live) {
            RC_ASSERT(h.GetIndex() < table.GetSize());
            RC_ASSERT(indices.insert(h.GetIndex()).second);
        }
    });

    rc::prop("a new session reuses released slots before growing", [] {
        GenerationTable<EntityKind::Face> table;
        const auto count = *rc::gen::inRange(1, 20);
        std::vector<BuilderId<EntityKind::Face>> handles;
        for (int i = 0; i < count; ++i) {
            handles.push_back(table.Allocate());
        }
        for (const auto& h : handles) {
            table.Release(h);
        }
        RC_ASSERT(table.GetAliveCount() == 0);
        table.BeginSession();
        for (int i = 0; i < count; ++i) {
            const auto reused = table.Allocate();
            RC_ASSERT(reused.GetIndex() < static_cast<std::uint32_t>(count));
            RC_ASSERT(reused.GetGeneration() == 1);
        }
        RC_ASSERT(table.GetSize() == static_cast<std::uint32_t>(count));
    });
}
