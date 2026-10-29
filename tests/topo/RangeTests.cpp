// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The read-only ranges of ADR-0036 to ADR-0038: the concepts they model and
// those they do not, their order and multiplicity, the snapshot they keep
// alive, and their behaviour on handles that name nothing.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::CoedgeRange;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::EntityRange;
using Stereon::ErrorCode;
using Stereon::FaceId;
using Stereon::LoopId;
using Stereon::LoopRange;
using Stereon::RadialCoedgeRange;
using Stereon::Shape;
using Stereon::ShellFaceRange;
using Stereon::SolidShellRange;
using Stereon::VertexId;

namespace {

/// What ADR-0036 asks of every local range: a forward range of handles by
/// value, not writable, not borrowed, and nothing it does not promise.
template <typename R>
constexpr bool IsTopologyRange =
    std::ranges::forward_range<R> && !std::ranges::bidirectional_range<R> &&
    !std::ranges::random_access_range<R> && !std::ranges::contiguous_range<R> &&
    !std::ranges::borrowed_range<R> && std::copyable<std::ranges::iterator_t<R>> &&
    !std::indirectly_writable<std::ranges::iterator_t<R>, std::ranges::range_value_t<R>> &&
    !std::is_reference_v<std::ranges::range_reference_t<R>>;

static_assert(IsTopologyRange<EntityRange<EntityKind::Vertex>>);
static_assert(IsTopologyRange<EntityRange<EntityKind::Solid>>);
static_assert(IsTopologyRange<LoopRange>);
static_assert(IsTopologyRange<CoedgeRange>);
static_assert(IsTopologyRange<RadialCoedgeRange>);
static_assert(IsTopologyRange<ShellFaceRange>);
static_assert(IsTopologyRange<SolidShellRange>);
static_assert(std::same_as<std::ranges::range_value_t<CoedgeRange>, CoedgeId>);
static_assert(std::same_as<std::ranges::range_value_t<LoopRange>, LoopId>);
// AllEntities knows its size: the table keeps the live count.
static_assert(std::ranges::sized_range<EntityRange<EntityKind::Edge>>);
static_assert(!std::ranges::sized_range<CoedgeRange>);
static_assert(!std::ranges::sized_range<RadialCoedgeRange>);

} // namespace

TEST_CASE("AllEntities gives every live entity once, in slot order", "[topo][ranges]")
{
    Body body = MakeBox();
    const Shape shape = std::move(body).Freeze();
    const auto vertices = Collect(Must(Stereon::AllEntities<EntityKind::Vertex>(shape)));
    REQUIRE(vertices.size() == 8);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        REQUIRE(vertices[i] == shape.Resolve(body.Vertices[i]));
        REQUIRE(vertices[i].GetIndex() == i);
    }
    const auto coedges = Must(Stereon::AllEntities<EntityKind::Coedge>(shape));
    REQUIRE(coedges.size() == 24);
    REQUIRE(static_cast<std::size_t>(std::ranges::distance(coedges)) == 24);
    REQUIRE(std::ranges::is_sorted(Collect(coedges)));
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Solid>(shape)).size() == 1);
}

TEST_CASE("AllEntities includes what no root reaches", "[topo][ranges]")
{
    // Free vertices, a wire and a face in no shell: none is below a solid.
    Body body = MakeWire(3);
    body.AddVertex({9, 9, 9});
    const BuilderFaceId f = body.AddFace(XyPlane());
    AddBoundaryPolygon(body, f, std::vector<Point3>{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});
    const Shape shape = std::move(body).Freeze();
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Vertex>(shape)).size() == 4 + 1 + 3);
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Edge>(shape)).size() == 3 + 3);
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Face>(shape)).size() == 1);
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Shell>(shape)).empty());
}

TEST_CASE("a loop walk emits its cycle from the anchor, each coedge once", "[topo][ranges]")
{
    SECTION("a triangle")
    {
        Body body = MakeFan(1);
        const Shape shape = std::move(body).Freeze();
        const auto loop = Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(body.Loops[0]))));
        REQUIRE(loop == std::vector{shape.Resolve(body.Coedges[0]), shape.Resolve(body.Coedges[1]),
                                    shape.Resolve(body.Coedges[2])});
    }
    SECTION("a one-coedge loop")
    {
        Body body = MakeCylinder();
        const Shape shape = std::move(body).Freeze();
        const auto disk = Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(body.Loops[1]))));
        REQUIRE(disk == std::vector{shape.Resolve(body.Coedges[4])});
    }
    SECTION("a long loop")
    {
        Body body = MakePrism(500);
        const Shape shape = std::move(body).Freeze();
        // The last loop is the top face: 500 coedges.
        const auto top = Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(body.Loops.back()))));
        REQUIRE(top.size() == 500);
        for (std::size_t i = 0; i < top.size(); ++i) {
            REQUIRE(EndOfUse(shape, top[i]) == StartOfUse(shape, top[(i + 1) % top.size()]));
            REQUIRE(Stereon::LoopOf(shape, top[i]) == shape.Resolve(body.Loops.back()));
        }
    }
}

TEST_CASE("a face's loops come in the order they were added", "[topo][ranges]")
{
    Body body = MakeHoleSheet();
    const Shape shape = std::move(body).Freeze();
    const auto loops = Collect(Must(Stereon::LoopsOf(shape, shape.Resolve(body.Faces[0]))));
    REQUIRE(loops == std::vector{shape.Resolve(body.Loops[0]), shape.Resolve(body.Loops[1])});
    for (const LoopId l : loops) {
        REQUIRE(Stereon::FaceOf(shape, l) == shape.Resolve(body.Faces[0]));
        REQUIRE(Collect(Must(Stereon::CoedgesOf(shape, l))).size() == 4);
    }
}

TEST_CASE("a radial walk emits every use of the edge, with multiplicity", "[topo][ranges][radial]")
{
    SECTION("one use: a sheet boundary")
    {
        Body body = MakeSquareSheet();
        const Shape shape = std::move(body).Freeze();
        for (const BuilderEdgeId e : body.Edges) {
            REQUIRE(Collect(Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(e)))).size() == 1);
        }
    }
    SECTION("two uses of opposite sense: a closed manifold")
    {
        Body body = MakeBox();
        const Shape shape = std::move(body).Freeze();
        for (const BuilderEdgeId e : body.Edges) {
            const auto ring = Collect(Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(e))));
            REQUIRE(ring.size() == 2);
            REQUIRE(Stereon::StoredEdgeUseSense(shape, ring[0]) !=
                    Stereon::StoredEdgeUseSense(shape, ring[1]));
            REQUIRE(Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[0])) !=
                    Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[1])));
        }
    }
    SECTION("k uses: a non-manifold fan")
    {
        for (const std::size_t k : {1u, 2u, 3u, 7u}) {
            Body body = MakeFan(k);
            const Shape shape = std::move(body).Freeze();
            const auto ring =
                Collect(Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(body.Edges[0]))));
            REQUIRE(ring.size() == k);
            std::vector<FaceId> faces;
            for (const CoedgeId c : ring) {
                REQUIRE(Stereon::EdgeOf(shape, c) == shape.Resolve(body.Edges[0]));
                faces.push_back(Stereon::FaceOf(shape, Stereon::LoopOf(shape, c)));
            }
            std::ranges::sort(faces);
            REQUIRE(std::ranges::adjacent_find(faces) == faces.end());
        }
    }
    SECTION("no uses: a wire edge")
    {
        Body body = MakeWire(2);
        const Shape shape = std::move(body).Freeze();
        const auto ring = Must(Stereon::RadialCoedgesOf(shape, shape.Resolve(body.Edges[0])));
        REQUIRE(ring.empty());
        REQUIRE(ring.begin() == ring.end());
    }
    SECTION("a seam: two uses on one face, each with its own pcurve")
    {
        Body body = MakeCylinder();
        const Shape shape = std::move(body).Freeze();
        const EdgeId seam = shape.Resolve(body.Edges[2]);
        REQUIRE(Stereon::KindOf(shape, seam) == EdgeCurveKind::Boundary);
        const auto ring = Collect(Must(Stereon::RadialCoedgesOf(shape, seam)));
        REQUIRE(ring.size() == 2);
        const FaceId lateral = shape.Resolve(body.Faces[0]);
        REQUIRE(Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[0])) == lateral);
        REQUIRE(Stereon::FaceOf(shape, Stereon::LoopOf(shape, ring[1])) == lateral);
        REQUIRE(Stereon::StoredEdgeUseSense(shape, ring[0]) !=
                Stereon::StoredEdgeUseSense(shape, ring[1]));
        const Point2 a = Eval(PCurveOfUse(shape, ring[0]), 1.0);
        const Point2 b = Eval(PCurveOfUse(shape, ring[1]), 1.0);
        REQUIRE(a.X == 2.0 * std::numbers::pi);
        REQUIRE(b.X == 0.0);
        REQUIRE(a.Y == b.Y);
        // The lateral loop uses the seam twice.
        const auto loop = Collect(Must(Stereon::CoedgesOf(shape, shape.Resolve(body.Loops[0]))));
        REQUIRE(std::ranges::count_if(loop, [&](const CoedgeId c) {
                    return Stereon::EdgeOf(shape, c) == seam;
                }) == 2);
    }
}

TEST_CASE("iterators are independent and traversals repeat", "[topo][ranges]")
{
    const Shape shape = MakePrism(6).Freeze();
    const auto range = Must(Stereon::AllEntities<EntityKind::Edge>(shape));
    auto a = range.begin();
    auto b = a;
    ++a;
    ++a;
    REQUIRE(*b == EdgeId::FromIndex(0));
    REQUIRE(*a == EdgeId::FromIndex(2));
    REQUIRE(*b++ == EdgeId::FromIndex(0));
    REQUIRE(*b == EdgeId::FromIndex(1));
    REQUIRE(Collect(range) == Collect(range));

    // Nested: every edge's ring inside a walk over every edge.
    std::size_t uses = 0;
    for (const EdgeId e : range) {
        for (const CoedgeId c : Must(Stereon::RadialCoedgesOf(shape, e))) {
            REQUIRE(Stereon::EdgeOf(shape, c) == e);
            ++uses;
        }
    }
    REQUIRE(uses == shape.Count<EntityKind::Coedge>());
}

TEST_CASE("a range keeps its snapshot alive after the Shape is gone", "[topo][ranges][lifetime]")
{
    std::optional<Shape> shape = MakeBox().Freeze();
    const std::weak_ptr<const Stereon::TopologySnapshot> weak = shape->GetSnapshot();
    auto edges = Must(Stereon::AllEntities<EntityKind::Edge>(*shape));
    const auto ring = Must(Stereon::RadialCoedgesOf(*shape, *edges.begin()));
    const auto expected = Collect(ring);
    // One token per range, none per iterator.
    REQUIRE(weak.use_count() == 3);
    {
        auto it = ring.begin();
        ++it;
        REQUIRE(weak.use_count() == 3);
    }
    shape.reset();
    REQUIRE_FALSE(weak.expired());
    REQUIRE(Collect(edges).size() == 12);
    REQUIRE(Collect(ring) == expected);

    // Moving a range moves its token.
    auto moved = std::move(edges);
    REQUIRE(weak.use_count() == 2);
    REQUIRE(Collect(moved).size() == 12);
}

TEST_CASE("the last owner releases the snapshot", "[topo][ranges][lifetime]")
{
    std::weak_ptr<const Stereon::TopologySnapshot> weak;
    {
        const Shape shape = MakeSquareSheet().Freeze();
        weak = shape.GetSnapshot();
        {
            // The range owning a snapshot is scoped, so the count is back
            // to one once the id is out (ADR-0046: no iterator from a temporary).
            const auto face = [&] {
                const auto faces = Must(Stereon::AllEntities<EntityKind::Face>(shape));
                return *faces.begin();
            }();
            REQUIRE(weak.use_count() == 1);
            const auto loops = Must(Stereon::LoopsOf(shape, face));
            REQUIRE(weak.use_count() == 2);
        }
        REQUIRE(weak.use_count() == 1);
    }
    REQUIRE(weak.expired());
}

TEST_CASE("a new snapshot leaves ranges over an older one unchanged", "[topo][ranges][lifetime]")
{
    const Shape first = MakePrism(4).Freeze();
    const auto old = Must(Stereon::AllEntities<EntityKind::Face>(first));
    const auto before = Collect(old);
    const Shape second = MakePrism(9).Freeze();
    REQUIRE(Collect(old) == before);
    REQUIRE(before.size() == 6);
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Face>(second)).size() == 11);
}

TEST_CASE("a query refuses a handle that names no live entity", "[topo][ranges]")
{
    const Shape shape = MakeBox().Freeze();
    const auto check = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    };
    check(Stereon::LoopsOf(shape, FaceId{}));
    check(Stereon::LoopsOf(shape, FaceId::FromIndex(6)));
    check(Stereon::CoedgesOf(shape, LoopId::FromIndex(1000)));
    check(Stereon::RadialCoedgesOf(shape, EdgeId::FromIndex(12)));
    check(Stereon::FacesOf(shape, Stereon::ShellId::FromIndex(1)));
    check(Stereon::ShellsOf(shape, Stereon::SolidId{}));
    if constexpr (STEREON_CHECKED_HANDLES) {
        // A handle of another shape, at a slot that exists in both.
        const Shape other = MakeBox().Freeze();
        const auto edges = Must(Stereon::AllEntities<EntityKind::Edge>(other));
        const auto edge = *edges.begin();
        check(Stereon::RadialCoedgesOf(shape, edge));
        REQUIRE(Stereon::Describe(shape, edge) == "Edge#0 [not in this shape]");
    }
}

TEST_CASE("handles from a range carry their snapshot in checked builds", "[topo][ranges]")
{
    const Shape shape = MakeBox().Freeze();
    for (const VertexId v : Must(Stereon::AllEntities<EntityKind::Vertex>(shape))) {
        REQUIRE(shape.GetSnapshot()->Contains(v));
        REQUIRE(v.GetStore() == (STEREON_CHECKED_HANDLES ? shape.GetSnapshot()->GetStore() : 0));
    }
}

TEST_CASE("separate ranges traverse one snapshot concurrently", "[topo][ranges][threads]")
{
    const Shape shape = MakePrism(64).Freeze();
    std::vector<CoedgeId> expected;
    for (const LoopId l : Must(Stereon::AllEntities<EntityKind::Loop>(shape))) {
        for (const CoedgeId c : Must(Stereon::CoedgesOf(shape, l))) {
            expected.push_back(c);
        }
    }
    std::vector<std::vector<CoedgeId>> seen(8);
    {
        std::vector<std::jthread> threads;
        for (std::size_t t = 0; t < seen.size(); ++t) {
            threads.emplace_back([&shape, &out = seen[t]] {
                for (int repeat = 0; repeat < 20; ++repeat) {
                    out.clear();
                    for (const LoopId l : *Stereon::AllEntities<EntityKind::Loop>(shape)) {
                        for (const CoedgeId c : *Stereon::CoedgesOf(shape, l)) {
                            out.push_back(c);
                        }
                    }
                }
            });
        }
    }
    for (const auto& out : seen) {
        REQUIRE(out == expected);
    }
}

// Each borrow as a template over the range or shape, so that a rejected call
// is a false requirement rather than an error.
template <typename R>
constexpr bool IteratorsFromLvalue = requires(const R& r) {
    r.begin();
    r.end();
};
template <typename R>
constexpr bool BeginFromRvalue = requires(R&& r) { std::move(r).begin(); };
template <typename R>
constexpr bool EndFromRvalue = requires(R&& r) { std::move(r).end(); };
template <typename S>
constexpr bool SnapshotFromLvalue = requires(const S& s) { s.GetSnapshot(); };
template <typename S>
constexpr bool SnapshotFromRvalue = requires(S&& s) { std::move(s).GetSnapshot(); };

TEST_CASE("topology ranges own their snapshot and lend iterators to lvalues only",
          "[topo][range][borrowing]")
{
    // ADR-0038 and ADR-0046, section 3 (AC-0046-02): a range is no borrowed
    // range, its iterators are taken from a named range only, and a
    // range-based for over a temporary binds the temporary first.
    using Faces = Stereon::EntityRange<EntityKind::Face>;
    STATIC_REQUIRE_FALSE(std::ranges::borrowed_range<Faces>);
    STATIC_REQUIRE_FALSE(std::ranges::borrowed_range<Stereon::LoopRange>);
    STATIC_REQUIRE(IteratorsFromLvalue<Faces>);
    STATIC_REQUIRE_FALSE(BeginFromRvalue<Faces>);
    STATIC_REQUIRE_FALSE(EndFromRvalue<Faces>);
    STATIC_REQUIRE(SnapshotFromLvalue<Shape>);
    STATIC_REQUIRE_FALSE(SnapshotFromRvalue<Shape>);

    const Shape shape = MakeSquareSheet().Freeze();
    std::size_t counted = 0;
    for ([[maybe_unused]] const Stereon::FaceId face :
         Must(Stereon::AllEntities<EntityKind::Face>(shape))) {
        ++counted;
    }
    REQUIRE(counted == Must(Stereon::AllEntities<EntityKind::Face>(shape)).size());
}

extern "C" const char* stereon_describe_handle(int kind, unsigned int index, unsigned int store);
extern "C" std::size_t stereon_describe_handle_into(int kind, unsigned int index,
                                                    unsigned int store, char* buffer,
                                                    std::size_t capacity);

TEST_CASE("a handle's description is copied into the caller's buffer", "[topo][describe]")
{
    // ADR-0046, section 7 (AC-0046-07): the registry's text leaves it only
    // as a copy made under its lock; the legacy debugger entry point keeps
    // its process-wide buffer for serialised use alone.
    if constexpr (STEREON_CHECKED_HANDLES) {
        const Shape shape = MakeBox().Freeze();
        const auto edges = Must(Stereon::AllEntities<EntityKind::Edge>(shape));
        const auto edge = *edges.begin();
        const int kind = static_cast<int>(EntityKind::Edge);
        const unsigned int index = edge.GetIndex();
        const unsigned int store = shape.GetSnapshot()->GetStore();
        const char* legacy = stereon_describe_handle(kind, index, store);
        REQUIRE(legacy != nullptr);
        const std::string expected = legacy;
        REQUIRE_FALSE(expected.empty());

        std::array<char, 256> buffer{};
        REQUIRE(stereon_describe_handle_into(kind, index, store, buffer.data(), buffer.size()) ==
                expected.size());
        REQUIRE(std::string(buffer.data()) == expected);
        std::array<char, 8> small{};
        REQUIRE(stereon_describe_handle_into(kind, index, store, small.data(), small.size()) ==
                expected.size());
        REQUIRE(std::string(small.data()) == expected.substr(0, 7));
        REQUIRE(stereon_describe_handle_into(kind, index, store, nullptr, 0) == expected.size());
        REQUIRE(stereon_describe_handle_into(kind, index, 0xFFFFFFFFu, buffer.data(),
                                             buffer.size()) == 0);
        REQUIRE(stereon_describe_handle_into(-1, index, store, buffer.data(), buffer.size()) == 0);

        // Concurrent callers each get a whole copy of their own.
        std::atomic<int> mismatches{0};
        std::vector<std::thread> threads;
        for (int t = 0; t < 8; ++t) {
            threads.emplace_back([&] {
                for (int i = 0; i < 200; ++i) {
                    std::array<char, 256> own{};
                    stereon_describe_handle_into(kind, index, store, own.data(), own.size());
                    if (std::string(own.data()) != expected) {
                        ++mismatches;
                    }
                }
            });
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
        REQUIRE(mismatches == 0);
    }
}
