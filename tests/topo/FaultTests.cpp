// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The seeded-fault suite of ADR-0004: stale handles, handles from another
// shape, off-by-one indices, wrong entity kinds and corrupted adjacency,
// each caught at the first access that returns a Result, or by Freeze,
// with a Describe() message naming the entity. The O(1) accessors
// (PositionOf, EdgeOf, ...) report the same faults as a contract violation
// in checked builds, which ends the process and so has no in-process test;
// every query and builder operation makes the same check through its
// Result, which is what this file exercises.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::FaceId;
using Stereon::LoopId;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::ShellId;
using Stereon::SolidId;
using Stereon::VertexId;

namespace {

template <typename T>
std::string Refused(const Stereon::Result<T>& result)
{
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    return std::string(result.error().GetMessage());
}

/// Every query on `shape` with handles made by `make` must be refused.
void RequireAllQueriesRefused(const Shape& shape, const std::uint32_t index)
{
    const std::string expected = "not a live entity of this Shape";
    REQUIRE(Refused(Stereon::LoopsOf(shape, FaceId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::CoedgesOf(shape, LoopId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::RadialCoedgesOf(shape, EdgeId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::FacesOf(shape, EdgeId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::EdgesOf(shape, VertexId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::FacesOf(shape, ShellId::FromIndex(index))).find(expected) !=
            std::string::npos);
    REQUIRE(Refused(Stereon::ShellsOf(shape, SolidId::FromIndex(index))).find(expected) !=
            std::string::npos);
}

} // namespace

TEST_CASE("stale builder handles are refused by every operation", "[topo][fault]")
{
    Body box = MakeBox();
    ShapeBuilder& b = box.Builder;
    const BuilderFaceId face = box.Faces[0];
    const BuilderCoedgeId coedge = box.Coedges[0];
    Must(b.RemoveFace(face));
    const std::string stale = "handle is not a live entity of this builder";
    REQUIRE(Refused(b.RemoveFace(face)).find("RemoveFace: Face " + stale) != std::string::npos);
    REQUIRE(Refused(b.AddLoop(face)).find("AddLoop: Face " + stale) != std::string::npos);
    REQUIRE(Refused(b.AddFaceToShell(box.Shells[0], face)).find(stale) != std::string::npos);
    REQUIRE(Refused(b.SetSurface(face, Stereon::SurfaceRef::FromIndex(0))).find(stale) !=
            std::string::npos);
    REQUIRE(Refused(b.RemoveCoedge(coedge)).find("Coedge " + stale) != std::string::npos);
    REQUIRE(Refused(b.SetPrimaryCoedge(box.Edges[0], coedge)).find(stale) != std::string::npos);
    REQUIRE(Refused(b.RemoveLoop(box.Loops[0])).find("Loop " + stale) != std::string::npos);
    // A handle with the right index and a wrong generation.
    const BuilderFaceId wrongGeneration =
        BuilderFaceId::Make(box.Faces[1].GetIndex(), box.Faces[1].GetGeneration() + 1, 0);
    REQUIRE(Refused(b.RemoveFace(wrongGeneration)).find(stale) != std::string::npos);
    REQUIRE(b.Count<EntityKind::Face>() == 5);
    // The frozen shape does not resolve it either.
    const Shape shape = std::move(box).Freeze();
    REQUIRE_FALSE(shape.Resolve(face).IsValid());
    REQUIRE_FALSE(shape.Resolve(coedge).IsValid());
    REQUIRE(shape.Resolve(box.Faces[1]).IsValid());
}

TEST_CASE("handles from another builder are refused", "[topo][fault]")
{
    Body a = MakeBox();
    Body b = MakeBox();
    if constexpr (STEREON_CHECKED_HANDLES) {
        // Live in both, by index; the store tells them apart.
        REQUIRE(Refused(a.Builder.AddLoop(b.Faces[0])).find("not a live entity of this builder") !=
                std::string::npos);
        REQUIRE_FALSE(a.Builder.RemoveFace(b.Faces[0]).has_value());
        REQUIRE(a.Builder.Count<EntityKind::Face>() == 6);
    }
    // An index past the other builder's table is refused in every build.
    Body wire = MakeWire(2);
    REQUIRE_FALSE(wire.Builder.RemoveFace(a.Faces[0]).has_value());
    REQUIRE_FALSE(wire.Builder
                      .AppendCoedge(a.Loops[0], wire.Edges[0], EdgeUseSense::Same,
                                    Must(wire.Builder.AddPCurve(AnyPCurve())))
                      .has_value());
}

TEST_CASE("handles from another shape are refused by every query", "[topo][fault]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    const Shape b = MakeBox().Freeze();
    if constexpr (STEREON_CHECKED_HANDLES) {
        const FaceId face = a.Resolve(box.Faces[0]);
        REQUIRE(Stereon::LoopsOf(a, face).has_value());
        REQUIRE(Refused(Stereon::LoopsOf(b, face)).find("not a live entity of this Shape") !=
                std::string::npos);
        REQUIRE_FALSE(Stereon::EdgesOf(b, a.Resolve(box.Vertices[0])).has_value());
        REQUIRE_FALSE(Stereon::FacesOf(b, a.Resolve(box.Edges[0])).has_value());
        REQUIRE_FALSE(Stereon::RadialCoedgesOf(b, a.Resolve(box.Edges[0])).has_value());
        REQUIRE_FALSE(Stereon::CoedgesOf(b, a.Resolve(box.Loops[0])).has_value());
        REQUIRE_FALSE(Stereon::FacesOf(b, a.Resolve(box.Shells[0])).has_value());
        REQUIRE_FALSE(Stereon::ShellsOf(b, a.Resolve(box.Solids[0])).has_value());
        REQUIRE(Stereon::Describe(b, face) == "Face#0 [not in this shape]");
        REQUIRE_FALSE(b.Resolve(box.Faces[0]).IsValid());
        // A handle of a shape that is gone keeps its store: no new shape
        // accepts it.
        FaceId orphan;
        {
            const Shape gone = MakeSquareSheet().Freeze();
            const auto faces = Must(Stereon::AllEntities<EntityKind::Face>(gone));
            orphan = *faces.begin();
        }
        REQUIRE_FALSE(Stereon::LoopsOf(b, orphan).has_value());
    }
    // A handle to a slot the other shape does not have is refused in every
    // build.
    const Shape wire = MakeWire(2).Freeze();
    REQUIRE_FALSE(Stereon::LoopsOf(wire, FaceId::FromIndex(0)).has_value());
    REQUIRE_FALSE(Stereon::EdgesOf(wire, VertexId::FromIndex(3)).has_value());
}

TEST_CASE("off-by-one indices and tombstones are refused", "[topo][fault]")
{
    Body box = MakeBox();
    Must(box.Builder.RemoveFace(box.Faces[3]));
    const Shape shape = std::move(box).Freeze();
    // One past the end of each table.
    RequireAllQueriesRefused(shape, 24);
    RequireAllQueriesRefused(shape, 0xFFFFFFFEu);
    // The tombstone of the removed face and its loop, as handles of this
    // very shape, and a handle of no store to the same slots.
    const Stereon::TopologySnapshot& snapshot = *shape.GetSnapshot();
    REQUIRE_FALSE(Stereon::LoopsOf(shape, snapshot.MakeId<EntityKind::Face>(3)).has_value());
    REQUIRE_FALSE(Stereon::CoedgesOf(shape, snapshot.MakeId<EntityKind::Loop>(3)).has_value());
    REQUIRE_FALSE(Stereon::LoopsOf(shape, FaceId::FromIndex(3)).has_value());
    REQUIRE_FALSE(Stereon::RadialCoedgesOf(shape, EdgeId::FromIndex(12)).has_value());
    REQUIRE(Stereon::LoopsOf(shape, snapshot.MakeId<EntityKind::Face>(2)).has_value());
    REQUIRE(Stereon::Describe(shape, FaceId::FromIndex(3)) == "Face#3 [not in this shape]");
    REQUIRE(Stereon::Describe(shape, FaceId::FromIndex(6)) == "Face#6 [not in this shape]");
    REQUIRE(Stereon::Describe(shape, FaceId{}) == "Face#- [not in this shape]");
}

namespace {

// Substitution into a concept fails quietly where a plain requires-expression
// outside a template is a hard error.
template <typename H>
concept LoopsOfTakes = requires(const Shape& s, H h) { Stereon::LoopsOf(s, h); };
template <typename H>
concept EdgesOfTakes = requires(const Shape& s, H h) { Stereon::EdgesOf(s, h); };
template <typename H>
concept FacesOfTakes = requires(const Shape& s, H h) { Stereon::FacesOf(s, h); };
template <typename H>
concept RadialCoedgesOfTakes = requires(const Shape& s, H h) { Stereon::RadialCoedgesOf(s, h); };
template <typename H>
concept RemoveFaceTakes = requires(ShapeBuilder& b, H h) { b.RemoveFace(h); };

} // namespace

TEST_CASE("mixing entity kinds or raw integers does not compile", "[topo][fault]")
{
    STATIC_REQUIRE(!std::is_convertible_v<std::uint32_t, FaceId>);
    STATIC_REQUIRE(!std::is_constructible_v<FaceId, int>);
    STATIC_REQUIRE(!std::is_convertible_v<EdgeId, FaceId>);
    STATIC_REQUIRE(!std::is_convertible_v<FaceId, std::uint32_t>);
    STATIC_REQUIRE(!std::is_convertible_v<BuilderFaceId, FaceId>);
    STATIC_REQUIRE(!std::is_convertible_v<FaceId, BuilderFaceId>);
    STATIC_REQUIRE(LoopsOfTakes<FaceId>);
    STATIC_REQUIRE(!LoopsOfTakes<EdgeId>);
    STATIC_REQUIRE(!LoopsOfTakes<std::uint32_t>);
    STATIC_REQUIRE(EdgesOfTakes<VertexId>);
    STATIC_REQUIRE(!EdgesOfTakes<FaceId>);
    STATIC_REQUIRE(FacesOfTakes<EdgeId>);
    STATIC_REQUIRE(FacesOfTakes<ShellId>);
    STATIC_REQUIRE(!FacesOfTakes<VertexId>);
    STATIC_REQUIRE(!RadialCoedgesOfTakes<LoopId>);
    STATIC_REQUIRE(RemoveFaceTakes<BuilderFaceId>);
    STATIC_REQUIRE(!RemoveFaceTakes<BuilderEdgeId>);
    STATIC_REQUIRE(!RemoveFaceTakes<FaceId>);
    STATIC_REQUIRE(!RemoveFaceTakes<std::uint32_t>);
    STATIC_REQUIRE(sizeof(FaceId) == (STEREON_CHECKED_HANDLES ? 12 : 4));
}

TEST_CASE("corrupted adjacency is refused by Freeze with the entity named", "[topo][fault]")
{
    using Stereon::Testing::TopologyFaults;
    const auto refusal = [](Body body) {
        const auto frozen = Stereon::Freeze(std::move(body.Builder));
        REQUIRE_FALSE(frozen.has_value());
        REQUIRE(frozen.error().GetCode() == ErrorCode::InvalidInput);
        return std::string(frozen.error().GetMessage());
    };
    SECTION("a loop link into another loop")
    {
        Body box = MakeBox();
        TopologyFaults::SetNext(box.Builder, box.Coedges[0], box.Coedges[4]);
        const std::string message = refusal(std::move(box));
        REQUIRE(message.find("S02 LoopLinksBroken: Coedge#0") != std::string::npos);
        REQUIRE(message.find("S03 LoopOpen: Loop#0") != std::string::npos);
    }
    SECTION("a radial link into another ring")
    {
        Body box = MakeBox();
        TopologyFaults::SetRadial(box.Builder, box.Coedges[0], box.Coedges[1]);
        const std::string message = refusal(std::move(box));
        REQUIRE(message.find("S05 RadialRingOpen: Edge#") != std::string::npos);
    }
    SECTION("a dead slot still linked")
    {
        Body box = MakeBox();
        TopologyFaults::KillCoedge(box.Builder, box.Coedges[5]);
        const std::string message = refusal(std::move(box));
        REQUIRE(message.find("S01 DanglingLink: Coedge#") != std::string::npos);
        REQUIRE(message.find("out of range or dead") != std::string::npos);
    }
    SECTION("a link past the end of a table")
    {
        Body box = MakeBox();
        TopologyFaults::SetNextOutOfRange(box.Builder, box.Coedges[2]);
        const std::string message = refusal(std::move(box));
        REQUIRE(message.find("S01 DanglingLink: Coedge#2") != std::string::npos);
    }
    SECTION("the incidence index")
    {
        Body box = MakeBox();
        TopologyFaults::SetFirstEdge(box.Builder, box.Vertices[0], BuilderEdgeId{});
        const std::string message = refusal(std::move(box));
        REQUIRE(message.find("S16 VertexIncidence: Edge#") != std::string::npos);
        REQUIRE(message.find("Vertex#0 [") != std::string::npos);
    }
    SECTION("a primary that is another edge's use")
    {
        Body cylinder = MakeCylinder();
        TopologyFaults::SetPrimary(cylinder.Builder, cylinder.Edges[2], cylinder.Coedges[0]);
        const std::string message = refusal(std::move(cylinder));
        REQUIRE(message.find("S08 PrimaryCoedge: Edge#2") != std::string::npos);
    }
}
