// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Compact() (ADR-0004): the explicit renumbering, its remap table, and the
// compacted shape being the same shape under new indices.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::Compacted;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::FaceId;
using Stereon::LoopId;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::VertexId;

namespace {

/// A box with its face 1 removed, frozen, and the box's handles.
struct Edited {
    Body Box = MakeBox();
    Shape Frozen;
};

Edited BoxMinusFace()
{
    Edited e;
    Must(e.Box.Builder.RemoveFace(e.Box.Faces[1]));
    e.Frozen = Must(Stereon::Freeze(std::move(e.Box.Builder)));
    return e;
}

} // namespace

TEST_CASE("Compact renumbers the live entities in their old order and drops the dead",
          "[topo][compact]")
{
    const Edited e = BoxMinusFace();
    const Compacted c = Must(Stereon::Compact(e.Frozen));
    REQUIRE(c.Shape.Count<EntityKind::Face>() == 5);
    REQUIRE(c.Shape.Count<EntityKind::Loop>() == 5);
    REQUIRE(c.Shape.Count<EntityKind::Coedge>() == 20);
    REQUIRE(c.Shape.Count<EntityKind::Edge>() == 12);
    REQUIRE(c.Shape.Count<EntityKind::Vertex>() == 8);
    REQUIRE(c.SourceStore == e.Frozen.GetSnapshot()->GetStore());
    REQUIRE(c.Shape.GetSnapshot()->GetStore() != c.SourceStore);
    // Dense, in the old order.
    std::uint32_t expected = 0;
    for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(c.Shape))) {
        REQUIRE(f.GetIndex() == expected++);
    }
    REQUIRE(c.Shape.GetSnapshot()->GetTables().Faces.GetSize() == 5);
    REQUIRE(c.Shape.GetSnapshot()->GetTables().Coedges.GetSize() == 20);
    // The remap: face 1 is gone, the rest shift down by one.
    REQUIRE(c.Remap.OldCount<EntityKind::Face>() == 6);
    REQUIRE(c.Remap.NewIndex<EntityKind::Face>(0) == 0);
    REQUIRE(c.Remap.NewIndex<EntityKind::Face>(1) == FaceId::Invalid);
    REQUIRE(c.Remap.NewIndex<EntityKind::Face>(2) == 1);
    REQUIRE(c.Remap.NewIndex<EntityKind::Face>(5) == 4);
    REQUIRE(c.Remap.NewIndex<EntityKind::Face>(6) == FaceId::Invalid);
    REQUIRE(c.Remap.NewIndex<EntityKind::Vertex>(7) == 7);
    for (std::uint32_t i = 0; i < 24; ++i) {
        const std::uint32_t mapped = c.Remap.NewIndex<EntityKind::Coedge>(i);
        REQUIRE((mapped == CoedgeId::Invalid) == (i >= 4 && i < 8));
    }
    // Unreferenced geometry goes too: face 1's surface and its four pcurves.
    const Stereon::Detail::TopologyTables& t = c.Shape.GetSnapshot()->GetTables();
    REQUIRE(t.Surfaces.GetSize() == 5);
    REQUIRE(t.PCurves.GetSize() == 20);
    REQUIRE(t.Curves.GetSize() == 12);
    REQUIRE(c.Remap.NewSurface(1) == Stereon::SurfaceRef::Invalid);
    REQUIRE(c.Remap.NewSurface(2) == 1);
    REQUIRE(c.Remap.NewCurve(11) == 11);
    REQUIRE(Stereon::ValidateTopology(c.Shape).IsValid());
    REQUIRE(c.Shape.GetSnapshot()->GetCertificate().ValidatorVersion ==
            Stereon::Detail::StructuralValidatorVersion);
}

TEST_CASE("Translate takes the source's handles to the compacted shape's", "[topo][compact]")
{
    const Edited e = BoxMinusFace();
    const Compacted c = Must(Stereon::Compact(e.Frozen));
    const FaceId dead = e.Frozen.GetSnapshot()->MakeId<EntityKind::Face>(1);
    REQUIRE_FALSE(e.Frozen.GetSnapshot()->Contains(dead));
    REQUIRE_FALSE(c.Translate(dead).IsValid());
    REQUIRE_FALSE(c.Translate(FaceId{}).IsValid());
    for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(e.Frozen))) {
        const FaceId g = c.Translate(f);
        REQUIRE(g.IsValid());
        REQUIRE(g.GetIndex() == c.Remap.NewIndex<EntityKind::Face>(f.GetIndex()));
        REQUIRE(Stereon::LoopsOf(c.Shape, g).has_value());
    }
    if constexpr (STEREON_CHECKED_HANDLES) {
        const Shape other = MakeBox().Freeze();
        REQUIRE_FALSE(c.Translate(other.Resolve(e.Box.Faces[0])).IsValid());
    }
}

TEST_CASE("the compacted shape is the same shape under new indices", "[topo][compact]")
{
    const Edited e = BoxMinusFace();
    const Shape& a = e.Frozen;
    const Compacted c = Must(Stereon::Compact(a));
    const Shape& b = c.Shape;
    for (const VertexId v : Must(Stereon::AllEntities<EntityKind::Vertex>(a))) {
        const VertexId w = c.Translate(v);
        REQUIRE(PositionOf(b, w) == PositionOf(a, v));
        REQUIRE(ToleranceOf(b, w) == ToleranceOf(a, v));
        std::vector<EdgeId> translated;
        for (const EdgeId x : Must(Stereon::EdgesOf(a, v))) {
            translated.push_back(c.Translate(x));
        }
        REQUIRE(Collect(Must(Stereon::EdgesOf(b, w))) == translated);
    }
    for (const EdgeId x : Must(Stereon::AllEntities<EntityKind::Edge>(a))) {
        const EdgeId y = c.Translate(x);
        REQUIRE(KindOf(b, y) == KindOf(a, x));
        REQUIRE(RangeOf(b, y) == RangeOf(a, x));
        REQUIRE(StartOf(b, y) == c.Translate(StartOf(a, x)));
        REQUIRE(EndOf(b, y) == c.Translate(EndOf(a, x)));
        REQUIRE(Eval(Deref(CurveOf(b, y)), 0.5) == Eval(Deref(CurveOf(a, x)), 0.5));
        std::vector<FaceId> faces;
        for (const FaceId f : Must(Stereon::FacesOf(a, x))) {
            faces.push_back(c.Translate(f));
        }
        REQUIRE(Must(Stereon::FacesOf(b, y)) == faces);
        std::vector<CoedgeId> ring;
        for (const CoedgeId u : Must(Stereon::RadialCoedgesOf(a, x))) {
            ring.push_back(c.Translate(u));
        }
        REQUIRE(Collect(Must(Stereon::RadialCoedgesOf(b, y))) == ring);
    }
    for (const FaceId f : Must(Stereon::AllEntities<EntityKind::Face>(a))) {
        const FaceId g = c.Translate(f);
        REQUIRE(IsReversed(b, g) == IsReversed(a, f));
        REQUIRE(ShellOf(b, g) == c.Translate(ShellOf(a, f)));
        REQUIRE(Deref(SurfaceOf(b, g).Target<Stereon::Plane>()).Origin ==
                Deref(SurfaceOf(a, f).Target<Stereon::Plane>()).Origin);
        for (const LoopId l : Must(Stereon::LoopsOf(a, f))) {
            std::vector<CoedgeId> cycle;
            for (const CoedgeId u : Must(Stereon::CoedgesOf(a, l))) {
                cycle.push_back(c.Translate(u));
                REQUIRE(StoredEdgeUseSense(b, c.Translate(u)) == StoredEdgeUseSense(a, u));
                REQUIRE(Eval(Deref(PCurveOf(b, c.Translate(u))), 0.5) ==
                        Eval(Deref(PCurveOf(a, u)), 0.5));
            }
            REQUIRE(Collect(Must(Stereon::CoedgesOf(b, c.Translate(l)))) == cycle);
        }
    }
}

TEST_CASE("Compact of a shape without tombstones renumbers by identity", "[topo][compact]")
{
    const Shape a = MakeCylinder().Freeze();
    const Compacted c = Must(Stereon::Compact(a));
    for (std::uint32_t i = 0; i < a.Count<EntityKind::Coedge>(); ++i) {
        REQUIRE(c.Remap.NewIndex<EntityKind::Coedge>(i) == i);
    }
    REQUIRE(c.Shape.Count<EntityKind::Face>() == 3);
    REQUIRE(c.Shape.GetSnapshot()->GetTables().Surfaces.GetSize() == 3);
    REQUIRE(Stereon::ValidateTopology(c.Shape).IsValid());
    // A shape of nothing compacts to a shape of nothing.
    const Compacted nothing = Must(Stereon::Compact(Must(Stereon::Freeze(ShapeBuilder{}))));
    REQUIRE(nothing.Shape.Count<EntityKind::Vertex>() == 0);
    REQUIRE(Stereon::Compact(Shape{}).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("the source's handles do not reach the compacted shape", "[topo][compact][fault]")
{
    const Edited e = BoxMinusFace();
    const Compacted c = Must(Stereon::Compact(e.Frozen));
    // The source's last face has no slot in the compacted shape, which has
    // five; a handle made for that slot is refused in every build.
    const FaceId last = e.Frozen.GetSnapshot()->MakeId<EntityKind::Face>(5);
    REQUIRE(Stereon::LoopsOf(e.Frozen, last).has_value());
    REQUIRE_FALSE(Stereon::LoopsOf(c.Shape, last).has_value());
    REQUIRE_FALSE(Stereon::LoopsOf(c.Shape, FaceId::FromIndex(5)).has_value());
    if constexpr (STEREON_CHECKED_HANDLES) {
        const FaceId ofSource = e.Frozen.Resolve(e.Box.Faces[0]);
        REQUIRE_FALSE(Stereon::LoopsOf(c.Shape, ofSource).has_value());
        REQUIRE(Stereon::Describe(c.Shape, ofSource) == "Face#0 [not in this shape]");
        REQUIRE_FALSE(c.Shape.Resolve(e.Box.Faces[0]).IsValid());
    }
}
