// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Imprinting (the Boolean prototype's M6): a clipped branch written into a
// face through the Euler operators, inside a transaction, and the result
// committed valid. On a box: a line through two edges (two splits, one face
// split), through two corners (no splits), a circle inside the top face (a
// ring and a new disc face), a circle through all four edges (four splits
// of the face). On a cylinder: an ellipse wholly inside the lateral face
// as an Exact ring, and the skew cylinders' quartic as a Trim ring with its
// cached curve.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::Box;
using BooleanTest::BranchOn;
using BooleanTest::CountOf;
using BooleanTest::ctx;
using BooleanTest::Drill;
using BooleanTest::faceBox;
using BooleanTest::PlaneFaceAlong;
using BooleanTest::VerticalPlane;
using Stereon::Box2;
using Stereon::BranchSide;
using Stereon::EdgeCurveKind;
using Stereon::EntityKind;
using Stereon::FaceClip;
using Stereon::FaceId;
using Stereon::Imprint;
using Stereon::IntersectionBranch;
using Stereon::OperationResult;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

/// The branch clipped to the face and imprinted into a builder opened on
/// the shape, inside a transaction, then committed.
struct Imprinted {
    Imprint Made;
    OperationResult Result;
};

Imprinted ImprintOn(const Shape& shape, const FaceId face, const IntersectionBranch& branch,
                    const BranchSide side)
{
    const auto clip = ClipToFace(shape, face, branch, side, ctx);
    REQUIRE(clip.has_value());
    auto builder = ShapeBuilder::Open(shape);
    REQUIRE(builder.has_value());
    const auto made = builder->Transaction([&](ShapeBuilder& candidate) {
        return ImprintClip(candidate, shape, face, branch, side, *clip, ctx);
    });
    if (!made) {
        INFO(std::string(made.error().GetMessage()));
    }
    REQUIRE(made.has_value());
    auto committed = Commit(std::move(*builder), ctx);
    if (!committed) {
        INFO(std::string(committed.error().GetMessage()));
    }
    REQUIRE(committed.has_value());
    return {*made, std::move(*committed)};
}

} // namespace

TEST_CASE("a line through a box face splits two edges and the face", "[boolean][imprint]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const IntersectionBranch branch = BranchOn(box, top, VerticalPlane(3.0), faceBox);
    const Imprinted done = ImprintOn(box, top, branch, BranchSide::A);
    REQUIRE(done.Made.SplitVertices.size() == 2);
    REQUIRE(done.Made.Pieces.size() == 1);
    REQUIRE(done.Made.Pieces.front().Face.IsValid());
    const Shape& shape = done.Result.Shape;
    REQUIRE(CountOf<EntityKind::Vertex>(shape) == 10);
    REQUIRE(CountOf<EntityKind::Edge>(shape) == 15);
    REQUIRE(CountOf<EntityKind::Face>(shape) == 7);
    REQUIRE(CountOf<EntityKind::Loop>(shape) == 7);
    REQUIRE_FALSE(done.Result.TopologyEvents.empty());
    // The new edge is the exact line, with the branch's tolerance.
    const Stereon::EdgeId edge = shape.Resolve(done.Made.Pieces.front().Edge);
    REQUIRE(edge.IsValid());
    REQUIRE(KindOf(shape, edge) == EdgeCurveKind::Exact);
    REQUIRE(CurveOf(shape, edge).has_value());
    REQUIRE(CurveOf(shape, edge)->Target<Stereon::Line>() != nullptr);
    for (const Stereon::VertexId v : {StartOf(shape, edge), EndOf(shape, edge)}) {
        const Point3 p = PositionOf(shape, v);
        REQUIRE(std::abs(p.X - 3.0) <= 1e-12);
        REQUIRE(std::abs(p.Z - 10.0) <= 1e-12);
    }
}

TEST_CASE("a line through two corners splits the face without splitting an edge",
          "[boolean][imprint]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const Surface diagonal{Plane{Point3{}, Normalized(Vec3{1.0, 1.0, 0.0}), Vec3{0.0, 0.0, 1.0}}};
    const IntersectionBranch branch = BranchOn(box, top, diagonal, faceBox);
    const Imprinted done = ImprintOn(box, top, branch, BranchSide::A);
    REQUIRE(done.Made.SplitVertices.empty());
    REQUIRE(done.Made.Pieces.size() == 1);
    const Shape& shape = done.Result.Shape;
    REQUIRE(CountOf<EntityKind::Vertex>(shape) == 8);
    REQUIRE(CountOf<EntityKind::Edge>(shape) == 13);
    REQUIRE(CountOf<EntityKind::Face>(shape) == 7);
}

TEST_CASE("a circle inside a face becomes a ring and a new face", "[boolean][imprint]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const IntersectionBranch branch = BranchOn(box, top, Drill(2.0), faceBox);
    const Imprinted done = ImprintOn(box, top, branch, BranchSide::A);
    REQUIRE(done.Made.Pieces.size() == 1);
    REQUIRE(done.Made.Pieces.front().Start == done.Made.Pieces.front().End);
    REQUIRE(done.Made.Pieces.front().Face.IsValid());
    const Shape& shape = done.Result.Shape;
    REQUIRE(CountOf<EntityKind::Vertex>(shape) == 9);
    REQUIRE(CountOf<EntityKind::Edge>(shape) == 13);
    REQUIRE(CountOf<EntityKind::Face>(shape) == 7);
    REQUIRE(CountOf<EntityKind::Loop>(shape) == 8);
    // The top face keeps its outer loop and gains the ring; the disc has one.
    const auto loops = LoopsOf(shape, shape.Resolve(done.Made.Pieces.front().Face));
    REQUIRE(loops.has_value());
    std::size_t discLoops = 0;
    for ([[maybe_unused]] const Stereon::LoopId loop : *loops) {
        ++discLoops;
    }
    REQUIRE(discLoops == 1);
    const Stereon::EdgeId edge = shape.Resolve(done.Made.Pieces.front().Edge);
    REQUIRE(KindOf(shape, edge) == EdgeCurveKind::Exact);
    REQUIRE(StartOf(shape, edge) == EndOf(shape, edge));
}

TEST_CASE("a circle through all four edges splits the face four times", "[boolean][imprint]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const IntersectionBranch branch = BranchOn(box, top, Drill(6.0), faceBox);
    const Imprinted done = ImprintOn(box, top, branch, BranchSide::A);
    REQUIRE(done.Made.SplitVertices.size() == 8);
    REQUIRE(done.Made.Pieces.size() == 4);
    const Shape& shape = done.Result.Shape;
    REQUIRE(CountOf<EntityKind::Vertex>(shape) == 16);
    REQUIRE(CountOf<EntityKind::Edge>(shape) == 24);
    REQUIRE(CountOf<EntityKind::Face>(shape) == 10);
    for (const Stereon::ImprintedPiece& piece : done.Made.Pieces) {
        REQUIRE(piece.Face.IsValid());
        const Stereon::EdgeId edge = shape.Resolve(piece.Edge);
        REQUIRE(CurveOf(shape, edge)->Target<Stereon::Circle>() != nullptr);
    }
}

TEST_CASE("an ellipse inside a cylinder's lateral face is an exact ring, a quartic a cached one",
          "[boolean][imprint]")
{
    const Shape drum = BooleanTest::Drum();
    const FaceId lateral = BooleanTest::CylinderFace(drum);
    // The plane first: the cylinder is side B.
    const Surface oblique{
        Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Normalized(Vec3{0.0, 1.0, -1.0})}};
    const auto ellipse =
        Intersect(oblique, Box2::Everything(), SurfaceOf(drum, lateral), Box2::Everything(), ctx);
    REQUIRE(ellipse.has_value());
    const Imprinted exact = ImprintOn(drum, lateral, ellipse->Branches.front(), BranchSide::B);
    REQUIRE(exact.Made.Pieces.size() == 1);
    const Shape& cut = exact.Result.Shape;
    REQUIRE(CountOf<EntityKind::Face>(cut) == 4);
    REQUIRE(CountOf<EntityKind::Edge>(cut) == 4);
    REQUIRE(CountOf<EntityKind::Loop>(cut) == 5);
    const Stereon::EdgeId ring = cut.Resolve(exact.Made.Pieces.front().Edge);
    REQUIRE(KindOf(cut, ring) == EdgeCurveKind::Exact);
    REQUIRE(CurveOf(cut, ring)->Target<Stereon::Ellipse>() != nullptr);

    // A cylinder of radius 2 along x, its axis 4 from the drum's: one loop
    // of the quartic, a B-spline, between heights 3 and 7.
    const Surface across{Stereon::Cylinder{Point3{0.0, 4.0, 5.0}, 2.0, Vec3{0.0, 1.0, 0.0},
                                           Vec3{0.0, 0.0, 1.0}, true}};
    const auto quartic =
        Intersect(SurfaceOf(drum, lateral), Box2::Everything(), across, Box2::Everything(), ctx);
    REQUIRE(quartic.has_value());
    REQUIRE(quartic->Branches.size() == 1);
    const Imprinted trim = ImprintOn(drum, lateral, quartic->Branches.front(), BranchSide::A);
    REQUIRE(trim.Made.Pieces.size() == 1);
    const Shape& drilled = trim.Result.Shape;
    REQUIRE(CountOf<EntityKind::Face>(drilled) == 4);
    const Stereon::EdgeId edge = drilled.Resolve(trim.Made.Pieces.front().Edge);
    REQUIRE(KindOf(drilled, edge) == EdgeCurveKind::Trim);
    REQUIRE(CachedCurveOf(drilled, edge).has_value());
}
