// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The two-body driver (the Boolean prototype's pipeline, steps 1 to 4):
// every face pair intersected, every branch clipped to both faces, both
// bodies imprinted with the curve networks their faces receive, inside one
// transaction each, and committed valid. Two boxes overlapping at a corner
// meet in a closed chain of six lines through three junctions in each body;
// a cylinder through a box's top and bottom cuts two discs out of the box
// and two bands out of the cylinder.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <iterator>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::Box;
using BooleanTest::CountOf;
using BooleanTest::ctx;
using Stereon::Box3;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ImprintedBodies;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::Vec3;

namespace {

ImprintedBodies Imprinted(const Shape& a, const Shape& b)
{
    auto result = ImprintBodies(a, b, ctx);
    const std::string message = result ? std::string() : std::string(result.error().GetMessage());
    INFO(message);
    REQUIRE(result.has_value());
    return std::move(*result);
}

/// Every shared pair names live edges of the same geometry: ends within
/// tolerance of each other and the same kind.
void RequireShared(const ImprintedBodies& bodies, const std::size_t count)
{
    REQUIRE(bodies.Shared.size() == count);
    for (const auto& [ea, eb, transition] : bodies.Shared) {
        REQUIRE(ea.IsValid());
        REQUIRE(eb.IsValid());
        const Shape& a = bodies.A.Shape;
        const Shape& b = bodies.B.Shape;
        const Point3 startA = PositionOf(a, StartOf(a, ea));
        const Point3 endA = PositionOf(a, EndOf(a, ea));
        const Point3 startB = PositionOf(b, StartOf(b, eb));
        const Point3 endB = PositionOf(b, EndOf(b, eb));
        const double tolerance = 1e-9;
        REQUIRE(((Distance(startA, startB) <= tolerance && Distance(endA, endB) <= tolerance) ||
                 (Distance(startA, endB) <= tolerance && Distance(endA, startB) <= tolerance)));
    }
}

/// The edges with an end at the point, an edge with both ends there counted
/// twice. Within 1e-6: the pieces of a traced curve end where the tracing
/// put them, well inside the vertices' tolerance.
std::size_t EdgesAt(const Shape& shape, const Point3& p)
{
    std::size_t n = 0;
    for (const EdgeId e : *AllEntities<EntityKind::Edge>(shape)) {
        for (const Stereon::VertexId v : {StartOf(shape, e), EndOf(shape, e)}) {
            if (Distance(PositionOf(shape, v), p) <= 1e-6) {
                ++n;
            }
        }
    }
    return n;
}

/// The faces on a cylinder.
std::size_t LateralFacesOf(const Shape& shape)
{
    std::size_t n = 0;
    for (const Stereon::FaceId face : *AllEntities<EntityKind::Face>(shape)) {
        const Stereon::Surface surface = SurfaceOf(shape, face);
        if (surface.Target<Stereon::Cylinder>() != nullptr) {
            ++n;
        }
    }
    return n;
}

} // namespace

TEST_CASE("two boxes overlapping at a corner are imprinted with a six-piece chain",
          "[boolean][bodies]")
{
    const Shape a = Box();
    const auto other =
        MakeBoxSolid(Box3::FromCorners(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}), ctx);
    REQUIRE(other.has_value());
    const ImprintedBodies bodies = Imprinted(a, other->Shape);
    // Three faces of each box are cut, each by two pieces meeting at a
    // junction where an edge of the other box pierces it.
    for (const Stereon::OperationResult& body : {bodies.A, bodies.B}) {
        REQUIRE(CountOf<EntityKind::Face>(body.Shape) == 9);
        REQUIRE(CountOf<EntityKind::Vertex>(body.Shape) == 14);
        REQUIRE(CountOf<EntityKind::Edge>(body.Shape) == 21);
        REQUIRE_FALSE(body.TopologyEvents.empty());
    }
    RequireShared(bodies, 6);
    // Every piece lies on the cube [5, 10]^3 shared by both boxes.
    for (const auto& [ea, eb, transition] : bodies.Shared) {
        for (const Stereon::VertexId v : {StartOf(bodies.A.Shape, ea), EndOf(bodies.A.Shape, ea)}) {
            const Point3 p = PositionOf(bodies.A.Shape, v);
            REQUIRE(p.X >= 5.0 - 1e-9);
            REQUIRE(p.X <= 10.0 + 1e-9);
            REQUIRE(p.Y >= 5.0 - 1e-9);
            REQUIRE(p.Y <= 10.0 + 1e-9);
            REQUIRE(p.Z >= 5.0 - 1e-9);
            REQUIRE(p.Z <= 10.0 + 1e-9);
        }
    }
}

TEST_CASE("a cylinder through a box cuts two discs and two bands", "[boolean][bodies]")
{
    const Shape box = Box();
    const auto drum = MakeCylinderSolid(Point3{5.0, 5.0, -2.0}, Vec3{1.0, 0.0, 0.0},
                                        Vec3{0.0, 1.0, 0.0}, 3.0, 15.0, ctx);
    REQUIRE(drum.has_value());
    const ImprintedBodies bodies = Imprinted(box, drum->Shape);
    // The box: a ring and a disc on the top and on the bottom face.
    REQUIRE(CountOf<EntityKind::Face>(bodies.A.Shape) == 8);
    REQUIRE(CountOf<EntityKind::Edge>(bodies.A.Shape) == 14);
    REQUIRE(CountOf<EntityKind::Loop>(bodies.A.Shape) == 10);
    // The cylinder: its lateral face cut into three bands by two circles,
    // each a piece from the seam round to the seam.
    REQUIRE(CountOf<EntityKind::Face>(bodies.B.Shape) == 5);
    REQUIRE(CountOf<EntityKind::Edge>(bodies.B.Shape) == 7);
    RequireShared(bodies, 2);
    for (const auto& [ea, eb, transition] : bodies.Shared) {
        REQUIRE(CurveOf(bodies.A.Shape, ea)->Target<Stereon::Circle>() != nullptr);
        REQUIRE(KindOf(bodies.B.Shape, eb) == Stereon::EdgeCurveKind::Exact);
    }
}

TEST_CASE("a cylinder pierced across its seam is cut into faces that close in its parameters",
          "[boolean][bodies][seam][regression]")
{
    // The upright cylinder's seam is at u = 0, on +x. A cylinder along x,
    // through its axis or 1 off it, pierces it across the seam, where the
    // hole's curve meets the seam twice, and opposite it at u = pi, where a
    // closed curve meets no edge. Before #214 the second half of the curve
    // at the seam was joined to the seam's use on the far side, a jump of
    // 2 pi in a loop (G03; S07 while the seam was a Boundary edge), and the
    // closed curve was hung on whichever part kept the face's handle, a
    // half disc at u = 2 pi (docs/imprinting.rst, *Seams*).
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    const Vec3 ez{0.0, 0.0, 1.0};
    const auto upright = MakeCylinderSolid(Point3{0.0, 0.0, -5.0}, ex, ey, 3.0, 10.0, ctx);
    REQUIRE(upright.has_value());
    for (const double offset : {0.0, 1.0}) {
        INFO("offset " << offset);
        const auto across = MakeCylinderSolid(Point3{-6.0, offset, 0.0}, ey, ez, 1.5, 12.0, ctx);
        REQUIRE(across.has_value());
        // Committed at the Geometric level, so every loop closes in the
        // surface's parameters.
        const ImprintedBodies bodies = Imprinted(upright->Shape, across->Shape);
        const Shape& a = bodies.A.Shape;
        // The lateral face: the outer part, the two halves of the disc
        // across the seam, which the seam's middle piece divides, and the
        // disc at u = pi.
        std::size_t lateral = 0;
        std::size_t holed = 0;
        for (const Stereon::FaceId face : *AllEntities<EntityKind::Face>(a)) {
            const Stereon::Surface surface = SurfaceOf(a, face);
            if (surface.Target<Stereon::Cylinder>() == nullptr) {
                continue;
            }
            ++lateral;
            const auto loops = LoopsOf(a, face);
            REQUIRE(loops.has_value());
            const auto count = std::ranges::distance(loops->begin(), loops->end());
            if (count < 2) {
                continue;
            }
            // The face with a hole is the outer part, the one running round
            // the whole cylinder from its bottom circle to its top.
            ++holed;
            const auto outer = CoedgesOf(a, *loops->begin());
            REQUIRE(outer.has_value());
            REQUIRE(std::ranges::distance(outer->begin(), outer->end()) >= 6);
        }
        REQUIRE(lateral == 4);
        REQUIRE(holed == 1);
    }
}

TEST_CASE("boxes sharing their top and bottom planes are imprinted from the other's corner",
          "[boolean][bodies][coincident][regression]")
{
    // B = [5, 15]^2 x [0, 10] shares A's top and bottom planes. In each, the
    // overlap [5, 10]^2 is bounded inside A's face by two parts of B's
    // edges, which meet at B's corner (5, 5) and run out to A's boundary:
    // two pieces that both start at a junction no other piece reaches. A
    // strut used to be grown only from the start of a piece, so neither
    // could be written and the imprint refused ("a piece runs into the face
    // from a junction no other piece reaches"). Now the first is grown from
    // its far end, towards the junction, and the second closes the face
    // (docs/imprinting.rst, *Writing the pieces*).
    const Shape a = Box();
    const auto other =
        MakeBoxSolid(Box3::FromCorners(Point3{5.0, 5.0, 0.0}, Point3{15.0, 15.0, 10.0}), ctx);
    REQUIRE(other.has_value());
    const ImprintedBodies bodies = Imprinted(a, other->Shape);
    // Each box: its top and its bottom cut in two by two pieces, two of its
    // sides by one line each, and four of its edges split where the other
    // box's sides cross them.
    for (const Stereon::OperationResult& body : {bodies.A, bodies.B}) {
        REQUIRE(CountOf<EntityKind::Face>(body.Shape) == 10);
        REQUIRE(CountOf<EntityKind::Edge>(body.Shape) == 22);
        REQUIRE(CountOf<EntityKind::Vertex>(body.Shape) == 14);
    }
    // The square's four sides in the top and in the bottom, and the two
    // lines down the sides.
    RequireShared(bodies, 10);
    for (const double z : {0.0, 10.0}) {
        INFO("z = " << z);
        // B's corner in A's face and A's in B's: a junction of two pieces.
        REQUIRE(EdgesAt(bodies.A.Shape, Point3{5.0, 5.0, z}) == 2);
        REQUIRE(EdgesAt(bodies.B.Shape, Point3{10.0, 10.0, z}) == 2);
        // The same points in the box whose corners they are: its three edges.
        REQUIRE(EdgesAt(bodies.B.Shape, Point3{5.0, 5.0, z}) == 3);
        REQUIRE(EdgesAt(bodies.A.Shape, Point3{10.0, 10.0, z}) == 3);
    }
}

TEST_CASE("perpendicular cylinders of equal radii are imprinted through the nodes of their curves",
          "[boolean][bodies][seam][regression]")
{
    // The upright cylinder about z and one of the same radius about x meet
    // in two ellipses, in the planes x = z and x = -z, which cross at
    // (0, 3, 0) and (0, -3, 0), where the cylinders are tangent: four
    // pieces at one vertex in the middle of a face. A piece joining such a
    // vertex has a use to choose for each part of the face round it, all at
    // one place in the face's parameters. It used to take one by its sense
    // alone, which could be a use of the part across the crossing, and the
    // imprint refused ("a piece's ends are not both on the boundary of one
    // face"). Now it takes the use whose corner of the face holds the
    // piece's direction (docs/imprinting.rst, *Where a piece joins a
    // vertex*).
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    const Vec3 ez{0.0, 0.0, 1.0};
    const auto upright = MakeCylinderSolid(Point3{0.0, 0.0, -5.0}, ex, ey, 3.0, 10.0, ctx);
    const auto across = MakeCylinderSolid(Point3{-6.0, 0.0, 0.0}, ey, ez, 3.0, 12.0, ctx);
    REQUIRE(upright.has_value());
    REQUIRE(across.has_value());
    // Committed at the Geometric level, so every loop closes in its
    // surface's parameters.
    const ImprintedBodies bodies = Imprinted(upright->Shape, across->Shape);
    const Shape& a = bodies.A.Shape;
    const Shape& b = bodies.B.Shape;
    // The upright one, its seam at (3, 0, z): the parts above and below
    // both ellipses, the lens about -x between them, and the lens about +x,
    // which the seam divides. The ellipses cross the seam at z = 3 and
    // z = -3: six pieces, and the seam in three parts.
    REQUIRE(LateralFacesOf(a) == 5);
    REQUIRE(CountOf<EntityKind::Face>(a) == 7);
    REQUIRE(CountOf<EntityKind::Edge>(a) == 11);
    REQUIRE(CountOf<EntityKind::Vertex>(a) == 6);
    // The one about x, its seam at (x, 3, 0), through a node: the parts
    // beyond both ellipses towards either cap, each holding a half of the
    // seam, and the lenses above and below.
    REQUIRE(LateralFacesOf(b) == 4);
    REQUIRE(CountOf<EntityKind::Face>(b) == 6);
    REQUIRE(CountOf<EntityKind::Edge>(b) == 10);
    REQUIRE(CountOf<EntityKind::Vertex>(b) == 6);
    RequireShared(bodies, 6);
    // Four pieces at each node; at the one on B's seam, its two halves too.
    REQUIRE(EdgesAt(a, Point3{0.0, 3.0, 0.0}) == 4);
    REQUIRE(EdgesAt(a, Point3{0.0, -3.0, 0.0}) == 4);
    REQUIRE(EdgesAt(b, Point3{0.0, 3.0, 0.0}) == 6);
    REQUIRE(EdgesAt(b, Point3{0.0, -3.0, 0.0}) == 4);
    // Where the ellipses cross A's seam: two pieces and two parts of the
    // seam in A, two pieces in B.
    for (const double z : {3.0, -3.0}) {
        INFO("z = " << z);
        REQUIRE(EdgesAt(a, Point3{3.0, 0.0, z}) == 4);
        REQUIRE(EdgesAt(b, Point3{3.0, 0.0, z}) == 2);
    }
}

TEST_CASE("bodies apart are imprinted with nothing", "[boolean][bodies]")
{
    const Shape a = Box();
    const auto far =
        MakeBoxSolid(Box3::FromCorners(Point3{20.0, 20.0, 20.0}, Point3{30.0, 30.0, 30.0}), ctx);
    REQUIRE(far.has_value());
    const ImprintedBodies bodies = Imprinted(a, far->Shape);
    REQUIRE(bodies.Shared.empty());
    REQUIRE(CountOf<EntityKind::Face>(bodies.A.Shape) == 6);
    REQUIRE(CountOf<EntityKind::Face>(bodies.B.Shape) == 6);
}
