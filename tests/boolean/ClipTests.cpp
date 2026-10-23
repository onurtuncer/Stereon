// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Face clipping (the Boolean prototype's M6): a branch of two surfaces
// restricted to a face on one of them, by its pcurve against the face's
// loops. On a box's top face: a line through two edges, through two
// corners (crossings at the vertices, ADR-0001), beside the face, a
// circle inside it and a circle through all four edges. On a cylinder's
// lateral face: an ellipse's sinusoid pcurve wholly inside, and one that
// leaves through the top circle.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Boolean;

using Catch::Matchers::WithinAbs;
using Stereon::Box2;
using Stereon::Box3;
using Stereon::BranchSide;
using Stereon::ContactType;
using Stereon::Context;
using Stereon::Cylinder;
using Stereon::EntityKind;
using Stereon::FaceClip;
using Stereon::FaceId;
using Stereon::IntersectionBranch;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Context ctx = Context::Default();

/// The unit box scaled to 10, from the origin.
Shape Box()
{
    const auto box = MakeBoxSolid(Box3::FromCorners(Point3{}, Point3{10.0, 10.0, 10.0}), ctx);
    REQUIRE(box.has_value());
    return box->Shape;
}

/// A cylinder of radius 3 and height 10 about z from the origin.
Shape Drum()
{
    const auto drum =
        MakeCylinderSolid(Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, 3.0, 10.0, ctx);
    REQUIRE(drum.has_value());
    return drum->Shape;
}

/// The face whose plane's normal points along n.
FaceId PlaneFaceAlong(const Shape& shape, const Vec3& n)
{
    const auto faces = AllEntities<EntityKind::Face>(shape);
    REQUIRE(faces.has_value());
    for (const FaceId face : *faces) {
        const Surface surface = SurfaceOf(shape, face);
        if (const Plane* plane = surface.Target<Plane>();
            plane != nullptr && Dot(Normalized(Cross(plane->XAxis, plane->YAxis)), n) > 0.9) {
            return face;
        }
    }
    FAIL("no face along the direction");
    return {};
}

/// The face on a cylinder.
FaceId CylinderFace(const Shape& shape)
{
    const auto faces = AllEntities<EntityKind::Face>(shape);
    REQUIRE(faces.has_value());
    for (const FaceId face : *faces) {
        if (SurfaceOf(shape, face).Target<Cylinder>() != nullptr) {
            return face;
        }
    }
    FAIL("no cylindrical face");
    return {};
}

/// The one branch of the face's surface against `other`, the face's surface
/// as side A, over the face's parameter box (planes are unbounded).
IntersectionBranch BranchOn(const Shape& shape, const FaceId face, const Surface& other,
                            const Box2& box)
{
    const auto result = Intersect(SurfaceOf(shape, face), box, other, Box2::Everything(), ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    return result->Branches.front();
}

Surface VerticalPlane(const double x)
{
    return Surface{Plane{Point3{x, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};
}

/// The box's face parameters span [0, 10] in both directions; a margin lets
/// the lines be clipped beside the face too.
const Box2 faceBox{Interval{-1.0, 11.0}, Interval{-1.0, 11.0}};

} // namespace

TEST_CASE("a line through a face crosses two edges and lies inside between them", "[boolean][clip]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const IntersectionBranch branch = BranchOn(box, top, VerticalPlane(3.0), faceBox);
    const auto clip = ClipToFace(box, top, branch, BranchSide::A, ctx);
    REQUIRE(clip.has_value());
    REQUIRE(clip->Decided);
    REQUIRE(clip->OnBoundary.empty());
    REQUIRE(clip->Crossings.size() == 2);
    for (const Stereon::LoopCrossing& c : clip->Crossings) {
        REQUIRE_FALSE(c.Vertex.IsValid());
        REQUIRE(c.Contact == ContactType::Crossing);
        REQUIRE(c.Decided);
        const Point3 p = Eval(branch.Curve3d, c.OnBranch);
        REQUIRE_THAT(p.X, WithinAbs(3.0, 1e-12));
        REQUIRE_THAT(p.Z, WithinAbs(10.0, 1e-12));
        REQUIRE_THAT(std::abs(p.Y - 5.0), WithinAbs(5.0, 1e-12));
    }
    REQUIRE(clip->Inside.size() == 1);
    const Stereon::FacePiece& piece = clip->Inside.front();
    REQUIRE(piece.Start.has_value());
    REQUIRE(piece.End.has_value());
    REQUIRE(piece.Range.GetLower() == clip->Crossings[*piece.Start].OnBranch);
    REQUIRE(piece.Range.GetUpper() == clip->Crossings[*piece.End].OnBranch);
    const Point3 middle =
        Eval(branch.Curve3d, 0.5 * (piece.Range.GetLower() + piece.Range.GetUpper()));
    REQUIRE_THAT(middle.Y, WithinAbs(5.0, 1e-12));
}

TEST_CASE("a line through two corners crosses at the vertices", "[boolean][clip]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const Surface diagonal{Plane{Point3{}, Normalized(Vec3{1.0, 1.0, 0.0}), Vec3{0.0, 0.0, 1.0}}};
    const IntersectionBranch branch = BranchOn(box, top, diagonal, faceBox);
    const auto clip = ClipToFace(box, top, branch, BranchSide::A, ctx);
    REQUIRE(clip.has_value());
    REQUIRE(clip->Crossings.size() == 2);
    for (const Stereon::LoopCrossing& c : clip->Crossings) {
        REQUIRE(c.Vertex.IsValid());
        const Point3 p = Eval(branch.Curve3d, c.OnBranch);
        REQUIRE(Distance(p, PositionOf(box, c.Vertex)) <= ToleranceOf(box, c.Vertex));
        REQUIRE_THAT(std::abs(p.X - 5.0), WithinAbs(5.0, 1e-12));
    }
    REQUIRE(clip->Inside.size() == 1);
    REQUIRE(clip->Inside.front().Start.has_value());
    REQUIRE(clip->Inside.front().End.has_value());
}

TEST_CASE("a line beside a face has no crossing and nothing inside", "[boolean][clip]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const IntersectionBranch branch = BranchOn(box, top, VerticalPlane(-0.5), faceBox);
    const auto clip = ClipToFace(box, top, branch, BranchSide::A, ctx);
    REQUIRE(clip.has_value());
    REQUIRE(clip->Crossings.empty());
    REQUIRE(clip->Inside.empty());
    REQUIRE(clip->Decided);
}

TEST_CASE("a circle inside a face is one closed piece; through its edges, four arcs",
          "[boolean][clip]")
{
    const Shape box = Box();
    const FaceId top = PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const auto drill = [](const double radius) {
        return Surface{Cylinder{Point3{5.0, 5.0, 0.0}, radius, Vec3{1.0, 0.0, 0.0},
                                Vec3{0.0, 1.0, 0.0}, true}};
    };
    const IntersectionBranch small = BranchOn(box, top, drill(2.0), faceBox);
    const auto inside = ClipToFace(box, top, small, BranchSide::A, ctx);
    REQUIRE(inside.has_value());
    REQUIRE(inside->Crossings.empty());
    REQUIRE(inside->Inside.size() == 1);
    REQUIRE(inside->Inside.front().Range == small.Range);
    REQUIRE_FALSE(inside->Inside.front().Start.has_value());
    REQUIRE_FALSE(inside->Inside.front().End.has_value());

    // Radius 8 reaches past the corners: the circle runs around the face.
    const IntersectionBranch around = BranchOn(box, top, drill(8.0), faceBox);
    const auto outside = ClipToFace(box, top, around, BranchSide::A, ctx);
    REQUIRE(outside.has_value());
    REQUIRE(outside->Crossings.empty());
    REQUIRE(outside->Inside.empty());

    // Radius 6 crosses every edge twice, 5 from the centre, and misses the corners.
    const IntersectionBranch large = BranchOn(box, top, drill(6.0), faceBox);
    const auto arcs = ClipToFace(box, top, large, BranchSide::A, ctx);
    REQUIRE(arcs.has_value());
    REQUIRE(arcs->Decided);
    REQUIRE(arcs->Crossings.size() == 8);
    REQUIRE(arcs->Inside.size() == 4);
    double covered = 0.0;
    for (const Stereon::FacePiece& piece : arcs->Inside) {
        REQUIRE(piece.Start.has_value());
        REQUIRE(piece.End.has_value());
        const Point3 middle =
            Eval(large.Curve3d, 0.5 * (piece.Range.GetLower() + piece.Range.GetUpper()));
        REQUIRE(std::abs(middle.X - 5.0) < 5.0);
        REQUIRE(std::abs(middle.Y - 5.0) < 5.0);
        covered += piece.Range.GetUpper() - piece.Range.GetLower();
    }
    // The arcs inside the square: the whole turn less four arcs of 2 acos(5/6).
    REQUIRE_THAT(covered, WithinAbs(2.0 * pi - 8.0 * std::acos(5.0 / 6.0), 1e-9));
}

TEST_CASE("an ellipse's sinusoid pcurve is clipped on a cylinder's lateral face", "[boolean][clip]")
{
    const Shape drum = Drum();
    const FaceId lateral = CylinderFace(drum);
    const auto oblique = [](const double z) {
        return Surface{
            Plane{Point3{0.0, 0.0, z}, Vec3{1.0, 0.0, 0.0}, Normalized(Vec3{0.0, 1.0, -1.0})}};
    };
    // The plane first: the cylinder is side B. Heights 5 ± 3: wholly inside.
    const auto within = Intersect(oblique(5.0), Box2::Everything(), SurfaceOf(drum, lateral),
                                  Box2::Everything(), ctx);
    REQUIRE(within.has_value());
    REQUIRE(within->Branches.size() == 1);
    const auto whole = ClipToFace(drum, lateral, within->Branches.front(), BranchSide::B, ctx);
    REQUIRE(whole.has_value());
    REQUIRE_FALSE(whole->Decided); // A B-spline pcurve: bracketed to rounding.
    REQUIRE(whole->Inside.size() == 1);
    REQUIRE(whole->Inside.front().Range == within->Branches.front().Range);

    // Heights 9 ± 3: it leaves through the top circle and comes back.
    const auto leaving = Intersect(oblique(9.0), Box2::Everything(), SurfaceOf(drum, lateral),
                                   Box2::Everything(), ctx);
    REQUIRE(leaving.has_value());
    const IntersectionBranch& branch = leaving->Branches.front();
    const auto part = ClipToFace(drum, lateral, branch, BranchSide::B, ctx);
    REQUIRE(part.has_value());
    std::size_t topCrossings = 0;
    for (const Stereon::LoopCrossing& c : part->Crossings) {
        const Point3 p = Eval(branch.Curve3d, c.OnBranch);
        if (std::abs(p.Z - 10.0) <= 1e-6) {
            ++topCrossings;
        }
    }
    REQUIRE(topCrossings == 2);
    double covered = 0.0;
    for (const Stereon::FacePiece& piece : part->Inside) {
        const Point3 middle =
            Eval(branch.Curve3d, 0.5 * (piece.Range.GetLower() + piece.Range.GetUpper()));
        REQUIRE(middle.Z < 10.0);
        covered += piece.Range.GetUpper() - piece.Range.GetLower();
    }
    // Inside where 9 + 3 cos t < 10: all but 2 acos(1/3) of the turn.
    REQUIRE_THAT(covered, WithinAbs(2.0 * pi - 2.0 * std::acos(1.0 / 3.0), 1e-6));
}
