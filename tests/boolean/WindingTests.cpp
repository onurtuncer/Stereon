// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The bounded generalised winding number (M7, ADR-0013): the polyhedron a
// body's faces are triangulated into, accepted only where the point is
// separated from every triangle by more than the face can stray. A box and
// a drum give 1 inside and 0 outside, refine near a face and stay undecided
// on one; a drum cut by an oblique plane, whose lateral face then has a
// B-spline pcurve, is classified through that pcurve's derivative bounds.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>
#include <utility>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::Box;
using BooleanTest::ctx;
using BooleanTest::Drum;
using Catch::Matchers::WithinAbs;
using Stereon::Box2;
using Stereon::BranchSide;
using Stereon::ErrorCode;
using Stereon::FaceId;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::Surface;
using Stereon::Vec3;
using Stereon::WindingNumber;

namespace {

constexpr int Budget = 8;

WindingNumber Winding(const Shape& body, const Point3& p, const int refinements = Budget)
{
    const auto winding = WindingNumberOf(body, p, refinements, ctx);
    const std::string error = winding ? std::string() : std::string(winding.error().GetMessage());
    INFO(error);
    REQUIRE(winding.has_value());
    return *winding;
}

void RequireInside(const WindingNumber& w)
{
    REQUIRE(w.IsDecided());
    REQUIRE(w.IsInside());
    REQUIRE_THAT(w.Value, WithinAbs(1.0, 1e-9));
}

void RequireOutside(const WindingNumber& w)
{
    REQUIRE(w.IsDecided());
    REQUIRE_FALSE(w.IsInside());
    REQUIRE_THAT(w.Value, WithinAbs(0.0, 1e-9));
}

/// The drum with the oblique plane's ellipse imprinted on its lateral face,
/// as ImprintTests makes it: the lateral face splits in two bands along a
/// B-spline pcurve.
Shape CutDrum()
{
    const Shape drum = Drum();
    const FaceId lateral = BooleanTest::CylinderFace(drum);
    const Surface oblique{
        Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Normalized(Vec3{0.0, 1.0, -1.0})}};
    const auto ellipse =
        Intersect(oblique, Box2::Everything(), SurfaceOf(drum, lateral), Box2::Everything(), ctx);
    REQUIRE(ellipse.has_value());
    REQUIRE(ellipse->Branches.size() == 1);
    const auto& branch = ellipse->Branches.front();
    const auto clip = ClipToFace(drum, lateral, branch, BranchSide::B, ctx);
    REQUIRE(clip.has_value());
    auto builder = ShapeBuilder::Open(drum);
    REQUIRE(builder.has_value());
    const auto made = builder->Transaction([&](ShapeBuilder& candidate) {
        return ImprintClip(candidate, drum, lateral, branch, BranchSide::B, *clip, ctx);
    });
    REQUIRE(made.has_value());
    auto committed = Commit(std::move(*builder), ctx);
    REQUIRE(committed.has_value());
    REQUIRE(BooleanTest::CountOf<Stereon::EntityKind::Face>(committed->Shape) == 4);
    return committed->Shape;
}

} // namespace

TEST_CASE("the winding number is 1 inside a box and 0 outside", "[boolean][winding]")
{
    const Shape box = Box();
    const WindingNumber in = Winding(box, Point3{5.0, 5.0, 5.0});
    RequireInside(in);
    REQUIRE(in.Refinements == 0); // the first round separates a point 5 from every face
    REQUIRE(in.Triangles == 12);  // two per planar face, unsplit
    for (const Point3& p : {Point3{15.0, 5.0, 5.0}, Point3{5.0, -1.0, 5.0}, Point3{5.0, 5.0, 10.5},
                            Point3{20.0, 20.0, 20.0}, Point3{12.0, 5.0, 10.0}}) {
        INFO(p.X << ", " << p.Y << ", " << p.Z);
        RequireOutside(Winding(box, p));
    }
}

TEST_CASE("the winding number is exact near a box's face and undecided on it", "[boolean][winding]")
{
    const Shape box = Box();
    // Planar faces with line edges are triangulated exactly: their chords
    // are the edges, so the only deviation is the tolerances, four times
    // 1e-5, and a point 1e-4 inside is decided in the first round.
    const WindingNumber near = Winding(box, Point3{5.0, 5.0, 10.0 - 1e-4}, 0);
    RequireInside(near);
    REQUIRE(near.Refinements == 0);
    // Within those tolerances of the face no round can separate the point:
    // undecided, never a guess.
    const WindingNumber within = Winding(box, Point3{5.0, 5.0, 10.0 - 2e-5});
    REQUIRE_FALSE(within.IsDecided());
    REQUIRE(within.Bound >= 1.0);
    // On the face every round leaves a cell unseparated: undecided, never a
    // guess, and the interval says so.
    const WindingNumber on = Winding(box, Point3{5.0, 5.0, 10.0});
    REQUIRE_FALSE(on.IsDecided());
    REQUIRE(on.Bound >= 1.0);
    REQUIRE(on.Refinements >= 1);
}

TEST_CASE("the winding number refuses a negative budget", "[boolean][winding]")
{
    const auto winding = WindingNumberOf(Box(), Point3{5.0, 5.0, 5.0}, -1, ctx);
    REQUIRE_FALSE(winding.has_value());
    REQUIRE(winding.error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("the winding number of a drum comes from its tessellation", "[boolean][winding]")
{
    const Shape drum = Drum(); // radius 3, height 10, about z from the origin
    RequireInside(Winding(drum, Point3{0.0, 0.0, 5.0}));
    RequireInside(Winding(drum, Point3{2.5, 1.0, 0.5}));
    const WindingNumber nearWall = Winding(drum, Point3{2.999, 0.0, 5.0});
    RequireInside(nearWall);
    REQUIRE(nearWall.Refinements > 0);
    RequireOutside(Winding(drum, Point3{3.001, 0.0, 5.0}));
    RequireOutside(Winding(drum, Point3{0.0, 0.0, 10.5}));
    RequireOutside(Winding(drum, Point3{4.0, 4.0, 5.0}));
    RequireOutside(Winding(drum, Point3{0.0, 0.0, -0.001}));
    // On the wall.
    const WindingNumber on = Winding(drum, Point3{3.0, 0.0, 5.0});
    REQUIRE_FALSE(on.IsDecided());
    REQUIRE(on.Bound >= 1.0);
}

TEST_CASE("the winding number is 1 all round just inside a drum's wall", "[boolean][winding]")
{
    // A regression: corners a split placed on the surface left gaps against
    // the unsplit neighbour, and corners placed on long chords strayed past
    // the bound; either gave 0 here with a tiny bound.
    const Shape drum = Drum();
    const Shape cut = CutDrum();
    for (int k = 0; k < 12; ++k) {
        const double u = 0.5 * k + 0.1;
        const Point3 p{2.9 * std::cos(u), 2.9 * std::sin(u), 5.0};
        INFO("u = " << u);
        RequireInside(Winding(drum, p));
        RequireInside(Winding(cut, p));
        RequireOutside(Winding(cut, Point3{3.1 * std::cos(u), 3.1 * std::sin(u), 5.0}));
    }
}

TEST_CASE("the winding number of a drum cut along a B-spline pcurve", "[boolean][winding]")
{
    const Shape cut = CutDrum();
    // The ellipse lies in the plane z = 5 + y: points just inside and just
    // outside the wall near it, at angle 1 about the axis.
    const double c = std::cos(1.0);
    const double s = std::sin(1.0);
    RequireInside(Winding(cut, Point3{0.0, 0.0, 2.0}));
    RequireInside(Winding(cut, Point3{0.0, 0.0, 8.0}));
    const WindingNumber inside = Winding(cut, Point3{2.99 * c, 2.99 * s, 5.0 + 3.0 * s});
    RequireInside(inside);
    const WindingNumber outside = Winding(cut, Point3{3.01 * c, 3.01 * s, 5.0 + 3.0 * s});
    RequireOutside(outside);
    RequireOutside(Winding(cut, Point3{0.0, 0.0, 11.0}));
}
