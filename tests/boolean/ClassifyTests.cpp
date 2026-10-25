// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Classification (the Boolean prototype's M7, ADR-0013): the regions of two
// imprinted bodies, each seeded by deterministic rays against the other
// body and cross-checked across the imprinted edges by the branches' side
// transitions. Two boxes overlapping at a corner: one region of each inside
// the other, one outside. A cylinder through a box: the box's two discs
// inside the cylinder, the cylinder's middle band inside the box.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::Box;
using BooleanTest::ctx;
using Stereon::Box3;
using Stereon::Classification;
using Stereon::ImprintedBodies;
using Stereon::Point3;
using Stereon::Region;
using Stereon::Shape;
using Stereon::Side;
using Stereon::Vec3;

namespace {

Classification Classified(const Shape& a, const Shape& b)
{
    auto bodies = ImprintBodies(a, b, ctx);
    const std::string imprintError =
        bodies ? std::string() : std::string(bodies.error().GetMessage());
    INFO(imprintError);
    REQUIRE(bodies.has_value());
    auto result = Classify(*bodies, ctx);
    const std::string classifyError =
        result ? std::string() : std::string(result.error().GetMessage());
    INFO(classifyError);
    REQUIRE(result.has_value());
    return std::move(*result);
}

/// How many regions of the body lie on each side, and the faces they hold.
struct Tally {
    std::size_t Inside = 0;
    std::size_t Outside = 0;
    std::size_t FacesInside = 0;
    std::size_t FacesOutside = 0;
};

Tally TallyOf(const std::vector<Region>& regions)
{
    Tally t;
    for (const Region& region : regions) {
        REQUIRE(region.Side != Side::Unresolved);
        REQUIRE(region.Rays >= 1);
        REQUIRE(region.SeedFace.IsValid());
        if (region.Side == Side::Inside) {
            ++t.Inside;
            t.FacesInside += region.Faces.size();
        } else {
            ++t.Outside;
            t.FacesOutside += region.Faces.size();
        }
    }
    return t;
}

} // namespace

TEST_CASE("two boxes overlapping at a corner each have one region inside the other",
          "[boolean][classify]")
{
    const auto other =
        MakeBoxSolid(Box3::FromCorners(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}), ctx);
    REQUIRE(other.has_value());
    const Classification c = Classified(Box(), other->Shape);
    for (const std::vector<Region>& regions : {c.A, c.B}) {
        REQUIRE(regions.size() == 2);
        const Tally t = TallyOf(regions);
        REQUIRE(t.Inside == 1);
        REQUIRE(t.Outside == 1);
        // The three corner squares inside, the three L-shaped faces and the
        // three untouched faces outside.
        REQUIRE(t.FacesInside == 3);
        REQUIRE(t.FacesOutside == 6);
    }
    // The inside region's seed lies in the shared cube.
    for (const Region& region : c.A) {
        if (region.Side == Side::Inside) {
            REQUIRE(region.Seed.X > 5.0);
            REQUIRE(region.Seed.Y > 5.0);
            REQUIRE(region.Seed.Z > 5.0);
        }
    }
}

TEST_CASE("a cylinder through a box: the discs inside the cylinder, the band inside the box",
          "[boolean][classify]")
{
    const auto drum = MakeCylinderSolid(Point3{5.0, 5.0, -2.0}, Vec3{1.0, 0.0, 0.0},
                                        Vec3{0.0, 1.0, 0.0}, 3.0, 15.0, ctx);
    REQUIRE(drum.has_value());
    const Classification c = Classified(Box(), drum->Shape);
    // The box: the two discs, each its own region inside; everything else one
    // region outside.
    const Tally box = TallyOf(c.A);
    REQUIRE(c.A.size() == 3);
    REQUIRE(box.Inside == 2);
    REQUIRE(box.FacesInside == 2);
    REQUIRE(box.Outside == 1);
    REQUIRE(box.FacesOutside == 6);
    // The cylinder: the middle band inside the box; each end band with its
    // cap outside.
    const Tally drumTally = TallyOf(c.B);
    REQUIRE(c.B.size() == 3);
    REQUIRE(drumTally.Inside == 1);
    REQUIRE(drumTally.FacesInside == 1);
    REQUIRE(drumTally.Outside == 2);
    REQUIRE(drumTally.FacesOutside == 4);
}

TEST_CASE("bodies apart are each one region outside the other", "[boolean][classify]")
{
    const auto far =
        MakeBoxSolid(Box3::FromCorners(Point3{20.0, 20.0, 20.0}, Point3{30.0, 30.0, 30.0}), ctx);
    REQUIRE(far.has_value());
    const Classification c = Classified(Box(), far->Shape);
    REQUIRE(c.A.size() == 1);
    REQUIRE(c.A.front().Side == Side::Outside);
    REQUIRE(c.B.size() == 1);
    REQUIRE(c.B.front().Side == Side::Outside);
}
