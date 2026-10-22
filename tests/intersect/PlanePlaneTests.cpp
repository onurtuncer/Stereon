// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The plane–plane handler (the Boolean prototype's M3, ADR-0010, ADR-0029):
// a line for planes that are not parallel, with both pcurves on the line's
// parameter and the range clipped to the requested boxes; a proved
// separation for parallel planes beyond the pair tolerance; an overlap,
// exact or declared, within it, with the parameter map between the planes.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using IntersectTest::PcurveDeviation;
using Stereon::Box2;
using Stereon::CoincidenceEvidence;
using Stereon::ContactType;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Surface xy{Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
const Surface yz{Plane{Point3{}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};

Surface Lifted(const double z)
{
    return Surface{Plane{Point3{0.0, 0.0, z}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
}

} // namespace

TEST_CASE("planes that are not parallel meet in a line on both pcurves", "[intersect][plane]")
{
    const auto result = Intersect(xy, yz, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.empty());
    REQUIRE(result->Overlaps.empty());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& line = result->Branches.front();
    REQUIRE(line.Curve3d.Target<Stereon::Line>() != nullptr);
    REQUIRE(line.Contact == ContactType::Crossing);
    REQUIRE(line.Transition == SideTransition::OppositeSide);
    REQUIRE_FALSE(line.IsCertified());
    REQUIRE(line.Range == Interval::Everything());
    // Along the y axis, at unit speed, through the origin.
    const Stereon::Line& l = *line.Curve3d.Target<Stereon::Line>();
    REQUIRE_THAT(std::abs(l.Direction.Y), WithinAbs(1.0, 1e-15));
    REQUIRE(Distance(l.Origin, Point3{}) == 0.0);
    REQUIRE(PcurveDeviation(line, xy, yz) <= 1e-14);
}

TEST_CASE("the line's range is clipped to the requested boxes", "[intersect][plane]")
{
    // The line is x = z = 0, which both pcurves give as a line along their
    // first axis (xy: y; yz: its X axis is world y).
    const Box2 boxA{Interval{-1.0, 1.0}, Interval{0.0, 2.0}};
    const Box2 boxB{Interval{0.0, 5.0}, Interval{-7.0, 7.0}};
    const auto result = Intersect(xy, boxA, yz, boxB, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Interval range = result->Branches.front().Range;
    const double sign = result->Branches.front().Curve3d.Target<Stereon::Line>()->Direction.Y;
    // Within y in [0, 2], either way along the line.
    // Outward-rounded, so within a few ulps.
    REQUIRE_THAT(range.GetLower(), WithinAbs(sign > 0.0 ? 0.0 : -2.0, 1e-14));
    REQUIRE_THAT(range.GetUpper(), WithinAbs(sign > 0.0 ? 2.0 : 0.0, 1e-14));

    // A box that the line misses: a proved separation over the domains.
    const Box2 aside{Interval{3.0, 4.0}, Interval{0.0, 2.0}};
    const auto missed = Intersect(xy, aside, yz, boxB, ctx);
    REQUIRE(missed.has_value());
    REQUIRE(missed->IsEmpty());
}

TEST_CASE("nearly parallel planes still give their line", "[intersect][plane]")
{
    // Tilted by a nanoradian about x: not exactly parallel, so a line, far
    // away but on both planes.
    const double e = 1e-9;
    const Surface tilted{
        Plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, std::cos(e), std::sin(e)}}};
    const auto result = Intersect(xy, tilted, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::Line& l = *result->Branches.front().Curve3d.Target<Stereon::Line>();
    REQUIRE_THAT(std::abs(l.Direction.X), WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(std::abs(l.Origin.Y), WithinAbs(1.0 / std::tan(e), 1e-3));
    REQUIRE(PcurveDeviation(result->Branches.front(), xy, tilted) <= 1e-6);
}

TEST_CASE("parallel planes beyond the pair tolerance are proved separate", "[intersect][plane]")
{
    const auto apart =
        Intersect(xy, Box2::Everything(), Lifted(1e-3), Box2::Everything(), 1e-4, ctx);
    REQUIRE(apart.has_value());
    REQUIRE(apart->IsEmpty());
    // At the Context's resolution by default: 1e-5 < 1e-3.
    REQUIRE(Intersect(xy, Lifted(1e-3), ctx)->IsEmpty());
    // Opposite normals change nothing.
    const Surface down{Plane{Point3{0.0, 0.0, 1e-3}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}}};
    REQUIRE(Intersect(xy, down, ctx)->IsEmpty());
}

TEST_CASE("parallel planes within the pair tolerance are declared coincident", "[intersect][plane]")
{
    std::vector<std::string> log;
    const Stereon::Context logging = IntersectTest::Logging(log);
    const auto result =
        Intersect(xy, Box2::Everything(), Lifted(1e-3), Box2::Everything(), 2e-3, logging);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Overlaps.size() == 1);
    const Stereon::OverlapRegion& overlap = result->Overlaps.front();
    REQUIRE(overlap.Evidence == CoincidenceEvidence::Tolerant);
    REQUIRE_THAT(overlap.Distance, WithinAbs(1e-3, 1e-18));
    REQUIRE(overlap.SameOrientation);
    REQUIRE(overlap.FullDomain);
    REQUIRE(overlap.Boundary.empty());
    REQUIRE(log.size() == 1);
    REQUIRE(log.front().find("declared coincident") != std::string::npos);
    // Exactly at the tolerance counts as within it.
    REQUIRE(Intersect(xy, Box2::Everything(), Lifted(1e-3), Box2::Everything(), 1e-3, ctx)
                ->Overlaps.size() == 1);
}

TEST_CASE("exactly coincident planes map one frame onto the other", "[intersect][plane]")
{
    // The same plane z = 5 with a turned, shifted frame and the same normal.
    const Surface a{Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Surface b{Plane{Point3{3.0, 4.0, 5.0}, Vec3{0.0, 1.0, 0.0}, Vec3{-1.0, 0.0, 0.0}}};
    const auto result = Intersect(a, b, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.size() == 1);
    const Stereon::OverlapRegion& overlap = result->Overlaps.front();
    REQUIRE(overlap.Evidence == CoincidenceEvidence::Exact);
    REQUIRE(overlap.Distance == 0.0);
    REQUIRE(overlap.SameOrientation);
    REQUIRE(overlap.FullDomain);
    for (const Point2 uv : {Point2{0.0, 0.0}, Point2{2.0, -3.0}, Point2{-7.5, 1.25}}) {
        REQUIRE(Distance(Eval(b, Mapped(overlap.AToB, uv)), Eval(a, uv)) <= 1e-14);
    }

    // The same plane with the normal turned round.
    const Surface reversed{Plane{Point3{1.0, 1.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}}};
    const auto opposite = Intersect(a, reversed, ctx);
    REQUIRE(opposite.has_value());
    REQUIRE(opposite->Overlaps.size() == 1);
    REQUIRE_FALSE(opposite->Overlaps.front().SameOrientation);
    REQUIRE(opposite->Overlaps.front().Evidence == CoincidenceEvidence::Exact);
}

TEST_CASE("an overlap is restricted to the requested boxes", "[intersect][plane]")
{
    const Surface a{Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Box2 big{Interval{0.0, 2.0}, Interval{0.0, 2.0}};
    const Box2 small{Interval{0.0, 1.0}, Interval{0.0, 1.0}};
    const auto whole = Intersect(a, big, a, big, ctx);
    REQUIRE(whole->Overlaps.front().FullDomain);
    REQUIRE(whole->Overlaps.front().OnA == big);

    const auto part = Intersect(a, big, a, small, ctx);
    REQUIRE(part->Overlaps.size() == 1);
    REQUIRE_FALSE(part->Overlaps.front().FullDomain);
    // The image of a box is outward-rounded: the small box, within ulps.
    for (const Box2& on : {part->Overlaps.front().OnA, part->Overlaps.front().OnB}) {
        REQUIRE(Contains(on, small));
        REQUIRE_THAT(Width(on.X), WithinAbs(1.0, 1e-14));
        REQUIRE_THAT(Width(on.Y), WithinAbs(1.0, 1e-14));
    }

    // Boxes that do not meet: coincident supports, no overlap in the domains.
    const Box2 elsewhere{Interval{5.0, 6.0}, Interval{0.0, 1.0}};
    const auto none = Intersect(a, big, a, elsewhere, ctx);
    REQUIRE(none.has_value());
    REQUIRE(none->IsEmpty());
}

TEST_CASE("the pair tolerance must be finite and not negative", "[intersect][plane]")
{
    REQUIRE(
        Intersect(xy, Box2::Everything(), yz, Box2::Everything(), -1e-3, ctx).error().GetCode() ==
        ErrorCode::InvalidInput);
    REQUIRE(Intersect(xy, Box2::Everything(), yz, Box2::Everything(), 0.0, ctx).has_value());
}
