// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The cylinder–cylinder handler for parallel axes (the Boolean prototype's
// M3, ADR-0010, ADR-0029): a coincidence, exact or declared within the pair
// tolerance, with the map between the frames; an external or internal
// tangent ruling; two rulings; a proved separation apart or nested. Skew
// axes are M3b and report NotImplemented.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
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
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

Surface AboutZ(const Point3& origin, const double radius)
{
    return Surface{Cylinder{origin, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
}

/// About -z: the frame turns the other way.
Surface AboutMinusZ(const Point3& origin, const double radius)
{
    return Surface{Cylinder{origin, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}, true}};
}

const Surface a = AboutZ(Point3{}, 3.0);

} // namespace

TEST_CASE("coaxial cylinders of one radius coincide exactly, frames mapped",
          "[intersect][cylinder]")
{
    const Surface b = AboutMinusZ(Point3{0.0, 0.0, 9.0}, 3.0);
    const auto result = Intersect(a, b, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Overlaps.size() == 1);
    const Stereon::OverlapRegion& overlap = result->Overlaps.front();
    REQUIRE(overlap.Evidence == CoincidenceEvidence::Exact);
    REQUIRE(overlap.Distance == 0.0);
    REQUIRE(overlap.SameOrientation);
    REQUIRE(overlap.FullDomain);
    // B's angle runs the other way and its height is measured downwards.
    REQUIRE(overlap.AToB.U == Stereon::Vec2{-1.0, 0.0});
    REQUIRE(overlap.AToB.V == Stereon::Vec2{0.0, -1.0});
    REQUIRE(overlap.AToB.Offset == Stereon::Vec2{0.0, 9.0});
    for (const Point2 uv : {Point2{0.0, 0.0}, Point2{1.0, -3.0}, Point2{4.5, 12.0}}) {
        REQUIRE(Distance(Eval(b, Mapped(overlap.AToB, uv)), Eval(a, uv)) <= 1e-14);
    }
    // Turned inside out: the normals disagree.
    const Surface inward{Cylinder{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, false}};
    REQUIRE_FALSE(Intersect(a, inward, ctx)->Overlaps.front().SameOrientation);
}

TEST_CASE("nearly coaxial cylinders within the pair tolerance are declared coincident",
          "[intersect][cylinder]")
{
    // Axes 3e-4 apart and radii 5e-4 apart: 8e-4 in all.
    const Surface b = AboutZ(Point3{3e-4, 0.0, 0.0}, 3.0005);
    std::vector<std::string> log;
    const auto declared =
        Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-3, IntersectTest::Logging(log));
    REQUIRE(declared.has_value());
    REQUIRE(declared->Overlaps.size() == 1);
    REQUIRE(declared->Overlaps.front().Evidence == CoincidenceEvidence::Tolerant);
    REQUIRE_THAT(declared->Overlaps.front().Distance, WithinAbs(8e-4, 1e-12));
    REQUIRE(log.size() == 1);
    // Beyond it, B's axis is inside A and B is the larger: nested, no contact.
    const auto nested = Intersect(a, Box2::Everything(), b, Box2::Everything(), 5e-4, ctx);
    REQUIRE(nested.has_value());
    REQUIRE(nested->IsEmpty());
}

TEST_CASE("coincident cylinders intersect their requested height domains", "[intersect][cylinder]")
{
    const Interval turn{0.0, Stereon::CirclePeriod};
    const Box2 first{turn, Interval{0.0, 1.0}};
    const auto separate = Intersect(a, first, a, Box2{turn, Interval{2.0, 3.0}}, ctx);
    REQUIRE(separate.has_value());
    REQUIRE(separate->IsEmpty());
    const auto same = Intersect(a, first, a, first, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE(same->Overlaps.front().FullDomain);
    REQUIRE(same->Overlaps.front().OnA == first);
    REQUIRE(same->Overlaps.front().OnB == first);

    const auto partial = Intersect(a, first, a, Box2{turn, Interval{0.5, 2.0}}, ctx);
    REQUIRE(partial.has_value());
    REQUIRE(partial->Overlaps.size() == 1);
    const auto& overlap = partial->Overlaps.front();
    REQUIRE_FALSE(overlap.FullDomain);
    REQUIRE_THAT(overlap.OnA.Y.GetLower(), WithinAbs(0.5, 1e-14));
    REQUIRE(overlap.OnA.Y.GetUpper() == 1.0);
    REQUIRE(overlap.OnB.Y.GetLower() == 0.5);
    REQUIRE_THAT(overlap.OnB.Y.GetUpper(), WithinAbs(1.0, 1e-14));

    const Surface down = AboutMinusZ(Point3{0.0, 0.0, 9.0}, 3.0);
    const Box2 reversed{turn, Interval{8.0, 9.0}};
    const auto mapped = Intersect(a, first, down, reversed, ctx);
    REQUIRE(mapped.has_value());
    REQUIRE(mapped->Overlaps.size() == 1);
    REQUIRE(mapped->Overlaps.front().FullDomain);
    REQUIRE(mapped->Overlaps.front().OnB == reversed);
    REQUIRE(Intersect(a, first, down, Box2{turn, Interval{4.0, 5.0}}, ctx)->IsEmpty());
}

TEST_CASE("coincident cylinder arcs split at the periodic seam", "[intersect][cylinder]")
{
    const Interval height{0.0, 1.0};
    const Box2 first{Interval{0.0, 1.0}, height};
    REQUIRE(Intersect(a, first, a, Box2{Interval{2.0, 3.0}, height}, ctx)->IsEmpty());
    const auto same = Intersect(a, first, a, first, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE(same->Overlaps.front().FullDomain);

    const Box2 whole{Interval{0.0, Stereon::CirclePeriod}, height};
    const Box2 across{Interval{-0.5, 0.5}, height};
    const auto split = Intersect(a, whole, a, across, ctx);
    REQUIRE(split.has_value());
    REQUIRE(split->Overlaps.size() == 2);
    double covered = 0.0;
    for (const auto& overlap : split->Overlaps) {
        REQUIRE_FALSE(overlap.FullDomain);
        REQUIRE(Contains(whole.X, overlap.OnA.X));
        REQUIRE(Contains(across.X, overlap.OnB.X));
        covered += Width(overlap.OnA.X);
    }
    REQUIRE_THAT(covered, WithinAbs(1.0, 1e-14));

    // Half-unbounded charts still cover a turn, but their returned chart
    // must stay inside the requested interval.
    const Box2 halfUnbounded{Interval{10.0, std::numeric_limits<double>::infinity()}, height};
    const auto unbounded = Intersect(a, first, a, halfUnbounded, ctx);
    REQUIRE(unbounded.has_value());
    REQUIRE_FALSE(unbounded->IsEmpty());
    for (const auto& overlap : unbounded->Overlaps) {
        REQUIRE_FALSE(overlap.OnB.X.IsEmpty());
        REQUIRE(Contains(halfUnbounded.X, overlap.OnB.X));
    }
}

TEST_CASE("coincident cylinder overlap coverage agrees with an angular membership oracle",
          "[intersect][cylinder][property]")
{
    // Exact quarter-turn frames and many independently sampled chart points.
    // The oracle tests each requested interval modulo the period directly.
    const Surface turned{Cylinder{Point3{}, 3.0, Vec3{0.0, 1.0, 0.0}, Vec3{-1.0, 0.0, 0.0}, true}};
    const Surface reversed = AboutMinusZ(Point3{}, 3.0);
    const auto inside = [](const Interval& arc, const double u) {
        for (int k = -3; k <= 3; ++k) {
            if (Contains(arc, u + k * Stereon::CirclePeriod)) {
                return true;
            }
        }
        return false;
    };
    for (const Surface& other : {a, turned, reversed}) {
        const double sense = other.Target<Cylinder>()->YAxis.Y < 0.0 ? -1.0 : 1.0;
        const double offset = other.Target<Cylinder>()->XAxis.Y > 0.0 ? -pi / 2.0 : 0.0;
        for (const double start : {-1.0, 0.25, 2.0, 5.5}) {
            const Box2 boxA{Interval{-0.75, 5.0}, Interval{-1.0, 1.0}};
            const Box2 boxB{Interval{start, start + 1.5}, Interval{-1.0, 1.0}};
            const auto result = Intersect(a, boxA, other, boxB, ctx);
            REQUIRE(result.has_value());
            for (int sample = 0; sample < 100; ++sample) {
                const double u = -0.75 + 5.75 * (sample + 0.5) / 100.0;
                const bool expected = inside(boxB.X, offset + sense * u);
                const bool found = std::ranges::any_of(result->Overlaps, [&](const auto& overlap) {
                    return Contains(overlap.OnA.X, u);
                });
                REQUIRE(found == expected);
            }
        }
    }
}

TEST_CASE("cylinders touching from outside share one tangent ruling, same side",
          "[intersect][cylinder]")
{
    const Surface b = AboutZ(Point3{6.0, 0.0, 0.0}, 3.0);
    const auto result = Intersect(a, b, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    REQUIRE(branch.Contact == ContactType::Tangential);
    REQUIRE(branch.Transition == SideTransition::SameSide);
    REQUIRE(Distance(branch.Curve3d.Target<Stereon::Line>()->Origin, Point3{3.0, 0.0, 0.0}) == 0.0);
    REQUIRE(branch.OnA.Target<Stereon::Line2d>()->Origin.X == 0.0);
    REQUIRE_THAT(branch.OnB.Target<Stereon::Line2d>()->Origin.X, WithinAbs(pi, 1e-15));
    REQUIRE(PcurveDeviation(branch, a, b) <= 1e-14);
}

TEST_CASE("cylinders touching from inside share one tangent ruling", "[intersect][cylinder]")
{
    // B inside A, touching at A's +x side.
    const Surface smaller = AboutZ(Point3{1.0, 0.0, 0.0}, 2.0);
    const auto inside = Intersect(a, smaller, ctx);
    REQUIRE(inside->Branches.size() == 1);
    REQUIRE(inside->Branches.front().Contact == ContactType::Tangential);
    REQUIRE(Distance(inside->Branches.front().Curve3d.Target<Stereon::Line>()->Origin,
                     Point3{3.0, 0.0, 0.0}) == 0.0);
    REQUIRE(PcurveDeviation(inside->Branches.front(), a, smaller) <= 1e-14);
    // A inside B, touching at A's -x side.
    const Surface larger = AboutZ(Point3{1.0, 0.0, 0.0}, 4.0);
    const auto outside = Intersect(a, larger, ctx);
    REQUIRE(outside->Branches.size() == 1);
    REQUIRE(Distance(outside->Branches.front().Curve3d.Target<Stereon::Line>()->Origin,
                     Point3{-3.0, 0.0, 0.0}) == 0.0);
    REQUIRE(PcurveDeviation(outside->Branches.front(), a, larger) <= 1e-14);
}

TEST_CASE("crossing cylinders share two rulings, in order of angle on the first",
          "[intersect][cylinder]")
{
    const Surface b = AboutZ(Point3{5.0, 0.0, 0.0}, 3.0);
    const auto result = Intersect(a, b, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 2);
    const double y = std::sqrt(9.0 - 6.25);
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.Transition == SideTransition::OppositeSide);
        const Point3 origin = branch.Curve3d.Target<Stereon::Line>()->Origin;
        REQUIRE_THAT(origin.X, WithinAbs(2.5, 1e-14));
        REQUIRE_THAT(std::abs(origin.Y), WithinAbs(y, 1e-14));
        REQUIRE(PcurveDeviation(branch, a, b) <= 1e-14);
    }
    REQUIRE(result->Branches[0].Curve3d.Target<Stereon::Line>()->Origin.Y > 0.0);
    REQUIRE(result->Branches[1].Curve3d.Target<Stereon::Line>()->Origin.Y < 0.0);

    // Heights: A's box from 0 to 2, B's from 1 to 5 (B's origin at the same z).
    const Box2 boxA{Interval{0.0, Stereon::CirclePeriod}, Interval{0.0, 2.0}};
    const Box2 boxB{Interval{0.0, Stereon::CirclePeriod}, Interval{1.0, 5.0}};
    const auto clipped = Intersect(a, boxA, b, boxB, ctx);
    REQUIRE(clipped->Branches.size() == 2);
    REQUIRE(clipped->Branches.front().IsCertified());
    REQUIRE(clipped->Branches.front().Bound < 1e-12);
    REQUIRE_THAT(clipped->Branches.front().Range.GetLower(), WithinAbs(1.0, 1e-14));
    REQUIRE_THAT(clipped->Branches.front().Range.GetUpper(), WithinAbs(2.0, 1e-14));

    // B's height runs the other way: its pcurve descends.
    const Surface down = AboutMinusZ(Point3{5.0, 0.0, 4.0}, 3.0);
    const auto opposite = Intersect(a, down, ctx);
    REQUIRE(opposite->Branches.size() == 2);
    REQUIRE(opposite->Branches.front().OnB.Target<Stereon::Line2d>()->Direction.Y == -1.0);
    REQUIRE(PcurveDeviation(opposite->Branches.front(), a, down) <= 1e-14);
}

TEST_CASE("cylinders apart or nested are proved separate", "[intersect][cylinder]")
{
    REQUIRE(Intersect(a, AboutZ(Point3{10.0, 0.0, 0.0}, 3.0), ctx)->IsEmpty());
    REQUIRE(Intersect(a, AboutZ(Point3{0.5, 0.0, 0.0}, 1.0), ctx)->IsEmpty());
    REQUIRE(Intersect(a, AboutZ(Point3{0.0, 0.0, 7.0}, 1.0), ctx)->IsEmpty());
    // The corpus's near-tangent pair: 5.0005 apart, radii 2 and 3, as numbers disjoint.
    const Surface two = AboutZ(Point3{}, 2.0);
    REQUIRE(Intersect(two, AboutZ(Point3{5.0005, 0.0, 0.0}, 3.0), ctx)->IsEmpty());
}

TEST_CASE("cylinders with skew axes are not implemented yet", "[intersect][cylinder]")
{
    const Surface aboutX{
        Cylinder{Point3{5.0, 0.0, 0.0}, 3.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true}};
    const auto result = Intersect(a, aboutX, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::NotImplemented);
    // Nearly parallel is skew too: the decision is exact (RS-0001).
    const Surface nearly{
        IntersectTest::CylinderAbout(Point3{5.0, 0.0, 0.0}, Vec3{1e-9, 0.0, 1.0}, 3.0)};
    REQUIRE(Intersect(a, nearly, ctx).error().GetCode() == ErrorCode::NotImplemented);
}
