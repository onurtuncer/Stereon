// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <format>

import Stereon.Core;

#include "Formatting.hpp"

using Stereon::Box2;
using Stereon::Box3;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

rc::Gen<double> SmallDyadic()
{
    return rc::gen::map(rc::gen::inRange(-1024, 1025),
                        [](const int value) { return static_cast<double>(value) / 64.0; });
}

rc::Gen<Point3> SmallPoint3()
{
    return rc::gen::construct<Point3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
}

rc::Gen<Box3> SmallBox3()
{
    return rc::gen::map(rc::gen::pair(SmallPoint3(), SmallPoint3()), [](const auto& corners) {
        return Box3::FromCorners(corners.first, corners.second);
    });
}

} // namespace

TEST_CASE("boxes are products of intervals", "[core][box]")
{
    STATIC_REQUIRE(Box3{}.IsEmpty());
    STATIC_REQUIRE(Box3::Empty() == Box3{});
    STATIC_REQUIRE(Box3{Interval{0, 1}, Interval{}, Interval{0, 1}}.IsEmpty());
    STATIC_REQUIRE_FALSE(Box3::Everything().IsEmpty());

    // Corners may be given in any order.
    constexpr Box3 box = Box3::FromCorners(Point3{1, 2, 3}, Point3{0, 5, 2});
    STATIC_REQUIRE(box == Box3{Interval{0, 1}, Interval{2, 5}, Interval{2, 3}});
    STATIC_REQUIRE(box == Box3::FromCorners(Point3{0, 5, 2}, Point3{1, 2, 3}));
    STATIC_REQUIRE(box[0] == box.X);
    STATIC_REQUIRE(box[2] == box.Z);

    STATIC_REQUIRE(MinCorner(box) == Point3{0, 2, 2});
    STATIC_REQUIRE(MaxCorner(box) == Point3{1, 5, 3});
    STATIC_REQUIRE(Center(box) == Point3{0.5, 3.5, 2.5});
    STATIC_REQUIRE(Diagonal(box) == Vec3{1, 3, 1});

    // Boxes are closed: corners and faces are inside.
    STATIC_REQUIRE(Contains(box, Point3{0, 2, 2}));
    STATIC_REQUIRE(Contains(box, Point3{1, 5, 3}));
    STATIC_REQUIRE(Contains(box, Point3{0.5, 3, 2.5}));
    STATIC_REQUIRE_FALSE(Contains(box, Point3{1.5, 3, 2.5}));
    STATIC_REQUIRE_FALSE(Contains(Box3{}, Point3{}));
    STATIC_REQUIRE(Contains(Box3::Everything(), Point3{1e300, -1e300, 0}));

    Box3 mutableBox = box;
    mutableBox[1] = Interval{0, 1};
    REQUIRE(mutableBox.Y == Interval{0, 1});
}

TEST_CASE("hull grows a box around points and boxes", "[core][box]")
{
    constexpr Box3 fromPoints =
        Hull(Hull(Hull(Box3{}, Point3{1, 1, 1}), Point3{-1, 2, 0}), Point3{0, 0, 3});
    STATIC_REQUIRE(fromPoints == Box3::FromCorners(Point3{-1, 0, 0}, Point3{1, 2, 3}));

    constexpr Box3 a = Box3::FromCorners(Point3{0, 0, 0}, Point3{1, 1, 1});
    constexpr Box3 b = Box3::FromCorners(Point3{2, 2, 2}, Point3{3, 3, 3});
    STATIC_REQUIRE(Hull(a, b) == Box3::FromCorners(Point3{0, 0, 0}, Point3{3, 3, 3}));
    STATIC_REQUIRE(Hull(a, Box3{}) == a);
    STATIC_REQUIRE(Hull(Box3{}, b) == b);

    STATIC_REQUIRE(Contains(Hull(a, b), a));
    STATIC_REQUIRE(Contains(Hull(a, b), b));
    STATIC_REQUIRE_FALSE(Contains(a, b));
    STATIC_REQUIRE(Contains(a, Box3{}));
}

TEST_CASE("intersection and overlap tests", "[core][box]")
{
    constexpr Box3 a = Box3::FromCorners(Point3{0, 0, 0}, Point3{2, 2, 2});
    constexpr Box3 b = Box3::FromCorners(Point3{1, 1, 1}, Point3{3, 3, 3});
    STATIC_REQUIRE(Intersection(a, b) == Box3::FromCorners(Point3{1, 1, 1}, Point3{2, 2, 2}));
    STATIC_REQUIRE(Intersects(a, b));

    // Touching faces still intersect; a gap along one axis is enough not to.
    constexpr Box3 touching = Box3::FromCorners(Point3{2, 0, 0}, Point3{3, 2, 2});
    STATIC_REQUIRE(Intersects(a, touching));
    STATIC_REQUIRE(Intersection(a, touching) ==
                   Box3{Interval{2.0}, Interval{0, 2}, Interval{0, 2}});
    constexpr Box3 apart = Box3::FromCorners(Point3{0, 0, 2.5}, Point3{2, 2, 3});
    STATIC_REQUIRE_FALSE(Intersects(a, apart));
    // The result is the canonical empty box, whichever axis separated them.
    STATIC_REQUIRE(Intersection(a, apart) == Box3{});
    STATIC_REQUIRE(Intersection(a, Box3{}) == Box3{});
    STATIC_REQUIRE_FALSE(Intersects(a, Box3{}));
}

TEST_CASE("inflation moves every face", "[core][box]")
{
    constexpr Box3 unit = Box3::FromCorners(Point3{0, 0, 0}, Point3{1, 1, 1});
    STATIC_REQUIRE(Inflated(unit, 0.5) ==
                   Box3::FromCorners(Point3{-0.5, -0.5, -0.5}, Point3{1.5, 1.5, 1.5}));
    STATIC_REQUIRE(Inflated(unit, -0.25) ==
                   Box3::FromCorners(Point3{0.25, 0.25, 0.25}, Point3{0.75, 0.75, 0.75}));
    STATIC_REQUIRE(Inflated(unit, -1.0) == Box3{});
    STATIC_REQUIRE(Inflated(Box3{}, 1.0) == Box3{});
}

TEST_CASE("two-dimensional boxes work the same way", "[core][box]")
{
    constexpr Box2 box = Box2::FromCorners(Point2{3, 0}, Point2{1, 2});
    STATIC_REQUIRE(box == Box2{Interval{1, 3}, Interval{0, 2}});
    STATIC_REQUIRE(Contains(box, Point2{2, 1}));
    STATIC_REQUIRE_FALSE(Contains(box, Point2{0, 1}));
    STATIC_REQUIRE(Center(box) == Point2{2, 1});
    STATIC_REQUIRE(Diagonal(box) == Vec2{2, 2});
    STATIC_REQUIRE(Intersects(box, Box2::FromCorners(Point2{3, 2}, Point2{4, 4})));
    STATIC_REQUIRE(Hull(Box2{}, Point2{1, 1}) == Box2{Interval{1.0}, Interval{1.0}});
    STATIC_REQUIRE(Box2{}.IsEmpty());
}

TEST_CASE("boxes format as corner pairs", "[core][box]")
{
    REQUIRE(std::format("{}", Box3::FromCorners(Point3{0, 0, 0}, Point3{1, 2, 3})) ==
            "[(0, 0, 0), (1, 2, 3)]");
    REQUIRE(std::format("{}", Box2::FromCorners(Point2{1, 2}, Point2{0, 0})) == "[(0, 0), (1, 2)]");
    REQUIRE(std::format("{}", Box3{}) == "[]");
}

TEST_CASE("box set identities hold on dyadic inputs", "[core][box][property]")
{
    rc::prop("corner order does not matter", [] {
        const Point3 p = *SmallPoint3();
        const Point3 q = *SmallPoint3();
        RC_ASSERT(Box3::FromCorners(p, q) == Box3::FromCorners(q, p));
        RC_ASSERT(Contains(Box3::FromCorners(p, q), p));
        RC_ASSERT(Contains(Box3::FromCorners(p, q), q));
    });

    rc::prop("the hull contains its parts and a hulled point", [] {
        const Box3 a = *SmallBox3();
        const Box3 b = *SmallBox3();
        const Point3 p = *SmallPoint3();
        RC_ASSERT(Contains(Hull(a, b), a));
        RC_ASSERT(Contains(Hull(a, b), b));
        RC_ASSERT(Contains(Hull(a, p), p));
        RC_ASSERT(Contains(Hull(a, p), a));
    });

    rc::prop("the intersection lies in both and is empty exactly when they are apart", [] {
        const Box3 a = *SmallBox3();
        const Box3 b = *SmallBox3();
        const Box3 both = Intersection(a, b);
        RC_ASSERT(Contains(a, both));
        RC_ASSERT(Contains(b, both));
        RC_ASSERT(Intersects(a, b) == !both.IsEmpty());
        const Point3 p = *SmallPoint3();
        RC_ASSERT(Contains(both, p) == (Contains(a, p) && Contains(b, p)));
    });
}
