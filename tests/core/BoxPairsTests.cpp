// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The broad phase must list exactly the pairs that testing every pair lists,
// in the same order: callers keep a fixed order of their pairs (ADR-0007)
// and take this list in its place.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

import Stereon.Core;

using Stereon::Box3;
using Stereon::BoxPair;
using Stereon::Interval;
using Stereon::Point3;

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();

/// Every pair, tested in order.
std::vector<BoxPair> AllPairs(const std::vector<Box3>& a, const std::vector<Box3>& b)
{
    std::vector<BoxPair> pairs;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            if (Intersects(a[i], b[j])) {
                pairs.push_back({i, j});
            }
        }
    }
    return pairs;
}

/// An interval on a coarse lattice, so that boxes often touch exactly; now
/// and then a point, empty, or unbounded on one side or both.
Interval RandomInterval()
{
    const int kind = *rc::gen::inRange(0, 20);
    const double lower = *rc::gen::inRange(-8, 8) * 0.5;
    const double width = *rc::gen::inRange(0, 6) * 0.5;
    switch (kind) {
    case 0:
        return Interval::Empty();
    case 1:
        return {-Infinity, lower};
    case 2:
        return {lower, Infinity};
    case 3:
        return Interval::Everything();
    case 4:
        return Interval{lower};
    default:
        return {lower, lower + width};
    }
}

std::vector<Box3> RandomBoxes()
{
    const int count = *rc::gen::inRange(0, 60);
    std::vector<Box3> boxes;
    for (int k = 0; k < count; ++k) {
        boxes.push_back({RandomInterval(), RandomInterval(), RandomInterval()});
        // Repeats, which split no hierarchy.
        if (*rc::gen::inRange(0, 8) == 0) {
            boxes.push_back(boxes.back());
        }
    }
    return boxes;
}

/// An n x n grid of unit cubes spaced `pitch` apart in x and y.
std::vector<Box3> Grid(const int n, const double pitch, const double shift)
{
    std::vector<Box3> boxes;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const Point3 lo{pitch * i + shift, pitch * j, 0.0};
            boxes.push_back(Box3::FromCorners(lo, Point3{lo.X + 1.0, lo.Y + 1.0, 1.0}));
        }
    }
    return boxes;
}

} // namespace

TEST_CASE("the broad phase lists exactly the pairs that meet, in order", "[core][boxpairs]")
{
    rc::prop("it agrees with testing every pair", [] {
        const std::vector<Box3> a = RandomBoxes();
        const std::vector<Box3> b = RandomBoxes();
        RC_ASSERT(OverlappingPairs(a, b) == AllPairs(a, b));
    });
}

TEST_CASE("the broad phase on lists with nothing to pair", "[core][boxpairs]")
{
    const std::vector<Box3> none;
    const std::vector<Box3> some{Box3::FromCorners(Point3{}, Point3{1.0, 1.0, 1.0})};
    REQUIRE(OverlappingPairs(none, none).empty());
    REQUIRE(OverlappingPairs(none, some).empty());
    REQUIRE(OverlappingPairs(some, none).empty());
    const std::vector<Box3> empty{Box3::Empty(), Box3::Empty()};
    REQUIRE(OverlappingPairs(empty, some).empty());
    REQUIRE(OverlappingPairs(some, empty).empty());
}

TEST_CASE("boxes that touch at a corner meet; boxes a step apart do not", "[core][boxpairs]")
{
    const std::vector<Box3> a{Box3::FromCorners(Point3{}, Point3{1.0, 1.0, 1.0})};
    const std::vector<Box3> corner{Box3::FromCorners(Point3{1.0, 1.0, 1.0}, Point3{2.0, 2.0, 2.0})};
    REQUIRE(OverlappingPairs(a, corner) == std::vector<BoxPair>{{0, 0}});
    const double apart = std::nextafter(1.0, 2.0);
    const std::vector<Box3> gap{Box3::FromCorners(Point3{apart, 0.0, 0.0}, Point3{2.0, 1.0, 1.0})};
    REQUIRE(OverlappingPairs(a, gap).empty());
}

TEST_CASE("the broad phase on large, mostly disjoint grids", "[core][boxpairs]")
{
    // Two 60 x 60 grids of unit cubes two apart, the second shifted by half
    // a unit in x: each cube meets only its counterpart, 3 600 pairs of the
    // 13 million.
    const std::vector<Box3> a = Grid(60, 2.0, 0.0);
    const std::vector<Box3> b = Grid(60, 2.0, 0.5);
    const std::vector<BoxPair> pairs = OverlappingPairs(a, b);
    REQUIRE(pairs == AllPairs(a, b));
    REQUIRE(pairs.size() == 60 * 60);
}
