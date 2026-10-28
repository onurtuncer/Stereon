// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// EvalGrid against Eval, bit for bit, on random polynomial and rational
// surfaces of degrees 1 to 6 and on grids of sorted, shuffled and
// degenerate (empty, single, repeated) parameter lists.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;

namespace {

/// A random non-decreasing knot vector for the degree and pole count,
/// clamped or not, interior multiplicities up to the degree.
std::vector<double> RandomKnots(const int p, const std::size_t poles)
{
    const std::size_t count = poles + static_cast<std::size_t>(p) + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (*rc::gen::arbitrary<bool>()) {
        const auto q = static_cast<std::size_t>(p);
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    return knots;
}

BSplineSurface RandomSurface(const bool rational)
{
    const int p = *rc::gen::inRange(1, 7);
    const int q = *rc::gen::inRange(1, 7);
    const auto countU = static_cast<std::size_t>(p + 1 + *rc::gen::inRange(0, 4));
    const auto countV = static_cast<std::size_t>(q + 1 + *rc::gen::inRange(0, 4));
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t k = 0; k < countU * countV; ++k) {
        poles.push_back({static_cast<double>(k / countV) + *rc::gen::inRange(-5, 5) / 7.0,
                         static_cast<double>(k % countV) + *rc::gen::inRange(-5, 5) / 7.0,
                         *rc::gen::inRange(-30, 30) / 7.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(3, 30) / 7.0);
        }
    }
    const auto made = Stereon::MakeBSplineSurface(p, q, RandomKnots(p, countU),
                                                  RandomKnots(q, countV), poles, weights);
    RC_PRE(made.has_value());
    return *made;
}

/// Parameters over an interval, sorted or shuffled, with its ends.
std::vector<double> RandomParameters(const Interval& d)
{
    std::vector<double> t{d.GetLower(), d.GetUpper()};
    const int count = *rc::gen::inRange(0, 24);
    for (int k = 0; k < count; ++k) {
        const double s = *rc::gen::inRange(0, 10000) / 10000.0;
        t.push_back(std::min(d.GetLower() + s * (d.GetUpper() - d.GetLower()), d.GetUpper()));
    }
    std::ranges::sort(t);
    if (*rc::gen::inRange(0, 3) == 0) {
        std::ranges::reverse(t); // descending: every run still shares spans
    } else if (*rc::gen::inRange(0, 2) == 0) {
        std::ranges::rotate(t, t.begin() + static_cast<std::ptrdiff_t>(t.size() / 2));
    }
    return t;
}

bool SameBits(const Point3& a, const Point3& b)
{
    return std::memcmp(&a, &b, sizeof(Point3)) == 0;
}

void CheckGrid(const BSplineSurface& s, const std::vector<double>& u, const std::vector<double>& v)
{
    std::vector<Point3> points(u.size() * v.size());
    Stereon::EvalGrid(s, u, v, points);
    for (std::size_t i = 0; i < u.size(); ++i) {
        for (std::size_t j = 0; j < v.size(); ++j) {
            RC_ASSERT(SameBits(points[i * v.size() + j], Eval(s, Point2{u[i], v[j]})));
        }
    }
}

} // namespace

TEST_CASE("EvalGrid returns Eval's points bit for bit", "[surface][bspline][grid][property]")
{
    rc::prop("random surfaces and grids", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Box2 d = Domain(s);
        CheckGrid(s, RandomParameters(d.X), RandomParameters(d.Y));
    });
}

TEST_CASE("EvalGrid handles empty, single and repeated parameter lists", "[surface][bspline][grid]")
{
    rc::prop("degenerate grids", [] {
        const BSplineSurface s = RandomSurface(true);
        const Box2 d = Domain(s);
        const std::vector<double> none;
        const std::vector<double> one{d.X.GetLower()};
        const std::vector<double> repeated(9, d.Y.GetUpper());
        CheckGrid(s, none, RandomParameters(d.Y));
        CheckGrid(s, RandomParameters(d.X), none);
        CheckGrid(s, one, repeated);
        CheckGrid(s, std::vector<double>(5, d.X.GetUpper()), RandomParameters(d.Y));
    });
}
