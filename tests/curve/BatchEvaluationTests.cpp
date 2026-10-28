// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// EvalAll against Eval, bit for bit: the batch evaluator promises the very
// doubles Eval returns, whatever the SIMD width, the order of the
// parameters or how they fall into knot spans. Equality is therefore
// tested with ==, never a tolerance, on sorted parameters (the vector
// path), shuffled ones (mostly the scalar path), parameters on knots and at
// the domain ends, and lists of every length around the lane count.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <utility>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;

namespace {

/// A random valid curve: degree 1 to 7, clamped or not, interior knots with
/// multiplicities up to the degree, optionally rational.
BSplineCurve RandomCurve(const bool rational)
{
    const int p = *rc::gen::inRange(1, 8);
    const auto q = static_cast<std::size_t>(p);
    const std::size_t poles = q + 1 + static_cast<std::size_t>(*rc::gen::inRange(0, 8));
    const bool clamped = *rc::gen::arbitrary<bool>();
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (clamped) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        points.push_back({*rc::gen::inRange(-50, 50) / 7.0, *rc::gen::inRange(-50, 50) / 7.0,
                          *rc::gen::inRange(-50, 50) / 7.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(2, 40) / 7.0);
        }
    }
    const auto made = Stereon::MakeBSplineCurve(p, knots, points, weights);
    RC_PRE(made.has_value());
    return *made;
}

/// Parameters over the domain: random interior values, the knots inside
/// it and its two ends, sorted or shuffled.
std::vector<double> RandomParameters(const BSplineCurve& c, const bool sorted)
{
    const Interval d = Domain(c);
    std::vector<double> t;
    const int count = *rc::gen::inRange(0, 60);
    for (int k = 0; k < count; ++k) {
        const double s = *rc::gen::inRange(0, 100000) / 100000.0;
        t.push_back(std::min(d.GetLower() + s * (d.GetUpper() - d.GetLower()), d.GetUpper()));
    }
    for (const double u : c.Knots) {
        if (Contains(d, u) && *rc::gen::inRange(0, 3) == 0) {
            t.push_back(u);
        }
    }
    t.push_back(d.GetLower());
    t.push_back(d.GetUpper());
    if (sorted) {
        std::ranges::sort(t);
    } else {
        const auto order =
            *rc::gen::container<std::vector<int>>(t.size(), rc::gen::inRange(0, 1 << 20));
        std::vector<std::pair<int, double>> keyed;
        for (std::size_t k = 0; k < t.size(); ++k) {
            keyed.emplace_back(order[k], t[k]);
        }
        std::ranges::sort(keyed);
        for (std::size_t k = 0; k < t.size(); ++k) {
            t[k] = keyed[k].second;
        }
    }
    return t;
}

/// Bitwise equality, so that -0.0 and 0.0 differ and no tolerance hides a
/// changed rounding.
template <typename P>
bool SameBits(const P& a, const P& b)
{
    return std::memcmp(&a, &b, sizeof(P)) == 0;
}

BSplineCurve2d Flattened(const BSplineCurve& c)
{
    BSplineCurve2d flat{c.Degree, c.Knots, {}, c.Weights};
    for (const Point3& pole : c.Poles) {
        flat.Poles.push_back({pole.X, pole.Y});
    }
    return flat;
}

} // namespace

TEST_CASE("EvalAll returns Eval's points bit for bit", "[curve][bspline][batch][property]")
{
    rc::prop("in space, sorted or shuffled", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c, *rc::gen::arbitrary<bool>());
        std::vector<Point3> points(t.size());
        Stereon::EvalAll(c, t, points);
        for (std::size_t k = 0; k < t.size(); ++k) {
            RC_ASSERT(SameBits(points[k], Eval(c, t[k])));
        }
    });
    rc::prop("in the plane", [] {
        const BSplineCurve2d c = Flattened(RandomCurve(*rc::gen::arbitrary<bool>()));
        const Interval d = Domain(c);
        std::vector<double> t;
        const int count = *rc::gen::inRange(0, 40);
        for (int k = 0; k <= count; ++k) {
            t.push_back(d.GetLower() + (d.GetUpper() - d.GetLower()) * k / (count + 1));
        }
        std::vector<Point2> points(t.size());
        Stereon::EvalAll(c, t, points);
        for (std::size_t k = 0; k < t.size(); ++k) {
            RC_ASSERT(SameBits(points[k], Eval(c, t[k])));
        }
    });
}

TEST_CASE("EvalAll handles every length around the lane count", "[curve][bspline][batch]")
{
    // A cubic with uneven knots and weights, so that runs cross spans.
    const BSplineCurve c{3,
                         {0.0, 0.0, 0.0, 0.0, 0.3, 0.35, 1.0, 1.0, 1.0, 1.0},
                         {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, -1.0}, Point3{2.5, -1.0, 0.5},
                          Point3{3.0, 1.5, 2.0}, Point3{4.5, 0.0, 1.0}, Point3{5.0, 2.0, -0.5}},
                         {1.0, 2.0, 0.5, 1.5, 3.0, 1.0}};
    for (std::size_t n = 0; n <= 4 * Stereon::Simd::Width<double> + 3; ++n) {
        std::vector<double> t;
        for (std::size_t k = 0; k < n; ++k) {
            t.push_back(static_cast<double>(k) / static_cast<double>(n > 1 ? n - 1 : 1));
        }
        std::vector<Point3> points(n);
        Stereon::EvalAll(c, t, points);
        for (std::size_t k = 0; k < n; ++k) {
            REQUIRE(SameBits(points[k], Eval(c, t[k])));
        }
    }
}
