// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// EvalAllFast against Eval: not bit for bit, but within the bound
// docs/nurbs.rst, *Fast batch evaluation*, documents,
//
//     |EvalAllFast - Eval| <= FastBound(p) * eps * S,
//     S = max_i |P_i| * max_i w_i / min_i w_i,
//
// and, like EvalAll, independent of the SIMD width, the parameters' order
// and how they fall into spans, which is tested with ==.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>

import Stereon.Curve;

using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;

namespace {

constexpr double Epsilon = std::numeric_limits<double>::epsilon();

/// The documented bound on |EvalAllFast - Eval| in units of eps * S, by
/// degree (docs/nurbs.rst, *Fast batch evaluation*).
double FastBound(const int degree)
{
    return 8.0 * (degree + 1) + std::ldexp(1.0, degree - 1);
}

/// S = max_i |P_i| * max_i w_i / min_i w_i.
double Scale(const BSplineCurve& c)
{
    double poles = 0.0;
    for (const Point3& pole : c.Poles) {
        poles = std::max(poles, Length(Stereon::ToVec(pole)));
    }
    if (c.Weights.empty()) {
        return poles;
    }
    const auto [low, high] = std::ranges::minmax(c.Weights);
    return poles * high / low;
}

/// A curve of the given degree from `next`, clamped or not, interior knots
/// with multiplicities up to the degree, optionally rational.
template <typename Next>
std::pair<BSplineCurve, bool> MakeCurve(const int p, const bool rational, Next&& next)
{
    const auto q = static_cast<std::size_t>(p);
    const std::size_t poles = q + 1 + static_cast<std::size_t>(next(0, 8));
    const bool clamped = next(0, 2) == 1;
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += next(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(next(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (clamped) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        points.push_back({next(-50, 50) / 7.0, next(-50, 50) / 7.0, next(-50, 50) / 7.0});
        if (rational) {
            weights.push_back(next(2, 40) / 7.0);
        }
    }
    const auto made = Stereon::MakeBSplineCurve(p, knots, points, weights);
    return {made ? *made : BSplineCurve{}, made.has_value()};
}

std::vector<double> SortedParameters(const BSplineCurve& c, const int count)
{
    const Interval d = Domain(c);
    std::vector<double> t;
    for (int k = 0; k <= count; ++k) {
        t.push_back(d.GetLower() + (d.GetUpper() - d.GetLower()) * k / count);
    }
    for (const double u : c.Knots) {
        if (Contains(d, u)) {
            t.push_back(u);
        }
    }
    std::ranges::sort(t);
    return t;
}

/// max_k |fast_k - Eval(c, t_k)| / (eps * S).
double Deviation(const BSplineCurve& c, const std::vector<double>& t,
                 const std::vector<Point3>& fast)
{
    double worst = 0.0;
    for (std::size_t k = 0; k < t.size(); ++k) {
        worst = std::max(worst, Length(fast[k] - Eval(c, t[k])));
    }
    return worst / (Epsilon * Scale(c));
}

} // namespace

TEST_CASE("EvalAllFast stays within its documented bound of Eval",
          "[curve][bspline][batch][property]")
{
    rc::prop("for random curves of every degree, polynomial and rational", [] {
        const int p = *rc::gen::inRange(1, Stereon::BSplineMaxDegree + 1);
        const auto [c, made] = MakeCurve(p, *rc::gen::arbitrary<bool>(), [](int low, int high) {
            return *rc::gen::inRange(low, high);
        });
        RC_PRE(made);
        const std::vector<double> t = SortedParameters(c, *rc::gen::inRange(1, 80));
        std::vector<Point3> fast(t.size());
        Stereon::EvalAllFast(c, t, fast);
        RC_ASSERT(Deviation(c, t, fast) <= FastBound(p));
    });
}

TEST_CASE("EvalAllFast does not depend on the order of its parameters",
          "[curve][bspline][batch][property]")
{
    // Lanes run the scalar Horner steps, so a point is the same whether it
    // was evaluated in a vector (sorted runs) or alone (shuffled ones).
    rc::prop("sorted and shuffled lists give the same bits per parameter", [] {
        const auto [c, made] =
            MakeCurve(*rc::gen::inRange(1, 8), *rc::gen::arbitrary<bool>(),
                      [](int low, int high) { return *rc::gen::inRange(low, high); });
        RC_PRE(made);
        const std::vector<double> sorted = SortedParameters(c, *rc::gen::inRange(1, 60));
        std::vector<std::size_t> order(sorted.size());
        for (std::size_t k = 0; k < order.size(); ++k) {
            order[k] = k;
        }
        std::ranges::reverse(order);
        std::ranges::rotate(order, order.begin() + static_cast<std::ptrdiff_t>(order.size() / 3));
        std::vector<double> shuffled;
        for (const std::size_t k : order) {
            shuffled.push_back(sorted[k]);
        }
        const Stereon::SpanPolynomials polynomials = Stereon::ToSpanPolynomials(c);
        std::vector<Point3> a(sorted.size());
        std::vector<Point3> b(sorted.size());
        Stereon::EvalAllFast(polynomials, sorted, a);
        Stereon::EvalAllFast(polynomials, shuffled, b);
        for (std::size_t k = 0; k < order.size(); ++k) {
            RC_ASSERT(b[k] == a[order[k]]);
        }
        std::vector<Point3> once(sorted.size());
        Stereon::EvalAllFast(c, sorted, once);
        RC_ASSERT(once == a);
    });
}

TEST_CASE("EvalAllFast handles every length and the plane", "[curve][bspline][batch]")
{
    const BSplineCurve c{3,
                         {0.0, 0.0, 0.0, 0.0, 0.3, 0.35, 1.0, 1.0, 1.0, 1.0},
                         {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, -1.0}, Point3{2.5, -1.0, 0.5},
                          Point3{3.0, 1.5, 2.0}, Point3{4.5, 0.0, 1.0}, Point3{5.0, 2.0, -0.5}},
                         {1.0, 2.0, 0.5, 1.5, 3.0, 1.0}};
    BSplineCurve2d flat{c.Degree, c.Knots, {}, c.Weights};
    for (const Point3& pole : c.Poles) {
        flat.Poles.push_back({pole.X, pole.Y});
    }
    for (std::size_t n = 0; n <= 4 * Stereon::Simd::Width<double> + 3; ++n) {
        std::vector<double> t;
        for (std::size_t k = 0; k < n; ++k) {
            t.push_back(static_cast<double>(k) / static_cast<double>(n > 1 ? n - 1 : 1));
        }
        std::vector<Point3> points(n);
        Stereon::EvalAllFast(c, t, points);
        REQUIRE(Deviation(c, t, points) <= FastBound(3));
        std::vector<Point2> planar(n);
        Stereon::EvalAllFast(flat, t, planar);
        for (std::size_t k = 0; k < n; ++k) {
            REQUIRE(planar[k].X == points[k].X);
            REQUIRE(planar[k].Y == points[k].Y);
        }
    }
}

TEST_CASE("EvalAllFast deviation from Eval, measured by degree", "[.][measure][batch]")
{
    // The figures docs/nurbs.rst quotes: run with
    // stereon_curve_tests "[measure]" -s.
    std::mt19937 engine(20261028);
    const auto next = [&engine](const int low, const int high) {
        return std::uniform_int_distribution<int>(low, high - 1)(engine);
    };
    for (int p = 1; p <= Stereon::BSplineMaxDegree; ++p) {
        for (const bool rational : {false, true}) {
            double worst = 0.0;
            double total = 0.0;
            int curves = 0;
            while (curves < (p <= 12 ? 2000 : 300)) {
                const auto [c, made] = MakeCurve(p, rational, next);
                if (!made) {
                    continue;
                }
                const std::vector<double> t = SortedParameters(c, 200);
                std::vector<Point3> fast(t.size());
                Stereon::EvalAllFast(c, t, fast);
                const double deviation = Deviation(c, t, fast);
                worst = std::max(worst, deviation);
                total += deviation;
                ++curves;
            }
            WARN("degree " << p << (rational ? " rational" : " polynomial") << ": max " << worst
                           << ", mean " << total / curves << " (eps * S), bound " << FastBound(p));
        }
    }
}
