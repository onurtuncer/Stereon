// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Two plane B-splines crossing at, or within rounding of, an interior knot
// (v0.2 ledger, 2D curve-curve). Each segment's search brackets such a root
// across the segment's end, where the segment's polynomial continued is no
// longer the curve; the handler proves it over the whole spline, piecewise,
// and merges the brackets from either side only when one lies in the box
// the other was proved unique in. The expected answers are closed forms:
// the cubic A below draws x = t exactly (linear precision, its x poles at
// the Greville abscissae) and is antisymmetric about its knot t = 3, so
// A(3) = (3, 0); off the knot, the oracle is an independent de Boor
// evaluation in long double, not the kernel's. Metamorphic relations
// (ADR-0053): exchange, reversal, shape-preserving knot insertion, rigid
// motion and scaling by powers of two.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::BSplineCurve2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Result;

namespace {

BSplineCurve2d Spline(const int degree, std::vector<double> knots, std::vector<Point2> poles,
                      std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// The cubic on [0, 6] with one interior knot, 3: x = t exactly, y odd
/// about t = 3, y'(3) = -1/2. Its x poles are the Greville abscissae.
BSplineCurve2d Wave(std::vector<double> weights = {})
{
    return Spline(3, {0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 6.0, 6.0, 6.0},
                  {{0.0, 0.0}, {1.0, 1.0}, {3.0, 0.0}, {5.0, -1.0}, {6.0, 0.0}},
                  std::move(weights));
}

/// The wave with x and y exchanged about (3, 0): B(u) = (3 + y_A(u), u - 3),
/// through (3, 0) at its own knot u = 3, across A there (A' = (1, -1/2),
/// B' = (-1/2, 1)).
BSplineCurve2d Across()
{
    return Spline(3, {0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 6.0, 6.0, 6.0},
                  {{3.0, -3.0}, {4.0, -2.0}, {3.0, 0.0}, {2.0, 2.0}, {3.0, 3.0}});
}

/// The vertical line x = c from (c, -4) to (c, 4), a degree-1 B-spline
/// with a knot at its middle, (c, 0).
BSplineCurve2d Post(const double c)
{
    return Spline(1, {0.0, 0.0, 0.5, 1.0, 1.0}, {{c, -4.0}, {c, 0.0}, {c, 4.0}});
}

/// The spline with every pole mapped by f.
template <typename F>
BSplineCurve2d Mapped(const BSplineCurve2d& s, F f)
{
    std::vector<Point2> poles;
    for (const Point2& p : s.Poles) {
        poles.push_back(f(p));
    }
    return Spline(s.Degree, s.Knots, std::move(poles), s.Weights);
}

// --- The oracle: de Boor in long double, independent of the kernel -------

struct Exact {
    long double X;
    long double Y;
};

/// The point of a (rational) B-spline at t by de Boor's algorithm on
/// homogeneous poles, in long double.
Exact DeBoor(const BSplineCurve2d& s, const long double t)
{
    const auto p = static_cast<std::size_t>(s.Degree);
    const std::vector<double>& u = s.Knots;
    std::size_t k = p;
    while (k + 1 < s.Poles.size() && static_cast<long double>(u[k + 1]) <= t) {
        ++k;
    }
    struct H {
        long double X;
        long double Y;
        long double W;
    };
    std::vector<H> d;
    for (std::size_t j = 0; j <= p; ++j) {
        const std::size_t i = k - p + j;
        const long double w = s.Weights.empty() ? 1.0L : s.Weights[i];
        d.push_back({w * s.Poles[i].X, w * s.Poles[i].Y, w});
    }
    for (std::size_t r = 1; r <= p; ++r) {
        for (std::size_t j = p; j >= r; --j) {
            const std::size_t i = k - p + j;
            const long double alpha =
                (t - u[i]) / (static_cast<long double>(u[i + p + 1 - r]) - u[i]);
            d[j] = {(1 - alpha) * d[j - 1].X + alpha * d[j].X,
                    (1 - alpha) * d[j - 1].Y + alpha * d[j].Y,
                    (1 - alpha) * d[j - 1].W + alpha * d[j].W};
        }
    }
    return {d[p].X / d[p].W, d[p].Y / d[p].W};
}

// --- Checks ----------------------------------------------------------------

/// The one point of the intersection, required to succeed with exactly one
/// decided crossing.
PlanarPoint OneCrossing(const Curve2d& a, const Curve2d& b)
{
    const Result<PlanarIntersection> result = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    REQUIRE(result->Points.size() == 1);
    const PlanarPoint& p = result->Points.front();
    REQUIRE(p.Contact == ContactType::Crossing);
    REQUIRE(p.Decided);
    return p;
}

void RequireAt(const PlanarPoint& p, const double ta, const double tb, const Point2& at,
               const double within = 1e-13)
{
    REQUIRE_THAT(p.ParameterA, WithinAbs(ta, within));
    REQUIRE_THAT(p.ParameterB, WithinAbs(tb, within));
    REQUIRE_THAT(p.Position.X, WithinAbs(at.X, within * std::max(1.0, std::abs(at.X))));
    REQUIRE_THAT(p.Position.Y, WithinAbs(at.Y, within * std::max(1.0, std::abs(at.Y))));
}

/// The crossing at both knots and its images under the metamorphic
/// relations: the parameters kept (or reflected), the position moved.
void RequireKnotCrossing(const BSplineCurve2d& a, const BSplineCurve2d& b)
{
    const Curve2d ca{a};
    const Curve2d cb{b};
    RequireAt(OneCrossing(ca, cb), 3.0, 3.0, Point2{3.0, 0.0});
    // MR-02: exchanged.
    RequireAt(OneCrossing(cb, ca), 3.0, 3.0, Point2{3.0, 0.0});
    // MR-04: reversed, t -> 6 - t on either or both.
    RequireAt(OneCrossing(Curve2d{Reversed(a)}, cb), 3.0, 3.0, Point2{3.0, 0.0});
    RequireAt(OneCrossing(Curve2d{Reversed(a)}, Curve2d{Reversed(b)}), 3.0, 3.0, Point2{3.0, 0.0});
    // MR-04: knots inserted, the shape kept exactly (dyadic ratios); the
    // crossing knot raised to multiplicity three as well.
    for (const double knot : {1.5, 3.0, 4.5}) {
        auto refined = InsertKnot(a, knot, knot == 3.0 ? 2 : 1);
        REQUIRE(refined.has_value());
        RequireAt(OneCrossing(Curve2d{*refined}, cb), 3.0, 3.0, Point2{3.0, 0.0});
    }
    // MR-01: translated by exact doubles, and turned a quarter.
    const auto shift = [](const Point2& p) { return Point2{p.X + 1024.0, p.Y - 2048.0}; };
    RequireAt(OneCrossing(Curve2d{Mapped(a, shift)}, Curve2d{Mapped(b, shift)}), 3.0, 3.0,
              Point2{1027.0, -2048.0}, 1e-12);
    const auto turn = [](const Point2& p) { return Point2{-p.Y, p.X}; };
    RequireAt(OneCrossing(Curve2d{Mapped(a, turn)}, Curve2d{Mapped(b, turn)}), 3.0, 3.0,
              Point2{0.0, 3.0});
    // MR-03: scaled by 2^20 and 2^-20, the parameters kept.
    for (const double k : {0x1p20, 0x1p-20}) {
        const auto scale = [k](const Point2& p) { return Point2{k * p.X, k * p.Y}; };
        const PlanarPoint p = OneCrossing(Curve2d{Mapped(a, scale)}, Curve2d{Mapped(b, scale)});
        REQUIRE_THAT(p.ParameterA, WithinAbs(3.0, 1e-13));
        REQUIRE_THAT(p.ParameterB, WithinAbs(3.0, 1e-13));
        REQUIRE_THAT(p.Position.X / k, WithinAbs(3.0, 1e-13));
        REQUIRE_THAT(p.Position.Y / k, WithinAbs(0.0, 1e-13));
    }
}

} // namespace

TEST_CASE("two B-splines crossing at both their interior knots are decided",
          "[intersect][planar][spline][knot]")
{
    SECTION("polynomial cubics")
    {
        RequireKnotCrossing(Wave(), Across());
    }
    SECTION("a rational cubic against a polynomial one")
    {
        // Symmetric weights keep A(3) = (3, 0); B is the polynomial cubic
        // crossing there, and the rational wave's slope is still across it.
        const BSplineCurve2d rational = Wave({1.0, 2.0, 1.0, 2.0, 1.0});
        const Exact at = DeBoor(rational, 3.0L);
        REQUIRE(std::abs(at.X - 3.0L) <= 1e-18L);
        REQUIRE(std::abs(at.Y) <= 1e-18L);
        RequireKnotCrossing(rational, Across());
    }
}

TEST_CASE("a B-spline crossing a polyline's vertex at its own knot is decided",
          "[intersect][planar][spline][knot]")
{
    // The post x = 3 has its knot, a vertex, at (3, 0): the crossing is at
    // t = 3 on the wave and u = 1/2 on the post.
    const Curve2d wave{Wave()};
    const Curve2d post{Post(3.0)};
    RequireAt(OneCrossing(wave, post), 3.0, 0.5, Point2{3.0, 0.0});
    RequireAt(OneCrossing(post, wave), 0.5, 3.0, Point2{3.0, 0.0});
    RequireAt(OneCrossing(wave, Curve2d{Reversed(Post(3.0))}), 3.0, 0.5, Point2{3.0, 0.0});
}

TEST_CASE("crossings a rounding and 2^-40 either side of a knot are decided where the oracle "
          "puts them",
          "[intersect][planar][spline][knot]")
{
    const BSplineCurve2d a = Wave();
    const double ulp = std::nextafter(3.0, 4.0) - 3.0;
    for (const double delta : {0.0, ulp, -ulp, 0x1p-40, -0x1p-40, 3.0 * 0x1p-40}) {
        const double c = 3.0 + delta;
        DYNAMIC_SECTION("x = 3 + " << delta)
        {
            // x = t exactly, so the wave meets x = c at t = c; the post's
            // parameter is (y + 4) / 8 there.
            const Exact at = DeBoor(a, static_cast<long double>(c));
            REQUIRE(std::abs(at.X - static_cast<long double>(c)) <= 1e-18L);
            const auto y = static_cast<double>(at.Y);
            const double u = (y + 4.0) / 8.0;
            // A tolerance far below 2^-40: the side of the knot is the oracle's.
            RequireAt(OneCrossing(Curve2d{a}, Curve2d{Post(c)}), c, u, Point2{c, y}, 1e-15);
        }
    }
}

TEST_CASE("a crossing at a corner of a B-spline is reported, never dropped",
          "[intersect][planar][spline][knot]")
{
    // The wave with its knot raised to full multiplicity and its middle
    // pole moved: a corner at (3, 0), entering with slope -1 and leaving
    // with slope -1/4, both across the post.
    const BSplineCurve2d corner = Spline(
        3, {0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 6.0, 6.0, 6.0, 6.0},
        {{0.0, 3.0}, {1.0, 2.0}, {2.0, 1.0}, {3.0, 0.0}, {4.0, -0.25}, {5.0, -0.5}, {6.0, -0.75}});
    const Result<PlanarIntersection> result = Intersect(Curve2d{corner}, Curve2d{Post(3.0)}, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 1);
    RequireAt(result->Points.front(), 3.0, 0.5, Point2{3.0, 0.0});
    REQUIRE(result->Points.front().Decided);
}
