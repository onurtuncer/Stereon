// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every ordered pair of built-in plane curves (v0.2's 2D curve–curve gate):
// lines against hyperbolas and parabolas decided exactly, conic pairs by the
// exact quartic of a rational chart, B-splines, polynomial and rational,
// against conics and against each other by proved subdivision. The safety
// properties under test: an empty result is a proved separation (exact
// configurations and their neighbours a rounding either side), a tangency
// is reported, never dropped, and is Decided only when settled exactly,
// what cannot be proved fails with IntersectionFailed, and exchanging the
// curves or reversing one gives the same points.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
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
using Stereon::BSplineCurve2d;
using Stereon::Circle2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::Hyperbola2d;
using Stereon::Interval;
using Stereon::Line2d;
using Stereon::Parabola2d;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Result;
using Stereon::Vec2;

namespace {

constexpr Vec2 ex{1.0, 0.0};
constexpr Vec2 ey{0.0, 1.0};

/// x² - y² = 1, the branch x > 0.
const Hyperbola2d unitHyperbola{Point2{}, 1.0, 1.0, ex, ey};
/// y² = 4x.
const Parabola2d unitParabola{Point2{}, 1.0, ex, ey};

Curve2d LineThrough(const Point2& origin, const Vec2& direction)
{
    return Curve2d{Line2d{origin, direction}};
}

Curve2d Vertical(const double x)
{
    return LineThrough(Point2{x, 0.0}, ey);
}

Curve2d CircleAt(const Point2& center, const double radius)
{
    return Curve2d{Circle2d{center, radius, ex, ey}};
}

BSplineCurve2d Spline(const int degree, std::vector<double> knots, std::vector<Point2> poles,
                      std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// A Bezier curve of degree poles - 1 on [0, 1].
BSplineCurve2d Bezier(std::vector<Point2> poles, std::vector<double> weights = {})
{
    const std::size_t n = poles.size();
    std::vector<double> knots(n, 0.0);
    knots.resize(2 * n, 1.0);
    return Spline(static_cast<int>(n) - 1, std::move(knots), std::move(poles), std::move(weights));
}

/// A quarter of the circle of radius 2 about the origin, from (2, 0) to
/// (0, 2), as a rational quadratic.
BSplineCurve2d QuarterCircle()
{
    return Bezier({{2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}, {1.0, std::numbers::sqrt2 / 2.0, 1.0});
}

/// Every point lies on both curves at its parameters, within `within`, and
/// the points come in order of the first parameter.
void RequireOnBoth(const PlanarIntersection& result, const Curve2d& a, const Curve2d& b,
                   const double within = 1e-9)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const PlanarPoint& p : result.Points) {
        REQUIRE(Distance(Eval(a, p.ParameterA), p.Position) <= within);
        REQUIRE(Distance(Eval(b, p.ParameterB), p.Position) <= within);
        REQUIRE(p.ParameterA >= previous);
        previous = p.ParameterA;
    }
}

/// The intersection, required to succeed, its points on both curves.
PlanarIntersection Meet(const Curve2d& a, const Curve2d& b, const double within = 1e-9)
{
    const Result<PlanarIntersection> result = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireOnBoth(*result, a, b, within);
    return *result;
}

/// The same points from the other side: the parameters exchanged, the
/// contacts and decisions kept.
void RequireSwapped(const PlanarIntersection& ab, const Curve2d& a, const Curve2d& b)
{
    const PlanarIntersection ba = Meet(b, a);
    REQUIRE(ba.Points.size() == ab.Points.size());
    REQUIRE(ba.Overlaps.size() == ab.Overlaps.size());
    for (const PlanarPoint& p : ab.Points) {
        REQUIRE(std::ranges::any_of(ba.Points, [&](const PlanarPoint& q) {
            return Distance(p.Position, q.Position) <= 1e-9 && q.Contact == p.Contact &&
                   q.Decided == p.Decided && std::abs(q.ParameterA - p.ParameterB) <= 1e-9 &&
                   std::abs(q.ParameterB - p.ParameterA) <= 1e-9;
        }));
    }
}

void RequireAll(const PlanarIntersection& result, const ContactType contact, const bool decided)
{
    for (const PlanarPoint& p : result.Points) {
        REQUIRE(p.Contact == contact);
        REQUIRE(p.Decided == decided);
    }
}

} // namespace

// --- Lines against hyperbolas and parabolas -------------------------------

TEST_CASE("a line meets a hyperbola's branch in two points, touches it or misses it, decided "
          "exactly",
          "[intersect][planar][hyperbola]")
{
    const Curve2d hyperbola{unitHyperbola};
    const PlanarIntersection two = Meet(Vertical(2.0), hyperbola);
    REQUIRE(two.Points.size() == 2);
    RequireAll(two, ContactType::Crossing, true);
    REQUIRE_THAT(two.Points[0].ParameterA, WithinAbs(-std::sqrt(3.0), 1e-15));
    REQUIRE_THAT(two.Points[1].ParameterB, WithinAbs(std::acosh(2.0), 1e-15));
    RequireSwapped(two, Vertical(2.0), hyperbola);

    // At the vertex, a tangency; a rounding either side, two points or none.
    const PlanarIntersection touch = Meet(Vertical(1.0), hyperbola);
    REQUIRE(touch.Points.size() == 1);
    RequireAll(touch, ContactType::Tangential, true);
    REQUIRE(touch.Points.front().ParameterB == 0.0);
    REQUIRE(Meet(Vertical(std::nextafter(1.0, 2.0)), hyperbola).Points.size() == 2);
    REQUIRE(Meet(Vertical(std::nextafter(1.0, 0.0)), hyperbola).IsEmpty());

    // The other branch is not the curve: through it, or touching it.
    REQUIRE(Meet(Vertical(-2.0), hyperbola).IsEmpty());
    REQUIRE(Meet(Vertical(-1.0), hyperbola).IsEmpty());

    // Parallel to an asymptote: one crossing (the linear case), at the
    // vertex for y = x - 1; y = x + 1 meets only the other branch; the
    // asymptote itself meets neither.
    const PlanarIntersection once = Meet(LineThrough(Point2{0.0, -1.0}, Vec2{1.0, 1.0}), hyperbola);
    REQUIRE(once.Points.size() == 1);
    RequireAll(once, ContactType::Crossing, true);
    REQUIRE_THAT(once.Points.front().ParameterB, WithinAbs(0.0, 1e-15));
    REQUIRE(Meet(LineThrough(Point2{0.0, 1.0}, Vec2{1.0, 1.0}), hyperbola).IsEmpty());
    REQUIRE(Meet(LineThrough(Point2{}, Vec2{1.0, 1.0}), hyperbola).IsEmpty());

    // The requested domains keep the upper point only.
    const auto upper =
        Intersect(Vertical(2.0), Interval::Everything(), hyperbola, Interval{0.0, 10.0}, ctx);
    REQUIRE(upper.has_value());
    REQUIRE(upper->Points.size() == 1);
    REQUIRE(upper->Points.front().Position.Y > 0.0);
}

TEST_CASE("a line meets a parabola in two points, one, a tangency or none, decided exactly",
          "[intersect][planar][parabola]")
{
    const Curve2d parabola{unitParabola};
    const PlanarIntersection two = Meet(Vertical(1.0), parabola);
    REQUIRE(two.Points.size() == 2);
    RequireAll(two, ContactType::Crossing, true);
    REQUIRE_THAT(two.Points[0].ParameterB, WithinAbs(-2.0, 1e-15));
    REQUIRE_THAT(two.Points[1].ParameterB, WithinAbs(2.0, 1e-15));
    RequireSwapped(two, Vertical(1.0), parabola);

    // Tangent at the vertex and at (1, 2), along y = x + 1; a rounding
    // aside, nothing.
    const PlanarIntersection vertex = Meet(Vertical(0.0), parabola);
    REQUIRE(vertex.Points.size() == 1);
    RequireAll(vertex, ContactType::Tangential, true);
    REQUIRE(Meet(Vertical(-0x1p-60), parabola).IsEmpty());
    const PlanarIntersection slope = Meet(LineThrough(Point2{0.0, 1.0}, Vec2{1.0, 1.0}), parabola);
    REQUIRE(slope.Points.size() == 1);
    RequireAll(slope, ContactType::Tangential, true);
    REQUIRE_THAT(slope.Points.front().ParameterB, WithinAbs(2.0, 1e-15));
    REQUIRE(Meet(LineThrough(Point2{0.0, std::nextafter(1.0, 2.0)}, Vec2{1.0, 1.0}), parabola)
                .IsEmpty());

    // Parallel to the axis: exactly one crossing.
    const PlanarIntersection along = Meet(LineThrough(Point2{0.0, 1.0}, ex), parabola);
    REQUIRE(along.Points.size() == 1);
    RequireAll(along, ContactType::Crossing, true);
    REQUIRE_THAT(along.Points.front().Position.X, WithinAbs(0.25, 1e-15));
}

TEST_CASE("reversing a hyperbola or a parabola negates its parameters",
          "[intersect][planar][metamorphic]")
{
    const Curve2d line = LineThrough(Point2{3.0, -5.0}, Vec2{-0.25, 1.0});
    for (const Curve2d& conic : {Curve2d{unitHyperbola}, Curve2d{unitParabola}}) {
        const PlanarIntersection forward = Meet(line, conic);
        const PlanarIntersection backward = Meet(line, Reversed(conic));
        REQUIRE_FALSE(forward.Points.empty());
        REQUIRE(backward.Points.size() == forward.Points.size());
        for (std::size_t i = 0; i < forward.Points.size(); ++i) {
            REQUIRE_THAT(backward.Points[i].ParameterB,
                         WithinAbs(-forward.Points[i].ParameterB, 1e-12));
        }
    }
}

// --- Two conics -----------------------------------------------------------

TEST_CASE("a hyperbola meets circles: crossings proved, a tangency at a double decided",
          "[intersect][planar][hyperbola]")
{
    const Curve2d hyperbola{unitHyperbola};
    // x² + y² = 7 meets x² - y² = 1 at (2, ±sqrt 3), and the other branch.
    const PlanarIntersection two = Meet(hyperbola, CircleAt(Point2{}, std::sqrt(7.0)));
    REQUIRE(two.Points.size() == 2);
    RequireAll(two, ContactType::Crossing, true);
    for (const PlanarPoint& p : two.Points) {
        REQUIRE_THAT(p.Position.X, WithinAbs(2.0, 1e-12));
    }
    RequireSwapped(two, hyperbola, CircleAt(Point2{}, std::sqrt(7.0)));

    // The unit circle touches the vertex (u = 0 exactly): Decided.
    const PlanarIntersection touch = Meet(hyperbola, CircleAt(Point2{}, 1.0));
    REQUIRE(touch.Points.size() == 1);
    RequireAll(touch, ContactType::Tangential, true);
    REQUIRE(touch.Points.front().ParameterA == 0.0);
    // A rounding smaller, apart; larger, two crossings.
    REQUIRE(Meet(hyperbola, CircleAt(Point2{}, std::nextafter(1.0, 0.0))).IsEmpty());
    const PlanarIntersection near = Meet(hyperbola, CircleAt(Point2{}, std::nextafter(1.0, 2.0)));
    REQUIRE_FALSE(near.Points.empty());
    for (const PlanarPoint& p : near.Points) {
        REQUIRE(Distance(p.Position, Point2{1.0, 0.0}) <= 1e-6);
    }
    // The circle about the other branch's vertex: apart from this one.
    REQUIRE(Meet(hyperbola, CircleAt(Point2{-2.0, 0.0}, 1.0)).IsEmpty());
}

TEST_CASE("conic pairs with hyperbolas and parabolas meet where the algebra says",
          "[intersect][planar][hyperbola][parabola]")
{
    const Curve2d hyperbola{unitHyperbola};
    const Curve2d parabola{unitParabola};
    const Curve2d ellipse{Ellipse2d{Point2{}, 2.0, 1.0, ex, ey}};

    // x²/4 + y² = 1 and y² = 4x: x² + 16 x - 4 = 0, x = -8 + sqrt 68.
    const PlanarIntersection ep = Meet(ellipse, parabola);
    REQUIRE(ep.Points.size() == 2);
    RequireAll(ep, ContactType::Crossing, true);
    for (const PlanarPoint& p : ep.Points) {
        REQUIRE_THAT(p.Position.X, WithinAbs(-8.0 + std::sqrt(68.0), 1e-12));
    }
    RequireSwapped(ep, ellipse, parabola);

    // y² = 4x and x² = 4y (the same parabola turned a quarter) cross at the
    // common vertex, exactly, and at (4, 4).
    const Curve2d turned{Parabola2d{Point2{}, 1.0, ey, Vec2{-1.0, 0.0}}};
    const PlanarIntersection pp = Meet(parabola, turned);
    REQUIRE(pp.Points.size() == 2);
    RequireAll(pp, ContactType::Crossing, true);
    REQUIRE(Distance(pp.Points.front().Position, Point2{}) == 0.0);
    REQUIRE(Distance(pp.Points.back().Position, Point2{4.0, 4.0}) <= 1e-12);
    RequireSwapped(pp, parabola, turned);

    // x² - y² = 1 and x² - y²/4 = 1 touch at the common vertex only.
    const Curve2d wide{Hyperbola2d{Point2{}, 1.0, 2.0, ex, ey}};
    const PlanarIntersection hh = Meet(hyperbola, wide);
    REQUIRE(hh.Points.size() == 1);
    RequireAll(hh, ContactType::Tangential, true);
    // The conjugate y² - x² = 1 never meets it: a proved separation.
    REQUIRE(
        Meet(hyperbola, Curve2d{Hyperbola2d{Point2{}, 1.0, 1.0, ey, Vec2{-1.0, 0.0}}}).IsEmpty());

    // A hyperbola and a parabola: y² = 4x meets x² - y² = 1 where
    // x² - 4x - 1 = 0, x = 2 + sqrt 5.
    const PlanarIntersection hp = Meet(hyperbola, parabola);
    REQUIRE(hp.Points.size() == 2);
    RequireAll(hp, ContactType::Crossing, true);
    for (const PlanarPoint& p : hp.Points) {
        REQUIRE_THAT(p.Position.X, WithinAbs(2.0 + std::sqrt(5.0), 1e-12));
    }
    RequireSwapped(hp, hyperbola, parabola);
    RequireSwapped(Meet(ellipse, hyperbola), ellipse, hyperbola);
}

TEST_CASE("an osculating circle's contact of order four is settled exactly", "[intersect][planar]")
{
    // The ellipse x²/4 + y² = 1 curves with radius 1/2 at (2, 0): the circle
    // about (3/2, 0) through it meets it there only, with multiplicity 4.
    const Curve2d ellipse{Ellipse2d{Point2{}, 2.0, 1.0, ex, ey}};
    const Curve2d osculating = CircleAt(Point2{1.5, 0.0}, 0.5);
    const PlanarIntersection touch = Meet(osculating, ellipse);
    REQUIRE(touch.Points.size() == 1);
    RequireAll(touch, ContactType::Tangential, true);
    REQUIRE(Distance(touch.Points.front().Position, Point2{2.0, 0.0}) == 0.0);
    RequireSwapped(touch, osculating, ellipse);
}

TEST_CASE("coincident conics overlap under their parameters' correspondence", "[intersect][planar]")
{
    // The same hyperbola, the second traversed the other way: t' = -t.
    const Curve2d hyperbola{unitHyperbola};
    const auto same = Intersect(hyperbola, Interval{-1.0, 2.0}, Curve2d{Reversed(unitHyperbola)},
                                Interval{0.0, 5.0}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Points.empty());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE_FALSE(same->Overlaps.front().SameDirection);
    REQUIRE_THAT(same->Overlaps.front().OnA.GetLower(), WithinAbs(-1.0, 1e-12));
    REQUIRE_THAT(same->Overlaps.front().OnA.GetUpper(), WithinAbs(0.0, 1e-12));
    // The other branch is another curve: apart.
    const Curve2d left{Hyperbola2d{Point2{}, 1.0, 1.0, Vec2{-1.0, 0.0}, ey}};
    REQUIRE(Meet(hyperbola, left).IsEmpty());

    // An ellipse and the same ellipse half a turn on: B's angle is s + pi.
    const Curve2d ellipse{Ellipse2d{Point2{}, 2.0, 1.0, ex, ey}};
    const Curve2d half{Ellipse2d{Point2{}, 2.0, 1.0, Vec2{-1.0, 0.0}, Vec2{0.0, -1.0}}};
    const PlanarIntersection turned = Meet(ellipse, half);
    REQUIRE(turned.Points.empty());
    REQUIRE_FALSE(turned.Overlaps.empty());
    for (const Stereon::PlanarOverlap& overlap : turned.Overlaps) {
        REQUIRE(overlap.SameDirection);
    }
    // A parabola with itself.
    const PlanarIntersection parabola = Meet(Curve2d{unitParabola}, Curve2d{unitParabola});
    REQUIRE(parabola.Overlaps.size() == 1);
    REQUIRE(parabola.Overlaps.front().SameDirection);
}

// --- B-splines against conics ---------------------------------------------

TEST_CASE("a rational B-spline meets lines and conics where the circle it draws does",
          "[intersect][planar][spline][rational]")
{
    const Curve2d arc{QuarterCircle()};
    // The diagonal: one proved crossing at (sqrt 2, sqrt 2).
    const PlanarIntersection diagonal = Meet(arc, LineThrough(Point2{}, Vec2{1.0, 1.0}));
    REQUIRE(diagonal.Points.size() == 1);
    RequireAll(diagonal, ContactType::Crossing, true);
    REQUIRE_THAT(diagonal.Points.front().Position.X, WithinAbs(std::numbers::sqrt2, 1e-12));
    RequireSwapped(diagonal, arc, LineThrough(Point2{}, Vec2{1.0, 1.0}));

    // x² - y² = 1 on the circle of radius 2: x² = 5/2, the upper point.
    const PlanarIntersection hyperbola = Meet(arc, Curve2d{unitHyperbola});
    REQUIRE(hyperbola.Points.size() == 1);
    REQUIRE_THAT(hyperbola.Points.front().Position.X, WithinAbs(std::sqrt(2.5), 1e-12));
    RequireSwapped(hyperbola, arc, Curve2d{unitHyperbola});

    // y² = 4x on x² + y² = 4: x = -2 + 2 sqrt 2.
    const PlanarIntersection parabola = Meet(arc, Curve2d{unitParabola});
    REQUIRE(parabola.Points.size() == 1);
    REQUIRE_THAT(parabola.Points.front().Position.X,
                 WithinAbs(-2.0 + 2.0 * std::numbers::sqrt2, 1e-12));
    // The ellipse x²/9 + y² = 1 meets the circle where 8 x² / 9 = 3.
    const Curve2d ellipse{Ellipse2d{Point2{}, 3.0, 1.0, ex, ey}};
    const PlanarIntersection inside = Meet(arc, ellipse);
    REQUIRE(inside.Points.size() == 1);
    REQUIRE_THAT(inside.Points.front().Position.X, WithinAbs(std::sqrt(27.0 / 8.0), 1e-12));

    // Along the circle it draws, nothing can be proved: it fails, never
    // an empty success.
    const auto along = Intersect(arc, CircleAt(Point2{}, 2.0), ctx);
    REQUIRE_FALSE(along.has_value());
    REQUIRE(along.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a B-spline's touch of a line is reported, bracketed, never dropped",
          "[intersect][planar][spline]")
{
    // y = (1 - 2s)², x = 2s: touching the x axis at (1, 0).
    const Curve2d arch{Bezier({{0.0, 1.0}, {1.0, -1.0}, {2.0, 1.0}})};
    const PlanarIntersection touch = Meet(arch, LineThrough(Point2{}, ex), 1e-6);
    REQUIRE(touch.Points.size() == 1);
    REQUIRE(touch.Points.front().Contact == ContactType::Tangential);
    REQUIRE_FALSE(touch.Points.front().Decided);
    REQUIRE(Distance(touch.Points.front().Position, Point2{1.0, 0.0}) <= 1e-6);
    // 3 2^-40 below the axis it is apart; above it, two proved crossings.
    REQUIRE(Meet(arch, LineThrough(Point2{0.0, -3.0 * 0x1p-40}, ex)).IsEmpty());
    const PlanarIntersection two = Meet(arch, LineThrough(Point2{0.0, 3.0 * 0x1p-40}, ex), 1e-6);
    REQUIRE(two.Points.size() == 2);
    RequireAll(two, ContactType::Crossing, true);
    // The parabola y = (x - 1)² / 2 touches it there too.
    const Curve2d cup{Parabola2d{Point2{1.0, 0.0}, 0.5, ey, Vec2{-1.0, 0.0}}};
    const auto both = Intersect(arch, cup, ctx);
    REQUIRE(both.has_value());
    REQUIRE_FALSE(both->Points.empty());
}

TEST_CASE("a B-spline along a line is an overlap", "[intersect][planar][spline]")
{
    const Curve2d along{Bezier({{1.0, 0.0}, {2.0, 0.0}, {4.0, 0.0}})};
    const PlanarIntersection overlap = Meet(along, LineThrough(Point2{}, ex));
    REQUIRE(overlap.Points.empty());
    REQUIRE(overlap.Overlaps.size() == 1);
    REQUIRE(overlap.Overlaps.front().OnA == Interval{0.0, 1.0});
    REQUIRE_THAT(overlap.Overlaps.front().OnB.GetLower(), WithinAbs(1.0, 1e-14));
    REQUIRE_THAT(overlap.Overlaps.front().OnB.GetUpper(), WithinAbs(4.0, 1e-14));
    // Beyond the line's requested stretch, apart; cut by it, not clipped yet.
    REQUIRE(Intersect(along, Interval{0.0, 1.0}, LineThrough(Point2{}, ex), Interval{5.0, 6.0}, ctx)
                ->IsEmpty());
    REQUIRE(Intersect(along, Interval{0.0, 1.0}, LineThrough(Point2{}, ex), Interval{2.0, 6.0}, ctx)
                .error()
                .GetCode() == ErrorCode::NotImplemented);
}

// --- Two B-splines --------------------------------------------------------

TEST_CASE("two B-splines cross where a spline and the line it draws do",
          "[intersect][planar][spline]")
{
    // A cubic wave and a line through it, the line also as a degree-1
    // B-spline with an interior knot: the same three crossings.
    const Curve2d wave{Bezier({{0.0, 0.0}, {1.0, 2.0}, {2.0, -2.0}, {3.0, 0.0}})};
    const Curve2d line = LineThrough(Point2{-1.0, 0.125}, Vec2{1.0, -0.0625});
    const Curve2d segment{
        Spline(1, {0.0, 0.0, 0.5, 1.0, 1.0}, {{-1.0, 0.125}, {1.5, -0.03125}, {4.0, -0.1875}})};
    const PlanarIntersection reference = Meet(wave, line);
    const PlanarIntersection pair = Meet(wave, segment);
    REQUIRE(reference.Points.size() == 3);
    REQUIRE(pair.Points.size() == 3);
    for (std::size_t i = 0; i < pair.Points.size(); ++i) {
        REQUIRE(pair.Points[i].Decided);
        REQUIRE(pair.Points[i].Contact == ContactType::Crossing);
        REQUIRE(Distance(pair.Points[i].Position, reference.Points[i].Position) <= 1e-9);
    }
    RequireSwapped(pair, wave, segment);
    // Reversing the wave keeps the points, its parameters reflected.
    const PlanarIntersection backward = Meet(Reversed(wave), segment);
    REQUIRE(backward.Points.size() == 3);
    REQUIRE_THAT(backward.Points.front().ParameterA,
                 WithinAbs(1.0 - pair.Points.back().ParameterA, 1e-9));
}

TEST_CASE("two rational and polynomial B-splines cross at a piece's middle",
          "[intersect][planar][spline][rational]")
{
    // The diagonal from (0, 0) to (2 sqrt 2, 2 sqrt 2) crosses the quarter
    // circle at its own middle, s = 1/2, where the subdivision cuts.
    const double r = 2.0 * std::numbers::sqrt2;
    const Curve2d diagonal{Bezier({{0.0, 0.0}, {r, r}})};
    const PlanarIntersection one = Meet(Curve2d{QuarterCircle()}, diagonal);
    REQUIRE(one.Points.size() == 1);
    RequireAll(one, ContactType::Crossing, true);
    REQUIRE_THAT(one.Points.front().ParameterB, WithinAbs(0.5, 1e-12));
    RequireSwapped(one, Curve2d{QuarterCircle()}, diagonal);
}

TEST_CASE("two B-splines that touch are a bracketed touch; 3 2^-40 apart they are proved apart",
          "[intersect][planar][spline]")
{
    // The arch y = (x - 1)^2 touches the axis at (1, 0), s = t = 1/2
    // (PlanarSplineTangencyTests.cpp has the full set).
    const Curve2d arch{Bezier({{0.0, 1.0}, {1.0, -1.0}, {2.0, 1.0}})};
    const Curve2d axis{Bezier({{-1.0, 0.0}, {3.0, 0.0}})};
    const PlanarIntersection touch = Meet(arch, axis);
    REQUIRE(touch.Points.size() == 1);
    RequireAll(touch, ContactType::Tangential, false);
    REQUIRE_THAT(touch.Points.front().ParameterA, WithinAbs(0.5, 1e-6));
    REQUIRE_THAT(touch.Points.front().ParameterB, WithinAbs(0.5, 1e-6));
    const Curve2d below{Bezier({{-1.0, -3.0 * 0x1p-40}, {3.0, -3.0 * 0x1p-40}})};
    REQUIRE(Meet(arch, below).IsEmpty());
    // Identical splines overlap on their common requested stretch.
    const auto same = Intersect(arch, Interval{0.25, 1.0}, arch, Interval{0.0, 0.5}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE(same->Overlaps.front().OnA == Interval{0.25, 0.5});
}

// --- Properties -----------------------------------------------------------

TEST_CASE("every sign change of a conic's equation along a random B-spline has a point",
          "[intersect][planar][spline][property]")
{
    // Completeness against sampling, as for curves against surfaces: where
    // the other curve's equation changes sign between two samples of the
    // spline, a reported point lies between them, and an empty success has
    // no sign change at all. A hyperbola's equation changes sign on its
    // other branch too, so only samples on the near side count there.
    constexpr double noise = 1e-9;
    constexpr double slack = 1e-9;
    const auto check = [](const BSplineCurve2d& spline, const Curve2d& other, const auto& level) {
        const auto result = Intersect(Curve2d{spline}, other, ctx);
        if (!result) {
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        const std::vector<double> t = IntersectTest::SamplesOf(Domain(spline), 256);
        for (std::size_t i = 0; i + 1 < t.size(); ++i) {
            const double a = level(Eval(spline, t[i]));
            const double b = level(Eval(spline, t[i + 1]));
            if ((a < -noise && b > noise) || (a > noise && b < -noise)) {
                RC_ASSERT(std::ranges::any_of(result->Points, [&](const PlanarPoint& p) {
                    return t[i] - slack <= p.ParameterA && p.ParameterA <= t[i + 1] + slack;
                }));
            }
        }
        for (const PlanarPoint& p : result->Points) {
            RC_ASSERT(Distance(Eval(spline, p.ParameterA), p.Position) <= 1e-9);
            RC_ASSERT(Distance(Eval(other, p.ParameterB), p.Position) <= 1e-9);
        }
    };
    rc::prop("for random cubic B-splines against a circle, a parabola and a hyperbola", [&] {
        const auto coordinate = [] { return *rc::gen::inRange(-24, 25) / 8.0; };
        std::vector<Point2> poles;
        for (int i = 0; i < 6; ++i) {
            poles.push_back(Point2{coordinate(), coordinate()});
        }
        const BSplineCurve2d spline = *Stereon::MakeBSplineCurve(
            3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0}, poles);
        check(spline, CircleAt(Point2{}, 1.5),
              [](const Point2& p) { return p.X * p.X + p.Y * p.Y - 2.25; });
        check(spline, Curve2d{unitParabola}, [](const Point2& p) { return p.Y * p.Y - 4.0 * p.X; });
        check(spline, Curve2d{unitHyperbola},
              [](const Point2& p) { return p.X > 0.5 ? p.X * p.X - p.Y * p.Y - 1.0 : -1.0; });
    });
}

TEST_CASE("random conic pairs agree with themselves exchanged",
          "[intersect][planar][property][metamorphic]")
{
    // Points lie on both curves, and exchanging the curves exchanges the
    // parameters and keeps the contacts and the decisions.
    rc::prop("for random hyperbolas, parabolas, ellipses and circles", [] {
        const auto value = [](const int lo, const int hi) {
            return *rc::gen::inRange(lo, hi) / 4.0;
        };
        const double angle = value(0, 26);
        const Vec2 x{std::cos(angle), std::sin(angle)};
        const Vec2 y{-x.Y, x.X};
        const Point2 c{value(-8, 9), value(-8, 9)};
        const std::array<Curve2d, 4> curves{
            Curve2d{Hyperbola2d{c, value(1, 9), value(1, 9), x, y}},
            Curve2d{Parabola2d{Point2{value(-8, 9), value(-8, 9)}, value(1, 9), y, -x}},
            Curve2d{
                Ellipse2d{Point2{value(-8, 9), value(-8, 9)}, value(5, 13), value(1, 5), ex, ey}},
            CircleAt(Point2{value(-8, 9), value(-8, 9)}, value(1, 13))};
        for (std::size_t i = 0; i < curves.size(); ++i) {
            for (std::size_t j = i + 1; j < curves.size(); ++j) {
                const auto ab = Intersect(curves[i], curves[j], ctx);
                const auto ba = Intersect(curves[j], curves[i], ctx);
                RC_ASSERT(ab.has_value() == ba.has_value());
                if (!ab) {
                    RC_ASSERT(ab.error().GetCode() == ErrorCode::IntersectionFailed);
                    continue;
                }
                if (ab->Points.size() != ba->Points.size()) {
                    RC_LOG() << "pair " << i << " " << j << ": " << ab->Points.size() << " vs "
                             << ba->Points.size() << "\n";
                    for (const PlanarPoint& p : ab->Points) {
                        RC_LOG() << p.ParameterA << " " << p.ParameterB << " (" << p.Position.X
                                 << ", " << p.Position.Y << ") " << p.Decided << "\n";
                    }
                    for (const PlanarPoint& p : ba->Points) {
                        RC_LOG() << p.ParameterA << " " << p.ParameterB << " (" << p.Position.X
                                 << ", " << p.Position.Y << ") " << p.Decided << "\n";
                    }
                }
                RC_ASSERT(ab->Points.size() == ba->Points.size());
                for (const PlanarPoint& p : ab->Points) {
                    RC_ASSERT(Distance(Eval(curves[i], p.ParameterA), p.Position) <= 1e-7);
                    RC_ASSERT(Distance(Eval(curves[j], p.ParameterB), p.Position) <= 1e-7);
                }
            }
        }
    });
}

TEST_CASE("every ordered pair of built-in plane curves has a handler",
          "[intersect][planar][dispatch]")
{
    const std::array<Curve2d, 6> samples{
        LineThrough(Point2{0.5, -3.0}, Vec2{0.125, 1.0}),
        CircleAt(Point2{0.25, 0.0}, 1.0),
        Curve2d{Ellipse2d{Point2{}, 2.0, 0.75, ex, ey}},
        Curve2d{Bezier({{-2.0, -1.0}, {0.0, 3.0}, {1.0, -2.0}, {3.0, 1.0}})},
        Curve2d{unitHyperbola},
        Curve2d{unitParabola},
    };
    for (const Curve2d& a : samples) {
        for (const Curve2d& b : samples) {
            const auto result = Intersect(a, b, ctx);
            if (!result) {
                // Only a coincidence a B-spline cannot prove may fail.
                REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
                continue;
            }
            RequireOnBoth(*result, a, b);
        }
    }
}
