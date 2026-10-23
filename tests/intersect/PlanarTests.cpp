// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Plane curve intersection (the Boolean prototype's M6, face clipping):
// lines, circles and ellipses against each other with every configuration
// decided exactly, and B-splines against lines and circles bracketed to
// rounding; parameters placed in the requested domains, points in order of
// the first curve's parameter, overlaps of coincident lines and circles.

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
using Stereon::Circle2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line2d;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Vec2;

namespace {

constexpr double pi = std::numbers::pi;

const Curve2d xAxis{Line2d{Point2{}, Vec2{1.0, 0.0}}};
const Curve2d unitCircle{Circle2d{Point2{}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
/// Semi-axes 2 and 1 along x and y.
const Curve2d ellipse{Ellipse2d{Point2{}, 2.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};

Curve2d Horizontal(const double y)
{
    return Curve2d{Line2d{Point2{0.0, y}, Vec2{1.0, 0.0}}};
}

Curve2d Vertical(const double x)
{
    return Curve2d{Line2d{Point2{x, 0.0}, Vec2{0.0, 1.0}}};
}

/// Every point lies on both curves at its parameters, to rounding or to
/// `within`, and the points come in order of the first parameter.
void RequireOnBoth(const PlanarIntersection& result, const Curve2d& a, const Curve2d& b,
                   const double within = 1e-12)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const PlanarPoint& p : result.Points) {
        REQUIRE(Distance(Eval(a, p.ParameterA), p.Position) <= within);
        REQUIRE(Distance(Eval(b, p.ParameterB), p.Position) <= within);
        REQUIRE(p.ParameterA >= previous);
        previous = p.ParameterA;
    }
}

} // namespace

TEST_CASE("lines cross, pass or coincide, decided exactly", "[intersect][planar]")
{
    const Curve2d upright{Line2d{Point2{1.0, -1.0}, Vec2{0.0, 1.0}}};
    const auto cross = Intersect(xAxis, upright, ctx);
    REQUIRE(cross.has_value());
    REQUIRE(cross->Points.size() == 1);
    REQUIRE(cross->Points.front().ParameterA == 1.0);
    REQUIRE(cross->Points.front().ParameterB == 1.0);
    REQUIRE(cross->Points.front().Contact == ContactType::Crossing);
    REQUIRE(cross->Points.front().Decided);
    RequireOnBoth(*cross, xAxis, upright);

    REQUIRE(Intersect(xAxis, Horizontal(1.0), ctx)->IsEmpty());

    // Coincident: B's parameter s maps to A's t = 2 + 2 s.
    const Curve2d along{Line2d{Point2{2.0, 0.0}, Vec2{2.0, 0.0}}};
    const auto same = Intersect(xAxis, Interval{0.0, 10.0}, along, Interval{0.0, 1.0}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Points.empty());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE(same->Overlaps.front().OnB == Interval{0.0, 1.0});
    REQUIRE_THAT(same->Overlaps.front().OnA.GetLower(), WithinAbs(2.0, 1e-14));
    REQUIRE_THAT(same->Overlaps.front().OnA.GetUpper(), WithinAbs(4.0, 1e-14));
    REQUIRE(same->Overlaps.front().SameDirection);
    const Curve2d back{Line2d{Point2{2.0, 0.0}, Vec2{-1.0, 0.0}}};
    REQUIRE_FALSE(Intersect(xAxis, back, ctx)->Overlaps.front().SameDirection);
    // Beyond A's domain, nothing.
    REQUIRE(Intersect(xAxis, Interval{0.0, 1.0}, along, Interval{0.0, 1.0}, ctx)->IsEmpty());

    // (3, 1) and (0.3, 0.1) are parallel in the reals and not in the
    // doubles: the decision is on the doubles, a crossing far away.
    const Curve2d tilted{Line2d{Point2{0.0, 1.0}, Vec2{3.0, 1.0}}};
    const Curve2d nearly{Line2d{Point2{0.0, 0.0}, Vec2{0.3, 0.1}}};
    const auto far = Intersect(tilted, nearly, ctx);
    REQUIRE(far.has_value());
    REQUIRE(far->Overlaps.empty());
    REQUIRE(far->Points.size() == 1);
    REQUIRE(far->Points.front().Decided);
}

TEST_CASE("a line cuts a circle in two points, touches it in one or misses it",
          "[intersect][planar]")
{
    const auto two = Intersect(Horizontal(0.5), unitCircle, ctx);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    RequireOnBoth(*two, Horizontal(0.5), unitCircle);
    const double x = std::sqrt(0.75);
    REQUIRE_THAT(two->Points[0].ParameterA, WithinAbs(-x, 1e-15));
    REQUIRE_THAT(two->Points[1].ParameterA, WithinAbs(x, 1e-15));
    REQUIRE_THAT(two->Points[1].ParameterB, WithinAbs(pi / 6.0, 1e-15));
    for (const PlanarPoint& p : two->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
        REQUIRE(p.Decided);
    }

    const auto touch = Intersect(unitCircle, Horizontal(1.0), ctx);
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(touch->Points.front().Contact == ContactType::Tangential);
    REQUIRE(touch->Points.front().Decided);
    REQUIRE_THAT(touch->Points.front().ParameterA, WithinAbs(pi / 2.0, 1e-15));
    REQUIRE_THAT(touch->Points.front().ParameterB, WithinAbs(0.0, 1e-15));

    REQUIRE(Intersect(Horizontal(1.5), unitCircle, ctx)->IsEmpty());
    // Just inside and just outside the radius are decided on the doubles.
    REQUIRE(Intersect(Horizontal(std::nextafter(1.0, 0.0)), unitCircle, ctx)->Points.size() == 2);
    REQUIRE(Intersect(Horizontal(std::nextafter(1.0, 2.0)), unitCircle, ctx)->IsEmpty());

    // The circle's requested domain keeps the upper point only, and its
    // parameter is placed in that domain.
    const auto upper =
        Intersect(Vertical(0.0), Interval::Everything(), unitCircle, Interval{0.0, pi}, ctx);
    REQUIRE(upper->Points.size() == 1);
    REQUIRE_THAT(upper->Points.front().ParameterB, WithinAbs(pi / 2.0, 1e-15));
    const auto shifted = Intersect(Vertical(0.0), Interval::Everything(), unitCircle,
                                   Interval{2.0 * pi, 3.0 * pi}, ctx);
    REQUIRE(shifted->Points.size() == 1);
    REQUIRE_THAT(shifted->Points.front().ParameterB, WithinAbs(2.5 * pi, 1e-14));
}

TEST_CASE("circles cross, touch, nest or coincide, decided exactly", "[intersect][planar]")
{
    const auto circleAt = [](const double x, const double r) {
        return Curve2d{Circle2d{Point2{x, 0.0}, r, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    };
    const auto two = Intersect(unitCircle, circleAt(1.0, 1.0), ctx);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    RequireOnBoth(*two, unitCircle, circleAt(1.0, 1.0));
    for (const PlanarPoint& p : two->Points) {
        REQUIRE_THAT(p.Position.X, WithinAbs(0.5, 1e-15));
        REQUIRE(p.Contact == ContactType::Crossing);
    }
    const auto outside = Intersect(unitCircle, circleAt(2.0, 1.0), ctx);
    REQUIRE(outside->Points.size() == 1);
    REQUIRE(outside->Points.front().Contact == ContactType::Tangential);
    REQUIRE(Distance(outside->Points.front().Position, Point2{1.0, 0.0}) == 0.0);
    const auto inside = Intersect(unitCircle, circleAt(0.5, 0.5), ctx);
    REQUIRE(inside->Points.size() == 1);
    REQUIRE(inside->Points.front().Contact == ContactType::Tangential);
    REQUIRE(Distance(inside->Points.front().Position, Point2{1.0, 0.0}) == 0.0);
    REQUIRE(Intersect(unitCircle, circleAt(0.2, 0.5), ctx)->IsEmpty());
    REQUIRE(Intersect(unitCircle, circleAt(3.0, 1.0), ctx)->IsEmpty());
    REQUIRE(Intersect(unitCircle, circleAt(0.0, 0.5), ctx)->IsEmpty());

    // Coincident, B's frame a quarter turn on: B's angle is A's less pi/2,
    // so A's [0, pi] meets B's turn in two stretches.
    const Curve2d turned{Circle2d{Point2{}, 1.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}}};
    const auto same =
        Intersect(unitCircle, Interval{0.0, pi}, turned, Interval{0.0, 2.0 * pi}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Points.empty());
    REQUIRE(same->Overlaps.size() == 2);
    for (const Stereon::PlanarOverlap& overlap : same->Overlaps) {
        REQUIRE(overlap.SameDirection);
        // Not Midpoint: GCC 16.2 rejects its contract here as not constant.
        const double middle = 0.5 * (overlap.OnA.GetLower() + overlap.OnA.GetUpper());
        for (const double t : {overlap.OnA.GetLower(), middle, overlap.OnA.GetUpper()}) {
            const double s = t - pi / 2.0;
            // The stretches are outward-rounded, so an end maps a rounding
            // outside the other stretch.
            const Interval onB = Inflated(overlap.OnB, 1e-12);
            REQUIRE(
                (Contains(onB, s) || Contains(onB, s + 2.0 * pi) || Contains(onB, s - 2.0 * pi)));
            REQUIRE(Distance(Eval(unitCircle, t), Eval(turned, s)) <= 1e-15);
        }
    }
    // Reflected frames run the other way.
    const Curve2d reflected{Circle2d{Point2{}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, -1.0}}};
    REQUIRE_FALSE(Intersect(unitCircle, reflected, ctx)->Overlaps.front().SameDirection);
}

TEST_CASE("a line cuts an ellipse in two points or touches it in one", "[intersect][planar]")
{
    const auto two = Intersect(Vertical(1.0), ellipse, ctx);
    REQUIRE(two.has_value());
    REQUIRE(two->Points.size() == 2);
    RequireOnBoth(*two, Vertical(1.0), ellipse);
    const double y = std::sqrt(0.75);
    REQUIRE_THAT(two->Points[0].ParameterA, WithinAbs(-y, 1e-15));
    REQUIRE_THAT(two->Points[1].ParameterA, WithinAbs(y, 1e-15));
    REQUIRE_THAT(two->Points[1].ParameterB, WithinAbs(pi / 3.0, 1e-14));
    const auto touch = Intersect(ellipse, Vertical(2.0), ctx);
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(touch->Points.front().Contact == ContactType::Tangential);
    REQUIRE(touch->Points.front().Decided);
    REQUIRE_THAT(touch->Points.front().ParameterA, WithinAbs(0.0, 1e-15));
    REQUIRE(Intersect(Vertical(3.0), ellipse, ctx)->IsEmpty());
    REQUIRE(Intersect(Vertical(std::nextafter(2.0, 3.0)), ellipse, ctx)->IsEmpty());
}

TEST_CASE("two conics meet in up to four points, crossings certified exactly",
          "[intersect][planar]")
{
    // The unit circle about (1, 0) cuts the ellipse at x = 2/3 and touches it
    // at (2, 0), where both tangents are vertical.
    const Curve2d circle{Circle2d{Point2{1.0, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    const auto mixed = Intersect(circle, ellipse, ctx);
    INFO(IntersectTest::ErrorText(mixed));
    REQUIRE(mixed.has_value());
    RequireOnBoth(*mixed, circle, ellipse, 1e-9);
    std::size_t crossings = 0;
    for (const PlanarPoint& p : mixed->Points) {
        if (p.Contact == ContactType::Crossing) {
            REQUIRE(p.Decided);
            REQUIRE_THAT(p.Position.X, WithinAbs(2.0 / 3.0, 1e-9));
            ++crossings;
        } else {
            // The touch is bracketed to rounding, not decided; it may also be
            // missed, as a root of even multiplicity.
            REQUIRE_FALSE(p.Decided);
            REQUIRE(Distance(p.Position, Point2{2.0, 0.0}) <= 1e-6);
        }
    }
    REQUIRE(crossings == 2);
    REQUIRE(Intersect(ellipse, circle, ctx)->Points.size() == mixed->Points.size());

    // The same ellipse turned a quarter: four crossings at x = ±y.
    const Curve2d turned{Ellipse2d{Point2{}, 2.0, 1.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}}};
    const auto four = Intersect(ellipse, turned, ctx);
    INFO(IntersectTest::ErrorText(four));
    REQUIRE(four.has_value());
    REQUIRE(four->Points.size() == 4);
    RequireOnBoth(*four, ellipse, turned, 1e-9);
    for (const PlanarPoint& p : four->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
        REQUIRE(p.Decided);
        REQUIRE_THAT(std::abs(p.Position.X), WithinAbs(2.0 / std::sqrt(5.0), 1e-9));
        REQUIRE_THAT(std::abs(p.Position.Y), WithinAbs(2.0 / std::sqrt(5.0), 1e-9));
    }
    // Apart, nothing; nested, nothing.
    const Curve2d far{Circle2d{Point2{5.0, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    REQUIRE(Intersect(far, ellipse, ctx)->IsEmpty());
    const Curve2d small{Circle2d{Point2{}, 0.5, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    REQUIRE(Intersect(small, ellipse, ctx)->IsEmpty());
}

TEST_CASE("a B-spline meets lines and circles to rounding", "[intersect][planar]")
{
    // A cubic through samples of (t, sin t) over a turn.
    std::vector<Point2> samples;
    for (int k = 0; k <= 64; ++k) {
        const double t = 2.0 * pi * k / 64.0;
        samples.push_back(Point2{t, std::sin(t)});
    }
    const auto fit = Interpolate(samples, 3);
    REQUIRE(fit.has_value());
    const auto wave = Reparametrized(fit->Curve, Interval{0.0, 2.0 * pi});
    REQUIRE(wave.has_value());
    const Curve2d spline{*wave};

    const auto half = Intersect(spline, Horizontal(0.5), ctx);
    INFO(IntersectTest::ErrorText(half));
    REQUIRE(half.has_value());
    REQUIRE(half->Points.size() == 2);
    RequireOnBoth(*half, spline, Horizontal(0.5), 1e-9);
    // The fit is parametrised by chord length, so the positions, not the
    // parameters, are where sin t = 1/2.
    REQUIRE_THAT(half->Points[0].Position.X, WithinAbs(pi / 6.0, 1e-4));
    REQUIRE_THAT(half->Points[1].Position.X, WithinAbs(5.0 * pi / 6.0, 1e-4));
    for (const PlanarPoint& p : half->Points) {
        REQUIRE(p.Contact == ContactType::Crossing);
        REQUIRE_FALSE(p.Decided);
    }
    // Through the middle, at a segment join: one point, not doubled. The
    // crossings at the very ends of the spline's domain are roots at the
    // ends of the isolation range, which may be missed (a face's vertices
    // are met by their tolerances, not by roots).
    const auto axis = Intersect(spline, Horizontal(0.0), ctx);
    REQUIRE_FALSE(axis->Points.empty());
    REQUIRE(axis->Points.size() <= 3);
    bool middle = false;
    for (const PlanarPoint& p : axis->Points) {
        middle = middle || std::abs(p.Position.X - pi) <= 1e-6;
    }
    REQUIRE(middle);
    // The line first gives the same points, the sides exchanged.
    const auto swapped = Intersect(Horizontal(0.5), spline, ctx);
    REQUIRE(swapped->Points.size() == 2);
    REQUIRE(swapped->Points[0].ParameterB == half->Points[0].ParameterA);
    REQUIRE(swapped->Points[0].ParameterA == half->Points[0].ParameterB);
    REQUIRE(Intersect(spline, Horizontal(1.5), ctx)->IsEmpty());

    // A circle about (pi, 0): two crossings of the wave.
    const Curve2d circle{Circle2d{Point2{pi, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    const auto round = Intersect(spline, circle, ctx);
    INFO(IntersectTest::ErrorText(round));
    REQUIRE(round.has_value());
    REQUIRE(round->Points.size() == 2);
    RequireOnBoth(*round, spline, circle, 1e-9);
    REQUIRE(Intersect(circle, spline, ctx)->Points.size() == 2);
    // Only the part of the spline before pi.
    REQUIRE(
        Intersect(spline, Interval{0.0, pi}, circle, Interval{0.0, 2.0 * pi}, ctx)->Points.size() ==
        1);
}

TEST_CASE("plane curve pairs dispatch like the others", "[intersect][planar][dispatch]")
{
    STATIC_REQUIRE(Stereon::Detail::PlanarPairTable.size() == Stereon::BuiltinCurve2dTypes::Size);
    const Curve2d hyperbola{
        Stereon::Hyperbola2d{Point2{}, 1.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    const auto stub = Intersect(hyperbola, xAxis, ctx);
    REQUIRE_FALSE(stub.has_value());
    REQUIRE(stub.error().GetCode() == ErrorCode::NotImplemented);
    REQUIRE(std::string(stub.error().GetMessage()).starts_with("Hyperbola2d x Line2d"));
    REQUIRE(Intersect(Curve2d{}, xAxis, ctx).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(
        Intersect(xAxis, Interval::Empty(), xAxis, Interval::Everything(), ctx).error().GetCode() ==
        ErrorCode::InvalidInput);
}
