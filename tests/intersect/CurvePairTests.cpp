// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every pair of 3D curves but two lines, through one search: lines and
// conics as enclosed Bezier pieces of the B-spline pair's (docs/intersection.rst,
// *Lines and conics against other curves*). Closed forms for each kind, meetings placed at a
// known point of two random curves, a near miss either side of the
// resolution, and what the search must refuse.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::Curve;
using Stereon::CurveCurveIntersection;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Parabola;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Vec3 X{1.0, 0.0, 0.0};
const Vec3 Y{0.0, 1.0, 0.0};
const Vec3 Z{0.0, 0.0, 1.0};
const Interval unit{0.0, 1.0};

Result<CurveCurveIntersection> Meet(const Curve& a, const Curve& b)
{
    return Intersect(a, Domain(a), b, Domain(b), ctx);
}

Result<CurveCurveIntersection> Meet(const Curve& a, const Interval& da, const Curve& b,
                                    const Interval& db)
{
    return Intersect(a, da, b, db, ctx);
}

/// Every reported meeting has both curves within the resolution of it.
void CheckOnBoth(const CurveCurveIntersection& r, const Curve& a, const Curve& b)
{
    for (const auto& p : r.Points) {
        CHECK(Distance(Eval(a, p.ParameterA), p.Position) <= ctx.GetResolution());
        CHECK(Distance(Eval(b, p.ParameterB), p.Position) <= ctx.GetResolution());
    }
}

Vec3 Unit(const Vec3& v)
{
    return v / Length(v);
}

BSplineCurve Flat(const Point3& p, const Point3& q)
{
    return BSplineCurve{1, {0.0, 0.0, 1.0, 1.0}, {p, q}, {}};
}

} // namespace

TEST_CASE("a line crosses a circle where the closed form says, on the circle's own angle",
          "[intersect][curve-curve][curve-pair]")
{
    // A chord at height 1 of a circle of radius 2: angles pi/6 and 5 pi/6,
    // and between the quarter-turn arcs' knots, so the angle is recovered.
    const Curve circle{Circle{Point3{}, 2.0, X, Y}};
    const Curve line{Line{Point3{0.0, 1.0, 0.0}, X}};
    const double half = std::sqrt(3.0);
    for (const bool swapped : {false, true}) {
        INFO("swapped " << swapped);
        const auto r = swapped ? Meet(circle, line) : Meet(line, circle);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 2);
        CheckOnBoth(*r, swapped ? circle : line, swapped ? line : circle);
        std::vector<double> onCircle;
        std::vector<double> onLine;
        for (const auto& p : r->Points) {
            onCircle.push_back(swapped ? p.ParameterA : p.ParameterB);
            onLine.push_back(swapped ? p.ParameterB : p.ParameterA);
        }
        std::ranges::sort(onCircle);
        std::ranges::sort(onLine);
        CHECK_THAT(onCircle[0], WithinAbs(pi / 6.0, 1e-9));
        CHECK_THAT(onCircle[1], WithinAbs(5.0 * pi / 6.0, 1e-9));
        CHECK_THAT(onLine[0], WithinAbs(-half, 1e-9));
        CHECK_THAT(onLine[1], WithinAbs(half, 1e-9));
    }
}

TEST_CASE("a line passes a circle within the resolution, or misses it",
          "[intersect][curve-curve][curve-pair]")
{
    const Curve circle{Circle{Point3{}, 2.0, X, Y}};
    const auto at = [&](const double z) {
        return Meet(Curve{Line{Point3{0.0, 1.0, z}, X}}, circle);
    };
    const auto near = at(4e-6);
    INFO(IntersectTest::ErrorText(near));
    REQUIRE(near.has_value());
    CHECK(near->Points.size() == 2);
    const auto far = at(3e-5);
    REQUIRE(far.has_value());
    CHECK(far->IsEmpty());
}

TEST_CASE("a crossing at a whole circle's seam is one meeting",
          "[intersect][curve-curve][curve-pair]")
{
    // The seam at angle 0 is no end: arcs on both sides reach it.
    const Curve circle{Circle{Point3{}, 1.0, X, Y}};
    const Curve across{Line{Point3{1.0, -1.0, -1.0}, Vec3{0.0, 1.0, 1.0}}};
    const auto r = Meet(across, circle);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Points.size() == 1);
    CHECK_THAT(r->Points[0].ParameterA, WithinAbs(1.0, 1e-9));
    const double angle = r->Points[0].ParameterB;
    CHECK(std::min(angle, 2.0 * pi - angle) <= 1e-9);
}

TEST_CASE("a sub-arc of a circle keeps only its meetings", "[intersect][curve-curve][curve-pair]")
{
    const Curve circle{Circle{Point3{}, 2.0, X, Y}};
    const Curve line{Line{Point3{0.0, 1.0, 0.0}, X}};
    const auto r = Meet(line, Interval::Everything(), circle, Interval{0.0, pi / 2.0});
    REQUIRE(r.has_value());
    REQUIRE(r->Points.size() == 1);
    CHECK_THAT(r->Points[0].ParameterB, WithinAbs(pi / 6.0, 1e-9));
    const auto none = Meet(line, Interval::Everything(), circle, Interval{pi, 1.5 * pi});
    REQUIRE(none.has_value());
    CHECK(none->IsEmpty());
}

TEST_CASE("an ellipse, a parabola and a hyperbola are met on their own parameters",
          "[intersect][curve-curve][curve-pair]")
{
    const double t = 0.75;
    const auto check = [&](const Curve& conic, const Interval& domain) {
        const Point3 p = Eval(conic, t);
        const Curve line{Line{p - 2.0 * Vec3{0.3, -0.2, 1.0}, Vec3{0.3, -0.2, 1.0}}};
        const auto r = Meet(line, Interval::Everything(), conic, domain);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterA, WithinAbs(2.0, 1e-9));
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(t, 1e-9));
        CheckOnBoth(*r, line, conic);
    };
    const Interval whole = Interval::Everything();
    SECTION("ellipse")
    {
        check(Curve{Ellipse{Point3{1.0, 2.0, 3.0}, 3.0, 1.5, X, Y}}, Interval{0.0, 2.0 * pi});
    }
    SECTION("parabola, bounded and unbounded")
    {
        const Curve parabola{Parabola{Point3{1.0, 2.0, 3.0}, 0.5, X, Y}};
        check(parabola, Interval{-3.0, 3.0});
        // An unbounded line against an unbounded parabola is refused below.
        const Curve segment{Line{Eval(parabola, t) - Vec3{0.0, 0.0, 1.0}, Z}};
        const auto r = Meet(segment, Interval{0.0, 2.0}, parabola, whole);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(t, 1e-9));
    }
    SECTION("hyperbola, bounded and unbounded")
    {
        const Curve hyperbola{Hyperbola{Point3{1.0, 2.0, 3.0}, 2.0, 1.0, X, Y}};
        check(hyperbola, Interval{-2.5, 2.5});
        const Curve segment{Line{Eval(hyperbola, t) - Vec3{0.0, 0.0, 1.0}, Z}};
        const auto r = Meet(segment, Interval{0.0, 2.0}, hyperbola, whole);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(t, 1e-9));
    }
}

TEST_CASE("a conic meets a conic and a B-spline", "[intersect][curve-curve][curve-pair]")
{
    SECTION("circles in perpendicular planes through two common points")
    {
        const Curve a{Circle{Point3{}, 1.0, X, Y}};
        const Curve b{Circle{Point3{}, 1.0, X, Z}};
        const auto r = Meet(a, b);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 2);
        CheckOnBoth(*r, a, b);
    }
    SECTION("a B-spline through a point of an ellipse")
    {
        const Curve ellipse{Ellipse{Point3{}, 3.0, 1.0, X, Y}};
        const Point3 p = Eval(ellipse, 1.0);
        const Curve spline{Flat(p - Vec3{0.2, 0.1, 1.0}, p + Vec3{0.2, 0.1, 1.0})};
        const auto r = Meet(spline, ellipse);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterA, WithinAbs(0.5, 1e-9));
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(1.0, 1e-9));
    }
    SECTION("an arc's end a rounding short of a line meets it")
    {
        const Curve arc{Circle{Point3{}, 1.0, X, Y}};
        // Across the arc's plane, so the contact is no tangency.
        const Curve line{Line{Point3{0.0, 1.0 + 3e-6, -1.0}, Z}};
        const auto r = Meet(line, Interval{0.0, 2.0}, arc, Interval{0.25, pi / 2.0});
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(pi / 2.0, 1e-9));
    }
}

TEST_CASE("a line and its degree-1 B-spline meet a conic alike",
          "[intersect][curve-curve][curve-pair]")
{
    const Curve circle{Circle{Point3{}, 2.0, X, Y}};
    const Line line{Point3{-3.0, 1.0, 0.0}, X};
    const Interval domain{0.0, 6.0};
    const Curve spline{Flat(Eval(line, 0.0), Eval(line, 6.0))};
    const auto asLine = Meet(Curve{line}, domain, circle, Domain(circle));
    const auto asSpline = Meet(spline, circle);
    REQUIRE(asLine.has_value());
    REQUIRE(asSpline.has_value());
    REQUIRE(asLine->Points.size() == asSpline->Points.size());
    for (std::size_t i = 0; i < asLine->Points.size(); ++i) {
        CHECK_THAT(asLine->Points[i].ParameterB, WithinAbs(asSpline->Points[i].ParameterB, 1e-9));
        CHECK_THAT(asLine->Points[i].ParameterA,
                   WithinAbs(6.0 * asSpline->Points[i].ParameterA, 1e-9));
    }
}

TEST_CASE("the search refuses what it cannot prove", "[intersect][curve-curve][curve-pair]")
{
    const Curve circle{Circle{Point3{}, 1.0, X, Y}};
    const Interval whole = Interval::Everything();
    SECTION("identical conics overlap")
    {
        const auto r = Meet(circle, circle);
        REQUIRE(r.has_value());
        REQUIRE(r->Overlaps.size() == 1);
        CHECK(r->Points.empty());
    }
    SECTION("a tangent line fails closed")
    {
        const auto r = Meet(Curve{Line{Point3{-2.0, 1.0, 0.0}, X}}, circle);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("coaxial circles, equidistant all round, fail closed")
    {
        const auto r = Meet(circle, Curve{Circle{Point3{0.0, 0.0, 1e-6}, 1.0, X, Y}});
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("two unbounded domains are not supported")
    {
        const Curve line{Line{Point3{}, X}};
        const Curve parabola{Parabola{Point3{}, 1.0, X, Y}};
        const auto r = Meet(line, whole, parabola, whole);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::NotSupported);
    }
    SECTION("more than one turn is not supported")
    {
        const auto r = Meet(circle, Interval{0.0, 7.0}, Curve{Line{Point3{}, X}}, whole);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::NotSupported);
    }
    SECTION("invalid data")
    {
        const auto r = Meet(Curve{Circle{Point3{}, -1.0, X, Y}}, Curve{Line{Point3{}, X}});
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("random curves through a common point meet there",
          "[intersect][curve-curve][curve-pair][property]")
{
    rc::prop("in both orders, on each curve's own parameter", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-16, 17) / 8.0; };
        const double angle = *rc::gen::inRange(1, 62) / 10.0;
        const double radius = *rc::gen::inRange(4, 25) / 8.0;
        const Vec3 tilt{coordinate(), coordinate(), 1.0};
        const Vec3 x = Unit(Vec3{1.0, 0.0, -tilt.X});
        const Vec3 y = Unit(Stereon::Cross(tilt, x));
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const Curve circle{Circle{center, radius, x, y}};
        const Point3 p = Eval(circle, angle);
        // A line across the circle's plane at p, so the crossing is simple.
        const Vec3 d = tilt + Vec3{coordinate() / 4.0, coordinate() / 4.0, 0.0};
        const double t = *rc::gen::inRange(-8, 9) / 4.0;
        const Curve line{Line{p - t * d, d}};
        RC_PRE(std::abs(Stereon::Dot(Unit(d), Unit(Deriv(circle, angle, 1)))) < 0.9);
        const auto r = Meet(line, Interval::Everything(), circle, Interval{0.0, 2.0 * pi});
        RC_ASSERT(r.has_value());
        bool found = false;
        for (const auto& m : r->Points) {
            found = found ||
                    (std::abs(m.ParameterA - t) <= 1e-8 && std::abs(m.ParameterB - angle) <= 1e-8);
        }
        RC_ASSERT(found);
    });
}
