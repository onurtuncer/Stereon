// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string_view>
#include <type_traits>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::CirclePeriod;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Ellipse;
using Stereon::Ellipse2d;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

/// A curve type that is not built in: a line under another name.
struct UserLine {
    Line Base;
};
Point3 Eval(const UserLine& c, const double t)
{
    return Eval(c.Base, t);
}
Vec3 Deriv(const UserLine& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const UserLine& c)
{
    return Domain(c.Base);
}
Stereon::Box3 BoundingBox(const UserLine& c)
{
    return BoundingBox(c.Base);
}
double Project(const UserLine& c, const Point3& p)
{
    return Project(c.Base, p);
}
UserLine Reversed(const UserLine& c)
{
    return {Reversed(c.Base)};
}
UserLine Transformed(const UserLine& c, const Transform& x)
{
    return {Transformed(c.Base, x)};
}

/// Which path a visitor was given.
enum class Path { Line, Circle, Ellipse, Erased };

Path PathOf(const Curve& curve)
{
    return VisitBuiltin(curve, []<typename G>(const G&) {
        if constexpr (std::is_same_v<G, Line>) {
            return Path::Line;
        } else if constexpr (std::is_same_v<G, Circle>) {
            return Path::Circle;
        } else if constexpr (std::is_same_v<G, Ellipse>) {
            return Path::Ellipse;
        } else {
            static_assert(std::is_same_v<G, Curve>);
            return Path::Erased;
        }
    });
}

Circle SomeCircle(const double radius)
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, -2.0, 0.5}, 0.3);
    return {Point3{0.5, 1.0, -1.0}, radius, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Ellipse SomeEllipse(const double minorRadius)
{
    const Transform tilt = Stereon::MakeRotation(Vec3{0.5, 1.0, -2.0}, 1.1);
    return {Point3{-0.5, 2.0, 1.0}, minorRadius + 1.25, minorRadius,
            Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

} // namespace

TEST_CASE("VisitBuiltin hands built-in curves to the visitor as themselves", "[curve][visit]")
{
    const Line line{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 1.0, 0.0}};
    const Curve lineCurve(line);
    const Curve circleCurve(SomeCircle(2.0));
    const Curve ellipseCurve(SomeEllipse(0.5));
    const Curve userCurve(UserLine{line});

    REQUIRE(PathOf(lineCurve) == Path::Line);
    REQUIRE(PathOf(circleCurve) == Path::Circle);
    REQUIRE(PathOf(ellipseCurve) == Path::Ellipse);
    REQUIRE(PathOf(userCurve) == Path::Erased);

    REQUIRE(BuiltinIndexOf(lineCurve) == 0);
    REQUIRE(BuiltinIndexOf(circleCurve) == 1);
    REQUIRE(BuiltinIndexOf(ellipseCurve) == 2);
    REQUIRE(BuiltinIndexOf(userCurve) == Stereon::BuiltinCurveTypes::Size);
    STATIC_REQUIRE(Stereon::BuiltinCurveTypes::Size == 3);
    STATIC_REQUIRE(Stereon::GeometryTypeName<Line> == std::string_view{"Line"});
    STATIC_REQUIRE(Stereon::GeometryTypeName<Ellipse> == std::string_view{"Ellipse"});
    STATIC_REQUIRE(Stereon::GeometryTypeName<Ellipse2d> == std::string_view{"Ellipse2d"});
    STATIC_REQUIRE(Stereon::GeometryTypeName<UserLine> == std::string_view{"geometry"});

    // The visitor's result comes back; a loop inside it runs on the type.
    const double sum = VisitBuiltin(circleCurve, [](const auto& g) {
        double total = 0.0;
        for (int i = 0; i < 8; ++i) {
            total += Eval(g, i * CirclePeriod / 8.0).Z;
        }
        return total;
    });
    double expected = 0.0;
    for (int i = 0; i < 8; ++i) {
        expected += Eval(circleCurve, i * CirclePeriod / 8.0).Z;
    }
    REQUIRE(sum == expected);

    const Curve2d planar(Circle2d{Point2{}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}});
    REQUIRE(BuiltinIndexOf(planar) == 1);
    REQUIRE(VisitBuiltin(planar, [](const auto& g) { return Eval(g, 0.0); }) == Point2{1.0, 0.0});
    REQUIRE(BuiltinIndexOf(Curve2d(Line2d{Point2{}, Vec2{1.0, 0.0}})) == 0);
    REQUIRE(BuiltinIndexOf(
                Curve2d(Ellipse2d{Point2{}, 2.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}})) == 2);
}

TEST_CASE("the fast path and the erased path give identical results", "[curve][visit][property]")
{
    // ADR-0003 keeps two dispatch paths and requires them to agree; they
    // run the same code, so they agree to the bit.
    rc::prop("for every built-in curve type and every operation", [] {
        const double t = *rc::gen::inRange(-6283, 6283) / 1000.0;
        const int n = *rc::gen::inRange(1, 5);
        const Point3 p{*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                       *rc::gen::inRange(-50, 50) / 10.0};
        for (const Curve& curve : {Curve(Line{Point3{1.0, 2.0, 3.0}, Vec3{0.5, -1.0, 2.0}}),
                                   Curve(SomeCircle(*rc::gen::inRange(1, 40) / 8.0)),
                                   Curve(SomeEllipse(*rc::gen::inRange(1, 40) / 8.0))}) {
            RC_ASSERT(VisitBuiltin(curve, [&](const auto& g) { return Eval(g, t); }) ==
                      Eval(curve, t));
            RC_ASSERT(VisitBuiltin(curve, [&](const auto& g) { return Deriv(g, t, n); }) ==
                      Deriv(curve, t, n));
            RC_ASSERT(VisitBuiltin(curve, [&](const auto& g) { return Project(g, p); }) ==
                      Project(curve, p));
            RC_ASSERT(VisitBuiltin(curve, [](const auto& g) { return BoundingBox(g); }) ==
                      BoundingBox(curve));
        }
    });
}

TEST_CASE("periodic geometry reports its period through the wrappers", "[curve][periods]")
{
    STATIC_REQUIRE(Stereon::PeriodicCurve<Circle>);
    STATIC_REQUIRE(Stereon::PeriodicCurve<Curve>);
    STATIC_REQUIRE(Stereon::PeriodicCurve<Ellipse>);
    STATIC_REQUIRE_FALSE(Stereon::PeriodicCurve<Line>);
    REQUIRE(Period(Curve(SomeEllipse(1.0))) == CirclePeriod);
    REQUIRE(Period(Curve(SomeCircle(1.0))) == CirclePeriod);
    REQUIRE(Period(Curve(Line{Point3{}, Vec3{1.0, 0.0, 0.0}})) == 0.0);
    // A periodic curve accepts every finite parameter, also through Curve.
    const Curve circle(SomeCircle(2.0));
    REQUIRE(Eval(circle, 3.0 * CirclePeriod).X == Eval(circle, 3.0 * CirclePeriod).X);
    REQUIRE(Stereon::Detail::EvaluatorChecks == static_cast<bool>(STEREON_HOT_CONTRACTS));
}
