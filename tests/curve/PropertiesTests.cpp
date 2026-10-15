// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local differential properties and arc length of curves, against closed
// forms: the circle, the line, the helix (a user type with torsion), the
// parabola, and the ellipse, whose perimeter is a complete elliptic
// integral of the second kind.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::ArcLength;
using Stereon::BSplineCurve;
using Stereon::CenterOfCurvature;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::Curvature;
using Stereon::Curve;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Parabola;
using Stereon::ParameterAtArcLength;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::PrincipalNormal;
using Stereon::Torsion;
using Stereon::UnitTangent;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-11;
constexpr double Pi = std::numbers::pi;

bool Near(const double a, const double b, const double tolerance = Tolerance)
{
    return std::abs(a - b) <= tolerance * (1.0 + std::abs(b));
}

bool Near(const Vec3& a, const Vec3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(b));
}

/// The circular helix (r cos t, r sin t, c t) on [0, 4 pi], a user curve
/// type: curvature r / (r^2 + c^2), torsion c / (r^2 + c^2).
struct Helix {
    double Radius;
    double Pitch;
};

Point3 Eval(const Helix& h, const double t)
{
    return {h.Radius * std::cos(t), h.Radius * std::sin(t), h.Pitch * t};
}

Vec3 Deriv(const Helix& h, const double t, const int n)
{
    // The n-th derivative of (cos t, sin t) is (cos, sin) of t + n pi/2.
    const double turned = t + n * Pi / 2.0;
    return {h.Radius * std::cos(turned), h.Radius * std::sin(turned), n == 1 ? h.Pitch : 0.0};
}

Interval Domain(const Helix& /*helix*/)
{
    return {0.0, 4.0 * Pi};
}

Circle TiltedCircle()
{
    const double s = std::numbers::sqrt2 / 2.0;
    return {Point3{1.0, -2.0, 0.5}, 3.0, Vec3{s, s, 0.0}, Vec3{0.0, 0.0, 1.0}};
}

} // namespace

TEST_CASE("a circle's properties are its radius and centre", "[curve][properties]")
{
    const Circle circle = TiltedCircle();
    for (const double t : {0.0, 0.3, 2.0, 5.5}) {
        const Vec3 tangent = *UnitTangent(circle, t);
        CHECK(Near(Length(tangent), 1.0));
        CHECK(Near(Dot(tangent, Eval(circle, t) - circle.Center), 0.0));
        CHECK(Near(*Curvature(circle, t), 1.0 / circle.Radius));
        CHECK(Near(*PrincipalNormal(circle, t), (circle.Center - Eval(circle, t)) / circle.Radius));
        CHECK(Near(*Torsion(circle, t), 0.0));
        CHECK(Near(ToVec(*CenterOfCurvature(circle, t)), ToVec(circle.Center)));
    }
}

TEST_CASE("plane curvature is signed by the sense of turning", "[curve][properties]")
{
    const Circle2d left{Point2{1.0, 2.0}, 2.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}};
    const Circle2d right{Point2{1.0, 2.0}, 2.0, Vec2{1.0, 0.0}, Vec2{0.0, -1.0}};
    CHECK(Near(*Curvature(left, 0.7), 0.5));
    CHECK(Near(*Curvature(right, 0.7), -0.5));
}

TEST_CASE("a line has no curvature and no osculating plane", "[curve][properties]")
{
    const Line line{Point3{1.0, 2.0, 3.0}, Vec3{0.0, 3.0, 4.0}};
    CHECK(Near(*UnitTangent(line, 2.0), Vec3{0.0, 0.6, 0.8}));
    CHECK(*Curvature(line, 2.0) == 0.0);
    CHECK(PrincipalNormal(line, 2.0).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(Torsion(line, 2.0).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(CenterOfCurvature(line, 2.0).error().GetCode() == ErrorCode::DegenerateResult);
}

TEST_CASE("a helix has constant curvature and torsion", "[curve][properties]")
{
    const Helix helix{2.0, 0.5};
    const double denominator = helix.Radius * helix.Radius + helix.Pitch * helix.Pitch;
    for (const double t : {0.0, 1.0, 7.0}) {
        CHECK(Near(*Curvature(helix, t), helix.Radius / denominator));
        CHECK(Near(*Torsion(helix, t), helix.Pitch / denominator));
        // The principal normal points horizontally at the axis.
        const Point3 p = Eval(helix, t);
        CHECK(Near(*PrincipalNormal(helix, t), Vec3{-p.X, -p.Y, 0.0} / helix.Radius));
    }
    // Its length is that of the hypotenuse of the unrolled cylinder.
    const auto length = *ArcLength(helix, Domain(helix));
    CHECK(Near(length.Length, 4.0 * Pi * std::sqrt(denominator)));
    CHECK(length.ErrorEstimate <= Tolerance);
}

TEST_CASE("a curve that stops has no tangent there", "[curve][properties]")
{
    // A repeated pole: the speed is zero on the whole first span.
    const BSplineCurve stops{1,
                             {0.0, 0.0, 1.0, 2.0, 2.0},
                             {Point3{0.0, 0.0, 0.0}, Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}},
                             {}};
    CHECK(UnitTangent(stops, 0.5).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(Curvature(stops, 0.5).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(PrincipalNormal(stops, 0.5).error().GetCode() == ErrorCode::DegenerateResult);
    CHECK(Near(*UnitTangent(stops, 1.5), Vec3{1.0, 0.0, 0.0}));
    CHECK(Near(ArcLength(stops, Domain(stops))->Length, 1.0));
    // Any parameter of the stopped stretch has length 0 from its start.
    CHECK(Near(*ParameterAtArcLength(stops, 0.0, 0.5), 1.5));
}

TEST_CASE("arc length matches closed forms", "[curve][properties][arclength]")
{
    SECTION("a parabola, from its vertex")
    {
        const Parabola parabola{Point3{0.0, 0.0, 0.0}, 0.5, Vec3{1.0, 0.0, 0.0},
                                Vec3{0.0, 1.0, 0.0}};
        // The length of (t^2 / 4f, t) over [0, t] is f (u sqrt(1 + u^2) +
        // asinh u) with u = t / 2f.
        const double t = 3.0;
        const double u = t / (2.0 * parabola.FocalLength);
        const double expected = parabola.FocalLength * (u * std::sqrt(1.0 + u * u) + std::asinh(u));
        CHECK(Near(ArcLength(parabola, Interval{0.0, t})->Length, expected));
        CHECK(Near(ArcLength(parabola, Interval{-t, 0.0})->Length, expected));
        CHECK(Near(*Curvature(parabola, 0.0), 1.0 / (2.0 * parabola.FocalLength)));
    }
    SECTION("an ellipse, round its perimeter")
    {
        const Ellipse ellipse{Point3{0.0, 0.0, 0.0}, 5.0, 2.0, Vec3{1.0, 0.0, 0.0},
                              Vec3{0.0, 1.0, 0.0}};
        // 4 a E(e), e the eccentricity, E the complete elliptic integral.
        const double e = std::sqrt(1.0 - 4.0 / 25.0);
        const double perimeter = 4.0 * 5.0 * std::comp_ellint_2(e);
        CHECK(Near(ArcLength(ellipse, Domain(ellipse))->Length, perimeter));
        // Its exact NURBS form has the same length, over its own parameter.
        const BSplineCurve nurbs = *ToBSpline(ellipse, Domain(ellipse));
        CHECK(Near(ArcLength(nurbs, Domain(nurbs))->Length, perimeter));
        // A periodic curve's length counts every turn of the range.
        CHECK(Near(ArcLength(ellipse, Interval{-1.0, 4.0 * Pi - 1.0})->Length, 2.0 * perimeter));
    }
}

TEST_CASE("arc length checks its range", "[curve][properties][arclength]")
{
    const Ellipse ellipse{Point3{0.0, 0.0, 0.0}, 5.0, 2.0, Vec3{1.0, 0.0, 0.0},
                          Vec3{0.0, 1.0, 0.0}};
    const BSplineCurve nurbs = *ToBSpline(ellipse, Interval{0.0, 1.0});
    CHECK(ArcLength(nurbs, Interval{0.0, 2.0}).error().GetCode() == ErrorCode::InvalidInput);
    const Line line{Point3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(ArcLength(line, Interval{0.0, inf}).error().GetCode() == ErrorCode::InvalidInput);
    CHECK(ParameterAtArcLength(nurbs, 0.5, 100.0).error().GetCode() == ErrorCode::InvalidInput);
    CHECK(ParameterAtArcLength(nurbs, 2.0, 0.1).error().GetCode() == ErrorCode::InvalidInput);
    CHECK(ParameterAtArcLength(nurbs, 0.5, inf).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("the parameter at a length inverts arc length", "[curve][properties][arclength]")
{
    SECTION("on a circle, past the end of its domain")
    {
        const Circle circle = TiltedCircle();
        // Two and a half turns forwards from 6 lands half a turn on.
        const double t = *ParameterAtArcLength(circle, 6.0, 2.5 * 2.0 * Pi * circle.Radius);
        CHECK(Near(t, 6.0 + 5.0 * Pi, 1e-10));
        CHECK(Near(*ParameterAtArcLength(circle, 1.0, -circle.Radius), 0.0, 1e-10));
    }
    SECTION("on an unbounded line")
    {
        const Line line{Point3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 2.0}};
        CHECK(Near(*ParameterAtArcLength(line, 1.0, -1000.0), -499.0));
    }
    SECTION("on a rational spline, erased")
    {
        const Ellipse ellipse{Point3{0.0, 0.0, 0.0}, 5.0, 2.0, Vec3{1.0, 0.0, 0.0},
                              Vec3{0.0, 1.0, 0.0}};
        const Curve curve(*ToBSpline(ellipse, Interval{0.2, 5.0}));
        rc::prop("for every pair of parameters", [&curve] {
            const Interval domain = Domain(curve);
            const double span = domain.GetUpper() - domain.GetLower();
            const double a = domain.GetLower() + span * *rc::gen::inRange(0, 1000) / 1000.0;
            const double b = domain.GetLower() + span * *rc::gen::inRange(0, 1000) / 1000.0;
            const double length =
                ArcLength(curve, Interval{std::min(a, b), std::max(a, b)})->Length;
            const double t = *ParameterAtArcLength(curve, a, b >= a ? length : -length);
            RC_ASSERT(std::abs(t - b) <= 1e-9);
        });
    }
}
