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
using Stereon::FrenetFrame;
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

// A straight curve with nonconstant speed; its second derivative is
// nonzero but parallel to its tangent. Analytic derivatives avoid making
// this regression depend on the spline evaluator's rounding.
struct StraightQuadratic {};

Point3 Eval(const StraightQuadratic&, const double t)
{
    const double x = 2.0 * t + t * t;
    return {x, x, x};
}

Vec3 Deriv(const StraightQuadratic&, const double t, const int n)
{
    const double d = n == 1 ? 2.0 + 2.0 * t : (n == 2 ? 2.0 : 0.0);
    return {d, d, d};
}

Circle TiltedCircle()
{
    const double s = std::numbers::sqrt2 / 2.0;
    return {Point3{1.0, -2.0, 0.5}, 3.0, Vec3{s, s, 0.0}, Vec3{0.0, 0.0, 1.0}};
}

/// The width of the Blip's moving stretch, 2^-20: far narrower than the
/// gap between the Gauss nodes of a piece a period wide.
constexpr double BlipWidth = 0x1p-20;

/// A degree-1 curve in the periodic form over [0, 1], period 1: it rises
/// one unit over [0, e] and falls back over [e, 2e], e = BlipWidth, then
/// stands still to the seam, so its length over any period is exactly 2.
/// Its knots are stored for one period only, as the form requires; the
/// release audit of 2026-10-30 (F1) found the quadrature over a later
/// period looking for kinks at the stored knots, finding none in range
/// and returning 0 with an error estimate of 0.
BSplineCurve Blip()
{
    constexpr double e = BlipWidth;
    return *Stereon::MakeBSplineCurve(
        1, std::vector<double>{-1.0 + 2.0 * e, 0.0, e, 2.0 * e, 1.0, 1.0 + e},
        std::vector<Point3>{{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        std::vector<double>{});
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
    CHECK(FrenetFrame(line, 2.0).error().GetCode() == ErrorCode::DegenerateResult);
}

TEST_CASE("a nonuniform straight curve has no principal normal", "[curve][properties]")
{
    const StraightQuadratic curve;
    for (const double t : {0.1, 0.2, 0.7}) {
        const auto normal = PrincipalNormal(curve, t);
        REQUIRE_FALSE(normal.has_value());
        CHECK(normal.error().GetCode() == ErrorCode::DegenerateResult);
        CHECK(*Curvature(curve, t) == 0.0);
        const auto frame = FrenetFrame(curve, t);
        REQUIRE_FALSE(frame.has_value());
        CHECK(frame.error().GetCode() == ErrorCode::DegenerateResult);
    }
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
    // The frame from one evaluation agrees with the separate functions and
    // is a right-handed orthonormal triple.
    const auto frame = *FrenetFrame(helix, 1.0);
    CHECK(Near(frame.Tangent, *UnitTangent(helix, 1.0)));
    CHECK(Near(frame.Normal, *PrincipalNormal(helix, 1.0)));
    CHECK(Near(frame.Binormal, Cross(frame.Tangent, frame.Normal)));
    CHECK(Near(frame.Curvature, *Curvature(helix, 1.0)));
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
        // Also when erased: the quadrature still finds the knots, where the
        // speed has kinks that would fool its error estimate.
        const Curve erased(nurbs);
        CHECK(Near(ArcLength(erased, Domain(erased))->Length, perimeter));
        // A periodic curve's length counts every turn of the range.
        CHECK(Near(ArcLength(ellipse, Interval{-1.0, 4.0 * Pi - 1.0})->Length, 2.0 * perimeter));
    }
}

TEST_CASE("arc length of a periodic B-spline finds the knots of every period",
          "[curve][properties][arclength][periodic]")
{
    const BSplineCurve blip = Blip();
    REQUIRE(Stereon::IsPeriodic(blip));
    REQUIRE(Stereon::Period(blip) == 1.0);
    constexpr double e = BlipWidth;
    // The speed is constant between the knots' images, so the quadrature
    // is exact on every piece once it splits there.
    constexpr double Exact = 1e-9;
    SECTION("one period, from anywhere")
    {
        for (const double start : {0.0, 2.0, -3.0, 0.5, 1.0 - e, 1.0e6 + 0.25}) {
            const auto result = *ArcLength(blip, Interval{start, start + 1.0});
            CHECK(Near(result.Length, 2.0, Exact));
            CHECK(result.ErrorEstimate <= Exact);
        }
    }
    SECTION("many periods")
    {
        CHECK(Near(ArcLength(blip, Interval{0.25, 1000.25})->Length, 2000.0, Exact));
        CHECK(Near(ArcLength(blip, Interval{-7.0, 3.0})->Length, 20.0, Exact));
    }
    SECTION("a narrow range of another period")
    {
        // Half the rise and half the fall.
        CHECK(Near(ArcLength(blip, Interval{2.0 + 0.5 * e, 2.0 + 1.5 * e})->Length, 1.0, Exact));
        CHECK(Near(ArcLength(blip, Interval{-1.0 + 0.5 * e, -1.0 + 1.5 * e})->Length, 1.0, Exact));
        // The still stretch.
        CHECK(ArcLength(blip, Interval{5.0 + 2.0 * e, 6.0})->Length == 0.0);
    }
    SECTION("the same through an erased wrapper")
    {
        const Curve erased(blip);
        CHECK(Near(ArcLength(erased, Interval{2.0, 3.0})->Length, 2.0, Exact));
        CHECK(Near(ArcLength(erased, Interval{2.0 + 0.5 * e, 2.0 + 1.5 * e})->Length, 1.0, Exact));
    }
    SECTION("a trimmed arc across the seam")
    {
        // The blip of the next period, at [1, 1 + 2e], lies inside the arc.
        const auto across = Stereon::MakeTrimmedCurve(Curve(blip), Interval{0.5, 1.5});
        REQUIRE(across.has_value());
        CHECK(Near(ArcLength(*across, Domain(*across))->Length, 2.0, Exact));
        CHECK(Near(ArcLength(Stereon::Reversed(*across), Domain(*across))->Length, 2.0, Exact));
        // Half the rise alone, on an arc that starts in the period before
        // the seam and ends just past it (an arc spans at most one period).
        const auto rise = Stereon::MakeTrimmedCurve(Curve(blip), Interval{-0.5, 0.5 * e});
        REQUIRE(rise.has_value());
        CHECK(Near(ArcLength(*rise, Domain(*rise))->Length, 0.5, Exact));
    }
    SECTION("the parameter at a length follows later periods")
    {
        // The speed is positive where these lengths are reached, so the
        // parameters are unique.
        CHECK(Near(*ParameterAtArcLength(blip, 2.0, 1.0), 2.0 + e, Exact));
        CHECK(Near(*ParameterAtArcLength(blip, 2.0, 2.5), 3.0 + 0.5 * e, Exact));
        CHECK(Near(*ParameterAtArcLength(blip, 0.5, -1.0), e, Exact));
        // Three whole turns and the whole of the next rise.
        CHECK(Near(*ParameterAtArcLength(blip, -4.0, 7.0), -1.0 + e, Exact));
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

TEST_CASE("arc length inversion reports unmet residual targets", "[curve][properties][arclength]")
{
    const Line line{Point3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    SECTION("adjacent parameters cannot represent the requested displacement")
    {
        for (const double length : {0x1p-60, -0x1p-60}) {
            const auto result = ParameterAtArcLength(line, 1.0, length);
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error().GetCode() == ErrorCode::Unresolved);
        }
    }
    SECTION("the iteration budget runs out before the bracket stagnates")
    {
        // The first unit bracket needs more than 60 halvings to reach
        // this start after the rounded Newton step fails to advance.
        const auto result = ParameterAtArcLength(line, 0x1p-100, 0x1p-160);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error().GetCode() == ErrorCode::Unresolved);
    }
    SECTION("small representable displacements and zero still succeed")
    {
        for (const double length : {0.0, 0x1p-40, -0x1p-40}) {
            const auto result = ParameterAtArcLength(line, 1.0, length);
            REQUIRE(result.has_value());
            CHECK(*result == 1.0 + length);
        }
    }
}

TEST_CASE("periodic arc length survives an overflowing parameter width",
          "[curve][properties][regression]")
{
    const Circle circle{Point3{0.0, 0.0, 0.0}, 0.1, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    for (const Interval range : {Interval{-1e308, 1e308}, Interval{-1.7e308, 1e308}}) {
        const double expected = 0.1 * range.GetUpper() - 0.1 * range.GetLower();
        REQUIRE(std::isfinite(expected));
        const auto length = ArcLength(circle, range);
        REQUIRE(length.has_value());
        REQUIRE(std::isfinite(length->Length));
        REQUIRE(std::isfinite(length->ErrorEstimate));
        CHECK(std::abs(length->Length / expected - 1.0) <= Tolerance);
        const auto erased = ArcLength(Curve{circle}, range);
        REQUIRE(erased.has_value());
        CHECK(erased->Length == length->Length);
    }
}

TEST_CASE("arc length refuses an unrepresentable result", "[curve][properties][regression]")
{
    const Circle circle{Point3{0.0, 0.0, 0.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // Exercise both the finite-width and split-width paths.
    for (const Interval range : {Interval{0.0, 1e308}, Interval{-1e308, 1e308}}) {
        const auto length = ArcLength(circle, range);
        REQUIRE_FALSE(length.has_value());
        CHECK(length.error().GetCode() == ErrorCode::InvalidResult);
    }
}
