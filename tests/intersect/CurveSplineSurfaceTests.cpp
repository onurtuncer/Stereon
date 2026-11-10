// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Curves against B-spline surfaces (docs/intersection.rst, *Curves against
// B-spline surfaces*). The saddle z = u^2 - v^2 is a degree-2 B-spline
// exactly, so a line through a known point of it is a closed form; a NURBS
// cylinder is met where the analytic cylinder's own handler meets it; and
// tangent contact, a line along the surface and roots at a domain's end fail
// closed.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double lo = -1.0;
constexpr double hi = 2.0;

/// z = u^2 - v^2 over [lo, hi]^2, x = u, y = v: degree 2 in each
/// direction, its poles the Bernstein coefficients of u, v, u^2 and v^2.
BSplineSurface Saddle()
{
    const std::vector<double> knots{lo, lo, lo, hi, hi, hi};
    const std::array<double, 3> line{lo, 0.5 * (lo + hi), hi};
    const std::array<double, 3> square{lo * lo, lo * hi, hi * hi};
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({line[i], line[j], square[i] - square[j]});
        }
    }
    return *Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
}

Result<CurveSurfaceIntersection> Meet(const Curve& c, const Surface& s)
{
    return Intersect(c, Domain(c), s, Domain(s), ctx);
}

/// Every point lies on both, within the resolution.
void CheckOnBoth(const CurveSurfaceIntersection& r, const Curve& c, const Surface& s)
{
    for (const auto& p : r.Points) {
        CHECK(Distance(Eval(c, p.Parameter), p.Position) <= ctx.GetResolution());
        CHECK(Distance(Eval(s, p.OnSurface), p.Position) <= ctx.GetResolution());
        CHECK(p.Decided);
        CHECK(p.Contact == Stereon::ContactType::Crossing);
    }
}

bool HasPoint(const CurveSurfaceIntersection& r, const double t, const Point2& uv,
              const double within)
{
    return std::ranges::any_of(r.Points, [&](const auto& p) {
        return std::abs(p.Parameter - t) <= within && std::abs(p.OnSurface.X - uv.X) <= within &&
               std::abs(p.OnSurface.Y - uv.Y) <= within;
    });
}

} // namespace

TEST_CASE("a line crosses a B-spline saddle where the closed form says",
          "[intersect][curve-surface][spline-surface]")
{
    const Surface saddle{Saddle()};
    const double p = 0.75;
    const double q = 1.25;
    // Vertical: z = p^2 - q^2 at t = that height plus 3.
    const Curve line{Line{Point3{p, q, -3.0}, Vec3{0.0, 0.0, 1.0}}};
    const auto r = Meet(line, saddle);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Points.size() == 1);
    CHECK(r->Overlaps.empty());
    CheckOnBoth(*r, line, saddle);
    CHECK_THAT(r->Points[0].Parameter, WithinAbs(p * p - q * q + 3.0, 1e-10));
    CHECK_THAT(r->Points[0].OnSurface.X, WithinAbs(p, 1e-10));
    CHECK_THAT(r->Points[0].OnSurface.Y, WithinAbs(q, 1e-10));
}

TEST_CASE("lines through random points of the saddle meet it there",
          "[intersect][curve-surface][spline-surface][property]")
{
    const Surface saddle{Saddle()};
    rc::prop("at the point's parameters, in the line's own", [&] {
        const double p = *rc::gen::inRange(-15, 32) / 16.0;
        const double q = *rc::gen::inRange(-15, 32) / 16.0;
        const Vec3 d{*rc::gen::inRange(-4, 5) / 8.0, *rc::gen::inRange(-4, 5) / 8.0, 1.0};
        // Transverse: the direction well off the tangent plane, whose normal
        // is (-2p, 2q, 1).
        const double across = Dot(d, Vec3{-2.0 * p, 2.0 * q, 1.0}) /
                              (Length(d) * Length(Vec3{-2.0 * p, 2.0 * q, 1.0}));
        RC_PRE(std::abs(across) > 0.2);
        const double t = *rc::gen::inRange(-8, 9) / 4.0;
        const Point3 at{p, q, p * p - q * q};
        const Curve line{Line{at - t * d, d}};
        const auto r = Meet(line, saddle);
        RC_ASSERT(r.has_value());
        RC_ASSERT(HasPoint(*r, t, Point2{p, q}, 1e-8));
        for (const auto& m : r->Points) {
            RC_ASSERT(Distance(Eval(saddle, m.OnSurface), m.Position) <= ctx.GetResolution());
        }
    });
}

TEST_CASE("a NURBS cylinder is crossed where the analytic cylinder is",
          "[intersect][curve-surface][spline-surface]")
{
    const Cylinder cylinder{Point3{1.0, -2.0, 0.5}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const auto nurbs =
        Stereon::ToBSpline(cylinder, Box2{Interval{0.0, 2.0 * pi}, Interval{-20.0, 20.0}});
    REQUIRE(nurbs.has_value());
    const Surface rational{*nurbs};
    for (const Line& line : {Line{Point3{-5.0, -1.5, 0.0}, Vec3{1.0, 0.1, 0.3}},
                             Line{Point3{1.0, -2.0, -3.0}, Vec3{0.2, 1.0, 0.7}},
                             Line{Point3{3.5, 0.0, 1.0}, Vec3{-0.3, -1.0, 0.0}}}) {
        const auto exact = Intersect(Curve{line}, Surface{cylinder}, ctx);
        const auto ours = Meet(Curve{line}, rational);
        INFO(IntersectTest::ErrorText(ours));
        REQUIRE(exact.has_value());
        REQUIRE(ours.has_value());
        REQUIRE(ours->Points.size() == exact->Points.size());
        CheckOnBoth(*ours, Curve{line}, rational);
        for (std::size_t i = 0; i < ours->Points.size(); ++i) {
            CHECK_THAT(ours->Points[i].Parameter, WithinAbs(exact->Points[i].Parameter, 1e-9));
            CHECK(Distance(ours->Points[i].Position, exact->Points[i].Position) <= 1e-9);
        }
    }
}

TEST_CASE("conics and B-spline curves cross a B-spline surface",
          "[intersect][curve-surface][spline-surface]")
{
    const Surface saddle{Saddle()};
    SECTION("a circle standing across the saddle, crossing it twice")
    {
        // In the plane y = 0.5, about (0.5, 0.5, 0): it meets z = x^2 - 0.25.
        const Curve circle{
            Circle{Point3{0.5, 0.5, 0.0}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};
        const auto r = Meet(circle, saddle);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 2);
        CheckOnBoth(*r, circle, saddle);
    }
    SECTION("a cubic through a known point")
    {
        const Point3 at{0.5, -0.25, 0.25 - 0.0625};
        const Vec3 d{0.2, 0.1, 1.0};
        const Curve cubic{BSplineCurve{3,
                                       {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
                                       {at - d, at - d / 3.0 + Vec3{0.05, 0.0, 0.0},
                                        at + d / 3.0 - Vec3{0.05, 0.0, 0.0}, at + d},
                                       {}}};
        const auto r = Meet(cubic, saddle);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CheckOnBoth(*r, cubic, saddle);
        CHECK_THAT(r->Points[0].Parameter, WithinAbs(0.5, 1e-9));
        CHECK_THAT(r->Points[0].OnSurface.X, WithinAbs(0.5, 1e-9));
        CHECK_THAT(r->Points[0].OnSurface.Y, WithinAbs(-0.25, 1e-9));
    }
}

TEST_CASE("requested domains keep only their crossings",
          "[intersect][curve-surface][spline-surface]")
{
    const Surface saddle{Saddle()};
    const Curve line{Line{Point3{0.75, 1.25, -3.0}, Vec3{0.0, 0.0, 1.0}}};
    const auto outside = Intersect(line, Interval::Everything(), saddle,
                                   Box2{Interval{lo, 0.5}, Interval{lo, hi}}, ctx);
    REQUIRE(outside.has_value());
    CHECK(outside->Points.empty());
    const auto shortLine = Intersect(line, Interval{0.0, 1.0}, saddle, Domain(saddle), ctx);
    REQUIRE(shortLine.has_value());
    CHECK(shortLine->Points.empty());
    const auto far = Meet(Curve{Line{Point3{0.0, 0.0, 10.0}, Vec3{1.0, 0.0, 0.0}}}, saddle);
    REQUIRE(far.has_value());
    CHECK(far->Points.empty());
}

TEST_CASE("what cannot be isolated fails closed", "[intersect][curve-surface][spline-surface]")
{
    const Surface saddle{Saddle()};
    SECTION("a line tangent to the saddle")
    {
        // Along u at v = 0, z = 0: it touches z = u^2 at the origin.
        const auto r = Meet(Curve{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}}, saddle);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("a line along the saddle")
    {
        // u = v gives z = 0: the diagonal lies in the surface.
        const auto r = Meet(Curve{Line{Point3{}, Vec3{1.0, 1.0, 0.0}}}, saddle);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("invalid surface data")
    {
        BSplineSurface bad = Saddle();
        bad.Poles.pop_back();
        const auto r = Intersect(Curve{Line{Point3{}, Vec3{0.0, 0.0, 1.0}}}, Interval::Everything(),
                                 Surface{bad}, Box2{Interval{lo, hi}, Interval{lo, hi}}, ctx);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("lines across a wavy bicubic and a rational cylinder are settled",
          "[intersect][curve-surface][spline-surface]")
{
    // The cases tests/occt/IntersectionDifferentialTests.cpp compares with
    // GeomAPI_IntCS, run here without OCCT: each settles, every point on both.
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<Point3> poles;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double z = 0.4 * std::sin(1.3 * i + 0.7 * j) - 0.1 * i * j;
            poles.push_back({i / 3.0 * 4.0, j / 3.0 * 3.0, z});
        }
    }
    const Surface patch{*Stereon::MakeBSplineSurface(3, 3, knots, knots, poles, {})};
    const Cylinder cylinder{Point3{1.0, -2.0, 0.5}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Surface rational{
        *Stereon::ToBSpline(cylinder, Box2{Interval{0.0, 2.0 * pi}, Interval{-20.0, 20.0}})};
    for (int k = 0; k < 24; ++k) {
        INFO("case " << k);
        const double u = 0.05 + 0.9 * std::fmod(0.37 * (k + 1), 1.0);
        const double v = 0.05 + 0.9 * std::fmod(0.61 * (k + 1), 1.0);
        const Vec3 d{0.15 * std::cos(k), 0.15 * std::sin(k), 1.0};
        const Curve down{Line{Eval(patch, Point2{u, v}) - 2.0 * d, d}};
        const auto onPatch = Meet(down, patch);
        INFO(IntersectTest::ErrorText(onPatch));
        REQUIRE(onPatch.has_value());
        CHECK(HasPoint(*onPatch, 2.0, Point2{u, v}, 1e-8));
        CheckOnBoth(*onPatch, down, patch);
        const Vec3 across{std::cos(0.7 * k), std::sin(0.7 * k), 0.3 * std::cos(1.1 * k)};
        const Curve through{Line{Point3{1.0, -2.0, 0.5 + 0.2 * k - 2.0}, across}};
        const auto onCylinder = Meet(through, rational);
        INFO(IntersectTest::ErrorText(onCylinder));
        REQUIRE(onCylinder.has_value());
        CHECK(onCylinder->Points.size() == 2);
        CheckOnBoth(*onCylinder, through, rational);
    }
}

TEST_CASE("a crossing at a closed surface's seam is one point at the domain's end",
          "[intersect][curve-surface][spline-surface]")
{
    // The NURBS cylinder's u domain ends at the seam, angle 0 and 2 pi at
    // once. A line along X through the axis crosses it there and at pi.
    // The crossing at the seam is found from the patches on both sides and
    // lies within rounding of the domain's end: it is placed at an end,
    // once, not refused (the first run of the bicubic-and-cylinder cases
    // found it refused).
    const Cylinder cylinder{Point3{1.0, -2.0, 0.5}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Surface rational{
        *Stereon::ToBSpline(cylinder, Box2{Interval{0.0, 2.0 * pi}, Interval{-20.0, 20.0}})};
    const Curve line{Line{Point3{-3.0, -2.0, 0.0}, Vec3{1.0, 0.0, 0.0}}};
    const auto r = Meet(line, rational);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Points.size() == 2);
    CheckOnBoth(*r, line, rational);
    CHECK_THAT(r->Points[0].Parameter, WithinAbs(2.0, 1e-9));
    CHECK_THAT(r->Points[1].Parameter, WithinAbs(6.0, 1e-9));
    const double seam = r->Points[1].OnSurface.X;
    CHECK((seam == 0.0 || seam == 2.0 * pi));
    // A box short of the whole turn has ends, not a seam: refused there.
    const auto part = Intersect(line, Interval::Everything(), rational,
                                Box2{Interval{0.0, pi}, Interval{-20.0, 20.0}}, ctx);
    REQUIRE_FALSE(part.has_value());
    CHECK(part.error().GetCode() == ErrorCode::IntersectionFailed);
}
