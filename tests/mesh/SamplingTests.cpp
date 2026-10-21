// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Display sampling (ADR-0035 verification): the chord deviation stays within
// the display tolerance on every built-in type, B-spline polylines pass
// through every distinct knot, unbounded geometry is clipped to its range,
// and bad tolerances and ranges are refused.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

import Stereon.Mesh;

#include "FaultGeometry.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::DisplayTolerance;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::IsoDirection;
using Stereon::Line;
using Stereon::Parabola;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Polyline;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr Vec3 X{1.0, 0.0, 0.0};
constexpr Vec3 Y{0.0, 1.0, 0.0};
constexpr Point3 O{1.0, -2.0, 0.5};

/// The chord deviation is checked at seven interior points of a segment;
/// between them it may exceed the tolerance slightly on a curve whose
/// curvature varies, which the documented contract allows. The tests check
/// at fifteen points against this much more.
constexpr double Slack = 1.05;

DisplayTolerance Tolerance(const double chord, const double angle,
                           const std::size_t maxPoints = 65536)
{
    DisplayTolerance tolerance;
    tolerance.Chord = chord;
    tolerance.Angle = angle;
    tolerance.MaxPoints = maxPoints;
    return tolerance;
}

double DistanceToSegment(const Point3& p, const Point3& a, const Point3& b)
{
    const Vec3 v = b - a;
    const Vec3 w = p - a;
    const double length2 = Stereon::Dot(v, v);
    if (length2 == 0.0) {
        return Stereon::Length(w);
    }
    const double s = std::clamp(Stereon::Dot(w, v) / length2, 0.0, 1.0);
    return Stereon::Length(w - s * v);
}

/// The largest distance from the curve to the polyline's segments, at
/// fifteen points inside each.
double MaxDeviation(const Curve& curve, const Polyline& line)
{
    double worst = 0.0;
    for (std::size_t i = 0; i + 1 < line.Points.size(); ++i) {
        const double t0 = line.Parameters[i];
        const double t1 = line.Parameters[i + 1];
        for (int k = 1; k < 16; ++k) {
            const Point3 p = Stereon::Eval(curve, t0 + (t1 - t0) * k / 16);
            worst = std::max(worst, DistanceToSegment(p, line.Points[i], line.Points[i + 1]));
        }
    }
    return worst;
}

void RequireConsistent(const Curve& curve, const Polyline& line)
{
    REQUIRE(line.Points.size() >= 2);
    REQUIRE(line.Points.size() == line.Parameters.size());
    REQUIRE(std::ranges::is_sorted(line.Parameters));
    REQUIRE(std::ranges::adjacent_find(line.Parameters) == line.Parameters.end());
    for (std::size_t i = 0; i < line.Points.size(); ++i) {
        REQUIRE(line.Points[i] == Stereon::Eval(curve, line.Parameters[i]));
    }
}

bool Contains(const std::vector<double>& values, const double x)
{
    return std::ranges::find(values, x) != values.end();
}

/// A cubic with an interior knot of multiplicity 3: a kink at t = 0.5.
BSplineCurve Kinked()
{
    return BSplineCurve{3,
                        {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0},
                        {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
                         Point3{3.0, 1.0, 2.0}, Point3{4.0, 3.0, 0.0}, Point3{5.0, 0.0, 0.0},
                         Point3{6.0, 2.0, 1.0}, Point3{7.0, 0.0, 0.0}},
                        {1.0, 0.5, 2.0, 0.75, 1.0, 1.0, 1.5, 1.0}};
}

/// A discontinuity at the start of a finite domain. Each leftmost segment
/// still turns through nearly pi/2, however often it is subdivided.
struct StepCurve {
    double Height = 1.0;
};
Point3 Eval(const StepCurve& c, const double t)
{
    return {t, t == 0.0 ? 0.0 : c.Height, 0.0};
}
Vec3 Deriv(const StepCurve&, const double, const int n)
{
    return n == 1 ? X : Vec3{};
}
Interval Domain(const StepCurve&)
{
    return {0.0, 1.0};
}
Stereon::Box3 BoundingBox(const StepCurve& c)
{
    return {{0, 1}, {0, c.Height}, {0, 0}};
}
double Project(const StepCurve&, const Point3& p)
{
    return std::clamp(p.X, 0.0, 1.0);
}
StepCurve Reversed(const StepCurve& c)
{
    return c;
}
StepCurve Transformed(const StepCurve& c, const Stereon::Transform&)
{
    return c;
}

} // namespace

TEST_CASE("every built-in curve is sampled within the chord tolerance", "[mesh][sampling]")
{
    const DisplayTolerance tolerance = Tolerance(1e-3, 0.2);
    const Interval range{0.25, 2.0};
    const std::vector<std::pair<std::string, Curve>> curves{
        {"line", Curve(Line{O, X})},
        {"circle", Curve(Circle{O, 2.5, X, Y})},
        {"ellipse", Curve(Ellipse{O, 3.0, 1.25, X, Y})},
        {"hyperbola", Curve(Hyperbola{O, 2.0, 0.75, X, Y})},
        {"parabola", Curve(Parabola{O, 0.5, X, Y})},
        {"rational B-spline", Curve(Kinked())},
    };
    for (const auto& [name, curve] : curves) {
        INFO(name);
        const auto line = Stereon::SampleCurve(curve, tolerance, range);
        REQUIRE(line);
        RequireConsistent(curve, *line);
        REQUIRE(MaxDeviation(curve, *line) <= Slack * tolerance.Chord);
    }
}

TEST_CASE("a circle's polyline deviates by at most the chord tolerance", "[mesh][sampling]")
{
    // On an arc the deviation peaks at the segment's middle, a point the
    // sampler checks: the bound holds without slack.
    const Curve circle(Circle{O, 10.0, X, Y});
    const DisplayTolerance tolerance = Tolerance(1e-2, 1.0);
    const auto line = Stereon::SampleCurve(circle, tolerance);
    REQUIRE(line);
    REQUIRE(line->Parameters.front() == Stereon::Domain(circle).GetLower());
    REQUIRE(line->Parameters.back() == Stereon::Domain(circle).GetUpper());
    REQUIRE(MaxDeviation(circle, *line) <= tolerance.Chord * (1.0 + 1e-9));
}

TEST_CASE("the angle tolerance bounds the turn between half-segments", "[mesh][sampling]")
{
    // A chord tolerance this loose would draw a circle as a few segments;
    // the angle keeps it round: at most 0.1 rad of turn per half segment.
    const Curve circle(Circle{O, 1.0, X, Y});
    const DisplayTolerance tolerance = Tolerance(1e3, 0.1);
    const auto line = Stereon::SampleCurve(circle, tolerance);
    REQUIRE(line);
    REQUIRE(line->Points.size() >= static_cast<std::size_t>(2 * std::numbers::pi / 0.2));
}

TEST_CASE("a B-spline's polyline passes through every distinct knot", "[mesh][sampling]")
{
    const Curve curve(Kinked());
    const auto line = Stereon::SampleCurve(curve, DisplayTolerance{});
    REQUIRE(line);
    for (const double knot : {0.0, 0.25, 0.5, 1.0}) {
        INFO(knot);
        REQUIRE(Contains(line->Parameters, knot));
    }
    REQUIRE(Stereon::Breakpoints(curve, Interval{0.1, 0.9}) ==
            std::vector<double>{0.1, 0.25, 0.5, 0.9});
}

TEST_CASE("a trimmed B-spline keeps its basis's knots, reversed or not", "[mesh][sampling]")
{
    const Curve basis(Kinked());
    const Interval trim{0.1, 0.8};
    const auto forward = Stereon::MakeTrimmedCurve(basis, trim);
    REQUIRE(forward);
    const Curve reversed(forward->Reversed());
    const auto line = Stereon::SampleCurve(reversed, DisplayTolerance{});
    REQUIRE(line);
    RequireConsistent(reversed, *line);
    for (const double knot : {0.25, 0.5}) {
        INFO(knot);
        REQUIRE(Contains(line->Parameters, forward->Reflect(knot)));
    }
}

TEST_CASE("unbounded geometry needs a range, and a range is clipped to the domain",
          "[mesh][sampling]")
{
    const Curve line(Line{O, X});
    const auto unbounded = Stereon::SampleCurve(line, DisplayTolerance{});
    REQUIRE_FALSE(unbounded);
    REQUIRE(unbounded.error().GetCode() == ErrorCode::InvalidInput);

    const auto drawn = Stereon::SampleCurve(line, DisplayTolerance{}, Interval{-1.0, 4.0});
    REQUIRE(drawn);
    REQUIRE(drawn->Parameters.front() == -1.0);
    REQUIRE(drawn->Parameters.back() == 4.0);
    REQUIRE(drawn->Points.size() == 5); // one span, its four first segments flat

    const Curve bspline(Kinked());
    const auto clipped = Stereon::SampleCurve(bspline, DisplayTolerance{}, Interval{-1.0, 0.4});
    REQUIRE(clipped);
    REQUIRE(clipped->Parameters.front() == 0.0);
    REQUIRE(clipped->Parameters.back() == 0.4);

    const auto outside = Stereon::SampleCurve(bspline, DisplayTolerance{}, Interval{2.0, 3.0});
    REQUIRE_FALSE(outside);
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(Stereon::SampleCurve(line, DisplayTolerance{}, Interval{0.0, inf}));
}

TEST_CASE("bad display tolerances are refused, and a budget overrun is reported",
          "[mesh][sampling]")
{
    const Curve circle(Circle{O, 1.0, X, Y});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const DisplayTolerance bad :
         {Tolerance(0.0, 0.1), Tolerance(-1.0, 0.1), Tolerance(nan, 0.1), Tolerance(0.01, 0.0),
          Tolerance(0.01, 0.1, 1)}) {
        const auto refused = Stereon::SampleCurve(circle, bad);
        REQUIRE_FALSE(refused);
        REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
    }
    const auto overrun = Stereon::SampleCurve(circle, Tolerance(1e-9, 0.1, 50));
    REQUIRE_FALSE(overrun);
    REQUIRE(overrun.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE_FALSE(Stereon::SampleCurve(Curve{}, DisplayTolerance{}));
}

TEST_CASE("plane curves are sampled like space curves", "[mesh][sampling]")
{
    const Curve2d circle(Circle2d{Point2{0.5, 0.5}, 0.25, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}});
    const auto line = Stereon::SampleCurve(circle, Tolerance(1e-4, 0.1));
    REQUIRE(line);
    REQUIRE(line->Points.size() == line->Parameters.size());
    for (std::size_t i = 0; i < line->Points.size(); ++i) {
        REQUIRE(line->Points[i] == Stereon::Eval(circle, line->Parameters[i]));
    }
}

TEST_CASE("a surface is drawn as iso-lines in both directions", "[mesh][sampling]")
{
    const Surface sphere(Sphere{O, 2.5, X, Y});
    const auto lines = Stereon::SampleSurface(sphere, DisplayTolerance{}, std::nullopt, 5);
    REQUIRE(lines);
    const auto alongU =
        std::ranges::count(*lines, IsoDirection::AlongU, &Stereon::IsoLine::Direction);
    const auto alongV =
        std::ranges::count(*lines, IsoDirection::AlongV, &Stereon::IsoLine::Direction);
    REQUIRE(alongU == 5);
    REQUIRE(alongV == 5);
    for (const Stereon::IsoLine& iso : *lines) {
        for (std::size_t i = 0; i < iso.Line.Points.size(); ++i) {
            const double s = iso.Line.Parameters[i];
            const Point2 uv =
                iso.Direction == IsoDirection::AlongU ? Point2{s, iso.Fixed} : Point2{iso.Fixed, s};
            REQUIRE(iso.Line.Points[i] == Stereon::Eval(sphere, uv));
        }
    }
}

TEST_CASE("a B-spline surface gets an iso-line at every interior knot", "[mesh][sampling]")
{
    const Surface patch(
        BSplineSurface{1,
                       1,
                       {0.0, 0.0, 0.3, 1.0, 1.0},
                       {0.0, 0.0, 1.0, 1.0},
                       {Point3{0.0, 0.0, 0.0}, Point3{0.0, 1.0, 0.0}, Point3{1.0, 0.0, 1.0},
                        Point3{1.0, 1.0, 1.0}, Point3{2.0, 0.0, 0.0}, Point3{2.0, 1.0, 0.0}},
                       {}});
    const auto lines = Stereon::SampleSurface(patch, DisplayTolerance{}, std::nullopt, 2);
    REQUIRE(lines);
    std::vector<double> fixedU;
    for (const Stereon::IsoLine& iso : *lines) {
        if (iso.Direction == IsoDirection::AlongV) {
            fixedU.push_back(iso.Fixed);
        }
        // Each AlongU line passes through the kink at u = 0.3.
        if (iso.Direction == IsoDirection::AlongU) {
            REQUIRE(Contains(iso.Line.Parameters, 0.3));
        }
    }
    REQUIRE(fixedU == std::vector<double>{0.0, 0.3, 1.0});
}

TEST_CASE("an unbounded surface needs a range", "[mesh][sampling]")
{
    const Surface plane(Plane{O, X, Y});
    REQUIRE_FALSE(Stereon::SampleSurface(plane, DisplayTolerance{}));
    const auto lines = Stereon::SampleSurface(
        plane, DisplayTolerance{}, Stereon::Box2{Interval{0.0, 2.0}, Interval{-1.0, 1.0}}, 3);
    REQUIRE(lines);
    REQUIRE(lines->size() == 6);
    REQUIRE_FALSE(Stereon::SampleSurface(
        plane, DisplayTolerance{}, Stereon::Box2{Interval{0.0, 2.0}, Interval{-1.0, 1.0}}, 1));
}

TEST_CASE("sampling refuses unresolved segments at the subdivision limit",
          "[mesh][sampling][property]")
{
    for (const double height : {0.5, 1.0, 2.0, 4.0}) {
        const Curve step(StepCurve{height});
        const auto sampled = Stereon::SampleCurve(step, Tolerance(1e-6, 0.1));
        REQUIRE_FALSE(sampled);
        REQUIRE(sampled.error().GetCode() == ErrorCode::NotSupported);
        REQUIRE(std::string(sampled.error().GetMessage()).find("subdivision") != std::string::npos);
    }
}

TEST_CASE("sampling refuses non-finite evaluations instead of publishing a polyline",
          "[mesh][sampling]")
{
    for (const double bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        FaultGeometry::Curve faulty;
        faulty.Value = bad;
        const auto sampled = Stereon::SampleCurve(Curve(faulty), Tolerance(1e-3, 0.2));
        REQUIRE_FALSE(sampled);
        REQUIRE(sampled.error().GetCode() == ErrorCode::NotSupported);
    }
}
