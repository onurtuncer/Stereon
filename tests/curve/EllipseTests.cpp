// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

import Stereon.Curve;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicEllipse;
using Stereon::Circle;
using Stereon::CirclePeriod;
using Stereon::Ellipse;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;

Ellipse TiltedEllipse(const double major = 3.0, const double minor = 1.25)
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, -1.0, 2.0}, 0.9);
    return {Point3{0.5, -1.0, 2.0}, major, minor, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Transform SomeSimilarity()
{
    return Stereon::MakeTranslation(Vec3{3.0, -1.0, 2.0}) *
           Stereon::MakeRotation(Vec3{0.0, 1.0, 1.0}, 1.1) * Stereon::MakeScaling(1.5);
}

bool Near(const Point3& a, const Point3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(b));
}

BasicEllipse<Dual, 3> Lifted(const Ellipse& e)
{
    const auto lift = [](const Vec3& v) { return Stereon::BasicVec3<Dual>{v.X, v.Y, v.Z}; };
    return {Stereon::BasicPoint3<Dual>{e.Center.X, e.Center.Y, e.Center.Z}, e.MajorRadius,
            e.MinorRadius, lift(e.XAxis), lift(e.YAxis)};
}

template <typename P>
Vec3 DerivPart(const P& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

/// The smallest distance from p to the ellipse over a dense sample of its
/// parameter: an upper bound on the true distance, which a correct
/// projection must meet.
double SampledDistance(const Ellipse& e, const Point3& p)
{
    constexpr int samples = 20000;
    double best = Distance(Eval(e, 0.0), p);
    for (int i = 1; i < samples; ++i) {
        best = std::min(best, Distance(Eval(e, i * CirclePeriod / samples), p));
    }
    return best;
}

Point3 RandomPoint()
{
    return {*rc::gen::inRange(-6000, 6000) / 1000.0, *rc::gen::inRange(-6000, 6000) / 1000.0,
            *rc::gen::inRange(-6000, 6000) / 1000.0};
}

} // namespace

TEST_CASE("ellipses satisfy the curve concepts and fit in the buffer", "[curve][ellipse]")
{
    STATIC_REQUIRE(Stereon::CurveGeometry<Ellipse>);
    STATIC_REQUIRE(Stereon::Curve2dGeometry<Ellipse2d>);
    STATIC_REQUIRE(Stereon::PeriodicCurve<Ellipse>);
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Ellipse, Stereon::GeometryBufferSize>);
}

TEST_CASE("an ellipse is evaluated at its vertices", "[curve][ellipse]")
{
    const Ellipse e{Point3{1.0, 1.0, 0.0}, 3.0, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Eval(e, 0.0) == Point3{4.0, 1.0, 0.0});
    REQUIRE(Near(Eval(e, std::numbers::pi / 2), Point3{1.0, 3.0, 0.0}));
    REQUIRE(Near(Eval(e, std::numbers::pi), Point3{-2.0, 1.0, 0.0}));
    REQUIRE(Near(Eval(e, 1.5 * std::numbers::pi), Point3{1.0, -1.0, 0.0}));
    REQUIRE(Domain(e) == Stereon::Interval{0.0, CirclePeriod});
    REQUIRE(Period(e) == CirclePeriod);
    // Quarter turns: the first derivative at 0 is the minor axis, the second
    // points back along the major axis, the fourth equals the radius vector.
    REQUIRE(Deriv(e, 0.0, 1) == Vec3{0.0, 2.0, 0.0});
    REQUIRE(Deriv(e, 0.0, 2) == Vec3{-3.0, 0.0, 0.0});
    REQUIRE(Deriv(e, 0.0, 4) == Vec3{3.0, 0.0, 0.0});
}

TEST_CASE("ellipse derivatives are the derivatives of its evaluation", "[curve][ellipse][property]")
{
    rc::prop("Deriv(e, t, n + 1) is d/dt of Deriv(e, t, n), and n = 1 of Eval", [] {
        const Ellipse e = TiltedEllipse();
        const double t = *rc::gen::inRange(-6000, 6000) / 1000.0;
        const int n = *rc::gen::inRange(1, 7);
        const Dual dt{t, 1.0};
        RC_ASSERT(Near(Deriv(e, t, 1), DerivPart(Eval(Lifted(e), dt))));
        RC_ASSERT(Near(Deriv(e, t, n + 1), DerivPart(Deriv(Lifted(e), dt, n))));
    });
}

TEST_CASE("Project inverts Eval and finds the closest point", "[curve][ellipse][property]")
{
    rc::prop("a point of the ellipse projects to its own parameter", [] {
        const Ellipse e = TiltedEllipse(*rc::gen::inRange(1000, 4000) / 1000.0, 1.0);
        const double t = *rc::gen::inRange(0, 6283) / 1000.0;
        const double back = Project(e, Eval(e, t));
        RC_ASSERT(back >= 0.0 && back < CirclePeriod);
        RC_ASSERT(Distance(Eval(e, back), Eval(e, t)) <= 1e-9);
    });

    rc::prop("no sampled point of the ellipse is closer than the projection", [] {
        const Ellipse e = TiltedEllipse(*rc::gen::inRange(1000, 5000) / 1000.0, 1.0);
        const Point3 p = RandomPoint();
        const double t = Project(e, p);
        const double distance = Distance(Eval(e, t), p);
        RC_ASSERT(distance <= SampledDistance(e, p) + 1e-12);
        // The foot is a stationary point: p - Eval is normal to the tangent.
        const Vec3 tangent = Deriv(e, t, 1);
        RC_ASSERT(std::abs(Dot(p - Eval(e, t), tangent)) <= 1e-9 * (1.0 + Length(tangent)));
    });
}

TEST_CASE("points on the axes of symmetry project as documented", "[curve][ellipse]")
{
    const Ellipse e{Point3{}, 4.0, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // The centre is closest to the minor vertices, and gets the one at pi/2.
    REQUIRE_THAT(Project(e, Point3{}), WithinAbs(std::numbers::pi / 2, Tolerance));
    // On the major axis outside the evolute's cusp at (a^2 - b^2) / a = 3:
    // the vertex. Inside it: two closest points off the axis, the one with
    // t in (0, pi).
    REQUIRE(Project(e, Point3{3.5, 0.0, 0.0}) == 0.0);
    REQUIRE_THAT(Project(e, Point3{-3.5, 0.0, 0.0}), WithinAbs(std::numbers::pi, Tolerance));
    const double inside = Project(e, Point3{1.5, 0.0, 0.0});
    REQUIRE(inside > 0.0);
    REQUIRE(inside < std::numbers::pi / 2);
    // Known closed form inside the cusp: x = a^2 x0 / (a^2 - b^2).
    REQUIRE_THAT(Eval(e, inside).X, WithinAbs(16.0 * 1.5 / 12.0, Tolerance));
    // On the minor axis: always the minor vertex on p's side.
    REQUIRE_THAT(Project(e, Point3{0.0, -0.5, 7.0}), WithinAbs(1.5 * std::numbers::pi, Tolerance));
    // Off the plane: the same as the projection onto the plane.
    REQUIRE(Project(e, Point3{1.0, 1.0, 9.0}) == Project(e, Point3{1.0, 1.0, 0.0}));
}

TEST_CASE("an ellipse with equal radii is the circle", "[curve][ellipse]")
{
    const Circle c{Point3{1.0, 2.0, 3.0}, 2.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const Ellipse e{c.Center, 2.0, 2.0, c.XAxis, c.YAxis};
    for (int i = 0; i <= 16; ++i) {
        const double t = i * CirclePeriod / 16.0;
        REQUIRE(Near(Eval(e, t), Eval(c, t)));
        REQUIRE(Near(Deriv(e, t, 3), Deriv(c, t, 3)));
        const Point3 p{1.0 + std::cos(t), 4.0, 3.0 + std::sin(t)};
        const double gap = std::abs(Project(e, p) - Project(c, p));
        REQUIRE(std::min(gap, CirclePeriod - gap) <= 1e-12);
    }
    REQUIRE(Project(e, c.Center) == Project(c, c.Center));
}

TEST_CASE("ellipses reverse, move and bound consistently", "[curve][ellipse][property]")
{
    rc::prop("Reversed runs backwards, Transformed commutes, the box contains the curve", [] {
        const Ellipse e = TiltedEllipse();
        const Transform x = SomeSimilarity();
        const double t = *rc::gen::inRange(0, 6283) / 1000.0;
        RC_ASSERT(Near(Eval(Reversed(e), CirclePeriod - t), Eval(e, t)));
        RC_ASSERT(Near(Eval(Transformed(e, x), t), Transformed(Eval(e, t), x)));
        RC_ASSERT(Contains(BoundingBox(e), Eval(e, t)));
    });

    // The box is tight on an axis-aligned ellipse.
    const Ellipse flat{Point3{}, 3.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Stereon::Box3 box = BoundingBox(flat);
    REQUIRE_THAT(box.X.GetUpper(), WithinAbs(3.0, Tolerance));
    REQUIRE_THAT(box.Y.GetLower(), WithinAbs(-1.0, Tolerance));
    REQUIRE(Contains(box.Z, 0.0));
    REQUIRE(box.Z.GetUpper() < Tolerance);
    REQUIRE(Reversed(Reversed(flat)) == flat);
}

TEST_CASE("planar ellipses evaluate and project in the plane", "[curve][ellipse]")
{
    const Ellipse2d e{Point2{1.0, 0.0}, 2.0, 1.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}};
    REQUIRE(Eval(e, 0.0) == Point2{1.0, 2.0});
    REQUIRE(Near(Point3{Eval(e, std::numbers::pi / 2).X, Eval(e, std::numbers::pi / 2).Y, 0.0},
                 Point3{0.0, 0.0, 0.0}));
    REQUIRE_THAT(Project(e, Point2{1.0, 5.0}), WithinAbs(0.0, Tolerance));
    REQUIRE(Contains(BoundingBox(e), Point2{1.0, -2.0}));
    const Stereon::Curve2d wrapped(e);
    REQUIRE(Eval(wrapped, 1.0) == Eval(e, 1.0));
}

TEST_CASE("MakeEllipse validates data from outside the kernel", "[curve][ellipse][factory]")
{
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    const auto made = Stereon::MakeEllipse(Point3{}, 2.0, 1.0, x, y);
    REQUIRE(made.has_value());
    REQUIRE(*made == Ellipse{Point3{}, 2.0, 1.0, x, y});
    REQUIRE(Stereon::MakeEllipse(Point2{}, 1.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}).has_value());

    const auto codeOf = [](const auto& result) { return result.error().GetCode(); };
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, 1.0, 2.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, 2.0, 0.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, 2.0, -1.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, std::nan(""), 1.0, x, y)) ==
            ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{std::nan(""), 0.0, 0.0}, 2.0, 1.0, x, y)) ==
            ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, 2.0, 1.0, x, Vec3{0.1, 1.0, 0.0})) ==
            ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeEllipse(Point3{}, 2.0, 1.0, x, 2.0 * y)) ==
            ErrorCode::InvalidInput);
    // Within FrameAllowance.
    REQUIRE(Stereon::MakeEllipse(Point3{}, 2.0, 1.0, x, Vec3{1e-14, 1.0, 0.0}).has_value());
}

TEST_CASE("the frame allowance is a closed bound", "[curve][factory]")
{
    // ADR-0028: a frame whose error equals FrameAllowance is accepted, one
    // just above it is not. With y = (a, 1, 0), x . y is exactly a, and
    // y . y - 1 = a^2 = 2^-80 rounds away. The test above, at 1e-14, runs the
    // same check without pinning its boundary: the mutation pilot of
    // ADR-0032 found that "<=" could become "<" unnoticed.
    const Vec3 x{1.0, 0.0, 0.0};
    const double a = Stereon::FrameAllowance;
    REQUIRE(Stereon::MakeCircle(Point3{}, 1.0, x, Vec3{a, 1.0, 0.0}).has_value());
    const double above = std::nextafter(a, 1.0);
    REQUIRE_FALSE(Stereon::MakeCircle(Point3{}, 1.0, x, Vec3{above, 1.0, 0.0}).has_value());
}

TEST_CASE("the other analytic curves validate their data too", "[curve][factory]")
{
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    REQUIRE(Stereon::MakeCircle(Point3{}, 1.0, x, y).has_value());
    REQUIRE_FALSE(Stereon::MakeCircle(Point3{}, 0.0, x, y).has_value());
    REQUIRE_FALSE(Stereon::MakeCircle(Point3{}, 1.0, x, x).has_value());
    REQUIRE(Stereon::MakeLine(Point3{}, x).has_value());
    REQUIRE_FALSE(Stereon::MakeLine(Point3{}, Vec3{}).has_value());
    REQUIRE_FALSE(Stereon::MakeLine(Point3{}, Vec3{INFINITY, 0.0, 0.0}).has_value());
    // The message names the factory and the value.
    const auto bad = Stereon::MakeCircle(Point3{}, -2.0, x, y);
    REQUIRE(bad.error().GetMessage() == "MakeCircle: Radius must be finite and positive");
}
