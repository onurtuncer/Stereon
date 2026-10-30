// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The hyperbola and the parabola, and the polynomial root finder their
// projections rest on: implicit equations, derivatives against duals,
// projections against dense sampling.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

import Stereon.Curve;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicHyperbola;
using Stereon::BasicParabola;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Hyperbola2d;
using Stereon::Interval;
using Stereon::Parabola;
using Stereon::Parabola2d;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;

bool Near(const Point3& a, const Point3& b, const double tolerance = Tolerance)
{
    return Length(a - b) <= tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b, const double tolerance = Tolerance)
{
    return Length(a - b) <= tolerance * (1.0 + Length(b));
}

struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{2.0, -1.0, 0.5}, 0.8);
    return {Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Hyperbola TiltedHyperbola(const double a = 2.0, const double b = 1.25)
{
    const Frame f = TiltedFrame();
    return {Point3{1.0, -0.5, 2.0}, a, b, f.X, f.Y};
}

Parabola TiltedParabola(const double focal = 0.75)
{
    const Frame f = TiltedFrame();
    return {Point3{-1.0, 0.5, 1.0}, focal, f.X, f.Y};
}

Transform SomeSimilarity()
{
    return Stereon::MakeTranslation(Vec3{3.0, -1.0, 2.0}) *
           Stereon::MakeRotation(Vec3{0.0, 1.0, 1.0}, 1.1) * Stereon::MakeScaling(1.5);
}

Stereon::BasicVec3<Dual> Lift(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

template <typename P>
Vec3 DerivPart(const P& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

/// The smallest distance from p over samples of the curve on [-limit, limit].
template <typename C>
double SampledDistance(const C& c, const Point3& p, const double limit)
{
    constexpr int samples = 40000;
    double best = Distance(Eval(c, -limit), p);
    for (int i = 1; i <= samples; ++i) {
        best = std::min(best, Distance(Eval(c, -limit + 2.0 * limit * i / samples), p));
    }
    return best;
}

Point3 RandomPoint()
{
    return {*rc::gen::inRange(-6000, 6000) / 1000.0, *rc::gen::inRange(-6000, 6000) / 1000.0,
            *rc::gen::inRange(-6000, 6000) / 1000.0};
}

} // namespace

TEST_CASE("the open conics satisfy the curve concepts", "[curve][conic]")
{
    STATIC_REQUIRE(Stereon::CurveGeometry<Hyperbola>);
    STATIC_REQUIRE(Stereon::Curve2dGeometry<Hyperbola2d>);
    STATIC_REQUIRE(Stereon::CurveGeometry<Parabola>);
    STATIC_REQUIRE(Stereon::Curve2dGeometry<Parabola2d>);
    STATIC_REQUIRE_FALSE(Stereon::PeriodicCurve<Hyperbola>);
    REQUIRE(Domain(TiltedHyperbola()) == Interval::Everything());
    REQUIRE(BoundingBox(TiltedParabola()) == Stereon::Box3::Everything());
}

TEST_CASE("PolynomialRoots isolates every simple root", "[curve][conic][roots]")
{
    // (t - 1)(t - 2)(t - 3) = t^3 - 6 t^2 + 11 t - 6.
    const std::vector<double> cubic{-6.0, 11.0, -6.0, 1.0};
    const std::vector<double> roots = Stereon::Detail::PolynomialRoots(cubic, Interval{0.0, 4.0});
    REQUIRE(roots.size() == 3);
    REQUIRE_THAT(roots[0], WithinAbs(1.0, 1e-14));
    REQUIRE_THAT(roots[1], WithinAbs(2.0, 1e-14));
    REQUIRE_THAT(roots[2], WithinAbs(3.0, 1e-14));
    REQUIRE(Stereon::Detail::PolynomialRoots(cubic, Interval{3.5, 9.0}).empty());
    REQUIRE(Stereon::Detail::PolynomialRoots(std::vector{1.0, 0.0, 1.0}, Interval{-5.0, 5.0})
                .empty()); // t^2 + 1
    // A double root touches zero; whatever is reported lies next to it.
    for (const double r :
         Stereon::Detail::PolynomialRoots(std::vector{1.0, -2.0, 1.0}, Interval{-2.0, 3.0})) {
        REQUIRE_THAT(r, WithinAbs(1.0, 1e-6));
    }
}

TEST_CASE("a hyperbola lies on its implicit equation", "[curve][conic][property]")
{
    rc::prop("(x/a)^2 - (y/b)^2 = 1 with x > 0, and derivatives are the duals'", [] {
        const Hyperbola h = TiltedHyperbola();
        const double t = *rc::gen::inRange(-3000, 3000) / 1000.0;
        const Vec3 offset = Eval(h, t) - h.Center;
        const double x = Dot(offset, h.XAxis) / h.MajorRadius;
        const double y = Dot(offset, h.YAxis) / h.MinorRadius;
        RC_ASSERT(std::abs(x * x - y * y - 1.0) <= 1e-12 * (1.0 + x * x));
        RC_ASSERT(x > 0.0);
        const BasicHyperbola<Dual, 3> lifted{
            Stereon::BasicPoint3<Dual>{h.Center.X, h.Center.Y, h.Center.Z}, h.MajorRadius,
            h.MinorRadius, Lift(h.XAxis), Lift(h.YAxis)};
        const int n = *rc::gen::inRange(1, 6);
        RC_ASSERT(Near(Deriv(h, t, 1), DerivPart(Eval(lifted, Dual{t, 1.0}))));
        RC_ASSERT(Near(Deriv(h, t, n + 1), DerivPart(Deriv(lifted, Dual{t, 1.0}, n))));
    });
}

TEST_CASE("hyperbola projection finds the closest point", "[curve][conic][property]")
{
    rc::prop("no sample is closer; points of the curve project to themselves", [] {
        const Hyperbola h = TiltedHyperbola(*rc::gen::inRange(500, 3000) / 1000.0,
                                            *rc::gen::inRange(500, 3000) / 1000.0);
        const Point3 p = RandomPoint();
        const double t = Project(h, p);
        RC_ASSERT(Distance(Eval(h, t), p) <= SampledDistance(h, p, 6.0) + 1e-12);
        const double on = *rc::gen::inRange(-3000, 3000) / 1000.0;
        RC_ASSERT(std::abs(Project(h, Eval(h, on)) - on) <= 1e-9);
    });
}

TEST_CASE("a parabola lies on its implicit equation", "[curve][conic][property]")
{
    rc::prop("y^2 = 4 f x, and derivatives are the duals'", [] {
        const Parabola c = TiltedParabola();
        const double t = *rc::gen::inRange(-10000, 10000) / 1000.0;
        const Vec3 offset = Eval(c, t) - c.Vertex;
        const double x = Dot(offset, c.XAxis);
        const double y = Dot(offset, c.YAxis);
        RC_ASSERT(std::abs(y * y - 4.0 * c.FocalLength * x) <= 1e-12 * (1.0 + y * y));
        const BasicParabola<Dual, 3> lifted{
            Stereon::BasicPoint3<Dual>{c.Vertex.X, c.Vertex.Y, c.Vertex.Z}, c.FocalLength,
            Lift(c.XAxis), Lift(c.YAxis)};
        RC_ASSERT(Near(Deriv(c, t, 1), DerivPart(Eval(lifted, Dual{t, 1.0}))));
        RC_ASSERT(Near(Deriv(c, t, 2), DerivPart(Deriv(lifted, Dual{t, 1.0}, 1))));
        RC_ASSERT(Deriv(c, t, 3) == Vec3{});
    });
}

TEST_CASE("parabola projection finds the closest point", "[curve][conic][property]")
{
    rc::prop("no sample is closer; points of the curve project to themselves", [] {
        const Parabola c = TiltedParabola(*rc::gen::inRange(100, 3000) / 1000.0);
        const Point3 p = RandomPoint();
        const double t = Project(c, p);
        RC_ASSERT(Distance(Eval(c, t), p) <= SampledDistance(c, p, 20.0) + 1e-12);
        const double on = *rc::gen::inRange(-5000, 5000) / 1000.0;
        RC_ASSERT(std::abs(Project(c, Eval(c, on)) - on) <= 1e-9);
    });
}

TEST_CASE("points on the axes project as documented", "[curve][conic]")
{
    const Hyperbola h{Point3{}, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // On the symmetry axis the vertex is closest up to the centre of
    // curvature, x = a + b^2/a = 2.5 here; beyond it two mirror points tie.
    REQUIRE_THAT(Project(h, Point3{2.2, 0.0, 0.0}), WithinAbs(0.0, 1e-12));
    const Point3 beyond{6.0, 0.0, 0.0};
    const double far = Project(h, beyond);
    REQUIRE(std::abs(far) > 0.5);
    REQUIRE_THAT(Distance(Eval(h, far), beyond), WithinAbs(Distance(Eval(h, -far), beyond), 1e-12));
    // The parabola's centre of curvature at the vertex is at x = 2f = 2.
    const Parabola c{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Project(c, Point3{}) == 0.0);
    REQUIRE_THAT(Project(c, Point3{1.5, 0.0, 0.0}), WithinAbs(0.0, 1e-12));
    REQUIRE(std::abs(Project(c, Point3{5.0, 0.0, 0.0})) > 1.0);
}

TEST_CASE("the open conics reverse, move and work in the plane", "[curve][conic][property]")
{
    rc::prop("Reversed runs backwards and Transformed commutes with Eval", [] {
        const Hyperbola h = TiltedHyperbola();
        const Parabola c = TiltedParabola();
        const Transform x = SomeSimilarity();
        const double t = *rc::gen::inRange(-2000, 2000) / 1000.0;
        RC_ASSERT(Near(Eval(Reversed(h), -t), Eval(h, t)));
        RC_ASSERT(Near(Eval(Reversed(c), -t), Eval(c, t)));
        RC_ASSERT(Near(Eval(Transformed(h, x), t), Transformed(Eval(h, t), x), 1e-11));
        RC_ASSERT(Near(Eval(Transformed(c, x), 1.5 * t), Transformed(Eval(c, t), x), 1e-11));
    });
    const Parabola2d planar{Point2{}, 0.5, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}};
    REQUIRE(Eval(planar, 2.0) == Point2{-2.0, 2.0});
    REQUIRE_THAT(Project(planar, Point2{-2.0, 2.0}), WithinAbs(2.0, 1e-12));
    const Hyperbola2d planarH{Point2{}, 1.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}};
    REQUIRE_THAT(Project(planarH, Point2{1.0, 0.0}), WithinAbs(0.0, 1e-12));
    const Stereon::Curve tilted(TiltedHyperbola());
    REQUIRE(tilted.Target<Hyperbola>() != nullptr);
}

TEST_CASE("the open conics validate data from outside the kernel", "[curve][conic][factory]")
{
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    REQUIRE(Stereon::MakeHyperbola(Point3{}, 1.0, 3.0, x, y).has_value()); // any order of radii
    REQUIRE(Stereon::MakeHyperbola(Point2{}, 1.0, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}).has_value());
    REQUIRE(Stereon::MakeHyperbola(Point3{}, 0.0, 1.0, x, y).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(Stereon::MakeHyperbola(Point3{}, 1.0, 1.0, x, x).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(Stereon::MakeParabola(Point3{}, 1.0, x, y).has_value());
    REQUIRE(Stereon::MakeParabola(Point3{}, -1.0, x, y).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(Stereon::MakeParabola(Point2{NAN, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0})
                .error()
                .GetCode() == ErrorCode::InvalidInput);
}
