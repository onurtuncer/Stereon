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
#include <type_traits>

import Stereon.Surface;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicCylinder;
using Stereon::CirclePeriod;
using Stereon::Cylinder;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::SurfaceGeometry;
using Stereon::Transform;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;

Cylinder TiltedCylinder()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{2.0, -1.0, 1.0}, 0.4);
    return {Point3{0.5, 1.0, -2.0}, 1.75, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Plane TiltedPlane()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 1.0, 0.0}, 0.9);
    return {Point3{1.0, 2.0, 3.0}, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

constexpr double Scale = 1.5;

Transform SomeSimilarity()
{
    return Stereon::MakeTranslation(Vec3{-1.0, 4.0, 0.5}) *
           Stereon::MakeRotation(Vec3{0.0, 0.0, 1.0}, 2.2) * Stereon::MakeScaling(Scale);
}

bool Near(const Point3& a, const Point3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(b));
}

BasicCylinder<Dual> Lifted(const Cylinder& c)
{
    const auto lift = [](const Vec3& v) { return Stereon::BasicVec3<Dual>{v.X, v.Y, v.Z}; };
    return {Stereon::BasicPoint3<Dual>{c.Origin.X, c.Origin.Y, c.Origin.Z}, c.Radius, lift(c.XAxis),
            lift(c.YAxis), c.RightHanded};
}

template <typename P>
Vec3 DerivPart(const P& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

Point2 RandomUv()
{
    return {*rc::gen::inRange(0, 6283) / 1000.0, *rc::gen::inRange(-5000, 5000) / 1000.0};
}

} // namespace

TEST_CASE("surface types satisfy the concept", "[surface][concepts]")
{
    STATIC_REQUIRE(SurfaceGeometry<Plane>);
    STATIC_REQUIRE(SurfaceGeometry<Cylinder>);
    STATIC_REQUIRE(SurfaceGeometry<Surface>);
    STATIC_REQUIRE_FALSE(SurfaceGeometry<Stereon::Circle>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<Surface, Stereon::Line>);
}

TEST_CASE("a plane is evaluated, projected and reversed at compile time", "[surface][plane]")
{
    constexpr Plane plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    STATIC_REQUIRE(Eval(plane, Point2{2.0, 3.0}) == Point3{2.0, 3.0, 1.0});
    STATIC_REQUIRE(Deriv(plane, Point2{}, 1, 0) == Vec3{1.0, 0.0, 0.0});
    STATIC_REQUIRE(Deriv(plane, Point2{}, 0, 1) == Vec3{0.0, 1.0, 0.0});
    STATIC_REQUIRE(Deriv(plane, Point2{}, 1, 1) == Vec3{});
    STATIC_REQUIRE(Project(plane, Point3{4.0, -2.0, 9.0}) == Point2{4.0, -2.0});
    STATIC_REQUIRE(Domain(plane) == Stereon::Box2::Everything());
    STATIC_REQUIRE(BoundingBox(plane) == Stereon::Box3::Everything());
    STATIC_REQUIRE(Eval(Reversed(plane), Point2{2.0, -3.0}) == Eval(plane, Point2{2.0, 3.0}));
}

TEST_CASE("plane operations agree with each other", "[surface][plane][property]")
{
    rc::prop("Project inverts Eval; Transformed scales the parameter", [] {
        const Plane plane = TiltedPlane();
        const Transform x = SomeSimilarity();
        const Point2 uv = RandomUv();
        const Point2 back = Project(plane, Eval(plane, uv));
        RC_ASSERT(std::abs(back.X - uv.X) <= 1e-9 && std::abs(back.Y - uv.Y) <= 1e-9);
        const Point2 scaled{Scale * uv.X, Scale * uv.Y};
        RC_ASSERT(Near(Eval(Transformed(plane, x), scaled), Transformed(Eval(plane, uv), x)));
    });
}

TEST_CASE("cylinder derivatives are the derivatives of its evaluation",
          "[surface][cylinder][property]")
{
    rc::prop("partials in u and v, of every order", [] {
        const Cylinder c = TiltedCylinder();
        const Point2 uv = RandomUv();
        const int n = *rc::gen::inRange(1, 6);
        const Stereon::BasicPoint2<Dual> du{Dual{uv.X, 1.0}, Dual{uv.Y}};
        const Stereon::BasicPoint2<Dual> dv{Dual{uv.X}, Dual{uv.Y, 1.0}};
        RC_ASSERT(Near(Deriv(c, uv, 1, 0), DerivPart(Eval(Lifted(c), du))));
        RC_ASSERT(Near(Deriv(c, uv, 0, 1), DerivPart(Eval(Lifted(c), dv))));
        RC_ASSERT(Near(Deriv(c, uv, n + 1, 0), DerivPart(Deriv(Lifted(c), du, n, 0))));
        RC_ASSERT(Near(Deriv(c, uv, n, 1), DerivPart(Deriv(Lifted(c), dv, n, 0))));
        RC_ASSERT(Deriv(c, uv, 0, n + 1) == Vec3{});
    });
}

TEST_CASE("cylinder operations agree with each other", "[surface][cylinder][property]")
{
    rc::prop("normal points outward, Reversed flips it, Project and Transformed invert", [] {
        const Cylinder c = TiltedCylinder();
        const Transform x = SomeSimilarity();
        const Point2 uv = RandomUv();
        const Point3 p = Eval(c, uv);
        const Vec3 normal = Cross(Deriv(c, uv, 1, 0), Deriv(c, uv, 0, 1));
        const Point3 foot = c.Origin + uv.Y * Stereon::Axis(c);
        RC_ASSERT(Dot(normal, p - foot) > 0.0);

        const Point2 mirrored{CirclePeriod - uv.X, uv.Y};
        const Cylinder r = Reversed(c);
        RC_ASSERT(Near(Eval(r, mirrored), p));
        RC_ASSERT(Dot(Cross(Deriv(r, mirrored, 1, 0), Deriv(r, mirrored, 0, 1)), normal) < 0.0);

        const Point2 back = Project(c, p);
        const double gap = std::abs(back.X - uv.X);
        RC_ASSERT(std::min(gap, CirclePeriod - gap) <= 1e-9);
        RC_ASSERT(std::abs(back.Y - uv.Y) <= 1e-9);
        RC_ASSERT(Near(Eval(Transformed(c, x), Point2{uv.X, Scale * uv.Y}), Transformed(p, x)));
    });

    // A reflection keeps the axis the image of the old axis and the normal
    // the image of the old normal.
    {
        const Cylinder c = TiltedCylinder();
        const Transform mirror = Stereon::MakeScaling(Vec3{-1.0, 1.0, 1.0});
        const Cylinder m = Transformed(c, mirror);
        REQUIRE(Near(Stereon::Axis(m), Transformed(Stereon::Axis(c), mirror)));
        const Point2 uv{0.7, 1.5};
        REQUIRE(Near(Eval(m, uv), Transformed(Eval(c, uv), mirror)));
        REQUIRE(m.RightHanded != c.RightHanded);
    }
    const Cylinder c = TiltedCylinder();
    REQUIRE(Domain(c) ==
            Stereon::Box2{Stereon::Interval{0.0, CirclePeriod}, Stereon::Interval::Everything()});
    REQUIRE_THAT(Length(Eval(c, Point2{1.0, 0.0}) - c.Origin), WithinAbs(c.Radius, Tolerance));
    REQUIRE(Project(c, c.Origin).X == 0.0); // a point on the axis
}

TEST_CASE("Surface wraps any surface and answers as the surface does", "[surface][erased]")
{
    const Cylinder cylinder = TiltedCylinder();
    const Surface surface(cylinder);
    const Point2 uv{0.8, -1.25};
    REQUIRE(Eval(surface, uv) == Eval(cylinder, uv));
    REQUIRE(Deriv(surface, uv, 2, 0) == Deriv(cylinder, uv, 2, 0));
    REQUIRE(Project(surface, Eval(cylinder, uv)) == Project(cylinder, Eval(cylinder, uv)));
    REQUIRE(Domain(surface) == Domain(cylinder));
    REQUIRE(BoundingBox(surface) == BoundingBox(cylinder));
    REQUIRE(*surface.Target<Cylinder>() == cylinder);
    REQUIRE(surface.Target<Plane>() == nullptr);
    REQUIRE(*Reversed(surface).Target<Cylinder>() == Reversed(cylinder));
    REQUIRE(*Transformed(surface, SomeSimilarity()).Target<Cylinder>() ==
            Transformed(cylinder, SomeSimilarity()));

    Surface copy = surface;
    REQUIRE(copy.Target<Cylinder>() != surface.Target<Cylinder>()); // small: copied
    copy = Surface(TiltedPlane());
    REQUIRE(copy.Target<Plane>() != nullptr);
    REQUIRE(Surface{}.IsEmpty());
}
