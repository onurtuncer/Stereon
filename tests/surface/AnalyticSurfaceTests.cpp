// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The sphere, cone and torus. Each is checked against what it must be
// independently of its own code: derivatives against forward-mode duals of
// its evaluation, projections against closed-form distances computed a
// different way, points against the implicit equation of the surface.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

import Stereon.Surface;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicCone;
using Stereon::BasicSphere;
using Stereon::BasicTorus;
using Stereon::CirclePeriod;
using Stereon::Cone;
using Stereon::ErrorCode;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Transform;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;
constexpr double Scale = 1.5;

struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 2.0, -0.5}, 0.65);
    return {Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

Sphere TiltedSphere()
{
    const Frame f = TiltedFrame();
    return {Point3{1.0, -2.0, 0.5}, 2.5, f.X, f.Y};
}

Cone TiltedCone(const double semiAngle = 0.5)
{
    const Frame f = TiltedFrame();
    return {Point3{-1.0, 0.5, 1.0}, 1.25, semiAngle, f.X, f.Y};
}

Torus TiltedTorus()
{
    const Frame f = TiltedFrame();
    return {Point3{0.5, 1.0, -1.5}, 3.0, 0.75, f.X, f.Y};
}

Transform SomeSimilarity()
{
    return Stereon::MakeTranslation(Vec3{-1.0, 4.0, 0.5}) *
           Stereon::MakeRotation(Vec3{0.0, 1.0, 1.0}, 2.2) * Stereon::MakeScaling(Scale);
}

bool Near(const Point3& a, const Point3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(b));
}

/// Equal angles modulo 2 pi.
bool SameAngle(const double a, const double b)
{
    const double gap = std::fmod(std::abs(a - b), CirclePeriod);
    return std::min(gap, CirclePeriod - gap) <= 1e-9;
}

Stereon::BasicVec3<Dual> Lift(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

Stereon::BasicPoint3<Dual> Lift(const Point3& p)
{
    return {p.X, p.Y, p.Z};
}

BasicSphere<Dual> Lifted(const Sphere& s)
{
    return {Lift(s.Center), s.Radius, Lift(s.XAxis), Lift(s.YAxis), s.RightHanded};
}

BasicCone<Dual> Lifted(const Cone& c)
{
    return {Lift(c.Origin), c.Radius, c.SemiAngle, Lift(c.XAxis), Lift(c.YAxis), c.RightHanded};
}

BasicTorus<Dual> Lifted(const Torus& t)
{
    return {Lift(t.Center), t.MajorRadius, t.MinorRadius,
            Lift(t.XAxis),  Lift(t.YAxis), t.RightHanded};
}

template <typename P>
Vec3 DerivPart(const P& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

/// Partials in u and v of every order agree with the duals of the evaluator
/// one order lower.
template <typename S>
void CheckDerivatives(const S& s, const Point2& uv, const int n)
{
    const Stereon::BasicPoint2<Dual> du{Dual{uv.X, 1.0}, Dual{uv.Y}};
    const Stereon::BasicPoint2<Dual> dv{Dual{uv.X}, Dual{uv.Y, 1.0}};
    RC_ASSERT(Near(Deriv(s, uv, 1, 0), DerivPart(Eval(Lifted(s), du))));
    RC_ASSERT(Near(Deriv(s, uv, 0, 1), DerivPart(Eval(Lifted(s), dv))));
    RC_ASSERT(Near(Deriv(s, uv, n + 1, 0), DerivPart(Deriv(Lifted(s), du, n, 0))));
    RC_ASSERT(Near(Deriv(s, uv, n, 1), DerivPart(Deriv(Lifted(s), dv, n, 0))));
    RC_ASSERT(Near(Deriv(s, uv, 0, n + 1), DerivPart(Deriv(Lifted(s), dv, 0, n))));
    RC_ASSERT(Near(Deriv(s, uv, 1, n), DerivPart(Deriv(Lifted(s), du, 0, n))));
}

/// The foot of the projection is a stationary point of the distance:
/// p - S(uv) is normal to both partials.
template <typename S>
bool IsStationary(const S& s, const Point3& p, const Point2& uv)
{
    const Vec3 gap = p - Eval(s, uv);
    const Vec3 su = Deriv(s, uv, 1, 0);
    const Vec3 sv = Deriv(s, uv, 0, 1);
    const double scale = 1e-9 * (1.0 + Length(gap));
    return std::abs(Dot(gap, su)) <= scale * (1.0 + Length(su)) &&
           std::abs(Dot(gap, sv)) <= scale * (1.0 + Length(sv));
}

Point3 RandomPoint()
{
    return {*rc::gen::inRange(-8000, 8000) / 1000.0, *rc::gen::inRange(-8000, 8000) / 1000.0,
            *rc::gen::inRange(-8000, 8000) / 1000.0};
}

double RandomAngle()
{
    return *rc::gen::inRange(0, 6283) / 1000.0;
}

/// p in the frame of a surface of revolution: radial distance from the axis
/// and height along it.
struct Meridian {
    double Rho;
    double H;
};

template <typename S>
Meridian MeridianOf(const S& s, const Point3& origin, const Point3& p)
{
    const Vec3 offset = p - origin;
    const Vec3 axis = Stereon::Axis(s);
    const double h = Dot(offset, axis);
    return {Length(offset - h * axis), h};
}

} // namespace

TEST_CASE("analytic surfaces satisfy the surface concept", "[surface][concepts]")
{
    STATIC_REQUIRE(Stereon::SurfaceGeometry<Sphere>);
    STATIC_REQUIRE(Stereon::SurfaceGeometry<Cone>);
    STATIC_REQUIRE(Stereon::SurfaceGeometry<Torus>);
    STATIC_REQUIRE(Stereon::PeriodicSurface<Sphere>);
    STATIC_REQUIRE(Stereon::PeriodicSurface<Cone>);
    STATIC_REQUIRE(Stereon::PeriodicSurface<Torus>);
    // Every built-in surface is stored in the wrapper's own buffer.
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Sphere, Stereon::SurfaceBufferSize>);
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Cone, Stereon::SurfaceBufferSize>);
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Torus, Stereon::SurfaceBufferSize>);
}

TEST_CASE("a sphere is evaluated at its poles and equator", "[surface][sphere]")
{
    const Sphere s{Point3{1.0, 2.0, 3.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Eval(s, Point2{0.0, 0.0}) == Point3{3.0, 2.0, 3.0});
    REQUIRE(Near(Eval(s, Point2{std::numbers::pi / 2, 0.0}), Point3{1.0, 4.0, 3.0}));
    // Every longitude meets at a pole.
    for (int i = 0; i < 8; ++i) {
        REQUIRE(Near(Eval(s, Point2{i * 0.8, Stereon::SpherePoleLatitude}), Point3{1.0, 2.0, 5.0}));
        REQUIRE(
            Near(Eval(s, Point2{i * 0.8, -Stereon::SpherePoleLatitude}), Point3{1.0, 2.0, 1.0}));
    }
    REQUIRE(Domain(s) == Stereon::Box2{Stereon::Interval{0.0, CirclePeriod},
                                       Stereon::Interval{-Stereon::SpherePoleLatitude,
                                                         Stereon::SpherePoleLatitude}});
    REQUIRE(Periods(s) == Stereon::Vec2{CirclePeriod, 0.0});
    // The centre and points on the axis project as documented.
    REQUIRE(Project(s, s.Center) == Point2{0.0, 0.0});
    REQUIRE(Project(s, Point3{1.0, 2.0, 9.0}).X == 0.0);
    REQUIRE_THAT(Project(s, Point3{1.0, 2.0, 9.0}).Y, WithinAbs(std::numbers::pi / 2, Tolerance));
    REQUIRE_THAT(Project(s, Point3{1.0, 2.0, -9.0}).Y, WithinAbs(-std::numbers::pi / 2, Tolerance));
    // The bounding box is tight.
    const Stereon::Box3 box = BoundingBox(s);
    REQUIRE_THAT(box.X.GetLower(), WithinAbs(-1.0, Tolerance));
    REQUIRE_THAT(box.Z.GetUpper(), WithinAbs(5.0, Tolerance));
}

TEST_CASE("sphere operations agree with each other", "[surface][sphere][property]")
{
    rc::prop("derivatives are the derivatives of the evaluation", [] {
        const Point2 uv{RandomAngle(), *rc::gen::inRange(-1570, 1570) / 1000.0};
        CheckDerivatives(TiltedSphere(), uv, *rc::gen::inRange(1, 6));
    });

    rc::prop("points lie on the sphere and the normal points outward", [] {
        const Sphere s = TiltedSphere();
        const Point2 uv{RandomAngle(), *rc::gen::inRange(-1500, 1500) / 1000.0};
        const Point3 p = Eval(s, uv);
        RC_ASSERT(std::abs(Distance(p, s.Center) - s.Radius) <= Tolerance * s.Radius);
        RC_ASSERT(Contains(BoundingBox(s), p));
        RC_ASSERT(Dot(Cross(Deriv(s, uv, 1, 0), Deriv(s, uv, 0, 1)), p - s.Center) > 0.0);
        // Reversed runs u backwards and turns the normal round.
        const Sphere r = Reversed(s);
        const Point2 mirrored{CirclePeriod - uv.X, uv.Y};
        RC_ASSERT(Near(Eval(r, mirrored), p));
        RC_ASSERT(Dot(Cross(Deriv(r, mirrored, 1, 0), Deriv(r, mirrored, 0, 1)), p - s.Center) <
                  0.0);
        const Transform x = SomeSimilarity();
        RC_ASSERT(Near(Eval(Transformed(s, x), uv), Transformed(p, x)));
    });

    rc::prop("Project inverts Eval and finds the closest point", [] {
        const Sphere s = TiltedSphere();
        const Point2 uv{RandomAngle(), *rc::gen::inRange(-1500, 1500) / 1000.0};
        const Point2 back = Project(s, Eval(s, uv));
        RC_ASSERT(SameAngle(back.X, uv.X));
        RC_ASSERT(std::abs(back.Y - uv.Y) <= 1e-9);

        const Point3 p = RandomPoint();
        const Point2 foot = Project(s, p);
        RC_ASSERT(Contains(Domain(s), foot));
        const double expected = std::abs(Distance(p, s.Center) - s.Radius);
        RC_ASSERT(std::abs(Distance(Eval(s, foot), p) - expected) <= 1e-12 * (1.0 + expected));
    });
}

TEST_CASE("a cone is evaluated along its generators", "[surface][cone]")
{
    const double angle = std::numbers::pi / 4;
    const Cone c{Point3{}, 1.0, angle, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Eval(c, Point2{0.0, 0.0}) == Point3{1.0, 0.0, 0.0});
    // v is arc length along the generator.
    REQUIRE_THAT(Distance(Eval(c, Point2{0.3, 2.0}), Eval(c, Point2{0.3, 0.0})),
                 WithinAbs(2.0, Tolerance));
    // Every u meets at the apex, at v = -Radius / sin a.
    const double apexV = -1.0 / std::sin(angle);
    REQUIRE(Near(Stereon::Apex(c), Point3{0.0, 0.0, -1.0}));
    for (int i = 0; i < 8; ++i) {
        REQUIRE(Near(Eval(c, Point2{i * 0.8, apexV}), Stereon::Apex(c)));
    }
    // Higher derivatives in v vanish; the cone is ruled.
    REQUIRE(Deriv(c, Point2{0.4, 1.0}, 0, 2) == Vec3{});
    REQUIRE(Deriv(c, Point2{0.4, 1.0}, 2, 3) == Vec3{});
    REQUIRE(Domain(c) ==
            Stereon::Box2{Stereon::Interval{0.0, CirclePeriod}, Stereon::Interval::Everything()});
    REQUIRE(BoundingBox(c) == Stereon::Box3::Everything());
    // A point on the axis gets u = 0; beyond the apex it is nearest the
    // second nappe.
    REQUIRE(Project(c, Point3{0.0, 0.0, 5.0}).X == 0.0);
    REQUIRE(Project(c, Point3{0.0, 0.0, -3.0}).Y < apexV);
}

TEST_CASE("cone operations agree with each other", "[surface][cone][property]")
{
    rc::prop("derivatives are the derivatives of the evaluation", [] {
        const Point2 uv{RandomAngle(), *rc::gen::inRange(-5000, 5000) / 1000.0};
        CheckDerivatives(TiltedCone(*rc::gen::inRange(-1500, 1500) / 1000.0 + 0.01), uv,
                         *rc::gen::inRange(1, 6));
    });

    rc::prop("points satisfy the implicit equation, Reversed and Transformed commute", [] {
        const Cone c = TiltedCone();
        const Point2 uv{RandomAngle(), *rc::gen::inRange(-5000, 5000) / 1000.0};
        const Point3 p = Eval(c, uv);
        // rho = |Radius + (h / cos a) sin a|: the distance from the axis grows
        // linearly with the height.
        const Meridian m = MeridianOf(c, c.Origin, p);
        RC_ASSERT(std::abs(m.Rho - std::abs(c.Radius + m.H * std::tan(c.SemiAngle))) <= 1e-12 * 10);
        const Point2 mirrored{CirclePeriod - uv.X, uv.Y};
        RC_ASSERT(Near(Eval(Reversed(c), mirrored), p));
        const Transform x = SomeSimilarity();
        RC_ASSERT(Near(Eval(Transformed(c, x), Point2{uv.X, Scale * uv.Y}), Transformed(p, x)));
    });

    rc::prop("Project finds the closest point of either nappe", [] {
        const Cone c = TiltedCone();
        const Point3 p = RandomPoint();
        const Point2 foot = Project(c, p);
        // In the meridian plane the cone is two lines through the apex, at
        // angles +-a from the axis; the distance to a line through the apex
        // is |rho cos a -+ (h - hApex) sin a|.
        const Meridian m = MeridianOf(c, Stereon::Apex(c), p);
        const double s = std::sin(c.SemiAngle);
        const double k = std::cos(c.SemiAngle);
        const double expected =
            std::min(std::abs(m.Rho * k - m.H * s), std::abs(m.Rho * k + m.H * s));
        RC_ASSERT(std::abs(Distance(Eval(c, foot), p) - expected) <= 1e-11 * (1.0 + expected));
        RC_ASSERT(IsStationary(c, p, foot));
        RC_ASSERT(foot.X >= 0.0 && foot.X < CirclePeriod);
    });
}

TEST_CASE("a torus is evaluated around its tube", "[surface][torus]")
{
    const Torus t{Point3{}, 3.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Eval(t, Point2{0.0, 0.0}) == Point3{4.0, 0.0, 0.0});
    REQUIRE(Near(Eval(t, Point2{0.0, std::numbers::pi}), Point3{2.0, 0.0, 0.0}));
    REQUIRE(Near(Eval(t, Point2{0.0, std::numbers::pi / 2}), Point3{3.0, 0.0, 1.0}));
    REQUIRE(Periods(t) == Stereon::Vec2{CirclePeriod, CirclePeriod});
    // The box is tight: outer radius 4 in the plane, tube radius 1 along z.
    const Stereon::Box3 box = BoundingBox(t);
    REQUIRE_THAT(box.X.GetUpper(), WithinAbs(4.0, Tolerance));
    REQUIRE_THAT(box.Z.GetLower(), WithinAbs(-1.0, Tolerance));
    // On the axis u = 0; on the centre circle v = 0.
    REQUIRE(Project(t, Point3{0.0, 0.0, 2.0}).X == 0.0);
    REQUIRE(Project(t, Point3{0.0, 3.0, 0.0}).Y == 0.0);
}

TEST_CASE("torus operations agree with each other", "[surface][torus][property]")
{
    rc::prop("derivatives are the derivatives of the evaluation", [] {
        CheckDerivatives(TiltedTorus(), Point2{RandomAngle(), RandomAngle()},
                         *rc::gen::inRange(1, 6));
    });

    rc::prop("points satisfy the implicit equation; the normal leaves the tube", [] {
        const Torus t = TiltedTorus();
        const Point2 uv{RandomAngle(), RandomAngle()};
        const Point3 p = Eval(t, uv);
        const Meridian m = MeridianOf(t, t.Center, p);
        RC_ASSERT(std::abs(std::hypot(m.Rho - t.MajorRadius, m.H) - t.MinorRadius) <= 1e-12 * 10);
        RC_ASSERT(Contains(BoundingBox(t), p));
        const Point3 core =
            t.Center + t.MajorRadius * Stereon::RadialDirection(t.XAxis, t.YAxis, uv.X, 0);
        RC_ASSERT(Dot(Cross(Deriv(t, uv, 1, 0), Deriv(t, uv, 0, 1)), p - core) > 0.0);
        RC_ASSERT(Near(Eval(Reversed(t), Point2{CirclePeriod - uv.X, uv.Y}), p));
        const Transform x = SomeSimilarity();
        RC_ASSERT(Near(Eval(Transformed(t, x), uv), Transformed(p, x)));
    });

    rc::prop("Project inverts Eval and finds the closest point", [] {
        const Torus t = TiltedTorus();
        const Point2 uv{RandomAngle(), RandomAngle()};
        const Point2 back = Project(t, Eval(t, uv));
        RC_ASSERT(SameAngle(back.X, uv.X) && SameAngle(back.Y, uv.Y));

        const Point3 p = RandomPoint();
        const Point2 foot = Project(t, p);
        const Meridian m = MeridianOf(t, t.Center, p);
        const double expected = std::abs(std::hypot(m.Rho - t.MajorRadius, m.H) - t.MinorRadius);
        RC_ASSERT(std::abs(Distance(Eval(t, foot), p) - expected) <= 1e-11 * (1.0 + expected));
        RC_ASSERT(Contains(Domain(t), foot));
    });
}

TEST_CASE("reflections keep the axis and the normal the images of the old ones",
          "[surface][property]")
{
    const Transform mirror = Stereon::MakeScaling(Vec3{-1.0, 1.0, 1.0});
    const Point2 uv{0.7, 0.4};
    const Sphere s = Transformed(TiltedSphere(), mirror);
    const Cone c = Transformed(TiltedCone(), mirror);
    const Torus t = Transformed(TiltedTorus(), mirror);
    REQUIRE(Near(Eval(s, uv), Transformed(Eval(TiltedSphere(), uv), mirror)));
    REQUIRE(Near(Eval(c, uv), Transformed(Eval(TiltedCone(), uv), mirror)));
    REQUIRE(Near(Eval(t, uv), Transformed(Eval(TiltedTorus(), uv), mirror)));
    REQUIRE(Near(Stereon::Axis(t), Transformed(Stereon::Axis(TiltedTorus()), mirror)));
    REQUIRE_FALSE(s.RightHanded);
}

TEST_CASE("the wrapper holds every analytic surface", "[surface][erased]")
{
    const Point2 uv{1.0, 0.5};
    const Point3 p{1.0, 1.0, 1.0};
    const Surface sphere(TiltedSphere());
    const Surface cone(TiltedCone());
    const Surface torus(TiltedTorus());
    REQUIRE(Eval(sphere, uv) == Eval(TiltedSphere(), uv));
    REQUIRE(Deriv(cone, uv, 1, 1) == Deriv(TiltedCone(), uv, 1, 1));
    REQUIRE(Project(torus, p) == Project(TiltedTorus(), p));
    REQUIRE(Periods(torus) == Stereon::Vec2{CirclePeriod, CirclePeriod});
    const Surface reversedCone = Reversed(cone);
    const Cone* reversed = reversedCone.Target<Cone>();
    REQUIRE(reversed != nullptr);
    REQUIRE(*reversed == Reversed(TiltedCone()));
    // A torus accepts any finite (u, v) through the wrapper; a sphere is
    // periodic in u only.
    REQUIRE(Eval(torus, Point2{-7.0, 20.0}).X == Eval(TiltedTorus(), Point2{-7.0, 20.0}).X);
}

TEST_CASE("analytic surfaces validate data from outside the kernel", "[surface][factory]")
{
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    const Vec3 skew{0.5, 1.0, 0.0};
    const auto codeOf = [](const auto& result) { return result.error().GetCode(); };

    REQUIRE(Stereon::MakeSphere(Point3{}, 1.0, x, y).has_value());
    REQUIRE(codeOf(Stereon::MakeSphere(Point3{}, 0.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeSphere(Point3{}, 1.0, x, skew)) == ErrorCode::InvalidInput);

    REQUIRE(Stereon::MakeCone(Point3{}, 0.0, 0.3, x, y).has_value()); // apex at the origin
    REQUIRE(Stereon::MakeCone(Point3{}, 1.0, -0.3, x, y).has_value());
    REQUIRE(codeOf(Stereon::MakeCone(Point3{}, -1.0, 0.3, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeCone(Point3{}, 1.0, 0.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeCone(Point3{}, 1.0, std::numbers::pi / 2, x, y)) ==
            ErrorCode::InvalidInput);

    REQUIRE(Stereon::MakeTorus(Point3{}, 2.0, 1.0, x, y).has_value());
    REQUIRE(codeOf(Stereon::MakeTorus(Point3{}, 1.0, 1.0, x, y)) == ErrorCode::InvalidInput);
    REQUIRE(codeOf(Stereon::MakeTorus(Point3{}, 2.0, 0.0, x, y)) == ErrorCode::InvalidInput);

    REQUIRE(Stereon::MakePlane(Point3{}, x, y).has_value());
    REQUIRE(codeOf(Stereon::MakePlane(Point3{}, x, x)) == ErrorCode::InvalidInput);
    const auto leftHanded = Stereon::MakeCylinder(Point3{}, 1.0, x, y, false);
    REQUIRE(leftHanded.has_value());
    REQUIRE_FALSE(leftHanded.value().RightHanded);
    REQUIRE(codeOf(Stereon::MakeCylinder(Point3{INFINITY, 0.0, 0.0}, 1.0, x, y)) ==
            ErrorCode::InvalidInput);
}
