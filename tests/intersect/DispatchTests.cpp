// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <string>

import Stereon.Intersect;

using Stereon::Circle;
using Stereon::Context;
using Stereon::Curve;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Context ctx = Context::Default();

Curve SomeLine()
{
    return Curve(Line{Point3{}, Vec3{1.0, 0.0, 0.0}});
}

Curve SomeCircle()
{
    return Curve(Circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

Surface SomePlane()
{
    return Surface(Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

Surface SomeCylinder()
{
    return Surface(Cylinder{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

Surface SomeSphere()
{
    return Surface(Stereon::Sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

/// A curve type that is not built in.
struct UserCurve {
    Line Base;
};
Point3 Eval(const UserCurve& c, const double t)
{
    return Eval(c.Base, t);
}
Vec3 Deriv(const UserCurve& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const UserCurve& c)
{
    return Domain(c.Base);
}
Stereon::Box3 BoundingBox(const UserCurve& c)
{
    return BoundingBox(c.Base);
}
double Project(const UserCurve& c, const Point3& p)
{
    return Project(c.Base, p);
}
UserCurve Reversed(const UserCurve& c)
{
    return {Reversed(c.Base)};
}
UserCurve Transformed(const UserCurve& c, const Stereon::Transform& x)
{
    return {Transformed(c.Base, x)};
}

std::string MessageOf(const Stereon::KernelError& error)
{
    return std::string(error.GetMessage());
}

} // namespace

TEST_CASE("every built-in pair has a handler that names both types", "[intersect][dispatch]")
{
    STATIC_REQUIRE(Stereon::Detail::CurvePairTable.size() == Stereon::BuiltinCurveTypes::Size);
    STATIC_REQUIRE(Stereon::Detail::SurfacePairTable[1].size() ==
                   Stereon::BuiltinSurfaceTypes::Size);
    for (const auto& row : Stereon::Detail::CurvePairTable) {
        for (const auto handler : row) {
            // Compared outside REQUIRE: Catch2's wrapper around a pointer to a
            // function returning std::expected sends libc++ into a recursive
            // constraint check.
            const bool present = handler != nullptr;
            REQUIRE(present);
        }
    }

    // Every 3D curve pair has a handler: it answers, or fails closed.
    REQUIRE(Intersect(SomeLine(), SomeLine(), ctx).has_value());
    for (const auto& result :
         {Intersect(SomeLine(), SomeCircle(), ctx), Intersect(SomeCircle(), SomeLine(), ctx)}) {
        REQUIRE(result.has_value());
    }
    REQUIRE(Intersect(SomeCircle(), SomeCircle(), ctx).has_value());

    // Every pair of planes, cylinders, cones, spheres and tori has its
    // analytic handler (M3, v0.2); a B-spline surface against an analytic
    // one still names both types.
    REQUIRE(Intersect(SomePlane(), SomeSphere(), ctx).has_value());
    REQUIRE(Intersect(SomeSphere(), SomeCylinder(), ctx).has_value());
    constexpr std::size_t spline = Stereon::BuiltinSurfaceTypes::Size - 1;
    const auto splinePlane = Stereon::Detail::SurfacePairTable[spline][0](
        SomePlane(), Stereon::Box2::Everything(), SomePlane(), Stereon::Box2::Everything(), 0.0,
        ctx);
    REQUIRE(splinePlane.error().GetCode() == ErrorCode::NotImplemented);
    REQUIRE(MessageOf(splinePlane.error()).starts_with("BSplineSurface x Plane"));
    REQUIRE(Intersect(SomePlane(), SomeCylinder(), ctx).has_value());
    REQUIRE(Intersect(SomeCylinder(), SomePlane(), ctx).has_value());
    REQUIRE(Intersect(SomePlane(), SomePlane(), ctx).has_value());
    REQUIRE(Intersect(SomeCylinder(), SomeCylinder(), ctx).has_value());
}

TEST_CASE("every curve-surface pair has a handler; lines against planes, spheres and cylinders "
          "have theirs",
          "[intersect][dispatch]")
{
    STATIC_REQUIRE(Stereon::Detail::CurveSurfaceTable.size() == Stereon::BuiltinCurveTypes::Size);
    STATIC_REQUIRE(Stereon::Detail::CurveSurfaceTable[0].size() ==
                   Stereon::BuiltinSurfaceTypes::Size);
    for (const auto& row : Stereon::Detail::CurveSurfaceTable) {
        for (const auto handler : row) {
            const bool present = handler != nullptr;
            REQUIRE(present);
        }
    }
    REQUIRE(Intersect(SomeLine(), SomePlane(), ctx).has_value());
    REQUIRE(Intersect(SomeLine(), SomeSphere(), ctx).has_value());
    REQUIRE(Intersect(SomeLine(), SomeCylinder(), ctx).has_value());
    REQUIRE(Intersect(SomeCircle(), SomePlane(), ctx).has_value());
    const Surface torus(
        Stereon::Torus{Point3{}, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    REQUIRE(Intersect(SomeLine(), torus, ctx).has_value());
    // The ellipsoid (ADR-0048) has the handlers the elementary surfaces
    // have: a line, a circle, an ellipse and a B-spline.
    const Surface ellipsoid(
        Stereon::Ellipsoid{Point3{}, 3.0, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    REQUIRE(Intersect(SomeLine(), ellipsoid, ctx).has_value());
    REQUIRE(Intersect(SomeCircle(), ellipsoid, ctx).has_value());
    // A parabola has its handler too, over a bounded domain; over its whole,
    // unbounded one it is NotSupported, and the message names both types,
    // the curve first.
    const Curve parabola(
        Stereon::Parabola{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    REQUIRE(Intersect(parabola, Stereon::Interval{-2.0, 2.0}, ellipsoid, Domain(ellipsoid), ctx)
                .has_value());
    const auto whole = Intersect(parabola, SomePlane(), ctx);
    REQUIRE(whole.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(MessageOf(whole.error()).starts_with("Parabola x Plane"));
    REQUIRE(
        MessageOf(Intersect(parabola, ellipsoid, ctx).error()).starts_with("Parabola x Ellipsoid"));
}

TEST_CASE("empty inputs are invalid input", "[intersect][dispatch]")
{
    REQUIRE(Intersect(Curve{}, SomeLine(), ctx).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Intersect(SomePlane(), Surface{}, ctx).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Intersect(Curve{}, SomePlane(), ctx).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Intersect(SomeLine(), Surface{}, ctx).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("types that are not built in are not supported yet", "[intersect][dispatch]")
{
    const Curve user(UserCurve{Line{Point3{}, Vec3{0.0, 1.0, 0.0}}});
    const auto result = Intersect(user, SomeLine(), ctx);
    REQUIRE(result.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(Intersect(SomeCircle(), user, ctx).error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(Intersect(user, SomePlane(), ctx).error().GetCode() == ErrorCode::NotSupported);
}
