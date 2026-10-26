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

    // Line-line is implemented; the remaining 3D curve pairs still name
    // their types in NotImplemented diagnostics, ordered as called.
    REQUIRE(Intersect(SomeLine(), SomeLine(), ctx).has_value());
    const auto lineCircle = Intersect(SomeLine(), SomeCircle(), ctx);
    REQUIRE_FALSE(lineCircle.has_value());
    REQUIRE(lineCircle.error().GetCode() == ErrorCode::NotImplemented);
    REQUIRE(MessageOf(lineCircle.error()).starts_with("Line x Circle"));
    REQUIRE(
        MessageOf(Intersect(SomeCircle(), SomeLine(), ctx).error()).starts_with("Circle x Line"));
    REQUIRE(MessageOf(Intersect(SomeCircle(), SomeCircle(), ctx).error())
                .starts_with("Circle x Circle"));

    // Planes and cylinders have their analytic handlers (M3); spheres wait.
    const auto planeSphere = Intersect(SomePlane(), SomeSphere(), ctx);
    REQUIRE(planeSphere.error().GetCode() == ErrorCode::NotImplemented);
    REQUIRE(MessageOf(planeSphere.error()).starts_with("Plane x Sphere"));
    REQUIRE(MessageOf(Intersect(SomeSphere(), SomeCylinder(), ctx).error())
                .starts_with("Sphere x Cylinder"));
    REQUIRE(Intersect(SomePlane(), SomeCylinder(), ctx).has_value());
    REQUIRE(Intersect(SomeCylinder(), SomePlane(), ctx).has_value());
    REQUIRE(Intersect(SomePlane(), SomePlane(), ctx).has_value());
    REQUIRE(Intersect(SomeCylinder(), SomeCylinder(), ctx).has_value());
}

TEST_CASE("empty inputs are invalid input", "[intersect][dispatch]")
{
    REQUIRE(Intersect(Curve{}, SomeLine(), ctx).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Intersect(SomePlane(), Surface{}, ctx).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("types that are not built in are not supported yet", "[intersect][dispatch]")
{
    const Curve user(UserCurve{Line{Point3{}, Vec3{0.0, 1.0, 0.0}}});
    const auto result = Intersect(user, SomeLine(), ctx);
    REQUIRE(result.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(Intersect(SomeCircle(), user, ctx).error().GetCode() == ErrorCode::NotSupported);
}
