// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Rectangular trimmed surfaces: a basis restricted to a box of its
// parameter plane, evaluated as the basis, projected onto the patch rather
// than the basis, and converted, moved and bounded as the patch.

#include "RapidCheck.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string_view>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::Box2;
using Stereon::Box3;
using Stereon::BSplineSurface;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Transform;
using Stereon::TrimmedSurface;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-10;
constexpr double Tau = 2.0 * std::numbers::pi;

Cylinder SomeCylinder()
{
    return {Point3{1.0, -2.0, 0.5}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Plane SomePlane()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 0.0, 1.0}, 0.6);
    return {Point3{0.0, -1.0, 0.5}, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

/// A wavy bicubic over [0, 1]^2.
BSplineSurface SomeBSpline()
{
    std::vector<Point3> poles;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double z = ((i + j) % 2 == 0 ? 0.8 : -0.8) * (i == 0 || i == 3 ? 0.5 : 1.0);
            poles.push_back({i * 1.0, j * 1.0, z});
        }
    }
    return {3, 3, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}, poles, {}};
}

/// A surface type that is not built in: a plane under another name, with
/// no NURBS form the kernel knows.
struct UserPlane {
    Plane Base;
};
Point3 Eval(const UserPlane& s, const Point2& uv)
{
    return Eval(s.Base, uv);
}
Vec3 Deriv(const UserPlane& s, const Point2& uv, const int du, const int dv)
{
    return Deriv(s.Base, uv, du, dv);
}
Box2 Domain(const UserPlane& s)
{
    return Domain(s.Base);
}
Box3 BoundingBox(const UserPlane& s)
{
    return BoundingBox(s.Base);
}
Point2 Project(const UserPlane& s, const Point3& p)
{
    return Project(s.Base, p);
}
UserPlane Reversed(const UserPlane& s)
{
    return {Reversed(s.Base)};
}
UserPlane Transformed(const UserPlane& s, const Transform& x)
{
    return {Transformed(s.Base, x)};
}

bool Near(const Point3& a, const Point3& b, const double tolerance = Tolerance)
{
    return Distance(a, b) <= tolerance;
}

/// The least distance from p to the surface over a grid of the box, the
/// referee for projections.
double SampledDistance(const Surface& s, const Box2& box, const Point3& p, const int n = 60)
{
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            // std::lerp reaches the upper bound exactly at the last sample.
            const Point2 uv{std::lerp(box.X.GetLower(), box.X.GetUpper(), 1.0 * i / n),
                            std::lerp(box.Y.GetLower(), box.Y.GetUpper(), 1.0 * j / n)};
            best = std::min(best, Distance(Eval(s, uv), p));
        }
    }
    return best;
}

Point2 Clamped(const Point2& uv, const Box2& box)
{
    return {std::clamp(uv.X, box.X.GetLower(), box.X.GetUpper()),
            std::clamp(uv.Y, box.Y.GetLower(), box.Y.GetUpper())};
}

} // namespace

TEST_CASE("a trimmed surface evaluates as its basis on its box", "[surface][trimmed]")
{
    const Surface basis(SomeCylinder());
    const Box2 box{Interval{0.2, 1.5}, Interval{-1.0, 3.0}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    REQUIRE(Domain(t) == box);
    for (const double u : {0.2, 0.7, 1.5}) {
        for (const double v : {-1.0, 0.25, 3.0}) {
            const Point2 uv{u, v};
            REQUIRE(Eval(t, uv) == Eval(basis, uv));
            REQUIRE(Deriv(t, uv, 1, 0) == Deriv(basis, uv, 1, 0));
            REQUIRE(Deriv(t, uv, 0, 1) == Deriv(basis, uv, 0, 1));
            REQUIRE(Deriv(t, uv, 2, 1) == Deriv(basis, uv, 2, 1));
        }
    }
    // Not periodic, whatever the basis repeats.
    const Surface erased(t);
    REQUIRE(Periods(erased) == Vec2{});
    REQUIRE(Domain(erased) == box);
    REQUIRE(erased.Target<TrimmedSurface>() != nullptr);
    STATIC_REQUIRE(Stereon::GeometryTypeName<TrimmedSurface> == std::string_view{"TrimmedSurface"});
}

TEST_CASE("a trimmed surface projects onto its patch, not its basis", "[surface][trimmed]")
{
    // A wavy B-spline restricted to a corner of its domain: the basis's
    // closest point to a point above its middle lies outside the box, and
    // clamping its parameters into the box gives a point farther than the
    // patch's own closest point.
    const Surface basis(SomeBSpline());
    const Box2 box{Interval{0.0, 0.35}, Interval{0.0, 0.35}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    REQUIRE(t.GetPatch() != nullptr);
    const Point3 p{1.1, 1.3, 2.0};
    const Point2 global = Project(basis, p);
    REQUIRE_FALSE(Contains(box, global));
    const Point2 found = Project(t, p);
    REQUIRE(Contains(box, found));
    const double distance = Distance(Eval(t, found), p);
    REQUIRE(distance <= SampledDistance(Surface(t), box, p) + 1e-6);
    REQUIRE(distance <= Distance(Eval(basis, Clamped(global, box)), p));
}

TEST_CASE("trimmed surfaces project to the closest point of the patch",
          "[surface][trimmed][property]")
{
    rc::prop("for cylinders, planes and B-splines, points near and far", [] {
        const int kind = *rc::gen::inRange(0, 3);
        const Surface basis = kind == 0   ? Surface(SomeCylinder())
                              : kind == 1 ? Surface(SomePlane())
                                          : Surface(SomeBSpline());
        const Box2 domain = Domain(basis);
        const auto sub = [&](const Interval& d, const double period) {
            const double low = period > 0.0 ? *rc::gen::inRange(-200, 400) / 100.0
                                            : (std::isfinite(d.GetLower()) ? d.GetLower() : -3.0);
            const double top = period > 0.0 ? low + *rc::gen::inRange(20, 300) / 100.0
                                            : (std::isfinite(d.GetUpper()) ? d.GetUpper() : 3.0);
            const double a = low + (top - low) * *rc::gen::inRange(0, 50) / 100.0;
            const double b = a + (top - a) * *rc::gen::inRange(20, 100) / 100.0;
            return Interval{a, b};
        };
        const Vec2 periods = Periods(basis);
        const Box2 box{sub(domain.X, periods.X), sub(domain.Y, periods.Y)};
        const auto made = Stereon::MakeTrimmedSurface(basis, box);
        RC_ASSERT(made.has_value());
        const TrimmedSurface t = *made;
        const Point3 p{*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                       *rc::gen::inRange(-50, 50) / 10.0};
        const Point2 found = Project(t, p);
        RC_ASSERT(Contains(box, found));
        const double distance = Distance(Eval(t, found), p);
        RC_ASSERT(distance <= SampledDistance(Surface(t), box, p) + 1e-5);
        // Never worse than clamping the basis's projection into the box.
        const Point2 clamped = Clamped(Project(basis, p), box);
        RC_ASSERT(distance <= Distance(Eval(basis, clamped), p) + Tolerance);
    });
}

TEST_CASE("a point of a patch's edge projects to the edge, not a corner", "[surface][trimmed]")
{
    // The cylinder's own projection of the patch's edge point comes back a
    // rounding error above the box's upper u. Wrapping that down by a
    // period and clamping would land on the lower u, the far end of the
    // arc, and leave the corner as the best candidate.
    const Surface basis(SomeCylinder());
    const Box2 box{Interval{0.2175, 0.5643}, Interval{-3.0, -0.48}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    const Point3 p{-5.0, -4.4, -2.4};
    const Point2 found = Project(t, p);
    REQUIRE(found.X == 0.5643);
    REQUIRE(found.Y == Catch::Approx(-2.9).margin(1e-9));
    REQUIRE(Distance(Eval(t, found), p) <= SampledDistance(Surface(t), box, p) + 1e-9);
}

TEST_CASE("a cylinder trimmed across its seam is one patch", "[surface][trimmed]")
{
    const Surface basis(SomeCylinder());
    const Box2 box{Interval{5.0, 7.0}, Interval{0.0, 1.0}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    REQUIRE(t.GetPatch() != nullptr);
    const Point2 uv{6.5, 0.4};
    REQUIRE(Eval(t, uv) == Eval(basis, uv));
    // The point at u = 6.5 projects to a parameter in the box, not to the
    // basis's 6.5 - 2 pi.
    const Point2 found = Project(t, Eval(basis, uv) + Vec3{0.0, 0.0, 0.0});
    REQUIRE(Contains(box, found));
    REQUIRE(Near(Eval(t, found), Eval(basis, uv), 1e-8));
    REQUIRE(Stereon::MakeTrimmedSurface(basis, Box2{Interval{5.0, 5.0 + Tau + 0.1}, box.Y})
                .error()
                .GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("a reversed trimmed surface runs u backwards", "[surface][trimmed]")
{
    const Surface basis(SomeBSpline());
    const Box2 box{Interval{0.1, 0.9}, Interval{0.2, 0.7}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    const TrimmedSurface r = Reversed(t);
    REQUIRE(r.IsReversed());
    REQUIRE(Domain(r) == box);
    for (const double u : {0.1, 0.4, 0.9}) {
        for (const double v : {0.2, 0.55, 0.7}) {
            const Point2 uv{u, v};
            // 1 - 0.9 is an ulp below 0.1, outside the box: clamp it in.
            const Point2 reflected = Clamped(Point2{1.0 - u, v}, box);
            REQUIRE(Near(Eval(r, reflected), Eval(t, uv)));
            REQUIRE(Length(Deriv(r, reflected, 1, 0) + Deriv(t, uv, 1, 0)) <= Tolerance);
            REQUIRE(Length(Deriv(r, reflected, 0, 1) - Deriv(t, uv, 0, 1)) <= Tolerance);
            REQUIRE(Length(Deriv(r, reflected, 2, 0) - Deriv(t, uv, 2, 0)) <= Tolerance);
            const Vec3 n = *Stereon::UnitNormal(t, uv);
            const Vec3 m = *Stereon::UnitNormal(r, reflected);
            REQUIRE(Dot(n, m) <= -1.0 + 1e-9);
        }
    }
    REQUIRE(Reversed(r).IsReversed() == false);
    REQUIRE(r.Reflect(0.1) == 0.9);
    REQUIRE(r.Reflect(0.9) == 0.1);
    // Projection comes back in the reversed parametrisation.
    const Point3 p{0.5, 0.5, 3.0};
    const Point2 onT = Project(t, p);
    const Point2 onR = Project(r, p);
    REQUIRE(Near(Eval(r, onR), Eval(t, onT), 1e-8));
}

TEST_CASE("a trimmed surface converts, moves and bounds as its patch", "[surface][trimmed]")
{
    SECTION("a B-spline basis keeps its parameters")
    {
        const Surface basis(SomeBSpline());
        const Box2 box{Interval{0.25, 0.75}, Interval{0.0, 0.5}};
        const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
        const BSplineSurface patch = *Stereon::ToBSpline(t);
        REQUIRE(Domain(patch) == box);
        const BSplineSurface part = *Stereon::ToBSpline(t, Box2{Interval{0.3, 0.6}, box.Y});
        for (const double u : {0.3, 0.45, 0.6}) {
            for (const double v : {0.0, 0.3, 0.5}) {
                REQUIRE(Near(Eval(patch, Point2{u, v}), Eval(t, Point2{u, v})));
                REQUIRE(Near(Eval(part, Point2{u, v}), Eval(t, Point2{u, v})));
            }
        }
        const BSplineSurface reversed = *Stereon::ToBSpline(Reversed(t));
        REQUIRE(Near(Eval(reversed, Point2{0.3, 0.1}), Eval(Reversed(t), Point2{0.3, 0.1})));
    }
    SECTION("a cylinder's patch holds the same points")
    {
        const Surface basis(SomeCylinder());
        const Box2 box{Interval{0.5, 2.5}, Interval{-1.0, 1.0}};
        const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
        const BSplineSurface patch = *Stereon::ToBSpline(t);
        const Box3 bound = BoundingBox(t);
        for (const double u : {0.5, 1.2, 2.5}) {
            for (const double v : {-1.0, 0.0, 1.0}) {
                const Point3 p = Eval(t, Point2{u, v});
                REQUIRE(Contains(bound, p));
                REQUIRE(Near(Eval(patch, Project(patch, p)), p, 1e-8));
            }
        }
    }
    SECTION("moved")
    {
        const Surface basis(SomeCylinder());
        const Box2 box{Interval{0.5, 2.5}, Interval{-1.0, 1.0}};
        const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
        const Transform x = Stereon::MakeRotation(Vec3{0.3, 1.0, -0.2}, 0.8);
        const TrimmedSurface moved = Transformed(t, x);
        REQUIRE(moved.GetRange() == box);
        for (const double u : {0.5, 1.7}) {
            for (const double v : {-1.0, 0.4}) {
                REQUIRE(Near(Eval(moved, Point2{u, v}), Transformed(Eval(t, Point2{u, v}), x)));
            }
        }
    }
    SECTION("scaled: the range follows the parameters the similarity rescales")
    {
        // A plane's parameters are distances along its axes, which its
        // Transformed renormalises: the unit patch scaled by two is the
        // patch over [0, 2]^2, its far corner at (2, 2, 0), not (1, 1, 0).
        // A cylinder's u is an angle and its v a distance along the axis: v
        // follows the scale, u does not. A sphere's are both angles.
        const Transform scaling = Stereon::MakeScaling(Point3{}, 2.0);
        const TrimmedSurface unit = *Stereon::MakeTrimmedSurface(
            Surface(Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}),
            Box2{Interval{0.0, 1.0}, Interval{0.0, 1.0}});
        const TrimmedSurface plane = Transformed(unit, scaling);
        REQUIRE(plane.GetRange() == Box2{Interval{0.0, 2.0}, Interval{0.0, 2.0}});
        REQUIRE(Near(Eval(plane, Point2{2.0, 2.0}), Point3{2.0, 2.0, 0.0}));
        REQUIRE(Near(Eval(plane, Point2{2.0, 2.0}),
                     Transformed(Eval(unit, Point2{1.0, 1.0}), scaling)));
        const Box2 box{Interval{0.5, 2.5}, Interval{-1.0, 1.0}};
        const TrimmedSurface t = *Stereon::MakeTrimmedSurface(Surface(SomeCylinder()), box);
        const TrimmedSurface drum = Transformed(t, scaling);
        REQUIRE(drum.GetRange() == Box2{box.X, Interval{-2.0, 2.0}});
        for (const double u : {0.5, 1.7, 2.5}) {
            for (const double v : {-1.0, 0.4, 1.0}) {
                REQUIRE(Near(Eval(drum, Point2{u, 2.0 * v}),
                             Transformed(Eval(t, Point2{u, v}), scaling)));
            }
        }
        const TrimmedSurface ball = Transformed(
            *Stereon::MakeTrimmedSurface(
                Surface(Sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}), box),
            scaling);
        REQUIRE(ball.GetRange() == box);
    }
}

TEST_CASE("a user's basis type projects approximately and does not convert", "[surface][trimmed]")
{
    const Surface basis(UserPlane{SomePlane()});
    const Box2 box{Interval{-1.0, 1.0}, Interval{0.0, 2.0}};
    const TrimmedSurface t = *Stereon::MakeTrimmedSurface(basis, box);
    REQUIRE(t.GetPatch() == nullptr);
    REQUIRE(Stereon::ToBSpline(t).error().GetCode() == ErrorCode::NotSupported);
    // A plane's closest point clamped into the box is the patch's closest
    // point, so the approximation is exact here.
    const Point3 p{4.0, 3.0, -2.0};
    const Point2 found = Project(t, p);
    REQUIRE(found == Clamped(Project(basis, p), box));
    REQUIRE(Distance(Eval(t, found), p) <= SampledDistance(Surface(t), box, p) + 1e-6);
    REQUIRE(BoundingBox(t) == BoundingBox(basis));
}

TEST_CASE("MakeTrimmedSurface validates the range", "[surface][trimmed][validation]")
{
    const auto invalid = [](const auto& result) {
        return !result.has_value() && result.error().GetCode() == ErrorCode::InvalidInput;
    };
    const Surface spline(SomeBSpline());
    const Surface cylinder(SomeCylinder());
    constexpr double inf = std::numeric_limits<double>::infinity();
    REQUIRE(invalid(Stereon::MakeTrimmedSurface(Surface{}, Box2{Interval{0, 1}, Interval{0, 1}})));
    REQUIRE(invalid(Stereon::MakeTrimmedSurface(spline, Box2{Interval{0, 2}, Interval{0, 1}})));
    REQUIRE(invalid(Stereon::MakeTrimmedSurface(spline, Box2{Interval{0, 1}, Interval{-0.5, 1}})));
    REQUIRE(invalid(Stereon::MakeTrimmedSurface(spline, Box2{Interval{0.5, 0.5}, Interval{0, 1}})));
    REQUIRE(invalid(Stereon::MakeTrimmedSurface(cylinder, Box2{Interval{0, 1}, Interval{0, inf}})));
    REQUIRE(invalid(
        Stereon::MakeTrimmedSurface(cylinder, Box2{Interval{0, Tau + 1e-3}, Interval{0, 1}})));
    REQUIRE(Stereon::MakeTrimmedSurface(cylinder, Box2{Interval{-3.0, 3.0}, Interval{-9.0, 9.0}})
                .has_value());
    REQUIRE(Stereon::MakeTrimmedSurface(spline, Box2{Interval{0, 1}, Interval{0, 1}}).has_value());
    const Surface sphere(Sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    const Box2 sphereDomain = Domain(sphere);
    REQUIRE(Stereon::MakeTrimmedSurface(sphere, Box2{Interval{0, 1}, sphereDomain.Y}).has_value());
}
