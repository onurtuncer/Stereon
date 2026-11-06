// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The ellipsoid (ADR-0048): its points satisfy its implicit equation and
// its derivatives are those of its evaluation; with equal semi-axes it is
// the sphere, point for point and projection for projection; Project
// inverts Eval on the surface and finds a closest point off it, checked
// against a dense sampling, with the medial cases (the centre, the axes,
// a point inside on a plane of symmetry) deterministic; and data from
// outside the kernel is validated.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

import Stereon.Surface;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicEllipsoid;
using Stereon::CirclePeriod;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Transform;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double Tolerance = 1e-12;

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};

struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 2.0, -0.5}, 0.65);
    return {Transformed(ex, tilt), Transformed(ey, tilt)};
}

Ellipsoid TiltedEllipsoid()
{
    const Frame f = TiltedFrame();
    return {Point3{1.0, -2.0, 0.5}, 3.0, 1.75, 1.25, f.X, f.Y};
}

Transform SomeSimilarity()
{
    return Stereon::MakeTranslation(Vec3{-1.0, 4.0, 0.5}) *
           Stereon::MakeRotation(Vec3{0.0, 1.0, 1.0}, 2.2) * Stereon::MakeScaling(1.5);
}

bool Near(const Point3& a, const Point3& b, const double within = Tolerance)
{
    return Length(a - b) <= within * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b)
{
    return Length(a - b) <= Tolerance * (1.0 + Length(b));
}

bool SameAngle(const double a, const double b)
{
    const double gap = std::fmod(std::abs(a - b), CirclePeriod);
    return std::min(gap, CirclePeriod - gap) <= 1e-9;
}

/// (x/A)² + (y/B)² + (z/C)² - 1 in the ellipsoid's frame.
double Implicit(const Ellipsoid& e, const Point3& p)
{
    const Vec3 d = p - e.Center;
    const double x = Dot(d, e.XAxis) / e.A;
    const double y = Dot(d, e.YAxis) / e.B;
    const double z = Dot(d, Axis(e)) / e.C;
    return x * x + y * y + z * z - 1.0;
}

Stereon::BasicVec3<Dual> Lift(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

BasicEllipsoid<Dual> Lifted(const Ellipsoid& e)
{
    return {{e.Center.X, e.Center.Y, e.Center.Z},
            e.A,
            e.B,
            e.C,
            Lift(e.XAxis),
            Lift(e.YAxis),
            e.RightHanded};
}

template <typename P>
Vec3 DerivPart(const P& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

Point3 RandomPoint()
{
    return {*rc::gen::inRange(-6000, 6000) / 1000.0, *rc::gen::inRange(-6000, 6000) / 1000.0,
            *rc::gen::inRange(-6000, 6000) / 1000.0};
}

Point2 RandomParameters()
{
    return {*rc::gen::inRange(0, 6283) / 1000.0, *rc::gen::inRange(-1500, 1500) / 1000.0};
}

/// The smallest distance from p to a dense grid of the ellipsoid's points:
/// no closer than the closest point, and within the grid's spacing of it.
double SampledDistance(const Ellipsoid& e, const Point3& p)
{
    double best = std::numeric_limits<double>::infinity();
    constexpr int nu = 360;
    constexpr int nv = 180;
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j <= nv; ++j) {
            const Point2 uv{CirclePeriod * i / nu, -pi / 2 + pi * j / nv};
            best = std::min(best, Distance(Eval(e, uv), p));
        }
    }
    return best;
}

} // namespace

TEST_CASE("the ellipsoid is a surface stored inline", "[surface][ellipsoid][concepts]")
{
    STATIC_REQUIRE(Stereon::SurfaceGeometry<Ellipsoid>);
    STATIC_REQUIRE(Stereon::PeriodicSurface<Ellipsoid>);
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Ellipsoid, Stereon::SurfaceBufferSize>);
    const Surface wrapped{TiltedEllipsoid()};
    REQUIRE(wrapped.Target<Ellipsoid>() != nullptr);
    REQUIRE(Near(Eval(wrapped, Point2{0.3, 0.2}), Eval(TiltedEllipsoid(), Point2{0.3, 0.2})));
    REQUIRE(Stereon::GeometryNameOf(static_cast<const Ellipsoid*>(nullptr)) == "Ellipsoid");
}

TEST_CASE("an ellipsoid is evaluated at its vertices and poles", "[surface][ellipsoid]")
{
    const Ellipsoid e{Point3{1.0, 2.0, 3.0}, 3.0, 2.0, 1.0, ex, ey};
    REQUIRE(Eval(e, Point2{0.0, 0.0}) == Point3{4.0, 2.0, 3.0});
    REQUIRE(Near(Eval(e, Point2{pi / 2, 0.0}), Point3{1.0, 4.0, 3.0}));
    REQUIRE(Near(Eval(e, Point2{pi, 0.0}), Point3{-2.0, 2.0, 3.0}));
    for (int i = 0; i < 8; ++i) {
        REQUIRE(Near(Eval(e, Point2{i * 0.8, Stereon::SpherePoleLatitude}), Point3{1.0, 2.0, 4.0}));
        REQUIRE(
            Near(Eval(e, Point2{i * 0.8, -Stereon::SpherePoleLatitude}), Point3{1.0, 2.0, 2.0}));
    }
    REQUIRE(Domain(e) == Stereon::Box2{Stereon::Interval{0.0, CirclePeriod},
                                       Stereon::Interval{-Stereon::SpherePoleLatitude,
                                                         Stereon::SpherePoleLatitude}});
    REQUIRE(Periods(e) == Stereon::Vec2{CirclePeriod, 0.0});
    // The bounding box is the support along each axis: tight.
    const Stereon::Box3 box = BoundingBox(e);
    REQUIRE_THAT(box.X.GetLower(), WithinAbs(-2.0, Tolerance));
    REQUIRE_THAT(box.X.GetUpper(), WithinAbs(4.0, Tolerance));
    REQUIRE_THAT(box.Y.GetLower(), WithinAbs(0.0, Tolerance));
    REQUIRE_THAT(box.Z.GetUpper(), WithinAbs(4.0, Tolerance));
}

TEST_CASE("ellipsoid operations agree with each other", "[surface][ellipsoid][property]")
{
    rc::prop("derivatives are the derivatives of the evaluation", [] {
        const Ellipsoid e = TiltedEllipsoid();
        const Point2 uv = RandomParameters();
        const int n = *rc::gen::inRange(1, 6);
        const Stereon::BasicPoint2<Dual> du{Dual{uv.X, 1.0}, Dual{uv.Y}};
        const Stereon::BasicPoint2<Dual> dv{Dual{uv.X}, Dual{uv.Y, 1.0}};
        RC_ASSERT(Near(Deriv(e, uv, 1, 0), DerivPart(Eval(Lifted(e), du))));
        RC_ASSERT(Near(Deriv(e, uv, 0, 1), DerivPart(Eval(Lifted(e), dv))));
        RC_ASSERT(Near(Deriv(e, uv, n + 1, 0), DerivPart(Deriv(Lifted(e), du, n, 0))));
        RC_ASSERT(Near(Deriv(e, uv, n, 1), DerivPart(Deriv(Lifted(e), dv, n, 0))));
        RC_ASSERT(Near(Deriv(e, uv, 0, n + 1), DerivPart(Deriv(Lifted(e), dv, 0, n))));
        RC_ASSERT(Near(Deriv(e, uv, 1, n), DerivPart(Deriv(Lifted(e), du, 0, n))));
    });

    rc::prop("points satisfy the implicit equation and the normal points outward", [] {
        const Ellipsoid e = TiltedEllipsoid();
        const Point2 uv = RandomParameters();
        const Point3 p = Eval(e, uv);
        RC_ASSERT(std::abs(Implicit(e, p)) <= 1e-13);
        RC_ASSERT(Contains(BoundingBox(e), p));
        const Vec3 normal = Cross(Deriv(e, uv, 1, 0), Deriv(e, uv, 0, 1));
        RC_ASSERT(Dot(normal, p - e.Center) > 0.0);
        const Ellipsoid r = Reversed(e);
        const Point2 mirrored{CirclePeriod - uv.X, uv.Y};
        RC_ASSERT(Near(Eval(r, mirrored), p));
        RC_ASSERT(Dot(Cross(Deriv(r, mirrored, 1, 0), Deriv(r, mirrored, 0, 1)), p - e.Center) <
                  0.0);
        const Transform x = SomeSimilarity();
        RC_ASSERT(Near(Eval(Transformed(e, x), uv), Transformed(p, x)));
    });

    rc::prop("Project inverts Eval", [] {
        const Ellipsoid e = TiltedEllipsoid();
        const Point2 uv = RandomParameters();
        const Point2 back = Project(e, Eval(e, uv));
        RC_ASSERT(SameAngle(back.X, uv.X));
        RC_ASSERT(std::abs(back.Y - uv.Y) <= 1e-9);
    });
}

TEST_CASE("an ellipsoid with equal semi-axes is the sphere", "[surface][ellipsoid][property]")
{
    rc::prop("point for point and projection for projection", [] {
        const Frame f = TiltedFrame();
        const double r = *rc::gen::inRange(1, 40) / 8.0;
        const Point3 c{0.5, -1.0, 2.0};
        const Ellipsoid e{c, r, r, r, f.X, f.Y};
        const Sphere s{c, r, f.X, f.Y};
        const Point2 uv = RandomParameters();
        RC_ASSERT(Near(Eval(e, uv), Eval(s, uv)));
        RC_ASSERT(Near(Deriv(e, uv, 1, 1), Deriv(s, uv, 1, 1)));
        const Point3 p = RandomPoint();
        const Point2 a = Project(e, p);
        const Point2 b = Project(s, p);
        RC_ASSERT(Near(Eval(e, a), Eval(s, b), 1e-10));
        RC_ASSERT(Project(e, c) == Project(s, c));
    });
}

TEST_CASE("Project finds a closest point of the ellipsoid", "[surface][ellipsoid][property]")
{
    rc::prop("no sample of the surface is closer", [] {
        const Ellipsoid e = TiltedEllipsoid();
        const Point3 p = RandomPoint();
        const Point2 foot = Project(e, p);
        RC_ASSERT(Contains(Domain(e), foot));
        const double found = Distance(Eval(e, foot), p);
        RC_ASSERT(found <= SampledDistance(e, p) + 1e-12 * (1.0 + found));
        // And it is the same every time.
        RC_ASSERT(Project(e, p) == foot);
    });
}

TEST_CASE("the medial cases of Project are deterministic", "[surface][ellipsoid]")
{
    const Ellipsoid e{Point3{}, 3.0, 2.0, 1.0, ex, ey};
    // The centre is nearest the two poles of the shortest axis: the north
    // one.
    const Point2 centre = Project(e, Point3{});
    REQUIRE_THAT(centre.Y, WithinAbs(pi / 2, Tolerance));
    // Points on the axes beyond the surface project to the vertices.
    REQUIRE(Near(Eval(e, Project(e, Point3{5.0, 0.0, 0.0})), Point3{3.0, 0.0, 0.0}));
    REQUIRE(Near(Eval(e, Project(e, Point3{0.0, -5.0, 0.0})), Point3{0.0, -2.0, 0.0}));
    REQUIRE(Near(Eval(e, Project(e, Point3{0.0, 0.0, -5.0})), Point3{0.0, 0.0, -1.0}));
    // A point inside on the plane z = 0, close to the centre, is nearer the
    // poles than the equator: its foot is off the plane, on the north side.
    const Point3 inside{0.1, 0.05, 0.0};
    const Point2 foot = Project(e, inside);
    REQUIRE(foot.Y > 0.0);
    REQUIRE(Distance(Eval(e, foot), inside) <= SampledDistance(e, inside) + Tolerance);
    // Points on the major axis inside the evolute of the x-z section have
    // two closest points, mirror images; Project picks the north one.
    const Point3 onAxis{1.0, 0.0, 0.0};
    const Point2 axisFoot = Project(e, onAxis);
    REQUIRE(axisFoot.Y > 0.0);
    REQUIRE(Distance(Eval(e, axisFoot), onAxis) <= SampledDistance(e, onAxis) + Tolerance);
}

TEST_CASE("an ellipsoid converts exactly to a rational B-spline surface",
          "[surface][ellipsoid][conversion]")
{
    const Ellipsoid e = TiltedEllipsoid();
    const Stereon::Box2 whole = Domain(e);
    const Stereon::Box2 patch{Stereon::Interval{0.5, 2.0}, Stereon::Interval{-0.75, 1.0}};
    for (const Stereon::Box2& range : {whole, patch}) {
        const auto converted = ToBSpline(e, range);
        REQUIRE(converted.has_value());
        const Stereon::Box2 d = Domain(*converted);
        for (int i = 0; i <= 16; ++i) {
            for (int j = 0; j <= 16; ++j) {
                const Point2 st{d.X.GetLower() + Width(d.X) * i / 16.0,
                                d.Y.GetLower() + Width(d.Y) * j / 16.0};
                REQUIRE(std::abs(Implicit(e, Eval(*converted, st))) <= 1e-13);
            }
        }
        // The corners are the ellipsoid's at the corners of the range.
        REQUIRE(Near(Eval(*converted, Point2{d.X.GetLower(), d.Y.GetLower()}),
                     Eval(e, Point2{range.X.GetLower(), range.Y.GetLower()})));
        REQUIRE(Near(Eval(*converted, Point2{d.X.GetUpper(), d.Y.GetUpper()}),
                     Eval(e, Point2{range.X.GetUpper(), range.Y.GetUpper()})));
    }
}

TEST_CASE("an ellipsoid validates data from outside the kernel", "[surface][ellipsoid][factory]")
{
    REQUIRE(Stereon::MakeEllipsoid(Point3{}, 3.0, 2.0, 1.0, ex, ey).has_value());
    REQUIRE(Stereon::MakeEllipsoid(Point3{}, 1.0, 1.0, 1.0, ex, ey).has_value());
    const auto code = [](const auto& made) { return made.error().GetCode(); };
    REQUIRE(code(Stereon::MakeEllipsoid(Point3{}, 0.0, 2.0, 1.0, ex, ey)) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoid(Point3{}, 3.0, -2.0, 1.0, ex, ey)) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoid(Point3{}, 3.0, 2.0, std::nan(""), ex, ey)) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoid(Point3{std::numeric_limits<double>::infinity(), 0.0, 0.0},
                                        3.0, 2.0, 1.0, ex, ey)) == ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeEllipsoid(Point3{}, 3.0, 2.0, 1.0, ex, Vec3{1.0, 1.0, 0.0})) ==
            ErrorCode::InvalidInput);
}
