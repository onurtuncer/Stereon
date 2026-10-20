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
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <type_traits>
#include <utility>

import Stereon.Curve;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicCircle;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::CirclePeriod;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Curve2dGeometry;
using Stereon::CurveGeometry;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;

/// A tilted unit-radius-free circle: centre off the origin, axes rotated.
Circle TiltedCircle()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 2.0, 3.0}, 0.7);
    return {Point3{1.0, -2.0, 0.5}, 2.5, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

/// A similarity: uniform scaling, rotation, translation.
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

/// The circle with Dual coordinates, for differentiating Eval and Deriv.
BasicCircle<Dual, 3> Lifted(const Circle& c)
{
    const auto lift = [](const Vec3& v) { return Stereon::BasicVec3<Dual>{v.X, v.Y, v.Z}; };
    return {Stereon::BasicPoint3<Dual>{c.Center.X, c.Center.Y, c.Center.Z}, c.Radius, lift(c.XAxis),
            lift(c.YAxis)};
}

Vec3 DerivPart(const Stereon::BasicPoint3<Dual>& p)
{
    return {p.X.Deriv, p.Y.Deriv, p.Z.Deriv};
}

Vec3 DerivPart(const Stereon::BasicVec3<Dual>& v)
{
    return {v.X.Deriv, v.Y.Deriv, v.Z.Deriv};
}

/// Not a curve: it has Eval and nothing else.
struct OnlyEval {};
[[maybe_unused]] Point3 Eval(const OnlyEval&, double)
{
    return {};
}

/// A curve too large for the wrapper's buffer: a line with a payload.
struct BigLine {
    Line Base;
    std::array<double, 32> Payload{};
};
Point3 Eval(const BigLine& c, const double t)
{
    return Eval(c.Base, t);
}
Vec3 Deriv(const BigLine& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const BigLine& c)
{
    return Domain(c.Base);
}
Stereon::Box3 BoundingBox(const BigLine& c)
{
    return BoundingBox(c.Base);
}
double Project(const BigLine& c, const Point3& p)
{
    return Project(c.Base, p);
}
BigLine Reversed(const BigLine& c)
{
    return {Reversed(c.Base), c.Payload};
}
BigLine Transformed(const BigLine& c, const Transform& x)
{
    return {Transformed(c.Base, x), c.Payload};
}

} // namespace

TEST_CASE("curve types satisfy their concepts and nothing else does", "[curve][concepts]")
{
    STATIC_REQUIRE(CurveGeometry<Line>);
    STATIC_REQUIRE(CurveGeometry<Circle>);
    STATIC_REQUIRE(CurveGeometry<Curve>);
    STATIC_REQUIRE(CurveGeometry<BigLine>);
    STATIC_REQUIRE(Curve2dGeometry<Line2d>);
    STATIC_REQUIRE(Curve2dGeometry<Circle2d>);
    STATIC_REQUIRE(Curve2dGeometry<Curve2d>);

    STATIC_REQUIRE_FALSE(CurveGeometry<Line2d>);
    STATIC_REQUIRE_FALSE(Curve2dGeometry<Line>);
    STATIC_REQUIRE_FALSE(CurveGeometry<OnlyEval>);
    STATIC_REQUIRE_FALSE(CurveGeometry<int>);
    // A Curve cannot be made from a planar curve or from a non-curve.
    STATIC_REQUIRE_FALSE(std::is_constructible_v<Curve, Line2d>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<Curve, OnlyEval>);
    STATIC_REQUIRE_FALSE(std::is_convertible_v<Line, Curve>); // explicit only
}

TEST_CASE("a line is evaluated, projected, reversed and moved at compile time", "[curve][line]")
{
    constexpr Line line{Point3{1.0, 2.0, 3.0}, Vec3{2.0, 0.0, 0.0}};
    STATIC_REQUIRE(Eval(line, 1.5) == Point3{4.0, 2.0, 3.0});
    STATIC_REQUIRE(Deriv(line, 0.0, 1) == Vec3{2.0, 0.0, 0.0});
    STATIC_REQUIRE(Deriv(line, 0.0, 2) == Vec3{});
    STATIC_REQUIRE(Domain(line) == Stereon::Interval::Everything());
    STATIC_REQUIRE(BoundingBox(line) == Stereon::Box3::Everything());
    STATIC_REQUIRE(Project(line, Point3{5.0, 9.0, -4.0}) == 2.0);
    STATIC_REQUIRE(Eval(Reversed(line), -1.5) == Eval(line, 1.5));
    STATIC_REQUIRE(Reversed(Reversed(line)) == line);

    constexpr Transform shift = Stereon::MakeTranslation(Vec3{0.0, 0.0, 1.0});
    STATIC_REQUIRE(Eval(Transformed(line, shift), 1.5) == Point3{4.0, 2.0, 4.0});

    constexpr Line2d planar{Point2{0.0, 1.0}, Vec2{1.0, 1.0}};
    STATIC_REQUIRE(Eval(planar, 2.0) == Point2{2.0, 3.0});
    STATIC_REQUIRE(Project(planar, Point2{2.0, 1.0}) == 1.0);
    STATIC_REQUIRE(BoundingBox(planar) == Stereon::Box2::Everything());
}

TEST_CASE("a circle is evaluated at its quarter points", "[curve][circle]")
{
    const Circle c{Point3{1.0, 1.0, 0.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Eval(c, 0.0) == Point3{3.0, 1.0, 0.0});
    REQUIRE(Near(Eval(c, std::numbers::pi / 2), Point3{1.0, 3.0, 0.0}));
    REQUIRE(Near(Eval(c, std::numbers::pi), Point3{-1.0, 1.0, 0.0}));
    REQUIRE(Domain(c) == Stereon::Interval{0.0, CirclePeriod});
    REQUIRE(Deriv(c, 0.0, 1) == Vec3{0.0, 2.0, 0.0});
    REQUIRE(Deriv(c, 0.0, 2) == Vec3{-2.0, 0.0, 0.0});
    REQUIRE(Deriv(c, 0.0, 4) == Vec3{2.0, 0.0, 0.0});

    // The centre projects to parameter 0; points project into [0, 2 pi).
    REQUIRE(Project(c, c.Center) == 0.0);
    REQUIRE_THAT(Project(c, Point3{1.0, -5.0, 7.0}), WithinAbs(1.5 * std::numbers::pi, Tolerance));

    const Circle2d planar{Point2{0.0, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, -1.0}}; // clockwise
    REQUIRE(Eval(planar, std::numbers::pi / 2).Y < 0.0);
    REQUIRE(Contains(BoundingBox(planar), Point2{0.0, -1.0}));
}

TEST_CASE("circle derivatives are the derivatives of its evaluation", "[curve][circle][property]")
{
    rc::prop("Deriv(c, t, n + 1) is d/dt of Deriv(c, t, n), and n = 1 of Eval", [] {
        const Circle c = TiltedCircle();
        const double t = *rc::gen::inRange(-6000, 6000) / 1000.0;
        const int n = *rc::gen::inRange(1, 7);
        const Dual dt{t, 1.0};
        RC_ASSERT(Near(Deriv(c, t, 1), DerivPart(Eval(Lifted(c), dt))));
        RC_ASSERT(Near(Deriv(c, t, n + 1), DerivPart(Deriv(Lifted(c), dt, n))));
    });
}

TEST_CASE("circles project, reverse, move and bound consistently", "[curve][circle][property]")
{
    rc::prop("Project inverts Eval on the domain", [] {
        const Circle c = TiltedCircle();
        const double t = *rc::gen::inRange(0, 6283) / 1000.0;
        const double back = Project(c, Eval(c, t));
        RC_ASSERT(back >= 0.0 && back < CirclePeriod);
        // 0 and 2 pi are the same point.
        const double gap = std::abs(back - t);
        RC_ASSERT(std::min(gap, CirclePeriod - gap) <= 1e-9);
    });

    rc::prop("Reversed traverses backwards and Transformed commutes with Eval", [] {
        const Circle c = TiltedCircle();
        const Transform x = SomeSimilarity();
        const double t = *rc::gen::inRange(0, 6283) / 1000.0;
        RC_ASSERT(Near(Eval(Reversed(c), CirclePeriod - t), Eval(c, t)));
        RC_ASSERT(Near(Eval(Transformed(c, x), t), Transformed(Eval(c, t), x)));
        RC_ASSERT(Contains(BoundingBox(c), Eval(c, t)));
    });

    // The box is tight: it reaches the radius on an axis the circle spans.
    const Circle flat{Point3{0.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Stereon::Box3 box = BoundingBox(flat);
    REQUIRE_THAT(box.X.GetUpper(), WithinAbs(3.0, Tolerance));
    REQUIRE(Contains(box, Point3{3.0, 0.0, 0.0}));
    REQUIRE(Contains(box, Point3{0.0, -3.0, 0.0}));
    REQUIRE(Contains(box.Z, 0.0));
}

TEST_CASE("Curve wraps any curve and answers as the curve does", "[curve][erased]")
{
    const Circle circle = TiltedCircle();
    const Curve curve(circle);
    REQUIRE_FALSE(curve.IsEmpty());
    for (int i = 0; i <= 8; ++i) {
        const double t = i * CirclePeriod / 8.0;
        REQUIRE(Eval(curve, t) == Eval(circle, t));
        REQUIRE(Deriv(curve, t, 2) == Deriv(circle, t, 2));
        REQUIRE(Project(curve, Eval(circle, t)) == Project(circle, Eval(circle, t)));
    }
    REQUIRE(Domain(curve) == Domain(circle));
    REQUIRE(BoundingBox(curve) == BoundingBox(circle));

    REQUIRE(curve.Target<Circle>() != nullptr);
    REQUIRE(*curve.Target<Circle>() == circle);
    REQUIRE(curve.Target<Line>() == nullptr);

    const Curve reversed = Reversed(curve);
    REQUIRE(*reversed.Target<Circle>() == Reversed(circle));
    const Curve moved = Transformed(curve, SomeSimilarity());
    REQUIRE(*moved.Target<Circle>() == Transformed(circle, SomeSimilarity()));

    const Curve empty;
    REQUIRE(empty.IsEmpty());
    REQUIRE(empty.Target<Circle>() == nullptr);

    const Curve2d planar(Circle2d{Point2{0.0, 0.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}});
    REQUIRE(Eval(planar, 0.0) == Point2{1.0, 0.0});
    REQUIRE(Reversed(planar).Target<Circle2d>() != nullptr);
}

TEST_CASE("Curve stores small curves in place and shares large ones", "[curve][erased]")
{
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Circle, Stereon::CurveBufferSize>);
    STATIC_REQUIRE_FALSE(Stereon::Detail::FitsInline<BigLine, Stereon::CurveBufferSize>);

    const Curve small(Line{Point3{}, Vec3{1.0, 0.0, 0.0}});
    Curve smallCopy = small;
    REQUIRE(smallCopy.Target<Line>() != small.Target<Line>()); // its own copy
    REQUIRE(*smallCopy.Target<Line>() == *small.Target<Line>());

    const Curve big(BigLine{Line{Point3{}, Vec3{0.0, 1.0, 0.0}}, {}});
    const Curve bigCopy = big;
    REQUIRE(bigCopy.Target<BigLine>() == big.Target<BigLine>()); // one shared value
    REQUIRE(Eval(bigCopy, 2.0) == Point3{0.0, 2.0, 0.0});

    // Assignment and moves replace the value; a moved-from curve is empty.
    smallCopy = big;
    REQUIRE(smallCopy.Target<BigLine>() == big.Target<BigLine>());
    Curve moved = std::move(smallCopy);
    REQUIRE(smallCopy.IsEmpty()); // NOLINT(bugprone-use-after-move)
    REQUIRE(Eval(moved, 1.0) == Point3{0.0, 1.0, 0.0});
    moved = small;
    REQUIRE(moved.Target<Line>() != nullptr);

    // Curves are semiregular, so the paged stores of ADR-0004 can hold them.
    Stereon::PagedArena<Curve, 4> pool;
    for (int i = 0; i < 6; ++i) {
        pool.PushBack(i % 2 == 0 ? small : big);
    }
    const Stereon::PagedArena<Curve, 4> snapshot = pool;
    pool.Set(5, Curve(TiltedCircle()));
    REQUIRE(snapshot[5].Target<BigLine>() != nullptr);
    REQUIRE(pool[5].Target<Circle>() != nullptr);
}
