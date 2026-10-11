// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Dual.hpp"
#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <format>
#include <type_traits>

import Stereon.Core;

#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicPoint3;
using Stereon::BasicVec3;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Scalar;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::Dual;
using StereonTest::Variable;

namespace {

// Multiples of 1/64 up to 16 in magnitude: sums, differences and the products
// the tests below form stay exactly representable, so they can compare for
// bit identity.
rc::Gen<double> SmallDyadic()
{
    return rc::gen::map(rc::gen::inRange(-1024, 1025),
                        [](const int value) { return static_cast<double>(value) / 64.0; });
}

rc::Gen<Vec3> SmallVec3()
{
    return rc::gen::construct<Vec3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
}

rc::Gen<Point3> SmallPoint3()
{
    return rc::gen::construct<Point3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
}

// A requires-expression outside a template is ill-formed when its expression
// is, so "does not compile" checks go through a concept.
template <typename A, typename B>
concept Addable = requires(A a, B b) { a + b; };

} // namespace

TEST_CASE("Scalar accepts real-number types only", "[core][math]")
{
    STATIC_REQUIRE(Scalar<double>);
    STATIC_REQUIRE(Scalar<float>);
    STATIC_REQUIRE(Scalar<Dual>);
    STATIC_REQUIRE_FALSE(Scalar<int>);
    STATIC_REQUIRE_FALSE(Scalar<bool>);
}

TEST_CASE("vectors add, subtract and scale", "[core][math]")
{
    STATIC_REQUIRE(Vec3{1, 2, 3} + Vec3{4, 5, 6} == Vec3{5, 7, 9});
    STATIC_REQUIRE(Vec3{4, 5, 6} - Vec3{1, 2, 3} == Vec3{3, 3, 3});
    STATIC_REQUIRE(-Vec3{1, -2, 3} == Vec3{-1, 2, -3});
    STATIC_REQUIRE(2.0 * Vec3{1, 2, 3} == Vec3{2, 4, 6});
    STATIC_REQUIRE(Vec3{1, 2, 3} * 2.0 == Vec3{2, 4, 6});
    STATIC_REQUIRE(Vec3{2, 4, 6} / 2.0 == Vec3{1, 2, 3});

    STATIC_REQUIRE(Vec2{1, 2} + Vec2{3, 4} == Vec2{4, 6});
    STATIC_REQUIRE(Vec2{1, 2} * 3.0 == Vec2{3, 6});

    Vec3 v{1, 2, 3};
    v += Vec3{1, 1, 1};
    v *= 2.0;
    v -= Vec3{0, 0, 2};
    v /= 2.0;
    REQUIRE(v == Vec3{2, 3, 3});
}

TEST_CASE("dot and cross products", "[core][math]")
{
    STATIC_REQUIRE(Dot(Vec3{1, 2, 3}, Vec3{4, 5, 6}) == 32.0);
    STATIC_REQUIRE(Dot(Vec2{1, 2}, Vec2{3, 4}) == 11.0);

    constexpr Vec3 x{1, 0, 0};
    constexpr Vec3 y{0, 1, 0};
    constexpr Vec3 z{0, 0, 1};
    STATIC_REQUIRE(Cross(x, y) == z);
    STATIC_REQUIRE(Cross(y, z) == x);
    STATIC_REQUIRE(Cross(z, x) == y);

    STATIC_REQUIRE(Cross(Vec2{1, 0}, Vec2{0, 1}) == 1.0);
    STATIC_REQUIRE(Cross(Vec2{0, 1}, Vec2{1, 0}) == -1.0);
    STATIC_REQUIRE(Perpendicular(Vec2{1, 0}) == Vec2{0, 1});
}

TEST_CASE("lengths and normalisation", "[core][math]")
{
    STATIC_REQUIRE(LengthSquared(Vec3{1, 2, 2}) == 9.0);
    STATIC_REQUIRE(LengthSquared(Vec2{3, 4}) == 25.0);
    REQUIRE(Length(Vec3{3, 4, 0}) == 5.0);
    REQUIRE(Length(Vec2{3, 4}) == 5.0);
    REQUIRE(Normalized(Vec3{0, 0, 2}) == Vec3{0, 0, 1});
    REQUIRE(Normalized(Vec2{-3, 0}) == Vec2{-1, 0});
}

TEST_CASE("points and vectors are distinct types", "[core][math]")
{
    STATIC_REQUIRE(Point3{1, 2, 3} - Point3{0, 1, 1} == Vec3{1, 1, 2});
    STATIC_REQUIRE(Point3{1, 2, 3} + Vec3{1, 1, 1} == Point3{2, 3, 4});
    STATIC_REQUIRE(Vec3{1, 1, 1} + Point3{1, 2, 3} == Point3{2, 3, 4});
    STATIC_REQUIRE(Point3{1, 2, 3} - Vec3{1, 1, 1} == Point3{0, 1, 2});
    STATIC_REQUIRE(Point2{1, 2} - Point2{3, 5} == Vec2{-2, -3});

    // A position plus a position has no meaning, so it does not compile.
    STATIC_REQUIRE((Addable<Point3, Vec3>));
    STATIC_REQUIRE_FALSE((Addable<Point3, Point3>));
    STATIC_REQUIRE_FALSE((Addable<Point2, Point2>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<Vec3, Point3>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<Point3, Vec3>));
    STATIC_REQUIRE_FALSE((std::is_constructible_v<Point3, Vec3>));

    STATIC_REQUIRE(ToVec(Point3{1, 2, 3}) == Vec3{1, 2, 3});
    STATIC_REQUIRE(ToPoint(Vec3{1, 2, 3}) == Point3{1, 2, 3});
    STATIC_REQUIRE(ToVec(Point2{1, 2}) == Vec2{1, 2});
    STATIC_REQUIRE(ToPoint(Vec2{1, 2}) == Point2{1, 2});
}

TEST_CASE("distances, interpolation and midpoints", "[core][math]")
{
    STATIC_REQUIRE(DistanceSquared(Point3{0, 0, 0}, Point3{3, 4, 0}) == 25.0);
    REQUIRE(Distance(Point3{0, 0, 0}, Point3{3, 4, 0}) == 5.0);
    REQUIRE(Distance(Point2{1, 1}, Point2{4, 5}) == 5.0);

    STATIC_REQUIRE(Lerp(Point3{0, 0, 0}, Point3{4, 8, 12}, 0.25) == Point3{1, 2, 3});
    STATIC_REQUIRE(Lerp(Vec3{0, 0, 0}, Vec3{4, 8, 12}, 0.5) == Vec3{2, 4, 6});
    STATIC_REQUIRE(Lerp(Vec2{0, 0}, Vec2{4, 8}, 0.5) == Vec2{2, 4});
    STATIC_REQUIRE(Midpoint(Point3{0, 0, 0}, Point3{2, 4, 6}) == Point3{1, 2, 3});
    STATIC_REQUIRE(Midpoint(Point2{0, 0}, Point2{2, 4}) == Point2{1, 2});
}

TEST_CASE("components are indexable", "[core][math]")
{
    constexpr Vec3 v{1, 2, 3};
    STATIC_REQUIRE(v[0] == 1.0);
    STATIC_REQUIRE(v[1] == 2.0);
    STATIC_REQUIRE(v[2] == 3.0);

    Point3 p{1, 2, 3};
    p[1] = 7.0;
    REQUIRE(p == Point3{1, 7, 3});

    Vec2 w{1, 2};
    w[0] = 5.0;
    REQUIRE(w == Vec2{5, 2});
    REQUIRE(Point2{3, 4}[1] == 4.0);
}

TEST_CASE("vectors and points format as tuples", "[core][math]")
{
    REQUIRE(std::format("{}", Vec3{1, 2.5, -3}) == "(1, 2.5, -3)");
    REQUIRE(std::format("{}", Vec2{0.5, 0}) == "(0.5, 0)");
    REQUIRE(std::format("{}", Point3{1, 2, 3}) == "(1, 2, 3)");
    REQUIRE(std::format("{}", Point2{-1, 2}) == "(-1, 2)");

    STATIC_REQUIRE(std::formattable<Vec3, char>);
    // A scalar that cannot be formatted gives no formatter for its vectors.
    STATIC_REQUIRE_FALSE(std::formattable<BasicVec3<Dual>, char>);
}

TEST_CASE("vector algebra is scalar-generic", "[core][math][adr-0021]")
{
    // v(p) = (p, 2p, 2p) has |v| = 3p.
    const Dual p = Variable(2.0);
    const BasicVec3<Dual> v{p, Dual{2} * p, Dual{2} * p};

    const Dual length = Length(v);
    REQUIRE(length.Value == 6.0);
    REQUIRE(length.Deriv == 3.0);

    // d/dp (9 p^2) = 18 p.
    const Dual squared = LengthSquared(v);
    REQUIRE(squared.Value == 36.0);
    REQUIRE(squared.Deriv == 36.0);

    // The direction of v does not change with p.
    const BasicVec3<Dual> unit = Normalized(v);
    REQUIRE_THAT(unit.X.Value, WithinAbs(1.0 / 3.0, 1e-15));
    REQUIRE_THAT(unit.X.Deriv, WithinAbs(0.0, 1e-15));

    // (p, 0, 0) x (0, p, 0) = (0, 0, p^2), with derivative 2p.
    const BasicVec3<Dual> a{p, Dual{}, Dual{}};
    const BasicVec3<Dual> b{Dual{}, p, Dual{}};
    const BasicVec3<Dual> cross = Cross(a, b);
    REQUIRE(cross.Z.Value == 4.0);
    REQUIRE(cross.Z.Deriv == 4.0);

    // Points move with their coordinates.
    const BasicPoint3<Dual> q{p, Dual{1}, Dual{}};
    const BasicVec3<Dual> fromOrigin = q - BasicPoint3<Dual>{};
    REQUIRE(fromOrigin.X.Deriv == 1.0);
    REQUIRE(fromOrigin.Y.Deriv == 0.0);
}

TEST_CASE("vector identities hold bit for bit on dyadic inputs", "[core][math][property]")
{
    rc::prop("addition commutes", [] {
        const Vec3 a = *SmallVec3();
        const Vec3 b = *SmallVec3();
        RC_ASSERT(a + b == b + a);
    });

    rc::prop("the dot product commutes", [] {
        const Vec3 a = *SmallVec3();
        const Vec3 b = *SmallVec3();
        RC_ASSERT(Dot(a, b) == Dot(b, a));
    });

    rc::prop("the cross product is antisymmetric and orthogonal to its factors", [] {
        const Vec3 a = *SmallVec3();
        const Vec3 b = *SmallVec3();
        RC_ASSERT(Cross(a, b) == -Cross(b, a));
        RC_ASSERT(Dot(Cross(a, b), a) == 0.0);
        RC_ASSERT(Dot(Cross(a, b), b) == 0.0);
    });

    rc::prop("a displacement applied to a point is recovered by subtraction", [] {
        const Point3 p = *SmallPoint3();
        const Vec3 v = *SmallVec3();
        RC_ASSERT((p + v) - p == v);
        RC_ASSERT((p + v) - v == p);
    });

    rc::prop("interpolation reaches its endpoints", [] {
        const Point3 p = *SmallPoint3();
        const Point3 q = *SmallPoint3();
        RC_ASSERT(Lerp(p, q, 0.0) == p);
        RC_ASSERT(Lerp(p, q, 1.0) == q);
    });
}

TEST_CASE("dual derivatives agree with central differences", "[core][math][adr-0021][property]")
{
    // f(p) = |(p, p^2, 1)|, differentiated two ways.
    const auto length = [](const double p) { return Length(Vec3{p, p * p, 1.0}); };

    rc::prop("the derivative of a length matches a finite difference", [&] {
        const double p = static_cast<double>(*rc::gen::inRange(-200, 201)) / 100.0;
        const Dual dp = Variable(p);
        const Dual f = Length(BasicVec3<Dual>{dp, dp * dp, Dual{1}});

        constexpr double h = 1e-6;
        const double centralDifference = (length(p + h) - length(p - h)) / (2.0 * h);
        RC_ASSERT(f.Value == length(p));
        RC_ASSERT(std::abs(f.Deriv - centralDifference) < 1e-7);
    });
}
