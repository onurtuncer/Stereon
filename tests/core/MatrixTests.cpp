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
#include <cstddef>
#include <format>
#include <numbers>

import Stereon.Core;

#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BasicMat;
using Stereon::BasicPoint3;
using Stereon::BasicTransform;
using Stereon::BasicVec3;
using Stereon::MakeScaling;
using Stereon::Mat2;
using Stereon::Mat3;
using Stereon::Mat4;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec3;
using StereonTest::Dual;
using StereonTest::Variable;

namespace {

constexpr double Pi = std::numbers::pi;

rc::Gen<double> SmallDyadic()
{
    return rc::gen::map(rc::gen::inRange(-1024, 1025),
                        [](const int value) { return static_cast<double>(value) / 64.0; });
}

rc::Gen<Mat3> SmallMat3()
{
    return rc::gen::map(rc::gen::container<std::array<double, 9>>(SmallDyadic()),
                        [](const std::array<double, 9>& elements) { return Mat3{elements}; });
}

rc::Gen<Point3> SmallPoint3()
{
    return rc::gen::construct<Point3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
}

template <std::size_t N>
bool ApproxEqual(const BasicMat<double, N>& a, const BasicMat<double, N>& b, const double tolerance)
{
    for (std::size_t i = 0; i < N * N; ++i) {
        if (std::abs(a.Elements[i] - b.Elements[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

bool ApproxEqual(const Point3& a, const Point3& b, const double tolerance)
{
    return Distance(a, b) <= tolerance;
}

} // namespace

TEST_CASE("matrices are row-major and indexable", "[core][math]")
{
    constexpr Mat2 m{1, 2, 3, 4};
    STATIC_REQUIRE(m(0, 0) == 1.0);
    STATIC_REQUIRE(m(0, 1) == 2.0);
    STATIC_REQUIRE(m(1, 0) == 3.0);
    STATIC_REQUIRE(m(1, 1) == 4.0);

    STATIC_REQUIRE(Mat3::Identity() == Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
    STATIC_REQUIRE(Mat3::Diagonal({1, 2, 3}) == Mat3{1, 0, 0, 0, 2, 0, 0, 0, 3});
    STATIC_REQUIRE(Mat3{} == Mat3{0, 0, 0, 0, 0, 0, 0, 0, 0});

    STATIC_REQUIRE(Transposed(m) == Mat2{1, 3, 2, 4});
    STATIC_REQUIRE(Trace(m) == 5.0);

    Mat3 n;
    n(2, 0) = 7.0;
    REQUIRE(n.Elements[6] == 7.0);
}

TEST_CASE("matrix arithmetic", "[core][math]")
{
    constexpr Mat2 a{1, 2, 3, 4};
    constexpr Mat2 b{5, 6, 7, 8};
    STATIC_REQUIRE(a + b == Mat2{6, 8, 10, 12});
    STATIC_REQUIRE(b - a == Mat2{4, 4, 4, 4});
    STATIC_REQUIRE(-a == Mat2{-1, -2, -3, -4});
    STATIC_REQUIRE(2.0 * a == Mat2{2, 4, 6, 8});
    STATIC_REQUIRE(a * 2.0 == Mat2{2, 4, 6, 8});
    STATIC_REQUIRE(a * b == Mat2{19, 22, 43, 50});
    STATIC_REQUIRE(Mat2::Identity() * a == a);

    constexpr Mat3 m{1, 2, 3, 4, 5, 6, 7, 8, 9};
    STATIC_REQUIRE(m * Vec3{1, 0, 0} == Vec3{1, 4, 7});
    STATIC_REQUIRE(m * Vec3{0, 1, 0} == Vec3{2, 5, 8});
    STATIC_REQUIRE(a * Stereon::Vec2{1, 1} == Stereon::Vec2{3, 7});
}

TEST_CASE("determinants, cofactors and inverses", "[core][math]")
{
    STATIC_REQUIRE(Determinant(Mat2{1, 2, 3, 4}) == -2.0);
    STATIC_REQUIRE(Determinant(BasicMat<double, 1>{5}) == 5.0);

    // Unimodular, so the inverse is an integer matrix and the test is exact.
    constexpr Mat3 m{1, 2, 3, 0, 1, 4, 5, 6, 0};
    STATIC_REQUIRE(Determinant(m) == 1.0);
    STATIC_REQUIRE(Cofactor(m, 0, 0) == -24.0);
    STATIC_REQUIRE(Cofactor(m, 0, 1) == 20.0);
    STATIC_REQUIRE(Inverse(m) == Mat3{-24, 18, 5, 20, -15, -4, -5, 4, 1});
    STATIC_REQUIRE(m * Inverse(m) == Mat3::Identity());
    STATIC_REQUIRE(Inverse(m) * m == Mat3::Identity());

    constexpr Mat4 u{2, 1, 0, 0, 0, 3, 0, 0, 0, 0, 4, 1, 0, 0, 0, 5};
    STATIC_REQUIRE(Determinant(u) == 120.0);
    REQUIRE(ApproxEqual(u * Inverse(u), Mat4::Identity(), 1e-15));
    STATIC_REQUIRE(Inverse(Mat4::Identity()) == Mat4::Identity());
}

TEST_CASE("transforms compose and invert", "[core][math]")
{
    constexpr Vec3 d{1, 2, 3};
    constexpr Point3 p{1, 1, 1};

    STATIC_REQUIRE(Transform{} == Transform::Identity());
    STATIC_REQUIRE(Transformed(p, Transform::Identity()) == p);

    constexpr Transform translation = MakeTranslation(d);
    STATIC_REQUIRE(Transformed(Point3{}, translation) == Point3{1, 2, 3});
    // A displacement has no position, so translation leaves it alone.
    STATIC_REQUIRE(Transformed(Vec3{1, 1, 1}, translation) == Vec3{1, 1, 1});

    constexpr Transform scaling = MakeScaling(2.0);
    STATIC_REQUIRE(Transformed(p, scaling) == Point3{2, 2, 2});
    STATIC_REQUIRE(Transformed(Vec3{1, 1, 1}, scaling) == Vec3{2, 2, 2});
    STATIC_REQUIRE(Transformed(p, MakeScaling(Vec3{1, 2, 3})) == Point3{1, 2, 3});

    // (a * b) applies b first.
    constexpr Transform scaleThenTranslate = translation * scaling;
    STATIC_REQUIRE(Transformed(p, scaleThenTranslate) == Point3{3, 4, 5});
    STATIC_REQUIRE(Transformed(p, scaling * translation) == Point3{4, 6, 8});

    STATIC_REQUIRE(Inverse(translation) == MakeTranslation(-d));
    STATIC_REQUIRE(Transformed(Transformed(p, scaleThenTranslate), Inverse(scaleThenTranslate)) ==
                   p);
    STATIC_REQUIRE(Inverse(scaleThenTranslate) * scaleThenTranslate == Transform::Identity());
}

TEST_CASE("rotations are right-handed and length-preserving", "[core][math]")
{
    const Transform quarterTurn = MakeRotation(Vec3{0, 0, 1}, Pi / 2);
    REQUIRE(ApproxEqual(Transformed(Point3{1, 0, 0}, quarterTurn), Point3{0, 1, 0}, 1e-15));
    REQUIRE(ApproxEqual(Transformed(Point3{0, 1, 0}, quarterTurn), Point3{-1, 0, 0}, 1e-15));
    REQUIRE(ApproxEqual(Transformed(Point3{0, 0, 1}, quarterTurn), Point3{0, 0, 1}, 1e-15));

    // The axis is normalised, so its length does not matter.
    const Transform scaledAxis = MakeRotation(Vec3{0, 0, 5}, Pi / 2);
    REQUIRE(scaledAxis == quarterTurn);

    const Transform r = MakeRotation(Vec3{1, 2, 3}, 0.6);
    REQUIRE_THAT(Determinant(r.Linear), WithinAbs(1.0, 1e-15));
    REQUIRE(ApproxEqual(r.Linear * Transposed(r.Linear), Mat3::Identity(), 1e-15));
    REQUIRE(ApproxEqual(Inverse(r).Linear, Transposed(r.Linear), 1e-15));
    REQUIRE(ApproxEqual(MakeRotation(Vec3{1, 2, 3}, -0.6).Linear, Transposed(r.Linear), 1e-15));
}

TEST_CASE("ToMat4 embeds a transform in homogeneous coordinates", "[core][math]")
{
    constexpr Transform t{Mat3{1, 2, 3, 4, 5, 6, 7, 8, 9}, Vec3{10, 11, 12}};
    STATIC_REQUIRE(ToMat4(t) == Mat4{1, 2, 3, 10, 4, 5, 6, 11, 7, 8, 9, 12, 0, 0, 0, 1});
}

TEST_CASE("matrices and transforms format as nested rows", "[core][math]")
{
    REQUIRE(std::format("{}", Mat2{1, 2, 3, 4}) == "[[1, 2], [3, 4]]");
    REQUIRE(std::format("{}", MakeTranslation(Vec3{1, 2, 3})) ==
            "{Linear: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], Translation: (1, 2, 3)}");
}

TEST_CASE("transforms are scalar-generic", "[core][math][adr-0021]")
{
    // Rotating (1, 0, 0) about Z by theta gives (cos theta, sin theta, 0), whose
    // derivative in theta is (-sin theta, cos theta, 0).
    const Dual theta = Variable(0.7);
    const BasicTransform<Dual> rotation =
        MakeRotation(BasicVec3<Dual>{Dual{}, Dual{}, Dual{1}}, theta);
    const BasicPoint3<Dual> moved =
        Transformed(BasicPoint3<Dual>{Dual{1}, Dual{}, Dual{}}, rotation);

    REQUIRE_THAT(moved.X.Value, WithinAbs(std::cos(0.7), 1e-15));
    REQUIRE_THAT(moved.Y.Value, WithinAbs(std::sin(0.7), 1e-15));
    REQUIRE_THAT(moved.X.Deriv, WithinAbs(-std::sin(0.7), 1e-15));
    REQUIRE_THAT(moved.Y.Deriv, WithinAbs(std::cos(0.7), 1e-15));
    REQUIRE(moved.Z.Deriv == 0.0);

    // A rotation keeps unit determinant for every theta.
    const Dual det = Determinant(rotation.Linear);
    REQUIRE_THAT(det.Value, WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(det.Deriv, WithinAbs(0.0, 1e-15));

    // A translation by a parameter-dependent offset passes the derivative on.
    const Dual s = Variable(3.0);
    const BasicTransform<Dual> slide = MakeTranslation(BasicVec3<Dual>{s, s * s, Dual{}});
    const BasicPoint3<Dual> slid = Transformed(BasicPoint3<Dual>{}, slide);
    REQUIRE(slid.X.Deriv == 1.0);
    REQUIRE(slid.Y.Deriv == 6.0);
}

TEST_CASE("the Taylor residual of a rotation decays at second order", "[core][math][adr-0021]")
{
    // ADR-0021 verification: |x(t + h) - x(t) - h J| shrinks like h^2 over a
    // step range above round-off.
    const Vec3 axis{1, 2, 3};
    const Point3 p{1, -2, 0.5};
    constexpr double theta = 0.4;

    const Dual dualTheta = Variable(theta);
    const BasicPoint3<Dual> x =
        Transformed(BasicPoint3<Dual>{p.X, p.Y, p.Z},
                    MakeRotation(BasicVec3<Dual>{axis.X, axis.Y, axis.Z}, dualTheta));
    const Vec3 derivative{x.X.Deriv, x.Y.Deriv, x.Z.Deriv};

    const auto residual = [&](const double h) {
        const Point3 moved = Transformed(p, MakeRotation(axis, theta + h));
        const Point3 predicted = Transformed(p, MakeRotation(axis, theta)) + h * derivative;
        return Distance(moved, predicted);
    };

    const double coarse = residual(1e-2);
    const double fine = residual(1e-3);
    REQUIRE(coarse > 1e-6);
    REQUIRE_THAT(coarse / fine, WithinAbs(100.0, 5.0));
}

TEST_CASE("matrix identities hold bit for bit on dyadic inputs", "[core][math][property]")
{
    rc::prop("transposition is an involution", [] {
        const Mat3 m = *SmallMat3();
        RC_ASSERT(Transposed(Transposed(m)) == m);
    });

    rc::prop("the identity is neutral", [] {
        const Mat3 m = *SmallMat3();
        RC_ASSERT(m * Mat3::Identity() == m);
        RC_ASSERT(Mat3::Identity() * m == m);
    });

    rc::prop("transposition reverses products", [] {
        const Mat3 a = *SmallMat3();
        const Mat3 b = *SmallMat3();
        RC_ASSERT(Transposed(a * b) == Transposed(b) * Transposed(a));
    });

    rc::prop("translations compose by adding offsets", [] {
        const Point3 p = *SmallPoint3();
        const Vec3 d = *rc::gen::construct<Vec3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
        const Vec3 e = *rc::gen::construct<Vec3>(SmallDyadic(), SmallDyadic(), SmallDyadic());
        RC_ASSERT(Transformed(p, MakeTranslation(d) * MakeTranslation(e)) ==
                  Transformed(p, MakeTranslation(d + e)));
    });
}

TEST_CASE("rotation derivatives agree with central differences", "[core][math][adr-0021][property]")
{
    rc::prop("for random angles and points", [] {
        const double theta = static_cast<double>(*rc::gen::inRange(-300, 301)) / 100.0;
        const Point3 p = *SmallPoint3();
        const Vec3 axis{0.3, -0.4, 0.5};

        const Dual dualTheta = Variable(theta);
        const BasicPoint3<Dual> x =
            Transformed(BasicPoint3<Dual>{p.X, p.Y, p.Z},
                        MakeRotation(BasicVec3<Dual>{axis.X, axis.Y, axis.Z}, dualTheta));

        constexpr double h = 1e-5;
        const Vec3 centralDifference = (Transformed(p, MakeRotation(axis, theta + h)) -
                                        Transformed(p, MakeRotation(axis, theta - h))) /
                                       (2.0 * h);
        const Vec3 derivative{x.X.Deriv, x.Y.Deriv, x.Z.Deriv};
        RC_ASSERT(Length(derivative - centralDifference) < 1e-8 * (1.0 + Length(ToVec(p))));
    });
}
