// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <limits>

import Stereon.Robust;

#include "ExactDeterminant.hpp"
#include "RobustFormatting.hpp"

using Stereon::InCircle;
using Stereon::InSphere;
using Stereon::Orient2D;
using Stereon::Orient2DExact;
using Stereon::Orient3D;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sign;
using Stereon::Detail::CheckDomain;
using Stereon::Detail::DomainStatus;
using Stereon::Detail::HighestBitExponent;
using Stereon::Detail::LowestBitExponent;
using Stereon::Detail::ScaledByPowerOfTwo;

namespace {

Point2 Scaled(const Point2& p, const int k)
{
    return {std::ldexp(p.X, k), std::ldexp(p.Y, k)};
}

Point3 Scaled(const Point3& p, const int k)
{
    return {std::ldexp(p.X, k), std::ldexp(p.Y, k), std::ldexp(p.Z, k)};
}

Point2 SmallPoint2()
{
    return {static_cast<double>(*rc::gen::inRange<std::int64_t>(-1000, 1001)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-1000, 1001))};
}

Point3 SmallPoint3()
{
    return {static_cast<double>(*rc::gen::inRange<std::int64_t>(-1000, 1001)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-1000, 1001)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-1000, 1001))};
}

} // namespace

TEST_CASE("bit exponents are read off the representation", "[robust][domain]")
{
    STATIC_REQUIRE(LowestBitExponent(1.0) == 0);
    STATIC_REQUIRE(HighestBitExponent(1.0) == 0);
    STATIC_REQUIRE(LowestBitExponent(0.75) == -2);
    STATIC_REQUIRE(HighestBitExponent(0.75) == -1);
    STATIC_REQUIRE(LowestBitExponent(-6.0) == 1);
    STATIC_REQUIRE(HighestBitExponent(-6.0) == 2);
    STATIC_REQUIRE(LowestBitExponent(3.0 * 0x1p-600) == -600);
    STATIC_REQUIRE(HighestBitExponent(3.0 * 0x1p-600) == -599);
    constexpr double smallest = std::numeric_limits<double>::denorm_min();
    STATIC_REQUIRE(LowestBitExponent(smallest) == -1074);
    STATIC_REQUIRE(HighestBitExponent(smallest) == -1074);
    STATIC_REQUIRE(LowestBitExponent(5.0 * smallest) == -1074);
    STATIC_REQUIRE(HighestBitExponent(5.0 * smallest) == -1072);

    STATIC_REQUIRE(ScaledByPowerOfTwo(3.0, 10) == 3072.0);
    STATIC_REQUIRE(ScaledByPowerOfTwo(0x1p-1000, 1500) == 0x1p500);
    STATIC_REQUIRE(ScaledByPowerOfTwo(0x1p900, -1800) == 0x1p-900);
}

TEST_CASE("the domain check accepts, rescales or rejects", "[robust][domain]")
{
    // Model coordinates are inside without scaling.
    STATIC_REQUIRE(CheckDomain({Point2{1.0, 2.0}, Point2{500'000.0, 0.25}}, 5).Status ==
                   DomainStatus::Inside);
    STATIC_REQUIRE(CheckDomain({Point2{}, Point2{}}, 2).Status == DomainStatus::Inside);

    // The review's example: products of 2^-600 underflow, a shift fixes it.
    constexpr auto tiny = CheckDomain({Point2{0x1p-600, 0.0}, Point2{0.0, 0x1p-600}}, 2);
    STATIC_REQUIRE(tiny.Status == DomainStatus::Scalable);
    STATIC_REQUIRE(tiny.Shift > 0);
    // Huge coordinates are scaled down.
    constexpr auto huge = CheckDomain({Point3{0x1p400, 1.0 + 0x1p-40, 0.0}}, 3);
    STATIC_REQUIRE(huge.Status == DomainStatus::Scalable);
    STATIC_REQUIRE(huge.Shift < 0);

    STATIC_REQUIRE(CheckDomain({Point2{0x1p-1000, 0.0}, Point2{0x1p500, 0.0}}, 2).Status ==
                   DomainStatus::TooWide);
    REQUIRE(CheckDomain({Point2{std::numeric_limits<double>::infinity(), 0.0}}, 2).Status ==
            DomainStatus::NotFinite);
    REQUIRE(CheckDomain({Point2{std::numeric_limits<double>::quiet_NaN(), 0.0}}, 2).Status ==
            DomainStatus::NotFinite);
}

TEST_CASE("predicates are exact on tiny and huge coordinates", "[robust][domain]")
{
    // Run time, not compile time: GCC 16.2 rejects a hot-path contract met in
    // constant evaluation on the rescaling path as "not constant"
    // (docs/toolchain/gcc16-feature-probe.md).
    // The review's counterexample: the determinant is 2^-1200, positive, but
    // both terms of the plain formula underflow to zero.
    const double t = 0x1p-600;
    REQUIRE(Orient2D(Point2{t, 0.0}, Point2{0.0, t}, Point2{}) == Sign::Positive);
    REQUIRE(Orient2DExact(Point2{t, 0.0}, Point2{0.0, t}, Point2{}) == Sign::Positive);
    REQUIRE(Orient2D(Point2{0.0, t}, Point2{t, 0.0}, Point2{}) == Sign::Negative);
    REQUIRE(Orient2D(Point2{t, t}, Point2{2.0 * t, 2.0 * t}, Point2{}) == Sign::Zero);

    // The same configurations at the bottom of the subnormal range.
    const double d = std::numeric_limits<double>::denorm_min();
    REQUIRE(Orient2D(Point2{d, 0.0}, Point2{0.0, d}, Point2{}) == Sign::Positive);
    REQUIRE(Orient3D(Point3{d, 0.0, 0.0}, Point3{0.0, d, 0.0}, Point3{0.0, 0.0, d}, Point3{}) ==
            Sign::Positive);

    // A unit-circle test shrunk to 2^-500 and blown up to 2^200.
    for (const int k : {-500, -300, 0, 150, 200}) {
        const auto s = [k](const double x, const double y) { return Scaled(Point2{x, y}, k); };
        REQUIRE(InCircle(s(1, 0), s(0, 1), s(-1, 0), s(0, 0)) == Sign::Positive);
        REQUIRE(InCircle(s(1, 0), s(0, 1), s(-1, 0), s(0, -1)) == Sign::Zero);
        REQUIRE(InCircle(s(1, 0), s(0, 1), s(-1, 0), s(0, 2)) == Sign::Negative);
    }
}

#if STEREON_HAS_EXACT_ORACLE

TEST_CASE("scaling by any power of two keeps every predicate's answer",
          "[robust][domain][property]")
{
    // The determinants are homogeneous, so 2^k * (all coordinates) has the
    // sign of the unscaled input, which the integer oracle decides. k ranges
    // over every exponent that keeps the coordinates finite and exact.
    rc::prop("Orient2D and Orient3D", [] {
        const int k = *rc::gen::inRange(-1074, 1013);
        const Point2 a = SmallPoint2();
        const Point2 b = SmallPoint2();
        const Point2 c = SmallPoint2();
        RC_ASSERT(Orient2D(Scaled(a, k), Scaled(b, k), Scaled(c, k)) ==
                  Stereon::Tests::ExactOrient2D(a, b, c));
        const Point3 p = SmallPoint3();
        const Point3 q = SmallPoint3();
        const Point3 r = SmallPoint3();
        const Point3 s = SmallPoint3();
        RC_ASSERT(Orient3D(Scaled(p, k), Scaled(q, k), Scaled(r, k), Scaled(s, k)) ==
                  Stereon::Tests::ExactOrient3D(p, q, r, s));
    });

    rc::prop("InCircle and InSphere", [] {
        const int k = *rc::gen::inRange(-1074, 1013);
        const Point2 a = SmallPoint2();
        const Point2 b = SmallPoint2();
        const Point2 c = SmallPoint2();
        const Point2 d = SmallPoint2();
        RC_ASSERT(InCircle(Scaled(a, k), Scaled(b, k), Scaled(c, k), Scaled(d, k)) ==
                  Stereon::Tests::ExactInCircle(a, b, c, d));
        const Point3 p = SmallPoint3();
        const Point3 q = SmallPoint3();
        const Point3 r = SmallPoint3();
        const Point3 s = SmallPoint3();
        const Point3 e = SmallPoint3();
        RC_ASSERT(InSphere(Scaled(p, k), Scaled(q, k), Scaled(r, k), Scaled(s, k), Scaled(e, k)) ==
                  Stereon::Tests::ExactInSphere(p, q, r, s, e));
    });

    rc::prop("degenerate configurations stay degenerate at every scale", [] {
        const int k = *rc::gen::inRange(-1074, 1013);
        const Point2 a = SmallPoint2();
        const Point2 b = SmallPoint2();
        const Point2 c{2.0 * b.X - a.X, 2.0 * b.Y - a.Y}; // b is the midpoint of ac
        RC_ASSERT(Orient2D(Scaled(a, k), Scaled(b, k), Scaled(c, k)) == Sign::Zero);
    });
}

#endif // STEREON_HAS_EXACT_ORACLE
