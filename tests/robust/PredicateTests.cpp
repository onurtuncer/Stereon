// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <numeric>
#include <random>

import Stereon.Robust;

#include "ExactDeterminant.hpp"
#include "RobustFormatting.hpp"

using Stereon::InCircle;
using Stereon::InCircleExact;
using Stereon::InSphere;
using Stereon::InSphereExact;
using Stereon::Orient2D;
using Stereon::Orient2DExact;
using Stereon::Orient3D;
using Stereon::Orient3DExact;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sign;
using Stereon::SignOf;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

// --- Plain floating point: what the predicates fix --------------------------

Sign NaiveOrient2D(const Point2& a, const Point2& b, const Point2& c)
{
    return SignOf((a.X - c.X) * (b.Y - c.Y) - (a.Y - c.Y) * (b.X - c.X));
}

double NaiveDet3(const Vec3& p, const Vec3& q, const Vec3& r)
{
    return p.X * (q.Y * r.Z - q.Z * r.Y) - p.Y * (q.X * r.Z - q.Z * r.X) +
           p.Z * (q.X * r.Y - q.Y * r.X);
}

Sign NaiveOrient3D(const Point3& a, const Point3& b, const Point3& c, const Point3& d)
{
    return SignOf(NaiveDet3(a - d, b - d, c - d));
}

Sign NaiveInCircle(const Point2& a, const Point2& b, const Point2& c, const Point2& d)
{
    const Vec2 ad = a - d;
    const Vec2 bd = b - d;
    const Vec2 cd = c - d;
    const double det = (ad.X * ad.X + ad.Y * ad.Y) * (bd.X * cd.Y - bd.Y * cd.X) +
                       (bd.X * bd.X + bd.Y * bd.Y) * (cd.X * ad.Y - cd.Y * ad.X) +
                       (cd.X * cd.X + cd.Y * cd.Y) * (ad.X * bd.Y - ad.Y * bd.X);
    return SignOf(det);
}

Sign NaiveInSphere(const Point3& a, const Point3& b, const Point3& c, const Point3& d,
                   const Point3& e)
{
    const Vec3 ae = a - e;
    const Vec3 be = b - e;
    const Vec3 ce = c - e;
    const Vec3 de = d - e;
    const auto lift = [](const Vec3& v) { return v.X * v.X + v.Y * v.Y + v.Z * v.Z; };
    const double det = lift(de) * NaiveDet3(ae, be, ce) - lift(ce) * NaiveDet3(ae, be, de) +
                       lift(be) * NaiveDet3(ae, ce, de) - lift(ae) * NaiveDet3(be, ce, de);
    return SignOf(det);
}

// --- Deterministic integer cases --------------------------------------------

Point2 ScaledPoint(const Point2& p, const double factor)
{
    return {p.X * factor, p.Y * factor};
}

/// 1 for a wrong answer, 0 for a right one, as a count.
std::size_t Count(const bool wrong)
{
    return wrong ? 1 : 0;
}

/// Integer points, exactly degenerate configurations built from integer
/// parametrisations, and one-unit perturbations of them. Everything stays
/// within the magnitude limits of the oracle in ExactDeterminant.hpp.
class Cases {
public:
    explicit Cases(const std::uint64_t seed)
        : m_Engine(seed)
    {}

    std::int64_t Int(const std::int64_t lo, const std::int64_t hi)
    {
        return std::uniform_int_distribution<std::int64_t>{lo, hi}(m_Engine);
    }
    double Coordinate(const std::int64_t bound) { return static_cast<double>(Int(-bound, bound)); }
    Point2 RandomPoint2(const std::int64_t bound) { return {Coordinate(bound), Coordinate(bound)}; }
    Point3 RandomPoint3(const std::int64_t bound)
    {
        return {Coordinate(bound), Coordinate(bound), Coordinate(bound)};
    }

    /// Moves one coordinate by -1, 0 or +1.
    void Perturb(Point2& p)
    {
        const double delta = Coordinate(1);
        (Int(0, 1) == 0 ? p.X : p.Y) += delta;
    }
    void Perturb(Point3& p)
    {
        const double delta = Coordinate(1);
        const std::int64_t axis = Int(0, 2);
        (axis == 0 ? p.X : (axis == 1 ? p.Y : p.Z)) += delta;
    }

    /// a + s d near a and a + t d far along the line, so that the products
    /// of the plain formula exceed 2^53 and round while a one-unit move of c
    /// changes the determinant by only |b - a|.
    std::array<Point2, 3> Collinear()
    {
        const Point2 a = RandomPoint2(std::int64_t{1} << 30);
        Vec2 d{Coordinate(1000), Coordinate(1000)};
        if (d.X == 0.0 && d.Y == 0.0) {
            d.X = 1.0;
        }
        const double s = Coordinate(1024);
        const double t = Coordinate(std::int64_t{1} << 30);
        return {a, {a.X + s * d.X, a.Y + s * d.Y}, {a.X + t * d.X, a.Y + t * d.Y}};
    }

    /// Three points a + s u + t v in the plane through a spanned by u and v.
    std::array<Point3, 4> Coplanar()
    {
        const Point3 a = RandomPoint3(std::int64_t{1} << 28);
        const Vec3 u{Coordinate(1000), Coordinate(1000), Coordinate(1000)};
        const Vec3 v{Coordinate(1000), Coordinate(1000), Coordinate(1000)};
        std::array<Point3, 4> points{a, a, a, a};
        for (std::size_t i = 1; i < 4; ++i) {
            const double s = Coordinate(512);
            const double t = Coordinate(512);
            points[i] = {a.X + s * u.X + t * v.X, a.Y + s * u.Y + t * v.Y, a.Z + s * u.Z + t * v.Z};
        }
        return points;
    }

    /// Four integer points on one circle: (R (1 - t^2), 2 R t) / (1 + t^2) for
    /// distinct integers t, with R a multiple of every 1 + t^2.
    std::array<Point2, 4> Cocircular()
    {
        std::array<std::int64_t, 13> parameters{};
        std::iota(parameters.begin(), parameters.end(), std::int64_t{-6});
        std::shuffle(parameters.begin(), parameters.end(), m_Engine);
        std::int64_t radius = Int(1, 4);
        for (std::size_t i = 0; i < 4; ++i) {
            radius *= 1 + parameters[i] * parameters[i];
        }
        const Point2 center = RandomPoint2(std::int64_t{1} << 22);
        std::array<Point2, 4> points{};
        for (std::size_t i = 0; i < 4; ++i) {
            const std::int64_t t = parameters[i];
            const std::int64_t scale = radius / (1 + t * t);
            points[i] = {center.X + static_cast<double>(scale * (1 - t * t)),
                         center.Y + static_cast<double>(scale * 2 * t)};
        }
        return points;
    }

    /// Five integer points on one sphere: (2t, 2u, 1 - t^2 - u^2) R / (1 + t^2
    /// + u^2) for distinct integer pairs (t, u), with R a multiple of every
    /// 1 + t^2 + u^2.
    std::array<Point3, 5> Cospherical()
    {
        std::array<std::int64_t, 25> pairs{};
        std::iota(pairs.begin(), pairs.end(), std::int64_t{0});
        std::shuffle(pairs.begin(), pairs.end(), m_Engine);
        std::int64_t radius = Int(1, 2);
        for (std::size_t i = 0; i < 5; ++i) {
            const std::int64_t t = pairs[i] / 5 - 2;
            const std::int64_t u = pairs[i] % 5 - 2;
            radius *= 1 + t * t + u * u;
        }
        const Point3 center = RandomPoint3(std::int64_t{1} << 18);
        std::array<Point3, 5> points{};
        for (std::size_t i = 0; i < 5; ++i) {
            const std::int64_t t = pairs[i] / 5 - 2;
            const std::int64_t u = pairs[i] % 5 - 2;
            const std::int64_t scale = radius / (1 + t * t + u * u);
            points[i] = {center.X + static_cast<double>(scale * 2 * t),
                         center.Y + static_cast<double>(scale * 2 * u),
                         center.Z + static_cast<double>(scale * (1 - t * t - u * u))};
        }
        return points;
    }

private:
    std::mt19937_64 m_Engine;
};

#if STEREON_HAS_EXACT_ORACLE

/// One adversarial sweep: degenerate configurations, two thirds of them
/// moved by one unit in one coordinate, against the integer oracle.
struct Sweep {
    std::size_t Cases = 0;
    std::size_t Degenerate = 0;
    std::size_t NaiveWrong = 0;
    std::size_t AdaptiveWrong = 0;
    std::size_t ExactWrong = 0;

    void Record(const Sign expected, const Sign naive, const Sign adaptive, const Sign exact)
    {
        ++Cases;
        Degenerate += Count(expected == Sign::Zero);
        NaiveWrong += Count(naive != expected);
        AdaptiveWrong += Count(adaptive != expected);
        ExactWrong += Count(exact != expected);
    }
};

Sweep SweepOrient2D(Cases& cases, const std::size_t count)
{
    Sweep sweep;
    for (std::size_t i = 0; i < count; ++i) {
        auto [a, b, c] = cases.Collinear();
        if (cases.Int(0, 2) != 0) {
            cases.Perturb(c);
        }
        sweep.Record(Stereon::Tests::ExactOrient2D(a, b, c), NaiveOrient2D(a, b, c),
                     Orient2D(a, b, c), Orient2DExact(a, b, c));
    }
    return sweep;
}

Sweep SweepOrient3D(Cases& cases, const std::size_t count)
{
    Sweep sweep;
    for (std::size_t i = 0; i < count; ++i) {
        auto [a, b, c, d] = cases.Coplanar();
        if (cases.Int(0, 2) != 0) {
            cases.Perturb(d);
        }
        sweep.Record(Stereon::Tests::ExactOrient3D(a, b, c, d), NaiveOrient3D(a, b, c, d),
                     Orient3D(a, b, c, d), Orient3DExact(a, b, c, d));
    }
    return sweep;
}

Sweep SweepInCircle(Cases& cases, const std::size_t count)
{
    Sweep sweep;
    for (std::size_t i = 0; i < count; ++i) {
        auto [a, b, c, d] = cases.Cocircular();
        if (cases.Int(0, 2) != 0) {
            cases.Perturb(d);
        }
        sweep.Record(Stereon::Tests::ExactInCircle(a, b, c, d), NaiveInCircle(a, b, c, d),
                     InCircle(a, b, c, d), InCircleExact(a, b, c, d));
    }
    return sweep;
}

Sweep SweepInSphere(Cases& cases, const std::size_t count)
{
    Sweep sweep;
    for (std::size_t i = 0; i < count; ++i) {
        auto [a, b, c, d, e] = cases.Cospherical();
        if (cases.Int(0, 2) != 0) {
            cases.Perturb(e);
        }
        sweep.Record(Stereon::Tests::ExactInSphere(a, b, c, d, e), NaiveInSphere(a, b, c, d, e),
                     InSphere(a, b, c, d, e), InSphereExact(a, b, c, d, e));
    }
    return sweep;
}

void RequireExact(const Sweep& sweep, const std::size_t count)
{
    INFO("cases " << sweep.Cases << ", degenerate " << sweep.Degenerate << ", naive wrong "
                  << sweep.NaiveWrong);
    REQUIRE(sweep.Cases == count);
    REQUIRE(sweep.AdaptiveWrong == 0);
    REQUIRE(sweep.ExactWrong == 0);
    // The sweep has teeth: plain floating point fails on it.
    REQUIRE(sweep.NaiveWrong > 0);
    REQUIRE(sweep.Degenerate > 0);
    REQUIRE(sweep.Degenerate < sweep.Cases);
}

#endif // STEREON_HAS_EXACT_ORACLE

// --- RapidCheck generators ---------------------------------------------------

Point2 ArbitraryPoint2(const std::int64_t bound)
{
    return {static_cast<double>(*rc::gen::inRange<std::int64_t>(-bound, bound + 1)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-bound, bound + 1))};
}

Point3 ArbitraryPoint3(const std::int64_t bound)
{
    return {static_cast<double>(*rc::gen::inRange<std::int64_t>(-bound, bound + 1)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-bound, bound + 1)),
            static_cast<double>(*rc::gen::inRange<std::int64_t>(-bound, bound + 1))};
}

} // namespace

TEST_CASE("predicates answer the textbook cases at compile time", "[robust][predicates]")
{
    constexpr Point2 o{0.0, 0.0};
    constexpr Point2 x{1.0, 0.0};
    constexpr Point2 y{0.0, 1.0};
    STATIC_REQUIRE(Orient2D(o, x, y) == Sign::Positive); // counterclockwise
    STATIC_REQUIRE(Orient2D(o, y, x) == Sign::Negative);
    STATIC_REQUIRE(Orient2D(o, x, Point2{2.0, 0.0}) == Sign::Zero);
    STATIC_REQUIRE(Orient2D(o, o, o) == Sign::Zero);
    STATIC_REQUIRE(Orient2DExact(o, x, y) == Sign::Positive);

    constexpr Point3 o3{0.0, 0.0, 0.0};
    constexpr Point3 x3{1.0, 0.0, 0.0};
    constexpr Point3 y3{0.0, 1.0, 0.0};
    constexpr Point3 z3{0.0, 0.0, 1.0};
    STATIC_REQUIRE(Orient3D(o3, x3, y3, z3) == Sign::Negative); // d above the plane
    STATIC_REQUIRE(Orient3D(o3, x3, y3, Point3{0.0, 0.0, -1.0}) == Sign::Positive);
    STATIC_REQUIRE(Orient3D(o3, x3, y3, Point3{3.0, 5.0, 0.0}) == Sign::Zero);
    STATIC_REQUIRE(Orient3DExact(o3, x3, y3, z3) == Sign::Negative);

    // The unit circle through (1, 0), (0, 1), (-1, 0), counterclockwise.
    constexpr Point2 a{1.0, 0.0};
    constexpr Point2 b{0.0, 1.0};
    constexpr Point2 c{-1.0, 0.0};
    STATIC_REQUIRE(Orient2D(a, b, c) == Sign::Positive);
    STATIC_REQUIRE(InCircle(a, b, c, o) == Sign::Positive);                // inside
    STATIC_REQUIRE(InCircle(a, b, c, Point2{0.0, -1.0}) == Sign::Zero);    // on
    STATIC_REQUIRE(InCircle(a, b, c, Point2{2.0, 0.0}) == Sign::Negative); // outside
    STATIC_REQUIRE(InCircle(a, c, b, o) == Sign::Negative);                // clockwise flips
    STATIC_REQUIRE(InCircleExact(a, b, c, Point2{0.0, -1.0}) == Sign::Zero);

    // The unit sphere through four points, positively oriented.
    constexpr Point3 d3{-1.0, 0.0, 0.0};
    STATIC_REQUIRE(Orient3D(x3, y3, z3, d3) == Sign::Positive);
    STATIC_REQUIRE(InSphere(x3, y3, z3, d3, o3) == Sign::Positive);                 // inside
    STATIC_REQUIRE(InSphere(x3, y3, z3, d3, Point3{0.0, -1.0, 0.0}) == Sign::Zero); // on
    STATIC_REQUIRE(InSphere(x3, y3, z3, d3, Point3{3.0, 0.0, 0.0}) == Sign::Negative);
    STATIC_REQUIRE(InSphere(y3, x3, z3, d3, o3) == Sign::Negative); // orientation flips
    STATIC_REQUIRE(InSphereExact(x3, y3, z3, d3, Point3{0.0, -1.0, 0.0}) == Sign::Zero);
}

TEST_CASE("the orientation grid of Kettner et al. is decided exactly", "[robust][predicates]")
{
    // a = (0.5 + i u, 0.5 + j u) with u one unit in the last place of 0.5,
    // b = (12, 12), c = (24, 24): b and c lie on the diagonal, so the exact
    // answer is the sign of j - i. Plain floating point gets hundreds of
    // these wrong (Kettner, Mehlhorn, Pion, Schirra and Yap, 2008).
    const Point2 b{12.0, 12.0};
    const Point2 c{24.0, 24.0};
    std::size_t adaptiveWrong = 0;
    std::size_t exactWrong = 0;
    std::size_t naiveWrong = 0;
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) {
            const Point2 a{0.5 + i * 0x1p-53, 0.5 + j * 0x1p-53};
            const Sign expected = i == j ? Sign::Zero : (i < j ? Sign::Positive : Sign::Negative);
#if STEREON_HAS_EXACT_ORACLE
            REQUIRE(Stereon::Tests::ExactOrient2D(ScaledPoint(a, 0x1p53), ScaledPoint(b, 0x1p53),
                                                  ScaledPoint(c, 0x1p53)) == expected);
#endif
            adaptiveWrong += Count(Orient2D(a, b, c) != expected);
            exactWrong += Count(Orient2DExact(a, b, c) != expected);
            naiveWrong += Count(NaiveOrient2D(a, b, c) != expected);
        }
    }
    REQUIRE(adaptiveWrong == 0);
    REQUIRE(exactWrong == 0);
    REQUIRE(naiveWrong > 100);
}

TEST_CASE("adaptive and exact stages agree on nearly degenerate floating-point input",
          "[robust][predicates]")
{
    // Constructions in plain double arithmetic: the inputs are not dyadic, so
    // there is no integer oracle, but the exact stage is exact for any
    // doubles and the adaptive one must agree with it, and every predicate
    // must be antisymmetric under swapping two points.
    std::mt19937_64 engine{20261012};
    std::uniform_real_distribution<double> unit{-1.0, 1.0};
    const auto point2 = [&] { return Point2{unit(engine), unit(engine)}; };
    const auto point3 = [&] { return Point3{unit(engine), unit(engine), unit(engine)}; };
    const auto onCircle = [&](const Point2& center, const double radius) {
        const double angle = unit(engine) * std::numbers::pi;
        return Point2{center.X + radius * std::cos(angle), center.Y + radius * std::sin(angle)};
    };
    const auto onSphere = [&](const Point3& center, const double radius) {
        Vec3 v{unit(engine), unit(engine), unit(engine)};
        const double length = std::sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
        return Point3{center.X + radius * v.X / length, center.Y + radius * v.Y / length,
                      center.Z + radius * v.Z / length};
    };

    for (int i = 0; i < 2000; ++i) {
        const Point2 a = point2();
        const Point2 b = point2();
        const double t = unit(engine) * 2.0;
        const Point2 c{a.X + t * (b.X - a.X), a.Y + t * (b.Y - a.Y)}; // nearly collinear
        REQUIRE(Orient2D(a, b, c) == Orient2DExact(a, b, c));
        REQUIRE(Orient2D(b, a, c) == -Orient2D(a, b, c));

        const Point3 p = point3();
        const Point3 q = point3();
        const Point3 r = point3();
        const double s = unit(engine);
        const Point3 d{p.X + s * (q.X - p.X) + t * (r.X - p.X),
                       p.Y + s * (q.Y - p.Y) + t * (r.Y - p.Y),
                       p.Z + s * (q.Z - p.Z) + t * (r.Z - p.Z)}; // nearly coplanar
        REQUIRE(Orient3D(p, q, r, d) == Orient3DExact(p, q, r, d));
        REQUIRE(Orient3D(q, p, r, d) == -Orient3D(p, q, r, d));

        const Point2 center = point2();
        const double radius = 0.25 + std::fabs(unit(engine));
        const Point2 ca = onCircle(center, radius);
        const Point2 cb = onCircle(center, radius);
        const Point2 cc = onCircle(center, radius);
        const Point2 cd = onCircle(center, radius); // nearly cocircular
        REQUIRE(InCircle(ca, cb, cc, cd) == InCircleExact(ca, cb, cc, cd));
        REQUIRE(InCircle(cb, ca, cc, cd) == -InCircle(ca, cb, cc, cd));

        const Point3 center3 = point3();
        const Point3 sa = onSphere(center3, radius);
        const Point3 sb = onSphere(center3, radius);
        const Point3 sc = onSphere(center3, radius);
        const Point3 sd = onSphere(center3, radius);
        const Point3 se = onSphere(center3, radius); // nearly cospherical
        REQUIRE(InSphere(sa, sb, sc, sd, se) == InSphereExact(sa, sb, sc, sd, se));
        REQUIRE(InSphere(sb, sa, sc, sd, se) == -InSphere(sa, sb, sc, sd, se));
    }
}

#if STEREON_HAS_EXACT_ORACLE

TEST_CASE("predicates agree with exact integer arithmetic on random points",
          "[robust][predicates][property]")
{
    rc::prop("Orient2D", [] {
        const Point2 a = ArbitraryPoint2(std::int64_t{1} << 40);
        const Point2 b = ArbitraryPoint2(std::int64_t{1} << 40);
        const Point2 c = ArbitraryPoint2(std::int64_t{1} << 40);
        const Sign expected = Stereon::Tests::ExactOrient2D(a, b, c);
        RC_ASSERT(Orient2D(a, b, c) == expected);
        RC_ASSERT(Orient2DExact(a, b, c) == expected);
        RC_ASSERT(Orient2D(b, c, a) == expected); // cyclic
        RC_ASSERT(Orient2D(b, a, c) == -expected);
    });

    rc::prop("Orient3D", [] {
        const Point3 a = ArbitraryPoint3(std::int64_t{1} << 30);
        const Point3 b = ArbitraryPoint3(std::int64_t{1} << 30);
        const Point3 c = ArbitraryPoint3(std::int64_t{1} << 30);
        const Point3 d = ArbitraryPoint3(std::int64_t{1} << 30);
        const Sign expected = Stereon::Tests::ExactOrient3D(a, b, c, d);
        RC_ASSERT(Orient3D(a, b, c, d) == expected);
        RC_ASSERT(Orient3DExact(a, b, c, d) == expected);
        RC_ASSERT(Orient3D(b, a, c, d) == -expected);
        RC_ASSERT(Orient3D(b, c, a, d) == expected); // even permutation
    });

    rc::prop("InCircle", [] {
        const Point2 a = ArbitraryPoint2(std::int64_t{1} << 26);
        const Point2 b = ArbitraryPoint2(std::int64_t{1} << 26);
        const Point2 c = ArbitraryPoint2(std::int64_t{1} << 26);
        const Point2 d = ArbitraryPoint2(std::int64_t{1} << 26);
        const Sign expected = Stereon::Tests::ExactInCircle(a, b, c, d);
        RC_ASSERT(InCircle(a, b, c, d) == expected);
        RC_ASSERT(InCircleExact(a, b, c, d) == expected);
        RC_ASSERT(InCircle(b, a, c, d) == -expected);
    });

    rc::prop("InSphere", [] {
        const Point3 a = ArbitraryPoint3(std::int64_t{1} << 20);
        const Point3 b = ArbitraryPoint3(std::int64_t{1} << 20);
        const Point3 c = ArbitraryPoint3(std::int64_t{1} << 20);
        const Point3 d = ArbitraryPoint3(std::int64_t{1} << 20);
        const Point3 e = ArbitraryPoint3(std::int64_t{1} << 20);
        const Sign expected = Stereon::Tests::ExactInSphere(a, b, c, d, e);
        RC_ASSERT(InSphere(a, b, c, d, e) == expected);
        RC_ASSERT(InSphereExact(a, b, c, d, e) == expected);
        RC_ASSERT(InSphere(b, a, c, d, e) == -expected);
    });

    rc::prop("translating every point changes nothing", [] {
        const Point3 a = ArbitraryPoint3(std::int64_t{1} << 19);
        const Point3 b = ArbitraryPoint3(std::int64_t{1} << 19);
        const Point3 c = ArbitraryPoint3(std::int64_t{1} << 19);
        const Point3 d = ArbitraryPoint3(std::int64_t{1} << 19);
        const Point3 e = ArbitraryPoint3(std::int64_t{1} << 19);
        const Vec3 shift = ArbitraryPoint3(std::int64_t{1} << 19) - Point3{};
        RC_ASSERT(Orient3D(a + shift, b + shift, c + shift, d + shift) == Orient3D(a, b, c, d));
        RC_ASSERT(InSphere(a + shift, b + shift, c + shift, d + shift, e + shift) ==
                  InSphere(a, b, c, d, e));
    });
}

TEST_CASE("predicates survive the adversarial sweep", "[robust][predicates]")
{
    constexpr std::size_t count = 20'000;
    Cases cases{1};
    RequireExact(SweepOrient2D(cases, count), count);
    RequireExact(SweepOrient3D(cases, count), count);
    RequireExact(SweepInCircle(cases, count), count);
    RequireExact(SweepInSphere(cases, count), count);
}

// The ADR-0002 verification run: a million near-degenerate cases per
// predicate. Hidden from the default run; invoke with
// `stereon_robust_tests "[exhaustive]"`.
TEST_CASE("predicates survive a million near-degenerate cases per predicate",
          "[.][robust][predicates][exhaustive]")
{
    constexpr std::size_t count = 1'000'000;
    Cases cases{20261012};
    RequireExact(SweepOrient2D(cases, count), count);
    RequireExact(SweepOrient3D(cases, count), count);
    RequireExact(SweepInCircle(cases, count), count);
    RequireExact(SweepInSphere(cases, count), count);
}

#endif // STEREON_HAS_EXACT_ORACLE
