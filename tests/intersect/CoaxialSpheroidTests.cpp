// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Coaxial spheroids (ADR-0010, ADR-0029): an ellipsoid with two equal
// semi-axes is the turn of its meridian ellipse about the third, so against
// a cylinder, cone, sphere, torus or another spheroid on the same axis it
// meets in the circles turned by the points where the profiles meet, and
// touches along those where they touch. Every expected circle below is a
// closed form of the profiles, worked by hand in its comment; none comes
// from the kernel. Each configuration is checked in both operand orders,
// placed at the origin, turned onto another axis by an exact permutation of
// the frame, moved far out by an exact translation and scaled by powers of
// two, under which the closed forms move exactly.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::ContactType;
using Stereon::Cylinder;
using Stereon::Ellipsoid;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceSurfaceIntersection;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

/// Where a configuration is put: its frame, a point the local origin goes
/// to, and a scale. Every placement here is exact in doubles: the frames
/// are permutations of the axes, the shifts and scales dyadic.
struct Placement {
    const char* Name;
    Vec3 X;
    Vec3 Y;
    Vec3 Shift;
    double Scale;

    [[nodiscard]] Vec3 Z() const { return Cross(X, Y); }
    [[nodiscard]] Point3 At(const double z) const { return Point3{} + Shift + (Scale * z) * Z(); }
};

const std::array<Placement, 5> Placements{{
    {"origin", {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {}, 1.0},
    {"axis along x", {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {}, 1.0},
    {"far", {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {480000.5, -300000.25, 210000.75}, 1.0},
    {"scaled 2^-3", {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-3.5, 2.25, 1.0}, 0x1p-3},
    {"scaled 2^5", {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {64.0, 0.0, -32.0}, 0x1p5},
}};

/// The surfaces of a configuration, built in a placement: local heights
/// along the common axis, local lengths.
Surface Spheroid(const Placement& p, const double z, const double across, const double along)
{
    return Surface{
        Ellipsoid{p.At(z), p.Scale * across, p.Scale * across, p.Scale * along, p.X, p.Y, true}};
}

Surface CylinderOf(const Placement& p, const double radius)
{
    return Surface{Cylinder{p.At(0.0), p.Scale * radius, p.X, p.Y, true}};
}

Surface ConeApexAt(const Placement& p, const double z, const double angle)
{
    return Surface{Cone{p.At(z), 0.0, angle, p.X, p.Y, true}};
}

Surface SphereOf(const Placement& p, const double z, const double radius)
{
    return Surface{Sphere{p.At(z), p.Scale * radius, p.X, p.Y, true}};
}

Surface TorusOf(const Placement& p, const double major, const double minor)
{
    return Surface{Torus{p.At(0.0), p.Scale * major, p.Scale * minor, p.X, p.Y, true}};
}

/// A circle about the common axis, at local height z and radius rho, and
/// whether the surfaces touch along it.
struct Expected {
    double Rho;
    double Z;
    bool Tangent;
};

/// The result matches the closed form: every branch an exact circle about
/// the axis, one of the expected ones, with its contact and side
/// transition, certified, on both surfaces at its pcurves; each expected
/// circle covered over a whole turn; and a point at every expected local
/// height on the axis.
void CheckAgainst(const SurfaceSurfaceIntersection& result, const Surface& a, const Surface& b,
                  const Placement& p, const std::vector<Expected>& circles,
                  const std::vector<double>& axisPoints, const double tolerance)
{
    const double close = 1e-9 * std::max(1.0, p.Scale);
    std::vector<double> turn(circles.size(), 0.0);
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const Circle* c = branch.Curve3d.Target<Circle>();
        REQUIRE(c != nullptr);
        const Vec3 offset = c->Center - (Point3{} + p.Shift);
        const double z = Dot(offset, p.Z()) / p.Scale;
        CAPTURE(c->Radius / p.Scale, z);
        REQUIRE(Length(Cross(offset, p.Z())) <= close);
        REQUIRE(std::abs(Dot(Cross(c->XAxis, c->YAxis), p.Z())) >= 1.0 - 1e-12);
        const auto found = std::ranges::find_if(circles, [&](const Expected& e) {
            return std::abs(e.Rho * p.Scale - c->Radius) <= close &&
                   std::abs(e.Z - z) * p.Scale <= close;
        });
        REQUIRE(found != circles.end());
        const bool tangent = found->Tangent;
        REQUIRE(branch.Contact == (tangent ? ContactType::Tangential : ContactType::Crossing));
        REQUIRE(branch.Transition ==
                (tangent ? SideTransition::SameSide : SideTransition::OppositeSide));
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 128) <= ctx.GetResolution());
        turn[static_cast<std::size_t>(found - circles.begin())] +=
            branch.Range.GetUpper() - branch.Range.GetLower();
    }
    for (const double t : turn) {
        REQUIRE(std::abs(t - 2.0 * pi) <= 1e-9);
    }
    REQUIRE(result.Points.size() == axisPoints.size());
    for (const double z : axisPoints) {
        const bool seen = std::ranges::any_of(result.Points, [&](const auto& point) {
            return Distance(point.Position, p.At(z)) <= close;
        });
        REQUIRE(seen);
    }
    REQUIRE(result.Overlaps.empty());
    IntersectTest::CheckToleranceEvidence(result, tolerance);
}

/// Every placement, both operand orders.
template <typename Build>
void CheckEverywhere(const Build& build, const std::vector<Expected>& circles,
                     const std::vector<double>& axisPoints = {})
{
    for (const Placement& p : Placements) {
        CAPTURE(p.Name);
        const auto [a, b] = build(p);
        const double tolerance = 1e-5 * p.Scale;
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const Surface& first = swapped ? b : a;
            const Surface& second = swapped ? a : b;
            const auto result =
                Intersect(first, Box2::Everything(), second, Box2::Everything(), tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            CheckAgainst(*result, first, second, p, circles, axisPoints, tolerance);
        }
    }
}

struct Pair {
    Surface A;
    Surface B;
};

} // namespace

TEST_CASE("a cylinder about a spheroid's axis cuts two circles, touches its equator, or misses",
          "[intersect][elementary][spheroid]")
{
    // The oblate spheroid rho^2 / 9 + z^2 / 4 = 1 and the cylinder rho = r:
    // z = +-2 sqrt(1 - r^2 / 9), one double root at r = 3, none past it.
    const double z = 2.0 * std::sqrt(5.0) / 3.0;
    CheckEverywhere(
        [](const Placement& p) { return Pair{CylinderOf(p, 2.0), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {{2.0, -z, false}, {2.0, z, false}});
    CheckEverywhere(
        [](const Placement& p) { return Pair{CylinderOf(p, 3.0), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {{3.0, 0.0, true}});
    CheckEverywhere(
        [](const Placement& p) { return Pair{CylinderOf(p, 4.0), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {});
}

TEST_CASE("a sphere on a spheroid's axis cuts circles or touches along one or at the poles",
          "[intersect][elementary][spheroid]")
{
    // rho^2 = 9 - 9 w^2 / 4 on the spheroid, rho^2 = R^2 - (w - h)^2 on the
    // sphere centred h up the axis. Their difference, -5 w^2 / 4 - 2 h w +
    // h^2 + 9 - R^2, vanishes at
    // - h = 0, R = 2.5: w^2 = 2.2, two circles of radius sqrt(4.05);
    // - h = 0, R = 3: w = 0 twice, the equator, along which they touch;
    // - h = 0, R = 2: w = +-2, the poles, where rho = 0: two points;
    // - h = 1, R^2 = 9 (1 + h^2 / 5) = 10.8: a double root w = -0.8, the
    //   circle of radius sqrt(9 - 2.25 * 0.64) = sqrt(7.56), touching.
    const double w = std::sqrt(2.2);
    CheckEverywhere(
        [](const Placement& p) { return Pair{SphereOf(p, 0.0, 2.5), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {{std::sqrt(4.05), -w, false}, {std::sqrt(4.05), w, false}});
    CheckEverywhere(
        [](const Placement& p) { return Pair{SphereOf(p, 0.0, 3.0), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {{3.0, 0.0, true}});
    CheckEverywhere(
        [](const Placement& p) { return Pair{SphereOf(p, 0.0, 2.0), Spheroid(p, 0.0, 3.0, 2.0)}; },
        {}, {-2.0, 2.0});
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{SphereOf(p, 1.0, std::sqrt(10.8)), Spheroid(p, 0.0, 3.0, 2.0)};
        },
        {{std::sqrt(7.56), -0.8, true}});
}

TEST_CASE("a cone about a spheroid's axis cuts two circles or touches along one",
          "[intersect][elementary][spheroid]")
{
    // A 45 degree cone with its apex h below the centre: its upper generator
    // rho = w + h. Against rho^2 / 9 + w^2 / 4 = 1:
    // - h = 3: 13 w^2 + 24 w = 0, w = 0 (rho = 3) and w = -24 / 13
    //   (rho = 15 / 13); the lower generator, rho = -(w + 3), misses;
    // - h = sqrt(13): the generator's distance from the centre, h / sqrt(2),
    //   equals the ellipse's support sqrt((9 + 4) / 2) along its normal
    //   (1, -1) / sqrt(2): a touch at (9, -4) / sqrt(13).
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{ConeApexAt(p, -3.0, pi / 4.0), Spheroid(p, 0.0, 3.0, 2.0)};
        },
        {{15.0 / 13.0, -24.0 / 13.0, false}, {3.0, 0.0, false}});
    const double s = std::sqrt(13.0);
    CheckEverywhere(
        [&](const Placement& p) {
            return Pair{ConeApexAt(p, -s, pi / 4.0), Spheroid(p, 0.0, 3.0, 2.0)};
        },
        {{9.0 / s, -4.0 / s, true}});
}

TEST_CASE("a torus about a spheroid's axis cuts two circles or touches along its inner equator",
          "[intersect][elementary][spheroid]")
{
    // The tube's circle (rho - 3)^2 + w^2 = 1:
    // - the spheroid rho^2 / 4 + w^2 / 2.25 = 1 touches it from inside the
    //   hole at (2, 0), both curves' extreme point in rho;
    // - the spheroid rho^2 / 9 + w^2 / 2.25 = 1 passes through the tube's
    //   centre (3, 0) and leaves it at one point above and one below: on the
    //   tube's circle (3 + cos v, sin v), ((3 + cos v) / 3)^2 + (sin v /
    //   1.5)^2 - 1 is positive for v in [0, pi / 2], where 3 + cos v >= 3,
    //   and negative at pi; bisection over (pi / 2, pi) finds the crossing.
    CheckEverywhere(
        [](const Placement& p) { return Pair{TorusOf(p, 3.0, 1.0), Spheroid(p, 0.0, 2.0, 1.5)}; },
        {{2.0, 0.0, true}});
    const auto f = [](const double v) {
        return std::pow((3.0 + std::cos(v)) / 3.0, 2.0) + std::pow(std::sin(v) / 1.5, 2.0) - 1.0;
    };
    double lo = 0.5 * pi;
    double hi = pi;
    REQUIRE(f(lo) > 0.0);
    REQUIRE(f(hi) < 0.0);
    for (int k = 0; k < 200 && lo < 0.5 * (lo + hi) && 0.5 * (lo + hi) < hi; ++k) {
        const double m = 0.5 * (lo + hi);
        (f(m) > 0.0 ? lo : hi) = m;
    }
    const double v = 0.5 * (lo + hi);
    const double rho = 3.0 + std::cos(v);
    CheckEverywhere(
        [](const Placement& p) { return Pair{TorusOf(p, 3.0, 1.0), Spheroid(p, 0.0, 3.0, 1.5)}; },
        {{rho, -std::sin(v), false}, {rho, std::sin(v), false}});
}

TEST_CASE("coaxial spheroids cross, touch along a circle or at a pole, or coincide",
          "[intersect][elementary][spheroid]")
{
    // Against rho^2 / 9 + w^2 / 4 = 1, centred at the origin:
    // - rho^2 / 9 + w^2 = 1: subtracting, 3 w^2 / 4 = 0, the equator twice;
    // - rho^2 / 4 + w^2 / 9 = 1: by symmetry rho^2 = w^2 = 36 / 13;
    // - the same spheroid 2 up: w = 1, rho = 3 sqrt(3) / 2;
    // - the same spheroid 4 up: only the poles meet, at w = 2;
    // - rho^2 / a^2 + (w - 2)^2 / (4 a^2) = 1 with a^2 = 81 / 8: the
    //   difference is -2 w^2 - w - 1 / 8 = -2 (w + 1 / 4)^2, a touch at
    //   w = -1 / 4 along the circle of radius 3 sqrt(63) / 8.
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{Spheroid(p, 0.0, 3.0, 2.0), Spheroid(p, 0.0, 3.0, 1.0)};
        },
        {{3.0, 0.0, true}});
    const double m = 6.0 / std::sqrt(13.0);
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{Spheroid(p, 0.0, 3.0, 2.0), Spheroid(p, 0.0, 2.0, 3.0)};
        },
        {{m, -m, false}, {m, m, false}});
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{Spheroid(p, 0.0, 3.0, 2.0), Spheroid(p, 2.0, 3.0, 2.0)};
        },
        {{1.5 * std::sqrt(3.0), 1.0, false}});
    CheckEverywhere(
        [](const Placement& p) {
            return Pair{Spheroid(p, 0.0, 3.0, 2.0), Spheroid(p, 4.0, 3.0, 2.0)};
        },
        {}, {2.0});
    const double a = std::sqrt(81.0 / 8.0);
    CheckEverywhere(
        [&](const Placement& p) {
            return Pair{Spheroid(p, 0.0, 3.0, 2.0), Spheroid(p, 2.0, a, 2.0 * a)};
        },
        {{3.0 * std::sqrt(63.0) / 8.0, -0.25, true}});

    for (const Placement& p : Placements) {
        CAPTURE(p.Name);
        const Surface s = Spheroid(p, 0.0, 3.0, 2.0);
        const auto same = Intersect(s, Box2::Everything(), Spheroid(p, 0.0, 3.0, 2.0),
                                    Box2::Everything(), 1e-5 * p.Scale, ctx);
        INFO(IntersectTest::ErrorText(same));
        REQUIRE(same.has_value());
        REQUIRE(same->Overlaps.size() == 1);
        REQUIRE(same->Overlaps.front().FullDomain);
        REQUIRE(same->Overlaps.front().Evidence == Stereon::CoincidenceEvidence::Exact);
    }
}

TEST_CASE("a spheroid touching within the tolerance is declared tangent, past it is not",
          "[intersect][elementary][spheroid]")
{
    // The cylinder 5e-4 inside or outside the equator of rho^2 / 9 +
    // w^2 / 4 = 1. Within the pair tolerance 1e-3 the touch is declared
    // along the equator, its circle certified within the tolerance and
    // marked as declared (ADR-0050). Under 1e-4 the inner cylinder cuts two
    // circles at w = +-2 sqrt(1 - r^2 / 9) and the outer misses.
    const Placement& p = Placements[0];
    const Surface spheroid = Spheroid(p, 0.0, 3.0, 2.0);
    for (const double r : {3.0 - 5e-4, 3.0 + 5e-4}) {
        CAPTURE(r);
        const Surface cylinder = CylinderOf(p, r);
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto touch = swapped ? Intersect(spheroid, Box2::Everything(), cylinder,
                                                   Box2::Everything(), 1e-3, ctx)
                                       : Intersect(cylinder, Box2::Everything(), spheroid,
                                                   Box2::Everything(), 1e-3, ctx);
            INFO(IntersectTest::ErrorText(touch));
            REQUIRE(touch.has_value());
            REQUIRE(touch->Points.empty());
            REQUIRE_FALSE(touch->Branches.empty());
            for (const Stereon::IntersectionBranch& arc : touch->Branches) {
                REQUIRE(arc.Contact == ContactType::Tangential);
                REQUIRE(arc.IsCertified());
                REQUIRE(arc.Bound >= 2e-4);
                REQUIRE(arc.Bound <= 1e-3);
                REQUIRE(arc.BoundCause == Stereon::ToleranceCause::TolerantCoincidence);
            }
            IntersectTest::CheckToleranceEvidence(*touch, 1e-3);
        }
        const auto strict =
            Intersect(cylinder, Box2::Everything(), spheroid, Box2::Everything(), 1e-4, ctx);
        INFO(IntersectTest::ErrorText(strict));
        REQUIRE(strict.has_value());
        if (r < 3.0) {
            const double z = 2.0 * std::sqrt(1.0 - r * r / 9.0);
            CheckAgainst(*strict, cylinder, spheroid, p, {{r, -z, false}, {r, z, false}}, {}, 1e-4);
        } else {
            REQUIRE(strict->IsEmpty());
        }
    }
}

TEST_CASE("a spheroid's circles are clipped to the requested boxes",
          "[intersect][elementary][spheroid]")
{
    // The equator of rho^2 / 9 + w^2 / 4 = 1 touched by the cylinder
    // rho = 3: a box of half a turn on the cylinder keeps half the circle;
    // a box on the spheroid below its equator keeps none of it.
    const Placement& p = Placements[0];
    const Surface spheroid = Spheroid(p, 0.0, 3.0, 2.0);
    const Surface cylinder = CylinderOf(p, 3.0);
    const auto half = Intersect(cylinder, Box2{Interval{0.0, pi}, Interval::Everything()}, spheroid,
                                Box2::Everything(), 1e-5, ctx);
    INFO(IntersectTest::ErrorText(half));
    REQUIRE(half.has_value());
    double turn = 0.0;
    for (const Stereon::IntersectionBranch& arc : half->Branches) {
        REQUIRE(arc.Contact == ContactType::Tangential);
        turn += arc.Range.GetUpper() - arc.Range.GetLower();
    }
    REQUIRE(std::abs(turn - pi) <= 1e-9);
    const auto below =
        Intersect(cylinder, Box2::Everything(), spheroid,
                  Box2{Interval{0.0, 2.0 * pi}, Interval{-0.5 * pi, -0.25}}, 1e-5, ctx);
    INFO(IntersectTest::ErrorText(below));
    REQUIRE(below.has_value());
    REQUIRE(below->IsEmpty());
}
