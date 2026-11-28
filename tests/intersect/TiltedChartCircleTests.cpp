// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Circles about a common axis on a chart that does not turn about it: a
// sphere whose poles are off the axis, and a spheroid turning about its
// frame's X axis (ADR-0010, ADR-0029). The circle is the coaxial closed
// form's, exact; its pcurve on that chart is no line of constant v but the
// certified interpolant of OvalPcurve.hpp, cut where the circle crosses the
// chart's seam. Tangent circles, which the trace refused, are answered so.
// The expected circles are the profiles' closed forms worked in each
// comment, in both operand orders, the configuration placed at the origin,
// turned by an exact permutation, moved far out and scaled by 2^3.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
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
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

/// An exact placement: a permutation of the axes, a dyadic shift and scale.
struct Placement {
    const char* Name;
    Vec3 X;
    Vec3 Y;
    Vec3 Shift;
    double Scale;

    /// The local z direction, the common axis.
    [[nodiscard]] Vec3 Z() const { return Cross(X, Y); }
    /// The point at local coordinates.
    [[nodiscard]] Point3 At(const double x, const double y, const double z) const
    {
        return Point3{} + Shift + Scale * (x * X + y * Y + z * Z());
    }
    /// A local direction.
    [[nodiscard]] Vec3 Along(const Vec3& v) const { return v.X * X + v.Y * Y + v.Z * Z(); }
};

const std::array<Placement, 4> Placements{{
    {"origin", {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {}, 1.0},
    {"permuted", {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {}, 1.0},
    {"far", {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {480000.5, -300000.25, 210000.75}, 1.0},
    {"scaled 2^3", {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-16.0, 8.0, 4.0}, 0x1p3},
}};

/// A sphere's frame turned about local x by a rational rotation: its poles,
/// +-(0, -0.8, 0.6), lie off the common axis and off every circle below.
const Vec3 TiltedY{0.0, 0.6, 0.8};

/// A sphere centred z up the axis, its poles tilted.
Surface TiltedSphere(const Placement& p, const double z, const double r)
{
    return Surface{Sphere{p.At(0, 0, z), p.Scale * r, p.X, p.Along(TiltedY), true}};
}

/// An expected circle about the axis.
struct Expected {
    double Rho;
    double Z;
    bool Tangent;
};

/// A configuration in local coordinates.
struct Configuration {
    const char* Name;
    std::function<std::array<Surface, 2>(const Placement&)> Build;
    std::array<Box2, 2> Boxes;
    std::vector<Expected> Circles;
};

const Box2 All = Box2::Everything();
const Box2 Ruled{Interval{0.0, 2.0 * pi}, Interval{-100.0, 100.0}};

/// The configurations: profiles' closed forms.
std::vector<Configuration> Configurations()
{
    const double q = std::sqrt(0.5);
    return {
        // The cylinder rho = 1 and the unit sphere: the equator, touching.
        {"cylinder round the equator",
         [](const Placement& p) {
             return std::array{Surface{Cylinder{p.At(0, 0, 0), p.Scale, p.X, p.Y, true}},
                               TiltedSphere(p, 0.0, 1.0)};
         },
         {Ruled, All},
         {{1.0, 0.0, true}}},
        // The cylinder rho = 1 and the sphere of radius 2: rho^2 + z^2 = 4,
        // two circles at z = +-sqrt(3).
        {"cylinder through the sphere",
         [](const Placement& p) {
             return std::array{Surface{Cylinder{p.At(0, 0, 0), p.Scale, p.X, p.Y, true}},
                               TiltedSphere(p, 0.0, 2.0)};
         },
         {Ruled, All},
         {{1.0, -std::sqrt(3.0), false}, {1.0, std::sqrt(3.0), false}}},
        // A 45 degree cone from the origin and a sphere of radius sqrt(2)
        // centred 2 up: the generator rho = z is sqrt(2) from (0, 2), the
        // foot (1, 1).
        {"sphere seated in a cone",
         [](const Placement& p) {
             return std::array{Surface{Cone{p.At(0, 0, 0), 0.0, pi / 4.0, p.X, p.Y, true}},
                               TiltedSphere(p, 2.0, std::sqrt(2.0))};
         },
         {Ruled, All},
         {{1.0, 1.0, true}}},
        // The torus of radii 3 and 1 and the concentric sphere of radius 4:
        // the outer equator, touching; of radius 2, the inner one.
        {"sphere round the torus",
         [](const Placement& p) {
             return std::array{
                 Surface{Torus{p.At(0, 0, 0), p.Scale * 3.0, p.Scale, p.X, p.Y, true}},
                 TiltedSphere(p, 0.0, 4.0)};
         },
         {All, All},
         {{4.0, 0.0, true}}},
        {"sphere in the torus's hole",
         [](const Placement& p) {
             return std::array{
                 Surface{Torus{p.At(0, 0, 0), p.Scale * 3.0, p.Scale, p.X, p.Y, true}},
                 TiltedSphere(p, 0.0, 2.0)};
         },
         {All, All},
         {{2.0, 0.0, true}}},
        // The spheroid w^2 / 9 + rho^2 / 4 = 1 turns about its frame's x
        // axis, here the common z. A sphere centred h = 1 up it, of radius
        // R: rho^2 differs by 5 w^2 / 9 - 2 h w + h^2 + 4 - R^2, a double
        // root where R^2 = 4 - 4 h^2 / 5 = 3.2, at w = 9 h / 5 = 1.8,
        // rho^2 = 4 (1 - 1.8^2 / 9) = 2.56.
        {"sphere touching a spheroid inside",
         [](const Placement& p) {
             return std::array{TiltedSphere(p, 1.0, std::sqrt(3.2)),
                               Surface{Ellipsoid{p.At(0, 0, 0), 3.0 * p.Scale, 2.0 * p.Scale,
                                                 2.0 * p.Scale, p.Z(), p.X, true}}};
         },
         {All, All},
         {{1.6, 1.8, true}}},
        // The cylinder of radius sqrt(2): w^2 = 9 (1 - 2 / 4), w = +-3 q.
        {"cylinder through a spheroid",
         [](const Placement& p) {
             return std::array{
                 Surface{Cylinder{p.At(0, 0, 0), std::sqrt(2.0) * p.Scale, p.X, p.Y, true}},
                 Surface{Ellipsoid{p.At(0, 0, 0), 3.0 * p.Scale, 2.0 * p.Scale, 2.0 * p.Scale,
                                   p.Z(), p.X, true}}};
         },
         {Ruled, All},
         {{std::sqrt(2.0), -3.0 * q, false}, {std::sqrt(2.0), 3.0 * q, false}}},
    };
}

/// Every branch an exact circle about the axis, one of the expected ones,
/// with its contact, certified, on both surfaces at its pcurves; every
/// expected circle covered over a whole turn.
void CheckCircles(const Stereon::SurfaceSurfaceIntersection& result, const Surface& a,
                  const Surface& b, const Placement& p, const std::vector<Expected>& circles)
{
    REQUIRE(result.Points.empty());
    REQUIRE(result.Overlaps.empty());
    const double close = 1e-9 * p.Scale;
    std::vector<double> turn(circles.size(), 0.0);
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const Circle* c = branch.Curve3d.Target<Circle>();
        REQUIRE(c != nullptr);
        const Vec3 offset = c->Center - p.At(0, 0, 0);
        const double z = Dot(offset, p.Z()) / p.Scale;
        CAPTURE(c->Radius / p.Scale, z);
        REQUIRE(Length(Cross(offset, p.Z())) <= close);
        const auto found = std::ranges::find_if(circles, [&](const Expected& e) {
            return std::abs(e.Rho * p.Scale - c->Radius) <= close &&
                   std::abs(e.Z - z) * p.Scale <= close;
        });
        REQUIRE(found != circles.end());
        REQUIRE(branch.Contact ==
                (found->Tangent ? ContactType::Tangential : ContactType::Crossing));
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 256) <= ctx.GetResolution());
        turn[static_cast<std::size_t>(found - circles.begin())] +=
            branch.Range.GetUpper() - branch.Range.GetLower();
    }
    for (const double t : turn) {
        REQUIRE(std::abs(t - 2.0 * pi) <= 1e-9);
    }
}

/// Either a refusal, IntersectionFailed, or every branch certified on both
/// surfaces at its pcurves, in both orders.
void CheckRefusedOrCertified(const Surface& a, const Surface& b)
{
    for (const bool swapped : {false, true}) {
        CAPTURE(swapped);
        const Surface& first = swapped ? b : a;
        const Surface& second = swapped ? a : b;
        const Box2 boxFirst = first.Target<Cylinder>() != nullptr ? Ruled : All;
        const Box2 boxSecond = second.Target<Cylinder>() != nullptr ? Ruled : All;
        const auto result = Intersect(first, boxFirst, second, boxSecond, 1e-5, ctx);
        if (!result) {
            REQUIRE(result.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
            continue;
        }
        for (const Stereon::IntersectionBranch& branch : result->Branches) {
            REQUIRE(branch.IsCertified());
            REQUIRE(IntersectTest::PcurveDeviation(branch, first, second, 512) <=
                    ctx.GetResolution());
        }
    }
}

} // namespace

/// Closed forms: circles on a chart that does not turn about the axis are
/// answered exactly, tangent ones included, in both orders.
TEST_CASE("circles on a chart that does not turn about the common axis are exact",
          "[intersect][elementary][tilted-chart]")
{
    for (const Configuration& c : Configurations()) {
        for (const Placement& p : Placements) {
            CAPTURE(c.Name, p.Name);
            const auto [a, b] = c.Build(p);
            for (const bool swapped : {false, true}) {
                CAPTURE(swapped);
                const auto result =
                    swapped ? Intersect(b, c.Boxes[1], a, c.Boxes[0], 1e-5 * p.Scale, ctx)
                            : Intersect(a, c.Boxes[0], b, c.Boxes[1], 1e-5 * p.Scale, ctx);
                INFO(IntersectTest::ErrorText(result));
                REQUIRE(result.has_value());
                CheckCircles(*result, swapped ? b : a, swapped ? a : b, p, c.Circles);
                IntersectTest::CheckToleranceEvidence(*result, 1e-5 * p.Scale);
            }
        }
    }
}

/// A circle through a pole of the tilted chart has no longitude there: it
/// fails closed, never answered with a pcurve that skips the pole.
TEST_CASE("a circle through a pole of the tilted chart fails closed",
          "[intersect][elementary][tilted-chart]")
{
    // The unit sphere with poles +-(0, -0.8, 0.6): the cylinder of radius
    // 0.8 cuts it in the circles at z = +-0.6, the upper through the pole
    // (0, -0.8, 0.6), the lower through the other.
    const Placement& p = Placements[0];
    const Surface sphere = TiltedSphere(p, 0.0, 1.0);
    const Surface cylinder{Cylinder{p.At(0, 0, 0), 0.8, p.X, p.Y, true}};
    // The spheroid about its frame's x axis: its equator, the circle the
    // cylinder of radius 2 touches, holds both poles of its chart, +-2 y'.
    const Surface spheroid{Ellipsoid{p.At(0, 0, 0), 3.0, 2.0, 2.0, p.Z(), p.X, true}};
    const Surface round{Cylinder{p.At(0, 0, 0), 2.0, p.X, p.Y, true}};
    CheckRefusedOrCertified(cylinder, sphere);
    CheckRefusedOrCertified(round, spheroid);
}
