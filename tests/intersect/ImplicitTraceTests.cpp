// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The parametric-implicit trace (ADR-0010, the pairs of the analytic table
// without a closed form): every branch certified within the resolution of
// both surfaces, the number of curves where the configuration's topology
// is known, completeness against a dense-sampling oracle, and a refusal,
// never a guess, where a node makes the topology undecidable.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
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
using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Vec3 X{1.0, 0.0, 0.0};
const Vec3 Y{0.0, 1.0, 0.0};

/// A signed function that vanishes exactly on a surface.
using Implicit = std::function<double(const Point3&)>;

Implicit SphereImplicit(const Point3& c, const double r)
{
    return [=](const Point3& p) { return LengthSquared(p - c) - r * r; };
}

/// About the z axis through `c`.
Implicit TorusImplicit(const Point3& c, const double big, const double small)
{
    return [=](const Point3& p) {
        const Vec3 d = p - c;
        const double rho = std::hypot(d.X, d.Y);
        return (rho - big) * (rho - big) + d.Z * d.Z - small * small;
    };
}

/// Every point of a branch on both surfaces within the resolution, every
/// branch certified, its pcurves agreeing with its curve.
void RequireCertified(const Stereon::SurfaceSurfaceIntersection& result, const Surface& a,
                      const Surface& b)
{
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= ctx.GetResolution());
        REQUIRE(branch.Contact == Stereon::ContactType::Crossing);
        REQUIRE(branch.Transition == Stereon::SideTransition::OppositeSide);
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 200) <= ctx.GetResolution());
    }
}

/// The dense-sampling oracle: on a grid over a's parameter box, every sign
/// change of b's implicit function along a grid line is a point of the
/// intersection, bracketed and bisected; each must lie near a branch.
void RequireComplete(const Stereon::SurfaceSurfaceIntersection& result, const Surface& a,
                     const Box2& grid, const Implicit& onB, const std::size_t n = 120)
{
    std::vector<Point3> traced;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 2000)) {
            traced.push_back(Eval(branch.Curve3d, t));
        }
    }
    const double cells = static_cast<double>(n);
    const auto at = [&](const double i, const double j) {
        const double u = grid.X.GetLower() + (grid.X.GetUpper() - grid.X.GetLower()) * i / cells;
        const double v = grid.Y.GetLower() + (grid.Y.GetUpper() - grid.Y.GetLower()) * j / cells;
        return Eval(a, Point2{u, v});
    };
    std::size_t roots = 0;
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double fi = static_cast<double>(i);
            const double fj = static_cast<double>(j);
            if ((onB(at(fi, fj)) < 0.0) == (onB(at(fi, fj + 1.0)) < 0.0)) {
                continue;
            }
            double lo = fj;
            double hi = fj + 1.0;
            for (int k = 0; k < 60; ++k) {
                const double m = 0.5 * (lo + hi);
                ((onB(at(fi, m)) < 0.0) == (onB(at(fi, lo)) < 0.0) ? lo : hi) = m;
            }
            const Point3 root = at(fi, lo);
            double nearest = 1e300;
            for (const Point3& p : traced) {
                nearest = std::min(nearest, Distance(p, root));
            }
            INFO("an oracle root at " << root.X << ", " << root.Y << ", " << root.Z);
            REQUIRE(nearest <= 1e-2);
            ++roots;
        }
    }
    REQUIRE(roots > 0);
}

/// Independent circle oracle for a plane section, including total angular
/// coverage. It does not assume a fixed branch count or parameterisation.
void RequireCircleSection(const Stereon::SurfaceSurfaceIntersection& result, const Plane& plane,
                          const Sphere& sphere)
{
    REQUIRE(result.Points.empty());
    REQUIRE(result.Overlaps.empty());
    REQUIRE_FALSE(result.Branches.empty());
    REQUIRE(IntersectTest::ConnectivityOf(result).Curves == 1);
    const Vec3 normal = Cross(plane.XAxis, plane.YAxis);
    const Vec3 unit = normal / Length(normal);
    const double height = Dot(sphere.Center - plane.Origin, unit);
    const Point3 center = sphere.Center - height * unit;
    const double radius = std::sqrt(sphere.Radius * sphere.Radius - height * height);
    double sweep = 0.0;
    for (const auto& branch : result.Branches) {
        double previous = 0.0;
        double angleSum = 0.0;
        std::size_t index = 0;
        for (const double t : IntersectTest::SamplesOf(branch.Range, 128)) {
            const Point3 point = Eval(branch.Curve3d, t);
            CHECK(std::abs(Dot(point - plane.Origin, unit)) <= ctx.GetResolution());
            CHECK(std::abs(Distance(point, center) - radius) <= ctx.GetResolution());
            const Vec3 delta = point - center;
            const double angle = std::atan2(Dot(delta, plane.YAxis), Dot(delta, plane.XAxis));
            if (index++ > 0) {
                angleSum += std::remainder(angle - previous, 2.0 * pi);
            }
            previous = angle;
        }
        sweep += std::abs(angleSum);
    }
    CHECK(std::abs(sweep - 2.0 * pi) <=
          8.0 * static_cast<double>(result.Branches.size()) * ctx.GetResolution() / radius);
}

const Surface torus{Torus{Point3{}, 3.0, 1.0, X, Y, true}};
const Box2 torusChart{Interval{0.0, 2.0 * pi}, Interval{0.0, 2.0 * pi}};

} // namespace

TEST_CASE("a cylinder through a sphere's surface traces one closed curve", "[intersect][trace]")
{
    const Surface cylinder{IntersectTest::CylinderAbout(Point3{}, Vec3{0.0, 0.0, 1.0}, 2.0)};
    const Surface sphere{Sphere{Point3{1.0, 0.0, 0.0}, 2.0, X, Y, true}};
    const auto result = Intersect(cylinder, sphere, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_FALSE(result->Branches.empty());
    RequireCertified(*result, cylinder, sphere);
    REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == 1);
    RequireComplete(*result, cylinder, Box2{Interval{0.0, 2.0 * pi}, Interval{-3.0, 3.0}},
                    SphereImplicit(Point3{1.0, 0.0, 0.0}, 2.0));

    // The sphere first: the same curve, its sides exchanged.
    const auto swapped = Intersect(sphere, cylinder, ctx);
    INFO(IntersectTest::ErrorText(swapped));
    REQUIRE(swapped.has_value());
    RequireCertified(*swapped, sphere, cylinder);
    REQUIRE(IntersectTest::ConnectivityOf(*swapped).Curves == 1);
}

TEST_CASE("a thin cylinder inside a sphere leaves through both caps: two turns",
          "[intersect][trace]")
{
    const Surface cylinder{
        IntersectTest::CylinderAbout(Point3{0.2, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.0)};
    const Surface sphere{Sphere{Point3{}, 2.0, X, Y, true}};
    const auto result = Intersect(cylinder, sphere, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, cylinder, sphere);
    REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == 2);
}

TEST_CASE("planes parallel to a torus's axis cut its spiric sections", "[intersect][trace]")
{
    // At 3.5 from the axis one oval, at 2.5 one waisted curve, at 1.5 two.
    for (const auto& [d, curves] : {std::pair{3.5, 1}, std::pair{2.5, 1}, std::pair{1.5, 2}}) {
        const Surface plane{Plane{Point3{d, 0.0, 0.0}, Y, Vec3{0.0, 0.0, 1.0}}};
        const auto result = Intersect(plane, torus, ctx);
        INFO("d = " << d << ": " << IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        RequireCertified(*result, plane, torus);
        REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == static_cast<std::size_t>(curves));
        RequireComplete(*result, torus, torusChart, [d](const Point3& p) { return p.X - d; });
    }
}

TEST_CASE("an oblique plane through a torus", "[intersect][trace]")
{
    const Surface plane{
        IntersectTest::PlaneWithNormal(Point3{0.3, -0.2, 0.1}, Vec3{0.3, 0.2, 1.0})};
    const auto result = Intersect(torus, plane, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, torus, plane);
    const Vec3 n = Normalized(Vec3{0.3, 0.2, 1.0});
    RequireComplete(*result, torus, torusChart,
                    [n](const Point3& p) { return Dot(p - Point3{0.3, -0.2, 0.1}, n); });
}

TEST_CASE("a sphere against a torus, off the axis", "[intersect][trace]")
{
    const Surface sphere{Sphere{Point3{3.0, 0.5, 0.5}, 1.2, X, Y, true}};
    const auto result = Intersect(sphere, torus, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, sphere, torus);
    RequireComplete(*result, torus, torusChart, SphereImplicit(Point3{3.0, 0.5, 0.5}, 1.2));
}

TEST_CASE("two tori crossing at an angle", "[intersect][trace]")
{
    const Surface other{Torus{Point3{2.5, 0.0, 0.0}, 2.0, 0.7, X, Vec3{0.0, 0.0, 1.0}, true}};
    const auto result = Intersect(torus, other, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, torus, other);
    REQUIRE_FALSE(result->Branches.empty());
    // The oracle on the second torus, about z through its centre after a
    // quarter turn about x: its local frame is (x, z, -y).
    const auto onOther = TorusImplicit(Point3{}, 2.0, 0.7);
    RequireComplete(*result, torus, torusChart,
                    [&](const Point3& p) { return onOther(Point3{p.X - 2.5, p.Z, -p.Y}); });
}

TEST_CASE("a cylinder and a cone with skew axes, the cylinder's box bounded", "[intersect][trace]")
{
    const Surface cylinder{
        IntersectTest::CylinderAbout(Point3{0.5, 0.0, 0.0}, Vec3{1.0, 0.0, 1.0}, 1.0)};
    const Surface cone{Cone{Point3{0.0, 0.0, 3.0}, 1.5, pi / 6.0, X, Y, true}};
    const Box2 cylinderBox{Interval{0.0, 2.0 * pi}, Interval{-4.0, 6.0}};
    const auto result = Intersect(cylinder, cylinderBox, cone, Box2::Everything(), ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, cylinder, cone);
    REQUIRE_FALSE(result->Branches.empty());
}

TEST_CASE("a cone against a sphere off its axis", "[intersect][trace]")
{
    const Surface cone{Cone{Point3{}, 1.0, pi / 6.0, X, Y, true}};
    const Surface sphere{Sphere{Point3{1.0, 0.3, 2.0}, 1.5, X, Y, true}};
    const auto result = Intersect(cone, sphere, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, cone, sphere);
    REQUIRE_FALSE(result->Branches.empty());
}

TEST_CASE("a sphere far from a torus is proved separate", "[intersect][trace]")
{
    const Surface sphere{Sphere{Point3{20.0, 0.0, 0.0}, 1.0, X, Y, true}};
    const auto result = Intersect(sphere, torus, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->IsEmpty());
}

TEST_CASE("a sliver of a branch at a corner of the box is a point, not a separation",
          "[intersect][trace][regression]")
{
    // The audit of 2026-11-05: a chain whose samples lie within 2^-40 of one
    // another was dropped without a proof, an empty result claiming a
    // separation. The cylinder of radius 2 about z meets the sphere of
    // radius 2 about (1, 0, 0) where z² = 4 cos u - 1, falling with u. The
    // cylinder's box has its corner 1e-13 below that curve at u0 = pi / 6,
    // so the curve crosses the box over about 1e-13 near the corner.
    const Stereon::Cylinder cylinder{Point3{}, 2.0, X, Y, true};
    const Surface sphere{Sphere{Point3{1.0, 0.0, 0.0}, 2.0, X, Y, true}};
    const double u0 = pi / 6.0;
    const double z0 = std::sqrt(4.0 * std::cos(u0) - 1.0);
    const Box2 corner{Interval{u0, u0 + 1.0}, Interval{z0 - 1e-13, z0 + 1.0}};
    const auto result = Intersect(Surface{cylinder}, corner, sphere, Box2::Everything(), ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_FALSE(result->IsEmpty());
    REQUIRE(result->Branches.empty());
    const Point3 sliver = Eval(cylinder, Point2{u0, z0});
    for (const Stereon::IntersectionPoint& p : result->Points) {
        REQUIRE(p.Bound <= ctx.GetResolution());
        REQUIRE(Distance(p.Position, sliver) <= p.Bound + 1e-12);
        REQUIRE(Distance(Eval(cylinder, p.OnA), p.Position) <= p.Bound);
        REQUIRE(Distance(Eval(sphere, p.OnB), p.Position) <= p.Bound);
    }
}

TEST_CASE("a node of the intersection fails closed", "[intersect][trace]")
{
    // The plane at R - r from the axis cuts the lemniscate of Bernoulli,
    // whose node no box can certify.
    const Surface plane{Plane{Point3{2.0, 0.0, 0.0}, Y, Vec3{0.0, 0.0, 1.0}}};
    const auto result = Intersect(plane, torus, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
}

// --- The ellipsoid (ADR-0048 §3.4) -----------------------------------------

namespace {

/// The ellipsoid about the origin with semi-axes 3, 2 and 1 along x, y, z.
const Surface stretched{Stereon::Ellipsoid{Point3{}, 3.0, 2.0, 1.0, X, Y, true}};
const Box2 ellipsoidChart{Interval{0.0, 2.0 * pi}, Interval{-pi / 2.0, pi / 2.0}};

/// (x / a)² + (y / b)² + (z / c)² - 1 about the origin, the frame of space.
Implicit EllipsoidImplicit(const double a, const double b, const double c)
{
    return [=](const Point3& p) {
        return (p.X / a) * (p.X / a) + (p.Y / b) * (p.Y / b) + (p.Z / c) * (p.Z / c) - 1.0;
    };
}

} // namespace

TEST_CASE("a cylinder through an ellipsoid leaves through its top and its bottom",
          "[intersect][trace][ellipsoid]")
{
    // The cylinder x² + y² = 9/4 lies inside the equator (x / 3)² + (y / 2)² =
    // 1, so it cuts the ellipsoid in two closed curves, one above and one
    // below.
    const Surface cylinder{IntersectTest::CylinderAbout(Point3{}, Vec3{0.0, 0.0, 1.0}, 1.5)};
    const auto result = Intersect(cylinder, stretched, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, cylinder, stretched);
    REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == 2);
    RequireComplete(*result, cylinder, Box2{Interval{0.0, 2.0 * pi}, Interval{-2.0, 2.0}},
                    EllipsoidImplicit(3.0, 2.0, 1.0));
    // The other way round: the same curves.
    const auto swapped = Intersect(stretched, cylinder, ctx);
    REQUIRE(swapped.has_value());
    RequireCertified(*swapped, stretched, cylinder);
    REQUIRE(IntersectTest::ConnectivityOf(*swapped).Curves == 2);
}

TEST_CASE("a sphere off the centre of an ellipsoid, and a cone through it",
          "[intersect][trace][ellipsoid]")
{
    const Surface sphere{Sphere{Point3{2.5, 0.5, 0.25}, 1.0, X, Y, true}};
    const auto onSphere = Intersect(sphere, stretched, ctx);
    INFO(IntersectTest::ErrorText(onSphere));
    REQUIRE(onSphere.has_value());
    RequireCertified(*onSphere, sphere, stretched);
    REQUIRE(IntersectTest::ConnectivityOf(*onSphere).Curves == 1);
    RequireComplete(*onSphere, stretched, ellipsoidChart,
                    SphereImplicit(Point3{2.5, 0.5, 0.25}, 1.0));

    const Surface cone{Cone{Point3{0.5, 0.0, 0.0}, 0.5, pi / 6.0, X, Y, true}};
    const auto onCone = Intersect(cone, stretched, ctx);
    INFO(IntersectTest::ErrorText(onCone));
    REQUIRE(onCone.has_value());
    RequireCertified(*onCone, cone, stretched);
    REQUIRE_FALSE(onCone->Branches.empty());
}

TEST_CASE("two ellipsoids crossing, and a torus about one", "[intersect][trace][ellipsoid]")
{
    // Semi-axes (3, 2, 1) and (1, 1.5, 3) along x, y, z: each pokes out of
    // the other along its own long axis, in two closed curves.
    const Surface tall{Stereon::Ellipsoid{Point3{}, 1.0, 1.5, 3.0, X, Y, true}};
    const auto result = Intersect(stretched, tall, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    RequireCertified(*result, stretched, tall);
    REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == 2);
    RequireComplete(*result, stretched, ellipsoidChart, EllipsoidImplicit(1.0, 1.5, 3.0));

    // The tube about the circle of radius 2.4, its outer equator at 2.9: inside
    // the ellipsoid along x, outside it along y.
    const Surface ring{Torus{Point3{}, 2.4, 0.5, X, Y, true}};
    const auto onTorus = Intersect(ring, stretched, ctx);
    INFO(IntersectTest::ErrorText(onTorus));
    REQUIRE(onTorus.has_value());
    RequireCertified(*onTorus, ring, stretched);
    REQUIRE_FALSE(onTorus->Branches.empty());
    RequireComplete(*onTorus, ring, torusChart, EllipsoidImplicit(3.0, 2.0, 1.0));
}

TEST_CASE("an ellipsoid with equal semi-axes traces the sphere's curves",
          "[intersect][trace][ellipsoid]")
{
    // ADR-0048's oracle on the trace: a cylinder through a sphere and
    // through the ellipsoid that is the same sphere.
    const Surface cylinder{IntersectTest::CylinderAbout(Point3{}, Vec3{0.0, 0.0, 1.0}, 2.0)};
    const Surface sphere{Sphere{Point3{1.0, 0.0, 0.0}, 2.0, X, Y, true}};
    const Surface round{Stereon::Ellipsoid{Point3{1.0, 0.0, 0.0}, 2.0, 2.0, 2.0, X, Y, true}};
    const auto onSphere = Intersect(cylinder, sphere, ctx);
    const auto onEllipsoid = Intersect(cylinder, round, ctx);
    INFO(IntersectTest::ErrorText(onEllipsoid));
    REQUIRE(onSphere.has_value());
    REQUIRE(onEllipsoid.has_value());
    RequireCertified(*onEllipsoid, cylinder, round);
    REQUIRE(IntersectTest::ConnectivityOf(*onEllipsoid).Curves ==
            IntersectTest::ConnectivityOf(*onSphere).Curves);
    for (const Stereon::IntersectionBranch& branch : onEllipsoid->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            REQUIRE(std::abs(Distance(Eval(branch.Curve3d, t), Point3{1.0, 0.0, 0.0}) - 2.0) <=
                    ctx.GetResolution());
        }
    }
}

TEST_CASE("a sphere far from an ellipsoid is proved separate; a tangent one fails closed",
          "[intersect][trace][ellipsoid]")
{
    const Surface far{Sphere{Point3{20.0, 0.0, 0.0}, 1.0, X, Y, true}};
    const auto apart = Intersect(far, stretched, ctx);
    REQUIRE(apart.has_value());
    REQUIRE(apart->IsEmpty());
    // A sphere resting on the end of the longest semi-axis, touching it at
    // (3, 0, 0) from outside: a tangency the trace cannot resolve.
    const Surface resting{Sphere{Point3{4.0, 0.0, 0.0}, 1.0, X, Y, true}};
    const auto touch = Intersect(resting, stretched, ctx);
    REQUIRE_FALSE(touch.has_value());
    REQUIRE(touch.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
}

TEST_CASE("MR audit narrow sphere seam excursions remain complete certified circles",
          "[intersect][trace][regression][metamorphic]")
{
    // Coverage seed 13593231413355092879. The circle crosses the sphere's
    // seam twice between adjacent samples of the initial 64-step scan.
    // MR-01 translation, MR-02 exchange, MR-03 dyadic scale and MR-04 a
    // quarter-turn of the sphere's chart all preserve this circle oracle.
    const Vec3 normal{-2.5, -4.875, -3.0};
    for (const double scale : {0.5, 1.0, 2.0}) {
        for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
            const Point3 center = Point3{} + scale * Vec3{-5.0, -5.0, -5.0} + shift;
            const Plane plane = IntersectTest::PlaneWithNormal(
                center + (-5.0 / 8.0 * scale) * (normal / Length(normal)), normal);
            for (const bool turned : {false, true}) {
                const Sphere sphere{center, scale, turned ? Y : X, turned ? -X : Y};
                for (const bool exchanged : {false, true}) {
                    INFO("scale=" << scale << " shift=" << shift.X << " turned=" << turned
                                  << " exchanged=" << exchanged);
                    const Surface a = exchanged ? Surface{sphere} : Surface{plane};
                    const Surface b = exchanged ? Surface{plane} : Surface{sphere};
                    const auto result = Intersect(a, b, ctx);
                    INFO(IntersectTest::ErrorText(result));
                    REQUIRE(result.has_value());
                    RequireCertified(*result, a, b);
                    RequireCircleSection(*result, plane, sphere);
                }
            }
        }
    }
}
