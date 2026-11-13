// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A B-spline surface against every analytic surface (ADR-0010): the
// parametric-implicit trace with the B-spline as the carrier. Each answer is
// held to a closed form written here, never derived from the kernel's path:
// the curves a plane z = h cuts from a plane, a cylinder, a cone, a sphere,
// a torus and an ellipsoid, and a curved and a rational B-spline against a
// plane. The same answers under exchange of the surfaces, a translation,
// scaling by powers of two and reversed parameters (ADR-0053, MR-01 to
// MR-04). A tangency fails closed; a separation is a proved empty result.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Cone;
using Stereon::Context;
using Stereon::Cylinder;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double Height = 0.5;

/// A translation and a power-of-two scale of space, and the context scaled
/// with it.
struct Motion {
    double Scale = 1.0;
    Vec3 Shift{0.0, 0.0, 0.0};
    [[nodiscard]] Point3 Of(const Point3& p) const
    {
        return {Scale * p.X + Shift.X, Scale * p.Y + Shift.Y, Scale * p.Z + Shift.Z};
    }
    [[nodiscard]] Point3 Back(const Point3& p) const
    {
        return {(p.X - Shift.X) / Scale, (p.Y - Shift.Y) / Scale, (p.Z - Shift.Z) / Scale};
    }
    [[nodiscard]] Context Policy() const
    {
        const Context base = Context::Default();
        const double resolution = base.GetResolution() * Scale;
        const double ceiling = base.GetToleranceCeiling() * Scale;
        return Scale >= 1.0 ? base.WithToleranceCeiling(ceiling).WithResolution(resolution)
                            : base.WithResolution(resolution).WithToleranceCeiling(ceiling);
    }
};

/// The bilinear B-spline that is exactly z = Height over [-6, 6]^2, its u
/// reversed when asked.
BSplineSurface Flat(const Motion& m, const bool reversed = false)
{
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    std::vector<Point3> poles;
    for (const double x : reversed ? std::array{6.0, -6.0} : std::array{-6.0, 6.0}) {
        for (const double y : {-6.0, 6.0}) {
            poles.push_back(m.Of(Point3{x, y, Height}));
        }
    }
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// The analytic partner: the surface, how many closed curves or lines the
/// plane z = Height cuts from it, points on them, and how far a point is
/// from the partner in closed form.
struct Partner {
    std::string Name;
    std::function<Surface(const Motion&)> Make;
    std::size_t Curves;
    std::vector<Point3> Witnesses;
    std::function<double(const Point3&)> Off;
};

std::vector<Partner> Partners()
{
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    const Vec3 ez{0.0, 0.0, 1.0};
    const double h = Height;
    std::vector<Partner> out;
    out.push_back({"plane x = 1/2",
                   [=](const Motion& m) { return Surface{Plane{m.Of({0.5, 0.0, 0.0}), ey, ez}}; },
                   1,
                   {Point3{0.5, 0.0, h}, Point3{0.5, 5.0, h}},
                   [](const Point3& q) { return std::abs(q.X - 0.5); }});
    out.push_back({"cylinder of radius 2",
                   [=](const Motion& m) {
                       return Surface{Cylinder{m.Of({0.0, 0.0, 0.0}), 2.0 * m.Scale, ex, ey, true}};
                   },
                   1,
                   {Point3{2.0, 0.0, h}, Point3{0.0, -2.0, h}},
                   [](const Point3& q) { return std::abs(std::hypot(q.X, q.Y) - 2.0); }});
    {
        const double a = 0.5;
        const double rho = 1.0 + h * std::tan(a);
        out.push_back({"cone",
                       [=](const Motion& m) {
                           return Surface{Cone{m.Of({0.0, 0.0, 0.0}), m.Scale, a, ex, ey, true}};
                       },
                       1,
                       {Point3{rho, 0.0, h}, Point3{0.0, rho, h}},
                       [=](const Point3& q) {
                           return std::abs(std::hypot(q.X, q.Y) - (1.0 + q.Z * std::tan(a))) *
                                  std::cos(a);
                       }});
    }
    {
        const double rho = std::sqrt(4.0 - h * h);
        out.push_back({"sphere of radius 2",
                       [=](const Motion& m) {
                           return Surface{
                               Sphere{m.Of({0.0, 0.0, 0.0}), 2.0 * m.Scale, ex, ey, true}};
                       },
                       1,
                       {Point3{rho, 0.0, h}, Point3{-rho, 0.0, h}},
                       [](const Point3& q) { return std::abs(std::hypot(q.X, q.Y, q.Z) - 2.0); }});
    }
    {
        const double w = std::sqrt(1.0 - h * h);
        out.push_back({"torus 3, 1",
                       [=](const Motion& m) {
                           return Surface{
                               Torus{m.Of({0.0, 0.0, 0.0}), 3.0 * m.Scale, m.Scale, ex, ey, true}};
                       },
                       2,
                       {Point3{3.0 - w, 0.0, h}, Point3{0.0, 3.0 + w, h}},
                       [](const Point3& q) {
                           return std::abs(std::hypot(std::hypot(q.X, q.Y) - 3.0, q.Z) - 1.0);
                       }});
    }
    {
        const double k = std::sqrt(1.0 - (h / 1.5) * (h / 1.5));
        out.push_back({"ellipsoid 3, 2, 3/2",
                       [=](const Motion& m) {
                           return Surface{Ellipsoid{m.Of({0.0, 0.0, 0.0}), 3.0 * m.Scale,
                                                    2.0 * m.Scale, 1.5 * m.Scale, ex, ey, true}};
                       },
                       1,
                       {Point3{3.0 * k, 0.0, h}, Point3{0.0, 2.0 * k, h}},
                       [](const Point3& q) {
                           const double f =
                               q.X * q.X / 9.0 + q.Y * q.Y / 4.0 + q.Z * q.Z / 2.25 - 1.0;
                           return std::abs(f) * 0.75; // f's gradient is at least 4/3 near it
                       }});
    }
    return out;
}

/// How far a point is from the result's branches, by their sampled polylines.
double DistanceToBranches(const Stereon::SurfaceSurfaceIntersection& r, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 1024)) {
            nearest = std::min(nearest, Distance(p, Eval(branch.Curve3d, t)));
        }
    }
    return nearest;
}

/// The closed form, moved: the curves' count, the witnesses on them, and
/// every sample of every certified branch on z = Height and on the partner.
void Check(const Partner& partner, const Motion& m, const bool swapped, const bool reversed)
{
    INFO("B-spline x " << partner.Name << " scale " << m.Scale << " shift (" << m.Shift.X << ", "
                       << m.Shift.Y << ", " << m.Shift.Z << ")" << (swapped ? " swapped" : "")
                       << (reversed ? " reversed" : ""));
    const Surface spline{Flat(m, reversed)};
    const Surface other = partner.Make(m);
    const Context policy = m.Policy();
    const double tolerance = policy.GetResolution();
    const Box2 all = Box2::Everything();
    const auto r = swapped ? Intersect(other, all, spline, all, tolerance, policy)
                           : Intersect(spline, all, other, all, tolerance, policy);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Overlaps.empty());
    REQUIRE(r->Points.empty());
    REQUIRE(IntersectTest::ConnectivityOf(*r, 1e-6 * m.Scale).Curves == partner.Curves);
    for (const Point3& w : partner.Witnesses) {
        REQUIRE(DistanceToBranches(*r, m.Of(w)) <= 1e-2 * m.Scale);
    }
    for (const Stereon::IntersectionBranch& branch : r->Branches) {
        REQUIRE(branch.IsCertified());
        const double reach = 4.0 * (branch.Bound + tolerance);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            const Point3 local = m.Back(Eval(branch.Curve3d, t));
            REQUIRE(std::abs(local.Z - Height) * m.Scale <= reach);
            REQUIRE(partner.Off(local) * m.Scale <= reach);
        }
    }
}

} // namespace

TEST_CASE("a B-spline plane cuts every analytic surface as its closed form says",
          "[intersect][spline-analytic]")
{
    // Before this handler every pair returned NotImplemented.
    for (const Partner& p : Partners()) {
        Check(p, {}, false, false);
    }
}

TEST_CASE("B-spline against analytic answers move with the surfaces",
          "[intersect][spline-analytic][metamorphic]")
{
    for (const Partner& p : Partners()) {
        SECTION("MR-02: either order, " + p.Name)
        {
            Check(p, {}, true, false);
        }
        SECTION("MR-01: a translation exact in binary, " + p.Name)
        {
            Check(p, {1.0, Vec3{16.0, -32.0, 8.0}}, false, false);
            Check(p, {1.0, Vec3{16.0, -32.0, 8.0}}, true, false);
        }
        SECTION("MR-03: scaling by powers of two, " + p.Name)
        {
            for (const int k : {-4, 6}) {
                Check(p, {std::ldexp(1.0, k), Vec3{0.0, 0.0, 0.0}}, false, false);
            }
        }
        SECTION("MR-04: the B-spline's u reversed, " + p.Name)
        {
            Check(p, {}, false, true);
            Check(p, {}, true, true);
        }
    }
}

namespace {

/// z = (x - 1/2)^2 + (y - 1/2)^2 + lift over the unit square, biquadratic.
BSplineSurface Paraboloid(const double lift)
{
    const std::array<double, 3> a{0.25, -0.25, 0.25};
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({static_cast<double>(i) / 2.0, static_cast<double>(j) / 2.0,
                             lift + a.at(i) + a.at(j)});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// A unit quarter cylinder about the y axis, exact: rational quadratic in u.
BSplineSurface QuarterCylinder()
{
    const double s = std::sqrt(0.5);
    const std::vector<Point3> poles{{1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0},
                                    {1.0, 1.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}};
    const std::vector<double> weights{1.0, 1.0, s, s, 1.0, 1.0};
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 1, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};

} // namespace

TEST_CASE("curved and rational B-splines against analytic surfaces", "[intersect][spline-analytic]")
{
    const Context ctx = Context::Default();
    const Box2 all = Box2::Everything();
    SECTION("a paraboloid cut by the plane z = 0.1: the circle of radius sqrt(0.1)")
    {
        const Surface bowl{Paraboloid(0.0)};
        const Surface plane{Plane{Point3{0.0, 0.0, 0.1}, ex, ey}};
        for (const bool swapped : {false, true}) {
            const auto r = swapped ? Intersect(plane, all, bowl, all, 1e-6, ctx)
                                   : Intersect(bowl, all, plane, all, 1e-6, ctx);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            REQUIRE(IntersectTest::ConnectivityOf(*r, 1e-6).Curves == 1);
            for (const Stereon::IntersectionBranch& branch : r->Branches) {
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Point3 p = Eval(branch.Curve3d, t);
                    REQUIRE(std::abs(std::hypot(p.X - 0.5, p.Y - 0.5) - std::sqrt(0.1)) <= 1e-4);
                    REQUIRE(std::abs(p.Z - 0.1) <= 1e-4);
                }
            }
        }
    }
    SECTION("a rational quarter cylinder cut by the plane z = 1/2: the line x = sqrt(3)/2")
    {
        const Surface quarter{QuarterCylinder()};
        const Surface plane{Plane{Point3{0.0, 0.0, 0.5}, ex, ey}};
        for (const bool swapped : {false, true}) {
            const auto r = swapped ? Intersect(plane, all, quarter, all, 1e-6, ctx)
                                   : Intersect(quarter, all, plane, all, 1e-6, ctx);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            REQUIRE(r->Branches.size() == 1);
            REQUIRE(r->Branches[0].IsCertified());
            for (const double t : IntersectTest::SamplesOf(r->Branches[0].Range, 40)) {
                const Point3 p = Eval(r->Branches[0].Curve3d, t);
                REQUIRE(std::abs(p.X - std::sqrt(0.75)) <= 1e-4);
                REQUIRE(std::abs(p.Z - 0.5) <= 1e-4);
            }
        }
    }
    SECTION("a sphere apart from the B-spline plane: a proved empty result")
    {
        const Surface flat{Flat({})};
        const Surface ball{Sphere{Point3{0.0, 0.0, 10.0}, 2.0, ex, ey, true}};
        const auto r = Intersect(flat, all, ball, all, 1e-6, ctx);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->IsEmpty());
    }
    SECTION("tangencies fail closed or answer only the point of contact")
    {
        // The paraboloid's lowest point on the plane z = 0, and a sphere of
        // radius 2 resting on the plane z = Height.
        const Surface bowl{Paraboloid(0.0)};
        const Surface plane{Plane{Point3{0.0, 0.0, 0.0}, ex, ey}};
        const Surface flat{Flat({})};
        const Surface ball{Sphere{Point3{0.0, 0.0, Height + 2.0}, 2.0, ex, ey, true}};
        for (const auto& r : {Intersect(bowl, all, plane, all, 1e-6, ctx),
                              Intersect(flat, all, ball, all, 1e-6, ctx)}) {
            if (!r) {
                REQUIRE(r.error().GetCode() == ErrorCode::IntersectionFailed);
                continue;
            }
            REQUIRE(r->Branches.empty());
            REQUIRE(r->Points.size() == 1);
        }
    }
    SECTION("a requested box off the B-spline's domain leaves nothing to intersect")
    {
        const Surface flat{Flat({})};
        const Surface ball{Sphere{Point3{0.0, 0.0, 0.0}, 2.0, ex, ey, true}};
        const Box2 outside{Stereon::Interval{2.0, 3.0}, Stereon::Interval{0.0, 1.0}};
        const auto r = Intersect(flat, outside, ball, all, 1e-6, ctx);
        REQUIRE(r.has_value());
        REQUIRE(r->IsEmpty());
    }
}
