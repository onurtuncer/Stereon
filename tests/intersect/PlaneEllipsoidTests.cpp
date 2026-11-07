// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A plane against an ellipsoid (ADR-0048 §3.3). The configuration is
// decided exactly on the input doubles: a plane 2^-40 beyond a vertex is
// apart, one at it touches, one 2^-40 inside cuts a section. A section is
// an exact ellipse with its exact pcurve on the plane and a fitted pcurve
// on the ellipsoid, certified within the resolution and cut at the chart's
// seam; a plane through a pole of the chart fails closed. With equal
// semi-axes the sections lie on the sphere.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::ContactType;
using Stereon::Ellipse;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::IntersectionBranch;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceSurfaceIntersection;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Vec3 X{1.0, 0.0, 0.0};
const Vec3 Y{0.0, 1.0, 0.0};
const Vec3 Z{0.0, 0.0, 1.0};

/// The ellipsoid about the origin with semi-axes 3, 2 and 1 along x, y, z.
const Ellipsoid stretched{Point3{}, 3.0, 2.0, 1.0, X, Y};

Surface Horizontal(const double z)
{
    return Surface{Plane{Point3{0.0, 0.0, z}, X, Y}};
}

/// (x / A)² + (y / B)² + (z / C)² - 1 in the ellipsoid's frame.
double Level(const Ellipsoid& e, const Point3& p)
{
    const Vec3 w = p - e.Center;
    const double x = Dot(w, e.XAxis) / e.A;
    const double y = Dot(w, e.YAxis) / e.B;
    const double z = Dot(w, Axis(e)) / e.C;
    return x * x + y * y + z * z - 1.0;
}

/// Every branch certified, crossing, its 3D curve an ellipse, and its
/// curve and pcurves agreeing with both surfaces along it; its points on
/// the plane and on the ellipsoid's implicit equation. Returns the total
/// angle the branches cover.
double RequireSections(const SurfaceSurfaceIntersection& result, const Surface& a, const Surface& b,
                       const Ellipsoid& e, const Plane& plane)
{
    double covered = 0.0;
    for (const IntersectionBranch& branch : result.Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= ctx.GetResolution());
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.Transition == SideTransition::OppositeSide);
        REQUIRE(branch.Curve3d.Target<Ellipse>() != nullptr);
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b) <= ctx.GetResolution());
        const Vec3 n = Normalized(Cross(plane.XAxis, plane.YAxis));
        for (const double t : IntersectTest::SamplesOf(branch.Range, 32)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(std::abs(Dot(p - plane.Origin, n)) <= 1e-12 * (1.0 + Length(p - Point3{})));
            REQUIRE(std::abs(Level(e, p)) <= 1e-12);
        }
        covered += branch.Range.GetUpper() - branch.Range.GetLower();
    }
    return covered;
}

} // namespace

TEST_CASE("a horizontal plane cuts an ellipsoid in a scaled copy of its equator",
          "[intersect][ellipsoid][plane-ellipsoid]")
{
    // At height z the section is the ellipse of semi-axes 3 k and 2 k,
    // k = sqrt(1 - z²), about (0, 0, z). It goes round the axis, so it
    // crosses the chart's seam once: one branch, the whole turn.
    for (const double z : {0.0, 0.5, -0.75}) {
        INFO("z = " << z);
        const Surface plane = Horizontal(z);
        const auto result = Intersect(plane, Surface{stretched}, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        REQUIRE(result->Points.empty());
        REQUIRE(result->Branches.size() == 1);
        const double covered =
            RequireSections(*result, plane, Surface{stretched}, stretched, *plane.Target<Plane>());
        REQUIRE(std::abs(covered - 2.0 * pi) <= 1e-12);
        const Ellipse& section = *result->Branches.front().Curve3d.Target<Ellipse>();
        const double k = std::sqrt(1.0 - z * z);
        REQUIRE(std::abs(section.MajorRadius - 3.0 * k) <= 1e-14);
        REQUIRE(std::abs(section.MinorRadius - 2.0 * k) <= 1e-14);
        REQUIRE(Distance(section.Center, Point3{0.0, 0.0, z}) <= 1e-14);
        // The pcurve on the ellipsoid keeps the latitude asin(z).
        const IntersectionBranch& branch = result->Branches.front();
        for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
            const Point2 uv = Eval(branch.OnB, t);
            REQUIRE(std::abs(uv.Y - std::asin(z)) <= 1e-9);
            REQUIRE(uv.X >= 0.0);
            REQUIRE(uv.X <= 2.0 * pi);
        }
    }
}

TEST_CASE("a plane is apart from, touches or cuts an ellipsoid, decided exactly",
          "[intersect][ellipsoid][plane-ellipsoid]")
{
    const double gap = 0x1p-40;
    // Over the pole (0, 0, 1).
    const auto apart = Intersect(Horizontal(1.0 + gap), Surface{stretched}, ctx);
    REQUIRE(apart.has_value());
    REQUIRE(apart->IsEmpty());
    const auto touch = Intersect(Horizontal(1.0), Surface{stretched}, ctx);
    INFO(IntersectTest::ErrorText(touch));
    REQUIRE(touch.has_value());
    REQUIRE(touch->Branches.empty());
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(Distance(touch->Points.front().Position, Point3{0.0, 0.0, 1.0}) <= 1e-15);
    REQUIRE(std::abs(touch->Points.front().OnB.Y - pi / 2.0) <= 1e-15);
    const auto inside = Intersect(Horizontal(1.0 - gap), Surface{stretched}, ctx);
    INFO(IntersectTest::ErrorText(inside));
    REQUIRE(inside.has_value());
    REQUIRE(inside->Branches.size() == 1);
    const Surface near = Horizontal(1.0 - gap);
    RequireSections(*inside, near, Surface{stretched}, stretched, *near.Target<Plane>());

    // Beside the end of the longest semi-axis, the plane x = 3.
    const auto beside = [](const double x) { return Surface{Plane{Point3{x, 0.0, 0.0}, Y, Z}}; };
    REQUIRE(Intersect(beside(3.0 + gap), Surface{stretched}, ctx)->IsEmpty());
    const auto vertex = Intersect(beside(3.0), Surface{stretched}, ctx);
    REQUIRE(vertex.has_value());
    REQUIRE(vertex->Points.size() == 1);
    REQUIRE(Distance(vertex->Points.front().Position, Point3{3.0, 0.0, 0.0}) <= 1e-15);
}

TEST_CASE("an oblique plane through an ellipsoid is cut at the chart's seam",
          "[intersect][ellipsoid][plane-ellipsoid]")
{
    // The plane x + y + 2 z = 1/2, tilted across every axis. Its section
    // misses both poles; it crosses the seam, the half-plane y = 0, x > 0,
    // wherever it does, and the branches between cover the whole turn.
    const Plane tilted = IntersectTest::PlaneWithNormal(Point3{0.5, 0.0, 0.0}, Vec3{1.0, 1.0, 2.0});
    const auto result = Intersect(Surface{tilted}, Surface{stretched}, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_FALSE(result->Branches.empty());
    const double covered =
        RequireSections(*result, Surface{tilted}, Surface{stretched}, stretched, tilted);
    REQUIRE(std::abs(covered - 2.0 * pi) <= 1e-12);

    // The other way round: the same branches, the pcurves exchanged.
    const auto swapped = Intersect(Surface{stretched}, Surface{tilted}, ctx);
    REQUIRE(swapped.has_value());
    REQUIRE(swapped->Branches.size() == result->Branches.size());
    for (std::size_t i = 0; i < result->Branches.size(); ++i) {
        const IntersectionBranch& p = result->Branches[i];
        const IntersectionBranch& q = swapped->Branches[i];
        for (const double t : IntersectTest::SamplesOf(p.Range, 8)) {
            REQUIRE(Distance(Eval(p.OnA, t), Eval(q.OnB, t)) <= 1e-15);
            REQUIRE(Distance(Eval(p.OnB, t), Eval(q.OnA, t)) <= 1e-15);
        }
    }
}

TEST_CASE("a plane through a pole of the ellipsoid's chart fails closed",
          "[intersect][ellipsoid][plane-ellipsoid]")
{
    // The plane x = y holds the axis, so both poles: there the section's
    // longitude has no value, and no pcurve in the chart can follow it.
    const Plane through = IntersectTest::PlaneWithNormal(Point3{}, Vec3{1.0, -1.0, 0.0});
    const auto result = Intersect(Surface{through}, Surface{stretched}, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a plane within rounding of a pole fails closed as one through it",
          "[intersect][ellipsoid][plane-ellipsoid][regression]")
{
    // CI's counterexample of 2026-11-07 (gcc16-release, gcc16-asan and
    // macOS): the plane was made through the south pole of the sphere of
    // radius 3/2 about (-2, -2, -2), but its frame, normalised from
    // (-1, 4, 2) / 8, is rounded, so it misses the pole by a unit in the
    // last place. The exact pole test on the rounded plane then says no,
    // and the section's pcurve cannot be fitted within the budget so near
    // the pole. That failed with CertificationFailed, which the oracle did
    // not expect; it is IntersectionFailed, as for a plane exactly through
    // a pole.
    const Point3 center{-2.0, -2.0, -2.0};
    const Plane plane = IntersectTest::PlaneWithNormal(center + Vec3{0.25, -0.375, -0.625},
                                                       Vec3{-0.125, 0.5, 0.25});
    const auto result =
        Intersect(Surface{plane}, Surface{Ellipsoid{center, 1.5, 1.5, 1.5, X, Y}}, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a bounded plane box keeps the section's stretches inside it",
          "[intersect][ellipsoid][plane-ellipsoid]")
{
    // The equator, x² / 9 + y² / 4 = 1, against the plane's box [0, 4] x
    // [-4, 4]: the half with x >= 0, from (0, -2) round to (0, 2).
    const Surface plane = Horizontal(0.0);
    const Box2 half{Interval{0.0, 4.0}, Interval{-4.0, 4.0}};
    const auto result = Intersect(plane, half, Surface{stretched}, Domain(Surface{stretched}), ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_FALSE(result->Branches.empty());
    double covered = 0.0;
    for (const IntersectionBranch& branch : result->Branches) {
        REQUIRE(branch.IsCertified());
        for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
            REQUIRE(Eval(branch.OnA, t).X >= -1e-12);
        }
        covered += branch.Range.GetUpper() - branch.Range.GetLower();
    }
    REQUIRE(std::abs(covered - pi) <= 1e-9);
}

TEST_CASE("with equal semi-axes a plane's section lies on the sphere",
          "[intersect][ellipsoid][plane-ellipsoid][property]")
{
    // ADR-0048's oracle: an ellipsoid of equal semi-axes cut by a plane
    // gives the circle the sphere gives, as an ellipse of equal radii; its
    // points lie on the sphere and the plane. Either may fail closed where
    // the other decides (a plane through a pole of the ellipsoid's chart).
    rc::prop("for planes and spheres on dyadic coordinates", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-16, 17) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const double r = *rc::gen::inRange(4, 25) / 8.0;
        const Vec3 normal{coordinate(), coordinate(), coordinate()};
        RC_PRE(normal != Vec3{});
        const Plane plane = IntersectTest::PlaneWithNormal(
            center + Vec3{coordinate() / 2.0, coordinate() / 2.0, coordinate() / 2.0}, normal);
        const Ellipsoid ellipsoid{center, r, r, r, X, Y};
        const auto result = Intersect(Surface{plane}, Surface{ellipsoid}, ctx);
        if (!result) {
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        const auto onSphere = Intersect(Surface{plane}, Surface{Sphere{center, r, X, Y}}, ctx);
        RC_ASSERT(onSphere.has_value());
        RC_ASSERT(result->Points.size() == onSphere->Points.size());
        RC_ASSERT(result->Branches.empty() == onSphere->Branches.empty());
        for (const IntersectionBranch& branch : result->Branches) {
            RC_ASSERT(branch.IsCertified());
            for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
                const Point3 p = Eval(branch.Curve3d, t);
                RC_ASSERT(std::abs(Distance(p, center) - r) <= 1e-12 * (1.0 + r));
            }
        }
    });
}

TEST_CASE("every plane through an ellipsoid turned at random gives certified sections or fails "
          "closed",
          "[intersect][ellipsoid][plane-ellipsoid][property]")
{
    rc::prop("for planes and ellipsoids on dyadic coordinates", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-16, 17) / 8.0; };
        const auto semi = [] { return *rc::gen::inRange(4, 25) / 8.0; };
        const Point3 center{coordinate(), coordinate(), coordinate()};
        const Vec3 axis{coordinate(), coordinate(), coordinate()};
        const Vec3 normal{coordinate(), coordinate(), coordinate()};
        RC_PRE(axis != Vec3{} && normal != Vec3{});
        const auto frame = IntersectTest::CylinderAbout(center, axis, 1.0);
        const Ellipsoid ellipsoid{center, semi(), semi(), semi(), frame.XAxis, frame.YAxis};
        const Plane plane = IntersectTest::PlaneWithNormal(
            center + Vec3{coordinate() / 2.0, coordinate() / 2.0, coordinate() / 2.0}, normal);
        const auto result = Intersect(Surface{plane}, Surface{ellipsoid}, ctx);
        if (!result) {
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            RC_CLASSIFY(true, "failed closed");
            return;
        }
        RC_CLASSIFY(!result->Branches.empty(), "section");
        RC_CLASSIFY(result->IsEmpty(), "apart");
        double covered = 0.0;
        for (const IntersectionBranch& branch : result->Branches) {
            RC_ASSERT(branch.IsCertified());
            RC_ASSERT(IntersectTest::PcurveDeviation(branch, Surface{plane}, Surface{ellipsoid}) <=
                      ctx.GetResolution());
            for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
                RC_ASSERT(std::abs(Level(ellipsoid, Eval(branch.Curve3d, t))) <= 1e-11);
                const Point2 uv = Eval(branch.OnB, t);
                RC_ASSERT(uv.X >= 0.0 && uv.X <= 2.0 * pi);
            }
            covered += branch.Range.GetUpper() - branch.Range.GetLower();
        }
        if (!result->Branches.empty()) {
            RC_ASSERT(std::abs(covered - 2.0 * pi) <= 1e-12);
        }
    });
}
