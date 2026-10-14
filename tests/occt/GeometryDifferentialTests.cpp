// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests: Stereon's built-in geometry against OpenCASCADE's
// (ADR-0015's independent referee). The same line, circle, ellipse, plane,
// cylinder, cone, sphere and torus are built in both kernels and evaluated,
// with first and second derivatives, over a sweep of parameters; the two
// must agree to a few units in the last place. Both kernels use the same
// parametrisations here: a line by arc length along a unit direction, the
// conics and surfaces of revolution by angle from the X axis around the
// normal, a plane by its X and Y axes, a cone by distance along a
// generator. Projections are compared by distance with OCCT's extrema,
// which reach the closest point by a different algorithm.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

import Stereon.Surface;

using Stereon::Circle;
using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

/// Agreement to a relative 1e-13: a few units in the last place for the
/// magnitudes below, which both kernels compute by different expressions.
constexpr double RelativeTolerance = 1e-13;

bool Agree(const Point3& a, const gp_Pnt& b)
{
    const double scale = 1.0 + std::abs(b.X()) + std::abs(b.Y()) + std::abs(b.Z());
    return std::abs(a.X - b.X()) <= RelativeTolerance * scale &&
           std::abs(a.Y - b.Y()) <= RelativeTolerance * scale &&
           std::abs(a.Z - b.Z()) <= RelativeTolerance * scale;
}

bool Agree(const Vec3& a, const gp_Vec& b)
{
    const double scale = 1.0 + std::abs(b.X()) + std::abs(b.Y()) + std::abs(b.Z());
    return std::abs(a.X - b.X()) <= RelativeTolerance * scale &&
           std::abs(a.Y - b.Y()) <= RelativeTolerance * scale &&
           std::abs(a.Z - b.Z()) <= RelativeTolerance * scale;
}

gp_Pnt ToOcct(const Point3& p)
{
    return {p.X, p.Y, p.Z};
}

gp_Dir ToOcctDir(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

/// A right-handed orthonormal frame, tilted away from the coordinate axes so
/// that every component of every result is non-trivial.
struct Frame {
    Point3 Origin;
    Vec3 X;
    Vec3 Y;
    Vec3 Z;
};

Frame TiltedFrame()
{
    const Stereon::Transform tilt = Stereon::MakeRotation(Vec3{1.0, -2.0, 0.5}, 0.7);
    return {Point3{1.25, -0.5, 3.0}, Transformed(Vec3{1.0, 0.0, 0.0}, tilt),
            Transformed(Vec3{0.0, 1.0, 0.0}, tilt), Transformed(Vec3{0.0, 0.0, 1.0}, tilt)};
}

gp_Ax3 ToOcctFrame(const Frame& f)
{
    return {ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)};
}

/// Parameters across two periods and both signs.
constexpr int SampleCount = 41;

double Sample(const int i)
{
    return -2.0 * std::numbers::pi + i * (4.0 * std::numbers::pi / (SampleCount - 1));
}

/// Points around the geometry, inside and outside, on a lattice tilted
/// against the frame: the probes for the projection comparisons.
constexpr int ProbeCount = 7;

Point3 Probe(const Frame& f, const int i, const int j, const int k)
{
    const auto at = [](const int n) { return -4.5 + 1.5 * n + 0.1 * n * n; };
    return f.Origin + at(i) * f.X + at(j) * f.Y + at(k) * f.Z;
}

/// Our projection must find a point as close as OCCT's closest one: no
/// farther than OCCT's distance plus round-off, and, since OCCT's is a true
/// point of the geometry too, no closer than it minus round-off.
bool SameDistance(const double ours, const double occt)
{
    return std::abs(ours - occt) <= 1e-9 * (1.0 + occt);
}

/// Every first and second partial of a Stereon surface against OCCT's D2,
/// and one third-order partial against DN.
template <typename S>
void CheckSecondOrder(const S& surface, const Geom_Surface& occt, const Point2& uv)
{
    gp_Pnt p;
    gp_Vec du;
    gp_Vec dv;
    gp_Vec duu;
    gp_Vec dvv;
    gp_Vec duv;
    occt.D2(uv.X, uv.Y, p, du, dv, duu, dvv, duv);
    REQUIRE(Agree(Eval(surface, uv), p));
    REQUIRE(Agree(Deriv(surface, uv, 1, 0), du));
    REQUIRE(Agree(Deriv(surface, uv, 0, 1), dv));
    REQUIRE(Agree(Deriv(surface, uv, 2, 0), duu));
    REQUIRE(Agree(Deriv(surface, uv, 0, 2), dvv));
    REQUIRE(Agree(Deriv(surface, uv, 1, 1), duv));
    REQUIRE(Agree(Deriv(surface, uv, 2, 1), occt.DN(uv.X, uv.Y, 2, 1)));
}

/// Projections of the probe lattice agree with OCCT's extrema by distance.
template <typename S>
void CheckProjections(const S& surface, const occ::handle<Geom_Surface>& occt, const Frame& f)
{
    for (int i = 0; i < ProbeCount; ++i) {
        for (int j = 0; j < ProbeCount; ++j) {
            for (int k = 0; k < ProbeCount; ++k) {
                const Point3 p = Probe(f, i, j, k);
                GeomAPI_ProjectPointOnSurf extrema(ToOcct(p), occt);
                REQUIRE(extrema.NbPoints() > 0);
                const double ours = Distance(Eval(surface, Project(surface, p)), p);
                REQUIRE(SameDistance(ours, extrema.LowerDistance()));
            }
        }
    }
}

} // namespace

TEST_CASE("lines agree with Geom_Line", "[occt][curve]")
{
    const Frame f = TiltedFrame();
    const Line line{f.Origin, f.X}; // unit direction: OCCT parametrises by arc length
    const occ::handle<Geom_Line> occt = new Geom_Line(ToOcct(f.Origin), ToOcctDir(f.X));
    for (int i = 0; i < SampleCount; ++i) {
        const double t = Sample(i);
        gp_Pnt p;
        gp_Vec d1;
        occt->D1(t, p, d1);
        REQUIRE(Agree(Eval(line, t), p));
        REQUIRE(Agree(Deriv(line, t, 1), d1));
        REQUIRE(Agree(Deriv(line, t, 2), occt->DN(t, 2)));
    }
}

TEST_CASE("circles agree with Geom_Circle", "[occt][curve]")
{
    const Frame f = TiltedFrame();
    const double radius = 2.75;
    const Circle circle{f.Origin, radius, f.X, f.Y};
    const occ::handle<Geom_Circle> occt =
        new Geom_Circle(gp_Ax2(ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)), radius);
    for (int i = 0; i < SampleCount; ++i) {
        const double t = Sample(i);
        gp_Pnt p;
        gp_Vec d1;
        gp_Vec d2;
        occt->D2(t, p, d1, d2);
        REQUIRE(Agree(Eval(circle, t), p));
        REQUIRE(Agree(Deriv(circle, t, 1), d1));
        REQUIRE(Agree(Deriv(circle, t, 2), d2));
        REQUIRE(Agree(Deriv(circle, t, 3), occt->DN(t, 3)));
    }
    REQUIRE(occt->Period() == Stereon::CirclePeriod);
}

TEST_CASE("ellipses agree with Geom_Ellipse", "[occt][curve]")
{
    const Frame f = TiltedFrame();
    const Ellipse ellipse{f.Origin, 3.5, 1.25, f.X, f.Y};
    const occ::handle<Geom_Ellipse> occt =
        new Geom_Ellipse(gp_Ax2(ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)), 3.5, 1.25);
    for (int i = 0; i < SampleCount; ++i) {
        const double t = Sample(i);
        gp_Pnt p;
        gp_Vec d1;
        gp_Vec d2;
        occt->D2(t, p, d1, d2);
        REQUIRE(Agree(Eval(ellipse, t), p));
        REQUIRE(Agree(Deriv(ellipse, t, 1), d1));
        REQUIRE(Agree(Deriv(ellipse, t, 2), d2));
        REQUIRE(Agree(Deriv(ellipse, t, 3), occt->DN(t, 3)));
    }
    for (int i = 0; i < ProbeCount; ++i) {
        for (int j = 0; j < ProbeCount; ++j) {
            for (int k = 0; k < ProbeCount; ++k) {
                const Point3 p = Probe(f, i, j, k);
                GeomAPI_ProjectPointOnCurve extrema(ToOcct(p), occt);
                REQUIRE(extrema.NbPoints() > 0);
                const double ours = Distance(Eval(ellipse, Project(ellipse, p)), p);
                REQUIRE(SameDistance(ours, extrema.LowerDistance()));
            }
        }
    }
}

TEST_CASE("planes agree with Geom_Plane", "[occt][surface]")
{
    const Frame f = TiltedFrame();
    const Plane plane{f.Origin, f.X, f.Y};
    const occ::handle<Geom_Plane> occt = new Geom_Plane(ToOcctFrame(f));
    for (int i = 0; i < SampleCount; ++i) {
        const Point2 uv{Sample(i), 0.5 * Sample(SampleCount - 1 - i)};
        gp_Pnt p;
        gp_Vec du;
        gp_Vec dv;
        occt->D1(uv.X, uv.Y, p, du, dv);
        REQUIRE(Agree(Eval(plane, uv), p));
        REQUIRE(Agree(Deriv(plane, uv, 1, 0), du));
        REQUIRE(Agree(Deriv(plane, uv, 0, 1), dv));
    }
}

TEST_CASE("cylinders agree with Geom_CylindricalSurface", "[occt][surface]")
{
    const Frame f = TiltedFrame();
    const double radius = 1.5;
    const Cylinder cylinder{f.Origin, radius, f.X, f.Y};
    const occ::handle<Geom_CylindricalSurface> occt =
        new Geom_CylindricalSurface(ToOcctFrame(f), radius);
    for (int i = 0; i < SampleCount; ++i) {
        const Point2 uv{Sample(i), 0.75 * Sample(SampleCount - 1 - i)};
        gp_Pnt p;
        gp_Vec du;
        gp_Vec dv;
        gp_Vec duu;
        gp_Vec dvv;
        gp_Vec duv;
        occt->D2(uv.X, uv.Y, p, du, dv, duu, dvv, duv);
        REQUIRE(Agree(Eval(cylinder, uv), p));
        REQUIRE(Agree(Deriv(cylinder, uv, 1, 0), du));
        REQUIRE(Agree(Deriv(cylinder, uv, 0, 1), dv));
        REQUIRE(Agree(Deriv(cylinder, uv, 2, 0), duu));
        REQUIRE(Agree(Deriv(cylinder, uv, 0, 2), dvv));
        REQUIRE(Agree(Deriv(cylinder, uv, 1, 1), duv));
    }
    REQUIRE(occt->UPeriod() == Stereon::CirclePeriod);
}

TEST_CASE("cones agree with Geom_ConicalSurface", "[occt][surface]")
{
    const Frame f = TiltedFrame();
    const Cone cone{f.Origin, 1.25, 0.6, f.X, f.Y};
    const occ::handle<Geom_ConicalSurface> occt =
        new Geom_ConicalSurface(ToOcctFrame(f), 0.6, 1.25);
    for (int i = 0; i < SampleCount; ++i) {
        CheckSecondOrder(cone, *occt, Point2{Sample(i), 0.75 * Sample(SampleCount - 1 - i)});
    }
    REQUIRE(Agree(Stereon::Apex(cone), occt->Apex()));
    CheckProjections(cone, occt, f);
}

TEST_CASE("spheres agree with Geom_SphericalSurface", "[occt][surface]")
{
    const Frame f = TiltedFrame();
    const Sphere sphere{f.Origin, 2.25, f.X, f.Y};
    const occ::handle<Geom_SphericalSurface> occt = new Geom_SphericalSurface(ToOcctFrame(f), 2.25);
    for (int i = 0; i < SampleCount; ++i) {
        // Latitudes inside [-pi/2, pi/2].
        CheckSecondOrder(sphere, *occt, Point2{Sample(i), 0.24 * Sample(SampleCount - 1 - i)});
    }
    CheckProjections(sphere, occt, f);
}

TEST_CASE("tori agree with Geom_ToroidalSurface", "[occt][surface]")
{
    const Frame f = TiltedFrame();
    const Torus torus{f.Origin, 3.0, 0.8, f.X, f.Y};
    const occ::handle<Geom_ToroidalSurface> occt =
        new Geom_ToroidalSurface(ToOcctFrame(f), 3.0, 0.8);
    for (int i = 0; i < SampleCount; ++i) {
        CheckSecondOrder(torus, *occt, Point2{Sample(i), Sample(SampleCount - 1 - i)});
    }
    REQUIRE(occt->UPeriod() == Stereon::CirclePeriod);
    REQUIRE(occt->VPeriod() == Stereon::CirclePeriod);
    CheckProjections(torus, occt, f);
}

TEST_CASE("the erased wrappers agree with OCCT too", "[occt][erased]")
{
    const Frame f = TiltedFrame();
    const Stereon::Curve circle(Circle{f.Origin, 2.0, f.X, f.Y});
    const occ::handle<Geom_Circle> occt =
        new Geom_Circle(gp_Ax2(ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)), 2.0);
    const Stereon::Surface cylinder(Cylinder{f.Origin, 2.0, f.X, f.Y});
    const occ::handle<Geom_CylindricalSurface> occtCylinder =
        new Geom_CylindricalSurface(ToOcctFrame(f), 2.0);
    for (int i = 0; i < SampleCount; ++i) {
        const double t = Sample(i);
        REQUIRE(Agree(Eval(circle, t), occt->Value(t)));
        REQUIRE(Agree(Eval(cylinder, Point2{t, t}), occtCylinder->Value(t, t)));
    }
}
