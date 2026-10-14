// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests: Stereon's built-in geometry against OpenCASCADE's
// (ADR-0015's independent referee). The same line, circle, plane and
// cylinder are built in both kernels and evaluated, with first and second
// derivatives, over a sweep of parameters; the two must agree to a few units
// in the last place. Both kernels use the same parametrisations here: a
// line by arc length along a unit direction, a circle and cylinder by angle
// from the X axis around the normal, a plane by its X and Y axes.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
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
using Stereon::Cylinder;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
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
