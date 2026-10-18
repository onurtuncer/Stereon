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
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_Line.hxx>
#include <Geom_Parabola.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

import Stereon.Surface;

#include "OcctConversion.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Hyperbola;
using Stereon::Line;
using Stereon::Parabola;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Torus;
using Stereon::Vec3;
using StereonTest::ToOcct;

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

/// Clamped and unclamped, polynomial and rational, low and high degree,
/// with a double interior knot: the cases the algorithms branch on.
std::vector<BSplineCurve> SomeBSplineCurves()
{
    const std::vector<Point3> poles{Point3{0.0, 0.0, 0.0},  Point3{1.0, 2.0, -1.0},
                                    Point3{2.5, -1.0, 0.5}, Point3{3.0, 1.5, 2.0},
                                    Point3{4.5, 0.0, 1.0},  Point3{5.0, 2.0, -0.5},
                                    Point3{6.0, -1.5, 0.0}};
    const std::vector<double> weights{1.0, 2.0, 0.5, 1.5, 3.0, 0.75, 1.0};
    return {
        {3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0}, poles, {}},
        {3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0}, poles, weights},
        {2, {-1.0, -0.5, 0.0, 0.75, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0}, poles, weights},
        {5, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, poles, weights},
    };
}

TEST_CASE("B-spline curves agree with Geom_BSplineCurve", "[occt][curve][bspline]")
{
    for (const BSplineCurve& curve : SomeBSplineCurves()) {
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(curve);
        const Stereon::Interval domain = Domain(curve);
        REQUIRE(occt->FirstParameter() == domain.GetLower());
        REQUIRE(occt->LastParameter() == domain.GetUpper());
        for (int i = 0; i < SampleCount; ++i) {
            const double t =
                domain.GetLower() + i * (domain.GetUpper() - domain.GetLower()) / (SampleCount - 1);
            gp_Pnt p;
            gp_Vec d1;
            gp_Vec d2;
            gp_Vec d3;
            occt->D3(t, p, d1, d2, d3);
            REQUIRE(Agree(Eval(curve, t), p));
            // Derivatives grow with the knot spacing; compare relative to
            // their size, a few ulps of the largest term.
            REQUIRE(Agree(Deriv(curve, t, 1) / (1.0 + Length(Deriv(curve, t, 1))),
                          d1 / (1.0 + d1.Magnitude())));
            REQUIRE(Agree(Deriv(curve, t, 2) / (1.0 + Length(Deriv(curve, t, 2))),
                          d2 / (1.0 + d2.Magnitude())));
            REQUIRE(Agree(Deriv(curve, t, 3) / (1.0 + Length(Deriv(curve, t, 3))),
                          d3 / (1.0 + d3.Magnitude())));
        }
    }
}

TEST_CASE("B-spline projections agree with GeomAPI_ProjectPointOnCurve", "[occt][curve][bspline]")
{
    const Frame f{Point3{3.0, 0.0, 0.5}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
                  Vec3{0.0, 0.0, 1.0}};
    for (const BSplineCurve& curve : SomeBSplineCurves()) {
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(curve);
        for (int i = 0; i < ProbeCount; ++i) {
            for (int j = 0; j < ProbeCount; ++j) {
                for (int k = 0; k < ProbeCount; ++k) {
                    const Point3 p = Probe(f, i, j, k);
                    GeomAPI_ProjectPointOnCurve extrema(ToOcct(p), occt);
                    const double ours = Distance(Eval(curve, Project(curve, p)), p);
                    // OCCT reports interior extrema only; the ends count too.
                    double reference = std::min(Distance(Eval(curve, Domain(curve).GetLower()), p),
                                                Distance(Eval(curve, Domain(curve).GetUpper()), p));
                    if (extrema.NbPoints() > 0) {
                        reference = std::min(reference, extrema.LowerDistance());
                    }
                    REQUIRE(ours <= reference + 1e-9 * (1.0 + reference));
                    REQUIRE(SameDistance(ours, reference));
                }
            }
        }
    }
}

/// Polynomial and rational, clamped and unclamped, with a double knot.
std::vector<BSplineSurface> SomeBSplineSurfaces()
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) {
            poles.push_back({1.5 * i, 1.25 * j, 0.5 * ((i * 7 + j * 3) % 5) - 1.0});
            weights.push_back(1.0 + 0.25 * ((i + 2 * j) % 4));
        }
    }
    const std::vector<double> clampedU{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0};
    const std::vector<double> clampedV{0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0};
    const std::vector<double> openU{-1.0, -0.5, 0.0, 0.25, 0.25, 1.0, 1.5, 2.0};
    return {
        {3, 2, clampedU, clampedV, poles, {}},
        {3, 2, clampedU, clampedV, poles, weights},
        {2, 2, openU, clampedV, poles, weights},
    };
}

TEST_CASE("B-spline surfaces agree with Geom_BSplineSurface", "[occt][surface][bspline]")
{
    for (const BSplineSurface& surface : SomeBSplineSurfaces()) {
        const occ::handle<Geom_BSplineSurface> occt = ToOcct(surface);
        const Stereon::Box2 d = Domain(surface);
        double u0 = 0.0;
        double u1 = 0.0;
        double v0 = 0.0;
        double v1 = 0.0;
        occt->Bounds(u0, u1, v0, v1);
        REQUIRE(u0 == d.X.GetLower());
        REQUIRE(u1 == d.X.GetUpper());
        REQUIRE(v0 == d.Y.GetLower());
        REQUIRE(v1 == d.Y.GetUpper());
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                const Point2 uv{u0 + i * (u1 - u0) / 10.0, v0 + j * (v1 - v0) / 10.0};
                CheckSecondOrder(surface, *occt, uv);
            }
        }
    }
}

TEST_CASE("B-spline surface projections are as close as GeomAPI_ProjectPointOnSurf's",
          "[occt][surface][bspline]")
{
    const Frame f{Point3{3.0, 2.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
                  Vec3{0.0, 0.0, 1.0}};
    for (const BSplineSurface& surface : SomeBSplineSurfaces()) {
        const occ::handle<Geom_BSplineSurface> occt = ToOcct(surface);
        for (int i = 0; i < ProbeCount; ++i) {
            for (int j = 0; j < ProbeCount; ++j) {
                for (int k = 0; k < ProbeCount; ++k) {
                    const Point3 p = Probe(f, i, j, k);
                    const double ours = Distance(Eval(surface, Project(surface, p)), p);
                    // OCCT's extrema are interior stationary points and can
                    // miss a closest point on the boundary, so ours must be
                    // at least as close, not equal.
                    GeomAPI_ProjectPointOnSurf extrema(ToOcct(p), occt);
                    if (extrema.NbPoints() > 0) {
                        REQUIRE(ours <= extrema.LowerDistance() + 1e-9 * (1.0 + ours));
                    }
                }
            }
        }
    }
}

/// Values, D1 to D3, and projections of the probe lattice by distance.
template <typename C>
void CheckOpenConic(const C& curve, const occ::handle<Geom_Curve>& occt, const Frame& f,
                    const double range)
{
    for (int i = 0; i < SampleCount; ++i) {
        const double t = -range + i * (2.0 * range / (SampleCount - 1));
        gp_Pnt p;
        gp_Vec d1;
        gp_Vec d2;
        gp_Vec d3;
        occt->D3(t, p, d1, d2, d3);
        REQUIRE(Agree(Eval(curve, t), p));
        REQUIRE(Agree(Deriv(curve, t, 1), d1));
        REQUIRE(Agree(Deriv(curve, t, 2), d2));
        REQUIRE(Agree(Deriv(curve, t, 3), d3));
    }
    for (int i = 0; i < ProbeCount; ++i) {
        for (int j = 0; j < ProbeCount; ++j) {
            for (int k = 0; k < ProbeCount; ++k) {
                const Point3 p = Probe(f, i, j, k);
                GeomAPI_ProjectPointOnCurve extrema(ToOcct(p), occt);
                REQUIRE(extrema.NbPoints() > 0);
                const double ours = Distance(Eval(curve, Project(curve, p)), p);
                REQUIRE(SameDistance(ours, extrema.LowerDistance()));
            }
        }
    }
}

TEST_CASE("hyperbolas agree with Geom_Hyperbola", "[occt][curve]")
{
    const Frame f = TiltedFrame();
    const Hyperbola hyperbola{f.Origin, 2.5, 1.5, f.X, f.Y};
    const occ::handle<Geom_Hyperbola> occt =
        new Geom_Hyperbola(gp_Ax2(ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)), 2.5, 1.5);
    CheckOpenConic(hyperbola, occt, f, 3.0);
}

TEST_CASE("parabolas agree with Geom_Parabola", "[occt][curve]")
{
    const Frame f = TiltedFrame();
    const Parabola parabola{f.Origin, 0.75, f.X, f.Y};
    const occ::handle<Geom_Parabola> occt =
        new Geom_Parabola(gp_Ax2(ToOcct(f.Origin), ToOcctDir(f.Z), ToOcctDir(f.X)), 0.75);
    CheckOpenConic(parabola, occt, f, 8.0);
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

TEST_CASE("swept surfaces agree with OCCT construction evaluators", "[occt][surface][swept]")
{
    const BSplineCurve profile{
        2, {0, 0, 0, 1, 1, 1}, {{2, 0, -1}, {3, 1, 1}, {2, 0, 3}}, {1, 2, 1}};
    const Vec3 direction = Normalized(Vec3{1, -2, 3});
    const auto extruded =
        Stereon::MakeExtrusionSurface(profile, direction, Stereon::Interval{-2, 3});
    REQUIRE(extruded.has_value());
    const Geom_SurfaceOfLinearExtrusion occtExtrusion(
        ToOcct(profile), gp_Dir(direction.X, direction.Y, direction.Z));
    for (const double u : {0.0, 0.17, 0.6, 1.0}) {
        for (const double v : {-2.0, 0.3, 3.0}) {
            REQUIRE(Agree(Eval(*extruded, Point2{u, v}), occtExtrusion.Value(u, v)));
        }
    }
    const Point3 origin{1, -2, 0.5};
    const Vec3 axis{2, -1, 3};
    const double start = 0.2;
    const double span = 1.0;
    const auto revolved = Stereon::MakeRevolutionSurface(profile, origin, axis,
                                                         Stereon::Interval{start, start + span});
    REQUIRE(revolved.has_value());
    const Geom_SurfaceOfRevolution occtRevolution(
        ToOcct(profile),
        gp_Ax1(gp_Pnt(origin.X, origin.Y, origin.Z), gp_Dir(axis.X, axis.Y, axis.Z)));
    for (const double t : {0.0, 0.13, 0.5, 0.79, 1.0}) {
        const double angle =
            start + 2.0 * std::atan2(t * std::sin(span / 2.0), 1.0 - t + t * std::cos(span / 2.0));
        for (const double v : {0.0, 0.17, 0.6, 1.0}) {
            REQUIRE(Agree(Eval(*revolved, Point2{start + t * span, v}),
                          occtRevolution.Value(angle, v)));
        }
    }
}
