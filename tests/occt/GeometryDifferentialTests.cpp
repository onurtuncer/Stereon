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
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
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
#include <Geom_RectangularTrimmedSurface.hxx>
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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

/// Agreement after an operation that rebuilds the poles (elevation,
/// concatenation): each kernel rounds its own sequence of convex
/// combinations, so a relative 1e-11 rather than Agree's few ulps.
bool AgreeRebuilt(const Point3& a, const gp_Pnt& b)
{
    return Distance(a, Point3{b.X(), b.Y(), b.Z()}) <= 1e-11 * (1.0 + Length(Stereon::ToVec(a)));
}

TEST_CASE("degree elevation agrees with Geom_BSplineCurve::IncreaseDegree",
          "[occt][curve][bspline][operations]")
{
    for (const BSplineCurve& curve : SomeBSplineCurves()) {
        for (const int degree : {curve.Degree + 1, curve.Degree + 3}) {
            const BSplineCurve ours = *Stereon::ElevatedDegree(curve, degree);
            // OCCT is given the clamped form, the same curve on the same
            // domain: on the unclamped fixture its IncreaseDegree disagrees
            // with ours, while ours reproduces the input (checked below and
            // by the property tests of tests/curve/BSplineOperationsTests.cpp).
            const occ::handle<Geom_BSplineCurve> occt = ToOcct(Stereon::Clamped(curve));
            occt->IncreaseDegree(degree);
            REQUIRE(occt->Degree() == ours.Degree);
            const Stereon::Interval domain = Domain(ours);
            for (int i = 0; i < SampleCount; ++i) {
                const double t = domain.GetLower() +
                                 i * (domain.GetUpper() - domain.GetLower()) / (SampleCount - 1);
                REQUIRE(AgreeRebuilt(Eval(ours, t), occt->Value(t)));
                REQUIRE(Distance(Eval(ours, t), Eval(curve, t)) <=
                        1e-11 * (1.0 + Length(Stereon::ToVec(Eval(curve, t)))));
            }
            // Both keep the elevated curve minimal: the same number of poles,
            // since both start from the clamped form.
            REQUIRE(occt->NbPoles() == static_cast<int>(ours.Poles.size()));
        }
    }
}

TEST_CASE("degree elevation over a tiny knot interval agrees with IncreaseDegree",
          "[occt][curve][bspline][operations]")
{
    // The audit fixture of 2026-10-28 (A1): degree 9, a first interior knot
    // at 1e-8 and poles inside a 1 x 2 box, which the unchecked knot
    // removals of the elevation once moved by 3.8e7. The elevated curve must
    // agree with OCCT's and with the input. OCCT's IncreaseDegree is a global
    // algorithm and reaches the minimal pole count; ours may keep a knot
    // copy where a removal would not be exact, so pole counts are not
    // compared here.
    std::vector<double> knots(10, 0.0);
    knots.insert(knots.end(), {1e-8, 0.5, 0.75});
    knots.insert(knots.end(), 10, 1.0);
    std::vector<Point3> poles;
    for (int i = 0; i < 13; ++i) {
        poles.push_back({static_cast<double>(i % 2), static_cast<double>(i % 3), 0.0});
    }
    const BSplineCurve curve = *Stereon::MakeBSplineCurve(9, knots, poles);
    for (const int degree : {10, 12}) {
        const auto elevated = Stereon::ElevatedDegree(curve, degree);
        REQUIRE(elevated.has_value());
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(curve);
        occt->IncreaseDegree(degree);
        REQUIRE(occt->Degree() == elevated->Degree);
        for (int i = 0; i < SampleCount; ++i) {
            const double t = i / (SampleCount - 1.0);
            REQUIRE(AgreeRebuilt(Eval(*elevated, t), occt->Value(t)));
            REQUIRE(Distance(Eval(*elevated, t), Eval(curve, t)) <= 1e-12);
        }
        REQUIRE(Distance(Eval(*elevated, 0.155), Eval(curve, 0.155)) <= 1e-12);
    }
}

TEST_CASE("knot removal agrees with Geom_BSplineCurve::RemoveKnot",
          "[occt][curve][bspline][operations]")
{
    for (const BSplineCurve& curve : SomeBSplineCurves()) {
        // A knot inserted 43% of the way along the domain, on no existing
        // knot of any fixture, is exactly removable; both kernels must remove
        // it and give the same poles back.
        const Stereon::Interval domain = Domain(curve);
        const double inserted = domain.GetLower() + 0.43 * (domain.GetUpper() - domain.GetLower());
        const BSplineCurve refined = *Stereon::InsertKnot(curve, inserted);
        const Stereon::KnotRemoval ours = *Stereon::RemoveKnot(refined, inserted, 1, 1e-9);
        REQUIRE(ours.Removed == 1);
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(refined);
        int index = 0;
        for (int k = 1; k <= occt->NbKnots(); ++k) {
            if (occt->Knot(k) == inserted) {
                index = k;
            }
        }
        REQUIRE(index > 0);
        REQUIRE(occt->RemoveKnot(index, 0, 1e-9));
        REQUIRE(occt->NbPoles() == static_cast<int>(ours.Curve.Poles.size()));
        for (std::size_t i = 0; i < ours.Curve.Poles.size(); ++i) {
            const gp_Pnt pole = occt->Pole(static_cast<int>(i) + 1);
            REQUIRE(Distance(ours.Curve.Poles[i], Point3{pole.X(), pole.Y(), pole.Z()}) <= 1e-9);
        }
        // A knot that is not removable: both refuse it within the tolerance.
        const std::vector<double> inside = [&] {
            std::vector<double> values;
            for (const double u : curve.Knots) {
                if (u > domain.GetLower() && u < domain.GetUpper() &&
                    (values.empty() || values.back() != u)) {
                    values.push_back(u);
                }
            }
            return values;
        }();
        for (const double u : inside) {
            const auto times = static_cast<int>(Stereon::KnotMultiplicity(curve.Knots, u));
            const Stereon::KnotRemoval tried = *Stereon::RemoveKnot(curve, u, times, 1e-9);
            const occ::handle<Geom_BSplineCurve> reference = ToOcct(curve);
            int at = 0;
            for (int k = 1; k <= reference->NbKnots(); ++k) {
                if (reference->Knot(k) == u) {
                    at = k;
                }
            }
            const bool occtRemoved = reference->RemoveKnot(at, times - 1, 1e-9);
            REQUIRE((tried.Removed > 0) == occtRemoved);
        }
    }
}

TEST_CASE("concatenation agrees with GeomConvert_CompCurveToBSplineCurve",
          "[occt][curve][bspline][operations]")
{
    const std::vector<BSplineCurve> curves = SomeBSplineCurves();
    for (std::size_t i = 0; i + 1 < curves.size(); ++i) {
        // Move the second curve so that it starts where the first ends: a corner.
        const BSplineCurve& first = curves[i];
        const Point3 end = Eval(first, Domain(first).GetUpper());
        BSplineCurve second = curves[i + 1];
        const Point3 start = Eval(second, Domain(second).GetLower());
        for (Point3& pole : second.Poles) {
            pole = pole + (end - start);
        }
        const BSplineCurve ours = *Stereon::Concatenated(first, second, 1e-9);
        // OCCT's composite builder throws on an unclamped knot vector ("BSpline
        // curve: # Poles and degree mismatch"); give it the clamped forms.
        GeomConvert_CompCurveToBSplineCurve compound(ToOcct(Stereon::Clamped(first)));
        REQUIRE(compound.Add(ToOcct(Stereon::Clamped(second)), 1e-9, true, false, 0));
        const occ::handle<Geom_BSplineCurve> occt = compound.BSplineCurve();
        REQUIRE(occt->Degree() == ours.Degree);
        // The same point set: every sample of ours lies on OCCT's curve, and
        // the ends agree. (OCCT keeps its own parametrisation of the second
        // piece, so the parameters are not compared.)
        const Stereon::Interval domain = Domain(ours);
        REQUIRE(AgreeRebuilt(Eval(ours, domain.GetLower()), occt->Value(occt->FirstParameter())));
        REQUIRE(AgreeRebuilt(Eval(ours, domain.GetUpper()), occt->Value(occt->LastParameter())));
        for (int k = 0; k < SampleCount; ++k) {
            const double t =
                domain.GetLower() + k * (domain.GetUpper() - domain.GetLower()) / (SampleCount - 1);
            const Point3 p = Eval(ours, t);
            GeomAPI_ProjectPointOnCurve projection(ToOcct(p), occt);
            double distance = std::numeric_limits<double>::infinity();
            if (projection.NbPoints() > 0) {
                distance = projection.LowerDistance();
            }
            const gp_Pnt occtStart = occt->StartPoint();
            const gp_Pnt occtEnd = occt->EndPoint();
            distance = std::min({distance,
                                 Distance(p, Point3{occtStart.X(), occtStart.Y(), occtStart.Z()}),
                                 Distance(p, Point3{occtEnd.X(), occtEnd.Y(), occtEnd.Z()})});
            // The projection converges to about 1e-9 itself; allow it that.
            REQUIRE(distance <= 1e-8 * (1.0 + Length(Stereon::ToVec(p))));
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

namespace {

/// Our projection onto `surface` must find a point at least as close as
/// OCCT's closest one on `occt`, a surface holding no point ours does not:
/// OCCT's extrema are interior stationary points and can miss a closest
/// point on a trimmed surface's boundary.
template <typename S>
void CheckProjectionsAtLeastAsClose(const S& surface, const occ::handle<Geom_Surface>& occt,
                                    const Frame& f)
{
    for (int i = 0; i < ProbeCount; ++i) {
        for (int j = 0; j < ProbeCount; ++j) {
            for (int k = 0; k < ProbeCount; ++k) {
                const Point3 p = Probe(f, i, j, k);
                const double ours = Distance(Eval(surface, Project(surface, p)), p);
                GeomAPI_ProjectPointOnSurf extrema(ToOcct(p), occt);
                if (extrema.NbPoints() > 0) {
                    REQUIRE(ours <= extrema.LowerDistance() + 1e-9 * (1.0 + ours));
                }
            }
        }
    }
}

} // namespace

TEST_CASE("the sweep wrappers agree with OCCT's unbounded sweeps", "[occt][surface][swept]")
{
    const BSplineCurve profile{
        2, {0, 0, 0, 1, 1, 1}, {{2, 0, -1}, {3, 1, 1}, {2, 0, 3}}, {1, 2, 1}};
    const occ::handle<Geom_BSplineCurve> occtProfile = ToOcct(profile);
    const Stereon::Curve erasedProfile(profile);

    const Vec3 direction = Normalized(Vec3{1, -2, 3});
    const Stereon::ExtrusionSurface extrusion = *Stereon::MakeExtrusion(erasedProfile, direction);
    const occ::handle<Geom_SurfaceOfLinearExtrusion> occtExtrusion =
        new Geom_SurfaceOfLinearExtrusion(occtProfile, ToOcctDir(direction));
    // Unbounded along the direction: v well beyond any finite patch.
    for (const double u : {0.0, 0.17, 0.6, 1.0}) {
        for (const double v : {-40.0, -2.0, 0.3, 3.0, 55.0}) {
            CheckSecondOrder(extrusion, *occtExtrusion, Point2{u, v});
        }
    }
    // OCCT's extrema need finite bounds: a patch of the same surface.
    const occ::handle<Geom_RectangularTrimmedSurface> finiteExtrusion =
        new Geom_RectangularTrimmedSurface(occtExtrusion, 0.0, 1.0, -30.0, 30.0);
    CheckProjectionsAtLeastAsClose(extrusion, finiteExtrusion, TiltedFrame());

    const Point3 origin{1, -2, 0.5};
    const Vec3 axis = Normalized(Vec3{2, -1, 3});
    const Stereon::RevolutionSurface revolution =
        *Stereon::MakeRevolution(erasedProfile, origin, axis);
    const occ::handle<Geom_SurfaceOfRevolution> occtRevolution =
        new Geom_SurfaceOfRevolution(occtProfile, gp_Ax1(ToOcct(origin), ToOcctDir(axis)));
    // The wrapper's u is the angle, as OCCT's is, right-handed about the
    // axis and periodic: two turns of samples either way.
    for (int i = 0; i < SampleCount; ++i) {
        for (const double v : {0.0, 0.17, 0.6, 1.0}) {
            CheckSecondOrder(revolution, *occtRevolution, Point2{Sample(i), v});
        }
    }
    CheckProjectionsAtLeastAsClose(revolution, occtRevolution, TiltedFrame());
}

TEST_CASE("trimmed surfaces agree with Geom_RectangularTrimmedSurface", "[occt][surface][trimmed]")
{
    SECTION("a cylinder across its seam")
    {
        const Frame f = TiltedFrame();
        const Stereon::Surface cylinder(Cylinder{f.Origin, 2.0, f.X, f.Y});
        const occ::handle<Geom_CylindricalSurface> occtCylinder =
            new Geom_CylindricalSurface(ToOcctFrame(f), 2.0);
        const Stereon::Box2 box{Stereon::Interval{5.0, 7.5}, Stereon::Interval{-1.0, 2.0}};
        const Stereon::TrimmedSurface trimmed = *Stereon::MakeTrimmedSurface(cylinder, box);
        const occ::handle<Geom_RectangularTrimmedSurface> occt =
            new Geom_RectangularTrimmedSurface(occtCylinder, 5.0, 7.5, -1.0, 2.0);
        for (const double u : {5.0, 5.9, 2.0 * std::numbers::pi, 6.7, 7.5}) {
            for (const double v : {-1.0, 0.4, 2.0}) {
                CheckSecondOrder(trimmed, *occt, Point2{u, v});
            }
        }
        CheckProjectionsAtLeastAsClose(trimmed, occt, f);
    }
    SECTION("a corner of each B-spline surface")
    {
        const Frame f{Point3{3.0, 2.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
                      Vec3{0.0, 0.0, 1.0}};
        for (const BSplineSurface& surface : SomeBSplineSurfaces()) {
            const Stereon::Box2 domain = Domain(surface);
            const auto part = [](const Stereon::Interval& d) {
                return Stereon::Interval{std::lerp(d.GetLower(), d.GetUpper(), 0.2),
                                         std::lerp(d.GetLower(), d.GetUpper(), 0.7)};
            };
            const Stereon::Box2 box{part(domain.X), part(domain.Y)};
            const Stereon::TrimmedSurface trimmed =
                *Stereon::MakeTrimmedSurface(Stereon::Surface(surface), box);
            const occ::handle<Geom_RectangularTrimmedSurface> occt =
                new Geom_RectangularTrimmedSurface(ToOcct(surface), box.X.GetLower(),
                                                   box.X.GetUpper(), box.Y.GetLower(),
                                                   box.Y.GetUpper());
            for (const double s : {0.0, 0.35, 1.0}) {
                for (const double t : {0.0, 0.6, 1.0}) {
                    CheckSecondOrder(trimmed, *occt,
                                     Point2{std::lerp(box.X.GetLower(), box.X.GetUpper(), s),
                                            std::lerp(box.Y.GetLower(), box.Y.GetUpper(), t)});
                }
            }
            CheckProjectionsAtLeastAsClose(trimmed, occt, f);
        }
    }
}

TEST_CASE("horn torus derivatives agree with OCCT", "[occt][surface][regression]")
{
    const Frame f = TiltedFrame();
    const Torus torus{f.Origin, 1.0, 1.0, f.X, f.Y};
    const occ::handle<Geom_ToroidalSurface> occt =
        new Geom_ToroidalSurface(ToOcctFrame(f), 1.0, 1.0);
    for (int i = 0; i < SampleCount; ++i) {
        CheckSecondOrder(torus, *occt, Point2{Sample(i), Sample(SampleCount - 1 - i)});
    }
}
