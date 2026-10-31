// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests of the analytic surface intersections (the Boolean
// prototype's M3) against OpenCASCADE's IntAna_QuadQuadGeo, the same closed
// forms by another hand: for each plane–plane, plane–cylinder and parallel
// cylinder–cylinder configuration the two kernels must agree on the kind of
// result and, where there are curves, every point of Stereon's branches
// must lie on one of OCCT's. Cylinders with skew axes (M3b) are compared
// with IntAna_IntQuadQuad's traces of the quartic, both ways: every sample
// of Stereon's branches lies on OCCT's curves and every sample of OCCT's
// curves lies on Stereon's branches.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <IntAna_Curve.hxx>
#include <IntAna_IntConicQuad.hxx>
#include <IntAna_IntQuadQuad.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_Quadric.hxx>
#include <IntAna_ResultType.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Sphere.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "OcctConversion.hpp"

using Stereon::Curve;
using Stereon::Cylinder;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

using StereonTest::ToOcct;

const Stereon::Context ctx = Stereon::Context::Default();

gp_Dir ToDir(const Vec3& v)
{
    return gp_Dir(v.X, v.Y, v.Z);
}

gp_Pln ToOcct(const Plane& plane)
{
    return gp_Pln(
        gp_Ax3(ToOcct(plane.Origin), ToDir(Cross(plane.XAxis, plane.YAxis)), ToDir(plane.XAxis)));
}

gp_Cylinder ToOcct(const Cylinder& cylinder)
{
    return gp_Cylinder(
        gp_Ax3(ToOcct(cylinder.Origin), ToDir(Axis(cylinder)), ToDir(cylinder.XAxis)),
        cylinder.Radius);
}

/// The distance from p to the nearest of OCCT's result curves.
gp_Lin ToOcct(const Line& line)
{
    return gp_Lin(ToOcct(line.Origin), ToDir(line.Direction));
}

gp_Sphere ToOcct(const Sphere& sphere)
{
    return gp_Sphere(gp_Ax3(ToOcct(sphere.Center), ToDir(Axis(sphere)), ToDir(sphere.XAxis)),
                     sphere.Radius);
}

/// Stereon's line against the surface, and OCCT's IntAna_IntConicQuad on
/// the same line and quadric: the same count, every point within `within`
/// of one of OCCT's, and the parameter along the unit direction the same.
template <typename S>
void RequireAgreesWithConicQuad(const Line& line, const S& surface, const double within)
{
    const auto ours = Intersect(Curve{line}, Surface{surface}, ctx);
    REQUIRE(ours.has_value());
    const IntAna_Quadric quadric(ToOcct(surface));
    const IntAna_IntConicQuad occt(ToOcct(line), quadric);
    REQUIRE(occt.IsDone());
    REQUIRE_FALSE(occt.IsInQuadric());
    REQUIRE(ours->Overlaps.empty());
    REQUIRE(static_cast<int>(ours->Points.size()) == occt.NbPoints());
    const double speed = Length(line.Direction);
    for (const Stereon::CurveSurfacePoint& p : ours->Points) {
        double nearest = std::numeric_limits<double>::infinity();
        for (int i = 1; i <= occt.NbPoints(); ++i) {
            const double apart = occt.Point(i).Distance(ToOcct(p.Position));
            if (apart < nearest) {
                nearest = apart;
                REQUIRE(std::abs(p.Parameter * speed - occt.ParamOnConic(i)) <= within);
            }
        }
        REQUIRE(nearest <= within);
    }
}

double DistanceToResult(const IntAna_QuadQuadGeo& occt, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (int i = 1; i <= occt.NbSolutions(); ++i) {
        switch (occt.TypeInter()) {
        case IntAna_Line:
            nearest = std::min(nearest, occt.Line(i).Distance(ToOcct(p)));
            break;
        case IntAna_Circle:
            nearest = std::min(nearest, occt.Circle(i).Distance(ToOcct(p)));
            break;
        case IntAna_Ellipse: {
            GeomAPI_ProjectPointOnCurve projection(ToOcct(p), new Geom_Ellipse(occt.Ellipse(i)));
            nearest = std::min(nearest, projection.LowerDistance());
            break;
        }
        default:
            break;
        }
    }
    return nearest;
}

/// Every sample of every branch lies on one of OCCT's curves.
void RequireOnOcct(const Stereon::SurfaceSurfaceIntersection& result,
                   const IntAna_QuadQuadGeo& occt, const double tolerance)
{
    REQUIRE(occt.IsDone());
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const double lower = std::isinf(branch.Range.GetLower()) ? -7.0 : branch.Range.GetLower();
        const double upper = std::isinf(branch.Range.GetUpper()) ? 7.0 : branch.Range.GetUpper();
        for (int k = 0; k <= 40; ++k) {
            const double t = lower + (upper - lower) * k / 40.0;
            REQUIRE(DistanceToResult(occt, Eval(branch.Curve3d, t)) <= tolerance);
        }
    }
}

const Plane xy{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Cylinder aboutZ{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
const Cylinder aboutX{Point3{}, 5.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true};

Plane Tilted(const Point3& origin, const Vec3& normal)
{
    const Vec3 n = Normalized(normal);
    const Vec3 seed = std::abs(n.X) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
    const Vec3 x = Normalized(Cross(seed, n));
    return Plane{origin, x, Cross(n, x)};
}

} // namespace

TEST_CASE("lines against planes, spheres and cylinders agree with IntAna_IntConicQuad",
          "[occt][intersect]")
{
    // Generic lines, away from tangency and from parallel, where OCCT's
    // tolerances cannot change the count; the surfaces off the axes.
    const Plane plane{Point3{1.0, -2.0, 0.5}, Normalized(Vec3{1.0, 1.0, 0.0}), Vec3{0.0, 0.0, 1.0}};
    const Sphere sphere{Point3{2.0, 1.0, -1.0}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Cylinder cylinder{Point3{-1.0, 0.5, 2.0}, 2.0, Normalized(Vec3{1.0, 1.0, 0.0}),
                            Vec3{0.0, 0.0, 1.0}};
    for (int k = 0; k < 24; ++k) {
        const double a = static_cast<double>(k) / 24.0 * 2.0 * std::numbers::pi;
        // Through a point near each surface's centre, so the quadrics are
        // crossed twice and the plane once.
        const Vec3 direction{std::cos(a) + 0.3, std::sin(a) - 0.2, 0.7 - 0.05 * k};
        const Point3 through{1.5 + 0.25 * std::sin(3.0 * a), 0.5 + 0.25 * std::cos(2.0 * a),
                             0.5 + 0.1 * k / 24.0};
        const Line line{through - 4.0 * direction, direction};
        RequireAgreesWithConicQuad(line, plane, 1e-9);
        RequireAgreesWithConicQuad(line, sphere, 1e-9);
        RequireAgreesWithConicQuad(line, cylinder, 1e-9);
    }
}

TEST_CASE("plane-plane intersections agree with OCCT", "[occt][intersect]")
{
    const Plane tilted = Tilted(Point3{1.0, 2.0, 3.0}, Vec3{1.0, 1.0, 2.0});
    const auto line = Intersect(Surface(xy), Surface(tilted), ctx);
    REQUIRE(line->Branches.size() == 1);
    const IntAna_QuadQuadGeo occtLine(ToOcct(xy), ToOcct(tilted), Precision::Angular(),
                                      Precision::Confusion());
    REQUIRE(occtLine.TypeInter() == IntAna_Line);
    RequireOnOcct(*line, occtLine, 1e-9);

    const Plane lifted{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Intersect(Surface(xy), Surface(lifted), ctx)->IsEmpty());
    REQUIRE(
        IntAna_QuadQuadGeo(ToOcct(xy), ToOcct(lifted), Precision::Angular(), Precision::Confusion())
            .TypeInter() == IntAna_Empty);

    const Plane same{Point3{3.0, 4.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{-1.0, 0.0, 0.0}};
    REQUIRE(Intersect(Surface(xy), Surface(same), ctx)->Overlaps.size() == 1);
    REQUIRE(
        IntAna_QuadQuadGeo(ToOcct(xy), ToOcct(same), Precision::Angular(), Precision::Confusion())
            .TypeInter() == IntAna_Same);
}

TEST_CASE("plane-cylinder intersections agree with OCCT", "[occt][intersect]")
{
    struct Case {
        Plane Cut;
        Cylinder Body;
        IntAna_ResultType Type;
        std::size_t Branches;
    };
    const std::vector<Case> cases{
        {Plane{Point3{0.0, 0.0, 3.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutZ,
         IntAna_Circle, 1},
        {Tilted(Point3{0.0, 1.0, 0.0}, Vec3{0.0, 1.0, 1.0}), aboutZ, IntAna_Ellipse, 1},
        {Tilted(Point3{2.0, 0.0, 1.0}, Vec3{1.0, 0.0, 3.0}), aboutZ, IntAna_Ellipse, 1},
        {Plane{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Line, 2},
        {Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Line, 1},
        {Plane{Point3{0.0, 0.0, 6.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Empty, 0},
    };
    for (const Case& c : cases) {
        const auto result = Intersect(Surface(c.Cut), Surface(c.Body), ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.size() == c.Branches);
        const IntAna_QuadQuadGeo occt(ToOcct(c.Cut), ToOcct(c.Body), Precision::Angular(),
                                      Precision::Confusion());
        REQUIRE(occt.TypeInter() == c.Type);
        REQUIRE(static_cast<std::size_t>(std::max(occt.NbSolutions(), 0)) == c.Branches);
        RequireOnOcct(*result, occt, 1e-9);
    }
}

TEST_CASE("parallel cylinder-cylinder intersections agree with OCCT", "[occt][intersect]")
{
    struct Case {
        Cylinder Other;
        IntAna_ResultType Type;
        std::size_t Branches;
    };
    const Cylinder a{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
    const std::vector<Case> cases{
        {Cylinder{Point3{5.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 2},
        {Cylinder{Point3{6.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 1},
        {Cylinder{Point3{1.0, 0.0, 0.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 1},
        {Cylinder{Point3{10.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Empty, 0},
        {Cylinder{Point3{0.0, 0.0, 9.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}, true},
         IntAna_Same, 0},
    };
    for (const Case& c : cases) {
        const auto result = Intersect(Surface(a), Surface(c.Other), ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.size() == c.Branches);
        REQUIRE(result->Overlaps.empty() == (c.Type != IntAna_Same));
        const IntAna_QuadQuadGeo occt(ToOcct(a), ToOcct(c.Other), Precision::Confusion());
        REQUIRE(occt.TypeInter() == c.Type);
        RequireOnOcct(*result, occt, 1e-9);
    }
}

namespace {

/// The nearest of the points to p.
double DistanceToNearest(const std::vector<Point3>& points, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Point3& q : points) {
        nearest = std::min(nearest, Distance(p, q));
    }
    return nearest;
}

/// Samples of every branch of Stereon's result.
std::vector<Point3> SamplesOf(const Stereon::SurfaceSurfaceIntersection& result, const int count)
{
    std::vector<Point3> points;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        for (int k = 0; k <= count; ++k) {
            const double t = branch.Range.GetLower() +
                             (branch.Range.GetUpper() - branch.Range.GetLower()) * k / count;
            points.push_back(Eval(branch.Curve3d, t));
        }
    }
    return points;
}

/// Samples of every curve OCCT traces for the quartic of two cylinders,
/// with the largest gap between neighbours along a curve: a point within
/// that gap of the nearest sample is on the curve to that accuracy.
struct OcctTrace {
    std::vector<Point3> Points;
    double Gap = 0.0;
};

OcctTrace OcctQuarticSamples(const Cylinder& a, const Cylinder& b, const int count)
{
    const IntAna_IntQuadQuad occt(ToOcct(a), IntAna_Quadric(ToOcct(b)), Precision::Confusion());
    REQUIRE(occt.IsDone());
    REQUIRE_FALSE(occt.IdenticalElements());
    OcctTrace trace;
    for (int i = 1; i <= occt.NbCurve(); ++i) {
        IntAna_Curve curve = occt.Curve(i);
        double first = 0.0;
        double last = 0.0;
        curve.Domain(first, last);
        for (int k = 0; k <= count; ++k) {
            const gp_Pnt p = curve.Value(first + (last - first) * k / count);
            const Point3 point{p.X(), p.Y(), p.Z()};
            if (k > 0) {
                trace.Gap = std::max(trace.Gap, Distance(point, trace.Points.back()));
            }
            trace.Points.push_back(point);
        }
    }
    return trace;
}

/// Refine the reference sampling itself until its spatial gap is small
/// enough. Uniform parameter spacing need not imply uniform spatial spacing,
/// particularly next to a quartic's folds. Keep a finite work budget.
OcctTrace RefinedOcctQuarticSamples(const Cylinder& a, const Cylinder& b)
{
    constexpr int InitialCount = 4000;
    constexpr int MaxCount = 64000;
    constexpr double MaxGap = 0.1;
    for (int count = InitialCount;; count *= 2) {
        OcctTrace trace = OcctQuarticSamples(a, b, count);
        REQUIRE_FALSE(trace.Points.empty());
        if (trace.Gap < MaxGap) {
            return trace;
        }
        INFO("OCCT reference sampling count " << count << ", gap " << trace.Gap);
        REQUIRE(count < MaxCount);
    }
}

/// The distance from p to the nearest of Stereon's branches, by the global
/// projection onto their B-spline curves.
double DistanceToBranches(const Stereon::SurfaceSurfaceIntersection& result, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const auto* curve = branch.Curve3d.Target<Stereon::BSplineCurve>();
        REQUIRE(curve != nullptr);
        nearest = std::min(nearest, Distance(p, Eval(*curve, Project(*curve, p))));
    }
    return nearest;
}

} // namespace

TEST_CASE("skew cylinder-cylinder intersections agree with OCCT", "[occt][intersect]")
{
    const Cylinder a{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
    const auto alongX = [](const Point3& through, const double radius) {
        return Cylinder{through, radius, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true};
    };
    const std::vector<Cylinder> others{
        alongX(Point3{0.0, 4.0, 0.0}, 2.0), // one loop
        alongX(Point3{0.0, 1.0, 0.0}, 1.0), // two loops, B through A
        alongX(Point3{}, 2.0),              // meeting axes, two loops
        Cylinder{Point3{1.0, 2.0, 0.5}, 2.5, Vec3{0.0, 0.6, 0.8}, Vec3{0.0, -0.8, 0.6}, true},
    };
    for (const Cylinder& other : others) {
        for (const auto& [first, second] : {std::pair{a, other}, std::pair{other, a}}) {
            const auto result = Intersect(Surface(first), Surface(second), ctx);
            REQUIRE(result.has_value());
            REQUIRE_FALSE(result->Branches.empty());
            // Stereon's samples lie within a sample gap of OCCT's traces,
            // and OCCT's samples lie on Stereon's branches to the certified
            // bound (within the resolution) plus OCCT's own accuracy.
            const OcctTrace occt = RefinedOcctQuarticSamples(first, second);
            REQUIRE_FALSE(occt.Points.empty());
            REQUIRE(occt.Gap < 0.1);
            for (const Point3& p : SamplesOf(*result, 40)) {
                REQUIRE(DistanceToNearest(occt.Points, p) <= occt.Gap + 1e-6);
            }
            for (std::size_t i = 0; i < occt.Points.size(); i += 50) {
                REQUIRE(DistanceToBranches(*result, occt.Points[i]) <= 1e-4);
            }
        }
    }
}

TEST_CASE("the tangent point and the two ellipses of skew cylinders agree with OCCT",
          "[occt][intersect]")
{
    const Cylinder a{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
    // Externally tangent: axes the sum of the radii apart.
    const Cylinder touching{Point3{0.0, 5.0, 0.0}, 2.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0},
                            true};
    const auto point = Intersect(Surface(a), Surface(touching), ctx);
    REQUIRE(point.has_value());
    REQUIRE(point->Points.size() == 1);
    const IntAna_QuadQuadGeo occtPoint(ToOcct(a), ToOcct(touching), Precision::Confusion());
    REQUIRE(occtPoint.IsDone());
    REQUIRE(occtPoint.TypeInter() == IntAna_Point);
    REQUIRE(occtPoint.Point(1).Distance(ToOcct(point->Points.front().Position)) <= 1e-9);

    // Meeting axes and equal radii: two ellipses.
    const Cylinder crossing{Point3{}, 3.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true};
    const auto ellipses = Intersect(Surface(a), Surface(crossing), ctx);
    REQUIRE(ellipses.has_value());
    const IntAna_QuadQuadGeo occtEllipses(ToOcct(a), ToOcct(crossing), Precision::Confusion());
    REQUIRE(occtEllipses.TypeInter() == IntAna_Ellipse);
    REQUIRE(occtEllipses.NbSolutions() == 2);
    RequireOnOcct(*ellipses, occtEllipses, 1e-4);
}
