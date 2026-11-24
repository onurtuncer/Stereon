// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests of the curve-surface handlers (v0.2 ledger,
// curve-surface) against OpenCASCADE's GeomAPI_IntCS: every curve kind
// beyond the line (B-spline curves, polynomial and rational, circles,
// ellipses, hyperbolas and parabolas) against every surface kind OCCT shares
// with Stereon (planes, cylinders, cones, spheres, tori and B-spline
// surfaces, polynomial and rational), on the same curves over the same
// parameter domains. Lines are compared in IntersectionDifferentialTests.cpp
// and circles and ellipses against the quadrics with IntAna_IntConicQuad
// there; the B-spline handlers and the conics against the torus and against
// B-spline surfaces are compared only here.
//
// Each curve is placed in eight tilted frames about a point of the surface,
// where it crosses the surface transversally and away from the domains'
// ends, so that neither kernel's tolerance can change the count. The kernels
// must find the same points: as many of Stereon's as distinct points of
// OCCT's (OCCT reports a crossing on a closed surface's seam once per side,
// Stereon once), each of Stereon's within 1e-7 of one of OCCT's and each of
// OCCT's within 1e-7 of one of Stereon's. A configuration Stereon fails
// closed on (IntersectionFailed, or a point it reports undecided or as a
// touch) is not a disagreement, but it is counted, and at most one in
// twenty may be.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <GeomAPI_IntCS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_Parabola.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "OcctConversion.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::Curve;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Parabola;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

using StereonTest::ToOcct;

const Stereon::Context ctx = Stereon::Context::Default();
constexpr double pi = std::numbers::pi;

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// How close a point of one kernel must be to one of the other's.
constexpr double Reach = 1e-7;

gp_Dir ToDir(const Vec3& v)
{
    return gp_Dir(v.X, v.Y, v.Z);
}

gp_Ax3 FrameOf(const Point3& origin, const Vec3& axis, const Vec3& x)
{
    return gp_Ax3(ToOcct(origin), ToDir(axis), ToDir(x));
}

gp_Ax2 ConicFrame(const Point3& center, const Vec3& x, const Vec3& y)
{
    return gp_Ax2(ToOcct(center), ToDir(Cross(x, y)), ToDir(x));
}

// --- Surfaces -----------------------------------------------------------------

occ::handle<Geom_Surface> GeomOf(const Plane& p)
{
    return new Geom_Plane(FrameOf(p.Origin, Cross(p.XAxis, p.YAxis), p.XAxis));
}
occ::handle<Geom_Surface> GeomOf(const Cylinder& c)
{
    return new Geom_CylindricalSurface(FrameOf(c.Origin, Axis(c), c.XAxis), c.Radius);
}
occ::handle<Geom_Surface> GeomOf(const Cone& c)
{
    return new Geom_ConicalSurface(FrameOf(c.Origin, Axis(c), c.XAxis), c.SemiAngle, c.Radius);
}
occ::handle<Geom_Surface> GeomOf(const Sphere& s)
{
    return new Geom_SphericalSurface(FrameOf(s.Center, Axis(s), s.XAxis), s.Radius);
}
occ::handle<Geom_Surface> GeomOf(const Torus& t)
{
    return new Geom_ToroidalSurface(FrameOf(t.Center, Axis(t), t.XAxis), t.MajorRadius,
                                    t.MinorRadius);
}
occ::handle<Geom_Surface> GeomOf(const BSplineSurface& s)
{
    return ToOcct(s);
}

/// A surface in both kernels, a point of it the test curves are placed
/// about, and the size they are scaled to.
struct TestSurface {
    std::string Name;
    Surface Ours;
    occ::handle<Geom_Surface> Theirs;
    Point3 Anchor;
    double Scale;
};

template <typename S>
TestSurface MakeSurface(std::string name, const S& s, const Point3& anchor, const double scale)
{
    return {std::move(name), Surface{s}, GeomOf(s), anchor, scale};
}

/// A bicubic patch over [0, 4] x [0, 3], wavy in z.
BSplineSurface WavyPatch()
{
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<Point3> poles;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double z = 0.4 * std::sin(1.3 * i + 0.7 * j) - 0.1 * i * j;
            poles.push_back({i / 3.0 * 4.0, j / 3.0 * 3.0, z});
        }
    }
    auto made = Stereon::MakeBSplineSurface(3, 3, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

std::vector<TestSurface> Surfaces()
{
    const Cylinder round{Point3{1.0, -2.0, 0.5}, 2.0, ex, ey};
    const auto rational =
        Stereon::ToBSpline(round, Stereon::Box2{Interval{0.0, 2.0 * pi}, Interval{-20.0, 20.0}});
    REQUIRE(rational.has_value());
    return {
        MakeSurface("plane", Plane{Point3{1.0, -0.5, 0.5}, Normalized(Vec3{1.0, 1.0, 0.0}), ez},
                    Point3{1.0, -0.5, 0.5}, 1.0),
        MakeSurface("cylinder",
                    Cylinder{Point3{-0.5, 0.5, 1.0}, 1.25, Normalized(Vec3{1.0, 1.0, 0.0}), ez},
                    Point3{-0.5, 0.5, 1.0}, 1.0),
        MakeSurface("cone", Cone{Point3{0.5, 0.0, -0.5}, 0.75, 0.4, ex, ey, true},
                    Point3{0.5, 0.0, 0.0}, 0.8),
        MakeSurface("sphere", Sphere{Point3{1.0, 0.5, -0.5}, 1.5, ex, ey, true},
                    Point3{1.0, 0.5, -0.5}, 1.0),
        MakeSurface("torus", Torus{Point3{0.3, -0.2, 0.1}, 2.0, 0.6, ex, ey, true},
                    Point3{0.3, -0.2, 0.1}, 1.2),
        MakeSurface("bicubic patch", WavyPatch(), Point3{2.0, 1.5, 0.0}, 0.55),
        MakeSurface("rational cylinder", *rational, Point3{1.0, -2.0, 0.5}, 1.4),
    };
}

// --- Curves -------------------------------------------------------------------

/// A placement: a point and an orthonormal frame, x and y spanning the
/// plane the conics lie in, z out of it; and a scale.
struct Placement {
    Point3 At;
    Vec3 X;
    Vec3 Y;
    Vec3 Z;
    double Scale;

    [[nodiscard]] Point3 Map(const double lx, const double ly, const double lz) const
    {
        return At + Scale * lx * X + Scale * ly * Y + Scale * lz * Z;
    }
};

/// The k-th of eight frames turned about a tilted axis, near the anchor.
Placement PlacementOf(const TestSurface& s, const int k)
{
    const double a = static_cast<double>(k) / 8.0 * 2.0 * pi + 0.3;
    const Vec3 x = Normalized(Vec3{std::cos(a), std::sin(a), 0.3});
    const Vec3 y = Normalized(Cross(Normalized(Vec3{0.2, -0.3 + 0.1 * k, 1.0}), x));
    const Vec3 offset{0.1 * std::sin(2.0 * a), 0.12 * std::cos(a), 0.05 * k / 8.0};
    return {s.Anchor + s.Scale * offset, x, y, Cross(x, y), s.Scale};
}

/// A test curve: Stereon's, OCCT's over the domain, and the domain.
struct TestCurve {
    std::string Name;
    Curve Ours;
    occ::handle<Geom_Curve> Theirs;
    Interval Domain;
};

TestCurve Trimmed(std::string name, Curve ours, const occ::handle<Geom_Curve>& theirs,
                  const Interval& domain)
{
    return {std::move(name), std::move(ours),
            new Geom_TrimmedCurve(theirs, domain.GetLower(), domain.GetUpper()), domain};
}

BSplineCurve Spline(const int degree, std::vector<double> knots, std::vector<Point3> poles,
                    std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

TestCurve SplineCurve(std::string name, const BSplineCurve& c)
{
    return {std::move(name), Curve{c}, ToOcct(c), Domain(c)};
}

/// The six curves in placement p: a cubic wave out of the plane, a rational
/// quadratic quarter circle, a circle, an ellipse, a hyperbola and a
/// parabola, each over a bounded domain.
std::vector<TestCurve> Curves(const Placement& p)
{
    const double w = std::sqrt(0.5);
    const BSplineCurve cubic =
        Spline(3, {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
               {p.Map(-1.8, -0.6, 0.2), p.Map(-0.6, 1.2, -0.3), p.Map(0.6, -1.2, 0.25),
                p.Map(1.8, 0.6, -0.2), p.Map(2.1, -0.3, 0.1)});
    const BSplineCurve quarter =
        Spline(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
               {p.Map(1.7, -0.2, 0.0), p.Map(1.7, 1.5, 0.0), p.Map(0.0, 1.5, 0.0)}, {1.0, w, 1.0});
    const Circle circle{p.At, 1.75 * p.Scale, p.X, p.Y};
    const Ellipse ellipse{p.At, 2.25 * p.Scale, 1.0 * p.Scale, p.X, p.Y};
    const Hyperbola hyperbola{p.Map(-0.6, 0.0, 0.0), 0.6 * p.Scale, 0.5 * p.Scale, p.X, p.Y};
    const Parabola parabola{p.Map(-1.0, 0.0, 0.0), 0.4 * p.Scale, p.X, p.Y};
    return {
        SplineCurve("cubic B-spline", cubic),
        SplineCurve("rational quarter", quarter),
        Trimmed(
            "circle", Curve{circle},
            new Geom_Circle(ConicFrame(circle.Center, circle.XAxis, circle.YAxis), circle.Radius),
            Interval{0.0, 2.0 * pi}),
        Trimmed("ellipse", Curve{ellipse},
                new Geom_Ellipse(ConicFrame(ellipse.Center, ellipse.XAxis, ellipse.YAxis),
                                 ellipse.MajorRadius, ellipse.MinorRadius),
                Interval{0.0, 2.0 * pi}),
        Trimmed("hyperbola", Curve{hyperbola},
                new Geom_Hyperbola(ConicFrame(hyperbola.Center, hyperbola.XAxis, hyperbola.YAxis),
                                   hyperbola.MajorRadius, hyperbola.MinorRadius),
                Interval{-1.5, 1.5}),
        Trimmed("parabola", Curve{parabola},
                new Geom_Parabola(ConicFrame(parabola.Vertex, parabola.XAxis, parabola.YAxis),
                                  parabola.FocalLength),
                Interval{-1.5 * p.Scale, 1.5 * p.Scale}),
    };
}

// --- The comparison -----------------------------------------------------------

/// OCCT's points of the curve on the surface, those a rounding apart on a
/// closed surface's seam counted once.
std::vector<gp_Pnt> OcctPoints(const TestCurve& c, const TestSurface& s)
{
    GeomAPI_IntCS occt(c.Theirs, s.Theirs);
    REQUIRE(occt.IsDone());
    REQUIRE(occt.NbSegments() == 0);
    std::vector<gp_Pnt> distinct;
    for (int i = 1; i <= occt.NbPoints(); ++i) {
        if (std::ranges::none_of(
                distinct, [&](const gp_Pnt& q) { return q.Distance(occt.Point(i)) <= Reach; })) {
            distinct.push_back(occt.Point(i));
        }
    }
    return distinct;
}

double Nearest(const std::vector<gp_Pnt>& points, const gp_Pnt& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const gp_Pnt& q : points) {
        nearest = std::min(nearest, q.Distance(p));
    }
    return nearest;
}

/// The comparison's tally: pairs compared, and pairs Stereon failed closed
/// on.
struct Tally {
    int Compared = 0;
    int Closed = 0;
};

void Compare(const TestCurve& c, const TestSurface& s, Tally& tally)
{
    INFO(c.Name << " x " << s.Name);
    const auto ours = Intersect(c.Ours, c.Domain, s.Ours, Domain(s.Ours), ctx);
    if (!ours) {
        INFO(std::string(ours.error().GetMessage()));
        REQUIRE(ours.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
        ++tally.Closed;
        return;
    }
    REQUIRE(ours->Overlaps.empty());
    if (std::ranges::any_of(ours->Points, [](const Stereon::CurveSurfacePoint& p) {
            return !p.Decided || p.Contact != Stereon::ContactType::Crossing;
        })) {
        ++tally.Closed;
        return;
    }
    const std::vector<gp_Pnt> theirs = OcctPoints(c, s);
    std::vector<gp_Pnt> mine;
    for (const Stereon::CurveSurfacePoint& p : ours->Points) {
        mine.push_back(ToOcct(p.Position));
    }
    REQUIRE(mine.size() == theirs.size());
    for (const gp_Pnt& p : mine) {
        REQUIRE(Nearest(theirs, p) <= Reach);
    }
    for (const gp_Pnt& p : theirs) {
        REQUIRE(Nearest(mine, p) <= Reach);
    }
    ++tally.Compared;
}

} // namespace

TEST_CASE("curves against surfaces agree with GeomAPI_IntCS", "[occt][intersect][curve-surface]")
{
    Tally tally;
    int meetings = 0;
    for (const TestSurface& s : Surfaces()) {
        for (int k = 0; k < 8; ++k) {
            INFO("placement " << k);
            for (const TestCurve& c : Curves(PlacementOf(s, k))) {
                const int before = tally.Compared;
                Compare(c, s, tally);
                if (tally.Compared > before) {
                    meetings += static_cast<int>(OcctPoints(c, s).size());
                }
            }
        }
    }
    INFO(tally.Compared << " pairs compared, " << tally.Closed << " failed closed, " << meetings
                        << " points");
    // Seven surfaces, six curves, eight placements: 336 pairs.
    REQUIRE(tally.Compared + tally.Closed == 336);
    REQUIRE(tally.Closed * 20 <= 336);
    // The placements cross: the comparison is not of empty results.
    REQUIRE(meetings >= 100);
}
