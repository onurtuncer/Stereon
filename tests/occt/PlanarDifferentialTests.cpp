// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests of the plane curve intersections (v0.2, every pair of
// lines, circles, ellipses, hyperbolas, parabolas and B-splines, rational
// ones included) against OpenCASCADE's Geom2dAPI_InterCurveCurve, on the
// same curves over the same parameter domains: the same number of points,
// each of Stereon's within reach of one of OCCT's and each of OCCT's within
// reach of one of Stereon's. The curves are placed so that every meeting is
// a transverse crossing away from the domains' ends, where neither kernel's
// tolerance can change the count.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_Hyperbola.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_Parabola.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Ax22d.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>

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

using Stereon::BSplineCurve2d;
using Stereon::Circle2d;
using Stereon::Curve2d;
using Stereon::Ellipse2d;
using Stereon::Hyperbola2d;
using Stereon::Interval;
using Stereon::Line2d;
using Stereon::Parabola2d;
using Stereon::Point2;
using Stereon::Vec2;

namespace {

const Stereon::Context ctx = Stereon::Context::Default();
constexpr double pi = std::numbers::pi;

gp_Pnt2d ToOcct(const Point2& p)
{
    return {p.X, p.Y};
}

gp_Dir2d ToDir(const Vec2& v)
{
    return {v.X, v.Y};
}

/// The frame of a conic: its centre or vertex and its axes, the sense from
/// the Y axis, as both kernels parametrise it.
gp_Ax22d FrameOf(const Point2& origin, const Vec2& x, const Vec2& y)
{
    return gp_Ax22d(ToOcct(origin), ToDir(x), ToDir(y));
}

occ::handle<Geom2d_Curve> ToOcct(const Line2d& l)
{
    // OCCT's line runs at unit speed: the test's lines have unit directions.
    return new Geom2d_Line(ToOcct(l.Origin), ToDir(l.Direction));
}
occ::handle<Geom2d_Curve> ToOcct(const Circle2d& c)
{
    return new Geom2d_Circle(FrameOf(c.Center, c.XAxis, c.YAxis), c.Radius);
}
occ::handle<Geom2d_Curve> ToOcct(const Ellipse2d& e)
{
    return new Geom2d_Ellipse(FrameOf(e.Center, e.XAxis, e.YAxis), e.MajorRadius, e.MinorRadius);
}
occ::handle<Geom2d_Curve> ToOcct(const Hyperbola2d& h)
{
    return new Geom2d_Hyperbola(FrameOf(h.Center, h.XAxis, h.YAxis), h.MajorRadius, h.MinorRadius);
}
occ::handle<Geom2d_Curve> ToOcct(const Parabola2d& p)
{
    return new Geom2d_Parabola(FrameOf(p.Vertex, p.XAxis, p.YAxis), p.FocalLength);
}
occ::handle<Geom2d_Curve> ToOcct(const BSplineCurve2d& c)
{
    const int count = static_cast<int>(c.Poles.size());
    NCollection_Array1<gp_Pnt2d> poles(1, count);
    NCollection_Array1<double> weights(1, count);
    for (int i = 0; i < count; ++i) {
        const auto k = static_cast<std::size_t>(i);
        poles(i + 1) = ToOcct(c.Poles[k]);
        weights(i + 1) = c.Weights.empty() ? 1.0 : c.Weights[k];
    }
    std::vector<double> values;
    std::vector<int> multiplicities;
    for (const double u : c.Knots) {
        if (values.empty() || values.back() != u) {
            values.push_back(u);
            multiplicities.push_back(0);
        }
        ++multiplicities.back();
    }
    const int n = static_cast<int>(values.size());
    NCollection_Array1<double> knots(1, n);
    NCollection_Array1<int> mults(1, n);
    for (int i = 0; i < n; ++i) {
        knots(i + 1) = values[static_cast<std::size_t>(i)];
        mults(i + 1) = multiplicities[static_cast<std::size_t>(i)];
    }
    return new Geom2d_BSplineCurve(poles, weights, knots, mults, c.Degree);
}

/// A test curve: Stereon's, OCCT's, and the domain both are asked over.
struct TestCurve {
    std::string Name;
    Curve2d Ours;
    occ::handle<Geom2d_Curve> Theirs;
    Interval Domain;
};

template <typename T>
TestCurve Make(std::string name, const T& curve, const Interval& domain)
{
    return {std::move(name), Curve2d{curve}, ToOcct(curve), domain};
}

BSplineCurve2d Spline(const int degree, std::vector<double> knots, std::vector<Point2> poles,
                      std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

/// Seven curves about the origin, crossing one another transversally.
std::vector<TestCurve> Curves()
{
    const double c = std::cos(0.4);
    const double s = std::sin(0.4);
    const Vec2 ex{1.0, 0.0};
    const Vec2 ey{0.0, 1.0};
    const Vec2 slope = Normalized(Vec2{1.0, 0.3});
    const double w = std::sqrt(0.5);
    return {
        Make("line", Line2d{Point2{-5.0, -1.2}, slope}, Interval{0.0, 11.0}),
        Make("circle", Circle2d{Point2{0.3, 0.2}, 2.0, ex, ey}, Interval{0.0, 2.0 * pi}),
        Make("ellipse", Ellipse2d{Point2{-0.4, 0.1}, 3.0, 1.2, Vec2{c, s}, Vec2{-s, c}},
             Interval{0.0, 2.0 * pi}),
        Make("hyperbola", Hyperbola2d{Point2{-1.5, 0.2}, 1.0, 0.8, ex, ey}, Interval{-2.0, 2.0}),
        Make("parabola", Parabola2d{Point2{-2.2, -0.5}, 0.6, ex, ey}, Interval{-3.0, 3.0}),
        Make("cubic B-spline",
             Spline(3, {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
                    {Point2{-3.0, -2.0}, Point2{-1.0, 3.0}, Point2{1.0, -3.0}, Point2{3.0, 2.0},
                     Point2{4.0, -1.0}}),
             Interval{0.0, 1.0}),
        Make("rational arc",
             Spline(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                    {Point2{3.0, -0.3}, Point2{3.0, 2.2}, Point2{0.5, 2.2}}, {1.0, w, 1.0}),
             Interval{0.0, 1.0}),
    };
}

/// OCCT's points of the pair over the domains, from the curves trimmed to
/// them.
std::vector<Point2> OcctPoints(const TestCurve& a, const TestCurve& b)
{
    const occ::handle<Geom2d_Curve> ta =
        new Geom2d_TrimmedCurve(a.Theirs, a.Domain.GetLower(), a.Domain.GetUpper());
    const occ::handle<Geom2d_Curve> tb =
        new Geom2d_TrimmedCurve(b.Theirs, b.Domain.GetLower(), b.Domain.GetUpper());
    const Geom2dAPI_InterCurveCurve occt(ta, tb, 1e-9);
    REQUIRE(occt.NbSegments() == 0);
    std::vector<Point2> points;
    for (int i = 1; i <= occt.NbPoints(); ++i) {
        const gp_Pnt2d p = occt.Point(i);
        points.push_back(Point2{p.X(), p.Y()});
    }
    return points;
}

double NearestTo(const std::vector<Point2>& points, const Point2& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Point2& q : points) {
        nearest = std::min(nearest, Distance(p, q));
    }
    return nearest;
}

} // namespace

TEST_CASE("every plane curve pair agrees with Geom2dAPI_InterCurveCurve", "[occt][planar]")
{
    const std::vector<TestCurve> curves = Curves();
    for (const TestCurve& a : curves) {
        for (const TestCurve& b : curves) {
            if (&a == &b) {
                continue;
            }
            INFO(a.Name << " x " << b.Name);
            const auto ours = Intersect(a.Ours, a.Domain, b.Ours, b.Domain, ctx);
            INFO((ours ? std::string() : std::string(ours.error().GetMessage())));
            REQUIRE(ours.has_value());
            REQUIRE(ours->Overlaps.empty());
            const std::vector<Point2> theirs = OcctPoints(a, b);
            std::vector<Point2> mine;
            for (const Stereon::PlanarPoint& p : ours->Points) {
                REQUIRE(p.Contact == Stereon::ContactType::Crossing);
                mine.push_back(p.Position);
            }
            REQUIRE(mine.size() == theirs.size());
            for (const Point2& p : mine) {
                REQUIRE(NearestTo(theirs, p) <= 1e-7);
            }
            for (const Point2& p : theirs) {
                REQUIRE(NearestTo(mine, p) <= 1e-7);
            }
        }
    }
}
