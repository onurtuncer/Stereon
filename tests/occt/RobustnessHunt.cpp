// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A hunt for inputs on which OCCT gives a wrong answer, or none, and
// Stereon a right one: adversarial geometry for projection, arc length and
// point at a length, each answer checked against an independent reference
// (exact values where they exist, dense sampling otherwise). The cases are
// hidden ("[.hunt]"): they report what they find rather than pass or fail,
// and the findings worth keeping become ordinary tests in
// OcctRobustnessTests.cpp. Run with
//
//   stereon_occt_tests "[.hunt]"

#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Standard_Failure.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <optional>
#include <vector>

import Stereon.Surface;

#include "OcctConversion.hpp"
#include "ReferenceLengths.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Cylinder;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Torus;
using Stereon::Vec3;
using StereonTest::PolylineLength;
using StereonTest::PolylineReference;
using StereonTest::ToOcct;

namespace {

/// A fixed linear congruential generator, so every hunt finds the same
/// cases on every machine.
class Random {
public:
    explicit Random(const std::uint64_t seed)
        : m_State(seed)
    {}
    /// Uniform in [low, high).
    double Next(const double low, const double high)
    {
        m_State = m_State * 6364136223846793005ULL + 1442695040888963407ULL;
        const double unit = static_cast<double>(m_State >> 11) * 0x1p-53;
        return low + (high - low) * unit;
    }

private:
    std::uint64_t m_State;
};

/// What a hunt found: how many cases, and how often each kernel was wrong
/// (farther or less accurate than the reference allows) or gave no answer.
struct Tally {
    int Cases = 0;
    int OcctWrong = 0;
    int OcctNoAnswer = 0;
    int StereonWrong = 0;
    double WorstOcct = 0.0;
    double WorstStereon = 0.0;
};

void Report(const char* name, const Tally& t)
{
    std::printf("HUNT %s: %d cases; OCCT wrong %d, no answer %d (worst %.3g); "
                "Stereon wrong %d (worst %.3g)\n",
                name, t.Cases, t.OcctWrong, t.OcctNoAnswer, t.WorstOcct, t.StereonWrong,
                t.WorstStereon);
    std::fflush(stdout);
}

/// OCCT's closest distance from p to the curve, empty when it finds none
/// or throws.
std::optional<double> OcctDistance(const Point3& p, const occ::handle<Geom_Curve>& c)
{
    try {
        GeomAPI_ProjectPointOnCurve extrema(ToOcct(p), c);
        if (extrema.NbPoints() == 0) {
            return std::nullopt;
        }
        return extrema.LowerDistance();
    } catch (const Standard_Failure&) {
        return std::nullopt;
    }
}

std::optional<double> OcctDistance(const Point3& p, const occ::handle<Geom_Surface>& s)
{
    try {
        GeomAPI_ProjectPointOnSurf extrema(ToOcct(p), s);
        if (extrema.NbPoints() == 0) {
            return std::nullopt;
        }
        return extrema.LowerDistance();
    } catch (const Standard_Failure&) {
        return std::nullopt;
    }
}

/// The smallest distance from p to `samples` + 1 evenly spaced points of
/// the curve, an upper bound on the true distance.
double SampledDistance(const BSplineCurve& c, const Point3& p, const int samples)
{
    const Interval d = Domain(c);
    double best = Distance(Eval(c, d.GetLower()), p);
    for (int i = 1; i <= samples; ++i) {
        const double t = d.GetLower() + (d.GetUpper() - d.GetLower()) * i / samples;
        best = std::min(best, Distance(Eval(c, std::min(t, d.GetUpper())), p));
    }
    return best;
}

/// Records one projection: a kernel is wrong when it is farther than the
/// sampled bound by more than `slack`.
void Score(Tally& t, const double stereon, const std::optional<double> occt, const double sampled,
           const double slack)
{
    ++t.Cases;
    if (stereon > sampled + slack) {
        ++t.StereonWrong;
        t.WorstStereon = std::max(t.WorstStereon, stereon - sampled);
    }
    if (!occt) {
        ++t.OcctNoAnswer;
    } else if (*occt > std::min(stereon, sampled) + slack) {
        ++t.OcctWrong;
        t.WorstOcct = std::max(t.WorstOcct, *occt - std::min(stereon, sampled));
    }
}

/// A clamped curve along the x axis with one pole raised to height h and
/// weighted w: a narrow spike, narrower as w grows.
BSplineCurve SpikeCurve(const int degree, const int count, const int spike, const double h,
                        const double w)
{
    BSplineCurve c{degree, {}, {}, {}};
    c.Knots.assign(static_cast<std::size_t>(degree) + 1, 0.0);
    for (int i = 1; i < count - degree; ++i) {
        c.Knots.push_back(static_cast<double>(i));
    }
    c.Knots.insert(c.Knots.end(), static_cast<std::size_t>(degree) + 1,
                   static_cast<double>(count - degree));
    for (int i = 0; i < count; ++i) {
        c.Poles.push_back({static_cast<double>(i), i == spike ? h : 0.0, 0.0});
        c.Weights.push_back(i == spike ? w : 1.0);
    }
    return c;
}

/// A clamped polyline (degree 1) through the points with the given knots:
/// its length between two parameters is exact arithmetic on the points.
BSplineCurve Polyline(const std::vector<Point3>& points, const std::vector<double>& breaks)
{
    BSplineCurve c{1, {}, points, {}};
    c.Knots.push_back(breaks.front());
    c.Knots.insert(c.Knots.end(), breaks.begin(), breaks.end());
    c.Knots.push_back(breaks.back());
    return c;
}

/// A random polyline with stalls (repeated points) and knots spread over
/// many orders of magnitude.
BSplineCurve RandomPolyline(Random& random)
{
    const int count = static_cast<int>(random.Next(4.0, 40.0));
    std::vector<Point3> points;
    std::vector<double> breaks{0.0};
    for (int i = 0; i < count; ++i) {
        if (i > 0 && random.Next(0.0, 1.0) < 0.2) {
            points.push_back(points.back()); // a stall: zero speed on a span
        } else {
            points.push_back(
                {random.Next(-10.0, 10.0), random.Next(-10.0, 10.0), random.Next(-10.0, 10.0)});
        }
        if (i > 0) {
            breaks.push_back(breaks.back() + std::pow(10.0, random.Next(-6.0, 2.0)));
        }
    }
    return Polyline(points, breaks);
}

/// Prints a polyline and a request so a finding can be reproduced.
void PrintCase(const char* what, const BSplineCurve& c, const double target, const double error)
{
    std::printf("HUNT CASE %s: error %.3g, target length %.17g\n  knots:", what, error, target);
    for (const double u : c.Knots) {
        std::printf(" %.17g", u);
    }
    std::printf("\n  points:");
    for (const Point3& p : c.Poles) {
        std::printf(" (%.17g, %.17g, %.17g)", p.X, p.Y, p.Z);
    }
    std::printf("\n");
}

/// The smallest distance from p to a (samples + 1)^2 grid of the surface,
/// an upper bound on the true distance.
double SampledDistance(const BSplineSurface& s, const Point3& p, const int samples)
{
    const Stereon::Box2 d = Domain(s);
    double best = Distance(Eval(s, Point2{d.X.GetLower(), d.Y.GetLower()}), p);
    for (int i = 0; i <= samples; ++i) {
        const double u = std::min(d.X.GetUpper(), d.X.GetLower() + Width(d.X) * i / samples);
        for (int j = 0; j <= samples; ++j) {
            const double v = std::min(d.Y.GetUpper(), d.Y.GetLower() + Width(d.Y) * j / samples);
            best = std::min(best, Distance(Eval(s, Point2{u, v}), p));
        }
    }
    return best;
}

/// A clamped bicubic sheet over an n x n net in the plane z = 0 with one
/// pole raised to height h: a narrow bump, a few knot spans wide.
BSplineSurface BumpSurface(const int n, const int bumpU, const int bumpV, const double h)
{
    constexpr int degree = 3;
    std::vector<double> knots(degree + 1, 0.0);
    for (int i = 1; i < n - degree; ++i) {
        knots.push_back(static_cast<double>(i));
    }
    knots.insert(knots.end(), degree + 1, static_cast<double>(n - degree));
    BSplineSurface s{degree, degree, knots, knots, {}, {}};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const bool bump = i == bumpU && j == bumpV;
            s.Poles.push_back({static_cast<double>(i), static_cast<double>(j), bump ? h : 0.0});
        }
    }
    return s;
}

/// A clamped cubic hairpin: out along y = +gap/2, round a tight turn, back
/// along y = -gap/2, so two branches lie close together.
BSplineCurve Hairpin(const double length, const double gap, const int count)
{
    BSplineCurve c{3, {}, {}, {}};
    for (int i = 0; i < count; ++i) {
        c.Poles.push_back({length * i / (count - 1), 0.5 * gap, 0.0});
    }
    for (int i = count - 1; i >= 0; --i) {
        c.Poles.push_back({length * i / (count - 1), -0.5 * gap, 0.0});
    }
    const int poles = static_cast<int>(c.Poles.size());
    c.Knots.assign(4, 0.0);
    for (int i = 1; i < poles - 3; ++i) {
        c.Knots.push_back(static_cast<double>(i));
    }
    c.Knots.insert(c.Knots.end(), 4, static_cast<double>(poles - 3));
    return c;
}

/// A smooth rational cubic whose weights span orders of magnitude, so its
/// speed does too.
BSplineCurve WildWeights(Random& random)
{
    constexpr int count = 9;
    BSplineCurve c{3, {}, {}, {}};
    c.Knots.assign(4, 0.0);
    for (int i = 1; i < count - 3; ++i) {
        c.Knots.push_back(static_cast<double>(i));
    }
    c.Knots.insert(c.Knots.end(), 4, static_cast<double>(count - 3));
    for (int i = 0; i < count; ++i) {
        c.Poles.push_back(
            {random.Next(-10.0, 10.0), random.Next(-10.0, 10.0), random.Next(-10.0, 10.0)});
        c.Weights.push_back(std::pow(10.0, random.Next(-3.0, 3.0)));
    }
    return c;
}

/// A tilted orthonormal frame with exact rational axes:
/// X = (1, 2, 2) / 3, Y = (2, 1, -2) / 3, Z = X x Y = (-2, 2, -1) / 3.
constexpr Point3 Origin{1.0, -2.0, 0.5};
constexpr Vec3 AxisX{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
constexpr Vec3 AxisY{2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0};
constexpr Vec3 AxisZ{-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0};

gp_Ax3 OcctFrame()
{
    return {ToOcct(Origin), gp_Dir(AxisZ.X, AxisZ.Y, AxisZ.Z), gp_Dir(AxisX.X, AxisX.Y, AxisX.Z)};
}

/// The point with frame coordinates (x, y, z).
Point3 InFrame(const double x, const double y, const double z)
{
    return Origin + x * AxisX + y * AxisY + z * AxisZ;
}

/// Records a projection judged against an exact distance: a kernel is
/// wrong when it is farther than the exact distance by more than slack.
void ScoreExact(Tally& t, const double stereon, const std::optional<double> occt,
                const double exact, const double slack)
{
    ++t.Cases;
    if (stereon > exact + slack) {
        ++t.StereonWrong;
        t.WorstStereon = std::max(t.WorstStereon, stereon - exact);
    }
    if (!occt) {
        ++t.OcctNoAnswer;
    } else if (*occt > exact + slack) {
        ++t.OcctWrong;
        t.WorstOcct = std::max(t.WorstOcct, *occt - exact);
    }
}

/// A point near the axis of the frame: at distance 10^e from it, e in
/// [-12, -1], at height z, in a random direction around it.
Point3 NearAxis(Random& random, const double z)
{
    const double r = std::pow(10.0, random.Next(-12.0, -1.0));
    const double a = random.Next(0.0, 2.0 * std::numbers::pi);
    return InFrame(r * std::cos(a), r * std::sin(a), z);
}

/// An angle within 10^e of the seam at 0, e in [-15, -3], on either side.
double NearSeam(Random& random)
{
    const double offset = std::pow(10.0, random.Next(-15.0, -3.0));
    return random.Next(0.0, 1.0) < 0.5 ? offset : 2.0 * std::numbers::pi - offset;
}

} // namespace

TEST_CASE("hunt: projection onto narrow spikes", "[.hunt][curve]")
{
    Tally t;
    Random random(1);
    for (int trial = 0; trial < 400; ++trial) {
        const int degree = static_cast<int>(random.Next(2.0, 6.0));
        const int count = static_cast<int>(random.Next(12.0, 40.0));
        const int spike = static_cast<int>(random.Next(degree, count - degree));
        const double h = random.Next(2.0, 40.0);
        const double w = std::pow(10.0, random.Next(0.0, 3.0));
        const BSplineCurve curve = SpikeCurve(degree, count, spike, h, w);
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(curve);
        // Points above the spike's tip, where the closest point is on it.
        const Point3 p{spike + random.Next(-0.5, 0.5), h + random.Next(0.0, 5.0),
                       random.Next(-1.0, 1.0)};
        const double ours = Distance(Eval(curve, Project(curve, p)), p);
        Score(t, ours, OcctDistance(p, occt), SampledDistance(curve, p, 200000), 1e-7);
    }
    Report("curve projection, spikes", t);
}

TEST_CASE("hunt: arc length of polylines with stalls and uneven knots", "[.hunt][curve]")
{
    Tally t;
    Random random(2);
    constexpr double tolerance = 1e-10;
    for (int trial = 0; trial < 400; ++trial) {
        const BSplineCurve curve = RandomPolyline(random);
        const Interval d = Domain(curve);
        const double exact = PolylineLength(curve, d.GetUpper());
        const double ours = Stereon::ArcLength(curve, d)->Length;
        ++t.Cases;
        const double ourError = std::abs(ours - exact) / exact;
        if (ourError > 1e-9) {
            ++t.StereonWrong;
            t.WorstStereon = std::max(t.WorstStereon, ourError);
        }
        try {
            const GeomAdaptor_Curve adaptor(ToOcct(curve));
            const double occt =
                GCPnts_AbscissaPoint::Length(adaptor, d.GetLower(), d.GetUpper(), tolerance);
            const double occtError = std::abs(occt - exact) / exact;
            if (!std::isfinite(occt) || occtError > 1e-9) {
                ++t.OcctWrong;
                t.WorstOcct = std::max(t.WorstOcct, occtError);
            }
        } catch (const Standard_Failure&) {
            ++t.OcctNoAnswer;
        }
    }
    Report("arc length, polylines", t);
}

TEST_CASE("hunt: point at a length on polylines with stalls and uneven knots", "[.hunt][curve]")
{
    Tally t;
    Random random(3);
    constexpr double tolerance = 1e-10;
    for (int trial = 0; trial < 400; ++trial) {
        const BSplineCurve curve = RandomPolyline(random);
        const Interval d = Domain(curve);
        const double total = PolylineLength(curve, d.GetUpper());
        const double target = total * random.Next(0.05, 0.95);
        ++t.Cases;
        const auto ours = Stereon::ParameterAtArcLength(curve, d.GetLower(), target);
        const double ourError =
            ours ? std::abs(PolylineLength(curve, *ours) - target) / total : 1.0;
        if (ourError > 1e-9) {
            ++t.StereonWrong;
            if (ourError > t.WorstStereon) {
                PrintCase("Stereon point at a length", curve, target, ourError);
            }
            t.WorstStereon = std::max(t.WorstStereon, ourError);
        }
        try {
            const GeomAdaptor_Curve adaptor(ToOcct(curve));
            const GCPnts_AbscissaPoint occt(tolerance, adaptor, target, d.GetLower());
            if (!occt.IsDone()) {
                ++t.OcctNoAnswer;
                continue;
            }
            const double occtError =
                std::abs(PolylineLength(curve, occt.Parameter()) - target) / total;
            if (occtError > 1e-9) {
                ++t.OcctWrong;
                if (occtError > t.WorstOcct) {
                    PrintCase("OCCT point at a length", curve, target, occtError);
                }
                t.WorstOcct = std::max(t.WorstOcct, occtError);
            }
        } catch (const Standard_Failure&) {
            ++t.OcctNoAnswer;
        }
    }
    Report("point at a length, polylines", t);
}

TEST_CASE("hunt: projection onto narrow bumps", "[.hunt][surface]")
{
    Tally t;
    Random random(4);
    for (int trial = 0; trial < 120; ++trial) {
        const int n = static_cast<int>(random.Next(12.0, 40.0));
        const int bu = static_cast<int>(random.Next(3.0, n - 3.0));
        const int bv = static_cast<int>(random.Next(3.0, n - 3.0));
        const double h = random.Next(1.0, 30.0);
        const BSplineSurface surface = BumpSurface(n, bu, bv, h);
        const occ::handle<Geom_BSplineSurface> occt = ToOcct(surface);
        // Points beside and above the bump, where its flank is closest.
        const Point3 p{bu + random.Next(-2.0, 2.0), bv + random.Next(-2.0, 2.0),
                       random.Next(0.3, 1.2) * h * 0.25};
        const double ours = Distance(Eval(surface, Project(surface, p)), p);
        const std::optional<double> theirs = OcctDistance(p, occt);
        if (theirs && *theirs - ours > t.WorstOcct) {
            std::printf("HUNT CASE bump: n %d, pole (%d, %d), h %.17g, p (%.17g, %.17g, %.17g), "
                        "Stereon %.17g, OCCT %.17g\n",
                        n, bu, bv, h, p.X, p.Y, p.Z, ours, *theirs);
        }
        Score(t, ours, theirs, SampledDistance(surface, p, 600), 1e-6);
    }
    Report("surface projection, bumps", t);
}

TEST_CASE("hunt: projection between the branches of a hairpin", "[.hunt][curve]")
{
    Tally t;
    Random random(5);
    for (int trial = 0; trial < 400; ++trial) {
        const double gap = std::pow(10.0, random.Next(-4.0, 0.0));
        const BSplineCurve curve =
            Hairpin(random.Next(5.0, 50.0), gap, static_cast<int>(random.Next(4.0, 20.0)));
        const occ::handle<Geom_BSplineCurve> occt = ToOcct(curve);
        const Interval d = Domain(curve);
        // A point just off the curve, near one branch, anywhere along it.
        const Point3 on = Eval(curve, random.Next(d.GetLower(), d.GetUpper()));
        const Point3 p{on.X, on.Y + random.Next(-0.25, 0.25) * gap, random.Next(-1.0, 1.0) * gap};
        const double ours = Distance(Eval(curve, Project(curve, p)), p);
        Score(t, ours, OcctDistance(p, occt), SampledDistance(curve, p, 400000), 1e-9 + 1e-6 * gap);
    }
    Report("curve projection, hairpins", t);
}

TEST_CASE("hunt: arc length and point at a length with wild weights", "[.hunt][curve]")
{
    Tally length;
    Tally inverse;
    Random random(6);
    constexpr double tolerance = 1e-10;
    for (int trial = 0; trial < 60; ++trial) {
        const BSplineCurve curve = WildWeights(random);
        const Interval d = Domain(curve);
        const double exact = PolylineReference(curve, d.GetLower(), d.GetUpper());
        const GeomAdaptor_Curve adaptor(ToOcct(curve));
        ++length.Cases;
        const double ours = std::abs(Stereon::ArcLength(curve, d)->Length - exact) / exact;
        if (ours > 1e-8) {
            ++length.StereonWrong;
            length.WorstStereon = std::max(length.WorstStereon, ours);
        }
        try {
            const double occt = std::abs(GCPnts_AbscissaPoint::Length(adaptor, d.GetLower(),
                                                                      d.GetUpper(), tolerance) -
                                         exact) /
                                exact;
            if (!(occt <= 1e-8)) {
                std::printf("HUNT CASE weights trial %d: OCCT length error %.3g\n", trial, occt);
                ++length.OcctWrong;
                length.WorstOcct = std::max(length.WorstOcct, occt);
            }
        } catch (const Standard_Failure&) {
            ++length.OcctNoAnswer;
        }
        // The inverse: the parameter at a random fraction of the length,
        // checked by the reference length to it.
        const double target = exact * random.Next(0.05, 0.95);
        ++inverse.Cases;
        const auto t = Stereon::ParameterAtArcLength(curve, d.GetLower(), target);
        const double ourError =
            t ? std::abs(PolylineReference(curve, d.GetLower(), *t) - target) / exact : 1.0;
        if (ourError > 1e-8) {
            ++inverse.StereonWrong;
            inverse.WorstStereon = std::max(inverse.WorstStereon, ourError);
        }
        try {
            const GCPnts_AbscissaPoint occt(tolerance, adaptor, target, d.GetLower());
            if (!occt.IsDone()) {
                ++inverse.OcctNoAnswer;
                continue;
            }
            const double occtError =
                std::abs(PolylineReference(curve, d.GetLower(), occt.Parameter()) - target) / exact;
            if (occtError > 1e-8) {
                std::printf("HUNT CASE weights trial %d: OCCT inverse error %.3g, target %.17g\n",
                            trial, occtError, target);
                ++inverse.OcctWrong;
                inverse.WorstOcct = std::max(inverse.WorstOcct, occtError);
            }
        } catch (const Standard_Failure&) {
            ++inverse.OcctNoAnswer;
        }
    }
    Report("arc length, wild weights", length);
    Report("point at a length, wild weights", inverse);
}

TEST_CASE("hunt: projection near the poles of a sphere", "[.hunt][surface][singular]")
{
    Tally analytic;
    Tally nurbs;
    Random random(7);
    constexpr double radius = 2.5;
    const Sphere sphere{Origin, radius, AxisX, AxisY};
    const BSplineSurface converted = *ToBSpline(sphere, Domain(sphere));
    const occ::handle<Geom_Surface> occtSphere = new Geom_SphericalSurface(OcctFrame(), radius);
    const occ::handle<Geom_Surface> occtNurbs = ToOcct(converted);
    for (int trial = 0; trial < 300; ++trial) {
        // Inside, on or outside the sphere, near the north or south pole.
        const double side = random.Next(0.0, 1.0) < 0.5 ? 1.0 : -1.0;
        const Point3 p = NearAxis(random, side * radius * random.Next(0.5, 1.5));
        const double exact = std::abs(Distance(p, Origin) - radius);
        ScoreExact(analytic, Distance(Eval(sphere, Project(sphere, p)), p),
                   OcctDistance(p, occtSphere), exact, 1e-12);
        const std::optional<double> theirs = OcctDistance(p, occtNurbs);
        if (theirs && *theirs - exact > std::max(1e-9, nurbs.WorstOcct)) {
            std::printf("HUNT CASE NURBS pole: p (%.17g, %.17g, %.17g), exact %.17g, OCCT %.17g\n",
                        p.X, p.Y, p.Z, exact, *theirs);
        }
        ScoreExact(nurbs, Distance(Eval(converted, Project(converted, p)), p), theirs, exact, 1e-9);
    }
    Report("sphere projection near the poles, analytic", analytic);
    Report("sphere projection near the poles, NURBS", nurbs);
}

TEST_CASE("hunt: curvature next to the collapsed edge of a NURBS sphere",
          "[.hunt][surface][singular]")
{
    Tally t;
    Random random(8);
    constexpr double radius = 2.5;
    const Sphere sphere{Origin, radius, AxisX, AxisY};
    const BSplineSurface converted = *ToBSpline(sphere, Domain(sphere));
    const occ::handle<Geom_Surface> occt = ToOcct(converted);
    const Stereon::Box2 d = Domain(converted);
    for (int trial = 0; trial < 300; ++trial) {
        // A parameter 10^e of the v range away from either collapsed edge.
        const double offset = std::pow(10.0, random.Next(-14.0, -2.0)) * Width(d.Y);
        const double v =
            random.Next(0.0, 1.0) < 0.5 ? d.Y.GetLower() + offset : d.Y.GetUpper() - offset;
        const Point2 uv{random.Next(d.X.GetLower(), d.X.GetUpper()), v};
        ++t.Cases;
        // Both principal curvatures of a sphere are 1/R in magnitude.
        const auto ours = Stereon::PrincipalCurvatures(converted, uv);
        const double ourError = ours ? std::max(std::abs(std::abs(ours->Min) * radius - 1.0),
                                                std::abs(std::abs(ours->Max) * radius - 1.0))
                                     : 1.0;
        if (ourError > 1e-6) {
            ++t.StereonWrong;
            t.WorstStereon = std::max(t.WorstStereon, ourError);
        }
        try {
            GeomLProp_SLProps props(occt, uv.X, uv.Y, 2, 1e-9);
            if (!props.IsCurvatureDefined()) {
                ++t.OcctNoAnswer;
                continue;
            }
            const double occtError =
                std::max(std::abs(std::abs(props.MinCurvature()) * radius - 1.0),
                         std::abs(std::abs(props.MaxCurvature()) * radius - 1.0));
            if (!(occtError <= 1e-6)) {
                ++t.OcctWrong;
                t.WorstOcct = std::max(t.WorstOcct, occtError);
            }
        } catch (const Standard_Failure&) {
            ++t.OcctNoAnswer;
        }
    }
    Report("sphere curvature next to a collapsed edge, NURBS", t);
}

TEST_CASE("hunt: projection near periodic seams", "[.hunt][curve][surface][seam]")
{
    Tally circles;
    Tally closed;
    Tally cylinders;
    Tally tori;
    Random random(9);
    constexpr double big = 3.0;
    constexpr double small = 0.75;
    const Circle circle{Origin, big, AxisX, AxisY};
    const BSplineCurve loop = *ToBSpline(circle, Interval{0.0, 2.0 * std::numbers::pi});
    const Cylinder cylinder{Origin, big, AxisX, AxisY};
    const Torus torus{Origin, big, small, AxisX, AxisY};
    const occ::handle<Geom_Curve> occtCircle =
        new Geom_Circle(gp_Ax2(ToOcct(Origin), gp_Dir(AxisZ.X, AxisZ.Y, AxisZ.Z),
                               gp_Dir(AxisX.X, AxisX.Y, AxisX.Z)),
                        big);
    const occ::handle<Geom_Curve> occtLoop = ToOcct(loop);
    const occ::handle<Geom_Surface> occtCylinder = new Geom_CylindricalSurface(OcctFrame(), big);
    const occ::handle<Geom_Surface> occtTorus = new Geom_ToroidalSurface(OcctFrame(), big, small);
    for (int trial = 0; trial < 300; ++trial) {
        const double a = NearSeam(random);
        const double rho = random.Next(0.2, 2.0) * big;
        const double z = random.Next(-1.0, 1.0);
        const Point3 p = InFrame(rho * std::cos(a), rho * std::sin(a), z);
        const double toCircle = std::hypot(rho - big, z);
        ScoreExact(circles, Distance(Eval(circle, Project(circle, p)), p),
                   OcctDistance(p, occtCircle), toCircle, 1e-12);
        ScoreExact(closed, Distance(Eval(loop, Project(loop, p)), p), OcctDistance(p, occtLoop),
                   toCircle, 1e-9);
        ScoreExact(cylinders, Distance(Eval(cylinder, Project(cylinder, p)), p),
                   OcctDistance(p, occtCylinder), std::abs(rho - big), 1e-12);
        // Near the torus's seam in v too: just off the outer equator.
        const double b = NearSeam(random);
        const double off = random.Next(0.1, 1.9) * small;
        const Point3 q = InFrame((big + off * std::cos(b)) * std::cos(a),
                                 (big + off * std::cos(b)) * std::sin(a), off * std::sin(b));
        ScoreExact(tori, Distance(Eval(torus, Project(torus, q)), q), OcctDistance(q, occtTorus),
                   std::abs(off - small), 1e-12);
    }
    Report("projection near the seam, circle", circles);
    Report("projection near the seam, closed NURBS circle", closed);
    Report("projection near the seam, cylinder", cylinders);
    Report("projection near the seams, torus", tori);
}

TEST_CASE("hunt: lengths across the seam of a circle", "[.hunt][curve][seam]")
{
    Tally length;
    Tally inverse;
    Random random(10);
    constexpr double radius = 3.0;
    const Circle circle{Origin, radius, AxisX, AxisY};
    const occ::handle<Geom_Curve> occtCircle =
        new Geom_Circle(gp_Ax2(ToOcct(Origin), gp_Dir(AxisZ.X, AxisZ.Y, AxisZ.Z),
                               gp_Dir(AxisX.X, AxisX.Y, AxisX.Z)),
                        radius);
    const GeomAdaptor_Curve adaptor(occtCircle);
    constexpr double tolerance = 1e-10;
    constexpr double turn = 2.0 * std::numbers::pi;
    for (int trial = 0; trial < 300; ++trial) {
        // From just before the seam, forwards over up to five turns.
        const double start = turn - std::pow(10.0, random.Next(-12.0, 0.0));
        const double angle = random.Next(1e-6, 5.0 * turn);
        const double exact = radius * angle;
        ++length.Cases;
        const double ours = Stereon::ArcLength(circle, Interval{start, start + angle})->Length;
        if (std::abs(ours - exact) > 1e-9 * exact) {
            ++length.StereonWrong;
            length.WorstStereon = std::max(length.WorstStereon, std::abs(ours - exact) / exact);
        }
        try {
            const double occt =
                GCPnts_AbscissaPoint::Length(adaptor, start, start + angle, tolerance);
            if (!(std::abs(occt - exact) <= 1e-9 * exact)) {
                ++length.OcctWrong;
                length.WorstOcct = std::max(length.WorstOcct, std::abs(occt - exact) / exact);
            }
        } catch (const Standard_Failure&) {
            ++length.OcctNoAnswer;
        }
        // The inverse: the parameter at that length, compared by the arc
        // between it and the exact one, modulo whole turns.
        const auto arcError = [&](const double t) {
            const double gap = std::remainder(t - (start + angle), turn);
            return std::abs(gap) * radius / exact;
        };
        ++inverse.Cases;
        const auto t = Stereon::ParameterAtArcLength(circle, start, exact);
        const double ourError = t ? arcError(*t) : 1.0;
        if (ourError > 1e-9) {
            ++inverse.StereonWrong;
            inverse.WorstStereon = std::max(inverse.WorstStereon, ourError);
        }
        try {
            const GCPnts_AbscissaPoint occt(tolerance, adaptor, exact, start);
            if (!occt.IsDone()) {
                ++inverse.OcctNoAnswer;
                continue;
            }
            if (arcError(occt.Parameter()) > 1e-9) {
                ++inverse.OcctWrong;
                inverse.WorstOcct = std::max(inverse.WorstOcct, arcError(occt.Parameter()));
            }
        } catch (const Standard_Failure&) {
            ++inverse.OcctNoAnswer;
        }
    }
    Report("arc length across the seam, circle", length);
    Report("point at a length across the seam, circle", inverse);
}
