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
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Standard_Failure.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <vector>

import Stereon.Surface;

#include "OcctConversion.hpp"
#include "ReferenceLengths.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
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
std::optional<double> OcctDistance(const Point3& p, const occ::handle<Geom_BSplineCurve>& c)
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

std::optional<double> OcctDistance(const Point3& p, const occ::handle<Geom_BSplineSurface>& s)
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
        std::printf("HUNT CASE weights trial %d: length %.17g\n  points:", trial, exact);
        for (std::size_t i = 0; i < curve.Poles.size(); ++i) {
            std::printf(" (%.17g, %.17g, %.17g; %.17g)", curve.Poles[i].X, curve.Poles[i].Y,
                        curve.Poles[i].Z, curve.Weights[i]);
        }
        std::printf("\n");
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
