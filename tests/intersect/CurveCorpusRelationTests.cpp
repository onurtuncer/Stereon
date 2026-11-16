// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0053 relations MR-01, MR-03 and MR-04 on every case of
// corpus/curve-curve/, for every ordered pair of 3D curves. Each case is a
// pair whose closest approach near stated parameters is a stated gap,
// constructed in closed form by tools/generate_curve_corpus.py: the curves
// meet near those parameters exactly when the gap is within the
// resolution. The relations move the case exactly and map the stated
// parameters in closed form, and the claim must still hold, in both orders:
//
// - MR-01: cyclic axis permutations and a quarter turn about z (exact
//   rotations), with and without the translation (8, -4, 2);
// - MR-03: scales 1/4 and 4, the resolution and the tolerance ceiling
//   scaled with them (ADR-0053: policy transformed separately); a
//   parabola's parameter is a length and scales, the others do not;
// - MR-04: one operand at a time in an equivalent chart built from the
//   data, a line's origin moved two directions on with its direction
//   doubled (t' = (t - 2) / 2), a circle's frame turned a quarter
//   (t' = t - pi/2), an ellipse's turned half (t' = t - pi), a hyperbola's
//   or a parabola's y axis negated (t' = -t), a B-spline reversed
//   (t' = a + b - t).
//
// MR-02 is the corpus runner's own both-orders check (CurveCorpusTests).

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"
#include "CurveCorpus.hpp"
#include "CurveWitnessChecks.hpp"

using Stereon::Context;
using Stereon::Curve;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

/// How far a reported parameter may be from the mapped one, as the corpus
/// runner allows.
constexpr double Near = 1e-4;
constexpr double pi = std::numbers::pi;

/// p -> Offset + Scale Q p, with Q an exact rotation given by its columns.
struct Motion {
    std::array<Vec3, 3> Columns{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    Vec3 Offset{};
    double Scale = 1.0;

    [[nodiscard]] Vec3 Turn(const Vec3& v) const
    {
        return v.X * Columns[0] + v.Y * Columns[1] + v.Z * Columns[2];
    }
    [[nodiscard]] Point3 Of(const Point3& p) const
    {
        return Point3{} + Offset + Scale * Turn(p - Point3{});
    }
    [[nodiscard]] Context Policy() const
    {
        return Context::Default()
            .WithResolution(Context::DefaultResolution * Scale)
            .WithToleranceCeiling(Context::DefaultToleranceCeiling * Scale)
            .WithScheduler(Stereon::Scheduler::Serial);
    }
};

/// A curve of the case moved, with its domain and the map of its parameter.
struct Moved {
    Curve Value;
    Interval Domain;
    double (*Forward)(double, double) = nullptr; ///< t' from t and an argument
    double Argument = 0.0;
};

double Same(const double t, double /*unused*/)
{
    return t;
}
double Scaled(const double t, const double s)
{
    return s * t;
}
double LineChart(const double t, double /*unused*/)
{
    return (t - 2.0) / 2.0;
}
double QuarterChart(const double t, double /*unused*/)
{
    return t - pi / 2.0;
}
double HalfChart(const double t, double /*unused*/)
{
    return t - pi;
}
double Negated(const double t, double /*unused*/)
{
    return -t;
}
double Reflected(const double t, const double sum)
{
    return sum - t;
}

Interval MapDomain(const Interval& d, double (*f)(double, double), const double arg)
{
    const double a = f(d.GetLower(), arg);
    const double b = f(d.GetUpper(), arg);
    return Interval{std::min(a, b), std::max(a, b)};
}

/// The curve under a motion: points moved, unit frames turned, lengths and
/// a line's direction scaled. Parameters keep their meaning, but a
/// parabola's, a length along its y axis, scales.
Moved Under(const JsonValue& c, const Motion& m)
{
    using namespace Stereon;
    const std::string type = c["type"].AsString();
    const Interval domain = StereonTest::DomainFrom(c);
    const auto point = [&](const char* key) { return m.Of(StereonTest::PointFrom(c[key])); };
    const auto unit = [&](const char* key) { return m.Turn(StereonTest::VectorFrom(c[key])); };
    if (type == "line") {
        return {Curve{Line{point("origin"), m.Scale * unit("direction")}}, domain, &Same};
    }
    if (type == "bspline") {
        BSplineCurve s = StereonTest::SplineFrom(c);
        for (Point3& p : s.Poles) {
            p = m.Of(p);
        }
        return {Curve{std::move(s)}, domain, &Same};
    }
    const Vec3 x = unit("x_axis");
    const Vec3 y = unit("y_axis");
    if (type == "circle") {
        return {Curve{Circle{point("center"), m.Scale * c["radius"].AsNumber(), x, y}}, domain,
                &Same};
    }
    if (type == "ellipse") {
        return {Curve{Ellipse{point("center"), m.Scale * c["major"].AsNumber(),
                              m.Scale * c["minor"].AsNumber(), x, y}},
                domain, &Same};
    }
    if (type == "parabola") {
        return {Curve{Parabola{point("vertex"), m.Scale * c["focal"].AsNumber(), x, y}},
                MapDomain(domain, &Scaled, m.Scale), &Scaled, m.Scale};
    }
    return {Curve{Hyperbola{point("center"), m.Scale * c["major"].AsNumber(),
                            m.Scale * c["minor"].AsNumber(), x, y}},
            domain, &Same};
}

/// The curve in an equivalent chart, built from the data (MR-04).
Moved Rechart(const JsonValue& c)
{
    using namespace Stereon;
    const std::string type = c["type"].AsString();
    const Interval domain = StereonTest::DomainFrom(c);
    const auto with = [&](Curve curve, double (*f)(double, double), const double arg = 0.0) {
        return Moved{std::move(curve), MapDomain(domain, f, arg), f, arg};
    };
    if (type == "line") {
        const Point3 o = StereonTest::PointFrom(c["origin"]);
        const Vec3 d = StereonTest::VectorFrom(c["direction"]);
        return with(Curve{Line{o + 2.0 * d, 2.0 * d}}, &LineChart);
    }
    if (type == "bspline") {
        BSplineCurve s = StereonTest::SplineFrom(c);
        const double sum = s.Knots.front() + s.Knots.back();
        std::vector<double> knots;
        for (auto it = s.Knots.rbegin(); it != s.Knots.rend(); ++it) {
            knots.push_back(sum - *it);
        }
        std::ranges::reverse(s.Poles);
        s.Knots = std::move(knots);
        return with(Curve{std::move(s)}, &Reflected, sum);
    }
    const Vec3 x = StereonTest::VectorFrom(c["x_axis"]);
    const Vec3 y = StereonTest::VectorFrom(c["y_axis"]);
    if (type == "circle") {
        return with(
            Curve{Circle{StereonTest::PointFrom(c["center"]), c["radius"].AsNumber(), y, -x}},
            &QuarterChart);
    }
    if (type == "ellipse") {
        return with(Curve{Ellipse{StereonTest::PointFrom(c["center"]), c["major"].AsNumber(),
                                  c["minor"].AsNumber(), -x, -y}},
                    &HalfChart);
    }
    if (type == "parabola") {
        return with(
            Curve{Parabola{StereonTest::PointFrom(c["vertex"]), c["focal"].AsNumber(), x, -y}},
            &Negated);
    }
    return with(Curve{Hyperbola{StereonTest::PointFrom(c["center"]), c["major"].AsNumber(),
                                c["minor"].AsNumber(), x, -y}},
                &Negated);
}

bool MeetsNear(const Stereon::CurveCurveIntersection& r, const double a, const double b,
               const bool swapped)
{
    for (const Stereon::CurveCurvePoint& p : r.Points) {
        const double onA = swapped ? p.ParameterB : p.ParameterA;
        const double onB = swapped ? p.ParameterA : p.ParameterB;
        if (std::abs(onA - a) <= Near && std::abs(onB - b) <= Near) {
            return true;
        }
    }
    return false;
}

/// The corpus claim for the moved pair, in both orders.
void Claim(const JsonValue& c, const Moved& a, const Moved& b, const Context& policy,
           const std::string& what)
{
    const double ta = a.Forward(c["expected"]["parameters"][0].AsNumber(), a.Argument);
    const double tb = b.Forward(c["expected"]["parameters"][1].AsNumber(), b.Argument);
    const bool meets = c["expected"]["meets"].AsBool();
    for (const bool swapped : {false, true}) {
        INFO(c["id"].AsString() << ' ' << what << (swapped ? " swapped" : ""));
        const auto r = swapped ? Intersect(b.Value, b.Domain, a.Value, a.Domain, policy)
                               : Intersect(a.Value, a.Domain, b.Value, b.Domain, policy);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        CHECK(MeetsNear(*r, ta, tb, swapped) == meets);
        for (const auto& p : r->Points) {
            CHECK(IntersectTest::ValidCurveWitness(
                p, swapped ? b.Value : a.Value, swapped ? b.Domain : a.Domain,
                swapped ? a.Value : b.Value, swapped ? a.Domain : b.Domain, policy));
        }
    }
}

std::vector<JsonValue> Cases()
{
    std::vector<JsonValue> cases;
    for (const auto& file : StereonTest::CurveCorpusFiles()) {
        cases.push_back(StereonTest::ReadJsonFile(file));
    }
    REQUIRE(cases.size() >= 168);
    return cases;
}

} // namespace

TEST_CASE("MR-01: curve-curve corpus claims hold under exact rigid motions",
          "[intersect][curve-curve][corpus][metamorphic]")
{
    const std::array<Motion, 4> turns{
        Motion{{Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}}},
        Motion{{Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}},
        Motion{{Vec3{0, 1, 0}, Vec3{-1, 0, 0}, Vec3{0, 0, 1}}}, // a quarter turn about z
        Motion{},
    };
    for (const JsonValue& c : Cases()) {
        const auto& curves = c["inputs"]["curves"];
        for (std::size_t i = 0; i < turns.size(); ++i) {
            for (const Vec3 offset : {Vec3{}, Vec3{8, -4, 2}}) {
                if (i + 1 == turns.size() && offset == Vec3{}) {
                    continue; // the identity is the corpus runner's own case
                }
                Motion m = turns[i];
                m.Offset = offset;
                Claim(c, Under(curves[0], m), Under(curves[1], m), m.Policy(),
                      "turn " + std::to_string(i) + " offset " + std::to_string(offset.X));
            }
        }
    }
}

TEST_CASE("MR-03: curve-curve corpus claims hold under power-of-two scales",
          "[intersect][curve-curve][corpus][metamorphic]")
{
    for (const JsonValue& c : Cases()) {
        const auto& curves = c["inputs"]["curves"];
        for (const double scale : {0.25, 4.0}) {
            Motion m;
            m.Scale = scale;
            Claim(c, Under(curves[0], m), Under(curves[1], m), m.Policy(),
                  "scale " + std::to_string(scale));
        }
    }
}

TEST_CASE("MR-04: curve-curve corpus claims hold in equivalent charts",
          "[intersect][curve-curve][corpus][metamorphic]")
{
    const Motion identity;
    for (const JsonValue& c : Cases()) {
        const auto& curves = c["inputs"]["curves"];
        Claim(c, Rechart(curves[0]), Under(curves[1], identity), identity.Policy(), "A recharted");
        Claim(c, Under(curves[0], identity), Rechart(curves[1]), identity.Policy(), "B recharted");
    }
}
