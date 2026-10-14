// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Exact NURBS conversion of the analytic curves. The converted curve must
// trace exactly the original's arc over the range: every point of it lies on
// the original inside the range, every point of the original's arc lies on
// it, the two agree at every knot, and they run the same way.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::CirclePeriod;
using Stereon::Curve;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Parabola;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-11;

struct Frame {
    Vec3 X;
    Vec3 Y;
};

Frame TiltedFrame()
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 1.0, -2.0}, 0.9);
    return {Transformed(Vec3{1.0, 0.0, 0.0}, tilt), Transformed(Vec3{0.0, 1.0, 0.0}, tilt)};
}

/// Parameters equal modulo the period, or plainly equal for period 0.
bool SameParameter(const double a, const double b, const double period)
{
    if (period == 0.0) {
        return std::abs(a - b) <= 1e-9 * (1.0 + std::abs(b));
    }
    const double gap = std::fmod(std::abs(a - b), period);
    return std::min(gap, period - gap) <= 1e-9;
}

/// True when t lies in [low, high] modulo the period.
bool InRange(const double t, const Interval range, const double period)
{
    if (period == 0.0) {
        return t >= range.GetLower() - 1e-9 && t <= range.GetUpper() + 1e-9;
    }
    const double shifted = std::fmod(t - range.GetLower(), period);
    const double offset = shifted < 0.0 ? shifted + period : shifted;
    return offset <= range.GetUpper() - range.GetLower() + 1e-9 || offset >= period - 1e-9;
}

template <typename G>
double PeriodOf(const G& g)
{
    if constexpr (Stereon::PeriodicCurve<G>) {
        return Period(g);
    } else {
        static_cast<void>(g);
        return 0.0;
    }
}

/// Every property of an exact conversion, at a sampled parameter.
template <typename G>
void CheckConversion(const G& g, const Interval range, const bool sameParametrisation)
{
    const auto converted = Stereon::ToBSpline(g, range);
    RC_ASSERT(converted.has_value());
    const BSplineCurve& c = *converted;
    RC_ASSERT(Domain(c) == range);
    const double period = PeriodOf(g);
    for (const double knot : c.Knots) {
        RC_ASSERT(Distance(Eval(c, knot), Eval(g, knot)) <= Tolerance);
    }
    const double s = range.GetLower() +
                     (*rc::gen::inRange(0, 1000) / 1000.0) * (range.GetUpper() - range.GetLower());
    const Point3 p = Eval(c, s);
    const double t = Project(g, p);
    RC_ASSERT(Distance(Eval(g, t), p) <= Tolerance);
    RC_ASSERT(InRange(t, range, period));
    RC_ASSERT(Dot(Deriv(c, s, 1), Deriv(g, t, 1)) > 0.0);
    if (sameParametrisation) {
        RC_ASSERT(Distance(Eval(c, s), Eval(g, s)) <= Tolerance);
    }
    // Coverage: the original's point at s lies on the converted curve.
    const Point3 q = Eval(g, s);
    RC_ASSERT(Distance(Eval(c, Project(c, q)), q) <= 1e-9);
    RC_ASSERT(SameParameter(Project(g, Eval(c, Project(c, q))), s, period));
}

Interval RandomRange(const double start, const double maxLength)
{
    const double low = start * *rc::gen::inRange(-1000, 1000) / 1000.0;
    const double length = maxLength * *rc::gen::inRange(1, 1000) / 1000.0;
    return {low, low + length};
}

} // namespace

TEST_CASE("analytic curves convert exactly to NURBS", "[curve][conversion][property]")
{
    const Frame f = TiltedFrame();
    const Point3 c{1.0, -2.0, 0.5};
    rc::prop("circles and ellipses over up to one period", [&] {
        CheckConversion(Circle{c, 2.5, f.X, f.Y}, RandomRange(10.0, CirclePeriod), false);
        CheckConversion(Ellipse{c, 3.0, 1.25, f.X, f.Y}, RandomRange(10.0, CirclePeriod), false);
    });
    rc::prop("hyperbolas, parabolas and lines over finite ranges", [&] {
        CheckConversion(Hyperbola{c, 2.0, 1.5, f.X, f.Y}, RandomRange(3.0, 4.0), false);
        CheckConversion(Parabola{c, 0.75, f.X, f.Y}, RandomRange(5.0, 10.0), true);
        CheckConversion(Line{c, f.X}, RandomRange(5.0, 10.0), true);
    });
}

TEST_CASE("a full circle is four rational quadratic arcs", "[curve][conversion]")
{
    const Circle circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const BSplineCurve c = *Stereon::ToBSpline(circle, Interval{0.0, CirclePeriod});
    REQUIRE(c.Degree == 2);
    REQUIRE(c.Poles.size() == 9);
    REQUIRE(Stereon::IsRational(c));
    REQUIRE(std::abs(c.Weights[1] - std::numbers::sqrt2 / 2.0) <= 1e-15);
    REQUIRE(Eval(c, CirclePeriod) == Eval(c, 0.0)); // closed, both ends exact
}

TEST_CASE("a B-spline converts to its restriction", "[curve][conversion]")
{
    const BSplineCurve curve{3,
                             {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0},
                             {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
                              Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}},
                             {}};
    const BSplineCurve piece = *Stereon::ToBSpline(curve, Interval{0.25, 1.5});
    REQUIRE(Domain(piece) == Interval{0.25, 1.5});
    for (int i = 0; i <= 10; ++i) {
        const double t = 0.25 + i * 0.125;
        REQUIRE(Distance(Eval(piece, t), Eval(curve, t)) <= 1e-13);
    }
    REQUIRE(*Stereon::ToBSpline(curve, Domain(curve)) == curve);
}

TEST_CASE("conversion through the wrapper and its errors", "[curve][conversion]")
{
    const Circle circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(*Stereon::ToBSpline(Curve(circle), Interval{0.0, 1.0}) ==
            *Stereon::ToBSpline(circle, Interval{0.0, 1.0}));
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::ToBSpline(circle, Interval{0.0, 7.0}).error().GetCode() == invalid);
    REQUIRE(Stereon::ToBSpline(circle, Interval{1.0, 1.0}).error().GetCode() == invalid);
    REQUIRE(Stereon::ToBSpline(Line{Point3{}, Vec3{1.0, 0.0, 0.0}}, Interval::Everything())
                .error()
                .GetCode() == invalid);
    // Ranges of a periodic curve need not lie in [0, 2 pi].
    REQUIRE(Stereon::ToBSpline(circle, Interval{-10.0, -9.0}).has_value());
}
