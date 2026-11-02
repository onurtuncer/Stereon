// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A B-spline curve in the periodic form is periodic: its period is the width
// of its domain, it evaluates at every finite parameter, knot insertion keeps
// the form, and its pieces may cross the seam. The fixtures are both makers
// of the form: a periodic interpolation (seam of multiplicity 1, polynomial)
// and AsPeriodic of a full circle (seam of multiplicity p, rational).

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::ArcLength;
using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::Curve;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

/// A wavy closed loop, interpolated periodically: a seam of multiplicity 1.
BSplineCurve Loop(const int degree)
{
    std::vector<Point3> points;
    for (int k = 0; k < 9; ++k) {
        const double a = 2.0 * pi * k / 9.0;
        points.push_back({3.0 * std::cos(a), 2.0 * std::sin(a), 0.5 * std::sin(3.0 * a)});
    }
    return Stereon::InterpolatePeriodic(points, degree)->Curve;
}

/// The unit circle's exact NURBS form over one turn, put in the periodic
/// form: rational, with a seam of multiplicity p.
BSplineCurve Circle()
{
    const Stereon::Circle circle{Point3{1.0, 2.0, 0.0}, 2.0, Vec3{1.0, 0.0, 0.0},
                                 Vec3{0.0, 1.0, 0.0}};
    const BSplineCurve closed = *Stereon::ToBSpline(circle, Interval{0.0, 2.0 * pi});
    return *Stereon::AsPeriodic(closed, 1e-12);
}

BSplineCurve2d Loop2d()
{
    std::vector<Point2> points;
    for (int k = 0; k < 7; ++k) {
        const double a = 2.0 * pi * k / 7.0;
        points.push_back({2.0 * std::cos(a) + 0.3 * std::cos(2.0 * a), std::sin(a)});
    }
    return Stereon::InterpolatePeriodic(std::span<const Point2>{points}, 3)->Curve;
}

std::vector<BSplineCurve> Fixtures()
{
    return {Loop(2), Loop(3), Loop(5), Circle()};
}

/// The largest distance between the two curves over `samples` parameters
/// spread over three periods around the first curve's domain.
double Apart(const BSplineCurve& a, const BSplineCurve& b, const int samples = 300)
{
    const Interval domain = Stereon::Domain(a);
    const double width = Width(domain);
    double worst = 0.0;
    for (int k = 0; k <= samples; ++k) {
        const double t = domain.GetLower() - width + 3.0 * width * k / samples;
        worst = std::max(worst, Distance(Eval(a, t), Eval(b, t)));
    }
    return worst;
}

/// The multiplicity of u among one period's knots, [a, b), counting knots
/// within 1e-12 of it: a knot inserted at an image of u is u carried back by
/// whole periods, which rounds.
std::size_t MultiplicityInPeriod(const BSplineCurve& c, const double u)
{
    const Interval domain = Stereon::Domain(c);
    return static_cast<std::size_t>(std::ranges::count_if(c.Knots, [&](const double k) {
        return std::abs(k - u) <= 1e-12 && k >= domain.GetLower() && k < domain.GetUpper();
    }));
}

} // namespace

static_assert(Stereon::PeriodicCurve<BSplineCurve>);
static_assert(Stereon::PeriodicCurve<BSplineCurve2d>);

TEST_CASE("a B-spline curve in the periodic form has its domain's width as its period",
          "[curve][periodic]")
{
    for (const BSplineCurve& c : Fixtures()) {
        REQUIRE(Stereon::IsPeriodic(c));
        REQUIRE(Stereon::Period(c) == Width(Stereon::Domain(c)));
        REQUIRE(Stereon::Period(Curve{c}) == Stereon::Period(c));
        // Clamped, it is the same closed curve but no longer in the form.
        REQUIRE(Stereon::Period(Stereon::Clamped(c)) == 0.0);
    }
    const BSplineCurve open{2,
                            {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                            {Point3{}, Point3{1.0, 1.0, 0.0}, Point3{2.0, 0.0, 0.0}},
                            {}};
    REQUIRE(Stereon::Period(open) == 0.0);
    REQUIRE(Stereon::Period(Curve{open}) == 0.0);
    REQUIRE(Stereon::Period(Loop2d()) == Width(Stereon::Domain(Loop2d())));
}

TEST_CASE("a periodic B-spline's seam is right-continuous: its end evaluates as its start",
          "[curve][periodic]")
{
    // The property below found it with its 89th sample: at t = a and
    // a + T = b, the quadratic loop's second derivative came from the first
    // span at a and from the last at b, 246 apart, since the seam of a
    // periodic interpolation is only C^(p-1). Every order now repeats.
    for (const BSplineCurve& c : Fixtures()) {
        const Interval domain = Stereon::Domain(c);
        const double a = domain.GetLower();
        const double b = domain.GetUpper();
        INFO("degree " << c.Degree);
        REQUIRE(Eval(c, b) == Eval(c, a));
        for (int n = 1; n <= c.Degree + 1; ++n) {
            INFO("order " << n);
            REQUIRE(Deriv(c, b, n) == Deriv(c, a, n));
        }
        // The batch evaluators carry b the same way, so they keep Eval's bits.
        const std::vector<double> t{a, b, b + Stereon::Period(c)};
        std::vector<Point3> batch(t.size());
        Stereon::EvalAll(c, t, batch);
        for (std::size_t i = 0; i < t.size(); ++i) {
            REQUIRE(batch[i] == Eval(c, t[i]));
        }
        std::vector<Point3> fast(t.size());
        Stereon::EvalAllFast(Stereon::ToSpanPolynomials(c), t, fast);
        Stereon::EvalAllFast(c, t, batch);
        REQUIRE(fast == batch);
    }
}

TEST_CASE("a periodic B-spline evaluates at every parameter, a period apart alike",
          "[curve][periodic][property]")
{
    const std::vector<BSplineCurve> fixtures = Fixtures();
    rc::prop("points and derivatives repeat with the period", [&] {
        const BSplineCurve& c = fixtures[*rc::gen::inRange<std::size_t>(0, fixtures.size())];
        const Interval domain = Stereon::Domain(c);
        const double period = Stereon::Period(c);
        const double t = domain.GetLower() + period * *rc::gen::inRange(0, 1000) / 1000.0;
        const int k = *rc::gen::inRange(-4, 5);
        const double shifted = t + k * period;
        // The image is carried back by whole periods, so it differs from t
        // by the rounding of t + k T: the points by that much times the speed.
        const double slack = 1e-12 * (1.0 + std::abs(shifted)) * 100.0;
        RC_ASSERT(Distance(Eval(c, shifted), Eval(c, t)) <= slack);
        RC_ASSERT(Length(Deriv(c, shifted, 1) - Deriv(c, t, 1)) <= slack * 10.0);
        RC_ASSERT(Length(Deriv(c, shifted, 2) - Deriv(c, t, 2)) <= slack * 100.0);
        const Curve erased{c};
        RC_ASSERT(Distance(Eval(erased, shifted), Eval(c, t)) <= slack);
    });
}

TEST_CASE("the batch evaluators carry parameters into a periodic curve's domain",
          "[curve][periodic]")
{
    const BSplineCurve c = Loop(3);
    const Interval domain = Stereon::Domain(c);
    const double period = Stereon::Period(c);
    std::vector<double> parameters;
    for (int k = 0; k <= 64; ++k) {
        parameters.push_back(domain.GetLower() - 2.0 * period + 5.0 * period * k / 64.0);
    }
    std::vector<Point3> all(parameters.size());
    std::vector<Point3> fast(parameters.size());
    Stereon::EvalAll(c, parameters, all);
    Stereon::EvalAllFast(c, parameters, fast);
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        INFO("t = " << parameters[k]);
        REQUIRE(Distance(all[k], Eval(c, parameters[k])) == 0.0);
        REQUIRE(Distance(fast[k], Eval(c, parameters[k])) <= 1e-12);
    }
}

TEST_CASE("knot insertion keeps a curve in the periodic form", "[curve][periodic][property]")
{
    const std::vector<BSplineCurve> fixtures = Fixtures();
    rc::prop("the knot goes in at every image; the curve and its period stay", [&] {
        const BSplineCurve& c = fixtures[*rc::gen::inRange<std::size_t>(0, fixtures.size())];
        const Interval domain = Stereon::Domain(c);
        const double period = Stereon::Period(c);
        // Anywhere, in any period, the seam itself among the choices.
        const int pick = *rc::gen::inRange(0, 5);
        const double inside = pick == 0 ? domain.GetLower()
                              : pick == 1
                                  ? c.Knots[static_cast<std::size_t>(c.Degree) + 1]
                                  : domain.GetLower() + period * *rc::gen::inRange(1, 999) / 1000.0;
        const double u = inside + *rc::gen::inRange(-2, 3) * period;
        const std::size_t before = MultiplicityInPeriod(c, inside);
        RC_PRE(before < static_cast<std::size_t>(c.Degree));
        const int times = *rc::gen::inRange(1, c.Degree - static_cast<int>(before) + 1);
        const auto inserted = Stereon::InsertKnot(c, u, times);
        RC_ASSERT(inserted.has_value());
        RC_ASSERT(Stereon::IsPeriodic(*inserted));
        RC_ASSERT(Stereon::Domain(*inserted) == domain);
        RC_ASSERT(Stereon::Period(*inserted) == period);
        RC_ASSERT(inserted->Poles.size() == c.Poles.size() + static_cast<std::size_t>(times));
        RC_ASSERT(MultiplicityInPeriod(*inserted, inside) ==
                  before + static_cast<std::size_t>(times));
        RC_ASSERT(Apart(*inserted, c) <= 1e-12 * 10.0 * (1.0 + period));
    });
}

TEST_CASE("knot insertion refuses a multiplicity above the degree at the seam", "[curve][periodic]")
{
    // AsPeriodic's seam has multiplicity p: nothing more fits there, at any
    // image of it.
    const BSplineCurve c = Circle();
    const double a = Stereon::Domain(c).GetLower();
    const double period = Stereon::Period(c);
    for (const double u : {a, a + period, a - 3.0 * period}) {
        const auto refused = Stereon::InsertKnot(c, u);
        REQUIRE_FALSE(refused.has_value());
        REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("the piece of a periodic B-spline may cross its seam", "[curve][periodic]")
{
    for (const BSplineCurve& c : Fixtures()) {
        const Interval domain = Stereon::Domain(c);
        const double a = domain.GetLower();
        const double t = Stereon::Period(c);
        for (const Interval range :
             {Interval{a + 0.7 * t, a + 1.4 * t}, Interval{a - 0.5 * t, a + 0.2 * t},
              Interval{a + 5.25 * t, a + 6.1 * t}, Interval{a + 0.3 * t, a + 1.3 * t},
              Interval{a, a + t}, Interval{a + 2.0 * t, a + 2.5 * t}}) {
            INFO("range [" << range.GetLower() << ", " << range.GetUpper() << "]");
            const auto piece = Stereon::ToBSpline(c, range);
            REQUIRE(piece.has_value());
            REQUIRE(Stereon::Domain(*piece) == range);
            // A clamped piece, except the whole domain, which is the curve.
            REQUIRE(Stereon::Period(*piece) == (range == domain ? t : 0.0));
            double worst = 0.0;
            for (int k = 0; k <= 200; ++k) {
                const double u = range.GetLower() + Width(range) * k / 200.0;
                worst = std::max(worst, Distance(Eval(*piece, u), Eval(c, u)));
            }
            REQUIRE(worst <= 1e-11);
        }
        // Inside the domain the piece is the restriction, as before; the
        // whole domain gives the curve back unchanged.
        REQUIRE(*Stereon::ToBSpline(c, domain) == c);
        // More than a period is still refused.
        const auto longer = Stereon::ToBSpline(c, Interval{a, a + 1.5 * t});
        REQUIRE_FALSE(longer.has_value());
        REQUIRE(longer.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("a trimmed periodic B-spline crosses the seam", "[curve][periodic][trimmed]")
{
    const BSplineCurve c = Loop(3);
    const double b = Stereon::Domain(c).GetUpper();
    const double t = Stereon::Period(c);
    const Interval range{b - 0.25 * t, b + 0.3 * t};
    const auto trimmed = Stereon::MakeTrimmedCurve(Curve{c}, range);
    REQUIRE(trimmed.has_value());
    const auto exact = Stereon::ToBSpline(Curve{*trimmed}, range);
    REQUIRE(exact.has_value());
    for (int k = 0; k <= 100; ++k) {
        const double u = range.GetLower() + Width(range) * k / 100.0;
        REQUIRE(Distance(Eval(*trimmed, u), Eval(c, u)) == 0.0);
        REQUIRE(Distance(Eval(*exact, u), Eval(c, u)) <= 1e-11);
    }
}

TEST_CASE("a curve not in the periodic form keeps to its domain", "[curve][periodic]")
{
    const BSplineCurve open{2,
                            {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                            {Point3{}, Point3{1.0, 1.0, 0.0}, Point3{2.0, 0.0, 0.0}},
                            {}};
    REQUIRE(Stereon::InsertKnot(open, 1.5).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Stereon::ToBSpline(open, Interval{0.5, 1.2}).error().GetCode() ==
            ErrorCode::InvalidInput);
    // A closed curve clamped is not in the form either.
    const BSplineCurve clamped = Stereon::Clamped(Loop(3));
    const Interval domain = Stereon::Domain(clamped);
    REQUIRE(Stereon::InsertKnot(clamped, domain.GetUpper() + 0.1).error().GetCode() ==
            ErrorCode::InvalidInput);
}

TEST_CASE("a periodic plane B-spline is periodic as in space", "[curve][periodic]")
{
    const BSplineCurve2d c = Loop2d();
    const Interval domain = Stereon::Domain(c);
    const double period = Stereon::Period(c);
    REQUIRE(period > 0.0);
    for (int k = 0; k <= 50; ++k) {
        const double t = domain.GetLower() + period * k / 50.0;
        REQUIRE(Distance(Eval(c, t + 3.0 * period), Eval(c, t)) <= 1e-11);
    }
    const auto inserted = Stereon::InsertKnot(c, domain.GetLower() - 0.4 * period);
    REQUIRE(inserted.has_value());
    REQUIRE(Stereon::IsPeriodic(*inserted));
    for (int k = 0; k <= 50; ++k) {
        const double t = domain.GetLower() + period * k / 50.0;
        REQUIRE(Distance(Eval(*inserted, t), Eval(c, t)) <= 1e-12);
    }
}

TEST_CASE("a periodic curve opened at any parameter traces it once from there", "[curve][periodic]")
{
    for (const BSplineCurve& c : Fixtures()) {
        const double a = Stereon::Domain(c).GetLower();
        const double period = Stereon::Period(c);
        for (const double t : {a, a + 0.37 * period, a - 2.6 * period, a + 7.0 * period}) {
            INFO("t = " << t);
            const auto opened = Stereon::Opened(c, t);
            REQUIRE(opened.has_value());
            REQUIRE(Stereon::Domain(*opened) == Interval{t, t + period});
            REQUIRE(Stereon::Period(*opened) == 0.0);
            REQUIRE(Stereon::IsClosed(*opened, 1e-11));
            for (int k = 0; k <= 100; ++k) {
                const double u = t + period * k / 100.0;
                REQUIRE(Distance(Eval(*opened, u), Eval(c, u)) <= 1e-11);
            }
            const auto closed = Stereon::AsPeriodic(*opened, 1e-11);
            INFO((closed ? "" : std::string(closed.error().GetMessage())));
            REQUIRE(closed.has_value());
            REQUIRE(Stereon::IsPeriodic(*closed));
        }
    }
    const BSplineCurve open = Stereon::Clamped(Loop(3));
    REQUIRE(Stereon::Opened(open, 0.5).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("a periodic curve splits at the image of any parameter", "[curve][periodic]")
{
    const BSplineCurve c = Loop(3);
    const Interval domain = Stereon::Domain(c);
    const double period = Stereon::Period(c);
    const double inside = domain.GetLower() + 0.4 * period;
    const auto near = Stereon::Split(c, inside);
    const auto far = Stereon::Split(c, inside + 3.0 * period);
    REQUIRE(near.has_value());
    REQUIRE(far.has_value());
    REQUIRE(Stereon::Domain(far->first).GetLower() == domain.GetLower());
    REQUIRE(std::abs(Stereon::Domain(far->first).GetUpper() - inside) <= 1e-14);
    REQUIRE(Stereon::Domain(far->second).GetUpper() == domain.GetUpper());
    // The seam has nothing on one side.
    REQUIRE(Stereon::Split(c, domain.GetLower() + 2.0 * period).error().GetCode() ==
            ErrorCode::InvalidInput);
}

TEST_CASE("span polynomials of a periodic curve carry parameters in", "[curve][periodic]")
{
    const BSplineCurve c = Circle();
    const Stereon::SpanPolynomials polynomials = Stereon::ToSpanPolynomials(c);
    REQUIRE(polynomials.Period == Stereon::Period(c));
    REQUIRE(Stereon::ToSpanPolynomials(Stereon::Clamped(c)).Period == 0.0);
    const Interval domain = Stereon::Domain(c);
    const double period = Stereon::Period(c);
    std::vector<double> parameters;
    for (int k = 0; k <= 80; ++k) {
        parameters.push_back(domain.GetLower() - 3.0 * period + 7.0 * period * k / 80.0);
    }
    std::vector<Point3> points(parameters.size());
    Stereon::EvalAllFast(polynomials, parameters, points);
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        INFO("t = " << parameters[k]);
        REQUIRE(Distance(points[k], Eval(c, parameters[k])) <= 1e-12);
    }
}

TEST_CASE("the arc length of a periodic B-spline counts its turns from anywhere",
          "[curve][periodic][arclength]")
{
    // The quadrature aims at a relative 2^-40 on each range; two ranges a
    // whole number of periods apart integrate the same speeds between the
    // same kinks, so they agree to a few units of that target.
    constexpr double Agreement = 1e-10;
    for (const BSplineCurve& c : Fixtures()) {
        const Interval domain = Stereon::Domain(c);
        const double period = Stereon::Period(c);
        const double turn = ArcLength(c, domain)->Length;
        INFO("degree " << c.Degree << ", one turn " << turn);
        SECTION("any range one period long covers the curve once")
        {
            for (const double offset : {-2.5, -0.3, 0.37, 1.0 + 1e-3, 7.0}) {
                const double start = domain.GetLower() + offset * period;
                const auto result = *ArcLength(c, Interval{start, start + period});
                INFO("from " << start);
                CHECK(std::abs(result.Length - turn) <= Agreement * turn);
                CHECK(result.ErrorEstimate <= Agreement * turn);
            }
        }
        SECTION("lengths add up across the seam and over several turns")
        {
            const double a = domain.GetLower() - 0.4 * period;
            const double m = domain.GetLower() + 0.6 * period;
            const double b = domain.GetUpper() + 2.3 * period;
            const double whole = ArcLength(c, Interval{a, b})->Length;
            const double parts =
                ArcLength(c, Interval{a, m})->Length + ArcLength(c, Interval{m, b})->Length;
            CHECK(std::abs(whole - parts) <= Agreement * whole);
            CHECK(std::abs(ArcLength(c, Interval{a, a + 3.0 * period})->Length - 3.0 * turn) <=
                  Agreement * 3.0 * turn);
        }
        SECTION("the erased curve and the parameter at a length agree")
        {
            const Curve erased(c);
            const double start = domain.GetUpper() + 0.2 * period;
            const Interval range{start, start + 1.7 * period};
            CHECK(std::abs(ArcLength(erased, range)->Length - ArcLength(c, range)->Length) <=
                  Agreement * turn);
            const double t = *Stereon::ParameterAtArcLength(c, start, 2.5 * turn);
            CHECK(std::abs(ArcLength(c, Interval{start, t})->Length - 2.5 * turn) <=
                  Agreement * turn);
        }
    }
}

TEST_CASE("the arc length of a periodic B-spline is the same a period on",
          "[curve][periodic][arclength][property]")
{
    const std::vector<BSplineCurve> fixtures = Fixtures();
    rc::prop("for every range of up to five periods, anywhere", [&fixtures] {
        const BSplineCurve& c = fixtures[static_cast<std::size_t>(
            *rc::gen::inRange(0, static_cast<int>(fixtures.size())))];
        const Interval domain = Stereon::Domain(c);
        const double period = Stereon::Period(c);
        const double start = domain.GetLower() + period * *rc::gen::inRange(-4000, 4000) / 1000.0;
        const double width = period * *rc::gen::inRange(0, 5000) / 1000.0;
        const double shift = period * *rc::gen::inRange(-6, 7);
        const double here = ArcLength(c, Interval{start, start + width})->Length;
        const double there = ArcLength(c, Interval{start + shift, start + shift + width})->Length;
        RC_ASSERT(std::abs(here - there) <= 1e-10 * (1.0 + here));
    });
}

namespace {

/// The knots of a curve of degree p in the periodic form over [0, 1], period
/// 1, from the distinct knots of one period, u[0] = 0 < ... < u[m-1] < 1: one
/// period's knots with p of their images on each side. Every knot is dyadic,
/// so a parameter moved by whole periods stays exact in doubles and the seam
/// can be checked bit for bit, not to a tolerance.
std::vector<double> PeriodicKnots(const int p, const std::vector<double>& u)
{
    const auto m = u.size();
    const auto q = static_cast<std::size_t>(p);
    std::vector<double> knots;
    for (std::size_t i = q; i >= 1; --i) {
        knots.push_back(u[m - i] - 1.0);
    }
    knots.insert(knots.end(), u.begin(), u.end());
    knots.push_back(1.0);
    for (std::size_t i = 1; i <= q; ++i) {
        knots.push_back(u[i] + 1.0);
    }
    return knots;
}

/// m distinct poles and their first p repeated: the periodic form, whose seam
/// at 0 = 1 has multiplicity 1, so only its derivatives up to order p - 1 are
/// continuous there.
template <typename P>
std::vector<P> PeriodicPoles(const int p, const std::vector<P>& distinct)
{
    std::vector<P> poles = distinct;
    for (int i = 0; i < p; ++i) {
        poles.push_back(distinct[static_cast<std::size_t>(i)]);
    }
    return poles;
}

/// Nonuniform knots of one period: spans of 1/8, 3/8, 1/4, 1/4.
const std::vector<double> OnePeriod{0.0, 0.125, 0.5, 0.75};

/// A curve in the periodic form over [0, 1] of degree p, rational or not.
BSplineCurve SeamCurve(const int p, const bool rational)
{
    const std::vector<Point3> distinct{
        {0.0, 0.0, 0.0}, {4.0, 1.0, -1.0}, {5.0, 6.0, 2.0}, {-1.0, 4.0, 0.5}};
    std::vector<double> weights;
    if (rational) {
        weights = PeriodicPoles<double>(p, {1.0, 2.5, 0.75, 1.5});
    }
    return *Stereon::MakeBSplineCurve(p, PeriodicKnots(p, OnePeriod), PeriodicPoles(p, distinct),
                                      weights);
}

BSplineCurve2d SeamCurve2d(const int p)
{
    const std::vector<Point2> distinct{{0.0, 0.0}, {4.0, 1.0}, {5.0, 6.0}, {-1.0, 4.0}};
    return *Stereon::MakeBSplineCurve(p, PeriodicKnots(p, OnePeriod), PeriodicPoles(p, distinct),
                                      std::vector<double>{});
}

/// The seam fixtures: degrees 1, 2 and 3, polynomial and rational. Degree 5
/// is covered by Loop(5) in the seam test above.
std::vector<BSplineCurve> SeamCurves()
{
    std::vector<BSplineCurve> curves;
    for (const int p : {1, 2, 3}) {
        curves.push_back(SeamCurve(p, false));
        curves.push_back(SeamCurve(p, true));
    }
    return curves;
}

} // namespace

TEST_CASE("the seam of a periodic B-spline is right-continuous at every degree, bit for bit",
          "[curve][periodic][seam]")
{
    for (const BSplineCurve& c : SeamCurves()) {
        INFO("degree " << c.Degree << (c.Weights.empty() ? "" : ", rational"));
        REQUIRE(Stereon::IsPeriodic(c));
        REQUIRE(Stereon::Period(c) == 1.0);
        // The seam's own derivative of order p differs from side to side:
        // the fixture is one where the convention matters.
        if (c.Degree >= 2) {
            REQUIRE(Deriv(c, std::nextafter(1.0, 0.0), c.Degree) != Deriv(c, 0.0, c.Degree));
        }
        // At the end, at the start and a whole number of periods away, every
        // order is the start's, exactly: the parameters are dyadic, so moving
        // them by periods rounds nothing.
        for (const double t : {0.0, 0.375, 0.875}) {
            for (const double k : {-2.0, -1.0, 1.0, 2.0}) {
                INFO("t = " << t << ", shifted by " << k);
                REQUIRE(Eval(c, t + k) == Eval(c, t));
                for (int n = 1; n <= c.Degree + 1; ++n) {
                    REQUIRE(Deriv(c, t + k, n) == Deriv(c, t, n));
                }
            }
        }
        REQUIRE(Eval(c, 1.0) == Eval(c, 0.0));
        for (int n = 1; n <= c.Degree + 1; ++n) {
            REQUIRE(Deriv(c, 1.0, n) == Deriv(c, 0.0, n));
        }
        // Next to the seam on either side: the point is continuous, so the
        // neighbours of b and of a lie within rounding of each other.
        const double scale = 1e-12 * 8.0;
        REQUIRE(Distance(Eval(c, std::nextafter(1.0, 2.0)), Eval(c, std::nextafter(0.0, 1.0))) <=
                scale);
        REQUIRE(Distance(Eval(c, std::nextafter(0.0, -1.0)), Eval(c, std::nextafter(1.0, 0.0))) <=
                scale);
    }
}

TEST_CASE("the batch evaluators follow the seam's convention on mixed batches",
          "[curve][periodic][seam][batch]")
{
    // Inside, outside, on both ends and beside them, in no particular order.
    const std::vector<double> t{0.375, 1.0,    -2.0, 0.0, std::nextafter(1.0, 0.0), 3.0,
                                1.875, -0.125, 0.5,  2.0, std::nextafter(0.0, 1.0), 1.0,
                                0.75};
    for (const BSplineCurve& c : SeamCurves()) {
        INFO("degree " << c.Degree << (c.Weights.empty() ? "" : ", rational"));
        std::vector<Point3> all(t.size());
        Stereon::EvalAll(c, t, all);
        for (std::size_t i = 0; i < t.size(); ++i) {
            INFO("t = " << t[i]);
            REQUIRE(all[i] == Eval(c, t[i]));
        }
        // The fast overloads agree with each other bit for bit, and with Eval
        // within their documented bound (docs/nurbs.rst, *Fast batch
        // evaluation*), with S the largest pole times the weight ratio.
        std::vector<Point3> converting(t.size());
        std::vector<Point3> prepared(t.size());
        Stereon::EvalAllFast(c, t, converting);
        Stereon::EvalAllFast(Stereon::ToSpanPolynomials(c), t, prepared);
        REQUIRE(converting == prepared);
        double s = 0.0;
        for (const Point3& pole : c.Poles) {
            s = std::max(s, Length(Stereon::ToVec(pole)));
        }
        if (!c.Weights.empty()) {
            const auto [low, high] = std::ranges::minmax(c.Weights);
            s *= high / low;
        }
        const double bound = (8.0 * (c.Degree + 1) + std::ldexp(1.0, c.Degree - 1)) *
                             std::numeric_limits<double>::epsilon() * s;
        for (std::size_t i = 0; i < t.size(); ++i) {
            INFO("t = " << t[i]);
            REQUIRE(Distance(converting[i], Eval(c, t[i])) <= bound);
        }
    }
}

TEST_CASE("a periodic plane B-spline's seam follows the same convention", "[curve][periodic][seam]")
{
    for (const int p : {1, 2, 3}) {
        const BSplineCurve2d c = SeamCurve2d(p);
        INFO("degree " << p);
        REQUIRE(Stereon::IsPeriodic(c));
        REQUIRE(Eval(c, 1.0) == Eval(c, 0.0));
        REQUIRE(Eval(c, -2.0) == Eval(c, 0.0));
        for (int n = 1; n <= p + 1; ++n) {
            REQUIRE(Deriv(c, 1.0, n) == Deriv(c, 0.0, n));
            REQUIRE(Deriv(c, 2.375, n) == Deriv(c, 0.375, n));
        }
        const std::vector<double> t{1.0, 0.0, -1.0, 0.375};
        std::vector<Point2> all(t.size());
        Stereon::EvalAll(c, t, all);
        for (std::size_t i = 0; i < t.size(); ++i) {
            REQUIRE(all[i] == Eval(c, t[i]));
        }
    }
}

TEST_CASE("splitting at the seam leaves no sliver; elsewhere it splits at the image",
          "[curve][periodic][seam]")
{
    for (const BSplineCurve& c : SeamCurves()) {
        INFO("degree " << c.Degree << (c.Weights.empty() ? "" : ", rational"));
        // The seam and its images: nothing on one side, so a refusal, never a
        // piece of zero length.
        for (const double t : {0.0, 1.0, -1.0, 3.0}) {
            const auto refused = Stereon::Split(c, t);
            REQUIRE_FALSE(refused.has_value());
            REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
        }
        // An image of an interior parameter splits there, two pieces covering the domain.
        const auto pieces = Stereon::Split(c, 2.375);
        REQUIRE(pieces.has_value());
        REQUIRE(Stereon::Domain(pieces->first) == Interval{0.0, 0.375});
        REQUIRE(Stereon::Domain(pieces->second) == Interval{0.375, 1.0});
        REQUIRE(Distance(Eval(pieces->first, 0.375), Eval(c, 0.375)) <= 1e-12);
    }
}

TEST_CASE("a curve not in the periodic form keeps its domain's end as it is",
          "[curve][periodic][seam]")
{
    for (const BSplineCurve& periodic : SeamCurves()) {
        // The same closed curve, clamped: not periodic, so its end is the end
        // of a bounded domain, evaluated in the last span, and no parameter is
        // carried anywhere.
        const BSplineCurve c = Stereon::Clamped(periodic);
        INFO("degree " << c.Degree << (c.Weights.empty() ? "" : ", rational"));
        REQUIRE(Stereon::Period(c) == 0.0);
        const Interval domain = Stereon::Domain(c);
        const double b = domain.GetUpper();
        REQUIRE(Distance(Eval(c, b), Stereon::ToPoint(Stereon::ToVec(c.Poles.back()))) <= 1e-12);
        if (c.Degree >= 2 && c.Weights.empty()) {
            // The order-p derivative at b is the last span's, not the first's:
            // a polynomial's is constant on a span, so the bits agree.
            REQUIRE(Deriv(c, b, c.Degree) == Deriv(c, std::nextafter(b, 0.0), c.Degree));
        }
        const std::vector<double> t{domain.GetLower(), b};
        std::vector<Point3> all(t.size());
        Stereon::EvalAll(c, t, all);
        REQUIRE(all[0] == Eval(c, t[0]));
        REQUIRE(all[1] == Eval(c, t[1]));
    }
}
