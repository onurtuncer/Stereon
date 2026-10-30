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
#include <numbers>
#include <span>
#include <string>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

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
