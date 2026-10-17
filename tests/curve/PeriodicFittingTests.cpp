// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Fitting closed, periodic curves. A periodic interpolant must pass through
// every point and close with the smoothness of its degree at the seam, and
// reproduce a periodic spline of its own space exactly; periodic least
// squares must satisfy its normal equations over every pole; periodic
// approximation must meet its budget.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::ApproximationBudget;
using Stereon::BSplineCurve;
using Stereon::CurveFit;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Parametrization;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec3;

namespace {

/// Points around a bumpy loop, unevenly spaced and never coinciding, that
/// do not return to the first point: the closing chord is a gap like the
/// others.
std::vector<Point3> RandomLoop(const std::size_t count)
{
    std::vector<double> steps;
    double total = 0.0;
    for (std::size_t k = 0; k <= count; ++k) {
        steps.push_back(*rc::gen::inRange(5, 40) / 100.0);
        total += steps.back();
    }
    std::vector<Point3> points;
    double s = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        const double angle = 2.0 * std::numbers::pi * s / total;
        const double radius = 1.0 + 0.3 * std::cos(3.0 * angle) + *rc::gen::inRange(-5, 5) / 100.0;
        points.push_back(
            {radius * std::cos(angle), radius * std::sin(angle), 0.2 * std::sin(2.0 * angle)});
        s += steps[k];
    }
    return points;
}

/// `count` points evenly spaced on the unit circle in the plane z = 0.
std::vector<Point3> CirclePoints(const std::size_t count)
{
    std::vector<Point3> points;
    for (std::size_t k = 0; k < count; ++k) {
        const double a =
            2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(count);
        points.push_back({std::cos(a), std::sin(a), 0.0});
    }
    return points;
}

Parametrization RandomParametrization()
{
    const std::array choices{Parametrization::Uniform, Parametrization::ChordLength,
                             Parametrization::Centripetal};
    return choices[*rc::gen::inRange<std::size_t>(0, choices.size())];
}

double Scale(const std::span<const Point3> points)
{
    double scale = 1.0;
    for (const Point3& p : points) {
        scale = std::max({scale, std::abs(p.X), std::abs(p.Y), std::abs(p.Z)});
    }
    return scale;
}

/// The curve is one period of a periodic spline: it has the wrapped poles
/// and the periodic knots of Detail::WrappedPoles and
/// Detail::PeriodicKnotVector, and its derivatives up to order p - 1 agree
/// at the two ends of the domain.
void RequirePeriodic(const BSplineCurve& c, const std::size_t distinctPoles)
{
    const auto p = static_cast<std::size_t>(c.Degree);
    RC_ASSERT((Domain(c) == Interval{0.0, 1.0}));
    RC_ASSERT(c.Poles.size() == distinctPoles + p);
    RC_ASSERT(c.Knots.size() == distinctPoles + 2 * p + 1);
    for (std::size_t i = 0; i < p; ++i) {
        RC_ASSERT(c.Poles[distinctPoles + i] == c.Poles[i]);
    }
    for (std::size_t i = 0; i + distinctPoles < c.Knots.size(); ++i) {
        RC_ASSERT(std::abs(c.Knots[i + distinctPoles] - c.Knots[i] - 1.0) <= 1e-15);
    }
    RC_ASSERT(Stereon::MakeBSplineCurve(c.Degree, c.Knots, c.Poles).has_value());
    const int order = c.Degree - 1;
    std::vector<Vec3> start(p);
    std::vector<Vec3> end(p);
    Derivatives(c, 0.0, order, std::span{start});
    Derivatives(c, 1.0, order, std::span{end});
    for (std::size_t k = 0; k < p; ++k) {
        RC_ASSERT(Length(start[k] - end[k]) <= 1e-9 * (1.0 + Length(start[k])));
    }
}

} // namespace

TEST_CASE("periodic interpolation closes through every point", "[curve][fitting][property]")
{
    rc::prop("for every degree and parametrisation", [] {
        const int degree = *rc::gen::inRange(1, 7);
        const std::size_t count =
            static_cast<std::size_t>(degree) + 1 + *rc::gen::inRange<std::size_t>(0, 40);
        const std::vector<Point3> points = RandomLoop(count);
        const auto fit = Stereon::InterpolatePeriodic(points, degree, RandomParametrization());
        RC_ASSERT(fit.has_value());
        RequirePeriodic(fit->Curve, count);
        RC_ASSERT(fit->Parameters.size() == count);
        RC_ASSERT(std::ranges::is_sorted(fit->Parameters));
        RC_ASSERT(std::ranges::adjacent_find(fit->Parameters) == fit->Parameters.end());
        RC_ASSERT(fit->Parameters.front() >= 0.0 && fit->Parameters.back() < 1.0);
        RC_ASSERT((fit->Parameters.front() == 0.0) == (degree % 2 == 1));
        RC_ASSERT(fit->MaxDeviation <= 1e-10 * Scale(points));
        for (std::size_t k = 0; k < count; ++k) {
            RC_ASSERT(Distance(Eval(fit->Curve, fit->Parameters[k]), points[k]) <=
                      1e-10 * Scale(points));
        }
    });
}

TEST_CASE("periodic interpolation reproduces a periodic spline of its own space",
          "[curve][fitting][property]")
{
    // A periodic spline with uniform knots, sampled at its knots (odd
    // degree) or halfway between them (even degree), is the unique periodic
    // interpolant of those samples under uniform parametrisation.
    rc::prop("for every degree", [] {
        const int degree = *rc::gen::inRange(1, 7);
        const auto p = static_cast<std::size_t>(degree);
        const std::size_t n = p + 1 + *rc::gen::inRange<std::size_t>(0, 25);
        std::vector<Point3> poles;
        std::vector<double> periodKnots;
        for (std::size_t i = 0; i <= n; ++i) {
            periodKnots.push_back(static_cast<double>(i) / static_cast<double>(n));
            if (i < n) {
                poles.push_back({*rc::gen::inRange(-20, 20) / 10.0,
                                 *rc::gen::inRange(-20, 20) / 10.0,
                                 *rc::gen::inRange(-20, 20) / 10.0});
            }
        }
        const BSplineCurve reference{degree,
                                     Stereon::Detail::PeriodicKnotVector(periodKnots, degree),
                                     Stereon::Detail::WrappedPoles(poles, degree),
                                     {}};
        std::vector<Point3> samples;
        for (std::size_t k = 0; k < n; ++k) {
            const double t =
                (static_cast<double>(k) + (degree % 2 == 0 ? 0.5 : 0.0)) / static_cast<double>(n);
            samples.push_back(Eval(reference, t));
        }
        const auto fit = Stereon::InterpolatePeriodic(samples, degree, Parametrization::Uniform);
        RC_ASSERT(fit.has_value());
        const double scale = Scale(poles);
        for (int i = 0; i <= 200; ++i) {
            const double s = i / 200.0;
            RC_ASSERT(Distance(Eval(fit->Curve, s), Eval(reference, s)) <= 1e-9 * scale);
        }
    });
}

TEST_CASE("an even degree interpolates an even number of points", "[curve][fitting]")
{
    // With knots at the data, a periodic quadratic through eight evenly
    // spaced points would be singular; the midpoint knots make it a
    // well-conditioned fit whose poles stay near the circle.
    const std::vector<Point3> points = CirclePoints(8);
    for (const int degree : {2, 4}) {
        const CurveFit fit = *Stereon::InterpolatePeriodic(points, degree);
        REQUIRE(fit.MaxDeviation <= 1e-12);
        for (const Point3& pole : fit.Curve.Poles) {
            REQUIRE(Distance(pole, Point3{}) <= 1.5);
        }
    }
}

TEST_CASE("periodic least squares satisfies its normal equations", "[curve][fitting][property]")
{
    rc::prop("the residual is orthogonal to every periodic basis function", [] {
        const int degree = *rc::gen::inRange(1, 6);
        const auto p = static_cast<std::size_t>(degree);
        const std::vector<Point3> points =
            RandomLoop(static_cast<std::size_t>(*rc::gen::inRange(10, 60)));
        const std::size_t n = p + 1 + *rc::gen::inRange<std::size_t>(0, points.size() - p - 1);
        const auto fit =
            Stereon::FitLeastSquaresPeriodic(points, degree, n, RandomParametrization());
        RC_ASSERT(fit.has_value());
        const BSplineCurve& c = fit->Curve;
        RequirePeriodic(c, n);
        // Gradient of the squared residual with respect to distinct pole i:
        // the sum over the points of B_i(t_k) (C(t_k) - Q_k), B_i the
        // periodic basis function that multiplies pole i and its repeat.
        std::vector<Vec3> gradient(n, Vec3{});
        std::vector<double> basis(p + 1);
        for (std::size_t k = 0; k < points.size(); ++k) {
            const double t = fit->Parameters[k];
            const std::size_t span = Stereon::FindSpan(c.Knots, degree, c.Poles.size(), t);
            Stereon::BasisFunctions(c.Knots, degree, span, t, basis);
            for (std::size_t a = 0; a <= p; ++a) {
                gradient[(span - p + a) % n] += basis[a] * (Eval(c, t) - points[k]);
            }
        }
        const double scale = std::max(Scale(points), Scale(c.Poles));
        for (std::size_t i = 0; i < n; ++i) {
            RC_ASSERT(Length(gradient[i]) <= 1e-9 * scale);
        }
    });
}

TEST_CASE("periodic approximation meets its budget with few poles", "[curve][fitting][property]")
{
    rc::prop("the deviation is within the budget", [] {
        const std::vector<Point3> points =
            RandomLoop(static_cast<std::size_t>(*rc::gen::inRange(10, 80)));
        const double budget = *rc::gen::inRange(1, 300) / 1000.0;
        const auto fit = Stereon::ApproximatePeriodic(points, 3, ApproximationBudget{budget});
        RC_ASSERT(fit.has_value());
        RC_ASSERT(fit->MaxDeviation <= budget);
        RC_ASSERT(fit->Curve.Poles.size() <= points.size() + 3);
    });
    // A finely sampled circle needs a handful of poles for a thousandth of
    // a millimetre: the error of a cubic falls with the fourth power of the
    // span.
    const CurveFit fit =
        *Stereon::ApproximatePeriodic(CirclePoints(40), 3, ApproximationBudget{1e-3});
    REQUIRE(fit.MaxDeviation <= 1e-3);
    REQUIRE(fit.Curve.Poles.size() - 3 <= 12);
}

TEST_CASE("planar periodic fitting works in the plane", "[curve][fitting]")
{
    std::vector<Point2> points;
    for (int k = 0; k < 12; ++k) {
        const double a = 2.0 * std::numbers::pi * k / 12.0;
        points.push_back({2.0 * std::cos(a), std::sin(a)});
    }
    const auto fit = Stereon::InterpolatePeriodic(std::span<const Point2>{points}, 3);
    REQUIRE(fit->MaxDeviation <= 1e-12);
    REQUIRE(Distance(Eval(fit->Curve, 0.0), Eval(fit->Curve, 1.0)) <= 1e-15);
    REQUIRE(fit->Curve.Poles.size() == 15);
    const auto approx =
        Stereon::ApproximatePeriodic(std::span<const Point2>{points}, 3, ApproximationBudget{1e-2});
    REQUIRE(approx->MaxDeviation <= 1e-2);
    const auto lsq = Stereon::FitLeastSquaresPeriodic(std::span<const Point2>{points}, 2, 6);
    REQUIRE(lsq->Curve.Poles.size() == 8);
}

TEST_CASE("periodic fitting validates data from outside the kernel", "[curve][fitting][factory]")
{
    const std::vector<Point3> few{Point3{}, Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}};
    const std::vector<Point3> returning{Point3{}, Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0},
                                        Point3{0.0, 1.0, 0.0}, Point3{}};
    const std::vector<Point3> nan{Point3{}, Point3{NAN, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}};
    const std::vector<Point3> square = CirclePoints(4);
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::InterpolatePeriodic(few, 3).error().GetCode() == invalid);
    REQUIRE(Stereon::InterpolatePeriodic(few, 0).error().GetCode() == invalid);
    REQUIRE(Stereon::InterpolatePeriodic(nan, 1).error().GetCode() == invalid);
    // The first point given again at the end is a zero closing chord.
    REQUIRE(Stereon::InterpolatePeriodic(returning, 2).error().GetCode() == invalid);
    REQUIRE(Stereon::InterpolatePeriodic(returning, 2, Parametrization::Uniform).has_value());
    REQUIRE(Stereon::FitLeastSquaresPeriodic(square, 2, 2).error().GetCode() == invalid);
    REQUIRE(Stereon::FitLeastSquaresPeriodic(square, 2, 5).error().GetCode() == invalid);
    REQUIRE(Stereon::FitLeastSquaresPeriodic(square, 2, 4).has_value());
    REQUIRE(Stereon::ApproximatePeriodic(square, 1, ApproximationBudget{0.0}).error().GetCode() ==
            invalid);
    REQUIRE(Stereon::ApproximatePeriodic(square, 1, ApproximationBudget{NAN}).error().GetCode() ==
            invalid);
}
