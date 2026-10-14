// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Curve fitting. Interpolation must pass through every point and reproduce
// any curve of its own space exactly; least squares must satisfy its
// normal equations, the first-order optimality conditions; approximation
// must meet its budget.

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

/// Points on a wavy space curve, unevenly spaced.
std::vector<Point3> RandomPoints(const std::size_t count)
{
    std::vector<Point3> points;
    double s = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        s += *rc::gen::inRange(5, 40) / 100.0;
        points.push_back(
            {s, std::sin(2.0 * s) + *rc::gen::inRange(-10, 10) / 100.0, std::cos(3.0 * s)});
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

} // namespace

TEST_CASE("interpolation passes through every point", "[curve][fitting][property]")
{
    rc::prop("for every degree and parametrisation", [] {
        const int degree = *rc::gen::inRange(1, 6);
        const std::vector<Point3> points = RandomPoints(static_cast<std::size_t>(degree) + 1 +
                                                        *rc::gen::inRange<std::size_t>(0, 40));
        const auto fit = Stereon::Interpolate(points, degree, RandomParametrization());
        RC_ASSERT(fit.has_value());
        RC_ASSERT((Domain(fit->Curve) == Interval{0.0, 1.0}));
        RC_ASSERT(fit->Parameters.front() == 0.0 && fit->Parameters.back() == 1.0);
        RC_ASSERT(std::ranges::is_sorted(fit->Parameters));
        RC_ASSERT(fit->MaxDeviation <= 1e-11 * Scale(points));
        RC_ASSERT(Eval(fit->Curve, 0.0) == points.front());
        RC_ASSERT(Distance(Eval(fit->Curve, 1.0), points.back()) <= 1e-12 * Scale(points));
    });
}

TEST_CASE("interpolation reproduces a curve of its own space", "[curve][fitting]")
{
    // A cubic polynomial in t lies in every cubic spline space on [0, 1],
    // so interpolating its points at their own parameters returns it.
    const auto cubic = [](const double t) { return Point3{t, t * t - t, 2.0 * t * t * t}; };
    std::vector<Point3> points;
    for (int k = 0; k <= 12; ++k) {
        points.push_back(cubic(k / 12.0));
    }
    const CurveFit fit = *Stereon::Interpolate(points, 3, Parametrization::Uniform);
    for (int i = 0; i <= 100; ++i) {
        const double s = i / 100.0;
        REQUIRE(Distance(Eval(fit.Curve, s), cubic(s)) <= 1e-12);
    }
}

TEST_CASE("least squares satisfies its normal equations", "[curve][fitting][property]")
{
    rc::prop("the residual is orthogonal to every interior basis function", [] {
        const int degree = *rc::gen::inRange(1, 5);
        const std::vector<Point3> points =
            RandomPoints(static_cast<std::size_t>(*rc::gen::inRange(10, 60)));
        const std::size_t poles =
            static_cast<std::size_t>(degree) + 1 +
            *rc::gen::inRange<std::size_t>(0, points.size() - static_cast<std::size_t>(degree) - 1);
        const auto fit = Stereon::FitLeastSquares(points, degree, poles, RandomParametrization());
        RC_ASSERT(fit.has_value());
        const BSplineCurve& c = fit->Curve;
        RC_ASSERT(c.Poles.size() == poles);
        RC_ASSERT(c.Poles.front() == points.front());
        RC_ASSERT(c.Poles.back() == points.back());
        // Gradient of the squared residual with respect to interior pole i:
        // sum over interior points of N_i(t_k) (C(t_k) - Q_k).
        std::vector<Vec3> gradient(poles, Vec3{});
        const auto p = static_cast<std::size_t>(degree);
        std::vector<double> basis(p + 1);
        for (std::size_t k = 1; k + 1 < points.size(); ++k) {
            const double t = fit->Parameters[k];
            const std::size_t span = Stereon::FindSpan(c.Knots, degree, poles, t);
            Stereon::BasisFunctions(c.Knots, degree, span, t, basis);
            for (std::size_t a = 0; a <= p; ++a) {
                gradient[span - p + a] += basis[a] * (Eval(c, t) - points[k]);
            }
        }
        for (std::size_t i = 1; i + 1 < poles; ++i) {
            RC_ASSERT(Length(gradient[i]) <= 1e-9 * Scale(points));
        }
    });
}

TEST_CASE("least squares reproduces data from its own space", "[curve][fitting]")
{
    // Points on a straight line: any degree and pole count fits them exactly.
    std::vector<Point3> points;
    for (int k = 0; k <= 20; ++k) {
        const double s = k * k / 400.0;
        points.push_back({1.0 + 2.0 * s, -s, 0.5 * s});
    }
    for (const std::size_t poles : {std::size_t{4}, std::size_t{7}, std::size_t{12}}) {
        const CurveFit fit = *Stereon::FitLeastSquares(points, 3, poles);
        REQUIRE(fit.MaxDeviation <= 1e-12);
    }
}

TEST_CASE("approximation meets its budget with few poles", "[curve][fitting][property]")
{
    rc::prop("the deviation is within the budget", [] {
        const std::vector<Point3> points =
            RandomPoints(static_cast<std::size_t>(*rc::gen::inRange(10, 80)));
        const double budget = *rc::gen::inRange(1, 300) / 1000.0;
        const auto fit = Stereon::Approximate(points, 3, ApproximationBudget{budget});
        RC_ASSERT(fit.has_value());
        RC_ASSERT(fit->MaxDeviation <= budget);
        RC_ASSERT(fit->Curve.Poles.size() <= points.size());
    });
    // A generous budget takes the fewest poles, degree + 1.
    std::vector<Point3> line;
    for (int k = 0; k <= 30; ++k) {
        line.push_back({k * 0.1, 0.0, 0.0});
    }
    REQUIRE(Stereon::Approximate(line, 3, ApproximationBudget{1e-6})->Curve.Poles.size() == 4);
}

TEST_CASE("planar fitting works in the plane", "[curve][fitting]")
{
    std::vector<Point2> points;
    for (int k = 0; k <= 16; ++k) {
        const double a = k * std::numbers::pi / 16.0;
        points.push_back({std::cos(a), std::sin(a)});
    }
    const auto fit = Stereon::Interpolate(std::span<const Point2>{points}, 3);
    REQUIRE(fit->MaxDeviation <= 1e-12);
    const auto approx =
        Stereon::Approximate(std::span<const Point2>{points}, 3, ApproximationBudget{1e-3});
    REQUIRE(approx->MaxDeviation <= 1e-3);
    const auto lsq = Stereon::FitLeastSquares(std::span<const Point2>{points}, 2, 6);
    REQUIRE(lsq->Curve.Poles.front() == points.front());
}

TEST_CASE("fitting validates data from outside the kernel", "[curve][fitting][factory]")
{
    const std::vector<Point3> few{Point3{}, Point3{1.0, 0.0, 0.0}};
    const std::vector<Point3> repeated{Point3{}, Point3{1.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0},
                                       Point3{2.0, 1.0, 0.0}};
    const std::vector<Point3> nan{Point3{}, Point3{NAN, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}};
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::Interpolate(few, 3).error().GetCode() == invalid);
    REQUIRE(Stereon::Interpolate(few, 0).error().GetCode() == invalid);
    REQUIRE(Stereon::Interpolate(repeated, 2).error().GetCode() == invalid);
    REQUIRE(Stereon::Interpolate(repeated, 2, Parametrization::Uniform).has_value());
    REQUIRE(Stereon::Interpolate(nan, 1).error().GetCode() == invalid);
    REQUIRE(Stereon::FitLeastSquares(repeated, 1, 9, Parametrization::Uniform).error().GetCode() ==
            invalid);
    REQUIRE(Stereon::Approximate(repeated, 1, ApproximationBudget{0.0}, Parametrization::Uniform)
                .error()
                .GetCode() == invalid);
    REQUIRE(Stereon::Approximate(repeated, 1, ApproximationBudget{NAN}, Parametrization::Uniform)
                .error()
                .GetCode() == invalid);
}
