// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Conversions from the exchange formats' representations: IGES 112's
// piecewise polynomials must become B-splines that evaluate as the
// polynomials do, and IGES 104's implicit conics must become the conic whose
// points satisfy the equation, of the right type, with the right frame.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::Curve2d;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::Hyperbola2d;
using Stereon::ImplicitConic;
using Stereon::Interval;
using Stereon::Parabola2d;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double Pi = std::numbers::pi;

// --- Piecewise polynomials ---------------------------------------------------

/// p(t) = c0 + c1 t + c2 t^2 + c3 t^3.
struct Cubic {
    std::vector<Vec3> C;

    [[nodiscard]] Point3 At(const double t) const
    {
        return Point3{} + C[0] + t * C[1] + (t * t) * C[2] + (t * t * t) * C[3];
    }

    /// The coefficients of p in s = t - t0, by Taylor expansion at t0.
    [[nodiscard]] std::vector<Vec3> About(const double t0) const
    {
        return {Stereon::ToVec(At(t0)), C[1] + (2.0 * t0) * C[2] + (3.0 * t0 * t0) * C[3],
                C[2] + (3.0 * t0) * C[3], C[3]};
    }
};

// --- Implicit conics -----------------------------------------------------------

/// The implicit equation of alpha u^2 + beta v^2 + gamma u + delta v + eps = 0
/// in the frame u = X.(p - c), v = Y.(p - c), X at angle phi, Y a quarter
/// turn on.
ImplicitConic Equation(const Point2 c, const double phi, const double alpha, const double beta,
                       const double gamma, const double delta, const double eps)
{
    const Vec2 x{std::cos(phi), std::sin(phi)};
    const Vec2 y{-x.Y, x.X};
    const double u0 = -Dot(x, Stereon::ToVec(c));
    const double v0 = -Dot(y, Stereon::ToVec(c));
    ImplicitConic k;
    k.A = alpha * x.X * x.X + beta * y.X * y.X;
    k.B = 2.0 * (alpha * x.X * x.Y + beta * y.X * y.Y);
    k.C = alpha * x.Y * x.Y + beta * y.Y * y.Y;
    k.D = 2.0 * alpha * u0 * x.X + 2.0 * beta * v0 * y.X + gamma * x.X + delta * y.X;
    k.E = 2.0 * alpha * u0 * x.Y + 2.0 * beta * v0 * y.Y + gamma * x.Y + delta * y.Y;
    k.F = alpha * u0 * u0 + beta * v0 * v0 + gamma * u0 + delta * v0 + eps;
    return k;
}

ImplicitConic Scaled(ImplicitConic k, const double s)
{
    for (double* c : {&k.A, &k.B, &k.C, &k.D, &k.E, &k.F}) {
        *c *= s;
    }
    return k;
}

/// The largest residual of the equation over points of the conic.
double WorstResidual(const ImplicitConic& k, const Curve2d& conic, const Interval range)
{
    double worst = 0.0;
    for (int i = 0; i <= 64; ++i) {
        const double t = range.GetLower() + (range.GetUpper() - range.GetLower()) * i / 64.0;
        worst = std::max(worst, std::abs(Stereon::ImplicitResidual(k, Eval(conic, t))));
    }
    return worst;
}

ErrorCode CodeOf(const ImplicitConic& k)
{
    return Stereon::ConicFromImplicit(k).error().GetCode();
}

} // namespace

TEST_CASE("a piecewise polynomial converts to the B-spline it is", "[curve][interchange]")
{
    const Cubic p{
        {Vec3{1.0, -2.0, 0.5}, Vec3{0.5, 1.0, -1.0}, Vec3{-0.3, 0.2, 0.4}, Vec3{0.05, -0.1, 0.02}}};
    const std::vector<double> breaks{-1.0, 0.25, 1.5, 4.0};
    std::vector<std::vector<Vec3>> pieces;
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        pieces.push_back(p.About(breaks[i]));
    }
    const Stereon::PiecewiseConversion result =
        *Stereon::PiecewisePolynomialToBSpline(breaks, pieces);
    REQUIRE(result.Curve.Degree == 3);
    REQUIRE(result.Curve.Poles.size() == 10);
    REQUIRE(Domain(result.Curve) == Interval{-1.0, 4.0});
    REQUIRE(result.MaxJoinGap <= 1e-13);
    for (int i = 0; i <= 100; ++i) {
        const double t = -1.0 + 5.0 * i / 100.0;
        REQUIRE(Distance(Eval(result.Curve, t), p.At(t)) <= 1e-12);
    }
}

TEST_CASE("pieces of lower degree are raised and gaps are reported", "[curve][interchange]")
{
    // A line from (0, 0, 0) to (1, 0, 0) on [0, 1], then a parabola from
    // (1, 0.5, 0) on [1, 3]: the join misses by 0.5.
    const std::vector<double> breaks{0.0, 1.0, 3.0};
    const std::vector<std::vector<Vec3>> pieces{
        {Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}},
        {Vec3{1.0, 0.5, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}},
    };
    const Stereon::PiecewiseConversion result =
        *Stereon::PiecewisePolynomialToBSpline(breaks, pieces);
    REQUIRE(result.Curve.Degree == 2);
    REQUIRE(std::abs(result.MaxJoinGap - 0.5) <= 1e-15);
    REQUIRE(Distance(Eval(result.Curve, 1.0), Point3{1.0, 0.25, 0.0}) <= 1e-15);
    REQUIRE(Distance(Eval(result.Curve, 3.0), Point3{1.0, 2.5, 4.0}) <= 1e-14);
    // Halfway along each piece the moved join pole weighs 1/4: a quarter
    // of its 0.25 shift.
    REQUIRE(Distance(Eval(result.Curve, 0.5), Point3{0.5, 0.0625, 0.0}) <= 1e-15);
    REQUIRE(Distance(Eval(result.Curve, 2.0), Point3{1.0, 1.4375, 1.0}) <= 1e-14);
}

TEST_CASE("piecewise polynomials convert in the plane", "[curve][interchange][2d]")
{
    const std::vector<double> breaks{2.0, 5.0};
    const std::vector<std::vector<Vec2>> pieces{{Vec2{1.0, 1.0}, Vec2{0.0, 2.0}, Vec2{1.0, 0.0}}};
    const Stereon::PiecewiseConversion2d result =
        *Stereon::PiecewisePolynomialToBSpline(breaks, pieces);
    REQUIRE(Distance(Eval(result.Curve, 4.0), Point2{5.0, 5.0}) <= 1e-14);
}

TEST_CASE("PiecewisePolynomialToBSpline validates its input", "[curve][interchange][validation]")
{
    const auto code = [](const std::vector<double>& breaks,
                         const std::vector<std::vector<Vec3>>& pieces) {
        return Stereon::PiecewisePolynomialToBSpline(breaks, pieces).error().GetCode();
    };
    const std::vector<Vec3> line{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(code({0.0}, {}) == ErrorCode::InvalidInput);
    REQUIRE(code({0.0, 1.0, 2.0}, {line}) == ErrorCode::InvalidInput);
    REQUIRE(code({0.0, 0.0}, {line}) == ErrorCode::InvalidInput);
    REQUIRE(code({1.0, 0.0}, {line}) == ErrorCode::InvalidInput);
    REQUIRE(code({0.0, NaN}, {line}) == ErrorCode::InvalidInput);
    REQUIRE(code({0.0, 1.0}, {{}}) == ErrorCode::InvalidInput);
    REQUIRE(code({0.0, 1.0}, {{Vec3{NaN, 0.0, 0.0}}}) == ErrorCode::InvalidInput);
    // A constant piece is a degree-1 curve standing still.
    REQUIRE(Stereon::PiecewisePolynomialToBSpline(
                std::vector{0.0, 1.0}, std::vector<std::vector<Vec3>>{{Vec3{1.0, 2.0, 3.0}}})
                ->Curve.Degree == 1);
}

TEST_CASE("an implicit ellipse becomes its ellipse", "[curve][interchange][conic]")
{
    const Point2 c{2.0, -1.0};
    const double phi = 0.7;
    const ImplicitConic k = Equation(c, phi, 1.0 / 9.0, 1.0 / 4.0, 0.0, 0.0, -1.0);
    for (const double s : {1.0, -36.0, 1e-6, 1e7}) {
        const Curve2d conic = *Stereon::ConicFromImplicit(Scaled(k, s));
        const Ellipse2d* e = conic.Target<Ellipse2d>();
        REQUIRE(e != nullptr);
        REQUIRE(Distance(e->Center, c) <= 1e-12);
        REQUIRE(std::abs(e->MajorRadius - 3.0) <= 1e-12);
        REQUIRE(std::abs(e->MinorRadius - 2.0) <= 1e-12);
        REQUIRE(std::abs(std::abs(Dot(e->XAxis, Vec2{std::cos(phi), std::sin(phi)})) - 1.0) <=
                1e-12);
        // Counter-clockwise.
        REQUIRE(Cross(e->XAxis, e->YAxis) > 0.0);
        REQUIRE(WorstResidual(k, conic, Interval{0.0, 2.0 * Pi}) <= 1e-14);
    }
}

TEST_CASE("implicit conics round-trip through their equations",
          "[curve][interchange][conic][property]")
{
    rc::prop("the converted conic satisfies the equation",
             [](const double cx, const double cy, const double angle, const unsigned ra,
                const unsigned rb, const unsigned kind) {
                 const Point2 c{std::fmod(cx, 50.0), std::fmod(cy, 50.0)};
                 const double phi = std::fmod(angle, 2.0 * Pi);
                 const double a = 0.5 + (ra % 1000) / 100.0;
                 const double b = 0.5 + (rb % 1000) / 100.0;
                 ImplicitConic k;
                 Interval range{0.0, 2.0 * Pi};
                 switch (kind % 3) {
                 case 0:
                     k = Equation(c, phi, 1.0 / (a * a), 1.0 / (b * b), 0.0, 0.0, -1.0);
                     break;
                 case 1:
                     k = Equation(c, phi, 1.0 / (a * a), -1.0 / (b * b), 0.0, 0.0, -1.0);
                     range = Interval{-2.0, 2.0};
                     break;
                 default:
                     k = Equation(c, phi, 0.0, -1.0 / (4.0 * a), 1.0, 0.0, 0.0);
                     range = Interval{-3.0 * a, 3.0 * a};
                     break;
                 }
                 const auto conic = Stereon::ConicFromImplicit(k);
                 RC_ASSERT(conic.has_value());
                 RC_ASSERT((kind % 3 != 0 || conic->Target<Ellipse2d>() != nullptr));
                 RC_ASSERT((kind % 3 != 1 || conic->Target<Hyperbola2d>() != nullptr));
                 RC_ASSERT((kind % 3 != 2 || conic->Target<Parabola2d>() != nullptr));
                 RC_ASSERT(WorstResidual(k, *conic, range) <= 1e-9);
             });
}

TEST_CASE("an implicit hyperbola takes the branch it is asked for", "[curve][interchange][conic]")
{
    // x^2 / 4 - y^2 = 1 about (1, 1): branches through (3, 1) and (-1, 1).
    const ImplicitConic k = Equation(Point2{1.0, 1.0}, 0.0, 0.25, -1.0, 0.0, 0.0, -1.0);
    const Curve2d right = *Stereon::ConicFromImplicit(k, Point2{10.0, 1.0});
    const Curve2d left = *Stereon::ConicFromImplicit(k, Point2{-10.0, 3.0});
    REQUIRE(Distance(Eval(right, 0.0), Point2{3.0, 1.0}) <= 1e-14);
    REQUIRE(Distance(Eval(left, 0.0), Point2{-1.0, 1.0}) <= 1e-14);
    const Hyperbola2d* h = right.Target<Hyperbola2d>();
    REQUIRE(h != nullptr);
    REQUIRE(std::abs(h->MajorRadius - 2.0) <= 1e-14);
    REQUIRE(std::abs(h->MinorRadius - 1.0) <= 1e-14);

    // y^2 - x^2 / 4 = 1: the transverse axis is y.
    const ImplicitConic upright = Equation(Point2{0.0, 0.0}, 0.0, -0.25, 1.0, 0.0, 0.0, -1.0);
    const Curve2d up = *Stereon::ConicFromImplicit(upright, Point2{0.0, 5.0});
    REQUIRE(Distance(Eval(up, 0.0), Point2{0.0, 1.0}) <= 1e-14);
}

TEST_CASE("an implicit parabola becomes its parabola", "[curve][interchange][conic]")
{
    // y = x^2 / 8 + 1: vertex (0, 1), focal length 2, opening up.
    const ImplicitConic k{1.0, 0.0, 0.0, 0.0, -8.0, 8.0};
    const Curve2d conic = *Stereon::ConicFromImplicit(k);
    const Parabola2d* p = conic.Target<Parabola2d>();
    REQUIRE(p != nullptr);
    REQUIRE(Distance(p->Vertex, Point2{0.0, 1.0}) <= 1e-14);
    REQUIRE(std::abs(p->FocalLength - 2.0) <= 1e-14);
    REQUIRE(Distance(Point2{} + p->XAxis, Point2{0.0, 1.0}) <= 1e-14);

    // An exporter's rounding leaves B^2 - 4AC a little off zero; it is
    // still a parabola.
    const ImplicitConic rotated = Equation(Point2{3.0, -2.0}, 0.3, 0.0, -0.125, 1.0, 0.0, 0.0);
    ImplicitConic rounded = rotated;
    rounded.A *= 1.0 + 1e-13;
    const Curve2d nearly = *Stereon::ConicFromImplicit(rounded);
    REQUIRE(nearly.Target<Parabola2d>() != nullptr);
    REQUIRE(WorstResidual(rotated, nearly, Interval{-5.0, 5.0}) <= 1e-11);
}

TEST_CASE("ConicFromImplicit rejects equations that are not conics",
          "[curve][interchange][conic][validation]")
{
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(CodeOf({1.0, 0.0, 1.0, 0.0, 0.0, 1.0}) == ErrorCode::InvalidInput);  // no point
    REQUIRE(CodeOf({1.0, 0.0, 1.0, 0.0, 0.0, 0.0}) == ErrorCode::InvalidInput);  // one point
    REQUIRE(CodeOf({1.0, 0.0, -1.0, 0.0, 0.0, 0.0}) == ErrorCode::InvalidInput); // crossing
    REQUIRE(CodeOf({1.0, 0.0, 0.0, 0.0, 0.0, -1.0}) == ErrorCode::InvalidInput); // parallel
    REQUIRE(CodeOf({0.0, 0.0, 0.0, 1.0, 2.0, 3.0}) == ErrorCode::InvalidInput);  // a line
    REQUIRE(CodeOf({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) == ErrorCode::InvalidInput);  // the plane
    REQUIRE(CodeOf({NaN, 0.0, 1.0, 0.0, 0.0, -1.0}) == ErrorCode::InvalidInput);
    REQUIRE(Stereon::ConicFromImplicit({1.0, 0.0, 1.0, 0.0, 0.0, -1.0}, Point2{NaN, 0.0})
                .error()
                .GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Stereon::ImplicitResidual({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Point2{1.0, 1.0}) == 0.0);
}
