// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// B-spline and NURBS curves. The kernel's algorithms (A2.1-A2.3, A3.2,
// A4.2, A5.1, Bernstein root isolation) are checked against independent
// references written here the textbook way: the Cox-de Boor recursion for
// the basis, de Boor's algorithm for points, exact polynomials and circle
// identities for derivatives, and dense sampling for projection.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

import Stereon.Curve;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::Curve;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::Vec3;
using StereonTest::Dual;

namespace {

constexpr double Tolerance = 1e-12;

bool Near(const Point3& a, const Point3& b, const double tolerance = Tolerance)
{
    return Length(a - b) <= tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

bool Near(const Vec3& a, const Vec3& b, const double tolerance = Tolerance)
{
    return Length(a - b) <= tolerance * (1.0 + Length(b));
}

// ---------------------------------------------------------------------------
// Independent references

/// N_{i,p}(t) by the Cox-de Boor recursion, with 0/0 = 0, for t inside the
/// half-open span structure (not at the end of the domain).
double CoxDeBoor(const std::vector<double>& u, const std::size_t i, const int p, const double t)
{
    if (p == 0) {
        return (u[i] <= t && t < u[i + 1]) ? 1.0 : 0.0;
    }
    const auto q = static_cast<std::size_t>(p);
    const double a = u[i + q] - u[i];
    const double b = u[i + q + 1] - u[i + 1];
    const double left = a == 0.0 ? 0.0 : (t - u[i]) / a * CoxDeBoor(u, i, p - 1, t);
    const double right = b == 0.0 ? 0.0 : (u[i + q + 1] - t) / b * CoxDeBoor(u, i + 1, p - 1, t);
    return left + right;
}

/// The k-th derivative of N_{i,p} by the derivative recursion
/// N' = p/(u_{i+p}-u_i) N_{i,p-1} - p/(u_{i+p+1}-u_{i+1}) N_{i+1,p-1}.
double CoxDeBoorDerivative(const std::vector<double>& u, const std::size_t i, const int p,
                           const int k, const double t)
{
    if (k == 0) {
        return CoxDeBoor(u, i, p, t);
    }
    if (p == 0) {
        return 0.0;
    }
    const auto q = static_cast<std::size_t>(p);
    const double a = u[i + q] - u[i];
    const double b = u[i + q + 1] - u[i + 1];
    const double left = a == 0.0 ? 0.0 : p / a * CoxDeBoorDerivative(u, i, p - 1, k - 1, t);
    const double right = b == 0.0 ? 0.0 : p / b * CoxDeBoorDerivative(u, i + 1, p - 1, k - 1, t);
    return left - right;
}

/// The point by summing every basis function times its homogeneous pole.
Point3 ReferencePoint(const BSplineCurve& c, const double t)
{
    Vec3 numerator{};
    double denominator = 0.0;
    for (std::size_t i = 0; i < c.Poles.size(); ++i) {
        const double w = c.Weights.empty() ? 1.0 : c.Weights[i];
        const double n = CoxDeBoor(c.Knots, i, c.Degree, t);
        numerator += (n * w) * Stereon::ToVec(c.Poles[i]);
        denominator += n * w;
    }
    return Stereon::ToPoint(numerator / denominator);
}

// ---------------------------------------------------------------------------
// Generators

/// A random valid curve: degree 1 to 5, clamped or not, interior knots with
/// multiplicities up to the degree, optionally rational.
BSplineCurve RandomCurve(const bool rational)
{
    const int p = *rc::gen::inRange(1, 6);
    const auto q = static_cast<std::size_t>(p);
    const std::size_t poles = q + 1 + static_cast<std::size_t>(*rc::gen::inRange(0, 6));
    const bool clamped = *rc::gen::arbitrary<bool>();
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    // Clamped: the first and last p + 1 knots share their value.
    if (clamped) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        points.push_back({*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                          *rc::gen::inRange(-50, 50) / 10.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(2, 40) / 10.0);
        }
    }
    const auto made = Stereon::MakeBSplineCurve(p, knots, points, weights);
    RC_PRE(made.has_value()); // the generator may produce an empty domain
    return *made;
}

double RandomParameter(const BSplineCurve& c)
{
    const Interval d = Domain(c);
    const double s = *rc::gen::inRange(0, 1000) / 1000.0;
    return d.GetLower() + s * (d.GetUpper() - d.GetLower());
}

/// (t, t^2, t^3) on [0, 1] as a cubic Bezier curve: the monomial t^k has
/// Bernstein coefficients C(i, k) / C(3, k).
BSplineCurve MonomialCurve()
{
    return {3,
            {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
            {Point3{0.0, 0.0, 0.0}, Point3{1.0 / 3.0, 0.0, 0.0}, Point3{2.0 / 3.0, 1.0 / 3.0, 0.0},
             Point3{1.0, 1.0, 1.0}},
            {}};
}

/// The quarter of the unit circle in the XY plane as a rational quadratic
/// Bezier curve, the classic exact conic.
BSplineCurve QuarterCircle()
{
    return {2,
            {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
            {Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}, Point3{0.0, 1.0, 0.0}},
            {1.0, std::numbers::sqrt2 / 2.0, 1.0}};
}

} // namespace

TEST_CASE("B-spline curves satisfy the curve concepts and are shared when wrapped",
          "[curve][bspline]")
{
    STATIC_REQUIRE(Stereon::CurveGeometry<BSplineCurve>);
    STATIC_REQUIRE(Stereon::Curve2dGeometry<BSplineCurve2d>);
    STATIC_REQUIRE_FALSE(Stereon::Detail::FitsInline<BSplineCurve, Stereon::GeometryBufferSize>);
    const Curve curve(MonomialCurve());
    const Curve copy = curve; // NOLINT(performance-unnecessary-copy-initialization)
    REQUIRE(copy.Target<BSplineCurve>() == curve.Target<BSplineCurve>());
    REQUIRE(Eval(copy, 0.5) == Eval(MonomialCurve(), 0.5));
}

TEST_CASE("FindSpan finds the span of every parameter", "[curve][bspline][basis]")
{
    const std::vector<double> u{0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0};
    const std::size_t poles = u.size() - 3;
    REQUIRE(Stereon::FindSpan(u, 2, poles, 0.0) == 2);
    REQUIRE(Stereon::FindSpan(u, 2, poles, 0.5) == 2);
    REQUIRE(Stereon::FindSpan(u, 2, poles, 1.0) == 3);
    REQUIRE(Stereon::FindSpan(u, 2, poles, 2.0) == 5); // past the double knot
    REQUIRE(Stereon::FindSpan(u, 2, poles, 3.99) == 6);
    REQUIRE(Stereon::FindSpan(u, 2, poles, 4.0) == 6); // the end belongs to the last span
    REQUIRE(Stereon::KnotMultiplicity(u, 2.0) == 2);
    REQUIRE(Stereon::KnotMultiplicity(u, 4.0) == 3);
    REQUIRE(Stereon::KnotMultiplicity(u, 2.5) == 0);
}

TEST_CASE("basis functions agree with the Cox-de Boor recursion",
          "[curve][bspline][basis][property]")
{
    rc::prop("values, partition of unity and derivatives of every order", [] {
        const BSplineCurve c = RandomCurve(false);
        const double t = RandomParameter(c);
        RC_PRE(t < Domain(c).GetUpper());
        const auto p = static_cast<std::size_t>(c.Degree);
        const std::size_t span = Stereon::FindSpan(c.Knots, c.Degree, c.Poles.size(), t);
        const int order = c.Degree + 1;
        std::vector<double> ders((static_cast<std::size_t>(order) + 1) * (p + 1));
        Stereon::BasisDerivatives(c.Knots, c.Degree, span, t, order, ders);
        std::vector<double> values(p + 1);
        Stereon::BasisFunctions(c.Knots, c.Degree, span, t, values);
        double sum = 0.0;
        for (std::size_t j = 0; j <= p; ++j) {
            const std::size_t i = span - p + j;
            RC_ASSERT(values[j] >= 0.0);
            RC_ASSERT(std::abs(values[j] - CoxDeBoor(c.Knots, i, c.Degree, t)) <= 1e-13);
            RC_ASSERT(values[j] == ders[j]);
            for (int k = 1; k <= order; ++k) {
                const double expected = CoxDeBoorDerivative(c.Knots, i, c.Degree, k, t);
                const double actual = ders[static_cast<std::size_t>(k) * (p + 1) + j];
                RC_ASSERT(std::abs(actual - expected) <= 1e-9 * (1.0 + std::abs(expected)));
            }
            sum += values[j];
        }
        RC_ASSERT(std::abs(sum - 1.0) <= 1e-14);
    });
}

TEST_CASE("Eval agrees with the basis-function sum for polynomial and rational curves",
          "[curve][bspline][property]")
{
    rc::prop("every random curve", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double t = RandomParameter(c);
        RC_PRE(t < Domain(c).GetUpper());
        RC_ASSERT(Near(Eval(c, t), ReferencePoint(c, t), 1e-12));
        RC_ASSERT(Contains(BoundingBox(c), Eval(c, t)));
    });
}

TEST_CASE("clamped curves start and end at their end poles", "[curve][bspline]")
{
    const BSplineCurve c = QuarterCircle();
    REQUIRE(Eval(c, 0.0) == c.Poles.front());
    REQUIRE(Eval(c, 1.0) == c.Poles.back());
    REQUIRE(Domain(c) == Interval{0.0, 1.0});
    REQUIRE(Stereon::IsRational(c));
    REQUIRE_FALSE(Stereon::IsRational(MonomialCurve()));
}

TEST_CASE("derivatives of a polynomial curve are exact", "[curve][bspline]")
{
    // The monomial curve, also after knot insertion splits it into spans.
    const BSplineCurve one = MonomialCurve();
    const BSplineCurve many = *Stereon::InsertKnot(*Stereon::InsertKnot(one, 0.3), 0.7, 2);
    for (const BSplineCurve& c : {one, many}) {
        for (int i = 0; i <= 20; ++i) {
            const double t = i / 20.0;
            REQUIRE(Near(Eval(c, t), Point3{t, t * t, t * t * t}));
            REQUIRE(Near(Deriv(c, t, 1), Vec3{1.0, 2.0 * t, 3.0 * t * t}));
            REQUIRE(Near(Deriv(c, t, 2), Vec3{0.0, 2.0, 6.0 * t}));
            REQUIRE(Near(Deriv(c, t, 3), Vec3{0.0, 0.0, 6.0}));
            REQUIRE(Deriv(c, t, 4) == Vec3{});
        }
    }
}

TEST_CASE("rational derivatives satisfy the circle's identities", "[curve][bspline][nurbs]")
{
    // On the unit circle |C| = 1, so C.C' = 0, C.C'' + |C'|^2 = 0 and
    // C.C''' + 3 C'.C'' = 0: every order of the quotient rule is checked.
    const BSplineCurve c = QuarterCircle();
    for (int i = 0; i <= 20; ++i) {
        const double t = i / 20.0;
        const Vec3 p = Stereon::ToVec(Eval(c, t));
        const Vec3 d1 = Deriv(c, t, 1);
        const Vec3 d2 = Deriv(c, t, 2);
        const Vec3 d3 = Deriv(c, t, 3);
        REQUIRE_THAT(Length(p), WithinAbs(1.0, Tolerance));
        REQUIRE_THAT(Dot(p, d1), WithinAbs(0.0, Tolerance));
        REQUIRE_THAT(Dot(p, d2) + Dot(d1, d1), WithinAbs(0.0, 1e-11));
        REQUIRE_THAT(Dot(p, d3) + 3.0 * Dot(d1, d2), WithinAbs(0.0, 1e-10));
        // A rational curve's derivatives do not stop at its degree.
        REQUIRE(Length(Deriv(c, t, 4)) > 0.0);
    }
}

TEST_CASE("derivatives agree with finite differences on random curves",
          "[curve][bspline][property]")
{
    rc::prop("first and second derivatives, polynomial and rational", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const Interval d = Domain(c);
        const double h = 1e-4 * (d.GetUpper() - d.GetLower());
        const double t = RandomParameter(c);
        // Away from knots, where the curve may be only C^0.
        for (const double u : c.Knots) {
            RC_PRE(std::abs(u - t) > 3.0 * h);
        }
        // Fourth-order central differences.
        const auto at = [&](const double k) { return Stereon::ToVec(Eval(c, t + k * h)); };
        const Vec3 d1 = (at(-2.0) - 8.0 * at(-1.0) + 8.0 * at(1.0) - at(2.0)) / (12.0 * h);
        const Vec3 d2 =
            (-1.0 * at(-2.0) + 16.0 * at(-1.0) - 30.0 * at(0.0) + 16.0 * at(1.0) - at(2.0)) /
            (12.0 * h * h);
        const double scale = 1.0 + Length(Deriv(c, t, 1)) + Length(Deriv(c, t, 2));
        RC_ASSERT(Length(Deriv(c, t, 1) - d1) <= 1e-6 * scale);
        RC_ASSERT(Length(Deriv(c, t, 2) - d2) <= 1e-3 * scale);
    });
}

TEST_CASE("the derivative with respect to a pole is its basis function",
          "[curve][bspline][property]")
{
    // ADR-0021: Eval is generic over the poles' scalar, so a dual pole gives
    // dC/dP_i = N_i w_i / sum N_j w_j in that coordinate.
    rc::prop("for a pole picked at random", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double t = RandomParameter(c);
        RC_PRE(t < Domain(c).GetUpper());
        const std::size_t picked = *rc::gen::inRange<std::size_t>(0, c.Poles.size());
        Stereon::BasicBSplineCurve<Dual, 3> lifted{c.Degree, c.Knots, {}, {}};
        double denominator = 0.0;
        for (std::size_t i = 0; i < c.Poles.size(); ++i) {
            lifted.Poles.push_back(
                {Dual{c.Poles[i].X, i == picked ? 1.0 : 0.0}, c.Poles[i].Y, c.Poles[i].Z});
            const double w = c.Weights.empty() ? 1.0 : c.Weights[i];
            denominator += w * CoxDeBoor(c.Knots, i, c.Degree, t);
        }
        for (const double w : c.Weights) {
            lifted.Weights.push_back(Dual{w});
        }
        const double w = c.Weights.empty() ? 1.0 : c.Weights[picked];
        const double expected = w * CoxDeBoor(c.Knots, picked, c.Degree, t) / denominator;
        RC_ASSERT(std::abs(Eval(lifted, t).X.Deriv - expected) <= 1e-12);
    });
}

TEST_CASE("knot insertion, splitting and Bezier extraction keep the curve",
          "[curve][bspline][property]")
{
    rc::prop("the point set and the parametrisation are unchanged", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double t = RandomParameter(c);
        const double knot = RandomParameter(c);
        const auto inserted = Stereon::InsertKnot(c, knot);
        if (inserted) {
            RC_ASSERT(Near(Eval(*inserted, t), Eval(c, t), 1e-11));
            RC_ASSERT(inserted->Poles.size() == c.Poles.size() + 1);
        }
        for (const BSplineCurve& segment : Stereon::ToBezierSegments(c)) {
            const Interval d = Domain(segment);
            RC_ASSERT(segment.Knots.size() == 2 * segment.Poles.size());
            RC_ASSERT(Near(Eval(segment, d.GetLower()), Eval(c, d.GetLower()), 1e-11));
            RC_ASSERT(Near(Eval(segment, d.GetUpper()), Eval(c, d.GetUpper()), 1e-11));
            RC_ASSERT(
                Near(Eval(segment, Stereon::Midpoint(d)), Eval(c, Stereon::Midpoint(d)), 1e-11));
        }
        const auto halves = Stereon::Split(c, knot);
        if (halves) {
            const BSplineCurve& piece = t <= knot ? halves->first : halves->second;
            RC_ASSERT(Near(Eval(piece, t), Eval(c, t), 1e-11));
            RC_ASSERT(Domain(halves->first).GetUpper() == knot);
            RC_ASSERT(Domain(halves->second).GetLower() == knot);
        }
    });
}

TEST_CASE("reversed, moved and reparametrised curves trace the same points",
          "[curve][bspline][property]")
{
    rc::prop("Reversed, Transformed by an affine map, Reparametrized", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const Interval d = Domain(c);
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(Reversed(c), d.GetLower() + d.GetUpper() - t), Eval(c, t), 1e-11));
        // B-splines are affinely invariant: a non-uniform scaling is fine.
        const Transform x = Stereon::MakeTranslation(Vec3{1.0, -2.0, 0.5}) *
                            Stereon::MakeScaling(Vec3{2.0, 0.5, -1.5});
        RC_ASSERT(Near(Eval(Transformed(c, x), t), Transformed(Eval(c, t), x), 1e-11));
        const BSplineCurve unit = *Stereon::Reparametrized(c, Interval{0.0, 1.0});
        RC_ASSERT((Domain(unit) == Interval{0.0, 1.0}));
        const double s = (t - d.GetLower()) / (d.GetUpper() - d.GetLower());
        RC_ASSERT(Near(Eval(unit, std::min(s, 1.0)), Eval(c, t), 1e-10));
    });
}

TEST_CASE("Project finds the global closest point", "[curve][bspline][property]")
{
    rc::prop("no sampled point is closer, and points of the curve project to themselves", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const Point3 p{*rc::gen::inRange(-60, 60) / 10.0, *rc::gen::inRange(-60, 60) / 10.0,
                       *rc::gen::inRange(-60, 60) / 10.0};
        const double t = Project(c, p);
        RC_ASSERT(Contains(Domain(c), t));
        const double distance = Distance(Eval(c, t), p);
        const Interval d = Domain(c);
        constexpr int samples = 4000;
        double sampled = Distance(Eval(c, d.GetLower()), p);
        for (int i = 1; i <= samples; ++i) {
            const double u = d.GetLower() + i * (d.GetUpper() - d.GetLower()) / samples;
            sampled = std::min(sampled, Distance(Eval(c, std::min(u, d.GetUpper())), p));
        }
        RC_ASSERT(distance <= sampled + 1e-12);

        const Point3 on = Eval(c, RandomParameter(c));
        RC_ASSERT(Distance(Eval(c, Project(c, on)), on) <= 1e-9);
    });
}

TEST_CASE("Project handles degenerate and simple cases", "[curve][bspline]")
{
    // From the centre every point of the quarter circle is at distance 1,
    // up to rounding; any parameter of the domain is a valid answer.
    const double centre = Project(QuarterCircle(), Point3{});
    REQUIRE(Contains(Interval{0.0, 1.0}, centre));
    REQUIRE_THAT(Length(Stereon::ToVec(Eval(QuarterCircle(), centre))), WithinAbs(1.0, Tolerance));
    // Interior feet of the perpendicular, in space and in the plane.
    const BSplineCurve line{
        1, {0.0, 0.0, 1.0, 1.0}, {Point3{-1.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}}, {}};
    REQUIRE_THAT(Project(line, Point3{0.25, 3.0, 0.0}), WithinAbs(0.625, Tolerance));
    const BSplineCurve2d planar{1, {0.0, 0.0, 1.0, 1.0}, {Point2{0.0, 0.0}, Point2{2.0, 0.0}}, {}};
    REQUIRE_THAT(Project(planar, Point2{1.5, 1.0}), WithinAbs(0.75, Tolerance));
}

TEST_CASE("a projection tie across segments goes to the smaller parameter",
          "[curve][bspline][projection][audit]")
{
    // Audit 2026-10-20, F9. The polyline (20, 0) -> (5, 0) -> (12, -1) ->
    // (3, 4) on the knots 0, 1, 2, 3 is at distance 5 from the origin at
    // t = 1, the vertex (5, 0), and at t = 3, the vertex (3, 4), and farther
    // everywhere else; every distance and box bound here is exact. The last
    // segment's pole box comes within 3 of the origin, so the search meets
    // it first and holds t = 3. The first two segments' boxes are at exactly
    // 5, the best distance so far, and must still be searched: a tie goes to
    // the smaller parameter, so pruning a segment whose bound only equals
    // the best distance would answer 3.
    const BSplineCurve polyline{1,
                                {0.0, 0.0, 1.0, 2.0, 3.0, 3.0},
                                {Point3{20.0, 0.0, 0.0}, Point3{5.0, 0.0, 0.0},
                                 Point3{12.0, -1.0, 0.0}, Point3{3.0, 4.0, 0.0}},
                                {}};
    const Point3 origin{};
    REQUIRE(DistanceSquared(Eval(polyline, 1.0), origin) == 25.0);
    REQUIRE(DistanceSquared(Eval(polyline, 3.0), origin) == 25.0);
    REQUIRE(Project(polyline, origin) == 1.0);
    REQUIRE(Stereon::CurveProjector(polyline).Project(origin) == 1.0);
    // Reversed, the tie is between t = 0, the start, and t = 2.
    REQUIRE(Project(Reversed(polyline), origin) == 0.0);
}

TEST_CASE("MakeBSplineCurve validates data from outside the kernel", "[curve][bspline][factory]")
{
    const std::vector<Point3> poles{Point3{}, Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}};
    const std::vector<double> clamped{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    REQUIRE(Stereon::MakeBSplineCurve(2, clamped, poles).has_value());
    REQUIRE(Stereon::MakeBSplineCurve(2, clamped, poles, {1.0, 2.0, 1.0}).has_value());
    const auto codeOf = [](const auto& r) { return r.error().GetCode(); };
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(0, clamped, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(26, clamped, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(3, clamped, poles)) == invalid); // too few poles
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, {0.0, 0.0, 1.0, 1.0, 1.0}, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, {0.0, 0.0, 0.0, 1.0, 0.5, 1.0}, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0}, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, {0.0, 0.0, 0.0, 0.0, 1.0},
                                             {Point3{}, Point3{}, Point3{}})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, {0.0, 0.0, 0.0, NAN, 1.0, 1.0}, poles)) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, clamped, poles, {1.0, 0.0, 1.0})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(2, clamped, poles, {1.0, 1.0})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(
                2, clamped, {Point3{}, Point3{INFINITY, 0.0, 0.0}, Point3{}})) == invalid);
    // An interior knot of multiplicity p + 1 would break the curve.
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, {0.0, 0.0, 0.5, 0.5, 1.0, 1.0},
                                             {Point3{}, Point3{}, Point3{}, Point3{}})) == invalid);
    REQUIRE(Stereon::MakeBSplineCurve(2, clamped, poles, {1.0, -1.0, 1.0}).error().GetMessage() ==
            "MakeBSplineCurve: weights must be finite and positive");
}

TEST_CASE("weights outside the supported range are refused, and a common scale changes nothing",
          "[curve][bspline][factory][audit]")
{
    // Audit 2026-10-20, F6: weights of 1e308 were admitted, and the weighted
    // pole sum overflowed at t = 0.5 although the curve is the segment from
    // (10, 0, 0) to (20, 0, 0).
    constexpr double MinWeight = Stereon::Detail::MinWeight;
    constexpr double MaxWeight = Stereon::Detail::MaxWeight;
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const std::vector<Point3> poles{Point3{10.0, 0.0, 0.0}, Point3{20.0, 0.0, 0.0}};
    const auto codeOf = [](const auto& r) { return r.error().GetCode(); };
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, knots, poles, {1e308, 1e308})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, knots, poles, {1e-320, 1.0})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, knots, poles, {2.0 * MaxWeight, 1.0})) == invalid);
    REQUIRE(codeOf(Stereon::MakeBSplineCurve(1, knots, poles, {1.0, 0.5 * MinWeight})) == invalid);
    // At the ends of the range the sums are exact.
    for (const double w : {MinWeight, MaxWeight}) {
        const auto curve = Stereon::MakeBSplineCurve(1, knots, poles, {w, w});
        REQUIRE(curve.has_value());
        REQUIRE(Eval(*curve, 0.5) == Point3{15.0, 0.0, 0.0});
    }

    // A rational cubic with every weight scaled by one power of two: the
    // same points and derivatives, bit for bit, and the same projection.
    const std::vector<double> cubicKnots{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0};
    const std::vector<Point3> cubicPoles{Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0},
                                         Point3{2.0, -1.0, 1.0}, Point3{3.0, 1.0, 2.0},
                                         Point3{4.0, 0.0, 0.0}};
    const std::vector<double> weights{1.0, 0.5, 2.0, 0.75, 1.0};
    const BSplineCurve reference = *Stereon::MakeBSplineCurve(3, cubicKnots, cubicPoles, weights);
    const Point3 p{1.0, 2.0, -0.5};
    for (const int e : {-200, -60, 1, 60, 200}) {
        INFO("weights scaled by 2^" << e);
        std::vector<double> scaled;
        for (const double w : weights) {
            scaled.push_back(std::ldexp(w, e));
        }
        const BSplineCurve curve = *Stereon::MakeBSplineCurve(3, cubicKnots, cubicPoles, scaled);
        for (const double t : {0.0, 0.1, 0.37, 0.5, 0.81, 1.0}) {
            REQUIRE(Eval(curve, t) == Eval(reference, t));
            REQUIRE(Deriv(curve, t, 1) == Deriv(reference, t, 1));
            REQUIRE(Deriv(curve, t, 2) == Deriv(reference, t, 2));
        }
        REQUIRE_THAT(Project(curve, p), WithinAbs(Project(reference, p), Tolerance));
        const BSplineCurve refined = *Stereon::InsertKnot(curve, 0.3, 2);
        REQUIRE(Near(Eval(refined, 0.7), Eval(reference, 0.7)));
    }
}

TEST_CASE("a knot vector whose width overflows is refused", "[curve][bspline][factory][audit]")
{
    // Audit 2026-10-20, F6: knots of -1e308 and 1e308 are each finite, but
    // the basis functions divide by their difference, which is not.
    const std::vector<Point3> poles{Point3{}, Point3{1.0, 0.0, 0.0}};
    REQUIRE(Stereon::MakeBSplineCurve(1, {-1e308, -1e308, 1e308, 1e308}, poles).error().GetCode() ==
            ErrorCode::InvalidInput);
    const auto wide = Stereon::MakeBSplineCurve(1, {-1e307, -1e307, 1e307, 1e307}, poles);
    REQUIRE(wide.has_value());
    REQUIRE(Near(Eval(*wide, 0.0), Point3{0.5, 0.0, 0.0}));
}

TEST_CASE("the knot operations validate their parameters", "[curve][bspline][factory]")
{
    const BSplineCurve c = MonomialCurve();
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::InsertKnot(c, 2.0).error().GetCode() == invalid);
    REQUIRE(Stereon::InsertKnot(c, 0.5, 4).error().GetCode() == invalid); // beyond the degree
    REQUIRE(Stereon::InsertKnot(c, 0.0).error().GetCode() == invalid);    // already p + 1
    REQUIRE(Stereon::InsertKnot(c, 0.5, 3).has_value());
    REQUIRE(Stereon::Split(c, 0.0).error().GetCode() == invalid);
    REQUIRE(Stereon::Split(c, NAN).error().GetCode() == invalid);
    REQUIRE(Stereon::Reparametrized(c, Interval{1.0, 1.0}).error().GetCode() == invalid);
    REQUIRE(Stereon::ToBezierSegments(c).size() == 1);
    REQUIRE(Stereon::ToBezierSegments(*Stereon::InsertKnot(c, 0.5)).size() == 2);

    const auto expanded = Stereon::ExpandKnots(std::vector{0.0, 1.0, 2.0}, std::vector{3, 1, 3});
    REQUIRE(*expanded == std::vector{0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0});
    REQUIRE_FALSE(Stereon::ExpandKnots(std::vector{0.0, 0.0}, std::vector{1, 1}).has_value());
    REQUIRE_FALSE(Stereon::ExpandKnots(std::vector{0.0, 1.0}, std::vector{1, 0}).has_value());
    REQUIRE_FALSE(Stereon::ExpandKnots(std::vector{0.0}, std::vector{1, 2}).has_value());
}

TEST_CASE("Derivatives returns every order at once", "[curve][bspline][property]")
{
    rc::prop("the same values as Eval and Deriv, order by order", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double t = RandomParameter(c);
        const int n = *rc::gen::inRange(0, c.Degree + 3);
        std::vector<Vec3> all(static_cast<std::size_t>(n) + 1);
        Stereon::Derivatives(c, t, n, std::span{all});
        RC_ASSERT(Near(Stereon::ToPoint(all[0]), Eval(c, t), 1e-13));
        for (int k = 1; k <= n; ++k) {
            RC_ASSERT(all[static_cast<std::size_t>(k)] == Deriv(c, t, k));
        }
    });
}

TEST_CASE("a CurveProjector answers as Project does", "[curve][bspline][property]")
{
    rc::prop("bit for bit, for many points on one prepared curve", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const Stereon::CurveProjector projector(c);
        for (int i = 0; i < 8; ++i) {
            const Point3 p{*rc::gen::inRange(-60, 60) / 10.0, *rc::gen::inRange(-60, 60) / 10.0,
                           *rc::gen::inRange(-60, 60) / 10.0};
            RC_ASSERT(projector.Project(p) == Project(c, p));
        }
    });
    // A planar curve is projected in its plane.
    const BSplineCurve2d planar{1, {0.0, 0.0, 1.0, 1.0}, {Point2{0.0, 0.0}, Point2{2.0, 0.0}}, {}};
    REQUIRE(Stereon::CurveProjector(planar).Project(Point2{1.5, 1.0}) ==
            Project(planar, Point2{1.5, 1.0}));
}
