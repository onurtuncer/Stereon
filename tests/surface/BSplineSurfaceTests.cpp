// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// B-spline and NURBS surfaces, against independent references written here:
// the tensor product of the Cox-de Boor recursion for points, an exact
// polynomial surface and an extruded circle for derivatives, finite
// differences for random surfaces, and dense sampling for projection.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

import Stereon.Surface;

#include "Dual.hpp"
#include "Formatting.hpp"

using Catch::Matchers::WithinAbs;
using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::ParameterDirection;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
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

/// N_{i,p}(t) by the Cox-de Boor recursion, with 0/0 = 0.
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

/// The point as the tensor-product sum over every pole.
Point3 ReferencePoint(const BSplineSurface& s, const Point2& uv)
{
    const std::size_t countV = Stereon::PoleCountV(s);
    Vec3 numerator{};
    double denominator = 0.0;
    for (std::size_t k = 0; k < s.Poles.size(); ++k) {
        const double w = s.Weights.empty() ? 1.0 : s.Weights[k];
        const double n = CoxDeBoor(s.KnotsU, k / countV, s.DegreeU, uv.X) *
                         CoxDeBoor(s.KnotsV, k % countV, s.DegreeV, uv.Y);
        numerator += (n * w) * Stereon::ToVec(s.Poles[k]);
        denominator += n * w;
    }
    return Stereon::ToPoint(numerator / denominator);
}

/// A random knot vector for degree p and the given pole count, clamped or
/// not, with interior multiplicities up to p.
std::vector<double> RandomKnots(const int p, const std::size_t poles)
{
    const auto q = static_cast<std::size_t>(p);
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (*rc::gen::arbitrary<bool>()) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    return knots;
}

BSplineSurface RandomSurface(const bool rational)
{
    const int p = *rc::gen::inRange(1, 5);
    const int q = *rc::gen::inRange(1, 5);
    const auto countU = static_cast<std::size_t>(p + 1 + *rc::gen::inRange(0, 4));
    const auto countV = static_cast<std::size_t>(q + 1 + *rc::gen::inRange(0, 4));
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t k = 0; k < countU * countV; ++k) {
        poles.push_back({static_cast<double>(k / countV) + *rc::gen::inRange(-5, 5) / 10.0,
                         static_cast<double>(k % countV) + *rc::gen::inRange(-5, 5) / 10.0,
                         *rc::gen::inRange(-30, 30) / 10.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(3, 30) / 10.0);
        }
    }
    const auto made = Stereon::MakeBSplineSurface(p, q, RandomKnots(p, countU),
                                                  RandomKnots(q, countV), poles, weights);
    RC_PRE(made.has_value());
    return *made;
}

/// uv kept inside the box, against rounding in the grid's last step.
Point2 Clamped(const Point2& uv, const Box2& box)
{
    return {std::min(uv.X, box.X.GetUpper()), std::min(uv.Y, box.Y.GetUpper())};
}

Point2 RandomParameter(const BSplineSurface& s)
{
    const Box2 d = Domain(s);
    const double a = *rc::gen::inRange(0, 1000) / 1000.0;
    const double b = *rc::gen::inRange(0, 1000) / 1000.0;
    return {d.X.GetLower() + a * (d.X.GetUpper() - d.X.GetLower()),
            d.Y.GetLower() + b * (d.Y.GetUpper() - d.Y.GetLower())};
}

/// (u, v, u^2 v^3) on [0, 1]^2 as a Bezier surface of degree (2, 3): the
/// monomial u^a v^b has Bernstein coefficients C(i, a)/C(2, a) C(j, b)/C(3, b).
BSplineSurface MonomialSurface()
{
    std::vector<Point3> poles;
    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j <= 3; ++j) {
            poles.push_back({i / 2.0, j / 3.0, (i == 2 && j == 3) ? 1.0 : 0.0});
        }
    }
    return {2,     3, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
            poles, {}};
}

/// A quarter of the unit cylinder: the rational quarter circle in u
/// extruded linearly along z in v.
BSplineSurface QuarterCylinder()
{
    const double w = std::numbers::sqrt2 / 2.0;
    return {2,
            1,
            {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
            {0.0, 0.0, 1.0, 1.0},
            {Point3{1.0, 0.0, 0.0}, Point3{1.0, 0.0, 2.0}, Point3{1.0, 1.0, 0.0},
             Point3{1.0, 1.0, 2.0}, Point3{0.0, 1.0, 0.0}, Point3{0.0, 1.0, 2.0}},
            {1.0, 1.0, w, w, 1.0, 1.0}};
}

} // namespace

TEST_CASE("B-spline surfaces satisfy the surface concept and are shared when wrapped",
          "[surface][bspline]")
{
    STATIC_REQUIRE(Stereon::SurfaceGeometry<BSplineSurface>);
    STATIC_REQUIRE_FALSE(Stereon::Detail::FitsInline<BSplineSurface, Stereon::GeometryBufferSize>);
    const Surface surface(MonomialSurface());
    const Surface copy = surface; // NOLINT(performance-unnecessary-copy-initialization)
    REQUIRE(copy.Target<BSplineSurface>() == surface.Target<BSplineSurface>());
    REQUIRE(Stereon::PoleCountU(MonomialSurface()) == 3);
    REQUIRE(Stereon::PoleCountV(MonomialSurface()) == 4);
    REQUIRE(Domain(MonomialSurface()) == Box2{Interval{0.0, 1.0}, Interval{0.0, 1.0}});
}

TEST_CASE("Eval agrees with the tensor-product basis sum", "[surface][bspline][property]")
{
    rc::prop("polynomial and rational random surfaces", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const Box2 d = Domain(s);
        RC_PRE(uv.X < d.X.GetUpper() && uv.Y < d.Y.GetUpper());
        RC_ASSERT(Near(Eval(s, uv), ReferencePoint(s, uv), 1e-12));
        RC_ASSERT(Contains(BoundingBox(s), Eval(s, uv)));
    });
}

TEST_CASE("derivatives of a polynomial surface are exact", "[surface][bspline]")
{
    const BSplineSurface one = MonomialSurface();
    const auto refinedU = Stereon::InsertKnot(one, ParameterDirection::U, 0.4, 2);
    const BSplineSurface many = *Stereon::InsertKnot(*refinedU, ParameterDirection::V, 0.7);
    for (const BSplineSurface& s : {one, many}) {
        for (int i = 0; i <= 8; ++i) {
            for (int j = 0; j <= 8; ++j) {
                const double u = i / 8.0;
                const double v = j / 8.0;
                const Point2 uv{u, v};
                REQUIRE(Near(Eval(s, uv), Point3{u, v, u * u * v * v * v}));
                REQUIRE(Near(Deriv(s, uv, 1, 0), Vec3{1.0, 0.0, 2.0 * u * v * v * v}));
                REQUIRE(Near(Deriv(s, uv, 0, 1), Vec3{0.0, 1.0, 3.0 * u * u * v * v}));
                REQUIRE(Near(Deriv(s, uv, 1, 1), Vec3{0.0, 0.0, 6.0 * u * v * v}));
                REQUIRE(Near(Deriv(s, uv, 2, 1), Vec3{0.0, 0.0, 6.0 * v * v}));
                REQUIRE(Near(Deriv(s, uv, 1, 2), Vec3{0.0, 0.0, 12.0 * u * v}));
                REQUIRE(Near(Deriv(s, uv, 2, 3), Vec3{0.0, 0.0, 12.0}));
                REQUIRE(Deriv(s, uv, 3, 0) == Vec3{});
                REQUIRE(Deriv(s, uv, 0, 4) == Vec3{});
            }
        }
    }
}

TEST_CASE("rational surface derivatives satisfy the cylinder's identities",
          "[surface][bspline][nurbs]")
{
    // x^2 + y^2 = 1 in u (the circle), linear in z along v.
    const BSplineSurface s = QuarterCylinder();
    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 4; ++j) {
            const Point2 uv{i / 10.0, j / 4.0};
            const Vec3 p = Stereon::ToVec(Eval(s, uv));
            const Vec3 radial{p.X, p.Y, 0.0};
            const Vec3 su = Deriv(s, uv, 1, 0);
            const Vec3 suu = Deriv(s, uv, 2, 0);
            REQUIRE_THAT(Length(radial), WithinAbs(1.0, Tolerance));
            REQUIRE_THAT(p.Z, WithinAbs(2.0 * uv.Y, Tolerance));
            REQUIRE_THAT(Dot(radial, su), WithinAbs(0.0, Tolerance));
            REQUIRE_THAT(Dot(radial, suu) + Dot(su, su), WithinAbs(0.0, 1e-11));
            REQUIRE(Near(Deriv(s, uv, 0, 1), Vec3{0.0, 0.0, 2.0}));
            REQUIRE(Near(Deriv(s, uv, 1, 1), Vec3{}));
            REQUIRE(Near(Deriv(s, uv, 3, 1), Vec3{}, 1e-10));
        }
    }
}

TEST_CASE("surface derivatives agree with finite differences", "[surface][bspline][property]")
{
    rc::prop("first, second and mixed partials of random surfaces", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Box2 d = Domain(s);
        const double hu = 1e-4 * (d.X.GetUpper() - d.X.GetLower());
        const double hv = 1e-4 * (d.Y.GetUpper() - d.Y.GetLower());
        const Point2 uv = RandomParameter(s);
        for (const double k : s.KnotsU) {
            RC_PRE(std::abs(k - uv.X) > 3.0 * hu);
        }
        for (const double k : s.KnotsV) {
            RC_PRE(std::abs(k - uv.Y) > 3.0 * hv);
        }
        const auto at = [&](const double a, const double b) {
            return Stereon::ToVec(Eval(s, Point2{uv.X + a * hu, uv.Y + b * hv}));
        };
        const Vec3 su = (at(-2, 0) - 8.0 * at(-1, 0) + 8.0 * at(1, 0) - at(2, 0)) / (12.0 * hu);
        const Vec3 sv = (at(0, -2) - 8.0 * at(0, -1) + 8.0 * at(0, 1) - at(0, 2)) / (12.0 * hv);
        const Vec3 suv = (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4.0 * hu * hv);
        const double scale = 1.0 + Length(Deriv(s, uv, 1, 0)) + Length(Deriv(s, uv, 0, 1)) +
                             Length(Deriv(s, uv, 1, 1));
        RC_ASSERT(Length(Deriv(s, uv, 1, 0) - su) <= 1e-6 * scale);
        RC_ASSERT(Length(Deriv(s, uv, 0, 1) - sv) <= 1e-6 * scale);
        RC_ASSERT(Length(Deriv(s, uv, 1, 1) - suv) <= 1e-4 * scale);
    });
}

TEST_CASE("the derivative with respect to a pole is its basis product",
          "[surface][bspline][property]")
{
    rc::prop("for a pole picked at random", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const Box2 d = Domain(s);
        RC_PRE(uv.X < d.X.GetUpper() && uv.Y < d.Y.GetUpper());
        const std::size_t picked = *rc::gen::inRange<std::size_t>(0, s.Poles.size());
        const std::size_t countV = Stereon::PoleCountV(s);
        Stereon::BasicBSplineSurface<Dual> lifted{s.DegreeU, s.DegreeV, s.KnotsU, s.KnotsV, {}, {}};
        double denominator = 0.0;
        for (std::size_t k = 0; k < s.Poles.size(); ++k) {
            lifted.Poles.push_back(
                {s.Poles[k].X, Dual{s.Poles[k].Y, k == picked ? 1.0 : 0.0}, s.Poles[k].Z});
            const double w = s.Weights.empty() ? 1.0 : s.Weights[k];
            denominator += w * CoxDeBoor(s.KnotsU, k / countV, s.DegreeU, uv.X) *
                           CoxDeBoor(s.KnotsV, k % countV, s.DegreeV, uv.Y);
        }
        for (const double w : s.Weights) {
            lifted.Weights.push_back(Dual{w});
        }
        const double w = s.Weights.empty() ? 1.0 : s.Weights[picked];
        const double expected = w * CoxDeBoor(s.KnotsU, picked / countV, s.DegreeU, uv.X) *
                                CoxDeBoor(s.KnotsV, picked % countV, s.DegreeV, uv.Y) / denominator;
        RC_ASSERT(std::abs(Eval(lifted, uv).Y.Deriv - expected) <= 1e-12);
    });
}

TEST_CASE("knot insertion, Bezier extraction and iso-curves keep the surface",
          "[surface][bspline][property]")
{
    rc::prop("the same points through every derived representation", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const Point2 knot = RandomParameter(s);
        for (const ParameterDirection direction : {ParameterDirection::U, ParameterDirection::V}) {
            const double value = direction == ParameterDirection::U ? knot.X : knot.Y;
            const auto inserted = Stereon::InsertKnot(s, direction, value);
            if (inserted) {
                RC_ASSERT(Near(Eval(*inserted, uv), Eval(s, uv), 1e-11));
            }
        }
        const auto alongV = Stereon::IsoCurve(s, ParameterDirection::U, uv.X);
        const auto alongU = Stereon::IsoCurve(s, ParameterDirection::V, uv.Y);
        RC_ASSERT(Near(Eval(*alongV, uv.Y), Eval(s, uv), 1e-11));
        RC_ASSERT(Near(Eval(*alongU, uv.X), Eval(s, uv), 1e-11));
        for (const BSplineSurface& patch : Stereon::ToBezierPatches(s)) {
            const Box2 r = Domain(patch);
            const Point2 centre{Stereon::Midpoint(r.X), Stereon::Midpoint(r.Y)};
            RC_ASSERT(Near(Eval(patch, centre), Eval(s, centre), 1e-11));
        }
    });
}

TEST_CASE("reversed and moved surfaces trace the same points", "[surface][bspline][property]")
{
    rc::prop("Reversed runs u backwards, Transformed takes any affine map", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const Interval du = Domain(s).X;
        const Point2 mirrored{du.GetLower() + du.GetUpper() - uv.X, uv.Y};
        RC_ASSERT(Near(Eval(Reversed(s), mirrored), Eval(s, uv), 1e-11));
        const Transform x = Stereon::MakeTranslation(Vec3{1.0, -2.0, 0.5}) *
                            Stereon::MakeScaling(Vec3{2.0, 0.5, -1.5});
        RC_ASSERT(Near(Eval(Transformed(s, x), uv), Transformed(Eval(s, uv), x), 1e-11));
    });
}

TEST_CASE("Project finds the closest point of a B-spline surface", "[surface][bspline][property]")
{
    rc::prop("no point of a dense grid is closer, and surface points project to themselves", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point3 p{*rc::gen::inRange(-20, 60) / 10.0, *rc::gen::inRange(-20, 60) / 10.0,
                       *rc::gen::inRange(-40, 40) / 10.0};
        const Point2 foot = Project(s, p);
        const Box2 d = Domain(s);
        RC_ASSERT(Contains(d, foot));
        const double distance = Distance(Eval(s, foot), p);
        constexpr int samples = 80;
        double sampled = distance + 1.0;
        for (int i = 0; i <= samples; ++i) {
            for (int j = 0; j <= samples; ++j) {
                const double u = d.X.GetLower() + i * (d.X.GetUpper() - d.X.GetLower()) / samples;
                const double v = d.Y.GetLower() + j * (d.Y.GetUpper() - d.Y.GetLower()) / samples;
                sampled = std::min(sampled, Distance(Eval(s, Clamped(Point2{u, v}, d)), p));
            }
        }
        RC_ASSERT(distance <= sampled + 1e-9);
        const Point3 on = Eval(s, RandomParameter(s));
        RC_ASSERT(Distance(Eval(s, Project(s, on)), on) <= 1e-8);
    });
}

TEST_CASE("projection onto a NURBS cylinder matches the analytic one", "[surface][bspline]")
{
    const BSplineSurface s = QuarterCylinder();
    for (int i = 0; i <= 6; ++i) {
        const double angle = i * std::numbers::pi / 12.0;
        const Point3 p{3.0 * std::cos(angle), 3.0 * std::sin(angle), 0.75};
        const Point3 foot = Eval(s, Project(s, p));
        REQUIRE(Near(foot, Point3{std::cos(angle), std::sin(angle), 0.75}, 1e-10));
    }
}

TEST_CASE("MakeBSplineSurface validates data from outside the kernel",
          "[surface][bspline][factory]")
{
    const BSplineSurface good = MonomialSurface();
    const auto make = [&](const int p, const int q, std::vector<double> u, std::vector<double> v,
                          std::vector<Point3> poles, std::vector<double> w = {}) {
        return Stereon::MakeBSplineSurface(p, q, std::move(u), std::move(v), std::move(poles),
                                           std::move(w));
    };
    REQUIRE(make(2, 3, good.KnotsU, good.KnotsV, good.Poles).has_value());
    REQUIRE(make(2, 3, good.KnotsU, good.KnotsV, good.Poles,
                 std::vector<double>(good.Poles.size(), 2.0))
                .has_value());
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(make(0, 3, good.KnotsU, good.KnotsV, good.Poles).error().GetCode() == invalid);
    REQUIRE(make(2, 3, good.KnotsU, {0.0, 0.0, 1.0}, good.Poles).error().GetCode() == invalid);
    std::vector<Point3> fewer = good.Poles;
    fewer.pop_back();
    REQUIRE(make(2, 3, good.KnotsU, good.KnotsV, fewer).error().GetCode() == invalid);
    std::vector<Point3> nan = good.Poles;
    nan[3].Z = NAN;
    REQUIRE(make(2, 3, good.KnotsU, good.KnotsV, nan).error().GetCode() == invalid);
    REQUIRE(make(2, 3, good.KnotsU, good.KnotsV, good.Poles, {1.0}).error().GetCode() == invalid);
    REQUIRE(Stereon::InsertKnot(good, ParameterDirection::U, 2.0).error().GetCode() == invalid);
    REQUIRE(Stereon::InsertKnot(good, ParameterDirection::V, 0.5, 4).error().GetCode() == invalid);
    REQUIRE(Stereon::IsoCurve(good, ParameterDirection::V, -0.5).error().GetCode() == invalid);
    REQUIRE(Stereon::ToBezierPatches(good).size() == 1);
}
