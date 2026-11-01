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
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
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
    STATIC_REQUIRE_FALSE(Stereon::Detail::FitsInline<BSplineSurface, Stereon::SurfaceBufferSize>);
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

TEST_CASE("a point on a rational bilinear surface projects to itself",
          "[surface][bspline][regression]")
{
    // The shrunk counterexample of a property run (RapidCheck seed
    // 1556740052594640537, CI 2026-10-16): an unclamped bilinear patch whose
    // pole weights range from 0.3 to 2.0, and a surface point near its low-u
    // edge, which projected 3.2e-3 away.
    const BSplineSurface s = *Stereon::MakeBSplineSurface(
        1, 1, std::vector<double>{0.25, 0.5, 1.25, 1.5, 2.0, 2.25},
        std::vector<double>{0.25, 0.5, 0.75, 1.0},
        std::vector<Point3>{Point3{-0.5, -0.2, -3.0}, Point3{0.0, 0.5, -0.4},
                            Point3{0.5, -0.2, 0.0}, Point3{0.5, 0.7, 0.2}, Point3{1.5, -0.5, -3.0},
                            Point3{1.5, 0.5, -3.0}, Point3{2.5, -0.5, -3.0},
                            Point3{2.5, 0.5, -3.0}},
        std::vector<double>{1.6, 0.3, 0.3, 2.0, 0.3, 0.3, 0.3, 0.3});
    const Point3 on = Eval(s, Point2{0.557, 0.652});
    const Point2 foot = Project(s, on);
    CHECK(Distance(Eval(s, foot), on) <= 1e-8);
}

TEST_CASE("a minimum just inside a patch is not lost to its boundary",
          "[surface][bspline][regression]")
{
    // Case 9441 of the B1.3 convergence corpus b13-v1 (bench/convergence/
    // ProjectionCorpus.cpp, 2026-11-01): a rational bilinear surface and a
    // point 1.3e-4 from it, whose closest point lies 2.4e-4 inside the patch
    // below the knot line v = 0.97234686. The leaf's Newton run walked onto
    // that knot line, the patch's boundary, and was accepted there: the
    // projector answered with the boundary point, 1.37e-3 away.
    const std::vector<double> knotsU{
        0, 0, 0.055867720548746752, 0.40721838323166171, 0.8942613552601637, 1, 1};
    const std::vector<double> knotsV{
        0, 0, 0.09695669622675776, 0.10291401872957094, 0.14536816132113231, 0.97234686380282498,
        1, 1};
    const std::vector<Point3> poles{
        {-0.53702742898688527, -0.28196200361572538, 1.4578296119761234},
        {0.51427311924283392, 1.6241224317803389, -1.4339291592438115},
        {0.33645814525446405, 3.7549856736452325, 0.65766850113939634},
        {0.23100015374094029, 6.5941033045926458, -1.3246503937034673},
        {0.29103864446232985, 8.3181466359173353, -2.2226952993646143},
        {0.29075656318378063, 10.228892205893205, 0.27036552304624495},
        {1.8238210872794256, 0.036108953362221197, 1.2309808918650553},
        {2.4256859273257354, 2.3526642517137666, -0.81860835775450091},
        {1.9284018043686761, 3.931132618312009, -1.8823107802638457},
        {2.42780810815061, 6.1220080641908012, 0.91268063048847869},
        {1.9146232067427311, 7.7154911804683328, 2.3216753914996948},
        {2.5176315642630285, 9.6102436105804223, -0.35665577560027728},
        {3.445197473375488, -0.23485075244806075, 2.4790606734135796},
        {4.277172029177243, 2.1444894962209067, 2.9219298458250114},
        {3.8255655295594297, 3.5506310568521111, 0.49852520505326225},
        {3.9018347357374354, 6.5995999368982474, -1.1057038976126328},
        {3.7423449709216947, 8.5316089093180949, -0.12718303312314783},
        {4.4141376297947854, 10.400675085787363, -1.7315409082397686},
        {6.1186747503747894, -0.24040651578845207, 0.95280842798900522},
        {5.9199632688478117, 2.0853155571701296, 2.9856360576692893},
        {5.5477193816411532, 4.5761435496722882, 1.5354103159042047},
        {5.6175519703599264, 6.3067346241146662, -0.35123118778195206},
        {6.1834982851255074, 8.5581840438070298, 0.81083112366696497},
        {6.4291193775503279, 9.9946510573221481, 1.933952439397574},
        {8.3287266625238789, -0.35280039901697891, -0.69167529891162216},
        {7.9941753669034243, 1.7898631158588272, -2.509729075446562},
        {8.2579608693368414, 4.3696552505336941, -1.0397230961405257},
        {8.5809326826917491, 5.9705678381512044, 0.7288733316828151},
        {8.2321343842881518, 7.8980009501071002, -0.50173492557591093},
        {8.1186961115692107, 10.441360694801618, 2.4287467385210819}};
    const std::vector<double> weights{
        1.0575550228578328,  0.87575871249378134, 1.5177156145397412,  1.6869029486608607,
        0.45863900602563357, 2.551447380984746,   0.51163045506365745, 0.79963687374181414,
        1.8605254454432958,  1.4809800028032518,  0.96077131124524606, 0.73366510974815768,
        0.70119030927996662, 2.3888259159672613,  2.5483203690660363,  2.9988955791424239,
        0.94882044894051787, 0.77636448339163566, 0.7449750820182639,  2.9081370227010135,
        0.84367353530685829, 1.4611431354395943,  1.0333191477018169,  1.8066572423204028,
        2.387852984475114,   1.8828354836238828,  1.391735001863319,   0.58916933385262715,
        1.2085213198994931,  0.54522851614087053};
    const BSplineSurface s = *Stereon::MakeBSplineSurface(1, 1, knotsU, knotsV, poles, weights);
    const Point3 point{2.7338072080782165, 8.0800729599302805, 1.2233894280799094};
    // The closest point the corpus's independent reference found, and its distance.
    const Point2 reference{0.21437379822600636, 0.97210604249558241};
    const double closest = Distance(Eval(s, reference), point);
    REQUIRE(closest < 1.32e-4);
    const Point2 foot = Project(s, point);
    INFO("Project gave (" << foot.X << ", " << foot.Y << ")");
    CHECK(Distance(Eval(s, foot), point) <= closest * (1.0 + 1e-9));
    CHECK(Distance(Eval(s, Stereon::SurfaceProjector(s).Project(point)), point) <=
          closest * (1.0 + 1e-9));
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

TEST_CASE("surface weights outside the supported range are refused",
          "[surface][bspline][factory][audit]")
{
    // Audit 2026-10-20, F6: the surface factory shares the curve's weight
    // range, so no rational patch overflows its weighted sums either.
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const std::vector<Point3> poles{Point3{0.0, 0.0, 0.0}, Point3{0.0, 2.0, 0.0},
                                    Point3{2.0, 0.0, 0.0}, Point3{2.0, 2.0, 4.0}};
    const auto make = [&](const double w) {
        return Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, std::vector<double>(4, w));
    };
    REQUIRE(make(1e300).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(make(1e-300).error().GetCode() == ErrorCode::InvalidInput);
    for (const double w : {Stereon::Detail::MinWeight, 1.0, Stereon::Detail::MaxWeight}) {
        const auto patch = make(w);
        REQUIRE(patch.has_value());
        REQUIRE(Eval(*patch, Point2{0.5, 0.5}) == Point3{1.0, 1.0, 1.0});
    }
}

TEST_CASE("Derivatives returns every partial at once", "[surface][bspline][property]")
{
    rc::prop("the same values as Eval and Deriv, partial by partial", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const int du = *rc::gen::inRange(0, 4);
        const int dv = *rc::gen::inRange(0, 4);
        std::vector<Vec3> all(static_cast<std::size_t>((du + 1) * (dv + 1)));
        Stereon::Derivatives(s, uv, du, dv, std::span{all});
        RC_ASSERT(Near(Stereon::ToPoint(all[0]), Eval(s, uv), 1e-13));
        for (int k = 0; k <= du; ++k) {
            for (int l = 0; l <= dv; ++l) {
                const Vec3& partial = all[static_cast<std::size_t>(k * (dv + 1) + l)];
                if (k + l == 1) {
                    // Deriv of first order runs the first-derivative kernel
                    // (Detail::FirstDerivatives), the general path here
                    // differences the poles: the same values up to rounding.
                    const double scale = 1.0 + Length(partial) + Length(all[0]);
                    RC_ASSERT(Length(partial - Deriv(s, uv, k, l)) <= 1e-12 * scale);
                } else if (k + l > 1) {
                    RC_ASSERT(partial == Deriv(s, uv, k, l));
                }
            }
        }
    });
}

TEST_CASE("TotalDerivatives returns the partials of bounded total order",
          "[surface][bspline][property]")
{
    // First order has a kernel of its own (Detail::FirstDerivatives), which
    // sums basis derivatives where the general path differences poles: the
    // same values up to rounding. Every other order runs the general path
    // and agrees to the bit.
    rc::prop("those of Derivatives with k + l <= n", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        const int n = *rc::gen::inRange(0, 5);
        const auto size = static_cast<std::size_t>((n + 1) * (n + 1));
        std::vector<Vec3> square(size);
        std::vector<Vec3> total(size);
        Stereon::Derivatives(s, uv, n, n, std::span{square});
        Stereon::TotalDerivatives(s, uv, n, std::span{total});
        for (int k = 0; k <= n; ++k) {
            for (int l = 0; k + l <= n; ++l) {
                const auto i = static_cast<std::size_t>(k * (n + 1) + l);
                if (n == 1) {
                    const double scale = 1.0 + Length(square[i]) + Length(square[0]);
                    RC_ASSERT(Length(total[i] - square[i]) <= 1e-12 * scale);
                } else {
                    RC_ASSERT(total[i] == square[i]);
                }
            }
        }
    });
}

TEST_CASE("the first-derivative kernel agrees with the general path and with Deriv",
          "[surface][bspline][property]")
{
    rc::prop("S, S_u and S_v from one pass over the poles", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Point2 uv = RandomParameter(s);
        std::array<Vec3, 4> first{};
        std::array<Vec3, 4> general{};
        Stereon::TotalDerivatives(s, uv, 1, std::span<Vec3>{first});
        Stereon::Derivatives(s, uv, 1, 1, std::span<Vec3>{general});
        const double scale = 1.0 + Length(general[0]) + Length(general[1]) + Length(general[2]);
        RC_ASSERT(Length(first[0] - general[0]) <= 1e-12 * scale);
        RC_ASSERT(Length(first[1] - general[1]) <= 1e-12 * scale);
        RC_ASSERT(Length(first[2] - general[2]) <= 1e-12 * scale);
        // Deriv of first order takes the same kernel.
        RC_ASSERT(Deriv(s, uv, 1, 0) == first[2]);
        RC_ASSERT(Deriv(s, uv, 0, 1) == first[1]);
        RC_ASSERT(Length(ToVec(Eval(s, uv)) - first[0]) <= 1e-12 * scale);
    });
}

TEST_CASE("a SurfaceProjector answers as Project does", "[surface][bspline][property]")
{
    rc::prop("bit for bit, for several points on one prepared surface", [] {
        const BSplineSurface s = RandomSurface(*rc::gen::arbitrary<bool>());
        const Stereon::SurfaceProjector projector(s);
        for (int i = 0; i < 4; ++i) {
            const Point3 p{*rc::gen::inRange(-20, 60) / 10.0, *rc::gen::inRange(-20, 60) / 10.0,
                           *rc::gen::inRange(-40, 40) / 10.0};
            RC_ASSERT(projector.Project(p) == Project(s, p));
        }
    });
}

TEST_CASE("surface factories refuse reciprocal-overflow spans in either direction",
          "[surface][bspline][factory][regression]")
{
    const std::vector<double> regular{0.0, 0.0, 1.0, 1.0};
    const std::vector<double> narrow{0.0, 0.0, 1e-309, 1e-309};
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
    for (const bool narrowU : {false, true}) {
        const auto surface =
            Stereon::MakeBSplineSurface(1, 1, narrowU ? narrow : regular,
                                        narrowU ? regular : narrow, poles, std::vector<double>{});
        REQUIRE_FALSE(surface.has_value());
        CHECK(surface.error().GetCode() == ErrorCode::InvalidInput);
    }
}
