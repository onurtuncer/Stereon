// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Phase 2: a successful empty result requires a mathematical separation
// proof, never simply failure to find a seed or a small iteration budget;
// branches are held to closed forms, their certificates and a sampling
// oracle; contact that is not transverse is a point or fails closed.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Context;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();

/// A bilinear NURBS patch over a unit parameter square, with the optional
/// diagonal pole deliberately moved to stretch the convex hull.
BSplineSurface Patch(const double z, const bool rational = false)
{
    std::vector<Point3> poles;
    poles.push_back({0.0, 0.0, z});
    poles.push_back({0.0, 1.0, z});
    poles.push_back({1.0, 0.0, z});
    poles.push_back({1.0, 1.0, z});
    const std::vector<double> weights =
        rational ? std::vector<double>{1.0, 3.0, 2.0, 4.0} : std::vector<double>{};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

const Box2 middle{Interval{0.25, 0.75}, Interval{0.25, 0.75}};

} // namespace

TEST_CASE("separated polynomial and rational NURBS control nets prove empty SSI",
          "[intersect][generic][broadphase]")
{
    for (const bool rationalA : {false, true}) {
        for (const bool rationalB : {false, true}) {
            const Surface a{Patch(0.0, rationalA)};
            const Surface b{Patch(10.0, rationalB)};
            for (const Box2 box : {Box2::Everything(), middle}) {
                const auto forward = Intersect(a, box, b, box, 1e-6, ctx);
                const auto reverse = Intersect(b, box, a, box, 1e-6, ctx);
                REQUIRE(forward.has_value());
                REQUIRE(reverse.has_value());
                REQUIRE(forward->IsEmpty());
                REQUIRE(reverse->IsEmpty());
            }
        }
    }
}

TEST_CASE("touching or overlapping NURBS hulls never masquerade as proved separation",
          "[intersect][generic][broadphase]")
{
    const Surface base{Patch(0.0)};
    for (const double z : {0.0, 1e-7}) {
        const Surface other{Patch(z)};
        const auto result = Intersect(base, middle, other, middle, 1e-6, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    const auto rational = Intersect(base, middle, Surface{Patch(0.0, true)}, middle, 1e-6, ctx);
    REQUIRE_FALSE(rational.has_value());
    REQUIRE(rational.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("tolerance prevents premature disjoint classification",
          "[intersect][generic][broadphase]")
{
    const Surface a{Patch(0.0)};
    const Surface b{Patch(2.0)};
    const auto separate = Intersect(a, middle, b, middle, 0.01, ctx);
    REQUIRE(separate.has_value());
    REQUIRE(separate->IsEmpty());
    const auto uncertain = Intersect(a, middle, b, middle, 3.0, ctx);
    REQUIRE_FALSE(uncertain.has_value());
    REQUIRE(uncertain.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("malformed rational control nets fail rather than certify separation",
          "[intersect][generic][broadphase]")
{
    const Surface good{Patch(10.0)};
    BSplineSurface bad = Patch(0.0, true);
    bad.Weights[0] = -1.0;
    const auto negative = Intersect(Surface{bad}, good, ctx);
    REQUIRE_FALSE(negative.has_value());
    REQUIRE(negative.error().GetCode() == ErrorCode::InvalidInput);

    bad = Patch(0.0);
    bad.Poles[0].Z = std::numeric_limits<double>::quiet_NaN();
    // Use an explicit requested box: a malformed aggregate may not have
    // a valid computed domain, but must never reach the broad-phase proof.
    const auto nan = Intersect(Surface{bad}, middle, good, middle, 0.0, ctx);
    REQUIRE_FALSE(nan.has_value());
    REQUIRE(nan.error().GetCode() == ErrorCode::InvalidInput);
}

namespace {

/// A biquadratic bowl over the unit square, curved along u: its height is
/// lift + 1 at u = 0 and 1 and lift + 1/2 - dip/2 at u = 1/2, while its net
/// reaches down to lift - dip. The optional weights vary along v only, so
/// they do not move the height.
BSplineSurface Bowl(const double lift, const double dip = 2.0, const bool rational = false)
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            poles.push_back({i / 2.0, j / 2.0, lift + (i == 1 ? -dip : 1.0)});
            weights.push_back(j == 1 ? 2.0 : 1.0);
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles,
                                                  rational ? weights : std::vector<double>{});
    REQUIRE(made.has_value());
    return *made;
}

/// A bilinear sheet over the unit square, rising with u from low to high.
BSplineSurface Slope(const double low, const double high)
{
    const std::vector<Point3> poles{
        {0.0, 0.0, low}, {0.0, 1.0, low}, {1.0, 0.0, high}, {1.0, 1.0, high}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// A bicubic with the given heights, x and y running over [0, 3].
BSplineSurface Bicubic(const std::vector<double>& z, const double lift)
{
    std::vector<Point3> poles;
    for (std::size_t k = 0; k < z.size(); ++k) {
        poles.push_back({static_cast<double>(k / 4), static_cast<double>(k % 4), z[k] + lift});
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(3, 3, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = (x - 1/2)^2 + (y - 1/2)^2 over the unit square, biquadratic.
BSplineSurface Paraboloid()
{
    const std::vector<double> a{0.25, -0.25, 0.25};
    std::vector<Point3> poles;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            poles.push_back({i / 2.0, j / 2.0,
                             a[static_cast<std::size_t>(i)] + a[static_cast<std::size_t>(j)]});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = (u - 1/2)(v - 1/2) over the unit square, bilinear.
BSplineSurface Saddle()
{
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.25}, {0.0, 1.0, -0.25}, {1.0, 0.0, -0.25}, {1.0, 1.0, 0.25}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// A unit quarter cylinder about the y axis, exact: rational quadratic in u.
BSplineSurface QuarterCylinder()
{
    const double s = std::sqrt(0.5);
    const std::vector<Point3> poles{{1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0},
                                    {1.0, 1.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}};
    const std::vector<double> weights{1.0, 1.0, s, s, 1.0, 1.0};
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 1, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// How far a parameter point of A is from the nearest sampled point of any
/// branch's pcurve on A.
double NearestOnA(const Stereon::SurfaceSurfaceIntersection& r, const Stereon::Point2& uv)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 256)) {
            const Stereon::Point2 q = Eval(branch.OnA, t);
            nearest = std::min(nearest, std::hypot(q.X - uv.X, q.Y - uv.Y));
        }
    }
    return nearest;
}

/// Over a grid of a pair sharing x and y parameters, wherever the height
/// difference changes sign along a grid line a branch must pass nearby: an
/// oracle that owes nothing to the search.
bool EverySignChangeIsFound(const Surface& a, const Surface& b,
                            const Stereon::SurfaceSurfaceIntersection& r)
{
    constexpr int n = 8;
    const double step = 1.0 / n;
    const auto gap = [&](const int i, const int j) {
        const Stereon::Point2 uv{i * step, j * step};
        return Eval(a, uv).Z - Eval(b, uv).Z;
    };
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (gap(i, j) * gap(i, j + 1) < 0.0 &&
                NearestOnA(r, Stereon::Point2{i * step, (j + 0.5) * step}) > step) {
                return false;
            }
            if (gap(j, i) * gap(j + 1, i) < 0.0 &&
                NearestOnA(r, Stereon::Point2{(j + 0.5) * step, i * step}) > step) {
                return false;
            }
        }
    }
    return true;
}

} // namespace

TEST_CASE("subdivision proves apart surfaces whose whole nets overlap",
          "[intersect][generic][subdivision]")
{
    for (const bool rational : {false, true}) {
        const Surface a{Bowl(0.0, 0.5, rational)};
        const Surface b{Bowl(0.5, 0.5, rational)};
        // The nets' boxes share z in [0, 1]; the surfaces stay 0.5 apart in
        // height everywhere.
        const auto forward = Intersect(a, b, ctx);
        const auto reverse = Intersect(b, a, ctx);
        INFO((forward ? std::string() : std::string(forward.error().GetMessage())));
        REQUIRE(forward.has_value());
        REQUIRE(reverse.has_value());
        REQUIRE(forward->IsEmpty());
        REQUIRE(reverse->IsEmpty());
    }
}

TEST_CASE("the requested domains decide what must be separated",
          "[intersect][generic][subdivision]")
{
    const Surface flat{Patch(0.0)};
    const Surface slope{Slope(-0.5, 0.5)}; // Through z = 0 at u = 1/2.
    const Box2 highU{Interval{0.75, 1.0}, Interval{0.0, 1.0}};
    const auto away = Intersect(flat, Box2::Everything(), slope, highU, 1e-6, ctx);
    REQUIRE(away.has_value());
    REQUIRE(away->IsEmpty());
    const auto swapped = Intersect(slope, highU, flat, Box2::Everything(), 1e-6, ctx);
    REQUIRE(swapped.has_value());
    REQUIRE(swapped->IsEmpty());

    const Box2 acrossU{Interval{0.25, 0.75}, Interval{0.0, 1.0}};
    const auto crossing = Intersect(flat, Box2::Everything(), slope, acrossU, 1e-6, ctx);
    INFO(IntersectTest::ErrorText(crossing));
    REQUIRE(crossing.has_value());
    REQUIRE(crossing->Branches.size() == 1);
    for (const double t : IntersectTest::SamplesOf(crossing->Branches[0].Range, 40)) {
        const Point3 p = Eval(crossing->Branches[0].Curve3d, t);
        REQUIRE(std::abs(p.X - 0.5) <= 2e-5);
        REQUIRE(std::abs(p.Z) <= 2e-5);
    }

    // Domains that miss the surface altogether leave nothing to intersect.
    const Box2 outside{Interval{2.0, 3.0}, Interval{0.0, 1.0}};
    const auto none = Intersect(flat, Box2::Everything(), slope, outside, 1e-6, ctx);
    REQUIRE(none.has_value());
    REQUIRE(none->IsEmpty());
}

TEST_CASE("a surface dipping through another is found on both sides of the dip",
          "[intersect][generic][marching]")
{
    // Lowest height 1/2 - 1.002/2 = -0.001: only a sliver crosses z = 0, in
    // two lines either side of u = 1/2.
    const Surface bowl{Bowl(0.0, 1.002)};
    const Surface flat{Patch(0.0)};
    const auto result = Intersect(flat, bowl, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 2);
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, flat, bowl, 100) <= branch.Bound + 1e-9);
    }
}

TEST_CASE("a plane cuts a paraboloid in one closed loop", "[intersect][generic][marching]")
{
    constexpr double h = 0.1;
    const Surface bowl{Paraboloid()};
    const Surface flat{Patch(h)};
    for (const bool swapped : {false, true}) {
        const auto result = swapped ? Intersect(bowl, flat, ctx) : Intersect(flat, bowl, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.size() == 1);
        REQUIRE(IntersectTest::ConnectivityOf(*result, 1e-5).Curves == 1);
        const Stereon::IntersectionBranch& branch = result->Branches[0];
        REQUIRE(branch.IsCertified());
        REQUIRE(Distance(Eval(branch.Curve3d, branch.Range.GetLower()),
                         Eval(branch.Curve3d, branch.Range.GetUpper())) <= 1e-5);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 200)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(std::abs(std::hypot(p.X - 0.5, p.Y - 0.5) - std::sqrt(h)) <= 1e-4);
            REQUIRE(std::abs(p.Z - h) <= 2e-5);
        }
    }
}

TEST_CASE("a rational quarter cylinder is cut along its exact line",
          "[intersect][generic][marching]")
{
    const Surface quarter{QuarterCylinder()};
    const Surface flat{Patch(0.5)};
    const auto result = Intersect(quarter, flat, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    REQUIRE(result->Branches[0].IsCertified());
    for (const double t : IntersectTest::SamplesOf(result->Branches[0].Range, 40)) {
        const Point3 p = Eval(result->Branches[0].Curve3d, t);
        REQUIRE(std::abs(p.X - std::sqrt(0.75)) <= 2e-5);
        REQUIRE(std::abs(p.Z - 0.5) <= 2e-5);
    }
}

TEST_CASE("near extrema never replace closed loops with isolated touches",
          "[intersect][generic][regression]")
{
    const Surface bowl{Paraboloid()};
    for (const double height : {0.0, 5e-7}) {
        const auto result = Intersect(bowl, Box2::Everything(), Surface{Patch(height)},
                                      Box2::Everything(), 1e-6, ctx);
        if (result) {
            REQUIRE(height > 0.0);
            REQUIRE_FALSE(result->Branches.empty());
            REQUIRE(result->Points.empty());
        } else {
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
}

TEST_CASE("tangency along a curve and a branch point fail closed", "[intersect][generic][special]")
{
    const Surface flat{Patch(0.0)};
    // Height exactly 0 along u = 1/2: tangent to the plane along a line.
    const auto tangent = Intersect(flat, Surface{Bowl(0.0, 1.0)}, ctx);
    REQUIRE_FALSE(tangent.has_value());
    REQUIRE(tangent.error().GetCode() == ErrorCode::IntersectionFailed);
    // z = (u - 1/2)(v - 1/2): two lines crossing at the centre.
    const auto node = Intersect(flat, Surface{Saddle()}, ctx);
    REQUIRE_FALSE(node.has_value());
    REQUIRE(node.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a bilinear pair's transverse line certifies against B-spline surfaces",
          "[intersect][generic][certificate]")
{
    const Surface flat{Patch(0.0)};
    const Surface slope{Slope(-0.5, 0.5)};
    const auto result = Intersect(slope, flat, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = result->Branches[0];
    REQUIRE(branch.IsCertified());
    REQUIRE(branch.Bound <= ctx.GetResolution());
}

TEST_CASE("branches of random bicubic pairs are complete and on both surfaces",
          "[intersect][generic][marching][property]")
{
    rc::prop("sign changes of the height difference all have a branch beside them", [] {
        std::vector<double> z(16);
        std::vector<double> w(16);
        for (std::size_t k = 0; k < 16; ++k) {
            z[k] = *rc::gen::inRange(-16, 17) / 16.0;
            w[k] = *rc::gen::inRange(-16, 17) / 16.0;
        }
        const double lift = *rc::gen::inRange(0, 25) / 8.0;
        const Surface a{Bicubic(z, 0.0)};
        const Surface b{Bicubic(w, lift)};
        const auto result = Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-6, ctx);
        if (!result) {
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        for (const Stereon::IntersectionBranch& branch : result->Branches) {
            RC_ASSERT(branch.IsCertified());
            RC_ASSERT(IntersectTest::PcurveDeviation(branch, a, b, 32) <= branch.Bound + 1e-9);
        }
        RC_ASSERT(EverySignChangeIsFound(a, b, *result));
    });
}
TEST_CASE("direct NURBS construction cannot bypass factory weight bounds",
          "[intersect][generic][broadphase][regression]")
{
    BSplineSurface invalid = Patch(0.0, true);
    invalid.Weights[0] = 0x1p257;
    const auto result =
        Intersect(Surface{invalid}, middle, Surface{Patch(100.0)}, middle, 0.0, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("initial SSI patch pairs obey the search budget", "[intersect][generic][regression]")
{
    constexpr int count = 130;
    std::vector<double> knots{0.0};
    for (int i = 0; i < count; ++i) {
        knots.push_back(static_cast<double>(i) / (count - 1));
    }
    knots.push_back(1.0);
    std::vector<Point3> poles;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            poles.push_back(
                {static_cast<double>(i) / (count - 1), static_cast<double>(j) / (count - 1), 0.0});
        }
    }
    const auto many = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles);
    REQUIRE(many.has_value());
    const auto result = Intersect(Surface{*many}, Surface{Patch(0.0)}, ctx);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
}
