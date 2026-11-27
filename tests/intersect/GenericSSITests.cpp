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
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
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
        INFO(IntersectTest::ErrorText(result));
        REQUIRE((result.has_value() || height > 0.0));
        if (result && height == 0.0) {
            // An exact touch at the bowl's bottom: one proved point, no
            // branch (the owner's decision of 2026-11-26 asks it answered).
            REQUIRE(result->Branches.empty());
            REQUIRE(result->Points.size() == 1);
            REQUIRE(result->Points.front().Bound <= ctx.GetResolution());
        } else if (result) {
            // A loop some 1.4e-3 across, wider than the resolution: never a
            // point.
            REQUIRE_FALSE(result->Branches.empty());
            REQUIRE(result->Points.empty());
        } else {
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
}

TEST_CASE("tangency along a curve fails closed; a transverse node is split",
          "[intersect][generic][special]")
{
    const Surface flat{Patch(0.0)};
    // Height exactly 0 along u = 1/2: tangent to the plane along a line whose
    // ends lie where both domains end, which the trace does not settle
    // (GenericSSITangentTests.cpp answers the tangencies it does).
    const auto tangent = Intersect(flat, Surface{Bowl(0.0, 1.0)}, ctx);
    REQUIRE_FALSE(tangent.has_value());
    REQUIRE(tangent.error().GetCode() == ErrorCode::IntersectionFailed);
    // z = (u - 1/2)(v - 1/2): two lines crossing at the centre, four branches
    // ending there (GenericSSINodeTests.cpp holds them to the closed form).
    const auto node = Intersect(flat, Surface{Saddle()}, ctx);
    INFO(IntersectTest::ErrorText(node));
    REQUIRE(node.has_value());
    REQUIRE(node->Branches.size() == 4);
    REQUIRE(IntersectTest::ConnectivityOf(*node, 1e-9).Curves == 1);
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

namespace {

/// A rigid placement: local (x, y, z) to x E1 + y E2 + z E3 + O.
struct Placement {
    Stereon::Vec3 E1;
    Stereon::Vec3 E2;
    Stereon::Vec3 E3;
    Point3 O;
};

const Placement atOrigin{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}};
/// The corpus's tilted placement: a rational rotation, orthonormal to rounding.
const Placement tilted{{0.6, 0.8, 0.0}, {-0.48, 0.36, 0.8}, {0.64, -0.48, 0.6}, {12.0, -7.0, 5.0}};

Point3 Placed(const Placement& at, const double x, const double y, const double z)
{
    return at.O + x * at.E1 + y * at.E2 + z * at.E3;
}

/// A point's local coordinates under a placement.
Stereon::Vec3 Local(const Placement& at, const Point3& p)
{
    const Stereon::Vec3 d = p - at.O;
    return {Dot(d, at.E1), Dot(d, at.E2), Dot(d, at.E3)};
}

/// The unit sphere, exact: a rational quadratic full circle in u, a half
/// circle from the south pole to the north one in v, the first and last rows
/// of v collapsed to the poles.
BSplineSurface UnitSphere(const Placement& at)
{
    const double s = std::sqrt(0.5);
    const std::vector<std::array<double, 2>> circle{{1, 0},   {1, 1},  {0, 1},  {-1, 1}, {-1, 0},
                                                    {-1, -1}, {0, -1}, {1, -1}, {1, 0}};
    const std::vector<double> circleWeights{1.0, s, 1.0, s, 1.0, s, 1.0, s, 1.0};
    const std::vector<std::array<double, 2>> half{{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}};
    const std::vector<double> halfWeights{1.0, s, 1.0, s, 1.0};
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < circle.size(); ++i) {
        for (std::size_t j = 0; j < half.size(); ++j) {
            const double rho = half[j][0];
            poles.push_back(Placed(at, rho * circle[i][0], rho * circle[i][1], half[j][1]));
            weights.push_back(circleWeights[i] * halfWeights[j]);
        }
    }
    const std::vector<double> knotsU{0.0, 0.0,  0.0,  0.25, 0.25, 0.5,
                                     0.5, 0.75, 0.75, 1.0,  1.0,  1.0};
    const std::vector<double> knotsV{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// A plane in the placement's local coordinates: n . x = Offset for
/// n = (sin Tilt, 0, cos Tilt), so z = Offset when Tilt is zero.
struct Cut {
    double Tilt;
    double Offset;
};

/// The square of side 4 in the cut's plane about the foot of the origin,
/// bilinear, its sides along (cos Tilt, 0, -sin Tilt) and y.
BSplineSurface Square(const Cut& cut, const Placement& at)
{
    const double c = std::cos(cut.Tilt);
    const double s = std::sin(cut.Tilt);
    const auto corner = [&](const double x, const double y) {
        return Placed(at, x * c + cut.Offset * s, y, cut.Offset * c - x * s);
    };
    const std::vector<Point3> poles{corner(-2.0, -2.0), corner(-2.0, 2.0), corner(2.0, -2.0),
                                    corner(2.0, 2.0)};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// What a result must show of the cut's section of the unit sphere, the
/// circle of radius sqrt(1 - Offset^2) about Offset n: every sampled point
/// of every branch on the circle within the tolerance, and the branches
/// together going `turn` round it, so that no arc is missing.
void CheckSection(const Stereon::SurfaceSurfaceIntersection& r, const Cut& cut, const Placement& at,
                  const double tolerance, const double turn)
{
    const double c = std::cos(cut.Tilt);
    const double s = std::sin(cut.Tilt);
    const double radius = std::sqrt(1.0 - cut.Offset * cut.Offset);
    REQUIRE(r.Points.empty());
    REQUIRE(r.Overlaps.empty());
    REQUIRE(IntersectTest::ConnectivityOf(r, 2.0 * tolerance).Curves == 1);
    double turned = 0.0;
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        REQUIRE(branch.Contact == Stereon::ContactType::Crossing);
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= ctx.GetResolution());
        double previous = 0.0;
        bool first = true;
        for (const double t : IntersectTest::SamplesOf(branch.Range, 200)) {
            const Stereon::Vec3 p = Local(at, Eval(branch.Curve3d, t));
            const double across = p.X * c - p.Z * s;
            REQUIRE(std::abs(std::hypot(across, p.Y) - radius) <= tolerance);
            REQUIRE(std::abs(p.X * s + p.Z * c - cut.Offset) <= tolerance);
            const double angle = std::atan2(p.Y, across);
            if (!first) {
                turned += std::abs(std::remainder(angle - previous, 2.0 * std::numbers::pi));
            }
            previous = angle;
            first = false;
        }
    }
    // An arc's ends are placed within the resolution, an angle of that over
    // the radius each.
    REQUIRE(std::abs(turned - turn) <= 1e-6 + 2.0 * ctx.GetResolution() / radius);
}

} // namespace

// The plane z = 0.999 meets the unit sphere in the circle of radius
// sqrt(1 - 0.999^2) = 0.0447 about the north pole, transverse everywhere at
// acos(0.999) = 2.56 degrees, along the latitude v* = 0.98433 of the Bezier
// patch that ends at the collapsed pole row. The partials of the sphere's
// rational patches, enclosed by the quotient rule about the origin, were
// widened by the patch's distance from it, one unit, times the width of the
// weight's derivative: S_u . z, zero on the sphere, enclosed as +-0.04 on a
// cell 0.002 wide, where the normals' 2.56 degrees part only by 0.045. No
// advance direction and no boundary root was proved near the circle, and the
// search cut pairs down to the tolerance and refused them as a tangential
// contact (corpus/nurbs-ssi/sphere-plane/near-pole-*). Tilted off the world
// axes, the boxes of the partials wrapped them further and the search ran
// out of its budget. MR: the same pair under the corpus's rigid tilt, and
// with the operands exchanged.
TEST_CASE("a plane a thousandth below a rational sphere's pole cuts its small circle",
          "[intersect][generic][nurbs][regression]")
{
    constexpr Cut cut{0.0, 1.0 - 1e-3};
    constexpr double tolerance = 1e-5;
    for (const Placement* at : {&atOrigin, &tilted}) {
        const Surface sphere{UnitSphere(*at)};
        const Surface flat{Square(cut, *at)};
        for (const bool swapped : {false, true}) {
            INFO((at == &tilted ? "tilted, " : "at the origin, ")
                 << (swapped ? "plane first" : "sphere first"));
            const auto result = swapped ? Intersect(flat, Box2::Everything(), sphere,
                                                    Box2::Everything(), tolerance, ctx)
                                        : Intersect(sphere, Box2::Everything(), flat,
                                                    Box2::Everything(), tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            CheckSection(*result, cut, *at, tolerance, 2.0 * std::numbers::pi);
        }
    }
}

TEST_CASE("a plane a thousandth above a rational sphere's pole is proved apart",
          "[intersect][generic][nurbs][regression]")
{
    const Surface sphere{UnitSphere(atOrigin)};
    const Surface flat{Square({0.0, 1.0 + 1e-3}, atOrigin)};
    const auto result = Intersect(sphere, Box2::Everything(), flat, Box2::Everything(), 1e-5, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Points.empty());
}

// The whole sphere, not only the fixture's latitude: both caps, both
// hemispheres, a circle a thousandth off the equator's knot line and one in
// a plane a radian off the axis, which crosses every patch obliquely. The
// two a thousandth from the poles are refused without the enclosures about
// the patch; the others were answered before and must stay so.
TEST_CASE("planes cut a rational sphere in circles from pole to pole",
          "[intersect][generic][nurbs]")
{
    constexpr double tolerance = 1e-5;
    const Surface sphere{UnitSphere(atOrigin)};
    for (const Cut& cut :
         {Cut{0.0, -1.0 + 1e-3}, Cut{0.0, -0.9}, Cut{0.0, -0.5}, Cut{0.0, 1e-3}, Cut{0.0, 0.25},
          Cut{0.0, 0.5}, Cut{0.0, 0.9}, Cut{0.0, 0.99}, Cut{1.0, 0.2}}) {
        INFO("tilt " << cut.Tilt << ", offset " << cut.Offset);
        const Surface flat{Square(cut, atOrigin)};
        const auto result =
            Intersect(sphere, Box2::Everything(), flat, Box2::Everything(), tolerance, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        CheckSection(*result, cut, atOrigin, tolerance, 2.0 * std::numbers::pi);
    }
}

// Sections the search does not resolve yet, each a transverse circle: three
// nearer the pole than the fixture (radii 0.0141, 0.0045 and 0.00045, where
// the sphere's partial along u is as short), the equator, along a knot line,
// a circle through the corner of four patches at (1, 0, 0), one through the
// pole, where the chart is singular, and the fixture's circle tipped 0.05
// radians, which passes 0.0053 from the pole.
// Whatever is answered is the circle, whole; anything else is refused.
TEST_CASE("a sphere's section that is not resolved is refused, never answered in part",
          "[intersect][generic][nurbs][safety]")
{
    constexpr double tolerance = 1e-5;
    const double third = std::numbers::pi / 6.0;
    const Surface sphere{UnitSphere(atOrigin)};
    for (const Cut& cut :
         {Cut{0.0, 1.0 - 1e-4}, Cut{0.0, 1.0 - 1e-5}, Cut{0.0, 1.0 - 1e-7}, Cut{0.0, 0.0},
          Cut{third, 0.5}, Cut{third, std::cos(third)}, Cut{0.05, 0.999}}) {
        INFO("tilt " << cut.Tilt << ", offset " << cut.Offset);
        const Surface flat{Square(cut, atOrigin)};
        const auto result =
            Intersect(sphere, Box2::Everything(), flat, Box2::Everything(), tolerance, ctx);
        if (result) {
            CheckSection(*result, cut, atOrigin, tolerance, 2.0 * std::numbers::pi);
        } else {
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
}

// The requested domain decides how much of the small circle there is: half
// the sphere's turn, u from 0.1 to 0.6, holds half of it, between two
// antipodal points on the domain's edges (each quarter turn of u is
// parameterized alike), and the sphere below the latitude v = 0.9
// (z = 0.956) none. This exercises cuts away from the interior knots;
// the exact-knot boundary regression below exercises those separately.
TEST_CASE("the small circle near the pole follows the sphere's requested domain",
          "[intersect][generic][nurbs]")
{
    constexpr Cut cut{0.0, 1.0 - 1e-3};
    constexpr double tolerance = 1e-5;
    const Surface sphere{UnitSphere(atOrigin)};
    const Surface flat{Square(cut, atOrigin)};
    const Box2 half{Interval{0.1, 0.6}, Interval{0.0, 1.0}};
    const auto arc = Intersect(sphere, half, flat, Box2::Everything(), tolerance, ctx);
    INFO(IntersectTest::ErrorText(arc));
    REQUIRE(arc.has_value());
    CheckSection(*arc, cut, atOrigin, tolerance, std::numbers::pi);
    const Box2 below{Interval{0.0, 1.0}, Interval{0.0, 0.9}};
    const auto none = Intersect(sphere, below, flat, Box2::Everything(), tolerance, ctx);
    INFO(IntersectTest::ErrorText(none));
    REQUIRE(none.has_value());
    REQUIRE(none->Branches.empty());
    REQUIRE(none->Points.empty());
}

// The old PatchesOver treated the right-hand patch touching u=0.5 as an
// entire searched patch, so a requested hemisphere answered up to u=0.75.
// Every successful result must stay in its own requested chart box, under
// both operand orders. The closed-form circle must still be half complete.
TEST_CASE("a rational sphere's exact interior knot bounds the searched SSI domain",
          "[intersect][generic][nurbs][regression][metamorphic]")
{
    constexpr Cut cut{0.0, 1.0 - 1e-3};
    constexpr double tolerance = 1e-5;
    const Surface sphere{UnitSphere(atOrigin)};
    const Surface flat{Square(cut, atOrigin)};
    for (const Interval u : {Interval{0.0, 0.5}, Interval{0.5, 1.0}}) {
        const Box2 wanted{u, Interval{0.0, 1.0}};
        for (const bool swapped : {false, true}) {
            INFO("requested u [" << u.GetLower() << ", " << u.GetUpper() << "]"
                                 << (swapped ? " plane first" : " sphere first"));
            const auto result =
                swapped ? Intersect(flat, Box2::Everything(), sphere, wanted, tolerance, ctx)
                        : Intersect(sphere, wanted, flat, Box2::Everything(), tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            CheckSection(*result, cut, atOrigin, tolerance, std::numbers::pi);
            for (const Stereon::IntersectionBranch& branch : result->Branches) {
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Stereon::Point2 uv = Eval(swapped ? branch.OnB : branch.OnA, t);
                    CHECK(uv.X >= u.GetLower() - 1e-10);
                    CHECK(uv.X <= u.GetUpper() + 1e-10);
                    CHECK(uv.Y >= -1e-10);
                    CHECK(uv.Y <= 1.0 + 1e-10);
                    CHECK(Distance(Eval(branch.Curve3d, t), Eval(sphere, uv)) <= tolerance);
                }
            }
        }
    }
}
