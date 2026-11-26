// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The products of a surface–surface intersection at its nodes (ADR-0021,
// section 9; ADR-0054; the v0.2 linearisation gate). A node, where the
// surfaces are tangent and two branches cross, is where the handlers end
// four branches: two equal cylinders whose axes meet (two ellipses through
// two nodes) and a transverse node of two B-spline surfaces (#337). A
// sample there is a candidate event, not a point that moves: IntersectJvp,
// IntersectVjp and Linearize on the whole result answer Nondifferentiable
// and carry no values. The nodes are held to closed forms written here, not
// read from the kernel. Samples away from the node keep their derivatives,
// and two ends meeting at a seam or a cut are no node.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::BranchSample;
using Stereon::BSplineSurface;
using Stereon::Cylinder;
using Stereon::DerivativeStatus;
using Stereon::ErrorCode;
using Stereon::IntersectionBranch;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceBinding;
using Stereon::SurfaceComponent;
using Stereon::SurfaceSurfaceIntersection;
using Stereon::Vec3;

namespace {

/// One sample of a result, the status its products and Linearize give, and
/// whether they carry values.
struct Answer {
    DerivativeStatus Jvp;
    DerivativeStatus Vjp;
    DerivativeStatus Linearized;
    bool JvpValues;
    bool VjpValues;
};

/// The answers for the samples, with a parameter moving `bindingsA` on a.
Answer AnswerFor(const Surface& a, const std::vector<SurfaceBinding>& bindingsA, const Surface& b,
                 const SurfaceSurfaceIntersection& at, const std::vector<BranchSample>& samples,
                 SensitivityContext& sensitivity)
{
    const std::vector<SurfaceBinding> none;
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const auto jvp = Stereon::IntersectJvp(a, bindingsA, b, none, at, samples, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(sensitivity
                .SetAdjoint(std::vector<double>(Stereon::BranchPointValues * samples.size(), 1.0))
                .has_value());
    const auto vjp = Stereon::IntersectVjp(a, bindingsA, b, none, at, samples, sensitivity);
    INFO(IntersectTest::ErrorText(vjp));
    REQUIRE(vjp.has_value());
    DerivativeStatus worst = DerivativeStatus::Valid;
    for (const BranchSample& sample : samples) {
        const auto lin = Stereon::Linearize(a, b, at, sample);
        INFO(IntersectTest::ErrorText(lin));
        REQUIRE(lin.has_value());
        if (lin->Status != DerivativeStatus::Valid) {
            worst = lin->Status;
        }
    }
    return {jvp->GetStatus(), vjp->GetStatus(), worst, jvp->GetValues().has_value(),
            vjp->GetValues().has_value()};
}

/// Whether the branch's point at t is one of `nodes`, within `within`.
bool AtOneOf(const IntersectionBranch& branch, const double t, const std::vector<Point3>& nodes,
             const double within)
{
    const Point3 p = Eval(branch.Curve3d, t);
    for (const Point3& node : nodes) {
        if (Distance(p, node) <= within) {
            return true;
        }
    }
    return false;
}

/// Checks every end of every branch, and the middle of each: an end at one
/// of `nodes` (closed forms) is Nondifferentiable with no values, every
/// other sample keeps its derivative. Returns how many ends were at nodes.
std::size_t RequireNodesSingled(const Surface& a, const std::vector<SurfaceBinding>& bindingsA,
                                const Surface& b, const SurfaceSurfaceIntersection& at,
                                const std::vector<Point3>& nodes, const double within,
                                SensitivityContext& sensitivity)
{
    std::size_t atNodes = 0;
    for (std::size_t k = 0; k < at.Branches.size(); ++k) {
        const IntersectionBranch& branch = at.Branches[k];
        REQUIRE(branch.Range.IsBounded());
        const double lo = branch.Range.GetLower();
        const double hi = branch.Range.GetUpper();
        for (const double t : {lo, hi, 0.5 * (lo + hi)}) {
            INFO("branch " << k << " at " << t);
            const Answer answer = AnswerFor(a, bindingsA, b, at, {BranchSample{k, t}}, sensitivity);
            if (AtOneOf(branch, t, nodes, within)) {
                ++atNodes;
                CHECK(answer.Jvp == DerivativeStatus::Nondifferentiable);
                CHECK(answer.Vjp == DerivativeStatus::Nondifferentiable);
                CHECK(answer.Linearized == DerivativeStatus::Nondifferentiable);
                CHECK_FALSE(answer.JvpValues);
                CHECK_FALSE(answer.VjpValues);
            } else {
                CHECK(HasDerivative(answer.Jvp));
                CHECK(answer.Vjp == answer.Jvp);
                CHECK(answer.JvpValues);
                CHECK(answer.VjpValues);
            }
        }
    }
    return atNodes;
}

/// z = x' y' over the unit square, x' = x - 1/2 and y' = y - 1/2: bilinear,
/// crossing z = 0 along x' = 0 and y' = 0, its node at (1/2, 1/2, 0).
BSplineSurface Saddle()
{
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.25}, {0.0, 1.0, -0.25}, {1.0, 0.0, -0.25}, {1.0, 1.0, 0.25}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = 0, bilinear over the unit square.
BSplineSurface Floor()
{
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

} // namespace

TEST_CASE("a sample at a node of two equal cylinders has no derivative",
          "[intersect][linearization][branch][node]")
{
    // Radius 3 about z, and radius 3 about x: two ellipses in the planes
    // z = x and z = -x, through the nodes (0, +-3, 0), where the axes'
    // common normal meets both cylinders and they are tangent.
    const Surface upright{Cylinder{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
    const Surface lying{Cylinder{Point3{}, 3.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true}};
    const auto at = Stereon::Intersect(upright, Box2::Everything(), lying, Box2::Everything(), ctx);
    INFO(IntersectTest::ErrorText(at));
    REQUIRE(at.has_value());
    REQUIRE(IntersectTest::ConnectivityOf(*at).Nodes == 2);

    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    const std::vector<SurfaceBinding> radius{{p, SurfaceComponent::Radius}};
    const std::vector<Point3> nodes{Point3{0.0, 3.0, 0.0}, Point3{0.0, -3.0, 0.0}};
    // Four branch ends at each node.
    CHECK(RequireNodesSingled(upright, radius, lying, *at, nodes, 1e-6, sensitivity) == 8);

    SECTION("one sample at a node takes the whole result's derivative away")
    {
        std::size_t node = at->Branches.size();
        for (std::size_t k = 0; k < at->Branches.size(); ++k) {
            const IntersectionBranch& branch = at->Branches[k];
            if (AtOneOf(branch, branch.Range.GetLower(), nodes, 1e-6)) {
                node = k;
            }
        }
        REQUIRE(node < at->Branches.size());
        const IntersectionBranch& branch = at->Branches[node];
        const std::vector<BranchSample> samples{
            {node, 0.5 * (branch.Range.GetLower() + branch.Range.GetUpper())},
            {node, branch.Range.GetLower()}};
        const Answer answer = AnswerFor(upright, radius, lying, *at, samples, sensitivity);
        CHECK(answer.Jvp == DerivativeStatus::Nondifferentiable);
        CHECK(answer.Vjp == DerivativeStatus::Nondifferentiable);
        CHECK_FALSE(answer.JvpValues);
    }
    SECTION("the single-branch Linearize cannot see the node; the whole result's can")
    {
        for (std::size_t k = 0; k < at->Branches.size(); ++k) {
            const IntersectionBranch& branch = at->Branches[k];
            const double t = branch.Range.GetLower();
            if (!AtOneOf(branch, t, nodes, 1e-6)) {
                continue;
            }
            const auto whole = Stereon::Linearize(upright, lying, *at, BranchSample{k, t});
            REQUIRE(whole.has_value());
            CHECK(whole->Status == DerivativeStatus::Nondifferentiable);
            // The surfaces are tangent at the node: the system is singular,
            // or within rounding of it.
            CHECK(whole->Conditioning < 1e-6);
        }
    }
}

TEST_CASE("a sample at a transverse node of two B-spline surfaces has no derivative",
          "[intersect][linearization][branch][node]")
{
    // The saddle z = x' y' crosses the floor along x' = 0 and y' = 0: four
    // branches from the node (1/2, 1/2, 0) to the square's edges (#337).
    const Surface saddle{Saddle()};
    const Surface floor{Floor()};
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "floor first" : "saddle first"));
        const Surface& a = swapped ? floor : saddle;
        const Surface& b = swapped ? saddle : floor;
        const auto at = Stereon::Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-6, ctx);
        INFO(IntersectTest::ErrorText(at));
        REQUIRE(at.has_value());
        REQUIRE(at->Branches.size() == 4);
        ParameterSet set;
        const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
        SensitivityContext sensitivity(std::move(set));
        // The parameter lifts a's first pole.
        const std::vector<SurfaceBinding> lift{{p, SurfaceComponent::PoleZ, 0}};
        CHECK(RequireNodesSingled(a, lift, b, *at, {Point3{0.5, 0.5, 0.0}}, 1e-6, sensitivity) ==
              4);
    }
}

TEST_CASE("two branch ends meeting at a seam are no node",
          "[intersect][linearization][branch][node]")
{
    // The plane z = 1 cuts the sphere of radius 2 in a circle that its
    // handlers split at the sphere's seam: every end is a regular point.
    const Surface ball{Sphere{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Surface cut{Plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const auto at = Stereon::Intersect(cut, ball, ctx);
    INFO(IntersectTest::ErrorText(at));
    REQUIRE(at.has_value());
    REQUIRE_FALSE(at->Branches.empty());
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    const std::vector<SurfaceBinding> rise{{p, SurfaceComponent::OriginZ}};
    CHECK(RequireNodesSingled(cut, rise, ball, *at, {}, 0.0, sensitivity) == 0);
}

TEST_CASE("Linearize on a whole result refuses what it cannot read",
          "[intersect][linearization][branch][node]")
{
    const Surface ball{Sphere{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Surface cut{Plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const auto at = Stereon::Intersect(cut, ball, ctx);
    REQUIRE(at.has_value());
    REQUIRE_FALSE(at->Branches.empty());
    const auto beyond = Stereon::Linearize(cut, ball, *at, BranchSample{at->Branches.size(), 0.0});
    REQUIRE_FALSE(beyond.has_value());
    CHECK(beyond.error().GetCode() == ErrorCode::InvalidInput);
    const auto empty = Stereon::Linearize(Surface{}, ball, *at, BranchSample{0, 0.0});
    REQUIRE_FALSE(empty.has_value());
    CHECK(empty.error().GetCode() == ErrorCode::InvalidInput);
    const IntersectionBranch& branch = at->Branches.front();
    const auto off =
        Stereon::Linearize(cut, ball, *at, BranchSample{0, branch.Range.GetUpper() + 1.0});
    REQUIRE_FALSE(off.has_value());
    CHECK(off.error().GetCode() == ErrorCode::InvalidInput);
    // A regular sample agrees with the single-branch overload.
    const double t = 0.5 * (branch.Range.GetLower() + branch.Range.GetUpper());
    const auto whole = Stereon::Linearize(cut, ball, *at, BranchSample{0, t});
    const auto single = Stereon::Linearize(cut, ball, branch, t);
    REQUIRE(whole.has_value());
    REQUIRE(single.has_value());
    CHECK(whole->Status == single->Status);
    CHECK(whole->Position == single->Position);
    CHECK(whole->Conditioning == single->Conditioning);
}
