// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The certificates of fitted branches (issue #155): the bisection of a
// B-spline pcurve starts at the knots and bounds each span by the smaller of
// a mean-value form and a Taylor form whose remainder alone is enclosed over
// the span, from the fourth derivative. Whatever the bound's form, it
// must never fall below the deviation actually sampled, and moving a pcurve
// off its curve by more than the target must still be refuted.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::BSplineCurve2d;
using Stereon::Curve2d;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Surface aboutZ{Cylinder{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
const Surface tilted{Plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.8, 0.6}}};
const Surface skewCrossing{
    Cylinder{Point3{1.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.6, 0.8}, true}};

/// The largest |c(t) - S(p(t))| at n + 1 parameters of the range.
double Sampled(const Stereon::IntersectionBranch& b, const Surface& s, const Curve2d& p,
               const int n = 20000)
{
    double worst = 0.0;
    for (int k = 0; k <= n; ++k) {
        const double t = b.Range.GetLower() + Width(b.Range) * k / n;
        worst = std::max(worst, Distance(Eval(b.Curve3d, t), Eval(s, Eval(p, t))));
    }
    return worst;
}

/// The pcurve's B-spline. A throw, not a dereference of null, where it is
/// not one: the throw is a path the optimiser sees end, so -O3 has no null
/// dereference to warn of, and the test fails with a message.
const BSplineCurve2d& SplineOf(const Curve2d& pcurve)
{
    const BSplineCurve2d* spline = pcurve.Target<BSplineCurve2d>();
    if (spline == nullptr) {
        throw std::runtime_error("the pcurve is not a B-spline");
    }
    return *spline;
}

/// The B-spline pcurve with every pole moved by `dv` along v.
Curve2d MovedAlongV(const Curve2d& pcurve, const double dv)
{
    BSplineCurve2d spline = SplineOf(pcurve);
    for (Point2& pole : spline.Poles) {
        pole.Y += dv;
    }
    return Curve2d{spline};
}

} // namespace

TEST_CASE("fitted branches are certified above their sampled deviation", "[intersect][certify]")
{
    for (const auto& [a, b] : {std::pair{&tilted, &aboutZ}, std::pair{&aboutZ, &skewCrossing}}) {
        const auto result = Intersect(*a, *b, ctx);
        REQUIRE(result.has_value());
        REQUIRE_FALSE(result->Branches.empty());
        for (const Stereon::IntersectionBranch& branch : result->Branches) {
            REQUIRE(branch.Range.IsBounded());
            for (const auto& [surface, pcurve] :
                 {std::pair{a, &branch.OnA}, std::pair{b, &branch.OnB}}) {
                const auto cert = CertifyDeviation(branch.Curve3d, *surface, *pcurve, branch.Range,
                                                   ctx.GetResolution());
                REQUIRE(cert.has_value());
                REQUIRE(cert->Bound <= ctx.GetResolution());
                REQUIRE(cert->Bound >= Sampled(branch, *surface, *pcurve));
            }
        }
    }
}

TEST_CASE("a fitted pcurve moved off its curve is refuted or bounded above the move",
          "[intersect][certify]")
{
    // The ellipse's pcurve on the cylinder is a cubic B-spline in (u, v);
    // moving it along v moves its image along the axis by as much.
    const auto result = Intersect(tilted, aboutZ, ctx);
    REQUIRE(result.has_value());
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    REQUIRE(branch.OnB.Target<BSplineCurve2d>() != nullptr);
    const double target = ctx.GetResolution();
    for (const double dv : {2e-6, 5e-6, 2e-5, 1e-3}) {
        INFO("moved by " << dv);
        const Curve2d moved = MovedAlongV(branch.OnB, dv);
        const double sampled = Sampled(branch, aboutZ, moved);
        const auto cert = CertifyDeviation(branch.Curve3d, aboutZ, moved, branch.Range, target);
        if (sampled > target) {
            REQUIRE_FALSE(cert.has_value());
            REQUIRE(cert.error().GetCode() == ErrorCode::CertificationFailed);
        } else if (cert) {
            REQUIRE(cert->Bound >= sampled);
            REQUIRE(cert->Bound <= target);
        }
        // A loose target always holds and is never below the sample.
        const auto loose = CertifyDeviation(branch.Curve3d, aboutZ, moved, branch.Range, 1e-2);
        REQUIRE(loose.has_value());
        REQUIRE(loose->Bound >= sampled);
    }
}

TEST_CASE("each knot span of a fitted branch is certified whole", "[intersect][certify]")
{
    // The mean-value form alone needed some forty spans per knot span, and
    // the Taylor form with its third derivative enclosed over the span two;
    // with the fourth enclosed, one (docs/performance.rst).
    const auto result = Intersect(tilted, aboutZ, ctx);
    REQUIRE(result.has_value());
    const Stereon::IntersectionBranch& branch = result->Branches.front();
    const BSplineCurve2d& spline = SplineOf(branch.OnB);
    std::vector<double> knots(spline.Knots.begin(), spline.Knots.end());
    knots.erase(std::ranges::unique(knots).begin(), knots.end());
    const auto cert =
        CertifyDeviation(branch.Curve3d, aboutZ, branch.OnB, branch.Range, ctx.GetResolution());
    REQUIRE(cert.has_value());
    REQUIRE(cert->Spans <= knots.size() - 1);
}
