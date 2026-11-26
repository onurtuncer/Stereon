// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The products of a crossing where a B-spline is not smooth (ADR-0021,
// section 9; the v0.2 linearisation gate). Since #334 the planar handlers
// decide crossings at a spline's own end and at a corner exactly. There the
// point leaves the curve, or passes from one arm to the other, as the
// geometry moves: Nondifferentiable, never Valid, whatever the one-sided
// tangent the evaluation happens to use says. The same holds for a B-spline
// curve's end or corner against a surface, and for a point on a B-spline
// surface's boundary; a point on an interior knot line of a surface where
// its partials may jump is Indeterminate. A smooth joint (a simple knot, or
// a joint whose one-sided tangents agree) keeps its derivative.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::BSplineSurface;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::DerivativeStatus;
using Stereon::Line;
using Stereon::LineBinding;
using Stereon::LineComponent;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Plane;
using Stereon::PlaneCurveBinding;
using Stereon::PlaneCurveComponent;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::PoleBinding;
using Stereon::PoleComponent;
using Stereon::SensitivityContext;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

BSplineCurve2d Spline(const int degree, std::vector<double> knots, std::vector<Point2> poles)
{
    auto made = Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), {});
    REQUIRE(made.has_value());
    return *made;
}

BSplineCurve Spline3(const int degree, std::vector<double> knots, std::vector<Point3> poles)
{
    auto made = Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), {});
    REQUIRE(made.has_value());
    return *made;
}

const std::vector<double> cubic{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

/// The fixtures of tests/intersect/PlanarSplineEndTests.cpp: a cubic from
/// the origin to (1, 1), one from (1, 1) to (2, 0), and the segment x = 1
/// from (1, -1) to (1, 3).
BSplineCurve2d Arch()
{
    return Spline(3, cubic, {{0, 0}, {0.3, 0.6}, {0.7, 1.2}, {1, 1}});
}
BSplineCurve2d Descent()
{
    return Spline(3, cubic, {{1, 1}, {1.4, 0.7}, {1.7, 0.2}, {2, 0}});
}
BSplineCurve2d Post()
{
    return Spline(1, {0.0, 0.0, 1.0, 1.0}, {{1, -1}, {1, 3}});
}

/// y = |x| over [-1, 1], its corner at t = 1, and y = 0 over [-1, 1] as
/// a polyline jointed at the origin, smooth there (equal speeds).
BSplineCurve2d Vee()
{
    return Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 1}, {0, 0}, {1, 1}});
}
BSplineCurve2d Bar()
{
    return Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 0}, {0, 0}, {1, 0}});
}

/// The segment x = 0 from (0, -1) to (0, 3).
BSplineCurve2d Pole()
{
    return Spline(1, {0.0, 0.0, 1.0, 1.0}, {{0, -1}, {0, 3}});
}

/// The statuses of a planar result: the forward and reverse products, with
/// the first pole of a moving in x, and Linearize at each point (the worst).
struct Statuses {
    DerivativeStatus Jvp;
    DerivativeStatus Vjp;
    DerivativeStatus Point;
    bool Values;
};

Statuses PlanarStatuses(const BSplineCurve2d& a, const BSplineCurve2d& b,
                        const PlanarIntersection& at)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    REQUIRE(sensitivity
                .SetAdjoint(std::vector<double>(Stereon::PlanarPointValues * at.Points.size(), 1.0))
                .has_value());
    const std::vector<PlaneCurveBinding> bindings{{p, PlaneCurveComponent::PoleX, 0}};
    const std::vector<PlaneCurveBinding> none;
    const auto jvp = Stereon::IntersectJvp(Curve2d{a}, bindings, Curve2d{b}, none, at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    const auto vjp = Stereon::IntersectVjp(Curve2d{a}, bindings, Curve2d{b}, none, at, sensitivity);
    REQUIRE(vjp.has_value());
    DerivativeStatus point = DerivativeStatus::Valid;
    for (const PlanarPoint& m : at.Points) {
        const auto lin = Stereon::Linearize(Curve2d{a}, Curve2d{b}, m);
        INFO(IntersectTest::ErrorText(lin));
        REQUIRE(lin.has_value());
        if (lin->Status != DerivativeStatus::Valid) {
            point = lin->Status;
        }
    }
    return {jvp->GetStatus(), vjp->GetStatus(), point, jvp->GetValues().has_value()};
}

PlanarIntersection Meet(const BSplineCurve2d& a, const BSplineCurve2d& b)
{
    const auto r = Stereon::Intersect(Curve2d{a}, Curve2d{b}, ctx);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Overlaps.empty());
    return *r;
}

/// A decided crossing at p, as #334 reports it, in both orders, is
/// Nondifferentiable with no values.
void RequireNoDerivative(const BSplineCurve2d& a, const BSplineCurve2d& b, const Point2& p)
{
    for (const bool exchanged : {false, true}) {
        INFO((exchanged ? "exchanged" : "as given"));
        const BSplineCurve2d& first = exchanged ? b : a;
        const BSplineCurve2d& second = exchanged ? a : b;
        const PlanarIntersection at = Meet(first, second);
        REQUIRE(at.Points.size() == 1);
        REQUIRE(at.Points.front().Decided);
        REQUIRE(at.Points.front().Contact == ContactType::Crossing);
        CHECK(Distance(at.Points.front().Position, p) <= 1e-12);
        const Statuses s = PlanarStatuses(first, second, at);
        CHECK(s.Jvp == DerivativeStatus::Nondifferentiable);
        CHECK(s.Vjp == DerivativeStatus::Nondifferentiable);
        CHECK(s.Point == DerivativeStatus::Nondifferentiable);
        CHECK_FALSE(s.Values);
    }
}

/// The status of IntersectJvp for a line, its origin's z moving, against
/// `surface` at one hand-made crossing (t, uv, x), and of Linearize there.
std::pair<DerivativeStatus, DerivativeStatus> LineStatus(const Line& line, const Surface& surface,
                                                         const CurveSurfacePoint& point)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<LineBinding> bindings{{p, LineComponent::OriginX}};
    const CurveSurfaceIntersection at{{point}, {}};
    const auto jvp = Stereon::IntersectJvp(line, bindings, surface, {}, at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    CHECK(jvp->GetValues().has_value() == HasDerivative(jvp->GetStatus()));
    const auto lin = Stereon::Linearize(line, surface, point);
    INFO(IntersectTest::ErrorText(lin));
    REQUIRE(lin.has_value());
    return {jvp->GetStatus(), lin->Status};
}

/// The same for a B-spline curve, its first pole moving in x.
std::pair<DerivativeStatus, DerivativeStatus>
SplineStatus(const BSplineCurve& curve, const Surface& surface, const CurveSurfacePoint& point)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<PoleBinding> bindings{{p, 0, PoleComponent::X}};
    const CurveSurfaceIntersection at{{point}, {}};
    const auto jvp = Stereon::IntersectJvp(curve, bindings, surface, {}, at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    CHECK(jvp->GetValues().has_value() == HasDerivative(jvp->GetStatus()));
    const auto lin = Stereon::Linearize(curve, surface, point);
    INFO(IntersectTest::ErrorText(lin));
    REQUIRE(lin.has_value());
    return {jvp->GetStatus(), lin->Status};
}

/// z = 0 over [0, 1]^2, bilinear, with an interior knot line u = 1/2 when
/// `split`: a C0 joint of the patch, flat across it here.
BSplineSurface Floor(const bool split)
{
    std::vector<double> knotsU{0.0, 0.0, 1.0, 1.0};
    std::vector<Point3> poles{{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
    if (split) {
        knotsU = {0.0, 0.0, 0.5, 1.0, 1.0};
        poles = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.5, 0.0, 0.0},
                 {0.5, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
    }
    const auto made =
        Stereon::MakeBSplineSurface(1, 1, knotsU, {0.0, 0.0, 1.0, 1.0}, std::move(poles), {});
    REQUIRE(made.has_value());
    return *made;
}

/// The vertical line through (x, y, 0), at t = 0 there.
Line Upright(const double x, const double y)
{
    return Line{Point3{x, y, 0.0}, Vec3{0.0, 0.0, 1.0}};
}

CurveSurfacePoint Crossing(const double t, const Point2& uv, const Point3& x)
{
    return CurveSurfacePoint{t, uv, x, ContactType::Crossing, true};
}

const Surface ground{Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};

} // namespace

TEST_CASE("plane B-splines crossing at a spline's own end have no derivative",
          "[intersect][linearization][planar][spline-end]")
{
    // Arch ends where Descent starts, both ends at (1, 1); Arch ends on the
    // post's middle. Before, both were Valid, on the end segment's
    // polynomial continued past the end.
    RequireNoDerivative(Arch(), Descent(), Point2{1, 1});
    RequireNoDerivative(Arch(), Post(), Point2{1, 1});
}

TEST_CASE("a plane B-spline crossed at its corner has no derivative",
          "[intersect][linearization][planar][spline-end]")
{
    // The pole x = 0 passes through the vee's corner outside its cone: both
    // arms cross it the same way, so #334 decides one crossing at t = 1,
    // where the arms' tangents (1, -1) and (1, 1) differ.
    RequireNoDerivative(Vee(), Pole(), Point2{0, 0});
}

TEST_CASE("a plane B-spline crossed at a smooth joint keeps its derivative",
          "[intersect][linearization][planar][spline-end]")
{
    // The bar's joint at the origin has the tangent (1, 0) on both sides.
    for (const bool exchanged : {false, true}) {
        INFO((exchanged ? "exchanged" : "as given"));
        const BSplineCurve2d first = exchanged ? Pole() : Bar();
        const BSplineCurve2d second = exchanged ? Bar() : Pole();
        const PlanarIntersection at = Meet(first, second);
        REQUIRE(at.Points.size() == 1);
        REQUIRE(at.Points.front().Decided);
        const Statuses s = PlanarStatuses(first, second, at);
        CHECK(s.Jvp == DerivativeStatus::Valid);
        CHECK(s.Vjp == DerivativeStatus::Valid);
        CHECK(s.Point == DerivativeStatus::Valid);
        CHECK(s.Values);
    }
    // The joint itself, given exactly: smooth, so Valid.
    const PlanarIntersection exact{
        {PlanarPoint{1.0, 0.25, Point2{0, 0}, ContactType::Crossing, true}}, {}};
    CHECK(PlanarStatuses(Bar(), Pole(), exact).Jvp == DerivativeStatus::Valid);
    // The vee's corner given exactly: Nondifferentiable.
    CHECK(PlanarStatuses(Vee(), Pole(), exact).Jvp == DerivativeStatus::Nondifferentiable);
}

TEST_CASE("a space B-spline's end or corner on a surface has no derivative",
          "[intersect][linearization][spline-end]")
{
    // A polyline in the xz-plane cornered at the origin, its arms along
    // (1, 0, 1) and (1, 0, 2): it crosses z = 0 there, at t = 1.
    const BSplineCurve bent =
        Spline3(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 0, -1}, {0, 0, 0}, {1, 0, 2}});
    const auto [cornerJvp, cornerLin] =
        SplineStatus(bent, ground, Crossing(1.0, Point2{0, 0}, Point3{}));
    CHECK(cornerJvp == DerivativeStatus::Nondifferentiable);
    CHECK(cornerLin == DerivativeStatus::Nondifferentiable);
    // Off the corner, on one arm: Valid.
    const BSplineCurve straight =
        Spline3(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 0, -1}, {0, 0, 0}, {1, 0, 1}});
    const auto [smoothJvp, smoothLin] =
        SplineStatus(straight, ground, Crossing(1.0, Point2{0, 0}, Point3{}));
    CHECK(smoothJvp == DerivativeStatus::Valid);
    CHECK(smoothLin == DerivativeStatus::Valid);
    // The first arm alone, ending on the plane at t = 1.
    const BSplineCurve arm = Spline3(1, {0.0, 0.0, 1.0, 1.0}, {{-1, 0, -1}, {0, 0, 0}});
    const auto [endJvp, endLin] = SplineStatus(arm, ground, Crossing(1.0, Point2{0, 0}, Point3{}));
    CHECK(endJvp == DerivativeStatus::Nondifferentiable);
    CHECK(endLin == DerivativeStatus::Nondifferentiable);
    // The primal reports the corner's crossing; whatever it says, never Valid.
    const auto met = Stereon::Intersect(Stereon::Curve{bent}, ground, ctx);
    INFO(IntersectTest::ErrorText(met));
    if (met && met->Points.size() == 1) {
        CHECK(SplineStatus(bent, ground, met->Points.front()).first != DerivativeStatus::Valid);
    }
}

TEST_CASE("a point on a B-spline surface's boundary or knot line keeps no derivative",
          "[intersect][linearization][spline-end]")
{
    const Surface whole{Floor(false)};
    const Surface split{Floor(true)};
    // Inside: Valid.
    CHECK(LineStatus(Upright(0.25, 0.5), whole,
                     Crossing(0.0, Point2{0.25, 0.5}, Point3{0.25, 0.5, 0}))
              .first == DerivativeStatus::Valid);
    // On the boundary u = 1 and on the corner (0, 0): Nondifferentiable.
    const auto [edgeJvp, edgeLin] =
        LineStatus(Upright(1.0, 0.5), whole, Crossing(0.0, Point2{1.0, 0.5}, Point3{1.0, 0.5, 0}));
    CHECK(edgeJvp == DerivativeStatus::Nondifferentiable);
    CHECK(edgeLin == DerivativeStatus::Nondifferentiable);
    CHECK(LineStatus(Upright(0.0, 0.0), whole, Crossing(0.0, Point2{0.0, 0.0}, Point3{})).first ==
          DerivativeStatus::Nondifferentiable);
    // On the interior knot line u = 1/2: Indeterminate; beside it, Valid.
    const auto [knotJvp, knotLin] =
        LineStatus(Upright(0.5, 0.5), split, Crossing(0.0, Point2{0.5, 0.5}, Point3{0.5, 0.5, 0}));
    CHECK(knotJvp == DerivativeStatus::Indeterminate);
    CHECK(knotLin == DerivativeStatus::Indeterminate);
    CHECK(LineStatus(Upright(0.25, 0.5), split,
                     Crossing(0.0, Point2{0.25, 0.5}, Point3{0.25, 0.5, 0}))
              .first == DerivativeStatus::Valid);
}
