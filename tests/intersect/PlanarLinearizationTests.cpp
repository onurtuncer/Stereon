// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The linearisation of points where two plane curves meet (ADR-0021,
// section 9): closed forms for a line against a circle, an ellipse and a
// segment; the second-order decay of the Taylor residual and the adjoint
// identity on random pairs of lines, circles, ellipses and cubic B-splines,
// rational half the time; every status; and the bindings refused.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::BSplineCurve2d;
using Stereon::Circle2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::DerivativeStatus;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line2d;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::PlaneCurveBinding;
using Stereon::PlaneCurveComponent;
using Stereon::Point2;
using Stereon::SensitivityContext;
using Stereon::Vec2;

namespace {

constexpr std::size_t Values = Stereon::PlanarPointValues;
constexpr double tau = 2.0 * std::numbers::pi;

const Curve2d XAxis{Line2d{Point2{0.0, 0.0}, Vec2{1.0, 0.0}}};

PlanarIntersection Meet(const Curve2d& a, const Curve2d& b)
{
    auto result = Stereon::Intersect(a, b, ctx);
    if (!result) {
        RC_FAIL(IntersectTest::ErrorText(result));
    }
    return *std::move(result);
}

/// The forward product of one binding of `a` moved at unit rate.
std::vector<double> RateOf(const Curve2d& a, const PlaneCurveComponent component, const Curve2d& b,
                           const std::size_t pole = 0)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<PlaneCurveBinding> bindings{{p, component, pole}};
    const auto at = Stereon::Intersect(a, b, ctx);
    REQUIRE(at.has_value());
    const auto jvp = Stereon::IntersectJvp(a, bindings, b, {}, *at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

// --- Random geometry for the properties ---------------------------------------

double Coordinate()
{
    return *rc::gen::inRange(-16, 17) / 8.0;
}

/// A random unit direction.
Vec2 RandomDirection()
{
    const double angle = *rc::gen::inRange(0, 64) * tau / 64.0;
    return {std::cos(angle), std::sin(angle)};
}

/// One side: the curve, and the components its bindings move.
struct RandomSide {
    Curve2d Curve;
    std::vector<PlaneCurveBinding> Shape; // parameters assigned later
};

RandomSide RandomLine()
{
    const Vec2 d = RandomDirection();
    return {Curve2d{Line2d{Point2{Coordinate(), Coordinate()}, d}},
            {{{}, PlaneCurveComponent::OriginX},
             {{}, PlaneCurveComponent::OriginY},
             {{}, PlaneCurveComponent::DirectionX},
             {{}, PlaneCurveComponent::DirectionY}}};
}

RandomSide RandomCircle()
{
    const Vec2 x = RandomDirection();
    return {Curve2d{Circle2d{Point2{Coordinate() / 2.0, Coordinate() / 2.0},
                             *rc::gen::inRange(10, 31) / 10.0, x, Vec2{-x.Y, x.X}}},
            {{{}, PlaneCurveComponent::CenterX},
             {{}, PlaneCurveComponent::CenterY},
             {{}, PlaneCurveComponent::Radius}}};
}

RandomSide RandomEllipse()
{
    const Vec2 x = RandomDirection();
    const double a = *rc::gen::inRange(15, 31) / 10.0;
    return {Curve2d{Ellipse2d{Point2{Coordinate() / 2.0, Coordinate() / 2.0}, a,
                              a * *rc::gen::inRange(4, 10) / 10.0, x, Vec2{-x.Y, x.X}}},
            {{{}, PlaneCurveComponent::CenterX},
             {{}, PlaneCurveComponent::CenterY},
             {{}, PlaneCurveComponent::MajorRadius},
             {{}, PlaneCurveComponent::MinorRadius}}};
}

/// A cubic Bezier across the origin, rational half the time, two of its
/// poles bound in x, y and, when rational, the weight.
RandomSide RandomSpline()
{
    const bool rational = *rc::gen::arbitrary<bool>();
    std::vector<Point2> poles;
    std::vector<double> weights;
    for (int i = 0; i < 4; ++i) {
        poles.push_back({-3.0 + 2.0 * i + Coordinate() / 4.0, Coordinate()});
        if (rational) {
            weights.push_back(*rc::gen::inRange(5, 21) / 10.0);
        }
    }
    const auto made = Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
                                                std::move(poles), std::move(weights));
    RC_ASSERT(made.has_value());
    RandomSide side{Curve2d{*made}, {}};
    for (const std::size_t pole : {std::size_t{1}, std::size_t{2}}) {
        side.Shape.push_back({{}, PlaneCurveComponent::PoleX, pole});
        side.Shape.push_back({{}, PlaneCurveComponent::PoleY, pole});
        if (rational) {
            side.Shape.push_back({{}, PlaneCurveComponent::PoleWeight, pole});
        }
    }
    return side;
}

RandomSide RandomCurve()
{
    switch (*rc::gen::inRange(0, 4)) {
    case 0:
        return RandomLine();
    case 1:
        return RandomCircle();
    case 2:
        return RandomEllipse();
    default:
        return RandomSpline();
    }
}

/// The curve with each binding's value moved by d(k), k counting from first.
template <typename Step>
Curve2d Moved(const RandomSide& side, const Step& d, const std::size_t first)
{
    using C = PlaneCurveComponent;
    const auto delta = [&](const C component, const std::size_t pole = 0) {
        double sum = 0.0;
        for (std::size_t k = 0; k < side.Shape.size(); ++k) {
            if (side.Shape[k].Component == component && side.Shape[k].Pole == pole) {
                sum += d(first + k);
            }
        }
        return sum;
    };
    if (const Line2d* line = side.Curve.Target<Line2d>()) {
        Line2d moved = *line;
        moved.Origin = moved.Origin + Vec2{delta(C::OriginX), delta(C::OriginY)};
        moved.Direction = moved.Direction + Vec2{delta(C::DirectionX), delta(C::DirectionY)};
        return Curve2d{moved};
    }
    if (const Circle2d* circle = side.Curve.Target<Circle2d>()) {
        Circle2d moved = *circle;
        moved.Center = moved.Center + Vec2{delta(C::CenterX), delta(C::CenterY)};
        moved.Radius += delta(C::Radius);
        return Curve2d{moved};
    }
    if (const Ellipse2d* ellipse = side.Curve.Target<Ellipse2d>()) {
        Ellipse2d moved = *ellipse;
        moved.Center = moved.Center + Vec2{delta(C::CenterX), delta(C::CenterY)};
        moved.MajorRadius += delta(C::MajorRadius);
        moved.MinorRadius += delta(C::MinorRadius);
        return Curve2d{moved};
    }
    BSplineCurve2d moved = *side.Curve.Target<BSplineCurve2d>();
    for (std::size_t i = 0; i < moved.Poles.size(); ++i) {
        moved.Poles[i] = moved.Poles[i] + Vec2{delta(C::PoleX, i), delta(C::PoleY, i)};
        if (!moved.Weights.empty()) {
            moved.Weights[i] += delta(C::PoleWeight, i);
        }
    }
    return Curve2d{moved};
}

/// The points as the products lay them out: x, y, s, t.
std::vector<double> Points(const Curve2d& a, const Curve2d& b)
{
    std::vector<double> values;
    for (const PlanarPoint& p : Meet(a, b).Points) {
        values.insert(values.end(), {p.Position.X, p.Position.Y, p.ParameterA, p.ParameterB});
    }
    return values;
}

double Norm(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

/// Every point a decided crossing well away from a tangency.
bool WellConditioned(const Curve2d& a, const Curve2d& b)
{
    const PlanarIntersection at = Meet(a, b);
    if (at.Points.empty() || !at.Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(at.Points, [&](const PlanarPoint& p) {
        const auto lin = Stereon::Linearize(a, b, p);
        return lin && lin->Status == DerivativeStatus::Valid && lin->Conditioning >= 0.25;
    });
}

/// Second-order Taylor decay and <J v, w> = <v, J^T w> for one pair. A step
/// that changes the number of points is a new regime and discarded; points
/// are paired by position, since a circle's angle may pass its seam; and
/// parameters are compared modulo a turn, which leaves a line's and a
/// B-spline's small differences as they are.
void RequireProducts(RandomSide a, RandomSide b)
{
    const std::size_t count = a.Shape.size() + b.Shape.size();
    ParameterSet set;
    std::vector<double> references;
    std::vector<double> seed;
    for (std::size_t k = 0; k < count; ++k) {
        const PlaneCurveBinding& shape =
            k < a.Shape.size() ? a.Shape[k] : b.Shape[k - a.Shape.size()];
        const bool weight = shape.Component == PlaneCurveComponent::PoleWeight;
        references.push_back(*rc::gen::inRange(1, 20) / 10.0);
        seed.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        const ParameterId id =
            *set.Add("p" + std::to_string(k),
                     weight ? ParameterQuantity::Dimensionless : ParameterQuantity::Length,
                     references.back());
        (k < a.Shape.size() ? a.Shape[k] : b.Shape[k - a.Shape.size()]).Parameter = id;
    }
    SensitivityContext sensitivity(std::move(set));
    RC_ASSERT(sensitivity.SetSeed(seed).has_value());
    const PlanarIntersection at = Meet(a.Curve, b.Curve);
    const auto jvp = Stereon::IntersectJvp(a.Curve, a.Shape, b.Curve, b.Shape, at, sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const std::vector<double> jv(jvp->GetValues()->begin(), jvp->GetValues()->end());
    const std::vector<double> base = Points(a.Curve, b.Curve);
    const auto moved = [&](const double h) {
        const auto d = [&](const std::size_t k) { return h * seed[k] * references[k]; };
        return Points(Moved(a, d, 0), Moved(b, d, a.Shape.size()));
    };
    for (const double h : {1e-2, 1e-3}) {
        RC_PRE(moved(h).size() == base.size());
    }
    const std::size_t points = base.size() / Values;
    const auto residual = [&](const double h) {
        const std::vector<double> after = moved(h);
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < points; ++i) {
            const auto gap = [&](const std::size_t j) {
                return std::hypot(after[Values * j] - base[Values * i],
                                  after[Values * j + 1] - base[Values * i + 1]);
            };
            std::size_t nearest = 0;
            for (std::size_t j = 1; j < points; ++j) {
                nearest = gap(j) < gap(nearest) ? j : nearest;
            }
            for (std::size_t k = 0; k < Values; ++k) {
                const double difference = after[Values * nearest + k] - base[Values * i + k];
                r[Values * i + k] = (k >= 2 ? std::remainder(difference, tau) : difference) -
                                    h * jv[Values * i + k];
            }
        }
        return Norm(r);
    };
    const double coarse = residual(1e-2);
    const double fine = residual(1e-3);
    RC_ASSERT(fine <= coarse / 30.0 + 1e-10 * (1.0 + Norm(base)));

    std::vector<double> w;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    RC_ASSERT(sensitivity.SetAdjoint(w).has_value());
    const auto vjp = Stereon::IntersectVjp(a.Curve, a.Shape, b.Curve, b.Shape, at, sensitivity);
    RC_ASSERT(vjp.has_value());
    RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
    const auto jtw = *vjp->GetValues();
    double left = 0.0;
    double right = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        left += jv[i] * w[i];
        scale += std::abs(jv[i] * w[i]);
    }
    for (std::size_t k = 0; k < seed.size(); ++k) {
        right += seed[k] * jtw[k];
        scale += std::abs(seed[k] * jtw[k]);
    }
    RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
}

} // namespace

TEST_CASE("plane crossings move as their closed forms say", "[intersect][linearization][planar]")
{
    SECTION("the x-axis through a circle of radius 2")
    {
        // The line meets the circle at s = -2 (angle pi), then s = 2 (angle
        // 0). A growing radius moves each point outward by one, s with it,
        // and leaves the angles. Raising the circle by h moves the points
        // along the line by -h / x, a second-order zero at h = 0, and the
        // angles by -1 / (r cos t): +1/2 at pi, -1/2 at 0.
        const Curve2d circle{Circle2d{Point2{}, 2.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
        const std::vector<double> radius = RateOf(circle, PlaneCurveComponent::Radius, XAxis);
        REQUIRE(radius.size() == 2 * Values);
        // A is the circle here: the points come in order of its angle.
        CHECK(std::abs(radius[0] - 1.0) <= 1e-15); // angle 0, x = 2
        CHECK(std::abs(radius[1]) <= 1e-15);
        CHECK(std::abs(radius[2]) <= 1e-15);
        CHECK(std::abs(radius[3] - 1.0) <= 1e-15);
        CHECK(std::abs(radius[Values] + 1.0) <= 1e-15); // angle pi, x = -2
        CHECK(std::abs(radius[Values + 2]) <= 1e-15);
        CHECK(std::abs(radius[Values + 3] + 1.0) <= 1e-15);
        const std::vector<double> rising = RateOf(circle, PlaneCurveComponent::CenterY, XAxis);
        CHECK(std::abs(rising[0]) <= 1e-15);
        CHECK(std::abs(rising[1]) <= 1e-15);
        CHECK(std::abs(rising[2] + 0.5) <= 1e-15);
        CHECK(std::abs(rising[Values + 2] - 0.5) <= 1e-15);
    }
    SECTION("the x-axis through an ellipse of semi-axes 3 and 2")
    {
        const Curve2d ellipse{Ellipse2d{Point2{}, 3.0, 2.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
        const std::vector<double> major = RateOf(ellipse, PlaneCurveComponent::MajorRadius, XAxis);
        REQUIRE(major.size() == 2 * Values);
        CHECK(std::abs(major[0] - 1.0) <= 1e-15);
        CHECK(std::abs(major[2]) <= 1e-15);
        CHECK(std::abs(major[Values] + 1.0) <= 1e-15);
        const std::vector<double> minor = RateOf(ellipse, PlaneCurveComponent::MinorRadius, XAxis);
        CHECK(std::abs(minor[0]) <= 1e-15);
        CHECK(std::abs(minor[2]) <= 1e-15);
    }
    SECTION("the x-axis through a segment from (0, -1) to (0, 1)")
    {
        // The degree-1 B-spline meets the axis at u = 1/2. Moving its upper
        // pole along x by one moves the point there by u = 1/2, and the
        // line's s with it; u stays.
        const auto segment = Stereon::MakeBSplineCurve(1, {0.0, 0.0, 1.0, 1.0},
                                                       {Point2{0.0, -1.0}, Point2{0.0, 1.0}});
        REQUIRE(segment.has_value());
        const std::vector<double> slant =
            RateOf(Curve2d{*segment}, PlaneCurveComponent::PoleX, XAxis, 1);
        REQUIRE(slant.size() == Values);
        CHECK(std::abs(slant[0] - 0.5) <= 1e-15);
        CHECK(std::abs(slant[1]) <= 1e-15);
        CHECK(std::abs(slant[2]) <= 1e-15);
        CHECK(std::abs(slant[3] - 0.5) <= 1e-15);
    }
}

TEST_CASE("plane crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][planar][property]")
{
    rc::prop("for random pairs of lines, circles, ellipses and cubic B-splines", [] {
        RandomSide a = RandomCurve();
        RandomSide b = RandomCurve();
        RC_PRE(WellConditioned(a.Curve, b.Curve));
        RequireProducts(std::move(a), std::move(b));
    });
}

TEST_CASE("a plane tangency, an overlap and unbound curves carry no derivative",
          "[intersect][linearization][planar]")
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<PlaneCurveBinding> radius{{p, PlaneCurveComponent::Radius}};
    // The circle of radius 1 about (0, 1) touches the x-axis at the origin.
    const Curve2d touching{Circle2d{Point2{0.0, 1.0}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    const auto touch = Stereon::Intersect(touching, XAxis, ctx);
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1u);
    REQUIRE(touch->Points[0].Contact == ContactType::Tangential);
    const auto tangent = Stereon::IntersectJvp(touching, radius, XAxis, {}, *touch, sensitivity);
    REQUIRE(tangent.has_value());
    CHECK(tangent->GetStatus() == DerivativeStatus::Nondifferentiable);
    CHECK_FALSE(tangent->GetValues().has_value());
    // A circle against itself: an overlap.
    const auto same = Stereon::Intersect(touching, touching, ctx);
    REQUIRE(same.has_value());
    REQUIRE_FALSE(same->Overlaps.empty());
    const auto overlap = Stereon::IntersectJvp(touching, radius, touching, {}, *same, sensitivity);
    REQUIRE(overlap.has_value());
    CHECK(overlap->GetStatus() == DerivativeStatus::Nondifferentiable);
    // A parabola binds nothing yet.
    const Curve2d parabola{
        Stereon::Parabola2d{Point2{0.0, -1.0}, 1.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}}};
    const auto unsupported =
        Stereon::IntersectJvp(parabola, {}, XAxis, {}, PlanarIntersection{}, sensitivity);
    REQUIRE(unsupported.has_value());
    CHECK(unsupported->GetStatus() == DerivativeStatus::NotSupported);
}

TEST_CASE("plane curve bindings are refused where the curve has no such value",
          "[intersect][linearization][planar]")
{
    ParameterSet set;
    const ParameterId length = *set.Add("p", ParameterQuantity::Length, 1.0);
    const ParameterId weight = *set.Add("w", ParameterQuantity::Dimensionless, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    const Curve2d circle{Circle2d{Point2{}, 2.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    const auto at = Stereon::Intersect(circle, XAxis, ctx);
    REQUIRE(at.has_value());
    const auto code = [&](const std::vector<PlaneCurveBinding>& bindings) {
        return Stereon::IntersectJvp(circle, bindings, XAxis, {}, *at, sensitivity)
            .error()
            .GetCode();
    };
    // A circle has no direction, no major radius, no poles; a radius is a
    // length, never dimensionless.
    for (const PlaneCurveComponent component :
         {PlaneCurveComponent::DirectionX, PlaneCurveComponent::MajorRadius,
          PlaneCurveComponent::PoleX}) {
        CHECK(code({{length, component}}) == ErrorCode::InvalidInput);
    }
    CHECK(code({{weight, PlaneCurveComponent::Radius}}) == ErrorCode::InvalidInput);
    // A stale result: the circle grown since.
    const Curve2d grown{Circle2d{Point2{}, 2.5, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    CHECK(Stereon::IntersectJvp(grown, {}, XAxis, {}, *at, sensitivity).error().GetCode() ==
          ErrorCode::InvalidInput);
}
