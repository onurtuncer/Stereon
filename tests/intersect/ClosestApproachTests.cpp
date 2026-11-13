// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The closest approach of two space curves (ADR-0052 and its addendum of
// 2026-11-12): closed forms for two skew lines, a line and a circle and a
// minimum at a domain end; an independent oracle (dense sampling, then
// Newton's method in long double) on random lines, circles and ellipses;
// the metamorphic relations of exchange, rigid motion, reversal and
// scaling; second-order Taylor decay and the adjoint identity of the
// products; every status; and the input refused.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::ApproachEnd;
using Stereon::ApproachUniqueness;
using Stereon::Circle;
using Stereon::ClosestApproachResult;
using Stereon::ConicBinding;
using Stereon::ConicComponent;
using Stereon::Curve;
using Stereon::CurveBindings;
using Stereon::DerivativeStatus;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::LineBinding;
using Stereon::LineComponent;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Vec3;

namespace {

constexpr std::size_t Values = Stereon::ClosestApproachValues;
constexpr double tau = 2.0 * std::numbers::pi;
const Interval Whole = Interval::Everything();

ClosestApproachResult Approach(const Curve& a, const Interval& da, const Curve& b,
                               const Interval& db)
{
    const auto r = Stereon::ClosestApproach(a, da, b, db, ctx);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    return *r;
}

/// The forward product of one binding at unit rate, every reference 1.
std::vector<double> Rate(const Curve& a, const Interval& da, const Curve& b, const Interval& db,
                         const bool onA, const std::optional<LineComponent> line,
                         const std::optional<ConicComponent> conic)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    std::vector<LineBinding> lines;
    std::vector<ConicBinding> conics;
    if (line) {
        lines.push_back({p, *line});
    } else {
        conics.push_back({p, *conic});
    }
    const CurveBindings bound{lines, conics};
    const CurveBindings none{};
    const ClosestApproachResult r = Approach(a, da, b, db);
    const auto jvp =
        Stereon::ClosestApproachJvp(a, onA ? bound : none, b, onA ? none : bound, r, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    INFO(jvp->GetReason());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    return {jvp->GetValues()->begin(), jvp->GetValues()->end()};
}

/// The twelve values against A(s)', B(t)', s', t' and |d|' (d' is A' - B').
void RequireRates(const std::vector<double>& v, const Vec3& xa, const Vec3& xb, const double s,
                  const double t, const double distance)
{
    REQUIRE(v.size() == Values);
    const std::array<double, Values> expected{
        xa.X, xa.Y, xa.Z, xb.X, xb.Y, xb.Z, s, t, xa.X - xb.X, xa.Y - xb.Y, xa.Z - xb.Z, distance};
    for (std::size_t i = 0; i < Values; ++i) {
        INFO("value " << i);
        CHECK(std::abs(v[i] - expected[i]) <= 1e-12);
    }
}

// --- The independent oracle ---------------------------------------------------

using Real = long double;

struct P3 {
    Real X, Y, Z;
};

P3 Of(const Point3& p)
{
    return {p.X, p.Y, p.Z};
}
P3 Of(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

/// A curve in long double, its point and first two derivatives.
struct OCurve {
    bool IsLine = true;
    P3 O{}, D{}, X{}, Y{};
    Real Ra = 0, Rb = 0;
    Real Lo = 0, Hi = 0;
    bool Periodic = false;

    [[nodiscard]] std::array<P3, 3> At(const Real u) const
    {
        if (IsLine) {
            return {P3{O.X + u * D.X, O.Y + u * D.Y, O.Z + u * D.Z}, D, P3{0, 0, 0}};
        }
        const Real c = std::cos(u);
        const Real s = std::sin(u);
        const auto mix = [&](const Real p, const Real q) {
            return P3{p * X.X + q * Y.X, p * X.Y + q * Y.Y, p * X.Z + q * Y.Z};
        };
        const P3 off = mix(Ra * c, Rb * s);
        return {P3{O.X + off.X, O.Y + off.Y, O.Z + off.Z}, mix(-Ra * s, Rb * c),
                mix(-Ra * c, -Rb * s)};
    }
    [[nodiscard]] Real Clamp(const Real u) const { return Periodic ? u : std::clamp(u, Lo, Hi); }
};

Real Dot(const P3& a, const P3& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}
P3 Minus(const P3& a, const P3& b)
{
    return {a.X - b.X, a.Y - b.Y, a.Z - b.Z};
}

OCurve OracleOf(const Curve& c, const Interval& domain)
{
    OCurve o;
    if (const Line* line = c.Target<Line>()) {
        o.O = Of(line->Origin);
        o.D = Of(line->Direction);
        o.Lo = domain.GetLower();
        o.Hi = domain.GetUpper();
        return o;
    }
    o.IsLine = false;
    if (const Circle* circle = c.Target<Circle>()) {
        o.O = Of(circle->Center);
        o.X = Of(circle->XAxis);
        o.Y = Of(circle->YAxis);
        o.Ra = o.Rb = circle->Radius;
    } else {
        const Ellipse& e = *c.Target<Ellipse>();
        o.O = Of(e.Center);
        o.X = Of(e.XAxis);
        o.Y = Of(e.YAxis);
        o.Ra = e.MajorRadius;
        o.Rb = e.MinorRadius;
    }
    o.Periodic = !domain.IsBounded() || domain.GetUpper() - domain.GetLower() >= tau;
    o.Lo = o.Periodic ? 0.0L : static_cast<Real>(domain.GetLower());
    o.Hi = o.Periodic ? static_cast<Real>(tau) : static_cast<Real>(domain.GetUpper());
    return o;
}

/// The nearest point of a line to a point: its parameter, held in domain.
Real Foot(const OCurve& line, const P3& p)
{
    return line.Clamp(Dot(Minus(p, line.O), line.D) / Dot(line.D, line.D));
}

Real Phi(const OCurve& a, const OCurve& b, const Real s, const Real t)
{
    const P3 d = Minus(a.At(s)[0], b.At(t)[0]);
    return Dot(d, d) / 2;
}

/// Projected Newton's method in long double on phi: a parameter resting
/// on an end its gradient pushes against stays there while the other takes
/// its own Newton step.
std::pair<Real, Real> Polish(const OCurve& a, const OCurve& b, Real s, Real t)
{
    const auto pinned = [](const OCurve& c, const Real u, const Real g) {
        return !c.Periodic && ((u <= c.Lo && g > 0) || (u >= c.Hi && g < 0));
    };
    const auto newton = [](const Real h, const Real g) {
        return h > 0 ? -g / h : -g / (std::abs(h) + 1);
    };
    for (int k = 0; k < 100; ++k) {
        const auto pa = a.At(s);
        const auto pb = b.At(t);
        const P3 d = Minus(pa[0], pb[0]);
        const Real gs = Dot(d, pa[1]);
        const Real gt = -Dot(d, pb[1]);
        const Real hss = Dot(pa[1], pa[1]) + Dot(d, pa[2]);
        const Real hst = -Dot(pa[1], pb[1]);
        const Real htt = Dot(pb[1], pb[1]) - Dot(d, pb[2]);
        const bool fixS = pinned(a, s, gs);
        const bool fixT = pinned(b, t, gt);
        Real ds = 0;
        Real dt = 0;
        const Real det = hss * htt - hst * hst;
        if (!fixS && !fixT && hss > 0 && det > 0) {
            ds = (-htt * gs + hst * gt) / det;
            dt = (hst * gs - hss * gt) / det;
        } else {
            ds = fixS ? 0 : newton(hss, gs);
            dt = fixT ? 0 : newton(htt, gt);
        }
        Real step = 1;
        const Real now = Phi(a, b, s, t);
        for (int h = 0; h < 60; ++h, step /= 2) {
            if (Phi(a, b, a.Clamp(s + step * ds), b.Clamp(t + step * dt)) <= now) {
                break;
            }
        }
        s = a.Clamp(s + step * ds);
        t = b.Clamp(t + step * dt);
    }
    return {s, t};
}

using Sample = std::pair<Real, std::pair<Real, Real>>;

/// Two lines, a convex quadratic: feet alternated from both ends and zero.
std::vector<Sample> LineSamples(const OCurve& a, const OCurve& b)
{
    std::vector<Sample> samples;
    for (const Real s0 : {a.Lo, a.Hi, Real{0}}) {
        Real s = std::isfinite(s0) ? a.Clamp(s0) : Real{0};
        Real t = 0;
        for (int k = 0; k < 400; ++k) {
            t = Foot(b, a.At(s)[0]);
            s = Foot(a, b.At(t)[0]);
        }
        samples.push_back({Phi(a, b, s, t), {s, t}});
    }
    return samples;
}

/// A line and a conic: the conic sampled finely, the line at each foot.
std::vector<Sample> FootSamples(const OCurve& a, const OCurve& b)
{
    constexpr int n = 5760;
    std::vector<Sample> samples;
    const bool lineFirst = a.IsLine;
    const OCurve& line = lineFirst ? a : b;
    const OCurve& conic = lineFirst ? b : a;
    for (int i = 0; i <= n; ++i) {
        const Real u = conic.Lo + (conic.Hi - conic.Lo) * static_cast<Real>(i) / n;
        const Real s = Foot(line, conic.At(u)[0]);
        const auto st = lineFirst ? std::pair{s, u} : std::pair{u, s};
        samples.push_back({Phi(a, b, st.first, st.second), st});
    }
    return samples;
}

/// Two conics: a grid over both ranges.
std::vector<Sample> GridSamples(const OCurve& a, const OCurve& b)
{
    constexpr int n = 720;
    const auto at = [](const OCurve& c, const int i) {
        return c.Lo + (c.Hi - c.Lo) * static_cast<Real>(i) / n;
    };
    std::vector<Sample> samples;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            samples.push_back({Phi(a, b, at(a, i), at(b, j)), {at(a, i), at(b, j)}});
        }
    }
    return samples;
}

/// The oracle's minimum: the best samples polished.
Sample OracleMinimum(const OCurve& a, const OCurve& b)
{
    std::vector<Sample> samples;
    if (a.IsLine && b.IsLine) {
        samples = LineSamples(a, b);
    } else if (a.IsLine || b.IsLine) {
        samples = FootSamples(a, b);
    } else {
        samples = GridSamples(a, b);
    }
    std::ranges::sort(samples, {}, [](const Sample& x) { return x.first; });
    Sample best{std::numeric_limits<Real>::infinity(), {}};
    for (std::size_t k = 0; k < std::min<std::size_t>(16, samples.size()); ++k) {
        const auto [s, t] = Polish(a, b, samples[k].second.first, samples[k].second.second);
        const Real value = Phi(a, b, s, t);
        if (value < best.first) {
            best = {value, {s, t}};
        }
    }
    return best;
}

// --- Random geometry ---------------------------------------------------------

double Coordinate()
{
    return *rc::gen::inRange(-16, 17) / 8.0;
}

Vec3 RandomUnit()
{
    const double theta = *rc::gen::inRange(0, 64) * tau / 64.0;
    const double z = *rc::gen::inRange(-15, 16) / 16.0;
    const double r = std::sqrt(1.0 - z * z);
    return {r * std::cos(theta), r * std::sin(theta), z};
}

std::pair<Vec3, Vec3> RandomFrame()
{
    const Vec3 x = RandomUnit();
    Vec3 helper = std::abs(x.Z) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
    Vec3 y = Cross(helper, x);
    y = y / Length(y);
    const double phi = *rc::gen::inRange(0, 32) * tau / 32.0;
    const Vec3 z = Cross(x, y);
    return {x, std::cos(phi) * y + std::sin(phi) * z};
}

/// A random line, circle or ellipse, and a domain for it: the whole curve
/// half the time, a segment or an arc otherwise.
std::pair<Curve, Interval> RandomCurve()
{
    const bool whole = *rc::gen::arbitrary<bool>();
    const Point3 at{Coordinate(), Coordinate(), Coordinate()};
    switch (*rc::gen::inRange(0, 3)) {
    case 0: {
        const double lo = Coordinate();
        return {Curve{Line{at, RandomUnit()}},
                whole ? Whole : Interval{lo, lo + 0.5 + *rc::gen::inRange(1, 17) / 4.0}};
    }
    case 1: {
        const auto [x, y] = RandomFrame();
        const double lo = *rc::gen::inRange(0, 32) * tau / 32.0;
        return {Curve{Circle{at, *rc::gen::inRange(4, 17) / 8.0, x, y}},
                whole ? Interval{0.0, tau} : Interval{lo, lo + *rc::gen::inRange(4, 40) / 8.0}};
    }
    default: {
        const auto [x, y] = RandomFrame();
        const double a = *rc::gen::inRange(6, 17) / 8.0;
        const double lo = *rc::gen::inRange(0, 32) * tau / 32.0;
        return {Curve{Ellipse{at, a, a * *rc::gen::inRange(4, 10) / 10.0, x, y}},
                whole ? Interval{0.0, tau} : Interval{lo, lo + *rc::gen::inRange(4, 40) / 8.0}};
    }
    }
}

double Scale(const ClosestApproachResult& r)
{
    return 1.0 + std::max({std::abs(r.OnA.X), std::abs(r.OnA.Y), std::abs(r.OnA.Z),
                           std::abs(r.OnB.X), std::abs(r.OnB.Y), std::abs(r.OnB.Z)});
}

/// Parameters compared on their curve: an angle modulo a turn.
double ParameterGap(const Curve& c, const double x, const double y)
{
    return c.Target<Line>() != nullptr ? std::abs(x - y) : std::abs(std::remainder(x - y, tau));
}

} // namespace

// --- Closed forms --------------------------------------------------------------

TEST_CASE("the closest approach of two skew lines is its closed form",
          "[intersect][closest-approach]")
{
    // A = (s, 0, 0), B = (0, t, 2): s = t = 0, d = (0, 0, -2). Raising B
    // lengthens the gap by one; sliding B along x moves s with it; tilting
    // A's direction up by e gives s = 2 e / (1 + e^2), so s' = 2, and A's
    // point moves along x by 2 while the gap keeps its length.
    const Curve a{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve b{Line{Point3{0.0, 0.0, 2.0}, Vec3{0.0, 1.0, 0.0}}};
    const ClosestApproachResult r = Approach(a, Whole, b, Whole);
    CHECK(r.Uniqueness == ApproachUniqueness::Unique);
    CHECK(std::abs(r.Distance - 2.0) <= 1e-15);
    CHECK(std::abs(r.ParameterA) <= 1e-15);
    CHECK(std::abs(r.ParameterB) <= 1e-15);
    CHECK(std::abs(r.Gap.Z + 2.0) <= 1e-15);
    CHECK_FALSE(r.Meets);
    RequireRates(Rate(a, Whole, b, Whole, false, LineComponent::OriginZ, {}), Vec3{},
                 Vec3{0.0, 0.0, 1.0}, 0.0, 0.0, 1.0);
    RequireRates(Rate(a, Whole, b, Whole, false, LineComponent::OriginX, {}), Vec3{1.0, 0.0, 0.0},
                 Vec3{1.0, 0.0, 0.0}, 1.0, 0.0, 0.0);
    RequireRates(Rate(a, Whole, b, Whole, true, LineComponent::DirectionZ, {}), Vec3{2.0, 0.0, 0.0},
                 Vec3{}, 2.0, 0.0, 0.0);
}

TEST_CASE("the closest approach of a line and a circle is its closed form",
          "[intersect][closest-approach]")
{
    // The unit circle about the origin in z = 0 and the vertical line
    // through (3, 0, 0): |d|^2 = 10 - 6 cos t, least at t = 0, s = 0,
    // distance 2. A growing radius shortens the gap by one. Moving the line
    // along y by e: |d|^2 = (3 - cos t)^2 + (e - sin t)^2, least at t = e /
    // 3, so t' = 1/3 and the circle's point moves along y by 1/3. Raising
    // the line moves s by -1 and no point.
    const Curve line{Line{Point3{3.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};
    const Curve circle{Circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Interval turn{0.0, tau};
    const ClosestApproachResult r = Approach(line, Whole, circle, turn);
    CHECK(r.Uniqueness == ApproachUniqueness::Unique);
    CHECK(std::abs(r.Distance - 2.0) <= 1e-14);
    CHECK(std::abs(std::remainder(r.ParameterB, tau)) <= 1e-7);
    RequireRates(Rate(line, Whole, circle, turn, false, {}, ConicComponent::Radius), Vec3{},
                 Vec3{1.0, 0.0, 0.0}, 0.0, 0.0, -1.0);
    RequireRates(Rate(line, Whole, circle, turn, true, LineComponent::OriginY, {}),
                 Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 1.0 / 3.0, 0.0}, 0.0, 1.0 / 3.0, 0.0);
    RequireRates(Rate(line, Whole, circle, turn, true, LineComponent::OriginZ, {}), Vec3{}, Vec3{},
                 -1.0, 0.0, 0.0);
}

TEST_CASE("a closest approach at a domain end holds that parameter",
          "[intersect][closest-approach]")
{
    // The skew lines above with A's domain [1, 3]: the unconstrained s = 0
    // lies outside, so s = 1, t = 0, d = (1, 0, -2), distance sqrt 5. Moving
    // A along x moves its end point and the gap's length by 1 / sqrt 5;
    // moving B along y moves t by -1 and no point.
    const Curve a{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve b{Line{Point3{0.0, 0.0, 2.0}, Vec3{0.0, 1.0, 0.0}}};
    const Interval end{1.0, 3.0};
    const ClosestApproachResult r = Approach(a, end, b, Whole);
    CHECK(r.EndA == ApproachEnd::Lower);
    CHECK(r.EndB == ApproachEnd::Interior);
    CHECK(r.ParameterA == 1.0);
    CHECK(std::abs(r.Distance - std::sqrt(5.0)) <= 1e-15);
    RequireRates(Rate(a, end, b, Whole, true, LineComponent::OriginX, {}), Vec3{1.0, 0.0, 0.0},
                 Vec3{}, 0.0, 0.0, 1.0 / std::sqrt(5.0));
    RequireRates(Rate(a, end, b, Whole, false, LineComponent::OriginY, {}), Vec3{}, Vec3{}, 0.0,
                 -1.0, 0.0);
    // An arc whose end is nearest: the quarter arc [pi/2, pi] of the unit
    // circle against the vertical line through (3, 0, 0) rests on t = pi/2,
    // (0, 1, 0), s = 0, distance sqrt 10. Growing the radius moves the point
    // along y by one, away from the line: the gap grows by 1 / sqrt 10.
    const Curve line{Line{Point3{3.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};
    const Curve circle{Circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Interval arc{std::numbers::pi / 2.0, std::numbers::pi};
    const ClosestApproachResult onArc = Approach(line, Whole, circle, arc);
    CHECK(onArc.EndB == ApproachEnd::Lower);
    CHECK(std::abs(onArc.Distance - std::sqrt(10.0)) <= 1e-14);
    RequireRates(Rate(line, Whole, circle, arc, false, {}, ConicComponent::Radius), Vec3{},
                 Vec3{std::cos(std::numbers::pi / 2.0), 1.0, 0.0}, 0.0, 0.0, 1.0 / std::sqrt(10.0));
}

// --- Statuses and refusals ---------------------------------------------------

namespace {

DerivativeStatus StatusOf(const Curve& a, const Interval& da, const Curve& b, const Interval& db)
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const auto r = Stereon::ClosestApproach(a, da, b, db, ctx);
    REQUIRE(r.has_value());
    const auto jvp = Stereon::ClosestApproachJvp(a, {}, b, {}, *r, sensitivity);
    REQUIRE(jvp.has_value());
    REQUIRE(sensitivity.SetAdjoint(std::vector<double>(Values, 1.0)).has_value());
    const auto vjp = Stereon::ClosestApproachVjp(a, {}, b, {}, *r, sensitivity);
    REQUIRE(vjp.has_value());
    CHECK(vjp->GetStatus() == jvp->GetStatus());
    CHECK(jvp->GetValues().has_value() == Stereon::HasDerivative(jvp->GetStatus()));
    static_cast<void>(p);
    return jvp->GetStatus();
}

} // namespace

TEST_CASE("a closest approach's status comes from what the primal proved",
          "[intersect][closest-approach]")
{
    const Curve xAxis{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve circle{Circle{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Interval turn{0.0, tau};
    SECTION("parallel lines: a set of closest pairs")
    {
        const Curve above{Line{Point3{0.0, 0.0, 1.0}, Vec3{2.0, 0.0, 0.0}}};
        CHECK(Approach(xAxis, Whole, above, Whole).Uniqueness == ApproachUniqueness::NotUnique);
        CHECK(StatusOf(xAxis, Whole, above, Whole) == DerivativeStatus::Nondifferentiable);
        // Apart along their direction, segments have one closest pair.
        const ClosestApproachResult apart =
            Approach(xAxis, Interval{0.0, 1.0}, above, Interval{2.0, 3.0});
        CHECK(apart.Uniqueness == ApproachUniqueness::Unique);
        CHECK(std::abs(apart.Distance - std::sqrt(9.0 + 1.0)) <= 1e-14);
    }
    SECTION("a line along a circle's axis, and two coaxial circles: a set")
    {
        const Curve axis{Line{Point3{0.0, 0.0, 5.0}, Vec3{0.0, 0.0, -1.0}}};
        CHECK(Approach(axis, Whole, circle, turn).Uniqueness == ApproachUniqueness::NotUnique);
        CHECK(StatusOf(axis, Whole, circle, turn) == DerivativeStatus::Nondifferentiable);
        const Curve higher{
            Circle{Point3{0.0, 0.0, 2.0}, 3.0, Vec3{0.0, 1.0, 0.0}, Vec3{-1.0, 0.0, 0.0}}};
        CHECK(Approach(circle, turn, higher, turn).Uniqueness == ApproachUniqueness::NotUnique);
    }
    SECTION("curves that meet")
    {
        const Curve yAxis{Line{Point3{}, Vec3{0.0, 1.0, 0.0}}};
        const ClosestApproachResult meet = Approach(xAxis, Whole, yAxis, Whole);
        CHECK(meet.Meets);
        CHECK(StatusOf(xAxis, Whole, yAxis, Whole) == DerivativeStatus::Nondifferentiable);
    }
    SECTION("an end the minimum may leave")
    {
        // A's domain starts where the unconstrained minimum is: s = 0.
        const Curve b{Line{Point3{0.0, 0.0, 2.0}, Vec3{0.0, 1.0, 0.0}}};
        const ClosestApproachResult r = Approach(xAxis, Interval{0.0, 3.0}, b, Whole);
        CHECK(r.EndA == ApproachEnd::Lower);
        CHECK(StatusOf(xAxis, Interval{0.0, 3.0}, b, Whole) == DerivativeStatus::Nondifferentiable);
    }
    SECTION("two equal minima the search cannot separate")
    {
        // The line (0, s, 2) is nearest the unit circle at t = pi/2 and 3pi/2.
        const Curve across{Line{Point3{0.0, 0.0, 2.0}, Vec3{0.0, 1.0, 0.0}}};
        const ClosestApproachResult r = Approach(across, Whole, circle, turn);
        CHECK(r.Uniqueness == ApproachUniqueness::Undecided);
        CHECK(std::abs(r.Distance - 2.0) <= 1e-12);
        CHECK(r.LowerBound <= r.Distance);
        CHECK(StatusOf(across, Whole, circle, turn) == DerivativeStatus::Indeterminate);
    }
    SECTION("nearly parallel skew lines")
    {
        const Curve tilted{Line{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 1e-3, 0.0}}};
        CHECK(Approach(xAxis, Whole, tilted, Whole).Uniqueness == ApproachUniqueness::Unique);
        CHECK(StatusOf(xAxis, Whole, tilted, Whole) == DerivativeStatus::NearEvent);
    }
    SECTION("a B-spline is not supported yet")
    {
        const auto spline = Stereon::MakeBSplineCurve(
            1, {0.0, 0.0, 1.0, 1.0}, {Point3{0.0, 0.0, 1.0}, Point3{1.0, 1.0, 1.0}});
        REQUIRE(spline.has_value());
        const auto r = Stereon::ClosestApproach(Curve{*spline}, xAxis, ctx);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::NotSupported);
        ParameterSet set;
        static_cast<void>(*set.Add("p", ParameterQuantity::Length, 1.0));
        SensitivityContext sensitivity(std::move(set));
        REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
        const auto jvp = Stereon::ClosestApproachJvp(Curve{*spline}, {}, xAxis, {},
                                                     ClosestApproachResult{}, sensitivity);
        REQUIRE(jvp.has_value());
        CHECK(jvp->GetStatus() == DerivativeStatus::NotSupported);
    }
}

TEST_CASE("closest approach input is refused where it is wrong", "[intersect][closest-approach]")
{
    const Curve xAxis{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve b{Line{Point3{0.0, 0.0, 2.0}, Vec3{0.0, 1.0, 0.0}}};
    const Curve ellipse{
        Ellipse{Point3{0.0, 0.0, 3.0}, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const auto code = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        return result.error().GetCode();
    };
    CHECK(code(Stereon::ClosestApproach(Curve{Line{Point3{}, Vec3{}}}, b, ctx)) ==
          ErrorCode::InvalidInput);
    CHECK(code(Stereon::ClosestApproach(xAxis, Interval::Empty(), b, Whole, ctx)) ==
          ErrorCode::InvalidInput);
    CHECK(code(Stereon::ClosestApproach(Curve{}, b, ctx)) == ErrorCode::InvalidInput);

    ParameterSet set;
    const ParameterId length = *set.Add("p", ParameterQuantity::Length, 1.0);
    const ParameterId angle = *set.Add("a", ParameterQuantity::Angle, 1.0);
    SensitivityContext sensitivity(std::move(set));
    const ClosestApproachResult r = Approach(xAxis, Whole, ellipse, Interval{0.0, tau});
    const auto jvp = [&](const CurveBindings& onA, const CurveBindings& onB) {
        return Stereon::ClosestApproachJvp(xAxis, onA, ellipse, onB, r, sensitivity);
    };
    // No seed yet.
    CHECK(code(jvp({}, {})) == ErrorCode::InvalidInput);
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    const std::vector<ConicBinding> conicOnLine{{length, ConicComponent::CenterX}};
    const std::vector<ConicBinding> radius{{length, ConicComponent::Radius}};
    const std::vector<ConicBinding> asAngle{{angle, ConicComponent::MajorRadius}};
    const std::vector<LineBinding> lineOnConic{{length, LineComponent::OriginX}};
    CHECK(code(jvp({{}, conicOnLine}, {})) == ErrorCode::InvalidInput);
    CHECK(code(jvp({}, {{}, radius})) == ErrorCode::InvalidInput);
    CHECK(code(jvp({}, {{}, asAngle})) == ErrorCode::InvalidInput);
    CHECK(code(jvp({}, {lineOnConic, {}})) == ErrorCode::InvalidInput);
    // A stale result: the ellipse moved since.
    const Curve moved{
        Ellipse{Point3{0.0, 0.0, 3.5}, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    CHECK(code(Stereon::ClosestApproachJvp(xAxis, {}, moved, {}, r, sensitivity)) ==
          ErrorCode::InvalidInput);
    // An adjoint of another length.
    REQUIRE(sensitivity.SetAdjoint(std::vector<double>(3, 1.0)).has_value());
    CHECK(code(Stereon::ClosestApproachVjp(xAxis, {}, ellipse, {}, r, sensitivity)) ==
          ErrorCode::InvalidInput);
}

// --- The oracle and the metamorphic relations --------------------------------

namespace {

/// A random pair whose closest approach the primal proved unique and
/// which does not meet.
struct Pair {
    Curve A;
    Interval DomainA;
    Curve B;
    Interval DomainB;
    ClosestApproachResult R;
};

/// The normalised Hessian's smaller eigenvalue at a proved approach, as
/// the products report it, and whether they are Valid.
std::pair<bool, double> Regularity(const Curve& a, const Curve& b, const ClosestApproachResult& r)
{
    ParameterSet set;
    static_cast<void>(*set.Add("p", ParameterQuantity::Length, 1.0));
    SensitivityContext sensitivity(std::move(set));
    RC_ASSERT(sensitivity.SetSeed(std::vector<double>{0.0}).has_value());
    const auto jvp = Stereon::ClosestApproachJvp(a, {}, b, {}, r, sensitivity);
    RC_ASSERT(jvp.has_value());
    return {jvp->GetStatus() == DerivativeStatus::Valid, jvp->GetConditioning()};
}

/// A random pair whose closest approach the primal proved unique, apart and
/// well conditioned, so that its parameters are well defined to rounding.
Pair RandomPair()
{
    auto [a, da] = RandomCurve();
    auto [b, db] = RandomCurve();
    const auto r = Stereon::ClosestApproach(a, da, b, db, ctx);
    RC_ASSERT(r.has_value());
    RC_PRE(r->Uniqueness == ApproachUniqueness::Unique && !r->Meets);
    const auto [valid, conditioning] = Regularity(a, b, *r);
    RC_PRE(valid && conditioning >= 0.01);
    return {a, da, b, db, *r};
}

/// The approaches agree on the distance; whether the transformed one was
/// proved unique too, so that its pair can be compared. The proof's bounds
/// are not invariant under the transformations (a scaled line's parameter
/// scales, its curves' bounds do not), so a borderline proof may fail on
/// one side: that one is then undecided, never wrong, and its distance
/// still brackets the minimum.
bool RequireSame(const ClosestApproachResult& x, const ClosestApproachResult& y,
                 const double distanceScale)
{
    RC_LOG() << "uniqueness " << static_cast<int>(x.Uniqueness) << " -> "
             << static_cast<int>(y.Uniqueness) << ", distance " << x.Distance << " -> "
             << y.Distance << ", parameters (" << x.ParameterA << ", " << x.ParameterB << ") -> ("
             << y.ParameterA << ", " << y.ParameterB << "), ends " << static_cast<int>(x.EndA)
             << static_cast<int>(x.EndB) << " -> " << static_cast<int>(y.EndA)
             << static_cast<int>(y.EndB) << '\n';
    const double tolerance = 1e-12 * (1.0 + y.Distance) * Scale(y);
    if (y.Uniqueness == ApproachUniqueness::Undecided) {
        RC_TAG("transformed not proved unique");
        RC_ASSERT(y.LowerBound <= x.Distance * distanceScale + tolerance);
        RC_ASSERT(x.Distance * distanceScale <= y.Distance + tolerance);
        return false;
    }
    RC_ASSERT(y.Uniqueness == ApproachUniqueness::Unique);
    RC_ASSERT(std::abs(x.Distance * distanceScale - y.Distance) <= tolerance);
    return true;
}

Vec3 Turned(const Vec3& v, const Vec3& axis, const double angle)
{
    return std::cos(angle) * v + std::sin(angle) * Cross(axis, v) +
           ((1.0 - std::cos(angle)) * Dot(axis, v)) * axis;
}

/// The rigid motion x -> R x + shift, R a turn about `axis`.
struct Motion {
    Vec3 Axis;
    double Angle;
    Vec3 Shift;

    [[nodiscard]] Point3 Of(const Point3& p) const
    {
        return Point3{} + Turned(p - Point3{}, Axis, Angle) + Shift;
    }
    [[nodiscard]] Curve Of(const Curve& c) const
    {
        if (const Line* l = c.Target<Line>()) {
            return Curve{Line{Of(l->Origin), Turned(l->Direction, Axis, Angle)}};
        }
        if (const Circle* k = c.Target<Circle>()) {
            return Curve{Circle{Of(k->Center), k->Radius, Turned(k->XAxis, Axis, Angle),
                                Turned(k->YAxis, Axis, Angle)}};
        }
        const Ellipse& e = *c.Target<Ellipse>();
        return Curve{Ellipse{Of(e.Center), e.MajorRadius, e.MinorRadius,
                             Turned(e.XAxis, Axis, Angle), Turned(e.YAxis, Axis, Angle)}};
    }
};

/// The curve traversed the other way, its domain mirrored: the point at
/// -u is the old one at u.
std::pair<Curve, Interval> Reversed(const Curve& c, const Interval& d)
{
    const Interval mirrored = d.IsBounded() ? Interval{-d.GetUpper(), -d.GetLower()} : d;
    if (const Line* l = c.Target<Line>()) {
        return {Curve{Line{l->Origin, -1.0 * l->Direction}}, mirrored};
    }
    if (const Circle* k = c.Target<Circle>()) {
        return {Curve{Circle{k->Center, k->Radius, k->XAxis, -1.0 * k->YAxis}}, mirrored};
    }
    const Ellipse& e = *c.Target<Ellipse>();
    return {Curve{Ellipse{e.Center, e.MajorRadius, e.MinorRadius, e.XAxis, -1.0 * e.YAxis}},
            mirrored};
}

/// The curve scaled by k about the origin; a line's direction kept, so
/// its parameter scales by k too.
std::pair<Curve, Interval> Scaled(const Curve& c, const Interval& d, const double k)
{
    const auto at = [k](const Point3& p) { return Point3{} + k * (p - Point3{}); };
    if (const Line* l = c.Target<Line>()) {
        const Interval scaled = d.IsBounded() ? Interval{k * d.GetLower(), k * d.GetUpper()} : d;
        return {Curve{Line{at(l->Origin), l->Direction}}, scaled};
    }
    if (const Circle* q = c.Target<Circle>()) {
        return {Curve{Circle{at(q->Center), k * q->Radius, q->XAxis, q->YAxis}}, d};
    }
    const Ellipse& e = *c.Target<Ellipse>();
    return {Curve{Ellipse{at(e.Center), k * e.MajorRadius, k * e.MinorRadius, e.XAxis, e.YAxis}},
            d};
}

} // namespace

TEST_CASE("the closest approach agrees with an independent long double oracle",
          "[intersect][closest-approach][property]")
{
    rc::prop("random lines, circles and ellipses, whole or cut", [] {
        auto [a, da] = RandomCurve();
        auto [b, db] = RandomCurve();
        const auto r = Stereon::ClosestApproach(a, da, b, db, ctx);
        RC_ASSERT(r.has_value());
        const auto [phi, st] = OracleMinimum(OracleOf(a, da), OracleOf(b, db));
        const double oracle = static_cast<double>(std::sqrt(2 * phi));
        // The oracle's rounding is long double's; the primal's double's.
        const double tolerance = 64.0 *
                                 (std::numeric_limits<double>::epsilon() +
                                  static_cast<double>(std::numeric_limits<Real>::epsilon())) *
                                 Scale(*r) * Scale(*r);
        RC_TAG(r->Uniqueness == ApproachUniqueness::Unique ? "unique" : "not proved unique");
        // Never worse than the oracle, never below its own lower bound.
        RC_ASSERT(r->Distance <= oracle + tolerance + 1e-9);
        RC_ASSERT(r->LowerBound <= oracle + tolerance);
        if (r->Uniqueness == ApproachUniqueness::Unique) {
            RC_ASSERT(std::abs(r->Distance - oracle) <= tolerance + 1e-12 * oracle);
            // Unique, apart and well conditioned: the oracle's pair is the
            // same, to the square root of rounding at worst.
            if (!r->Meets && Regularity(a, b, *r).second >= 0.01) {
                RC_ASSERT(ParameterGap(a, r->ParameterA, static_cast<double>(st.first)) <=
                          1e-6 * Scale(*r));
                RC_ASSERT(ParameterGap(b, r->ParameterB, static_cast<double>(st.second)) <=
                          1e-6 * Scale(*r));
            }
        }
    });
}

TEST_CASE("the closest approach obeys exchange, rigid motion, reversal and scaling",
          "[intersect][closest-approach][property][metamorphic]")
{
    rc::prop("exchanging the curves exchanges the pair and negates the gap", [] {
        const Pair p = RandomPair();
        const auto x = Stereon::ClosestApproach(p.B, p.DomainB, p.A, p.DomainA, ctx);
        RC_ASSERT(x.has_value());
        if (!RequireSame(p.R, *x, 1.0)) {
            return;
        }
        RC_ASSERT(ParameterGap(p.A, x->ParameterB, p.R.ParameterA) <= 1e-7);
        RC_ASSERT(ParameterGap(p.B, x->ParameterA, p.R.ParameterB) <= 1e-7);
        RC_ASSERT(Length(x->Gap + p.R.Gap) <= 1e-9 * Scale(p.R));
        RC_ASSERT(x->EndA == p.R.EndB && x->EndB == p.R.EndA);
    });
    rc::prop("a rigid motion moves the points and keeps the parameters", [] {
        const Pair p = RandomPair();
        const Motion m{RandomUnit(), *rc::gen::inRange(1, 64) * tau / 64.0,
                       Vec3{Coordinate(), Coordinate(), Coordinate()}};
        const auto x = Stereon::ClosestApproach(m.Of(p.A), p.DomainA, m.Of(p.B), p.DomainB, ctx);
        RC_ASSERT(x.has_value());
        if (!RequireSame(p.R, *x, 1.0)) {
            return;
        }
        RC_ASSERT(ParameterGap(p.A, x->ParameterA, p.R.ParameterA) <= 1e-7);
        RC_ASSERT(ParameterGap(p.B, x->ParameterB, p.R.ParameterB) <= 1e-7);
        RC_ASSERT(Distance(x->OnA, m.Of(p.R.OnA)) <= 1e-7 * Scale(p.R));
    });
    rc::prop("reversing a curve negates its parameter", [] {
        const Pair p = RandomPair();
        const auto [ra, rda] = Reversed(p.A, p.DomainA);
        const auto x = Stereon::ClosestApproach(ra, rda, p.B, p.DomainB, ctx);
        RC_ASSERT(x.has_value());
        if (!RequireSame(p.R, *x, 1.0)) {
            return;
        }
        RC_ASSERT(ParameterGap(p.A, x->ParameterA, -p.R.ParameterA) <= 1e-7);
        RC_ASSERT(Distance(x->OnA, p.R.OnA) <= 1e-7 * Scale(p.R));
    });
    rc::prop("scaling both curves by k scales the distance and a line's parameter", [] {
        const Pair p = RandomPair();
        const double k = *rc::gen::inRange(2, 9) / 4.0;
        const auto [sa, sda] = Scaled(p.A, p.DomainA, k);
        const auto [sb, sdb] = Scaled(p.B, p.DomainB, k);
        const auto x = Stereon::ClosestApproach(sa, sda, sb, sdb, ctx);
        RC_ASSERT(x.has_value());
        if (!RequireSame(p.R, *x, k)) {
            return;
        }
        const double expected = p.A.Target<Line>() != nullptr ? k * p.R.ParameterA : p.R.ParameterA;
        RC_ASSERT(ParameterGap(p.A, x->ParameterA, expected) <= 1e-7 * k * Scale(p.R));
    });
}

// --- Taylor decay and the transpose --------------------------------------------

namespace {

/// Each curve's bindings: every coordinate of a line's origin and direction,
/// a conic's centre and its radius or radii.
std::pair<std::vector<LineBinding>, std::vector<ConicBinding>>
BindingsOf(const Curve& c, std::vector<ParameterId>::const_iterator& next)
{
    std::vector<LineBinding> lines;
    std::vector<ConicBinding> conics;
    if (c.Target<Line>() != nullptr) {
        for (int k = 0; k < 6; ++k) {
            lines.push_back({*next++, static_cast<LineComponent>(k)});
        }
        return {lines, conics};
    }
    for (const ConicComponent k :
         {ConicComponent::CenterX, ConicComponent::CenterY, ConicComponent::CenterZ}) {
        conics.push_back({*next++, k});
    }
    if (c.Target<Circle>() != nullptr) {
        conics.push_back({*next++, ConicComponent::Radius});
    } else {
        conics.push_back({*next++, ConicComponent::MajorRadius});
        conics.push_back({*next++, ConicComponent::MinorRadius});
    }
    return {lines, conics};
}

std::size_t CountOf(const Curve& c)
{
    if (c.Target<Line>() != nullptr) {
        return 6;
    }
    return c.Target<Circle>() != nullptr ? 4 : 5;
}

/// The curve moved by d(first + k) along its k-th binding.
template <typename Step>
Curve MovedCurve(const Curve& c, const Step& d, const std::size_t first)
{
    const Vec3 shift{d(first), d(first + 1), d(first + 2)};
    if (const Line* l = c.Target<Line>()) {
        return Curve{
            Line{l->Origin + shift, l->Direction + Vec3{d(first + 3), d(first + 4), d(first + 5)}}};
    }
    if (const Circle* k = c.Target<Circle>()) {
        return Curve{Circle{k->Center + shift, k->Radius + d(first + 3), k->XAxis, k->YAxis}};
    }
    const Ellipse& e = *c.Target<Ellipse>();
    return Curve{Ellipse{e.Center + shift, e.MajorRadius + d(first + 3),
                         e.MinorRadius + d(first + 4), e.XAxis, e.YAxis}};
}

std::vector<double> Layout(const ClosestApproachResult& r)
{
    return {r.OnA.X,      r.OnA.Y,      r.OnA.Z, r.OnB.X, r.OnB.Y, r.OnB.Z,
            r.ParameterA, r.ParameterB, r.Gap.X, r.Gap.Y, r.Gap.Z, r.Distance};
}

double Norm(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

} // namespace

TEST_CASE("closest approach products: second-order Taylor decay and the transpose",
          "[intersect][closest-approach][property]")
{
    rc::prop("random lines, circles and ellipses with every defining value bound", [] {
        const Pair p = RandomPair();
        const std::size_t count = CountOf(p.A) + CountOf(p.B);
        ParameterSet set;
        std::vector<ParameterId> ids;
        std::vector<double> references;
        std::vector<double> seed;
        for (std::size_t k = 0; k < count; ++k) {
            references.push_back(*rc::gen::inRange(1, 20) / 10.0);
            seed.push_back(*rc::gen::inRange(-10, 11) / 10.0);
            ids.push_back(
                *set.Add("p" + std::to_string(k), ParameterQuantity::Length, references.back()));
        }
        auto next = ids.cbegin();
        const auto [linesA, conicsA] = BindingsOf(p.A, next);
        const auto [linesB, conicsB] = BindingsOf(p.B, next);
        const CurveBindings onA{linesA, conicsA};
        const CurveBindings onB{linesB, conicsB};
        SensitivityContext sensitivity(std::move(set));
        RC_ASSERT(sensitivity.SetSeed(seed).has_value());
        const auto jvp = Stereon::ClosestApproachJvp(p.A, onA, p.B, onB, p.R, sensitivity);
        RC_ASSERT(jvp.has_value());
        // Away from events: a firm Hessian and ends held firmly.
        RC_PRE(jvp->GetStatus() == DerivativeStatus::Valid && jvp->GetConditioning() >= 0.05);
        const std::vector<double> jv(jvp->GetValues()->begin(), jvp->GetValues()->end());
        const std::vector<double> base = Layout(p.R);
        const auto residual = [&](const double h) {
            const auto d = [&](const std::size_t k) { return h * seed[k] * references[k]; };
            const Curve a = MovedCurve(p.A, d, 0);
            const Curve b = MovedCurve(p.B, d, CountOf(p.A));
            const auto moved = Stereon::ClosestApproach(a, p.DomainA, b, p.DomainB, ctx);
            RC_ASSERT(moved.has_value());
            // A step that moves the minimiser onto or off an end, or to
            // another minimum, is a new regime.
            RC_PRE(moved->Uniqueness == ApproachUniqueness::Unique && moved->EndA == p.R.EndA &&
                   moved->EndB == p.R.EndB);
            const std::vector<double> after = Layout(*moved);
            std::vector<double> r(base.size());
            for (std::size_t i = 0; i < base.size(); ++i) {
                double difference = after[i] - base[i];
                if ((i == 6 && p.A.Target<Line>() == nullptr) ||
                    (i == 7 && p.B.Target<Line>() == nullptr)) {
                    difference = std::remainder(difference, tau);
                }
                r[i] = difference - h * jv[i];
            }
            return Norm(r);
        };
        const double coarse = residual(1e-3);
        const double fine = residual(1e-4);
        RC_ASSERT(fine <= coarse / 30.0 + 1e-9 * (1.0 + Norm(base)));

        std::vector<double> w;
        for (std::size_t i = 0; i < Values; ++i) {
            w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        }
        RC_ASSERT(sensitivity.SetAdjoint(w).has_value());
        const auto vjp = Stereon::ClosestApproachVjp(p.A, onA, p.B, onB, p.R, sensitivity);
        RC_ASSERT(vjp.has_value());
        RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
        const auto jtw = *vjp->GetValues();
        double left = 0.0;
        double right = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < Values; ++i) {
            left += jv[i] * w[i];
            scale += std::abs(jv[i] * w[i]);
        }
        for (std::size_t k = 0; k < seed.size(); ++k) {
            right += seed[k] * jtw[k];
            scale += std::abs(seed[k] * jtw[k]);
        }
        RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
    });
}
