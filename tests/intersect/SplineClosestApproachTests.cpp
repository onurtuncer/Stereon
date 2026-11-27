// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ClosestApproach with polynomial B-spline curves (ADR-0052, addendum of
// 2026-11-12, decision 5): the minimum proved global by subdivision, or
// Undecided with a certified lower bound. Closed forms: a parabola written
// as a quadratic Bezier against a line, a mirrored parabola, a circle (an
// end minimum, and a tie between both ends), a corner, a restricted domain.
// An independent oracle (dense sampling, then projected Newton's method on
// the curves' own evaluation) on random cubic B-splines against lines,
// circles and other B-splines: the lower bound never exceeds the oracle's
// minimum, and a Unique answer attains it. Rational B-splines answer
// NotSupported, and so do the products of any B-spline's approach.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::ApproachEnd;
using Stereon::ApproachUniqueness;
using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::ClosestApproachResult;
using Stereon::Curve;
using Stereon::DerivativeStatus;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Vec3;

namespace {

/// A B-spline curve, required valid.
BSplineCurve Spline(const int degree, std::vector<double> knots, std::vector<Point3> poles,
                    std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

const std::vector<double> quadratic{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

/// y = x^2 in the plane z = 0, x = 2 s - 1 over s in [0, 1]: the quadratic
/// Bezier with poles (-1, 1), (0, -1), (1, 1).
BSplineCurve Parabola()
{
    return Spline(2, quadratic, {{-1, 1, 0}, {0, -1, 0}, {1, 1, 0}});
}

/// y = -x^2 - 2, x = 2 t - 1.
BSplineCurve Mirrored()
{
    return Spline(2, quadratic, {{-1, -3, 0}, {0, -1, 0}, {1, -3, 0}});
}

/// y = |x| over [-1, 1], cornered at the origin, t = 1.
BSplineCurve Vee()
{
    return Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 1, 0}, {0, 0, 0}, {1, 1, 0}});
}

/// The line x = 0, y = -1, along z.
const Curve Below{Line{Point3{0.0, -1.0, 0.0}, Vec3{0.0, 0.0, 1.0}}};

/// ClosestApproach, required to succeed, with its points on the curves and its
/// gap and bounds consistent.
ClosestApproachResult Approach(const Curve& a, const Interval& da, const Curve& b,
                               const Interval& db)
{
    const auto r = Stereon::ClosestApproach(a, da, b, db, ctx);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    // The reported points are on the curves, and the gap is theirs.
    CHECK(Distance(r->OnA, Eval(a, r->ParameterA)) <= 1e-14 * (1.0 + Length(r->OnA - Point3{})));
    CHECK(Distance(r->OnB, Eval(b, r->ParameterB)) <= 1e-14 * (1.0 + Length(r->OnB - Point3{})));
    CHECK(std::abs(r->Distance - Length(r->Gap)) <= 1e-15 * (1.0 + r->Distance));
    CHECK(r->LowerBound <= r->Distance);
    return *r;
}

const Interval Whole = Interval::Everything();

// --- The oracle -------------------------------------------------------------

/// A curve over a bounded range, evaluated by the kernel's own Eval and
/// Deriv: the oracle shares the curves' arithmetic, not the search.
struct Bounded {
    Curve Geometry;
    double Lo = 0.0;
    double Hi = 0.0;
};

/// phi = |A(s) - B(t)|^2 / 2 on the oracle's curves.
double Phi(const Bounded& a, const Bounded& b, const double s, const double t)
{
    const Vec3 d = Eval(a.Geometry, s) - Eval(b.Geometry, t);
    return 0.5 * Dot(d, d);
}

/// Projected Newton's method on phi with backtracking, each parameter held
/// in its range; from a grid sample it reaches the nearby minimum.
std::pair<double, double> Polish(const Bounded& a, const Bounded& b, double s, double t)
{
    for (int k = 0; k < 60; ++k) {
        const Vec3 d = Eval(a.Geometry, s) - Eval(b.Geometry, t);
        const Vec3 a1 = Deriv(a.Geometry, s, 1);
        const Vec3 b1 = Deriv(b.Geometry, t, 1);
        const double gs = Dot(d, a1);
        const double gt = -Dot(d, b1);
        const double hss = Dot(a1, a1) + Dot(d, Deriv(a.Geometry, s, 2));
        const double hst = -Dot(a1, b1);
        const double htt = Dot(b1, b1) - Dot(d, Deriv(b.Geometry, t, 2));
        const bool atS = (s <= a.Lo && gs > 0.0) || (s >= a.Hi && gs < 0.0);
        const bool atT = (t <= b.Lo && gt > 0.0) || (t >= b.Hi && gt < 0.0);
        double ds = 0.0;
        double dt = 0.0;
        const double det = hss * htt - hst * hst;
        if (!atS && !atT && hss > 0.0 && det > 0.0) {
            ds = -(htt * gs - hst * gt) / det;
            dt = -(hss * gt - hst * gs) / det;
        } else {
            ds = atS ? 0.0 : (hss > 0.0 ? -gs / hss : -gs);
            dt = atT ? 0.0 : (htt > 0.0 ? -gt / htt : -gt);
        }
        const double before = Phi(a, b, s, t);
        double step = 1.0;
        for (int halving = 0; halving < 40; ++halving, step *= 0.5) {
            const double s1 = std::clamp(s + step * ds, a.Lo, a.Hi);
            const double t1 = std::clamp(t + step * dt, b.Lo, b.Hi);
            if (Phi(a, b, s1, t1) <= before) {
                s = s1;
                t = t1;
                break;
            }
        }
    }
    return {s, t};
}

/// The oracle's minimum distance: a 160 x 160 grid, its best 24 samples
/// polished.
double OracleDistance(const Bounded& a, const Bounded& b)
{
    constexpr int n = 160;
    std::vector<std::pair<double, std::pair<double, double>>> samples;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            const double s = a.Lo + (a.Hi - a.Lo) * i / n;
            const double t = b.Lo + (b.Hi - b.Lo) * j / n;
            samples.push_back({Phi(a, b, s, t), {s, t}});
        }
    }
    std::ranges::sort(samples, {}, [](const auto& x) { return x.first; });
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < 24; ++k) {
        const auto [s, t] = Polish(a, b, samples[k].second.first, samples[k].second.second);
        best = std::min(best, Phi(a, b, s, t));
    }
    return std::sqrt(2.0 * best);
}

/// A random coordinate in [-2, 2], a multiple of 1/8.
double Coordinate()
{
    return *rc::gen::inRange(-16, 17) / 8.0;
}

/// A random point with coordinates from Coordinate().
Point3 RandomPoint()
{
    return Point3{Coordinate(), Coordinate(), Coordinate()};
}

/// A clamped cubic with 4 to 7 random poles, uniform interior knots, one of
/// them doubled (a C1 joint) at random.
BSplineCurve RandomSpline()
{
    const auto count = static_cast<std::size_t>(*rc::gen::inRange(4, 8));
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < count; ++i) {
        poles.push_back(RandomPoint());
    }
    std::vector<double> knots{0.0, 0.0, 0.0, 0.0};
    const std::size_t inner = count - 4;
    std::vector<double> interior;
    for (std::size_t k = 1; k <= inner; ++k) {
        interior.push_back(static_cast<double>(k) / static_cast<double>(inner + 1));
    }
    if (inner >= 2 && *rc::gen::inRange(0, 2) == 1) {
        interior[1] = interior[0]; // a double knot: C1 there
    }
    knots.insert(knots.end(), interior.begin(), interior.end());
    knots.insert(knots.end(), {1.0, 1.0, 1.0, 1.0});
    return Spline(3, std::move(knots), std::move(poles));
}

/// A random partner: a bounded line, a circle or another cubic, with its
/// range.
Bounded RandomPartner()
{
    switch (*rc::gen::inRange(0, 3)) {
    case 0: {
        Vec3 d{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(d) > 0.1);
        const double lo = *rc::gen::inRange(-16, 0) / 8.0;
        return {Curve{Line{RandomPoint(), d}}, lo, lo + *rc::gen::inRange(1, 25) / 8.0};
    }
    case 1: {
        const Vec3 x{1.0, 0.0, 0.0};
        const double phi = *rc::gen::inRange(0, 16) * std::numbers::pi / 16.0;
        const Vec3 y{0.0, std::cos(phi), std::sin(phi)};
        return {Curve{Circle{RandomPoint(), *rc::gen::inRange(2, 17) / 8.0, x, y}}, 0.0,
                2.0 * std::numbers::pi};
    }
    default:
        return {Curve{RandomSpline()}, 0.0, 1.0};
    }
}

} // namespace

/// y = x^2 against a line, both orders, and over a restricted domain.
TEST_CASE("a B-spline's closest approach to a line, in closed form",
          "[intersect][closest-approach][spline]")
{
    // The nearest point of y = x^2 to (0, -1) is its vertex, at distance 1,
    // and the line's nearest point is its own at z = 0.
    for (const bool exchanged : {false, true}) {
        INFO((exchanged ? "line first" : "spline first"));
        const Curve parabola{Parabola()};
        const ClosestApproachResult r = exchanged ? Approach(Below, Whole, parabola, Whole)
                                                  : Approach(parabola, Whole, Below, Whole);
        CHECK(r.Uniqueness == ApproachUniqueness::Unique);
        const double s = exchanged ? r.ParameterB : r.ParameterA;
        const double t = exchanged ? r.ParameterA : r.ParameterB;
        CHECK(std::abs(s - 0.5) <= 1e-10);
        CHECK(std::abs(t) <= 1e-10);
        CHECK(std::abs(r.Distance - 1.0) <= 1e-14);
        CHECK(r.LowerBound <= 1.0);
        CHECK(r.LowerBound >= 1.0 - 1e-6);
        CHECK((exchanged ? r.EndB : r.EndA) == ApproachEnd::Interior);
        CHECK_FALSE(r.Meets);
    }
    // Over s in [0, 1/4], the nearest point is the end (-1/2, 1/4, 0).
    const ClosestApproachResult end =
        Approach(Curve{Parabola()}, Interval{0.0, 0.25}, Below, Whole);
    CHECK(end.Uniqueness == ApproachUniqueness::Unique);
    CHECK(end.ParameterA == 0.25);
    CHECK(end.EndA == ApproachEnd::Upper);
    CHECK(std::abs(end.Distance - std::sqrt(0.25 + 1.5625)) <= 1e-14);
}

/// y = x^2 against y = -x^2 - 2: 2 apart at their vertices.
TEST_CASE("two B-splines' closest approach, in closed form",
          "[intersect][closest-approach][spline]")
{
    // y = x^2 and y = -x^2 - 2 are nearest at their vertices, 2 apart.
    const ClosestApproachResult r = Approach(Curve{Parabola()}, Whole, Curve{Mirrored()}, Whole);
    CHECK(r.Uniqueness == ApproachUniqueness::Unique);
    CHECK(std::abs(r.ParameterA - 0.5) <= 1e-10);
    CHECK(std::abs(r.ParameterB - 0.5) <= 1e-10);
    CHECK(std::abs(r.Distance - 2.0) <= 1e-14);
    CHECK(r.LowerBound >= 2.0 - 1e-6);
}

/// A circle off the axis (an end minimum) and on it (a tie).
TEST_CASE("a B-spline against a circle: an end minimum, and a tie left undecided",
          "[intersect][closest-approach][spline]")
{
    // The unit circle about (0.3, 3, 0) in z = 0: over x in [-1, 1],
    // (x - 0.3)^2 + (x^2 - 3)^2 is least at x = 1 (4.49; 5.69 at x = -1,
    // and the interior critical points are maxima), so the distance is
    // sqrt(4.49) - 1, at the parabola's upper end.
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    const Curve off{Circle{Point3{0.3, 3.0, 0.0}, 1.0, x, y}};
    const ClosestApproachResult r = Approach(Curve{Parabola()}, Whole, off, Whole);
    CHECK(r.Uniqueness == ApproachUniqueness::Unique);
    CHECK(r.ParameterA == 1.0);
    CHECK(r.EndA == ApproachEnd::Upper);
    CHECK(std::abs(r.Distance - (std::sqrt(4.49) - 1.0)) <= 1e-12);
    // Centred on the axis, both ends are equally near: two minimisers the
    // search cannot separate. The best is reported with a lower bound.
    const Curve centred{Circle{Point3{0.0, 3.0, 0.0}, 1.0, x, y}};
    const ClosestApproachResult tie = Approach(Curve{Parabola()}, Whole, centred, Whole);
    CHECK(tie.Uniqueness == ApproachUniqueness::Undecided);
    CHECK(std::abs(tie.Distance - (std::sqrt(5.0) - 1.0)) <= 1e-12);
    CHECK(tie.LowerBound <= std::sqrt(5.0) - 1.0);
}

/// A vee's corner nearest a line: the distance right, no false proof.
TEST_CASE("a B-spline's corner nearest: the right distance, never a false proof",
          "[intersect][closest-approach][spline]")
{
    // The vee's corner is at distance 1 from the line; phi has no
    // derivative there, and a cluster across the corner is not proved.
    const ClosestApproachResult r = Approach(Curve{Vee()}, Whole, Below, Whole);
    CHECK(std::abs(r.Distance - 1.0) <= 1e-9);
    CHECK(r.LowerBound <= 1.0);
    if (r.Uniqueness == ApproachUniqueness::Unique) {
        CHECK(std::abs(r.ParameterA - 1.0) <= 1e-9);
    }
}

/// Domains beside the spline, rational splines and the products.
TEST_CASE("B-spline closest approach input outside what is proved",
          "[intersect][closest-approach][spline]")
{
    // A domain beside the spline's own.
    const auto beside =
        Stereon::ClosestApproach(Curve{Parabola()}, Interval{2.0, 3.0}, Below, Whole, ctx);
    REQUIRE_FALSE(beside.has_value());
    CHECK(beside.error().GetCode() == ErrorCode::InvalidInput);
    // A rational B-spline: its derivatives are not bounded by its poles.
    const Curve rational{
        Spline(2, quadratic, {{-1, 1, 0}, {0, -1, 0}, {1, 1, 0}}, {1.0, 2.0, 1.0})};
    const auto r = Stereon::ClosestApproach(rational, Below, ctx);
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().GetCode() == ErrorCode::NotSupported);
    // Equal weights are a polynomial.
    const Curve even{Spline(2, quadratic, {{-1, 1, 0}, {0, -1, 0}, {1, 1, 0}}, {2.0, 2.0, 2.0})};
    CHECK(Approach(even, Whole, Below, Whole).Uniqueness == ApproachUniqueness::Unique);
    // The products of a B-spline's approach are not supported yet.
    const ClosestApproachResult at = Approach(Curve{Parabola()}, Whole, Below, Whole);
    ParameterSet set;
    static_cast<void>(*set.Add("p", ParameterQuantity::Length, 1.0));
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const auto jvp = Stereon::ClosestApproachJvp(Curve{Parabola()}, {}, Below, {}, at, sensitivity);
    REQUIRE(jvp.has_value());
    CHECK(jvp->GetStatus() == DerivativeStatus::NotSupported);
    CHECK_FALSE(jvp->GetValues().has_value());
}

/// The lower bound and Unique answers against the oracle on random cubics.
TEST_CASE("a B-spline's closest approach agrees with an independent oracle",
          "[intersect][closest-approach][spline][property]")
{
    rc::prop("random cubics against lines, circles and cubics", [] {
        const Bounded a{Curve{RandomSpline()}, 0.0, 1.0};
        const Bounded b = RandomPartner();
        const bool exchanged = *rc::gen::inRange(0, 2) == 1;
        const Bounded& first = exchanged ? b : a;
        const Bounded& second = exchanged ? a : b;
        const auto r =
            Stereon::ClosestApproach(first.Geometry, Interval{first.Lo, first.Hi}, second.Geometry,
                                     Interval{second.Lo, second.Hi}, ctx);
        RC_ASSERT(r.has_value());
        const double oracle = OracleDistance(first, second);
        const double scale = 1e-9 * (1.0 + oracle);
        // A lower bound is never above the minimum, and a proved minimum is
        // the oracle's.
        RC_ASSERT(r->LowerBound <= oracle + scale);
        if (r->Uniqueness == ApproachUniqueness::Unique) {
            RC_ASSERT(std::abs(r->Distance - oracle) <= scale);
        }
        RC_TAG(r->Uniqueness == ApproachUniqueness::Unique ? "unique" : "undecided");
    });
}
