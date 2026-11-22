// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Two plane B-splines meeting at a spline's own end or at a corner, a knot
// of full multiplicity (v0.2 ledger, 2D curve-curve). Krawczyk's test
// proves such a root single only for the polynomials continued past the
// end, or for one arm of the corner, and rounding leaves its bracket across
// the end. The handler decides where the root lies against the end
// exactly, on the input doubles, when the end's point (a pole) is a pole
// the other spline passes through, or the other spline is straight there.
// The expected answers are closed forms: the meetings are at poles.
// Metamorphic relations (ADR-0053): exchange, reversal, translation and
// scaling by powers of two, all exact. Before the change every case here
// was reported undecided, the corner touches as crossings, and a spline
// ending a rounding short of a line was reported meeting it.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::BSplineCurve2d;
using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;

namespace {

BSplineCurve2d Spline(const int degree, std::vector<double> knots, std::vector<Point2> poles,
                      std::vector<double> weights = {})
{
    auto made =
        Stereon::MakeBSplineCurve(degree, std::move(knots), std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    return *made;
}

const std::vector<double> cubic{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

/// A cubic from the origin to (endX, 1).
BSplineCurve2d Arch(const double endX = 1.0)
{
    return Spline(3, cubic, {{0, 0}, {0.3, 0.6}, {0.7, 1.2}, {endX, 1}});
}

/// A cubic from (1, 1) to (2, 0), across Arch at their shared pole.
BSplineCurve2d Descent()
{
    return Spline(3, cubic, {{1, 1}, {1.4, 0.7}, {1.7, 0.2}, {2, 0}});
}

/// The vertical segment x = 1 from (1, -1) to (1, 3).
BSplineCurve2d Post()
{
    return Spline(1, {0.0, 0.0, 1.0, 1.0}, {{1, -1}, {1, 3}});
}

/// The corner y = |x| over [-1, 1], a polyline with its corner at t = 1,
/// and a quadratic with a double knot there, cornered at the same pole.
BSplineCurve2d Vee()
{
    return Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, 1}, {0, 0}, {1, 1}});
}
BSplineCurve2d QuadraticVee()
{
    return Spline(2, {0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0},
                  {{-1, 1}, {-0.5, 0.4}, {0, 0}, {0.5, 0.4}, {1, 1}});
}

/// A cubic along y = x / 2 from (-1, -0.5) to (1, 0.5), its poles on the
/// line: it meets the vees at their corner only, inside their cone, so it
/// does not cross them.
BSplineCurve2d Slant()
{
    return Spline(3, cubic, {{-1, -0.5}, {-0.3, -0.15}, {0.3, 0.15}, {1, 0.5}});
}

/// The spline with every pole mapped by f.
template <typename F>
BSplineCurve2d Mapped(const BSplineCurve2d& s, F f)
{
    std::vector<Point2> poles;
    for (const Point2& p : s.Poles) {
        poles.push_back(f(p));
    }
    return Spline(s.Degree, s.Knots, std::move(poles), s.Weights);
}

/// The spline reversed: t -> lo + hi - t.
BSplineCurve2d Reversed(const BSplineCurve2d& s)
{
    const double lo = s.Knots.front();
    const double hi = s.Knots.back();
    std::vector<double> knots;
    for (auto k = s.Knots.rbegin(); k != s.Knots.rend(); ++k) {
        knots.push_back(lo + hi - *k);
    }
    std::vector<Point2> poles(s.Poles.rbegin(), s.Poles.rend());
    std::vector<double> weights(s.Weights.rbegin(), s.Weights.rend());
    return Spline(s.Degree, std::move(knots), std::move(poles), std::move(weights));
}

PlanarIntersection Meet(const BSplineCurve2d& a, const BSplineCurve2d& b)
{
    const auto r = Intersect(Curve2d{a}, Curve2d{b}, ctx);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    REQUIRE(r->Overlaps.empty());
    return *r;
}

/// One decided meeting of the contact at (ta, tb) and p: exactly at the
/// parameters the exact decisions give, to rounding elsewhere.
void RequireOnly(const PlanarIntersection& r, const double ta, const double tb, const Point2& p,
                 const ContactType contact)
{
    REQUIRE(r.Points.size() == 1);
    const PlanarPoint& m = r.Points.front();
    CHECK(m.Decided);
    CHECK(m.Contact == contact);
    CHECK(std::abs(m.ParameterA - ta) <= 1e-15 * (1.0 + std::abs(ta)));
    CHECK(std::abs(m.ParameterB - tb) <= 1e-15 * (1.0 + std::abs(tb)));
    const double scale = 1e-15 * (1.0 + std::abs(p.X) + std::abs(p.Y));
    CHECK(std::abs(m.Position.X - p.X) <= scale);
    CHECK(std::abs(m.Position.Y - p.Y) <= scale);
}

/// The case under exchange, reversal of either spline, a translation and
/// scaling by 2^10 and 2^-10, all exact on these poles, with where the
/// meeting goes.
void UnderRelations(const BSplineCurve2d& a, const BSplineCurve2d& b, const double ta,
                    const double tb, const Point2& p, const ContactType contact)
{
    const double flipA = a.Knots.front() + a.Knots.back() - ta;
    const double flipB = b.Knots.front() + b.Knots.back() - tb;
    {
        INFO("as given");
        RequireOnly(Meet(a, b), ta, tb, p, contact);
    }
    {
        INFO("exchanged");
        RequireOnly(Meet(b, a), tb, ta, p, contact);
    }
    {
        INFO("first reversed");
        RequireOnly(Meet(Reversed(a), b), flipA, tb, p, contact);
    }
    {
        INFO("second reversed");
        RequireOnly(Meet(a, Reversed(b)), ta, flipB, p, contact);
    }
    const auto moved = [](const Point2& q) { return Point2{q.X + 1024.0, q.Y - 512.0}; };
    {
        INFO("translated");
        RequireOnly(Meet(Mapped(a, moved), Mapped(b, moved)), ta, tb, moved(p), contact);
    }
    for (const double s : {0x1p10, 0x1p-10}) {
        INFO("scaled by " << s);
        const auto scaled = [s](const Point2& q) { return Point2{s * q.X, s * q.Y}; };
        RequireOnly(Meet(Mapped(a, scaled), Mapped(b, scaled)), ta, tb, scaled(p), contact);
    }
}

} // namespace

TEST_CASE("two B-splines crossing at a shared end pole are decided there",
          "[intersect][planar][spline-end]")
{
    // Arch ends where Descent starts; the curves cross there.
    UnderRelations(Arch(), Descent(), 1.0, 0.0, Point2{1, 1}, ContactType::Crossing);
}

TEST_CASE("a B-spline ending on a straight B-spline is decided exactly",
          "[intersect][planar][spline-end]")
{
    // On the post at its middle.
    UnderRelations(Arch(), Post(), 1.0, 0.5, Point2{1, 1}, ContactType::Crossing);
    // On a quadratic whose poles are on x = 1, where its y,
    // 2 t^2 + 2 t - 1, is 1: t = (sqrt 5 - 1) / 2, not a pole. The arch's
    // parameter is decided exactly from the line; the quadratic's is placed.
    const BSplineCurve2d straight =
        Spline(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {{1, -1}, {1, 0}, {1, 3}});
    UnderRelations(Arch(), straight, 1.0, (std::sqrt(5.0) - 1.0) / 2.0, Point2{1, 1},
                   ContactType::Crossing);
}

TEST_CASE("a B-spline ending a rounding either side of a line is decided",
          "[intersect][planar][spline-end]")
{
    const double eps = 0x1p-52;
    // Short of the post: no meeting (before the change, an undecided one).
    CHECK(Meet(Arch(1.0 - eps), Post()).IsEmpty());
    CHECK(Meet(Post(), Arch(1.0 - eps)).IsEmpty());
    CHECK(Meet(Reversed(Arch(1.0 - eps)), Post()).IsEmpty());
    // Past it: one decided crossing inside the arch, at its end to
    // rounding.
    for (const bool exchanged : {false, true}) {
        INFO((exchanged ? "exchanged" : "as given"));
        const PlanarIntersection r =
            exchanged ? Meet(Post(), Arch(1.0 + eps)) : Meet(Arch(1.0 + eps), Post());
        REQUIRE(r.Points.size() == 1);
        const PlanarPoint& m = r.Points.front();
        CHECK(m.Decided);
        CHECK(m.Contact == ContactType::Crossing);
        const double onArch = exchanged ? m.ParameterB : m.ParameterA;
        CHECK(onArch <= 1.0);
        CHECK(onArch > 1.0 - 1e-13);
    }
}

TEST_CASE("a touch at a corner of a B-spline is decided as a touch",
          "[intersect][planar][spline-end]")
{
    // The slant passes through the vees' corner inside their cone: it meets
    // them there only and does not cross them.
    for (const BSplineCurve2d& vee : {Vee(), QuadraticVee()}) {
        INFO("degree " << vee.Degree);
        UnderRelations(vee, Slant(), 1.0, 0.5, Point2{0, 0}, ContactType::Tangential);
    }
}

TEST_CASE("a corner on another spline's corner is never dropped", "[intersect][planar][spline-end]")
{
    // Both splines have a corner at the meeting: neither is straight across
    // the other's, and the contact is not decided. It fails closed today; a
    // result, if one is given, holds the corner.
    const BSplineCurve2d lambda = Spline(1, {0.0, 0.0, 1.0, 2.0, 2.0}, {{-1, -1}, {0, 0}, {1, -1}});
    const auto r = Intersect(Curve2d{Vee()}, Curve2d{lambda}, ctx);
    if (r) {
        REQUIRE(r->Points.size() == 1);
        CHECK(Distance(r->Points.front().Position, Point2{0, 0}) <= 1e-12);
    } else {
        CHECK(r.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
    }
}

TEST_CASE("a closed B-spline's ends on another's start are both decided",
          "[intersect][planar][spline-end]")
{
    const BSplineCurve2d loop = Spline(3, {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
                                       {{1, 1}, {0, 2}, {-1, 0}, {0, -1}, {1, 1}});
    const PlanarIntersection r = Meet(loop, Descent());
    REQUIRE(r.Points.size() == 2);
    for (const PlanarPoint& m : r.Points) {
        CHECK(m.Decided);
        CHECK(m.ParameterB == 0.0);
        CHECK(m.Position.X == 1.0);
        CHECK(m.Position.Y == 1.0);
    }
    CHECK(r.Points[0].ParameterA == 0.0);
    CHECK(r.Points[1].ParameterA == 1.0);
}

TEST_CASE("tangent B-splines at a shared end stay undecided", "[intersect][planar][spline-end]")
{
    // A G1 join: the end tangents are parallel and no exact decision
    // applies. Reported, never dropped, never decided.
    const BSplineCurve2d in = Spline(3, cubic, {{0, 0}, {0.3, 0.6}, {0.6, 1.2}, {1, 1}});
    const BSplineCurve2d out = Spline(3, cubic, {{1, 1}, {1.4, 0.8}, {1.7, 0.2}, {2, 0}});
    const PlanarIntersection r = Meet(in, out);
    REQUIRE(r.Points.size() == 1);
    CHECK_FALSE(r.Points.front().Decided);
    CHECK(Distance(r.Points.front().Position, Point2{1, 1}) <= 1e-12);
}
