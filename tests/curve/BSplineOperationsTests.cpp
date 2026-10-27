// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Clamping, the periodic form, degree elevation, knot removal and
// concatenation of B-spline curves. Each operation is checked against the
// evaluator on random curves (the point set must not move, except by the
// bound knot removal reports), against closed forms where one exists (the
// Bernstein coefficients of the monomials, the unit circle) and for the
// structure of the knots and poles it promises.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;

namespace {

constexpr double Tolerance = 1e-11;

bool Near(const Point3& a, const Point3& b, const double tolerance = Tolerance)
{
    return Distance(a, b) <= tolerance * (1.0 + Length(Stereon::ToVec(b)));
}

/// A random valid curve: degree 1 to 5, clamped or not, interior knots with
/// multiplicities up to the degree, optionally rational (as in
/// BSplineTests.cpp).
BSplineCurve RandomCurve(const bool rational)
{
    const int p = *rc::gen::inRange(1, 6);
    const auto q = static_cast<std::size_t>(p);
    const std::size_t poles = q + 1 + static_cast<std::size_t>(*rc::gen::inRange(0, 6));
    const bool clamped = *rc::gen::arbitrary<bool>();
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (clamped) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        points.push_back({*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                          *rc::gen::inRange(-50, 50) / 10.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(2, 40) / 10.0);
        }
    }
    // Clamping may have merged an interior knot into an end beyond p + 1
    // copies: not a curve.
    RC_PRE(Stereon::MakeBSplineCurve(p, knots, points, weights).has_value());
    return {p, knots, points, weights};
}

double RandomParameter(const BSplineCurve& c)
{
    const Interval d = Domain(c);
    const double s = *rc::gen::inRange(0, 1000) / 1000.0;
    return d.GetLower() + s * (d.GetUpper() - d.GetLower());
}

bool IsClamped(const BSplineCurve& c)
{
    const auto p = static_cast<std::size_t>(c.Degree);
    const std::size_t m = c.Knots.size() - 1;
    return c.Knots[0] == c.Knots[p] && c.Knots[m - p] == c.Knots[m];
}

/// The distinct knots strictly inside the domain.
std::vector<double> InteriorKnots(const BSplineCurve& c)
{
    const Interval d = Domain(c);
    std::vector<double> inside;
    for (const double u : c.Knots) {
        if (u > d.GetLower() && u < d.GetUpper() && (inside.empty() || inside.back() != u)) {
            inside.push_back(u);
        }
    }
    return inside;
}

/// The largest distance between the two curves' points at `count` parameters
/// spread over the first curve's domain.
double SampledDeviation(const BSplineCurve& a, const BSplineCurve& b, const int count = 64)
{
    const Interval d = Domain(a);
    double largest = 0.0;
    for (int i = 0; i <= count; ++i) {
        const double t = d.GetLower() + i * (d.GetUpper() - d.GetLower()) / count;
        largest =
            std::max(largest, Distance(Eval(a, t), Eval(b, std::min(t, Domain(b).GetUpper()))));
    }
    return largest;
}

/// (t, t^2, t^3) on [0, 1] as a cubic Bezier curve.
BSplineCurve MonomialCurve()
{
    return {3,
            {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
            {Point3{0.0, 0.0, 0.0}, Point3{1.0 / 3.0, 0.0, 0.0}, Point3{2.0 / 3.0, 1.0 / 3.0, 0.0},
             Point3{1.0, 1.0, 1.0}},
            {}};
}

/// The unit circle as the classical nine-pole quadratic NURBS: four quarter
/// arcs on [0, 4], the corner poles with weight sqrt(2)/2.
BSplineCurve UnitCircle()
{
    const double w = std::numbers::sqrt2 / 2.0;
    return {2,
            {0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0},
            {Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}, Point3{0.0, 1.0, 0.0},
             Point3{-1.0, 1.0, 0.0}, Point3{-1.0, 0.0, 0.0}, Point3{-1.0, -1.0, 0.0},
             Point3{0.0, -1.0, 0.0}, Point3{1.0, -1.0, 0.0}, Point3{1.0, 0.0, 0.0}},
            {1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0}};
}

/// Ten points round a bumpy loop, for the periodic fits.
std::vector<Point3> LoopPoints()
{
    std::vector<Point3> points;
    for (int k = 0; k < 10; ++k) {
        const double angle = 2.0 * std::numbers::pi * k / 10.0;
        const double radius = 1.0 + 0.3 * std::cos(3.0 * angle);
        points.push_back(
            {radius * std::cos(angle), radius * std::sin(angle), 0.2 * std::sin(2.0 * angle)});
    }
    return points;
}

BSplineCurve2d Flattened(const BSplineCurve& c)
{
    BSplineCurve2d flat{c.Degree, c.Knots, {}, c.Weights};
    for (const Point3& pole : c.Poles) {
        flat.Poles.push_back({pole.X, pole.Y});
    }
    return flat;
}

} // namespace

// ---------------------------------------------------------------------------
// Clamping

TEST_CASE("Clamped keeps the curve on its domain and clamps both ends",
          "[curve][bspline][operations][property]")
{
    rc::prop("the clamped form evaluates as the original", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const BSplineCurve clamped = Stereon::Clamped(c);
        RC_ASSERT(IsClamped(clamped));
        RC_ASSERT(Domain(clamped) == Domain(c));
        RC_ASSERT(clamped.Degree == c.Degree);
        RC_ASSERT(IsRational(clamped) == IsRational(c));
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(clamped, t), Eval(c, t)));
        RC_ASSERT(Stereon::Clamped(clamped) == clamped);
        RC_ASSERT(InteriorKnots(clamped) == InteriorKnots(c));
        if (IsClamped(c)) {
            RC_ASSERT(clamped == c);
        }
    });
}

// ---------------------------------------------------------------------------
// Degree elevation

TEST_CASE("ElevatedDegree traces the same curve with every multiplicity raised",
          "[curve][bspline][operations][property]")
{
    rc::prop("the elevated curve evaluates as the original", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const int raise = *rc::gen::inRange(1, 4);
        const int degree = c.Degree + raise;
        const BSplineCurve e = *Stereon::ElevatedDegree(c, degree);
        RC_ASSERT(e.Degree == degree);
        RC_ASSERT(IsClamped(e));
        RC_ASSERT(Domain(e) == Domain(c));
        RC_ASSERT(IsRational(e) == IsRational(c));
        RC_ASSERT(e.Poles.size() + static_cast<std::size_t>(degree) + 1 == e.Knots.size());
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(e, t), Eval(c, t), 1e-10));
        RC_ASSERT(InteriorKnots(e) == InteriorKnots(c));
        for (const double u : InteriorKnots(c)) {
            RC_ASSERT(Stereon::KnotMultiplicity(e.Knots, u) ==
                      Stereon::KnotMultiplicity(c.Knots, u) + static_cast<std::size_t>(raise));
        }
        // The continuity the knots guarantee is unchanged; a single span
        // reports its degree, which rose.
        RC_ASSERT(Stereon::Continuity(e) ==
                  (InteriorKnots(c).empty() ? degree : Stereon::Continuity(c)));
        RC_ASSERT(*Stereon::ElevatedDegree(c, c.Degree) == c);
    });
}

TEST_CASE("elevating the monomial cubic gives the Bernstein coefficients of degree five",
          "[curve][bspline][operations]")
{
    const BSplineCurve e = *Stereon::ElevatedDegree(MonomialCurve(), 5);
    REQUIRE(e.Degree == 5);
    REQUIRE(e.Poles.size() == 6);
    // t^k in the degree-5 Bernstein basis has coefficients C(i, k) / C(5, k).
    for (std::size_t i = 0; i < 6; ++i) {
        const auto d = static_cast<double>(i);
        const Point3 expected{d / 5.0, d * (d - 1.0) / 20.0, d * (d - 1.0) * (d - 2.0) / 60.0};
        REQUIRE(Distance(e.Poles[i], expected) <= 1e-15);
    }
    for (std::size_t i = 0; i < 6; ++i) {
        REQUIRE(e.Knots[i] == 0.0);
        REQUIRE(e.Knots[i + 6] == 1.0);
    }
}

TEST_CASE("elevating the rational circle keeps it on the circle", "[curve][bspline][operations]")
{
    const BSplineCurve circle = UnitCircle();
    const BSplineCurve e = *Stereon::ElevatedDegree(circle, 4);
    REQUIRE(e.Degree == 4);
    REQUIRE(IsRational(e));
    REQUIRE(e.Poles.size() == 17); // one quarter of degree 4 (5 poles) each, sharing corners
    for (const double u : {1.0, 2.0, 3.0}) {
        REQUIRE(Stereon::KnotMultiplicity(e.Knots, u) == 4);
    }
    for (int i = 0; i <= 40; ++i) {
        const double t = i * 0.1;
        const Point3 p = Eval(e, t);
        REQUIRE(std::abs(Length(Stereon::ToVec(p)) - 1.0) <= 1e-14);
        REQUIRE(Near(p, Eval(circle, t), 1e-14));
    }
    REQUIRE(Stereon::ElevatedDegree(circle, 1).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(Stereon::ElevatedDegree(circle, 26).error().GetCode() == ErrorCode::InvalidInput);
}

// ---------------------------------------------------------------------------
// Knot removal

TEST_CASE("RemoveKnot undoes an insertion", "[curve][bspline][operations][property]")
{
    rc::prop("inserting then removing a knot gives the knot vector back", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double u = RandomParameter(c);
        const Interval d = Domain(c);
        RC_PRE(u > d.GetLower() && u < d.GetUpper());
        const auto inserted = Stereon::InsertKnot(c, u);
        RC_PRE(inserted.has_value());
        const auto removal = *Stereon::RemoveKnot(*inserted, u, 1, 1e-8);
        RC_ASSERT(removal.Removed == 1);
        RC_ASSERT(removal.Deviation <= 1e-8);
        RC_ASSERT(removal.Curve.Knots == c.Knots);
        RC_ASSERT(removal.Curve.Poles.size() == c.Poles.size());
        for (std::size_t i = 0; i < c.Poles.size(); ++i) {
            RC_ASSERT(Near(removal.Curve.Poles[i], c.Poles[i], 1e-8));
        }
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(removal.Curve, t), Eval(c, t), 1e-8));
    });
}

TEST_CASE("RemoveKnot moves the curve by no more than the bound it reports",
          "[curve][bspline][operations][property]")
{
    rc::prop("the sampled deviation is within Deviation, which is within the tolerance", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> inside = InteriorKnots(c);
        RC_PRE(!inside.empty());
        const double u =
            inside[static_cast<std::size_t>(*rc::gen::inRange(0, static_cast<int>(inside.size())))];
        const auto times = static_cast<int>(Stereon::KnotMultiplicity(c.Knots, u));
        const double tolerance = *rc::gen::inRange(0, 200) / 100.0;
        const auto removal = *Stereon::RemoveKnot(c, u, times, tolerance);
        RC_ASSERT(removal.Removed >= 0 && removal.Removed <= times);
        RC_ASSERT(removal.Deviation <= tolerance);
        RC_ASSERT(removal.Curve.Knots.size() + static_cast<std::size_t>(removal.Removed) ==
                  c.Knots.size());
        RC_ASSERT(Stereon::KnotMultiplicity(removal.Curve.Knots, u) ==
                  static_cast<std::size_t>(times - removal.Removed));
        if (removal.Removed == 0) {
            RC_ASSERT(removal.Curve == c);
            RC_ASSERT(removal.Deviation == 0.0);
        }
        RC_ASSERT(SampledDeviation(c, removal.Curve) <= removal.Deviation + 1e-9);
    });
}

TEST_CASE("RemoveKnot refuses a corner within a small tolerance and takes it within a large one",
          "[curve][bspline][operations]")
{
    const BSplineCurve corner{1,
                              {0.0, 0.0, 1.0, 2.0, 2.0},
                              {Point3{0.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}, Point3{2.0, 0.0, 0.0}},
                              {}};
    const auto kept = *Stereon::RemoveKnot(corner, 1.0, 1, 0.5);
    REQUIRE(kept.Removed == 0);
    REQUIRE(kept.Deviation == 0.0);
    REQUIRE(kept.Curve == corner);
    const auto removed = *Stereon::RemoveKnot(corner, 1.0, 1, 1.0);
    REQUIRE(removed.Removed == 1);
    REQUIRE(removed.Deviation == 1.0);
    REQUIRE(removed.Curve.Knots == std::vector{0.0, 0.0, 2.0, 2.0});
    REQUIRE(Near(Eval(removed.Curve, 1.0), Point3{1.0, 0.0, 0.0}));

    const BSplineCurve c = MonomialCurve();
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::RemoveKnot(c, 0.5, 1, 1e-9).error().GetCode() == invalid); // not a knot
    REQUIRE(Stereon::RemoveKnot(c, 0.0, 1, 1e-9).error().GetCode() == invalid); // a domain end
    const BSplineCurve inserted = *Stereon::InsertKnot(c, 0.5, 2);
    REQUIRE(Stereon::RemoveKnot(inserted, 0.5, 0, 1e-9).error().GetCode() == invalid);
    REQUIRE(Stereon::RemoveKnot(inserted, 0.5, 3, 1e-9).error().GetCode() == invalid);
    REQUIRE(Stereon::RemoveKnot(inserted, 0.5, 1, -1.0).error().GetCode() == invalid);
    REQUIRE(Stereon::RemoveKnot(inserted, NAN, 1, 1e-9).error().GetCode() == invalid);
    const auto both = *Stereon::RemoveKnot(inserted, 0.5, 2, 1e-12);
    REQUIRE(both.Removed == 2);
    REQUIRE(both.Curve.Knots == c.Knots);
}

// ---------------------------------------------------------------------------
// Queries

TEST_CASE("Continuity and IsClosed read the knots and the ends", "[curve][bspline][operations]")
{
    const std::vector<Point3> poles(7, Point3{0.0, 0.0, 0.0});
    const BSplineCurve c{3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0}, poles, {}};
    REQUIRE(Stereon::Continuity(c) == 1);
    REQUIRE(Stereon::Continuity(c, 1.0) == 2);
    REQUIRE(Stereon::Continuity(c, 2.0) == 1);
    REQUIRE(Stereon::Continuity(c, 0.5) == 3);
    REQUIRE(Stereon::Continuity(c, 0.0) == 3);
    REQUIRE(Stereon::Continuity(c, 3.0) == 3);
    REQUIRE(Stereon::Continuity(MonomialCurve()) == 3);
    REQUIRE(Stereon::Continuity(UnitCircle()) == 0);
    const BSplineCurve unclamped{
        2, {-1.0, -0.5, 0.0, 0.75, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0}, poles, {}};
    REQUIRE(Stereon::Continuity(unclamped) == 1);       // 0.75, 1, 2, 2.5 are simple
    REQUIRE(Stereon::Continuity(unclamped, -0.5) == 2); // outside the domain [0, 3]

    REQUIRE(Stereon::IsClosed(UnitCircle(), 0.0));
    REQUIRE_FALSE(Stereon::IsClosed(MonomialCurve(), 1.0));
    REQUIRE(Stereon::IsClosed(MonomialCurve(), 2.0));
    REQUIRE_FALSE(Stereon::IsClosed(MonomialCurve(), NAN));
}

// ---------------------------------------------------------------------------
// The periodic form

TEST_CASE("IsPeriodic recognises the periodic fits, and AsPeriodic rebuilds the form from the "
          "clamped curve",
          "[curve][bspline][operations]")
{
    const std::vector<Point3> points = LoopPoints();
    for (int degree = 1; degree <= 5; ++degree) {
        const BSplineCurve fit = Stereon::InterpolatePeriodic(points, degree)->Curve;
        REQUIRE(Stereon::IsPeriodic(fit));
        REQUIRE(Stereon::IsClosed(fit, 1e-12));
        const BSplineCurve clamped = Stereon::Clamped(fit);
        if (degree >= 2) {
            REQUIRE_FALSE(Stereon::IsPeriodic(clamped));
        }
        const BSplineCurve periodic = *Stereon::AsPeriodic(clamped, 1e-9);
        REQUIRE(Stereon::IsPeriodic(periodic));
        REQUIRE(Domain(periodic) == Domain(fit));
        REQUIRE(Stereon::KnotMultiplicity(periodic.Knots, 0.0) == static_cast<std::size_t>(degree));
        for (int i = 0; i <= 50; ++i) {
            const double t = i / 50.0;
            REQUIRE(Near(Eval(periodic, t), Eval(fit, t), 1e-9));
        }
        REQUIRE(*Stereon::AsPeriodic(periodic, 0.0) == periodic);
    }
}

TEST_CASE("AsPeriodic closes the unit circle exactly and refuses what does not close",
          "[curve][bspline][operations]")
{
    const BSplineCurve circle = UnitCircle();
    const BSplineCurve periodic = *Stereon::AsPeriodic(circle, 0.0);
    REQUIRE(Stereon::IsPeriodic(periodic));
    REQUIRE(periodic.Poles.size() == 10); // 8 distinct poles and the first 2 again
    REQUIRE(periodic.Knots.size() == 13);
    REQUIRE(Domain(periodic) == Interval{0.0, 4.0});
    for (int i = 0; i <= 40; ++i) {
        const double t = i * 0.1;
        REQUIRE(Near(Eval(periodic, t), Eval(circle, t), 1e-14));
    }
    // The seam has multiplicity 2 = p, and the knots repeat with the period.
    REQUIRE(Stereon::KnotMultiplicity(periodic.Knots, 0.0) == 2);
    REQUIRE(Stereon::KnotMultiplicity(periodic.Knots, 4.0) == 2);
    REQUIRE(periodic.Knots.front() == -1.0);
    REQUIRE(periodic.Knots.back() == 5.0);

    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::AsPeriodic(MonomialCurve(), 0.5).error().GetCode() == invalid);
    REQUIRE(Stereon::AsPeriodic(circle, -1.0).error().GetCode() == invalid);
    BSplineCurve reweighted = circle;
    reweighted.Weights.back() = 2.0; // the same point, another weight
    REQUIRE(Stereon::AsPeriodic(reweighted, 0.0).error().GetCode() == invalid);
    // Within tolerance, the end poles are merged at their midpoint.
    BSplineCurve open = MonomialCurve();
    open.Poles.back() = {0.1, 0.0, 0.0};
    const BSplineCurve merged = *Stereon::AsPeriodic(open, 0.1);
    REQUIRE(Stereon::IsPeriodic(merged));
    REQUIRE(Near(Eval(merged, 0.0), Point3{0.05, 0.0, 0.0}));
    REQUIRE(Near(Eval(merged, 1.0), Point3{0.05, 0.0, 0.0}));
}

TEST_CASE("the periodic form of any closed curve traces it",
          "[curve][bspline][operations][property]")
{
    rc::prop("a clamped curve whose ends are made to coincide", [] {
        BSplineCurve c = Stereon::Clamped(RandomCurve(*rc::gen::arbitrary<bool>()));
        c.Poles.back() = c.Poles.front();
        if (IsRational(c)) {
            c.Weights.back() = c.Weights.front();
        }
        const BSplineCurve periodic = *Stereon::AsPeriodic(c, 0.0);
        RC_ASSERT(Stereon::IsPeriodic(periodic));
        RC_ASSERT(Domain(periodic) == Domain(c));
        RC_ASSERT(periodic.Poles.size() == c.Poles.size() - 1 + static_cast<std::size_t>(c.Degree));
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(periodic, t), Eval(c, t)));
        RC_ASSERT(Near(Eval(periodic, Domain(c).GetLower()), c.Poles.front()));
        RC_ASSERT(Near(Eval(periodic, Domain(c).GetUpper()), c.Poles.front()));
    });
}

// ---------------------------------------------------------------------------
// Concatenation

TEST_CASE("Concatenated joins the two halves of a split curve, and knot removal at the junction "
          "restores it",
          "[curve][bspline][operations][property]")
{
    rc::prop("Split then Concatenated then RemoveKnot", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const double u = RandomParameter(c);
        const auto halves = Stereon::Split(c, u);
        RC_PRE(halves.has_value());
        const BSplineCurve joined = *Stereon::Concatenated(halves->first, halves->second, 0.0);
        RC_ASSERT(joined.Degree == c.Degree);
        RC_ASSERT(IsClamped(joined));
        RC_ASSERT(Domain(joined).GetLower() == Domain(c).GetLower());
        RC_ASSERT(std::abs(Domain(joined).GetUpper() - Domain(c).GetUpper()) <= 1e-12);
        const double t = RandomParameter(c);
        RC_ASSERT(Near(Eval(joined, std::min(t, Domain(joined).GetUpper())), Eval(c, t), 1e-10));
        const auto p = static_cast<std::size_t>(c.Degree);
        RC_ASSERT(Stereon::KnotMultiplicity(joined.Knots, u) == p);
        // The original is as smooth at u as its own multiplicity there says,
        // so that many copies of the junction knot come out exactly.
        const std::size_t removable = p - Stereon::KnotMultiplicity(c.Knots, u);
        RC_PRE(removable > 0);
        const auto removal = *Stereon::RemoveKnot(joined, u, static_cast<int>(removable), 1e-7);
        RC_ASSERT(removal.Removed == static_cast<int>(removable));
        const BSplineCurve clamped = Stereon::Clamped(c);
        RC_ASSERT(removal.Curve.Poles.size() == clamped.Poles.size());
        for (std::size_t i = 0; i < clamped.Poles.size(); ++i) {
            RC_ASSERT(Near(removal.Curve.Poles[i], clamped.Poles[i], 1e-7));
        }
    });
}

TEST_CASE("Concatenated elevates the lower degree, shifts the parameter and merges the junction",
          "[curve][bspline][operations]")
{
    const BSplineCurve line{
        1, {0.0, 0.0, 1.0, 1.0}, {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}}, {}};
    BSplineCurve cubic = MonomialCurve();
    for (Point3& pole : cubic.Poles) {
        pole.X += 1.0;
    }
    const BSplineCurve joined = *Stereon::Concatenated(line, cubic, 0.0);
    REQUIRE(joined.Degree == 3);
    REQUIRE(Domain(joined) == Interval{0.0, 2.0});
    REQUIRE(Stereon::Continuity(joined) == 0);
    REQUIRE(Stereon::KnotMultiplicity(joined.Knots, 1.0) == 3);
    REQUIRE(Near(Eval(joined, 0.5), Point3{0.5, 0.0, 0.0}));
    REQUIRE(Near(Eval(joined, 1.0), Point3{1.0, 0.0, 0.0}));
    REQUIRE(Near(Eval(joined, 1.5), Eval(cubic, 0.5)));
    REQUIRE(Near(Eval(joined, 2.0), Point3{2.0, 1.0, 1.0}));
    // The line and the cubic meet tangentially: the junction knot comes out once.
    const auto removal = *Stereon::RemoveKnot(joined, 1.0, 3, 1e-12);
    REQUIRE(removal.Removed == 1);

    // A gap: refused within a small tolerance, bridged at the midpoint within a large one.
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::Concatenated(line, MonomialCurve(), 0.5).error().GetCode() == invalid);
    REQUIRE(Stereon::Concatenated(line, MonomialCurve(), -1.0).error().GetCode() == invalid);
    const BSplineCurve bridged = *Stereon::Concatenated(line, MonomialCurve(), 1.0);
    REQUIRE(Near(Eval(bridged, 1.0), Point3{0.5, 0.0, 0.0}));

    // A rational quarter circle followed by a polynomial line: one rational curve.
    const BSplineCurve quarter = Stereon::ToBezierSegments(UnitCircle()).front();
    const BSplineCurve up{
        1, {0.0, 0.0, 1.0, 1.0}, {Point3{0.0, 1.0, 0.0}, Point3{0.0, 2.0, 0.0}}, {}};
    const BSplineCurve arcAndLine = *Stereon::Concatenated(quarter, up, 0.0);
    REQUIRE(IsRational(arcAndLine));
    REQUIRE(arcAndLine.Degree == 2);
    REQUIRE(Near(Eval(arcAndLine, 0.5), Eval(quarter, 0.5), 1e-14));
    REQUIRE(Near(Eval(arcAndLine, 1.5), Point3{0.0, 1.5, 0.0}, 1e-14));
    const BSplineCurve lineAndArc = *Stereon::Concatenated(Reversed(up), Reversed(quarter), 0.0);
    REQUIRE(Near(Eval(lineAndArc, 0.5), Point3{0.0, 1.5, 0.0}, 1e-14));
    REQUIRE(Near(Eval(lineAndArc, 1.5), Eval(quarter, 0.5), 1e-14));
}

// ---------------------------------------------------------------------------
// The plane

TEST_CASE("the operations work in the parameter plane too", "[curve][bspline][operations]")
{
    const BSplineCurve2d circle = Flattened(UnitCircle());
    REQUIRE(Stereon::IsClosed(circle, 0.0));
    REQUIRE(Stereon::Continuity(circle) == 0);
    REQUIRE(Stereon::Continuity(circle, 0.5) == 2);
    const BSplineCurve2d periodic = *Stereon::AsPeriodic(circle, 0.0);
    REQUIRE(Stereon::IsPeriodic(periodic));
    REQUIRE_FALSE(Stereon::IsPeriodic(circle));
    REQUIRE(Stereon::Clamped(periodic).Knots.front() == 0.0);
    const BSplineCurve2d elevated = *Stereon::ElevatedDegree(circle, 3);
    REQUIRE(elevated.Degree == 3);
    REQUIRE(Distance(Eval(elevated, 0.5), Eval(circle, 0.5)) <= 1e-14);
    const BSplineCurve2d inserted = *Stereon::InsertKnot(circle, 0.5);
    const auto removal = *Stereon::RemoveKnot(inserted, 0.5, 1, 1e-12);
    REQUIRE(removal.Removed == 1);
    REQUIRE(removal.Curve.Knots == circle.Knots);
    const BSplineCurve2d line{1, {0.0, 0.0, 1.0, 1.0}, {Point2{1.0, 0.0}, Point2{2.0, 0.0}}, {}};
    const BSplineCurve2d joined = *Stereon::Concatenated(circle, line, 0.0);
    REQUIRE(Domain(joined) == Interval{0.0, 5.0});
    REQUIRE(Distance(Eval(joined, 4.5), Point2{1.5, 0.0}) <= 1e-14);
}
