// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Two B-spline curves in space (SplineCurvePair.cpp): a meeting is a
// critical point of the distance within the resolution. Straight splines
// crossing, passing apart and passing within the resolution; edges sharing
// a vertex, exactly and a rounding apart; a rational quarter circle crossed
// at 45 degrees; random cubics crossed by a straight spline through one of
// their points, the point found and every point reported within the
// resolution, in either order; identical splines an overlap; and tangent
// curves failing closed, never dropped.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::BSplineCurve;
using Stereon::Curve;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Vec3;

namespace {

/// The straight degree-1 spline from p to q over [0, 1].
BSplineCurve Segment(const Point3& p, const Point3& q)
{
    return BSplineCurve{1, {0.0, 0.0, 1.0, 1.0}, {p, q}, {}};
}

/// A cubic Bezier over [0, 1].
BSplineCurve Cubic(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3)
{
    return BSplineCurve{3, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0}, {p0, p1, p2, p3}, {}};
}

Stereon::Result<Stereon::CurveCurveIntersection> Meet(const BSplineCurve& a, const BSplineCurve& b)
{
    return Stereon::Intersect(Curve{a}, Curve{b}, ctx);
}

double Coordinate()
{
    return *rc::gen::inRange(-40, 41) / 8.0;
}

} // namespace

TEST_CASE("straight splines in space meet where they cross, or within the resolution",
          "[intersect][curve-curve][spline]")
{
    const BSplineCurve a = Segment(Point3{-1.0, 0.0, 0.0}, Point3{3.0, 0.0, 0.0});
    SECTION("crossing")
    {
        const auto r = Meet(a, Segment(Point3{0.5, -2.0, 0.0}, Point3{0.5, 2.0, 0.0}));
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterA, WithinAbs(0.375, 1e-15));
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(0.5, 1e-15));
        CHECK(Distance(r->Points[0].Position, Point3{0.5, 0.0, 0.0}) <= 1e-15);
    }
    SECTION("skew, a millimetre apart: no meeting")
    {
        const auto r = Meet(a, Segment(Point3{0.5, -2.0, 1.0}, Point3{0.5, 2.0, 1.0}));
        REQUIRE(r.has_value());
        CHECK(r->IsEmpty());
    }
    SECTION("skew, within the resolution: a meeting at the closest approach")
    {
        const auto r = Meet(a, Segment(Point3{0.5, -2.0, 4e-6}, Point3{0.5, 2.0, 4e-6}));
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        CHECK_THAT(r->Points[0].ParameterA, WithinAbs(0.375, 1e-12));
    }
    SECTION("skew, beyond the resolution: no meeting")
    {
        const auto r = Meet(a, Segment(Point3{0.5, -2.0, 3e-5}, Point3{0.5, 2.0, 3e-5}));
        REQUIRE(r.has_value());
        CHECK(r->IsEmpty());
    }
    SECTION("crossing outside the requested domains: no meeting")
    {
        const auto r = Stereon::Intersect(
            Curve{a}, Interval{0.5, 1.0},
            Curve{Segment(Point3{0.5, -2.0, 0.0}, Point3{0.5, 2.0, 0.0})}, Interval{0.0, 1.0}, ctx);
        REQUIRE(r.has_value());
        CHECK(r->IsEmpty());
    }
}

TEST_CASE("edges sharing a vertex meet at it, also a rounding apart",
          "[intersect][curve-curve][spline]")
{
    const Point3 vertex{1.25, -0.5, 2.0};
    const BSplineCurve a =
        Cubic(Point3{-2.0, 0.0, 0.0}, Point3{-1.0, 1.0, 0.5}, Point3{0.0, 0.5, 1.5}, vertex);
    for (const double offset : {0.0, 1e-13, 1e-9}) {
        INFO("offset " << offset);
        const Point3 start = vertex + Vec3{offset, -offset, offset};
        const BSplineCurve b =
            Cubic(start, Point3{2.0, -1.0, 2.5}, Point3{3.0, 0.0, 2.0}, Point3{4.0, 1.0, 1.0});
        const auto r = Meet(a, b);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        // The curves cross at 12.5 degrees there, so the closest approach
        // is placed to the rounding of the distance's gradient over the
        // square of that angle's sine: about 1e-11, far within the
        // resolution. With the ends apart, it is the closest approach,
        // within the offset of the vertex.
        CHECK_THAT(r->Points[0].ParameterA, WithinAbs(1.0, 1e-10));
        CHECK_THAT(r->Points[0].ParameterB, WithinAbs(0.0, 1e-10));
        CHECK(Distance(r->Points[0].Position, vertex) <= 1e-10 + 4.0 * offset);
    }
}

TEST_CASE("a rational quarter circle is crossed where the closed form says",
          "[intersect][curve-curve][spline]")
{
    // The quarter of the circle of radius 2 about the origin, from (2, 0, 1)
    // to (0, 2, 1), as a rational quadratic; a straight spline along the
    // diagonal at height 1 crosses it at (sqrt 2, sqrt 2, 1).
    const double w = std::numbers::sqrt2 / 2.0;
    const BSplineCurve arc{2,
                           {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                           {Point3{2.0, 0.0, 1.0}, Point3{2.0, 2.0, 1.0}, Point3{0.0, 2.0, 1.0}},
                           {1.0, w, 1.0}};
    const BSplineCurve diagonal = Segment(Point3{0.0, 0.0, 1.0}, Point3{3.0, 3.0, 1.0});
    for (const bool swapped : {false, true}) {
        const auto r = swapped ? Meet(diagonal, arc) : Meet(arc, diagonal);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Points.size() == 1);
        const Point3 expected{std::numbers::sqrt2, std::numbers::sqrt2, 1.0};
        CHECK(Distance(r->Points[0].Position, expected) <= 1e-14);
        CHECK_THAT(swapped ? r->Points[0].ParameterB : r->Points[0].ParameterA,
                   WithinAbs(0.5, 1e-14));
    }
}

TEST_CASE("random cubics crossed by a straight spline through one of their points",
          "[intersect][curve-curve][spline][property]")
{
    rc::prop("the point is found, and every point reported is a meeting", [] {
        const BSplineCurve a = Cubic(Point3{Coordinate(), Coordinate(), Coordinate()},
                                     Point3{Coordinate(), Coordinate(), Coordinate()},
                                     Point3{Coordinate(), Coordinate(), Coordinate()},
                                     Point3{Coordinate(), Coordinate(), Coordinate()});
        const double s0 = *rc::gen::inRange(1, 16) / 16.0;
        const Point3 on = Eval(a, s0);
        const Vec3 tangent = Deriv(a, s0, 1);
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        // Steep to the cubic, so the crossing is well conditioned.
        RC_PRE(Length(tangent) >= 0.5 && Length(direction) >= 0.5);
        RC_PRE(Length(Cross(tangent, direction)) >= 0.5 * Length(tangent) * Length(direction));
        const BSplineCurve b = Segment(on - 2.0 * direction, on + 2.0 * direction);
        for (const bool swapped : {false, true}) {
            const auto r = swapped ? Meet(b, a) : Meet(a, b);
            if (!r) {
                RC_FAIL(IntersectTest::ErrorText(r));
            }
            bool found = false;
            for (const Stereon::CurveCurvePoint& p : r->Points) {
                const double sa = swapped ? p.ParameterB : p.ParameterA;
                const double sb = swapped ? p.ParameterA : p.ParameterB;
                RC_ASSERT(Distance(Eval(a, sa), Eval(b, sb)) <= 1e-5);
                found = found || (std::abs(sa - s0) <= 1e-9 && std::abs(sb - 0.5) <= 1e-9);
            }
            RC_ASSERT(found);
        }
    });
}

TEST_CASE("identical splines overlap, tangent ones fail closed", "[intersect][curve-curve][spline]")
{
    const BSplineCurve a = Cubic(Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0},
                                 Point3{2.0, -1.0, 1.0}, Point3{3.0, 0.0, 1.0});
    SECTION("identical")
    {
        const auto r =
            Stereon::Intersect(Curve{a}, Interval{0.25, 1.0}, Curve{a}, Interval{0.0, 0.75}, ctx);
        REQUIRE(r.has_value());
        REQUIRE(r->Overlaps.size() == 1);
        CHECK(r->Overlaps[0].OnA == Interval{0.25, 0.75});
        CHECK(r->Points.empty());
    }
    SECTION("a parabola touching a line: undecidable at a double root")
    {
        // y = x² on [-1, 1] in the plane z = 0 touches the x axis at the
        // origin: the distance's critical point there is degenerate.
        const BSplineCurve parabola{
            2,
            {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
            {Point3{-1.0, 1.0, 0.0}, Point3{0.0, -1.0, 0.0}, Point3{1.0, 1.0, 0.0}},
            {}};
        const auto r = Meet(parabola, Segment(Point3{-2.0, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}
