// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

namespace {

using namespace Stereon;
constexpr double pi = std::numbers::pi;

// MR-01: proper cyclic axis permutations, plus exactly represented offsets.
// MR-03: powers-of-two scale, with dimensional policy scaled separately.
struct Placement {
    double Scale;
    Vec3 Offset;
    Vec3 X;
    Vec3 Y;

    [[nodiscard]] Point3 Point(const double x, const double y, const double z) const
    {
        return Point3{} + Offset + Scale * (x * X + y * Y + z * Cross(X, Y));
    }

    [[nodiscard]] Context Policy() const
    {
        return Context::Default()
            .WithResolution(Context::DefaultResolution * Scale)
            .WithToleranceCeiling(Context::DefaultToleranceCeiling * Scale)
            .WithScheduler(Scheduler::Serial);
    }
};

std::vector<Placement> Placements()
{
    const std::array<Vec3, 3> axes{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    std::vector<Placement> result;
    for (const double scale : {0.25, 1.0, 4.0}) {
        for (std::size_t i = 0; i < axes.size(); ++i) {
            for (const Vec3 offset : {Vec3{}, Vec3{8, -4, 2}}) {
                result.push_back({scale, offset, axes[i], axes[(i + 1) % axes.size()]});
            }
        }
    }
    return result;
}

void Replay(const char* defect, const Placement& p)
{
    // INFO in the caller retains the full fixture; this checks the effective
    // numerical policy rather than letting a scale test change its meaning.
    CAPTURE(defect, p.Scale);
    CHECK(p.Policy().GetResolution() == Context::DefaultResolution * p.Scale);
    CHECK(p.Policy().GetToleranceCeiling() == Context::DefaultToleranceCeiling * p.Scale);
    CHECK(p.Policy().GetAngularResolution() == Context::DefaultAngularResolution);
}

CurveSurfaceIntersection Crossings(const Curve& line, const Surface& surface, const Placement& p,
                                   const double radius, const double height)
{
    const auto result = Intersect(line, Domain(line), surface, Domain(surface), p.Policy());
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    REQUIRE(result->Points.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        const auto& point = result->Points[i];
        const Point3 expected = p.Point(i == 0 ? -radius : radius, 0, height);
        CHECK(Distance(point.Position, expected) <= p.Policy().GetResolution());
        CHECK(Distance(Eval(line, point.Parameter), expected) <= p.Policy().GetResolution());
        CHECK(Distance(Eval(surface, point.OnSurface), expected) <= p.Policy().GetResolution());
        CHECK(point.Contact == ContactType::Crossing);
        CHECK(point.Decided);
    }
    return *result;
}

void FarOrigin(const bool cone)
{
    for (const auto& p : Placements()) {
        INFO("defect=" << (cone ? "PR-218" : "PR-205") << " scale=" << p.Scale
                       << " offset=" << p.Offset.X << ',' << p.Offset.Y << ',' << p.Offset.Z
                       << " axis=" << p.X.X << ',' << p.X.Y << ',' << p.X.Z);
        Replay(cone ? "PR-218" : "PR-205", p);
        const double height = cone ? 0.5 : 0.0;
        const double radius = cone ? 1.5 : 1.0;
        const Surface surface = cone ? Surface{Cone{p.Point(0, 0, 0), p.Scale, pi / 4.0, p.X, p.Y}}
                                     : Surface{Sphere{p.Point(0, 0, 0), p.Scale, p.X, p.Y}};
        const Curve near{Line{p.Point(-2, 0, height), p.Scale * p.X}};
        const auto base = Crossings(near, surface, p, radius, height);
        for (const double run : {0x1p20, 0x1p27}) {
            CAPTURE(run);
            const Curve far{Line{p.Point(-2 - run, 0, height), p.Scale * p.X}};
            const auto moved = Crossings(far, surface, p, radius, height);
            for (std::size_t i = 0; i < 2; ++i) {
                CHECK(Distance(base.Points[i].Position, moved.Points[i].Position) <=
                      p.Policy().GetResolution());
                // MR-04: line-origin shift changes t by run, not by scale*run.
                CHECK(std::abs(base.Points[i].Parameter + run - moved.Points[i].Parameter) <=
                      Context::DefaultResolution);
            }
        }
    }
}

} // namespace

TEST_CASE("MR PR-205 far-origin sphere crossings", "[intersect][metamorphic][backfill][PR-205]")
{
    FarOrigin(false);
}

TEST_CASE("MR PR-218 far-origin cone crossings", "[intersect][metamorphic][backfill][PR-218]")
{
    FarOrigin(true);
}

TEST_CASE("MR PR-206 singleton circle plane domain", "[intersect][metamorphic][backfill][PR-206]")
{
    for (const auto& p : Placements()) {
        INFO("PR-206 scale=" << p.Scale << " axis=" << p.X.X << ',' << p.X.Y << ',' << p.X.Z
                             << " offset=" << p.Offset.X << ',' << p.Offset.Y << ',' << p.Offset.Z);
        Replay("PR-206", p);
        const Surface plane{Plane{p.Point(0, 0, 0), p.X, p.Y}};
        const Box2 band{Interval{-2 * p.Scale, 2 * p.Scale}, Interval{-p.Scale, p.Scale}};
        for (const bool reversedChart : {false, true}) {
            CAPTURE(reversedChart);
            // MR-04: rotate the circle's chart by pi/2; subtract pi/2 from t.
            const Curve circle{Circle{p.Point(0, 0, 0), 1.5 * p.Scale, reversedChart ? p.Y : p.X,
                                      reversedChart ? -p.X : p.Y}};
            for (const bool inside : {true, false}) {
                const double angle = (inside ? 0.0 : pi / 2) - (reversedChart ? pi / 2 : 0.0);
                CAPTURE(inside, angle);
                const auto result = Intersect(circle, Interval{angle}, plane, band, p.Policy());
                REQUIRE(result.has_value());
                REQUIRE(result->Points.empty());
                if (!inside) {
                    REQUIRE(result->IsEmpty());
                    continue;
                }
                REQUIRE(result->Overlaps.size() == 1);
                const auto& overlap = result->Overlaps.front();
                CHECK(overlap.Range == Interval{angle});
                const Point3 expected = p.Point(1.5, 0, 0);
                CHECK(Distance(Eval(circle, angle), expected) <= p.Policy().GetResolution());
                CHECK(Distance(Eval(plane, Eval(overlap.OnSurface, angle)), expected) <=
                      p.Policy().GetResolution());
            }
        }
    }
}

TEST_CASE("MR PR-247 unequal cone angles remain unresolved",
          "[intersect][metamorphic][backfill][PR-247]")
{
    for (const auto& p : Placements()) {
        INFO("PR-247 scale=" << p.Scale << " axis=" << p.X.X << ',' << p.X.Y << ',' << p.X.Z
                             << " offset=" << p.Offset.X << ',' << p.Offset.Y << ',' << p.Offset.Z);
        Replay("PR-247", p);
        const Surface a{Cone{p.Point(0, 0, 0), 0, pi / 4, p.X, p.Y}};
        const Surface b{Cone{p.Point(0, 0, 0), 0, pi / 4 + 0x1p-42, p.X, p.Y}};
        // Independent premise: unequal angles mean unequal support sets.
        // This is a fail-closed contract regression, not resolved geometry.
        REQUIRE(pi / 4 + 0x1p-42 != pi / 4);
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto result = swapped ? Intersect(b, a, p.Policy()) : Intersect(a, b, p.Policy());
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
}

TEST_CASE("MR PR-244 near-axis sphere cone is not isolated points",
          "[intersect][metamorphic][backfill][PR-244]")
{
    // No translation: adding a large offset could round away the tiny
    // off-axis displacement and would not preserve this fixture's premise.
    for (const double offset : {1e-15, 1e-7, 5e-6}) {
        CAPTURE(offset);
        const Surface cone{Cone{Point3{}, 0, pi / 6, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
        const Surface sphere{Sphere{Point3{offset, 0, 4}, 2, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
        const Box2 box{Interval{0, 2 * pi}, Interval{-8, 8}};
        const Context policy = Context::Default().WithScheduler(Scheduler::Serial);
        const auto a =
            Intersect(sphere, Box2::Everything(), cone, box, policy.GetResolution(), policy);
        const auto b =
            Intersect(cone, box, sphere, Box2::Everything(), policy.GetResolution(), policy);
        REQUIRE(a.has_value() == b.has_value());
        if (!a) {
            CHECK(a.error().GetCode() == ErrorCode::IntersectionFailed);
            CHECK(b.error().GetCode() == a.error().GetCode());
            continue;
        }
        // Tolerant near-axis contact is a curve or unresolved, never two
        // meridian points. Check both supports independently on both sides.
        REQUIRE_FALSE(a->Branches.empty());
        REQUIRE_FALSE(b->Branches.empty());
        CHECK(a->Points.empty());
        CHECK(b->Points.empty());
        CHECK(a->Branches.size() == b->Branches.size());
        for (const auto& branch : a->Branches) {
            CHECK(IntersectTest::PcurveDeviation(branch, sphere, cone) <= policy.GetResolution());
        }
        for (const auto& branch : b->Branches) {
            CHECK(IntersectTest::PcurveDeviation(branch, cone, sphere) <= policy.GetResolution());
        }
    }
}

namespace {

void EllipsoidPoint(const Surface& plane, const Box2& planeBox, const Surface& ellipsoid,
                    const Box2& ellipsoidBox, const Point3& expected, const Context& policy)
{
    for (const bool swapped : {false, true}) {
        CAPTURE(swapped);
        const auto result = swapped ? Intersect(ellipsoid, ellipsoidBox, plane, planeBox, policy)
                                    : Intersect(plane, planeBox, ellipsoid, ellipsoidBox, policy);
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.empty());
        REQUIRE(result->Points.size() == 1);
        const auto& point = result->Points.front();
        const auto uvPlane = swapped ? point.OnB : point.OnA;
        const auto uvEllipsoid = swapped ? point.OnA : point.OnB;
        CHECK(Contains(planeBox.X, uvPlane.X));
        CHECK(Contains(planeBox.Y, uvPlane.Y));
        CHECK(Contains(ellipsoidBox.X, uvEllipsoid.X));
        CHECK(Contains(ellipsoidBox.Y, uvEllipsoid.Y));
        CHECK(Distance(point.Position, expected) <= policy.GetResolution());
        CHECK(Distance(Eval(plane, uvPlane), expected) <= policy.GetResolution());
        CHECK(Distance(Eval(ellipsoid, uvEllipsoid), expected) <= policy.GetResolution());
    }
}

void RoundedTangency(const Placement& p, const bool handed, const double side)
{
    const double a = std::sqrt(0.5);
    const Ellipsoid ellipsoid{Point3{},        3 * p.Scale,      2 * p.Scale, p.Scale,
                              a * (p.X + p.Y), a * (-p.X + p.Y), handed};
    const Point3 pole = Point3{} + side * ellipsoid.C * Axis(ellipsoid);
    const Surface surface{ellipsoid};
    const Surface plane{Plane{pole, p.X, p.Y}};
    EllipsoidPoint(plane, Domain(plane), surface, Domain(surface), pole, p.Policy());
    // One representable step outside the evaluator's pole must be separate.
    // The axis permutation is exact; no translation perturbs a tangent input.
    const Vec3 normal = Cross(p.X, p.Y);
    const double height = Dot(pole - Point3{}, normal);
    const double outside =
        std::nextafter(height, std::copysign(std::numeric_limits<double>::infinity(), height));
    const Surface apart{Plane{Point3{} + outside * normal, p.X, p.Y}};
    for (const bool swapped : {false, true}) {
        CAPTURE(swapped);
        const auto result =
            swapped ? Intersect(surface, apart, p.Policy()) : Intersect(apart, surface, p.Policy());
        REQUIRE(result.has_value());
        CHECK(result->IsEmpty());
    }
}

BSplineCurve EndSegment(const Placement& p, const bool second, const double gap,
                        const bool reversed, const bool remapped)
{
    const double low = remapped ? 2.0 : 0.0;
    const double high = remapped ? 6.0 : 1.0;
    BSplineCurve curve{
        1,
        {low, low, high, high},
        second ? std::vector<Point3>{p.Point(1 + gap, gap, 0), p.Point(2 + gap, gap + 0x1p-23, 0)}
               : std::vector<Point3>{p.Point(0, 0, 0), p.Point(1, 0, 0)},
        {}};
    return reversed ? Reversed(curve) : curve;
}

void SplineEndMeeting(const Placement& p, const bool reversed, const bool remapped,
                      const bool apart)
{
    const double gap = apart ? 0x1p-15 : 0x1p-20;
    const BSplineCurve a = EndSegment(p, false, gap, reversed, remapped);
    const BSplineCurve b = EndSegment(p, true, gap, reversed, remapped);
    const double low = remapped ? 2.0 : 0.0;
    const double high = remapped ? 6.0 : 1.0;
    const double ta = reversed ? low : high;
    const double tb = reversed ? high : low;
    // The two segments extend away from these endpoints, so this is the
    // constrained closest pair, regardless of their unbounded supports.
    const Point3 pa = p.Point(1, 0, 0);
    const Point3 pb = p.Point(1 + gap, gap, 0);
    REQUIRE((Distance(pa, pb) > p.Policy().GetResolution()) == apart);
    for (const bool swapped : {false, true}) {
        CAPTURE(swapped);
        const Curve first{swapped ? b : a};
        const Curve second{swapped ? a : b};
        const auto result = Intersect(first, second, p.Policy());
        REQUIRE(result.has_value());
        REQUIRE(result->Overlaps.empty());
        if (apart) {
            REQUIRE(result->IsEmpty());
            continue;
        }
        REQUIRE(result->Points.size() == 1);
        const auto& point = result->Points.front();
        CHECK(point.ParameterA == (swapped ? tb : ta));
        CHECK(point.ParameterB == (swapped ? ta : tb));
        CHECK(Distance(Eval(first, point.ParameterA), swapped ? pb : pa) <=
              p.Policy().GetResolution());
        CHECK(Distance(Eval(second, point.ParameterB), swapped ? pa : pb) <=
              p.Policy().GetResolution());
        // Space-curve results place the representative on operand A, so
        // exchanging operands does not promise bitwise-identical positions.
        CHECK(Distance(point.Position, swapped ? pb : pa) <= p.Policy().GetResolution());
    }
}

} // namespace

TEST_CASE("MR PR-230 ellipsoid rounded axis tangency", "[intersect][metamorphic][backfill][PR-230]")
{
    for (const auto& p : Placements()) {
        if (p.Offset.X != 0 || p.Offset.Y != 0 || p.Offset.Z != 0) {
            continue; // A rounded translation need not preserve exact tangency.
        }
        CAPTURE(p.Scale, p.X.X, p.X.Y, p.X.Z);
        Replay("PR-230", p);
        for (const bool handed : {false, true}) {
            for (const double side : {-1.0, 1.0}) {
                CAPTURE(handed, side);
                RoundedTangency(p, handed, side);
            }
        }
    }
}

TEST_CASE("MR PR-230-pole ellipsoid longitude boxes",
          "[intersect][metamorphic][backfill][PR-230-pole]")
{
    for (const auto& p : Placements()) {
        CAPTURE(p.Scale, p.X.X, p.X.Y, p.X.Z, p.Offset.X, p.Offset.Y, p.Offset.Z);
        Replay("PR-230-pole", p);
        const Surface surface{
            Ellipsoid{p.Point(0, 0, 0), 3 * p.Scale, 2 * p.Scale, p.Scale, p.X, p.Y}};
        const Box2 planeBox{Interval{-4 * p.Scale, 4 * p.Scale},
                            Interval{-4 * p.Scale, 4 * p.Scale}};
        for (const double side : {-1.0, 1.0}) {
            CAPTURE(side);
            const Point3 pole = p.Point(0, 0, side);
            const Surface plane{Plane{pole, p.X, p.Y}};
            for (const double turn : {0.0, 2 * pi}) {
                CAPTURE(turn);
                const Box2 box{Interval{1 + turn, 2 + turn}, Interval{-pi / 2, pi / 2}};
                EllipsoidPoint(plane, planeBox, surface, box, pole, p.Policy());
                // Latitude exclusion must remain an independently proved miss.
                const Box2 excluded{box.X, Interval{-pi / 4, pi / 4}};
                const auto result = Intersect(plane, planeBox, surface, excluded, p.Policy());
                REQUIRE(result.has_value());
                CHECK(result->IsEmpty());
            }
        }
    }
}

TEST_CASE("MR PR-241 spline boundary minimum", "[intersect][metamorphic][backfill][PR-240][PR-241]")
{
    for (const auto& p : Placements()) {
        CAPTURE(p.Scale, p.X.X, p.X.Y, p.X.Z, p.Offset.X, p.Offset.Y, p.Offset.Z);
        Replay("PR-241", p);
        for (const bool reversed : {false, true}) {
            for (const bool remapped : {false, true}) {
                for (const bool apart : {false, true}) {
                    CAPTURE(reversed, remapped, apart);
                    SplineEndMeeting(p, reversed, remapped, apart);
                }
            }
        }
    }
}
