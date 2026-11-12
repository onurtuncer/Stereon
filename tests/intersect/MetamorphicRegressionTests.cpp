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
            // An unresolved intersection or an uncertified pcurve must fail
            // closed; neither is a successful geometric answer.
            const auto code = a.error().GetCode();
            CHECK(
                (code == ErrorCode::IntersectionFailed || code == ErrorCode::CertificationFailed));
            CHECK(b.error().GetCode() == code);
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

namespace {

/// A proper cyclic permutation of coordinates: exact on every double, so a
/// fixture's premise (a rounding-sized gap, a far origin) survives it.
struct Cyclic {
    int Shift;

    [[nodiscard]] Vec3 operator()(const Vec3& v) const
    {
        const std::array<double, 3> c{v.X, v.Y, v.Z};
        const auto at = [&](const int i) {
            return c[static_cast<std::size_t>((i + 3 - Shift) % 3)];
        };
        return Vec3{at(0), at(1), at(2)};
    }

    [[nodiscard]] Point3 operator()(const Point3& p) const
    {
        return Point3{} + (*this)(p - Point3{});
    }
};

/// How far a point is from the cylinder of radius r about the line through
/// `on` along `axis`, from the closed form, not from the kernel.
double OffCylinder(const Point3& p, const Point3& on, const Vec3& axis, const double r)
{
    return std::abs(Length(Cross(p - on, axis)) / Length(axis) - r);
}

} // namespace

TEST_CASE("MR PR-200 far line origin crossings on the quadric",
          "[intersect][metamorphic][backfill][PR-200]")
{
    // PR #200's counterexample (RapidCheck seed 17734536904155436187): a
    // line about 24 units from a sphere and a cylinder of radius 1/2. Every
    // input is a dyadic rational, so moving the origin back along the line
    // by m directions (MR-04) and permuting the axes (MR-01) are exact. The
    // crossings must lie on each quadric to the point's own rounding, by the
    // closed-form distance, and the parameters must shift by exactly m.
    constexpr double radius = 0.5;
    const Context policy = Context::Default().WithScheduler(Scheduler::Serial);
    for (const int shift : {0, 1, 2}) {
        const Cyclic turn{shift};
        const Point3 center = turn(Point3{-11.0 / 8.0, -11.0 / 8.0, -12.0 / 8.0});
        const Vec3 direction = turn(Vec3{-30.0 / 8.0, -40.0 / 8.0, -39.0 / 8.0});
        const Vec3 axis = Cross(direction, turn(Vec3{-10.0 / 8.0, -20.0 / 8.0, -34.0 / 8.0}));
        const Point3 start = center + turn(Vec3{-0.125, -0.125, -0.1875}) - 3.0 * direction;
        const Surface sphere{Sphere{center, radius, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
        const Surface cylinder{IntersectTest::CylinderAbout(center, axis, radius)};
        for (const bool round : {true, false}) {
            const Surface& surface = round ? sphere : cylinder;
            std::vector<double> base;
            for (const double m : {0.0, 8.0, 64.0}) {
                INFO("PR-200 shift=" << shift << " sphere=" << round << " m=" << m);
                const Line line{start - m * direction, direction};
                const auto result = Intersect(Curve{line}, surface, policy);
                REQUIRE(result.has_value());
                REQUIRE(result->Points.size() == 2);
                for (std::size_t i = 0; i < 2; ++i) {
                    const auto& p = result->Points[i];
                    // The point's own rounding, and the line's: a point on
                    // it is placed no finer than a unit in the last place of
                    // t, times the direction's length.
                    const double eps = std::numeric_limits<double>::epsilon();
                    const double ulps = 32.0 * eps * std::max(1.0, Length(p.Position - Point3{})) +
                                        4.0 * eps * std::abs(p.Parameter) * Length(direction);
                    CHECK(p.Contact == ContactType::Crossing);
                    CHECK(Distance(Eval(line, p.Parameter), p.Position) <= ulps);
                    const double off = round ? std::abs(Distance(p.Position, center) - radius)
                                             : OffCylinder(p.Position, center, axis, radius);
                    CHECK(off <= ulps);
                    if (m == 0.0) {
                        base.push_back(p.Parameter);
                    } else {
                        CHECK(std::abs(p.Parameter - m - base[i]) <= 1e-12 * m);
                    }
                }
            }
        }
    }
}

TEST_CASE("MR PR-207 traced sliver at a box corner is a point",
          "[intersect][metamorphic][backfill][PR-207]")
{
    // PR #207's fixture: the cylinder of radius 2 about z meets the sphere
    // of radius 2 about (1, 0, 0) where z^2 = 4 cos u - 1. The cylinder's
    // box has its corner 1e-13 below that curve at u0 = pi / 6, so the
    // curve crosses the box over about 1e-13 there: a sliver, which must
    // be a point, never an empty result. MR-01 by exact axis permutations
    // (a translation or scale could round the 1e-13 premise away), MR-02 by
    // exchanging the operands. The oracle is the closed-form point and the
    // closed-form distances from both surfaces.
    const double u0 = pi / 6.0;
    const double z0 = std::sqrt(4.0 * std::cos(u0) - 1.0);
    const Box2 corner{Interval{u0, u0 + 1.0}, Interval{z0 - 1e-13, z0 + 1.0}};
    const Context policy = Context::Default().WithScheduler(Scheduler::Serial);
    for (const int shift : {0, 1, 2}) {
        const Cyclic turn{shift};
        const Vec3 x = turn(Vec3{1, 0, 0});
        const Vec3 y = turn(Vec3{0, 1, 0});
        const Vec3 z = turn(Vec3{0, 0, 1});
        const Point3 centre = turn(Point3{1, 0, 0});
        const Surface cylinder{Cylinder{Point3{}, 2.0, x, y, true}};
        const Surface sphere{Sphere{centre, 2.0, x, y, true}};
        const Point3 sliver = Point3{} + 2.0 * std::cos(u0) * x + 2.0 * std::sin(u0) * y + z0 * z;
        for (const bool swapped : {false, true}) {
            INFO("PR-207 shift=" << shift << " swapped=" << swapped);
            const auto result =
                swapped ? Intersect(sphere, Box2::Everything(), cylinder, corner, policy)
                        : Intersect(cylinder, corner, sphere, Box2::Everything(), policy);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE_FALSE(result->IsEmpty());
            CHECK(result->Branches.empty());
            for (const IntersectionPoint& p : result->Points) {
                CHECK(p.Bound <= policy.GetResolution());
                CHECK(Distance(p.Position, sliver) <= p.Bound + 1e-12);
                CHECK(std::abs(Distance(p.Position, centre) - 2.0) <= p.Bound + 1e-12);
                CHECK(OffCylinder(p.Position, Point3{}, z, 2.0) <= p.Bound + 1e-12);
            }
        }
    }
}

TEST_CASE("MR PR-238 rounded pole plane fails closed on ellipsoid and sphere",
          "[intersect][metamorphic][backfill][PR-238]")
{
    // PR #238's shrunk counterexample (CI seed 1794158452): a plane made
    // through the south pole (-2, -2, -2.5) of a round ellipsoid, whose
    // normalised frame misses the pole by a unit in the last place. The
    // section passes within rounding of the chart's pole, where the
    // pcurve's longitude has no rate: the answer must fail closed with
    // IntersectionFailed, never escape as InvalidInput. MR-02 exchanges the
    // operands; MR-04 replaces the ellipsoid with equal semi-axes by the
    // sphere, the same surface in an equivalent chart. No rigid motion or
    // scale: either would re-round the plane and change the premise.
    const Point3 center{-2.0, -2.0, -2.0};
    const double r = 0.5;
    const Vec3 x{1, 0, 0};
    const Vec3 y{0, 1, 0};
    const Plane support = IntersectTest::PlaneWithNormal(
        center + Vec3{-1.0 / 16.0, -0.5, -9.0 / 16.0}, Vec3{-15.0 / 8.0, 0.25, -0.125});
    const Surface plane{support};
    // Independent premise: the pole is off the rounded plane by rounding
    // only.
    const Point3 pole = center - r * Cross(x, y);
    REQUIRE(std::abs(Dot(pole - support.Origin, Cross(support.XAxis, support.YAxis))) <= 1e-15);
    const Context policy = Context::Default().WithScheduler(Scheduler::Serial);
    for (const bool round : {false, true}) {
        const Surface other =
            round ? Surface{Sphere{center, r, x, y}} : Surface{Ellipsoid{center, r, r, r, x, y}};
        for (const bool swapped : {false, true}) {
            INFO("PR-238 sphere=" << round << " swapped=" << swapped);
            const auto result =
                swapped ? Intersect(other, plane, policy) : Intersect(plane, other, policy);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
}

namespace {

/// A placement without translation: exact axis permutations and powers of
/// two, which keep a rounding-sized premise (a tangency on a seam) intact.
[[nodiscard]] bool Untranslated(const Placement& p)
{
    return p.Offset.X == 0 && p.Offset.Y == 0 && p.Offset.Z == 0;
}

[[nodiscard]] std::string Where(const Placement& p)
{
    return "scale=" + std::to_string(p.Scale) + " axis=" + std::to_string(p.X.X) + ',' +
           std::to_string(p.X.Y) + ',' + std::to_string(p.X.Z) +
           " offset=" + std::to_string(p.Offset.X) + ',' + std::to_string(p.Offset.Y) + ',' +
           std::to_string(p.Offset.Z);
}

/// The overlap between two heights' requested intervals, from the closed
/// form: their intersection.
void CheckClippedHeights(const SurfaceSurfaceIntersection& r, const Interval& expected,
                         const double scale)
{
    REQUIRE(r.Overlaps.size() == 1);
    const auto& o = r.Overlaps.front();
    CHECK(std::abs(o.OnA.Y.GetLower() - expected.GetLower()) <= 1e-14 * scale);
    CHECK(std::abs(o.OnA.Y.GetUpper() - expected.GetUpper()) <= 1e-14 * scale);
    CHECK(std::abs(o.OnB.Y.GetLower() - expected.GetLower()) <= 1e-14 * scale);
    CHECK(std::abs(o.OnB.Y.GetUpper() - expected.GetUpper()) <= 1e-14 * scale);
}

} // namespace

TEST_CASE("MR PR-112 coincident cylinders clipped to the requested domains",
          "[intersect][metamorphic][backfill][PR-112]")
{
    // PR #112: coincident cylinders returned the whole requested boxes as
    // their overlap, whatever the other's box, even boxes apart. The same
    // cylinder twice, under every placement (MR-01, MR-03: heights and the
    // resolution scale, angles do not) and in both orders (MR-02). The
    // oracle is the closed-form intersection of the requested intervals.
    for (const auto& p : Placements()) {
        INFO("PR-112 " << Where(p));
        Replay("PR-112", p);
        const Surface cylinder{Cylinder{p.Point(0, 0, 0), 3 * p.Scale, p.X, p.Y, true}};
        const Interval turn{0.0, 2 * pi};
        const Box2 first{turn, Interval{0, p.Scale}};
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto meet = [&](const Box2& a, const Box2& b) {
                return swapped ? Intersect(cylinder, b, cylinder, a, p.Policy().GetResolution(),
                                           p.Policy())
                               : Intersect(cylinder, a, cylinder, b, p.Policy().GetResolution(),
                                           p.Policy());
            };
            const auto apart = meet(first, Box2{turn, Interval{2 * p.Scale, 3 * p.Scale}});
            REQUIRE(apart.has_value());
            CHECK(apart->IsEmpty());
            const auto partial = meet(first, Box2{turn, Interval{0.5 * p.Scale, 2 * p.Scale}});
            REQUIRE(partial.has_value());
            CheckClippedHeights(*partial, Interval{0.5 * p.Scale, p.Scale}, p.Scale);
            CHECK_FALSE(partial->Overlaps.front().FullDomain);
            // Angular arcs [0, 1] and [2, 3] are apart, modulo a turn too.
            const Interval height{0, p.Scale};
            const auto arcs = meet(Box2{Interval{0, 1}, height}, Box2{Interval{2, 3}, height});
            REQUIRE(arcs.has_value());
            CHECK(arcs->IsEmpty());
        }
    }
}

TEST_CASE("MR PR-112-budget plane cylinder pcurve work is bounded",
          "[intersect][metamorphic][backfill][PR-112-budget]")
{
    // PR #112: the plane-cylinder ellipse pcurve converted and allocated a
    // sample count that a tiny resolution made unbounded. With the
    // resolution 1e-20 times the scale, the sample spacing, which goes as
    // (resolution / radius)^(1/4) times the radius, the count, angle over
    // spacing, is the same at every scale (MR-03) and placement (MR-01), in
    // both orders (MR-02): about 2e5 samples, beyond the budget of 2^16, so
    // the contract is CertificationFailed naming the sample budget.
    for (const auto& p : Placements()) {
        INFO("PR-112-budget " << Where(p));
        const Context policy = p.Policy().WithResolution(1e-20 * p.Scale);
        const Surface oblique{
            IntersectTest::PlaneWithNormal(p.Point(0, 0, 0), p.Y + Cross(p.X, p.Y))};
        const Surface cylinder{Cylinder{p.Point(0, 0, 0), 5 * p.Scale, p.X, p.Y, true}};
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto result = swapped ? Intersect(cylinder, oblique, policy)
                                        : Intersect(oblique, cylinder, policy);
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error().GetCode() == ErrorCode::CertificationFailed);
            CHECK(std::string(result.error().GetMessage()).find("sample budget") !=
                  std::string::npos);
        }
    }
}

TEST_CASE("MR PR-120 skew cylinder angular sectors of B",
          "[intersect][metamorphic][backfill][PR-120]")
{
    // PR #120's cases, repaired in #118 (3a2ba89): B's angular sector was
    // not clipped. A is the cylinder of radius 3 about z; B, of radius 2,
    // runs along x through (0, 4, 0), its angle from y towards z. They meet
    // only where cos(u_B) <= -1/2, the closed form. So the sector [-2, 2]
    // (wider than pi, a turn on too) and the singleton {0} must be proved
    // empty, and the sector [2.5, 3.8], inside [2 pi / 3, 4 pi / 3], must
    // not. Every placement (MR-01, MR-03) and both orders (MR-02).
    for (const auto& p : Placements()) {
        INFO("PR-120 " << Where(p));
        Replay("PR-120", p);
        const Vec3 z = Cross(p.X, p.Y);
        const Surface a{Cylinder{p.Point(0, 0, 0), 3 * p.Scale, p.X, p.Y, true}};
        const Surface b{Cylinder{p.Point(0, 4, 0), 2 * p.Scale, p.Y, z, true}};
        const auto meet = [&](const Box2& patch, const bool swapped) {
            return swapped ? Intersect(b, patch, a, Domain(a), p.Policy())
                           : Intersect(a, Domain(a), b, patch, p.Policy());
        };
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            for (const double shift : {0.0, 2 * pi}) {
                const auto wide =
                    meet(Box2{Interval{-2 + shift, 2 + shift}, Interval::Everything()}, swapped);
                INFO(IntersectTest::ErrorText(wide));
                REQUIRE(wide.has_value());
                CHECK(wide->IsEmpty());
            }
            const auto single = meet(Box2{Interval{0.0}, Interval::Everything()}, swapped);
            REQUIRE(single.has_value());
            CHECK(single->IsEmpty());
            const auto inside = meet(Box2{Interval{2.5, 3.8}, Interval::Everything()}, swapped);
            INFO(IntersectTest::ErrorText(inside));
            REQUIRE(inside.has_value());
            CHECK_FALSE(inside->IsEmpty());
        }
    }
}

TEST_CASE("MR PR-159 skew cylinders tangent on both seams meet at a node",
          "[intersect][metamorphic][backfill][PR-159]")
{
    // PR #159 (issue #155): radius 3 inside radius 5 with a tilted axis,
    // touching at (5, 0, 0), on the seam of each, was refused. Exact axis
    // permutations and powers of two (MR-01, MR-03; no translation, which
    // could round the tangency's premise), both orders (MR-02). The oracle:
    // one node at the closed-form point of tangency with four branch ends,
    // and every branch on both cylinders by the closed-form distance from
    // each axis.
    for (const auto& p : Placements()) {
        if (!Untranslated(p)) {
            continue;
        }
        INFO("PR-159 " << Where(p));
        Replay("PR-159", p);
        const Vec3 z = Cross(p.X, p.Y);
        const Vec3 tilted = 0.6 * p.Y + 0.8 * z;
        const Surface larger{Cylinder{p.Point(0, 0, 0), 5 * p.Scale, p.X, p.Y, true}};
        const Surface smaller{Cylinder{p.Point(2, 0, 0), 3 * p.Scale, p.X, tilted, true}};
        const Point3 node = p.Point(5, 0, 0);
        const double reach = 1e-7 * p.Scale;
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto result = swapped ? Intersect(smaller, larger, p.Policy())
                                        : Intersect(larger, smaller, p.Policy());
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            CHECK(IntersectTest::ConnectivityOf(*result, reach).Nodes == 1);
            std::size_t ends = 0;
            for (const IntersectionBranch& branch : result->Branches) {
                for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
                    ends += Distance(Eval(branch.Curve3d, t), node) <= reach ? 1 : 0;
                }
                for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
                    const Point3 q = Eval(branch.Curve3d, t);
                    const double slack = branch.Bound + p.Policy().GetResolution();
                    CHECK(OffCylinder(q, p.Point(0, 0, 0), z, 5 * p.Scale) <= slack);
                    CHECK(OffCylinder(q, p.Point(2, 0, 0), Cross(p.X, tilted), 3 * p.Scale) <=
                          slack);
                }
            }
            CHECK(ends == 4);
        }
    }
}

TEST_CASE("MR PR-236 concentric spheres within the tolerance are declared",
          "[intersect][metamorphic][backfill][PR-236]")
{
    // PR #236: concentric spheres whose radii differ by less than the pair
    // tolerance were an Exact overlap at distance 0. ADR-0029 section 6: a
    // coincidence declared within the tolerance is Tolerant and records how
    // far apart the surfaces are, here |R_b - R_a|, the closed form. Every
    // placement with the tolerance scaled (MR-01, MR-03), both orders
    // (MR-02); beyond the tolerance they are apart.
    for (const auto& p : Placements()) {
        INFO("PR-236 " << Where(p));
        Replay("PR-236", p);
        const double ra = 2 * p.Scale;
        const double rb = 2.0005 * p.Scale;
        const Surface a{Sphere{p.Point(0, 0, 0), ra, p.X, p.Y, true}};
        const Surface b{Sphere{p.Point(0, 0, 0), rb, p.X, p.Y, true}};
        for (const bool swapped : {false, true}) {
            CAPTURE(swapped);
            const auto meet = [&](const double tolerance) {
                return swapped ? Intersect(b, Box2::Everything(), a, Box2::Everything(), tolerance,
                                           p.Policy())
                               : Intersect(a, Box2::Everything(), b, Box2::Everything(), tolerance,
                                           p.Policy());
            };
            const auto declared = meet(1e-3 * p.Scale);
            REQUIRE(declared.has_value());
            REQUIRE(declared->Overlaps.size() == 1);
            CHECK(declared->Overlaps.front().Evidence == CoincidenceEvidence::Tolerant);
            CHECK(std::abs(declared->Overlaps.front().Distance - (rb - ra)) <= 1e-15 * p.Scale);
            const auto apart = meet(1e-4 * p.Scale);
            REQUIRE(apart.has_value());
            CHECK(apart->IsEmpty());
        }
    }
}

TEST_CASE("MR PR-236-tangent sphere on a cylinder's axis within the tolerance",
          "[intersect][metamorphic][backfill][PR-236-tangent]")
{
    // PR #236: a sphere centred on a cylinder's axis, radii 5e-4 apart under
    // a pair tolerance of 1e-3, is declared tangent along a circle; its
    // certification against the resolution failed the whole intersection.
    // Every placement with lengths and the tolerance scaled (MR-01, MR-03),
    // both orders (MR-02). Oracle: one tangential curve, certified, its
    // bound at least 0.8 of the closed-form gap |R_c - R_s| and within the
    // tolerance.
    for (const auto& p : Placements()) {
        INFO("PR-236-tangent " << Where(p));
        Replay("PR-236-tangent", p);
        const Surface sphere{Sphere{p.Point(0, 0, 0), 2 * p.Scale, p.X, p.Y, true}};
        for (const double radius : {2.0 - 5e-4, 2.0 + 5e-4}) {
            const Surface cylinder{Cylinder{p.Point(0, 0, 0), radius * p.Scale, p.X, p.Y, true}};
            const double gap = std::abs(radius - 2.0) * p.Scale;
            const double tolerance = 1e-3 * p.Scale;
            for (const bool swapped : {false, true}) {
                CAPTURE(radius, swapped);
                const auto touch = swapped ? Intersect(cylinder, Box2::Everything(), sphere,
                                                       Box2::Everything(), tolerance, p.Policy())
                                           : Intersect(sphere, Box2::Everything(), cylinder,
                                                       Box2::Everything(), tolerance, p.Policy());
                INFO(IntersectTest::ErrorText(touch));
                REQUIRE(touch.has_value());
                CHECK(IntersectTest::ConnectivityOf(*touch, 1e-7 * p.Scale).Curves == 1);
                for (const IntersectionBranch& arc : touch->Branches) {
                    CHECK(arc.Contact == ContactType::Tangential);
                    CHECK(arc.IsCertified());
                    CHECK(arc.Bound >= 0.8 * gap);
                    CHECK(arc.Bound <= tolerance);
                }
            }
        }
    }
}

TEST_CASE("MR PR-263-translated-boundary spline end meeting away from the origin",
          "[intersect][metamorphic][backfill][PR-263-translated-boundary]")
{
    // PR #263: the PR-241 fixture failed natively at scale 1/4 translated
    // by (8, -4, 2) (Actions run 43052745077, job 130633998895), resolved
    // at the origin: the pair search inverted an ill-conditioned Krawczyk
    // matrix where no stationary point could exist. The translated
    // placements, every scale (MR-01, MR-03), with reversal and knot maps
    // (MR-04) and both orders (MR-02), as SplineEndMeeting checks.
    for (const auto& p : Placements()) {
        if (Untranslated(p)) {
            continue;
        }
        CAPTURE(p.Scale, p.X.X, p.X.Y, p.X.Z);
        Replay("PR-263-translated-boundary", p);
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

// --- PR-194: one ledger entry, one test per root cause ---------------------

namespace {

/// The proper cyclic permutations of the axes and the identity: exact on
/// every double.
const std::array<std::array<Vec3, 3>, 3> Turns{{
    {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}},
    {Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}},
    {Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}},
}};

Vec3 Turned(const std::array<Vec3, 3>& t, const Vec3& v)
{
    return v.X * t[0] + v.Y * t[1] + v.Z * t[2];
}

Point3 Turned(const std::array<Vec3, 3>& t, const Point3& p)
{
    return Point3{} + Turned(t, p - Point3{});
}

/// A Bezier curve of degree poles - 1 on [0, 1].
BSplineCurve Bezier3(std::vector<Point3> poles)
{
    const std::size_t n = poles.size();
    std::vector<double> knots(n, 0.0);
    knots.resize(2 * n, 1.0);
    return BSplineCurve{static_cast<int>(n) - 1, std::move(knots), std::move(poles), {}};
}

BSplineCurve Reversed3(BSplineCurve c)
{
    const double sum = c.Knots.front() + c.Knots.back();
    std::vector<double> knots;
    for (auto it = c.Knots.rbegin(); it != c.Knots.rend(); ++it) {
        knots.push_back(sum - *it);
    }
    c.Knots = std::move(knots);
    std::ranges::reverse(c.Poles);
    return c;
}

} // namespace

TEST_CASE("MR PR-194 rotated bounded plane coincidence", "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (PlanePlane): a rotated, bounded coincidence was encoded as
    // the AABB of its image, which invents overlaps. The unit square of the
    // xy plane, and the square of half-side 0.35 turned 45 degrees about a
    // centre (c, c, 0): the diamond |x - c| + |y - c| <= 0.35 sqrt 2. For
    // c = 1.4 its least x + y is 2.8 - 0.495 > 2, the unit square's
    // greatest, so they are apart (closed form) though the AABBs overlap;
    // for c = 0.7 they truly overlap partly, which a rectangle cannot
    // represent: NotImplemented. MR-01 by exact axis permutations and the
    // translation (8, -4, 2), MR-02 by exchanging the planes.
    const double h = std::sqrt(0.5);
    const Box2 a{Interval{0.0, 1.0}, Interval{0.0, 1.0}};
    const Box2 b{Interval{-0.35, 0.35}, Interval{-0.35, 0.35}};
    for (const auto& t : Turns) {
        for (const Vec3 offset : {Vec3{}, Vec3{8, -4, 2}}) {
            const Surface xy{Plane{Turned(t, Point3{}) + offset, t[0], t[1]}};
            for (const double centre : {1.4, 0.7}) {
                const Surface turned{Plane{Turned(t, Point3{centre, centre, 0.0}) + offset,
                                           Turned(t, Vec3{h, h, 0.0}), Turned(t, Vec3{-h, h, 0.0})}};
                for (const bool swapped : {false, true}) {
                    INFO("PR-194 plane-plane centre=" << centre << " offset=" << offset.X
                                                      << " swapped=" << swapped);
                    const auto r = swapped ? Intersect(turned, b, xy, a, Context::Default())
                                           : Intersect(xy, a, turned, b, Context::Default());
                    if (centre > 1.0) {
                        REQUIRE(r.has_value());
                        CHECK(r->IsEmpty());
                    } else {
                        REQUIRE_FALSE(r.has_value());
                        CHECK(r.error().GetCode() == ErrorCode::NotImplemented);
                    }
                }
            }
        }
    }
}

TEST_CASE("MR PR-194 near-parallel planes underflowed cross norm",
          "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (PlanePlane): the squared norm of an exactly nonzero cross
    // product of normals underflowed, and normalising it gave no direction.
    // z = 0 and the plane through the origin spanned by x and (0, 1, 2^-600)
    // meet along the x axis (closed form): one branch, finite, along x,
    // through the origin. MR-01 by exact axis permutations (a translation
    // would not change the tilt), MR-02 by exchanging the planes.
    const double tilt = std::ldexp(1.0, -600);
    for (const auto& t : Turns) {
        const Surface xy{Plane{Point3{}, t[0], t[1]}};
        const Surface tilted{Plane{Point3{}, t[0], Turned(t, Vec3{0.0, 1.0, tilt})}};
        for (const bool swapped : {false, true}) {
            INFO("PR-194 underflow turn=" << t[0].X << t[0].Y << t[0].Z << " swapped=" << swapped);
            const auto r = swapped ? Intersect(tilted, xy, Context::Default())
                                   : Intersect(xy, tilted, Context::Default());
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            REQUIRE(r->Branches.size() == 1);
            const auto* line = r->Branches.front().Curve3d.Target<Line>();
            REQUIRE(line != nullptr);
            const Vec3 d = line->Direction;
            REQUIRE(std::isfinite(d.X));
            REQUIRE(std::isfinite(d.Y));
            REQUIRE(std::isfinite(d.Z));
            CHECK(Length(Cross(d, t[0])) <= 1e-12 * Length(d));
            CHECK(Length(Cross(line->Origin - Point3{}, t[0])) <= 1e-12);
        }
    }
}

TEST_CASE("MR PR-194 bounded plane box gets no uncut conic branch",
          "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (PlaneCylinder): a section was kept whole whatever the
    // plane's requested box. The cylinder of radius 5 about z; the plane
    // through the origin with normal (0, 1, 1), whose section is an ellipse
    // with semi-axes 5 sqrt 2 and 5 (closed form). Over the plane box
    // [-0.1, 0.1]^2, deep inside the ellipse, there is nothing; over
    // [-1, 1] x [4, 6], which the ellipse crosses, either a certified piece
    // whose every point lies in the box, or a refusal: never the whole
    // ellipse. MR-01 by exact axis permutations, MR-02 by exchanging the
    // surfaces.
    for (const auto& t : Turns) {
        const Surface cylinder{Cylinder{Point3{}, 5.0, t[0], t[1], true}};
        const Plane support =
            IntersectTest::PlaneWithNormal(Point3{}, Turned(t, Vec3{0.0, 1.0, 1.0}));
        const Surface plane{support};
        for (const bool partial : {false, true}) {
            const Box2 box = partial ? Box2{Interval{-1.0, 1.0}, Interval{4.0, 6.0}}
                                     : Box2{Interval{-0.1, 0.1}, Interval{-0.1, 0.1}};
            for (const bool swapped : {false, true}) {
                INFO("PR-194 plane-cylinder partial=" << partial << " swapped=" << swapped);
                const auto r = swapped
                                   ? Intersect(cylinder, Box2::Everything(), plane, box, Context::Default())
                                   : Intersect(plane, box, cylinder, Box2::Everything(), Context::Default());
                if (!partial) {
                    REQUIRE(r.has_value());
                    CHECK(r->IsEmpty());
                    continue;
                }
                if (!r) {
                    CHECK(r.error().GetCode() == ErrorCode::NotImplemented);
                    continue;
                }
                for (const IntersectionBranch& branch : r->Branches) {
                    for (const double s : IntersectTest::SamplesOf(branch.Range, 64)) {
                        const Point3 p = Eval(branch.Curve3d, s);
                        const double u = Dot(p - support.Origin, support.XAxis);
                        const double v = Dot(p - support.Origin, support.YAxis);
                        CHECK(Contains(Interval{-1.0 - 1e-9, 1.0 + 1e-9}, u));
                        CHECK(Contains(Interval{4.0 - 1e-9, 6.0 + 1e-9}, v));
                    }
                }
            }
        }
    }
}

TEST_CASE("MR PR-194 undecided spline root fails closed", "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (SplineSurface): a root enclosure left undecided (a touch
    // or a near miss, or a root on a knot) was reported as a contact. The
    // Bezier (0,0,1), (1,0,-1), (2,0,1) is z = (1 - 2s)^2: it touches z = 0
    // at s = 1/2, a double root (closed form); the polyline through
    // (0,0,-1), (1,0,0), (2,0,1) crosses z = 0 exactly at its knot. Both
    // must fail with IntersectionFailed. MR-01 by exact axis permutations,
    // MR-04 by reversing the curve (knots reflected, poles reversed).
    for (const auto& t : Turns) {
        const Surface ground{Plane{Point3{}, t[0], t[1]}};
        const auto at = [&](const double x, const double z) { return Turned(t, Point3{x, 0.0, z}); };
        const BSplineCurve touch = Bezier3({at(0, 1), at(1, -1), at(2, 1)});
        const BSplineCurve knotted{1, {0.0, 0.0, 1.0, 2.0, 2.0}, {at(0, -1), at(1, 0), at(2, 1)}, {}};
        for (const bool reversed : {false, true}) {
            for (const BSplineCurve& curve : {touch, knotted}) {
                INFO("PR-194 spline-surface reversed=" << reversed << " degree=" << curve.Degree);
                const Curve c{reversed ? Reversed3(curve) : curve};
                const auto r = Intersect(c, Domain(c), ground, Domain(ground), Context::Default());
                REQUIRE_FALSE(r.has_value());
                CHECK(r.error().GetCode() == ErrorCode::IntersectionFailed);
            }
        }
    }
}

TEST_CASE("MR PR-194 stale geometry is refused by JVP and VJP",
          "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (Linearization): a result computed for other geometry was
    // differentiated as if current. The line x = y = 0 meets z = 0 at the
    // origin; the plane moved 0.01 along its normal no longer holds that
    // point (closed form: the stored point is 0.01 off it), so JVP and VJP
    // must refuse with InvalidInput, and accept the current pair. MR-01 by
    // exact axis permutations and the translation (8, -4, 2), MR-03 by the
    // scales 1/4 and 4 (the move scales with the geometry).
    for (const auto& t : Turns) {
        for (const Vec3 offset : {Vec3{}, Vec3{8, -4, 2}}) {
            for (const double scale : {0.25, 1.0, 4.0}) {
                INFO("PR-194 stale offset=" << offset.X << " scale=" << scale);
                const Point3 origin = Point3{} + offset;
                const Vec3 normal = Cross(t[0], t[1]);
                const Plane ground{origin, t[0], t[1]};
                const Plane moved{origin + 0.01 * scale * normal, t[0], t[1]};
                const Line up{origin - scale * normal, scale * normal};
                const auto at = Intersect(Curve{up}, Surface{ground}, Context::Default());
                REQUIRE(at.has_value());
                REQUIRE(at->Points.size() == 1);
                SensitivityContext sensitivity{ParameterSet{}};
                const auto stale = IntersectJvp(up, {}, Surface{moved}, {}, *at, sensitivity);
                REQUIRE_FALSE(stale.has_value());
                CHECK(stale.error().GetCode() == ErrorCode::InvalidInput);
                const auto staleVjp = IntersectVjp(up, {}, Surface{moved}, {}, *at, sensitivity);
                REQUIRE_FALSE(staleVjp.has_value());
                CHECK(staleVjp.error().GetCode() == ErrorCode::InvalidInput);
                const auto current = IntersectJvp(up, {}, Surface{ground}, {}, *at, sensitivity);
                INFO(IntersectTest::ErrorText(current));
                CHECK(current.has_value());
            }
        }
    }
}

TEST_CASE("MR PR-194 NURBS hull applies the factory bounds",
          "[intersect][metamorphic][backfill][PR-194]")
{
    // Root cause (GenericSSI): a B-spline surface built directly, bypassing
    // MakeBSplineSurface, reached the convex-hull test with a weight beyond
    // the factory's bound [2^-256, 2^256], where the hull theorem's premise
    // fails. A bilinear patch with weight 2^257 must be refused with
    // InvalidInput (the factory's own verdict on the same data, checked
    // here). MR-01 by exact axis permutations of the poles, MR-02 by
    // exchanging the surfaces.
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const Box2 middle{Interval{0.25, 0.75}, Interval{0.25, 0.75}};
    for (const auto& t : Turns) {
        const auto patch = [&](const double z) {
            return std::vector<Point3>{Turned(t, Point3{0, 0, z}), Turned(t, Point3{0, 1, z}),
                                       Turned(t, Point3{1, 0, z}), Turned(t, Point3{1, 1, z})};
        };
        const std::vector<double> heavy{0x1p257, 3.0, 2.0, 4.0};
        REQUIRE_FALSE(MakeBSplineSurface(1, 1, knots, knots, patch(0.0), heavy).has_value());
        const BSplineSurface invalid{1, 1, knots, knots, patch(0.0), heavy};
        const auto valid = MakeBSplineSurface(1, 1, knots, knots, patch(100.0));
        REQUIRE(valid.has_value());
        for (const bool swapped : {false, true}) {
            INFO("PR-194 nurbs-hull swapped=" << swapped);
            const auto r = swapped
                               ? Intersect(Surface{*valid}, middle, Surface{invalid}, middle, 0.0, Context::Default())
                               : Intersect(Surface{invalid}, middle, Surface{*valid}, middle, 0.0, Context::Default());
            REQUIRE_FALSE(r.has_value());
            CHECK(r.error().GetCode() == ErrorCode::InvalidInput);
        }
    }
}
