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
    // each axis. Since #289 a refused parametrisation falls back to the
    // other cylinder's, which would hide the fault: the canonical one, by
    // the larger radius, must answer on its own, the meeting logged once.
    //
    // The fault was a seam crossing found a hair inside a branch's end: the
    // piece cut off there, far shorter than the resolution, was refused, or
    // would have been published as a branch of its own. Since #322 the seam
    // isolation puts the exact tangency's crossing on the end itself, where
    // nothing is cut, so the smaller cylinder's frame is also turned about
    // its axis by a power of two, either way. That leaves the tangency as
    // it is, exactly (the axis stays square to the common normal, the
    // offset and the radii keep their values), and moves its seam off the
    // node by 3 * turn * scale, 9e-8 and 7.2e-7 of a millimetre at scale 1
    // and under the resolution at every scale. The seam then crosses two
    // branches that far inside their ends, and the answer must not change:
    // the same node and four ends, and no branch that stays within the
    // resolution of its start.
    for (const auto& p : Placements()) {
        if (!Untranslated(p)) {
            continue;
        }
        INFO("PR-159 " << Where(p));
        Replay("PR-159", p);
        const Vec3 z = Cross(p.X, p.Y);
        const Vec3 tilted = 0.6 * p.Y + 0.8 * z;
        const Surface larger{Cylinder{p.Point(0, 0, 0), 5 * p.Scale, p.X, p.Y, true}};
        const Point3 node = p.Point(5, 0, 0);
        const double reach = 1e-7 * p.Scale;
        for (const double turn : {0.0, 0x1p-22, -0x1p-22, 0x1p-25, -0x1p-25}) {
            CAPTURE(turn);
            const Surface smaller{Cylinder{p.Point(2, 0, 0), 3 * p.Scale, p.X + turn * tilted,
                                           tilted - turn * p.X, true}};
            for (const bool swapped : {false, true}) {
                CAPTURE(swapped);
                std::vector<std::string> log;
                const Context policy = p.Policy().WithLogSink(
                    [&log](const LogRecord& record) { log.emplace_back(record.Message); });
                const auto result = swapped ? Intersect(smaller, larger, policy)
                                            : Intersect(larger, smaller, policy);
                INFO(IntersectTest::ErrorText(result));
                REQUIRE(result.has_value());
                const auto meetings = std::ranges::count_if(log, [](const std::string& line) {
                    return line.find("internally tangent") != std::string::npos;
                });
                CHECK(meetings == 1);
                CHECK(IntersectTest::ConnectivityOf(*result, reach).Nodes == 1);
                std::size_t ends = 0;
                for (const IntersectionBranch& branch : result->Branches) {
                    for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
                        if (Distance(Eval(branch.Curve3d, t), node) <= reach) {
                            ++ends;
                        }
                    }
                    const Point3 start = Eval(branch.Curve3d, branch.Range.GetLower());
                    double extent = 0.0;
                    for (const double t : IntersectTest::SamplesOf(branch.Range, 16)) {
                        const Point3 q = Eval(branch.Curve3d, t);
                        extent = std::max(extent, Distance(q, start));
                        const double slack = branch.Bound + p.Policy().GetResolution();
                        CHECK(OffCylinder(q, p.Point(0, 0, 0), z, 5 * p.Scale) <= slack);
                        CHECK(OffCylinder(q, p.Point(2, 0, 0), Cross(p.X, tilted), 3 * p.Scale) <=
                              slack);
                    }
                    CHECK(extent > p.Policy().GetResolution());
                }
                CHECK(ends == 4);
            }
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

TEST_CASE("MR PR-194 rotated bounded plane coincidence",
          "[intersect][metamorphic][backfill][PR-194]")
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
                                           Turned(t, Vec3{h, h, 0.0}),
                                           Turned(t, Vec3{-h, h, 0.0})}};
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
                const auto r =
                    swapped
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

TEST_CASE("MR PR-194 undecided spline root fails closed",
          "[intersect][metamorphic][backfill][PR-194]")
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
        const auto at = [&](const double x, const double z) {
            return Turned(t, Point3{x, 0.0, z});
        };
        const BSplineCurve touch = Bezier3({at(0, 1), at(1, -1), at(2, 1)});
        const BSplineCurve knotted{
            1, {0.0, 0.0, 1.0, 2.0, 2.0}, {at(0, -1), at(1, 0), at(2, 1)}, {}};
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
            const auto r = swapped ? Intersect(Surface{*valid}, middle, Surface{invalid}, middle,
                                               0.0, Context::Default())
                                   : Intersect(Surface{invalid}, middle, Surface{*valid}, middle,
                                               0.0, Context::Default());
            REQUIRE_FALSE(r.has_value());
            CHECK(r.error().GetCode() == ErrorCode::InvalidInput);
        }
    }
}

namespace {

/// The biquadratic z = alpha x'^2 + beta y'^2 + lift over [low, low +
/// width]^2, x' = x - 1/2, y' = y - 1/2, its poles row by row (a row per u,
/// along x): the Bernstein coefficients of (a + width t)^2, a = low - 1/2,
/// are a^2, a^2 + a width and (a + width)^2. Placed by an exact turn of the
/// axes, a power-of-two scale and a translation; u reversed when asked.
BSplineSurface PlacedQuadric(const double alpha, const double beta, const double lift,
                             const double low, const double width, const std::array<Vec3, 3>& turn,
                             const double scale, const Point3& shift, const bool reversed)
{
    const double a = low - 0.5;
    const std::array<double, 3> square{a * a, a * a + a * width, (a + width) * (a + width)};
    std::vector<Point3> poles;
    for (std::size_t r = 0; r < 3; ++r) {
        const std::size_t i = reversed ? 2 - r : r;
        for (std::size_t j = 0; j < 3; ++j) {
            const Point3 local{low + static_cast<double>(i) * width / 2.0,
                               low + static_cast<double>(j) * width / 2.0,
                               lift + alpha * square.at(i) + beta * square.at(j)};
            poles.push_back(shift + scale * Turned(turn, local - Point3{}));
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// The local point of a placed one: the inverse of the turn, scale and shift.
Point3 Unplaced(const std::array<Vec3, 3>& turn, const double scale, const Point3& shift,
                const Point3& p)
{
    const Vec3 w = (p - shift) / scale;
    return Point3{Dot(w, turn[0]), Dot(w, turn[1]), Dot(w, turn[2])};
}

} // namespace

TEST_CASE("MR SSI-translated-boundary far hyperbola branches end on the domain's sides",
          "[intersect][metamorphic][backfill][SSI-translated-boundary]")
{
    // Root causes (GenericSSI, repaired in #279): far from the origin the
    // patches' interval poles widen with the coordinates while Krawczyk's
    // acceptance is measured in parameters, and a root just outside the
    // other patch's cell left an edge undecided; both refused a pair that is
    // answered at the origin. The fixture of GenericSSITranslationTests:
    // 2x'^2 - 3y'^2 against the bowl x'^2 + y'^2 + 1/100 over a wider
    // square meet where x'^2 - 4y'^2 = 1/100: two branches across the unit
    // square, x' = +-sqrt(1/100 + 4y'^2), each with both ends on one side,
    // x = 0 or x = 1, at |y'| = sqrt(0.06). MR-01: translations far from the
    // origin and cyclic turns of the axes, exact; MR-02: both orders; MR-03:
    // a scale of 4 with the tolerance and the context; MR-04: the first
    // surface's u reversed (u' = 1 - u).
    struct Run {
        Point3 Shift;
        std::size_t Turn;
        double Scale;
        bool Reversed;
    };
    const std::array<Run, 6> runs{{
        {Point3{64.0, 0.0, 0.0}, 0, 1.0, false},
        {Point3{64.0, 0.0, 0.0}, 1, 1.0, false},
        {Point3{480000.0, -300000.0, 210000.0}, 0, 1.0, false},
        {Point3{480000.0, -300000.0, 210000.0}, 2, 1.0, false},
        {Point3{256.0, 0.0, 0.0}, 0, 4.0, false},
        {Point3{64.0, 0.0, 0.0}, 0, 1.0, true},
    }};
    const double half = std::sqrt(0.06);
    for (const Run& run : runs) {
        const auto& turn = Turns.at(run.Turn);
        const Context policy =
            Context::Default()
                .WithResolution(Context::DefaultResolution * run.Scale)
                .WithToleranceCeiling(Context::DefaultToleranceCeiling * run.Scale);
        const double tolerance = 1e-6 * run.Scale;
        const Surface a{
            PlacedQuadric(2.0, -3.0, 0.0, 0.0, 1.0, turn, run.Scale, run.Shift, run.Reversed)};
        const Surface b{
            PlacedQuadric(1.0, 1.0, 0.01, -0.5, 2.0, turn, run.Scale, run.Shift, false)};
        for (const bool swapped : {false, true}) {
            INFO("SSI-translated-boundary shift ("
                 << run.Shift.X << ", " << run.Shift.Y << ", " << run.Shift.Z << ") turn "
                 << run.Turn << " scale " << run.Scale << " reversed " << run.Reversed
                 << " swapped " << swapped);
            const auto r =
                swapped
                    ? Intersect(b, Box2::Everything(), a, Box2::Everything(), tolerance, policy)
                    : Intersect(a, Box2::Everything(), b, Box2::Everything(), tolerance, policy);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            REQUIRE(r->Overlaps.empty());
            REQUIRE(r->Points.empty());
            REQUIRE(r->Branches.size() == 2);
            REQUIRE(IntersectTest::ConnectivityOf(*r, 1e-9 * run.Scale).Curves == 2);
            std::array<std::size_t, 2> sides{0, 0};
            for (const IntersectionBranch& branch : r->Branches) {
                REQUIRE(branch.IsCertified());
                REQUIRE(branch.Range.IsBounded());
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Point3 p = Unplaced(turn, run.Scale, run.Shift, Eval(branch.Curve3d, t));
                    const double x = p.X - 0.5;
                    const double y = p.Y - 0.5;
                    CHECK(std::abs(x * x - 4.0 * y * y - 0.01) <= 8.0 * Context::DefaultResolution);
                    CHECK(std::abs(p.Z - (2.0 * x * x - 3.0 * y * y)) <=
                          4.0 * Context::DefaultResolution);
                }
                // Each branch runs across the square from one side back to
                // the same side: both ends on x = 0, or both on x = 1.
                for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
                    const Point3 p = Unplaced(turn, run.Scale, run.Shift, Eval(branch.Curve3d, t));
                    const double reach = 4.0 * Context::DefaultResolution;
                    CHECK(std::abs(std::abs(p.Y - 0.5) - half) <= reach);
                    if (std::abs(p.X) <= reach) {
                        ++sides[0];
                    } else if (std::abs(p.X - 1.0) <= reach) {
                        ++sides[1];
                    }
                }
            }
            CHECK(sides[0] == 2);
            CHECK(sides[1] == 2);
        }
    }
}

TEST_CASE("MR SSI-translated-boundary resolves both roots on one patch edge",
          "[intersect][generic][metamorphic][backfill][SSI-translated-boundary]")
{
    // In local coordinates the second patch is (s*u, s*v, s*(u + (v-1/2)^2 - 1/16)).
    // Its section by z=0 is one parabola. Both ends are on u=0 at v=1/4 and 3/4.
    // The wider plane contains the whole arc. With halving disabled the edge
    // search cannot isolate both roots before the pair reaches the tolerance.
    constexpr double s = 0x1p-12;
    constexpr double pairTolerance = 5e-5;
    constexpr double resolution = 1e-8;
    constexpr std::array<double, 3> square{0.25, -0.25, 0.25};
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const Context policy = Context::Default().WithResolution(resolution);
    for (const double shift : {0.0, 1.0}) {
        std::vector<Point3> planePoles;
        std::vector<Point3> curvedPoles;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double u = static_cast<double>(i) / 2.0;
                const double v = static_cast<double>(j) / 2.0;
                planePoles.push_back({shift + s * (2.0 * u - 0.5), s * (2.0 * v - 0.5), 0.0});
                curvedPoles.push_back({shift + s * u, s * v, s * (u + square.at(j) - 0.0625)});
            }
        }
        const auto plane = Stereon::MakeBSplineSurface(2, 2, knots, knots, planePoles, {});
        const auto curved = Stereon::MakeBSplineSurface(2, 2, knots, knots, curvedPoles, {});
        REQUIRE(plane.has_value());
        REQUIRE(curved.has_value());
        const Surface a{*plane};
        const Surface b{*curved};
        for (const bool swapped : {false, true}) {
            CAPTURE(shift, swapped);
            const auto result = Intersect(swapped ? b : a, Box2::Everything(), swapped ? a : b,
                                          Box2::Everything(), pairTolerance, policy);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->Branches.size() == 1);
            CHECK(result->Points.empty());
            CHECK(result->Overlaps.empty());
            const auto& branch = result->Branches.front();
            REQUIRE(branch.IsCertified());
            REQUIRE(branch.Range.IsBounded());
            const Point3 first = Eval(branch.Curve3d, branch.Range.GetLower());
            const Point3 last = Eval(branch.Curve3d, branch.Range.GetUpper());
            CHECK(std::abs(first.X - shift) <= resolution);
            CHECK(std::abs(last.X - shift) <= resolution);
            CHECK(std::abs(std::min(first.Y, last.Y) - 0.25 * s) <= resolution);
            CHECK(std::abs(std::max(first.Y, last.Y) - 0.75 * s) <= resolution);
            for (const double t : IntersectTest::SamplesOf(branch.Range, 32)) {
                const Point3 p = Eval(branch.Curve3d, t);
                const double y = p.Y / s - 0.5;
                CHECK(std::abs(p.Z) <= resolution);
                CHECK(std::abs((p.X - shift) - s * (0.0625 - y * y)) <= resolution);
            }
        }
    }
}

TEST_CASE("MR AUDIT-20261112-cone-polish far line origin keeps cone crossings",
          "[intersect][metamorphic][backfill][AUDIT-20261112-cone-polish]")
{
    // Root cause (SplineSurface.cpp, repaired in #286): a crossing of a line
    // and a cone was polished at origin + t direction rounded twice, which
    // loses the small offset when a far origin cancels t direction; near a
    // generator the root's correction magnifies it. CI seed
    // 17298607298449989786, shrunk: the line through (-1/4, -1/4, 21/8)
    // along (13, -23, 19)/8 and the right cone of semi-angle pi/4 about the
    // y axis with its apex at (23/8, -17/8, 3/2). For the ideal cone the
    // roots are t^2 + 382 t + 481 = 0, -191 +/- 60 sqrt(10), independently
    // of the kernel. MR-04: the line's origin moved back 2^19 directions,
    // t' = t + 2^19; MR-01: the three cyclic turns of the axes, exact, with
    // and without the translation (8, -4, 2).
    const Vec3 direction{13.0 / 8.0, -23.0 / 8.0, 19.0 / 8.0};
    const double run = 0x1p19;
    for (const auto& turn : Turns) {
        const Vec3 axis = Turned(turn, Vec3{0, 1, 0});
        const Vec3 across0 = Turned(turn, Vec3{0, 0, 1});
        const Vec3 across1 = Turned(turn, Vec3{1, 0, 0});
        const Vec3 d = Turned(turn, direction);
        for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
            INFO("cone-polish turn (" << turn[0].X << ", " << turn[0].Y << ", " << turn[0].Z
                                      << ") shift (" << shift.X << ", " << shift.Y << ", "
                                      << shift.Z << ")");
            const Line near{Turned(turn, Point3{-0.25, -0.25, 21.0 / 8.0}) + shift, d};
            const Line far{near.Origin - run * d, d};
            const Cone cone{Turned(turn, Point3{23.0 / 8.0, -17.0 / 8.0, 1.5}) + shift, 0.0,
                            pi / 4.0, across0, across1};
            const auto a = Intersect(Curve{near}, Surface{cone}, Context::Default());
            const auto b = Intersect(Curve{far}, Surface{cone}, Context::Default());
            INFO(IntersectTest::ErrorText(a));
            INFO(IntersectTest::ErrorText(b));
            REQUIRE(a.has_value());
            REQUIRE(b.has_value());
            REQUIRE(a->Points.size() == 2);
            REQUIRE(b->Points.size() == 2);
            const auto rounding = [&](const Line& line, const CurveSurfacePoint& point) {
                return 512.0 * std::numeric_limits<double>::epsilon() *
                       std::max({1.0, Length(point.Position - Point3{}),
                                 std::abs(point.Parameter) * Length(line.Direction)});
            };
            for (std::size_t k = 0; k < 2; ++k) {
                const CurveSurfacePoint& p = a->Points[k];
                const CurveSurfacePoint& q = b->Points[k];
                const double exact = -191.0 + (k == 0 ? -60.0 : 60.0) * std::sqrt(10.0);
                const double within = rounding(near, p) + rounding(far, q);
                // The cone's own support, in its frame: the distance from the
                // axis times cos(a) against the height times sin(a).
                for (const CurveSurfacePoint* point : {&p, &q}) {
                    const Vec3 w = point->Position - cone.Origin;
                    CHECK(std::abs(std::hypot(Dot(w, across0), Dot(w, across1)) *
                                       std::cos(cone.SemiAngle) -
                                   std::abs(Dot(w, axis)) * std::sin(cone.SemiAngle)) <= within);
                }
                CHECK(Distance(p.Position, Eval(near, p.Parameter)) <= rounding(near, p));
                CHECK(Distance(q.Position, Eval(far, q.Parameter)) <= rounding(far, q));
                CHECK(Distance(q.Position, Eval(near, exact)) <= within);
                CHECK(Distance(p.Position, q.Position) <= within);
                CHECK(std::abs(q.Parameter - p.Parameter - run) <= within / Length(direction));
                CHECK(p.Contact == ContactType::Crossing);
                CHECK(q.Contact == ContactType::Crossing);
            }
        }
    }
}

namespace {

/// The plane section of a sphere, in closed form: one closed curve, every
/// sample on the plane and at the circle's radius from its centre, and its
/// samples turning once about the centre.
void CheckSectionCircle(const SurfaceSurfaceIntersection& result, const Plane& plane,
                        const Sphere& sphere, const double resolution)
{
    CHECK(result.Points.empty());
    CHECK(result.Overlaps.empty());
    REQUIRE_FALSE(result.Branches.empty());
    CHECK(IntersectTest::ConnectivityOf(result, 1e-7 * sphere.Radius).Curves == 1);
    const Vec3 normal = Cross(plane.XAxis, plane.YAxis);
    const Vec3 unit = normal / Length(normal);
    const double height = Dot(sphere.Center - plane.Origin, unit);
    const Point3 center = sphere.Center - height * unit;
    const double radius = std::sqrt(sphere.Radius * sphere.Radius - height * height);
    double sweep = 0.0;
    for (const IntersectionBranch& branch : result.Branches) {
        CHECK(branch.IsCertified());
        CHECK(branch.Bound <= resolution);
        CHECK(branch.Contact == ContactType::Crossing);
        double previous = 0.0;
        double turned = 0.0;
        std::size_t index = 0;
        for (const double t : IntersectTest::SamplesOf(branch.Range, 128)) {
            const Point3 point = Eval(branch.Curve3d, t);
            CHECK(std::abs(Dot(point - plane.Origin, unit)) <= resolution);
            CHECK(std::abs(Distance(point, center) - radius) <= resolution);
            const Vec3 delta = point - center;
            const double angle = std::atan2(Dot(delta, plane.YAxis), Dot(delta, plane.XAxis));
            if (index++ > 0) {
                turned += std::remainder(angle - previous, 2.0 * pi);
            }
            previous = angle;
        }
        sweep += std::abs(turned);
    }
    CHECK(std::abs(sweep - 2.0 * pi) <=
          8.0 * static_cast<double>(result.Branches.size()) * resolution / radius);
}

} // namespace

TEST_CASE("MR AUDIT-20261113-narrow-sphere-seam plane section stays one whole circle",
          "[intersect][metamorphic][backfill][AUDIT-20261113-narrow-sphere-seam]")
{
    // Root cause (ImplicitTrace.cpp, repaired in #286): the circle a plane
    // cuts from a sphere crossed the sphere's seam twice between adjacent
    // samples of the trace's first scan; the fit clamped that excursion to
    // the seam and never certified, failing a supported section with the
    // sample budget. Coverage seed 13593231413355092879, shrunk: centre
    // (-5, -5, -5), radius 1, plane normal (-2.5, -4.875, -3) at the signed
    // offset -5/8. The oracle is the section's closed form: the circle in
    // the plane about the centre's foot, of radius sqrt(1 - (5/8)^2), whole.
    // MR-01: a cyclic turn of the axes and the translation (8, -4, 2),
    // exact; MR-02: both orders; MR-03: the scale 2 with the context;
    // MR-04: the sphere's chart turned a quarter about its axis.
    const Vec3 normal{-2.5, -4.875, -3.0};
    const Vec3 ex{1, 0, 0};
    const Vec3 ey{0, 1, 0};
    for (const std::size_t t : {std::size_t{0}, std::size_t{1}}) {
        const auto& turn = Turns.at(t);
        for (const double scale : {1.0, 2.0}) {
            for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
                if (scale != 1.0 && (t != 0 || shift.X != 0.0)) {
                    continue;
                }
                const Context policy =
                    Context::Default()
                        .WithResolution(Context::DefaultResolution * scale)
                        .WithToleranceCeiling(Context::DefaultToleranceCeiling * scale);
                const Vec3 n = Turned(turn, normal);
                const Point3 center =
                    Point3{} + scale * Turned(turn, Vec3{-5.0, -5.0, -5.0}) + shift;
                const Plane plane = IntersectTest::PlaneWithNormal(
                    center + (-5.0 / 8.0 * scale) * (n / Length(n)), n);
                const Vec3 x = Turned(turn, ex);
                const Vec3 y = Turned(turn, ey);
                for (const bool quarter : {false, true}) {
                    const Sphere sphere{center, scale, quarter ? y : x, quarter ? -x : y};
                    for (const bool exchanged : {false, true}) {
                        INFO("narrow-sphere-seam turn " << t << " scale " << scale << " shift "
                                                        << shift.X << " quarter " << quarter
                                                        << " exchanged " << exchanged);
                        const Surface a = exchanged ? Surface{sphere} : Surface{plane};
                        const Surface b = exchanged ? Surface{plane} : Surface{sphere};
                        const auto r = Intersect(a, b, policy);
                        INFO(IntersectTest::ErrorText(r));
                        REQUIRE(r.has_value());
                        for (const IntersectionBranch& branch : r->Branches) {
                            CHECK(IntersectTest::PcurveDeviation(branch, a, b, 200) <=
                                  policy.GetResolution());
                        }
                        CheckSectionCircle(*r, plane, sphere, policy.GetResolution());
                    }
                }
            }
        }
    }
}

TEST_CASE("MR CI-272-ellipse-plane crossings on the seam and a cut are decided",
          "[intersect][metamorphic][backfill][CI-272-ellipse-plane]")
{
    // Root cause (SplineSurface.cpp, repaired in #300): a conic is searched
    // in quarter-turn arcs, and a crossing exactly on a cut or on the seam
    // was seen undecided from both arcs and refused ("unresolved root
    // existence or multiplicity"). macOS run 43273260370, RapidCheck seed
    // 17041817369036985825, shrunk: the ellipse of semi-axes 3/2 and 3/5
    // about (-5/4, -5/4, -13/32), axes x and (0, cos phi, sin phi), phi the
    // rounded -pi/2, against the plane through its centre with axes x and
    // y. The plane's equation on the ellipse is b sin(phi) sin(t) = 0: two
    // transverse crossings, centre -+ a x, at t = pi (a cut) and t = 0 (the
    // seam), independently of the kernel. MR-01: the cyclic turns of the
    // axes and the translation (8, -4, 2), exact; MR-03: scales 1/4 and 4
    // with the context; MR-04: the ellipse reversed, and its frame turned
    // half (t' = t - pi).
    const double phi = -16.0 * (2.0 * pi) / 64.0;
    for (const auto& turn : Turns) {
        for (const double scale : {0.25, 1.0, 4.0}) {
            for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
                const Context policy =
                    Context::Default()
                        .WithResolution(Context::DefaultResolution * scale)
                        .WithToleranceCeiling(Context::DefaultToleranceCeiling * scale)
                        .WithScheduler(Scheduler::Serial);
                const double major = 1.5 * scale;
                const double minor = scale * (1.5 * 4.0 / 10.0);
                const Point3 center =
                    Point3{} +
                    scale * Turned(turn, Vec3{-40.0 / 32.0, -40.0 / 32.0, -13.0 / 32.0}) + shift;
                const Vec3 x = Turned(turn, Vec3{1, 0, 0});
                const Vec3 y = Turned(turn, Vec3{0, std::cos(phi), std::sin(phi)});
                const Surface plane{Plane{center, x, Turned(turn, Vec3{0, 1, 0})}};
                const Ellipse ellipse{center, major, minor, x, y};
                const std::array<Curve, 3> charts{Curve{ellipse}, Curve{Reversed(ellipse)},
                                                  Curve{Ellipse{center, major, minor, -x, -y}}};
                for (std::size_t c = 0; c < charts.size(); ++c) {
                    const Curve& curve = charts.at(c);
                    INFO("CI-272 turn (" << turn[0].X << ", " << turn[0].Y << ", " << turn[0].Z
                                         << ") scale " << scale << " shift " << shift.X << " chart "
                                         << c);
                    const auto r = Intersect(curve, plane, policy);
                    INFO(IntersectTest::ErrorText(r));
                    REQUIRE(r.has_value());
                    CHECK(r->Overlaps.empty());
                    // A closed turn may keep both representatives of the
                    // seam; there are only two spatial roots.
                    CHECK(r->Points.size() >= 2);
                    CHECK(r->Points.size() <= 3);
                    bool low = false;
                    bool high = false;
                    for (const CurveSurfacePoint& point : r->Points) {
                        const bool atLow =
                            Distance(point.Position, center - major * x) <= policy.GetResolution();
                        const bool atHigh =
                            Distance(point.Position, center + major * x) <= policy.GetResolution();
                        CHECK((atLow || atHigh));
                        low = low || atLow;
                        high = high || atHigh;
                        CHECK(point.Decided);
                        CHECK(point.Contact == ContactType::Crossing);
                        CHECK(Contains(Domain(curve), point.Parameter));
                        CHECK(Distance(Eval(curve, point.Parameter), point.Position) <=
                              policy.GetResolution());
                        CHECK(Distance(Eval(plane, point.OnSurface), point.Position) <=
                              policy.GetResolution());
                    }
                    CHECK(low);
                    CHECK(high);
                }
            }
        }
    }
}

TEST_CASE("MR AUDIT-20261122-plane-cylinder-periodic-box sections cover the requested turn",
          "[intersect][metamorphic][backfill][periodic-cylinder]")
{
    const Surface cylinder{Cylinder{Point3{}, 1.0, Vec3{-1, 0, 0}, Vec3{0, -1, 0}, true}};
    const Box2 planeBox{Interval{-4, 4}, Interval{-4, 4}};
    const Context policy = Context::Default();
    for (const bool tilted : {false, true}) {
        const Vec3 y = tilted ? Vec3{0, 0.6, 0.8} : Vec3{0, 1, 0};
        const Surface plane{Plane{Point3{0, 0, 0.5}, Vec3{1, 0, 0}, y}};
        for (const Interval arc :
             {Interval{-pi, pi}, Interval{-pi, -pi / 4}, Interval{2 * pi, 2.5 * pi}}) {
            const Box2 cylinderBox{arc, Interval{-4, 4}};
            for (const bool swapped : {false, true}) {
                DYNAMIC_SECTION("tilted " << tilted << " lower " << arc.GetLower() << " upper "
                                          << arc.GetUpper() << " swapped " << swapped)
                {
                    const auto result =
                        swapped ? Intersect(cylinder, cylinderBox, plane, planeBox, 1e-5, policy)
                                : Intersect(plane, planeBox, cylinder, cylinderBox, 1e-5, policy);
                    INFO(IntersectTest::ErrorText(result));
                    REQUIRE(result.has_value());
                    REQUIRE(result->Branches.size() == 1);
                    const IntersectionBranch& branch = result->Branches.front();
                    REQUIRE(branch.IsCertified());
                    const auto& onCylinder = swapped ? branch.OnA : branch.OnB;
                    // Check the entire requested arc including both endpoints. At the
                    // quarter turn of [-pi,pi] this includes the missing (0,1,11/6).
                    for (const double f : {0.0, 0.25, 0.5, 0.75, 1.0}) {
                        const double u = std::lerp(arc.GetLower(), arc.GetUpper(), f);
                        const double t =
                            std::lerp(branch.Range.GetLower(), branch.Range.GetUpper(), f);
                        const double worldY = -std::sin(u);
                        const Point3 expected{-std::cos(u), worldY,
                                              0.5 + (tilted ? 4.0 * worldY / 3.0 : 0.0)};
                        CHECK(Distance(Eval(branch.Curve3d, t), expected) <=
                              policy.GetResolution());
                        CHECK(std::abs(Eval(onCylinder, t).X - u) <= policy.GetResolution());
                    }
                }
            }
        }
    }
}

// The repairs merged after #263 up to #323 (ADR-0053; triage batch 4 in
// docs/planning/intersection-metamorphic-triage.md). Each fixture is the
// repair's own reproducer, deterministic, with no seed.

namespace {

/// The distance of p from the cylinder through o about the unit axis a.
double OffCylinderAxis(const Point3& p, const Point3& o, const Vec3& a, const double radius)
{
    const Vec3 w = p - o;
    return std::abs(Length(w - Dot(w, a) * a) - radius);
}

/// A placement's exact turn of a direction: the cyclic permutations of the
/// axes move no bit.
Vec3 TurnOf(const Placement& p, const Vec3& v)
{
    return v.X * p.X + v.Y * p.Y + v.Z * Cross(p.X, p.Y);
}

} // namespace

TEST_CASE("MR PR-289 skew cylinders are answered in either order",
          "[intersect][metamorphic][backfill][PR-289]")
{
    // Root cause (CylinderCylinderSkew.cpp, repaired in #289): the curve
    // was parametrised by the first operand's angle, and the radius-2
    // cylinder's parametrisation of analytic-ssi/cylinder-cylinder/
    // skew-tilted-non-unit exceeds the sample budget, so (A, B) answered and
    // (B, A) refused. Oracle, in closed form: the axes, z through
    // (0.1, -0.2) and (3, 4, 0) / 5 through (-3.1, 2.2, 1.3), are 4 apart,
    // between |3 - 2| and 3 + 2, so the section is one loop; every sample of
    // every certified branch lies on both cylinders by the distance from
    // each axis. MR-02: both orders. MR-01: the cyclic turns of the axes and
    // the translation (8, -4, 2). MR-03: scales 1/4 and 4, the radii, the
    // pair tolerance and the context scaled.
    for (const Placement& p : Placements()) {
        Replay("PR-289", p);
        const Point3 oa = p.Point(0.1, -0.2, 0.3);
        const Vec3 axisA = TurnOf(p, Vec3{0.0, 0.0, 1.0});
        const Point3 ob = p.Point(-3.1, 2.2, 1.3);
        const Vec3 axisB = TurnOf(p, Vec3{0.6, 0.8, 0.0});
        const Surface a{IntersectTest::CylinderAbout(oa, axisA, 3.0 * p.Scale)};
        const Surface b{IntersectTest::CylinderAbout(ob, axisB, 2.0 * p.Scale)};
        const Context policy = p.Policy();
        const double reach = 4.0 * policy.GetResolution();
        const Box2 all = Box2::Everything();
        for (const bool swapped : {false, true}) {
            INFO("PR-289 scale " << p.Scale << " offset " << p.Offset.X << " axis (" << p.X.X
                                 << ", " << p.X.Y << ", " << p.X.Z << ") "
                                 << (swapped ? "B, A" : "A, B"));
            const auto r = swapped ? Intersect(b, all, a, all, 1e-3 * p.Scale, policy)
                                   : Intersect(a, all, b, all, 1e-3 * p.Scale, policy);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            CHECK(r->Overlaps.empty());
            CHECK(IntersectTest::ConnectivityOf(*r, 1e-6 * p.Scale).Curves == 1);
            REQUIRE_FALSE(r->Branches.empty());
            for (const IntersectionBranch& branch : r->Branches) {
                CHECK(branch.IsCertified());
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Point3 q = Eval(branch.Curve3d, t);
                    CHECK(OffCylinderAxis(q, oa, axisA, 3.0 * p.Scale) <= reach);
                    CHECK(OffCylinderAxis(q, ob, axisB, 2.0 * p.Scale) <= reach);
                }
            }
        }
    }
}

namespace {

/// A circle's or an ellipse's point and tangent at t, in closed form.
Point3 ConicPoint(const Point3& c, const double a, const double b, const Vec3& x, const Vec3& y,
                  const double t)
{
    return c + a * std::cos(t) * x + b * std::sin(t) * y;
}

Vec3 ConicTangent(const double a, const double b, const Vec3& x, const Vec3& y, const double t)
{
    return -a * std::sin(t) * x + b * std::cos(t) * y;
}

/// A placed copy of a local point: turned exactly, scaled by a power of two
/// and moved by an exactly representable shift.
Point3 PlacedAt(const std::array<Vec3, 3>& turn, const double scale, const Vec3& shift,
                const Point3& q)
{
    return Point3{} + scale * Turned(turn, q - Point3{}) + shift;
}

Context ScaledPolicy(const double scale)
{
    return Context::Default()
        .WithResolution(Context::DefaultResolution * scale)
        .WithToleranceCeiling(Context::DefaultToleranceCeiling * scale)
        .WithScheduler(Scheduler::Serial);
}

/// PR-304's claim for one move of B, in either order: where A rests on its
/// upper end and B is free, B's parameter is the closed-form angle of A's end
/// point in B's frame. 1 when the case was of that kind and checked.
int FreeParameterPlaced(const Circle& a, const Interval& domainA, const Circle& b,
                        const Point3& end, const bool swapped, const Context& policy)
{
    const Interval domainB{0.0, 2.0 * pi};
    const auto r = swapped ? ClosestApproach(Curve{b}, domainB, Curve{a}, domainA, policy)
                           : ClosestApproach(Curve{a}, domainA, Curve{b}, domainB, policy);
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    const ApproachEnd onA = swapped ? r->EndB : r->EndA;
    const ApproachEnd onB = swapped ? r->EndA : r->EndB;
    if (r->Uniqueness != ApproachUniqueness::Unique || onA != ApproachEnd::Upper ||
        onB != ApproachEnd::Interior) {
        return 0;
    }
    const Vec3 v = end - b.Center;
    double expected = std::atan2(Dot(v, b.YAxis), Dot(v, b.XAxis));
    if (expected < 0.0) {
        expected += 2.0 * pi;
    }
    const double got = swapped ? r->ParameterA : r->ParameterB;
    INFO("parameter " << got << ", expected " << expected);
    CHECK(std::abs(got - expected) <= 1e-13);
    return 1;
}

} // namespace

TEST_CASE("MR PR-304 closest approach places the free parameter to rounding",
          "[intersect][metamorphic][backfill][PR-304]")
{
    // Root cause (ClosestApproach.cpp, repaired in #304): Refine took a
    // backtracking step only if phi = |A - B|^2 / 2 did not rise. Near the
    // minimiser a step of size h lowers phi by about h^2, below its rounding
    // once h is near sqrt(eps), so B's free parameter stopped about 1e-8
    // short. The fixture is the shrunk macOS seed 9376910471991398838: A
    // rests on its upper end, B is free, and B is moved by 17 x 17 small
    // steps. Oracle, in closed form: B's parameter is the angle of A's end
    // point projected onto B's plane, atan2 in B's frame. MR-01: the cyclic
    // turns of the axes and the translation (8, -4, 2), exact on these dyadic
    // centres. MR-02: both orders, the ends exchanged. MR-03: scales 1/4 and
    // 4 with the context; the angles do not scale.
    const Point3 ca{-0.125, -0.625, -1.25};
    const Vec3 xa{-0.46303632848027881, 0.8662800404646892, -0.1875};
    const Vec3 ya{0.78096480768274434, 0.49879441994300905, 0.37589638970701905};
    const Point3 cb{-1.625, -0.875, -0.875};
    const Vec3 xb{0.034108521304589856, 0.34630962847488705, -0.9375};
    const Vec3 yb{0.98021459587788196, -0.19459554818054772, -0.036220417671148579};
    const double radius = 1.125;
    const Interval domainA{0.98174770424681035, 1.9817477042468103};
    for (const auto& turn : Turns) {
        for (const double scale : {0.25, 1.0, 4.0}) {
            for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
                const Context policy = ScaledPolicy(scale);
                const Vec3 ax = Turned(turn, xa);
                const Vec3 ay = Turned(turn, ya);
                const Vec3 bx = Turned(turn, xb);
                const Vec3 by = Turned(turn, yb);
                const Circle a{PlacedAt(turn, scale, shift, ca), radius * scale, ax, ay};
                const Point3 end =
                    ConicPoint(a.Center, a.Radius, a.Radius, ax, ay, domainA.GetUpper());
                for (const bool swapped : {false, true}) {
                    int placed = 0;
                    for (int i = -8; i <= 8; ++i) {
                        for (int j = -8; j <= 8; ++j) {
                            const Vec3 step{i * 1.25e-4, j * 1.25e-4, (i - j) * 0.5e-4};
                            const Circle b{PlacedAt(turn, scale, shift, cb + step), radius * scale,
                                           bx, by};
                            INFO("PR-304 turn (" << turn[0].X << ", " << turn[0].Y << ", "
                                                 << turn[0].Z << ") scale " << scale << " shift "
                                                 << shift.X << " step (" << i << ", " << j << ") "
                                                 << (swapped ? "B, A" : "A, B"));
                            placed += FreeParameterPlaced(a, domainA, b, end, swapped, policy);
                        }
                    }
                    CHECK(placed > 200);
                }
            }
        }
    }
}

TEST_CASE("MR PR-312 closest approach steps on after one parameter lands on an end",
          "[intersect][metamorphic][backfill][PR-312]")
{
    // Root cause (ClosestApproach.cpp, repaired in #312): a coupled Newton
    // step that landed B's parameter on its end moved A's back by a
    // rounding, and the refinement took that as convergence, A 1.5e-7 short
    // of its own end and the distance 1.4e-7 too large. The fixture is the
    // Nightly counterexample (2026-11-17) of the long double oracle's
    // property: a circle arc and an ellipse arc whose minimum distance is
    // at the upper end of both, as that oracle found. Here the answer must
    // be that corner: both ends Upper, the distance the corner's in closed
    // form to 1e-12 relative, and the corner a first-order minimum on the
    // box (phi does not decrease into either domain). MR-01: the cyclic
    // turns of the axes and the translation (8, -4, 2). MR-02: both orders.
    // MR-03: scales 1/4 and 4 with the context.
    const double tau = 2.0 * pi;
    const auto unit = [&](const int theta, const int z) {
        const double h = z / 16.0;
        const double r = std::sqrt(1.0 - h * h);
        return Vec3{r * std::cos(theta * tau / 64.0), r * std::sin(theta * tau / 64.0), h};
    };
    const auto frame = [&](const int theta, const int z, const int turn) {
        const Vec3 x = unit(theta, z);
        const Vec3 helper = std::abs(x.Z) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
        const Vec3 y = Cross(helper, x) / Length(Cross(helper, x));
        const double phi = turn * tau / 32.0;
        return std::array<Vec3, 2>{x, std::cos(phi) * y + std::sin(phi) * Cross(x, y)};
    };
    const std::array<Vec3, 2> fa = frame(60, 11, 22);
    const std::array<Vec3, 2> fb = frame(50, -7, 4);
    const Interval da{26 * tau / 32.0, 26 * tau / 32.0 + 3.0};
    const Interval db{21 * tau / 32.0, 21 * tau / 32.0 + 0.5};
    for (const auto& turn : Turns) {
        for (const double scale : {0.25, 1.0, 4.0}) {
            for (const Vec3 shift : {Vec3{}, Vec3{8.0, -4.0, 2.0}}) {
                const Context policy = ScaledPolicy(scale);
                const Vec3 ax = Turned(turn, fa[0]);
                const Vec3 ay = Turned(turn, fa[1]);
                const Vec3 bx = Turned(turn, fb[0]);
                const Vec3 by = Turned(turn, fb[1]);
                const Point3 ca = PlacedAt(turn, scale, shift, Point3{-0.625, -0.375, -0.125});
                const Point3 cb = PlacedAt(turn, scale, shift, Point3{0.375, 1.5, -1.125});
                const double ra = 1.375 * scale;
                const double major = 1.875 * scale;
                const double minor = 1.875 * 0.9 * scale;
                const Curve a{Circle{ca, ra, ax, ay}};
                const Curve b{Ellipse{cb, major, minor, bx, by}};
                const Point3 endA = ConicPoint(ca, ra, ra, ax, ay, da.GetUpper());
                const Point3 endB = ConicPoint(cb, major, minor, bx, by, db.GetUpper());
                const Vec3 gap = endA - endB;
                const double corner = Length(gap);
                // phi's partials at the corner: into A's domain s decreases,
                // so phi must not fall there, nor into B's.
                CHECK(Dot(gap, ConicTangent(ra, ra, ax, ay, da.GetUpper())) <= 0.0);
                CHECK(Dot(gap, ConicTangent(major, minor, bx, by, db.GetUpper())) >= 0.0);
                for (const bool swapped : {false, true}) {
                    INFO("PR-312 turn (" << turn[0].X << ", " << turn[0].Y << ", " << turn[0].Z
                                         << ") scale " << scale << " shift " << shift.X << " "
                                         << (swapped ? "B, A" : "A, B"));
                    const auto r = swapped ? ClosestApproach(b, db, a, da, policy)
                                           : ClosestApproach(a, da, b, db, policy);
                    INFO(IntersectTest::ErrorText(r));
                    REQUIRE(r.has_value());
                    CHECK(r->EndA == ApproachEnd::Upper);
                    CHECK(r->EndB == ApproachEnd::Upper);
                    INFO("distance " << r->Distance << ", corner " << corner);
                    CHECK(std::abs(r->Distance - corner) <= 1e-12 * corner);
                }
            }
        }
    }
}

namespace {

/// A Bezier in the plane on [0, 1], of degree poles - 1.
BSplineCurve2d Bezier2(std::vector<Point2> poles)
{
    const std::size_t n = poles.size();
    std::vector<double> knots(n, 0.0);
    knots.resize(2 * n, 1.0);
    auto made = MakeBSplineCurve(static_cast<int>(n) - 1, std::move(knots), std::move(poles));
    REQUIRE(made.has_value());
    return *made;
}

/// A plane placement: a quarter turn applied `Turns` times, then a scale by
/// a power of two and a shift, all exact on dyadic data.
struct Placement2 {
    int Turns = 0;
    double Scale = 1.0;
    Vec2 Shift{};

    [[nodiscard]] Vec2 Turn(Vec2 v) const
    {
        for (int k = 0; k < Turns; ++k) {
            v = Vec2{-v.Y, v.X};
        }
        return v;
    }
    [[nodiscard]] Point2 Of(const Point2& p) const
    {
        const Vec2 v = Turn(Vec2{p.X, p.Y});
        return Point2{Shift.X + Scale * v.X, Shift.Y + Scale * v.Y};
    }
};

std::vector<Placement2> Placements2(const Vec2 far)
{
    std::vector<Placement2> out;
    for (const int turns : {0, 1}) {
        for (const double scale : {0.25, 1.0, 4.0}) {
            for (const Vec2 shift : {Vec2{}, far}) {
                out.push_back({turns, scale, shift});
            }
        }
    }
    return out;
}

/// A Bezier's point by de Casteljau in long double, independent of the
/// kernel's evaluation.
std::array<long double, 2> BezierAt(const std::vector<Point2>& poles, const double t)
{
    std::vector<std::array<long double, 2>> q;
    for (const Point2& p : poles) {
        q.push_back({p.X, p.Y});
    }
    const long double s = t;
    for (std::size_t n = q.size(); n > 1; --n) {
        for (std::size_t i = 0; i + 1 < n; ++i) {
            q[i] = {(1 - s) * q[i][0] + s * q[i + 1][0], (1 - s) * q[i][1] + s * q[i + 1][1]};
        }
    }
    return q.front();
}

/// The poles placed, and reversed when asked (t' = 1 - t).
std::vector<Point2> PlacedPoles(const std::vector<Point2>& poles, const Placement2& at,
                                const bool reversed)
{
    std::vector<Point2> placed;
    for (const Point2& p : poles) {
        placed.push_back(at.Of(p));
    }
    if (reversed) {
        std::ranges::reverse(placed);
    }
    return placed;
}

/// PR-316's oracle for one reported point: a decided crossing, both
/// operands at its parameters (de Casteljau in long double) and its
/// position within 64 units in the last place of the coordinates' size.
void CheckOnBothBeziers(const PlanarPoint& p, const std::vector<Point2>& qa,
                        const std::vector<Point2>& qb, const bool swapped)
{
    const double eps = std::numeric_limits<double>::epsilon();
    CHECK(p.Contact == ContactType::Crossing);
    CHECK(p.Decided);
    const double s = swapped ? p.ParameterB : p.ParameterA;
    const double t = swapped ? p.ParameterA : p.ParameterB;
    const std::array<long double, 2> onA = BezierAt(qa, s);
    const std::array<long double, 2> onB = BezierAt(qb, t);
    const double size = std::max({1.0, std::abs(p.Position.X), std::abs(p.Position.Y)});
    const auto gap = static_cast<double>(std::hypot(onA[0] - onB[0], onA[1] - onB[1]));
    const auto off = static_cast<double>(std::hypot(onA[0] - p.Position.X, onA[1] - p.Position.Y));
    INFO("s = " << s << ", t = " << t << ": |A(s) - B(t)| = " << gap);
    CHECK(gap <= 64.0 * eps * size);
    CHECK(off <= 64.0 * eps * size);
}

} // namespace

TEST_CASE("MR PR-316 proved crossings of two cubic Beziers lie on both curves",
          "[intersect][metamorphic][backfill][PR-316]")
{
    // Root cause (PlanarSplinePair.cpp, repaired in #316): a decided
    // crossing was placed at the middle of brackets that Krawczyk's
    // refinement on the subdivided pieces could not narrow below their
    // poles' enclosure width, 8e-13, so the curves' points at its
    // parameters were 1.5e-13 apart and IntersectJvp refused the primal's
    // own point as stale. The pair is the shrunk RapidCheck seed
    // 15971369584168069721. Oracle: each operand evaluated at the reported
    // parameter by de Casteljau in long double, independently of the kernel:
    // both points, and the reported position, within 64 units in the last
    // place of the coordinates' size. Two decided crossings. MR-01: a
    // quarter turn and the translation (8, -4), exact on these dyadic poles.
    // MR-02: both orders. MR-03: scales 1/4 and 4. MR-04: either cubic
    // reversed (poles reversed, t' = 1 - t). The positions must also be the
    // identity placement's, placed.
    const std::vector<Point2> pa{{-3.0, 1.625}, {-1.375, -0.5}, {0.8125, -1.375}, {2.59375, 1.5}};
    const std::vector<Point2> pb{
        {-3.3125, -1.0}, {-1.03125, -0.75}, {0.5, 1.25}, {2.65625, -1.625}};
    std::vector<Point2> reference;
    for (const Placement2& at : Placements2(Vec2{8.0, -4.0})) {
        for (const int reversed : {0, 1, 2}) {
            const std::vector<Point2> qa = PlacedPoles(pa, at, reversed == 1);
            const std::vector<Point2> qb = PlacedPoles(pb, at, reversed == 2);
            const Curve2d a{Bezier2(qa)};
            const Curve2d b{Bezier2(qb)};
            for (const bool swapped : {false, true}) {
                INFO("PR-316 turns " << at.Turns << " scale " << at.Scale << " shift " << at.Shift.X
                                     << " reversed " << reversed << " "
                                     << (swapped ? "B, A" : "A, B"));
                const auto r = swapped ? Intersect(b, a, IntersectTest::ctx)
                                       : Intersect(a, b, IntersectTest::ctx);
                INFO(IntersectTest::ErrorText(r));
                REQUIRE(r.has_value());
                CHECK(r->Overlaps.empty());
                REQUIRE(r->Points.size() == 2);
                const bool identity = at.Turns == 0 && at.Scale == 1.0 && at.Shift.X == 0.0 &&
                                      reversed == 0 && !swapped;
                for (const PlanarPoint& p : r->Points) {
                    CheckOnBothBeziers(p, qa, qb, swapped);
                    const double size =
                        std::max({1.0, std::abs(p.Position.X), std::abs(p.Position.Y)});
                    if (identity) {
                        reference.push_back(p.Position);
                        continue;
                    }
                    double nearest = std::numeric_limits<double>::infinity();
                    for (const Point2& q : reference) {
                        nearest = std::min(nearest, Distance(at.Of(q), p.Position));
                    }
                    CHECK(nearest <= 1e-12 * size);
                }
            }
        }
    }
    CHECK(reference.size() == 2);
}

TEST_CASE("MR PR-271 a B-spline contact of order four with a parabola is a touch",
          "[intersect][metamorphic][backfill][PR-271]")
{
    // Root cause (PlanarSpline.cpp, repaired in #271): along the spline the
    // parabola's implicit equation vanishes at an order-4 contact with its
    // first three derivatives, so no piece about the touch is apart,
    // monotone or convex in doubles, and the search halved until its budget
    // was spent and refused the pair, although the ledger promises touches
    // reported bracketed, never dropped. Found by the planar corpus
    // (osculating cases). Oracle, in closed form: the quartic Bezier
    // (3(2s - 1), 9(2s - 1)^2 + (2s - 1)^4) is y = x^2 + (x / 3)^4 with
    // exact poles, and touches y = x^2 to order four at s = 1/2, the
    // origin. The answer is one undecided tangential point there, within
    // the fourth root of the coefficients' rounding. MR-01: a quarter turn
    // and the translation (1024, -512). MR-02: both orders. MR-03: scales
    // 2^-10 and 2^10 (the reach scales with the curve). MR-04: the spline
    // reversed, which keeps the touch at s = 1/2.
    const std::vector<Point2> cup{
        {-3.0, 10.0}, {-1.5, -1.0}, {0.0, -2.0}, {1.5, -1.0}, {3.0, 10.0}};
    for (const int turns : {0, 1}) {
        for (const double k : {0x1p-10, 1.0, 0x1p10}) {
            for (const Vec2 shift : {Vec2{}, Vec2{1024.0, -512.0}}) {
                const Placement2 at{turns, k, shift};
                const Point2 touch = at.Of(Point2{0.0, 0.0});
                // y = x^2 / k about the touch: points touch + (-t, t^2 / k),
                // focal length k / 4, turned with the spline.
                const Curve2d parabola{
                    Parabola2d{touch, k / 4.0, at.Turn(Vec2{0.0, 1.0}), at.Turn(Vec2{-1.0, 0.0})}};
                for (const bool reversed : {false, true}) {
                    std::vector<Point2> poles;
                    for (const Point2& p : cup) {
                        poles.push_back(at.Of(p));
                    }
                    if (reversed) {
                        std::ranges::reverse(poles);
                    }
                    const Curve2d spline{Bezier2(poles)};
                    const double size = std::max(std::abs(touch.X), std::abs(touch.Y));
                    const double rounding = 0x1p-52 * (1.0 + size / k);
                    const double reach = 8.0 * std::pow(rounding, 0.25);
                    for (const bool swapped : {false, true}) {
                        INFO("PR-271 turns " << turns << " scale " << k << " shift " << shift.X
                                             << " reversed " << reversed << " "
                                             << (swapped ? "parabola first" : "spline first"));
                        const auto r = swapped ? Intersect(parabola, spline, IntersectTest::ctx)
                                               : Intersect(spline, parabola, IntersectTest::ctx);
                        INFO(IntersectTest::ErrorText(r));
                        REQUIRE(r.has_value());
                        CHECK(r->Overlaps.empty());
                        REQUIRE(r->Points.size() == 1);
                        const PlanarPoint& p = r->Points.front();
                        CHECK(p.Contact == ContactType::Tangential);
                        CHECK_FALSE(p.Decided);
                        const double s = swapped ? p.ParameterB : p.ParameterA;
                        CHECK(std::abs(s - 0.5) <= reach);
                        CHECK(Distance(p.Position, touch) <= 8.0 * k * reach);
                    }
                }
            }
        }
    }
}

TEST_CASE("MR PR-318 a request ending on a surface knot excludes crossings beyond it",
          "[intersect][metamorphic][backfill][PR-318]")
{
    // Root cause (BezierPatch.hpp, repaired in #318): a request ending
    // exactly on an interior knot kept the neighbouring Bezier patch, which
    // meets it only in that knot edge, whole, because the cut there could not
    // be split off; so a crossing beyond the request was answered, outside
    // the requested domain. The fixture is #318's two-span curve-NURBS
    // regression: the degree-1 spline plane S(u, v) = (u, v, 0) with a u
    // knot at 1/2, and the line through (x, 1/2, -1) along z, x = 1/4 or 3/4.
    // Oracle, in closed form: the crossing is S(x, 1/2) at the line's t = 1.
    // A request containing u = x has exactly that point; the other half,
    // ending on the knot, has none. MR-01: the cyclic turns of the axes and
    // the translation (8, -4, 2). MR-03: scales 1/4 and 4. MR-04: the
    // surface's u reversed (u' = 1 - u, the halves exchanged) and the line
    // run from (x, 1/2, 1) down (t' = 2 - t, still 1 at the crossing).
    // MR-02 does not apply: there is no exposed surface-curve operation.
    const std::vector<double> knotsU{0.0, 0.0, 0.5, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    for (const Placement& p : Placements()) {
        Replay("PR-318", p);
        const Context policy = p.Policy();
        const double reach = 4.0 * policy.GetResolution();
        for (const bool flipU : {false, true}) {
            std::vector<Point3> poles;
            for (const double u : {0.0, 0.5, 1.0}) {
                for (const double v : {0.0, 1.0}) {
                    poles.push_back(p.Point(flipU ? 1.0 - u : u, v, 0.0));
                }
            }
            const auto built = MakeBSplineSurface(1, 1, knotsU, knotsV, poles, {});
            REQUIRE(built.has_value());
            const Surface spline{*built};
            for (const double x : {0.25, 0.75}) {
                const double ux = flipU ? 1.0 - x : x;
                const Point3 crossing = p.Point(x, 0.5, 0.0);
                for (const bool down : {false, true}) {
                    const Curve line{
                        down ? Line{p.Point(x, 0.5, 1.0), TurnOf(p, Vec3{0.0, 0.0, -p.Scale})}
                             : Line{p.Point(x, 0.5, -1.0), TurnOf(p, Vec3{0.0, 0.0, p.Scale})}};
                    for (const bool lower : {false, true}) {
                        const Box2 half{lower ? Interval{0.0, 0.5} : Interval{0.5, 1.0},
                                        Interval{0.0, 1.0}};
                        const bool holds = lower == (ux < 0.5);
                        INFO("PR-318 scale " << p.Scale << " offset " << p.Offset.X << " axis ("
                                             << p.X.X << ", " << p.X.Y << ", " << p.X.Z
                                             << ") flipU " << flipU << " x " << x << " down "
                                             << down << " lower half " << lower);
                        const auto r = Intersect(line, Interval{0.0, 2.0}, spline, half, policy);
                        INFO(IntersectTest::ErrorText(r));
                        REQUIRE(r.has_value());
                        CHECK(r->Overlaps.empty());
                        if (!holds) {
                            CHECK(r->Points.empty());
                            continue;
                        }
                        REQUIRE(r->Points.size() == 1);
                        const CurveSurfacePoint& q = r->Points.front();
                        CHECK(Distance(q.Position, crossing) <= reach);
                        CHECK(std::abs(q.Parameter - 1.0) <= 1e-9);
                        CHECK(std::abs(q.OnSurface.X - ux) <= 1e-9);
                        CHECK(std::abs(q.OnSurface.Y - 0.5) <= 1e-9);
                    }
                }
            }
        }
    }
}

namespace {

/// A rigid placement of the unit sphere fixtures: local (x, y, z) to
/// x E1 + y E2 + z E3 + O.
struct Rigid {
    Vec3 E1;
    Vec3 E2;
    Vec3 E3;
    Point3 O;

    [[nodiscard]] Point3 Of(const double x, const double y, const double z) const
    {
        return O + x * E1 + y * E2 + z * E3;
    }
    [[nodiscard]] Vec3 Local(const Point3& p) const
    {
        const Vec3 d = p - O;
        return {Dot(d, E1), Dot(d, E2), Dot(d, E3)};
    }
};

/// The unit sphere, exact: a rational quadratic full circle in u, a half
/// circle from pole to pole in v, the end rows of v collapsed to the poles.
BSplineSurface RationalSphere(const Rigid& at)
{
    const double s = std::sqrt(0.5);
    const std::vector<std::array<double, 2>> circle{{1, 0},   {1, 1},  {0, 1},  {-1, 1}, {-1, 0},
                                                    {-1, -1}, {0, -1}, {1, -1}, {1, 0}};
    const std::vector<double> circleWeights{1.0, s, 1.0, s, 1.0, s, 1.0, s, 1.0};
    const std::vector<std::array<double, 2>> half{{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}};
    const std::vector<double> halfWeights{1.0, s, 1.0, s, 1.0};
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < circle.size(); ++i) {
        for (std::size_t j = 0; j < half.size(); ++j) {
            const double rho = half[j][0];
            poles.push_back(at.Of(rho * circle[i][0], rho * circle[i][1], half[j][1]));
            weights.push_back(circleWeights[i] * halfWeights[j]);
        }
    }
    const std::vector<double> knotsU{0.0, 0.0,  0.0,  0.25, 0.25, 0.5,
                                     0.5, 0.75, 0.75, 1.0,  1.0,  1.0};
    const std::vector<double> knotsV{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    const auto made = MakeBSplineSurface(2, 2, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// The bilinear square of side 4 in the plane z = height about the axis.
BSplineSurface Level(const double height, const Rigid& at)
{
    const std::vector<Point3> poles{at.Of(-2.0, -2.0, height), at.Of(-2.0, 2.0, height),
                                    at.Of(2.0, -2.0, height), at.Of(2.0, 2.0, height)};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// The section of the unit sphere by z = height, in closed form: the circle
/// of radius sqrt(1 - height^2). Every sample of every certified branch on it
/// within the tolerance, and the branches together one curve going once
/// round it.
void CheckLevelSection(const SurfaceSurfaceIntersection& r, const double height, const Rigid& at,
                       const double tolerance)
{
    const double radius = std::sqrt(1.0 - height * height);
    CHECK(r.Points.empty());
    CHECK(r.Overlaps.empty());
    CHECK(IntersectTest::ConnectivityOf(r, 2.0 * tolerance).Curves == 1);
    double turned = 0.0;
    for (const IntersectionBranch& branch : r.Branches) {
        CHECK(branch.Contact == ContactType::Crossing);
        CHECK(branch.IsCertified());
        CHECK(branch.Bound <= IntersectTest::ctx.GetResolution());
        double previous = 0.0;
        bool first = true;
        for (const double t : IntersectTest::SamplesOf(branch.Range, 200)) {
            const Vec3 q = at.Local(Eval(branch.Curve3d, t));
            CHECK(std::abs(std::hypot(q.X, q.Y) - radius) <= tolerance);
            CHECK(std::abs(q.Z - height) <= tolerance);
            const double angle = std::atan2(q.Y, q.X);
            if (!first) {
                turned += std::abs(std::remainder(angle - previous, 2.0 * pi));
            }
            previous = angle;
            first = false;
        }
    }
    CHECK(std::abs(turned - 2.0 * pi) <= 1e-6 + 2.0 * IntersectTest::ctx.GetResolution() / radius);
}

void RunLevelSection(const char* defect, const double height, const std::vector<Rigid>& placements)
{
    constexpr double tolerance = 1e-5;
    for (const Rigid& at : placements) {
        const Surface sphere{RationalSphere(at)};
        const Surface flat{Level(height, at)};
        for (const bool swapped : {false, true}) {
            INFO(defect << " origin (" << at.O.X << ", " << at.O.Y << ", " << at.O.Z << ") E1 ("
                        << at.E1.X << ", " << at.E1.Y << ", " << at.E1.Z << ") "
                        << (swapped ? "plane first" : "sphere first"));
            const Box2 all = Box2::Everything();
            const auto r = swapped
                               ? Intersect(flat, all, sphere, all, tolerance, IntersectTest::ctx)
                               : Intersect(sphere, all, flat, all, tolerance, IntersectTest::ctx);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            CheckLevelSection(*r, height, at, tolerance);
        }
    }
}

const Rigid identityAt{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
const Rigid turnedAt{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 0, 0}};
/// The corpus's far placement and tilt (corpus/nurbs-ssi): the tilt is a
/// rational rotation, orthonormal to rounding.
const Rigid farAt{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {480000.0, -300000.0, 210000.0}};
const Rigid farTurnedAt{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {480000.0, -300000.0, 210000.0}};
const Rigid tiltedAt{{0.6, 0.8, 0.0}, {-0.48, 0.36, 0.8}, {0.64, -0.48, 0.6}, {12.0, -7.0, 5.0}};

} // namespace

TEST_CASE("MR PR-315 far rational SSI branches certify against the input surfaces",
          "[intersect][metamorphic][backfill][PR-315]")
{
    // Root cause (GenericSSI.cpp, repaired in #315): far from the origin
    // the branches were marched about it but certified against the input
    // surfaces, whose rational enclosure widens every homogeneous component
    // by a margin scaled to the largest pole, a floor of about |P| 1.5e-8
    // (0.007 mm at |P| = 6e5) at any span width. Eight far nurbs-ssi
    // reference fixtures were refused by CertifyDeviation. The fixture: the
    // exact rational unit sphere and the plane z = 3/5 at the corpus's far
    // origin (480000, -300000, 210000). Oracle, in closed form: the circle of
    // radius 4/5, once round. MR-01: the same pair at the origin, and at both
    // places with the axes turned cyclically. MR-02: both orders.
    RunLevelSection("PR-315", 0.6, {identityAt, turnedAt, farAt, farTurnedAt});
}

TEST_CASE("MR PR-317 a plane near the pole of a rational sphere cuts its small circle",
          "[intersect][metamorphic][backfill][PR-317]")
{
    // Root causes (GenericSSI.cpp, repaired in #317): (1) a rational patch's
    // partials were enclosed by the quotient rule about the origin, widened
    // by the patch's distance from it times the width of the weight's
    // derivative, so near the pole no advance direction or boundary root was
    // proved and the search refused the pair as a tangential contact; (2)
    // tilted off the world axes, the boundary roots were not isolated in
    // world axes and the search ran out of its budget. The fixture is
    // corpus/nurbs-ssi/sphere-plane/near-pole: the exact rational unit
    // sphere and the plane z = 0.999. Oracle, in closed form: the circle of
    // radius sqrt(1 - 0.999^2), once round. MR-01: the identity, a cyclic
    // turn of the axes, the corpus's tilt and its far translation. MR-02:
    // both orders.
    RunLevelSection("PR-317", 1.0 - 1e-3, {identityAt, turnedAt, tiltedAt, farAt});
}
