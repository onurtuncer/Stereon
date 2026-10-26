// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

import Stereon.Intersect;

namespace {
using namespace Stereon;
const Context ctx = Context::Default();
const Interval whole = Interval::Everything();
const Line x{Point3{}, Vec3{1.0, 0.0, 0.0}};

Result<CurveCurveIntersection> Meet(const Line& a, const Line& b, const Interval& da = whole,
                                    const Interval& db = whole)
{
    return Intersect(Curve{a}, da, Curve{b}, db, ctx);
}

void CheckPoint(const CurveCurveIntersection& result, const Line& a, const Line& b,
                const Point3& expected)
{
    REQUIRE(result.Points.size() == 1);
    REQUIRE(result.Overlaps.empty());
    const auto& p = result.Points.front();
    REQUIRE(Distance(p.Position, expected) <= ctx.GetResolution());
    REQUIRE(Distance(Eval(a, p.ParameterA), p.Position) <= ctx.GetResolution());
    REQUIRE(Distance(Eval(b, p.ParameterB), p.Position) <= ctx.GetResolution());
}
} // namespace

TEST_CASE("3D lines retain non-unit parameters and operand order", "[intersect][line3d]")
{
    const Line a{Point3{1.0, 2.0, 3.0}, Vec3{2.0, 0.0, 0.0}};
    const Line b{Point3{5.0, -4.0, 3.0}, Vec3{0.0, 3.0, 0.0}};
    const auto result = Meet(a, b);
    REQUIRE(result.has_value());
    CheckPoint(*result, a, b, Point3{5.0, 2.0, 3.0});
    REQUIRE(result->Points[0].ParameterA == 2.0);
    REQUIRE(result->Points[0].ParameterB == 2.0);
    const auto swapped = Meet(b, a);
    REQUIRE(swapped.has_value());
    CheckPoint(*swapped, b, a, Point3{5.0, 2.0, 3.0});
}

TEST_CASE("3D skew and parallel lines are separate even below resolution", "[intersect][line3d]")
{
    const double gap = 0x1p-40;
    for (const Vec3 direction : {Vec3{0.0, 1.0, 0.0}, Vec3{2.0, 0.0, 0.0}}) {
        const auto result = Meet(x, Line{Point3{0.0, 0.0, gap}, direction});
        REQUIRE(result.has_value());
        REQUIRE(result->IsEmpty());
    }
}

TEST_CASE("3D line crossings clip against exact rational parameters", "[intersect][line3d]")
{
    const Line a{Point3{}, Vec3{3.0, 0.0, 0.0}};
    const Line b{Point3{1.0, -1.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // The double 1/3 rounds DOWN: accepting the rounded parameter would
    // invent a crossing just outside this requested segment.
    const double below = 1.0 / 3.0;
    const double above = std::nextafter(below, 1.0);
    const auto excluded = Meet(a, b, Interval{0.0, below});
    REQUIRE(excluded.has_value());
    REQUIRE(excluded->IsEmpty());
    const auto included = Meet(a, b, Interval{0.0, above});
    REQUIRE(included.has_value());
    CheckPoint(*included, a, b, Point3{1.0, 0.0, 0.0});
    const auto endpoint = Meet(x, b, Interval{1.0}, Interval{1.0});
    REQUIRE(endpoint.has_value());
    CheckPoint(*endpoint, x, b, Point3{1.0, 0.0, 0.0});
}

TEST_CASE("3D collinear segments distinguish overlap endpoint and gap", "[intersect][line3d]")
{
    const Line reversed{Point3{4.0, 0.0, 0.0}, Vec3{-2.0, 0.0, 0.0}};
    const auto overlap = Meet(x, reversed, Interval{0.0, 3.0}, Interval{0.0, 1.0});
    REQUIRE(overlap.has_value());
    REQUIRE(overlap->Points.empty());
    REQUIRE(overlap->Overlaps.size() == 1);
    const auto& o = overlap->Overlaps[0];
    REQUIRE_FALSE(o.SameDirection);
    REQUIRE(Contains(o.OnA, Interval{2.0, 3.0}));
    REQUIRE(Contains(o.OnB, Interval{0.5, 1.0}));
    REQUIRE(Width(o.OnA) < 1.0 + 1e-12);
    const auto endpoint = Meet(x, reversed, Interval{0.0, 2.0}, Interval{0.0, 1.0});
    REQUIRE(endpoint.has_value());
    CheckPoint(*endpoint, x, reversed, Point3{2.0, 0.0, 0.0});
    const auto gap = Meet(x, reversed, Interval{0.0, std::nextafter(2.0, 0.0)}, Interval{0.0, 1.0});
    REQUIRE(gap.has_value());
    REQUIRE(gap->IsEmpty());
}

TEST_CASE("3D coincident whole lines rays and singleton domains", "[intersect][line3d]")
{
    const double inf = std::numeric_limits<double>::infinity();
    const Line b{Point3{2.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}};
    const auto full = Meet(x, b);
    REQUIRE(full.has_value());
    REQUIRE(full->Overlaps.size() == 1);
    REQUIRE(full->Overlaps[0].OnA == whole);
    REQUIRE(full->Overlaps[0].OnB == whole);
    REQUIRE(full->Overlaps[0].SameDirection);
    const auto ray = Meet(x, b, Interval{0.0, inf}, Interval{0.0, inf});
    REQUIRE(ray.has_value());
    REQUIRE(ray->Overlaps.size() == 1);
    REQUIRE(Contains(ray->Overlaps[0].OnA, Interval{2.0, inf}));
    const auto point = Meet(x, b, Interval{4.0}, whole);
    REQUIRE(point.has_value());
    CheckPoint(*point, x, b, Point3{4.0, 0.0, 0.0});
}

TEST_CASE("3D near-parallel lines preserve an exact crossing", "[intersect][line3d]")
{
    const double tilt = 0x1p-40;
    const Line b{Point3{0.0, tilt, 0.0}, Vec3{1.0, -tilt, 0.0}};
    const auto result = Meet(x, b);
    REQUIRE(result.has_value());
    CheckPoint(*result, x, b, Point3{1.0, 0.0, 0.0});
}

TEST_CASE("3D line validation and numerical failures are explicit", "[intersect][line3d]")
{
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const Line bad : {Line{Point3{}, Vec3{}}, Line{Point3{inf, 0.0, 0.0}, x.Direction},
                           Line{Point3{}, Vec3{nan, 0.0, 0.0}}}) {
        const auto result = Meet(x, bad);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::InvalidInput);
    }
    const auto unsafe = Meet(x, Line{Point3{}, Vec3{0x1p-1000, 0.0, 0.0}});
    REQUIRE_FALSE(unsafe.has_value());
    REQUIRE(unsafe.error().GetCode() == ErrorCode::NotSupported);
    const auto infinityOnly = Meet(x, x, Interval{inf});
    REQUIRE_FALSE(infinityOnly.has_value());
    REQUIRE(infinityOnly.error().GetCode() == ErrorCode::InvalidInput);
    const Line a{Point3{0x1p50, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    const Line b{Point3{0x1p50, 1.0, 0.0}, Vec3{0.0, -1.0, 0.0}};
    const auto uncertified = Meet(a, b);
    REQUIRE_FALSE(uncertified.has_value());
    REQUIRE(uncertified.error().GetCode() == ErrorCode::CertificationFailed);
}

TEST_CASE("3D line intersections obey generated exact constructions",
          "[intersect][line3d][property]")
{
    rc::prop("known integer contacts survive translations and direction changes", [] {
        const auto integer = [] { return static_cast<double>(*rc::gen::inRange(-32, 33)); };
        const Point3 p{integer(), integer(), integer()};
        const Vec3 da{1.0, integer(), integer()};
        const Vec3 db{0.0, 1.0, integer()};
        const double t = integer();
        const double s = integer();
        const Line a{p - t * da, da};
        const Line b{p - s * db, db};
        const auto result = Meet(a, b);
        RC_ASSERT(result.has_value());
        RC_ASSERT(result->Points.size() == 1);
        RC_ASSERT(result->Points[0].ParameterA == t);
        RC_ASSERT(result->Points[0].ParameterB == s);
        RC_ASSERT(Distance(result->Points[0].Position, p) <= ctx.GetResolution());
        const auto swapped = Meet(b, a);
        RC_ASSERT(swapped.has_value());
        RC_ASSERT(swapped->Points.size() == 1);
        RC_ASSERT(swapped->Points[0].ParameterA == s);
        RC_ASSERT(swapped->Points[0].ParameterB == t);
    });
}

TEST_CASE("3D line determinant cancellation uses the exact fallback", "[intersect][line3d]")
{
    const double e = 0x1p-27;
    const Line a{Point3{}, Vec3{1.0, 1.0 + e, 0.0}};
    const Line b{Point3{}, Vec3{1.0 - e, 1.0, 0.0}};
    // Ordinary products round this determinant to zero; its exact value is e^2.
    const auto result = Meet(a, b);
    REQUIRE(result.has_value());
    CheckPoint(*result, a, b, Point3{});
}

TEST_CASE("3D collinear generated segments preserve overlap orientation",
          "[intersect][line3d][property]")
{
    rc::prop("reparameterized segments cover the same known interval", [] {
        const double start = static_cast<double>(*rc::gen::inRange(-32, 33));
        const double speed = static_cast<double>(*rc::gen::inRange(1, 17));
        const Line a{Point3{start, 2.0, 3.0}, Vec3{speed, 2.0 * speed, 0.0}};
        const Line b{Eval(a, 4.0), -2.0 * a.Direction};
        const auto result = Meet(a, b, Interval{0.0, 3.0}, Interval{0.0, 1.0});
        RC_ASSERT(result.has_value());
        RC_ASSERT(result->Points.empty());
        RC_ASSERT(result->Overlaps.size() == 1);
        RC_ASSERT(!result->Overlaps[0].SameDirection);
        RC_ASSERT(Contains(result->Overlaps[0].OnA, Interval{2.0, 3.0}));
        RC_ASSERT(Contains(result->Overlaps[0].OnB, Interval{0.5, 1.0}));
        RC_ASSERT(Width(result->Overlaps[0].OnA) < 1.0 + 1e-12);
    });
}
