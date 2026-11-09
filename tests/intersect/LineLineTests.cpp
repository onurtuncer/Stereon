// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
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

/// The segment of `line` over `domain` as a degree-1 B-spline on [0, 1].
BSplineCurve Flat(const Line& line, const Interval& domain)
{
    return BSplineCurve{1,
                        {0.0, 0.0, 1.0, 1.0},
                        {Eval(line, domain.GetLower()), Eval(line, domain.GetUpper())},
                        {}};
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

TEST_CASE("3D skew lines meet within the resolution, decided exactly", "[intersect][line3d]")
{
    // Lines meet where they come within the Context's resolution of each
    // other, as two B-splines do. The decision is exact: a gap of exactly
    // the resolution meets, the next double above it does not.
    const double resolution = ctx.GetResolution();
    for (const double gap : {0x1p-40, resolution}) {
        const Line b{Point3{0.0, 0.0, gap}, Vec3{0.0, 1.0, 0.0}};
        const auto result = Meet(x, b);
        REQUIRE(result.has_value());
        REQUIRE(result->Points.size() == 1);
        REQUIRE(result->Overlaps.empty());
        // Placed on the first line, at the closest approach.
        CHECK(result->Points[0].Position == Point3{});
        CHECK(result->Points[0].ParameterA == 0.0);
        CHECK(result->Points[0].ParameterB == 0.0);
    }
    const auto apart =
        Meet(x, Line{Point3{0.0, 0.0, std::nextafter(resolution, 1.0)}, Vec3{0.0, 1.0, 0.0}});
    REQUIRE(apart.has_value());
    CHECK(apart->IsEmpty());
}

TEST_CASE("3D parallel lines within the resolution overlap, beyond it are apart",
          "[intersect][line3d]")
{
    // Exactly parallel supports closer than the resolution cannot be told
    // from collinear ones, so they overlap over the stretch both cover;
    // ADR-0029: supports not exactly parallel meet at a point instead.
    const Line near{Point3{1.0, 0x1p-40, 0.0}, Vec3{2.0, 0.0, 0.0}};
    const auto overlap = Meet(x, near, Interval{0.0, 3.0}, Interval{0.0, 4.0});
    REQUIRE(overlap.has_value());
    REQUIRE(overlap->Points.empty());
    REQUIRE(overlap->Overlaps.size() == 1);
    CHECK(overlap->Overlaps[0].SameDirection);
    CHECK(Contains(overlap->Overlaps[0].OnA, Interval{1.0, 3.0}));
    CHECK(Contains(overlap->Overlaps[0].OnB, Interval{0.0, 1.0}));
    const Line far{Point3{1.0, 2e-5, 0.0}, Vec3{2.0, 0.0, 0.0}};
    const auto apart = Meet(x, far, Interval{0.0, 3.0}, Interval{0.0, 4.0});
    REQUIRE(apart.has_value());
    CHECK(apart->IsEmpty());
    // End to end across a gap within the resolution: a point at the ends.
    const auto ends = Meet(x, near, Interval{0.0, 1.0 - 0x1p-30}, Interval{0.0, 4.0});
    REQUIRE(ends.has_value());
    REQUIRE(ends->Points.size() == 1);
    CHECK(ends->Points[0].ParameterA == 1.0 - 0x1p-30);
    CHECK(ends->Points[0].ParameterB == 0.0);
}

TEST_CASE("3D segments a rounding apart meet at their ends", "[intersect][line3d]")
{
    // Two edges sharing a vertex, their ends a rounding apart, meet at it,
    // and so does a segment ending a rounding short of another's interior.
    const Line a{Point3{}, Vec3{1.0, 1.0, 0.0}};
    const Line b{Point3{1.0, 1.0 + 1e-12, 0.0}, Vec3{1.0, -1.0, 0.5}};
    const auto vertex = Meet(a, b, Interval{0.0, 1.0}, Interval{0.0, 1.0});
    REQUIRE(vertex.has_value());
    REQUIRE(vertex->Points.size() == 1);
    // A's end projects onto B just inside B's start: s = 1e-12 / |d|^2.
    CHECK(vertex->Points[0].ParameterA == 1.0);
    CHECK(std::abs(vertex->Points[0].ParameterB - 1e-12 / 2.25) <= 1e-15);
    const Line stem{Point3{0.5, 3e-6, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const auto tee = Meet(x, stem, Interval{0.0, 1.0}, Interval{0.0, 1.0});
    REQUIRE(tee.has_value());
    REQUIRE(tee->Points.size() == 1);
    CHECK(tee->Points[0].ParameterA == 0.5);
    CHECK(tee->Points[0].ParameterB == 0.0);
    const Line short_{Point3{0.5, 2e-5, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const auto apart = Meet(x, short_, Interval{0.0, 1.0}, Interval{0.0, 1.0});
    REQUIRE(apart.has_value());
    CHECK(apart->IsEmpty());
}

TEST_CASE("3D line crossings clip against exact rational parameters", "[intersect][line3d]")
{
    const Line a{Point3{}, Vec3{3.0, 0.0, 0.0}};
    const Line b{Point3{1.0, -1.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // The double 1/3 rounds DOWN, so the crossing lies just past this
    // requested segment's end. Whether it lies inside is decided on the
    // exact rational parameter; outside, the end is within the resolution
    // of the crossing, so the segments meet at the end.
    const double below = 1.0 / 3.0;
    const double above = std::nextafter(below, 1.0);
    const auto atEnd = Meet(a, b, Interval{0.0, below});
    REQUIRE(atEnd.has_value());
    CheckPoint(*atEnd, a, b, Point3{1.0, 0.0, 0.0});
    CHECK(atEnd->Points[0].ParameterA == below);
    const auto included = Meet(a, b, Interval{0.0, above});
    REQUIRE(included.has_value());
    CheckPoint(*included, a, b, Point3{1.0, 0.0, 0.0});
    CHECK(included->Points[0].ParameterA == below);
    const auto excluded = Meet(a, b, Interval{0.0, 0.3});
    REQUIRE(excluded.has_value());
    CHECK(excluded->IsEmpty());
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
    // A gap of a rounding is within the resolution: a point at the ends.
    const double short_ = std::nextafter(2.0, 0.0);
    const auto touching = Meet(x, reversed, Interval{0.0, short_}, Interval{0.0, 1.0});
    REQUIRE(touching.has_value());
    REQUIRE(touching->Points.size() == 1);
    CHECK(touching->Points[0].ParameterA == short_);
    CHECK(touching->Points[0].ParameterB == 1.0);
    const auto gap = Meet(x, reversed, Interval{0.0, 1.5}, Interval{0.0, 1.0});
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

TEST_CASE("a 3D line meets as the same segment as a degree-1 B-spline does",
          "[intersect][line3d][spline][property]")
{
    // The owner's rule: a Line and a flat B-spline are one curve, so they
    // meet alike: the same meeting, or none, or the spline fails closed
    // where its tangent-curve search cannot settle (parallel supports).
    // Segments on integer points, the second moved off by a gap around the
    // resolution, cross, touch at a vertex or a rounding apart, or miss.
    rc::prop("a line and its B-spline meet alike", [] {
        const auto integer = [] { return static_cast<double>(*rc::gen::inRange(-8, 9)); };
        const Point3 p0{integer(), integer(), integer()};
        const Point3 p1{integer(), integer(), integer()};
        const Point3 q0{integer(), integer(), integer()};
        const Point3 q1{integer(), integer(), integer()};
        RC_PRE(p0 != p1 && q0 != q1);
        const double gap = *rc::gen::element(0.0, 1e-13, 3e-6, 9.9e-6, 1.01e-5, 1e-4);
        const Vec3 off{gap, -gap, gap};
        const Line a{p0, p1 - p0};
        const Line b{q0 + off, q1 - q0};
        const Interval unit{0.0, 1.0};
        const auto lines = Meet(a, b, unit, unit);
        RC_ASSERT(lines.has_value());
        const auto splines = Intersect(Curve{Flat(a, unit)}, unit, Curve{Flat(b, unit)}, unit, ctx);
        if (!splines) {
            RC_ASSERT(splines.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        if (!lines->Overlaps.empty()) {
            return; // Parallel within the resolution; the splines found ends.
        }
        // One meeting, placed alike to within the resolution: the lines
        // place it at the segments' closest points, the splines at a domain's
        // end and the other curve's critical point continued past it, which
        // differ by about the gap at a vertex.
        RC_ASSERT(lines->Points.size() == splines->Points.size());
        for (std::size_t i = 0; i < lines->Points.size(); ++i) {
            const auto& l = lines->Points[i];
            const auto& s = splines->Points[i];
            RC_ASSERT(Distance(Eval(a, l.ParameterA), Eval(a, s.ParameterA)) <=
                      ctx.GetResolution());
            RC_ASSERT(Distance(Eval(b, l.ParameterB), Eval(b, s.ParameterB)) <=
                      ctx.GetResolution());
        }
    });
}

TEST_CASE("nearly parallel segments touching end to end meet as lines and as splines",
          "[intersect][line3d][spline]")
{
    // The supports' closest approach lies ten units behind both segments;
    // the end of one is 1.4e-6 from the start of the other.
    const Line a{Point3{}, Vec3{1.0, 0.0, 0.0}};
    const Line b{Point3{1.0 + 1e-6, 1e-6, 0.0}, Vec3{1.0, 1e-7, 0.0}};
    const Interval unit{0.0, 1.0};
    const auto lines = Meet(a, b, unit, unit);
    REQUIRE(lines.has_value());
    REQUIRE(lines->Points.size() == 1);
    const auto splines = Intersect(Curve{Flat(a, unit)}, unit, Curve{Flat(b, unit)}, unit, ctx);
    REQUIRE(splines.has_value());
    REQUIRE(splines->Points.size() == 1);
    CHECK(lines->Points[0].ParameterA == 1.0);
    CHECK(lines->Points[0].ParameterB == 0.0);
    CHECK(splines->Points[0].ParameterA == 1.0);
    CHECK(splines->Points[0].ParameterB == 0.0);
}
