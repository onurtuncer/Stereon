// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A line against the cone and the torus. As for every curve-surface pair, a
// successful empty result is a proof of separation, and a configuration the
// handler cannot resolve fails with IntersectionFailed: a line along a
// generator or through the apex over an unbounded domain, a line tangent to
// the tube. Every point lies on both geometries at its parameters, in order
// of the line's parameter, and every sign change of the surface's equation
// along a random line has a point in it.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Cone;
using Stereon::ContactType;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// The cone about z of radius 1 at its origin, opening at pi/4: its apex at
/// z = -1.
const Cone rightCone{Point3{}, 1.0, pi / 4.0, ex, ey};
/// The torus about z with radii 2 and 1.
const Torus ring{Point3{}, 2.0, 1.0, ex, ey};

template <typename S>
Result<CurveSurfaceIntersection> Meet(const Line& line, const S& surface)
{
    return Intersect(Curve{line}, Surface{surface}, ctx);
}

/// The cone's and the torus's equations at p, in their frames.
double ConeLevel(const Cone& cone, const Point3& p)
{
    const Vec3 w = p - cone.Origin;
    const double x = Dot(w, cone.XAxis);
    const double y = Dot(w, cone.YAxis);
    const double h = Dot(w, Axis(cone));
    const double c = std::cos(cone.SemiAngle);
    const double k = cone.Radius * c + h * std::sin(cone.SemiAngle);
    return c * c * (x * x + y * y) - k * k;
}

double TorusLevel(const Torus& torus, const Point3& p)
{
    const Vec3 w = p - torus.Center;
    const double x = Dot(w, torus.XAxis);
    const double y = Dot(w, torus.YAxis);
    const double big = torus.MajorRadius;
    const double small = torus.MinorRadius;
    const double k = Dot(w, w) + big * big - small * small;
    return k * k - 4.0 * big * big * (x * x + y * y);
}

/// Every point decided, a crossing, on both geometries to rounding of its
/// coordinates, in order of the line's parameter.
template <typename S>
void RequireOnBoth(const CurveSurfaceIntersection& result, const Line& line, const S& surface)
{
    double previous = -std::numeric_limits<double>::infinity();
    for (const CurveSurfacePoint& p : result.Points) {
        REQUIRE(p.Decided);
        REQUIRE(p.Contact == ContactType::Crossing);
        const double within = 1e-12 * std::max(1.0, Length(Stereon::ToVec(p.Position)));
        REQUIRE(Distance(Eval(line, p.Parameter), p.Position) <= within);
        REQUIRE(Distance(Eval(surface, p.OnSurface), p.Position) <= within);
        REQUIRE(p.Parameter >= previous);
        previous = p.Parameter;
    }
}

} // namespace

TEST_CASE("a line crosses a cone where its closed form says, on either nappe",
          "[intersect][line-surface][cone]")
{
    // Along x at height 1/2: the cone's radius there is 1 + 1/2 tan a = 3/2.
    const Line across{Point3{-5.0, 0.0, 0.5}, ex};
    const auto result = Meet(across, rightCone);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 2);
    REQUIRE(std::abs(result->Points[0].Parameter - 3.5) <= 1e-14);
    REQUIRE(std::abs(result->Points[1].Parameter - 6.5) <= 1e-14);
    RequireOnBoth(*result, across, rightCone);

    // Up x = 1/2: |1 + z| = 1/2 at z = -1/2 on the upper nappe and z = -3/2
    // on the lower, either side of the apex at z = -1.
    const Line up{Point3{0.5, 0.0, 0.0}, ez};
    const auto nappes = Meet(up, rightCone);
    REQUIRE(nappes.has_value());
    REQUIRE(nappes->Points.size() == 2);
    REQUIRE(std::abs(nappes->Points[0].Parameter + 1.5) <= 1e-14);
    REQUIRE(std::abs(nappes->Points[1].Parameter + 0.5) <= 1e-14);
    RequireOnBoth(*nappes, up, rightCone);

    // Over a half-line, only the crossing in it.
    const auto half =
        Intersect(Curve{across}, Interval{5.0, std::numeric_limits<double>::infinity()},
                  Surface{rightCone}, Domain(rightCone), ctx);
    REQUIRE(half.has_value());
    REQUIRE(half->Points.size() == 1);
    REQUIRE(std::abs(half->Points.front().Parameter - 6.5) <= 1e-14);

    // A horizontal line beside the circle the cone has at its height: apart.
    const Line beside{Point3{-5.0, 10.0, 0.5}, ex};
    const auto apart = Meet(beside, rightCone);
    REQUIRE(apart.has_value());
    REQUIRE(apart->IsEmpty());
}

TEST_CASE("a line along a generator or through the apex fails closed over an unbounded domain",
          "[intersect][line-surface][cone]")
{
    // (t, 0, t) runs at the cone's own angle: |t| = 1 + t only at t = -1/2,
    // the other root of the quadratic out at infinity, or as near it as the
    // rounding of pi/4 makes it. Over the whole line that is undecidable.
    const Line generator{Point3{}, Vec3{1.0, 0.0, 1.0}};
    const auto whole = Meet(generator, rightCone);
    REQUIRE_FALSE(whole.has_value());
    REQUIRE(whole.error().GetCode() == ErrorCode::IntersectionFailed);
    // Over a bounded domain the segment's proofs decide it: one crossing.
    const auto bounded = Intersect(Curve{generator}, Interval{-10.0, 10.0}, Surface{rightCone},
                                   Domain(rightCone), ctx);
    REQUIRE(bounded.has_value());
    REQUIRE(bounded->Points.size() == 1);
    REQUIRE(std::abs(bounded->Points.front().Parameter + 0.5) <= 1e-12);
    RequireOnBoth(*bounded, generator, rightCone);

    // Across the apex: a double root, never an empty success.
    const Line apex{Point3{0.0, 0.0, -1.0}, ex};
    const auto touch = Meet(apex, rightCone);
    REQUIRE_FALSE(touch.has_value());
    REQUIRE(touch.error().GetCode() == ErrorCode::IntersectionFailed);
}

TEST_CASE("a line crosses a torus where its closed form says", "[intersect][line-surface][torus]")
{
    // The x axis: through the tube at x = -3, -1, 1, 3.
    const Line axis{Point3{}, ex};
    const auto result = Meet(axis, ring);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 4);
    const double expected[] = {-3.0, -1.0, 1.0, 3.0};
    for (std::size_t k = 0; k < 4; ++k) {
        REQUIRE(std::abs(result->Points[k].Parameter - expected[k]) <= 1e-12);
    }
    RequireOnBoth(*result, axis, ring);

    // Down the axis, through the hole: apart. Far off: apart.
    REQUIRE(Meet(Line{Point3{}, ez}, ring)->IsEmpty());
    REQUIRE(Meet(Line{Point3{0.0, 0.0, 50.0}, ex}, ring)->IsEmpty());

    // Over a domain that holds only the outer crossing on the right.
    const auto right =
        Intersect(Curve{axis}, Interval{2.0, 10.0}, Surface{ring}, Domain(ring), ctx);
    REQUIRE(right.has_value());
    REQUIRE(right->Points.size() == 1);
    REQUIRE(std::abs(right->Points.front().Parameter - 3.0) <= 1e-12);
}

TEST_CASE("an empty torus result is a proof: 2^-30 above the tube is apart, below crosses",
          "[intersect][line-surface][torus][regression]")
{
    const double gap = std::ldexp(1.0, -30);
    // Along x at the tube's top, z = 1, the line touches it at x = -+2.
    const auto above = Meet(Line{Point3{0.0, 0.0, 1.0 + gap}, ex}, ring);
    REQUIRE(above.has_value());
    REQUIRE(above->IsEmpty());
    const auto below = Meet(Line{Point3{0.0, 0.0, 1.0 - gap}, ex}, ring);
    REQUIRE(below.has_value());
    REQUIRE(below->Points.size() == 4);
    const auto touching = Meet(Line{Point3{0.0, 0.0, 1.0}, ex}, ring);
    if (touching.has_value()) {
        REQUIRE_FALSE(touching->IsEmpty());
    } else {
        REQUIRE(touching.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}

TEST_CASE("every sign change along a random line has a point in it",
          "[intersect][line-surface][property]")
{
    // Signs are sampled in doubles over t in [-20, 20]: a sign counts beyond
    // `noise`, and a point may sit `slack` past its samples.
    constexpr double noise = 1e-9;
    constexpr double slack = 1e-9;
    const auto coordinate = [] { return *rc::gen::inRange(-24, 25) / 8.0; };
    const auto check = [](const Line& line, const auto& surface, const auto& level) {
        const auto result = Meet(line, surface);
        if (!result) {
            RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
            return;
        }
        const std::vector<double> t = IntersectTest::SamplesOf(Interval{-20.0, 20.0}, 512);
        for (std::size_t i = 0; i + 1 < t.size(); ++i) {
            const double a = level(Eval(line, t[i]));
            const double b = level(Eval(line, t[i + 1]));
            if ((a < -noise && b > noise) || (a > noise && b < -noise)) {
                RC_ASSERT(std::ranges::any_of(result->Points, [&](const CurveSurfacePoint& p) {
                    return t[i] - slack <= p.Parameter && p.Parameter <= t[i + 1] + slack;
                }));
            }
        }
        for (const CurveSurfacePoint& p : result->Points) {
            const double within = 1e-9 * (1.0 + Length(Stereon::ToVec(p.Position)));
            RC_ASSERT(Distance(Eval(line, p.Parameter), p.Position) <= within);
            RC_ASSERT(Distance(Eval(surface, p.OnSurface), p.Position) <= within);
        }
    };
    rc::prop("for random lines against a cone and a torus", [&] {
        const Point3 origin{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(Length(direction) >= 0.5);
        const Line line{origin, direction};
        const double sign = *rc::gen::arbitrary<bool>() ? 1.0 : -1.0;
        const Cone cone{Point3{0.25, -0.5, 0.125}, *rc::gen::inRange(4, 17) / 8.0,
                        sign * *rc::gen::inRange(2, 12) / 10.0, ex, ey};
        check(line, cone, [&](const Point3& p) { return ConeLevel(cone, p); });
        const double small = *rc::gen::inRange(4, 9) / 8.0;
        const Torus torus{Point3{-0.25, 0.5, 0.0}, small + *rc::gen::inRange(2, 13) / 8.0, small,
                          ex, ey};
        check(line, torus, [&](const Point3& p) { return TorusLevel(torus, p); });
    });
}
