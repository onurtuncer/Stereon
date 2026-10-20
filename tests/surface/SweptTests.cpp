// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string_view>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using namespace Stereon;

namespace {
constexpr double Tolerance = 1e-10;

BSplineCurve Profile()
{
    // A rational, non-planar cubic: neither a meridian nor a conic.
    return {
        3, {0, 0, 0, 0, 1, 1, 1, 1}, {{2, 0, -1}, {3, 1, 0}, {1, -1, 2}, {2, 0, 3}}, {1, 2, 3, 1}};
}

// Independent rational quadratic half-angle identity on a single arc.
double AngleAt(const double t, const double start, const double span)
{
    return start + 2.0 * std::atan2(t * std::sin(span / 2.0), 1.0 - t + t * std::cos(span / 2.0));
}
} // namespace

TEST_CASE("extrusion preserves rational profile parameters and mixed derivatives",
          "[surface][swept]")
{
    const auto profile = Profile();
    const Vec3 direction{1, -2, 3};
    const auto result = MakeExtrusionSurface(profile, direction, Interval{-2, 4});
    REQUIRE(result.has_value());
    const auto& surface = *result;
    REQUIRE(Domain(surface) == Box2{Interval{0, 1}, Interval{-2, 4}});
    for (const double u : {0.0, 0.17, 0.5, 0.83, 1.0}) {
        for (const double v : {-2.0, -0.3, 1.0, 4.0}) {
            REQUIRE(Distance(Eval(surface, Point2{u, v}), Eval(profile, u) + v * direction) <
                    Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 1, 0) - Deriv(profile, u, 1)) < Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 0, 1) - direction) < Tolerance);
            REQUIRE(Length(Deriv(surface, Point2{u, v}, 1, 1)) < Tolerance);
        }
    }
    const Surface erased{surface};
    REQUIRE(Distance(Eval(erased, Point2{0.4, 2}), Eval(surface, Point2{0.4, 2})) == 0.0);
}

TEST_CASE("extrusion commutes with translation", "[surface][swept][property]")
{
    rc::prop("translation of the profile translates every generated point", [] {
        const int seed = *rc::gen::inRange(-100, 101);
        auto shifted = Profile();
        const Vec3 shift{static_cast<double>(seed), -2.0, 3.0};
        for (auto& pole : shifted.Poles) {
            pole = pole + shift;
        }
        const Vec3 direction{2, 3, -1};
        const auto a = MakeExtrusionSurface(Profile(), direction, Interval{-1, 2});
        const auto b = MakeExtrusionSurface(shifted, direction, Interval{-1, 2});
        RC_ASSERT(a.has_value());
        RC_ASSERT(b.has_value());
        const double u = static_cast<double>(*rc::gen::inRange(0, 101)) / 100.0;
        RC_ASSERT(Distance(Eval(*a, Point2{u, 0.7}) + shift, Eval(*b, Point2{u, 0.7})) < Tolerance);
    });
}

TEST_CASE("revolution matches Rodrigues rotation between angular knots", "[surface][swept]")
{
    const auto profile = Profile();
    const Point3 origin{1, -2, 0.5};
    const Vec3 axis{2, -1, 3};
    const double start = 5.5; // The sweep crosses the conventional 2 pi seam.
    const double span = 1.0;
    const auto result = MakeRevolutionSurface(profile, origin, axis, Interval{start, start + span});
    REQUIRE(result.has_value());
    for (const double t : {0.0, 0.13, 0.5, 0.79, 1.0}) {
        const double angle = AngleAt(t, start, span);
        const auto rotation = MakeRotation(origin, axis, angle);
        for (const double v : {0.0, 0.23, 0.6, 1.0}) {
            REQUIRE(Distance(Eval(*result, Point2{start + t * span, v}),
                             Transformed(Eval(profile, v), rotation)) < Tolerance);
            REQUIRE(Length(Deriv(*result, Point2{start + t * span, v}, 0, 1) -
                           Transformed(Deriv(profile, v, 1), rotation)) < Tolerance);
        }
    }
}

TEST_CASE("full revolution closes and admits poles on the axis", "[surface][swept]")
{
    const BSplineCurve profile{1, {0, 0, 1, 1}, {{0, 0, 0}, {2, 0, 3}}, {}};
    const auto result =
        MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 2}, Interval{0, CirclePeriod});
    REQUIRE(result.has_value());
    for (const double v : {0.0, 0.3, 1.0}) {
        REQUIRE(Distance(Eval(*result, Point2{0, v}), Eval(*result, Point2{CirclePeriod, v})) <
                Tolerance);
    }
    for (const double u : {0.0, 0.7, 2.0, 4.0, CirclePeriod}) {
        REQUIRE(Distance(Eval(*result, Point2{u, 0}), Point3{}) < Tolerance);
        const auto p = Eval(*result, Point2{u, 0.5});
        REQUIRE(std::abs(std::hypot(p.X, p.Y) - 1.0) < Tolerance);
        REQUIRE(std::abs(p.Z - 1.5) < Tolerance);
    }
}

TEST_CASE("swept factories reject malformed and overflowing inputs", "[surface][swept]")
{
    const auto profile = Profile();
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const Vec3 direction : {Vec3{}, Vec3{nan, 0, 1}, Vec3{0, infinity, 1}}) {
        REQUIRE_FALSE(MakeExtrusionSurface(profile, direction, Interval{0, 1}));
        REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{}, direction, Interval{0, 1}));
    }
    for (const Interval range : {Interval{1, 1}, Interval{0, infinity}}) {
        REQUIRE_FALSE(MakeExtrusionSurface(profile, Vec3{0, 0, 1}, range));
        REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 1}, range));
    }
    REQUIRE_FALSE(
        MakeRevolutionSurface(profile, Point3{}, Vec3{0, 0, 1}, Interval{0, 2 * CirclePeriod}));
    REQUIRE_FALSE(MakeRevolutionSurface(profile, Point3{nan, 0, 0}, Vec3{0, 0, 1}, Interval{0, 1}));
    auto invalid = profile;
    invalid.Weights.pop_back();
    const auto rejected = MakeExtrusionSurface(invalid, Vec3{0, 0, 1}, Interval{0, 1});
    REQUIRE_FALSE(rejected);
    REQUIRE(rejected.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE_FALSE(MakeRevolutionSurface(invalid, Point3{}, Vec3{0, 0, 1}, Interval{0, 1}));
    const double huge = std::numeric_limits<double>::max();
    REQUIRE_FALSE(MakeExtrusionSurface(profile, Vec3{huge, 0, 0}, Interval{0, 2}));
}

// ---------------------------------------------------------------------------
// The wrapper types: unbounded sweeps evaluated from the profile itself

namespace {

/// A curve type that is not built in: a line under another name, with no
/// NURBS form the kernel knows.
struct UserLine {
    Line Base;
};
Point3 Eval(const UserLine& c, const double t)
{
    return Eval(c.Base, t);
}
Vec3 Deriv(const UserLine& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Interval Domain(const UserLine& c)
{
    return Domain(c.Base);
}
Box3 BoundingBox(const UserLine& c)
{
    return BoundingBox(c.Base);
}
double Project(const UserLine& c, const Point3& p)
{
    return Project(c.Base, p);
}
UserLine Reversed(const UserLine& c)
{
    return {Reversed(c.Base)};
}
UserLine Transformed(const UserLine& c, const Transform& x)
{
    return {Transformed(c.Base, x)};
}

Curve CircleProfile()
{
    // In the xz-plane, centred off the z axis: revolved about z it is a torus.
    return Curve(Circle{Point3{3.0, 0.0, 0.5}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}});
}

/// The least distance from p to the surface over a grid of the box.
double SampledDistance(const Surface& s, const Box2& box, const Point3& p, const int n = 80)
{
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            const Point2 uv{box.X.GetLower() + (box.X.GetUpper() - box.X.GetLower()) * i / n,
                            box.Y.GetLower() + (box.Y.GetUpper() - box.Y.GetLower()) * j / n};
            best = std::min(best, Distance(Eval(s, uv), p));
        }
    }
    return best;
}

/// Central differences of Eval, the referee for Deriv.
Vec3 NumericalPartial(const Surface& s, const Point2& uv, const bool alongU)
{
    const double h = 1e-5;
    const Point2 a = alongU ? Point2{uv.X - h, uv.Y} : Point2{uv.X, uv.Y - h};
    const Point2 b = alongU ? Point2{uv.X + h, uv.Y} : Point2{uv.X, uv.Y + h};
    return (Eval(s, b) - Eval(s, a)) / (2.0 * h);
}

constexpr double Infinite = std::numeric_limits<double>::infinity();

} // namespace

TEST_CASE("an extrusion evaluates its profile moved along the direction", "[surface][swept]")
{
    const Curve profile(Profile());
    const Vec3 direction{1.0, -2.0, 3.0};
    const ExtrusionSurface s = *MakeExtrusion(profile, direction);
    REQUIRE(Domain(s).X == Domain(profile));
    REQUIRE(Domain(s).Y.GetLower() == -Infinite);
    REQUIRE(Domain(s).Y.GetUpper() == Infinite);
    REQUIRE(Periods(s) == Vec2{});
    for (const double u : {0.0, 0.17, 0.5, 0.83, 1.0}) {
        for (const double v : {-1000.0, -0.3, 1.0, 4.0, 2500.0}) {
            const Point2 uv{u, v};
            REQUIRE(Distance(Eval(s, uv), Eval(profile, u) + v * direction) < Tolerance);
            REQUIRE(Deriv(s, uv, 1, 0) == Deriv(profile, u, 1));
            REQUIRE(Deriv(s, uv, 2, 0) == Deriv(profile, u, 2));
            REQUIRE(Deriv(s, uv, 0, 1) == direction);
            REQUIRE(Deriv(s, uv, 1, 1) == Vec3{});
            REQUIRE(Deriv(s, uv, 0, 2) == Vec3{});
        }
    }
    const Box3 box = BoundingBox(s);
    REQUIRE(box.X.GetLower() == -Infinite);
    REQUIRE(box.Z.GetUpper() == Infinite);
    const ExtrusionSurface circular = *MakeExtrusion(CircleProfile(), Vec3{0.0, 1.0, 0.0});
    REQUIRE(Periods(circular) == Vec2{CirclePeriod, 0.0});
    REQUIRE(BoundingBox(circular).X == BoundingBox(CircleProfile()).X);
    const Surface erased(s);
    REQUIRE(Eval(erased, Point2{0.4, -7.0}) == Eval(s, Point2{0.4, -7.0}));
    REQUIRE(Periods(erased) == Vec2{});
    STATIC_REQUIRE(GeometryTypeName<ExtrusionSurface> == std::string_view{"ExtrusionSurface"});
}

TEST_CASE("an extrusion projects onto the surface", "[surface][swept][property]")
{
    rc::prop("exactly for a profile with a NURBS form, as a plane for a line", [] {
        const int kind = *rc::gen::inRange(0, 3);
        const Curve profile = kind == 0 ? Curve(Profile())
                              : kind == 1
                                  ? CircleProfile()
                                  : Curve(Line{Point3{1.0, 2.0, 0.0}, Vec3{0.5, -1.0, 0.2}});
        const Vec3 direction{*rc::gen::inRange(-30, 30) / 10.0, *rc::gen::inRange(-30, 30) / 10.0,
                             1.0 + *rc::gen::inRange(0, 30) / 10.0};
        const ExtrusionSurface s = *MakeExtrusion(profile, direction);
        const Point3 p{*rc::gen::inRange(-60, 60) / 10.0, *rc::gen::inRange(-60, 60) / 10.0,
                       *rc::gen::inRange(-60, 60) / 10.0};
        const Point2 found = Project(s, p);
        RC_ASSERT(Contains(Domain(s), found));
        const double distance = Distance(Eval(s, found), p);
        const Interval uRange = kind == 2 ? Interval{-8.0, 8.0} : Domain(profile);
        const Box2 box{uRange, Interval{-8.0, 8.0}};
        RC_ASSERT(distance <= SampledDistance(Surface(s), box, p) + 2e-3);
        // The displacement to the closest point is normal to the direction.
        const Vec3 r = p - Eval(s, found);
        RC_ASSERT(std::abs(Dot(r, direction)) <= 1e-7 * (1.0 + Length(r)) * Length(direction));
    });
}

TEST_CASE("a line profile extrudes to a plane", "[surface][swept]")
{
    const Line line{Point3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Vec3 direction{0.0, 0.0, 2.0};
    const ExtrusionSurface s = *MakeExtrusion(Curve(line), direction);
    REQUIRE(s.GetProjection() == nullptr); // no finite domain, no NURBS form
    const Plane plane{Point3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    for (const Point3& p : {Point3{4.0, 3.0, -7.0}, Point3{-2.0, 100.0, 0.3}}) {
        const Point2 found = Project(s, p);
        REQUIRE(Distance(Eval(s, found), p) <= Distance(Eval(plane, Project(plane, p)), p) + 1e-12);
        REQUIRE(found == Point2{p.Y, p.Z / 2.0});
    }
    // A line along the direction extrudes to itself.
    const ExtrusionSurface degenerate = *MakeExtrusion(Curve(line), Vec3{0.0, 3.0, 0.0});
    const Point2 onLine = Project(degenerate, Point3{5.0, 2.0, 0.0});
    REQUIRE(Distance(Eval(degenerate, onLine), Point3{1.0, 2.0, 0.0}) <= 1e-12);
}

TEST_CASE("a profile without a NURBS form projects onto the surface, approximately",
          "[surface][swept]")
{
    const Curve profile(UserLine{Line{Point3{}, Vec3{1.0, 1.0, 0.0}}});
    const Vec3 up{0.0, 0.0, 1.0};
    const ExtrusionSurface s = *MakeExtrusion(profile, up);
    REQUIRE(s.GetProjection() == nullptr);
    const Point3 p{3.0, 1.0, 4.0};
    const Point2 found = Project(s, p);
    // On the surface, at p's height: the profile's own closest point moved
    // there, which is what the approximation is.
    REQUIRE(found.Y == 4.0);
    REQUIRE(Distance(Eval(s, found), Eval(profile, Project(profile, p)) + 4.0 * up) <= 1e-12);
}

TEST_CASE("a revolution turns its profile about the axis", "[surface][swept]")
{
    const Curve profile = CircleProfile();
    const RevolutionSurface s =
        *MakeRevolution(profile, Point3{0.0, 0.0, 0.5}, Vec3{0.0, 0.0, 2.0});
    REQUIRE(s.GetAxis() == Vec3{0.0, 0.0, 1.0});
    REQUIRE(Domain(s) == Box2{Interval{0.0, CirclePeriod}, Domain(profile)});
    REQUIRE(Periods(s) == Vec2{CirclePeriod, CirclePeriod});
    // The torus it is.
    const Torus torus{Point3{0.0, 0.0, 0.5}, 3.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Surface erased(s);
    for (const double u : {0.0, 0.7, 2.0, 4.5, CirclePeriod}) {
        for (const double v : {0.0, 1.1, 3.0, 5.9}) {
            const Point2 uv{u, v};
            const Point3 p = Eval(s, uv);
            REQUIRE(Distance(Eval(torus, Project(torus, p)), p) <= 1e-9);
            REQUIRE(Length(Deriv(s, uv, 1, 0) - NumericalPartial(erased, uv, true)) <= 1e-7);
            REQUIRE(Length(Deriv(s, uv, 0, 1) - NumericalPartial(erased, uv, false)) <= 1e-7);
            // S_uu points back to the axis: minus the radial part.
            REQUIRE(Length(Deriv(s, uv, 2, 0) + (p - Point3{0.0, 0.0, p.Z})) <= 1e-9);
            REQUIRE(Distance(Eval(s, Point2{u + CirclePeriod, v}), p) <= 1e-12);
        }
    }
    // Right-handed about the axis: a quarter turn takes +x towards +y.
    REQUIRE(Eval(s, Point2{std::numbers::pi / 2.0, 0.0}).Y > 3.9);
    STATIC_REQUIRE(GeometryTypeName<RevolutionSurface> == std::string_view{"RevolutionSurface"});
}

TEST_CASE("a revolution projects onto its full turn", "[surface][swept][property]")
{
    rc::prop("exactly for a profile with a NURBS form", [] {
        const Curve profile = *rc::gen::arbitrary<bool>() ? CircleProfile() : Curve(Profile());
        const Point3 origin{*rc::gen::inRange(-10, 10) / 10.0, 0.0,
                            *rc::gen::inRange(-10, 10) / 10.0};
        const Vec3 axis{*rc::gen::inRange(-5, 5) / 10.0, *rc::gen::inRange(-5, 5) / 10.0, 1.0};
        const RevolutionSurface s = *MakeRevolution(profile, origin, axis);
        RC_ASSERT(s.GetPiece() != nullptr);
        const Point3 p{*rc::gen::inRange(-60, 60) / 10.0, *rc::gen::inRange(-60, 60) / 10.0,
                       *rc::gen::inRange(-60, 60) / 10.0};
        const Point2 found = Project(s, p);
        RC_ASSERT(Contains(Domain(s), found));
        const double distance = Distance(Eval(s, found), p);
        RC_ASSERT(distance <= SampledDistance(Surface(s), Domain(s), p) + 5e-3);
    });
}

TEST_CASE("a profile without a NURBS form revolves and projects approximately", "[surface][swept]")
{
    const Curve profile(UserLine{Line{Point3{2.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}});
    const RevolutionSurface s = *MakeRevolution(profile, Point3{}, Vec3{0.0, 0.0, 1.0});
    REQUIRE(s.GetPiece() == nullptr);
    // A cylinder of radius 2 about z: the projection of (4, 0, 3) is (2, 0, 3).
    const Point2 found = Project(s, Point3{4.0, 0.0, 3.0});
    REQUIRE(Distance(Eval(s, found), Point3{2.0, 0.0, 3.0}) <= 1e-9);
    const Point3 inside{-1.0, -1.0, -2.0};
    const Point2 behind = Project(s, inside);
    REQUIRE(std::abs(Distance(Eval(s, behind), inside) - (2.0 - std::sqrt(2.0))) <= 1e-9);
    REQUIRE(BoundingBox(s).Z.GetLower() == -Infinite); // the line is unbounded
}

TEST_CASE("sweeps reverse, move and reflect consistently", "[surface][swept]")
{
    const Curve profile(Profile());
    const ExtrusionSurface e = *MakeExtrusion(profile, Vec3{1.0, -2.0, 3.0});
    const RevolutionSurface r =
        *MakeRevolution(profile, Point3{0.5, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0});
    const Transform rotation = MakeRotation(Vec3{1.0, 2.0, 0.5}, 0.7);
    const Transform mirror = MakePointMirror(Point3{1.0, -1.0, 2.0});
    for (const double t : {0.1, 0.6, 0.95}) {
        for (const double w : {0.3, 2.4, 5.0}) {
            // Extrusion: the profile runs backwards, so its parameter t maps
            // to 1 - t for this profile; the normal flips.
            const Point2 uv{t, w};
            const Point2 reflected{1.0 - t, w};
            REQUIRE(Distance(Eval(Reversed(e), reflected), Eval(e, uv)) <= 1e-9);
            REQUIRE(Dot(*UnitNormal(Reversed(e), reflected), *UnitNormal(e, uv)) <= -1.0 + 1e-9);
            for (const Transform& x : {rotation, mirror}) {
                REQUIRE(Distance(Eval(Transformed(e, x), uv), Transformed(Eval(e, uv), x)) <= 1e-9);
            }
            // Revolution: the angle w, the profile parameter t; reversed, the
            // turn runs the other way.
            const Point2 st{w, t};
            const Point2 turnedBack{CirclePeriod - w, t};
            REQUIRE(Distance(Eval(Reversed(r), turnedBack), Eval(r, st)) <= 1e-9);
            REQUIRE(Dot(*UnitNormal(Reversed(r), turnedBack), *UnitNormal(r, st)) <= -1.0 + 1e-9);
            for (const Transform& x : {rotation, mirror}) {
                REQUIRE(Distance(Eval(Transformed(r, x), st), Transformed(Eval(r, st), x)) <= 1e-9);
            }
        }
    }
}

TEST_CASE("ToBSpline of a sweep is the exact patch", "[surface][swept]")
{
    const Curve profile(Profile());
    const ExtrusionSurface e = *MakeExtrusion(profile, Vec3{1.0, -2.0, 3.0});
    const Box2 box{Interval{0.2, 0.8}, Interval{-2.0, 4.0}};
    const BSplineSurface patch = *ToBSpline(e, box);
    REQUIRE(Domain(patch) == box);
    for (const double u : {0.2, 0.5, 0.8}) {
        for (const double v : {-2.0, 1.0, 4.0}) {
            REQUIRE(Distance(Eval(patch, Point2{u, v}), Eval(e, Point2{u, v})) <= 1e-9);
        }
    }
    REQUIRE(ToBSpline(e, Box2{Interval{0.0, 1.5}, Interval{0.0, 1.0}}).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(ToBSpline(e, Box2{Interval{0.0, 1.0}, Interval{0.0, Infinite}}).error().GetCode() ==
            ErrorCode::InvalidInput);
    const RevolutionSurface r = *MakeRevolution(profile, Point3{}, Vec3{0.0, 0.0, 1.0});
    const Box2 part{Interval{0.0, 2.0}, Interval{0.0, 1.0}};
    const BSplineSurface turn = *ToBSpline(r, part);
    for (const double a : {0.0, 0.9, 2.0}) {
        for (const double t : {0.0, 0.4, 1.0}) {
            const Point3 p = Eval(turn, Point2{a, t});
            REQUIRE(Distance(Eval(r, Project(r, p)), p) <= 1e-8);
        }
    }
    REQUIRE(ToBSpline(r, Box2{Interval{0.0, CirclePeriod + 0.1}, Interval{0.0, 1.0}})
                .error()
                .GetCode() == ErrorCode::InvalidInput);
    // The erased conversion reaches the wrappers.
    REQUIRE(ToBSpline(Surface(e), box).has_value());
    REQUIRE(ToBSpline(Surface(r), part).has_value());
}

TEST_CASE("MakeExtrusion and MakeRevolution validate their data", "[surface][swept][validation]")
{
    const auto invalid = [](const auto& result) {
        return !result.has_value() && result.error().GetCode() == ErrorCode::InvalidInput;
    };
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const Curve profile(Profile());
    REQUIRE(invalid(MakeExtrusion(Curve{}, Vec3{0.0, 0.0, 1.0})));
    REQUIRE(invalid(MakeExtrusion(profile, Vec3{})));
    REQUIRE(invalid(MakeExtrusion(profile, Vec3{nan, 0.0, 1.0})));
    // Tiny but usable: 1e-300 would square to zero and is refused as zero.
    REQUIRE(MakeExtrusion(profile, Vec3{0.0, 0.0, 1e-150}).has_value());
    REQUIRE(invalid(MakeRevolution(Curve{}, Point3{}, Vec3{0.0, 0.0, 1.0})));
    REQUIRE(invalid(MakeRevolution(profile, Point3{nan, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0})));
    REQUIRE(invalid(MakeRevolution(profile, Point3{}, Vec3{})));
    REQUIRE(invalid(MakeRevolution(profile, Point3{}, Vec3{0.0, nan, 1.0})));
    REQUIRE(MakeRevolution(profile, Point3{}, Vec3{0.0, 0.0, -5.0}).has_value());
}

// ---------------------------------------------------------------------------
// Audit 2026-10-20: trimmed profiles (F2) and the numerical domain (F6)

TEST_CASE("an oblique extrusion of a trimmed line projects exactly", "[surface][swept][audit]")
{
    // F2: the profile C(u) = (u, 0, 0) trimmed to [0, 1], extruded along
    // d = (1, 1, 0). The point (0.5, 0.5, 0) lies on the surface at
    // (u, v) = (0, 0.5). Without the trimmed profile's exact form,
    // projection took the profile's own closest point, u = 0.5, and
    // returned (0.75, 0.25, 0): 0.35 from a point on the surface.
    const Curve profile(
        *MakeTrimmedCurve(Curve(Line{Point3{}, Vec3{1.0, 0.0, 0.0}}), Interval{0.0, 1.0}));
    const ExtrusionSurface s = *MakeExtrusion(profile, Vec3{1.0, 1.0, 0.0});
    REQUIRE(s.GetProjection() != nullptr);
    const Point3 p{0.5, 0.5, 0.0};
    const Point2 found = Project(s, p);
    REQUIRE(std::abs(found.X) <= 1e-12);
    REQUIRE(std::abs(found.Y - 0.5) <= 1e-12);
    REQUIRE(Distance(Eval(s, found), p) <= 1e-12);

    // Its exact patch, and a trimmed surface of it, which builds that patch.
    const Box2 box{Interval{0.0, 1.0}, Interval{-1.0, 1.0}};
    const BSplineSurface patch = *ToBSpline(s, box);
    REQUIRE(Distance(Eval(patch, Point2{0.0, 0.5}), p) <= 1e-12);
    const TrimmedSurface band = *MakeTrimmedSurface(Surface(s), box);
    REQUIRE(band.GetPatch() != nullptr);
    REQUIRE(ToBSpline(band).has_value());
    REQUIRE(Distance(Eval(band, Project(band, p)), p) <= 1e-12);
}

namespace {

/// Trimmed profiles of every kind of basis, forwards and reversed, and a
/// trimmed curve of a trimmed curve.
std::vector<Curve> TrimmedProfiles()
{
    const Curve line(Line{Point3{1.0, 0.0, 0.0}, Vec3{0.3, 0.2, 1.0}});
    const Curve spline(Profile());
    return {
        Curve(*MakeTrimmedCurve(line, Interval{-1.0, 2.0})),
        Curve(MakeTrimmedCurve(line, Interval{-1.0, 2.0})->Reversed()),
        Curve(*MakeTrimmedCurve(CircleProfile(), Interval{0.5, 2.5})),
        Curve(MakeTrimmedCurve(CircleProfile(), Interval{5.0, 7.0})->Reversed()),
        Curve(*MakeTrimmedCurve(spline, Interval{0.2, 0.7})),
        Curve(*MakeTrimmedCurve(Curve(*MakeTrimmedCurve(spline, Interval{0.1, 0.9})),
                                Interval{0.2, 0.6})),
    };
}

/// Points of the surface over a grid of the box project back onto it.
template <typename S>
void RequireProjectsBack(const S& s, const Box2& box, const double tolerance)
{
    for (int i = 0; i <= 6; ++i) {
        for (int j = 0; j <= 4; ++j) {
            const Point2 uv{box.X.GetLower() + Width(box.X) * i / 6.0,
                            box.Y.GetLower() + Width(box.Y) * j / 4.0};
            const Point3 p = Eval(s, uv);
            INFO("uv = (" << uv.X << ", " << uv.Y << ")");
            REQUIRE(Distance(Eval(s, Project(s, p)), p) <= tolerance);
        }
    }
}

} // namespace

TEST_CASE("sweeps of trimmed profiles have exact forms and project onto the surface",
          "[surface][swept][audit]")
{
    // F2: an extrusion or revolution of a trimmed profile keeps the
    // profile's exact NURBS form, so it projects exactly, converts to a
    // patch and can be trimmed in turn.
    int index = 0;
    for (const Curve& profile : TrimmedProfiles()) {
        INFO("profile " << index++);
        const Interval u = Domain(profile);
        const ExtrusionSurface e = *MakeExtrusion(profile, Vec3{0.4, -0.3, 1.0});
        REQUIRE(e.GetProjection() != nullptr);
        RequireProjectsBack(e, Box2{u, Interval{-1.5, 2.0}}, 1e-9);
        REQUIRE(ToBSpline(e, Box2{u, Interval{-1.0, 1.0}}).has_value());
        REQUIRE(MakeTrimmedSurface(Surface(e), Box2{u, Interval{0.0, 1.0}})->GetPatch() != nullptr);

        const RevolutionSurface r = *MakeRevolution(profile, Point3{}, Vec3{0.0, 0.0, 1.0});
        REQUIRE(r.GetPiece() != nullptr);
        RequireProjectsBack(r, Box2{Interval{0.1, 6.0}, u}, 1e-8);
        const Box2 part{Interval{0.0, 2.0}, u};
        REQUIRE(ToBSpline(r, part).has_value());
        REQUIRE(MakeTrimmedSurface(Surface(r), part)->GetPatch() != nullptr);
    }
}

TEST_CASE("a revolution's axis is normalised at any scale", "[surface][swept][audit]")
{
    // F6: Dot(axis, axis) overflowed for an axis of 1e200, which then
    // normalised to zero, and underflowed for one of 1e-200, which was
    // refused as zero although it has a direction. A power-of-two scale
    // changes nothing, bit for bit, and an ordinary axis normalises as
    // before.
    const Curve profile = CircleProfile();
    const Vec3 axis{2.0, -1.0, 3.0};
    const RevolutionSurface reference = *MakeRevolution(profile, Point3{}, axis);
    REQUIRE(reference.GetAxis() == Normalized(axis));
    for (const int e : {-1070, -1000, -600, -100, 1, 100, 600, 1021}) {
        INFO("scaled by 2^" << e);
        const Vec3 scaled{std::ldexp(axis.X, e), std::ldexp(axis.Y, e), std::ldexp(axis.Z, e)};
        const auto s = MakeRevolution(profile, Point3{}, scaled);
        REQUIRE(s.has_value());
        REQUIRE(s->GetAxis() == reference.GetAxis());
        REQUIRE(Eval(*s, Point2{0.7, 1.1}) == Eval(reference, Point2{0.7, 1.1}));
    }
    for (const double scale : {1e200, 1e-200, 1e-310, std::numeric_limits<double>::max()}) {
        INFO("scale " << scale);
        const auto up = MakeRevolution(profile, Point3{}, Vec3{0.0, 0.0, scale});
        REQUIRE(up.has_value());
        REQUIRE(up->GetAxis() == Vec3{0.0, 0.0, 1.0});
        const auto back = MakeRevolution(profile, Point3{}, Vec3{-scale, 0.0, 0.0});
        REQUIRE(back.has_value());
        REQUIRE(back->GetAxis() == Vec3{-1.0, 0.0, 0.0});
    }
}

TEST_CASE("an extrusion refuses a direction whose squared length is not a finite double",
          "[surface][swept][audit]")
{
    // F6: the extrusion keeps its direction as given and divides by its
    // squared length, which overflowed for a direction of 1e200 and turned
    // projection into nonsense; it is refused, as one of 1e-200 always was.
    const Curve profile(Profile());
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    REQUIRE(code(MakeExtrusion(profile, Vec3{1e200, 0.0, 0.0})) == ErrorCode::InvalidInput);
    REQUIRE(code(MakeExtrusion(profile, Vec3{0.0, -1e200, 1.0})) == ErrorCode::InvalidInput);
    REQUIRE(code(MakeExtrusion(profile, Vec3{0.0, 1e-200, 0.0})) == ErrorCode::InvalidInput);
    const ExtrusionSurface large = *MakeExtrusion(profile, Vec3{0.0, 0.0, 1e150});
    const Point3 p = Eval(large, Point2{0.4, 1e-150});
    REQUIRE(Distance(Eval(large, Project(large, p)), p) <= 1e-9);
}
