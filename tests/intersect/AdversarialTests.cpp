// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Adversarial properties shared by every curve–surface handler, over the
// whole grid of lines, circles, ellipses and polynomial and rational
// B-splines against planes, spheres, cylinders, cones, tori and ellipsoids
// (ADR-0048). A correct
// count can still carry wrong geometry marked certified (the audit of
// 2026-11-05: two crossings of a far line placed at one point), so every
// result is validated as geometry, not counted:
//
// - every point lies on both geometries to the rounding of its own size,
//   inside the requested domain and box, in order of the curve's
//   parameter; every overlap's pcurve traces the curve over its range;
// - a handler that cannot decide fails with IntersectionFailed, never with
//   another error and never with an unproved empty result.
//
// The inputs are multiples of 1/8 in axis-aligned frames, so translating
// by a multiple of 1/8 and scaling by a power of two transform them
// exactly: the exact problem is the transformed one, every exact decision
// must agree, and the points must move with the geometry. Reversing the
// curve keeps the points; moving a line's origin far along it shifts its
// parameters; a periodic surface's box a turn on gives the same points a
// turn on; and a domain of one parameter away from every point is a
// proved separation.

#include "RapidCheck.hpp"

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

using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::Curve;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfaceOverlap;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Result;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double eps = std::numeric_limits<double>::epsilon();

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// Axis-aligned right-handed frames, (X, Y) with X x Y the third axis.
const std::vector<std::pair<Vec3, Vec3>> frames{{ex, ey}, {ey, ez}, {ez, ex},
                                                {ey, ex}, {ez, ey}, {ex, ez}};

enum class CurveKind : int { Line, Circle, Ellipse, Cubic, Rational };
enum class SurfaceKind : int { Plane, Sphere, Cylinder, Cone, Torus, Ellipsoid };

const char* NameOf(const CurveKind k)
{
    constexpr std::array names{"Line", "Circle", "Ellipse", "Cubic", "Rational"};
    return names[static_cast<std::size_t>(k)];
}

const char* NameOf(const SurfaceKind k)
{
    constexpr std::array names{"Plane", "Sphere", "Cylinder", "Cone", "Torus", "Ellipsoid"};
    return names[static_cast<std::size_t>(k)];
}

/// An integer in [lo, hi), uniform whatever RapidCheck's size: inRange
/// alone favours the low end while the size is small, which left B-splines
/// and tori with one case in a hundred.
int Uniform(const int lo, const int hi)
{
    return *rc::gen::resize(rc::kNominalSize, rc::gen::inRange(lo, hi));
}

/// A coordinate, a multiple of 1/8 in [-3, 3].
double Coordinate()
{
    return Uniform(-24, 25) / 8.0;
}

Point3 SomePoint()
{
    return {Coordinate(), Coordinate(), Coordinate()};
}

/// A length, a multiple of 1/8 in [lo, hi] eighths.
double Length8(const int lo, const int hi)
{
    return Uniform(lo, hi + 1) / 8.0;
}

std::pair<Vec3, Vec3> SomeFrame()
{
    return frames[static_cast<std::size_t>(Uniform(0, 6))];
}

Curve SomeCurve(const CurveKind kind)
{
    switch (kind) {
    case CurveKind::Line: {
        const Vec3 d{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(d) > 0.0);
        return Curve{Line{SomePoint(), d}};
    }
    case CurveKind::Circle: {
        const auto [x, y] = SomeFrame();
        return Curve{Circle{SomePoint(), Length8(4, 20), x, y}};
    }
    case CurveKind::Ellipse: {
        const auto [x, y] = SomeFrame();
        const double minor = Length8(4, 16);
        return Curve{Ellipse{SomePoint(), minor + Length8(1, 8), minor, x, y}};
    }
    case CurveKind::Cubic:
    case CurveKind::Rational: {
        std::vector<Point3> poles;
        for (int i = 0; i < 5; ++i) {
            poles.push_back(SomePoint());
        }
        std::vector<double> weights;
        if (kind == CurveKind::Rational) {
            for (int i = 0; i < 5; ++i) {
                weights.push_back(Length8(4, 16));
            }
        }
        auto made = Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0},
                                              poles, weights);
        RC_PRE(made.has_value());
        return Curve{*made};
    }
    }
    return Curve{};
}

Surface SomeSurface(const SurfaceKind kind)
{
    const auto [x, y] = SomeFrame();
    const Point3 o = SomePoint();
    switch (kind) {
    case SurfaceKind::Plane:
        return Surface{Plane{o, x, y}};
    case SurfaceKind::Sphere:
        return Surface{Sphere{o, Length8(4, 20), x, y}};
    case SurfaceKind::Cylinder:
        return Surface{Cylinder{o, Length8(4, 20), x, y}};
    case SurfaceKind::Cone: {
        const double angle =
            std::array{pi / 6.0, pi / 4.0, pi / 3.0}[static_cast<std::size_t>(Uniform(0, 3))];
        return Surface{Cone{o, Length8(0, 12), angle, x, y}};
    }
    case SurfaceKind::Torus: {
        const double minor = Length8(4, 8);
        return Surface{Torus{o, minor + Length8(4, 16), minor, x, y}};
    }
    case SurfaceKind::Ellipsoid:
        return Surface{Ellipsoid{o, Length8(4, 20), Length8(4, 20), Length8(4, 20), x, y}};
    }
    return Surface{};
}

/// How many kinds of surface the grid has.
constexpr int SurfaceKinds = 6;

/// One configuration of the grid, named for the counterexample.
struct Pair {
    CurveKind CurveType;
    SurfaceKind SurfaceType;
    Curve C;
    Surface S;
};

Pair SomePair()
{
    const auto c = static_cast<CurveKind>(Uniform(0, 5));
    const auto s = static_cast<SurfaceKind>(Uniform(0, SurfaceKinds));
    RC_TAG(std::string(NameOf(c)) + " x " + NameOf(s));
    return {c, s, SomeCurve(c), SomeSurface(s)};
}

/// The size against which a point's rounding is measured: its distance
/// from the origin, or the run of the curve's parameter to it, whichever
/// is larger. A line's point carries the rounding of origin + t d.
double ScaleAt(const Curve& curve, const double t, const Point3& p)
{
    return std::max({1.0, Length(Stereon::ToVec(p)), std::abs(t) * Length(Deriv(curve, t, 1))});
}

/// A few hundred units in the last place of `scale`: the agreement a placed
/// point must reach.
double Rounding(const double scale)
{
    return 512.0 * eps * scale;
}

/// The result is either geometry that validates or IntersectionFailed.
/// Returns whether there is a result to compare further.
bool Validated(const Result<CurveSurfaceIntersection>& result, const Curve& curve,
               const Interval& domain, const Surface& surface, const Box2& box)
{
    if (!result) {
        RC_LOG() << "failed: " << IntersectTest::ErrorText(result) << "\n";
        RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
        return false;
    }
    double previous = -std::numeric_limits<double>::infinity();
    for (const CurveSurfacePoint& p : result->Points) {
        const double within = Rounding(ScaleAt(curve, p.Parameter, p.Position));
        RC_ASSERT(Contains(domain, p.Parameter));
        RC_ASSERT(Contains(box.X, p.OnSurface.X));
        RC_ASSERT(Contains(box.Y, p.OnSurface.Y));
        RC_ASSERT(p.Parameter >= previous);
        previous = p.Parameter;
        RC_ASSERT(Distance(Eval(curve, p.Parameter), p.Position) <= within);
        RC_ASSERT(Distance(Eval(surface, p.OnSurface), p.Position) <= within);
    }
    for (const CurveSurfaceOverlap& o : result->Overlaps) {
        RC_ASSERT(!o.Range.IsEmpty());
        RC_ASSERT(Contains(domain, o.Range));
        for (const double t : IntersectTest::SamplesOf(o.Range, 8)) {
            const Point3 on = Eval(curve, t);
            RC_ASSERT(Distance(Eval(surface, Eval(o.OnSurface, t)), on) <=
                      Rounding(ScaleAt(curve, t, on)));
        }
    }
    return true;
}

Result<CurveSurfaceIntersection> Meet(const Curve& c, const Surface& s)
{
    return Intersect(c, Domain(c), s, Domain(s), ctx);
}

/// Two results of transformed problems agree: the same points, contacts
/// and overlaps, each point's position mapped by `map` and its parameter
/// by `parameter`, within the rounding of both.
template <typename Map, typename Parameter>
void RequireSame(const CurveSurfaceIntersection& a, const CurveSurfaceIntersection& b,
                 const Curve& curveA, const Curve& curveB, Map&& map, Parameter&& parameter)
{
    RC_ASSERT(a.Points.size() == b.Points.size());
    RC_ASSERT(a.Overlaps.size() == b.Overlaps.size());
    for (std::size_t i = 0; i < a.Points.size(); ++i) {
        const CurveSurfacePoint& p = a.Points[i];
        const CurveSurfacePoint& q = b.Points[i];
        RC_ASSERT(p.Contact == q.Contact);
        RC_ASSERT(p.Decided == q.Decided);
        const Point3 mapped = map(p.Position);
        const double within = Rounding(ScaleAt(curveA, p.Parameter, p.Position)) +
                              Rounding(ScaleAt(curveB, q.Parameter, q.Position));
        RC_ASSERT(Distance(mapped, q.Position) <= within);
        RC_ASSERT(std::abs(parameter(p.Parameter) - q.Parameter) <=
                  1e-9 * std::max(1.0, std::abs(q.Parameter)));
    }
}

} // namespace

TEST_CASE("every curve-surface result validates as geometry or fails closed",
          "[intersect][adversarial][property]")
{
    rc::prop("over the grid of curves and surfaces", [] {
        const Pair p = SomePair();
        const auto result = Meet(p.C, p.S);
        if (Validated(result, p.C, Domain(p.C), p.S, Domain(p.S))) {
            RC_CLASSIFY(!result->Points.empty(), "points");
            RC_CLASSIFY(!result->Overlaps.empty(), "overlaps");
            RC_CLASSIFY(result->IsEmpty(), "apart");
        } else {
            RC_CLASSIFY(true, "failed closed");
        }
    });
}

TEST_CASE("translating the pair translates every point", "[intersect][adversarial][property]")
{
    rc::prop("by a multiple of 1/8", [] {
        const Pair p = SomePair();
        const Vec3 offset{Uniform(-256, 257) / 8.0, Uniform(-256, 257) / 8.0,
                          Uniform(-256, 257) / 8.0};
        const Stereon::Transform x = Stereon::MakeTranslation(offset);
        const Curve c = Transformed(p.C, x);
        const Surface s = Transformed(p.S, x);
        const auto before = Meet(p.C, p.S);
        const auto after = Meet(c, s);
        const bool a = Validated(before, p.C, Domain(p.C), p.S, Domain(p.S));
        const bool b = Validated(after, c, Domain(c), s, Domain(s));
        if (a && b) {
            RequireSame(
                *before, *after, p.C, c, [&](const Point3& q) { return q + offset; },
                [](const double t) { return t; });
        }
    });
}

TEST_CASE("scaling the pair by a power of two scales every point",
          "[intersect][adversarial][property]")
{
    rc::prop("by 2^k, |k| <= 8", [] {
        const Pair p = SomePair();
        const double factor = std::ldexp(1.0, Uniform(-8, 9));
        const Stereon::Transform x = Stereon::MakeScaling(factor);
        const Curve c = Transformed(p.C, x);
        const Surface s = Transformed(p.S, x);
        const auto before = Meet(p.C, p.S);
        const auto after = Meet(c, s);
        const bool a = Validated(before, p.C, Domain(p.C), p.S, Domain(p.S));
        const bool b = Validated(after, c, Domain(c), s, Domain(s));
        if (a && b) {
            RequireSame(
                *before, *after, p.C, c,
                [&](const Point3& q) { return Point3{} + factor * Stereon::ToVec(q); },
                [](const double t) { return t; });
        }
    });
}

TEST_CASE("reversing the curve keeps every point", "[intersect][adversarial][property]")
{
    rc::prop("over the grid of curves and surfaces", [] {
        const Pair p = SomePair();
        const Curve r = Reversed(p.C);
        const auto forward = Meet(p.C, p.S);
        const auto backward = Meet(r, p.S);
        const bool a = Validated(forward, p.C, Domain(p.C), p.S, Domain(p.S));
        const bool b = Validated(backward, r, Domain(r), p.S, Domain(p.S));
        if (!(a && b)) {
            return;
        }
        RC_ASSERT(forward->Points.size() == backward->Points.size());
        RC_ASSERT(forward->Overlaps.size() == backward->Overlaps.size());
        const std::size_t n = forward->Points.size();
        for (std::size_t i = 0; i < n; ++i) {
            // The points come in the other order along the reversed curve.
            const CurveSurfacePoint& f = forward->Points[i];
            const CurveSurfacePoint& g = backward->Points[n - 1 - i];
            RC_ASSERT(f.Contact == g.Contact);
            RC_ASSERT(Distance(f.Position, g.Position) <=
                      Rounding(ScaleAt(p.C, f.Parameter, f.Position)) +
                          Rounding(ScaleAt(r, g.Parameter, g.Position)));
        }
    });
}

TEST_CASE("a line's origin moved far along it shifts its parameters only",
          "[intersect][adversarial][property][regression]")
{
    // The audit of 2026-11-05: from 10^8 away a line's two crossings of a
    // circle were placed at one point. Here every surface, origins up to
    // 2^30 back along the line.
    rc::prop("for lines against every surface", [] {
        const auto s = static_cast<SurfaceKind>(Uniform(0, SurfaceKinds));
        RC_TAG(std::string("Line x ") + NameOf(s));
        const Curve near = SomeCurve(CurveKind::Line);
        const Surface surface = SomeSurface(s);
        const Line& line = *near.Target<Line>();
        const double run = std::ldexp(1.0, Uniform(10, 31));
        const Curve far{Line{line.Origin - run * line.Direction, line.Direction}};
        const auto a = Meet(near, surface);
        const auto b = Meet(far, surface);
        const bool va = Validated(a, near, Domain(near), surface, Domain(surface));
        const bool vb = Validated(b, far, Domain(far), surface, Domain(surface));
        if (va && vb) {
            RequireSame(
                *a, *b, near, far, [](const Point3& q) { return q; },
                [&](const double t) { return t + run; });
        }
    });
}

TEST_CASE("a periodic surface's box a turn on gives the same points a turn on",
          "[intersect][adversarial][property]")
{
    rc::prop("for spheres, cylinders, cones, tori and ellipsoids", [] {
        const auto c = static_cast<CurveKind>(Uniform(0, 5));
        const auto s = static_cast<SurfaceKind>(Uniform(1, SurfaceKinds));
        RC_TAG(std::string(NameOf(c)) + " x " + NameOf(s));
        const Curve curve = SomeCurve(c);
        const Surface surface = SomeSurface(s);
        const Box2 chart = Domain(surface);
        const bool torus = s == SurfaceKind::Torus;
        const Box2 base{Interval{0.0, 2.0 * pi}, torus ? Interval{0.0, 2.0 * pi} : chart.Y};
        const Box2 turned{Interval{2.0 * pi, 4.0 * pi}, torus ? Interval{-2.0 * pi, 0.0} : chart.Y};
        const auto a = Intersect(curve, Domain(curve), surface, base, ctx);
        const auto b = Intersect(curve, Domain(curve), surface, turned, ctx);
        const bool va = Validated(a, curve, Domain(curve), surface, base);
        const bool vb = Validated(b, curve, Domain(curve), surface, turned);
        if (!(va && vb)) {
            return;
        }
        RequireSame(
            *a, *b, curve, curve, [](const Point3& q) { return q; },
            [](const double t) { return t; });
    });
}

TEST_CASE("a domain of one parameter away from every point is a proved separation",
          "[intersect][adversarial][property]")
{
    rc::prop("over the grid of curves and surfaces", [] {
        const Pair p = SomePair();
        const Interval whole = Domain(p.C);
        const double s = Uniform(0, 1025) / 1024.0;
        const double t = whole.IsBounded() ? whole.GetLower() + s * Width(whole) : 16.0 * s - 8.0;
        const auto all = Meet(p.C, p.S);
        if (!Validated(all, p.C, whole, p.S, Domain(p.S))) {
            return;
        }
        const auto one = Intersect(p.C, Interval{t}, p.S, Domain(p.S), ctx);
        if (!Validated(one, p.C, Interval{t}, p.S, Domain(p.S))) {
            return;
        }
        const double margin = 1e-6 * std::max(1.0, std::abs(t));
        const bool nearPoint = std::ranges::any_of(all->Points, [&](const CurveSurfacePoint& q) {
            return std::abs(q.Parameter - t) <= margin;
        });
        const bool inOverlap =
            std::ranges::any_of(all->Overlaps, [&](const CurveSurfaceOverlap& o) {
                return o.Range.GetLower() - margin <= t && t <= o.Range.GetUpper() + margin;
            });
        if (!nearPoint && !inOverlap) {
            RC_ASSERT(one->IsEmpty());
        }
        const bool insideOverlap =
            std::ranges::any_of(all->Overlaps, [&](const CurveSurfaceOverlap& o) {
                return o.Range.GetLower() + margin < t && t < o.Range.GetUpper() - margin;
            });
        if (insideOverlap) {
            RC_ASSERT(one->Overlaps.size() == 1);
        }
    });
}

TEST_CASE("a curve in a plane is clipped to exactly the stretches in the box",
          "[intersect][adversarial][property]")
{
    // Overlaps are rare among random pairs. Here every curve lies in the
    // plane z = h, against a random box of the plane's parameters: a
    // parameter whose point is clearly inside the box lies in an overlap,
    // and one clearly outside lies in none.
    rc::prop("for circles, ellipses and B-splines in the plane", [] {
        const auto kind = static_cast<CurveKind>(Uniform(1, 5));
        RC_TAG(NameOf(kind));
        const double h = Coordinate();
        const Surface plane{Plane{Point3{0.0, 0.0, h}, ex, ey}};
        const Point3 c{Coordinate(), Coordinate(), h};
        Curve curve;
        if (kind == CurveKind::Circle) {
            curve = Curve{Circle{c, Length8(4, 20), ex, ey}};
        } else if (kind == CurveKind::Ellipse) {
            const double minor = Length8(4, 16);
            curve = Curve{Ellipse{c, minor + Length8(1, 8), minor, ey, Vec3{-1.0, 0.0, 0.0}}};
        } else {
            std::vector<Point3> poles;
            std::vector<double> weights;
            for (int i = 0; i < 5; ++i) {
                poles.push_back({Coordinate(), Coordinate(), h});
                if (kind == CurveKind::Rational) {
                    weights.push_back(Length8(4, 16));
                }
            }
            auto made = Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0},
                                                  poles, weights);
            RC_PRE(made.has_value());
            curve = Curve{*made};
        }
        const double x0 = Coordinate();
        const double y0 = Coordinate();
        const Box2 box{Interval{x0, x0 + Length8(1, 32)}, Interval{y0, y0 + Length8(1, 32)}};
        const Interval domain = Domain(curve);
        const auto result = Intersect(curve, domain, plane, box, ctx);
        if (!Validated(result, curve, domain, plane, box)) {
            RC_CLASSIFY(true, "failed closed");
            return;
        }
        RC_ASSERT(result->Points.empty());
        RC_CLASSIFY(!result->Overlaps.empty(), "clipped");
        for (const double t : IntersectTest::SamplesOf(domain, 400)) {
            const Point3 q = Eval(curve, t);
            const double margin = 1e-9 * std::max(1.0, Length(Stereon::ToVec(q)));
            const bool inside = q.X > box.X.GetLower() + margin &&
                                q.X < box.X.GetUpper() - margin &&
                                q.Y > box.Y.GetLower() + margin && q.Y < box.Y.GetUpper() - margin;
            const bool outside = q.X < box.X.GetLower() - margin ||
                                 q.X > box.X.GetUpper() + margin ||
                                 q.Y < box.Y.GetLower() - margin || q.Y > box.Y.GetUpper() + margin;
            const bool kept =
                std::ranges::any_of(result->Overlaps, [&](const CurveSurfaceOverlap& o) {
                    return Contains(o.Range, t);
                });
            if (inside) {
                RC_ASSERT(kept);
            }
            if (outside) {
                RC_ASSERT(!kept);
            }
        }
    });
}

TEST_CASE("a line 2^-k either side of a tangency gets the count its distance decides",
          "[intersect][adversarial][property]")
{
    // Lines along x at a height r + delta over a sphere's or a cylinder's
    // centre, and over a torus's tube top, with delta = ±2^-k for k from 20
    // to 50 or zero: the sphere and the cylinder decide the count exactly
    // (two crossings, a touch, none); the torus may fail closed within
    // rounding of its touch, but never reports a wrong count.
    rc::prop("for spheres, cylinders and tori", [] {
        const int which = Uniform(0, 3);
        const int sign = Uniform(-1, 2);
        const int k = Uniform(20, 51);
        const double delta = sign * std::ldexp(1.0, -k);
        const Point3 c{Coordinate(), Coordinate(), Coordinate()};
        const double r = Length8(4, 20);
        const auto expected = [&](const std::size_t below, const std::size_t at) {
            return sign < 0 ? below : sign == 0 ? at : std::size_t{0};
        };
        if (which < 2) {
            RC_TAG(which == 0 ? "Sphere" : "Cylinder");
            const Surface s =
                which == 0 ? Surface{Sphere{c, r, ex, ey}} : Surface{Cylinder{c, r, ex, ey}};
            const Curve line{Line{Point3{c.X - 8.0, c.Y + r + delta, c.Z}, ex}};
            const auto result = Meet(line, s);
            RC_ASSERT(Validated(result, line, Domain(line), s, Domain(s)));
            RC_ASSERT(result->Points.size() == expected(2, 1));
            return;
        }
        RC_TAG("Torus");
        const double minor = Length8(4, 8);
        const Surface s{Torus{c, minor + Length8(4, 16), minor, ex, ey}};
        const Curve line{Line{Point3{c.X - 8.0, c.Y, c.Z + minor + delta}, ex}};
        const auto result = Meet(line, s);
        if (!Validated(result, line, Domain(line), s, Domain(s))) {
            RC_CLASSIFY(true, "failed closed", sign, k / 5 * 5);
            return;
        }
        RC_ASSERT(result->Points.size() == expected(4, 2));
    });
}

TEST_CASE(
    "MR-01 and MR-04: far line retains two sphere crossings under translation and origin shift",
    "[intersect][adversarial][metamorphic][regression]")
{
    // PR #205: cancellation at a distant line origin collapsed two quadric
    // crossings into one. This fixed fixture complements the seeded property
    // campaign above and must detect the original bug deterministically.
    const Curve near{Line{Point3{-2.0, 0.0, 0.0}, ex}};
    const Surface sphere{Sphere{Point3{}, 1.0, ex, ey}};
    const double run = 0x1p27;
    const Curve far{Line{Point3{-2.0 - run, 0.0, 0.0}, ex}};

    const auto reference = Meet(near, sphere);
    const auto distant = Meet(far, sphere);
    REQUIRE(reference.has_value());
    REQUIRE(distant.has_value());
    REQUIRE(reference->Points.size() == 2);
    REQUIRE(distant->Points.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        const auto& p = reference->Points[i];
        const auto& q = distant->Points[i];
        const Point3 expected{i == 0 ? -1.0 : 1.0, 0.0, 0.0};
        CHECK(Distance(p.Position, expected) <= ctx.GetResolution());
        CHECK(Distance(q.Position, expected) <= ctx.GetResolution());
        CHECK(Distance(p.Position, q.Position) <= ctx.GetResolution());
        CHECK(std::abs((p.Parameter + run) - q.Parameter) <= ctx.GetResolution());
        CHECK(p.Contact == q.Contact);
        CHECK(p.Decided == q.Decided);
    }

    // MR-01: transform both operands, not just the answer. The original
    // contact count, parameters and the translated geometric points agree.
    const Vec3 offset{8.0, -4.0, 2.0};
    const Stereon::Transform transform = Stereon::MakeTranslation(offset);
    const Curve movedLine = Transformed(far, transform);
    const Surface movedSphere = Transformed(sphere, transform);
    const auto moved = Meet(movedLine, movedSphere);
    REQUIRE(moved.has_value());
    REQUIRE(moved->Points.size() == distant->Points.size());
    for (std::size_t i = 0; i < distant->Points.size(); ++i) {
        CHECK(Distance(moved->Points[i].Position, distant->Points[i].Position + offset) <=
              ctx.GetResolution());
        CHECK(std::abs(moved->Points[i].Parameter - distant->Points[i].Parameter) <=
              ctx.GetResolution());
        CHECK(moved->Points[i].Contact == distant->Points[i].Contact);
    }
}
