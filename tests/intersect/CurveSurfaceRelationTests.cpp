// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0053 relations MR-01, MR-03 and MR-04 for every ordered pair of a
// built-in curve and a built-in surface. MR-02 does not apply: there is no
// exposed surface-curve operation to exchange the operands into (ADR-0053
// section 1).
//
// The fixture is built in closed form. On each surface, a point P at chart
// parameters away from seams and poles, and its unit normal n, from the
// surface's defining formula evaluated here, not by the kernel. Each curve
// kind passes through P at a known parameter t*, its tangent there the
// transverse direction d = cos(b) n + sin(b) e (b = 0.4, e across n):
//
// - a line from P - 2d along d: t* = 2;
// - a circle, an ellipse and a hyperbola with their vertex at P, centre
//   P - a x, x across d, y = d: t* = 0;
// - a parabola with its vertex at P, x across d, y = d: t* = 0;
// - a cubic Bezier, polynomial and rational, symmetric about P along d:
//   u* = 1/2.
//
// The claim: the result has a point within the resolution of P whose curve
// parameter is t*, and every point reported lies on the curve and on the
// surface, the latter by the surface's closed-form distance. The relations
// move the fixture exactly and map P and t* in closed form:
//
// - MR-01: cyclic axis permutations and a quarter turn about z, exact, with
//   and without the translation (8, -4, 2);
// - MR-03: scales 1/4 and 4, the resolution and tolerance ceiling scaled;
//   a parabola's parameter is a length and scales (t* = 0 is fixed);
// - MR-04: the curve in an equivalent chart: a line's origin moved two
//   directions on with its direction doubled (t' = (t - 2) / 2), a
//   circle's frame turned a quarter (t' = t - pi/2), an ellipse's turned
//   half (t' = t - pi), a hyperbola's or a parabola's y axis negated
//   (t' = -t), a B-spline's knots mapped by u' = 3 - 2u (reversed and
//   stretched).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using namespace Stereon;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double tilt = 0.4;

/// p -> Offset + Scale Q p, with Q an exact rotation given by its columns.
struct Motion {
    std::array<Vec3, 3> Columns{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    Vec3 Offset{};
    double Scale = 1.0;

    [[nodiscard]] Vec3 Turn(const Vec3& v) const
    {
        return v.X * Columns[0] + v.Y * Columns[1] + v.Z * Columns[2];
    }
    [[nodiscard]] Point3 Of(const Point3& p) const
    {
        return Point3{} + Offset + Scale * Turn(p - Point3{});
    }
    /// The inverse: the local point of a placed one.
    [[nodiscard]] Point3 Back(const Point3& p) const
    {
        const Vec3 w = (p - Point3{} - Offset) / Scale;
        return Point3{Dot(w, Columns[0]), Dot(w, Columns[1]), Dot(w, Columns[2])};
    }
    [[nodiscard]] Context Policy() const
    {
        return Context::Default()
            .WithResolution(Context::DefaultResolution * Scale)
            .WithToleranceCeiling(Context::DefaultToleranceCeiling * Scale)
            .WithScheduler(Scheduler::Serial);
    }
};

/// A surface in its local frame, placed by a motion: the surface, a point
/// P on it and its unit normal, and the closed-form distance of a local
/// point from it (for the ellipsoid, its algebraic value times the least
/// semi-axis over two, an estimate near the surface).
struct SurfaceFixture {
    std::string Name;
    std::function<Surface(const Motion&)> Make;
    Point3 P;
    Vec3 Normal;
    std::function<double(const Point3&)> Off;
};

std::vector<SurfaceFixture> Surfaces()
{
    using std::cos;
    using std::sin;
    std::vector<SurfaceFixture> out;
    const Vec3 ex{1, 0, 0};
    const Vec3 ey{0, 1, 0};
    const double u = 1.0;
    {
        const Point3 p{0.3, 0.2, 0.0};
        out.push_back({"plane",
                       [](const Motion& m) {
                           return Surface{
                               Plane{m.Of(Point3{}), m.Turn(Vec3{1, 0, 0}), m.Turn(Vec3{0, 1, 0})}};
                       },
                       p, Vec3{0, 0, 1}, [](const Point3& q) { return std::abs(q.Z); }});
    }
    {
        const double r = 1.5;
        const Point3 p{r * cos(u), r * sin(u), 0.7};
        out.push_back({"cylinder",
                       [=](const Motion& m) {
                           return Surface{
                               Cylinder{m.Of(Point3{}), m.Scale * r, m.Turn(ex), m.Turn(ey), true}};
                       },
                       p, Vec3{cos(u), sin(u), 0},
                       [=](const Point3& q) { return std::abs(std::hypot(q.X, q.Y) - r); }});
    }
    {
        // Radius 1 at the origin, semi-angle 0.5: radial distance 1 + v sin a,
        // height v cos a; at v = 1.5. Both nappes: past the apex the radial
        // distance is |1 + v sin a|, which the distance below follows.
        const double r = 1.0;
        const double a = 0.5;
        const double v = 1.5;
        const double rho = r + v * sin(a);
        const Point3 p{rho * cos(u), rho * sin(u), v * cos(a)};
        out.push_back(
            {"cone",
             [=](const Motion& m) {
                 return Surface{Cone{m.Of(Point3{}), m.Scale * r, a, m.Turn(ex), m.Turn(ey), true}};
             },
             p, Vec3{cos(a) * cos(u), cos(a) * sin(u), -sin(a)},
             [=](const Point3& q) {
                 return std::abs(std::hypot(q.X, q.Y) - std::abs(r + q.Z * std::tan(a))) * cos(a);
             }});
    }
    {
        const double r = 2.0;
        const double v = 0.6;
        const Vec3 n{cos(v) * cos(u), cos(v) * sin(u), sin(v)};
        out.push_back({"sphere",
                       [=](const Motion& m) {
                           return Surface{
                               Sphere{m.Of(Point3{}), m.Scale * r, m.Turn(ex), m.Turn(ey), true}};
                       },
                       Point3{} + r * n, n,
                       [=](const Point3& q) { return std::abs(Length(q - Point3{}) - r); }});
    }
    {
        const double big = 3.0;
        const double small = 1.0;
        const double v = 0.7;
        const Vec3 radial{cos(u), sin(u), 0};
        const Vec3 n = cos(v) * radial + sin(v) * Vec3{0, 0, 1};
        out.push_back({"torus",
                       [=](const Motion& m) {
                           return Surface{Torus{m.Of(Point3{}), m.Scale * big, m.Scale * small,
                                                m.Turn(ex), m.Turn(ey), true}};
                       },
                       Point3{} + big * radial + small * n, n,
                       [=](const Point3& q) {
                           return std::abs(std::hypot(std::hypot(q.X, q.Y) - big, q.Z) - small);
                       }});
    }
    {
        const double a = 3.0;
        const double b = 2.0;
        const double c = 1.5;
        const double v = 0.5;
        const Point3 p{a * cos(v) * cos(u), b * cos(v) * sin(u), c * sin(v)};
        const Vec3 g{p.X / (a * a), p.Y / (b * b), p.Z / (c * c)};
        out.push_back({"ellipsoid",
                       [=](const Motion& m) {
                           return Surface{Ellipsoid{m.Of(Point3{}), m.Scale * a, m.Scale * b,
                                                    m.Scale * c, m.Turn(ex), m.Turn(ey), true}};
                       },
                       p, g / Length(g),
                       [=](const Point3& q) {
                           const double f = q.X * q.X / (a * a) + q.Y * q.Y / (b * b) +
                                            q.Z * q.Z / (c * c) - 1.0;
                           return std::abs(f) * c / 2.0;
                       }});
    }
    {
        // A bilinear B-spline patch that is exactly the plane z = 0 over
        // [-4, 4]^2.
        out.push_back({"bspline-plane",
                       [](const Motion& m) {
                           const std::vector<double> knots{0, 0, 1, 1};
                           std::vector<Point3> poles;
                           for (const double x : {-4.0, 4.0}) {
                               for (const double y : {-4.0, 4.0}) {
                                   poles.push_back(m.Of(Point3{x, y, 0}));
                               }
                           }
                           auto made = MakeBSplineSurface(1, 1, knots, knots, std::move(poles));
                           REQUIRE(made.has_value());
                           return Surface{*std::move(made)};
                       },
                       Point3{0.3, 0.2, 0.0}, Vec3{0, 0, 1},
                       [](const Point3& q) { return std::abs(q.Z); }});
    }
    return out;
}

/// A unit vector across n.
Vec3 Across(const Vec3& n)
{
    const Vec3 seed = std::abs(n.X) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const Vec3 e = Cross(n, seed);
    return e / Length(e);
}

/// A curve through the placed point along d, its domain and the parameter
/// at the point, each in the chart asked for.
struct CurveFixture {
    Curve Value;
    Interval Domain;
    double At = 0.0;
};

enum class Kind : std::uint8_t { Line, Circle, Ellipse, Hyperbola, Parabola, Spline, Rational };

constexpr std::array<Kind, 7> Kinds{Kind::Line,     Kind::Circle, Kind::Ellipse, Kind::Hyperbola,
                                    Kind::Parabola, Kind::Spline, Kind::Rational};

const char* NameOf(const Kind k)
{
    constexpr std::array<const char*, 7> names{"line",     "circle", "ellipse", "hyperbola",
                                               "parabola", "spline", "rational"};
    return names[static_cast<std::size_t>(k)];
}

/// The angle at which a circle or an ellipse passes through P: off the
/// quarter turns where the kernel cuts its arcs, so that the crossing is
/// interior to a piece.
constexpr double ConicAngle = 0.6;
/// The same for the B-splines' parameter: off the dyadic subdivision.
constexpr double SplineAt = 0.375;

/// The cubic B-spline through P along d at u = SplineAt: exactly the
/// quadratic c(u) = P + L (u - s) d + h (u - s)^2 x in cubic Bernstein
/// form. The rational one is the same curve with the linear weight
/// W(u) = 1 + u: numerator W c, weights W, each elevated to degree three, so
/// its points are c(u) to rounding and the claim does not change. Recharted,
/// u' = 3 - 2u: knots on [1, 3], poles and weights reversed.
BSplineCurve Cubic(const Point3& p, const Vec3& d, const Vec3& x, const double scale,
                   const bool rational, const bool recharted)
{
    const double length = 4.0 * scale;
    const double bend = 0.5 * scale;
    // c(u) in powers of u, as (along d, along x) coefficients.
    const double s0 = SplineAt;
    const std::array<double, 3> along{-length * s0, length, 0.0};
    const std::array<double, 3> across{bend * s0 * s0, -2.0 * bend * s0, bend};
    // Times W(u) = w0 + w1 u for the numerator; W = 1 for the polynomial.
    const double w0 = 1.0;
    const double w1 = rational ? 1.0 : 0.0;
    const auto times = [&](const std::array<double, 3>& q) {
        return std::array<double, 4>{w0 * q[0], w0 * q[1] + w1 * q[0], w0 * q[2] + w1 * q[1],
                                     w1 * q[2]};
    };
    // Power basis to cubic Bernstein: b_i = sum_j C(i, j) / C(3, j) a_j.
    const auto bernstein = [](const std::array<double, 4>& a) {
        return std::array<double, 4>{a[0], a[0] + a[1] / 3.0, a[0] + 2.0 * a[1] / 3.0 + a[2] / 3.0,
                                     a[0] + a[1] + a[2] + a[3]};
    };
    const std::array<double, 4> nd = bernstein(times(along));
    const std::array<double, 4> nx = bernstein(times(across));
    const std::array<double, 4> w = bernstein({w0, w1, 0.0, 0.0});
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < 4; ++i) {
        poles.push_back(p + (nd[i] / w[i]) * d + (nx[i] / w[i]) * x);
        weights.push_back(w[i]);
    }
    if (!rational) {
        weights.clear();
    }
    std::vector<double> knots{0, 0, 0, 0, 1, 1, 1, 1};
    if (recharted) {
        knots = {1, 1, 1, 1, 3, 3, 3, 3};
        std::ranges::reverse(poles);
        std::ranges::reverse(weights);
    }
    return BSplineCurve{3, std::move(knots), std::move(poles), std::move(weights)};
}

/// The frame of a circle or ellipse that passes through P at ConicAngle
/// with its tangent near d: X and Y turned so that cos(t) X + sin(t) Y is
/// x, across d, at t = ConicAngle.
std::pair<Vec3, Vec3> ConicFrame(const Vec3& d, const Vec3& x)
{
    const double c = std::cos(ConicAngle);
    const double s = std::sin(ConicAngle);
    return {c * x - s * d, s * x + c * d};
}

CurveFixture CurveThrough(const Kind kind, const Point3& p, const Vec3& d, const Vec3& x,
                          const double scale, const bool recharted)
{
    const double r = 3.0 * scale;
    const double minor = 2.0 * scale;
    const auto [fx, fy] = ConicFrame(d, x);
    const double c = std::cos(ConicAngle);
    const double s = std::sin(ConicAngle);
    switch (kind) {
    case Kind::Line:
        if (recharted) {
            return {Curve{Line{p, 2.0 * scale * d}}, Interval{-6.0, 4.0}, 0.0};
        }
        return {Curve{Line{p - 2.0 * scale * d, scale * d}}, Interval{-10.0, 10.0}, 2.0};
    case Kind::Circle: {
        const Point3 centre = p - r * x;
        if (recharted) {
            return {Curve{Circle{centre, r, fy, -fx}}, Interval{-pi, pi}, ConicAngle - pi / 2};
        }
        return {Curve{Circle{centre, r, fx, fy}}, Interval{-pi, pi}, ConicAngle};
    }
    case Kind::Ellipse: {
        const Point3 centre = p - r * c * fx - minor * s * fy;
        if (recharted) {
            return {Curve{Ellipse{centre, r, minor, -fx, -fy}}, Interval{-2 * pi, 0.0},
                    ConicAngle - pi};
        }
        return {Curve{Ellipse{centre, r, minor, fx, fy}}, Interval{-pi, pi}, ConicAngle};
    }
    case Kind::Hyperbola:
        return {Curve{Hyperbola{p - r * x, r, minor, x, recharted ? -d : d}}, Interval{-2.0, 2.0},
                0.0};
    case Kind::Parabola:
        return {Curve{Parabola{p, scale, x, recharted ? -d : d}},
                Interval{-4.0 * scale, 4.0 * scale}, 0.0};
    case Kind::Spline:
    case Kind::Rational:
        return {Curve{Cubic(p, d, x, scale, kind == Kind::Rational, recharted)},
                recharted ? Interval{1.0, 3.0} : Interval{0.0, 1.0},
                recharted ? 3.0 - 2.0 * SplineAt : SplineAt};
    }
    return {};
}

/// The parameters differ by at most this near t*: a few hundred roundings
/// of the curve's parameter scale, far below the next feature.
constexpr double ParameterReach = 1e-7;

struct Tally {
    std::size_t Answered = 0;
    std::size_t Refused = 0;
    std::size_t Unsupported = 0;
};

/// The claim for one pair under one motion and chart.
void Claim(const SurfaceFixture& s, const Kind kind, const Motion& m, const bool recharted,
           Tally& tally)
{
    INFO(NameOf(kind) << " x " << s.Name << " scale=" << m.Scale << " offset=" << m.Offset.X
                      << " turn=" << m.Columns[0].X << m.Columns[0].Y << m.Columns[0].Z
                      << " recharted=" << recharted);
    const Surface surface = s.Make(m);
    const Vec3 n = m.Turn(s.Normal);
    const Vec3 e = Across(n);
    const Vec3 d = std::cos(tilt) * n + std::sin(tilt) * e;
    const Vec3 x = Across(d);
    const Point3 p = m.Of(s.P);
    const CurveFixture c = CurveThrough(kind, p, d, x, m.Scale, recharted);
    const Context policy = m.Policy();
    const double reach = 4.0 * policy.GetResolution();
    const auto r = Intersect(c.Value, c.Domain, surface, Domain(surface), policy);
    if (!r) {
        INFO(IntersectTest::ErrorText(r));
        const ErrorCode code = r.error().GetCode();
        if (code == ErrorCode::NotImplemented) {
            // No handler for the pair yet (a hyperbola or a parabola against
            // a surface): its cells stay open, nothing is claimed.
            ++tally.Unsupported;
            return;
        }
        CHECK((code == ErrorCode::IntersectionFailed || code == ErrorCode::CertificationFailed));
        ++tally.Refused;
        return;
    }
    ++tally.Answered;
    CHECK(r->Overlaps.empty());
    bool found = false;
    for (const CurveSurfacePoint& q : r->Points) {
        CHECK(Distance(Eval(c.Value, q.Parameter), q.Position) <= reach);
        CHECK(s.Off(m.Back(q.Position)) * m.Scale <= reach);
        found = found || (Distance(q.Position, p) <= reach &&
                          std::abs(q.Parameter - c.At) <= ParameterReach * (1.0 + std::abs(c.At)));
    }
    CHECK(found);
}

const std::array<std::array<Vec3, 3>, 4> Turns{{
    {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}},
    {Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}},
    {Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}},
    {Vec3{0, 1, 0}, Vec3{-1, 0, 0}, Vec3{0, 0, 1}}, // a quarter turn about z
}};

void RunAll(const std::function<void(const SurfaceFixture&, Kind, Tally&)>& run)
{
    Tally tally;
    for (const SurfaceFixture& s : Surfaces()) {
        for (const Kind kind : Kinds) {
            run(s, kind, tally);
        }
    }
    WARN(tally.Answered << " answered, " << tally.Refused << " refused (fail-closed), "
                        << tally.Unsupported << " without a handler");
}

} // namespace

TEST_CASE("MR-01: curve-surface crossings move with exact rigid motions",
          "[intersect][curve-surface][metamorphic]")
{
    RunAll([](const SurfaceFixture& s, const Kind kind, Tally& tally) {
        for (const auto& turn : Turns) {
            for (const Vec3 offset : {Vec3{}, Vec3{8, -4, 2}}) {
                Motion m;
                m.Columns = turn;
                m.Offset = offset;
                Claim(s, kind, m, false, tally);
            }
        }
    });
}

TEST_CASE("MR-03: curve-surface crossings scale with the geometry and its resolution",
          "[intersect][curve-surface][metamorphic]")
{
    RunAll([](const SurfaceFixture& s, const Kind kind, Tally& tally) {
        for (const double scale : {0.25, 1.0, 4.0}) {
            Motion m;
            m.Scale = scale;
            Claim(s, kind, m, false, tally);
        }
    });
}

TEST_CASE("MR-04: curve-surface crossings keep their place in an equivalent curve chart",
          "[intersect][curve-surface][metamorphic]")
{
    RunAll([](const SurfaceFixture& s, const Kind kind, Tally& tally) {
        for (const auto& turn : Turns) {
            Motion m;
            m.Columns = turn;
            Claim(s, kind, m, false, tally);
            Claim(s, kind, m, true, tally);
        }
    });
}
