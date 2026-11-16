// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Intersection derivatives (ADR-0021, sections 6 and 9; the v0.2 gate): a
// line crossing a plane, a sphere or a cylinder, differentiated through
// F = C(t) - S(u, v) = 0. The forward product agrees with the closed form of
// each crossing; on random geometry its Taylor residual |y(p + h v) - y(p)
// - h J v| decays with the second power of h; the reverse product is its
// transpose; and every point that has no derivative says so: an overlap or
// a touch is Nondifferentiable, an undecided point or a pole of a sphere
// Indeterminate, a grazing crossing NearEvent, and none of them carries a
// usable value unless qualified.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Cone;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::DerivativeStatus;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::Line;
using Stereon::LineBinding;
using Stereon::LineComponent;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceBinding;
using Stereon::SurfaceComponent;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr std::size_t Values = Stereon::CurveSurfacePointValues;
constexpr double tau = 2.0 * std::numbers::pi;

/// The surface's origin, a sphere's centre, and whether it has a radius.
Point3& OriginOf(Plane& s)
{
    return s.Origin;
}
Point3& OriginOf(Sphere& s)
{
    return s.Center;
}
Point3& OriginOf(Cylinder& s)
{
    return s.Origin;
}
Point3& OriginOf(Cone& s)
{
    return s.Origin;
}
Point3& OriginOf(Torus& s)
{
    return s.Center;
}
Point3& OriginOf(Ellipsoid& s)
{
    return s.Center;
}

/// The frame's third axis: a plane's X x Y, any other surface's Axis.
template <typename S>
Vec3 ThirdAxis(const S& s)
{
    if constexpr (std::is_same_v<S, Plane>) {
        return Cross(s.XAxis, s.YAxis);
    } else {
        return Axis(s);
    }
}

/// v turned by the rotation vector w (Rodrigues).
Vec3 Rotated(const Vec3& v, const Vec3& w)
{
    const double angle = Length(w);
    if (angle == 0.0) {
        return v;
    }
    const Vec3 k = w / angle;
    return std::cos(angle) * v + std::sin(angle) * Cross(k, v) +
           ((1.0 - std::cos(angle)) * Dot(k, v)) * k;
}

/// A frame (x, y) turned by w0 x + w1 y + w2 z: an exact rotation, so the
/// frame stays orthonormal and the first-order motion is w x (P - origin).
void Turn(Vec3& x, Vec3& y, const Vec3& z, const double w0, const double w1, const double w2)
{
    const Vec3 w = w0 * x + w1 * y + w2 * z;
    x = Rotated(x, w);
    y = Rotated(y, w);
}

bool IsTurn(const SurfaceComponent c)
{
    return c == SurfaceComponent::RotationX || c == SurfaceComponent::RotationY ||
           c == SurfaceComponent::RotationZ;
}

/// What each surface's bindings move, in the order of their parameters:
/// the origin, then its radius, then a torus's minor radius or a cone's
/// semi-angle; an ellipsoid's three semi-axes in place of a radius; and,
/// when `turning`, the frame's turn about its three axes.
template <typename S>
std::vector<SurfaceComponent> ComponentsOf(const bool turning = false)
{
    std::vector<SurfaceComponent> components{SurfaceComponent::OriginX, SurfaceComponent::OriginY,
                                             SurfaceComponent::OriginZ};
    if constexpr (std::is_same_v<S, Ellipsoid>) {
        components.insert(components.end(),
                          {SurfaceComponent::SemiAxisA, SurfaceComponent::SemiAxisB,
                           SurfaceComponent::SemiAxisC});
    } else if constexpr (!std::is_same_v<S, Plane>) {
        components.push_back(SurfaceComponent::Radius);
    }
    if constexpr (std::is_same_v<S, Torus>) {
        components.push_back(SurfaceComponent::MinorRadius);
    }
    if constexpr (std::is_same_v<S, Cone>) {
        components.push_back(SurfaceComponent::SemiAngle);
    }
    if (turning) {
        components.insert(components.end(),
                          {SurfaceComponent::RotationX, SurfaceComponent::RotationY,
                           SurfaceComponent::RotationZ});
    }
    return components;
}

template <typename S>
std::size_t SurfaceParameters(const bool turning = false)
{
    return ComponentsOf<S>(turning).size();
}

/// The quantity of the parameter moving a surface component.
ParameterQuantity QuantityOf(const SurfaceComponent component)
{
    return component == SurfaceComponent::SemiAngle || IsTurn(component)
               ? ParameterQuantity::Angle
               : ParameterQuantity::Length;
}

/// The surface moved by d(first + k) along its k-th component; when
/// `turning`, its frame turned last.
template <typename S, typename Step>
void MoveSurface(S& surface, const Step& d, const std::size_t first, const bool turning = false)
{
    if (turning) {
        const std::size_t k = first + ComponentsOf<S>().size();
        Turn(surface.XAxis, surface.YAxis, ThirdAxis(surface), d(k), d(k + 1), d(k + 2));
    }
    OriginOf(surface) = OriginOf(surface) + Vec3{d(first), d(first + 1), d(first + 2)};
    if constexpr (std::is_same_v<S, Sphere> || std::is_same_v<S, Cylinder>) {
        surface.Radius += d(first + 3);
    } else if constexpr (std::is_same_v<S, Cone>) {
        surface.Radius += d(first + 3);
        surface.SemiAngle += d(first + 4);
    } else if constexpr (std::is_same_v<S, Torus>) {
        surface.MajorRadius += d(first + 3);
        surface.MinorRadius += d(first + 4);
    } else if constexpr (std::is_same_v<S, Ellipsoid>) {
        surface.A += d(first + 3);
        surface.B += d(first + 4);
        surface.C += d(first + 5);
    }
}

/// Parameters each bound to one coordinate: the line's origin and
/// direction, then the surface's origin, then its radius where it has one.
struct Model {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<LineBinding> LineBindings;
    std::vector<SurfaceBinding> SurfaceBindings;
};

/// The line's six bindings, then one per surface component, in order.
Model ModelOf(const std::vector<double>& references,
              const std::vector<SurfaceComponent>& components)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const ParameterQuantity quantity =
            k < 6 ? ParameterQuantity::Length : QuantityOf(components[k - 6]);
        ids.push_back(*set.Add("p" + std::to_string(k), quantity, references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < 6; ++k) {
        model.LineBindings.push_back({ids[k], static_cast<LineComponent>(k)});
    }
    for (std::size_t k = 6; k < references.size(); ++k) {
        model.SurfaceBindings.push_back({ids[k], components[k - 6]});
    }
    return model;
}

/// A plane's bindings for nine references, a sphere's or a cylinder's for
/// ten.
Model FullModel(const std::vector<double>& references)
{
    return ModelOf(references,
                   references.size() == 9 ? ComponentsOf<Plane>() : ComponentsOf<Sphere>());
}

template <typename S>
Model ModelFor(const std::vector<double>& references, const bool turning = false)
{
    return ModelOf(references, ComponentsOf<S>(turning));
}

/// The geometry moved by `step` times the seed, each parameter by its
/// reference per unit.
template <typename S>
std::pair<Line, S> Moved(Line line, S surface, const std::vector<double>& seed,
                         const std::vector<double>& references, const double step,
                         const bool turning = false)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    line.Origin = line.Origin + Vec3{d(0), d(1), d(2)};
    line.Direction = line.Direction + Vec3{d(3), d(4), d(5)};
    MoveSurface(surface, d, 6, turning);
    return {line, surface};
}

template <typename C, typename S>
CurveSurfaceIntersection Meet(const C& curve, const S& surface)
{
    auto result = Intersect(Stereon::Curve{curve}, Surface{surface}, ctx);
    if (!result) {
        RC_FAIL(IntersectTest::ErrorText(result));
    }
    return *std::move(result);
}

/// The points as the product's values are laid out: x, then t, u, v.
template <typename C, typename S>
std::vector<double> Points(const C& curve, const S& surface)
{
    std::vector<double> values;
    for (const CurveSurfacePoint& p : Meet(curve, surface).Points) {
        values.insert(values.end(), {p.Position.X, p.Position.Y, p.Position.Z, p.Parameter,
                                     p.OnSurface.X, p.OnSurface.Y});
    }
    return values;
}

template <typename S>
std::vector<double> Jvp(Model& model, const Line& line, const S& surface,
                        const std::vector<double>& seed)
{
    RC_ASSERT(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp =
        Stereon::IntersectJvp(line, model.LineBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(line, surface), model.Sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

double Norm(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

/// Random references and a random vector over `count` parameters.
std::pair<std::vector<double>, std::vector<double>> RandomSeed(const std::size_t count)
{
    std::vector<double> references;
    std::vector<double> seed;
    for (std::size_t k = 0; k < count; ++k) {
        references.push_back(*rc::gen::inRange(1, 20) / 10.0);
        seed.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    return {references, seed};
}

/// Every point crosses well away from a tangency and from the surface's
/// poles, so that the steps below stay in one regime.
template <typename S>
bool WellConditioned(const Line& line, const S& surface, const std::size_t points)
{
    const CurveSurfaceIntersection at = Meet(line, surface);
    if (at.Points.size() != points || !at.Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(at.Points, [&](const CurveSurfacePoint& p) {
        const auto lin = Stereon::Linearize(line, Surface{surface}, p);
        return lin && lin->Conditioning >= 0.25;
    });
}

/// The Taylor residual at steps 1e-2 and 1e-3 falls by second order. The
/// angle u is compared modulo its period: a step may carry it over a seam.
template <typename S>
void RequireTaylor(const Line& line, const S& surface, const bool turning = false)
{
    const auto [references, seed] = RandomSeed(6 + SurfaceParameters<S>(turning));
    Model model = ModelFor<S>(references, turning);
    const std::vector<double> jv = Jvp(model, line, surface, seed);
    const std::vector<double> base = Points(line, surface);
    const auto residual = [&](const double h) {
        const auto [l, s] = Moved(line, surface, seed, references, h, turning);
        const std::vector<double> moved = Points(l, s);
        RC_ASSERT(moved.size() == base.size());
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < base.size(); ++i) {
            const double difference = moved[i] - base[i];
            // u is an angle, and so is a torus's v; a step never moves
            // any other v by half a turn, so v is taken modulo one too.
            const bool angle = i % Values == 4 || i % Values == 5;
            r[i] = (angle ? std::remainder(difference, tau) : difference) - h * jv[i];
        }
        return Norm(r);
    };
    const double coarse = residual(1e-2);
    const double fine = residual(1e-3);
    // Second order: a tenth of the step leaves about a hundredth of the
    // residual; rounding sets a floor well below both.
    RC_ASSERT(fine <= coarse / 30.0 + 1e-10 * (1.0 + Norm(base)));
}

/// <J v, w> = <v, J^T w> over every point.
template <typename S>
void RequireTranspose(const Line& line, const S& surface, const bool turning = false)
{
    const auto [references, v] = RandomSeed(6 + SurfaceParameters<S>(turning));
    Model model = ModelFor<S>(references, turning);
    const std::vector<double> jv = Jvp(model, line, surface, v);
    std::vector<double> w;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
    const auto vjp =
        Stereon::IntersectVjp(line, model.LineBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(line, surface), model.Sensitivity);
    RC_ASSERT(vjp.has_value());
    RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
    const auto jtw = *vjp->GetValues();
    double left = 0.0;
    double right = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        left += jv[i] * w[i];
        scale += std::abs(jv[i] * w[i]);
    }
    for (std::size_t k = 0; k < v.size(); ++k) {
        right += v[k] * jtw[k];
        scale += std::abs(v[k] * jtw[k]);
    }
    RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
}

double Coordinate()
{
    return *rc::gen::inRange(-40, 41) / 8.0;
}

/// A random line through a random plane, well away from parallel.
std::pair<Line, Plane> RandomPlaneCrossing()
{
    const Vec3 normal{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
    RC_PRE(normal != Vec3{} && direction != Vec3{});
    // Not parallel: along the plane, the frame's rounding can leave an
    // exact slope that doubles cannot place, which Intersect refuses.
    RC_PRE(std::abs(Dot(normal, direction)) >= 0.25 * Length(normal) * Length(direction));
    const Plane plane =
        IntersectTest::PlaneWithNormal(Point3{Coordinate(), Coordinate(), Coordinate()}, normal);
    const Line line{Point3{Coordinate(), Coordinate(), Coordinate()}, direction};
    RC_PRE(WellConditioned(line, plane, 1));
    return {line, plane};
}

/// A random line through the inside of a random sphere and of a cylinder
/// about an axis across the line through the sphere's centre: within half
/// the radius of the centre on each coordinate, so every crossing is at a
/// steep angle. The sphere's crossings stay off its poles, where the angle
/// u turns too fast for a fixed step to stay in its regime.
std::pair<Line, std::pair<Sphere, Cylinder>> RandomQuadricCrossing()
{
    const Point3 center{Coordinate(), Coordinate(), Coordinate()};
    const double radius = *rc::gen::inRange(4, 33) / 8.0;
    const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 axis = Cross(direction, Vec3{Coordinate(), Coordinate(), Coordinate()});
    RC_PRE(direction != Vec3{} && Length(axis) >= 0.25 * Length(direction));
    const auto offset = [&] { return *rc::gen::inRange(-4, 5) / 8.0 * radius; };
    const Line line{center + Vec3{offset(), offset(), offset()} - 3.0 * direction, direction};
    const Sphere sphere{center, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Cylinder cylinder = IntersectTest::CylinderAbout(center, axis, radius);
    RC_PRE(WellConditioned(line, sphere, 2) && WellConditioned(line, cylinder, 2));
    RC_PRE(std::ranges::all_of(Meet(line, sphere).Points, [](const CurveSurfacePoint& p) {
        return std::abs(p.OnSurface.Y) <= 1.2;
    }));
    return {line, {sphere, cylinder}};
}

/// The crossings of a line at the parameters of their closed form, for the
/// cone and the torus, which Intersect has no line handler for yet.
template <typename S>
CurveSurfaceIntersection Crossings(const Line& line, const S& surface,
                                   const std::vector<double>& parameters)
{
    CurveSurfaceIntersection at;
    for (const double t : parameters) {
        const Point3 p = Eval(line, t);
        at.Points.push_back({t, Project(surface, p), p, Stereon::ContactType::Crossing, true});
    }
    return at;
}

/// The forward product at those crossings for the unit seed of parameter k,
/// every reference 1.
template <typename S>
std::vector<double> RateAt(const Line& line, const S& surface, const CurveSurfaceIntersection& at,
                           const std::size_t k)
{
    const std::size_t count = 6 + SurfaceParameters<S>();
    Model model = ModelFor<S>(std::vector<double>(count, 1.0));
    std::vector<double> seed(count, 0.0);
    seed[k] = 1.0;
    REQUIRE(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp = Stereon::IntersectJvp(line, model.LineBindings, Surface{surface},
                                           model.SurfaceBindings, at, model.Sensitivity);
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    return {jvp->GetValues()->begin(), jvp->GetValues()->end()};
}

} // namespace

TEST_CASE("the crossing of a plane moves as its closed form says", "[intersect][linearization]")
{
    // z = h meets O + t d at t = (h - Oz) / dz.
    const Line line{Point3{1.0, 2.0, -1.0}, Vec3{0.5, 0.25, 3.0}};
    const Plane plane{Point3{0.0, 0.0, 0.5}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const double t = (0.5 - -1.0) / 3.0;
    Model model = FullModel(std::vector<double>(9, 1.0));
    const auto rate = [&](const std::size_t k) {
        std::vector<double> seed(9, 0.0);
        seed[k] = 1.0;
        return Jvp(model, line, plane, seed);
    };
    // The plane's height: t' = 1 / dz, and the point moves along d.
    const std::vector<double> height = rate(8);
    REQUIRE(std::abs(height[3] - 1.0 / 3.0) <= 1e-15);
    REQUIRE(std::abs(height[0] - 0.5 / 3.0) <= 1e-15);
    REQUIRE(std::abs(height[2] - 1.0) <= 1e-15);
    REQUIRE(std::abs(height[4] - height[0]) <= 1e-15);
    // The line's height: t' = -1 / dz.
    REQUIRE(std::abs(rate(2)[3] + 1.0 / 3.0) <= 1e-15);
    // The line's vertical speed: t' = -t / dz, and the point stays put in z.
    const std::vector<double> speed = rate(5);
    REQUIRE(std::abs(speed[3] + t / 3.0) <= 1e-15);
    REQUIRE(std::abs(speed[2]) <= 1e-15);
    // Sliding the line along x moves the point along x, t unchanged.
    const std::vector<double> slide = rate(0);
    REQUIRE(std::abs(slide[0] - 1.0) <= 1e-15);
    REQUIRE(std::abs(slide[3]) <= 1e-15);
    // Sliding the plane in its own plane moves no point, only its (u, v).
    const std::vector<double> inPlane = rate(6);
    REQUIRE(std::abs(inPlane[0]) <= 1e-15);
    REQUIRE(std::abs(inPlane[4] + 1.0) <= 1e-15);
}

TEST_CASE("the crossings of a sphere and a cylinder move as their closed forms say",
          "[intersect][linearization]")
{
    // A line along x at height y0 crosses a circle of radius r about the z
    // axis at x = -+ sqrt(r² - y0²); growing r moves them by x' = r / x.
    constexpr double r = 1.5;
    constexpr double y0 = 0.5;
    const double x = std::sqrt(r * r - y0 * y0);
    const Line line{Point3{-5.0, y0, 0.25}, Vec3{1.0, 0.0, 0.0}};
    const auto rate = [&](const auto& surface, const std::size_t k) {
        Model model = FullModel(std::vector<double>(10, 1.0));
        std::vector<double> seed(10, 0.0);
        seed[k] = 1.0;
        return Jvp(model, line, surface, seed);
    };
    const Sphere sphere{Point3{0.0, 0.0, 0.25}, r, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Cylinder cylinder{Point3{}, r, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    for (const std::vector<double>& radius : {rate(sphere, 9), rate(cylinder, 9)}) {
        REQUIRE(radius.size() == 2 * Values);
        REQUIRE(std::abs(radius[0] + r / x) <= 1e-14);
        REQUIRE(std::abs(radius[3] + r / x) <= 1e-14);
        REQUIRE(std::abs(radius[Values] - r / x) <= 1e-14);
        REQUIRE(std::abs(radius[Values + 3] - r / x) <= 1e-14);
    }
    // Sliding either along x carries both crossings along x.
    for (const std::vector<double>& slide : {rate(sphere, 6), rate(cylinder, 6)}) {
        REQUIRE(std::abs(slide[0] - 1.0) <= 1e-14);
        REQUIRE(std::abs(slide[Values] - 1.0) <= 1e-14);
    }
    // Sliding the cylinder along its axis moves no point, only its height v.
    const std::vector<double> along = rate(cylinder, 8);
    REQUIRE(std::abs(along[2]) <= 1e-14);
    REQUIRE(std::abs(along[5] + 1.0) <= 1e-14);
}

TEST_CASE("the crossings of an ellipsoid move as their closed form says",
          "[intersect][linearization]")
{
    // A line along x through the centre of an ellipsoid with semi-axes
    // (A, B, C) along x, y and z crosses it at x = -+A, where u is pi or 0
    // and v is 0: growing A moves them by -+1, B and C not at all.
    constexpr double a = 2.0;
    const Ellipsoid egg{Point3{0.0, 0.0, 0.0}, a, 1.25, 0.75, Vec3{1.0, 0.0, 0.0},
                        Vec3{0.0, 1.0, 0.0}};
    const Line line{Point3{-5.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    const std::size_t count = 6 + SurfaceParameters<Ellipsoid>();
    const auto rate = [&](const std::size_t k) {
        Model model = ModelFor<Ellipsoid>(std::vector<double>(count, 1.0));
        std::vector<double> seed(count, 0.0);
        seed[k] = 1.0;
        return Jvp(model, line, egg, seed);
    };
    const std::vector<double> grow = rate(9);
    REQUIRE(grow.size() == 2 * Values);
    REQUIRE(std::abs(grow[0] + 1.0) <= 1e-14);
    REQUIRE(std::abs(grow[3] + 1.0) <= 1e-14);
    REQUIRE(std::abs(grow[Values] - 1.0) <= 1e-14);
    REQUIRE(std::abs(grow[Values + 3] - 1.0) <= 1e-14);
    for (const std::size_t k : {std::size_t{10}, std::size_t{11}}) {
        for (const double value : rate(k)) {
            REQUIRE(std::abs(value) <= 1e-14);
        }
    }
    // Sliding it along y moves no point, and turns u at the crossings by
    // -+1 / A as the axis of the line leaves the centre.
    const std::vector<double> slide = rate(7);
    REQUIRE(std::abs(slide[1]) <= 1e-14);
    REQUIRE(std::abs(std::abs(slide[4]) - 1.0 / 1.25) <= 1e-14);
}

TEST_CASE("the crossings of a cone and a torus move as their closed forms say",
          "[intersect][linearization]")
{
    constexpr double pi = std::numbers::pi;
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    // The cone about z with radius 1 at its origin and semi-angle pi/4, and
    // a line along x at height 1/2: crossings at x = -+(1 + 1/2 tan a) =
    // -+3/2. Growing R or a widens it: x' = -+1 and -+(1/2) sec² a = -+1.
    // Raising its origin narrows it at the line's height: x' = +-tan a.
    const Cone cone{Point3{}, 1.0, pi / 4.0, ex, ey};
    const Line line{Point3{-5.0, 0.0, 0.5}, ex};
    const CurveSurfaceIntersection onCone = Crossings(line, cone, {3.5, 6.5});
    const std::vector<double> radius = RateAt(line, cone, onCone, 9);
    const std::vector<double> angle = RateAt(line, cone, onCone, 10);
    const std::vector<double> raise = RateAt(line, cone, onCone, 8);
    for (const auto& [rate, first] :
         {std::pair{radius, -1.0}, std::pair{angle, -1.0}, std::pair{raise, 1.0}}) {
        REQUIRE(rate.size() == 2 * Values);
        REQUIRE(std::abs(rate[0] - first) <= 1e-14);
        REQUIRE(std::abs(rate[3] - first) <= 1e-14);
        REQUIRE(std::abs(rate[Values] + first) <= 1e-14);
        REQUIRE(std::abs(rate[Values + 3] + first) <= 1e-14);
    }

    // The torus about z with radii 2 and 1, and the x axis: crossings at
    // x = -3, -1, 1, 3. Growing R moves all four outward; growing r moves
    // the outer ones outward and the inner ones in.
    const Torus torus{Point3{}, 2.0, 1.0, ex, ey};
    const Line axis{Point3{-5.0, 0.0, 0.0}, ex};
    const CurveSurfaceIntersection onTorus = Crossings(axis, torus, {2.0, 4.0, 6.0, 8.0});
    const std::vector<double> major = RateAt(axis, torus, onTorus, 9);
    const std::vector<double> minor = RateAt(axis, torus, onTorus, 10);
    REQUIRE(major.size() == 4 * Values);
    const std::array<double, 4> outward{-1.0, -1.0, 1.0, 1.0};
    const std::array<double, 4> tube{-1.0, 1.0, -1.0, 1.0};
    for (std::size_t k = 0; k < 4; ++k) {
        REQUIRE(std::abs(major[k * Values] - outward[k]) <= 1e-14);
        REQUIRE(std::abs(minor[k * Values] - tube[k]) <= 1e-14);
        REQUIRE(std::abs(minor[k * Values + 3] - tube[k]) <= 1e-14);
    }
}

TEST_CASE("the Taylor residual of a crossing decays with the square of the step",
          "[intersect][linearization][property]")
{
    rc::prop("for random lines and planes, seeds and references", [] {
        const auto [line, plane] = RandomPlaneCrossing();
        RequireTaylor(line, plane);
    });
    rc::prop("for random lines through spheres and cylinders", [] {
        const auto [line, quadrics] = RandomQuadricCrossing();
        RequireTaylor(line, quadrics.first);
        RequireTaylor(line, quadrics.second);
    });
}

/// A random line through the inside of a random ellipsoid in a random
/// frame, steep to it and away from its poles, as for the sphere.
std::pair<Line, Ellipsoid> RandomEllipsoidCrossing()
{
    const Point3 center{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 axis{Coordinate(), Coordinate(), Coordinate()};
    RC_PRE(Length(axis) >= 0.5);
    const Cylinder frame = IntersectTest::CylinderAbout(center, axis, 1.0);
    const auto semi = [] { return *rc::gen::inRange(8, 33) / 8.0; };
    const Ellipsoid egg{center, semi(), semi(), semi(), frame.XAxis, frame.YAxis};
    const double reach = std::min({egg.A, egg.B, egg.C});
    // Mostly across the axis, so both crossings stay off the poles.
    const Vec3 direction =
        Coordinate() * frame.XAxis + Coordinate() * frame.YAxis + 0.25 * Coordinate() * Axis(egg);
    RC_PRE(Length(direction) >= 0.5);
    const auto offset = [&] { return *rc::gen::inRange(-2, 3) / 8.0 * reach; };
    const Line line{center + Vec3{offset(), offset(), offset()} - 6.0 * direction, direction};
    RC_PRE(WellConditioned(line, egg, 2));
    RC_PRE(std::ranges::all_of(Meet(line, egg).Points, [](const CurveSurfacePoint& p) {
        return std::abs(p.OnSurface.Y) <= 1.2;
    }));
    return {line, egg};
}

TEST_CASE("ellipsoid crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][property]")
{
    rc::prop("for random lines through ellipsoids", [] {
        const auto [line, egg] = RandomEllipsoidCrossing();
        RequireTaylor(line, egg);
    });
    rc::prop("<J v, w> = <v, J^T w> for ellipsoids", [] {
        const auto [line, egg] = RandomEllipsoidCrossing();
        RequireTranspose(line, egg);
    });
}

TEST_CASE("the reverse product is the forward one's transpose", "[intersect][linearization]")
{
    rc::prop("<J v, w> = <v, J^T w> for planes", [] {
        const auto [line, plane] = RandomPlaneCrossing();
        RequireTranspose(line, plane);
    });
    rc::prop("<J v, w> = <v, J^T w> for spheres and cylinders", [] {
        const auto [line, quadrics] = RandomQuadricCrossing();
        RequireTranspose(line, quadrics.first);
        RequireTranspose(line, quadrics.second);
    });
}

TEST_CASE("a point without a derivative says so and carries none", "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    Model model = FullModel(std::vector<double>(9, 1.0));
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(9, 1.0)).has_value());
    const auto jvp = [&](const Line& line, const CurveSurfaceIntersection& at) {
        return Stereon::IntersectJvp(line, model.LineBindings, Surface{ground},
                                     model.SurfaceBindings, at, model.Sensitivity);
    };
    SECTION("a line in the plane: an overlap, Nondifferentiable")
    {
        const Line along{Point3{}, Vec3{1.0, 1.0, 0.0}};
        const auto result = jvp(along, Meet(along, ground));
        REQUIRE(result.has_value());
        REQUIRE(result->GetStatus() == DerivativeStatus::Nondifferentiable);
        REQUIRE_FALSE(result->GetValues().has_value());
    }
    const Line up{Point3{0.0, 0.0, -1.0}, Vec3{0.0, 0.0, 1.0}};
    CurveSurfaceIntersection crossing = Meet(up, ground);
    SECTION("a touch: Nondifferentiable")
    {
        crossing.Points.front().Contact = Stereon::ContactType::Tangential;
        const auto result = jvp(up, crossing);
        REQUIRE(result->GetStatus() == DerivativeStatus::Nondifferentiable);
        REQUIRE_FALSE(result->GetValues().has_value());
    }
    SECTION("a point the solver left undecided: Indeterminate")
    {
        crossing.Points.front().Decided = false;
        const auto result = jvp(up, crossing);
        REQUIRE(result->GetStatus() == DerivativeStatus::Indeterminate);
        REQUIRE_FALSE(result->GetValues().has_value());
    }
    SECTION("a grazing crossing: NearEvent, qualified, with its conditioning")
    {
        const Line grazing{Point3{-1.0, 0.0, -1e-7}, Vec3{1.0, 0.0, 1e-7}};
        const auto result = jvp(grazing, Meet(grazing, ground));
        REQUIRE(result->GetStatus() == DerivativeStatus::NearEvent);
        REQUIRE(result->IsQualified());
        REQUIRE(result->GetValues().has_value());
        REQUIRE(result->GetConditioning() < Stereon::NearEventSine);
    }
    SECTION("a clean crossing: Valid, its residual and conditioning reported")
    {
        const auto result = jvp(up, crossing);
        REQUIRE(result->GetStatus() == DerivativeStatus::Valid);
        REQUIRE(result->GetResidual() <= 1e-15);
        REQUIRE(std::abs(result->GetConditioning() - 1.0) <= 1e-15);
    }
}

TEST_CASE("a quadric's touch and a sphere's poles carry no derivative",
          "[intersect][linearization]")
{
    const Sphere unit{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    Model model = FullModel(std::vector<double>(10, 1.0));
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(10, 1.0)).has_value());
    const auto status = [&](const Line& line) {
        const auto result =
            Stereon::IntersectJvp(line, model.LineBindings, Surface{unit}, model.SurfaceBindings,
                                  Meet(line, unit), model.Sensitivity);
        REQUIRE(result.has_value());
        // A status carries values exactly when it has a derivative.
        REQUIRE(HasDerivative(result->GetStatus()) == result->GetValues().has_value());
        return result->GetStatus();
    };
    // A line at exactly the radius touches: the two crossings merge there.
    REQUIRE(status(Line{Point3{-2.0, 1.0, 0.0}, Vec3{1.0, 0.0, 0.0}}) ==
            DerivativeStatus::Nondifferentiable);
    // Through both poles: the longitude is not defined there, so neither is
    // its derivative, though the geometry is regular.
    REQUIRE(status(Line{Point3{0.0, 0.0, -2.0}, Vec3{0.0, 0.0, 1.0}}) ==
            DerivativeStatus::Indeterminate);
    // Through the equator: Valid.
    REQUIRE(status(Line{Point3{-2.0, 0.25, 0.0}, Vec3{1.0, 0.0, 0.0}}) == DerivativeStatus::Valid);
}

TEST_CASE("Linearize exposes the residual, the system and the regularity",
          "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // At 30 degrees to the plane: the conditioning is the sine, 1/2.
    const Line slant{Point3{0.0, 0.0, -1.0}, Vec3{std::numbers::sqrt3, 0.0, 1.0}};
    const CurveSurfaceIntersection at = Meet(slant, ground);
    const auto lin = Stereon::Linearize(slant, Surface{ground}, at.Points.front());
    REQUIRE(lin.has_value());
    REQUIRE(Length(lin->Residual) <= 1e-15);
    REQUIRE(std::abs(lin->Conditioning - 0.5) <= 1e-15);
    REQUIRE(lin->Status == DerivativeStatus::Valid);
    REQUIRE(lin->Tangent == slant.Direction);
    REQUIRE(lin->Su == ground.XAxis);
    REQUIRE(lin->Sv == ground.YAxis);
    CurveSurfacePoint broken = at.Points.front();
    broken.Parameter = std::nan("");
    REQUIRE(Stereon::Linearize(slant, Surface{ground}, broken).error().GetCode() ==
            ErrorCode::InvalidInput);
    // On a cylinder, a crossing straight through the axis is perpendicular.
    const Cylinder unit{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Line across{Point3{-2.0, 0.0, 0.5}, Vec3{1.0, 0.0, 0.0}};
    for (const CurveSurfacePoint& p : Meet(across, unit).Points) {
        const auto onCylinder = Stereon::Linearize(across, Surface{unit}, p);
        REQUIRE(std::abs(onCylinder->Conditioning - 1.0) <= 1e-15);
        REQUIRE(Length(onCylinder->Residual) <= 1e-15);
    }
}

TEST_CASE("invalid sensitivity input is refused, and other surfaces are NotSupported",
          "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Line up{Point3{0.0, 0.0, -1.0}, Vec3{0.0, 0.0, 1.0}};
    const CurveSurfaceIntersection at = Meet(up, ground);
    Model model = FullModel(std::vector<double>(9, 1.0));
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    // No seed yet; no adjoint, or one of the wrong length.
    REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, Surface{ground},
                                       model.SurfaceBindings, at, model.Sensitivity)) ==
            ErrorCode::InvalidInput);
    REQUIRE(model.Sensitivity.SetAdjoint(std::vector<double>(5, 0.0)).has_value());
    REQUIRE(code(Stereon::IntersectVjp(up, model.LineBindings, Surface{ground},
                                       model.SurfaceBindings, at, model.Sensitivity)) ==
            ErrorCode::InvalidInput);
    // A parameter of another set, and an angle moving a coordinate.
    ParameterSet other;
    const ParameterId foreign = *other.Add("q", ParameterQuantity::Length);
    const ParameterId angle = *other.Add("a", ParameterQuantity::Angle);
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(9, 1.0)).has_value());
    for (const ParameterId id : {foreign, angle}) {
        const std::vector<LineBinding> bad{{id, LineComponent::OriginX}};
        REQUIRE(code(Stereon::IntersectJvp(up, bad, Surface{ground}, model.SurfaceBindings, at,
                                           model.Sensitivity)) == ErrorCode::InvalidInput);
    }
    // A plane has no radius to bind.
    const std::vector<SurfaceBinding> radius{
        {model.SurfaceBindings.front().Parameter, SurfaceComponent::Radius}};
    REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, Surface{ground}, radius, at,
                                       model.Sensitivity)) == ErrorCode::InvalidInput);
    // A cone has no minor radius, and its semi-angle is an Angle.
    const Surface cone(Cone{Point3{}, 1.0, 0.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    const ParameterId length = model.SurfaceBindings.front().Parameter;
    for (const SurfaceComponent component :
         {SurfaceComponent::MinorRadius, SurfaceComponent::SemiAngle}) {
        const std::vector<SurfaceBinding> bad{{length, component}};
        REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, cone, bad,
                                           CurveSurfaceIntersection{}, model.Sensitivity)) ==
                ErrorCode::InvalidInput);
    }
    // A B-spline surface's poles: none past its last, no radius on it, a
    // weight moved by a dimensionless parameter, and no pole on a plane.
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto patch = Stereon::MakeBSplineSurface(
        1, 1, knots, knots, {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}},
        {});
    REQUIRE(patch.has_value());
    for (const SurfaceBinding& bad : {SurfaceBinding{length, SurfaceComponent::PoleZ, 4},
                                      SurfaceBinding{length, SurfaceComponent::Radius},
                                      SurfaceBinding{length, SurfaceComponent::PoleWeight, 0}}) {
        const std::vector<SurfaceBinding> bindings{bad};
        REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, Surface{*patch}, bindings,
                                           CurveSurfaceIntersection{}, model.Sensitivity)) ==
                ErrorCode::InvalidInput);
    }
    const std::vector<SurfaceBinding> pole{{length, SurfaceComponent::PoleZ, 0}};
    REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, Surface{ground}, pole, at,
                                       model.Sensitivity)) == ErrorCode::InvalidInput);
    // A surface with no bindings, an extrusion: a NotSupported status,
    // never a value.
    const Surface extrusion{Stereon::ExtrusionSurface(
        Stereon::Curve{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}}, Vec3{0.0, 1.0, 0.0})};
    const auto unsupported = Stereon::IntersectJvp(up, model.LineBindings, extrusion, {},
                                                   CurveSurfaceIntersection{}, model.Sensitivity);
    REQUIRE(unsupported.has_value());
    REQUIRE(unsupported->GetStatus() == DerivativeStatus::NotSupported);
    REQUIRE_FALSE(unsupported->GetValues().has_value());
}

namespace {

using Stereon::BSplineSurface;

/// The bilinear patch over the unit square in the plane z = 0, pole (i, j)
/// at (i, j, 0), rational with weights (1, 2, 3, 4) when asked.
BSplineSurface FlatPatch(const bool rational)
{
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    return *Stereon::MakeBSplineSurface(
        1, 1, knots, knots, {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}},
        rational ? std::vector<double>{1.0, 2.0, 3.0, 4.0} : std::vector<double>{});
}

/// A bicubic patch over [0, 4] x [0, 3], its heights random, its weights
/// random in [0.5, 2] when rational.
BSplineSurface RandomPatch(const bool rational)
{
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            poles.push_back({i / 3.0 * 4.0, j / 3.0 * 3.0, *rc::gen::inRange(-5, 6) / 10.0});
            if (rational) {
                weights.push_back(*rc::gen::inRange(5, 21) / 10.0);
            }
        }
    }
    return *Stereon::MakeBSplineSurface(3, 3, knots, knots, std::move(poles), std::move(weights));
}

/// The line's six bindings, then each pole's components: x, y, z, and the
/// weight of a rational patch.
Model SplineSurfaceModel(const std::vector<double>& references,
                         const std::vector<std::size_t>& poles, const bool rational)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    const std::size_t perPole = rational ? 4 : 3;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const bool weight = k >= 6 && (k - 6) % perPole == 3;
        ids.push_back(*set.Add(
            "p" + std::to_string(k),
            weight ? ParameterQuantity::Dimensionless : ParameterQuantity::Length, references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < 6; ++k) {
        model.LineBindings.push_back({ids[k], static_cast<LineComponent>(k)});
    }
    for (std::size_t k = 6; k < references.size(); ++k) {
        const std::size_t c = (k - 6) % perPole;
        model.SurfaceBindings.push_back(
            {ids[k], static_cast<SurfaceComponent>(static_cast<int>(SurfaceComponent::PoleX) + c),
             poles[(k - 6) / perPole]});
    }
    return model;
}

/// The line and the patch moved by `step` times the seed.
std::pair<Line, BSplineSurface> MovedPatch(Line line, BSplineSurface patch, const Model& model,
                                           const std::vector<double>& seed,
                                           const std::vector<double>& references, const double step)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    line.Origin = line.Origin + Vec3{d(0), d(1), d(2)};
    line.Direction = line.Direction + Vec3{d(3), d(4), d(5)};
    for (std::size_t k = 0; k < model.SurfaceBindings.size(); ++k) {
        const SurfaceBinding& b = model.SurfaceBindings[k];
        Point3& pole = patch.Poles[b.Pole];
        switch (b.Component) {
        case SurfaceComponent::PoleX:
            pole.X += d(6 + k);
            break;
        case SurfaceComponent::PoleY:
            pole.Y += d(6 + k);
            break;
        case SurfaceComponent::PoleZ:
            pole.Z += d(6 + k);
            break;
        default:
            patch.Weights[b.Pole] += d(6 + k);
            break;
        }
    }
    return {line, patch};
}

} // namespace

TEST_CASE("the crossing of a B-spline surface moves with its poles as its closed form says",
          "[intersect][linearization][spline-surface]")
{
    // A vertical line through (0.3, 0.4) crosses the flat patch at u = 0.3,
    // v = 0.4. Raising pole (1, 1), index 3, by p makes the patch's height
    // p R_11(u, v): the point rises by R_11 = u v w_11 / W, t with it, and
    // (u, v) stay, the patch's x and y being unmoved.
    const Line up{Point3{0.3, 0.4, -1.0}, Vec3{0.0, 0.0, 1.0}};
    for (const bool rational : {false, true}) {
        INFO((rational ? "rational" : "polynomial"));
        const BSplineSurface patch = FlatPatch(rational);
        const auto at = Stereon::Intersect(Stereon::Curve{up}, Surface{patch}, ctx);
        REQUIRE(at.has_value());
        REQUIRE(at->Points.size() == 1u);
        ParameterSet set;
        const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
        SensitivityContext sensitivity(std::move(set));
        REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
        const std::vector<SurfaceBinding> lift{{p, SurfaceComponent::PoleZ, 3}};
        const auto jvp = Stereon::IntersectJvp(up, {}, Surface{patch}, lift, *at, sensitivity);
        REQUIRE(jvp.has_value());
        REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
        const auto values = *jvp->GetValues();
        REQUIRE(values.size() == Values);
        const Stereon::Point2 uv = at->Points[0].OnSurface;
        const Vec3 rate = Stereon::PoleVelocity(patch, 3, Stereon::PoleComponent::Z, uv);
        const double w = rational ? 4.0 / (1.0 * (1 - uv.X) * (1 - uv.Y) + 2.0 * (1 - uv.X) * uv.Y +
                                           3.0 * uv.X * (1 - uv.Y) + 4.0 * uv.X * uv.Y)
                                  : 1.0;
        CHECK(std::abs(rate.Z - uv.X * uv.Y * w) <= 1e-15);
        CHECK(std::abs(values[0]) <= 1e-15);
        CHECK(std::abs(values[1]) <= 1e-15);
        CHECK(std::abs(values[2] - rate.Z) <= 1e-15);
        CHECK(std::abs(values[3] - rate.Z) <= 1e-15);
        CHECK(std::abs(values[4]) <= 1e-15);
        CHECK(std::abs(values[5]) <= 1e-15);
    }
}

TEST_CASE("B-spline surface crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][spline-surface][property]")
{
    rc::prop("in the line, and in the poles and weights of the patch", [] {
        const bool rational = *rc::gen::arbitrary<bool>();
        const BSplineSurface patch = RandomPatch(rational);
        const double u = *rc::gen::inRange(2, 9) / 10.0;
        const double v = *rc::gen::inRange(2, 9) / 10.0;
        const int k = *rc::gen::inRange(0, 64);
        const Vec3 d{0.15 * std::cos(k), 0.15 * std::sin(k), 1.0};
        const Point3 on = Eval(Surface{patch}, Stereon::Point2{u, v});
        const Line line{on - 2.0 * d, d};
        const std::size_t count = Meet(line, patch).Points.size();
        RC_PRE(count >= 1u && WellConditioned(line, patch, count));
        // Every point inside the domain, so that a step leaves the count.
        for (const CurveSurfacePoint& p : Meet(line, patch).Points) {
            RC_PRE(p.OnSurface.X > 0.05 && p.OnSurface.X < 0.95 && p.OnSurface.Y > 0.05 &&
                   p.OnSurface.Y < 0.95);
        }
        std::vector<std::size_t> poles;
        for (int n = 0; n < 3; ++n) {
            poles.push_back(*rc::gen::inRange<std::size_t>(0, patch.Poles.size()));
        }
        const auto [references, seed] = RandomSeed(6 + poles.size() * (rational ? 4 : 3));
        Model model = SplineSurfaceModel(references, poles, rational);
        const std::vector<double> jv = Jvp(model, line, patch, seed);
        const std::vector<double> base = Points(line, patch);
        const auto residual = [&](const double h) {
            const auto [l, s] = MovedPatch(line, patch, model, seed, references, h);
            const std::vector<double> moved = Points(l, s);
            RC_ASSERT(moved.size() == base.size());
            std::vector<double> r(base.size());
            for (std::size_t i = 0; i < base.size(); ++i) {
                r[i] = moved[i] - base[i] - h * jv[i];
            }
            return Norm(r);
        };
        const double coarse = residual(1e-2);
        const double fine = residual(1e-3);
        RC_ASSERT(fine <= coarse / 30.0 + 1e-10 * (1.0 + Norm(base)));

        std::vector<double> w;
        for (std::size_t i = 0; i < jv.size(); ++i) {
            w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        }
        RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
        const auto vjp =
            Stereon::IntersectVjp(line, model.LineBindings, Surface{patch}, model.SurfaceBindings,
                                  Meet(line, patch), model.Sensitivity);
        RC_ASSERT(vjp.has_value());
        const auto jtw = *vjp->GetValues();
        double left = 0.0;
        double right = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < jv.size(); ++i) {
            left += jv[i] * w[i];
            scale += std::abs(jv[i] * w[i]);
        }
        for (std::size_t i = 0; i < seed.size(); ++i) {
            right += seed[i] * jtw[i];
            scale += std::abs(seed[i] * jtw[i]);
        }
        RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
    });
}

TEST_CASE("JVP and VJP reject intersection results from stale geometry",
          "[intersect][linearization][regression]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Point3 offset{0.0, 0.0, 0.01};
    const Plane shifted{offset, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Line up{Point3{0.0, 0.0, -1.0}, Vec3{0.0, 0.0, 1.0}};
    const CurveSurfaceIntersection at = Meet(up, ground);
    Model model = FullModel(std::vector<double>(9, 1.0));
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(9, 1.0)).has_value());
    REQUIRE(model.Sensitivity.SetAdjoint(std::vector<double>(6, 1.0)).has_value());

    const auto jvp = Stereon::IntersectJvp(up, model.LineBindings, Surface{shifted},
                                           model.SurfaceBindings, at, model.Sensitivity);
    REQUIRE_FALSE(jvp.has_value());
    REQUIRE(jvp.error().GetCode() == ErrorCode::InvalidInput);
    const auto vjp = Stereon::IntersectVjp(up, model.LineBindings, Surface{shifted},
                                           model.SurfaceBindings, at, model.Sensitivity);
    REQUIRE_FALSE(vjp.has_value());
    REQUIRE(vjp.error().GetCode() == ErrorCode::InvalidInput);

    CurveSurfaceIntersection altered = at;
    altered.Points.front().Position.Z = 0.01;
    const auto stored = Stereon::IntersectJvp(up, model.LineBindings, Surface{ground},
                                              model.SurfaceBindings, altered, model.Sensitivity);
    REQUIRE_FALSE(stored.has_value());
    REQUIRE(stored.error().GetCode() == ErrorCode::InvalidInput);

    const auto good = Stereon::IntersectJvp(up, model.LineBindings, Surface{ground},
                                            model.SurfaceBindings, at, model.Sensitivity);
    REQUIRE(good.has_value());
    REQUIRE(good->GetStatus() == DerivativeStatus::Valid);
}

namespace {

using Stereon::BSplineCurve;
using Stereon::PoleBinding;
using Stereon::PoleComponent;

/// A B-spline model: parameters bound to poles and weights, then to the
/// surface's origin and radius.
struct SplineModel {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<PoleBinding> PoleBindings;
    std::vector<SurfaceBinding> SurfaceBindings;
};

SplineModel MakeSplineModel(const std::vector<PoleBinding>& poles,
                            const std::vector<SurfaceComponent>& surface,
                            const std::vector<double>& references)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const bool weight = k < poles.size() && poles[k].Component == PoleComponent::Weight;
        const ParameterQuantity quantity = k >= poles.size() ? QuantityOf(surface[k - poles.size()])
                                           : weight          ? ParameterQuantity::Dimensionless
                                                             : ParameterQuantity::Length;
        ids.push_back(*set.Add("p" + std::to_string(k), quantity, references[k]));
    }
    SplineModel model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < poles.size(); ++k) {
        model.PoleBindings.push_back({ids[k], poles[k].Pole, poles[k].Component});
    }
    for (std::size_t k = 0; k < surface.size(); ++k) {
        model.SurfaceBindings.push_back({ids[poles.size() + k], surface[k]});
    }
    return model;
}

/// The curve and surface moved by `step` times the seed.
template <typename S>
std::pair<BSplineCurve, S> MovedSpline(BSplineCurve curve, S surface, const SplineModel& model,
                                       const std::vector<double>& seed,
                                       const std::vector<double>& references, const double step,
                                       const bool turning = false)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    for (std::size_t k = 0; k < model.PoleBindings.size(); ++k) {
        const PoleBinding& b = model.PoleBindings[k];
        Point3& pole = curve.Poles[b.Pole];
        switch (b.Component) {
        case PoleComponent::X:
            pole.X += d(k);
            break;
        case PoleComponent::Y:
            pole.Y += d(k);
            break;
        case PoleComponent::Z:
            pole.Z += d(k);
            break;
        case PoleComponent::Weight:
            curve.Weights[b.Pole] += d(k);
            break;
        }
    }
    MoveSurface(surface, d, model.PoleBindings.size(), turning);
    return {curve, surface};
}

template <typename S>
std::vector<double> SplineJvp(SplineModel& model, const BSplineCurve& curve, const S& surface,
                              const std::vector<double>& seed)
{
    RC_ASSERT(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp =
        Stereon::IntersectJvp(curve, model.PoleBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(curve, surface), model.Sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

/// A random cubic B-spline wandering through `center` along `direction`,
/// rational half the time, with random bindings to its poles and weights.
std::pair<BSplineCurve, std::vector<PoleBinding>>
RandomSplineThrough(const Point3& center, const Vec3& direction, const double spread)
{
    const bool rational = *rc::gen::arbitrary<bool>();
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 5; ++i) {
        const auto wiggle = [&] { return *rc::gen::inRange(-4, 5) / 16.0 * spread; };
        poles.push_back(center + static_cast<double>(i - 2) * direction +
                        Vec3{wiggle(), wiggle(), wiggle()});
        if (rational) {
            weights.push_back(*rc::gen::inRange(6, 15) / 10.0);
        }
    }
    const auto made =
        Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0}, poles, weights);
    RC_PRE(made.has_value());
    std::vector<PoleBinding> bindings;
    const int count = *rc::gen::inRange(1, 5);
    for (int k = 0; k < count; ++k) {
        const int top = rational ? 4 : 3;
        bindings.push_back({ParameterId{}, *rc::gen::inRange<std::size_t>(0, 5),
                            static_cast<PoleComponent>(*rc::gen::inRange(0, top))});
    }
    return {*made, bindings};
}

/// Every point of a B-spline's intersection decided and well conditioned,
/// and off a sphere's poles.
template <typename S>
bool SplineWellConditioned(const BSplineCurve& curve, const S& surface)
{
    // The primal's own completeness is tested elsewhere: a curve it cannot
    // resolve is no sample for the products.
    const auto result = Intersect(Stereon::Curve{curve}, Surface{surface}, ctx);
    if (!result || result->Points.empty() || !result->Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(result->Points, [&](const CurveSurfacePoint& p) {
        const auto lin = Stereon::Linearize(curve, Surface{surface}, p);
        bool offPole = true;
        if constexpr (std::is_same_v<S, Sphere>) {
            offPole = std::abs(p.OnSurface.Y) <= 1.2;
        } else if constexpr (std::is_same_v<S, Cone>) {
            offPole =
                std::abs(surface.Radius + p.OnSurface.Y * std::sin(surface.SemiAngle)) >= 0.25;
        }
        // Away from the curve's ends, so a step never carries a crossing
        // out of its domain: a cone, being unbounded, can lie near one.
        const Stereon::Interval domain = Domain(curve);
        const double margin = 0.025 * Width(domain);
        const bool inside =
            p.Parameter >= domain.GetLower() + margin && p.Parameter <= domain.GetUpper() - margin;
        return p.Decided && offPole && inside && lin && lin->Conditioning >= 0.25;
    });
}

/// Taylor decay and the transpose identity for a B-spline against a surface.
template <typename S>
void RequireSplineProducts(const BSplineCurve& curve, const std::vector<PoleBinding>& poles,
                           const S& surface, const bool turning = false)
{
    RC_PRE(SplineWellConditioned(curve, surface));
    const std::size_t count = poles.size() + SurfaceParameters<S>(turning);
    const auto [references, seed] = RandomSeed(count);
    SplineModel model = MakeSplineModel(poles, ComponentsOf<S>(turning), references);
    const std::vector<double> jv = SplineJvp(model, curve, surface, seed);
    const std::vector<double> base = Points(curve, surface);
    const auto residual = [&](const double h) {
        const auto [c, s] = MovedSpline(curve, surface, model, seed, references, h, turning);
        // A step that turns a near miss into crossings, or into a touch the
        // primal cannot resolve, crosses an event of the intersection: no
        // sample for its derivative. A cone, unbounded, offers near misses
        // far from the crossings the products are of.
        const bool resolved = Intersect(Stereon::Curve{c}, Surface{s}, ctx).has_value();
        RC_PRE(resolved);
        const std::vector<double> moved = Points(c, s);
        RC_PRE(moved.size() == base.size());
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < base.size(); ++i) {
            const double difference = moved[i] - base[i];
            // u is an angle, and so is a torus's v; a step never moves
            // any other v by half a turn, so v is taken modulo one too.
            const bool angle = i % Values == 4 || i % Values == 5;
            r[i] = (angle ? std::remainder(difference, tau) : difference) - h * jv[i];
        }
        return Norm(r);
    };
    RC_ASSERT(residual(1e-3) <= residual(1e-2) / 30.0 + 1e-9 * (1.0 + Norm(base)));
    // The transpose, on the same model.
    std::vector<double> w;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
    const auto vjp =
        Stereon::IntersectVjp(curve, model.PoleBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(curve, surface), model.Sensitivity);
    RC_ASSERT(vjp.has_value());
    const auto jtw = *vjp->GetValues();
    double left = 0.0;
    double right = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        left += jv[i] * w[i];
        scale += std::abs(jv[i] * w[i]);
    }
    for (std::size_t k = 0; k < seed.size(); ++k) {
        right += seed[k] * jtw[k];
        scale += std::abs(seed[k] * jtw[k]);
    }
    RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
}

} // namespace

TEST_CASE("a straight B-spline moves as the line it draws", "[intersect][linearization]")
{
    // Poles P0 and P1 on [0, 1] draw the line from P0 along P1 - P0: moving
    // P1 is moving the line's direction, through every surface.
    const Point3 p0{-2.0, 0.25, 0.125};
    const Point3 p1{2.0, -0.5, 0.375};
    const auto made = Stereon::MakeBSplineCurve(1, {0.0, 0.0, 1.0, 1.0}, {p0, p1});
    REQUIRE(made.has_value());
    const Line line{p0, p1 - p0};
    const auto compare = [&](const auto& surface) {
        for (std::size_t k = 0; k < 3; ++k) {
            ParameterSet set;
            const ParameterId id = *set.Add("p", ParameterQuantity::Length);
            SensitivityContext sensitivity(std::move(set));
            REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
            const std::vector<PoleBinding> pole{{id, 1, static_cast<PoleComponent>(k)}};
            const std::vector<LineBinding> direction{{id, static_cast<LineComponent>(3 + k)}};
            const auto fromSpline = Stereon::IntersectJvp(*made, pole, Surface{surface}, {},
                                                          Meet(*made, surface), sensitivity);
            const auto fromLine = Stereon::IntersectJvp(line, direction, Surface{surface}, {},
                                                        Meet(line, surface), sensitivity);
            REQUIRE(fromSpline->GetStatus() == DerivativeStatus::Valid);
            const auto a = *fromSpline->GetValues();
            const auto b = *fromLine->GetValues();
            REQUIRE(a.size() == b.size());
            for (std::size_t i = 0; i < a.size(); ++i) {
                REQUIRE(std::abs(a[i] - b[i]) <= 1e-12);
            }
        }
    };
    compare(Plane{Point3{0.0, 0.0, 0.25}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    compare(Sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    compare(Cylinder{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
}

TEST_CASE("B-spline crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][property]")
{
    rc::prop("for random cubic B-splines through spheres, cylinders and planes", [] {
        const Point3 center{Coordinate(), Coordinate(), Coordinate()};
        const double radius = *rc::gen::inRange(8, 33) / 8.0;
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        // The end poles lie two steps from the centre, well outside the
        // radius, so a small step never changes how many crossings there are.
        RC_PRE(Length(direction) >= radius);
        const auto [curve, poles] = RandomSplineThrough(center, direction, radius);
        const Vec3 axis = Cross(direction, Vec3{Coordinate(), Coordinate(), Coordinate()});
        RC_PRE(Length(axis) >= 0.25 * Length(direction));
        switch (*rc::gen::inRange(0, 3)) {
        case 0:
            RequireSplineProducts(curve, poles,
                                  Sphere{center, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
            break;
        case 1:
            RequireSplineProducts(curve, poles, IntersectTest::CylinderAbout(center, axis, radius));
            break;
        default:
            RequireSplineProducts(curve, poles, IntersectTest::PlaneWithNormal(center, direction));
            break;
        }
    });
}

TEST_CASE("B-spline crossings of cones and tori: Taylor decay and the transpose",
          "[intersect][linearization][property]")
{
    const auto frameAbout = [](const Point3& center) {
        const Vec3 axis{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(axis) >= 0.5);
        return IntersectTest::CylinderAbout(center, axis, 1.0);
    };
    rc::prop("for random cubic B-splines through cones", [&] {
        const Point3 center{Coordinate(), Coordinate(), Coordinate()};
        const Cylinder frame = frameAbout(center);
        const double radius = *rc::gen::inRange(8, 25) / 8.0;
        const double sign = *rc::gen::arbitrary<bool>() ? 1.0 : -1.0;
        const double angle = sign * *rc::gen::inRange(2, 11) / 10.0;
        const Cone cone{center, radius, angle, frame.XAxis, frame.YAxis};
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(direction) >= radius);
        const auto [curve, poles] = RandomSplineThrough(center, direction, radius);
        RequireSplineProducts(curve, poles, cone);
    });
    rc::prop("for random cubic B-splines through tori", [&] {
        const Point3 center{Coordinate(), Coordinate(), Coordinate()};
        const Cylinder frame = frameAbout(center);
        const double minor = *rc::gen::inRange(4, 9) / 8.0;
        const double major = minor + *rc::gen::inRange(4, 13) / 8.0;
        const Torus torus{center, major, minor, frame.XAxis, frame.YAxis};
        // Through a point of the tube's centre circle, so the tube is
        // crossed; the end poles, two steps away, lie outside the torus.
        const double u = *rc::gen::inRange(0, 16) / 16.0 * tau;
        const Point3 tube =
            center + major * (std::cos(u) * frame.XAxis + std::sin(u) * frame.YAxis);
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(direction) >= major + minor);
        const auto [curve, poles] = RandomSplineThrough(tube, direction, minor);
        RequireSplineProducts(curve, poles, torus);
    });
}

TEST_CASE("a B-spline's touch and bad pole bindings", "[intersect][linearization]")
{
    // z = (2 s - 1)²: the primal now fails closed on unresolved touches.
    // Check derivative status with an explicitly supplied candidate.
    const auto made = Stereon::MakeBSplineCurve(
        2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {{0.0, 0.0, 1.0}, {1.0, 0.0, -1.0}, {2.0, 0.0, 1.0}});
    REQUIRE(made.has_value());
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    ParameterSet set;
    const ParameterId lift = *set.Add("lift", ParameterQuantity::Length);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<PoleBinding> middle{{lift, 1, PoleComponent::Z}};
    CurveSurfaceIntersection primal;
    CurveSurfacePoint tangency{};
    tangency.Parameter = 0.5;
    tangency.OnSurface = {1.0, 0.0};
    tangency.Position = {1.0, 0.0, 0.0};
    tangency.Contact = Stereon::ContactType::Tangential;
    tangency.Decided = false;
    primal.Points.push_back(tangency);
    const auto touch =
        Stereon::IntersectJvp(*made, middle, Surface{ground}, {}, primal, sensitivity);
    REQUIRE(touch.has_value());
    REQUIRE(touch->GetStatus() == DerivativeStatus::Nondifferentiable);
    REQUIRE_FALSE(touch->GetValues().has_value());
    // A pole the curve does not have, and a weight moved by a length.
    for (const PoleBinding& bad :
         {PoleBinding{lift, 7, PoleComponent::Z}, PoleBinding{lift, 1, PoleComponent::Weight}}) {
        const std::vector<PoleBinding> one{bad};
        const auto invalid =
            Stereon::IntersectJvp(*made, one, Surface{ground}, {}, primal, sensitivity);
        REQUIRE_FALSE(invalid.has_value());
        REQUIRE(invalid.error().GetCode() == ErrorCode::InvalidInput);
    }
    // A point off the curve's domain.
    CurveSurfacePoint outside = primal.Points.front();
    outside.Parameter = 5.0;
    REQUIRE(Stereon::Linearize(*made, Surface{ground}, outside).error().GetCode() ==
            ErrorCode::InvalidInput);
}

TEST_CASE("Intersect declares its forward and reverse products", "[intersect][linearization]")
{
    const auto* capability = Stereon::FindCapability(Stereon::IntersectCapabilities(), "Intersect");
    REQUIRE(capability != nullptr);
    REQUIRE(capability->Supports(Stereon::DerivativeProduct::Jvp));
    REQUIRE(capability->Supports(Stereon::DerivativeProduct::Vjp));
    REQUIRE_FALSE(capability->Supports(Stereon::DerivativeProduct::DenseJacobian));
}

TEST_CASE("radius products are independent of world translation",
          "[intersect][linearization][regression]")
{
    const double angle = std::numbers::pi / 4.0;
    const Vec3 radial{std::cos(angle), std::sin(angle), 0.0};
    for (const double shift : {0.0, 0x1p30, 0x1p40, 0x1p54}) {
        const Point3 center{shift, 0.0, 0.0};
        const Sphere sphere{center, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const Cylinder cylinder{center, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        for (const Surface& surface : {Surface{sphere}, Surface{cylinder}}) {
            CAPTURE(shift);
            const Line line{center, radial};
            // Supply an analytic primal contact: avoid folding the intersection
            // solver's placement error into a derivative regression.
            CurveSurfaceIntersection at;
            at.Points.push_back(
                {1.0, {angle, 0.0}, Eval(line, 1.0), Stereon::ContactType::Crossing, true});
            Model model = FullModel(std::vector<double>(10, 1.0));
            std::vector<double> seed(10, 0.0);
            seed[9] = 1.0;
            REQUIRE(model.Sensitivity.SetSeed(seed).has_value());
            const auto jvp = Stereon::IntersectJvp(line, model.LineBindings, surface,
                                                   model.SurfaceBindings, at, model.Sensitivity);
            REQUIRE(jvp.has_value());
            REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
            const auto values = *jvp->GetValues();
            REQUIRE(std::abs(values[0] - radial.X) < 1e-12);
            REQUIRE(std::abs(values[1] - radial.Y) < 1e-12);
            REQUIRE(std::abs(values[2]) < 1e-12);
            REQUIRE(std::abs(values[3] - 1.0) < 1e-12);
            REQUIRE(std::abs(values[4]) < 1e-12);
            REQUIRE(std::abs(values[5]) < 1e-12);
            // Adjoint of t: its radius component must independently be one.
            const std::array<double, 6> adjoint{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            REQUIRE(model.Sensitivity.SetAdjoint(adjoint).has_value());
            const auto vjp = Stereon::IntersectVjp(line, model.LineBindings, surface,
                                                   model.SurfaceBindings, at, model.Sensitivity);
            REQUIRE(vjp.has_value());
            REQUIRE(vjp->GetStatus() == DerivativeStatus::Valid);
            REQUIRE(std::abs((*vjp->GetValues())[9] - 1.0) < 1e-12);
        }
    }
}

// --- Circles and ellipses ----------------------------------------------------

namespace {

using Stereon::Circle;
using Stereon::ConicBinding;
using Stereon::ConicComponent;
using Stereon::Ellipse;

/// The conic's bindings, its centre then its radii, in order; when
/// `turning`, then its frame's turn about X, Y and X x Y.
template <typename C>
std::vector<ConicComponent> ConicComponentsOf(const bool turning = false)
{
    std::vector<ConicComponent> components{ConicComponent::CenterX, ConicComponent::CenterY,
                                           ConicComponent::CenterZ};
    if constexpr (std::is_same_v<C, Circle>) {
        components.push_back(ConicComponent::Radius);
    } else {
        components.insert(components.end(),
                          {ConicComponent::MajorRadius, ConicComponent::MinorRadius});
    }
    if (turning) {
        components.insert(components.end(), {ConicComponent::RotationX, ConicComponent::RotationY,
                                             ConicComponent::RotationZ});
    }
    return components;
}

struct ConicModel {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<ConicBinding> ConicBindings;
    std::vector<SurfaceBinding> SurfaceBindings;
};

/// The conic's bindings, then the surface's, one parameter each.
template <typename C, typename S>
ConicModel ConicModelFor(const std::vector<double>& references, const bool turning = false)
{
    const std::vector<ConicComponent> conic = ConicComponentsOf<C>(turning);
    const std::vector<SurfaceComponent> surface = ComponentsOf<S>(turning);
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const bool turn = k < conic.size() && conic[k] >= ConicComponent::RotationX;
        const ParameterQuantity quantity = k >= conic.size() ? QuantityOf(surface[k - conic.size()])
                                           : turn            ? ParameterQuantity::Angle
                                                             : ParameterQuantity::Length;
        ids.push_back(*set.Add("p" + std::to_string(k), quantity, references[k]));
    }
    ConicModel model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < conic.size(); ++k) {
        model.ConicBindings.push_back({ids[k], conic[k]});
    }
    for (std::size_t k = 0; k < surface.size(); ++k) {
        model.SurfaceBindings.push_back({ids[conic.size() + k], surface[k]});
    }
    return model;
}

/// The conic and the surface moved by `step` times the seed.
template <typename C, typename S>
std::pair<C, S> MovedConic(C conic, S surface, const std::vector<double>& seed,
                           const std::vector<double>& references, const double step,
                           const bool turning = false)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    conic.Center = conic.Center + Vec3{d(0), d(1), d(2)};
    std::size_t next = 4;
    if constexpr (std::is_same_v<C, Circle>) {
        conic.Radius += d(3);
    } else {
        conic.MajorRadius += d(3);
        conic.MinorRadius += d(4);
        next = 5;
    }
    if (turning) {
        Turn(conic.XAxis, conic.YAxis, Cross(conic.XAxis, conic.YAxis), d(next), d(next + 1),
             d(next + 2));
        next += 3;
    }
    MoveSurface(surface, d, next, turning);
    return {conic, surface};
}

template <typename C, typename S>
std::vector<double> ConicJvp(ConicModel& model, const C& conic, const S& surface,
                             const std::vector<double>& seed)
{
    RC_ASSERT(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp =
        Stereon::IntersectJvp(conic, model.ConicBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(conic, surface), model.Sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

/// Every crossing well away from a tangency and from the surface's poles:
/// near a sphere's pole the angle u moves fast and bends sharply with the
/// geometry, a derivative too nonlinear for the steps below, though Valid.
/// The partials within a factor of ten of each other keep it away.
template <typename C, typename S>
bool ConicWellConditioned(const C& conic, const S& surface)
{
    const CurveSurfaceIntersection at = Meet(conic, surface);
    if (at.Points.empty() || !at.Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(at.Points, [&](const CurveSurfacePoint& p) {
        const auto lin = Stereon::Linearize(conic, Surface{surface}, p);
        return lin && lin->Conditioning >= 0.25 && Length(lin->Su) >= 0.1 * Length(lin->Sv) &&
               Length(lin->Sv) >= 0.1 * Length(lin->Su);
    });
}

/// The Taylor residual falls by second order, and <J v, w> = <v, J^T w>.
/// Angles (the conic's t, the surface's u and v) are compared modulo a
/// turn: a step may carry one over a seam.
template <typename C, typename S>
void RequireConicProducts(const C& conic, const S& surface, const bool turning = false)
{
    const auto [references, seed] =
        RandomSeed(ConicComponentsOf<C>(turning).size() + SurfaceParameters<S>(turning));
    ConicModel model = ConicModelFor<C, S>(references, turning);
    const std::vector<double> jv = ConicJvp(model, conic, surface, seed);
    const std::vector<double> base = Points(conic, surface);
    // A conic meets a quadric up to four times, a line twice: a conic that
    // passes near the surface away from its crossings gains two under a
    // step, a new regime, not an error of the derivative. Such steps are
    // discarded; the crossings that exist were checked transverse.
    for (const double h : {1e-2, 1e-3}) {
        const auto [c, s] = MovedConic(conic, surface, seed, references, h, turning);
        RC_PRE(Points(c, s).size() == base.size());
    }
    // Points come in order of the conic's angle, so a crossing near t = 0
    // that a step carries past the seam comes last: each base point is
    // paired with the moved point nearest it in space.
    const auto residual = [&](const double h) {
        const auto [c, s] = MovedConic(conic, surface, seed, references, h, turning);
        const std::vector<double> moved = Points(c, s);
        const std::size_t count = base.size() / Values;
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < count; ++i) {
            const auto gap = [&](const std::size_t j) {
                return Distance(
                    Point3{base[Values * i], base[Values * i + 1], base[Values * i + 2]},
                    Point3{moved[Values * j], moved[Values * j + 1], moved[Values * j + 2]});
            };
            std::size_t nearest = 0;
            for (std::size_t j = 1; j < count; ++j) {
                nearest = gap(j) < gap(nearest) ? j : nearest;
            }
            for (std::size_t k = 0; k < Values; ++k) {
                const double difference = moved[Values * nearest + k] - base[Values * i + k];
                const bool angle = k >= 3;
                r[Values * i + k] =
                    (angle ? std::remainder(difference, tau) : difference) - h * jv[Values * i + k];
            }
        }
        return Norm(r);
    };
    const double coarse = residual(1e-2);
    const double fine = residual(1e-3);
    RC_ASSERT(fine <= coarse / 30.0 + 1e-10 * (1.0 + Norm(base)));

    std::vector<double> w;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
    const auto vjp =
        Stereon::IntersectVjp(conic, model.ConicBindings, Surface{surface}, model.SurfaceBindings,
                              Meet(conic, surface), model.Sensitivity);
    RC_ASSERT(vjp.has_value());
    RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
    const auto jtw = *vjp->GetValues();
    double left = 0.0;
    double right = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        left += jv[i] * w[i];
        scale += std::abs(jv[i] * w[i]);
    }
    for (std::size_t k = 0; k < seed.size(); ++k) {
        right += seed[k] * jtw[k];
        scale += std::abs(seed[k] * jtw[k]);
    }
    RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
}

/// A random orthonormal frame: X in the xy-plane at angle theta, Y tilted
/// out of it by phi.
std::pair<Vec3, Vec3> RandomFrame()
{
    const double theta = *rc::gen::inRange(0, 64) * tau / 64.0;
    const double phi = *rc::gen::inRange(-16, 17) * tau / 64.0;
    return {Vec3{std::cos(theta), std::sin(theta), 0.0},
            Vec3{-std::sin(theta) * std::cos(phi), std::cos(theta) * std::cos(phi), std::sin(phi)}};
}

/// A circle or an ellipse near the origin, its radii in [1, 3].
template <typename C>
C RandomConic()
{
    const auto [x, y] = RandomFrame();
    const Point3 center{Coordinate() / 4.0, Coordinate() / 4.0, Coordinate() / 4.0};
    if constexpr (std::is_same_v<C, Circle>) {
        return Circle{center, *rc::gen::inRange(10, 31) / 10.0, x, y};
    } else {
        const double a = *rc::gen::inRange(15, 31) / 10.0;
        return Ellipse{center, a, a * *rc::gen::inRange(4, 10) / 10.0, x, y};
    }
}

/// A sphere, a cylinder about z and a plane through the origin, each
/// crossing a conic near the origin.
Sphere RandomSphere()
{
    return Sphere{Point3{Coordinate() / 2.0, Coordinate() / 2.0, Coordinate() / 2.0},
                  *rc::gen::inRange(10, 31) / 10.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Cylinder RandomCylinder()
{
    return Cylinder{Point3{Coordinate() / 2.0, Coordinate() / 2.0, 0.0},
                    *rc::gen::inRange(10, 31) / 10.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Plane RandomPlane()
{
    const auto [x, y] = RandomFrame();
    return Plane{Point3{Coordinate() / 4.0, Coordinate() / 4.0, Coordinate() / 4.0}, x, y};
}

template <typename C>
void RequireConicProperties(const bool turning = false)
{
    const C conic = RandomConic<C>();
    switch (*rc::gen::inRange(0, 3)) {
    case 0: {
        const Sphere sphere = RandomSphere();
        RC_PRE(ConicWellConditioned(conic, sphere));
        RequireConicProducts(conic, sphere, turning);
        break;
    }
    case 1: {
        const Cylinder cylinder = RandomCylinder();
        RC_PRE(ConicWellConditioned(conic, cylinder));
        RequireConicProducts(conic, cylinder, turning);
        break;
    }
    default: {
        const Plane plane = RandomPlane();
        RC_PRE(ConicWellConditioned(conic, plane));
        RequireConicProducts(conic, plane, turning);
        break;
    }
    }
}

/// The circle of radius 2 about the origin in the xz-plane, and the
/// ellipse of semi-axes 3 and 2 in it; the plane z = 1 cuts both.
const Circle UpCircle{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
const Ellipse UpEllipse{Point3{}, 3.0, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
const Plane Shelf{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

/// The forward product for one binding moved at unit rate.
template <typename C>
std::vector<double> RateOf(const C& conic, const ConicComponent component)
{
    ParameterSet set;
    const ParameterId p =
        *set.Add("p",
                 component >= ConicComponent::RotationX ? ParameterQuantity::Angle
                                                        : ParameterQuantity::Length,
                 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<ConicBinding> bindings{{p, component}};
    const auto at = Stereon::Intersect(Stereon::Curve{conic}, Surface{Shelf}, ctx);
    REQUIRE(at.has_value());
    const auto jvp = Stereon::IntersectJvp(conic, bindings, Surface{Shelf}, {}, *at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

} // namespace

TEST_CASE("a circle's and an ellipse's crossings of a plane move as their closed forms say",
          "[intersect][linearization][conic]")
{
    const double root3 = std::sqrt(3.0);
    SECTION("a circle: its radius and its centre's height")
    {
        // sin t = 1 / 2: t = pi / 6 at x = sqrt 3, then 5 pi / 6 at -sqrt 3.
        // A growing radius moves x by r / x and t by -(h / r^2) / cos t; a
        // rising centre moves x by (h - c) / x and t by -1 / (r cos t).
        const std::vector<double> radius = RateOf(UpCircle, ConicComponent::Radius);
        const std::vector<double> rising = RateOf(UpCircle, ConicComponent::CenterZ);
        REQUIRE(radius.size() == 2 * Values);
        for (const auto& [k, x, cosine] :
             {std::tuple{0, root3, root3 / 2.0}, std::tuple{1, -root3, -root3 / 2.0}}) {
            const std::size_t at = Values * static_cast<std::size_t>(k);
            CHECK(std::abs(radius[at] - 2.0 / x) <= 1e-14);
            CHECK(std::abs(radius[at + 2]) <= 1e-14);
            CHECK(std::abs(radius[at + 3] + 0.25 / cosine) <= 1e-14);
            CHECK(std::abs(radius[at + 4] - 2.0 / x) <= 1e-14);
            CHECK(std::abs(rising[at] - 1.0 / x) <= 1e-14);
            CHECK(std::abs(rising[at + 2]) <= 1e-14);
            CHECK(std::abs(rising[at + 3] + 0.5 / cosine) <= 1e-14);
        }
    }
    SECTION("an ellipse: its major and minor radii")
    {
        // 2 sin t = 1: t = pi / 6 at x = 3 sqrt 3 / 2. The major radius moves
        // x by cos t and leaves t; the minor one moves t by -(h / b^2) /
        // cos t = -1 / (2 sqrt 3), and x by -a sin t t' = 3 / (4 sqrt 3).
        const std::vector<double> major = RateOf(UpEllipse, ConicComponent::MajorRadius);
        const std::vector<double> minor = RateOf(UpEllipse, ConicComponent::MinorRadius);
        REQUIRE(major.size() == 2 * Values);
        CHECK(std::abs(major[0] - root3 / 2.0) <= 1e-14);
        CHECK(std::abs(major[3]) <= 1e-14);
        CHECK(std::abs(minor[0] - 3.0 / (4.0 * root3)) <= 1e-14);
        CHECK(std::abs(minor[2]) <= 1e-14);
        CHECK(std::abs(minor[3] + 1.0 / (2.0 * root3)) <= 1e-14);
    }
}

TEST_CASE("circle and ellipse crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][conic][property]")
{
    rc::prop("circles through spheres, cylinders and planes",
             [] { RequireConicProperties<Circle>(); });
    rc::prop("ellipses through spheres, cylinders and planes",
             [] { RequireConicProperties<Ellipse>(); });
}

TEST_CASE("a conic's bindings are refused where it has no such radius",
          "[intersect][linearization][conic]")
{
    ParameterSet set;
    const ParameterId length = *set.Add("p", ParameterQuantity::Length, 1.0);
    const ParameterId angle = *set.Add("a", ParameterQuantity::Angle, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    const auto circleAt = Stereon::Intersect(Stereon::Curve{UpCircle}, Surface{Shelf}, ctx);
    const auto ellipseAt = Stereon::Intersect(Stereon::Curve{UpEllipse}, Surface{Shelf}, ctx);
    REQUIRE(circleAt.has_value());
    REQUIRE(ellipseAt.has_value());
    for (const ConicComponent component :
         {ConicComponent::MajorRadius, ConicComponent::MinorRadius}) {
        const std::vector<ConicBinding> bad{{length, component}};
        REQUIRE(code(Stereon::IntersectJvp(UpCircle, bad, Surface{Shelf}, {}, *circleAt,
                                           sensitivity)) == ErrorCode::InvalidInput);
    }
    const std::vector<ConicBinding> radius{{length, ConicComponent::Radius}};
    REQUIRE(code(Stereon::IntersectJvp(UpEllipse, radius, Surface{Shelf}, {}, *ellipseAt,
                                       sensitivity)) == ErrorCode::InvalidInput);
    // A centre is moved by a length, never an angle.
    const std::vector<ConicBinding> turned{{angle, ConicComponent::CenterX}};
    REQUIRE(code(Stereon::IntersectJvp(UpCircle, turned, Surface{Shelf}, {}, *circleAt,
                                       sensitivity)) == ErrorCode::InvalidInput);
}

// --- Turning a frame ----------------------------------------------------------

namespace {

/// The forward product of one binding of the surface turned at unit rate.
template <typename S>
std::vector<double> TurnRate(const Line& line, const S& surface, const SurfaceComponent turn)
{
    ParameterSet set;
    const ParameterId w = *set.Add("w", ParameterQuantity::Angle, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<SurfaceBinding> bindings{{w, turn}};
    const auto at = Stereon::Intersect(Stereon::Curve{line}, Surface{surface}, ctx);
    REQUIRE(at.has_value());
    const auto jvp = Stereon::IntersectJvp(line, {}, Surface{surface}, bindings, *at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    return {jvp->GetValues()->begin(), jvp->GetValues()->end()};
}

/// The point's values at offset k against (x', t', u', v').
void RequireTurned(const std::vector<double>& v, const std::size_t k, const Vec3& x, const double t,
                   const double u, const double w)
{
    INFO("point at offset " << k);
    CHECK(std::abs(v[k] - x.X) <= 1e-14);
    CHECK(std::abs(v[k + 1] - x.Y) <= 1e-14);
    CHECK(std::abs(v[k + 2] - x.Z) <= 1e-14);
    CHECK(std::abs(v[k + 3] - t) <= 1e-14);
    CHECK(std::abs(v[k + 4] - u) <= 1e-14);
    CHECK(std::abs(v[k + 5] - w) <= 1e-14);
}

} // namespace

TEST_CASE("a turning frame moves crossings as their closed forms say",
          "[intersect][linearization][turn]")
{
    SECTION("the line y = 0, z = 2 across the cylinder of radius 1.5 about z")
    {
        // The points (-+1.5, 0, 2), t = x, u = pi then 0, v = 2. Turning the
        // axis about Y by w tilts it to (sin w, 0, cos w): the distance from
        // (x, 0, h) to it is R when x^2 + h^2 - (x sin w + h cos w)^2 = R^2,
        // so 2 x x' = 2 h x, x' = h, and v = (P - O) . axis moves by x.
        // About X, the line meets the tilted tube at second order only:
        // x' = 0, and u moves by h / x; about the axis, only u, by -1.
        const Line line{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}};
        const Cylinder tube{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const auto y = TurnRate(line, tube, SurfaceComponent::RotationY);
        REQUIRE(y.size() == 2 * Values);
        RequireTurned(y, 0, Vec3{2.0, 0.0, 0.0}, 2.0, 0.0, -1.5);
        RequireTurned(y, Values, Vec3{2.0, 0.0, 0.0}, 2.0, 0.0, 1.5);
        const auto x = TurnRate(line, tube, SurfaceComponent::RotationX);
        RequireTurned(x, 0, Vec3{}, 0.0, -2.0 / 1.5, 0.0);
        RequireTurned(x, Values, Vec3{}, 0.0, 2.0 / 1.5, 0.0);
        const auto z = TurnRate(line, tube, SurfaceComponent::RotationZ);
        RequireTurned(z, 0, Vec3{}, 0.0, -1.0, 0.0);
        RequireTurned(z, Values, Vec3{}, 0.0, -1.0, 0.0);
    }
    SECTION("the vertical line through (1, 0, 0) across the plane z = 0")
    {
        // About Y the plane's XAxis tips to (cos w, 0, -sin w): the line meets
        // it at z = -tan w, x' = (0, 0, -1), t' = -1. About X the point lies
        // on the axis and stays. About the normal the plane turns in itself:
        // no point moves, (u, v) = (x, y) moves by (y, -x).
        const Line line{Point3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
        const Plane floor{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        RequireTurned(TurnRate(line, floor, SurfaceComponent::RotationY), 0, Vec3{0.0, 0.0, -1.0},
                      -1.0, 0.0, 0.0);
        RequireTurned(TurnRate(line, floor, SurfaceComponent::RotationX), 0, Vec3{}, 0.0, 0.0, 0.0);
        RequireTurned(TurnRate(line, floor, SurfaceComponent::RotationZ), 0, Vec3{}, 0.0, 0.0,
                      -1.0);
    }
    SECTION("a sphere turned about its centre moves no point")
    {
        // The x-axis meets the sphere of radius 2 at (-+2, 0, 0) on its
        // equator. A turn about z moves u by -1; about y, v by -+1; about x,
        // nothing.
        const Sphere ball{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const Line line{Point3{}, Vec3{1.0, 0.0, 0.0}};
        RequireTurned(TurnRate(line, ball, SurfaceComponent::RotationZ), 0, Vec3{}, 0.0, -1.0, 0.0);
        const auto y = TurnRate(line, ball, SurfaceComponent::RotationY);
        RequireTurned(y, 0, Vec3{}, 0.0, 0.0, -1.0);
        RequireTurned(y, Values, Vec3{}, 0.0, 0.0, 1.0);
        RequireTurned(TurnRate(line, ball, SurfaceComponent::RotationX), Values, Vec3{}, 0.0, 0.0,
                      0.0);
    }
    SECTION("the circle of radius 2 in the xz-plane turning, cut by z = 1")
    {
        // The points (+-sqrt 3, 0, 1) at t = pi / 6 and 5 pi / 6; the plane's
        // (u, v) are their x and y. About X the circle tips out of y = 0:
        // x' = X x (C - O) = (0, -1, 0), v' = -1. About Y = z, x' =
        // (0, +-sqrt 3, 0). About its normal X x Y = -y it turns in itself:
        // no point moves and t' = -1.
        const double root3 = std::sqrt(3.0);
        const auto x = RateOf(UpCircle, ConicComponent::RotationX);
        const auto y = RateOf(UpCircle, ConicComponent::RotationY);
        const auto z = RateOf(UpCircle, ConicComponent::RotationZ);
        REQUIRE(x.size() == 2 * Values);
        for (const auto& [k, sign] : {std::pair{std::size_t{0}, 1.0}, std::pair{Values, -1.0}}) {
            RequireTurned(x, k, Vec3{0.0, -1.0, 0.0}, 0.0, 0.0, -1.0);
            RequireTurned(y, k, Vec3{0.0, sign * root3, 0.0}, 0.0, 0.0, sign * root3);
            RequireTurned(z, k, Vec3{}, -1.0, 0.0, 0.0);
        }
    }
}

TEST_CASE("turning frames: second-order Taylor decay and the transpose",
          "[intersect][linearization][turn][property]")
{
    rc::prop("for random lines and planes", [] {
        const auto [line, plane] = RandomPlaneCrossing();
        RequireTaylor(line, plane, true);
        RequireTranspose(line, plane, true);
    });
    rc::prop("for random lines through spheres and cylinders", [] {
        const auto [line, surfaces] = RandomQuadricCrossing();
        RequireTaylor(line, surfaces.first, true);
        RequireTaylor(line, surfaces.second, true);
        RequireTranspose(line, surfaces.first, true);
        RequireTranspose(line, surfaces.second, true);
    });
    rc::prop("for random lines through ellipsoids", [] {
        const auto [line, egg] = RandomEllipsoidCrossing();
        RequireTaylor(line, egg, true);
        RequireTranspose(line, egg, true);
    });
    const auto frameAbout = [](const Point3& center) {
        const Vec3 axis{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(axis) >= 0.5);
        return IntersectTest::CylinderAbout(center, axis, 1.0);
    };
    rc::prop("for random cubic B-splines through cones", [&] {
        const Point3 center{Coordinate(), Coordinate(), Coordinate()};
        const Cylinder frame = frameAbout(center);
        const double radius = *rc::gen::inRange(8, 25) / 8.0;
        const double sign = *rc::gen::arbitrary<bool>() ? 1.0 : -1.0;
        const double angle = sign * *rc::gen::inRange(2, 11) / 10.0;
        const Cone cone{center, radius, angle, frame.XAxis, frame.YAxis};
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(direction) >= radius);
        const auto [curve, poles] = RandomSplineThrough(center, direction, radius);
        RequireSplineProducts(curve, poles, cone, true);
    });
    rc::prop("for random cubic B-splines through tori", [&] {
        const Point3 center{Coordinate(), Coordinate(), Coordinate()};
        const Cylinder frame = frameAbout(center);
        const double minor = *rc::gen::inRange(4, 9) / 8.0;
        const double major = minor + *rc::gen::inRange(4, 13) / 8.0;
        const Torus torus{center, major, minor, frame.XAxis, frame.YAxis};
        const double u = *rc::gen::inRange(0, 16) / 16.0 * tau;
        const Point3 tube =
            center + major * (std::cos(u) * frame.XAxis + std::sin(u) * frame.YAxis);
        const Vec3 direction{Coordinate(), Coordinate(), Coordinate()};
        RC_PRE(Length(direction) >= major + minor);
        const auto [curve, poles] = RandomSplineThrough(tube, direction, minor);
        RequireSplineProducts(curve, poles, torus, true);
    });
    rc::prop("circles turning through turning spheres, cylinders and planes",
             [] { RequireConicProperties<Circle>(true); });
    rc::prop("ellipses turning through turning spheres, cylinders and planes",
             [] { RequireConicProperties<Ellipse>(true); });
}

TEST_CASE("a turn is refused where there is no frame or as a length",
          "[intersect][linearization][turn]")
{
    ParameterSet set;
    const ParameterId length = *set.Add("p", ParameterQuantity::Length, 1.0);
    const ParameterId angle = *set.Add("a", ParameterQuantity::Angle, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    const auto code = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        return result.error().GetCode();
    };
    // A turn is an angle.
    const Line line{Point3{1.0, 0.0, -1.0}, Vec3{0.0, 0.0, 1.0}};
    const Plane floor{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const auto onFloor = Stereon::Intersect(Stereon::Curve{line}, Surface{floor}, ctx);
    REQUIRE(onFloor.has_value());
    const std::vector<SurfaceBinding> asLength{{length, SurfaceComponent::RotationY}};
    CHECK(code(Stereon::IntersectJvp(line, {}, Surface{floor}, asLength, *onFloor, sensitivity)) ==
          ErrorCode::InvalidInput);
    const std::vector<ConicBinding> conicAsLength{{length, ConicComponent::RotationZ}};
    const auto circleAt = Stereon::Intersect(Stereon::Curve{UpCircle}, Surface{Shelf}, ctx);
    REQUIRE(circleAt.has_value());
    CHECK(code(Stereon::IntersectJvp(UpCircle, conicAsLength, Surface{Shelf}, {}, *circleAt,
                                     sensitivity)) == ErrorCode::InvalidInput);
    // A B-spline surface has poles, no frame.
    const BSplineSurface patch = FlatPatch(false);
    const auto onPatch = Stereon::Intersect(
        Stereon::Curve{Line{Point3{0.5, 0.5, -1.0}, Vec3{0.0, 0.0, 1.0}}}, Surface{patch}, ctx);
    REQUIRE(onPatch.has_value());
    const std::vector<SurfaceBinding> turned{{angle, SurfaceComponent::RotationX}};
    CHECK(code(Stereon::IntersectJvp(Line{Point3{0.5, 0.5, -1.0}, Vec3{0.0, 0.0, 1.0}}, {},
                                     Surface{patch}, turned, *onPatch, sensitivity)) ==
          ErrorCode::InvalidInput);
    // A turn of a tangency is still no derivative: the sphere of radius 1
    // about (0, 0, 1) touches the floor at the origin.
    const Sphere resting{Point3{0.0, 0.0, 1.0}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Line along{Point3{-2.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    const auto touch = Stereon::Intersect(Stereon::Curve{along}, Surface{resting}, ctx);
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1u);
    const std::vector<SurfaceBinding> spin{{angle, SurfaceComponent::RotationY}};
    const auto tangent =
        Stereon::IntersectJvp(along, {}, Surface{resting}, spin, *touch, sensitivity);
    REQUIRE(tangent.has_value());
    CHECK(tangent->GetStatus() == DerivativeStatus::Nondifferentiable);
}
