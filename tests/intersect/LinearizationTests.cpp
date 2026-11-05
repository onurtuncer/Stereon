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
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::DerivativeStatus;
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
template <typename S>
constexpr std::size_t SurfaceParameters = std::is_same_v<S, Plane> ? 3 : 4;

/// Parameters each bound to one coordinate: the line's origin and
/// direction, then the surface's origin, then its radius where it has one.
struct Model {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<LineBinding> LineBindings;
    std::vector<SurfaceBinding> SurfaceBindings;
};

Model FullModel(const std::vector<double>& references)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        ids.push_back(*set.Add("p" + std::to_string(k), ParameterQuantity::Length, references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < 6; ++k) {
        model.LineBindings.push_back({ids[k], static_cast<LineComponent>(k)});
    }
    for (std::size_t k = 6; k < references.size(); ++k) {
        model.SurfaceBindings.push_back({ids[k], static_cast<SurfaceComponent>(k - 6)});
    }
    return model;
}

/// The geometry moved by `step` times the seed, each parameter by its
/// reference per unit.
template <typename S>
std::pair<Line, S> Moved(Line line, S surface, const std::vector<double>& seed,
                         const std::vector<double>& references, const double step)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    line.Origin = line.Origin + Vec3{d(0), d(1), d(2)};
    line.Direction = line.Direction + Vec3{d(3), d(4), d(5)};
    OriginOf(surface) = OriginOf(surface) + Vec3{d(6), d(7), d(8)};
    if constexpr (SurfaceParameters<S> == 4) {
        surface.Radius += d(9);
    }
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
void RequireTaylor(const Line& line, const S& surface)
{
    const auto [references, seed] = RandomSeed(6 + SurfaceParameters<S>);
    Model model = FullModel(references);
    const std::vector<double> jv = Jvp(model, line, surface, seed);
    const std::vector<double> base = Points(line, surface);
    const auto residual = [&](const double h) {
        const auto [l, s] = Moved(line, surface, seed, references, h);
        const std::vector<double> moved = Points(l, s);
        RC_ASSERT(moved.size() == base.size());
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < base.size(); ++i) {
            const double difference = moved[i] - base[i];
            r[i] = (i % Values == 4 ? std::remainder(difference, tau) : difference) - h * jv[i];
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
void RequireTranspose(const Line& line, const S& surface)
{
    const auto [references, v] = RandomSeed(6 + SurfaceParameters<S>);
    Model model = FullModel(references);
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
    // A torus has no bindings yet: a NotSupported status, never a value.
    const Surface torus(
        Stereon::Torus{Point3{}, 2.0, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});
    const auto unsupported = Stereon::IntersectJvp(up, model.LineBindings, torus, {},
                                                   CurveSurfaceIntersection{}, model.Sensitivity);
    REQUIRE(unsupported.has_value());
    REQUIRE(unsupported->GetStatus() == DerivativeStatus::NotSupported);
    REQUIRE_FALSE(unsupported->GetValues().has_value());
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

SplineModel MakeSplineModel(const std::vector<PoleBinding>& poles, const std::size_t surface,
                            const std::vector<double>& references)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const bool weight = k < poles.size() && poles[k].Component == PoleComponent::Weight;
        ids.push_back(*set.Add(
            "p" + std::to_string(k),
            weight ? ParameterQuantity::Dimensionless : ParameterQuantity::Length, references[k]));
    }
    SplineModel model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < poles.size(); ++k) {
        model.PoleBindings.push_back({ids[k], poles[k].Pole, poles[k].Component});
    }
    for (std::size_t k = 0; k < surface; ++k) {
        model.SurfaceBindings.push_back({ids[poles.size() + k], static_cast<SurfaceComponent>(k)});
    }
    return model;
}

/// The curve and surface moved by `step` times the seed.
template <typename S>
std::pair<BSplineCurve, S> MovedSpline(BSplineCurve curve, S surface, const SplineModel& model,
                                       const std::vector<double>& seed,
                                       const std::vector<double>& references, const double step)
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
    const std::size_t s = model.PoleBindings.size();
    OriginOf(surface) = OriginOf(surface) + Vec3{d(s), d(s + 1), d(s + 2)};
    if constexpr (SurfaceParameters<S> == 4) {
        surface.Radius += d(s + 3);
    }
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
    RC_ASSERT(jvp.has_value());
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
    const CurveSurfaceIntersection at = Meet(curve, surface);
    if (at.Points.empty() || !at.Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(at.Points, [&](const CurveSurfacePoint& p) {
        const auto lin = Stereon::Linearize(curve, Surface{surface}, p);
        const bool offPole = !std::is_same_v<S, Sphere> || std::abs(p.OnSurface.Y) <= 1.2;
        return p.Decided && offPole && lin && lin->Conditioning >= 0.25;
    });
}

/// Taylor decay and the transpose identity for a B-spline against a surface.
template <typename S>
void RequireSplineProducts(const BSplineCurve& curve, const std::vector<PoleBinding>& poles,
                           const S& surface)
{
    RC_PRE(SplineWellConditioned(curve, surface));
    const std::size_t count = poles.size() + SurfaceParameters<S>;
    const auto [references, seed] = RandomSeed(count);
    SplineModel model = MakeSplineModel(poles, SurfaceParameters<S>, references);
    const std::vector<double> jv = SplineJvp(model, curve, surface, seed);
    const std::vector<double> base = Points(curve, surface);
    const auto residual = [&](const double h) {
        const auto [c, s] = MovedSpline(curve, surface, model, seed, references, h);
        const std::vector<double> moved = Points(c, s);
        RC_ASSERT(moved.size() == base.size());
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < base.size(); ++i) {
            const double difference = moved[i] - base[i];
            r[i] = (i % Values == 4 ? std::remainder(difference, tau) : difference) - h * jv[i];
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
