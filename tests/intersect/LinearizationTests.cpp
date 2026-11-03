// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The first intersection derivative (ADR-0021, sections 6 and 9; the v0.2
// gate): a line crossing a plane, differentiated through F = C(t) - S(u, v)
// = 0. The forward product agrees with the closed form of the crossing; on
// random geometry its Taylor residual |y(p + h v) - y(p) - h J v| decays
// with the second power of h; the reverse product is its transpose; and
// every point that has no derivative says so: an overlap or a touch is
// Nondifferentiable, an undecided point Indeterminate, a grazing crossing
// NearEvent, and none of them carries a usable value unless qualified.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::DerivativeStatus;
using Stereon::ErrorCode;
using Stereon::Line;
using Stereon::LineBinding;
using Stereon::LineComponent;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Plane;
using Stereon::PlaneBinding;
using Stereon::PlaneComponent;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Vec3;

namespace {

constexpr std::size_t Values = Stereon::CurveSurfacePointValues;

/// Nine parameters, each bound to one coordinate: the line's origin and
/// direction, then the plane's origin, with the given references.
struct Model {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<LineBinding> LineBindings;
    std::vector<PlaneBinding> PlaneBindings;
};

Model FullModel(const std::vector<double>& references)
{
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < 9; ++k) {
        ids.push_back(*set.Add("p" + std::to_string(k), ParameterQuantity::Length, references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < 6; ++k) {
        model.LineBindings.push_back({ids[k], static_cast<LineComponent>(k)});
    }
    for (std::size_t k = 0; k < 3; ++k) {
        model.PlaneBindings.push_back({ids[6 + k], static_cast<PlaneComponent>(k)});
    }
    return model;
}

/// The geometry moved by `step` times the seed, each parameter by its
/// reference per unit.
std::pair<Line, Plane> Moved(Line line, Plane plane, const std::vector<double>& seed,
                             const std::vector<double>& references, const double step)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    line.Origin = line.Origin + Vec3{d(0), d(1), d(2)};
    line.Direction = line.Direction + Vec3{d(3), d(4), d(5)};
    plane.Origin = plane.Origin + Vec3{d(6), d(7), d(8)};
    return {line, plane};
}

/// The one crossing as the product's values are laid out: x, then t, u, v.
std::vector<double> Crossing(const Line& line, const Plane& plane)
{
    const auto result = Intersect(Stereon::Curve{line}, Stereon::Surface{plane}, ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 1);
    const CurveSurfacePoint& p = result->Points.front();
    return {p.Position.X, p.Position.Y, p.Position.Z, p.Parameter, p.OnSurface.X, p.OnSurface.Y};
}

CurveSurfaceIntersection Meet(const Line& line, const Plane& plane)
{
    return *Intersect(Stereon::Curve{line}, Stereon::Surface{plane}, ctx);
}

std::vector<double> Jvp(Model& model, const Line& line, const Plane& plane,
                        const std::vector<double>& seed)
{
    REQUIRE(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp = Stereon::IntersectJvp(line, model.LineBindings, plane, model.PlaneBindings,
                                           Meet(line, plane), model.Sensitivity);
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
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

TEST_CASE("the Taylor residual of the crossing decays with the square of the step",
          "[intersect][linearization][property]")
{
    rc::prop("for random lines and planes, seeds and references", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Vec3 normal{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(normal != Vec3{} && direction != Vec3{});
        // Well away from grazing, so that h = 1e-3 stays in the regime.
        RC_PRE(std::abs(Dot(normal, direction)) >= 0.25 * Length(normal) * Length(direction));
        const Plane plane = IntersectTest::PlaneWithNormal(
            Point3{coordinate(), coordinate(), coordinate()}, normal);
        const Line line{Point3{coordinate(), coordinate(), coordinate()}, direction};
        std::vector<double> references;
        std::vector<double> seed;
        for (int k = 0; k < 9; ++k) {
            references.push_back(*rc::gen::inRange(1, 20) / 10.0);
            seed.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        }
        Model model = FullModel(references);
        const std::vector<double> jv = Jvp(model, line, plane, seed);
        const std::vector<double> base = Crossing(line, plane);
        const auto residual = [&](const double h) {
            const auto [l, p] = Moved(line, plane, seed, references, h);
            const std::vector<double> moved = Crossing(l, p);
            std::vector<double> r(Values);
            for (std::size_t i = 0; i < Values; ++i) {
                r[i] = moved[i] - base[i] - h * jv[i];
            }
            return Norm(r);
        };
        const double coarse = residual(1e-2);
        const double fine = residual(1e-3);
        // Second order: a tenth of the step leaves about a hundredth of the
        // residual; rounding sets a floor well below both.
        const double floor = 1e-10 * (1.0 + Norm(base));
        RC_ASSERT(fine <= coarse / 30.0 + floor);
    });
}

TEST_CASE("the reverse product is the forward one's transpose", "[intersect][linearization]")
{
    rc::prop("<J v, w> = <v, J^T w>", [] {
        const auto coordinate = [] { return *rc::gen::inRange(-40, 41) / 8.0; };
        const Vec3 normal{coordinate(), coordinate(), coordinate()};
        const Vec3 direction{coordinate(), coordinate(), coordinate()};
        RC_PRE(normal != Vec3{} && direction != Vec3{});
        RC_PRE(std::abs(Dot(normal, direction)) >= 0.25 * Length(normal) * Length(direction));
        const Plane plane = IntersectTest::PlaneWithNormal(
            Point3{coordinate(), coordinate(), coordinate()}, normal);
        const Line line{Point3{coordinate(), coordinate(), coordinate()}, direction};
        std::vector<double> references;
        std::vector<double> v;
        std::vector<double> w;
        for (std::size_t k = 0; k < 9; ++k) {
            references.push_back(*rc::gen::inRange(1, 20) / 10.0);
            v.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        }
        for (std::size_t k = 0; k < Values; ++k) {
            w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
        }
        Model model = FullModel(references);
        const std::vector<double> jv = Jvp(model, line, plane, v);
        RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
        const auto vjp = Stereon::IntersectVjp(line, model.LineBindings, plane, model.PlaneBindings,
                                               Meet(line, plane), model.Sensitivity);
        RC_ASSERT(vjp.has_value());
        RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
        const auto jtw = *vjp->GetValues();
        double left = 0.0;
        double right = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < Values; ++i) {
            left += jv[i] * w[i];
            scale += std::abs(jv[i] * w[i]);
        }
        for (std::size_t k = 0; k < 9; ++k) {
            right += v[k] * jtw[k];
            scale += std::abs(v[k] * jtw[k]);
        }
        RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
    });
}

TEST_CASE("a point without a derivative says so and carries none", "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    Model model = FullModel(std::vector<double>(9, 1.0));
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(9, 1.0)).has_value());
    const auto jvp = [&](const Line& line, const CurveSurfaceIntersection& at) {
        return Stereon::IntersectJvp(line, model.LineBindings, ground, model.PlaneBindings, at,
                                     model.Sensitivity);
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

TEST_CASE("Linearize exposes the residual, the system and the regularity",
          "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // At 30 degrees to the plane: the conditioning is the sine, 1/2.
    const Line slant{Point3{0.0, 0.0, -1.0}, Vec3{std::numbers::sqrt3, 0.0, 1.0}};
    const CurveSurfaceIntersection at = Meet(slant, ground);
    const auto lin = Stereon::Linearize(slant, ground, at.Points.front());
    REQUIRE(lin.has_value());
    REQUIRE(Length(lin->Residual) <= 1e-15);
    REQUIRE(std::abs(lin->Conditioning - 0.5) <= 1e-15);
    REQUIRE(lin->Status == DerivativeStatus::Valid);
    REQUIRE(lin->Tangent == slant.Direction);
    REQUIRE(lin->Su == ground.XAxis);
    REQUIRE(lin->Sv == ground.YAxis);
    CurveSurfacePoint broken = at.Points.front();
    broken.Parameter = std::nan("");
    REQUIRE(Stereon::Linearize(slant, ground, broken).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("invalid sensitivity input is refused", "[intersect][linearization]")
{
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Line up{Point3{0.0, 0.0, -1.0}, Vec3{0.0, 0.0, 1.0}};
    const CurveSurfaceIntersection at = Meet(up, ground);
    Model model = FullModel(std::vector<double>(9, 1.0));
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    // No seed yet; no adjoint, or one of the wrong length.
    REQUIRE(code(Stereon::IntersectJvp(up, model.LineBindings, ground, model.PlaneBindings, at,
                                       model.Sensitivity)) == ErrorCode::InvalidInput);
    REQUIRE(model.Sensitivity.SetAdjoint(std::vector<double>(5, 0.0)).has_value());
    REQUIRE(code(Stereon::IntersectVjp(up, model.LineBindings, ground, model.PlaneBindings, at,
                                       model.Sensitivity)) == ErrorCode::InvalidInput);
    // A parameter of another set, and an angle moving a coordinate.
    ParameterSet other;
    const ParameterId foreign = *other.Add("q", ParameterQuantity::Length);
    const ParameterId angle = *other.Add("a", ParameterQuantity::Angle);
    REQUIRE(model.Sensitivity.SetSeed(std::vector<double>(9, 1.0)).has_value());
    for (const ParameterId id : {foreign, angle}) {
        const std::vector<LineBinding> bad{{id, LineComponent::OriginX}};
        REQUIRE(code(Stereon::IntersectJvp(up, bad, ground, model.PlaneBindings, at,
                                           model.Sensitivity)) == ErrorCode::InvalidInput);
    }
}

TEST_CASE("Intersect declares its forward and reverse products", "[intersect][linearization]")
{
    const auto* capability = Stereon::FindCapability(Stereon::IntersectCapabilities(), "Intersect");
    REQUIRE(capability != nullptr);
    REQUIRE(capability->Supports(Stereon::DerivativeProduct::Jvp));
    REQUIRE(capability->Supports(Stereon::DerivativeProduct::Vjp));
    REQUIRE_FALSE(capability->Supports(Stereon::DerivativeProduct::DenseJacobian));
}
