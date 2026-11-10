// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Derivatives of surface–surface branch points (ADR-0021, section 9; the
// v0.2 gate), with the normal-plane tracking constraint: a sample of a
// branch moves within the plane through it normal to the branch. The
// forward product agrees with the closed form of a plane's circle on a
// sphere; on random planes, spheres and cylinders, exact and traced, its
// Taylor residual against the moved point in the same tracking plane
// decays with the second power of the step; the reverse product is its transpose; and a
// sample without a derivative says so: a tangential branch or an overlap
// is Nondifferentiable, an uncertified branch Indeterminate.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::BranchSample;
using Stereon::Cylinder;
using Stereon::DerivativeStatus;
using Stereon::Ellipsoid;
using Stereon::ErrorCode;
using Stereon::IntersectionBranch;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SensitivityContext;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceBinding;
using Stereon::SurfaceComponent;
using Stereon::SurfaceSurfaceIntersection;
using Stereon::Vec3;

namespace {

constexpr std::size_t Values = Stereon::BranchPointValues;
constexpr double tau = 2.0 * std::numbers::pi;

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
Point3& OriginOf(Ellipsoid& s)
{
    return s.Center;
}

/// The origin, then the radius where the surface has one, or an
/// ellipsoid's three semi-axes.
template <typename S>
std::vector<SurfaceComponent> ComponentsOf()
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
    return components;
}

/// The surface moved by d(first + k) along its k-th component.
template <typename S, typename Step>
void MoveSurface(S& surface, const Step& d, const std::size_t first)
{
    OriginOf(surface) = OriginOf(surface) + Vec3{d(first), d(first + 1), d(first + 2)};
    if constexpr (std::is_same_v<S, Ellipsoid>) {
        surface.A += d(first + 3);
        surface.B += d(first + 4);
        surface.C += d(first + 5);
    } else if constexpr (!std::is_same_v<S, Plane>) {
        surface.Radius += d(first + 3);
    }
}

/// One Length parameter per component of A, then of B.
struct Model {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<SurfaceBinding> BindingsA;
    std::vector<SurfaceBinding> BindingsB;
};

template <typename SA, typename SB>
Model ModelFor(const std::vector<double>& references)
{
    const std::vector<SurfaceComponent> a = ComponentsOf<SA>();
    const std::vector<SurfaceComponent> b = ComponentsOf<SB>();
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        ids.push_back(*set.Add("p" + std::to_string(k), ParameterQuantity::Length, references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < a.size(); ++k) {
        model.BindingsA.push_back({ids[k], a[k]});
    }
    for (std::size_t k = 0; k < b.size(); ++k) {
        model.BindingsB.push_back({ids[a.size() + k], b[k]});
    }
    return model;
}

template <typename SA, typename SB>
std::size_t ParameterCount()
{
    return ComponentsOf<SA>().size() + ComponentsOf<SB>().size();
}

SurfaceSurfaceIntersection Meet(const Surface& a, const Surface& b)
{
    auto result = Intersect(a, b, ctx);
    if (!result) {
        RC_FAIL(IntersectTest::ErrorText(result));
    }
    return *std::move(result);
}

/// Samples at the quarter points of every branch's range.
std::vector<BranchSample> SamplesOf(const SurfaceSurfaceIntersection& at)
{
    std::vector<BranchSample> samples;
    for (std::size_t k = 0; k < at.Branches.size(); ++k) {
        const Stereon::Interval range = at.Branches[k].Range;
        for (const double f : {0.2, 0.45, 0.7}) {
            samples.push_back({k, range.GetLower() + f * (range.GetUpper() - range.GetLower())});
        }
    }
    return samples;
}

double Norm(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

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

Stereon::Box2 Box2Of(const double lo, const double hi)
{
    return Stereon::Box2{Stereon::Interval{lo, hi}, Stereon::Interval{lo, hi}};
}

double Coordinate()
{
    return *rc::gen::inRange(-40, 41) / 8.0;
}

/// A pair of surfaces of the given types and their samples, for a check.
template <typename SA, typename SB>
struct Case {
    SA A;
    SB B;
    SurfaceSurfaceIntersection At;
    std::vector<BranchSample> Samples;
};

template <typename SA, typename SB>
std::vector<double> Jvp(Model& model, const Case<SA, SB>& c, const std::vector<double>& seed)
{
    RC_ASSERT(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp = Stereon::IntersectJvp(Surface{c.A}, model.BindingsA, Surface{c.B},
                                           model.BindingsB, c.At, c.Samples, model.Sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = jvp->GetValues();
    return {values->begin(), values->end()};
}

/// The samples' refined points and parameters, laid out as the product's
/// values: x, then u_A, v_A, u_B, v_B.
template <typename SA, typename SB>
std::vector<double> Base(const Case<SA, SB>& c)
{
    std::vector<double> values;
    for (const BranchSample& sample : c.Samples) {
        const auto lin = Stereon::Linearize(Surface{c.A}, Surface{c.B},
                                            c.At.Branches[sample.Branch], sample.Parameter);
        RC_ASSERT(lin.has_value());
        values.insert(values.end(), {lin->Position.X, lin->Position.Y, lin->Position.Z, lin->OnA.X,
                                     lin->OnA.Y, lin->OnB.X, lin->OnB.Y});
    }
    return values;
}

/// The same for the moved surfaces: the point in each sample's tracking
/// plane where they meet, which Newton's method reaches from the sample
/// (Linearize on the moved surfaces, the branch's bound widened to admit
/// the step), held to the moved equations here: on both moved surfaces and
/// in the unmoved tracking plane, to rounding.
template <typename SA, typename SB>
std::vector<double> Moved(const Case<SA, SB>& c, const SA& a, const SB& b)
{
    std::vector<double> values;
    for (const BranchSample& sample : c.Samples) {
        IntersectionBranch branch = c.At.Branches[sample.Branch];
        const Point3 x0 = Eval(branch.Curve3d, sample.Parameter);
        branch.Bound = 1.0;
        const auto lin = Stereon::Linearize(Surface{a}, Surface{b}, branch, sample.Parameter);
        RC_ASSERT(lin.has_value());
        RC_ASSERT(Distance(Eval(a, lin->OnA), Eval(b, lin->OnB)) <= 1e-13);
        RC_ASSERT(std::abs(Dot(lin->Tangent, lin->Position - x0)) <= 1e-13);
        values.insert(values.end(), {lin->Position.X, lin->Position.Y, lin->Position.Z, lin->OnA.X,
                                     lin->OnA.Y, lin->OnB.X, lin->OnB.Y});
    }
    return values;
}

/// The Taylor residual at steps 1e-2 and 1e-3 falls by second order;
/// parameters are compared modulo a turn, which a step never approaches
/// otherwise.
template <typename SA, typename SB>
void RequireTaylor(const Case<SA, SB>& c)
{
    const auto [references, seed] = RandomSeed(ParameterCount<SA, SB>());
    Model model = ModelFor<SA, SB>(references);
    const std::vector<double> jv = Jvp(model, c, seed);
    const std::vector<double> base = Base(c);
    RC_ASSERT(jv.size() == base.size());
    const auto residual = [&](const double h) {
        const auto d = [&](const std::size_t k) { return h * seed[k] * references[k]; };
        SA a = c.A;
        SB b = c.B;
        MoveSurface(a, d, 0);
        MoveSurface(b, d, ComponentsOf<SA>().size());
        const std::vector<double> moved = Moved(c, a, b);
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < base.size(); ++i) {
            const double difference = moved[i] - base[i];
            r[i] = (i % Values >= 3 ? std::remainder(difference, tau) : difference) - h * jv[i];
        }
        return Norm(r);
    };
    const double coarse = residual(1e-2);
    const double fine = residual(1e-3);
    RC_ASSERT(fine <= coarse / 30.0 + 1e-10 * (1.0 + Norm(base)));
}

/// <J v, w> = <v, J^T w> over every sample.
template <typename SA, typename SB>
void RequireTranspose(const Case<SA, SB>& c)
{
    const auto [references, v] = RandomSeed(ParameterCount<SA, SB>());
    Model model = ModelFor<SA, SB>(references);
    const std::vector<double> jv = Jvp(model, c, v);
    std::vector<double> w;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        w.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    RC_ASSERT(model.Sensitivity.SetAdjoint(w).has_value());
    const auto vjp = Stereon::IntersectVjp(Surface{c.A}, model.BindingsA, Surface{c.B},
                                           model.BindingsB, c.At, c.Samples, model.Sensitivity);
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

/// Every sample well away from a tangency and, for random geometry, from
/// the sphere's poles, where the angle u turns too fast for a fixed step.
template <typename SA, typename SB>
bool WellConditioned(const Case<SA, SB>& c, const double poleReach = 1.2)
{
    for (const BranchSample& sample : c.Samples) {
        const auto lin = Stereon::Linearize(Surface{c.A}, Surface{c.B},
                                            c.At.Branches[sample.Branch], sample.Parameter);
        if (!lin || lin->Status != DerivativeStatus::Valid || lin->Conditioning < 0.25) {
            return false;
        }
        if constexpr (std::is_same_v<SB, Sphere> || std::is_same_v<SB, Ellipsoid>) {
            if (std::abs(lin->OnB.Y) > poleReach) {
                return false;
            }
        }
    }
    return true;
}

/// A random plane cutting a random sphere well inside it: a circle.
Case<Plane, Sphere> RandomPlaneSphere()
{
    const Point3 center{Coordinate(), Coordinate(), Coordinate()};
    const double radius = *rc::gen::inRange(8, 33) / 8.0;
    const Vec3 normal{Coordinate(), Coordinate(), Coordinate()};
    RC_PRE(Length(normal) >= 0.5);
    const double offset = *rc::gen::inRange(-5, 6) / 8.0 * radius;
    const Plane plane =
        IntersectTest::PlaneWithNormal(center + offset * (normal / Length(normal)), normal);
    const Sphere sphere{center, radius, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    // A circle through a pole of the sphere fails closed (ADR-0009): its
    // pcurve has no chart there. Not this test's subject.
    const Point3 north = center + Vec3{0.0, 0.0, radius};
    const Point3 south = center - Vec3{0.0, 0.0, radius};
    const Vec3 unit = normal / Length(normal);
    RC_PRE(std::abs(Dot(unit, north - plane.Origin)) >= 0.125 * radius &&
           std::abs(Dot(unit, south - plane.Origin)) >= 0.125 * radius);
    Case<Plane, Sphere> c{plane, sphere, Meet(Surface{plane}, Surface{sphere}), {}};
    RC_PRE(!c.At.Branches.empty() && c.At.Overlaps.empty());
    c.Samples = SamplesOf(c.At);
    RC_PRE(WellConditioned(c));
    return c;
}

/// A random plane cutting a random ellipsoid in a random frame, well
/// inside it and away from its poles: an exact ellipse, its pcurve on the
/// ellipsoid a certified fit.
Case<Plane, Ellipsoid> RandomPlaneEllipsoid()
{
    const Point3 center{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 axis{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 normal{Coordinate(), Coordinate(), Coordinate()};
    RC_PRE(Length(axis) >= 0.5 && Length(normal) >= 0.5);
    const Cylinder frame = IntersectTest::CylinderAbout(center, axis, 1.0);
    const auto semi = [] { return *rc::gen::inRange(8, 33) / 8.0; };
    const Ellipsoid egg{center, semi(), semi(), semi(), frame.XAxis, frame.YAxis};
    const double reach = std::min({egg.A, egg.B, egg.C});
    const Vec3 unit = normal / Length(normal);
    const Plane plane = IntersectTest::PlaneWithNormal(
        center + (*rc::gen::inRange(-4, 5) / 8.0 * reach) * unit, normal);
    // Away from the poles, where the plane-ellipsoid handler fails closed.
    for (const double side : {1.0, -1.0}) {
        const Point3 pole = center + (side * egg.C) * Axis(egg);
        RC_PRE(std::abs(Dot(unit, pole - plane.Origin)) >= 0.125 * reach);
    }
    Case<Plane, Ellipsoid> c{plane, egg, Meet(Surface{plane}, Surface{egg}), {}};
    RC_PRE(!c.At.Branches.empty() && c.At.Overlaps.empty());
    c.Samples = SamplesOf(c.At);
    RC_PRE(WellConditioned(c));
    return c;
}

/// A random plane across a random cylinder, steep to its axis: an ellipse.
Case<Plane, Cylinder> RandomPlaneCylinder()
{
    const Point3 origin{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 axis{Coordinate(), Coordinate(), Coordinate()};
    const Vec3 normal{Coordinate(), Coordinate(), Coordinate()};
    RC_PRE(Length(axis) >= 0.5 && Length(normal) >= 0.5);
    RC_PRE(std::abs(Dot(axis, normal)) >= 0.5 * Length(axis) * Length(normal));
    const double radius = *rc::gen::inRange(8, 33) / 8.0;
    const Cylinder cylinder = IntersectTest::CylinderAbout(origin, axis, radius);
    const Plane plane = IntersectTest::PlaneWithNormal(
        origin + Vec3{Coordinate(), Coordinate(), Coordinate()}, normal);
    Case<Plane, Cylinder> c{plane, cylinder, Meet(Surface{plane}, Surface{cylinder}), {}};
    RC_PRE(!c.At.Branches.empty() && c.At.Overlaps.empty());
    c.Samples = SamplesOf(c.At);
    RC_PRE(WellConditioned(c));
    return c;
}

/// Cylinders inside a sphere, off its centre: two traced and certified
/// closed curves each. Traced once: a trace per random case would make
/// this test slow, and the products' properties are randomised over the
/// seeds and references instead.
const std::vector<Case<Cylinder, Sphere>>& TracedCases()
{
    static const std::vector<Case<Cylinder, Sphere>> cases = [] {
        std::vector<Case<Cylinder, Sphere>> made;
        const Sphere sphere{Point3{0.5, -1.0, 2.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        for (const auto& [offset, axis, tube] :
             {std::tuple{Vec3{0.5, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.0},
              std::tuple{Vec3{0.0, 0.4, 0.3}, Vec3{1.0, 0.5, -0.25}, 0.75}}) {
            const Cylinder cylinder =
                IntersectTest::CylinderAbout(sphere.Center + offset, axis, tube);
            auto at = Stereon::Intersect(Surface{cylinder}, Surface{sphere}, ctx);
            REQUIRE(at.has_value());
            REQUIRE(at->Branches.size() >= 2u);
            for (const IntersectionBranch& branch : at->Branches) {
                REQUIRE(branch.IsCertified());
                REQUIRE(branch.Bound > 0.0);
            }
            Case<Cylinder, Sphere> c{cylinder, sphere, *std::move(at), {}};
            c.Samples = SamplesOf(c.At);
            // Some samples come within 0.3 of a pole: a fixed step 1e-2 stays far
            // from it.
            REQUIRE(WellConditioned(c, 1.4));
            made.push_back(std::move(c));
        }
        return made;
    }();
    return cases;
}

/// One of them, at random.
const Case<Cylinder, Sphere>& RandomTracedCase()
{
    return TracedCases()[static_cast<std::size_t>(*rc::gen::inRange(0, 2))];
}

/// The closed form: the plane z = h cuts the sphere of radius R about the
/// origin in a circle of radius rho = sqrt(R^2 - h^2). Normal to the
/// circle, its point moves radially by rho' = (R R' - h h') / rho and up
/// by h'.
struct CircleCase {
    double Height = 0.6;
    double Radius = 1.5;
    Plane Cut{Point3{0.0, 0.0, 0.6}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    Sphere Ball{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
};

} // namespace

TEST_CASE("a plane's circle on a sphere moves as its closed form says",
          "[intersect][linearization][branch]")
{
    const CircleCase c;
    const auto at = Stereon::Intersect(Surface{c.Cut}, Surface{c.Ball}, ctx);
    REQUIRE(at.has_value());
    REQUIRE(at->Branches.size() == 1u);
    REQUIRE(at->Branches[0].IsCertified());
    const std::vector<BranchSample> samples = SamplesOf(*at);
    // The plane's height, then the sphere's radius.
    ParameterSet set;
    const ParameterId height = *set.Add("h", ParameterQuantity::Length, 1.0);
    const ParameterId radius = *set.Add("R", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    const std::vector<SurfaceBinding> onPlane{{height, SurfaceComponent::OriginZ}};
    const std::vector<SurfaceBinding> onSphere{{radius, SurfaceComponent::Radius}};
    const double rho = std::sqrt(c.Radius * c.Radius - c.Height * c.Height);
    for (const auto& [seed, radial, up] :
         {std::tuple{std::vector<double>{1.0, 0.0}, -c.Height / rho, 1.0},
          std::tuple{std::vector<double>{0.0, 1.0}, c.Radius / rho, 0.0}}) {
        REQUIRE(sensitivity.SetSeed(seed).has_value());
        const auto jvp = Stereon::IntersectJvp(Surface{c.Cut}, onPlane, Surface{c.Ball}, onSphere,
                                               *at, samples, sensitivity);
        REQUIRE(jvp.has_value());
        REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
        const auto values = *jvp->GetValues();
        REQUIRE(values.size() == Values * samples.size());
        for (std::size_t k = 0; k < samples.size(); ++k) {
            const auto lin = Stereon::Linearize(Surface{c.Cut}, Surface{c.Ball}, at->Branches[0],
                                                samples[k].Parameter);
            REQUIRE(lin.has_value());
            const Vec3 out{lin->Position.X / rho, lin->Position.Y / rho, 0.0};
            const Vec3 expected = radial * out + Vec3{0.0, 0.0, up};
            CHECK_THAT(values[Values * k], WithinAbs(expected.X, 1e-12));
            CHECK_THAT(values[Values * k + 1], WithinAbs(expected.Y, 1e-12));
            CHECK_THAT(values[Values * k + 2], WithinAbs(expected.Z, 1e-12));
            // On the plane, (u, v) are the point's x and y.
            CHECK_THAT(values[Values * k + 3], WithinAbs(expected.X, 1e-12));
            CHECK_THAT(values[Values * k + 4], WithinAbs(expected.Y, 1e-12));
        }
    }
}

TEST_CASE("Linearize exposes a branch point's system", "[intersect][linearization][branch]")
{
    const CircleCase c;
    const auto at = Stereon::Intersect(Surface{c.Cut}, Surface{c.Ball}, ctx);
    REQUIRE(at.has_value());
    const IntersectionBranch& branch = at->Branches[0];
    const double t =
        branch.Range.GetLower() + 0.3 * (branch.Range.GetUpper() - branch.Range.GetLower());
    const auto lin = Stereon::Linearize(Surface{c.Cut}, Surface{c.Ball}, branch, t);
    REQUIRE(lin.has_value());
    CHECK(lin->Status == DerivativeStatus::Valid);
    CHECK(Length(lin->Residual) <= 1e-14);
    CHECK_THAT(Length(lin->Tangent), WithinAbs(1.0, 1e-15));
    // The tracking constraint holds at the refined point.
    CHECK_THAT(Dot(lin->Tangent, lin->Position - Eval(branch.Curve3d, t)), WithinAbs(0.0, 1e-14));
    // The angle between a horizontal plane and the sphere at height h has
    // sine rho / R.
    const double rho = std::sqrt(c.Radius * c.Radius - c.Height * c.Height);
    CHECK_THAT(lin->Conditioning, WithinAbs(rho / c.Radius, 1e-12));
    CHECK(Stereon::Linearize(Surface{c.Cut}, Surface{c.Ball}, branch, branch.Range.GetUpper() + 1.0)
              .error()
              .GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("branch points: second-order Taylor decay and the transpose",
          "[intersect][linearization][branch][property]")
{
    rc::prop("for random planes through spheres", [] { RequireTaylor(RandomPlaneSphere()); });
    rc::prop("for random planes across cylinders", [] { RequireTaylor(RandomPlaneCylinder()); });
    rc::prop("for cylinders through a sphere, traced", [] { RequireTaylor(RandomTracedCase()); });
    rc::prop("for random planes through ellipsoids", [] { RequireTaylor(RandomPlaneEllipsoid()); });
    rc::prop("<J v, w> = <v, J^T w> for ellipsoids",
             [] { RequireTranspose(RandomPlaneEllipsoid()); });
    rc::prop("<J v, w> = <v, J^T w> for spheres", [] { RequireTranspose(RandomPlaneSphere()); });
    rc::prop("<J v, w> = <v, J^T w> for cylinders",
             [] { RequireTranspose(RandomPlaneCylinder()); });
    rc::prop("<J v, w> = <v, J^T w> for traced branches",
             [] { RequireTranspose(RandomTracedCase()); });
}

TEST_CASE("a branch point without a derivative says so and carries none",
          "[intersect][linearization][branch]")
{
    ParameterSet set;
    const ParameterId z = *set.Add("z", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<SurfaceBinding> bindings{{z, SurfaceComponent::OriginZ}};
    const std::vector<SurfaceBinding> none;
    const auto jvp = [&](const Surface& a, const Surface& b, const SurfaceSurfaceIntersection& at,
                         const std::vector<BranchSample>& samples) {
        return Stereon::IntersectJvp(a, bindings, b, none, at, samples, sensitivity);
    };
    SECTION("a plane tangent to a cylinder along a line: Nondifferentiable")
    {
        const Plane wall = IntersectTest::PlaneWithNormal(Point3{2.0, 0.0, 0.0}, Vec3{1, 0, 0});
        const Cylinder post{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const auto at = Stereon::Intersect(
            Surface{wall}, Box2Of(-5.0, 5.0), Surface{post},
            Stereon::Box2{Stereon::Interval{0.0, tau}, Stereon::Interval{-5.0, 5.0}}, ctx);
        REQUIRE(at.has_value());
        REQUIRE(at->Branches.size() == 1u);
        REQUIRE(at->Branches[0].Contact == Stereon::ContactType::Tangential);
        const auto result = jvp(Surface{wall}, Surface{post}, *at, SamplesOf(*at));
        REQUIRE(result.has_value());
        CHECK(result->GetStatus() == DerivativeStatus::Nondifferentiable);
        CHECK_FALSE(result->GetValues().has_value());
    }
    SECTION("coincident planes: an overlap, Nondifferentiable")
    {
        const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
        const auto at = Stereon::Intersect(Surface{ground}, Box2Of(-1.0, 1.0), Surface{ground},
                                           Box2Of(-1.0, 1.0), ctx);
        REQUIRE(at.has_value());
        REQUIRE_FALSE(at->Overlaps.empty());
        const auto result = jvp(Surface{ground}, Surface{ground}, *at, {});
        REQUIRE(result.has_value());
        CHECK(result->GetStatus() == DerivativeStatus::Nondifferentiable);
    }
    SECTION("an uncertified branch: Indeterminate")
    {
        const CircleCase c;
        auto at = Stereon::Intersect(Surface{c.Cut}, Surface{c.Ball}, ctx);
        REQUIRE(at.has_value());
        at->Branches[0].Bound = std::numeric_limits<double>::infinity();
        const auto result = jvp(Surface{c.Cut}, Surface{c.Ball}, *at, SamplesOf(*at));
        REQUIRE(result.has_value());
        CHECK(result->GetStatus() == DerivativeStatus::Indeterminate);
    }
}

TEST_CASE("branch products refuse bad input and other surfaces",
          "[intersect][linearization][branch]")
{
    const CircleCase c;
    const auto at = Stereon::Intersect(Surface{c.Cut}, Surface{c.Ball}, ctx);
    REQUIRE(at.has_value());
    ParameterSet set;
    const ParameterId z = *set.Add("z", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<SurfaceBinding> bindings{{z, SurfaceComponent::OriginZ}};
    const std::vector<SurfaceBinding> none;
    const std::vector<BranchSample> samples = SamplesOf(*at);
    SECTION("a sample of a branch the result does not have")
    {
        const std::vector<BranchSample> missing{{1, 0.0}};
        CHECK(Stereon::IntersectJvp(Surface{c.Cut}, bindings, Surface{c.Ball}, none, *at, missing,
                                    sensitivity)
                  .error()
                  .GetCode() == ErrorCode::InvalidInput);
    }
    SECTION("a result from other geometry: stale")
    {
        Sphere bigger = c.Ball;
        bigger.Radius = 1.6;
        CHECK(Stereon::IntersectJvp(Surface{c.Cut}, bindings, Surface{bigger}, none, *at, samples,
                                    sensitivity)
                  .error()
                  .GetCode() == ErrorCode::InvalidInput);
    }
    SECTION("a radius bound on a plane")
    {
        const std::vector<SurfaceBinding> radius{{z, SurfaceComponent::Radius}};
        CHECK(Stereon::IntersectJvp(Surface{c.Cut}, radius, Surface{c.Ball}, none, *at, samples,
                                    sensitivity)
                  .error()
                  .GetCode() == ErrorCode::InvalidInput);
    }
    SECTION("an adjoint of the wrong length")
    {
        REQUIRE(sensitivity.SetAdjoint(std::vector<double>(Values, 1.0)).has_value());
        CHECK(Stereon::IntersectVjp(Surface{c.Cut}, bindings, Surface{c.Ball}, none, *at, samples,
                                    sensitivity)
                  .error()
                  .GetCode() == ErrorCode::InvalidInput);
    }
    SECTION("a surface with no bindings, an extrusion: NotSupported")
    {
        const Surface extrusion{Stereon::ExtrusionSurface(
            Stereon::Curve{Stereon::Line{Point3{}, Vec3{1.0, 0.0, 0.0}}}, Vec3{0.0, 1.0, 0.0})};
        const auto result = Stereon::IntersectJvp(Surface{c.Cut}, bindings, extrusion, none, *at,
                                                  samples, sensitivity);
        REQUIRE(result.has_value());
        CHECK(result->GetStatus() == DerivativeStatus::NotSupported);
    }
}

TEST_CASE("two B-spline patches' branch moves with a pole as its closed form says",
          "[intersect][linearization][branch][spline-surface]")
{
    // A: the flat bilinear patch z = 0 over the unit square, pole (i, j) at
    // (i, j, 0). B: the vertical patch x = 0.5, (u, v) to (0.5, u, 2 v - 1).
    // They meet along x = 0.5, z = 0. Raising A's pole (1, 1), index 3, by
    // p makes A's height p u v: at the branch point (0.5, y, 0), which the
    // normal-plane constraint keeps at its y, the point rises by 0.5 y,
    // B's v by half that, and nothing else moves.
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto a = Stereon::MakeBSplineSurface(
        1, 1, knots, knots, {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}},
        {});
    const auto b = Stereon::MakeBSplineSurface(
        1, 1, knots, knots, {{0.5, 0.0, -1.0}, {0.5, 0.0, 1.0}, {0.5, 1.0, -1.0}, {0.5, 1.0, 1.0}},
        {});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const auto at = Stereon::Intersect(Surface{*a}, Surface{*b}, ctx);
    INFO(IntersectTest::ErrorText(at));
    REQUIRE(at.has_value());
    REQUIRE_FALSE(at->Branches.empty());
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<SurfaceBinding> lift{{p, SurfaceComponent::PoleZ, 3}};
    const std::vector<BranchSample> samples = SamplesOf(*at);
    const auto jvp =
        Stereon::IntersectJvp(Surface{*a}, lift, Surface{*b}, {}, *at, samples, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    const auto values = *jvp->GetValues();
    REQUIRE(values.size() == Values * samples.size());
    for (std::size_t k = 0; k < samples.size(); ++k) {
        const auto lin = Stereon::Linearize(Surface{*a}, Surface{*b},
                                            at->Branches[samples[k].Branch], samples[k].Parameter);
        REQUIRE(lin.has_value());
        const double y = lin->Position.Y;
        CHECK_THAT(values[Values * k], WithinAbs(0.0, 1e-12));
        CHECK_THAT(values[Values * k + 1], WithinAbs(0.0, 1e-12));
        CHECK_THAT(values[Values * k + 2], WithinAbs(0.5 * y, 1e-12));
        CHECK_THAT(values[Values * k + 3], WithinAbs(0.0, 1e-12));
        CHECK_THAT(values[Values * k + 4], WithinAbs(0.0, 1e-12));
        CHECK_THAT(values[Values * k + 5], WithinAbs(0.0, 1e-12));
        CHECK_THAT(values[Values * k + 6], WithinAbs(0.25 * y, 1e-12));
    }
}
