// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The derivative products of a hyperbola's and a parabola's crossings of a
// surface (ADR-0021, sections 6 and 9; the v0.2 linearisation gate): the
// handlers of #292 differentiated through F = C(t) - S(u, v) = 0, as a
// circle's are. Each partial is checked against a closed form: the hyperbola
// (3 cosh t, 0, 4 sinh t) cut by z = 3 at sinh t = 3/4, and the parabola
// (t^2, 0, t) cut by z = 2 at t = 2. On random conics through planes, spheres
// and cylinders the Taylor residual |y(p + h v) - y(p) - h J v| decays with
// the second power of h and the reverse product is the forward one's
// transpose. A touch and a curve in the surface are Nondifferentiable, an
// undecided point Indeterminate, a grazing crossing NearEvent, and a binding
// of a value the conic does not have is refused.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using Stereon::Circle;
using Stereon::ConicBinding;
using Stereon::ConicComponent;
using Stereon::ContactType;
using Stereon::CurveSurfaceIntersection;
using Stereon::CurveSurfacePoint;
using Stereon::Cylinder;
using Stereon::DerivativeStatus;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Parabola;
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
using Stereon::Vec3;

namespace {

constexpr std::size_t Values = Stereon::CurveSurfacePointValues;
constexpr double tau = 2.0 * std::numbers::pi;

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

/// The stretch of an open conic every case is asked over.
const Interval Stretch{-1.5, 1.5};

bool IsTurn(const ConicComponent c)
{
    return c == ConicComponent::RotationX || c == ConicComponent::RotationY ||
           c == ConicComponent::RotationZ;
}

/// The conic's bindings in the order of their parameters: its centre (a
/// parabola's vertex), its radii or focal length, then its frame's turn.
template <typename C>
std::vector<ConicComponent> ComponentsOf()
{
    std::vector<ConicComponent> components{ConicComponent::CenterX, ConicComponent::CenterY,
                                           ConicComponent::CenterZ};
    if constexpr (std::is_same_v<C, Hyperbola>) {
        components.insert(components.end(),
                          {ConicComponent::MajorRadius, ConicComponent::MinorRadius});
    } else {
        components.push_back(ConicComponent::FocalLength);
    }
    components.insert(components.end(), {ConicComponent::RotationX, ConicComponent::RotationY,
                                         ConicComponent::RotationZ});
    return components;
}

Point3& CentreOf(Hyperbola& h)
{
    return h.Center;
}
Point3& CentreOf(Parabola& p)
{
    return p.Vertex;
}

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

/// v turned by the rotation vector w (Rodrigues): an exact rotation, whose
/// first-order motion of a point is w x (P - centre).
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

template <typename C, typename S>
CurveSurfaceIntersection Meet(const C& conic, const S& surface)
{
    const Surface s{surface};
    auto result = Intersect(Stereon::Curve{conic}, Stretch, s, Domain(s), ctx);
    if (!result) {
        // A fail-closed refusal (a contact within rounding of a tangency or
        // of the stretch's end, ADR-0029) is no geometry a derivative applies
        // to: the case is discarded, never counted. Any other error fails.
        RC_PRE(result.error().GetCode() != ErrorCode::IntersectionFailed);
        RC_FAIL(IntersectTest::ErrorText(result));
    }
    return *std::move(result);
}

/// The points as the product's values are laid out: x, then t, u, v.
template <typename C, typename S>
std::vector<double> Points(const C& conic, const S& surface)
{
    std::vector<double> values;
    for (const CurveSurfacePoint& p : Meet(conic, surface).Points) {
        values.insert(values.end(), {p.Position.X, p.Position.Y, p.Position.Z, p.Parameter,
                                     p.OnSurface.X, p.OnSurface.Y});
    }
    return values;
}

struct Model {
    SensitivityContext Sensitivity{ParameterSet{}};
    std::vector<ConicBinding> ConicBindings;
    std::vector<SurfaceBinding> SurfaceBindings;
};

/// One parameter per binding: the conic's, then the surface's origin.
template <typename C>
Model ModelFor(const std::vector<double>& references)
{
    const std::vector<ConicComponent> conic = ComponentsOf<C>();
    ParameterSet set;
    std::vector<ParameterId> ids;
    for (std::size_t k = 0; k < references.size(); ++k) {
        const bool turn = k < conic.size() && IsTurn(conic[k]);
        ids.push_back(*set.Add("p" + std::to_string(k),
                               turn ? ParameterQuantity::Angle : ParameterQuantity::Length,
                               references[k]));
    }
    Model model{SensitivityContext(std::move(set)), {}, {}};
    for (std::size_t k = 0; k < conic.size(); ++k) {
        model.ConicBindings.push_back({ids[k], conic[k]});
    }
    const std::array origin{SurfaceComponent::OriginX, SurfaceComponent::OriginY,
                            SurfaceComponent::OriginZ};
    for (std::size_t k = 0; k < 3; ++k) {
        model.SurfaceBindings.push_back({ids[conic.size() + k], origin[k]});
    }
    return model;
}

/// The conic and the surface moved by `step` times the seed.
template <typename C, typename S>
std::pair<C, S> Moved(C conic, S surface, const std::vector<double>& seed,
                      const std::vector<double>& references, const double step)
{
    const auto d = [&](const std::size_t k) { return step * seed[k] * references[k]; };
    CentreOf(conic) = CentreOf(conic) + Vec3{d(0), d(1), d(2)};
    std::size_t next = 4;
    if constexpr (std::is_same_v<C, Hyperbola>) {
        conic.MajorRadius += d(3);
        conic.MinorRadius += d(4);
        next = 5;
    } else {
        conic.FocalLength += d(3);
    }
    const Vec3 w = d(next) * conic.XAxis + d(next + 1) * conic.YAxis +
                   d(next + 2) * Cross(conic.XAxis, conic.YAxis);
    conic.XAxis = Rotated(conic.XAxis, w);
    conic.YAxis = Rotated(conic.YAxis, w);
    OriginOf(surface) = OriginOf(surface) + Vec3{d(next + 3), d(next + 4), d(next + 5)};
    return {conic, surface};
}

double Norm(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

/// Every crossing well away from a tangency, from the surface's poles and
/// from the stretch's ends, so that the steps below stay in one regime.
template <typename C, typename S>
bool WellConditioned(const C& conic, const S& surface)
{
    const Surface s{surface};
    const auto result = Intersect(Stereon::Curve{conic}, Stretch, s, Domain(s), ctx);
    if (!result) {
        RC_ASSERT(result.error().GetCode() == ErrorCode::IntersectionFailed);
        return false;
    }
    if (result->Points.empty() || !result->Overlaps.empty()) {
        return false;
    }
    return std::ranges::all_of(result->Points, [&](const CurveSurfacePoint& p) {
        const auto lin = Stereon::Linearize(conic, s, p);
        return lin && lin->Status == DerivativeStatus::Valid && lin->Conditioning >= 0.25 &&
               Length(lin->Su) >= 0.1 * Length(lin->Sv) &&
               Length(lin->Sv) >= 0.1 * Length(lin->Su) && std::abs(p.Parameter) <= 1.25;
    });
}

/// The Taylor residual falls by second order, and <J v, w> = <v, J^T w>.
template <typename C, typename S>
void RequireProducts(const C& conic, const S& surface)
{
    const std::size_t count = ComponentsOf<C>().size() + 3;
    std::vector<double> references;
    std::vector<double> seed;
    for (std::size_t k = 0; k < count; ++k) {
        references.push_back(*rc::gen::inRange(1, 20) / 10.0);
        seed.push_back(*rc::gen::inRange(-10, 11) / 10.0);
    }
    Model model = ModelFor<C>(references);
    const CurveSurfaceIntersection at = Meet(conic, surface);
    RC_ASSERT(model.Sensitivity.SetSeed(seed).has_value());
    const auto jvp = Stereon::IntersectJvp(conic, model.ConicBindings, Surface{surface},
                                           model.SurfaceBindings, at, model.Sensitivity);
    if (!jvp) {
        RC_FAIL(IntersectTest::ErrorText(jvp));
    }
    RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
    const std::vector<double> jv(jvp->GetValues()->begin(), jvp->GetValues()->end());
    const std::vector<double> base = Points(conic, surface);
    // A conic that passes near the surface away from its crossings gains two
    // under a step, and one near the stretch's end may leave it: a new
    // regime, not an error of the derivative. Such steps are discarded.
    for (const double h : {1e-2, 1e-3}) {
        const auto [c, s] = Moved(conic, surface, seed, references, h);
        RC_PRE(Points(c, s).size() == base.size());
    }
    const auto residual = [&](const double h) {
        const auto [c, s] = Moved(conic, surface, seed, references, h);
        const std::vector<double> moved = Points(c, s);
        const std::size_t points = base.size() / Values;
        std::vector<double> r(base.size());
        for (std::size_t i = 0; i < points; ++i) {
            const auto gap = [&](const std::size_t j) {
                return Distance(
                    Point3{base[Values * i], base[Values * i + 1], base[Values * i + 2]},
                    Point3{moved[Values * j], moved[Values * j + 1], moved[Values * j + 2]});
            };
            std::size_t nearest = 0;
            for (std::size_t j = 1; j < points; ++j) {
                nearest = gap(j) < gap(nearest) ? j : nearest;
            }
            for (std::size_t k = 0; k < Values; ++k) {
                const double difference = moved[Values * nearest + k] - base[Values * i + k];
                // u and v may be angles, and a step may carry one over a seam;
                // the conic's own t is not one.
                const bool angle = k >= 4;
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
    const auto vjp = Stereon::IntersectVjp(conic, model.ConicBindings, Surface{surface},
                                           model.SurfaceBindings, at, model.Sensitivity);
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

double Coordinate()
{
    return *rc::gen::inRange(-40, 41) / 8.0;
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

template <typename C>
C RandomConic()
{
    const auto [x, y] = RandomFrame();
    const Point3 centre{Coordinate() / 4.0, Coordinate() / 4.0, Coordinate() / 4.0};
    if constexpr (std::is_same_v<C, Hyperbola>) {
        return Hyperbola{centre, *rc::gen::inRange(5, 16) / 10.0, *rc::gen::inRange(5, 21) / 10.0,
                         x, y};
    } else {
        return Parabola{centre, *rc::gen::inRange(2, 11) / 10.0, x, y};
    }
}

template <typename C>
void RequireProperties()
{
    const C conic = RandomConic<C>();
    // Every surface passes through the conic's point at t0, well inside the
    // stretch, so that most cases cross: placed at random, few did, and the
    // property gave up on its discards.
    const Point3 on = Eval(conic, *rc::gen::inRange(-8, 9) / 8.0);
    const double radius = *rc::gen::inRange(15, 41) / 10.0;
    switch (*rc::gen::inRange(0, 3)) {
    case 0: {
        const Vec3 away = RandomFrame().second;
        const Sphere sphere{on + radius * away, radius, ex, ey};
        RC_PRE(WellConditioned(conic, sphere));
        RequireProducts(conic, sphere);
        break;
    }
    case 1: {
        const Vec3 away = RandomFrame().first;
        const Point3 foot = on + radius * away;
        const Cylinder cylinder{Point3{foot.X, foot.Y, 0.0}, radius, ex, ey};
        RC_PRE(WellConditioned(conic, cylinder));
        RequireProducts(conic, cylinder);
        break;
    }
    default: {
        const auto [x, y] = RandomFrame();
        const Plane plane{on, x, y};
        RC_PRE(WellConditioned(conic, plane));
        RequireProducts(conic, plane);
        break;
    }
    }
}

/// The hyperbola (3 cosh t, 0, 4 sinh t) and the parabola (t^2, 0, t), both
/// in the xz-plane, their frame's normal X x Y = -y.
const Hyperbola UpHyperbola{Point3{}, 3.0, 4.0, ex, ez};
const Parabola UpParabola{Point3{}, 0.25, ex, ez};

Plane ShelfAt(const double height)
{
    return Plane{Point3{0.0, 0.0, height}, ex, ey};
}

template <typename C>
CurveSurfaceIntersection On(const C& conic, const Plane& plane)
{
    // Wider than the properties' stretch: the parabola meets z = 2 at t = 2.
    const Surface s{plane};
    const auto at = Intersect(Stereon::Curve{conic}, Interval{-3.0, 3.0}, s, Domain(s), ctx);
    INFO(IntersectTest::ErrorText(at));
    REQUIRE(at.has_value());
    return *at;
}

/// The forward product of the one crossing of `plane` for one binding moved
/// at unit rate: (x', t', u', v').
template <typename C>
std::vector<double> RateOf(const C& conic, const Plane& plane, const ConicComponent component)
{
    ParameterSet set;
    const ParameterId p = *set.Add(
        "p", IsTurn(component) ? ParameterQuantity::Angle : ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<ConicBinding> bindings{{p, component}};
    const CurveSurfaceIntersection at = On(conic, plane);
    REQUIRE(at.Points.size() == 1);
    REQUIRE(at.Points.front().Decided);
    const auto jvp = Stereon::IntersectJvp(conic, bindings, Surface{plane}, {}, at, sensitivity);
    INFO(IntersectTest::ErrorText(jvp));
    REQUIRE(jvp.has_value());
    REQUIRE(jvp->GetStatus() == DerivativeStatus::Valid);
    return {jvp->GetValues()->begin(), jvp->GetValues()->end()};
}

/// The values against x' and t'; on the shelf u' and v' are x' along x and y.
void RequireRate(const std::vector<double>& v, const Vec3& x, const double t)
{
    REQUIRE(v.size() == Values);
    CHECK(std::abs(v[0] - x.X) <= 1e-13);
    CHECK(std::abs(v[1] - x.Y) <= 1e-13);
    CHECK(std::abs(v[2] - x.Z) <= 1e-13);
    CHECK(std::abs(v[3] - t) <= 1e-13);
    CHECK(std::abs(v[4] - x.X) <= 1e-13);
    CHECK(std::abs(v[5] - x.Y) <= 1e-13);
}

} // namespace

TEST_CASE("a hyperbola's crossing of a plane moves as its closed form says",
          "[intersect][linearization][open-conic]")
{
    // 4 sinh t = 3: sinh t = 3/4, cosh t = 5/4, at (15/4, 0, 3), where
    // C' = (9/4, 0, 5). A binding moving the point by q at fixed t moves t by
    // -q.z / 5 and the crossing by q + C' t'.
    const Plane shelf = ShelfAt(3.0);
    const auto rate = [&](const ConicComponent c) { return RateOf(UpHyperbola, shelf, c); };
    SECTION("its centre")
    {
        RequireRate(rate(ConicComponent::CenterX), Vec3{1.0, 0.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::CenterY), Vec3{0.0, 1.0, 0.0}, 0.0);
        // Rising, it meets the shelf lower on itself: t' = -1/5.
        RequireRate(rate(ConicComponent::CenterZ), Vec3{-9.0 / 20.0, 0.0, 0.0}, -1.0 / 5.0);
    }
    SECTION("its radii")
    {
        // a moves the point by cosh t along x and leaves t; b moves it by
        // sinh t = 3/4 up, so t' = -3/20 and x' = 9/4 t' = -27/80.
        RequireRate(rate(ConicComponent::MajorRadius), Vec3{5.0 / 4.0, 0.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::MinorRadius), Vec3{-27.0 / 80.0, 0.0, 0.0}, -3.0 / 20.0);
    }
    SECTION("its frame turning")
    {
        // About X = x: x cross (15/4, 0, 3) = (0, -3, 0), out of the plane of
        // the curve, t' = 0. About Y = z: (0, 15/4, 0). About X x Y = -y:
        // (-3, 0, 15/4), so t' = -3/4 and x' = -3 - 27/16 = -75/16.
        RequireRate(rate(ConicComponent::RotationX), Vec3{0.0, -3.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::RotationY), Vec3{0.0, 15.0 / 4.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::RotationZ), Vec3{-75.0 / 16.0, 0.0, 0.0}, -3.0 / 4.0);
    }
}

TEST_CASE("a parabola's crossing of a plane moves as its closed form says",
          "[intersect][linearization][open-conic]")
{
    // t = 2 at (4, 0, 2), where C' = (t / (2 f), 0, 1) = (4, 0, 1).
    const Plane shelf = ShelfAt(2.0);
    const auto rate = [&](const ConicComponent c) { return RateOf(UpParabola, shelf, c); };
    SECTION("its vertex")
    {
        RequireRate(rate(ConicComponent::CenterX), Vec3{1.0, 0.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::CenterY), Vec3{0.0, 1.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::CenterZ), Vec3{-4.0, 0.0, 0.0}, -1.0);
    }
    SECTION("its focal length")
    {
        // x = t^2 / (4 f) at fixed t = 2: dx/df = -t^2 / (4 f^2) = -16.
        RequireRate(rate(ConicComponent::FocalLength), Vec3{-16.0, 0.0, 0.0}, 0.0);
    }
    SECTION("its frame turning")
    {
        // About X: x cross (4, 0, 2) = (0, -2, 0). About Y = z: (0, 4, 0).
        // About X x Y = -y: (-2, 0, 4), so t' = -4 and x' = -2 - 16 = -18.
        RequireRate(rate(ConicComponent::RotationX), Vec3{0.0, -2.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::RotationY), Vec3{0.0, 4.0, 0.0}, 0.0);
        RequireRate(rate(ConicComponent::RotationZ), Vec3{-18.0, 0.0, 0.0}, -4.0);
    }
}

TEST_CASE("hyperbola and parabola crossings: second-order Taylor decay and the transpose",
          "[intersect][linearization][open-conic][property]")
{
    rc::prop("hyperbolas through spheres, cylinders and planes",
             [] { RequireProperties<Hyperbola>(); });
    rc::prop("parabolas through spheres, cylinders and planes",
             [] { RequireProperties<Parabola>(); });
}

TEST_CASE("Linearize exposes an open conic's residual, system and regularity",
          "[intersect][linearization][open-conic]")
{
    const Plane shelf = ShelfAt(3.0);
    const CurveSurfaceIntersection at = On(UpHyperbola, shelf);
    REQUIRE(at.Points.size() == 1);
    const auto lin = Stereon::Linearize(UpHyperbola, Surface{shelf}, at.Points.front());
    REQUIRE(lin.has_value());
    CHECK(Length(lin->Residual) <= 1e-14);
    CHECK(Length(lin->Tangent - Vec3{9.0 / 4.0, 0.0, 5.0}) <= 1e-14);
    CHECK(lin->Su == ex);
    CHECK(lin->Sv == ey);
    // The sine of the angle between (9/4, 0, 5) and the shelf.
    CHECK(std::abs(lin->Conditioning - 20.0 / std::sqrt(481.0)) <= 1e-15);
    CHECK(lin->Status == DerivativeStatus::Valid);

    const CurveSurfaceIntersection onParabola = On(UpParabola, ShelfAt(2.0));
    const auto parabola =
        Stereon::Linearize(UpParabola, Surface{ShelfAt(2.0)}, onParabola.Points.front());
    REQUIRE(parabola.has_value());
    CHECK(Length(parabola->Tangent - Vec3{4.0, 0.0, 1.0}) <= 1e-14);
    CHECK(std::abs(parabola->Conditioning - 1.0 / std::sqrt(17.0)) <= 1e-15);
    CHECK(parabola->Status == DerivativeStatus::Valid);

    CurveSurfacePoint bad = at.Points.front();
    bad.Parameter = std::numeric_limits<double>::quiet_NaN();
    const auto refused = Stereon::Linearize(UpHyperbola, Surface{shelf}, bad);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().GetCode() == ErrorCode::InvalidInput);
    CHECK(Stereon::Linearize(UpParabola, Surface{}, at.Points.front()).error().GetCode() ==
          ErrorCode::InvalidInput);
}

TEST_CASE("an open conic's point without a derivative says so and carries none",
          "[intersect][linearization][open-conic]")
{
    ParameterSet set;
    const ParameterId p = *set.Add("p", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0}).has_value());
    const std::vector<ConicBinding> centre{{p, ConicComponent::CenterZ}};
    const auto status = [&](const auto& conic, const Plane& plane,
                            const CurveSurfaceIntersection& at) {
        const auto jvp = Stereon::IntersectJvp(conic, centre, Surface{plane}, {}, at, sensitivity);
        INFO(IntersectTest::ErrorText(jvp));
        REQUIRE(jvp.has_value());
        CHECK(jvp->GetValues().has_value() == HasDerivative(jvp->GetStatus()));
        REQUIRE(sensitivity.SetAdjoint(std::vector<double>(Values * at.Points.size(), 1.0))
                    .has_value());
        const auto vjp = Stereon::IntersectVjp(conic, centre, Surface{plane}, {}, at, sensitivity);
        REQUIRE(vjp.has_value());
        CHECK(vjp->GetStatus() == jvp->GetStatus());
        CHECK(vjp->GetValues().has_value() == HasDerivative(vjp->GetStatus()));
        return jvp->GetStatus();
    };
    const Plane shelf = ShelfAt(3.0);
    SECTION("a touch the solver reported: Nondifferentiable")
    {
        CurveSurfaceIntersection at = On(UpHyperbola, shelf);
        at.Points.front().Contact = ContactType::Tangential;
        CHECK(status(UpHyperbola, shelf, at) == DerivativeStatus::Nondifferentiable);
        CurveSurfaceIntersection onParabola = On(UpParabola, ShelfAt(1.0));
        onParabola.Points.front().Contact = ContactType::Tangential;
        CHECK(status(UpParabola, ShelfAt(1.0), onParabola) == DerivativeStatus::Nondifferentiable);
    }
    SECTION("a point the solver left undecided: Indeterminate")
    {
        CurveSurfaceIntersection at = On(UpHyperbola, shelf);
        at.Points.front().Decided = false;
        CHECK(status(UpHyperbola, shelf, at) == DerivativeStatus::Indeterminate);
        CurveSurfaceIntersection onParabola = On(UpParabola, ShelfAt(1.0));
        onParabola.Points.front().Decided = false;
        CHECK(status(UpParabola, ShelfAt(1.0), onParabola) == DerivativeStatus::Indeterminate);
    }
    SECTION("the curve in the plane: an overlap, Nondifferentiable")
    {
        const Plane own{Point3{}, ex, ez};
        const CurveSurfaceIntersection hyperbola = On(UpHyperbola, own);
        REQUIRE_FALSE(hyperbola.Overlaps.empty());
        CHECK(status(UpHyperbola, own, hyperbola) == DerivativeStatus::Nondifferentiable);
        const CurveSurfaceIntersection parabola = On(UpParabola, own);
        REQUIRE_FALSE(parabola.Overlaps.empty());
        CHECK(status(UpParabola, own, parabola) == DerivativeStatus::Nondifferentiable);
    }
    SECTION("a grazing crossing: NearEvent, qualified")
    {
        // A plane across x through the hyperbola's point a hair off its
        // vertex, where the curve runs within 3 sinh t / (4 cosh t) = 7.5e-8
        // of the plane's direction; the parabola likewise, t / (2 f) = 2e-7.
        for (const double t : {1e-7, -1e-7}) {
            const Point3 onHyperbola = Eval(UpHyperbola, t);
            const Plane wall{onHyperbola, ey, ez};
            const CurveSurfaceIntersection at{
                {CurveSurfacePoint{t, Point2{}, onHyperbola, ContactType::Crossing, true}}, {}};
            CHECK(status(UpHyperbola, wall, at) == DerivativeStatus::NearEvent);
            const auto lin = Stereon::Linearize(UpHyperbola, Surface{wall}, at.Points.front());
            REQUIRE(lin.has_value());
            CHECK(lin->Status == DerivativeStatus::NearEvent);
            CHECK(lin->Conditioning < Stereon::NearEventSine);
            const Point3 onParabola = Eval(UpParabola, t);
            const Plane other{onParabola, ey, ez};
            const CurveSurfaceIntersection grazing{
                {CurveSurfacePoint{t, Point2{}, onParabola, ContactType::Crossing, true}}, {}};
            CHECK(status(UpParabola, other, grazing) == DerivativeStatus::NearEvent);
        }
    }
    SECTION("the vertex on the plane across it: singular, Indeterminate")
    {
        // At t = 0 the hyperbola's tangent is along z, in the plane x = 3: a
        // crossing claimed there has a singular system, whatever it claims.
        const Plane wall{Point3{3.0, 0.0, 0.0}, ey, ez};
        const CurveSurfaceIntersection at{
            {CurveSurfacePoint{0.0, Point2{}, Point3{3.0, 0.0, 0.0}, ContactType::Crossing, true}},
            {}};
        CHECK(status(UpHyperbola, wall, at) == DerivativeStatus::Indeterminate);
    }
}

TEST_CASE("an open conic's bindings are refused where it has no such value",
          "[intersect][linearization][open-conic]")
{
    ParameterSet set;
    const ParameterId length = *set.Add("p", ParameterQuantity::Length, 1.0);
    const ParameterId angle = *set.Add("a", ParameterQuantity::Angle, 1.0);
    ParameterSet other;
    const ParameterId foreign = *other.Add("q", ParameterQuantity::Length, 1.0);
    SensitivityContext sensitivity(std::move(set));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    const Plane shelf = ShelfAt(1.0);
    const CurveSurfaceIntersection onHyperbola = On(UpHyperbola, shelf);
    const CurveSurfaceIntersection onParabola = On(UpParabola, shelf);
    const auto refused = [&](const auto& conic, const CurveSurfaceIntersection& at,
                             const ConicBinding binding) {
        const std::vector<ConicBinding> bindings{binding};
        const auto jvp =
            Stereon::IntersectJvp(conic, bindings, Surface{shelf}, {}, at, sensitivity);
        const auto vjp =
            Stereon::IntersectVjp(conic, bindings, Surface{shelf}, {}, at, sensitivity);
        return !jvp && jvp.error().GetCode() == ErrorCode::InvalidInput && !vjp &&
               vjp.error().GetCode() == ErrorCode::InvalidInput;
    };
    // A hyperbola has two radii, no single one and no focal length.
    CHECK(refused(UpHyperbola, onHyperbola, {length, ConicComponent::Radius}));
    CHECK(refused(UpHyperbola, onHyperbola, {length, ConicComponent::FocalLength}));
    // A parabola has a focal length and no radius.
    CHECK(refused(UpParabola, onParabola, {length, ConicComponent::Radius}));
    CHECK(refused(UpParabola, onParabola, {length, ConicComponent::MajorRadius}));
    CHECK(refused(UpParabola, onParabola, {length, ConicComponent::MinorRadius}));
    // A radius or a focal length is a Length, a turn an Angle.
    CHECK(refused(UpHyperbola, onHyperbola, {angle, ConicComponent::MajorRadius}));
    CHECK(refused(UpParabola, onParabola, {angle, ConicComponent::FocalLength}));
    CHECK(refused(UpHyperbola, onHyperbola, {length, ConicComponent::RotationZ}));
    CHECK(refused(UpParabola, onParabola, {length, ConicComponent::RotationX}));
    // A parameter of another set.
    CHECK(refused(UpParabola, onParabola, {foreign, ConicComponent::FocalLength}));

    // A circle and an ellipse have no focal length.
    const Circle circle{Point3{}, 2.0, ex, ez};
    const Ellipse ellipse{Point3{}, 3.0, 2.0, ex, ez};
    const std::vector<ConicBinding> focal{{length, ConicComponent::FocalLength}};
    const auto onCircle = Stereon::Intersect(Stereon::Curve{circle}, Surface{shelf}, ctx);
    const auto onEllipse = Stereon::Intersect(Stereon::Curve{ellipse}, Surface{shelf}, ctx);
    REQUIRE(onCircle.has_value());
    REQUIRE(onEllipse.has_value());
    CHECK(Stereon::IntersectJvp(circle, focal, Surface{shelf}, {}, *onCircle, sensitivity)
              .error()
              .GetCode() == ErrorCode::InvalidInput);
    CHECK(Stereon::IntersectJvp(ellipse, focal, Surface{shelf}, {}, *onEllipse, sensitivity)
              .error()
              .GetCode() == ErrorCode::InvalidInput);

    // A result of other geometry, and an adjoint of another length.
    const std::vector<ConicBinding> good{{length, ConicComponent::FocalLength}};
    const auto stale = Stereon::IntersectJvp(UpParabola, good, Surface{ShelfAt(1.25)}, {},
                                             onParabola, sensitivity);
    REQUIRE_FALSE(stale.has_value());
    CHECK(stale.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(sensitivity.SetAdjoint(std::vector<double>(Values + 1, 1.0)).has_value());
    const auto wrong =
        Stereon::IntersectVjp(UpParabola, good, Surface{shelf}, {}, onParabola, sensitivity);
    REQUIRE_FALSE(wrong.has_value());
    CHECK(wrong.error().GetCode() == ErrorCode::InvalidInput);
    // A surface no binding reaches is NotSupported, a status.
    const Surface swept{Stereon::ExtrusionSurface(Stereon::Curve{circle}, ey)};
    const auto unsupported =
        Stereon::IntersectJvp(UpParabola, good, swept, {}, CurveSurfaceIntersection{}, sensitivity);
    REQUIRE(unsupported.has_value());
    CHECK(unsupported->GetStatus() == DerivativeStatus::NotSupported);
    CHECK_FALSE(unsupported->GetValues().has_value());
}
