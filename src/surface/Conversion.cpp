// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Surface;

import :BSplineSurface;
import :Cone;
import :Conversion;
import :Cylinder;
import :Erased;
import :Plane;
import :Sphere;
import :Torus;
import :Visit;

namespace Stereon {

namespace {

/// InvalidInput unless the range is finite and non-empty in one parameter
/// and fits the domain, or one period when the surface is periodic in it.
[[nodiscard]] Result<void> CheckRange(const Interval range, const Interval domain,
                                      const double period, const std::string_view where)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    const bool finite = std::isfinite(low) && std::isfinite(high) && low < high;
    const bool fits =
        period > 0.0 ? high - low <= period : Contains(domain, low) && Contains(domain, high);
    if (!finite || !fits) {
        std::string message{where};
        message += ": the range must be finite, non-empty and inside the domain, or at most "
                   "one period of a periodic parameter";
        return Fail(ErrorCode::InvalidInput, std::move(message));
    }
    return {};
}

/// Periods(s) for a periodic surface, (0, 0) for the plane.
template <typename S>
[[nodiscard]] Vec2 PeriodsOf(const S& surface)
{
    if constexpr (PeriodicSurface<S>) {
        return Periods(surface);
    } else {
        static_cast<void>(surface);
        return {};
    }
}

template <typename S>
[[nodiscard]] Result<void> CheckBox(const S& surface, const Box2& range,
                                    const std::string_view where)
{
    const Box2 domain = Domain(surface);
    const Vec2 periods = PeriodsOf(surface);
    STEREON_TRY_VOID(CheckRange(range.X, domain.X, periods.X, where));
    return CheckRange(range.Y, domain.Y, periods.Y, where);
}

/// A meridian: a curve in the half-plane through the axis, with coordinates
/// (rho, h), the distance from the axis and the height along it.
struct Meridian {
    int Degree = 1;
    std::vector<double> Knots;
    std::vector<Point2> Poles;
    std::vector<double> Weights;
};

/// The meridian swept around the axis z through the angles of `around`
/// (algorithm A8.1): pole (i, j) is origin + rho_j U_i + h_j z, where U_i
/// are the poles of the unit circle arc in the (x, y) plane, and its weight
/// the product of the two curves' weights. The tensor product separates,
/// so the surface is origin + rho(v) e(u) + h(v) z, a surface of
/// revolution, including where rho is negative.
[[nodiscard]] BSplineSurface Revolve(const Point3& origin, const Vec3& x, const Vec3& y,
                                     const Vec3& z, const Meridian& meridian, const Interval around)
{
    const Detail::UnitArc arc = Detail::MakeUnitArc(around.GetLower(), around.GetUpper());
    BSplineSurface surface{2, meridian.Degree, arc.Knots, meridian.Knots, {}, {}};
    for (std::size_t i = 0; i < arc.Poles.size(); ++i) {
        const Vec3 radial = arc.Poles[i].X * x + arc.Poles[i].Y * y;
        for (std::size_t j = 0; j < meridian.Poles.size(); ++j) {
            surface.Poles.push_back(origin + meridian.Poles[j].X * radial +
                                    meridian.Poles[j].Y * z);
            surface.Weights.push_back(arc.Weights[i] * meridian.Weights[j]);
        }
    }
    return surface;
}

/// A meridian line from (rho0, h0) at v0 to (rho1, h1) at v1.
[[nodiscard]] Meridian LineMeridian(const Interval v, const Point2 start, const Point2 end)
{
    return {1, {v.GetLower(), v.GetLower(), v.GetUpper(), v.GetUpper()}, {start, end}, {1.0, 1.0}};
}

/// A meridian circle arc of the given radius about (centreRho, 0), over the
/// angles of v.
[[nodiscard]] Meridian ArcMeridian(const Interval v, const double centreRho, const double radius)
{
    const Detail::UnitArc arc = Detail::MakeUnitArc(v.GetLower(), v.GetUpper());
    Meridian meridian{2, arc.Knots, {}, arc.Weights};
    for (const Point2& p : arc.Poles) {
        meridian.Poles.push_back({centreRho + radius * p.X, radius * p.Y});
    }
    return meridian;
}

} // namespace

Result<BSplineSurface> ToBSpline(const Plane& plane, const Box2 range)
{
    STEREON_TRY_VOID(CheckBox(plane, range, "ToBSpline(Plane)"));
    const double u0 = range.X.GetLower();
    const double u1 = range.X.GetUpper();
    const double v0 = range.Y.GetLower();
    const double v1 = range.Y.GetUpper();
    return BSplineSurface{1,
                          1,
                          {u0, u0, u1, u1},
                          {v0, v0, v1, v1},
                          {Eval(plane, Point2{u0, v0}), Eval(plane, Point2{u0, v1}),
                           Eval(plane, Point2{u1, v0}), Eval(plane, Point2{u1, v1})},
                          {}};
}

Result<BSplineSurface> ToBSpline(const Cylinder& cylinder, const Box2 range)
{
    STEREON_TRY_VOID(CheckBox(cylinder, range, "ToBSpline(Cylinder)"));
    const Meridian meridian = LineMeridian(range.Y, Point2{cylinder.Radius, range.Y.GetLower()},
                                           Point2{cylinder.Radius, range.Y.GetUpper()});
    return Revolve(cylinder.Origin, cylinder.XAxis, cylinder.YAxis, Axis(cylinder), meridian,
                   range.X);
}

Result<BSplineSurface> ToBSpline(const Cone& cone, const Box2 range)
{
    STEREON_TRY_VOID(CheckBox(cone, range, "ToBSpline(Cone)"));
    const double s = Sin(cone.SemiAngle);
    const double c = Cos(cone.SemiAngle);
    const double v0 = range.Y.GetLower();
    const double v1 = range.Y.GetUpper();
    const Meridian meridian = LineMeridian(range.Y, Point2{cone.Radius + v0 * s, v0 * c},
                                           Point2{cone.Radius + v1 * s, v1 * c});
    return Revolve(cone.Origin, cone.XAxis, cone.YAxis, Axis(cone), meridian, range.X);
}

Result<BSplineSurface> ToBSpline(const Sphere& sphere, const Box2 range)
{
    STEREON_TRY_VOID(CheckBox(sphere, range, "ToBSpline(Sphere)"));
    return Revolve(sphere.Center, sphere.XAxis, sphere.YAxis, Axis(sphere),
                   ArcMeridian(range.Y, 0.0, sphere.Radius), range.X);
}

Result<BSplineSurface> ToBSpline(const Torus& torus, const Box2 range)
{
    STEREON_TRY_VOID(CheckBox(torus, range, "ToBSpline(Torus)"));
    return Revolve(torus.Center, torus.XAxis, torus.YAxis, Axis(torus),
                   ArcMeridian(range.Y, torus.MajorRadius, torus.MinorRadius), range.X);
}

Result<BSplineSurface> ToBSpline(const Surface& surface, const Box2 range)
{
    return VisitBuiltin(surface, [&]<typename G>(const G& geometry) -> Result<BSplineSurface> {
        if constexpr (std::is_same_v<G, Surface>) {
            return Fail(ErrorCode::NotSupported,
                        "ToBSpline: the surface's type has no exact NURBS form the kernel knows");
        } else {
            return ToBSpline(geometry, range);
        }
    });
}

} // namespace Stereon
