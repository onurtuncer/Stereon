// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md): every module
// implementation unit of a geometry library includes it, because one that
// calls core's inline functions with contracts and parses none itself
// leaves __tu_has_violation undefined, which the macOS linker rejects.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :Hyperbola;

/// Projection onto the open conics and their checked factories. Both
/// projections reduce the stationarity condition of the squared distance to
/// a polynomial and bound where the closest point can be: it is no farther
/// from p than the vertex, which confines its parameter to a finite range.

namespace Stereon {

namespace {

/// The best parameter so far and its squared in-plane distance.
struct Closest {
    double Parameter;
    double DistanceSquared;
};

void Keep(Closest& best, const double t, const double d)
{
    if (d < best.DistanceSquared || (d == best.DistanceSquared && t < best.Parameter)) {
        best = {t, d};
    }
}

template <std::size_t D>
[[nodiscard]] double ProjectHyperbola(const BasicHyperbola<double, D>& h,
                                      const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(Hyperbola)", "a non-finite point");
    const BasicVecN<double, D> offset = p - h.Center;
    const double x = Dot(offset, h.XAxis);
    const double y = Dot(offset, h.YAxis);
    const double a = h.MajorRadius;
    const double b = h.MinorRadius;
    const auto distance = [&](const double t) {
        const double dx = a * Cosh(t) - x;
        const double dy = b * Sinh(t) - y;
        return dx * dx + dy * dy;
    };
    Closest best{0.0, distance(0.0)};
    // The closest point is no farther than the vertex, so its y coordinate
    // b sinh t is within that distance of y.
    const double limit = Asinh((std::abs(y) + std::sqrt(best.DistanceSquared)) / b);
    const double sum = a * a + b * b;
    const std::array<double, 5> quartic{-sum, 2.0 * (a * x - b * y), 0.0, -2.0 * (a * x + b * y),
                                        sum};
    for (const double s : Detail::PolynomialRoots(quartic, Interval{Exp(-limit), Exp(limit)})) {
        const double t = Log(s);
        Keep(best, t, distance(t));
    }
    return best.Parameter;
}

template <std::size_t D>
[[nodiscard]] double ProjectParabola(const BasicParabola<double, D>& c,
                                     const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(Parabola)", "a non-finite point");
    const BasicVecN<double, D> offset = p - c.Vertex;
    const double x = Dot(offset, c.XAxis);
    const double y = Dot(offset, c.YAxis);
    const double f = c.FocalLength;
    const auto distance = [&](const double t) {
        const double dx = t * t / (4.0 * f) - x;
        const double dy = t - y;
        return dx * dx + dy * dy;
    };
    Closest best{0.0, distance(0.0)};
    if (best.DistanceSquared == 0.0) {
        return 0.0; // p is the vertex
    }
    // The closest point is no farther than the vertex, and its y coordinate
    // is t itself.
    const double reach = std::sqrt(best.DistanceSquared);
    const std::array<double, 4> cubic{-y, 1.0 - x / (2.0 * f), 0.0, 1.0 / (8.0 * f * f)};
    for (const double t : Detail::PolynomialRoots(cubic, Interval{y - reach, y + reach})) {
        Keep(best, t, distance(t));
    }
    return best.Parameter;
}

template <std::size_t D>
[[nodiscard]] Result<BasicHyperbola<double, D>>
MakeHyperbolaN(const BasicPointN<double, D>& center, const double majorRadius,
               const double minorRadius, const BasicVecN<double, D>& xAxis,
               const BasicVecN<double, D>& yAxis)
{
    constexpr std::string_view where = "MakeHyperbola";
    STEREON_TRY_VOID(Detail::CheckFinite(center, "Center", where));
    STEREON_TRY_VOID(Detail::CheckPositive(majorRadius, "MajorRadius", where));
    STEREON_TRY_VOID(Detail::CheckPositive(minorRadius, "MinorRadius", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return BasicHyperbola<double, D>{center, majorRadius, minorRadius, xAxis, yAxis};
}

template <std::size_t D>
[[nodiscard]] Result<BasicParabola<double, D>>
MakeParabolaN(const BasicPointN<double, D>& vertex, const double focalLength,
              const BasicVecN<double, D>& xAxis, const BasicVecN<double, D>& yAxis)
{
    constexpr std::string_view where = "MakeParabola";
    STEREON_TRY_VOID(Detail::CheckFinite(vertex, "Vertex", where));
    STEREON_TRY_VOID(Detail::CheckPositive(focalLength, "FocalLength", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return BasicParabola<double, D>{vertex, focalLength, xAxis, yAxis};
}

} // namespace

double Project(const Hyperbola& hyperbola, const Point3& p)
{
    return ProjectHyperbola(hyperbola, p);
}

double Project(const Hyperbola2d& hyperbola, const Point2& p)
{
    return ProjectHyperbola(hyperbola, p);
}

double Project(const Parabola& parabola, const Point3& p)
{
    return ProjectParabola(parabola, p);
}

double Project(const Parabola2d& parabola, const Point2& p)
{
    return ProjectParabola(parabola, p);
}

Result<Hyperbola> MakeHyperbola(const Point3& center, const double majorRadius,
                                const double minorRadius, const Vec3& xAxis, const Vec3& yAxis)
{
    return MakeHyperbolaN<3>(center, majorRadius, minorRadius, xAxis, yAxis);
}

Result<Hyperbola2d> MakeHyperbola(const Point2& center, const double majorRadius,
                                  const double minorRadius, const Vec2& xAxis, const Vec2& yAxis)
{
    return MakeHyperbolaN<2>(center, majorRadius, minorRadius, xAxis, yAxis);
}

Result<Parabola> MakeParabola(const Point3& vertex, const double focalLength, const Vec3& xAxis,
                              const Vec3& yAxis)
{
    return MakeParabolaN<3>(vertex, focalLength, xAxis, yAxis);
}

Result<Parabola2d> MakeParabola(const Point2& vertex, const double focalLength, const Vec2& xAxis,
                                const Vec2& yAxis)
{
    return MakeParabolaN<2>(vertex, focalLength, xAxis, yAxis);
}

} // namespace Stereon
