// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <numbers>

export module Stereon.Surface:Sphere;

import Stereon.Curve;
import :Concepts;

/// The sphere.

export namespace Stereon {

/// Center + Radius * (cos v * (cos u * XAxis + sin u * YAxis) + sin v * Axis(s))
/// with u in [0, 2 pi], the longitude, and v in [-pi/2, pi/2], the latitude.
/// XAxis and YAxis are orthonormal; the axis is Cross(XAxis, YAxis) when
/// RightHanded and its negative otherwise, as for BasicCylinder, so the
/// normal Su x Sv points outward in a right-handed frame and inward in a
/// left-handed one. Radius is positive.
///
/// The parametrisation is singular at the poles v = +-pi/2, where every u
/// gives the same point and Su vanishes: the normal there is the limit
/// +-Axis(s), not Su x Sv. These are the Degenerate edges of ADR-0009.
/// The conventions are STEP's spherical_surface and OCCT's
/// Geom_SphericalSurface.
template <Scalar T>
struct BasicSphere {
    BasicPoint3<T> Center;
    T Radius;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;
    bool RightHanded = true;

    constexpr bool operator==(const BasicSphere&) const = default;
};

using Sphere = BasicSphere<double>;

/// The latitude of the north pole, pi/2; the south pole is at its negative.
inline constexpr double SpherePoleLatitude = std::numbers::pi / 2;

/// The sphere after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, Radius > 0 and the axes are
/// orthonormal within FrameAllowance.
[[nodiscard]] Result<Sphere> MakeSphere(const Point3& center, const double radius,
                                        const Vec3& xAxis, const Vec3& yAxis,
                                        const bool rightHanded = true);

/// The direction of the axis, from the south pole to the north pole.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Axis(const BasicSphere<T>& sphere)
{
    return FrameAxis(sphere.XAxis, sphere.YAxis, sphere.RightHanded);
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicSphere<T>& sphere, const BasicPoint2<T>& uv)
{
    Detail::RequireFinite(uv, "Eval(Sphere)", "a non-finite parameter");
    const BasicVec3<T> radial = RadialDirection(sphere.XAxis, sphere.YAxis, uv.X, 0);
    return sphere.Center + sphere.Radius * (Cos(uv.Y) * radial + Sin(uv.Y) * Axis(sphere));
}

/// The point separates into a u part and a v part, R * cos v * radial(u) +
/// R * sin v * axis, so a mixed derivative is the product of derivatives of
/// the factors: R * cos(v + dv pi/2) * radial^(du)(u), plus
/// R * sin(v + dv pi/2) * axis when du is 0.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicSphere<T>& sphere, const BasicPoint2<T>& uv,
                                           const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(Sphere)", "a non-finite parameter");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(Sphere)",
                    "a derivative order that does not exist");
    const BasicCosSin<T> latitude = QuarterTurned(Cos(uv.Y), Sin(uv.Y), dv);
    const BasicVec3<T> radial = RadialDirection(sphere.XAxis, sphere.YAxis, uv.X, du);
    const BasicVec3<T> around = (sphere.Radius * latitude.Cos) * radial;
    if (du != 0) {
        return around;
    }
    return around + (sphere.Radius * latitude.Sin) * Axis(sphere);
}

/// The periods in u and v: 2 pi in u, none (0) in v.
template <Scalar T>
[[nodiscard]] constexpr Vec2 Periods([[maybe_unused]] const BasicSphere<T>& sphere) noexcept
{
    return {CirclePeriod, 0.0};
}

/// u in [0, 2 pi], v in [-pi/2, pi/2].
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicSphere<T>& sphere)
{
    return {Interval{0.0, CirclePeriod}, Interval{-SpherePoleLatitude, SpherePoleLatitude}};
}

/// The centre plus or minus the radius along every axis, in interval
/// arithmetic: tight, and it encloses the exact sphere.
[[nodiscard]] constexpr Box3 BoundingBox(const Sphere& sphere)
{
    const Interval reach{-sphere.Radius, sphere.Radius};
    return {Interval{sphere.Center.X} + reach, Interval{sphere.Center.Y} + reach,
            Interval{sphere.Center.Z} + reach};
}

/// The longitude and latitude of p seen from the centre: the closest point
/// of the sphere lies on the ray from the centre through p. u is in
/// [0, 2 pi). Points on the axis get u = 0, the centre gets (0, 0); every
/// point of the sphere is equally close to it.
[[nodiscard]] inline Point2 Project(const Sphere& sphere, const Point3& p)
{
    Detail::RequireFinite(p, "Project(Sphere)", "a non-finite point");
    const Vec3 offset = p - sphere.Center;
    const double x = Dot(offset, sphere.XAxis);
    const double y = Dot(offset, sphere.YAxis);
    const double z = Dot(offset, Axis(sphere));
    return {WrappedAngle(y, x), Atan2(z, Hypot(x, y))};
}

/// The same sphere with the opposite normal, by running u the other way:
/// Eval(Reversed(s), (u, v)) == Eval(s, (2 pi - u, v)). As for the
/// cylinder, YAxis is negated and the handedness flipped, which keeps the
/// axis.
template <Scalar T>
[[nodiscard]] constexpr BasicSphere<T> Reversed(const BasicSphere<T>& sphere)
{
    return {sphere.Center, sphere.Radius, sphere.XAxis, -sphere.YAxis, !sphere.RightHanded};
}

/// The image of the sphere under x, which must be a similarity; the
/// parameters are kept and the radius scales.
template <Scalar T>
[[nodiscard]] constexpr BasicSphere<T> Transformed(const BasicSphere<T>& sphere,
                                                   const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(sphere.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(sphere.YAxis, x);
    const T scale = Length(xAxis);
    // A reflection reverses Cross(XAxis, YAxis); the flag follows, so the
    // axis is the image of the old one.
    const bool mirrored = Determinant(x.Linear) < T{};
    return {Transformed(sphere.Center, x), sphere.Radius * scale, xAxis / scale, yAxis / scale,
            sphere.RightHanded != mirrored};
}

static_assert(SurfaceGeometry<Sphere>);
static_assert(Detail::FitsInline<Sphere, GeometryBufferSize>);

} // namespace Stereon
