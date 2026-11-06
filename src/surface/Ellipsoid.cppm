// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <cstddef>
#include <numbers>

export module Stereon.Surface:Ellipsoid;

import Stereon.Curve;
import :Concepts;

/// The ellipsoid (ADR-0048).

export namespace Stereon {

/// Center + cos v * (A * cos u * XAxis + B * sin u * YAxis) + C * sin v * Axis(e)
/// with u in [0, 2 pi], the longitude, and v in [-pi/2, pi/2], the latitude:
/// the sphere's chart with a semi-axis along each axis of the frame. XAxis
/// and YAxis are orthonormal; the axis is Cross(XAxis, YAxis) when
/// RightHanded and its negative otherwise, as for BasicSphere, so the
/// normal Su x Sv points outward in a right-handed frame and inward in a
/// left-handed one. A, B and C are positive; equal semi-axes are an
/// ellipsoid still, never converted to a sphere.
///
/// The parametrisation is singular at the poles v = +-pi/2, where every u
/// gives the same point and Su vanishes: the normal there is the limit
/// +-Axis(e). These are the Degenerate edges of ADR-0009.
template <Scalar T>
struct BasicEllipsoid {
    BasicPoint3<T> Center;
    T A;
    T B;
    T C;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;
    bool RightHanded = true;

    constexpr bool operator==(const BasicEllipsoid&) const = default;
};

using Ellipsoid = BasicEllipsoid<double>;

/// The ellipsoid after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, A, B and C are positive and
/// the axes are orthonormal within FrameAllowance.
[[nodiscard]] Result<Ellipsoid> MakeEllipsoid(const Point3& center, double a, double b, double c,
                                              const Vec3& xAxis, const Vec3& yAxis,
                                              bool rightHanded = true);

/// The direction of the axis, from the south pole to the north pole.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Axis(const BasicEllipsoid<T>& ellipsoid)
{
    return FrameAxis(ellipsoid.XAxis, ellipsoid.YAxis, ellipsoid.RightHanded);
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicEllipsoid<T>& ellipsoid,
                                            const BasicPoint2<T>& uv)
{
    Detail::RequireFinite(uv, "Eval(Ellipsoid)", "a non-finite parameter");
    const BasicVec3<T> around =
        RadialDirection(ellipsoid.A * ellipsoid.XAxis, ellipsoid.B * ellipsoid.YAxis, uv.X, 0);
    return ellipsoid.Center + (Cos(uv.Y) * around + (ellipsoid.C * Sin(uv.Y)) * Axis(ellipsoid));
}

/// The point separates into a u part and a v part, cos v * (A cos u X +
/// B sin u Y) + sin v * C axis, so a mixed derivative is the product of
/// derivatives of the factors, as for the sphere.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicEllipsoid<T>& ellipsoid,
                                           const BasicPoint2<T>& uv, const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(Ellipsoid)", "a non-finite parameter");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(Ellipsoid)",
                    "a derivative order that does not exist");
    const BasicCosSin<T> latitude = QuarterTurned(Cos(uv.Y), Sin(uv.Y), dv);
    const BasicVec3<T> around =
        latitude.Cos *
        RadialDirection(ellipsoid.A * ellipsoid.XAxis, ellipsoid.B * ellipsoid.YAxis, uv.X, du);
    if (du != 0) {
        return around;
    }
    return around + (ellipsoid.C * latitude.Sin) * Axis(ellipsoid);
}

/// The periods in u and v: 2 pi in u, none (0) in v.
template <Scalar T>
[[nodiscard]] constexpr Vec2 Periods([[maybe_unused]] const BasicEllipsoid<T>& ellipsoid) noexcept
{
    return {CirclePeriod, 0.0};
}

/// u in [0, 2 pi], v in [-pi/2, pi/2].
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicEllipsoid<T>& ellipsoid)
{
    return {Interval{0.0, CirclePeriod}, Interval{-std::numbers::pi / 2, std::numbers::pi / 2}};
}

} // namespace Stereon

namespace Stereon::Detail {

/// How far the ellipsoid reaches from its centre along one coordinate axis,
/// given the components x, y and z of XAxis, YAxis and the axis along it:
/// sqrt((A x)² + (B y)² + (C z)²), in interval arithmetic. The components
/// are passed one axis at a time, not indexed, for the reason TorusReach
/// gives: an inline non-template function that instantiates core's indexed
/// accessors corrupts GCC 16.2's module file.
[[nodiscard]] constexpr Interval EllipsoidReach(const Ellipsoid& ellipsoid, const double x,
                                                const double y, const double z)
{
    const Interval squares = Squared(Interval{ellipsoid.A} * Interval{x}) +
                             Squared(Interval{ellipsoid.B} * Interval{y}) +
                             Squared(Interval{ellipsoid.C} * Interval{z});
    // Clamped as in BoundingBox(Circle): outward rounding can step the lower
    // bound of a sum of squares below zero.
    const double lower = squares.GetLower() > 0.0 ? squares.GetLower() : 0.0;
    return Interval{-1.0, 1.0} * Sqrt(Interval{lower, squares.GetUpper()});
}

} // namespace Stereon::Detail

export namespace Stereon {

/// Along each axis of space the centre plus or minus the ellipsoid's
/// support, sqrt((A X_i)² + (B Y_i)² + (C Z_i)²), reached: tight, and in
/// interval arithmetic with the axis the evaluator uses, so it encloses the
/// exact ellipsoid.
[[nodiscard]] constexpr Box3 BoundingBox(const Ellipsoid& ellipsoid)
{
    const Vec3 axis = Axis(ellipsoid);
    const Vec3& x = ellipsoid.XAxis;
    const Vec3& y = ellipsoid.YAxis;
    return {Interval{ellipsoid.Center.X} + Detail::EllipsoidReach(ellipsoid, x.X, y.X, axis.X),
            Interval{ellipsoid.Center.Y} + Detail::EllipsoidReach(ellipsoid, x.Y, y.Y, axis.Y),
            Interval{ellipsoid.Center.Z} + Detail::EllipsoidReach(ellipsoid, x.Z, y.Z, axis.Z)};
}

} // namespace Stereon

export namespace Stereon {

/// The longitude and latitude of a point of the ellipsoid closest to p:
/// in the frame, the octant problem of ClosestOnEllipsoidOctant with the
/// semi-axes taken longest first, the point reflected back to p's octant
/// (ADR-0048 §2). p on the medial set, where two or more closest points
/// meet, gets one of them, the same one every time: the one in the
/// positive half of each coordinate p is zero in. With equal semi-axes
/// that is the sphere's Project, the centre included. Out of line, as
/// the factories are (Factories.cpp).
[[nodiscard]] Point2 Project(const Ellipsoid& ellipsoid, const Point3& p);

/// The same ellipsoid with the opposite normal, by running u the other
/// way: Eval(Reversed(e), (u, v)) == Eval(e, (2 pi - u, v)). As for the
/// sphere, YAxis is negated and the handedness flipped, which keeps the
/// axis.
template <Scalar T>
[[nodiscard]] constexpr BasicEllipsoid<T> Reversed(const BasicEllipsoid<T>& ellipsoid)
{
    return {ellipsoid.Center, ellipsoid.A,      ellipsoid.B,           ellipsoid.C,
            ellipsoid.XAxis,  -ellipsoid.YAxis, !ellipsoid.RightHanded};
}

/// The image of the ellipsoid under x, which must be a similarity; the
/// parameters are kept and the semi-axes scale. An uneven scaling is not
/// a similarity and gives an ellipsoid with another parametrisation, so
/// it is not this (ADR-0048 §1).
template <Scalar T>
[[nodiscard]] constexpr BasicEllipsoid<T> Transformed(const BasicEllipsoid<T>& ellipsoid,
                                                      const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(ellipsoid.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(ellipsoid.YAxis, x);
    const T scale = Length(xAxis);
    // A reflection reverses Cross(XAxis, YAxis); the flag follows, so the
    // axis is the image of the old one.
    const bool mirrored = Determinant(x.Linear) < T{};
    return {Transformed(ellipsoid.Center, x),
            ellipsoid.A * scale,
            ellipsoid.B * scale,
            ellipsoid.C * scale,
            xAxis / scale,
            yAxis / scale,
            ellipsoid.RightHanded != mirrored};
}

static_assert(SurfaceGeometry<Ellipsoid>);
static_assert(Detail::FitsInline<Ellipsoid, SurfaceBufferSize>);

} // namespace Stereon
