// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

export module Stereon.Surface:Cylinder;

import Stereon.Curve;
import :Concepts;

/// The circular cylinder.

export namespace Stereon {

/// Origin + Radius * (cos u * XAxis + sin u * YAxis) + v * Axis(c), with u in
/// [0, 2 pi] and v real. XAxis and YAxis are orthonormal. The axis is
/// Cross(XAxis, YAxis) when RightHanded and its negative otherwise, so the
/// frame may have either handedness: the normal Su x Sv points away from the
/// axis line in a right-handed frame and towards it in a left-handed one,
/// which is how Reversed() turns the normal round. A flag rather than a third
/// axis keeps the type within the wrappers' buffer (ADR-0003). Radius is
/// positive.
template <Scalar T>
struct BasicCylinder {
    BasicPoint3<T> Origin;
    T Radius;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;
    bool RightHanded = true;

    constexpr bool operator==(const BasicCylinder&) const = default;
};

using Cylinder = BasicCylinder<double>;

/// The direction of the cylinder's axis: the direction of growing v.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Axis(const BasicCylinder<T>& cylinder)
{
    const BasicVec3<T> z = Cross(cylinder.XAxis, cylinder.YAxis);
    return cylinder.RightHanded ? z : -z;
}

/// Radius * (cos u * XAxis + sin u * YAxis) turned by `quarterTurns` quarter
/// turns: the u-derivatives of the cylinder, by n mod 4, so that no multiple
/// of pi/2 is rounded into the argument.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> RadialTurned(const BasicCylinder<T>& cylinder, const T& u,
                                                  const int quarterTurns)
{
    // Not std::cos: see Stereon.Core:Trigonometry.
    const T c = Cos(u);
    const T s = Sin(u);
    switch (quarterTurns % 4) {
    case 1:
        return cylinder.Radius * (-s * cylinder.XAxis + c * cylinder.YAxis);
    case 2:
        return cylinder.Radius * (-c * cylinder.XAxis - s * cylinder.YAxis);
    case 3:
        return cylinder.Radius * (s * cylinder.XAxis - c * cylinder.YAxis);
    default:
        return cylinder.Radius * (c * cylinder.XAxis + s * cylinder.YAxis);
    }
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicCylinder<T>& cylinder,
                                            const BasicPoint2<T>& uv)
{
    return cylinder.Origin + RadialTurned(cylinder, uv.X, 0) + uv.Y * Axis(cylinder);
}

/// Derivatives in u turn the radius vector; the only non-zero derivative
/// involving v is the axis, for (0, 1).
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicCylinder<T>& cylinder,
                                           const BasicPoint2<T>& uv, const int du, const int dv)
    STEREON_HOT_PRE(du >= 0 && dv >= 0 && du + dv >= 1)
{
    if (dv == 0) {
        return RadialTurned(cylinder, uv.X, du);
    }
    if (du == 0 && dv == 1) {
        return Axis(cylinder);
    }
    return {};
}

/// u in [0, 2 pi], v unbounded.
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicCylinder<T>& cylinder)
{
    return {Interval{0.0, CirclePeriod}, Interval::Everything()};
}

/// Everything: the cylinder is unbounded along its axis.
[[nodiscard]] constexpr Box3 BoundingBox([[maybe_unused]] const Cylinder& cylinder)
{
    return Box3::Everything();
}

/// u from the angle around the axis, in [0, 2 pi), v from the height along
/// it. A point on the axis is equally close to every section point; it gets
/// u = 0.
[[nodiscard]] inline Point2 Project(const Cylinder& cylinder, const Point3& p)
{
    const Vec3 offset = p - cylinder.Origin;
    const double v = Dot(offset, Axis(cylinder));
    const double x = Dot(offset, cylinder.XAxis);
    const double y = Dot(offset, cylinder.YAxis);
    if (x == 0.0 && y == 0.0) {
        return {0.0, v};
    }
    const double angle = Atan2(y, x);
    if (angle >= 0.0) {
        return {angle, v};
    }
    // As for the circle: a tiny negative angle plus 2 pi rounds to 2 pi.
    const double wrapped = angle + CirclePeriod;
    return {wrapped < CirclePeriod ? wrapped : 0.0, v};
}

/// The same cylinder with the opposite normal, by running u the other way:
/// Eval(Reversed(s), (u, v)) == Eval(s, (2 pi - u, v)).
template <Scalar T>
[[nodiscard]] constexpr BasicCylinder<T> Reversed(const BasicCylinder<T>& cylinder)
{
    // Negating YAxis negates Cross(XAxis, YAxis); flipping the handedness
    // keeps the axis, so only u changes direction.
    return {cylinder.Origin, cylinder.Radius, cylinder.XAxis, -cylinder.YAxis,
            !cylinder.RightHanded};
}

/// The image of the cylinder under x, which must be a similarity. The radius
/// scales and the axes are renormalised, so v is rescaled by the scale
/// factor and u is kept.
template <Scalar T>
[[nodiscard]] constexpr BasicCylinder<T> Transformed(const BasicCylinder<T>& cylinder,
                                                     const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(cylinder.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(cylinder.YAxis, x);
    const T scale = Length(xAxis);
    // A reflection reverses Cross(XAxis, YAxis); the flag follows, so the
    // axis is the image of the old one.
    const bool mirrored = Determinant(x.Linear) < T{};
    return {Transformed(cylinder.Origin, x), cylinder.Radius * scale, xAxis / scale, yAxis / scale,
            cylinder.RightHanded != mirrored};
}

static_assert(SurfaceGeometry<Cylinder>);
static_assert(Detail::FitsInline<Cylinder, GeometryBufferSize>);

} // namespace Stereon
