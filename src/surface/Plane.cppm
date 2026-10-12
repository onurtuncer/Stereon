// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

export module Stereon.Surface:Plane;

import Stereon.Curve;
import :Concepts;

/// The plane.

export namespace Stereon {

/// Origin + u * XAxis + v * YAxis for every real u and v. XAxis and YAxis are
/// orthonormal; the normal is Cross(XAxis, YAxis).
template <Scalar T>
struct BasicPlane {
    BasicPoint3<T> Origin;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;

    constexpr bool operator==(const BasicPlane&) const = default;
};

using Plane = BasicPlane<double>;

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicPlane<T>& plane, const BasicPoint2<T>& uv)
{
    return plane.Origin + uv.X * plane.XAxis + uv.Y * plane.YAxis;
}

/// XAxis for (1, 0), YAxis for (0, 1), zero for every higher order.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicPlane<T>& plane,
                                           [[maybe_unused]] const BasicPoint2<T>& uv, const int du,
                                           const int dv)
    STEREON_HOT_PRE(du >= 0 && dv >= 0 && du + dv >= 1)
{
    if (du == 1 && dv == 0) {
        return plane.XAxis;
    }
    if (du == 0 && dv == 1) {
        return plane.YAxis;
    }
    return {};
}

/// The whole (u, v) plane.
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicPlane<T>& plane)
{
    return Box2::Everything();
}

/// Everything: a plane is unbounded.
[[nodiscard]] constexpr Box3 BoundingBox([[maybe_unused]] const Plane& plane)
{
    return Box3::Everything();
}

/// The parameter of the foot of the perpendicular from p.
[[nodiscard]] constexpr Point2 Project(const Plane& plane, const Point3& p)
{
    const Vec3 offset = p - plane.Origin;
    return {Dot(offset, plane.XAxis), Dot(offset, plane.YAxis)};
}

/// The same plane with the opposite normal:
/// Eval(Reversed(s), (u, v)) == Eval(s, (u, -v)).
template <Scalar T>
[[nodiscard]] constexpr BasicPlane<T> Reversed(const BasicPlane<T>& plane)
{
    return {plane.Origin, plane.XAxis, -plane.YAxis};
}

/// The image of the plane under x, which must be a similarity; the axes are
/// renormalised, so the parameter is rescaled by the scale factor.
template <Scalar T>
[[nodiscard]] constexpr BasicPlane<T> Transformed(const BasicPlane<T>& plane,
                                                  const BasicTransform<T>& x)
{
    return {Transformed(plane.Origin, x), Normalized(Transformed(plane.XAxis, x)),
            Normalized(Transformed(plane.YAxis, x))};
}

static_assert(SurfaceGeometry<Plane>);

} // namespace Stereon
