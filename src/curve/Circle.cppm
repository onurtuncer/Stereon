// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

export module Stereon.Curve:Circle;

import Stereon.Robust;
import :Concepts;

/// The circle, in space and in the parameter plane.

export namespace Stereon {

/// Center + Radius * (cos t * XAxis + sin t * YAxis) for t in [0, 2 pi].
/// XAxis and YAxis are orthonormal; in space the circle's normal is
/// Cross(XAxis, YAxis), in the plane YAxis is +-Perpendicular(XAxis) and
/// decides the sense. Radius is positive.
template <Scalar T, std::size_t D>
struct BasicCircle {
    BasicPointN<T, D> Center;
    T Radius;
    BasicVecN<T, D> XAxis;
    BasicVecN<T, D> YAxis;

    constexpr bool operator==(const BasicCircle&) const = default;
};

using Circle = BasicCircle<double, 3>;
using Circle2d = BasicCircle<double, 2>;

/// The parameter period, 2 pi.
inline constexpr double CirclePeriod = 2.0 * std::numbers::pi;

template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicPointN<T, D> Eval(const BasicCircle<T, D>& circle, const T& t)
{
    using std::cos;
    using std::sin;
    return circle.Center + circle.Radius * (cos(t) * circle.XAxis + sin(t) * circle.YAxis);
}

/// The n-th derivative: the radius vector turned by n quarter turns,
/// Radius * (cos(t + n pi/2) * XAxis + sin(t + n pi/2) * YAxis).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicVecN<T, D> Deriv(const BasicCircle<T, D>& circle, const T& t,
                                              const int n) STEREON_HOT_PRE(n >= 1)
{
    using std::cos;
    using std::sin;
    // cos and sin of t + n pi/2, by the quarter turn n mod 4, so that no
    // multiple of pi/2 is rounded into the argument.
    const T c = cos(t);
    const T s = sin(t);
    switch (n % 4) {
    case 1:
        return circle.Radius * (-s * circle.XAxis + c * circle.YAxis);
    case 2:
        return circle.Radius * (-c * circle.XAxis - s * circle.YAxis);
    case 3:
        return circle.Radius * (s * circle.XAxis - c * circle.YAxis);
    default:
        return circle.Radius * (c * circle.XAxis + s * circle.YAxis);
    }
}

/// [0, 2 pi].
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain([[maybe_unused]] const BasicCircle<T, D>& circle)
{
    return Interval{0.0, CirclePeriod};
}

/// A box containing the circle: along each axis i the circle reaches
/// Radius * sqrt(XAxis_i^2 + YAxis_i^2) from its centre. Computed in
/// interval arithmetic, so the box encloses the exact circle.
template <std::size_t D>
[[nodiscard]] constexpr BoxN<D> BoundingBox(const BasicCircle<double, D>& circle)
{
    BoxN<D> box;
    for (std::size_t i = 0; i < D; ++i) {
        // A sum of squares is never negative; outward rounding can step its
        // lower bound below zero, which Sqrt rejects, so it is clamped.
        const Interval squares =
            Squared(Interval{circle.XAxis[i]}) + Squared(Interval{circle.YAxis[i]});
        const Interval nonNegative{std::max(0.0, squares.GetLower()), squares.GetUpper()};
        const Interval reach = Interval{circle.Radius} * Sqrt(nonNegative);
        box[i] = Interval{circle.Center[i]} + Interval{-1.0, 1.0} * reach;
    }
    return box;
}

/// The parameter in [0, 2 pi) of the point of the circle closest to p: the
/// angle of p's projection onto the circle's plane. Every point of the
/// circle is equally close to a point on its axis; that case returns 0.
template <std::size_t D>
[[nodiscard]] double Project(const BasicCircle<double, D>& circle, const BasicPointN<double, D>& p)
{
    const BasicVecN<double, D> offset = p - circle.Center;
    const double x = Dot(offset, circle.XAxis);
    const double y = Dot(offset, circle.YAxis);
    if (x == 0.0 && y == 0.0) {
        return 0.0;
    }
    const double angle = std::atan2(y, x);
    if (angle >= 0.0) {
        return angle;
    }
    // A tiny negative angle plus 2 pi rounds to 2 pi itself, which is the
    // same point as 0 and outside [0, 2 pi).
    const double wrapped = angle + CirclePeriod;
    return wrapped < CirclePeriod ? wrapped : 0.0;
}

/// The same circle traversed the other way:
/// Eval(Reversed(c), t) == Eval(c, 2 pi - t).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicCircle<T, D> Reversed(const BasicCircle<T, D>& circle)
{
    return {circle.Center, circle.Radius, circle.XAxis, -circle.YAxis};
}

/// The image of the circle under x, which must be a similarity (a rigid
/// motion, possibly with uniform scaling); the parameter is kept. A
/// non-uniform scaling would make an ellipse, which this type cannot hold.
template <Scalar T>
[[nodiscard]] constexpr BasicCircle<T, 3> Transformed(const BasicCircle<T, 3>& circle,
                                                      const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(circle.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(circle.YAxis, x);
    const T scale = Length(xAxis);
    return {Transformed(circle.Center, x), circle.Radius * scale, xAxis / scale, yAxis / scale};
}

static_assert(CurveGeometry<Circle>);
static_assert(Curve2dGeometry<Circle2d>);

} // namespace Stereon
