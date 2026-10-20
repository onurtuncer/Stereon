// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>

export module Stereon.Curve:Ellipse;

import Stereon.Robust;
import :Circle;
import :Concepts;

/// The ellipse, in space and in the parameter plane.

export namespace Stereon {

/// Center + MajorRadius * cos t * XAxis + MinorRadius * sin t * YAxis for t
/// in [0, 2 pi]. XAxis and YAxis are orthonormal, XAxis along the major
/// axis; MajorRadius >= MinorRadius > 0. In space the normal is
/// Cross(XAxis, YAxis); in the plane YAxis decides the sense, as for
/// BasicCircle. The parameter t is the eccentric angle, not the polar angle
/// of the point: the two agree only at the four vertices.
///
/// The conventions are those of STEP's ellipse and OCCT's Geom_Ellipse,
/// so that geometry exchanged with either keeps its parameters.
template <Scalar T, std::size_t D>
struct BasicEllipse {
    BasicPointN<T, D> Center;
    T MajorRadius;
    T MinorRadius;
    BasicVecN<T, D> XAxis;
    BasicVecN<T, D> YAxis;

    constexpr bool operator==(const BasicEllipse&) const = default;
};

using Ellipse = BasicEllipse<double, 3>;
using Ellipse2d = BasicEllipse<double, 2>;

/// The ellipse after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, MajorRadius >= MinorRadius > 0
/// and the axes are orthonormal within FrameAllowance.
[[nodiscard]] Result<Ellipse> MakeEllipse(const Point3& center, const double majorRadius,
                                          const double minorRadius, const Vec3& xAxis,
                                          const Vec3& yAxis);
[[nodiscard]] Result<Ellipse2d> MakeEllipse(const Point2& center, const double majorRadius,
                                            const double minorRadius, const Vec2& xAxis,
                                            const Vec2& yAxis);

template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicPointN<T, D> Eval(const BasicEllipse<T, D>& ellipse, const T& t)
{
    Detail::RequireFinite(t, "Eval(Ellipse)", "a non-finite parameter");
    return ellipse.Center + (ellipse.MajorRadius * Cos(t)) * ellipse.XAxis +
           (ellipse.MinorRadius * Sin(t)) * ellipse.YAxis;
}

/// The n-th derivative, MajorRadius * cos(t + n pi/2) * XAxis +
/// MinorRadius * sin(t + n pi/2) * YAxis, by quarter turns (QuarterTurned).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicVecN<T, D> Deriv(const BasicEllipse<T, D>& ellipse, const T& t,
                                              const int n)
{
    Detail::RequireFinite(t, "Deriv(Ellipse)", "a non-finite parameter");
    Detail::Require(n >= 1, "Deriv(Ellipse)", "a derivative order below 1");
    const BasicCosSin<T> turned = QuarterTurned(Cos(t), Sin(t), n);
    return (ellipse.MajorRadius * turned.Cos) * ellipse.XAxis +
           (ellipse.MinorRadius * turned.Sin) * ellipse.YAxis;
}

/// [0, 2 pi].
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain([[maybe_unused]] const BasicEllipse<T, D>& ellipse)
{
    return Interval{0.0, CirclePeriod};
}

/// 2 pi: every finite parameter is valid.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr double Period([[maybe_unused]] const BasicEllipse<T, D>& ellipse) noexcept
{
    return CirclePeriod;
}

/// A box containing the ellipse: along axis i the coordinate
/// a cos t X_i + b sin t Y_i reaches sqrt((a X_i)^2 + (b Y_i)^2) from the
/// centre, the amplitude of a sum of a cosine and a sine. Computed in
/// interval arithmetic, so the box encloses the exact ellipse.
template <std::size_t D>
[[nodiscard]] constexpr BoxN<D> BoundingBox(const BasicEllipse<double, D>& ellipse)
{
    BoxN<D> box;
    for (std::size_t i = 0; i < D; ++i) {
        const Interval squares =
            Squared(Interval{ellipse.MajorRadius} * Interval{ellipse.XAxis[i]}) +
            Squared(Interval{ellipse.MinorRadius} * Interval{ellipse.YAxis[i]});
        // Clamped as in BoundingBox(Circle): outward rounding can step the
        // lower bound of a sum of squares below zero.
        const Interval nonNegative{std::max(0.0, squares.GetLower()), squares.GetUpper()};
        box[i] = Interval{ellipse.Center[i]} + Interval{-1.0, 1.0} * Sqrt(nonNegative);
    }
    return box;
}

} // namespace Stereon

namespace Stereon::Detail {

/// The largest number of bisection steps EllipseRoot needs: one per binary
/// digit of a double between the smallest subnormal and 2, so that the
/// bracket always shrinks to adjacent doubles first (Eberly).
inline constexpr int EllipseRootSteps =
    std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent;

/// The root s >= s0 of F(s) = (n0 / (s + r0))^2 + (z1 / (s + 1))^2 - 1,
/// bracketed by [s0, s1] and found by bisection, which converges for every
/// input where Newton's method can stall (Eberly, "Distance from a point to
/// an ellipse, an ellipsoid, or a hyperellipsoid", 2013).
[[nodiscard]] inline double EllipseRoot(const double r0, const double z0, const double z1,
                                        const double g)
{
    const double n0 = r0 * z0;
    double s0 = z1 - 1.0;
    double s1 = g < 0.0 ? 0.0 : Hypot(n0, z1) - 1.0;
    double s = 0.0;
    for (int i = 0; i < EllipseRootSteps; ++i) {
        s = 0.5 * (s0 + s1);
        if (s == s0 || s == s1) {
            break; // the bracket holds two adjacent doubles
        }
        const double ratio0 = n0 / (s + r0);
        const double ratio1 = z1 / (s + 1.0);
        const double value = ratio0 * ratio0 + ratio1 * ratio1 - 1.0;
        if (value > 0.0) {
            s0 = s;
        } else if (value < 0.0) {
            s1 = s;
        } else {
            break;
        }
    }
    return s;
}

/// The point (x0, x1) of the ellipse (x0/e0)^2 + (x1/e1)^2 = 1 closest to
/// (y0, y1), for e0 >= e1 > 0 and a point in the first quadrant, y0, y1 >= 0
/// (Eberly's case analysis). The closest point is then in the first quadrant
/// too, and the cases are: off both axes, on the minor axis, and on the
/// major axis, where a point inside the evolute has two closest points off
/// the axis and one outside it has the vertex.
[[nodiscard]] inline Point2 ClosestOnEllipseQuadrant(const double e0, const double e1,
                                                     const double y0, const double y1)
{
    if (y1 > 0.0 && y0 > 0.0) {
        const double z0 = y0 / e0;
        const double z1 = y1 / e1;
        const double g = z0 * z0 + z1 * z1 - 1.0;
        if (g == 0.0) {
            return {y0, y1}; // on the ellipse
        }
        const double r0 = (e0 / e1) * (e0 / e1);
        const double s = EllipseRoot(r0, z0, z1, g);
        return {r0 * y0 / (s + r0), y1 / (s + 1.0)};
    }
    if (y1 > 0.0) {
        return {0.0, e1}; // on the minor axis
    }
    const double numerator = e0 * y0;
    const double denominator = e0 * e0 - e1 * e1;
    if (numerator < denominator) {
        const double x = numerator / denominator;
        return {e0 * x, e1 * Sqrt(std::max(0.0, 1.0 - x * x))};
    }
    return {e0, 0.0};
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The parameter in [0, 2 pi) of a point of the ellipse closest to p. In
/// space the closest point is the closest point to p's projection onto the
/// ellipse's plane, since the distance splits into an in-plane part and a
/// normal part. A point on the major axis inside the evolute has two
/// closest points, mirror images across the axis; it gets the one with
/// t in (0, pi), so the result is deterministic. The centre is equally near
/// both vertices of the minor axis and gets the one at pi/2, unless the
/// ellipse is a circle, where it gets 0 as Project(Circle) does.
template <std::size_t D>
[[nodiscard]] double Project(const BasicEllipse<double, D>& ellipse,
                             const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(Ellipse)", "a non-finite point");
    const BasicVecN<double, D> offset = p - ellipse.Center;
    const double x = Dot(offset, ellipse.XAxis);
    const double y = Dot(offset, ellipse.YAxis);
    const Point2 q = Detail::ClosestOnEllipseQuadrant(ellipse.MajorRadius, ellipse.MinorRadius,
                                                      x < 0.0 ? -x : x, y < 0.0 ? -y : y);
    // Reflected back to p's quadrant; the eccentric angle of (q.X, q.Y) is
    // the polar angle of (q.X / a, q.Y / b).
    const double cosine = (x < 0.0 ? -q.X : q.X) / ellipse.MajorRadius;
    const double sine = (y < 0.0 ? -q.Y : q.Y) / ellipse.MinorRadius;
    return WrappedAngle(sine, cosine);
}

/// The same ellipse traversed the other way:
/// Eval(Reversed(e), t) == Eval(e, 2 pi - t).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicEllipse<T, D> Reversed(const BasicEllipse<T, D>& ellipse)
{
    return {ellipse.Center, ellipse.MajorRadius, ellipse.MinorRadius, ellipse.XAxis,
            -ellipse.YAxis};
}

/// The image of the ellipse under x, which must be a similarity; the
/// parameter is kept. A general affine map also takes an ellipse to an
/// ellipse, but with new axes and a shifted parameter, which this function
/// does not compute.
template <Scalar T>
[[nodiscard]] constexpr BasicEllipse<T, 3> Transformed(const BasicEllipse<T, 3>& ellipse,
                                                       const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(ellipse.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(ellipse.YAxis, x);
    const T scale = Length(xAxis);
    return {Transformed(ellipse.Center, x), ellipse.MajorRadius * scale,
            ellipse.MinorRadius * scale, xAxis / scale, yAxis / scale};
}

static_assert(CurveGeometry<Ellipse>);
static_assert(Curve2dGeometry<Ellipse2d>);
static_assert(Detail::FitsInline<Ellipse, CurveBufferSize>);
static_assert(Detail::FitsInline<Ellipse2d, Curve2dBufferSize>);

} // namespace Stereon
