// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>

export module Stereon.Surface:Torus;

import Stereon.Curve;
import :Concepts;

/// The ring torus.

export namespace Stereon {

/// Center + (MajorRadius + MinorRadius cos v) * (cos u * XAxis + sin u * YAxis)
/// + MinorRadius sin v * Axis(t), with u and v in [0, 2 pi]. u runs around
/// the axis, v around the tube; the tube's centre line is the circle of
/// radius MajorRadius in the (XAxis, YAxis) plane. XAxis and YAxis are
/// orthonormal; the axis follows the handedness flag as for BasicCylinder,
/// so the normal points out of the tube in a right-handed frame.
///
/// Only the ring torus, MajorRadius > MinorRadius > 0, is represented: the
/// horn and spindle tori, which touch or cross the axis, are
/// self-intersecting and singular, and STEP's degenerate_toroidal_surface
/// has its own entity. The conventions are STEP's toroidal_surface and
/// OCCT's Geom_ToroidalSurface.
template <Scalar T>
struct BasicTorus {
    BasicPoint3<T> Center;
    T MajorRadius;
    T MinorRadius;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;
    bool RightHanded = true;

    constexpr bool operator==(const BasicTorus&) const = default;
};

using Torus = BasicTorus<double>;

/// The torus after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, MajorRadius > MinorRadius > 0
/// and the axes are orthonormal within FrameAllowance.
[[nodiscard]] Result<Torus> MakeTorus(const Point3& center, const double majorRadius,
                                      const double minorRadius, const Vec3& xAxis,
                                      const Vec3& yAxis, const bool rightHanded = true);

/// The direction of the torus's axis.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Axis(const BasicTorus<T>& torus)
{
    return FrameAxis(torus.XAxis, torus.YAxis, torus.RightHanded);
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicTorus<T>& torus, const BasicPoint2<T>& uv)
{
    Detail::RequireFinite(uv, "Eval(Torus)", "a non-finite parameter");
    const T radius = torus.MajorRadius + torus.MinorRadius * Cos(uv.Y);
    return torus.Center + radius * RadialDirection(torus.XAxis, torus.YAxis, uv.X, 0) +
           (torus.MinorRadius * Sin(uv.Y)) * Axis(torus);
}

/// As for the sphere, the point separates into a u part and a v part:
/// (R + r cos v) * radial(u) + r sin v * axis. A derivative of order dv >= 1
/// in v drops R and turns (cos v, sin v) by dv quarter turns; the axis term
/// survives only when du is 0.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicTorus<T>& torus, const BasicPoint2<T>& uv,
                                           const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(Torus)", "a non-finite parameter");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(Torus)",
                    "a derivative order that does not exist");
    const BasicCosSin<T> tube = QuarterTurned(Cos(uv.Y), Sin(uv.Y), dv);
    const T radial = (dv == 0 ? torus.MajorRadius : T{}) + torus.MinorRadius * tube.Cos;
    const BasicVec3<T> around = radial * RadialDirection(torus.XAxis, torus.YAxis, uv.X, du);
    if (du != 0) {
        return around;
    }
    return around + (torus.MinorRadius * tube.Sin) * Axis(torus);
}

/// The periods in u and v: 2 pi in both.
template <Scalar T>
[[nodiscard]] constexpr Vec2 Periods([[maybe_unused]] const BasicTorus<T>& torus) noexcept
{
    return {CirclePeriod, CirclePeriod};
}

/// [0, 2 pi] in both u and v.
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicTorus<T>& torus)
{
    return {Interval{0.0, CirclePeriod}, Interval{0.0, CirclePeriod}};
}

} // namespace Stereon

namespace Stereon::Detail {

/// How far the torus reaches from its centre along one coordinate axis,
/// given the components x, y and z of XAxis, YAxis and the axis along it:
/// R * sqrt(x^2 + y^2) + r * sqrt(x^2 + y^2 + z^2), in interval arithmetic
/// (see BoundingBox(Torus)). The components are passed one axis at a time,
/// not indexed: an inline non-template function that instantiates core's
/// indexed accessors corrupts GCC 16.2's module file
/// (docs/toolchain/gcc16-feature-probe.md).
[[nodiscard]] constexpr Interval TorusReach(const Torus& torus, const double x, const double y,
                                            const double z)
{
    const Interval around = Squared(Interval{x}) + Squared(Interval{y});
    const Interval all = around + Squared(Interval{z});
    // Clamped as in BoundingBox(Circle).
    const Interval aroundNonNegative{std::max(0.0, around.GetLower()), around.GetUpper()};
    const Interval allNonNegative{std::max(0.0, all.GetLower()), all.GetUpper()};
    const Interval reach = Interval{torus.MajorRadius} * Sqrt(aroundNonNegative) +
                           Interval{torus.MinorRadius} * Sqrt(allNonNegative);
    return Interval{-1.0, 1.0} * reach;
}

} // namespace Stereon::Detail

export namespace Stereon {

/// A box containing the torus. Along axis i the coordinate is
/// (R + r cos v) * A cos(u - u_i) + r sin v * Z_i with A = sqrt(X_i^2 + Y_i^2),
/// which is at most R * A + r * (A cos v + Z_i sin v) because R + r cos v is
/// positive, and so at most R * A + r * sqrt(A^2 + Z_i^2). The bound is tight
/// and is computed in interval arithmetic with the axis the evaluator uses,
/// so the box encloses the exact torus.
[[nodiscard]] constexpr Box3 BoundingBox(const Torus& torus)
{
    const Vec3 axis = Axis(torus);
    const Vec3& x = torus.XAxis;
    const Vec3& y = torus.YAxis;
    return {Interval{torus.Center.X} + Detail::TorusReach(torus, x.X, y.X, axis.X),
            Interval{torus.Center.Y} + Detail::TorusReach(torus, x.Y, y.Y, axis.Y),
            Interval{torus.Center.Z} + Detail::TorusReach(torus, x.Z, y.Z, axis.Z)};
}

/// u and v, both in [0, 2 pi), of a point of the torus closest to p: u is
/// the angle of p around the axis, and v the angle of p around the tube's
/// centre circle in the half-plane through the axis and p, which holds the
/// nearer of the two tube sections. A point on the axis gets u = 0, and a
/// point on the centre circle v = 0; every candidate is equally close there.
[[nodiscard]] inline Point2 Project(const Torus& torus, const Point3& p)
{
    Detail::RequireFinite(p, "Project(Torus)", "a non-finite point");
    const Vec3 offset = p - torus.Center;
    const double x = Dot(offset, torus.XAxis);
    const double y = Dot(offset, torus.YAxis);
    const double h = Dot(offset, Axis(torus));
    return {WrappedAngle(y, x), WrappedAngle(h, Hypot(x, y) - torus.MajorRadius)};
}

/// The same torus with the opposite normal, by running u the other way:
/// Eval(Reversed(t), (u, v)) == Eval(t, (2 pi - u, v)).
template <Scalar T>
[[nodiscard]] constexpr BasicTorus<T> Reversed(const BasicTorus<T>& torus)
{
    return {torus.Center, torus.MajorRadius, torus.MinorRadius,
            torus.XAxis,  -torus.YAxis,      !torus.RightHanded};
}

/// The image of the torus under x, which must be a similarity; the
/// parameters are kept and both radii scale.
template <Scalar T>
[[nodiscard]] constexpr BasicTorus<T> Transformed(const BasicTorus<T>& torus,
                                                  const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(torus.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(torus.YAxis, x);
    const T scale = Length(xAxis);
    const bool mirrored = Determinant(x.Linear) < T{};
    return {Transformed(torus.Center, x),
            torus.MajorRadius * scale,
            torus.MinorRadius * scale,
            xAxis / scale,
            yAxis / scale,
            torus.RightHanded != mirrored};
}

static_assert(SurfaceGeometry<Torus>);
static_assert(Detail::FitsInline<Torus, GeometryBufferSize>);

} // namespace Stereon
