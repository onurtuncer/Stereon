// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

export module Stereon.Surface:Cone;

import Stereon.Curve;
import :Concepts;

/// The circular cone.

export namespace Stereon {

/// Origin + (Radius + v sin a) * (cos u * XAxis + sin u * YAxis) +
/// v cos a * Axis(c), with a = SemiAngle, u in [0, 2 pi] and v real. v is
/// the signed distance along a generator from the reference circle of
/// radius Radius through Origin, so the generators are unit-speed lines.
/// XAxis and YAxis are orthonormal; the axis follows the handedness flag as
/// for BasicCylinder. Radius >= 0 and 0 < |SemiAngle| < pi/2: a positive
/// angle opens the cone towards +Axis(c), a negative one towards -Axis(c).
///
/// The apex is at v = -Radius / sin a, where every u gives the same point:
/// the parametrisation is singular there (a Degenerate edge of ADR-0009)
/// and the surface continues through it as the second nappe, whose points
/// have Radius + v sin a < 0. The conventions are STEP's conical_surface and
/// OCCT's Geom_ConicalSurface.
template <Scalar T>
struct BasicCone {
    BasicPoint3<T> Origin;
    T Radius;
    T SemiAngle;
    BasicVec3<T> XAxis;
    BasicVec3<T> YAxis;
    bool RightHanded = true;

    constexpr bool operator==(const BasicCone&) const = default;
};

using Cone = BasicCone<double>;

/// The cone after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, Radius >= 0,
/// 0 < |SemiAngle| < pi/2 and the axes are orthonormal within
/// FrameAllowance.
[[nodiscard]] Result<Cone> MakeCone(const Point3& origin, const double radius,
                                    const double semiAngle, const Vec3& xAxis, const Vec3& yAxis,
                                    const bool rightHanded = true);

/// The direction of the cone's axis.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Axis(const BasicCone<T>& cone)
{
    return FrameAxis(cone.XAxis, cone.YAxis, cone.RightHanded);
}

/// The apex, Origin - (Radius / tan a) * Axis(c).
template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Apex(const BasicCone<T>& cone)
{
    return cone.Origin - (cone.Radius * Cos(cone.SemiAngle) / Sin(cone.SemiAngle)) * Axis(cone);
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Eval(const BasicCone<T>& cone, const BasicPoint2<T>& uv)
{
    Detail::RequireFinite(uv, "Eval(Cone)", "a non-finite parameter");
    const T radius = cone.Radius + uv.Y * Sin(cone.SemiAngle);
    return cone.Origin + radius * RadialDirection(cone.XAxis, cone.YAxis, uv.X, 0) +
           (uv.Y * Cos(cone.SemiAngle)) * Axis(cone);
}

/// The cone is linear in v: the v-derivative of the radius is sin a and of
/// the height cos a, and every derivative of order 2 or more in v is zero.
/// Derivatives in u turn the radial direction (RadialDirection).
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Deriv(const BasicCone<T>& cone, const BasicPoint2<T>& uv,
                                           const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(Cone)", "a non-finite parameter");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(Cone)",
                    "a derivative order that does not exist");
    if (dv >= 2) {
        return {};
    }
    const BasicVec3<T> radial = RadialDirection(cone.XAxis, cone.YAxis, uv.X, du);
    if (dv == 0) {
        return (cone.Radius + uv.Y * Sin(cone.SemiAngle)) * radial;
    }
    if (du != 0) {
        return Sin(cone.SemiAngle) * radial;
    }
    return Sin(cone.SemiAngle) * radial + Cos(cone.SemiAngle) * Axis(cone);
}

/// The periods in u and v: 2 pi in u, none (0) in v.
template <Scalar T>
[[nodiscard]] constexpr Vec2 Periods([[maybe_unused]] const BasicCone<T>& cone) noexcept
{
    return {CirclePeriod, 0.0};
}

/// u in [0, 2 pi], v unbounded.
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain([[maybe_unused]] const BasicCone<T>& cone)
{
    return {Interval{0.0, CirclePeriod}, Interval::Everything()};
}

/// Everything: the cone is unbounded along its axis.
[[nodiscard]] constexpr Box3 BoundingBox([[maybe_unused]] const Cone& cone)
{
    return Box3::Everything();
}

} // namespace Stereon

namespace Stereon::Detail {

/// The closest point of a cone in the meridian plane of p, in that plane's
/// coordinates: rho along the radial direction at angle u, h along the
/// axis. The plane cuts the cone in two generator lines through the apex:
/// one with parameter u, (Radius + v sin a, v cos a), and one with u + pi,
/// (-(Radius + v sin a), v cos a). Each is a unit-speed line, so the foot of
/// the perpendicular has a closed form; the nearer one wins, the line at u
/// on a tie, which is what a point on the axis meets.
struct ConeFoot {
    bool Opposite;
    double V;
};

[[nodiscard]] inline ConeFoot ConeFootInMeridian(const Cone& cone, const double rho, const double h)
{
    const double s = Sin(cone.SemiAngle);
    const double c = Cos(cone.SemiAngle);
    // Along (s, c) from (Radius, 0), and along (-s, c) from (-Radius, 0).
    const double near = (rho - cone.Radius) * s + h * c;
    const double far = -(rho + cone.Radius) * s + h * c;
    const double nearRho = rho - (cone.Radius + near * s);
    const double nearH = h - near * c;
    const double farRho = rho + (cone.Radius + far * s);
    const double farH = h - far * c;
    const bool opposite = farRho * farRho + farH * farH < nearRho * nearRho + nearH * nearH;
    return {opposite, opposite ? far : near};
}

} // namespace Stereon::Detail

export namespace Stereon {

/// u in [0, 2 pi) and v of a point of the cone closest to p. The closest
/// point lies in the plane through the axis and p, on one of the two
/// generators in it (see ConeFootInMeridian); a point on the axis gets
/// u = 0. Near the apex, the generator of the opposite nappe can be the
/// nearer one, and u is then the angle of p plus pi.
[[nodiscard]] inline Point2 Project(const Cone& cone, const Point3& p)
{
    Detail::RequireFinite(p, "Project(Cone)", "a non-finite point");
    const Vec3 offset = p - cone.Origin;
    const double x = Dot(offset, cone.XAxis);
    const double y = Dot(offset, cone.YAxis);
    const Detail::ConeFoot foot =
        Detail::ConeFootInMeridian(cone, Hypot(x, y), Dot(offset, Axis(cone)));
    // The opposite generator is the one at u + pi: the angle of (-x, -y).
    return {foot.Opposite ? WrappedAngle(-y, -x) : WrappedAngle(y, x), foot.V};
}

/// The same cone with the opposite normal, by running u the other way:
/// Eval(Reversed(c), (u, v)) == Eval(c, (2 pi - u, v)).
template <Scalar T>
[[nodiscard]] constexpr BasicCone<T> Reversed(const BasicCone<T>& cone)
{
    return {cone.Origin, cone.Radius, cone.SemiAngle, cone.XAxis, -cone.YAxis, !cone.RightHanded};
}

/// The image of the cone under x, which must be a similarity. Angles are
/// kept and lengths scale, so the radius scales and v is rescaled by the
/// scale factor; u is kept.
template <Scalar T>
[[nodiscard]] constexpr BasicCone<T> Transformed(const BasicCone<T>& cone,
                                                 const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(cone.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(cone.YAxis, x);
    const T scale = Length(xAxis);
    const bool mirrored = Determinant(x.Linear) < T{};
    return {Transformed(cone.Origin, x),
            cone.Radius * scale,
            cone.SemiAngle,
            xAxis / scale,
            yAxis / scale,
            cone.RightHanded != mirrored};
}

static_assert(SurfaceGeometry<Cone>);
static_assert(Detail::FitsInline<Cone, SurfaceBufferSize>);

} // namespace Stereon
