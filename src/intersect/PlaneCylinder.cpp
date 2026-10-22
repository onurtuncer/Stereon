// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Exact.hpp"

/// Plane–cylinder (ADR-0010's analytic table, RS-0001 Tier A). The plane's
/// normal n and the cylinder's axis a are cross products of the input axes,
/// so the configuration is decided exactly: n x a = 0 gives a circle; n . a
/// = 0 puts the axis along the plane, and the cylinder's height over the
/// plane against its radius gives two lines, one tangent line or a proved
/// separation; anything else gives an ellipse.
///
/// Every branch is parametrised by the cylinder's angle, so the pcurve on
/// the cylinder is a line for lines and circles. An ellipse's pcurve on the
/// cylinder is the sinusoid (u, v_c + B cos t), which no exact 2D type
/// represents: it is a cubic B-spline through samples, dense enough that
/// the interpolant stays within the Context's resolution of it (certified
/// by M4). The pcurve on the plane is the exact Ellipse2d.

namespace Stereon::Detail {

namespace {

/// The two surfaces with their normal and axis direction.
struct Sides {
    const Plane& P;
    Vec3 Normal;
    const Cylinder& C;
    Vec3 Axis;
};

[[nodiscard]] Sides SidesOf(const Plane& p, const Cylinder& c)
{
    return {p, Cross(p.XAxis, p.YAxis), c, Stereon::Axis(c)};
}

[[nodiscard]] ExactVector<4> ExactAxis(const Cylinder& c)
{
    const ExactVector<4> cross = ExactCross(Exactly(c.XAxis), Exactly(c.YAxis));
    return c.RightHanded ? cross : ExactNegated(cross);
}

[[nodiscard]] Vec2 InPlane(const Plane& plane, const Vec3& v)
{
    return {Dot(v, plane.XAxis), Dot(v, plane.YAxis)};
}

/// The axis is perpendicular to the plane: the section circle at the
/// plane's height, on the cylinder's angle.
[[nodiscard]] SurfaceSurfaceIntersection CircleSection(const Sides& s, const Box2& boxA,
                                                       const Box2& boxB)
{
    const double v = Dot(s.P.Origin - s.C.Origin, s.Axis);
    const Point3 center = s.C.Origin + v * s.Axis;
    const Circle circle{center, s.C.Radius, s.C.XAxis, s.C.YAxis};
    const Circle2d onPlane{Project(s.P, center), s.C.Radius, InPlane(s.P, s.C.XAxis),
                           InPlane(s.P, s.C.YAxis)};
    const Line2d onCylinder{Point2{0.0, v}, Vec2{1.0, 0.0}};
    const Box2 planeBox{Around(onPlane.Center.X, circle.Radius),
                        Around(onPlane.Center.Y, circle.Radius)};
    const Interval range = Intersection(Domain(circle), boxB.X);
    SurfaceSurfaceIntersection result;
    if (!range.IsEmpty() && Contains(boxB.Y, v) && Intersects(planeBox, boxA)) {
        result.Branches.push_back(
            Crossing(Curve(circle), Curve2d(onPlane), Curve2d(onCylinder), range));
    }
    return result;
}

/// The cylinder's ruling at angle u, with the line on each surface.
void AddRuling(SurfaceSurfaceIntersection& result, const Sides& s, const double u,
               const bool tangent, const Box2& boxA, const Box2& boxB)
{
    const Point3 point = Eval(s.C, Point2{u, 0.0});
    const Line2d onPlane{Project(s.P, point), InPlane(s.P, s.Axis)};
    const Line2d onCylinder{Point2{u, 0.0}, Vec2{0.0, 1.0}};
    const Interval range = Intersection(ClipLine(onPlane, boxA), ClipLine(onCylinder, boxB));
    if (range.IsEmpty()) {
        return;
    }
    const Curve line{Line{point, s.Axis}};
    result.Branches.push_back(tangent
                                  ? Tangency(line, Curve2d(onPlane), Curve2d(onCylinder), range)
                                  : Crossing(line, Curve2d(onPlane), Curve2d(onCylinder), range));
}

/// The axis lies along the plane: the rulings where the cylinder's height
/// over the plane, h + r cos(u - phi), vanishes. `heightSign` compares h
/// with r: zero for a tangency, negative for two rulings.
[[nodiscard]] SurfaceSurfaceIntersection Rulings(const Sides& s, const Box2& boxA, const Box2& boxB,
                                                 const Sign heightSign)
{
    const Vec3 offset = s.C.Origin - s.P.Origin;
    const double phi = Atan2(Dot(s.Normal, s.C.YAxis), Dot(s.Normal, s.C.XAxis));
    const double cosine =
        std::clamp(-Dot(offset, s.Normal) / (s.C.Radius * Length(s.Normal)), -1.0, 1.0);
    SurfaceSurfaceIntersection result;
    if (heightSign == Sign::Zero) {
        // Exactly tangent: the ruling on the plane's side of the axis.
        const double u = Dot(offset, s.Normal) > 0.0 ? phi + std::numbers::pi : phi;
        AddRuling(result, s, WrappedAngle(Sin(u), Cos(u)), true, boxA, boxB);
        return result;
    }
    const double delta = std::acos(cosine);
    std::array<double, 2> angles{WrappedAngle(Sin(phi - delta), Cos(phi - delta)),
                                 WrappedAngle(Sin(phi + delta), Cos(phi + delta))};
    std::ranges::sort(angles); // Deterministic order (ADR-0029 G1).
    for (const double u : angles) {
        AddRuling(result, s, u, false, boxA, boxB);
    }
    return result;
}

/// The cylinder's pcurve of an ellipse, (u0 + t, v + amplitude cos t) over
/// `range`: the clamped cubic interpolant through samples with the exact
/// end derivatives, dense enough for its interpolation error bound,
/// 5/384 h^4 max|f''''|, to stay within half the Context's resolution.
[[nodiscard]] Result<Curve2d> SinusoidPcurve(const double u0, const double v,
                                             const double amplitude, const Interval& range,
                                             const Context& ctx)
{
    const double budget = ctx.GetResolution() / 2.0;
    const double spacing =
        std::pow(384.0 * budget / (5.0 * std::max(std::abs(amplitude), budget)), 0.25);
    // Not Width(range): GCC 16.2 rejects its contract here as not constant.
    const double width = range.GetUpper() - range.GetLower();
    const auto n = static_cast<std::size_t>(std::max(16.0, std::ceil(width / spacing)));
    std::vector<Point2> samples(n + 1);
    for (std::size_t k = 0; k <= n; ++k) {
        const double t =
            range.GetLower() + width * (static_cast<double>(k) / static_cast<double>(n));
        samples[k] = Point2{u0 + t, v + amplitude * Cos(t)};
    }
    // d/dt (u0 + t, v + B cos t) = (1, -B sin t), with respect to the fit's
    // parameter on [0, 1], which runs `width` times faster.
    const auto slope = [&](const double t) { return Vec2{width, -amplitude * Sin(t) * width}; };
    const EndDerivatives2d ends{slope(range.GetLower()), slope(range.GetUpper())};
    STEREON_TRY(const CurveFit2d fit, Interpolate(samples, 3, ends, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d pcurve, Reparametrized(fit.Curve, range));
    return Curve2d(std::move(pcurve));
}

/// Neither parallel nor perpendicular: the ellipse, parametrised by the
/// cylinder's angle u = u0 + t, where u0 is the angle of the plane's tilt
/// around the axis. Its major axis lies along the tilt, its minor axis is
/// the cylinder's radius across it.
[[nodiscard]] Result<SurfaceSurfaceIntersection>
EllipseSection(const Sides& s, const Box2& boxA, const Box2& boxB, const Context& ctx)
{
    const double gamma = Dot(s.Normal, s.Axis);
    const Vec3 tilt = Normalized(s.Normal - gamma * s.Axis); // The normal across the axis.
    const double u0 = WrappedAngle(Dot(tilt, s.C.YAxis), Dot(tilt, s.C.XAxis));
    const double vCenter = Dot(s.P.Origin - s.C.Origin, s.Normal) / gamma;
    const Point3 center = s.C.Origin + vCenter * s.Axis;
    const double lambda = -Dot(tilt, s.Normal) / gamma; // tilt + lambda axis lies in the plane.
    const Vec3 major = tilt + lambda * s.Axis;
    const double majorRadius = s.C.Radius * Length(major);
    const Vec3 minor = Dot(tilt, s.C.XAxis) * s.C.YAxis - Dot(tilt, s.C.YAxis) * s.C.XAxis;
    const Ellipse ellipse{center, majorRadius, s.C.Radius, Normalized(major), minor};
    const Ellipse2d onPlane{Project(s.P, center), majorRadius, s.C.Radius,
                            InPlane(s.P, ellipse.XAxis), InPlane(s.P, minor)};
    const double amplitude = s.C.Radius * lambda; // v(t) = vCenter + amplitude cos t.
    // t runs from the seam round to it, so the pcurve stays within [0, 2 pi] in u.
    const Interval range = Intersection(Interval{-u0, CirclePeriod - u0},
                                        Interval{boxB.X.GetLower() - u0, boxB.X.GetUpper() - u0});
    const Box2 planeBox{Around(onPlane.Center.X, majorRadius),
                        Around(onPlane.Center.Y, majorRadius)};
    SurfaceSurfaceIntersection result;
    if (range.IsEmpty() || !Intersects(boxB.Y, Around(vCenter, std::abs(amplitude))) ||
        !Intersects(planeBox, boxA)) {
        return result;
    }
    STEREON_TRY(Curve2d onCylinder, SinusoidPcurve(u0, vCenter, amplitude, range, ctx));
    result.Branches.push_back(
        Crossing(Curve(ellipse), Curve2d(onPlane), std::move(onCylinder), range));
    return result;
}

} // namespace

Result<SurfaceSurfaceIntersection> PlaneCylinder(const Surface& a, const Box2& boxA,
                                                 const Surface& b, const Box2& boxB,
                                                 const double /*tolerance*/, const Context& ctx)
{
    const Plane* plane = a.Target<Plane>();
    const Cylinder* cylinder = b.Target<Cylinder>();
    if (plane == nullptr || cylinder == nullptr) {
        return Fail(ErrorCode::InvalidInput, "PlaneCylinder: dispatched for other surface types");
    }
    const Sides s = SidesOf(*plane, *cylinder);
    const ExactVector<4> normal = ExactCross(Exactly(s.P.XAxis), Exactly(s.P.YAxis));
    const ExactVector<4> axis = ExactAxis(s.C);
    if (IsExactlyZero(ExactCross(normal, axis))) {
        return CircleSection(s, boxA, boxB);
    }
    if (SignOf(ExactDot(normal, axis)) != Sign::Zero) {
        return EllipseSection(s, boxA, boxB, ctx);
    }
    // The axis runs along the plane: the axis's height over it against the radius.
    const auto height = ExactDot(ExactOffset(s.C.Origin, s.P.Origin), normal);
    const Sign sign = CompareSquared(height, Exactly(s.C.Radius), ExactDot(normal, normal));
    if (sign == Sign::Positive) {
        return SurfaceSurfaceIntersection{}; // Proved separate.
    }
    return Rulings(s, boxA, boxB, sign);
}

Result<SurfaceSurfaceIntersection> CylinderPlane(const Surface& a, const Box2& boxA,
                                                 const Surface& b, const Box2& boxB,
                                                 const double tolerance, const Context& ctx)
{
    STEREON_TRY(SurfaceSurfaceIntersection result, PlaneCylinder(b, boxB, a, boxA, tolerance, ctx));
    return Swapped(std::move(result));
}

} // namespace Stereon::Detail
