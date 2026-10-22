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
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Exact.hpp"

/// Cylinder–cylinder with parallel axes (ADR-0010's analytic table, RS-0001
/// Tier A). Whether the axes are parallel is decided exactly on the input
/// doubles; nearly parallel axes are the skew case, the quartic of M3b,
/// which reports NotImplemented until then. For parallel axes the distance
/// between them against the sum and difference of the radii decides
/// between a proved separation (apart or nested), an external or internal
/// tangent ruling, two rulings, and a coincidence, exact or declared within
/// the pair tolerance: the surfaces' Hausdorff distance is the axes'
/// distance plus the radii's difference.

namespace Stereon::Detail {

namespace {

/// A cylinder, erased and as itself, with its axis direction.
struct Frame {
    const Surface& Erased;
    const Cylinder& Geometry;
    Vec3 Axis;
    ExactVector<4> ExactAxis;
};

[[nodiscard]] Frame FrameOf(const Surface& erased, const Cylinder& c)
{
    const ExactVector<4> cross = ExactCross(Exactly(c.XAxis), Exactly(c.YAxis));
    return {erased, c, Stereon::Axis(c), c.RightHanded ? cross : ExactNegated(cross)};
}

/// The angle of a radial direction in the cylinder's frame.
[[nodiscard]] double AngleOn(const Cylinder& c, const Vec3& radial)
{
    return WrappedAngle(Dot(radial, c.YAxis), Dot(radial, c.XAxis));
}

/// The shared ruling through `point`, along A's axis, with a line on each
/// cylinder. B's v runs the same way or the opposite one.
[[nodiscard]] Result<void> AddRuling(SurfaceSurfaceIntersection& result, const Frame& a,
                                     const Frame& b, const Point3& point, const bool tangent,
                                     const Box2& boxA, const Box2& boxB, const Context& ctx)
{
    const Vec3 fromA = point - a.Geometry.Origin;
    const Vec3 fromB = point - b.Geometry.Origin;
    const double sense = Dot(a.Axis, b.Axis) > 0.0 ? 1.0 : -1.0;
    const Line2d onA{Point2{AngleOn(a.Geometry, fromA), 0.0}, Vec2{0.0, 1.0}};
    const Line2d onB{Point2{AngleOn(b.Geometry, fromB), Dot(fromB, b.Axis)}, Vec2{0.0, sense}};
    const Interval range = Intersection(ClipLine(onA, boxA), ClipLine(onB, boxB));
    if (range.IsEmpty()) {
        return {};
    }
    const Curve line{Line{point, a.Axis}};
    return Keep(result,
                tangent ? Tangency(line, Curve2d(onA), Curve2d(onB), range)
                        : Crossing(line, Curve2d(onA), Curve2d(onB), range),
                a.Erased, b.Erased, ctx);
}

/// The rulings of cylinders whose axes are `distance` apart across them,
/// B's in direction `across` from A's: the two where the sections cross
/// (`sumSign` negative, `differenceSign` positive), or the one where they
/// touch, externally or internally.
[[nodiscard]] Result<SurfaceSurfaceIntersection>
Rulings(const Frame& a, const Frame& b, const Vec3& across, const double distance,
        const Sign sumSign, const Sign differenceSign, const Box2& boxA, const Box2& boxB,
        const Context& ctx)
{
    const double rA = a.Geometry.Radius;
    const double rB = b.Geometry.Radius;
    SurfaceSurfaceIntersection result;
    if (sumSign == Sign::Zero) {
        STEREON_TRY_VOID(
            AddRuling(result, a, b, a.Geometry.Origin + rA * across, true, boxA, boxB, ctx));
    } else if (differenceSign == Sign::Zero) {
        // B's axis is inside A when A is the larger; the touch is on its side.
        const double side = rA > rB ? rA : -rA;
        STEREON_TRY_VOID(
            AddRuling(result, a, b, a.Geometry.Origin + side * across, true, boxA, boxB, ctx));
    } else {
        const double cosine = std::clamp(
            (distance * distance + rA * rA - rB * rB) / (2.0 * distance * rA), -1.0, 1.0);
        const double phi = std::acos(cosine);
        const Vec3 along = Cross(a.Axis, across);
        std::array<Point3, 2> points{
            a.Geometry.Origin + rA * (Cos(phi) * across - Sin(phi) * along),
            a.Geometry.Origin + rA * (Cos(phi) * across + Sin(phi) * along)};
        // Deterministic order, by the angle on A (ADR-0029 G1).
        std::ranges::sort(points, {}, [&](const Point3& p) {
            return AngleOn(a.Geometry, p - a.Geometry.Origin);
        });
        for (const Point3& p : points) {
            STEREON_TRY_VOID(AddRuling(result, a, b, p, false, boxA, boxB, ctx));
        }
    }
    return result;
}

/// Whether the surfaces are within the pair tolerance as a whole: their
/// Hausdorff distance, the axes' distance plus the radii's difference, is at
/// most the tolerance, decided exactly as |w x a|² <= (tol - |rA - rB|)² |a|²
/// with a non-negative margin.
[[nodiscard]] bool WithinTolerance(const ExactScalar<6144>& squaredDistance,
                                   const ExactScalar<2>& radiiDifference,
                                   const ExactScalar<96>& scale, const double tolerance)
{
    const ExactScalar<3> margin{
        Interval{tolerance} - Hull(radiiDifference.Hull, -radiiDifference.Hull),
        Difference(Expansion<1>(tolerance), SignOf(radiiDifference) == Sign::Negative
                                                ? Negated(radiiDifference.Value)
                                                : radiiDifference.Value)};
    return SignOf(margin) != Sign::Negative &&
           Compare(squaredDistance, margin, scale) != Sign::Positive;
}

/// Coincident cylinders: an overlap whose map turns A's angle into B's and
/// shifts the height; the angles run the same way when the frames turn the
/// same way about the axis.
[[nodiscard]] SurfaceSurfaceIntersection Coincident(const Frame& a, const Box2& boxA,
                                                    const Frame& b, const Box2& boxB,
                                                    const double distance, const double tolerance,
                                                    const Context& ctx)
{
    const Cylinder& ca = a.Geometry;
    const Cylinder& cb = b.Geometry;
    const bool sameTurn =
        SignOf(ExactDot(ExactCross(Exactly(ca.XAxis), Exactly(ca.YAxis)),
                        ExactCross(Exactly(cb.XAxis), Exactly(cb.YAxis)))) == Sign::Positive;
    OverlapRegion overlap;
    overlap.AToB.U = Vec2{sameTurn ? 1.0 : -1.0, 0.0};
    overlap.AToB.V = Vec2{0.0, Dot(a.Axis, b.Axis) > 0.0 ? 1.0 : -1.0};
    overlap.AToB.Offset = Vec2{AngleOn(cb, ca.XAxis), Dot(ca.Origin - cb.Origin, b.Axis)};
    overlap.OnA = boxA;
    overlap.OnB = boxB;
    overlap.FullDomain = boxA == Domain(ca) && boxB == Domain(cb);
    overlap.SameOrientation = ca.RightHanded == cb.RightHanded;
    if (distance > 0.0) {
        overlap.Evidence = CoincidenceEvidence::Tolerant;
        overlap.Distance = distance;
        ctx.Log(LogLevel::Info, "Intersect: exactly parallel cylinders " + Number(distance) +
                                    " mm apart, within the pair tolerance " + Number(tolerance) +
                                    " mm: declared coincident, nothing moved (ADR-0029)");
    }
    SurfaceSurfaceIntersection result;
    result.Overlaps.push_back(std::move(overlap));
    return result;
}

} // namespace

Result<SurfaceSurfaceIntersection> CylinderCylinder(const Surface& a, const Box2& boxA,
                                                    const Surface& b, const Box2& boxB,
                                                    const double tolerance, const Context& ctx)
{
    const Cylinder* ca = a.Target<Cylinder>();
    const Cylinder* cb = b.Target<Cylinder>();
    if (ca == nullptr || cb == nullptr) {
        return Fail(ErrorCode::InvalidInput,
                    "CylinderCylinder: dispatched for other surface types");
    }
    const Frame fa = FrameOf(a, *ca);
    const Frame fb = FrameOf(b, *cb);
    if (!IsExactlyZero(ExactCross(fa.ExactAxis, fb.ExactAxis))) {
        return Fail(ErrorCode::NotImplemented,
                    "Cylinder x Cylinder with skew or intersecting axes is not implemented yet "
                    "(the Boolean prototype's M3b)");
    }
    const double rA = fa.Geometry.Radius;
    const double rB = fb.Geometry.Radius;
    // The squared distance between the axes, in units of |a|²: |w x a|².
    const auto offset = ExactOffset(fb.Geometry.Origin, fa.Geometry.Origin);
    const auto squaredDistance = [&] {
        const auto cross = ExactCross(offset, fa.ExactAxis);
        return ExactDot(cross, cross);
    }();
    const auto scale = ExactDot(fa.ExactAxis, fa.ExactAxis);
    const Vec3 w = fb.Geometry.Origin - fa.Geometry.Origin;
    const Vec3 acrossVector = w - Dot(w, fa.Axis) * fa.Axis;
    const double distance = Length(acrossVector);
    const ExactScalar<2> radiiDifference = ExactDifference(rA, rB);
    if (WithinTolerance(squaredDistance, radiiDifference, scale, tolerance)) {
        const bool exact = SignOf(squaredDistance) == Sign::Zero && rA == rB;
        return Coincident(fa, boxA, fb, boxB, exact ? 0.0 : distance + std::abs(rA - rB), tolerance,
                          ctx);
    }
    const Sign sumSign = Compare(squaredDistance, ExactSum(rA, rB), scale);
    const Sign differenceSign = Compare(squaredDistance, radiiDifference, scale);
    if (sumSign == Sign::Positive || differenceSign == Sign::Negative) {
        return SurfaceSurfaceIntersection{}; // Proved separate: apart, or nested.
    }
    return Rulings(fa, fb, acrossVector / distance, distance, sumSign, differenceSign, boxA, boxB,
                   ctx);
}

} // namespace Stereon::Detail
