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

/// Plane–plane (ADR-0010's analytic table, RS-0001 Tier A). The normals are
/// the cross products of the input axes, so every decision is a polynomial
/// in the input doubles: exactly parallel normals (their cross product is
/// zero) give a coincidence or a proved separation, decided by the signed
/// height of B's origin over A against the pair tolerance; anything else
/// gives a line, which stays a Line (ADR-0009).

namespace Stereon::Detail {

namespace {

struct PlaneFrame {
    const Plane& Geometry;
    Vec3 Normal;
    ExactVector<4> ExactNormal;
};

[[nodiscard]] PlaneFrame FrameOf(const Plane& plane)
{
    return {plane, Cross(plane.XAxis, plane.YAxis),
            ExactCross(Exactly(plane.XAxis), Exactly(plane.YAxis))};
}

/// uv on A to uv on B: Project(B, Eval(A, uv)), which is affine.
[[nodiscard]] ParameterMap MapBetween(const Plane& a, const Plane& b)
{
    ParameterMap map;
    map.U = Vec2{Dot(a.XAxis, b.XAxis), Dot(a.XAxis, b.YAxis)};
    map.V = Vec2{Dot(a.YAxis, b.XAxis), Dot(a.YAxis, b.YAxis)};
    map.Offset = ToVec(Project(b, a.Origin));
    return map;
}

/// Parallel planes within the pair tolerance: an overlap, exact when the
/// height is exactly zero and declared otherwise (ADR-0029 section 6).
[[nodiscard]] Result<SurfaceSurfaceIntersection>
Coincident(const PlaneFrame& a, const Box2& boxA, const PlaneFrame& b, const Box2& boxB,
           const Sign heightSign, const double tolerance, const Context& ctx)
{
    OverlapRegion overlap;
    overlap.AToB = MapBetween(a.Geometry, b.Geometry);
    overlap.OnA = Intersection(boxA, ImageBox(Inverse(overlap.AToB), boxB));
    overlap.OnB = Intersection(boxB, ImageBox(overlap.AToB, overlap.OnA));
    SurfaceSurfaceIntersection result;
    if (overlap.OnA.X.IsEmpty() || overlap.OnA.Y.IsEmpty()) {
        return result; // Coincident supports, but not within the domains.
    }
    overlap.FullDomain = overlap.OnA == boxA && overlap.OnB == boxB;
    overlap.SameOrientation = SignOf(ExactDot(a.ExactNormal, b.ExactNormal)) == Sign::Positive;
    if (heightSign != Sign::Zero) {
        overlap.Evidence = CoincidenceEvidence::Tolerant;
        overlap.Distance =
            std::abs(Dot(b.Geometry.Origin - a.Geometry.Origin, a.Normal)) / Length(a.Normal);
        ctx.Log(LogLevel::Info, "Intersect: exactly parallel planes " + Number(overlap.Distance) +
                                    " mm apart, within the pair tolerance " + Number(tolerance) +
                                    " mm: declared coincident, nothing moved (ADR-0029)");
    }
    result.Overlaps.push_back(std::move(overlap));
    return result;
}

/// Planes that are not parallel: their line, parametrised by arc length,
/// with a Line2d on each plane on the same parameter.
[[nodiscard]] Result<SurfaceSurfaceIntersection> Crossing(const Surface& sa, const PlaneFrame& a,
                                                          const Box2& boxA, const Surface& sb,
                                                          const PlaneFrame& b, const Box2& boxB,
                                                          const Context& ctx)
{
    const Vec3 direction = Normalized(Cross(a.Normal, b.Normal));
    // In plane A and across the line: where B's plane equation vanishes.
    const Vec3 across = Cross(direction, a.Normal);
    const double s = Dot(b.Geometry.Origin - a.Geometry.Origin, b.Normal) / Dot(across, b.Normal);
    const Point3 point = a.Geometry.Origin + s * across;
    const Line2d onA{Project(a.Geometry, point),
                     Vec2{Dot(direction, a.Geometry.XAxis), Dot(direction, a.Geometry.YAxis)}};
    const Line2d onB{Project(b.Geometry, point),
                     Vec2{Dot(direction, b.Geometry.XAxis), Dot(direction, b.Geometry.YAxis)}};
    const Interval range = Intersection(ClipLine(onA, boxA), ClipLine(onB, boxB));
    SurfaceSurfaceIntersection result;
    if (!range.IsEmpty()) {
        STEREON_TRY_VOID(
            Keep(result,
                 Detail::Crossing(Curve(Line{point, direction}), Curve2d(onA), Curve2d(onB), range),
                 sa, sb, ctx));
    }
    return result;
}

} // namespace

Result<SurfaceSurfaceIntersection> PlanePlane(const Surface& a, const Box2& boxA, const Surface& b,
                                              const Box2& boxB, const double tolerance,
                                              const Context& ctx)
{
    const Plane* pa = a.Target<Plane>();
    const Plane* pb = b.Target<Plane>();
    if (pa == nullptr || pb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "PlanePlane: dispatched for other surface types");
    }
    const PlaneFrame fa = FrameOf(*pa);
    const PlaneFrame fb = FrameOf(*pb);
    if (!IsExactlyZero(ExactCross(fa.ExactNormal, fb.ExactNormal))) {
        return Crossing(a, fa, boxA, b, fb, boxB, ctx);
    }
    // The signed height of B's origin over A, in units of |nA|.
    const auto height =
        ExactDot(ExactOffset(fb.Geometry.Origin, fa.Geometry.Origin), fa.ExactNormal);
    const auto scale = ExactDot(fa.ExactNormal, fa.ExactNormal);
    if (CompareSquared(height, Exactly(tolerance), scale) == Sign::Positive) {
        return SurfaceSurfaceIntersection{}; // Proved separate.
    }
    return Coincident(fa, boxA, fb, boxB, SignOf(height), tolerance, ctx);
}

} // namespace Stereon::Detail
