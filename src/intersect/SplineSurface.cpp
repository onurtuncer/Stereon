// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
#include "Bernstein.hpp"
#include "Planar.hpp"

/// A B-spline curve, polynomial or rational, against a plane, a sphere or a
/// cylinder: the surface's implicit equation f = 0 composed with the curve.
/// On each Bezier segment, in homogeneous coordinates (X, W) = (w P, w),
/// W^k f(X / W) is a polynomial g of degree p for the plane and 2p for the
/// quadrics, with the sign of f since every weight is positive. Its
/// Bernstein coefficients are formed in outward-rounded interval
/// arithmetic from the segment's poles, which are themselves enclosed:
/// the segments come from Boehm's knot insertion done in intervals, not
/// from the rounded ToBezierSegments, so that every proof below holds for
/// the curve as given (ADR-0029: a result with no component is a proved
/// separation).
///
/// The segment's parameter is subdivided at halves. A piece is dropped only
/// when g's coefficients are all strictly of one sign, which proves g has no
/// root there (the Bernstein convex hull). A piece whose derivative
/// coefficients are strictly of one sign holds at most one root, a simple
/// one, so a crossing: proved present and bisected to rounding when g's end
/// values have strictly opposite signs, reported as undecided when an end
/// value straddles zero. A piece whose second differences are strictly of
/// one sign holds a single extremum; it is split on around it, and the
/// sub-piece still holding it at MaxDepth is reported as an undecided
/// touch: a tangency, a near miss or two crossings within rounding.
/// Anything else, a contact of higher order, a curve along a quadric, or a
/// search past PieceBudget, fails with IntersectionFailed: an unresolved
/// piece is never dropped (ADR-0010, completeness). The segments and the
/// subdivision are Bernstein.hpp's, shared with the plane curve handlers;
/// here every root that is not a proved crossing fails.
///
/// A curve whose every pole lies exactly in the plane lies in it, an
/// overlap; along a quadric there is no such exact decision yet, and the
/// search fails. The surface parameters of a point are those of its rounded
/// position, kept or dropped against the box as the line handlers do.

namespace Stereon::Detail {

namespace {

[[nodiscard]] Result<void> CheckSpline(const BSplineCurve& curve)
{
    constexpr std::string_view where = "Intersect(BSplineCurve, surface)";
    STEREON_TRY_VOID(CheckKnotVector(curve.Degree, curve.Knots, curve.Poles.size(), where));
    STEREON_TRY_VOID(CheckWeights(curve.Weights, curve.Poles.size(), where));
    for (const Point3& pole : curve.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return {};
}

// --- The implicit equations in Bernstein form ----------------------------

/// The coordinates of X - W c over a segment's poles: the curve's offset
/// from c, homogeneous, coordinate by coordinate.
struct Offsets {
    Coefficients X;
    Coefficients Y;
    Coefficients Z;
    Coefficients W;
};

[[nodiscard]] Offsets OffsetsFrom(const std::vector<HPoint>& poles, const Point3& c)
{
    Offsets o;
    for (const HPoint& h : poles) {
        o.X.push_back(h.X - h.W * Interval{c.X});
        o.Y.push_back(h.Y - h.W * Interval{c.Y});
        o.Z.push_back(h.Z - h.W * Interval{c.Z});
        o.W.push_back(h.W);
    }
    return o;
}

/// a x b, enclosed.
struct IVec3 {
    Interval X;
    Interval Y;
    Interval Z;
};

[[nodiscard]] IVec3 CrossOf(const Vec3& a, const Vec3& b)
{
    const IVec3 u{Interval{a.X}, Interval{a.Y}, Interval{a.Z}};
    const IVec3 v{Interval{b.X}, Interval{b.Y}, Interval{b.Z}};
    return {u.Y * v.Z - u.Z * v.Y, u.Z * v.X - u.X * v.Z, u.X * v.Y - u.Y * v.X};
}

/// x² + y² + z² - s w², s enclosing r² (and |n|² for a cylinder).
[[nodiscard]] Coefficients SumOfSquares(const Coefficients& x, const Coefficients& y,
                                        const Coefficients& z, const Coefficients& w,
                                        const Interval& s)
{
    Coefficients g = Product(x, x);
    const Coefficients yy = Product(y, y);
    const Coefficients zz = Product(z, z);
    const Coefficients ww = Product(w, w);
    for (std::size_t k = 0; k < g.size(); ++k) {
        g[k] = g[k] + yy[k] + zz[k] - s * ww[k];
    }
    return g;
}

/// The plane through O with normal n = X x Y: n . (X - W O), degree p.
[[nodiscard]] Coefficients PlaneForm(const std::vector<HPoint>& poles, const Plane& plane)
{
    const Offsets o = OffsetsFrom(poles, plane.Origin);
    const IVec3 n = CrossOf(plane.XAxis, plane.YAxis);
    Coefficients g;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        g.push_back(n.X * o.X[i] + n.Y * o.Y[i] + n.Z * o.Z[i]);
    }
    return g;
}

/// The sphere: |X - W C|² - r² W², degree 2p.
[[nodiscard]] Coefficients SphereForm(const std::vector<HPoint>& poles, const Sphere& sphere)
{
    const Offsets o = OffsetsFrom(poles, sphere.Center);
    const Interval r{sphere.Radius};
    return SumOfSquares(o.X, o.Y, o.Z, o.W, r * r);
}

/// The cylinder about n = X x Y, taken as it is as LineCylinder takes it:
/// |(X - W O) x n|² - r² |n|² W², degree 2p.
[[nodiscard]] Coefficients CylinderForm(const std::vector<HPoint>& poles, const Cylinder& cylinder)
{
    const Offsets o = OffsetsFrom(poles, cylinder.Origin);
    const IVec3 n = CrossOf(cylinder.XAxis, cylinder.YAxis);
    Offsets e;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        e.X.push_back(o.Y[i] * n.Z - o.Z[i] * n.Y);
        e.Y.push_back(o.Z[i] * n.X - o.X[i] * n.Z);
        e.Z.push_back(o.X[i] * n.Y - o.Y[i] * n.X);
    }
    const Interval r{cylinder.Radius};
    const Interval scale = r * r * (n.X * n.X + n.Y * n.Y + n.Z * n.Z);
    return SumOfSquares(e.X, e.Y, e.Z, o.W, scale);
}

/// Every root of g over the segments that reach the domain, merged and
/// placed in the domain.
template <typename Form>
[[nodiscard]] Result<std::vector<Found>> RootsOf(const BSplineCurve& spline, const Interval& domain,
                                                 Form&& form, const std::string_view handler)
{
    STEREON_TRY_VOID(CheckPeriodicDomain(spline, domain, handler));
    STEREON_TRY(const std::vector<Segment> segments, SegmentsOf(spline));
    std::vector<Found> found;
    for (const Segment& segment : segments) {
        if (Intersects(Interval{segment.T0, segment.T1}, domain)) {
            STEREON_TRY_VOID(RootsOfSegment(segment, form(segment.Poles), handler, found));
        }
    }
    Merge(found);
    std::vector<Found> kept;
    for (const Found& f : found) {
        if (const std::optional<Found> placed = InDomain(f, domain)) {
            // An uncertain root enclosure may be a near miss or multiple
            // crossings. It is not a certified intersection component.
            if (!placed->Decided) {
                return Fail(ErrorCode::IntersectionFailed,
                            std::string(handler) +
                                ": unresolved root existence or multiplicity in requested domain");
            }
            kept.push_back(*placed);
        }
    }
    return kept;
}

/// The points of the roots whose surface parameters lie in the box, the
/// angle u of a periodic surface placed there.
template <typename Parameters>
[[nodiscard]] CurveSurfaceIntersection PointsOf(const BSplineCurve& spline,
                                                const std::vector<Found>& found, const Box2& box,
                                                Parameters&& parameters, const bool periodicU)
{
    CurveSurfaceIntersection result;
    for (const Found& f : found) {
        const Point3 p = Eval(spline, f.T);
        Point2 uv = parameters(p);
        const std::optional<double> u =
            periodicU ? InPeriodicDomain(uv.X, box.X)
                      : (Contains(box.X, uv.X) ? std::optional<double>{uv.X} : std::nullopt);
        if (!u || !Contains(box.Y, uv.Y)) {
            continue;
        }
        uv.X = *u;
        result.Points.push_back({f.T, uv, p, f.Contact, f.Decided});
    }
    return result;
}

// --- A curve in a plane --------------------------------------------------

/// Whether every pole lies exactly in the plane, n . (P - O) = 0 decided on
/// the input doubles; the curve is then in the plane, in the hull of its
/// poles with positive weights.
[[nodiscard]] bool InPlane(const BSplineCurve& spline, const Plane& plane)
{
    return std::ranges::all_of(spline.Poles, [&](const Point3& pole) {
        return SignOfFormula([&](auto f) {
                   return Dot3(f, Cross3(f, Of3(f, plane.XAxis), Of3(f, plane.YAxis)),
                               Diff3(f, pole, plane.Origin));
               }) == Sign::Zero;
    });
}

/// The overlap of a curve in the plane: its pcurve is the spline of the
/// projected poles, a plane's parameters being affine. Clipping it to a box
/// its poles leave is not implemented.
[[nodiscard]] Result<CurveSurfaceIntersection>
AlongPlane(const BSplineCurve& spline, const Plane& plane, const Interval& domain, const Box2& box)
{
    BSplineCurve2d pcurve{spline.Degree, spline.Knots, {}, spline.Weights};
    for (const Point3& pole : spline.Poles) {
        const Point2 q = Project(plane, pole);
        if (!Contains(box.X, q.X) || !Contains(box.Y, q.Y)) {
            return Fail(ErrorCode::NotImplemented,
                        "BSplineCurve x Plane: a curve in the plane that leaves the requested box "
                        "is not clipped yet");
        }
        pcurve.Poles.push_back(q);
    }
    CurveSurfaceIntersection result;
    const Interval range = Intersection(domain, Domain(spline));
    if (!range.IsEmpty()) {
        result.Overlaps.push_back({range, Curve2d(std::move(pcurve))});
    }
    return result;
}

} // namespace

Result<CurveSurfaceIntersection> SplinePlane(const Curve& a, const Interval& domain,
                                             const Surface& b, const Box2& box,
                                             const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Plane* plane = b.Target<Plane>();
    if (spline == nullptr || plane == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplinePlane: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    if (InPlane(*spline, *plane)) {
        return AlongPlane(*spline, *plane, domain, box);
    }
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return PlaneForm(poles, *plane); },
                    "BSplineCurve x Plane"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*plane, p); }, false);
}

Result<CurveSurfaceIntersection> SplineSphere(const Curve& a, const Interval& domain,
                                              const Surface& b, const Box2& box,
                                              const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Sphere* sphere = b.Target<Sphere>();
    if (spline == nullptr || sphere == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineSphere: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return SphereForm(poles, *sphere); },
                    "BSplineCurve x Sphere"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*sphere, p); }, true);
}

Result<CurveSurfaceIntersection> SplineCylinder(const Curve& a, const Interval& domain,
                                                const Surface& b, const Box2& box,
                                                const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Cylinder* cylinder = b.Target<Cylinder>();
    if (spline == nullptr || cylinder == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineCylinder: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found, RootsOf(
                                                    *spline, domain,
                                                    [&](const std::vector<HPoint>& poles) {
                                                        return CylinderForm(poles, *cylinder);
                                                    },
                                                    "BSplineCurve x Cylinder"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*cylinder, p); }, true);
}

} // namespace Stereon::Detail
