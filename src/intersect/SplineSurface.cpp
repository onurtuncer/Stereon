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
#include "BoxClip.hpp"
#include "ConicArc.hpp"
#include "Planar.hpp"

/// A B-spline curve, polynomial or rational, or a circle or an ellipse,
/// against a plane, a cylinder, a cone, a sphere, a torus or an ellipsoid:
/// the surface's implicit equation f = 0 composed with the curve. On each
/// Bezier segment,
/// in homogeneous coordinates (X, W) = (w P, w), W^k f(X / W) is a
/// polynomial g of degree p for the plane, 2p for the quadrics and 4p for
/// the torus, with the sign of f since every weight is positive. Its
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
///
/// A circle or an ellipse is searched on its exact rational form: its
/// bounded domain cut into arcs of at most a quarter turn, each a rational
/// quadratic Bezier segment enclosed from the input angles, and every root
/// mapped back to the conic's angle (ConicSurface).

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

/// a + b and a - b, coefficient by coefficient, for polynomials of one degree.
[[nodiscard]] Coefficients Plus(Coefficients a, const Coefficients& b, const double sign = 1.0)
{
    for (std::size_t k = 0; k < a.size(); ++k) {
        a[k] = a[k] + Interval{sign} * b[k];
    }
    return a;
}

/// s a, coefficient by coefficient.
[[nodiscard]] Coefficients Scaled(Coefficients a, const Interval& s)
{
    for (Interval& c : a) {
        c = s * c;
    }
    return a;
}

/// The local coordinate of X - W O along a direction of the surface's
/// frame: (X - W O) . d, of the curve's degree.
[[nodiscard]] Coefficients Along(const Offsets& o, const Vec3& d)
{
    Coefficients g;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        g.push_back(Interval{d.X} * o.X[i] + Interval{d.Y} * o.Y[i] + Interval{d.Z} * o.Z[i]);
    }
    return g;
}

/// The ellipsoid of LineEllipsoid in homogeneous local coordinates (P, Q, H)
/// of X - W C and the weight W: (B C P)² + (A C Q)² + (A B H)² - (A B C W)²,
/// degree 2p, (A B C W)² times the value of (x / A)² + (y / B)² + (z / C)² -
/// 1 at X / W. Its frame is the evaluator's, Axis(ellipsoid) as rounded, as
/// the line's handler takes it.
[[nodiscard]] Coefficients EllipsoidForm(const std::vector<HPoint>& poles,
                                         const Ellipsoid& ellipsoid)
{
    const Offsets o = OffsetsFrom(poles, ellipsoid.Center);
    const Interval a{ellipsoid.A};
    const Interval b{ellipsoid.B};
    const Interval c{ellipsoid.C};
    const Coefficients p = Scaled(Along(o, ellipsoid.XAxis), b * c);
    const Coefficients q = Scaled(Along(o, ellipsoid.YAxis), a * c);
    const Coefficients h = Scaled(Along(o, Axis(ellipsoid)), a * b);
    const Interval abc = a * b * c;
    return SumOfSquares(p, q, h, o.W, abc * abc);
}

/// The cone of Implicit.hpp's ImplicitValue in homogeneous local
/// coordinates (P, Q, H) of X - W O and the weight W:
/// cos²a (P² + Q²) - (R cos a W + sin a H)², degree 2p, W² times the value
/// at X / W. Its frame is the evaluator's, Axis(cone) as rounded, and sin a
/// and cos a are enclosed.
[[nodiscard]] Coefficients ConeForm(const std::vector<HPoint>& poles, const Cone& cone)
{
    const Offsets o = OffsetsFrom(poles, cone.Origin);
    const Coefficients p = Along(o, cone.XAxis);
    const Coefficients q = Along(o, cone.YAxis);
    const Coefficients h = Along(o, Axis(cone));
    const IntervalScalar angle{Interval{cone.SemiAngle}};
    const Interval c = cos(angle).Value;
    const Interval s = sin(angle).Value;
    const Coefficients k = Plus(Scaled(o.W, Interval{cone.Radius} * c), Scaled(h, s));
    const Coefficients radial = Plus(Product(p, p), Product(q, q));
    return Plus(Scaled(radial, c * c), Product(k, k), -1.0);
}

/// The torus of ImplicitValue, homogeneous: with K = P² + Q² + H² +
/// (R² - r²) W², K² - 4 R² (P² + Q²) W², degree 4p, W⁴ times the value.
[[nodiscard]] Coefficients TorusForm(const std::vector<HPoint>& poles, const Torus& torus)
{
    const Offsets o = OffsetsFrom(poles, torus.Center);
    const Coefficients p = Along(o, torus.XAxis);
    const Coefficients q = Along(o, torus.YAxis);
    const Coefficients h = Along(o, Axis(torus));
    const Interval big{torus.MajorRadius};
    const Interval small{torus.MinorRadius};
    const Coefficients ww = Product(o.W, o.W);
    const Coefficients radial = Plus(Product(p, p), Product(q, q));
    const Coefficients k = Plus(Plus(radial, Product(h, h)), Scaled(ww, big * big - small * small));
    return Plus(Product(k, k), Scaled(Product(radial, ww), Interval{4.0} * big * big), -1.0);
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

/// A surface parameter placed in its range: a periodic one modulo its
/// period, any other kept only inside.
[[nodiscard]] std::optional<double> Placed(const double x, const Interval& range,
                                           const bool periodic)
{
    if (periodic) {
        return InPeriodicDomain(x, range);
    }
    return Contains(range, x) ? std::optional<double>{x} : std::nullopt;
}

/// The roots polished once each. The subdivision places a decided
/// crossing to the rounding of the enclosed form, whose coefficients carry
/// the scale of its degree: for a torus's quartic, a point 1.2e-13 off the
/// surface, hundreds of units in the last place of its coordinates. One
/// Newton step on the signed distance from the surface along the curve,
/// (C(t) - S(uv)) . n over C'(t) . n at the projection uv of the rounded
/// point, places it to the point's own rounding. The step is kept only
/// inside the root's proved bracket, and skipped where the normal vanishes,
/// at a pole of the surface's parameters.
template <typename C, typename S>
[[nodiscard]] std::vector<Found> Polished(std::vector<Found> found, const C& curve,
                                          const S& surface)
{
    for (Found& f : found) {
        if (!f.Decided || f.Contact != ContactType::Crossing) {
            continue;
        }
        const Point3 p = Eval(curve, f.T);
        const Point2 uv = Project(surface, p);
        const Vec3 n = Cross(Deriv(surface, uv, 1, 0), Deriv(surface, uv, 0, 1));
        const double slope = Dot(Deriv(curve, f.T, 1), n);
        if (!(std::abs(slope) > 0.0) || !std::isfinite(slope)) { // STEREON-LINT-ALLOW(epsilon)
            continue;
        }
        const double t = f.T - Dot(p - Eval(surface, uv), n) / slope;
        if (std::isfinite(t) && Contains(f.Bracket, t)) {
            f.T = t;
        }
    }
    return found;
}

/// The points of the roots whose surface parameters lie in the box, a
/// periodic surface's angles placed there: u for every surface of
/// revolution, v as well for the torus. `point` is the curve's point at a
/// root's parameter.
template <typename Point, typename Parameters>
[[nodiscard]] CurveSurfaceIntersection PointsOf(Point&& point, const std::vector<Found>& found,
                                                const Box2& box, Parameters&& parameters,
                                                const bool periodicU, const bool periodicV = false)
{
    CurveSurfaceIntersection result;
    for (const Found& f : found) {
        const Point3 p = point(f.T);
        const Point2 uv = parameters(p);
        const std::optional<double> u = Placed(uv.X, box.X, periodicU);
        const std::optional<double> v = Placed(uv.Y, box.Y, periodicV);
        if (!u || !v) {
            continue;
        }
        result.Points.push_back({f.T, Point2{*u, *v}, p, f.Contact, f.Decided});
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

/// Where a point lies against one bound of the box: inside, outside, or
/// within a few roundings of it, where its side is not trusted.
[[nodiscard]] Strict SideOf(const double x, const double bound, const double inward)
{
    const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
                          std::max({1.0, std::abs(x), std::abs(bound)});
    const double d = inward * (x - bound);
    if (d > margin) {
        return Strict::Positive;
    }
    return d < -margin ? Strict::Negative : Strict::Unknown;
}

/// The axis whose edge a stretch runs along: x = const or y = const.
enum class EdgeAxis : std::uint8_t { None, X, Y };

/// Whether a point is in the closed box, Unknown within rounding of an edge.
/// A point of a stretch along an edge of `onEdge`'s axis is on that edge's
/// line, so only the other axis decides it.
[[nodiscard]] Strict InBox(const Point2& p, const Box2& box, const EdgeAxis onEdge = EdgeAxis::None)
{
    Strict side = Strict::Positive;
    const auto bound = [&](const double x, const double b, const double inward,
                           const bool boundary) {
        if (!std::isfinite(b)) {
            return;
        }
        const Strict s = SideOf(x, b, inward);
        if (s == Strict::Negative) {
            side = Strict::Negative;
        } else if (s == Strict::Unknown && !boundary && side == Strict::Positive) {
            side = Strict::Unknown;
        }
    };
    bound(p.X, box.X.GetLower(), 1.0, onEdge == EdgeAxis::X);
    bound(p.X, box.X.GetUpper(), -1.0, onEdge == EdgeAxis::X);
    bound(p.Y, box.Y.GetLower(), 1.0, onEdge == EdgeAxis::Y);
    bound(p.Y, box.Y.GetUpper(), -1.0, onEdge == EdgeAxis::Y);
    return side;
}

/// Whether a point found on an edge's line lies in the closed box: the
/// rounded point, its distance outside within rounding counted in.
[[nodiscard]] bool OnBoxBoundary(const Point2& p, const Box2& box)
{
    const auto within = [](const double x, const Interval& range) {
        const double slack =
            64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(x)});
        return x >= range.GetLower() - slack && x <= range.GetUpper() + slack;
    };
    return within(p.X, box.X) && within(p.Y, box.Y);
}

/// Where the cuts and the boundary stretches of a curve against the box's
/// edges are collected.
struct BoxCuts {
    std::vector<double> At;
    std::vector<std::pair<Interval, EdgeAxis>> Along;
    std::vector<double> Touches;
};

/// The pcurve against the line of one finite edge, of the given axis: every
/// contact, from a complete planar intersection, is a cut; a stretch along
/// the line is on the boundary where the other axis keeps it in the box; a
/// contact in the box is kept even where both sides of it leave.
[[nodiscard]] Result<void> CutAtEdge(const Curve2d& pcurve, const Interval& range,
                                     const Line2d& edge, const EdgeAxis axis, const Box2& box,
                                     const Context& ctx, BoxCuts& cuts)
{
    STEREON_TRY(const PlanarIntersection meet,
                Intersect(pcurve, range, Curve2d(edge), Interval::Everything(), ctx));
    for (const PlanarPoint& p : meet.Points) {
        cuts.At.push_back(p.ParameterA);
        if (OnBoxBoundary(p.Position, box)) {
            cuts.Touches.push_back(p.ParameterA);
        }
    }
    for (const PlanarOverlap& o : meet.Overlaps) {
        cuts.At.push_back(o.OnA.GetLower());
        cuts.At.push_back(o.OnA.GetUpper());
        cuts.Along.emplace_back(o.OnA, axis);
    }
    return {};
}

/// The cuts of the pcurve over the range against every finite edge's line,
/// sorted, each once, the range's ends among them.
[[nodiscard]] Result<BoxCuts> CutsOf(const Curve2d& pcurve, const Interval& range, const Box2& box,
                                     const Context& ctx)
{
    BoxCuts cuts;
    cuts.At = {range.GetLower(), range.GetUpper()};
    const auto edge = [&](const double at, const EdgeAxis axis) -> Result<void> {
        if (!std::isfinite(at)) {
            return {};
        }
        const Line2d line = axis == EdgeAxis::X ? Line2d{Point2{at, 0.0}, Vec2{0.0, 1.0}}
                                                : Line2d{Point2{0.0, at}, Vec2{1.0, 0.0}};
        return CutAtEdge(pcurve, range, line, axis, box, ctx, cuts);
    };
    STEREON_TRY_VOID(edge(box.X.GetLower(), EdgeAxis::X));
    STEREON_TRY_VOID(edge(box.X.GetUpper(), EdgeAxis::X));
    STEREON_TRY_VOID(edge(box.Y.GetLower(), EdgeAxis::Y));
    STEREON_TRY_VOID(edge(box.Y.GetUpper(), EdgeAxis::Y));
    for (double& t : cuts.At) {
        t = std::clamp(t, range.GetLower(), range.GetUpper());
    }
    std::ranges::sort(cuts.At);
    const auto repeated = std::ranges::unique(cuts.At);
    cuts.At.erase(repeated.begin(), repeated.end());
    return cuts;
}

/// The axis of the edge whose line a piece runs along, if any.
[[nodiscard]] EdgeAxis AlongEdge(const BoxCuts& cuts, const Interval& piece)
{
    EdgeAxis onEdge = EdgeAxis::None;
    for (const auto& [stretch, axis] : cuts.Along) {
        if (Contains(stretch, piece)) {
            onEdge = axis;
        }
    }
    return onEdge;
}

/// The stretches, sorted, with the ones that meet or overlap joined.
[[nodiscard]] std::vector<Interval> Joined(std::vector<Interval> kept)
{
    std::ranges::sort(
        kept, [](const Interval& a, const Interval& b) { return a.GetLower() < b.GetLower(); });
    std::vector<Interval> joined;
    for (const Interval& k : kept) {
        if (!joined.empty() && k.GetLower() <= joined.back().GetUpper()) {
            joined.back() = Hull(joined.back(), k);
        } else {
            joined.push_back(k);
        }
    }
    return joined;
}

/// The parameter stretches over which a curve in the plane stays in the
/// closed box. The curve's range is cut wherever its pcurve meets the line
/// of a finite edge. Each planar intersection is complete: it reports every
/// contact or fails. Between two cuts, then, the pcurve stays strictly on
/// one side of every edge's line, and its middle point decides the piece. A
/// piece along an edge's line is on the boundary, and is kept where the
/// other axis keeps its middle point in the box. A contact in the box
/// between two pieces outside is kept as a stretch of one parameter, as the
/// line handler keeps a corner. A middle point within rounding of an edge
/// cannot be trusted, and the clip fails. A range of one parameter has no
/// piece between cuts: its point decides it, kept within rounding of an
/// edge as a contact on the boundary is (OnBoxBoundary).
[[nodiscard]] Result<std::vector<Interval>> InsideBox(const Curve2d& pcurve, const Interval& range,
                                                      const Box2& box, const Context& ctx,
                                                      const std::string_view handler)
{
    if (range.GetLower() == range.GetUpper()) {
        if (InBox(Eval(pcurve, range.GetLower()), box) == Strict::Negative) {
            return std::vector<Interval>{};
        }
        return std::vector<Interval>{range};
    }
    STEREON_TRY(const BoxCuts cuts, CutsOf(pcurve, range, box, ctx));
    std::vector<Interval> kept;
    for (std::size_t i = 0; i + 1 < cuts.At.size(); ++i) {
        const Interval piece{cuts.At[i], cuts.At[i + 1]};
        const Point2 middle = Eval(pcurve, 0.5 * (piece.GetLower() + piece.GetUpper()));
        const Strict side = InBox(middle, box, AlongEdge(cuts, piece));
        if (side == Strict::Unknown) {
            return Fail(ErrorCode::IntersectionFailed,
                        std::string(handler) +
                            ": a curve in the plane runs within rounding of the box's edge");
        }
        if (side == Strict::Positive) {
            kept.push_back(piece);
        }
    }
    for (const double t : cuts.Touches) {
        const double at = std::clamp(t, range.GetLower(), range.GetUpper());
        if (std::ranges::none_of(kept, [&](const Interval& k) { return Contains(k, at); })) {
            kept.push_back(Interval{at});
        }
    }
    return Joined(std::move(kept));
}

/// The overlaps of a curve in the plane over the stretches in the box, all
/// on the one pcurve.
[[nodiscard]] CurveSurfaceIntersection OverlapsOn(const Curve2d& pcurve,
                                                  const std::vector<Interval>& stretches)
{
    CurveSurfaceIntersection result;
    for (const Interval& stretch : stretches) {
        result.Overlaps.push_back({stretch, pcurve});
    }
    return result;
}

/// The overlap of a curve in the plane: its pcurve is the spline of the
/// projected poles, a plane's parameters being affine. When every pole is
/// in the box the curve is, in the hull of its poles; otherwise it is
/// clipped to the box.
[[nodiscard]] Result<CurveSurfaceIntersection> AlongPlane(const BSplineCurve& spline,
                                                          const Plane& plane,
                                                          const Interval& domain, const Box2& box,
                                                          const Context& ctx)
{
    BSplineCurve2d onPlane{spline.Degree, spline.Knots, {}, spline.Weights};
    bool inside = true;
    for (const Point3& pole : spline.Poles) {
        const Point2 q = Project(plane, pole);
        inside = inside && Contains(box.X, q.X) && Contains(box.Y, q.Y);
        onPlane.Poles.push_back(q);
    }
    const Interval range = Intersection(domain, Domain(spline));
    if (range.IsEmpty()) {
        return CurveSurfaceIntersection{};
    }
    const Curve2d pcurve(std::move(onPlane));
    if (inside) {
        return OverlapsOn(pcurve, {range});
    }
    STEREON_TRY(const std::vector<Interval> stretches,
                InsideBox(pcurve, range, box, ctx, "BSplineCurve x Plane"));
    return OverlapsOn(pcurve, stretches);
}

} // namespace

Result<std::vector<Interval>> ClipPcurveToBox(const Curve2d& pcurve, const Interval& range,
                                              const Box2& box, const Context& ctx,
                                              const std::string_view handler)
{
    if (range.IsEmpty() || box.IsEmpty()) {
        return std::vector<Interval>{};
    }
    return InsideBox(pcurve, range, box, ctx, handler);
}

Result<CurveSurfaceIntersection> SplinePlane(const Curve& a, const Interval& domain,
                                             const Surface& b, const Box2& box, const Context& ctx)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Plane* plane = b.Target<Plane>();
    if (spline == nullptr || plane == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplinePlane: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    if (InPlane(*spline, *plane)) {
        return AlongPlane(*spline, *plane, domain, box, ctx);
    }
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return PlaneForm(poles, *plane); },
                    "BSplineCurve x Plane"));
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *plane), box,
                    [&](const Point3& p) { return Project(*plane, p); }, false);
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
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *sphere), box,
                    [&](const Point3& p) { return Project(*sphere, p); }, true);
}

Result<CurveSurfaceIntersection> SplineEllipsoid(const Curve& a, const Interval& domain,
                                                 const Surface& b, const Box2& box,
                                                 const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Ellipsoid* ellipsoid = b.Target<Ellipsoid>();
    if (spline == nullptr || ellipsoid == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineEllipsoid: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found, RootsOf(
                                                    *spline, domain,
                                                    [&](const std::vector<HPoint>& poles) {
                                                        return EllipsoidForm(poles, *ellipsoid);
                                                    },
                                                    "BSplineCurve x Ellipsoid"));
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *ellipsoid), box,
                    [&](const Point3& p) { return Project(*ellipsoid, p); }, true);
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
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *cylinder), box,
                    [&](const Point3& p) { return Project(*cylinder, p); }, true);
}

Result<CurveSurfaceIntersection> SplineCone(const Curve& a, const Interval& domain,
                                            const Surface& b, const Box2& box,
                                            const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Cone* cone = b.Target<Cone>();
    if (spline == nullptr || cone == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineCone: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return ConeForm(poles, *cone); },
                    "BSplineCurve x Cone"));
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *cone), box,
                    [&](const Point3& p) { return Project(*cone, p); }, true);
}

Result<CurveSurfaceIntersection> SplineTorus(const Curve& a, const Interval& domain,
                                             const Surface& b, const Box2& box,
                                             const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Torus* torus = b.Target<Torus>();
    if (spline == nullptr || torus == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineTorus: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return TorusForm(poles, *torus); },
                    "BSplineCurve x Torus"));
    return PointsOf([&](const double t) { return Eval(*spline, t); },
                    Polished(found, *spline, *torus), box,
                    [&](const Point3& p) { return Project(*torus, p); }, true, true);
}

// --- A line against the cone and the torus -------------------------------

namespace {

/// The line over [t0, t1] as one Bezier segment of degree 1, its poles
/// enclosed from the input doubles.
[[nodiscard]] Segment LineSegment(const Line& line, const double t0, const double t1)
{
    const auto at = [&](const double t) {
        const Interval s{t};
        return HPoint{Interval{line.Origin.X} + s * Interval{line.Direction.X},
                      Interval{line.Origin.Y} + s * Interval{line.Direction.Y},
                      Interval{line.Origin.Z} + s * Interval{line.Direction.Z}, Interval{1.0}};
    };
    return {t0, t1, {at(t0), at(t1)}};
}

/// Every root of the surface's form along the line over a bounded range,
/// found and proved as a B-spline's are: a root left undecided fails.
template <typename Form>
[[nodiscard]] Result<std::vector<Found>> LineRootsOf(const Line& line, const Interval& range,
                                                     Form&& form, const std::string_view handler)
{
    std::vector<Found> found;
    if (range.IsEmpty()) {
        return found;
    }
    const Segment segment = LineSegment(line, range.GetLower(), range.GetUpper());
    STEREON_TRY_VOID(RootsOfSegment(segment, form(segment.Poles), handler, found));
    Merge(found);
    for (const Found& f : found) {
        if (!f.Decided) {
            return Fail(ErrorCode::IntersectionFailed,
                        std::string(handler) +
                            ": unresolved root existence or multiplicity in requested domain");
        }
    }
    return found;
}

/// The enclosed coordinates of `p - c`, and of the line's direction.
[[nodiscard]] IVec3 Offset(const Point3& p, const Point3& c)
{
    return {Interval{p.X} - Interval{c.X}, Interval{p.Y} - Interval{c.Y},
            Interval{p.Z} - Interval{c.Z}};
}

[[nodiscard]] IVec3 Enclosed(const Vec3& v)
{
    return {Interval{v.X}, Interval{v.Y}, Interval{v.Z}};
}

[[nodiscard]] Interval DotOf(const IVec3& a, const IVec3& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

/// The stretch of the line within distance `radius` of c, outward-rounded:
/// every enclosure taken on the side that widens it, so that no point of
/// the line within the radius is left out; Empty when the line provably
/// passes outside.
[[nodiscard]] Interval WithinSphere(const Line& line, const Point3& c, const double radius)
{
    const IVec3 w = Offset(line.Origin, c);
    const IVec3 d = Enclosed(line.Direction);
    const Interval dd = DotOf(d, d);
    const IVec3 n{w.Y * d.Z - w.Z * d.Y, w.Z * d.X - w.X * d.Z, w.X * d.Y - w.Y * d.X};
    // The squared distance of the line from c, |w x d|² / |d|².
    const Interval off = DotOf(n, n) / dd;
    const Interval r{radius};
    const double room = (r * r - off).GetUpper();
    if (!(room > 0.0)) {
        return Interval::Empty();
    }
    const double half = std::nextafter(std::sqrt((Interval{room} / dd).GetUpper()), Infinity);
    const Interval centre = (Interval{0.0} - DotOf(w, d)) / dd;
    return {std::nextafter(centre.GetLower() - half, -Infinity),
            std::nextafter(centre.GetUpper() + half, Infinity)};
}

/// A t² + 2 B t + C, enclosed.
struct EnclosedQuadratic {
    Interval A;
    Interval B;
    Interval C;
};

/// The quadratic the line through a point at the enclosed offset `w` from
/// the cone's origin, with the line's direction, cuts out of the cone,
/// enclosed. With (P, Q, H) the coordinates of w and (Pd, Qd, Hd) of the
/// direction in the cone's frame, and K = R cos a + H sin a, the cone's
/// cos² a (p² + q²) - (R cos a + h sin a)² gives A = cos² a (Pd² + Qd²) -
/// sin² a Hd², B = cos² a (P Pd + Q Qd) - K sin a Hd, and C = cos² a (P² +
/// Q²) - K².
[[nodiscard]] EnclosedQuadratic ConeQuadraticAt(const IVec3& w, const Line& line, const Cone& cone)
{
    const IVec3 d = Enclosed(line.Direction);
    const IVec3 x = Enclosed(cone.XAxis);
    const IVec3 y = Enclosed(cone.YAxis);
    const IVec3 h = Enclosed(Axis(cone));
    const IntervalScalar angle{Interval{cone.SemiAngle}};
    const Interval c = cos(angle).Value;
    const Interval s = sin(angle).Value;
    const Interval p = DotOf(w, x);
    const Interval q = DotOf(w, y);
    const Interval pd = DotOf(d, x);
    const Interval qd = DotOf(d, y);
    const Interval k = Interval{cone.Radius} * c + s * DotOf(w, h);
    const Interval kd = s * DotOf(d, h);
    const Interval cc = c * c;
    return {cc * (pd * pd + qd * qd) - kd * kd, cc * (p * pd + q * qd) - k * kd,
            cc * (p * p + q * q) - k * k};
}

/// The quadratic of the line through `origin`, a point as given.
[[nodiscard]] EnclosedQuadratic ConeQuadratic(const Line& line, const Point3& origin,
                                              const Cone& cone)
{
    return ConeQuadraticAt(Offset(origin, cone.Origin), line, cone);
}

/// The quadratic about the line's point at the parameter `t0`, that point's
/// offset from the cone's origin enclosed, so that what it proves holds for
/// the line as given. Its roots are the line's own less t0 and its
/// discriminant is the same number, but formed from the offset of a point
/// near the crossings. From an origin far from them the squares of that
/// distance cancel in the discriminant and leave its enclosure across zero,
/// however far apart the crossings are.
[[nodiscard]] EnclosedQuadratic ConeQuadraticAbout(const Line& line, const double t0,
                                                   const Cone& cone)
{
    const IVec3 w = Offset(line.Origin, cone.Origin);
    const IVec3 d = Enclosed(line.Direction);
    const Interval t{t0};
    return ConeQuadraticAt({w.X + t * d.X, w.Y + t * d.Y, w.Z + t * d.Z}, line, cone);
}

[[nodiscard]] double Middle(const Interval& x)
{
    return 0.5 * (x.GetLower() + x.GetUpper());
}

[[nodiscard]] LineQuadratic MiddleOf(const EnclosedQuadratic& q)
{
    return {Middle(q.A), Middle(q.B), Middle(q.C)};
}

/// A crossing of the cone placed again from the line re-centred at its
/// rounded point, where the quadratic is formed from a small offset and its
/// root nearest zero is placed to the point's own rounding; a step reaching
/// half way to the other crossing is refused.
[[nodiscard]] double PolishedOnCone(const Line& line, const Cone& cone, const double t,
                                    const double other)
{
    const LineQuadratic q = MiddleOf(ConeQuadratic(line, Eval(line, t), cone));
    if (q.A == 0.0) {
        return t;
    }
    const std::array<double, 2> roots = RootsOf(q, std::max(0.0, q.B * q.B - q.A * q.C));
    const double s = std::abs(roots[0]) <= std::abs(roots[1]) ? roots[0] : roots[1];
    const double polished = t + s;
    return std::abs(s) < 0.5 * std::abs(other - t) && std::isfinite(polished) ? polished : t;
}

/// What an unbounded domain gives no way to resolve.
[[nodiscard]] Result<std::vector<Found>> UnresolvedOnCone()
{
    return Fail(ErrorCode::IntersectionFailed,
                "Line x Cone: a line along a generator, through the apex or tangent to the cone "
                "within rounding cannot be resolved over an unbounded domain; a bounded one can "
                "be");
}

/// The crossings of the cone over an unbounded domain, from its quadratic.
/// The count is decided on enclosures, in two steps. From the line's
/// origin: a discriminant below zero is a proved separation, and a leading
/// coefficient of no strict sign, a line along a generator within rounding,
/// fails. Then about the middle of the two roots, -B / A as rounded
/// (ConeQuadraticAbout), so that the decision does not depend on where
/// along the line its origin lies: a discriminant below zero is a proved
/// separation; one above zero, two simple crossings, each placed in doubles
/// from the middle of that enclosure and polished. Else a line through the
/// apex or tangent within rounding fails.
[[nodiscard]] Result<std::vector<Found>> ConeQuadraticRoots(const Line& line, const Cone& cone,
                                                            const Interval& domain)
{
    const EnclosedQuadratic far = ConeQuadratic(line, line.Origin, cone);
    if ((far.B * far.B - far.A * far.C).GetUpper() < 0.0) {
        return std::vector<Found>{};
    }
    if (!(far.A.GetLower() > 0.0 || far.A.GetUpper() < 0.0)) {
        return UnresolvedOnCone();
    }
    const double middle = -Middle(far.B) / Middle(far.A);
    if (!std::isfinite(middle)) {
        return UnresolvedOnCone();
    }
    const EnclosedQuadratic q = ConeQuadraticAbout(line, middle, cone);
    const Interval discriminant = q.B * q.B - q.A * q.C;
    if (discriminant.GetUpper() < 0.0) {
        return std::vector<Found>{};
    }
    if (!(discriminant.GetLower() > 0.0)) {
        return UnresolvedOnCone();
    }
    const LineQuadratic rounded = MiddleOf(q);
    const std::array<double, 2> roots =
        RootsOf(rounded, std::max(0.0, rounded.B * rounded.B - rounded.A * rounded.C));
    std::vector<Found> found;
    for (std::size_t k = 0; k < 2; ++k) {
        const double t = PolishedOnCone(line, cone, middle + roots[k], middle + roots[1 - k]);
        if (std::isfinite(t) && Contains(domain, t)) {
            found.push_back({t, Interval{t}, ContactType::Crossing, true});
        }
    }
    std::ranges::sort(found, [](const Found& a, const Found& b) { return a.T < b.T; });
    return found;
}

} // namespace

Result<CurveSurfaceIntersection> LineCone(const Curve& a, const Interval& domain, const Surface& b,
                                          const Box2& box, const Context& /*ctx*/)
{
    const Line* line = a.Target<Line>();
    const Cone* cone = b.Target<Cone>();
    if (line == nullptr || cone == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineCone: dispatched for other types");
    }
    std::vector<Found> found;
    if (domain.IsBounded()) {
        STEREON_TRY(found,
                    LineRootsOf(
                        *line, domain,
                        [&](const std::vector<HPoint>& poles) { return ConeForm(poles, *cone); },
                        "Line x Cone"));
    } else {
        STEREON_TRY(found, ConeQuadraticRoots(*line, *cone, domain));
    }
    return PointsOf([&](const double t) { return Eval(*line, t); }, Polished(found, *line, *cone),
                    box, [&](const Point3& p) { return Project(*cone, p); }, true);
}

Result<CurveSurfaceIntersection> LineTorus(const Curve& a, const Interval& domain, const Surface& b,
                                           const Box2& box, const Context& /*ctx*/)
{
    const Line* line = a.Target<Line>();
    const Torus* torus = b.Target<Torus>();
    if (line == nullptr || torus == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineTorus: dispatched for other types");
    }
    // Every point of the torus lies within R + r of its centre: the line
    // need be searched only there, a bounded stretch whatever the domain.
    // The sphere is taken half a tube wider, so that the stretch's ends lie
    // well off the torus: a crossing on its outer equator, where the sphere
    // touches it, would otherwise sit within rounding of an end, whose sign
    // could not be decided.
    const double reach =
        (Interval{torus->MajorRadius} + Interval{1.5} * Interval{torus->MinorRadius}).GetUpper();
    const Interval range = Intersection(domain, WithinSphere(*line, torus->Center, reach));
    STEREON_TRY(const std::vector<Found> found,
                LineRootsOf(
                    *line, range,
                    [&](const std::vector<HPoint>& poles) { return TorusForm(poles, *torus); },
                    "Line x Torus"));
    return PointsOf([&](const double t) { return Eval(*line, t); }, Polished(found, *line, *torus),
                    box, [&](const Point3& p) { return Project(*torus, p); }, true, true);
}

// --- Circles and ellipses ------------------------------------------------

namespace {

/// The most quarter-turn arcs a conic's requested domain may be cut into:
/// 1 024 turns.
constexpr std::size_t MaxArcs = 4096;

/// The parabola's stretch [t0, t1] of V + t^2 A + t B (A = X / (4 f)), a
/// polynomial quadratic, as its Bezier segment on s in [0, 1]: the poles
/// P(t0), V + t0 t1 A + (t0 + t1) / 2 B and P(t1), weight 1, each enclosed
/// from the doubles t0 and t1, so the segment encloses the parabola as
/// given.
[[nodiscard]] Segment ParabolicSegment(const ConicData& parabola, const double t0, const double t1)
{
    const Interval a{t0};
    const Interval b{t1};
    const Interval one{1.0};
    return {0.0,
            1.0,
            {ConicPole(parabola, one, a * a, a),
             ConicPole(parabola, one, a * b, (a + b) / Interval{2.0}),
             ConicPole(parabola, one, b * b, b)}};
}

/// The parabola's parameter at the stretch's Bezier parameter s: linear,
/// exact at the stretch's ends.
[[nodiscard]] double ParabolicAt(const double t0, const double t1, const double s)
{
    if (s <= 0.0) {
        return t0;
    }
    if (s >= 1.0) {
        return t1;
    }
    return std::clamp(t0 + (t1 - t0) * s, t0, t1);
}

/// How a conic's requested domain becomes arcs: the conic as C + c(t) A +
/// s(t) B, the widest arc, an arc's rational quadratic Bezier segment, and
/// the conic's parameter at an arc's Bezier parameter (monotone, exact at
/// the arc's ends). Circles and ellipses take quarter turns (ArcSegment,
/// AngleAt), hyperbolas stretches of one unit of t (HyperbolicSegment,
/// HyperbolicAt; cosh h never vanishes), a parabola, a polynomial, one
/// stretch for the whole domain.
struct ConicChart {
    ConicData Data;
    double Span;
    Segment (*Arc)(const ConicData&, double, double);
    double (*At)(double, double, double);
};

/// The arc's root, from its Bezier parameter to the conic's parameter: the
/// stretch maps through the chart's monotone map, widened by a few units
/// of its rounding and kept in the arc.
[[nodiscard]] Found OnConic(const Found& f, const ConicChart& chart, const double t0,
                            const double t1)
{
    const auto widened = [&](const double x, const double toward) {
        const double step = 4.0 * std::numeric_limits<double>::epsilon() * (1.0 + std::abs(x));
        return std::clamp(toward < x ? x - step : x + step, t0, t1);
    };
    const double lower = f.Bracket.GetLower() <= 0.0
                             ? t0
                             : widened(chart.At(t0, t1, f.Bracket.GetLower()), t0 - 1.0);
    const double upper = f.Bracket.GetUpper() >= 1.0
                             ? t1
                             : widened(chart.At(t0, t1, f.Bracket.GetUpper()), t1 + 1.0);
    return {chart.At(t0, t1, f.T), Interval{lower, upper}, f.Contact, f.Decided, f.Isolated};
}

/// The roots of the surface's form on the arc [t0, t1], as angles.
template <typename Form>
[[nodiscard]] Result<void> ArcRoots(const ConicChart& chart, const double t0, const double t1,
                                    Form& form, const std::string_view handler,
                                    std::vector<Found>& found)
{
    const Segment arc = chart.Arc(chart.Data, t0, t1);
    std::vector<Found> local;
    STEREON_TRY_VOID(RootsOfSegment(arc, form(arc.Poles), handler, local));
    for (const Found& f : local) {
        found.push_back(OnConic(f, chart, t0, t1));
    }
    return {};
}

/// Whether a root found again across moved cuts settles its region: inside
/// it, and decided or clear of every moved cut.
[[nodiscard]] bool SettlesRegion(const Found& f, const Interval& region,
                                 const std::vector<double>& moved)
{
    return Contains(region, f.Bracket) &&
           (f.Decided ||
            std::ranges::none_of(moved, [&](const double c) { return Contains(f.Bracket, c); }));
}

/// The fraction of a lone arc at which AcrossCuts cuts it: the golden
/// section, 1 / phi = phi - 1 = (sqrt 5 - 1) / 2 (the subtraction exact in
/// doubles), far from every dyadic fraction of the arc it splits.
constexpr double GoldenCut = std::numbers::phi - 1.0;

/// A root undecided because it lies at or within rounding of where the
/// search cut the conic: a cut between arcs, where each arc brackets it
/// across its end (an arc monotone throughout brackets it with the whole
/// arc, and Merge joins such brackets), or a dyadic point of an arc's
/// subdivision, where each half brackets it across its end. It is searched
/// again over every arc its bracket reaches, the region [cuts.front(),
/// cuts.back()], cut elsewhere: each interior cut c moved to c + w/4, w the
/// narrowest arc, or a lone arc cut at its golden section. The conic is
/// smooth across every such point, which is the kernel's choice and no
/// feature of the curve, and in the new arcs the root lies off their
/// dyadic points. The new arcs are at most five quarters of the old (for a
/// circle, short of the half turn where a rational quadratic's middle
/// weight vanishes). The new search is complete over the region, so its
/// roots replace every root whose bracket lies in the region; when one of
/// them leaves the region or is undecided across a new cut, the region is
/// left as it was and the intersection fails closed. A tangency stays
/// undecided wherever the cuts are, and fails closed.
template <typename Form>
[[nodiscard]] Result<void> AcrossCuts(const ConicChart& chart, const std::vector<double>& cuts,
                                      Form& form, const std::string_view handler,
                                      std::vector<Found>& found)
{
    double w = cuts.back() - cuts.front();
    for (std::size_t k = 1; k < cuts.size(); ++k) {
        w = std::min(w, cuts[k] - cuts[k - 1]);
    }
    std::vector<double> moved(cuts.begin() + 1, cuts.end() - 1);
    for (double& c : moved) {
        c += 0.25 * w;
    }
    if (moved.empty()) {
        moved.push_back(cuts.front() + GoldenCut * w);
    }
    std::vector<double> ends{cuts.front()};
    ends.insert(ends.end(), moved.begin(), moved.end());
    ends.push_back(cuts.back());
    std::vector<Found> again;
    for (std::size_t k = 0; k + 1 < ends.size(); ++k) {
        STEREON_TRY_VOID(ArcRoots(chart, ends[k], ends[k + 1], form, handler, again));
    }
    Merge(again);
    const Interval region{cuts.front(), cuts.back()};
    if (!std::ranges::all_of(again,
                             [&](const Found& f) { return SettlesRegion(f, region, moved); })) {
        return {};
    }
    std::erase_if(found, [&](const Found& f) { return Contains(region, f.Bracket); });
    found.insert(found.end(), again.begin(), again.end());
    Merge(found);
    return {};
}

/// The first and last cut of the arcs an undecided root's bracket reaches;
/// nothing for a decided root.
[[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>>
StuckRegion(const Found& f, const std::vector<double>& cuts)
{
    if (f.Decided) {
        return std::nullopt;
    }
    std::size_t i = 0;
    while (i + 1 < cuts.size() && cuts[i + 1] <= f.Bracket.GetLower()) {
        ++i;
    }
    std::size_t j = cuts.size() - 1;
    while (j > i + 1 && cuts[j - 1] >= f.Bracket.GetUpper()) {
        --j;
    }
    return std::pair{i, j};
}

/// The regions of the stuck roots, those that share an arc or a cut joined
/// into one (a cut that ends one region would stay a cut of the other's
/// search, where its root would be stuck again), as runs of cuts.
[[nodiscard]] std::vector<std::vector<double>> StuckRegions(const std::vector<Found>& found,
                                                            const std::vector<double>& cuts)
{
    std::vector<std::pair<std::size_t, std::size_t>> spans;
    for (const Found& f : found) {
        if (const std::optional<std::pair<std::size_t, std::size_t>> s = StuckRegion(f, cuts)) {
            spans.push_back(*s);
        }
    }
    std::ranges::sort(spans);
    std::vector<std::pair<std::size_t, std::size_t>> joined;
    for (const auto& s : spans) {
        if (!joined.empty() && s.first <= joined.back().second) {
            joined.back().second = std::max(joined.back().second, s.second);
        } else {
            joined.push_back(s);
        }
    }
    std::vector<std::vector<double>> regions;
    for (const auto& [i, j] : joined) {
        regions.emplace_back(cuts.begin() + static_cast<std::ptrdiff_t>(i),
                             cuts.begin() + static_cast<std::ptrdiff_t>(j) + 1);
    }
    return regions;
}

/// Every root of the surface's form on the conic over the domain, cut into
/// the chart's arcs (quarter turns of a circle or an ellipse), as the
/// conic's parameters, merged and placed as RootsOf places a B-spline's; a
/// root undecided at a cut between arcs is searched again across it
/// (AcrossCuts).
template <typename Form>
[[nodiscard]] Result<std::vector<Found>> ConicRootsOf(const ConicChart& chart,
                                                      const Interval& domain, Form&& form,
                                                      const std::string_view handler)
{
    if (!domain.IsBounded()) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(handler) + ": a conic's domain must be bounded");
    }
    const double width = Width(domain);
    const double count = std::max(1.0, std::ceil(width / chart.Span));
    if (count > static_cast<double>(MaxArcs)) {
        return Fail(ErrorCode::NotSupported,
                    std::string(handler) + ": a domain of more than 4 096 arcs");
    }
    const auto arcs = static_cast<std::size_t>(count);
    std::vector<double> cuts{domain.GetLower()};
    for (std::size_t k = 1; k < arcs; ++k) {
        cuts.push_back(domain.GetLower() + width * static_cast<double>(k) / count);
    }
    cuts.push_back(domain.GetUpper());
    std::vector<Found> found;
    for (std::size_t k = 0; k < arcs; ++k) {
        STEREON_TRY_VOID(ArcRoots(chart, cuts[k], cuts[k + 1], form, handler, found));
    }
    Merge(found);
    for (const std::vector<double>& region : StuckRegions(found, cuts)) {
        STEREON_TRY_VOID(AcrossCuts(chart, region, form, handler, found));
    }
    for (const Found& f : found) {
        if (!f.Decided) {
            return Fail(ErrorCode::IntersectionFailed,
                        std::string(handler) +
                            ": unresolved root existence or multiplicity in requested domain");
        }
    }
    return found;
}

/// Whether the conic lies exactly in the plane: its centre at height zero
/// and both its axes across the normal, decided on the input doubles.
[[nodiscard]] bool ConicInPlane(const Point3& center, const Vec3& x, const Vec3& y,
                                const Plane& plane)
{
    const auto zero = [&](const auto& value) {
        return SignOfFormula([&](auto f) {
                   return Dot3(f, Cross3(f, Of3(f, plane.XAxis), Of3(f, plane.YAxis)), value(f));
               }) == Sign::Zero;
    };
    return zero([&](auto f) { return Diff3(f, center, plane.Origin); }) &&
           zero([&](auto f) { return Of3(f, x); }) && zero([&](auto f) { return Of3(f, y); });
}

[[nodiscard]] Vec2 InPlane2(const Plane& plane, const Vec3& v)
{
    return {Dot(v, plane.XAxis), Dot(v, plane.YAxis)};
}

/// The overlap of a conic in the plane over the requested domain: its
/// pcurve the same conic in the plane's parameters, on the same angle. A
/// conic whose bounding square is in the box is kept whole; otherwise it is
/// clipped to the box.
[[nodiscard]] Result<CurveSurfaceIntersection>
ConicAlongPlane(const Curve2d& pcurve, const Point2& center, const double reach,
                const Interval& domain, const Box2& box, const Context& ctx,
                const std::string_view handler)
{
    const Box2 square{Around(center.X, reach), Around(center.Y, reach)};
    if (Contains(box.X, square.X) && Contains(box.Y, square.Y)) {
        return OverlapsOn(pcurve, {domain});
    }
    STEREON_TRY(const std::vector<Interval> stretches,
                InsideBox(pcurve, domain, box, ctx, handler));
    return OverlapsOn(pcurve, stretches);
}

/// The conic against one elementary surface: its form, its parameters and
/// which of them are periodic.
template <typename Conic, typename S, typename Form>
[[nodiscard]] Result<CurveSurfaceIntersection>
ConicAgainst(const Conic& conic, const ConicChart& data, const Interval& domain, const S& surface,
             const Box2& box, Form&& form, const std::string_view handler, const bool periodicU,
             const bool periodicV)
{
    STEREON_TRY(const std::vector<Found> found,
                ConicRootsOf(
                    data, domain,
                    [&](const std::vector<HPoint>& poles) { return form(poles, surface); },
                    handler));
    return PointsOf([&](const double t) { return Eval(conic, t); }, Polished(found, conic, surface),
                    box, [&](const Point3& p) { return Project(surface, p); }, periodicU,
                    periodicV);
}

/// The conic in the plane's parameters: the same conic on the same angle,
/// its centre and axes projected.
[[nodiscard]] Curve2d InPlaneOf(const Circle& circle, const Plane& plane)
{
    return Curve2d(Circle2d{Project(plane, circle.Center), circle.Radius,
                            InPlane2(plane, circle.XAxis), InPlane2(plane, circle.YAxis)});
}

[[nodiscard]] Curve2d InPlaneOf(const Ellipse& ellipse, const Plane& plane)
{
    return Curve2d(Ellipse2d{Project(plane, ellipse.Center), ellipse.MajorRadius,
                             ellipse.MinorRadius, InPlane2(plane, ellipse.XAxis),
                             InPlane2(plane, ellipse.YAxis)});
}

[[nodiscard]] Curve2d InPlaneOf(const Hyperbola& hyperbola, const Plane& plane)
{
    return Curve2d(Hyperbola2d{Project(plane, hyperbola.Center), hyperbola.MajorRadius,
                               hyperbola.MinorRadius, InPlane2(plane, hyperbola.XAxis),
                               InPlane2(plane, hyperbola.YAxis)});
}

[[nodiscard]] Curve2d InPlaneOf(const Parabola& parabola, const Plane& plane)
{
    return Curve2d(Parabola2d{Project(plane, parabola.Vertex), parabola.FocalLength,
                              InPlane2(plane, parabola.XAxis), InPlane2(plane, parabola.YAxis)});
}

/// How far the conic reaches from its centre: a hyperbola or a parabola
/// has no bound, so in a plane it is always clipped to the box.
[[nodiscard]] double Reach(const Circle& circle)
{
    return circle.Radius;
}

[[nodiscard]] double Reach(const Ellipse& ellipse)
{
    return ellipse.MajorRadius;
}

[[nodiscard]] double Reach(const Hyperbola& /*hyperbola*/)
{
    return std::numeric_limits<double>::infinity();
}

[[nodiscard]] double Reach(const Parabola& /*parabola*/)
{
    return std::numeric_limits<double>::infinity();
}

template <typename Conic>
[[nodiscard]] Point3 CentreOf(const Conic& conic)
{
    return conic.Center;
}

[[nodiscard]] Point3 CentreOf(const Parabola& parabola)
{
    return parabola.Vertex;
}

/// Dispatch on the surface for one conic.
template <typename Conic>
[[nodiscard]] Result<CurveSurfaceIntersection>
ConicOnSurface(const Conic& conic, const ConicChart& data, const Interval& domain, const Surface& b,
               const Box2& box, const Context& ctx, const std::string_view name)
{
    const std::string prefix(name);
    if (!domain.IsBounded()) {
        return Fail(ErrorCode::InvalidInput, prefix + ": a conic's domain must be bounded");
    }
    if (const Plane* plane = b.Target<Plane>()) {
        if (ConicInPlane(CentreOf(conic), conic.XAxis, conic.YAxis, *plane)) {
            return ConicAlongPlane(InPlaneOf(conic, *plane), Project(*plane, CentreOf(conic)),
                                   Reach(conic), domain, box, ctx, prefix + " x Plane");
        }
        return ConicAgainst(conic, data, domain, *plane, box, PlaneForm, prefix + " x Plane", false,
                            false);
    }
    if (const Cylinder* cylinder = b.Target<Cylinder>()) {
        return ConicAgainst(conic, data, domain, *cylinder, box, CylinderForm,
                            prefix + " x Cylinder", true, false);
    }
    if (const Cone* cone = b.Target<Cone>()) {
        return ConicAgainst(conic, data, domain, *cone, box, ConeForm, prefix + " x Cone", true,
                            false);
    }
    if (const Sphere* sphere = b.Target<Sphere>()) {
        return ConicAgainst(conic, data, domain, *sphere, box, SphereForm, prefix + " x Sphere",
                            true, false);
    }
    if (const Torus* torus = b.Target<Torus>()) {
        return ConicAgainst(conic, data, domain, *torus, box, TorusForm, prefix + " x Torus", true,
                            true);
    }
    if (const Ellipsoid* ellipsoid = b.Target<Ellipsoid>()) {
        return ConicAgainst(conic, data, domain, *ellipsoid, box, EllipsoidForm,
                            prefix + " x Ellipsoid", true, false);
    }
    return Fail(ErrorCode::InvalidInput, "ConicSurface: dispatched for other types");
}

/// The elementary surface's type, for messages.
[[nodiscard]] std::string SurfaceName(const Surface& b)
{
    if (b.Target<Plane>() != nullptr) {
        return "Plane";
    }
    if (b.Target<Cylinder>() != nullptr) {
        return "Cylinder";
    }
    if (b.Target<Cone>() != nullptr) {
        return "Cone";
    }
    if (b.Target<Sphere>() != nullptr) {
        return "Sphere";
    }
    if (b.Target<Torus>() != nullptr) {
        return "Torus";
    }
    return b.Target<Ellipsoid>() != nullptr ? "Ellipsoid" : "Surface";
}

/// A hyperbola's branch C + a cosh t X + b sinh t Y, or a parabola
/// V + t^2 / (4 f) X + t Y, against a surface: as a circle, on its exact
/// rational (hyperbola) or polynomial (parabola) form, over a bounded
/// requested domain. Neither curve is bounded, and an unbounded domain is
/// NotSupported rather than clipped to some finite stretch.
[[nodiscard]] Result<CurveSurfaceIntersection> OpenConicSurface(const Curve& a,
                                                                const Interval& domain,
                                                                const Surface& b, const Box2& box,
                                                                const Context& ctx)
{
    const Hyperbola* hyperbola = a.Target<Hyperbola>();
    const Parabola* parabola = a.Target<Parabola>();
    if (hyperbola == nullptr && parabola == nullptr) {
        return Fail(ErrorCode::InvalidInput, "ConicSurface: dispatched for other types");
    }
    if (!domain.IsBounded()) {
        return Fail(ErrorCode::NotSupported,
                    std::string(hyperbola != nullptr ? "Hyperbola" : "Parabola") + " x " +
                        SurfaceName(b) +
                        ": the requested domain must be bounded (the curve is not, and is "
                        "never clipped to a guessed stretch)");
    }
    if (hyperbola != nullptr) {
        const ConicData data{hyperbola->Center, Times(hyperbola->MajorRadius, hyperbola->XAxis),
                             Times(hyperbola->MinorRadius, hyperbola->YAxis)};
        return ConicOnSurface(*hyperbola, {data, 1.0, &HyperbolicSegment, &HyperbolicAt}, domain, b,
                              box, ctx, "Hyperbola");
    }
    const double f = parabola->FocalLength;
    if (!(f > 0.0) || !std::isfinite(4.0 * f)) {
        return Fail(ErrorCode::InvalidInput, "Parabola: the focal length must be positive");
    }
    const Interval four{4.0 * f};
    const ConicData data{parabola->Vertex,
                         {Interval{parabola->XAxis.X} / four, Interval{parabola->XAxis.Y} / four,
                          Interval{parabola->XAxis.Z} / four},
                         Times(1.0, parabola->YAxis)};
    return ConicOnSurface(
        *parabola, {data, std::numeric_limits<double>::infinity(), &ParabolicSegment, &ParabolicAt},
        domain, b, box, ctx, "Parabola");
}

} // namespace

Result<CurveSurfaceIntersection> ConicSurface(const Curve& a, const Interval& domain,
                                              const Surface& b, const Box2& box, const Context& ctx)
{
    const double quarter = 0.5 * std::numbers::pi;
    if (const Circle* circle = a.Target<Circle>()) {
        const ConicData data{circle->Center, Times(circle->Radius, circle->XAxis),
                             Times(circle->Radius, circle->YAxis)};
        return ConicOnSurface(*circle, {data, quarter, &ArcSegment, &AngleAt}, domain, b, box, ctx,
                              "Circle");
    }
    if (const Ellipse* ellipse = a.Target<Ellipse>()) {
        const ConicData data{ellipse->Center, Times(ellipse->MajorRadius, ellipse->XAxis),
                             Times(ellipse->MinorRadius, ellipse->YAxis)};
        return ConicOnSurface(*ellipse, {data, quarter, &ArcSegment, &AngleAt}, domain, b, box, ctx,
                              "Ellipse");
    }
    return OpenConicSurface(a, domain, b, box, ctx);
}

} // namespace Stereon::Detail
