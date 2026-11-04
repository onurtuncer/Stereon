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
#include <cstdint>
#include <limits>
#include <memory>
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
#include "Bernstein.hpp"
#include "Planar.hpp"

/// A B-spline in the plane, polynomial or rational, against a line, a
/// circle, an ellipse, a hyperbola or a parabola (the Boolean prototype's
/// M6): the pcurve of a branch against a loop's edges. The other curve's
/// implicit equation f = 0 composed with the spline is, on each Bezier
/// segment and in homogeneous coordinates (X, Y, W) = (w P, w), a
/// polynomial g of degree p against a line and 2p against a conic, of the
/// sign of f since every weight is positive:
///
///   line       n . (X - W O),                n across the line;
///   circle     |X - W C|² - r² W²;
///   ellipse    b² xi² + a² eta² - a² b² W²,  (xi, eta) = (X - W C) . (X_B, Y_B);
///   hyperbola  b² xi² - a² eta² - a² b² W²;
///   parabola   eta² - 4 f xi W.
///
/// Its Bernstein coefficients are enclosed from the segments of
/// Bernstein.hpp, which come from Boehm's knot insertion done in intervals,
/// and its roots isolated by the same subdivision as the curve–surface
/// handlers' (SplineSurface.cpp): a piece is apart, a proved single simple
/// crossing (Decided), a touch or a crossing within rounding of a piece's
/// end (Decided false, bracketed), or unresolved. Unlike there, an
/// undecided root is a result here, reported with Decided false, as the
/// conic handlers report theirs; an unresolved one, a contact of higher
/// order or the spline along the other curve, fails with
/// IntersectionFailed, never dropped (ADR-0029). A hyperbola's equation
/// holds both branches: a root on the other branch, 2a away, is dropped by
/// the sign of its abscissa. A spline whose every pole lies exactly on a
/// line lies along it: an overlap, given when the line's requested domain
/// holds the whole spline, which is clamped and runs monotonically along
/// the line; any other such overlap is not clipped yet (NotImplemented).

namespace Stereon::Detail {

namespace {

// --- The other curve -----------------------------------------------------

enum class OtherKind : std::uint8_t { Line, Circle, Ellipse, Hyperbola, Parabola };

/// The other curve: a line (Origin, X its direction), a circle (Center,
/// radius A), or a conic in its frame (semi-axes A and B; a parabola's
/// vertex and focal length A).
struct Other {
    OtherKind Kind;
    Point2 Origin;
    Vec2 X;
    Vec2 Y;
    double A;
    double B;
};

[[nodiscard]] std::optional<Other> OtherOf(const Curve2d& curve)
{
    if (const Line2d* l = curve.Target<Line2d>()) {
        return Other{
            OtherKind::Line, l->Origin, l->Direction, Perpendicular(l->Direction), 0.0, 0.0};
    }
    if (const Circle2d* c = curve.Target<Circle2d>()) {
        return Other{OtherKind::Circle, c->Center, c->XAxis, c->YAxis, c->Radius, c->Radius};
    }
    if (const Ellipse2d* e = curve.Target<Ellipse2d>()) {
        return Other{OtherKind::Ellipse, e->Center,      e->XAxis,
                     e->YAxis,           e->MajorRadius, e->MinorRadius};
    }
    if (const Hyperbola2d* h = curve.Target<Hyperbola2d>()) {
        return Other{OtherKind::Hyperbola, h->Center,     h->XAxis, h->YAxis,
                     h->MajorRadius,       h->MinorRadius};
    }
    if (const Parabola2d* p = curve.Target<Parabola2d>()) {
        return Other{OtherKind::Parabola, p->Vertex, p->XAxis, p->YAxis, p->FocalLength, 0.0};
    }
    return std::nullopt;
}

/// The other curve's parameter of a point on it.
[[nodiscard]] double ParameterOf(const Other& o, const Point2& p)
{
    switch (o.Kind) {
    case OtherKind::Line:
        return Dot(p - o.Origin, o.X) / Dot(o.X, o.X);
    case OtherKind::Circle:
    case OtherKind::Ellipse:
        return AngleIn(o.Origin, o.X, o.Y, o.A, o.B, p);
    case OtherKind::Hyperbola:
        return std::asinh(Dot(p - o.Origin, o.Y) / o.B);
    case OtherKind::Parabola:
        return Dot(p - o.Origin, o.Y);
    }
    return 0.0;
}

[[nodiscard]] bool IsPeriodic(const Other& o)
{
    return o.Kind == OtherKind::Circle || o.Kind == OtherKind::Ellipse;
}

// --- The implicit equations in Bernstein form ----------------------------

/// X - W c over a segment's poles, coordinate by coordinate.
struct Offsets {
    Coefficients X;
    Coefficients Y;
    Coefficients W;
};

[[nodiscard]] Offsets OffsetsFrom(const std::vector<HPoint>& poles, const Point2& c)
{
    Offsets o;
    for (const HPoint& h : poles) {
        o.X.push_back(h.X - h.W * Interval{c.X});
        o.Y.push_back(h.Y - h.W * Interval{c.Y});
        o.W.push_back(h.W);
    }
    return o;
}

/// (X - W c) . e, coefficient by coefficient.
[[nodiscard]] Coefficients Along(const Offsets& o, const Vec2& e)
{
    Coefficients g;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        g.push_back(o.X[i] * Interval{e.X} + o.Y[i] * Interval{e.Y});
    }
    return g;
}

/// k1 p + k2 q + k3 r, coefficient by coefficient (all of one degree).
[[nodiscard]] Coefficients Sum3(const Interval& k1, const Coefficients& p, const Interval& k2,
                                const Coefficients& q, const Interval& k3, const Coefficients& r)
{
    Coefficients g;
    for (std::size_t i = 0; i < p.size(); ++i) {
        g.push_back(k1 * p[i] + k2 * q[i] + k3 * r[i]);
    }
    return g;
}

/// g on a segment's poles, as the table above gives it.
[[nodiscard]] Coefficients FormOf(const Other& o, const std::vector<HPoint>& poles)
{
    const Offsets offsets = OffsetsFrom(poles, o.Origin);
    const Coefficients ww = Product(offsets.W, offsets.W);
    const Interval a{o.A};
    const Interval b{o.B};
    switch (o.Kind) {
    case OtherKind::Line:
        return Along(offsets, o.Y);
    case OtherKind::Circle:
        return Sum3(Interval{1.0}, Product(offsets.X, offsets.X), Interval{1.0},
                    Product(offsets.Y, offsets.Y), -(a * a), ww);
    case OtherKind::Ellipse:
    case OtherKind::Hyperbola: {
        const Coefficients xi = Along(offsets, o.X);
        const Coefficients eta = Along(offsets, o.Y);
        const Interval sign{o.Kind == OtherKind::Ellipse ? 1.0 : -1.0};
        return Sum3(b * b, Product(xi, xi), sign * a * a, Product(eta, eta), -(a * a * b * b), ww);
    }
    case OtherKind::Parabola: {
        const Coefficients eta = Along(offsets, o.Y);
        const Coefficients xiW = Product(Along(offsets, o.X), offsets.W);
        const Coefficients zero(xiW.size(), Interval{0.0});
        return Sum3(Interval{1.0}, Product(eta, eta), Interval{-4.0} * a, xiW, Interval{0.0}, zero);
    }
    }
    return {};
}

// --- Roots ---------------------------------------------------------------

/// Every root of g over the segments that reach the domain, merged and
/// placed in it; an unresolved one fails.
[[nodiscard]] Result<std::vector<Found>> SplineRoots(const BSplineCurve2d& spline,
                                                     const Interval& domain, const Other& other,
                                                     const std::string_view handler)
{
    STEREON_TRY(const std::vector<Segment> segments, SegmentsOf(spline));
    std::vector<Found> found;
    for (const Segment& segment : segments) {
        if (Intersects(Interval{segment.T0, segment.T1}, domain)) {
            STEREON_TRY_VOID(RootsOfSegment(segment, FormOf(other, segment.Poles), handler, found));
        }
    }
    Merge(found);
    std::vector<Found> kept;
    for (const Found& f : found) {
        if (const std::optional<Found> placed = InDomain(f, domain)) {
            if (!placed->Isolated) {
                return Fail(ErrorCode::IntersectionFailed,
                            std::string(handler) +
                                ": neither a crossing nor a touch could be isolated near t = " +
                                Number(placed->T) +
                                " (a contact of higher order, or the spline along the other "
                                "curve; ADR-0010: it is never dropped)");
            }
            kept.push_back(*placed);
        }
    }
    return kept;
}

/// The point of a root, kept when the other curve's parameter lies in its
/// domain. A point on a hyperbola's other branch is not on it: its
/// abscissa is at most -a, 2a from the branch; one within a / 2 of the
/// axis would be on neither, and fails.
[[nodiscard]] Result<void> KeepRoot(PlanarIntersection& result, const BSplineCurve2d& spline,
                                    const Other& other, const Found& f, const Interval& domainA,
                                    const Interval& domainB)
{
    const Point2 p = Eval(spline, f.T);
    if (other.Kind == OtherKind::Hyperbola) {
        const double abscissa = Dot(p - other.Origin, other.X);
        if (std::abs(abscissa) < 0.5 * other.A) {
            return Fail(ErrorCode::IntersectionFailed,
                        "BSplineCurve2d x Hyperbola2d: a root off both branches");
        }
        if (abscissa < 0.0) {
            return {};
        }
    }
    KeepPoint(result, {f.T, ParameterOf(other, p), p, f.Contact, f.Decided}, domainA, false,
              domainB, IsPeriodic(other));
    return {};
}

// --- A spline along a line -----------------------------------------------

/// Whether every pole lies exactly on the line, d x (P - O) = 0 decided on
/// the input doubles: the spline then lies on it, in the hull of its poles.
[[nodiscard]] bool OnLine(const BSplineCurve2d& spline, const Other& line)
{
    return std::ranges::all_of(spline.Poles, [&](const Point2& pole) {
        return SignOfFormula([&](auto f) {
                   return Cross2(f, Of2(f, line.X), Diff2(f, pole, line.Origin));
               }) == Sign::Zero;
    });
}

/// The sign every step d . (P_i+1 - P_i) of the poles along the line
/// strictly has, exactly, or Zero: with positive weights the spline's
/// parameter on the line is then monotone (variation diminishing).
[[nodiscard]] Sign SenseAlong(const BSplineCurve2d& spline, const Vec2& d)
{
    std::optional<Sign> common;
    for (std::size_t i = 0; i + 1 < spline.Poles.size(); ++i) {
        const Sign s = SignOfFormula([&](auto f) {
            return Dot2(f, Of2(f, d), Diff2(f, spline.Poles[i + 1], spline.Poles[i]));
        });
        if (s == Sign::Zero || (common && *common != s)) {
            return Sign::Zero;
        }
        common = s;
    }
    return common.value_or(Sign::Zero);
}

/// Whether the knot vector is clamped: the spline starts at its first pole
/// and ends at its last.
[[nodiscard]] bool IsClamped(const BSplineCurve2d& spline)
{
    const auto p = static_cast<std::size_t>(spline.Degree);
    const std::vector<double>& k = spline.Knots;
    return k.size() > 2 * p + 1 && k[0] == k[p] && k[k.size() - 1] == k[k.size() - 1 - p];
}

/// The overlap of a spline along a line: proved apart when the hull of the
/// poles' parameters on the line misses its domain; the whole spline when
/// it is clamped, monotone along the line, requested over its whole domain
/// and inside the line's; NotImplemented otherwise.
[[nodiscard]] Result<PlanarIntersection> AlongLine(const BSplineCurve2d& spline, const Other& line,
                                                   const Interval& domainA, const Interval& domainB)
{
    PlanarIntersection result;
    const Interval scale =
        Interval{line.X.X} * Interval{line.X.X} + Interval{line.X.Y} * Interval{line.X.Y};
    Interval hull = Interval::Empty();
    for (const Point2& pole : spline.Poles) {
        const Interval along = (Interval{pole.X} - Interval{line.Origin.X}) * Interval{line.X.X} +
                               (Interval{pole.Y} - Interval{line.Origin.Y}) * Interval{line.X.Y};
        hull = Hull(hull, along / scale);
    }
    const Interval range = Intersection(domainA, Domain(spline));
    if (range.IsEmpty() || !Intersects(hull, domainB)) {
        return result;
    }
    const Sign sense = SenseAlong(spline, line.X);
    if (sense == Sign::Zero || !IsClamped(spline) || !Contains(domainB, hull) ||
        !Contains(range, Domain(spline))) {
        return Fail(ErrorCode::NotImplemented,
                    "BSplineCurve2d x Line2d: a B-spline along the line that is not monotone "
                    "along it, or that the requested domains cut, is not clipped yet");
    }
    result.Overlaps.push_back({range, hull, sense == Sign::Positive});
    return result;
}

} // namespace

Result<PlanarIntersection> SplineCurve2d(const Curve2d& a, const Interval& domainA,
                                         const Curve2d& b, const Interval& domainB,
                                         const Context& /*ctx*/)
{
    const BSplineCurve2d* spline = a.Target<BSplineCurve2d>();
    const std::optional<Other> other = OtherOf(b);
    if (spline == nullptr || !other) {
        return Fail(ErrorCode::InvalidInput, "SplineCurve2d: dispatched for other curve types");
    }
    const std::string handler = PairName(a, b);
    STEREON_TRY_VOID(CheckPlaneSpline(*spline, handler));
    STEREON_TRY_VOID(CheckPeriodicDomain(*spline, domainA, handler));
    if (other->Kind == OtherKind::Line && OnLine(*spline, *other)) {
        return AlongLine(*spline, *other, domainA, domainB);
    }
    STEREON_TRY(const std::vector<Found> found, SplineRoots(*spline, domainA, *other, handler));
    PlanarIntersection result;
    for (const Found& f : found) {
        STEREON_TRY_VOID(KeepRoot(result, *spline, *other, f, domainA, domainB));
    }
    Order(result);
    return result;
}

Result<PlanarIntersection> CurveSpline2d(const Curve2d& a, const Interval& domainA,
                                         const Curve2d& b, const Interval& domainB,
                                         const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, SplineCurve2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

} // namespace Stereon::Detail
