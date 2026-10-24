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
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Planar.hpp"

/// A non-rational B-spline in the plane against a line or a circle (the
/// Boolean prototype's M6): the pcurve of an ellipse on a cylinder, or of
/// the skew cylinders' quartic, against a loop's edges. Each Bezier
/// segment of the spline is a polynomial in its local parameter, so the
/// signed distance to the line, or the squared distance from the circle's
/// centre less its radius squared, is a polynomial too, whose roots the
/// Bernstein subdivision of PolynomialRoots isolates. Nothing here is
/// exact: the roots are bracketed to rounding, a root where the curve's
/// tangent runs along the other curve is a touch, and every point says
/// so with Decided = false. A root of even multiplicity may be missed, as
/// PolynomialRoots documents, and so may a root at the very end of the
/// spline's domain; a touch that is missed changes no side, and a face's
/// vertices are met by their tolerances, not by roots (ADR-0001).

namespace Stereon::Detail {

namespace {

/// The binomial coefficient C(n, k).
[[nodiscard]] double Binomial(const std::size_t n, const std::size_t k)
{
    double value = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        value = value * static_cast<double>(n + 1 - i) / static_cast<double>(i);
    }
    return value;
}

/// The power-basis coefficients on [0, 1] of the Bernstein polynomial with
/// the given coefficients: c_j = C(p, j) sum_i (-1)^(j - i) C(j, i) b_i.
[[nodiscard]] std::vector<double> PowerBasis(const std::span<const double> bernstein)
{
    const std::size_t p = bernstein.size() - 1;
    std::vector<double> power(bernstein.size(), 0.0);
    for (std::size_t j = 0; j <= p; ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i <= j; ++i) {
            const double sign = (j - i) % 2 == 0 ? 1.0 : -1.0;
            sum += sign * Binomial(j, i) * bernstein[i];
        }
        power[j] = Binomial(p, j) * sum;
    }
    return power;
}

[[nodiscard]] std::vector<double> Multiply(const std::span<const double> p,
                                           const std::span<const double> q)
{
    std::vector<double> r(p.size() + q.size() - 1, 0.0);
    for (std::size_t i = 0; i < p.size(); ++i) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            r[i + j] += p[i] * q[j];
        }
    }
    return r;
}

/// A Bezier segment of the spline as power-basis polynomials of its local
/// parameter in [0, 1], and the stretch of the spline's parameter it covers.
struct Segment {
    std::vector<double> X;
    std::vector<double> Y;
    Interval Range;
};

[[nodiscard]] Segment SegmentOf(const BSplineCurve2d& bezier)
{
    std::vector<double> x;
    std::vector<double> y;
    for (const Point2& pole : bezier.Poles) {
        x.push_back(pole.X);
        y.push_back(pole.Y);
    }
    return {PowerBasis(x), PowerBasis(y), Domain(bezier)};
}

/// The roots of the polynomial over the segment, as parameters of the
/// spline inside the requested domain, each with the point and its contact:
/// a crossing unless the spline's tangent there is along `grazing`'s
/// direction, decided in doubles.
template <typename Grazing, typename ParameterB>
void RootsOfSegment(PlanarIntersection& result, const BSplineCurve2d& spline,
                    const Segment& segment, const std::span<const double> polynomial,
                    const Interval& domainA, const Interval& domainB, const bool periodicB,
                    Grazing&& grazing, ParameterB&& parameterB)
{
    const double width = Width(segment.Range);
    for (const double local : PolynomialRoots(polynomial, Interval{0.0, 1.0})) {
        const double t = segment.Range.GetLower() + width * local;
        const Point2 p = Eval(spline, t);
        const bool touch = grazing(p, Deriv(spline, t, 1));
        KeepPoint(
            result,
            {t, parameterB(p), p, touch ? ContactType::Tangential : ContactType::Crossing, false},
            domainA, false, domainB, periodicB);
    }
}

/// Roots at a segment's end are found from both neighbours: one point.
void Dedupe(PlanarIntersection& result)
{
    Order(result);
    std::vector<PlanarPoint> distinct;
    for (const PlanarPoint& p : result.Points) {
        if (!distinct.empty() && std::abs(p.ParameterA - distinct.back().ParameterA) <=
                                     1e-12 * std::max(1.0, std::abs(p.ParameterA))) {
            continue; // STEREON-LINT-ALLOW(epsilon): a parameter gap of rounding
        }
        distinct.push_back(p);
    }
    result.Points = std::move(distinct);
}

/// The spline's segments that reach the requested domain.
[[nodiscard]] Result<std::vector<Segment>> SegmentsOf(const BSplineCurve2d& spline,
                                                      const Interval& domain)
{
    if (IsRational(spline)) {
        return Fail(ErrorCode::NotSupported,
                    "plane curve intersection covers non-rational B-splines only");
    }
    std::vector<Segment> segments;
    for (const BSplineCurve2d& bezier : ToBezierSegments(spline)) {
        Segment segment = SegmentOf(bezier);
        if (Intersects(segment.Range, domain)) {
            segments.push_back(std::move(segment));
        }
    }
    return segments;
}

} // namespace

Result<PlanarIntersection> SplineLine2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                        const Interval& domainB, const Context& /*ctx*/)
{
    const BSplineCurve2d* spline = a.Target<BSplineCurve2d>();
    const Line2d* line = b.Target<Line2d>();
    if (spline == nullptr || line == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineLine2d: dispatched for other curve types");
    }
    STEREON_TRY(const std::vector<Segment> segments, SegmentsOf(*spline, domainA));
    PlanarIntersection result;
    // The signed distance n . (B(t) - O), n across the line.
    const Vec2 normal = Perpendicular(line->Direction);
    const double scale = Dot(line->Direction, line->Direction);
    for (const Segment& segment : segments) {
        // The constant term carries the origin (no separate [0]: GCC's
        // null-dereference check cannot see the size).
        const double offset = Dot(normal, line->Origin - Point2{});
        std::vector<double> distance(segment.X.size());
        for (std::size_t i = 0; i < distance.size(); ++i) {
            distance[i] =
                normal.X * segment.X[i] + normal.Y * segment.Y[i] - (i == 0 ? offset : 0.0);
        }
        RootsOfSegment(
            result, *spline, segment, distance, domainA, domainB, false,
            [&](const Point2&, const Vec2& tangent) { return Dot(normal, tangent) == 0.0; },
            [&](const Point2& p) { return Dot(p - line->Origin, line->Direction) / scale; });
    }
    Dedupe(result);
    return result;
}

Result<PlanarIntersection> LineSpline2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                        const Interval& domainB, const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, SplineLine2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

Result<PlanarIntersection> SplineCircle2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& /*ctx*/)
{
    const BSplineCurve2d* spline = a.Target<BSplineCurve2d>();
    const Circle2d* circle = b.Target<Circle2d>();
    if (spline == nullptr || circle == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineCircle2d: dispatched for other curve types");
    }
    STEREON_TRY(const std::vector<Segment> segments, SegmentsOf(*spline, domainA));
    PlanarIntersection result;
    // |B(t) - C|² - r² = (x - Cx)² + (y - Cy)² - r².
    for (const Segment& segment : segments) {
        // The constant terms carry the centre and the radius (no separate
        // [0]: GCC's null-dereference check cannot see the size).
        std::vector<double> x(segment.X.size());
        std::vector<double> y(segment.Y.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = segment.X[i] - (i == 0 ? circle->Center.X : 0.0);
            y[i] = segment.Y[i] - (i == 0 ? circle->Center.Y : 0.0);
        }
        const std::vector<double> xx = Multiply(x, x);
        const std::vector<double> yy = Multiply(y, y);
        std::vector<double> squared(xx.size());
        for (std::size_t i = 0; i < squared.size(); ++i) {
            squared[i] = xx[i] + yy[i] - (i == 0 ? circle->Radius * circle->Radius : 0.0);
        }
        RootsOfSegment(
            result, *spline, segment, squared, domainA, domainB, true,
            [&](const Point2& p, const Vec2& tangent) {
                return Dot(p - circle->Center, tangent) == 0.0;
            },
            [&](const Point2& p) {
                return AngleIn(circle->Center, circle->XAxis, circle->YAxis, circle->Radius,
                               circle->Radius, p);
            });
    }
    Dedupe(result);
    return result;
}

Result<PlanarIntersection> CircleSpline2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, SplineCircle2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

} // namespace Stereon::Detail
