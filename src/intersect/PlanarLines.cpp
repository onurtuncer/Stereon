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
#include "Planar.hpp"

/// Lines and circles in the plane, and a line against an ellipse (the
/// Boolean prototype's M6): the pcurves of a face's loop against the
/// pcurve of a branch. Every decision is the sign of a polynomial in the
/// input doubles: the determinant of two directions, whether a parallel
/// line passes through a point, the discriminant of the quadratic a line
/// cuts out of a circle or an ellipse, and the distance between two
/// centres against the sum and the difference of the radii. A tangency is
/// a vanishing discriminant or a distance at a threshold, decided, never
/// "nearly".

namespace Stereon::Detail {

namespace {

/// The quadratic A t² + 2 B t + C of a line against a circle or an ellipse,
/// its coefficients as a formula over a field, and its roots by the form
/// that does not cancel.
struct Quadratic {
    double A;
    double B;
    double C;
};

/// The two roots, the smaller first; `discriminant` is B² - A C >= 0.
[[nodiscard]] std::array<double, 2> RootsOf(const Quadratic& q, const double discriminant)
{
    const double root = std::sqrt(discriminant);
    const double head = -(q.B + (q.B < 0.0 ? -root : root));
    const double t1 = head / q.A;
    const double t2 = head == 0.0 ? t1 : q.C / head;
    return {std::min(t1, t2), std::max(t1, t2)};
}

/// The points where a line cuts a conic, from the sign of the discriminant
/// (exact) and the roots (rounded): two crossings, one tangency or none.
template <typename Coefficients>
void CutConic(PlanarIntersection& result, const Line2d& line, const Interval& domainLine,
              const Interval& domainConic, Coefficients&& coefficients, const auto& conicParameter)
{
    const Sign sign = SignOfFormula([&](auto f) {
        using F = decltype(f);
        const auto [a, b, c] = coefficients(f);
        return F::Sub(F::Mul(b, b), F::Mul(a, c));
    });
    if (sign == Sign::Negative) {
        return;
    }
    const auto [a, b, c] = coefficients(Plain{});
    const Quadratic q{a, b, c};
    const auto keep = [&](const double t, const ContactType contact) {
        const Point2 p = Eval(line, t);
        KeepPoint(result, {t, conicParameter(p), p, contact, true}, domainLine, false, domainConic,
                  true);
    };
    if (sign == Sign::Zero) {
        keep(-q.B / q.A, ContactType::Tangential);
        return;
    }
    const double discriminant = std::max(0.0, q.B * q.B - q.A * q.C);
    for (const double t : RootsOf(q, discriminant)) {
        keep(t, ContactType::Crossing);
    }
}

// --- Line against line --------------------------------------------------

/// The stretch where coincident lines share their domains: B's parameter s
/// maps to A's t = t0 + k s.
void CoincidentLines(PlanarIntersection& result, const double t0, const double k,
                     const Interval& domainA, const Interval& domainB)
{
    const Interval onB = Intersection(domainB, AxisRange(t0, k, domainA));
    if (onB.IsEmpty()) {
        return;
    }
    const Interval onA = Intersection(domainA, Interval{t0} + Interval{k} * onB);
    result.Overlaps.push_back({onA, onB, k > 0.0});
}

} // namespace

Result<PlanarIntersection> LineLine2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                      const Interval& domainB, const Context& /*ctx*/)
{
    const Line2d* la = a.Target<Line2d>();
    const Line2d* lb = b.Target<Line2d>();
    if (la == nullptr || lb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineLine2d: dispatched for other curve types");
    }
    PlanarIntersection result;
    const Sign det = SignOfFormula(
        [&](auto f) { return Cross2(f, Of2(f, la->Direction), Of2(f, lb->Direction)); });
    if (det != Sign::Zero) {
        const Vec2 w = lb->Origin - la->Origin;
        const double determinant = Cross(la->Direction, lb->Direction);
        const double t = Cross(w, lb->Direction) / determinant;
        const double s = Cross(w, la->Direction) / determinant;
        KeepPoint(result, {t, s, Eval(*la, t), ContactType::Crossing, true}, domainA, false,
                  domainB, false);
        return result;
    }
    // Parallel: coincident when B's origin lies on A, decided exactly.
    const Sign apart = SignOfFormula(
        [&](auto f) { return Cross2(f, Diff2(f, lb->Origin, la->Origin), Of2(f, la->Direction)); });
    if (apart != Sign::Zero) {
        return result;
    }
    const double scale = Dot(la->Direction, la->Direction);
    CoincidentLines(result, Dot(lb->Origin - la->Origin, la->Direction) / scale,
                    Dot(lb->Direction, la->Direction) / scale, domainA, domainB);
    return result;
}

// --- Line against circle ------------------------------------------------

Result<PlanarIntersection> LineCircle2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                        const Interval& domainB, const Context& /*ctx*/)
{
    const Line2d* line = a.Target<Line2d>();
    const Circle2d* circle = b.Target<Circle2d>();
    if (line == nullptr || circle == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineCircle2d: dispatched for other curve types");
    }
    PlanarIntersection result;
    // |O + t d - C|² = r²: A = d.d, B = d.(O - C), C = |O - C|² - r².
    const auto coefficients = [&](auto f) {
        using F = decltype(f);
        const Pair2<F> d = Of2(f, line->Direction);
        const Pair2<F> w = Diff2(f, line->Origin, circle->Center);
        const auto r = F::Of(circle->Radius);
        return std::array{Dot2(f, d, d), Dot2(f, d, w), F::Sub(Dot2(f, w, w), F::Mul(r, r))};
    };
    CutConic(result, *line, domainA, domainB, coefficients, [&](const Point2& p) {
        return AngleIn(circle->Center, circle->XAxis, circle->YAxis, circle->Radius, circle->Radius,
                       p);
    });
    Order(result);
    return result;
}

Result<PlanarIntersection> CircleLine2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                        const Interval& domainB, const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, LineCircle2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

// --- Line against ellipse -----------------------------------------------

Result<PlanarIntersection> LineEllipse2d(const Curve2d& a, const Interval& domainA,
                                         const Curve2d& b, const Interval& domainB,
                                         const Context& /*ctx*/)
{
    const Line2d* line = a.Target<Line2d>();
    const Ellipse2d* ellipse = b.Target<Ellipse2d>();
    if (line == nullptr || ellipse == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineEllipse2d: dispatched for other curve types");
    }
    PlanarIntersection result;
    // In the ellipse's frame the line is (x0 + t dx, y0 + t dy), and the
    // ellipse b² x² + a² y² = a² b².
    const auto coefficients = [&](auto f) {
        using F = decltype(f);
        const Pair2<F> d = Of2(f, line->Direction);
        const Pair2<F> w = Diff2(f, line->Origin, ellipse->Center);
        const Pair2<F> x = Of2(f, ellipse->XAxis);
        const Pair2<F> y = Of2(f, ellipse->YAxis);
        const auto aa = F::Mul(F::Of(ellipse->MajorRadius), F::Of(ellipse->MajorRadius));
        const auto bb = F::Mul(F::Of(ellipse->MinorRadius), F::Of(ellipse->MinorRadius));
        const auto x0 = Dot2(f, w, x);
        const auto y0 = Dot2(f, w, y);
        const auto dx = Dot2(f, d, x);
        const auto dy = Dot2(f, d, y);
        return std::array{
            F::Add(F::Mul(bb, F::Mul(dx, dx)), F::Mul(aa, F::Mul(dy, dy))),
            F::Add(F::Mul(bb, F::Mul(x0, dx)), F::Mul(aa, F::Mul(y0, dy))),
            F::Sub(F::Add(F::Mul(bb, F::Mul(x0, x0)), F::Mul(aa, F::Mul(y0, y0))), F::Mul(aa, bb))};
    };
    CutConic(result, *line, domainA, domainB, coefficients, [&](const Point2& p) {
        return AngleIn(ellipse->Center, ellipse->XAxis, ellipse->YAxis, ellipse->MajorRadius,
                       ellipse->MinorRadius, p);
    });
    Order(result);
    return result;
}

Result<PlanarIntersection> EllipseLine2d(const Curve2d& a, const Interval& domainA,
                                         const Curve2d& b, const Interval& domainB,
                                         const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, LineEllipse2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

// --- Circle against circle ----------------------------------------------

namespace {

/// The sign of |C2 - C1|² - k² for k the sum or the difference of the radii.
[[nodiscard]] Sign CompareCentres(const Circle2d& a, const Circle2d& b, const bool sum)
{
    return SignOfFormula([&](auto f) {
        using F = decltype(f);
        const Pair2<F> w = Diff2(f, b.Center, a.Center);
        const auto k = sum ? F::Add(F::Of(a.Radius), F::Of(b.Radius))
                           : F::Sub(F::Of(a.Radius), F::Of(b.Radius));
        return F::Sub(Dot2(f, w, w), F::Mul(k, k));
    });
}

/// Coincident circles: B's angle is sense * t + offset for A's angle t,
/// the sense that of the frames' turns against each other, the offset the
/// angle of A's x axis in B's frame. The overlaps are the requested
/// domains' common stretches modulo the period, one per period shift.
void CoincidentCircles(PlanarIntersection& result, const Circle2d& a, const Circle2d& b,
                       const Interval& domainA, const Interval& domainB)
{
    const double sense = Cross(a.XAxis, a.YAxis) * Cross(b.XAxis, b.YAxis) > 0.0 ? 1.0 : -1.0;
    const double offset = Atan2(Dot(a.XAxis, b.YAxis), Dot(a.XAxis, b.XAxis));
    // One turn of each domain describes it; an unbounded one is a turn from 0.
    const auto turnOf = [](const Interval& domain) {
        const double start = std::isfinite(domain.GetLower()) ? domain.GetLower() : 0.0;
        return Intersection(domain, Interval{start, start + CirclePeriod});
    };
    const Interval arcA = turnOf(domainA);
    const Interval arcB = turnOf(domainB);
    // s = sense t + offset + 2 pi k: the shifts whose image can meet arcB.
    const Interval image = Interval{sense} * arcA + Interval{offset};
    const double first = std::floor((arcB.GetLower() - image.GetUpper()) / CirclePeriod);
    const double last = std::ceil((arcB.GetUpper() - image.GetLower()) / CirclePeriod);
    for (double k = first; k <= last; k += 1.0) {
        const Interval onB = Intersection(arcB, image + Interval{k} * Interval{CirclePeriod});
        if (onB.IsEmpty()) {
            continue;
        }
        const Interval onA =
            Intersection(arcA, Interval{sense} *
                                   (onB - Interval{offset} - Interval{k} * Interval{CirclePeriod}));
        if (!onA.IsEmpty()) {
            result.Overlaps.push_back({onA, onB, sense > 0.0});
        }
    }
}

/// The points of two circles that are not coincident: one tangency at the
/// sum or the difference of the radii, two crossings between, none beyond.
void MeetingCircles(PlanarIntersection& result, const Circle2d& a, const Circle2d& b,
                    const Sign outer, const Sign inner, const Interval& domainA,
                    const Interval& domainB)
{
    const Vec2 w = b.Center - a.Center;
    const double distance = Length(w);
    const Vec2 along = w / distance;
    const auto keep = [&](const Point2& p, const ContactType contact) {
        KeepPoint(result,
                  {AngleIn(a.Center, a.XAxis, a.YAxis, a.Radius, a.Radius, p),
                   AngleIn(b.Center, b.XAxis, b.YAxis, b.Radius, b.Radius, p), p, contact, true},
                  domainA, true, domainB, true);
    };
    if (outer == Sign::Zero) {
        keep(a.Center + a.Radius * along, ContactType::Tangential);
        return;
    }
    if (inner == Sign::Zero) {
        // B's centre is inside A when A is the larger: the touch is on its side.
        keep(a.Center + (a.Radius > b.Radius ? a.Radius : -a.Radius) * along,
             ContactType::Tangential);
        return;
    }
    const double x =
        (distance * distance + a.Radius * a.Radius - b.Radius * b.Radius) / (2.0 * distance);
    const double h = std::sqrt(std::max(0.0, a.Radius * a.Radius - x * x));
    const Vec2 across = Perpendicular(along);
    keep(a.Center + x * along + h * across, ContactType::Crossing);
    keep(a.Center + x * along - h * across, ContactType::Crossing);
}

} // namespace

Result<PlanarIntersection> CircleCircle2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& /*ctx*/)
{
    const Circle2d* ca = a.Target<Circle2d>();
    const Circle2d* cb = b.Target<Circle2d>();
    if (ca == nullptr || cb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "CircleCircle2d: dispatched for other curve types");
    }
    PlanarIntersection result;
    const Sign outer = CompareCentres(*ca, *cb, true);
    if (outer == Sign::Positive) {
        return result; // Apart.
    }
    const Sign inner = CompareCentres(*ca, *cb, false);
    if (inner == Sign::Negative) {
        return result; // Nested.
    }
    const Sign apart = SignOfFormula([&](auto f) {
        const auto w = Diff2(f, cb->Center, ca->Center);
        return Dot2(f, w, w);
    });
    if (apart == Sign::Zero && ca->Radius == cb->Radius) {
        CoincidentCircles(result, *ca, *cb, domainA, domainB);
        return result;
    }
    if (apart == Sign::Zero) {
        return result; // Concentric, different radii.
    }
    MeetingCircles(result, *ca, *cb, outer, inner, domainA, domainB);
    Order(result);
    return result;
}

} // namespace Stereon::Detail
