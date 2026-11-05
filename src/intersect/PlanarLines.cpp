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
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Planar.hpp"

/// Lines and circles in the plane, and a line against an ellipse, a
/// hyperbola or a parabola (the Boolean prototype's M6): the pcurves of a
/// face's loop against the pcurve of a branch. Every decision is the sign
/// of a polynomial in the input doubles: the determinant of two directions,
/// whether a parallel line passes through a point, the discriminant of the
/// quadratic a line cuts out of a conic, where its roots lie against the
/// ends of a hyperbola's rational chart, and the distance between two
/// centres against the sum and the difference of the radii. A tangency is
/// a vanishing discriminant or a distance at a threshold, decided, never
/// "nearly".

namespace Stereon::Detail {

namespace {

/// The points where a line cuts a conic, from the sign of the discriminant
/// (exact) and the roots (from coefficients to a few roundings, QuadraticOf):
/// two crossings, one tangency or none.
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
    const auto [q, discriminant] = QuadraticOf(coefficients);
    const auto keep = [&](const double t, const ContactType contact) {
        const Point2 p = Eval(line, t);
        KeepPoint(result, {t, conicParameter(p), p, contact, true}, domainLine, false, domainConic,
                  true);
    };
    if (sign == Sign::Zero) {
        keep(-q.B / q.A, ContactType::Tangential);
        return;
    }
    for (const double t : RootsOf(q, discriminant)) {
        keep(t, ContactType::Crossing);
    }
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
    AffineOverlap(result, Dot(lb->Origin - la->Origin, la->Direction) / scale,
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

// --- Line against hyperbola and parabola --------------------------------

namespace {

/// Where an end e of a hyperbola's chart lies against the two roots of a
/// quadratic f(u) = alpha u² + 2 beta u + gamma of positive discriminant.
enum class Side : std::uint8_t { Below, AtLower, Between, AtUpper, Above };

/// The side of e = ±1 from two exact signs: f(e) against alpha's (opposite
/// between the roots, alike outside them, zero on one), and alpha e + beta,
/// alpha times e's offset from the vertex -beta / alpha, which tells the
/// lower root from the upper. e is never the vertex: f is extreme there.
template <typename Quadratic>
[[nodiscard]] Side SideOf(const double e, Quadratic&& q, const Sign alpha)
{
    const Sign value = SignOfFormula([&](auto f) {
        using F = decltype(f);
        const auto [a, b, c] = q(f);
        return F::Add(F::Add(a, F::Mul(F::Of(2.0 * e), b)), c); // e² = 1
    });
    const Sign offset = SignOfFormula([&](auto f) {
        using F = decltype(f);
        const auto [a, b, c] = q(f);
        return F::Add(F::Mul(F::Of(e), a), b);
    });
    const bool right = offset == alpha;
    if (value == Sign::Zero) {
        return right ? Side::AtUpper : Side::AtLower;
    }
    if (value != alpha) {
        return Side::Between;
    }
    return right ? Side::Above : Side::Below;
}

/// Which of the two roots, lower and upper, lie strictly inside (-1, 1).
template <typename Quadratic>
[[nodiscard]] std::array<bool, 2> InsideChart(Quadratic&& q, const Sign alpha)
{
    const Side low = SideOf(-1.0, q, alpha);
    const Side high = SideOf(1.0, q, alpha);
    const bool lower = low == Side::Below && high != Side::Below && high != Side::AtLower;
    const bool upper = high == Side::Above && low != Side::Above && low != Side::AtUpper;
    return {lower, upper};
}

/// The single root of 2 beta u + gamma = 0 when alpha is exactly zero, on a
/// hyperbola's chart when |gamma| < 2 |beta|, decided exactly.
template <typename Quadratic, typename Place>
[[nodiscard]] Result<void> LinearRoot(Quadratic&& q, const bool chart, Place&& place)
{
    const Sign beta = SignOfFormula([&](auto f) { return q(f)[1]; });
    if (beta == Sign::Zero) {
        return {}; // A non-zero constant: a line along an asymptote, or none.
    }
    if (chart && SignOfFormula([&](auto f) {
                     using F = decltype(f);
                     const auto [a, b, c] = q(f);
                     return F::Sub(F::Mul(F::Of(4.0), F::Mul(b, b)), F::Mul(c, c));
                 }) != Sign::Positive) {
        return {};
    }
    const double b = ValueOfFormula([&](auto f) { return q(f)[1]; });
    const double c = ValueOfFormula([&](auto f) { return q(f)[2]; });
    return place(-c / (2.0 * b), ContactType::Crossing);
}

/// A line's meetings with a conic whose points C(u) it cuts as the roots of
/// f(u) = alpha u² + 2 beta u + gamma: by the exact sign of the
/// discriminant, two crossings, one tangency or none; one crossing when
/// alpha is exactly zero. On a hyperbola's chart (`chart`), only the roots
/// strictly inside (-1, 1) are points of the curve.
template <typename Quadratic, typename Place>
[[nodiscard]] Result<void> CutByLine(Quadratic&& q, const bool chart, Place&& place)
{
    const Sign alpha = SignOfFormula([&](auto f) { return q(f)[0]; });
    if (alpha == Sign::Zero) {
        return LinearRoot(q, chart, place);
    }
    const Sign sign = SignOfFormula([&](auto f) {
        using F = decltype(f);
        const auto [a, b, c] = q(f);
        return F::Sub(F::Mul(b, b), F::Mul(a, c));
    });
    if (sign == Sign::Negative) {
        return {};
    }
    const auto [quadratic, discriminant] = QuadraticOf(q);
    if (sign == Sign::Zero) {
        const Sign vertex = SignOfFormula([&](auto f) {
            using F = decltype(f);
            const auto [qa, qb, qc] = q(f);
            return F::Sub(F::Mul(qa, qa), F::Mul(qb, qb));
        });
        return !chart || vertex == Sign::Positive
                   ? place(-quadratic.B / quadratic.A, ContactType::Tangential)
                   : Result<void>{};
    }
    const std::array<bool, 2> inside =
        chart ? InsideChart(q, alpha) : std::array<bool, 2>{true, true};
    const std::array<double, 2> roots = RootsOf(quadratic, discriminant);
    for (std::size_t i = 0; i < roots.size(); ++i) {
        if (inside[i]) {
            STEREON_TRY_VOID(place(roots[i], ContactType::Crossing));
        }
    }
    return {};
}

/// The line's normal n and its offset k0 = n . (c - O) from a point c, and
/// n's components along a frame: the line's equation n . (P - O) = 0 in
/// the conic's frame.
template <typename F>
struct LineInFrame {
    typename F::Scalar K0;
    typename F::Scalar NX;
    typename F::Scalar NY;
};

template <typename F>
[[nodiscard]] LineInFrame<F> InFrame(F f, const Line2d& line, const Point2& c, const Vec2& x,
                                     const Vec2& y)
{
    const Pair2<F> n = Of2(f, Perpendicular(line.Direction));
    return {Dot2(f, n, Diff2(f, c, line.Origin)), Dot2(f, n, Of2(f, x)), Dot2(f, n, Of2(f, y))};
}

/// The line's parameter of a point on it.
[[nodiscard]] double OnLine(const Line2d& line, const Point2& p)
{
    return Dot(p - line.Origin, line.Direction) / Dot(line.Direction, line.Direction);
}

} // namespace

Result<PlanarIntersection> LineHyperbola2d(const Curve2d& a, const Interval& domainA,
                                           const Curve2d& b, const Interval& domainB,
                                           const Context& /*ctx*/)
{
    const Line2d* line = a.Target<Line2d>();
    const Hyperbola2d* hyperbola = b.Target<Hyperbola2d>();
    if (line == nullptr || hyperbola == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineHyperbola2d: dispatched for other curve types");
    }
    // The branch's rational chart, t = 2 artanh u for u in (-1, 1):
    // (1 - u²) C(t) = (1 - u²) c + a (1 + u²) X + 2 b u Y, so the line's
    // equation times (1 - u²) is (kx - k0) u² + 2 ky u + (k0 + kx) = 0 with
    // k0 = n . (c - O), kx = a n . X and ky = b n . Y.
    const auto quadratic = [&](auto f) {
        using F = decltype(f);
        const LineInFrame<F> k =
            InFrame(f, *line, hyperbola->Center, hyperbola->XAxis, hyperbola->YAxis);
        const auto kx = F::Mul(F::Of(hyperbola->MajorRadius), k.NX);
        const auto ky = F::Mul(F::Of(hyperbola->MinorRadius), k.NY);
        return std::array{F::Sub(kx, k.K0), ky, F::Add(k.K0, kx)};
    };
    PlanarIntersection result;
    const auto place = [&](const double u, const ContactType contact) -> Result<void> {
        const double t = 2.0 * std::atanh(u);
        if (!std::isfinite(t)) {
            return Fail(ErrorCode::IntersectionFailed,
                        "Line2d x Hyperbola2d: a meeting beyond the range of the hyperbola's "
                        "parameter in doubles");
        }
        const Point2 p = Eval(*hyperbola, t);
        KeepPoint(result, {OnLine(*line, p), t, p, contact, true}, domainA, false, domainB, false);
        return {};
    };
    STEREON_TRY_VOID(CutByLine(quadratic, true, place));
    Order(result);
    return result;
}

Result<PlanarIntersection> HyperbolaLine2d(const Curve2d& a, const Interval& domainA,
                                           const Curve2d& b, const Interval& domainB,
                                           const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, LineHyperbola2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

Result<PlanarIntersection> LineParabola2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& /*ctx*/)
{
    const Line2d* line = a.Target<Line2d>();
    const Parabola2d* parabola = b.Target<Parabola2d>();
    if (line == nullptr || parabola == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineParabola2d: dispatched for other curve types");
    }
    // 4 f C(t) = 4 f V + t² X + 4 f t Y, so the line's equation times 4 f
    // is (n . X) t² + 2 (2 f n . Y) t + 4 f n . (V - O) = 0.
    const auto quadratic = [&](auto f) {
        using F = decltype(f);
        const LineInFrame<F> k =
            InFrame(f, *line, parabola->Vertex, parabola->XAxis, parabola->YAxis);
        const auto focal = F::Of(parabola->FocalLength);
        return std::array{k.NX, F::Mul(F::Mul(F::Of(2.0), focal), k.NY),
                          F::Mul(F::Mul(F::Of(4.0), focal), k.K0)};
    };
    PlanarIntersection result;
    const auto place = [&](const double t, const ContactType contact) -> Result<void> {
        const Point2 p = Eval(*parabola, t);
        KeepPoint(result, {OnLine(*line, p), t, p, contact, true}, domainA, false, domainB, false);
        return {};
    };
    STEREON_TRY_VOID(CutByLine(quadratic, false, place));
    Order(result);
    return result;
}

Result<PlanarIntersection> ParabolaLine2d(const Curve2d& a, const Interval& domainA,
                                          const Curve2d& b, const Interval& domainB,
                                          const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, LineParabola2d(b, domainB, a, domainA, ctx));
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
/// angle of A's x axis in B's frame.
void CoincidentCircles(PlanarIntersection& result, const Circle2d& a, const Circle2d& b,
                       const Interval& domainA, const Interval& domainB)
{
    const double sense = Cross(a.XAxis, a.YAxis) * Cross(b.XAxis, b.YAxis) > 0.0 ? 1.0 : -1.0;
    const double offset = Atan2(Dot(a.XAxis, b.YAxis), Dot(a.XAxis, b.XAxis));
    PeriodicOverlaps(result, sense, offset, domainA, domainB);
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
