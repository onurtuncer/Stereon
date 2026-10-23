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
#include "Planar.hpp"

/// Two conics in the plane, a circle or an ellipse against an ellipse (the
/// Boolean prototype's M6): a quartic. The point of A at angle s, written
/// in B's frame, satisfies B's equation b² x² + a² y² = a² b², a
/// trigonometric polynomial of degree two in s. With u = tan(s/2) over
/// each half of the turn it is a quartic in u whose coefficients are
/// polynomials in the input doubles. Its simple roots are isolated by
/// Bernstein subdivision and each is certified a crossing by the exact
/// signs of the quartic on either side of it; a root where the exact
/// signs agree is a touch or a near miss that rounding cannot tell apart,
/// reported as a tangency bracketed to rounding, not decided. Coincident
/// conics are not handled: two ellipses given as one curve twice are
/// outside what a face's loops produce.

namespace Stereon::Detail {

namespace {

/// A conic as an ellipse: a circle has equal semi-axes.
struct Conic {
    Point2 Center;
    double A;
    double B;
    Vec2 X;
    Vec2 Y;
};

[[nodiscard]] Conic ConicOf(const Circle2d& c)
{
    return {c.Center, c.Radius, c.Radius, c.XAxis, c.YAxis};
}

[[nodiscard]] Conic ConicOf(const Ellipse2d& e)
{
    return {e.Center, e.MajorRadius, e.MinorRadius, e.XAxis, e.YAxis};
}

/// A polynomial over a field, coefficients in increasing degree.
template <typename F>
using Polynomial = std::vector<typename F::Scalar>;

template <typename F>
[[nodiscard]] Polynomial<F> Multiply(F, const Polynomial<F>& p, const Polynomial<F>& q)
{
    Polynomial<F> r(p.size() + q.size() - 1, F::Of(0.0));
    for (std::size_t i = 0; i < p.size(); ++i) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            r[i + j] = F::Add(r[i + j], F::Mul(p[i], q[j]));
        }
    }
    return r;
}

template <typename F>
[[nodiscard]] Polynomial<F> Scaled(F, const Polynomial<F>& p, const typename F::Scalar& k)
{
    Polynomial<F> r;
    r.reserve(p.size());
    for (const auto& c : p) {
        r.push_back(F::Mul(c, k));
    }
    return r;
}

template <typename F>
[[nodiscard]] Polynomial<F> Combine(F, const Polynomial<F>& p, const Polynomial<F>& q,
                                    const bool subtract)
{
    Polynomial<F> r(std::max(p.size(), q.size()), F::Of(0.0));
    for (std::size_t i = 0; i < r.size(); ++i) {
        const auto pi = i < p.size() ? p[i] : F::Of(0.0);
        const auto qi = i < q.size() ? q[i] : F::Of(0.0);
        r[i] = subtract ? F::Sub(pi, qi) : F::Add(pi, qi);
    }
    return r;
}

/// The quartic in u = tan(s/2) of A's point at angle s (or s + pi when
/// `otherHalf`) against B's equation, times (1 + u²)². A coordinate of A's
/// point in B's frame is p0 + p1 cos s + p2 sin s, which (1 + u²) turns into
/// the quadratic (p0 - p1) u² + 2 p2 u + (p0 + p1).
template <typename F>
[[nodiscard]] Polynomial<F> Quartic(F, const Conic& a, const Conic& b, const bool otherHalf)
{
    const auto half = otherHalf ? -1.0 : 1.0;
    const Pair2<F> w = Diff2(F{}, a.Center, b.Center);
    const auto coordinate = [&](const Vec2& axis) {
        const Pair2<F> e = Of2(F{}, axis);
        const auto p0 = Dot2(F{}, w, e);
        const auto p1 = F::Mul(F::Of(half * a.A), Dot2(F{}, Of2(F{}, a.X), e));
        const auto p2 = F::Mul(F::Of(half * a.B), Dot2(F{}, Of2(F{}, a.Y), e));
        return Polynomial<F>{F::Add(p0, p1), F::Mul(F::Of(2.0), p2), F::Sub(p0, p1)};
    };
    const Polynomial<F> x = coordinate(b.X);
    const Polynomial<F> y = coordinate(b.Y);
    const auto aa = F::Mul(F::Of(b.A), F::Of(b.A));
    const auto bb = F::Mul(F::Of(b.B), F::Of(b.B));
    const Polynomial<F> one{F::Of(1.0), F::Of(0.0), F::Of(1.0)};
    const Polynomial<F> lhs = Combine(F{}, Scaled(F{}, Multiply(F{}, x, x), bb),
                                      Scaled(F{}, Multiply(F{}, y, y), aa), false);
    return Combine(F{}, lhs, Scaled(F{}, Multiply(F{}, one, one), F::Mul(aa, bb)), true);
}

/// The sign of the polynomial at u, exactly.
[[nodiscard]] Sign ExactSignAt(const Polynomial<Exact>& p, const double u)
{
    DynamicExpansion value;
    const DynamicExpansion at{u};
    for (std::size_t i = p.size(); i-- > 0;) {
        value = Sum(Product(value, at), p[i]);
    }
    return value.GetSign();
}

/// Each isolated root u of the quartic, as a point: a crossing when the
/// exact signs on either side of it differ, otherwise a touch bracketed to
/// rounding. s = 2 atan u, plus pi for the other half.
void PointsOfHalf(PlanarIntersection& result, const Conic& a, const Conic& b, const bool otherHalf,
                  const Interval& domainA, const Interval& domainB)
{
    const Polynomial<Plain> rounded = Quartic(Plain{}, a, b, otherHalf);
    const Polynomial<Exact> exact = Quartic(Exact{}, a, b, otherHalf);
    for (const double u : PolynomialRoots(rounded, Interval{-1.0, 1.0})) {
        const double below = std::nextafter(std::nextafter(u, -Infinity), -Infinity);
        const double above = std::nextafter(std::nextafter(u, Infinity), Infinity);
        const Sign before = ExactSignAt(exact, below);
        const Sign after = ExactSignAt(exact, above);
        const bool crossing = before != Sign::Zero && after != Sign::Zero && before != after;
        const double s = 2.0 * std::atan(u) + (otherHalf ? std::numbers::pi : 0.0);
        const Point2 p = a.Center + a.A * Cos(s) * a.X + a.B * Sin(s) * a.Y;
        KeepPoint(result,
                  {s, AngleIn(b.Center, b.X, b.Y, b.A, b.B, p), p,
                   crossing ? ContactType::Crossing : ContactType::Tangential, crossing},
                  domainA, true, domainB, true);
    }
}

/// Points within rounding of each other, from the two halves meeting at
/// s = ±pi/2, are one point.
void Dedupe(PlanarIntersection& result)
{
    Order(result);
    std::vector<PlanarPoint> distinct;
    for (const PlanarPoint& p : result.Points) {
        if (!distinct.empty() &&
            std::abs(p.ParameterA - distinct.back().ParameterA) <= 1e-12 * CirclePeriod) {
            continue; // STEREON-LINT-ALLOW(epsilon): a parameter gap of rounding
        }
        distinct.push_back(p);
    }
    result.Points = std::move(distinct);
}

[[nodiscard]] Result<PlanarIntersection> ConicConic(const Conic& a, const Interval& domainA,
                                                    const Conic& b, const Interval& domainB)
{
    PlanarIntersection result;
    PointsOfHalf(result, a, b, false, domainA, domainB);
    PointsOfHalf(result, a, b, true, domainA, domainB);
    Dedupe(result);
    return result;
}

} // namespace

Result<PlanarIntersection> CircleEllipse2d(const Curve2d& a, const Interval& domainA,
                                           const Curve2d& b, const Interval& domainB,
                                           const Context& /*ctx*/)
{
    const Circle2d* circle = a.Target<Circle2d>();
    const Ellipse2d* ellipse = b.Target<Ellipse2d>();
    if (circle == nullptr || ellipse == nullptr) {
        return Fail(ErrorCode::InvalidInput, "CircleEllipse2d: dispatched for other curve types");
    }
    return ConicConic(ConicOf(*circle), domainA, ConicOf(*ellipse), domainB);
}

Result<PlanarIntersection> EllipseCircle2d(const Curve2d& a, const Interval& domainA,
                                           const Curve2d& b, const Interval& domainB,
                                           const Context& ctx)
{
    STEREON_TRY(PlanarIntersection result, CircleEllipse2d(b, domainB, a, domainA, ctx));
    return SwappedPlanar(std::move(result));
}

Result<PlanarIntersection> EllipseEllipse2d(const Curve2d& a, const Interval& domainA,
                                            const Curve2d& b, const Interval& domainB,
                                            const Context& /*ctx*/)
{
    const Ellipse2d* ea = a.Target<Ellipse2d>();
    const Ellipse2d* eb = b.Target<Ellipse2d>();
    if (ea == nullptr || eb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "EllipseEllipse2d: dispatched for other curve types");
    }
    return ConicConic(ConicOf(*ea), domainA, ConicOf(*eb), domainB);
}

} // namespace Stereon::Detail
