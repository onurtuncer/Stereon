// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

export module Stereon.Curve:Properties;

import Stereon.Robust;
import :BSplineCurve;
import :Concepts;

/// Local differential properties of curves (OCCT's GeomLProp_CLProps) and
/// arc length (GCPnts_AbscissaPoint), for any curve type: the built-in ones,
/// the erased Curve and Curve2d, and user types. Everything is computed
/// from Deriv, so a property is as exact as the derivatives it uses.
///
/// A property that is undefined at a parameter, such as the tangent where
/// the curve stops (zero speed) or the principal normal where it does not
/// bend, returns ErrorCode::DegenerateResult rather than an arbitrary
/// vector. The tests for zero are exact: near such points the property is
/// ill-conditioned, and callers that need a margin must apply their own
/// tolerance (ADR-0028), not inherit a hidden one.

export namespace Stereon {

/// The point type of a curve, Point2 or Point3, and its vector type.
template <typename C>
using CurvePointOf = decltype(Eval(std::declval<const C&>(), 0.0));

template <typename C>
using CurveVectorOf = decltype(Deriv(std::declval<const C&>(), 0.0, 1));

/// The unit tangent C'(t) / |C'(t)|. DegenerateResult where C'(t) = 0.
template <typename C>
[[nodiscard]] Result<CurveVectorOf<C>> UnitTangent(const C& curve, const double t)
{
    const CurveVectorOf<C> d1 = Deriv(curve, t, 1);
    const double speed = Length(d1);
    if (!(speed > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "UnitTangent: the curve has zero speed here");
    }
    return d1 / speed;
}

/// The curvature at t: |C' x C''| / |C'|^3 in space, which is never
/// negative, and the signed (C' x C'') / |C'|^3 in the plane, positive where
/// the curve turns counter-clockwise. DegenerateResult where C'(t) = 0.
template <typename C>
[[nodiscard]] Result<double> Curvature(const C& curve, const double t)
{
    const CurveVectorOf<C> d1 = Deriv(curve, t, 1);
    const CurveVectorOf<C> d2 = Deriv(curve, t, 2);
    const double speed = Length(d1);
    if (!(speed > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "Curvature: the curve has zero speed here");
    }
    const double cubed = speed * speed * speed;
    if constexpr (std::is_same_v<CurveVectorOf<C>, Vec2>) {
        return Cross(d1, d2) / cubed;
    } else {
        return Length(Cross(d1, d2)) / cubed;
    }
}

/// The unit principal normal of a curve in space, pointing towards the
/// centre of curvature: the component of C'' perpendicular to C',
/// normalised. DegenerateResult where the speed or the curvature is zero.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<Vec3> PrincipalNormal(const C& curve, const double t)
{
    const Vec3 d1 = Deriv(curve, t, 1);
    const Vec3 d2 = Deriv(curve, t, 2);
    const double speedSquared = LengthSquared(d1);
    if (!(speedSquared > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "PrincipalNormal: the curve has zero speed here");
    }
    const Vec3 normal = d2 - (Dot(d2, d1) / speedSquared) * d1;
    const double length = Length(normal);
    if (!(length > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "PrincipalNormal: the curve does not bend here (zero curvature)");
    }
    return normal / length;
}

/// The torsion of a curve in space, (C' x C'') . C''' / |C' x C''|^2:
/// how fast it leaves its osculating plane. DegenerateResult where the
/// curvature is zero, where the osculating plane is undefined.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<double> Torsion(const C& curve, const double t)
{
    const Vec3 binormal = Cross(Deriv(curve, t, 1), Deriv(curve, t, 2));
    const double squared = LengthSquared(binormal);
    if (!(squared > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "Torsion: the curve does not bend here (zero curvature)");
    }
    return Dot(binormal, Deriv(curve, t, 3)) / squared;
}

/// The centre of the osculating circle, C + N / kappa. DegenerateResult
/// where the speed or the curvature is zero.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<Point3> CenterOfCurvature(const C& curve, const double t)
{
    const Result<Vec3> normal = PrincipalNormal(curve, t);
    if (!normal) {
        return std::unexpected(normal.error());
    }
    return Eval(curve, t) + (1.0 / *Curvature(curve, t)) * *normal;
}

/// An arc length with the error estimate of the quadrature that produced
/// it (ADR-0028: the achieved bound is reported, not assumed).
struct ArcLengthResult {
    /// The length, in millimetres.
    double Length = 0.0;
    /// An estimate of its absolute error, in millimetres: the sum over the
    /// accepted pieces of the difference between the one-piece and
    /// two-piece Gauss-Legendre values. An estimate, not a certified bound.
    double ErrorEstimate = 0.0;
};

} // namespace Stereon

export namespace Stereon::Detail {

/// The points and weights of 16-point Gauss-Legendre quadrature on
/// [-1, 1], computed once by Newton's method on the Legendre polynomial
/// (defined in Properties.cpp), so that no 16-digit constants are written
/// out (ADR-0024).
[[nodiscard]] const std::array<std::pair<double, double>, 16>& GaussLegendre16() noexcept;

/// The relative accuracy the adaptive quadrature aims at, about 1e-12: a
/// piece is accepted when its one- and two-piece values agree this closely
/// relative to the length so far, which leaves several digits below the
/// precision of the evaluated derivatives. A solver target in the sense of
/// ADR-0028.
inline constexpr double ArcLengthRelativeTarget = 0x1p-40;

/// Halvings of a piece before the quadrature accepts it regardless; the
/// difference is still added to the error estimate.
inline constexpr int ArcLengthMaxDepth = 30;

template <typename C>
[[nodiscard]] double GaussLegendre(const C& curve, const double a, const double b)
{
    const double half = 0.5 * (b - a);
    const double middle = 0.5 * (a + b);
    double sum = 0.0;
    for (const auto& [node, weight] : GaussLegendre16()) {
        sum += weight * Length(Deriv(curve, middle + half * node, 1));
    }
    return half * sum;
}

/// Adaptive Gauss-Legendre on [a, b], given its one-piece value `whole`.
template <typename C>
void AdaptiveLength(const C& curve, const double a, const double b, const double whole,
                    const int depth, ArcLengthResult& result)
{
    const double mid = 0.5 * (a + b);
    const double left = GaussLegendre(curve, a, mid);
    const double right = GaussLegendre(curve, mid, b);
    const double difference = std::abs(left + right - whole);
    const double scale = std::max(result.Length + left + right, 0.0);
    if (difference <= ArcLengthRelativeTarget * scale || depth >= ArcLengthMaxDepth || mid == a ||
        mid == b) {
        result.Length += left + right;
        result.ErrorEstimate += difference;
        return;
    }
    AdaptiveLength(curve, a, mid, left, depth + 1, result);
    AdaptiveLength(curve, mid, b, right, depth + 1, result);
}

/// The break points of the quadrature in (low, high): a B-spline's
/// distinct knots, where its speed may have a kink; none for other types.
template <typename C>
[[nodiscard]] std::vector<double> LengthBreaks(const C& curve, const double low, const double high)
{
    std::vector<double> breaks{low};
    if constexpr (requires {
                      curve.Knots;
                      curve.Degree;
                  }) {
        for (const double u : curve.Knots) {
            if (u > breaks.back() && u < high) {
                breaks.push_back(u);
            }
        }
    }
    breaks.push_back(high);
    return breaks;
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The length of the curve between the parameters range.GetLower() and
/// range.GetUpper(), by adaptive 16-point Gauss-Legendre quadrature of the
/// speed |C'(t)|, split at a B-spline's knots. InvalidInput unless the
/// range is finite and inside the domain of a non-periodic curve; a range
/// of a periodic curve may be any finite range, and its length then counts
/// every turn it covers.
template <typename C>
[[nodiscard]] Result<ArcLengthResult> ArcLength(const C& curve, const Interval range)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    const bool finite = std::isfinite(low) && std::isfinite(high) && low <= high;
    bool inside = Contains(Domain(curve), range);
    if constexpr (PeriodicCurve<C>) {
        inside = inside || Period(curve) > 0.0;
    }
    if (!finite || !inside) {
        return Fail(ErrorCode::InvalidInput,
                    "ArcLength: the range must be finite and inside the curve's domain");
    }
    ArcLengthResult result;
    const std::vector<double> breaks = Detail::LengthBreaks(curve, low, high);
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double a = breaks[i];
        const double b = breaks[i + 1];
        if (a < b) {
            Detail::AdaptiveLength(curve, a, b, Detail::GaussLegendre(curve, a, b), 0, result);
        }
    }
    return result;
}

/// The parameter at which the arc length from `start` reaches `length`
/// (negative: backwards): the root of s(t) - length, where s is
/// non-decreasing, by Newton's method with ds/dt = |C'(t)|, kept inside a
/// bracket and falling back to bisection. A periodic curve may be followed
/// past the end of its domain, round as many turns as the length needs, and
/// so may a curve whose domain is unbounded. InvalidInput unless start lies
/// in the domain (any finite start for a periodic curve) and the curve is
/// long enough in that direction. Where the curve stops along a whole
/// stretch, the parameter is not unique and any one of the stretch may be
/// returned.
template <typename C>
[[nodiscard]] Result<double> ParameterAtArcLength(const C& curve, const double start,
                                                  const double length);

} // namespace Stereon

export namespace Stereon::Detail {

/// Newton steps for ParameterAtArcLength; the iteration converges
/// quadratically from inside the bracket, so the bound only limits the work
/// on degenerate input.
inline constexpr int ArcLengthNewtonSteps = 60;

/// The signed length from `start` to t.
template <typename C>
[[nodiscard]] double SignedLength(const C& curve, const double start, const double t)
{
    const Interval range{std::min(start, t), std::max(start, t)};
    const double length = ArcLength(curve, range)->Length;
    return t >= start ? length : -length;
}

template <typename C>
[[nodiscard]] double PeriodOrZero(const C& curve)
{
    if constexpr (PeriodicCurve<C>) {
        return Period(curve);
    } else {
        return 0.0;
    }
}

/// A parameter in the direction of `length` that the curve does not reach
/// before covering it: the domain end of a bounded curve, else a whole
/// number of periods or a doubling step. Empty if there is none.
template <typename C>
[[nodiscard]] std::optional<double> TravelEnd(const C& curve, const double start,
                                              const double length)
{
    const Interval domain = Domain(curve);
    const double bound = length > 0.0 ? domain.GetUpper() : domain.GetLower();
    const double period = PeriodOrZero(curve);
    if (period <= 0.0 && std::isfinite(bound)) {
        if (std::abs(SignedLength(curve, start, bound)) >= std::abs(length)) {
            return bound;
        }
        return std::nullopt;
    }
    const double sign = length > 0.0 ? 1.0 : -1.0;
    for (double step = period > 0.0 ? period : 1.0; std::isfinite(start + sign * step);
         step *= 2.0) {
        if (std::abs(SignedLength(curve, start, start + sign * step)) >= std::abs(length)) {
            return start + sign * step;
        }
    }
    return std::nullopt;
}

} // namespace Stereon::Detail

export namespace Stereon {

template <typename C>
Result<double> ParameterAtArcLength(const C& curve, const double start, const double length)
{
    const bool inside = Contains(Domain(curve), start) || Detail::PeriodOrZero(curve) > 0.0;
    if (!std::isfinite(start) || !std::isfinite(length) || !inside) {
        return Fail(ErrorCode::InvalidInput,
                    "ParameterAtArcLength: start must lie in the domain and length be finite");
    }
    if (length == 0.0) {
        return start;
    }
    const std::optional<double> end = Detail::TravelEnd(curve, start, length);
    if (!end) {
        return Fail(ErrorCode::InvalidInput,
                    "ParameterAtArcLength: the curve is not that long in that direction");
    }
    // s(t) - length changes sign on [low, high]; Newton's steps converge to
    // the root from inside it, and a step that would leave it bisects.
    double low = std::min(start, *end);
    double high = std::max(start, *end);
    const double settled = Detail::ArcLengthRelativeTarget * (high - low);
    double t = start + (*end - start) * (length / Detail::SignedLength(curve, start, *end));
    for (int step = 0; step < Detail::ArcLengthNewtonSteps; ++step) {
        const double residual = Detail::SignedLength(curve, start, t) - length;
        if (residual == 0.0) {
            return t;
        }
        (residual > 0.0 ? high : low) = t;
        const double speed = Length(Deriv(curve, t, 1));
        double next = speed > 0.0 ? t - residual / speed : 0.5 * (low + high);
        if (!(next > low && next < high)) {
            next = 0.5 * (low + high);
        }
        if (std::abs(next - t) <= settled) {
            return next;
        }
        t = next;
    }
    return t;
}

} // namespace Stereon
