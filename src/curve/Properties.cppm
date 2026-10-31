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
#include <limits>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

export module Stereon.Curve:Properties;

import Stereon.Robust;
import :BSplineCurve;
import :Concepts;
import :Visit;

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

} // namespace Stereon

export namespace Stereon::Detail {

/// The position vector and the derivatives up to order N at t, out[k] =
/// C^(k)(t): from one combined evaluation where the type has one
/// (Derivatives of a B-spline: one basis table for every order), else from
/// Eval and N calls of Deriv.
template <std::size_t N, typename C>
[[nodiscard]] std::array<CurveVectorOf<C>, N + 1> DerivativesTo(const C& curve, const double t)
{
    std::array<CurveVectorOf<C>, N + 1> out;
    if constexpr (requires(std::span<CurveVectorOf<C>> s) { Derivatives(curve, t, 2, s); }) {
        Derivatives(curve, t, static_cast<int>(N), std::span<CurveVectorOf<C>>{out});
    } else {
        out[0] = ToVec(Eval(curve, t));
        for (std::size_t k = 1; k <= N; ++k) {
            out[k] = Deriv(curve, t, static_cast<int>(k));
        }
    }
    return out;
}

/// The component of d2 perpendicular to d1, d1 != 0: the direction of the
/// principal normal, with length |d1 x d2| / |d1|.
template <typename V>
[[nodiscard]] V NormalPart(const V& d1, const V& d2)
{
    return d2 - (Dot(d2, d1) / LengthSquared(d1)) * d1;
}

} // namespace Stereon::Detail

export namespace Stereon {

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
    const auto d = Detail::DerivativesTo<2>(curve, t);
    const double speed = Length(d[1]);
    if (!(speed > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "Curvature: the curve has zero speed here");
    }
    const double cubed = speed * speed * speed;
    if constexpr (std::is_same_v<CurveVectorOf<C>, Vec2>) {
        return Cross(d[1], d[2]) / cubed;
    } else {
        return Length(Cross(d[1], d[2])) / cubed;
    }
}

/// The unit principal normal of a curve in space, pointing towards the
/// centre of curvature: the component of C'' perpendicular to C',
/// normalised. DegenerateResult where the speed or the curvature is zero.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<Vec3> PrincipalNormal(const C& curve, const double t)
{
    const auto d = Detail::DerivativesTo<2>(curve, t);
    if (!(LengthSquared(d[1]) > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "PrincipalNormal: the curve has zero speed here");
    }
    const Vec3 binormal = Cross(d[1], d[2]);
    const double bend = Length(binormal);
    if (!(bend > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "PrincipalNormal: the curve does not bend here (zero curvature)");
    }
    // The subtraction in NormalPart can leave roundoff parallel to the
    // tangent even when the derivatives are collinear. Use the same
    // curvature test and orthogonal construction as FrenetFrame.
    return Cross(binormal / bend, d[1] / Length(d[1]));
}

/// The torsion of a curve in space, (C' x C'') . C''' / |C' x C''|^2:
/// how fast it leaves its osculating plane. DegenerateResult where the
/// curvature is zero, where the osculating plane is undefined.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<double> Torsion(const C& curve, const double t)
{
    const auto d = Detail::DerivativesTo<3>(curve, t);
    const Vec3 binormal = Cross(d[1], d[2]);
    const double squared = LengthSquared(binormal);
    if (!(squared > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "Torsion: the curve does not bend here (zero curvature)");
    }
    return Dot(binormal, d[3]) / squared;
}

/// The Frenet frame of a curve in space and its curvature, all from one
/// evaluation of C, C' and C'': the counterpart of one SetParameter of
/// OCCT's GeomLProp_CLProps.
struct FrenetFrameResult {
    /// The unit tangent T.
    Vec3 Tangent;
    /// The unit principal normal N, towards the centre of curvature.
    Vec3 Normal;
    /// The unit binormal B = T x N, the normal of the osculating plane.
    Vec3 Binormal;
    /// The curvature, 1/mm, positive.
    double Curvature = 0.0;
};

/// The Frenet frame (T, N, B) and the curvature at t. DegenerateResult
/// where the speed or the curvature is zero, where N is undefined.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<FrenetFrameResult> FrenetFrame(const C& curve, const double t)
{
    const auto d = Detail::DerivativesTo<2>(curve, t);
    const double speed = Length(d[1]);
    if (!(speed > 0.0)) {
        return Fail(ErrorCode::DegenerateResult, "FrenetFrame: the curve has zero speed here");
    }
    const Vec3 binormal = Cross(d[1], d[2]);
    const double bend = Length(binormal);
    if (!(bend > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "FrenetFrame: the curve does not bend here (zero curvature)");
    }
    const Vec3 tangent = d[1] / speed;
    const Vec3 b = binormal / bend;
    return FrenetFrameResult{.Tangent = tangent,
                             .Normal = Cross(b, tangent),
                             .Binormal = b,
                             .Curvature = bend / (speed * speed * speed)};
}

/// The centre of the osculating circle, C + N / kappa. With n the part of
/// C'' perpendicular to C', N / kappa = n |C'|^4 / |C' x C''|^2, since
/// |n| = |C' x C''| / |C'|. DegenerateResult where the speed or the
/// curvature is zero.
template <typename C>
    requires std::is_same_v<CurveVectorOf<C>, Vec3>
[[nodiscard]] Result<Point3> CenterOfCurvature(const C& curve, const double t)
{
    const auto d = Detail::DerivativesTo<2>(curve, t);
    const double speedSquared = LengthSquared(d[1]);
    const double bend = LengthSquared(Cross(d[1], d[2]));
    if (!(speedSquared > 0.0) || !(bend > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "CenterOfCurvature: the curve has zero speed or zero curvature here");
    }
    return ToPoint(d[0]) + (speedSquared * speedSquared / bend) * Detail::NormalPart(d[1], d[2]);
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

template <typename C>
[[nodiscard]] double PeriodOrZero(const C& curve)
{
    if constexpr (PeriodicCurve<C>) {
        return Period(curve);
    } else {
        return 0.0;
    }
}

/// Appends to `knots` the images u + k T of the knot u, k an integer, that
/// lie strictly inside `range`; u itself when the period T is 0. The range
/// is at most one period wide, so the images to look at are the first at
/// or above its lower end and the next.
inline void AppendImages(const double u, const Interval range, const double period,
                         std::vector<double>& knots)
{
    const auto inside = [&](const double image) {
        return image > range.GetLower() && image < range.GetUpper();
    };
    if (!(period > 0.0)) {
        if (inside(u)) {
            knots.push_back(u);
        }
        return;
    }
    const double k = std::ceil((range.GetLower() - u) / period);
    for (const double image : {u + k * period, u + (k + 1.0) * period}) {
        if (inside(image)) {
            knots.push_back(image);
        }
    }
}

/// The distinct knots strictly inside `range`, sorted, where a B-spline's
/// speed may have a kink: a B-spline's own and, when it is periodic with
/// the given period, the images of one period's knots in the range (the
/// stored knots alone would leave a range outside the domain without its
/// kinks, and the quadrature with nothing to tell it where to look); those
/// of a B-spline inside an erased wrapper, found through VisitBuiltin; a
/// trimmed curve's basis's over the basis's range, in the trimmed
/// parameter; none for other types. The range is at most one period wide
/// (any width for a period of 0): PeriodicLength reduces a longer one.
template <typename C>
[[nodiscard]] std::vector<double> KnotsInside(const C& curve, const Interval range,
                                              const double period)
{
    std::vector<double> knots;
    if constexpr (requires {
                      curve.Knots;
                      curve.Degree;
                  }) {
        // In the periodic form, one period's knots [a, b) stand for all.
        std::span<const double> own = curve.Knots;
        if (period > 0.0) {
            const auto p = static_cast<std::size_t>(curve.Degree);
            own = own.subspan(p, curve.Poles.size() - p);
        }
        for (const double u : own) {
            AppendImages(u, range, period, knots);
        }
    } else if constexpr (requires {
                             curve.GetBasis();
                             curve.BasisParameter(0.0);
                         }) {
        // A trimmed curve: the basis's knots over the basis's range, which
        // may cross the seam of a periodic basis, taken to its own parameter.
        const double lower = curve.BasisParameter(range.GetLower());
        const double upper = curve.BasisParameter(range.GetUpper());
        for (const double u :
             KnotsInside(curve.GetBasis(), Interval{std::min(lower, upper), std::max(lower, upper)},
                         PeriodOrZero(curve.GetBasis()))) {
            const double t = curve.BasisParameter(u);
            if (t > range.GetLower() && t < range.GetUpper()) {
                knots.push_back(t);
            }
        }
    } else if constexpr (requires { curve.IsEmpty(); }) {
        VisitBuiltin(curve, [&](const auto& wrapped) {
            if constexpr (!std::is_same_v<std::remove_cvref_t<decltype(wrapped)>, C>) {
                knots = KnotsInside(wrapped, range, period);
            }
        });
    }
    std::ranges::sort(knots);
    const auto [first, last] = std::ranges::unique(knots);
    knots.erase(first, last);
    return knots;
}

/// The break points of the quadrature over the range: its ends and the
/// distinct knots between them, where a B-spline's speed may have a kink
/// that would slow a high-order rule down, or fool its error estimate.
template <typename C>
[[nodiscard]] std::vector<double> LengthBreaks(const C& curve, const Interval range,
                                               const double period)
{
    std::vector<double> breaks{range.GetLower()};
    for (const double u : KnotsInside(curve, range, period)) {
        if (u > breaks.back()) {
            breaks.push_back(u);
        }
    }
    breaks.push_back(range.GetUpper());
    return breaks;
}

/// The length over a range at most one period wide (any width for a period
/// of 0), by the adaptive quadrature between its break points.
template <typename C>
[[nodiscard]] ArcLengthResult QuadratureLength(const C& curve, const Interval range,
                                               const double period)
{
    ArcLengthResult result;
    const std::vector<double> breaks = LengthBreaks(curve, range, period);
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double a = breaks[i];
        const double b = breaks[i + 1];
        if (a < b) {
            AdaptiveLength(curve, a, b, GaussLegendre(curve, a, b), 0, result);
        }
    }
    return result;
}

/// The length over a range of a periodic curve that leaves the domain: the
/// whole periods it covers, each the length of one period over the domain,
/// and what remains, less than a period, from the range's start carried
/// into the domain, between the images of the knots it holds. Bounded work
/// for a range of any length, where a break at every knot image would not
/// be; the error estimate counts the period's estimate once per turn.
template <typename C>
[[nodiscard]] ArcLengthResult PeriodicLength(const C& curve, const Interval range,
                                             const double period)
{
    const double width = Width(range);
    double turns = std::floor(width / period);
    // width - turns * period in one rounding; the quotient may have rounded
    // up to a whole number.
    double rest = std::fma(-turns, period, width);
    if (rest < 0.0) {
        turns -= 1.0;
        rest += period;
    }
    ArcLengthResult result;
    if (turns > 0.0) {
        const ArcLengthResult turn = QuadratureLength(curve, Domain(curve), period);
        result.Length = turns * turn.Length;
        result.ErrorEstimate = turns * turn.ErrorEstimate;
    }
    if (rest > 0.0) {
        const double low = NearestInRange(range.GetLower(), Domain(curve), period);
        const ArcLengthResult tail = QuadratureLength(curve, Interval{low, low + rest}, period);
        result.Length += tail.Length;
        result.ErrorEstimate += tail.ErrorEstimate;
    }
    return result;
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The length of the curve between the parameters range.GetLower() and
/// range.GetUpper(), by adaptive 16-point Gauss-Legendre quadrature of the
/// speed |C'(t)|, split at a B-spline's knots. InvalidInput unless the
/// range is finite and inside the domain of a non-periodic curve; a range
/// of a periodic curve may be any finite range, and its length then counts
/// every turn it covers, split at the knots' images in every period.
template <typename C>
[[nodiscard]] Result<ArcLengthResult> ArcLength(const C& curve, const Interval range)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    const bool finite = std::isfinite(low) && std::isfinite(high) && low <= high;
    const bool inside = Contains(Domain(curve), range);
    const double period = Detail::PeriodOrZero(curve);
    if (!finite || !(inside || period > 0.0)) {
        return Fail(ErrorCode::InvalidInput,
                    "ArcLength: the range must be finite and inside the curve's domain");
    }
    if (!inside) {
        return Detail::PeriodicLength(curve, range, period);
    }
    return Detail::QuadratureLength(curve, range, period);
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
/// returned. Unresolved when floating-point stagnation or the iteration
/// budget prevents meeting the relative length residual target.
template <typename C>
[[nodiscard]] Result<double> ParameterAtArcLength(const C& curve, const double start,
                                                  const double length);

} // namespace Stereon

export namespace Stereon::Detail {

/// The relative residual at which ParameterAtArcLength stops: the length
/// to its parameter matches the request to about 1.5e-11, some 16 times
/// ArcLengthRelativeTarget, so the quadrature's own error cannot keep the
/// iteration going. The test is on the length, not the parameter: a
/// parameter step small for the whole range can still be large on a short
/// knot span where the curve is fast.
inline constexpr double ArcLengthResidualTarget = 0x1p-36;

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

/// A parameter and the signed length to it from start.
struct Station {
    double Parameter = 0.0;
    double Length = 0.0;
};

/// A parameter in the direction of `length` that the curve does not reach
/// before covering it: the domain end of a bounded curve, else a whole
/// number of periods or a doubling step. Empty if there is none. Each
/// doubling integrates only the stretch it adds, so the whole search
/// integrates its final range once.
template <typename C>
[[nodiscard]] std::optional<Station> TravelEnd(const C& curve, const double start,
                                               const double length)
{
    const Interval domain = Domain(curve);
    const double bound = length > 0.0 ? domain.GetUpper() : domain.GetLower();
    const double period = PeriodOrZero(curve);
    if (period <= 0.0 && std::isfinite(bound)) {
        const double reach = SignedLength(curve, start, bound);
        if (std::abs(reach) >= std::abs(length)) {
            return Station{bound, reach};
        }
        return std::nullopt;
    }
    const double sign = length > 0.0 ? 1.0 : -1.0;
    double reach = 0.0;
    double reached = start;
    for (double step = period > 0.0 ? period : 1.0; std::isfinite(start + sign * step);
         step *= 2.0) {
        const double end = start + sign * step;
        reach += SignedLength(curve, reached, end);
        reached = end;
        if (std::abs(reach) >= std::abs(length)) {
            return Station{end, reach};
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
    const std::optional<Detail::Station> end = Detail::TravelEnd(curve, start, length);
    if (!end) {
        return Fail(ErrorCode::InvalidInput,
                    "ParameterAtArcLength: the curve is not that long in that direction");
    }
    // s(t) - length changes sign on [low, high]; Newton's steps converge to
    // the root from inside it, and a step that would leave it bisects. The
    // length to each iterate is the length to the previous one plus the
    // short stretch between them, so no step integrates from start again.
    double low = std::min(start, end->Parameter);
    double high = std::max(start, end->Parameter);
    // The linear guess, kept inside the bracket against rounding.
    double t = std::clamp(start + (end->Parameter - start) * (length / end->Length), low, high);
    double reached = Detail::SignedLength(curve, start, t);
    for (int step = 0; step < Detail::ArcLengthNewtonSteps; ++step) {
        const double residual = reached - length;
        if (std::abs(residual) <= Detail::ArcLengthResidualTarget * std::abs(length)) {
            return t;
        }
        (residual > 0.0 ? high : low) = t;
        const double speed = Length(Deriv(curve, t, 1));
        double next = speed > 0.0 ? t - residual / speed : 0.5 * (low + high);
        if (!(next > low && next < high)) {
            next = 0.5 * (low + high);
        }
        if (next == t) {
            return Fail(
                ErrorCode::Unresolved,
                "ParameterAtArcLength: parameter stagnated before meeting the length target");
        }
        reached += Detail::SignedLength(curve, t, next);
        t = next;
    }
    // The final update has not yet had its residual tested.
    if (std::abs(reached - length) <= Detail::ArcLengthResidualTarget * std::abs(length)) {
        return t;
    }
    return Fail(
        ErrorCode::Unresolved,
        "ParameterAtArcLength: iteration budget exhausted before meeting the length target");
}

} // namespace Stereon
