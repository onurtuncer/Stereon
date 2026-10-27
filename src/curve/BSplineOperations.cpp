// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/BSplineCurve.cpp.
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;
import :BSplineOperations;
import :Homogeneous;

/// Clamping, the periodic form, degree elevation, knot removal and
/// concatenation of B-spline curves (The NURBS Book, chapter 5), on the
/// homogeneous form of Homogeneous.cppm. Knot removal is the one operation
/// here that can change the curve; everything else is exact up to rounding,
/// and the elevation uses the removal only where it is exact.

namespace Stereon {

namespace {

// ---------------------------------------------------------------------------
// Homogeneous poles

template <std::size_t D>
[[nodiscard]] std::size_t DegreeOf(const Homogeneous<D>& h)
{
    return static_cast<std::size_t>(h.Degree);
}

template <std::size_t D>
[[nodiscard]] BasicPointN<double, D> EuclideanOf(const HPoint<D>& q)
{
    BasicPointN<double, D> p{};
    for (std::size_t d = 0; d < D; ++d) {
        p[d] = q[d] / q[D];
    }
    return p;
}

/// |a - b| over all D + 1 coordinates.
template <std::size_t D>
[[nodiscard]] double HomogeneousDistance(const HPoint<D>& a, const HPoint<D>& b)
{
    double sum = 0.0;
    for (std::size_t d = 0; d <= D; ++d) {
        const double diff = a[d] - b[d];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

/// The pole q with Blend<D>(q, previous, alpha) == blended: the inverse of a
/// knot-insertion step from the front, (blended - (1 - alpha) previous) /
/// alpha. A polynomial curve keeps weight 1 exactly, as Blend does.
template <std::size_t D>
[[nodiscard]] HPoint<D> UnblendFront(const HPoint<D>& blended, const HPoint<D>& previous,
                                     const double alpha, const bool rational)
{
    HPoint<D> q{};
    for (std::size_t d = 0; d <= D; ++d) {
        q[d] = (blended[d] - (1.0 - alpha) * previous[d]) / alpha;
    }
    q[D] = rational ? q[D] : 1.0;
    return q;
}

/// The pole q with Blend<D>(next, q, alpha) == blended: (blended - alpha next)
/// / (1 - alpha).
template <std::size_t D>
[[nodiscard]] HPoint<D> UnblendBack(const HPoint<D>& blended, const HPoint<D>& next,
                                    const double alpha, const bool rational)
{
    HPoint<D> q{};
    for (std::size_t d = 0; d <= D; ++d) {
        q[d] = (blended[d] - alpha * next[d]) / (1.0 - alpha);
    }
    q[D] = rational ? q[D] : 1.0;
    return q;
}

/// The last index k with knots[k] <= t.
template <std::size_t D>
[[nodiscard]] std::size_t LastKnotIndex(const Homogeneous<D>& h, const double t)
{
    return static_cast<std::size_t>(std::ranges::upper_bound(h.Knots, t) - h.Knots.begin()) - 1;
}

/// The first index f with knots[f] >= t.
template <std::size_t D>
[[nodiscard]] std::size_t FirstKnotIndex(const Homogeneous<D>& h, const double t)
{
    return static_cast<std::size_t>(std::ranges::lower_bound(h.Knots, t) - h.Knots.begin());
}

template <std::size_t D>
[[nodiscard]] double LargestPoleNorm(const Homogeneous<D>& h)
{
    double largest = 0.0;
    for (const HPoint<D>& q : h.Poles) {
        largest = std::max(largest, Length(ToVec(EuclideanOf<D>(q))));
    }
    return largest;
}

template <std::size_t D>
[[nodiscard]] double SmallestWeight(const Homogeneous<D>& h)
{
    double smallest = std::numeric_limits<double>::infinity();
    for (const HPoint<D>& q : h.Poles) {
        smallest = std::min(smallest, q[D]);
    }
    return smallest;
}

// ---------------------------------------------------------------------------
// Clamping

/// Raises the start a of the domain to multiplicity p + 1 and drops the
/// knots and poles before it. With a at knots l - p + 1, ..., l, the basis
/// functions on the first non-empty span [a, u_{l+1}) depend on knots
/// l - p + 1 to l + p and poles l - p to l only, so knot l - p may become a
/// and everything before it may go.
template <std::size_t D>
void ClampStart(Homogeneous<D>& h)
{
    const std::size_t p = DegreeOf(h);
    const double a = h.Knots[p];
    const std::size_t s = KnotMultiplicity(h.Knots, a);
    if (s < p) {
        InsertKnotInPlace(h, a, p - s);
    }
    const auto drop = static_cast<std::ptrdiff_t>(LastKnotIndex(h, a) - p);
    h.Knots.erase(h.Knots.begin(), h.Knots.begin() + drop);
    h.Poles.erase(h.Poles.begin(), h.Poles.begin() + drop);
    h.Knots[0] = a;
}

/// The mirror of ClampStart at the end b of the domain: with b first at
/// knot f, the last non-empty span uses poles f - p - 1 to f - 1 and knots
/// up to f + p - 1, so the poles from f on go, knot f + p becomes b and the
/// knots after it go.
template <std::size_t D>
void ClampEnd(Homogeneous<D>& h)
{
    const std::size_t p = DegreeOf(h);
    const double b = h.Knots[h.Poles.size()];
    const std::size_t s = KnotMultiplicity(h.Knots, b);
    if (s < p) {
        InsertKnotInPlace(h, b, p - s);
    }
    const std::size_t f = FirstKnotIndex(h, b);
    h.Poles.resize(f);
    h.Knots.resize(f + p + 1);
    h.Knots[f + p] = b;
}

template <std::size_t D>
void Clamp(Homogeneous<D>& h)
{
    ClampStart(h);
    ClampEnd(h);
}

template <std::size_t D>
[[nodiscard]] bool IsClampedCurve(const BasicBSplineCurve<double, D>& curve)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    const std::size_t m = curve.Knots.size() - 1;
    return curve.Knots[0] == curve.Knots[p] && curve.Knots[m - p] == curve.Knots[m];
}

template <std::size_t D>
[[nodiscard]] BasicBSplineCurve<double, D> ClampedN(const BasicBSplineCurve<double, D>& curve)
{
    if (IsClampedCurve(curve)) {
        return curve;
    }
    Homogeneous<D> h = ToHomogeneous(curve);
    Clamp(h);
    return FromHomogeneous(std::move(h));
}

// ---------------------------------------------------------------------------
// Knot removal (A5.8), one copy at a time

/// Removes one copy of the knot at knots[r], its last index, of
/// multiplicity s: the poles P_{r-p}, ..., P_{r-s} that the insertion of
/// that copy would have produced are inverted from both ends of the run
/// towards its middle, where the two inversions meet, and one pole fewer
/// remains. Returns how far the two sides disagree at the meeting point,
/// in homogeneous coordinates: zero, up to rounding, when the knot is
/// removable exactly, and otherwise the distance by which one pole of the
/// curve with the knot re-inserted would miss the original.
template <std::size_t D>
[[nodiscard]] double RemoveOnce(Homogeneous<D>& h, const std::size_t r, const std::size_t s)
{
    const std::size_t p = DegreeOf(h);
    const double u = h.Knots[r];
    const std::size_t first = r - p;
    const std::size_t last = r - s;
    const std::size_t off = first - 1;
    std::vector<HPoint<D>> temp(last + 2 - off);
    temp[0] = h.Poles[off];
    temp[last + 1 - off] = h.Poles[last + 1];
    std::size_t i = first;
    std::size_t j = last;
    std::size_t ii = 1;
    std::size_t jj = last - off;
    while (j > i) {
        const double alfi = (u - h.Knots[i]) / (h.Knots[i + p + 1] - h.Knots[i]);
        const double alfj = (u - h.Knots[j]) / (h.Knots[j + p + 1] - h.Knots[j]);
        temp[ii] = UnblendFront<D>(h.Poles[i], temp[ii - 1], alfi, h.Rational);
        temp[jj] = UnblendBack<D>(h.Poles[j], temp[jj + 1], alfj, h.Rational);
        ++i;
        ++ii;
        --j;
        --jj;
    }
    double distance = 0.0;
    if (j < i) {
        distance = HomogeneousDistance<D>(temp[ii - 1], temp[jj + 1]);
    } else {
        const double alfi = (u - h.Knots[i]) / (h.Knots[i + p + 1] - h.Knots[i]);
        distance = HomogeneousDistance<D>(h.Poles[i],
                                          Blend<D>(temp[ii + 1], temp[ii - 1], alfi, h.Rational));
    }
    i = first;
    j = last;
    while (j > i) {
        h.Poles[i] = temp[i - off];
        h.Poles[j] = temp[j - off];
        ++i;
        --j;
    }
    h.Knots.erase(h.Knots.begin() + static_cast<std::ptrdiff_t>(r));
    h.Poles.erase(h.Poles.begin() + static_cast<std::ptrdiff_t>((first + last) / 2));
    return distance;
}

/// A bound on how far the curve moved when RemoveOnce reported `distance`:
/// one homogeneous pole of the re-inserted curve is off by at most that,
/// which moves a polynomial curve by at most the same, and a rational one,
/// C = A / w, by at most (|dA| + |C| |dw|) / w, bounded with the largest
/// pole norm of the curve before (|C| lies in the hull of its poles) and
/// the smallest weight after the removal.
template <std::size_t D>
[[nodiscard]] double DeviationBound(const Homogeneous<D>& after, const double largestNormBefore,
                                    const double distance)
{
    if (!after.Rational) {
        return distance;
    }
    return distance * (1.0 + largestNormBefore) / SmallestWeight(after);
}

template <std::size_t D>
[[nodiscard]] Result<BasicKnotRemoval<D>> RemoveKnotN(const BasicBSplineCurve<double, D>& curve,
                                                      const double t, const int times,
                                                      const double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        return Fail(ErrorCode::InvalidInput,
                    "RemoveKnot: the tolerance must be finite and non-negative");
    }
    const Interval domain = Domain(curve);
    if (!std::isfinite(t) || !(t > domain.GetLower() && t < domain.GetUpper())) {
        return Fail(ErrorCode::InvalidInput,
                    "RemoveKnot: the knot must lie strictly inside the domain");
    }
    const std::size_t s = KnotMultiplicity(curve.Knots, t);
    if (s == 0 || times < 1 || static_cast<std::size_t>(times) > s) {
        return Fail(ErrorCode::InvalidInput,
                    "RemoveKnot: times must be between 1 and the multiplicity of the knot");
    }
    Homogeneous<D> h = ToHomogeneous(curve);
    BasicKnotRemoval<D> result{};
    for (std::size_t k = 0; k < static_cast<std::size_t>(times); ++k) {
        Homogeneous<D> before = h;
        const double largest = LargestPoleNorm(before);
        const double distance = RemoveOnce(h, LastKnotIndex(h, t), s - k);
        const double bound = result.Deviation + DeviationBound(h, largest, distance);
        // Not removable within the tolerance, or not as a curve at all: the
        // inversion extrapolates, and can leave a rational curve a weight
        // that is not positive, with no bound to speak of.
        if (!(bound <= tolerance) || !(SmallestWeight(h) > 0.0)) {
            h = std::move(before);
            break;
        }
        result.Deviation = bound;
        ++result.Removed;
    }
    // Nothing removed: the input itself, not its round trip through the
    // homogeneous form, which rounds a rational curve's poles.
    result.Curve = result.Removed == 0 ? curve : FromHomogeneous(std::move(h));
    return result;
}

// ---------------------------------------------------------------------------
// Degree elevation (A5.9)

/// The Bezier poles raised by one degree, from n to n + 1:
/// Q_i = (i / (n + 1)) P_{i-1} + (1 - i / (n + 1)) P_i.
template <std::size_t D>
[[nodiscard]] std::vector<HPoint<D>> ElevatedOnce(const std::vector<HPoint<D>>& poles,
                                                  const bool rational)
{
    const std::size_t n = poles.size() - 1;
    std::vector<HPoint<D>> raised(n + 2);
    raised[0] = poles[0];
    raised[n + 1] = poles[n];
    for (std::size_t i = 1; i <= n; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(n + 1);
        raised[i] = Blend<D>(poles[i - 1], poles[i], alpha, rational);
    }
    return raised;
}

/// A distinct knot strictly inside the domain and its multiplicity.
struct InteriorKnot {
    double Value;
    std::size_t Multiplicity;
};

template <std::size_t D>
[[nodiscard]] std::vector<InteriorKnot> InteriorKnotsOf(const Homogeneous<D>& h)
{
    const double a = h.Knots[DegreeOf(h)];
    const double b = h.Knots[h.Poles.size()];
    std::vector<InteriorKnot> knots;
    for (const double u : h.Knots) {
        if (u > a && u < b) {
            if (knots.empty() || knots.back().Value != u) {
                knots.push_back({u, 0});
            }
            ++knots.back().Multiplicity;
        }
    }
    return knots;
}

/// The clamped curve h as one of `degree` >= its own: its Bezier segments
/// are raised one degree at a time and joined with every breakpoint at the
/// new degree's multiplicity q, then each interior knot of multiplicity s is
/// removed q - (s + q - p) = p - s times, which is exact because the
/// elevated curve is as smooth there as the original.
template <std::size_t D>
void ElevateClamped(Homogeneous<D>& h, const int degree)
{
    const std::size_t p = DegreeOf(h);
    const auto q = static_cast<std::size_t>(degree);
    if (q == p) {
        return;
    }
    const std::vector<InteriorKnot> interior = InteriorKnotsOf(h);
    const std::vector<Segment<D>> segments = ClampedBezierSegments(h);
    Homogeneous<D> raised{degree, {}, {}, h.Rational};
    raised.Knots.assign(q + 1, segments.front().Start);
    for (std::size_t k = 0; k < segments.size(); ++k) {
        std::vector<HPoint<D>> poles = segments[k].Poles;
        for (std::size_t e = p; e < q; ++e) {
            poles = ElevatedOnce<D>(poles, h.Rational);
        }
        if (k > 0) {
            raised.Knots.insert(raised.Knots.end(), q, segments[k].Start);
        }
        raised.Poles.insert(raised.Poles.end(), poles.begin() + (k > 0 ? 1 : 0), poles.end());
    }
    raised.Knots.insert(raised.Knots.end(), q + 1, segments.back().End);
    for (const InteriorKnot& knot : interior) {
        for (std::size_t s = knot.Multiplicity; s < p; ++s) {
            const std::size_t current = q - (s - knot.Multiplicity);
            static_cast<void>(RemoveOnce(raised, LastKnotIndex(raised, knot.Value), current));
        }
    }
    h = std::move(raised);
}

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
ElevatedDegreeN(const BasicBSplineCurve<double, D>& curve, const int degree)
{
    if (degree < curve.Degree || degree > BSplineMaxDegree) {
        return Fail(ErrorCode::InvalidInput,
                    "ElevatedDegree: the degree must be between the curve's and BSplineMaxDegree "
                    "(25)");
    }
    if (degree == curve.Degree) {
        return curve;
    }
    Homogeneous<D> h = ToHomogeneous(curve);
    Clamp(h);
    ElevateClamped(h, degree);
    return FromHomogeneous(std::move(h));
}

// ---------------------------------------------------------------------------
// Closure and the periodic form

template <std::size_t D>
[[nodiscard]] bool IsClosedN(const BasicBSplineCurve<double, D>& curve, const double tolerance)
{
    const Interval domain = Domain(curve);
    return Distance(Eval(curve, domain.GetLower()), Eval(curve, domain.GetUpper())) <= tolerance;
}

/// How far a knot of the periodic form may miss u_i + T, relative to the
/// size of the two knots compared: a few units of the rounding in forming
/// u_i +- T from knots of that size.
constexpr double PeriodicKnotSlack = 0x1p-48;

template <std::size_t D>
[[nodiscard]] bool IsPeriodicN(const BasicBSplineCurve<double, D>& curve)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    const std::size_t n = curve.Poles.size() - p;
    const double period = curve.Knots[n + p] - curve.Knots[p];
    for (std::size_t i = 0; i + n < curve.Knots.size(); ++i) {
        const double expected = curve.Knots[i] + period;
        const double slack = PeriodicKnotSlack * (std::abs(curve.Knots[i]) + std::abs(expected));
        if (std::abs(curve.Knots[i + n] - expected) > slack) {
            return false;
        }
    }
    for (std::size_t i = 0; i < p; ++i) {
        if (curve.Poles[n + i] != curve.Poles[i] ||
            Detail::WeightOf(curve, n + i) != Detail::WeightOf(curve, i)) {
            return false;
        }
    }
    return true;
}

/// The knots of the periodic form from the n + 1 knots v_0 <= ... <= v_n of
/// one period, v_n = v_0 + T: p copies below, shifted down by the period,
/// then the period, then p copies above. The shifts are formed as
/// a - (b - v) and b + (v - a), which keep the copies below at or below a
/// and those above at or above b whatever the rounding.
[[nodiscard]] std::vector<double> PeriodicKnotsOf(const std::span<const double> period,
                                                  const std::size_t p)
{
    const std::size_t n = period.size() - 1;
    const double a = period.front();
    const double b = period.back();
    std::vector<double> knots;
    knots.reserve(n + 2 * p + 1);
    for (std::size_t r = p; r >= 1; --r) {
        knots.push_back(a - (b - period[n - r]));
    }
    knots.insert(knots.end(), period.begin(), period.end());
    for (std::size_t r = 1; r <= p; ++r) {
        knots.push_back(b + (period[r] - a));
    }
    return knots;
}

/// The clamped closed curve h, with n + 1 poles whose ends have been merged
/// into `merged`, in the periodic form: the n distinct poles are
/// D_i = M_{(i - p + 1) mod n} with M_0 the merged pole, so that the seam
/// knot of multiplicity p puts D_{p-1} = M_0 at the start of the domain
/// and the curve on the domain keeps its poles in order; the last p repeat
/// the first p.
template <std::size_t D>
[[nodiscard]] Homogeneous<D> PeriodicFormOf(const Homogeneous<D>& h, const HPoint<D>& merged)
{
    const std::size_t p = DegreeOf(h);
    const std::size_t n = h.Poles.size() - 1;
    std::vector<double> period;
    period.reserve(n + 1);
    period.assign(p, h.Knots[p]);
    period.insert(period.end(), h.Knots.begin() + static_cast<std::ptrdiff_t>(p + 1),
                  h.Knots.begin() + static_cast<std::ptrdiff_t>(n + 1));
    period.push_back(h.Knots[n + 1]);
    Homogeneous<D> periodic{h.Degree, PeriodicKnotsOf(period, p), {}, h.Rational};
    periodic.Poles.reserve(n + p);
    for (std::size_t i = 0; i < n + p; ++i) {
        const std::size_t j = (i + n + 1 - p) % n;
        periodic.Poles.push_back(j == 0 ? merged : h.Poles[j]);
    }
    return periodic;
}

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
AsPeriodicN(const BasicBSplineCurve<double, D>& curve, const double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        return Fail(ErrorCode::InvalidInput,
                    "AsPeriodic: the tolerance must be finite and non-negative");
    }
    if (IsPeriodicN(curve)) {
        return curve;
    }
    Homogeneous<D> h = ToHomogeneous(curve);
    Clamp(h);
    const HPoint<D>& start = h.Poles.front();
    const HPoint<D>& end = h.Poles.back();
    if (h.Rational && start[D] != end[D]) {
        return Fail(ErrorCode::InvalidInput, "AsPeriodic: the weights of the two end poles differ");
    }
    if (Distance(EuclideanOf<D>(start), EuclideanOf<D>(end)) > tolerance) {
        return Fail(ErrorCode::InvalidInput,
                    "AsPeriodic: the curve does not close within the tolerance");
    }
    const HPoint<D> merged = Blend<D>(start, end, 0.5, h.Rational);
    return FromHomogeneous(PeriodicFormOf(h, merged));
}

// ---------------------------------------------------------------------------
// Concatenation

/// Scales every homogeneous pole of h by `factor`: the same curve with its
/// weights multiplied by the factor.
template <std::size_t D>
void ScaleWeights(Homogeneous<D>& h, const double factor)
{
    for (HPoint<D>& q : h.Poles) {
        for (double& c : q) {
            c *= factor;
        }
    }
}

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
ConcatenatedN(const BasicBSplineCurve<double, D>& first, const BasicBSplineCurve<double, D>& second,
              const double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        return Fail(ErrorCode::InvalidInput,
                    "Concatenated: the tolerance must be finite and non-negative");
    }
    Homogeneous<D> a = ToHomogeneous(first);
    Clamp(a);
    Homogeneous<D> b = ToHomogeneous(second);
    Clamp(b);
    if (Distance(EuclideanOf<D>(a.Poles.back()), EuclideanOf<D>(b.Poles.front())) > tolerance) {
        return Fail(ErrorCode::InvalidInput,
                    "Concatenated: the end of the first curve and the start of the second are "
                    "not within the tolerance");
    }
    const int degree = std::max(a.Degree, b.Degree);
    ElevateClamped(a, degree);
    ElevateClamped(b, degree);
    const bool rational = a.Rational || b.Rational;
    ScaleWeights(b, a.Poles.back()[D] / b.Poles.front()[D]);
    const auto q = static_cast<std::size_t>(degree);
    const double junction = a.Knots.back();
    const double start = b.Knots.front();
    Homogeneous<D> joined{degree, {}, {}, rational};
    joined.Knots.assign(a.Knots.begin(), a.Knots.end() - 1);
    for (std::size_t i = q + 1; i < b.Knots.size(); ++i) {
        joined.Knots.push_back(junction + (b.Knots[i] - start));
    }
    joined.Poles.assign(a.Poles.begin(), a.Poles.end() - 1);
    joined.Poles.push_back(Blend<D>(a.Poles.back(), b.Poles.front(), 0.5, rational));
    joined.Poles.insert(joined.Poles.end(), b.Poles.begin() + 1, b.Poles.end());
    return FromHomogeneous(std::move(joined));
}

} // namespace

// ---------------------------------------------------------------------------
// The exported overloads

bool IsClosed(const BSplineCurve& curve, const double tolerance)
{
    return IsClosedN(curve, tolerance);
}

bool IsClosed(const BSplineCurve2d& curve, const double tolerance)
{
    return IsClosedN(curve, tolerance);
}

bool IsPeriodic(const BSplineCurve& curve)
{
    return IsPeriodicN(curve);
}

bool IsPeriodic(const BSplineCurve2d& curve)
{
    return IsPeriodicN(curve);
}

BSplineCurve Clamped(const BSplineCurve& curve)
{
    return ClampedN(curve);
}

BSplineCurve2d Clamped(const BSplineCurve2d& curve)
{
    return ClampedN(curve);
}

Result<BSplineCurve> AsPeriodic(const BSplineCurve& curve, const double tolerance)
{
    return AsPeriodicN(curve, tolerance);
}

Result<BSplineCurve2d> AsPeriodic(const BSplineCurve2d& curve, const double tolerance)
{
    return AsPeriodicN(curve, tolerance);
}

Result<BSplineCurve> ElevatedDegree(const BSplineCurve& curve, const int degree)
{
    return ElevatedDegreeN(curve, degree);
}

Result<BSplineCurve2d> ElevatedDegree(const BSplineCurve2d& curve, const int degree)
{
    return ElevatedDegreeN(curve, degree);
}

Result<KnotRemoval> RemoveKnot(const BSplineCurve& curve, const double t, const int times,
                               const double tolerance)
{
    return RemoveKnotN(curve, t, times, tolerance);
}

Result<KnotRemoval2d> RemoveKnot(const BSplineCurve2d& curve, const double t, const int times,
                                 const double tolerance)
{
    return RemoveKnotN(curve, t, times, tolerance);
}

Result<BSplineCurve> Concatenated(const BSplineCurve& first, const BSplineCurve& second,
                                  const double tolerance)
{
    return ConcatenatedN(first, second, tolerance);
}

Result<BSplineCurve2d> Concatenated(const BSplineCurve2d& first, const BSplineCurve2d& second,
                                    const double tolerance)
{
    return ConcatenatedN(first, second, tolerance);
}

} // namespace Stereon
