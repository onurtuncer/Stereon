// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

module Stereon.Curve:Homogeneous;

import Stereon.Robust;
import :BSplineBasis;
import :BSplineCurve;

/// The homogeneous form of a B-spline curve, (w P, w), in which a rational
/// curve is a polynomial one and every algorithm of The NURBS Book chapter 5
/// applies unchanged, with the knot insertion and Bezier extraction built on
/// it. An implementation partition: shared by the implementation units of
/// the B-spline operations (BSplineCurve.cpp, BSplineOperations.cpp) and
/// exported from none of them.

namespace Stereon {

// ---------------------------------------------------------------------------
// Homogeneous form

/// A pole in homogeneous coordinates: (w x, w y[, w z], w).
template <std::size_t D>
using HPoint = std::array<double, D + 1>;

template <std::size_t D>
struct Homogeneous {
    int Degree;
    std::vector<double> Knots;
    std::vector<HPoint<D>> Poles;
    bool Rational;
};

template <std::size_t D>
[[nodiscard]] Homogeneous<D> ToHomogeneous(const BasicBSplineCurve<double, D>& curve)
{
    Homogeneous<D> h{curve.Degree, curve.Knots, {}, IsRational(curve)};
    h.Poles.reserve(curve.Poles.size());
    for (std::size_t i = 0; i < curve.Poles.size(); ++i) {
        const double w = Detail::WeightOf(curve, i);
        HPoint<D> q{};
        for (std::size_t d = 0; d < D; ++d) {
            q[d] = w * curve.Poles[i][d];
        }
        q[D] = w;
        h.Poles.push_back(q);
    }
    return h;
}

template <std::size_t D>
[[nodiscard]] BasicBSplineCurve<double, D> FromHomogeneous(Homogeneous<D> h)
{
    BasicBSplineCurve<double, D> curve{h.Degree, std::move(h.Knots), {}, {}};
    curve.Poles.reserve(h.Poles.size());
    for (const HPoint<D>& q : h.Poles) {
        BasicPointN<double, D> pole{};
        for (std::size_t d = 0; d < D; ++d) {
            pole[d] = h.Rational ? q[d] / q[D] : q[d];
        }
        curve.Poles.push_back(pole);
        if (h.Rational) {
            curve.Weights.push_back(q[D]);
        }
    }
    return curve;
}

/// alpha a + (1 - alpha) b. A polynomial curve keeps weight 1 exactly
/// rather than through rounding of alpha + (1 - alpha).
template <std::size_t D>
[[nodiscard]] HPoint<D> Blend(const HPoint<D>& a, const HPoint<D>& b, const double alpha,
                              const bool rational)
{
    HPoint<D> r{};
    for (std::size_t d = 0; d < D; ++d) {
        r[d] = alpha * a[d] + (1.0 - alpha) * b[d];
    }
    r[D] = rational ? alpha * a[D] + (1.0 - alpha) * b[D] : 1.0;
    return r;
}

/// The poles of the curve after inserting t `times` times, for knot index
/// k (the last with knots[k] <= t) and current multiplicity s; Boehm's
/// algorithm as A5.1 states it, with r = times <= p - s.
template <std::size_t D>
[[nodiscard]] std::vector<HPoint<D>> InsertedPoles(const Homogeneous<D>& h, const double t,
                                                   const std::size_t k, const std::size_t s,
                                                   const std::size_t r)
{
    const auto p = static_cast<std::size_t>(h.Degree);
    const std::size_t n = h.Poles.size() - 1;
    std::vector<HPoint<D>> q(n + 1 + r);
    for (std::size_t i = 0; i <= k - p; ++i) {
        q[i] = h.Poles[i];
    }
    for (std::size_t i = k - s; i <= n; ++i) {
        q[i + r] = h.Poles[i];
    }
    // The p - s + 1 poles the insertion touches, on the stack: at most
    // BSplineMaxDegree + 1 of them, 832 bytes in space, so an insertion
    // allocates its result and nothing else.
    std::array<HPoint<D>, static_cast<std::size_t>(BSplineMaxDegree) + 1>
        work; // written before read
    std::copy(h.Poles.begin() + static_cast<std::ptrdiff_t>(k - p),
              h.Poles.begin() + static_cast<std::ptrdiff_t>(k - s + 1), work.begin());
    for (std::size_t j = 1; j <= r; ++j) {
        const std::size_t l = k - p + j;
        for (std::size_t i = 0; i + j + s <= p; ++i) {
            const double alpha = (t - h.Knots[l + i]) / (h.Knots[i + k + 1] - h.Knots[l + i]);
            work[i] = Blend<D>(work[i + 1], work[i], alpha, h.Rational);
        }
        q[l] = work[0];
        q[k + r - j - s] = work[p - j - s];
    }
    const std::size_t l = k - p + r;
    for (std::size_t i = l + 1; i + s < k; ++i) {
        q[i] = work[i - l];
    }
    return q;
}

/// Inserts t `times` times; times + multiplicity(t) <= degree.
template <std::size_t D>
void InsertKnotInPlace(Homogeneous<D>& h, const double t, const std::size_t times)
{
    const auto upper = std::ranges::upper_bound(h.Knots, t);
    const auto k = static_cast<std::size_t>(upper - h.Knots.begin()) - 1;
    const std::size_t s = KnotMultiplicity(h.Knots, t);
    h.Poles = InsertedPoles(h, t, k, s, times);
    h.Knots.insert(h.Knots.begin() + static_cast<std::ptrdiff_t>(k + 1), times, t);
}

/// Raises every knot value in [low, high] to multiplicity `target`, at most
/// the degree, so that each span between them is a Bezier segment.
template <std::size_t D>
void RaiseMultiplicities(Homogeneous<D>& h, const double low, const double high,
                         const std::size_t target)
{
    std::vector<double> values;
    for (const double u : h.Knots) {
        if (u >= low && u <= high && (values.empty() || values.back() != u)) {
            values.push_back(u);
        }
    }
    for (const double u : values) {
        const std::size_t s = KnotMultiplicity(h.Knots, u);
        if (s < target) {
            InsertKnotInPlace(h, u, target - s);
        }
    }
}

/// One Bezier segment in homogeneous form, over [Start, End].
template <std::size_t D>
struct Segment {
    double Start;
    double End;
    std::vector<HPoint<D>> Poles;
};

/// The Bezier segments of a curve clamped at both ends in one pass
/// (algorithm A5.6): walking the distinct interior knots, each is raised to
/// multiplicity p by updating only the current segment's poles in place,
/// and the poles it shares with the next segment are copied forward. No
/// pole array is reallocated, unlike repeated knot insertion.
template <std::size_t D>
[[nodiscard]] std::vector<Segment<D>> ClampedBezierSegments(const Homogeneous<D>& h)
{
    const auto p = static_cast<std::size_t>(h.Degree);
    const std::vector<double>& u = h.Knots;
    const std::size_t m = u.size() - 1;
    std::vector<Segment<D>> segments;
    Segment<D> current{
        u[p], 0.0, {h.Poles.begin(), h.Poles.begin() + static_cast<std::ptrdiff_t>(p + 1)}};
    std::array<double, static_cast<std::size_t>(BSplineMaxDegree)> alphas; // written before read
    std::size_t a = p;
    std::size_t b = p + 1;
    while (b < m) {
        const std::size_t first = b;
        while (b < m && u[b + 1] == u[b]) {
            ++b;
        }
        const std::size_t mult = b - first + 1;
        if (b == m) {
            break; // the clamped end: the current segment is the last
        }
        current.End = u[b];
        Segment<D> next{u[b], 0.0, std::vector<HPoint<D>>(p + 1)};
        if (mult < p) {
            for (std::size_t j = p; j > mult; --j) {
                alphas[j - mult - 1] = (u[b] - u[a]) / (u[a + j] - u[a]);
            }
            const std::size_t r = p - mult;
            for (std::size_t j = 1; j <= r; ++j) {
                for (std::size_t k = p; k >= mult + j; --k) {
                    current.Poles[k] = Blend<D>(current.Poles[k], current.Poles[k - 1],
                                                alphas[k - mult - j], h.Rational);
                }
                next.Poles[r - j] = current.Poles[p];
            }
        }
        for (std::size_t i = p - std::min(mult, p); i <= p; ++i) {
            next.Poles[i] = h.Poles[b - p + i];
        }
        segments.push_back(std::move(current));
        current = std::move(next);
        a = b;
        ++b;
    }
    current.End = u[m];
    segments.push_back(std::move(current));
    return segments;
}

template <std::size_t D>
[[nodiscard]] std::vector<Segment<D>> BezierSegmentsOf(const BasicBSplineCurve<double, D>& curve)
{
    const auto degree = static_cast<std::size_t>(curve.Degree);
    const std::vector<double>& knots = curve.Knots;
    if (knots[0] == knots[degree] && knots[knots.size() - 1] == knots[knots.size() - 1 - degree]) {
        return ClampedBezierSegments(ToHomogeneous(curve));
    }
    const Interval domain = Domain(curve);
    const auto p = static_cast<std::size_t>(curve.Degree);
    Homogeneous<D> h = ToHomogeneous(curve);
    RaiseMultiplicities(h, domain.GetLower(), domain.GetUpper(), p);
    std::vector<Segment<D>> segments;
    for (std::size_t i = p; i + 1 < h.Poles.size() + 1; ++i) {
        const double a = h.Knots[i];
        const double b = h.Knots[i + 1];
        if (a < b && a >= domain.GetLower() && b <= domain.GetUpper()) {
            segments.push_back({a,
                                b,
                                {h.Poles.begin() + static_cast<std::ptrdiff_t>(i - p),
                                 h.Poles.begin() + static_cast<std::ptrdiff_t>(i + 1)}});
        }
    }
    return segments;
}

template <std::size_t D>
[[nodiscard]] BasicBSplineCurve<double, D> SegmentCurve(const Segment<D>& segment, const int degree,
                                                        const bool rational)
{
    const auto count = static_cast<std::size_t>(degree) + 1;
    Homogeneous<D> h{degree, {}, segment.Poles, rational};
    h.Knots.assign(count, segment.Start);
    h.Knots.insert(h.Knots.end(), count, segment.End);
    return FromHomogeneous(std::move(h));
}

// ---------------------------------------------------------------------------
// The periodic form

/// The knots of the periodic form from the n + 1 knots v_0 <= ... <= v_n of
/// one period, v_n = v_0 + T: p copies below, shifted down by the period,
/// then the period, then p copies above. The shifts are formed as
/// a - (b - v) and b + (v - a), which keep the copies below at or below a
/// and those above at or above b whatever the rounding.
[[nodiscard]] inline std::vector<double> PeriodicKnotsOf(const std::span<const double> period,
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

/// The knot of the periodic form's infinite sequence at place i, counted
/// from the start of the domain: v_r + q T for i = q n + r, 0 <= r < n, with
/// v_r the stored knots of one period. Shifted copies are formed from the
/// nearer end of the domain, a - (b - v) below and b + (v - a) above, as
/// PeriodicKnotsOf forms them.
template <std::size_t D>
[[nodiscard]] double PeriodicKnotAt(const Homogeneous<D>& h, const std::ptrdiff_t i)
{
    const auto p = static_cast<std::ptrdiff_t>(h.Degree);
    const auto n = static_cast<std::ptrdiff_t>(h.Poles.size()) - p;
    const double a = h.Knots[static_cast<std::size_t>(p)];
    const double b = h.Knots[static_cast<std::size_t>(n + p)];
    const std::ptrdiff_t q = i >= 0 ? i / n : -((-i + n - 1) / n);
    const std::ptrdiff_t r = i - q * n;
    const double v = h.Knots[static_cast<std::size_t>(p + r)];
    const double width = b - a;
    if (q == 0) {
        return v;
    }
    if (q > 0) {
        return b + (v - a) + static_cast<double>(q - 1) * width;
    }
    return a - (b - v) - static_cast<double>(-q - 1) * width;
}

/// h, in the periodic form, written out over `periods` more periods on each
/// side: the same curve on its domain and its periodic continuation on the
/// periods around it, an unclamped curve whose domain is [a - k T, b + k T].
/// The stored knots keep their values; the copies are made non-decreasing
/// against the rounding of their shifts.
template <std::size_t D>
[[nodiscard]] Homogeneous<D> PeriodicExtension(const Homogeneous<D>& h, const std::size_t periods)
{
    const auto p = static_cast<std::size_t>(h.Degree);
    const std::size_t n = h.Poles.size() - p;
    const std::size_t shift = periods * n;
    Homogeneous<D> e{h.Degree, {}, {}, h.Rational};
    const std::size_t poles = h.Poles.size() + 2 * shift;
    e.Poles.reserve(poles);
    for (std::size_t j = 0; j < poles; ++j) {
        e.Poles.push_back(h.Poles[j % n]); // the shift is a whole number of periods
    }
    e.Knots.reserve(poles + p + 1);
    for (std::size_t j = 0; j < poles + p + 1; ++j) {
        const auto w = static_cast<std::ptrdiff_t>(j) - static_cast<std::ptrdiff_t>(shift);
        const double knot = w >= 0 && static_cast<std::size_t>(w) < h.Knots.size()
                                ? h.Knots[static_cast<std::size_t>(w)]
                                : PeriodicKnotAt(h, w - static_cast<std::ptrdiff_t>(p));
        e.Knots.push_back(e.Knots.empty() ? knot : std::max(knot, e.Knots.back()));
    }
    return e;
}

/// h, in the periodic form, with t of its domain [a, b) inserted `times`
/// times at every one of its images a period apart, and so in the form
/// again with `times` more distinct poles. The form is written out over
/// enough periods that the window of the result lies well inside, every
/// image is inserted (Boehm's algorithm), and one period of knots and poles
/// is taken from the middle, the window built from it as PeriodicKnotsOf
/// builds one: the last p poles repeat the first p exactly, whatever the
/// rounding of the shifted copies. The images are formed from the nearer
/// end of the domain, so those below stay below a and those above at or
/// above b; t must be far enough below b that a - (b - t) < a, which the
/// caller ensures by taking t within rounding of b as a.
template <std::size_t D>
[[nodiscard]] Homogeneous<D> PeriodicInserted(const Homogeneous<D>& h, const double t,
                                              const std::size_t times)
{
    const auto p = static_cast<std::size_t>(h.Degree);
    const std::size_t n = h.Poles.size() - p;
    const double a = h.Knots[p];
    const double b = h.Knots[n + p];
    const double width = b - a;
    const std::size_t periods = p / n + 2;
    Homogeneous<D> e = PeriodicExtension(h, periods);
    for (std::size_t k = periods; k >= 1; --k) {
        InsertKnotInPlace(e, a - (b - t) - static_cast<double>(k - 1) * width, times);
    }
    InsertKnotInPlace(e, t, times);
    for (std::size_t k = 1; k <= periods; ++k) {
        InsertKnotInPlace(e, b + (t - a) + static_cast<double>(k - 1) * width, times);
    }
    // The first copy of a: its place in the extension, moved by the images
    // inserted below it; those at t and above go after it.
    const std::size_t count = n + times;
    const std::size_t start = periods * n + p + periods * times;
    const std::span<const double> period{e.Knots.data() + start, count + 1};
    Homogeneous<D> result{h.Degree, PeriodicKnotsOf(period, p), {}, h.Rational};
    result.Poles.reserve(count + p);
    for (std::size_t i = 0; i < count + p; ++i) {
        result.Poles.push_back(e.Poles[start - p + i % count]);
    }
    return result;
}

/// The restriction of a B-spline to a range inside its domain: split at the
/// range's ends where they are inside.
template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>> Restricted(BasicBSplineCurve<double, D> piece,
                                                              const Interval range)
{
    if (range.GetLower() > Domain(piece).GetLower()) {
        STEREON_TRY(auto halves, Split(piece, range.GetLower()));
        piece = std::move(halves.second);
    }
    if (range.GetUpper() < Domain(piece).GetUpper()) {
        STEREON_TRY(auto halves, Split(piece, range.GetUpper()));
        piece = std::move(halves.first);
    }
    return piece;
}

/// The piece of a curve in the periodic form over a range of at most one
/// period, anywhere: the range moved by whole periods to start
/// in the domain, cut from the form written out over one more period on
/// each side, and its knots moved back, the ends set to the range's own.
template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
PeriodicPiece(const BasicBSplineCurve<double, D>& curve, const Interval range, const double period)
{
    // The start in [a, b): b is the seam a, and taking a keeps the range's
    // end inside the extension, whose own ends are not clamped.
    double low = NearestInRange(range.GetLower(), Domain(curve), period);
    if (low >= Domain(curve).GetUpper()) {
        low = Domain(curve).GetLower();
    }
    const double shift = range.GetLower() - low;
    const BasicBSplineCurve<double, D> extended =
        FromHomogeneous(PeriodicExtension(ToHomogeneous(curve), 1));
    STEREON_TRY(auto piece, Restricted(extended, Interval{low, low + Width(range)}));
    const auto p = static_cast<std::size_t>(piece.Degree);
    for (std::size_t i = 0; i < piece.Knots.size(); ++i) {
        const bool start = i <= p;
        const bool end = i + p + 1 >= piece.Knots.size();
        piece.Knots[i] =
            start ? range.GetLower()
            : end ? range.GetUpper()
                  : std::clamp(piece.Knots[i] + shift, range.GetLower(), range.GetUpper());
    }
    return piece;
}

} // namespace Stereon
