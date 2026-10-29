// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <algorithm>
#include <array>
#include <cstddef>
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

} // namespace Stereon
