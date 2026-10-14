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
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;

/// The B-spline curve operations that are not evaluators: checked
/// construction, knot insertion, Bezier extraction, splitting,
/// reparametrisation and projection. They work on the homogeneous form of
/// the curve, (w P, w), in which a rational curve is a polynomial one and
/// every algorithm of The NURBS Book chapter 5 applies unchanged.

namespace Stereon {

namespace {

// ---------------------------------------------------------------------------
// Homogeneous form

/// A pole in homogeneous coordinates: (w x, w y[, w z], w).
template <std::size_t D>
using HPoint = std::array<double, D + 1>;

template <std::size_t D>
struct Homogeneous {
    int Degree = 1;
    std::vector<double> Knots;
    std::vector<HPoint<D>> Poles;
    bool Rational = false;
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
    std::vector<HPoint<D>> work(h.Poles.begin() + static_cast<std::ptrdiff_t>(k - p),
                                h.Poles.begin() + static_cast<std::ptrdiff_t>(k - s + 1));
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

template <std::size_t D>
[[nodiscard]] std::vector<Segment<D>> BezierSegmentsOf(const BasicBSplineCurve<double, D>& curve)
{
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
// Bernstein polynomials on [0, 1], for the projection

using Bernstein = std::vector<double>;

[[nodiscard]] double Binomial(const std::size_t n, const std::size_t k)
{
    double c = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        c = c * static_cast<double>(n - k + i) / static_cast<double>(i);
    }
    return c;
}

/// The product of two Bernstein polynomials of degrees m and n, in the
/// Bernstein basis of degree m + n: c_k = sum_{i+j=k} C(m,i) C(n,j) / C(m+n,k)
/// a_i b_j.
[[nodiscard]] Bernstein Product(const Bernstein& a, const Bernstein& b)
{
    const std::size_t m = a.size() - 1;
    const std::size_t n = b.size() - 1;
    Bernstein c(m + n + 1, 0.0);
    for (std::size_t i = 0; i <= m; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            c[i + j] += Binomial(m, i) * Binomial(n, j) * a[i] * b[j];
        }
    }
    for (std::size_t k = 0; k <= m + n; ++k) {
        c[k] /= Binomial(m + n, k);
    }
    return c;
}

/// The derivative, in the Bernstein basis one degree lower:
/// m (a_{i+1} - a_i).
[[nodiscard]] Bernstein Derivative(const Bernstein& a)
{
    const std::size_t m = a.size() - 1;
    Bernstein d(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        d[i] = static_cast<double>(m) * (a[i + 1] - a[i]);
    }
    return d;
}

/// The value at s by de Casteljau's algorithm.
[[nodiscard]] double BernsteinValue(Bernstein a, const double s)
{
    if (a.empty()) {
        return 0.0; // never called so; tells the optimiser the data exists
    }
    const std::size_t m = a.size() - 1;
    for (std::size_t r = 1; r <= m; ++r) {
        for (std::size_t i = 0; i + r <= m; ++i) {
            a[i] = (1.0 - s) * a[i] + s * a[i + 1];
        }
    }
    return a[0];
}

/// The coefficients of the two halves [0, 1/2] and [1/2, 1], rescaled to
/// [0, 1], read off the sides of de Casteljau's triangle; returns the value
/// at 1/2.
double Halves(Bernstein a, Bernstein& left, Bernstein& right)
{
    if (a.empty()) {
        return 0.0; // as in BernsteinValue
    }
    const std::size_t m = a.size() - 1;
    left.assign(m + 1, 0.0);
    right.assign(m + 1, 0.0);
    left[0] = a[0];
    right[m] = a[m];
    for (std::size_t r = 1; r <= m; ++r) {
        for (std::size_t i = 0; i + r <= m; ++i) {
            a[i] = 0.5 * (a[i] + a[i + 1]);
        }
        left[r] = a[0];
        right[m - r] = a[m - r];
    }
    return a[0];
}

/// Sign changes among the non-zero coefficients: an upper bound on the
/// number of roots in (0, 1), of the same parity (Descartes' rule of signs
/// in the Bernstein basis).
[[nodiscard]] std::size_t SignChanges(const Bernstein& a)
{
    std::size_t changes = 0;
    double previous = 0.0;
    for (const double c : a) {
        if (c != 0.0) {
            changes += (previous != 0.0 && (c < 0.0) != (previous < 0.0)) ? 1 : 0;
            previous = c;
        }
    }
    return changes;
}

/// Bisection steps that always reach two adjacent doubles in [0, 1], one
/// per binary digit between the smallest subnormal and 1.
constexpr int BisectionSteps =
    std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent;

/// Subdivisions of a parameter interval before it is as narrow as the
/// precision of a double allows: one per bit of the significand.
constexpr int SubdivisionDepth = std::numeric_limits<double>::digits;

/// The root in (0, 1) of a polynomial whose end values a(0), a(1) have
/// opposite signs, by bisection to adjacent doubles.
[[nodiscard]] double BisectRoot(const Bernstein& a)
{
    double low = 0.0;
    double high = 1.0;
    const bool lowNegative = a.front() < 0.0;
    for (int i = 0; i < BisectionSteps; ++i) {
        const double mid = 0.5 * (low + high);
        if (mid == low || mid == high) {
            break;
        }
        const double value = BernsteinValue(a, mid);
        if (value == 0.0) {
            return mid;
        }
        if ((value < 0.0) == lowNegative) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

/// Appends the roots in [low, high] of the polynomial with Bernstein
/// coefficients a on that interval. No sign change: no root. One sign
/// change with opposite end values: exactly one, found by bisection.
/// Otherwise the interval is halved; a root at the midpoint is recorded
/// there, since neither half would see it as an interior sign change.
void IsolateRoots(const Bernstein& a, const double low, const double high, const int depth,
                  std::vector<double>& roots)
{
    const std::size_t changes = SignChanges(a);
    if (changes == 0) {
        return;
    }
    if (changes == 1 && a.front() * a.back() < 0.0) {
        roots.push_back(low + BisectRoot(a) * (high - low));
        return;
    }
    const double mid = 0.5 * (low + high);
    if (depth >= SubdivisionDepth || mid == low || mid == high) {
        roots.push_back(mid);
        return;
    }
    Bernstein left;
    Bernstein right;
    if (Halves(a, left, right) == 0.0) {
        roots.push_back(mid);
    }
    IsolateRoots(left, low, mid, depth + 1, roots);
    IsolateRoots(right, mid, high, depth + 1, roots);
}

/// The numerator of the derivative of the squared distance from p to a
/// rational Bezier segment, in Bernstein form on [0, 1]. With
/// A(s) = sum B_i w_i P_i, w(s) = sum B_i w_i and g = |A - w p|^2,
/// d/ds |A/w - p|^2 = (g' w - 2 g w') / w^3, and w > 0, so the stationary
/// parameters are the roots of h = g' w - 2 g w', of degree 3p - 1.
template <std::size_t D>
[[nodiscard]] Bernstein StationarityPolynomial(const Segment<D>& segment,
                                               const BasicPointN<double, D>& p)
{
    const std::size_t count = segment.Poles.size();
    Bernstein w(count);
    for (std::size_t i = 0; i < count; ++i) {
        w[i] = segment.Poles[i][D];
    }
    Bernstein g(2 * count - 1, 0.0);
    for (std::size_t d = 0; d < D; ++d) {
        Bernstein offset(count);
        for (std::size_t i = 0; i < count; ++i) {
            offset[i] = segment.Poles[i][d] - w[i] * p[d];
        }
        const Bernstein square = Product(offset, offset);
        for (std::size_t k = 0; k < g.size(); ++k) {
            g[k] += square[k];
        }
    }
    Bernstein h = Product(Derivative(g), w);
    const Bernstein correction = Product(g, Derivative(w));
    for (std::size_t k = 0; k < h.size(); ++k) {
        h[k] -= 2.0 * correction[k];
    }
    return h;
}

/// The squared distance from p to the box of the segment's poles: a lower
/// bound on its squared distance to the segment, which lies in their
/// convex hull.
template <std::size_t D>
[[nodiscard]] double LowerBoundSquared(const Segment<D>& segment, const BasicPointN<double, D>& p)
{
    double total = 0.0;
    for (std::size_t d = 0; d < D; ++d) {
        double low = std::numeric_limits<double>::infinity();
        double high = -low;
        for (const HPoint<D>& q : segment.Poles) {
            low = std::min(low, q[d] / q[D]);
            high = std::max(high, q[d] / q[D]);
        }
        const double gap = std::max({low - p[d], p[d] - high, 0.0});
        total += gap * gap;
    }
    return total;
}

/// The best parameter found so far and its squared distance.
struct Closest {
    double Parameter;
    double DistanceSquared;
};

template <std::size_t D>
void Consider(const BasicBSplineCurve<double, D>& curve, const BasicPointN<double, D>& p,
              const double t, Closest& best)
{
    const double d = DistanceSquared(Eval(curve, t), p);
    if (d < best.DistanceSquared || (d == best.DistanceSquared && t < best.Parameter)) {
        best = {t, d};
    }
}

template <std::size_t D>
[[nodiscard]] double ProjectOnCurve(const BasicBSplineCurve<double, D>& curve,
                                    const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(BSplineCurve)", "a non-finite point");
    const double start = Domain(curve).GetLower();
    Closest best{start, DistanceSquared(Eval(curve, start), p)};
    for (const Segment<D>& segment : BezierSegmentsOf(curve)) {
        if (LowerBoundSquared(segment, p) > best.DistanceSquared) {
            continue;
        }
        Consider(curve, p, segment.End, best);
        std::vector<double> roots;
        IsolateRoots(StationarityPolynomial(segment, p), 0.0, 1.0, 0, roots);
        for (const double s : roots) {
            const double t = segment.Start + s * (segment.End - segment.Start);
            Consider(curve, p, std::clamp(t, segment.Start, segment.End), best);
        }
    }
    return best.Parameter;
}

// ---------------------------------------------------------------------------
// Checked construction and the knot operations

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
MakeCurve(const int degree, std::vector<double> knots, std::vector<BasicPointN<double, D>> poles,
          std::vector<double> weights)
{
    constexpr std::string_view where = "MakeBSplineCurve";
    STEREON_TRY_VOID(Detail::CheckKnotVector(degree, knots, poles.size(), where));
    for (const BasicPointN<double, D>& pole : poles) {
        STEREON_TRY_VOID(Detail::CheckFinite(pole, "every pole", where));
    }
    STEREON_TRY_VOID(Detail::CheckWeights(weights, poles.size(), where));
    return BasicBSplineCurve<double, D>{degree, std::move(knots), std::move(poles),
                                        std::move(weights)};
}

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
InsertKnotN(const BasicBSplineCurve<double, D>& curve, const double t, const int times)
{
    if (!std::isfinite(t) || !Contains(Domain(curve), t)) {
        return Fail(ErrorCode::InvalidInput, "InsertKnot: the knot must lie in the domain");
    }
    const std::size_t s = KnotMultiplicity(curve.Knots, t);
    if (times < 1 || s + static_cast<std::size_t>(times) > static_cast<std::size_t>(curve.Degree)) {
        return Fail(ErrorCode::InvalidInput,
                    "InsertKnot: times must be positive and keep the multiplicity at most the "
                    "degree");
    }
    Homogeneous<D> h = ToHomogeneous(curve);
    InsertKnotInPlace(h, t, static_cast<std::size_t>(times));
    return FromHomogeneous(std::move(h));
}

template <std::size_t D>
[[nodiscard]] std::vector<BasicBSplineCurve<double, D>>
BezierCurvesOf(const BasicBSplineCurve<double, D>& curve)
{
    std::vector<BasicBSplineCurve<double, D>> curves;
    for (const Segment<D>& segment : BezierSegmentsOf(curve)) {
        curves.push_back(SegmentCurve(segment, curve.Degree, IsRational(curve)));
    }
    return curves;
}

template <std::size_t D>
[[nodiscard]] Result<std::pair<BasicBSplineCurve<double, D>, BasicBSplineCurve<double, D>>>
SplitN(const BasicBSplineCurve<double, D>& curve, const double t)
{
    const Interval domain = Domain(curve);
    if (!std::isfinite(t) || !(t > domain.GetLower() && t < domain.GetUpper())) {
        return Fail(ErrorCode::InvalidInput, "Split: the parameter must lie inside the domain");
    }
    const auto p = static_cast<std::size_t>(curve.Degree);
    Homogeneous<D> h = ToHomogeneous(curve);
    const std::size_t s = KnotMultiplicity(h.Knots, t);
    if (s < p) {
        InsertKnotInPlace(h, t, p - s);
    }
    // t now occupies knots k - p + 1, ..., k.
    const auto k =
        static_cast<std::size_t>(std::ranges::upper_bound(h.Knots, t) - h.Knots.begin()) - 1;
    Homogeneous<D> left{h.Degree,
                        {h.Knots.begin(), h.Knots.begin() + static_cast<std::ptrdiff_t>(k + 1)},
                        {h.Poles.begin(), h.Poles.begin() + static_cast<std::ptrdiff_t>(k - p + 1)},
                        h.Rational};
    left.Knots.push_back(t);
    Homogeneous<D> right{h.Degree,
                         {t},
                         {h.Poles.begin() + static_cast<std::ptrdiff_t>(k - p), h.Poles.end()},
                         h.Rational};
    right.Knots.insert(right.Knots.end(), h.Knots.begin() + static_cast<std::ptrdiff_t>(k - p + 1),
                       h.Knots.end());
    return std::pair{FromHomogeneous(std::move(left)), FromHomogeneous(std::move(right))};
}

template <std::size_t D>
[[nodiscard]] Result<BasicBSplineCurve<double, D>>
ReparametrizedN(const BasicBSplineCurve<double, D>& curve, const Interval domain)
{
    const double low = domain.GetLower();
    const double high = domain.GetUpper();
    if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        return Fail(ErrorCode::InvalidInput,
                    "Reparametrized: the new domain must be finite and non-empty");
    }
    const Interval old = Domain(curve);
    const double scale = (high - low) / (old.GetUpper() - old.GetLower());
    BasicBSplineCurve<double, D> mapped = curve;
    for (double& u : mapped.Knots) {
        // The domain ends map exactly, so clamped ends stay clamped.
        if (u == old.GetUpper()) {
            u = high;
        } else if (u != old.GetLower()) {
            u = low + (u - old.GetLower()) * scale;
        } else {
            u = low;
        }
    }
    return mapped;
}

} // namespace

Result<BSplineCurve> MakeBSplineCurve(const int degree, std::vector<double> knots,
                                      std::vector<Point3> poles, std::vector<double> weights)
{
    return MakeCurve<3>(degree, std::move(knots), std::move(poles), std::move(weights));
}

Result<BSplineCurve2d> MakeBSplineCurve(const int degree, std::vector<double> knots,
                                        std::vector<Point2> poles, std::vector<double> weights)
{
    return MakeCurve<2>(degree, std::move(knots), std::move(poles), std::move(weights));
}

double Project(const BSplineCurve& curve, const Point3& p)
{
    return ProjectOnCurve(curve, p);
}

double Project(const BSplineCurve2d& curve, const Point2& p)
{
    return ProjectOnCurve(curve, p);
}

Result<BSplineCurve> InsertKnot(const BSplineCurve& curve, const double t, const int times)
{
    return InsertKnotN(curve, t, times);
}

Result<BSplineCurve2d> InsertKnot(const BSplineCurve2d& curve, const double t, const int times)
{
    return InsertKnotN(curve, t, times);
}

std::vector<BSplineCurve> ToBezierSegments(const BSplineCurve& curve)
{
    return BezierCurvesOf(curve);
}

std::vector<BSplineCurve2d> ToBezierSegments(const BSplineCurve2d& curve)
{
    return BezierCurvesOf(curve);
}

Result<std::pair<BSplineCurve, BSplineCurve>> Split(const BSplineCurve& curve, const double t)
{
    return SplitN(curve, t);
}

Result<std::pair<BSplineCurve2d, BSplineCurve2d>> Split(const BSplineCurve2d& curve, const double t)
{
    return SplitN(curve, t);
}

Result<BSplineCurve> Reparametrized(const BSplineCurve& curve, const Interval domain)
{
    return ReparametrizedN(curve, domain);
}

Result<BSplineCurve2d> Reparametrized(const BSplineCurve2d& curve, const Interval domain)
{
    return ReparametrizedN(curve, domain);
}

} // namespace Stereon
