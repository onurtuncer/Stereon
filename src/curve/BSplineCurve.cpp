// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md): every module
// implementation unit of a geometry library includes it, because one that
// calls core's inline functions with contracts and parses none itself
// leaves __tu_has_violation undefined, which the macOS linker rejects.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Parallel.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <memory>
#include <span>
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
    std::vector<double> alphas(p);
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
    std::vector<double> rm(m + 1, 1.0);
    std::vector<double> rn(n + 1, 1.0);
    for (std::size_t i = 1; i <= m; ++i) {
        rm[i] = rm[i - 1] * static_cast<double>(m - i + 1) / static_cast<double>(i);
    }
    for (std::size_t j = 1; j <= n; ++j) {
        rn[j] = rn[j - 1] * static_cast<double>(n - j + 1) / static_cast<double>(j);
    }
    for (std::size_t i = 0; i <= m; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            c[i + j] += rm[i] * rn[j] * a[i] * b[j];
        }
    }
    double rmn = 1.0; // C(m + n, k), stepped along k
    for (std::size_t k = 0; k <= m + n; ++k) {
        c[k] /= rmn;
        rmn = rmn * static_cast<double>(m + n - k) / static_cast<double>(k + 1);
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

/// L(X) = X' w - 2 X w', the operator that turns the squared distance's
/// numerator g into the numerator of its derivative: with A = sum B_i w_i P_i
/// and w = sum B_i w_i on a segment, d/ds |A/w - p|^2 = L(g) / w^3 for
/// g = |A - w p|^2. L is linear, and
/// g = sum_d A_d^2 - 2 sum_d p_d (A_d w) + |p|^2 w^2, so the stationarity
/// polynomial h = L(g), of degree 3p - 1, is a combination of four
/// point-independent polynomials per segment with the point's coordinates as
/// weights. CurveProjector computes those once.
[[nodiscard]] Bernstein Stationary(const Bernstein& x, const Bernstein& w, const Bernstein& dw)
{
    Bernstein h = Product(Derivative(x), w);
    const Bernstein correction = Product(x, dw);
    for (std::size_t k = 0; k < h.size(); ++k) {
        h[k] -= 2.0 * correction[k];
    }
    return h;
}

} // namespace

} // namespace Stereon

namespace Stereon::Detail {

/// One Bezier segment, prepared for projection: its parameter range, the box
/// of its Euclidean poles, and the point-independent parts of its
/// stationarity polynomial (see Stationary).
struct ProjectionSegment {
    double Start = 0.0;
    double End = 0.0;
    std::array<double, 3> Low{};
    std::array<double, 3> High{};
    Bernstein Constant;              // L(sum_d A_d^2)
    std::array<Bernstein, 3> Linear; // L(A_d w)
    Bernstein Quadratic;             // L(w^2)
};

struct CurveProjectorData {
    BSplineCurve Curve;
    std::vector<ProjectionSegment> Segments;
};

} // namespace Stereon::Detail

namespace Stereon {

namespace {

[[nodiscard]] Detail::ProjectionSegment Prepared(const Segment<3>& segment)
{
    const std::size_t count = segment.Poles.size();
    Bernstein w(count);
    std::array<Bernstein, 3> a{Bernstein(count), Bernstein(count), Bernstein(count)};
    Detail::ProjectionSegment prepared{segment.Start, segment.End, {}, {}, {}, {}, {}};
    prepared.Low.fill(std::numeric_limits<double>::infinity());
    prepared.High.fill(-std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < count; ++i) {
        const HPoint<3>& q = segment.Poles[i];
        w[i] = q[3];
        for (std::size_t d = 0; d < 3; ++d) {
            a[d][i] = q[d];
            prepared.Low[d] = std::min(prepared.Low[d], q[d] / q[3]);
            prepared.High[d] = std::max(prepared.High[d], q[d] / q[3]);
        }
    }
    const Bernstein dw = Derivative(w);
    Bernstein squares = Product(a[0], a[0]);
    for (std::size_t d = 1; d < 3; ++d) {
        const Bernstein square = Product(a[d], a[d]);
        for (std::size_t k = 0; k < squares.size(); ++k) {
            squares[k] += square[k];
        }
    }
    prepared.Constant = Stationary(squares, w, dw);
    for (std::size_t d = 0; d < 3; ++d) {
        prepared.Linear[d] = Stationary(Product(a[d], w), w, dw);
    }
    prepared.Quadratic = Stationary(Product(w, w), w, dw);
    return prepared;
}

/// The squared distance from p to the segment's pole box: a lower bound on
/// its squared distance to the segment, which lies in their convex hull.
[[nodiscard]] double LowerBoundSquared(const Detail::ProjectionSegment& segment, const Point3& p)
{
    double total = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double gap = std::max({segment.Low[d] - p[d], p[d] - segment.High[d], 0.0});
        total += gap * gap;
    }
    return total;
}

/// The stationarity polynomial of the segment for the point p.
[[nodiscard]] Bernstein StationarityFor(const Detail::ProjectionSegment& segment, const Point3& p)
{
    Bernstein h = segment.Constant;
    const double square = p.X * p.X + p.Y * p.Y + p.Z * p.Z;
    for (std::size_t k = 0; k < h.size(); ++k) {
        h[k] += square * segment.Quadratic[k] -
                2.0 * (p.X * segment.Linear[0][k] + p.Y * segment.Linear[1][k] +
                       p.Z * segment.Linear[2][k]);
    }
    return h;
}

/// The best parameter found so far and its squared distance.
struct Closest {
    double Parameter;
    double DistanceSquared;
};

void Consider(const BSplineCurve& curve, const Point3& p, const double t, Closest& best)
{
    const double d = DistanceSquared(Eval(curve, t), p);
    if (d < best.DistanceSquared || (d == best.DistanceSquared && t < best.Parameter)) {
        best = {t, d};
    }
}

/// The segments in increasing order of their lower bound for p, ties by
/// index, so the search meets close candidates first and can stop at the
/// first segment that cannot improve on them.
[[nodiscard]] std::vector<std::pair<double, std::size_t>>
SegmentsByBound(const std::vector<Detail::ProjectionSegment>& segments, const Point3& p)
{
    std::vector<std::pair<double, std::size_t>> order;
    order.reserve(segments.size());
    for (std::size_t i = 0; i < segments.size(); ++i) {
        order.emplace_back(LowerBoundSquared(segments[i], p), i);
    }
    std::ranges::sort(order);
    return order;
}

[[nodiscard]] BSplineCurve Lifted(const BSplineCurve2d& curve)
{
    BSplineCurve lifted{curve.Degree, curve.Knots, {}, curve.Weights};
    for (const Point2& pole : curve.Poles) {
        lifted.Poles.push_back({pole.X, pole.Y, 0.0});
    }
    return lifted;
}

[[nodiscard]] std::shared_ptr<const Detail::CurveProjectorData> Prepare(BSplineCurve curve)
{
    auto data = std::make_shared<Detail::CurveProjectorData>();
    for (const Segment<3>& segment : BezierSegmentsOf(curve)) {
        data->Segments.push_back(Prepared(segment));
    }
    data->Curve = std::move(curve);
    return data;
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

} // namespace Stereon

namespace Stereon::Detail {

std::vector<double> PolynomialRoots(const std::span<const double> coefficients,
                                    const Interval range)
{
    // Taylor shift to t = a + h s on s in [0, 1]: d_j = sum_k c_k C(k, j)
    // a^(k-j) h^j, then the Bernstein coefficients b_i = sum_{j<=i}
    // C(i, j) / C(n, j) d_j.
    const double a = range.GetLower();
    const double h = range.GetUpper() - range.GetLower();
    const std::size_t n = coefficients.empty() ? 0 : coefficients.size() - 1;
    std::vector<double> shifted(n + 1, 0.0);
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            shifted[j] += coefficients[k] * Binomial(k, j) *
                          std::pow(a, static_cast<double>(k - j)) *
                          std::pow(h, static_cast<double>(j));
        }
    }
    Bernstein bernstein(n + 1, 0.0);
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            bernstein[i] += Binomial(i, j) / Binomial(n, j) * shifted[j];
        }
    }
    std::vector<double> local;
    if (n > 0) {
        IsolateRoots(bernstein, 0.0, 1.0, 0, local);
    }
    std::vector<double> roots;
    for (const double s : local) {
        roots.push_back(std::clamp(a + s * h, range.GetLower(), range.GetUpper()));
    }
    std::ranges::sort(roots);
    return roots;
}

} // namespace Stereon::Detail

namespace Stereon {

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

CurveProjector::CurveProjector(const BSplineCurve& curve)
    : m_Data(Prepare(curve))
{}

CurveProjector::CurveProjector(const BSplineCurve2d& curve)
    : m_Data(Prepare(Lifted(curve)))
{}

double CurveProjector::Project(const Point3& p) const
{
    Detail::RequireFinite(p, "Project(BSplineCurve)", "a non-finite point");
    const BSplineCurve& curve = m_Data->Curve;
    const double start = Domain(curve).GetLower();
    Closest best{start, DistanceSquared(Eval(curve, start), p)};
    std::vector<double> roots; // reused by every segment: one allocation per query
    for (const auto& [bound, index] : SegmentsByBound(m_Data->Segments, p)) {
        if (bound > best.DistanceSquared) {
            break; // every later segment is at least as far
        }
        const Detail::ProjectionSegment& segment = m_Data->Segments[index];
        Consider(curve, p, segment.End, best);
        roots.clear();
        IsolateRoots(StationarityFor(segment, p), 0.0, 1.0, 0, roots);
        for (const double r : roots) {
            const double t = segment.Start + r * (segment.End - segment.Start);
            Consider(curve, p, std::clamp(t, segment.Start, segment.End), best);
        }
    }
    return best.Parameter;
}

double CurveProjector::Project(const Point2& p) const
{
    return Project(Point3{p.X, p.Y, 0.0});
}

std::vector<double> ProjectAll(const CurveProjector& projector,
                               const std::span<const Point3> points, const Context& context)
{
    std::vector<double> parameters(points.size());
    struct Work {
        const CurveProjector* Projector;
        const Point3* Points;
        double* Parameters;
    } work{&projector, points.data(), parameters.data()};
    Detail::ParallelFor(points.size(), context.GetThreadCount(),
                        context.GetScheduler() == Scheduler::Serial, &work,
                        [](void* state, const std::size_t begin, const std::size_t end) {
                            const Work& w = *static_cast<const Work*>(state);
                            for (std::size_t i = begin; i < end; ++i) {
                                w.Parameters[i] = w.Projector->Project(w.Points[i]);
                            }
                        });
    return parameters;
}

double Project(const BSplineCurve& curve, const Point3& p)
{
    return CurveProjector(curve).Project(p);
}

double Project(const BSplineCurve2d& curve, const Point2& p)
{
    return CurveProjector(curve).Project(p);
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
