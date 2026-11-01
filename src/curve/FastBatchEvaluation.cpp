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
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;
import :BatchEvaluation;

/// EvalAllFast: each non-empty span converted to a polynomial in a local
/// parameter centred on the span, then every point by Horner's rule
/// (docs/nurbs.rst, *Fast batch evaluation*). Its coefficients are the
/// homogeneous curve's Taylor coefficients at the span's midpoint in
/// u = (t - m) / h, which Detail::SpanCoefficients computes by de Boor
/// differencing in u itself, so that a narrow span cannot overflow them;
/// |u|^k <= 2^-k, so Horner's rule sums terms that shrink. A run of
/// Simd::Width<double> parameters in one span is
/// evaluated in the lanes of one vector by the scalar Horner steps, so lanes
/// and scalars agree bit for bit and the result does not depend on the SIMD
/// width. The overloads on a curve convert a span when a parameter first
/// lands in it, into a buffer on the stack, by the same function that
/// ToSpanPolynomials uses, so both give the same bits.

namespace Stereon {

namespace {

using Lanes = Simd::Vec<double>;
constexpr std::size_t LaneCount = Simd::Width<double>;
constexpr std::size_t MaxBasis = Detail::BSplineBasisCount;

/// One span's polynomial: its midpoint, 1 / width, and (Degree + 1) rows of
/// D + 1 coefficients.
struct SpanForm {
    double Midpoint = 0.0;
    double InverseWidth = 0.0;
    const double* Coefficients = nullptr;
};

/// The coefficients of the non-empty span [knots[span], knots[span + 1]],
/// midpoint m and width h, into `out`, (p + 1) (D + 1) values; work arrays
/// sized for the degree.
template <std::size_t Capacity, std::size_t D>
void CoefficientsWithCapacity(const BasicBSplineCurve<double, D>& curve, const std::size_t span,
                              const double m, const double h, double* out)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    std::array<BasicVecN<double, D>, Capacity> aw;
    std::array<double, Capacity> w;
    Detail::SpanCoefficients<Capacity>(curve, span, m, h, aw, w);
    for (std::size_t k = 0; k <= p; ++k) {
        for (std::size_t d = 0; d < D; ++d) {
            out[k * (D + 1) + d] = aw[k][d];
        }
        out[k * (D + 1) + D] = w[k];
    }
}

/// The polynomial of a non-empty span, its coefficients into `out`.
template <std::size_t D>
[[nodiscard]] SpanForm Convert(const BasicBSplineCurve<double, D>& curve, const std::size_t span,
                               double* out)
{
    const double a = curve.Knots[span];
    const double b = curve.Knots[span + 1];
    const double m = a + 0.5 * (b - a);
    const double h = b - a;
    if (static_cast<std::size_t>(curve.Degree) < Detail::SmallDegreeCapacity) {
        CoefficientsWithCapacity<Detail::SmallDegreeCapacity>(curve, span, m, h, out);
    } else {
        CoefficientsWithCapacity<MaxBasis>(curve, span, m, h, out);
    }
    return {m, 1.0 / h, out};
}

template <std::size_t D>
[[nodiscard]] BasicSpanPolynomials<D> SpanPolynomialsOf(const BasicBSplineCurve<double, D>& curve)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    const std::size_t n = curve.Poles.size();
    const std::size_t stride = (p + 1) * (D + 1);
    BasicSpanPolynomials<D> result;
    result.Degree = curve.Degree;
    result.Rational = IsRational(curve);
    result.Knots = curve.Knots;
    result.PoleCount = n;
    result.Slots.assign(curve.Knots.size(), n + 1);
    for (std::size_t span = p; span < n; ++span) {
        if (!(curve.Knots[span] < curve.Knots[span + 1])) {
            continue;
        }
        result.Slots[span] = result.Midpoints.size();
        result.Coefficients.resize(result.Coefficients.size() + stride);
        const SpanForm form =
            Convert(curve, span, result.Coefficients.data() + result.Coefficients.size() - stride);
        result.Midpoints.push_back(form.Midpoint);
        result.InverseWidths.push_back(form.InverseWidth);
    }
    return result;
}

/// The span FindSpan gives t, tried first against `guess`.
[[nodiscard]] std::size_t SpanOf(const std::span<const double> knots, const int degree,
                                 const std::size_t poles, const double t, const std::size_t guess)
{
    if (knots[guess] <= t && t < knots[guess + 1]) {
        return guess;
    }
    return FindSpan(knots, degree, poles, t);
}

/// The degree a kernel runs at: Fixed when it is not 0, so that the
/// compiler unrolls Horner's loop, otherwise the curve's.
template <std::size_t Fixed>
[[nodiscard]] constexpr std::size_t DegreeOf(const std::size_t runtime) noexcept
{
    return Fixed != 0 ? Fixed : runtime;
}

/// One point by Horner's rule.
template <std::size_t D, std::size_t Fixed, bool Rational>
[[nodiscard]] BasicPointN<double, D> PointAt(const std::size_t degree, const SpanForm& form,
                                             const double t)
{
    const std::size_t p = DegreeOf<Fixed>(degree);
    const double* c = form.Coefficients;
    const double u = (t - form.Midpoint) * form.InverseWidth;
    std::array<double, D + 1> sum;
    for (std::size_t d = 0; d <= D; ++d) {
        sum[d] = c[p * (D + 1) + d];
    }
    for (std::size_t k = p; k-- > 0;) {
        for (std::size_t d = 0; d <= D; ++d) {
            sum[d] = sum[d] * u + c[k * (D + 1) + d];
        }
    }
    BasicPointN<double, D> point;
    for (std::size_t d = 0; d < D; ++d) {
        point[d] = Rational ? sum[d] / sum[D] : sum[d];
    }
    return point;
}

/// The points at t[0, LaneCount), which share the span, in lanes, by
/// PointAt's operations in each.
template <std::size_t D, std::size_t Fixed, bool Rational>
void PointsInLanes(const std::size_t degree, const SpanForm& form, const double* t,
                   BasicPointN<double, D>* points)
{
    const std::size_t p = DegreeOf<Fixed>(degree);
    const double* c = form.Coefficients;
    const Lanes u =
        (Simd::Load(t) - Simd::Broadcast(form.Midpoint)) * Simd::Broadcast(form.InverseWidth);
    std::array<Lanes, D + 1> sum;
    for (std::size_t d = 0; d <= D; ++d) {
        sum[d] = Simd::Broadcast(c[p * (D + 1) + d]);
    }
    for (std::size_t k = p; k-- > 0;) {
        for (std::size_t d = 0; d <= D; ++d) {
            sum[d] = sum[d] * u + Simd::Broadcast(c[k * (D + 1) + d]);
        }
    }
    std::array<std::array<double, LaneCount>, D> coordinates; // filled by Store
    for (std::size_t d = 0; d < D; ++d) {
        const Lanes value = Rational ? sum[d] / sum[D] : sum[d];
        Simd::Store(value, coordinates[d].data());
    }
    for (std::size_t lane = 0; lane < LaneCount; ++lane) {
        for (std::size_t d = 0; d < D; ++d) {
            points[lane][d] = coordinates[d][lane];
        }
    }
}

/// The loop over the parameters, for one degree (Fixed, or the curve's
/// when Fixed is 0) and one kind of curve.
template <std::size_t D, std::size_t Fixed, bool Rational, typename Form>
void EvalKernel(const std::span<const double> knots, const int degree, const std::size_t poles,
                Form& form, const std::span<const double> parameters,
                const std::span<BasicPointN<double, D>> points)
{
    const auto p = static_cast<std::size_t>(degree);
    std::size_t guess = p; // the first span's index or below
    std::size_t k = 0;
    for (; k + LaneCount <= parameters.size(); k += LaneCount) {
        const double* t = parameters.data() + k;
        guess = SpanOf(knots, degree, poles, t[0], guess);
        bool shared = true;
        for (std::size_t lane = 1; lane < LaneCount && shared; ++lane) {
            shared = SpanOf(knots, degree, poles, t[lane], guess) == guess;
        }
        if (shared) {
            PointsInLanes<D, Fixed, Rational>(p, form(guess), t, points.data() + k);
            continue;
        }
        for (std::size_t lane = 0; lane < LaneCount; ++lane) {
            guess = SpanOf(knots, degree, poles, t[lane], guess);
            points[k + lane] = PointAt<D, Fixed, Rational>(p, form(guess), t[lane]);
        }
    }
    for (; k < parameters.size(); ++k) {
        guess = SpanOf(knots, degree, poles, parameters[k], guess);
        points[k] = PointAt<D, Fixed, Rational>(p, form(guess), parameters[k]);
    }
}

/// What both overloads share: the knots to find spans in, and `Form`, the
/// polynomial of a span by its knot index. Degrees below
/// Detail::SmallDegreeCapacity, which cover the curves met in practice, get
/// a kernel each, generated by an expansion statement (P1306); higher ones
/// share the kernel of a run-time degree. Every kernel performs the same
/// operations, so the choice changes no result.
template <std::size_t D, typename Form>
void EvalAllFastN(const std::span<const double> knots, const int degree, const std::size_t poles,
                  const bool rational, Form&& form, const std::span<const double> parameters,
                  const std::span<BasicPointN<double, D>> points)
{
    Detail::Require(points.size() == parameters.size(), "EvalAllFast(BSplineCurve)",
                    "as many points as parameters");
    const auto p = static_cast<std::size_t>(degree);
    if constexpr (Detail::EvaluatorChecks) {
        const Interval domain{knots[p], knots[poles]};
        for (const double t : parameters) {
            Detail::RequireFinite(t, "EvalAllFast(BSplineCurve)", "a non-finite parameter");
            Detail::Require(Contains(domain, t), "EvalAllFast(BSplineCurve)",
                            "a parameter outside the domain");
        }
    }
    bool done = false;
    template for (constexpr std::size_t fixed : Detail::Indices<Detail::SmallDegreeCapacity>)
    {
        if (fixed != 0 && !done && p == fixed) {
            if (rational) {
                EvalKernel<D, fixed, true>(knots, degree, poles, form, parameters, points);
            } else {
                EvalKernel<D, fixed, false>(knots, degree, poles, form, parameters, points);
            }
            done = true;
        }
    }
    if (done) {
        return;
    }
    if (rational) {
        EvalKernel<D, 0, true>(knots, degree, poles, form, parameters, points);
    } else {
        EvalKernel<D, 0, false>(knots, degree, poles, form, parameters, points);
    }
}

template <std::size_t D>
void EvalPrepared(const BasicSpanPolynomials<D>& polynomials,
                  const std::span<const double> parameters,
                  const std::span<BasicPointN<double, D>> points)
{
    const std::size_t stride = static_cast<std::size_t>(polynomials.Degree + 1) * (D + 1);
    const auto form = [&polynomials, stride](const std::size_t span) {
        const std::size_t slot = polynomials.Slots[span];
        return SpanForm{polynomials.Midpoints[slot], polynomials.InverseWidths[slot],
                        polynomials.Coefficients.data() + slot * stride};
    };
    // A periodic curve's parameters outside the domain or on its end, carried
    // in as Detail::IntoDomain carries them for the curve itself: by whole
    // periods, and the end b to the start a.
    const auto p = static_cast<std::size_t>(polynomials.Degree);
    const Interval domain{polynomials.Knots[p], polynomials.Knots[polynomials.PoleCount]};
    const auto offSeam = [&](const double t) {
        return !Contains(domain, t) || t == domain.GetUpper();
    };
    std::vector<double> storage;
    std::span<const double> inside = parameters;
    if (polynomials.Period > 0.0 && std::ranges::any_of(parameters, offSeam)) {
        storage.reserve(parameters.size());
        for (const double t : parameters) {
            double carried = t;
            if (std::isfinite(t) && offSeam(t)) {
                carried = NearestInRange(t, domain, polynomials.Period);
                carried = carried == domain.GetUpper() ? domain.GetLower() : carried;
            }
            storage.push_back(carried);
        }
        inside = storage;
    }
    EvalAllFastN<D>(polynomials.Knots, polynomials.Degree, polynomials.PoleCount,
                    polynomials.Rational, form, inside, points);
}

template <std::size_t D>
void EvalConverting(const BasicBSplineCurve<double, D>& curve,
                    const std::span<const double> parameters,
                    const std::span<BasicPointN<double, D>> points)
{
    std::array<double, MaxBasis*(D + 1)> buffer; // written by Convert before use
    std::size_t converted = curve.Poles.size();  // no span has this index
    SpanForm current;
    const auto form = [&](const std::size_t span) {
        if (span != converted) {
            current = Convert(curve, span, buffer.data());
            converted = span;
        }
        return current;
    };
    EvalAllFastN<D>(curve.Knots, curve.Degree, curve.Poles.size(), IsRational(curve), form,
                    parameters, points);
}

} // namespace

SpanPolynomials ToSpanPolynomials(const BSplineCurve& curve)
{
    SpanPolynomials polynomials = SpanPolynomialsOf(curve);
    polynomials.Period = Period(curve);
    return polynomials;
}

SpanPolynomials2d ToSpanPolynomials(const BSplineCurve2d& curve)
{
    SpanPolynomials2d polynomials = SpanPolynomialsOf(curve);
    polynomials.Period = Period(curve);
    return polynomials;
}

void EvalAllFast(const SpanPolynomials& polynomials, const std::span<const double> parameters,
                 const std::span<Point3> points)
{
    EvalPrepared(polynomials, parameters, points);
}

void EvalAllFast(const SpanPolynomials2d& polynomials, const std::span<const double> parameters,
                 const std::span<Point2> points)
{
    EvalPrepared(polynomials, parameters, points);
}

void EvalAllFast(const BSplineCurve& curve, const std::span<const double> parameters,
                 const std::span<Point3> points)
{
    std::vector<double> storage;
    EvalConverting(curve, Detail::IntoDomain(curve, parameters, storage), points);
}

void EvalAllFast(const BSplineCurve2d& curve, const std::span<const double> parameters,
                 const std::span<Point2> points)
{
    std::vector<double> storage;
    EvalConverting(curve, Detail::IntoDomain(curve, parameters, storage), points);
}

} // namespace Stereon
