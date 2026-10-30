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

#include <array>
#include <cstddef>
#include <span>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;
import :BatchEvaluation;

/// The lanes of a vector evaluate consecutive parameters of one knot span.
/// Everything that depends on the span alone (the knots of the basis
/// recurrence, the poles and weights) is the same in every lane and is
/// broadcast; the parameter, and so the basis functions, differ by lane.
/// Each lane performs exactly the operations of BasisFunctions and Eval, in
/// their order: IEEE arithmetic is the same in a lane as in a scalar, and the
/// build keeps contraction off (-ffp-contract=off), so no lane fuses what the
/// scalar code rounds twice.

namespace Stereon {

namespace {

using Lanes = Simd::Vec<double>;
constexpr std::size_t LaneCount = Simd::Width<double>;
constexpr std::size_t MaxBasis = static_cast<std::size_t>(BSplineMaxDegree) + 1;

/// The span FindSpan gives t, tried first against `guess`: a parameter with
/// knots[guess] <= t < knots[guess + 1] in a non-empty span of the domain is
/// in that span, by FindSpan's definition, and the binary search is left for
/// the others (a new span, or the end of the domain, which FindSpan gives
/// the last non-empty span).
template <std::size_t D>
[[nodiscard]] std::size_t SpanOf(const BasicBSplineCurve<double, D>& curve, const double t,
                                 const std::size_t guess)
{
    const std::vector<double>& u = curve.Knots;
    if (u[guess] <= t && t < u[guess + 1]) {
        return guess;
    }
    return FindSpan(u, curve.Degree, curve.Poles.size(), t);
}

/// The span of the run of LaneCount parameters at `t`, or `none` when they
/// do not share one; `guess` is updated to the first lane's span.
template <std::size_t D>
[[nodiscard]] std::size_t SharedSpan(const BasicBSplineCurve<double, D>& curve, const double* t,
                                     std::size_t& guess, const std::size_t none)
{
    guess = SpanOf(curve, t[0], guess);
    for (std::size_t lane = 1; lane < LaneCount; ++lane) {
        if (SpanOf(curve, t[lane], guess) != guess) {
            return none;
        }
    }
    return guess;
}

/// BasisFunctions (A2.2) in every lane: values[j] = N_{span-p+j,p}(t).
void LaneBasis(const std::span<const double> knots, const std::size_t p, const std::size_t span,
               const Lanes& t, std::array<Lanes, MaxBasis>& values)
{
    std::array<Lanes, MaxBasis> left;  // written before read, as in BasisFunctions
    std::array<Lanes, MaxBasis> right; // likewise
    values[0] = Simd::Broadcast(1.0);
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = t - Simd::Broadcast(knots[span + 1 - j]);
        right[j] = Simd::Broadcast(knots[span + j]) - t;
        Lanes saved = Simd::Broadcast(0.0);
        for (std::size_t r = 0; r < j; ++r) {
            const Lanes temp = values[r] / (right[r + 1] + left[j - r]);
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
    }
}

/// Eval's combination of the poles in every lane, into points[0, LaneCount).
template <std::size_t D>
void EvalLanes(const BasicBSplineCurve<double, D>& curve, const std::size_t span, const double* t,
               BasicPointN<double, D>* points)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    std::array<Lanes, MaxBasis> basis; // filled by LaneBasis
    LaneBasis(curve.Knots, p, span, Simd::Load(t), basis);
    std::array<Lanes, D> numerator;
    numerator.fill(Simd::Broadcast(0.0));
    Lanes denominator = Simd::Broadcast(0.0);
    for (std::size_t j = 0; j <= p; ++j) {
        const std::size_t i = span - p + j;
        const Lanes weighted = basis[j] * Simd::Broadcast(Detail::WeightOf(curve, i));
        for (std::size_t d = 0; d < D; ++d) {
            numerator[d] = numerator[d] + weighted * Simd::Broadcast(curve.Poles[i][d]);
        }
        denominator = denominator + weighted;
    }
    std::array<std::array<double, LaneCount>, D> coordinates; // filled by Store
    for (std::size_t d = 0; d < D; ++d) {
        const Lanes value = IsRational(curve) ? numerator[d] / denominator : numerator[d];
        Simd::Store(value, coordinates[d].data());
    }
    for (std::size_t lane = 0; lane < LaneCount; ++lane) {
        for (std::size_t d = 0; d < D; ++d) {
            points[lane][d] = coordinates[d][lane];
        }
    }
}

template <std::size_t D>
void EvalAllN(const BasicBSplineCurve<double, D>& curve, const std::span<const double> given,
              const std::span<BasicPointN<double, D>> points)
{
    std::vector<double> storage;
    const std::span<const double> parameters = Detail::IntoDomain(curve, given, storage);
    Detail::Require(points.size() == parameters.size(), "EvalAll(BSplineCurve)",
                    "as many points as parameters");
    for (const double t : parameters) {
        Detail::CheckEvaluation(curve, t, "EvalAll(BSplineCurve)");
    }
    const std::size_t none = curve.Poles.size();                // no span has this index
    std::size_t guess = static_cast<std::size_t>(curve.Degree); // the first span's index or below
    std::size_t k = 0;
    for (; k + LaneCount <= parameters.size(); k += LaneCount) {
        const std::size_t span = SharedSpan(curve, parameters.data() + k, guess, none);
        if (span != none) {
            EvalLanes(curve, span, parameters.data() + k, points.data() + k);
            continue;
        }
        for (std::size_t lane = 0; lane < LaneCount; ++lane) {
            points[k + lane] = Eval(curve, parameters[k + lane]);
        }
    }
    for (; k < parameters.size(); ++k) {
        points[k] = Eval(curve, parameters[k]);
    }
}

} // namespace

void EvalAll(const BSplineCurve& curve, const std::span<const double> parameters,
             const std::span<Point3> points)
{
    EvalAllN(curve, parameters, points);
}

void EvalAll(const BSplineCurve2d& curve, const std::span<const double> parameters,
             const std::span<Point2> points)
{
    EvalAllN(curve, parameters, points);
}

} // namespace Stereon
