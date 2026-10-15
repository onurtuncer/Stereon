// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;
import :Fitting;

namespace Stereon {

// Instantiates Vec3 before any arithmetic on it: GCC 16.2 crashes when the
// first instantiation of BasicVec3<double> in this unit is reached through
// its hidden-friend operators inside a function
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

// ---------------------------------------------------------------------------
// A banded linear solver

/// An n x n matrix with non-zeros only within `Half` of the diagonal, stored
/// row by row in 2 Half + 1 columns.
class BandMatrix {
public:
    BandMatrix(const std::size_t n, const std::size_t half)
        : m_N(n)
        , m_Half(half)
        , m_Entries(n * (2 * half + 1), 0.0)
    {}

    [[nodiscard]] double& At(const std::size_t i, const std::size_t j)
    {
        return m_Entries[i * (2 * m_Half + 1) + (j + m_Half - i)];
    }
    [[nodiscard]] std::size_t Size() const { return m_N; }
    [[nodiscard]] std::size_t Half() const { return m_Half; }

private:
    std::size_t m_N;
    std::size_t m_Half;
    std::vector<double> m_Entries;
};

/// Solves M x = b in place by Gaussian elimination without pivoting, for
/// three right-hand sides at once (one per coordinate). Without pivoting the
/// band structure is kept; it is stable for the matrices used here, B-spline
/// collocation matrices (totally positive) and normal equations (symmetric
/// positive definite).
[[nodiscard]] std::vector<Vec3> SolveBanded(BandMatrix& m, std::vector<Vec3> b)
{
    const std::size_t n = m.Size();
    const std::size_t w = m.Half();
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = k + 1; i < n && i <= k + w; ++i) {
            const double factor = m.At(i, k) / m.At(k, k);
            for (std::size_t j = k + 1; j < n && j <= k + w; ++j) {
                m.At(i, j) -= factor * m.At(k, j);
            }
            b[i] -= factor * b[k];
        }
    }
    for (std::size_t r = n; r > 0; --r) {
        const std::size_t i = r - 1;
        Vec3 sum = b[i];
        for (std::size_t j = i + 1; j < n && j <= i + w; ++j) {
            sum -= m.At(i, j) * b[j];
        }
        b[i] = sum / m.At(i, i);
    }
    return b;
}

// ---------------------------------------------------------------------------
// Validation and dimension handling

[[nodiscard]] Result<void> CheckPoints(const std::span<const Point3> points, const int degree,
                                       const std::string_view where)
{
    if (degree < 1 || degree > BSplineMaxDegree) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} + ": the degree must be between 1 and 25");
    }
    if (points.size() < static_cast<std::size_t>(degree) + 1) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} + ": fitting needs at least degree + 1 points");
    }
    for (const Point3& p : points) {
        STEREON_TRY_VOID(Detail::CheckFinite(p, "every point", where));
    }
    return {};
}

[[nodiscard]] Result<std::vector<double>> CheckedParameters(const std::span<const Point3> points,
                                                            const Parametrization parametrization,
                                                            const std::string_view where)
{
    std::vector<double> parameters = Detail::ParametersOf(points, parametrization);
    if (parameters.empty()) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} +
                        ": consecutive points coincide, so chord-length parameters are undefined");
    }
    return parameters;
}

[[nodiscard]] std::vector<Point3> Lifted(const std::span<const Point2> points)
{
    std::vector<Point3> lifted;
    lifted.reserve(points.size());
    for (const Point2& p : points) {
        lifted.push_back({p.X, p.Y, 0.0});
    }
    return lifted;
}

[[nodiscard]] CurveFit2d Flattened(CurveFit fit)
{
    CurveFit2d flat{{fit.Curve.Degree, std::move(fit.Curve.Knots), {}, {}},
                    std::move(fit.Parameters),
                    fit.MaxDeviation};
    for (const Point3& p : fit.Curve.Poles) {
        flat.Curve.Poles.push_back({p.X, p.Y});
    }
    return flat;
}

template <typename F>
[[nodiscard]] Result<CurveFit2d> InPlane(const std::span<const Point2> points, F&& fit)
{
    const std::vector<Point3> lifted = Lifted(points);
    STEREON_TRY(CurveFit fitted, fit(std::span<const Point3>{lifted}));
    return Flattened(std::move(fitted));
}

/// The fit's largest deviation at its points.
[[nodiscard]] double MaxDeviationOf(const BSplineCurve& curve, const std::span<const Point3> points,
                                    const std::span<const double> parameters)
{
    double worst = 0.0;
    for (std::size_t k = 0; k < points.size(); ++k) {
        worst = std::max(worst, Distance(Eval(curve, parameters[k]), points[k]));
    }
    return worst;
}

// ---------------------------------------------------------------------------
// Least squares

/// Knots for a least-squares fit with poleCount poles (eq. 9.69): with
/// d = (m + 1) / (n - p + 1), interior knot p + j blends the parameters
/// around j d, so that every knot span holds at least one parameter and the
/// normal equations are positive definite.
[[nodiscard]] std::vector<double> SpreadKnots(const std::span<const double> t, const int degree,
                                              const std::size_t poleCount)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    std::vector<double> knots(p + 1, 0.0);
    const double d = static_cast<double>(t.size()) / static_cast<double>(n - p + 1);
    for (std::size_t j = 1; j + p <= n; ++j) {
        const double position = static_cast<double>(j) * d;
        const auto i = static_cast<std::size_t>(position);
        const double alpha = position - static_cast<double>(i);
        knots.push_back((1.0 - alpha) * t[i - 1] + alpha * t[i]);
    }
    knots.insert(knots.end(), p + 1, 1.0);
    return knots;
}

/// The interior poles of the least-squares fit: the rows of N are the
/// basis functions at the interior parameters, restricted to the interior
/// poles, and R_k is Q_k minus the contributions of the fixed end poles.
[[nodiscard]] std::vector<Point3> LeastSquaresPoles(const std::span<const Point3> q,
                                                    const std::span<const double> t,
                                                    const std::span<const double> knots,
                                                    const int degree, const std::size_t poleCount)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    std::vector<Point3> poles(poleCount);
    poles.front() = q.front();
    poles.back() = q.back();
    if (n < 2) {
        return poles; // no interior pole
    }
    BandMatrix normal(n - 1, p);
    std::vector<Vec3> rhs(n - 1);
    std::array<double, BSplineMaxDegree + 1> basis; // filled by BasisFunctions
    for (std::size_t k = 1; k + 1 < q.size(); ++k) {
        const std::size_t span = FindSpan(knots, degree, poleCount, t[k]);
        BasisFunctions(knots, degree, span, t[k], std::span{basis}.first(p + 1));
        Vec3 r = ToVec(q[k]);
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = span - p + a;
            r -= (i == 0 ? basis[a] * ToVec(q.front()) : Vec3{}) +
                 (i == n ? basis[a] * ToVec(q.back()) : Vec3{});
        }
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = span - p + a;
            if (i == 0 || i == n) {
                continue;
            }
            rhs[i - 1] += basis[a] * r;
            for (std::size_t b = 0; b <= p; ++b) {
                const std::size_t j = span - p + b;
                if (j != 0 && j != n) {
                    normal.At(i - 1, j - 1) += basis[a] * basis[b];
                }
            }
        }
    }
    const std::vector<Vec3> interior = SolveBanded(normal, std::move(rhs));
    for (std::size_t i = 1; i < n; ++i) {
        poles[i] = ToPoint(interior[i - 1]);
    }
    return poles;
}

[[nodiscard]] Result<CurveFit> FitCount(const std::span<const Point3> points, const int degree,
                                        const std::size_t poleCount,
                                        const Parametrization parametrization);

/// The fit with `count` poles: interpolation when there are as many poles
/// as points, least squares otherwise.
[[nodiscard]] Result<CurveFit> FitWith(const std::span<const Point3> points, const int degree,
                                       const std::size_t count,
                                       const Parametrization parametrization)
{
    return count == points.size() ? Interpolate(points, degree, parametrization)
                                  : FitCount(points, degree, count, parametrization);
}

Result<CurveFit> FitCount(const std::span<const Point3> points, const int degree,
                          const std::size_t poleCount, const Parametrization parametrization)
{
    constexpr std::string_view where = "FitLeastSquares";
    STEREON_TRY_VOID(CheckPoints(points, degree, where));
    if (poleCount < static_cast<std::size_t>(degree) + 1 || poleCount > points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "FitLeastSquares: the pole count must be between degree + 1 and the number "
                    "of points");
    }
    STEREON_TRY(std::vector<double> t, CheckedParameters(points, parametrization, where));
    std::vector<double> knots = SpreadKnots(t, degree, poleCount);
    std::vector<Point3> poles = LeastSquaresPoles(points, t, knots, degree, poleCount);
    BSplineCurve curve{degree, std::move(knots), std::move(poles), {}};
    const double deviation = MaxDeviationOf(curve, points, t);
    return CurveFit{std::move(curve), std::move(t), deviation};
}

/// The smallest pole count among those tried whose fit is within the
/// budget: doubling from degree + 1 until one passes, then bisection
/// between the last failure and the first success.
[[nodiscard]] Result<CurveFit> ApproximateWithin(const std::span<const Point3> points,
                                                 const int degree, const double budget,
                                                 const Parametrization parametrization)
{
    std::size_t failing = 0;
    std::size_t count = static_cast<std::size_t>(degree) + 1;
    double best = 0.0;
    while (true) {
        STEREON_TRY(CurveFit fit, FitWith(points, degree, count, parametrization));
        if (fit.MaxDeviation <= budget) {
            break;
        }
        best = fit.MaxDeviation;
        failing = count;
        if (count == points.size()) {
            return std::unexpected(
                KernelError(ErrorCode::ToleranceExceeded,
                            "Approximate: even the interpolant exceeds the budget")
                    .WithTolerance(budget)
                    .WithTolerance(best));
        }
        count = std::min(2 * count, points.size());
    }
    while (count - failing > 1 && failing != 0) {
        const std::size_t mid = failing + (count - failing) / 2;
        STEREON_TRY(CurveFit fit, FitWith(points, degree, mid, parametrization));
        if (fit.MaxDeviation <= budget) {
            count = mid;
        } else {
            failing = mid;
        }
    }
    return FitWith(points, degree, count, parametrization);
}

} // namespace

Result<CurveFit> Interpolate(const std::span<const Point3> points, const int degree,
                             const Parametrization parametrization)
{
    constexpr std::string_view where = "Interpolate";
    STEREON_TRY_VOID(CheckPoints(points, degree, where));
    STEREON_TRY(std::vector<double> t, CheckedParameters(points, parametrization, where));
    std::vector<double> knots = Detail::AveragedKnots(t, degree);
    std::vector<Point3> poles = Detail::InterpolatingPoles(points, t, knots, degree);
    BSplineCurve curve{degree, std::move(knots), std::move(poles), {}};
    const double deviation = MaxDeviationOf(curve, points, t);
    return CurveFit{std::move(curve), std::move(t), deviation};
}

Result<CurveFit2d> Interpolate(const std::span<const Point2> points, const int degree,
                               const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Interpolate(lifted, degree, parametrization);
    });
}

Result<CurveFit> FitLeastSquares(const std::span<const Point3> points, const int degree,
                                 const std::size_t poleCount, const Parametrization parametrization)
{
    return FitCount(points, degree, poleCount, parametrization);
}

Result<CurveFit2d> FitLeastSquares(const std::span<const Point2> points, const int degree,
                                   const std::size_t poleCount,
                                   const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return FitCount(lifted, degree, poleCount, parametrization);
    });
}

Result<CurveFit> Approximate(const std::span<const Point3> points, const int degree,
                             const ApproximationBudget budget,
                             const Parametrization parametrization)
{
    if (!std::isfinite(budget.Deviation) || !(budget.Deviation > 0.0)) {
        return Fail(ErrorCode::InvalidInput, "Approximate: the budget must be finite and positive");
    }
    STEREON_TRY_VOID(CheckPoints(points, degree, "Approximate"));
    return ApproximateWithin(points, degree, budget.Deviation, parametrization);
}

Result<CurveFit2d> Approximate(const std::span<const Point2> points, const int degree,
                               const ApproximationBudget budget,
                               const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Approximate(lifted, degree, budget, parametrization);
    });
}

} // namespace Stereon

namespace Stereon::Detail {

std::vector<double> ParametersOf(const std::span<const Point3> points,
                                 const Parametrization parametrization)
{
    const std::size_t n = points.size() - 1;
    std::vector<double> steps(n + 1, 0.0);
    for (std::size_t k = 1; k <= n; ++k) {
        const double chord = Distance(points[k], points[k - 1]);
        if (parametrization != Parametrization::Uniform && chord == 0.0) {
            return {};
        }
        steps[k] = steps[k - 1] + (parametrization == Parametrization::Uniform       ? 1.0
                                   : parametrization == Parametrization::Centripetal ? Sqrt(chord)
                                                                                     : chord);
    }
    std::vector<double> t(n + 1);
    for (std::size_t k = 0; k <= n; ++k) {
        t[k] = steps[k] / steps[n];
    }
    t.back() = 1.0; // exactly
    return t;
}

std::vector<double> AveragedKnots(const std::span<const double> parameters, const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = parameters.size() - 1;
    std::vector<double> knots(p + 1, 0.0);
    for (std::size_t j = 1; j + p <= n; ++j) {
        double sum = 0.0;
        for (std::size_t i = j; i < j + p; ++i) {
            sum += parameters[i];
        }
        knots.push_back(sum / static_cast<double>(p));
    }
    knots.insert(knots.end(), p + 1, 1.0);
    return knots;
}

std::vector<Point3> InterpolatingPoles(const std::span<const Point3> points,
                                       const std::span<const double> parameters,
                                       const std::span<const double> knots, const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t count = points.size();
    BandMatrix matrix(count, p);
    std::vector<Vec3> rhs;
    rhs.reserve(count);
    std::array<double, BSplineMaxDegree + 1> basis; // filled by BasisFunctions
    for (std::size_t k = 0; k < count; ++k) {
        const std::size_t span = FindSpan(knots, degree, count, parameters[k]);
        BasisFunctions(knots, degree, span, parameters[k], std::span{basis}.first(p + 1));
        // Averaged knots satisfy the Schoenberg-Whitney conditions, which
        // puts every non-zero within p of the diagonal.
        Require(span >= k && span <= k + p, "Interpolate", "a collocation row outside the band");
        for (std::size_t a = 0; a <= p; ++a) {
            matrix.At(k, span - p + a) = basis[a];
        }
        rhs.push_back(ToVec(points[k]));
    }
    std::vector<Point3> poles;
    for (const Vec3& v : SolveBanded(matrix, std::move(rhs))) {
        poles.push_back(ToPoint(v));
    }
    return poles;
}

} // namespace Stereon::Detail
