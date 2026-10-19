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
#include <functional>
#include <optional>
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
// Collocation rows

/// The knot span of t for parameters that come in increasing order: scanned
/// forward from the previous row's span rather than searched afresh, so a
/// sweep over n sorted parameters costs O(n) rather than O(n log n). The
/// ends keep FindSpan's rules for the empty spans there, and checked builds
/// verify the agreement with FindSpan on every row.
[[nodiscard]] std::size_t SpanAfter(const std::span<const double> knots, const int degree,
                                    const std::size_t poleCount, const double t,
                                    const std::size_t previous)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    if (t <= knots[p] || t >= knots[n + 1]) {
        return FindSpan(knots, degree, poleCount, t);
    }
    std::size_t span = std::max(previous, p);
    while (knots[span + 1] <= t) {
        ++span;
    }
    Detail::Require(span == FindSpan(knots, degree, poleCount, t), "Interpolate",
                    "a scanned span that differs from the searched one");
    return span;
}

/// The collocation rows of increasing parameters: for row k the span of
/// t_k and the p + 1 basis functions there, Basis[k * (p + 1) + a] =
/// N_{span-p+a,p}(t_k). Built once per fit; the system is assembled from
/// them and the deviation measured with them.
struct CollocationRows {
    std::vector<std::size_t> Spans;
    std::vector<double> Basis;
};

[[nodiscard]] CollocationRows CollocationRowsOf(const std::span<const double> parameters,
                                                const std::span<const double> knots,
                                                const int degree, const std::size_t poleCount)
{
    const auto p = static_cast<std::size_t>(degree);
    CollocationRows rows;
    rows.Spans.reserve(parameters.size());
    rows.Basis.resize(parameters.size() * (p + 1));
    std::size_t span = p;
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        span = SpanAfter(knots, degree, poleCount, parameters[k], span);
        rows.Spans.push_back(span);
        BasisFunctions(knots, degree, span, parameters[k],
                       std::span{rows.Basis}.subspan(k * (p + 1), p + 1));
    }
    return rows;
}

/// The largest distance between the points and the curve with these poles
/// at the rows' parameters: what MaxDeviationOf measures, to the bit, since
/// Eval of a polynomial curve sums the same basis functions and poles in
/// the same order; found here without a second span search and basis
/// evaluation per point.
[[nodiscard]] double RowDeviation(const CollocationRows& rows, const std::span<const Point3> poles,
                                  const std::span<const Point3> points, const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    double worst = 0.0;
    for (std::size_t k = 0; k < points.size(); ++k) {
        Vec3 sum{};
        for (std::size_t a = 0; a <= p; ++a) {
            sum += rows.Basis[k * (p + 1) + a] * ToVec(poles[rows.Spans[k] - p + a]);
        }
        worst = std::max(worst, Distance(ToPoint(sum), points[k]));
    }
    return worst;
}

/// The poles of the unconstrained interpolant from its rows: the banded
/// collocation system, solved (see Detail::InterpolatingPoles).
[[nodiscard]] std::vector<Point3> InterpolatingPolesFrom(const std::span<const Point3> points,
                                                         const CollocationRows& rows,
                                                         const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t count = points.size();
    BandMatrix matrix(count, p);
    std::vector<Vec3> rhs;
    rhs.reserve(count);
    for (std::size_t k = 0; k < count; ++k) {
        const std::size_t span = rows.Spans[k];
        // Averaged knots satisfy the Schoenberg-Whitney conditions, which
        // puts every non-zero within p of the diagonal.
        Detail::Require(span >= k && span <= k + p, "Interpolate",
                        "a collocation row outside the band");
        for (std::size_t a = 0; a <= p; ++a) {
            matrix.At(k, span - p + a) = rows.Basis[k * (p + 1) + a];
        }
        rhs.push_back(ToVec(points[k]));
    }
    std::vector<Point3> poles;
    poles.reserve(count);
    for (const Vec3& v : SolveBanded(matrix, std::move(rhs))) {
        poles.push_back(ToPoint(v));
    }
    return poles;
}

// ---------------------------------------------------------------------------
// Validation and dimension handling

/// End derivatives when the caller prescribed them.
using Ends = std::optional<EndDerivatives>;

/// The lowest degree that can take end derivatives: with degree 1 the
/// shifted knots of the constrained interpolant would repeat the end knots,
/// u_{p+1} = t_0 = 0, and the fixed second pole would divide by zero.
constexpr int LowestConstrainedDegree = 2;

/// With end derivatives two poles are fixed at each end, so the fewest
/// poles any constrained fit can have.
constexpr std::size_t FewestConstrainedPoles = 4;

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

[[nodiscard]] Result<void> CheckEnds(const Ends& ends, const int degree,
                                     const std::string_view where)
{
    if (!ends) {
        return {};
    }
    if (degree < LowestConstrainedDegree) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} + ": end derivatives need degree 2 or more");
    }
    STEREON_TRY_VOID(Detail::CheckFinite(ends->Start, "the start derivative", where));
    STEREON_TRY_VOID(Detail::CheckFinite(ends->End, "the end derivative", where));
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

[[nodiscard]] EndDerivatives Lifted(const EndDerivatives2d& ends)
{
    return {{ends.Start.X, ends.Start.Y, 0.0}, {ends.End.X, ends.End.Y, 0.0}};
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
// Fixed poles

/// The poles that the data and the end derivatives fix before any system is
/// solved, in a vector of poleCount poles whose other entries are still to
/// be computed: the first and last pole are the first and last point, and
/// with end derivatives the second and the second to last follow from
/// C'(0) = p (P_1 - P_0) / u_{p+1} and C'(1) = p (P_n - P_{n-1}) / (1 - u_n)
/// (eq. 9.10). Returns the poles and how many are fixed at each end.
[[nodiscard]] std::pair<std::vector<Point3>, std::size_t>
FixedPoles(const std::span<const Point3> q, const std::span<const double> knots, const int degree,
           const std::size_t poleCount, const Ends& ends)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    // Filled with the first point rather than assigned through front() and
    // back(): GCC 16 at -O2 cannot prove the vector non-empty and reports
    // a potential null dereference (-Wnull-dereference). Every caller
    // passes at least two poles.
    std::vector<Point3> poles(poleCount, q.front());
    poles[n] = q.back();
    if (!ends) {
        return {std::move(poles), 1};
    }
    const double pd = static_cast<double>(degree);
    poles[1] = q.front() + (knots[p + 1] / pd) * ends->Start;
    poles[n - 1] = q.back() - ((1.0 - knots[n]) / pd) * ends->End;
    return {std::move(poles), 2};
}

// ---------------------------------------------------------------------------
// Interpolation with end derivatives

/// Knots for the interpolant with end derivatives (eq. 9.11, for any
/// degree): two more poles than points, so two more interior knots than
/// AveragedKnots gives, each the mean of p parameters starting one earlier,
/// interior knot j + p averaging parameters j - 1, ..., j + p - 2. The
/// Schoenberg-Whitney conditions then hold for point k against basis
/// function k + 1, the position it has once the fixed second pole is
/// counted.
[[nodiscard]] std::vector<double> ShiftedAveragedKnots(const std::span<const double> parameters,
                                                       const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t m = parameters.size() - 1;
    std::vector<double> knots(p + 1, 0.0);
    for (std::size_t j = 1; j + p <= m + 2; ++j) {
        double sum = 0.0;
        for (std::size_t i = j - 1; i + 2 <= j + p; ++i) {
            sum += parameters[i];
        }
        knots.push_back(sum / static_cast<double>(p));
    }
    knots.insert(knots.end(), p + 1, 1.0);
    return knots;
}

/// The poles of the interpolant with end derivatives: with the four end
/// poles fixed, the interior points' collocation rows restricted to the
/// free poles P_2, ..., P_{n-2} form a square banded system, still totally
/// positive, whose right-hand side is each point less the fixed poles'
/// contributions.
[[nodiscard]] std::vector<Point3> ConstrainedInterpolatingPoles(const std::span<const Point3> q,
                                                                const CollocationRows& rows,
                                                                const std::span<const double> knots,
                                                                const int degree,
                                                                const EndDerivatives& ends)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t poleCount = q.size() + 2;
    const std::size_t n = poleCount - 1;
    auto [poles, fixed] = FixedPoles(q, knots, degree, poleCount, ends);
    const std::size_t free = poleCount - 2 * fixed;
    BandMatrix matrix(free, p);
    std::vector<Vec3> rhs(free);
    for (std::size_t k = 1; k + 1 < q.size(); ++k) {
        const std::size_t span = rows.Spans[k];
        const std::span<const double> basis = std::span{rows.Basis}.subspan(k * (p + 1), p + 1);
        // Point k is row k - 1 and pole i column i - 2; the shifted knots put
        // basis function k + 1 on the diagonal, within p of every non-zero.
        Detail::Require(span >= k + 1 && span <= k + 1 + p, "Interpolate",
                        "a collocation row outside the band");
        Vec3 r = ToVec(q[k]);
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = span - p + a;
            if (i < fixed || i > n - fixed) {
                r -= basis[a] * ToVec(poles[i]);
            } else {
                matrix.At(k - 1, i - fixed) = basis[a];
            }
        }
        rhs[k - 1] = r;
    }
    const std::vector<Vec3> interior = SolveBanded(matrix, std::move(rhs));
    for (std::size_t i = 0; i < free; ++i) {
        poles[i + fixed] = ToPoint(interior[i]);
    }
    return poles;
}

/// Interpolation with or without end derivatives.
[[nodiscard]] Result<CurveFit> InterpolateWith(const std::span<const Point3> points,
                                               const int degree, const Ends& ends,
                                               const Parametrization parametrization)
{
    constexpr std::string_view where = "Interpolate";
    STEREON_TRY_VOID(CheckPoints(points, degree, where));
    STEREON_TRY_VOID(CheckEnds(ends, degree, where));
    STEREON_TRY(std::vector<double> t, CheckedParameters(points, parametrization, where));
    std::vector<double> knots =
        ends ? ShiftedAveragedKnots(t, degree) : Detail::AveragedKnots(t, degree);
    const std::size_t poleCount = ends ? points.size() + 2 : points.size();
    const CollocationRows rows = CollocationRowsOf(t, knots, degree, poleCount);
    std::vector<Point3> poles =
        ends ? ConstrainedInterpolatingPoles(points, rows, knots, degree, *ends)
             : InterpolatingPolesFrom(points, rows, degree);
    const double deviation = RowDeviation(rows, poles, points, degree);
    BSplineCurve curve{degree, std::move(knots), std::move(poles), {}};
    return CurveFit{std::move(curve), std::move(t), deviation};
}

// ---------------------------------------------------------------------------
// Least squares

/// Solves the normal equations for the free poles of a least-squares fit
/// whose `fixed` first and last poles are already in `poles`: the rows of
/// N are the basis functions at the interior parameters restricted to the
/// free poles, and R_k is Q_k less the contributions of the fixed poles
/// (algorithm A9.7, with two fixed poles per end when end derivatives are
/// prescribed).
void SolveFreePoles(const std::span<const Point3> q, const std::span<const double> t,
                    const std::span<const double> knots, const int degree,
                    std::vector<Point3>& poles, const std::size_t fixed)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poles.size() - 1;
    if (poles.size() <= 2 * fixed) {
        return; // no free pole
    }
    const std::size_t free = poles.size() - 2 * fixed;
    const auto isFixed = [&](const std::size_t i) { return i < fixed || i > n - fixed; };
    BandMatrix normal(free, p);
    std::vector<Vec3> rhs(free);
    std::array<double, BSplineMaxDegree + 1> basis; // filled by BasisFunctions
    for (std::size_t k = 1; k + 1 < q.size(); ++k) {
        const std::size_t span = FindSpan(knots, degree, poles.size(), t[k]);
        BasisFunctions(knots, degree, span, t[k], std::span{basis}.first(p + 1));
        Vec3 r = ToVec(q[k]);
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = span - p + a;
            if (isFixed(i)) {
                r -= basis[a] * ToVec(poles[i]);
            }
        }
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = span - p + a;
            if (isFixed(i)) {
                continue;
            }
            rhs[i - fixed] += basis[a] * r;
            for (std::size_t b = 0; b <= p; ++b) {
                const std::size_t j = span - p + b;
                if (!isFixed(j)) {
                    normal.At(i - fixed, j - fixed) += basis[a] * basis[b];
                }
            }
        }
    }
    const std::vector<Vec3> interior = SolveBanded(normal, std::move(rhs));
    for (std::size_t i = 0; i < free; ++i) {
        poles[i + fixed] = ToPoint(interior[i]);
    }
}

/// The least-squares fit with `poleCount` poles, with or without end
/// derivatives.
[[nodiscard]] Result<CurveFit> FitCount(const std::span<const Point3> points, const int degree,
                                        const std::size_t poleCount, const Ends& ends,
                                        const Parametrization parametrization)
{
    constexpr std::string_view where = "FitLeastSquares";
    STEREON_TRY_VOID(CheckPoints(points, degree, where));
    STEREON_TRY_VOID(CheckEnds(ends, degree, where));
    const std::size_t fewest =
        std::max(static_cast<std::size_t>(degree) + 1, ends ? FewestConstrainedPoles : 0);
    if (poleCount < fewest || poleCount > points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "FitLeastSquares: the pole count must be between degree + 1 (4 with end "
                    "derivatives) and the number of points");
    }
    STEREON_TRY(std::vector<double> t, CheckedParameters(points, parametrization, where));
    std::vector<double> knots = Detail::SpreadKnots(t, degree, poleCount);
    auto [poles, fixed] = FixedPoles(points, knots, degree, poleCount, ends);
    SolveFreePoles(points, t, knots, degree, poles, fixed);
    BSplineCurve curve{degree, std::move(knots), std::move(poles), {}};
    const double deviation = MaxDeviationOf(curve, points, t);
    return CurveFit{std::move(curve), std::move(t), deviation};
}

// ---------------------------------------------------------------------------
// Approximation to a budget

/// The fits Approximate searches, by level: least squares with `level`
/// poles up to as many as points, and above that the interpolant, which has
/// as many poles as points, or two more with end derivatives.
[[nodiscard]] Result<CurveFit> FitAtLevel(const std::span<const Point3> points, const int degree,
                                          const std::size_t level, const Ends& ends,
                                          const Parametrization parametrization)
{
    const std::size_t interpolant = points.size() + (ends ? 1 : 0);
    return level >= interpolant ? InterpolateWith(points, degree, ends, parametrization)
                                : FitCount(points, degree, level, ends, parametrization);
}

[[nodiscard]] Result<CurveFit> ApproximateWith(const std::span<const Point3> points,
                                               const int degree, const ApproximationBudget budget,
                                               const Ends& ends,
                                               const Parametrization parametrization)
{
    constexpr std::string_view where = "Approximate";
    if (!std::isfinite(budget.Deviation) || !(budget.Deviation > 0.0)) {
        return Fail(ErrorCode::InvalidInput, "Approximate: the budget must be finite and positive");
    }
    STEREON_TRY_VOID(CheckPoints(points, degree, where));
    STEREON_TRY_VOID(CheckEnds(ends, degree, where));
    const std::size_t fewest =
        std::max(static_cast<std::size_t>(degree) + 1, ends ? FewestConstrainedPoles : 0);
    const std::size_t interpolant = points.size() + (ends ? 1 : 0);
    STEREON_TRY(const std::size_t level,
                Detail::FewestPolesWithin(
                    fewest, interpolant, budget.Deviation,
                    [&](const std::size_t count) -> Result<double> {
                        STEREON_TRY(CurveFit fit,
                                    FitAtLevel(points, degree, count, ends, parametrization));
                        return fit.MaxDeviation;
                    },
                    where));
    return FitAtLevel(points, degree, level, ends, parametrization);
}

} // namespace

Result<CurveFit> Interpolate(const std::span<const Point3> points, const int degree,
                             const Parametrization parametrization)
{
    return InterpolateWith(points, degree, std::nullopt, parametrization);
}

Result<CurveFit2d> Interpolate(const std::span<const Point2> points, const int degree,
                               const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Interpolate(lifted, degree, parametrization);
    });
}

Result<CurveFit> Interpolate(const std::span<const Point3> points, const int degree,
                             const EndDerivatives ends, const Parametrization parametrization)
{
    return InterpolateWith(points, degree, ends, parametrization);
}

Result<CurveFit2d> Interpolate(const std::span<const Point2> points, const int degree,
                               const EndDerivatives2d ends, const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Interpolate(lifted, degree, Lifted(ends), parametrization);
    });
}

Result<CurveFit> FitLeastSquares(const std::span<const Point3> points, const int degree,
                                 const std::size_t poleCount, const Parametrization parametrization)
{
    return FitCount(points, degree, poleCount, std::nullopt, parametrization);
}

Result<CurveFit2d> FitLeastSquares(const std::span<const Point2> points, const int degree,
                                   const std::size_t poleCount,
                                   const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return FitLeastSquares(lifted, degree, poleCount, parametrization);
    });
}

Result<CurveFit> FitLeastSquares(const std::span<const Point3> points, const int degree,
                                 const std::size_t poleCount, const EndDerivatives ends,
                                 const Parametrization parametrization)
{
    return FitCount(points, degree, poleCount, ends, parametrization);
}

Result<CurveFit2d> FitLeastSquares(const std::span<const Point2> points, const int degree,
                                   const std::size_t poleCount, const EndDerivatives2d ends,
                                   const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return FitLeastSquares(lifted, degree, poleCount, Lifted(ends), parametrization);
    });
}

Result<CurveFit> Approximate(const std::span<const Point3> points, const int degree,
                             const ApproximationBudget budget,
                             const Parametrization parametrization)
{
    return ApproximateWith(points, degree, budget, std::nullopt, parametrization);
}

Result<CurveFit2d> Approximate(const std::span<const Point2> points, const int degree,
                               const ApproximationBudget budget,
                               const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Approximate(lifted, degree, budget, parametrization);
    });
}

Result<CurveFit> Approximate(const std::span<const Point3> points, const int degree,
                             const ApproximationBudget budget, const EndDerivatives ends,
                             const Parametrization parametrization)
{
    return ApproximateWith(points, degree, budget, ends, parametrization);
}

Result<CurveFit2d> Approximate(const std::span<const Point2> points, const int degree,
                               const ApproximationBudget budget, const EndDerivatives2d ends,
                               const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return Approximate(lifted, degree, budget, Lifted(ends), parametrization);
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
    return InterpolatingPolesFrom(
        points, CollocationRowsOf(parameters, knots, degree, points.size()), degree);
}

std::vector<double> SpreadKnots(const std::span<const double> parameters, const int degree,
                                const std::size_t poleCount)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = poleCount - 1;
    std::vector<double> knots(p + 1, 0.0);
    const double d = static_cast<double>(parameters.size()) / static_cast<double>(n - p + 1);
    for (std::size_t j = 1; j + p <= n; ++j) {
        const double position = static_cast<double>(j) * d;
        const auto i = static_cast<std::size_t>(position);
        const double alpha = position - static_cast<double>(i);
        knots.push_back((1.0 - alpha) * parameters[i - 1] + alpha * parameters[i]);
    }
    knots.insert(knots.end(), p + 1, 1.0);
    return knots;
}

std::vector<Point3> LeastSquaresPoles(const std::span<const Point3> points,
                                      const std::span<const double> parameters,
                                      const std::span<const double> knots, const int degree,
                                      const std::size_t poleCount)
{
    auto [poles, fixed] = FixedPoles(points, knots, degree, poleCount, std::nullopt);
    SolveFreePoles(points, parameters, knots, degree, poles, fixed);
    return poles;
}

Result<std::size_t> FewestPolesWithin(const std::size_t fewest, const std::size_t most,
                                      const double budget,
                                      const std::function<Result<double>(std::size_t)>& deviationOf,
                                      const std::string_view where)
{
    // Doubling from `fewest` until a count passes or `most` fails.
    std::size_t failing = 0;
    std::size_t count = fewest;
    double best = 0.0;
    while (true) {
        STEREON_TRY(const double deviation, deviationOf(count));
        if (deviation <= budget) {
            break;
        }
        best = deviation;
        failing = count;
        if (count >= most) {
            return std::unexpected(
                KernelError(ErrorCode::ToleranceExceeded,
                            std::string{where} + ": even the interpolant exceeds the budget")
                    .WithTolerance(budget)
                    .WithTolerance(best));
        }
        count = std::min(2 * count, most);
    }
    // Bisection between the last failure and the first success, keeping the
    // passing half.
    while (count - failing > 1 && failing != 0) {
        const std::size_t mid = failing + (count - failing) / 2;
        STEREON_TRY(const double deviation, deviationOf(mid));
        if (deviation <= budget) {
            count = mid;
        } else {
            failing = mid;
        }
    }
    return count;
}

} // namespace Stereon::Detail
