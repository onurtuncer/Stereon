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

// The fits of closed, periodic curves. They are a unit of their own because
// everything that distinguishes them from the open fits of Fitting.cpp is
// structural: the parameters run round a closing chord, the knots repeat
// with the period, every pole is free, and the linear systems are cyclic,
// with entries that wrap from the last poles to the first. A cyclic banded
// system is solved here by bordering (Golub and Van Loan, section 4.8): the
// last p unknowns are split off, the remaining band is factored once and
// applied to p + 1 right-hand sides, and a p x p Schur complement closes the
// loop. The banded factorisation is therefore kept apart from its
// substitution, which the single-use solver of Fitting.cpp has no need of.

namespace Stereon {

// See Fitting.cpp: instantiated before any arithmetic on it in this unit.
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

// ---------------------------------------------------------------------------
// Banded and cyclic banded solvers

/// An n x n matrix with non-zeros only within `half` of the diagonal, stored
/// row by row in 2 half + 1 columns, factored in place without pivoting and
/// then applied to any number of right-hand sides. Without pivoting the band
/// is kept; it is stable for the matrices used here, B-spline collocation
/// matrices (totally positive) and normal equations (symmetric positive
/// definite), and for their leading principal blocks.
class BandLU {
public:
    BandLU(const std::size_t n, const std::size_t half)
        : m_N(n)
        , m_Half(half)
        , m_Entries(n * (2 * half + 1), 0.0)
    {}

    [[nodiscard]] double& At(const std::size_t i, const std::size_t j)
    {
        return m_Entries[i * (2 * m_Half + 1) + (j + m_Half - i)];
    }

    /// Gaussian elimination without pivoting: the multipliers replace the
    /// entries below the diagonal.
    void Factor()
    {
        for (std::size_t k = 0; k < m_N; ++k) {
            for (std::size_t i = k + 1; i < m_N && i <= k + m_Half; ++i) {
                const double factor = At(i, k) / At(k, k);
                At(i, k) = factor;
                for (std::size_t j = k + 1; j < m_N && j <= k + m_Half; ++j) {
                    At(i, j) -= factor * At(k, j);
                }
            }
        }
    }

    /// Solves with the factored matrix, in place, for scalar or vector
    /// right-hand sides.
    template <typename V>
    void Solve(std::vector<V>& b) const
    {
        for (std::size_t k = 0; k < m_N; ++k) {
            for (std::size_t i = k + 1; i < m_N && i <= k + m_Half; ++i) {
                b[i] -= Entry(i, k) * b[k];
            }
        }
        for (std::size_t r = m_N; r > 0; --r) {
            const std::size_t i = r - 1;
            V sum = b[i];
            for (std::size_t j = i + 1; j < m_N && j <= i + m_Half; ++j) {
                sum -= Entry(i, j) * b[j];
            }
            b[i] = sum / Entry(i, i);
        }
    }

private:
    [[nodiscard]] double Entry(const std::size_t i, const std::size_t j) const
    {
        return m_Entries[i * (2 * m_Half + 1) + (j + m_Half - i)];
    }

    std::size_t m_N;
    std::size_t m_Half;
    std::vector<double> m_Entries;
};

/// Solves the dense system a x = b, a of size n x n row by row, by Gaussian
/// elimination with partial pivoting; false when a pivot is zero or not
/// finite, which means the matrix is singular. For the small Schur
/// complements of the cyclic solver, n <= BSplineMaxDegree.
[[nodiscard]] bool SolveDense(std::vector<double>& a, std::vector<Vec3>& b, const std::size_t n)
{
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(a[i * n + k]) > std::abs(a[pivot * n + k])) {
                pivot = i;
            }
        }
        if (!(std::abs(a[pivot * n + k]) > 0.0) || !std::isfinite(a[pivot * n + k])) {
            return false;
        }
        if (pivot != k) {
            std::swap_ranges(a.begin() + static_cast<std::ptrdiff_t>(k * n),
                             a.begin() + static_cast<std::ptrdiff_t>((k + 1) * n),
                             a.begin() + static_cast<std::ptrdiff_t>(pivot * n));
            std::swap(b[k], b[pivot]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = a[i * n + k] / a[k * n + k];
            for (std::size_t j = k + 1; j < n; ++j) {
                a[i * n + j] -= factor * a[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    for (std::size_t r = n; r > 0; --r) {
        const std::size_t i = r - 1;
        Vec3 sum = b[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= a[i * n + j] * b[j];
        }
        b[i] = sum / a[i * n + i];
    }
    return true;
}

/// An n x n matrix, n >= half + 1, whose entry (i, j) may be non-zero only
/// when i and j are within `half` of each other cyclically: the systems of
/// a periodic spline, whose last poles neighbour its first. It is stored
/// bordered. An entry that wraps round has its row or its column among the
/// last `half`, so the leading (n - half) x (n - half) block B is an
/// ordinary band matrix, and the rest is the border: the right block C, the
/// bottom block D and the corner E, all small.
class CyclicBandMatrix {
public:
    CyclicBandMatrix(const std::size_t n, const std::size_t half)
        : m_Inner(n - half)
        , m_Border(half)
        , m_Band(m_Inner, half)
        , m_Right(m_Inner * half, 0.0)
        , m_Bottom(half * m_Inner, 0.0)
        , m_Corner(half * half, 0.0)
    {}

    [[nodiscard]] double& At(const std::size_t i, const std::size_t j)
    {
        if (i < m_Inner && j < m_Inner) {
            Detail::Require(i + m_Border >= j && j + m_Border >= i, "CyclicBandMatrix",
                            "an entry outside the cyclic band");
            return m_Band.At(i, j);
        }
        if (i < m_Inner) {
            return m_Right[i * m_Border + (j - m_Inner)];
        }
        if (j < m_Inner) {
            return m_Bottom[(i - m_Inner) * m_Inner + j];
        }
        return m_Corner[(i - m_Inner) * m_Border + (j - m_Inner)];
    }

    /// The solution of M x = b, or nothing when the system is singular.
    /// Bordering: with x = (x_1, x_2) and b = (b_1, b_2) split at the
    /// border, B x_1 + C x_2 = b_1 and D x_1 + E x_2 = b_2, so x_2 solves
    /// the Schur complement system (E - D B^-1 C) x_2 = b_2 - D B^-1 b_1 and
    /// x_1 = B^-1 (b_1 - C x_2). B is factored once for its half + 1
    /// right-hand sides, the columns of C and b_1.
    [[nodiscard]] std::optional<std::vector<Vec3>> Solve(const std::vector<Vec3>& b)
    {
        m_Band.Factor();
        const std::vector<std::vector<double>> x = BorderColumns();
        std::vector<Vec3> y(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(m_Inner));
        m_Band.Solve(y);
        std::vector<double> schur = m_Corner;
        std::vector<Vec3> r(b.begin() + static_cast<std::ptrdiff_t>(m_Inner), b.end());
        for (std::size_t a = 0; a < m_Border; ++a) {
            for (std::size_t j = 0; j < m_Inner; ++j) {
                const double d = m_Bottom[a * m_Inner + j];
                if (d == 0.0) {
                    continue;
                }
                r[a] -= d * y[j];
                for (std::size_t c = 0; c < m_Border; ++c) {
                    schur[a * m_Border + c] -= d * x[c][j];
                }
            }
        }
        if (!SolveDense(schur, r, m_Border)) {
            return std::nullopt;
        }
        for (std::size_t i = 0; i < m_Inner; ++i) {
            for (std::size_t c = 0; c < m_Border; ++c) {
                y[i] -= x[c][i] * r[c];
            }
        }
        y.insert(y.end(), r.begin(), r.end());
        return y;
    }

private:
    /// B^-1 C, one column per border unknown.
    [[nodiscard]] std::vector<std::vector<double>> BorderColumns() const
    {
        std::vector<std::vector<double>> columns(m_Border, std::vector<double>(m_Inner, 0.0));
        for (std::size_t c = 0; c < m_Border; ++c) {
            for (std::size_t i = 0; i < m_Inner; ++i) {
                columns[c][i] = m_Right[i * m_Border + c];
            }
            m_Band.Solve(columns[c]);
        }
        return columns;
    }

    std::size_t m_Inner;
    std::size_t m_Border;
    BandLU m_Band;
    std::vector<double> m_Right;
    std::vector<double> m_Bottom;
    std::vector<double> m_Corner;
};

// ---------------------------------------------------------------------------
// Validation and dimension handling

[[nodiscard]] Result<void> CheckClosedPoints(const std::span<const Point3> points, const int degree,
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

/// The parameters of the closed polygon through the points, c_0 = 0, ...,
/// c_m, and c_{m+1} = 1 after the closing chord back to the first point:
/// ParametersOf over the points with the first repeated at the end.
[[nodiscard]] Result<std::vector<double>> ClosedParameters(const std::span<const Point3> points,
                                                           const Parametrization parametrization,
                                                           const std::string_view where)
{
    std::vector<Point3> closed(points.begin(), points.end());
    closed.push_back(points.front());
    std::vector<double> parameters = Detail::ParametersOf(closed, parametrization);
    if (parameters.empty()) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} +
                        ": consecutive points coincide, the last and the first included, so "
                        "chord-length parameters are undefined");
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

// ---------------------------------------------------------------------------
// Parameters and knots over one period

/// Where a periodic fit puts its points and its knots over the period
/// [0, 1].
struct PeriodLayout {
    /// The parameter of each point, increasing in [0, 1).
    std::vector<double> Parameters;
    /// The distinct knots v_0 = 0 < v_1 < ... < v_n = 1 of a fit with n
    /// distinct poles.
    std::vector<double> Knots;
};

/// The layout for the closed-polygon parameters c_0, ..., c_{m+1} and n
/// distinct poles, n <= m + 1. Think of the points as sitting at the
/// integer positions 0, ..., m + 1 of a line on which the parameter is
/// piecewise linear; the knots are then placed d = (m + 1) / n positions
/// apart, so that every knot span holds a point, as eq. 9.69 of The NURBS
/// Book does for open fits. With as many poles as points, d = 1, that puts
/// one knot at every point, which is the right rule for periodic
/// interpolation of odd degree (a periodic cubic spline with knots at its
/// data is the classical case, always uniquely solvable); for an even
/// degree the knots must instead lie halfway between the data
/// (docs/fitting.rst): with knots at the data, a periodic quadratic through an
/// even number of uniformly spaced points is exactly singular, because a
/// control polygon alternating up and down is invisible at the data. An
/// even degree therefore moves the knots back by half a position and moves
/// the parameters up by half the share of the closing chord, so that the
/// first knot still falls on 0 and the parameters still increase from the
/// first point and stay below 1.
[[nodiscard]] PeriodLayout LayOutPeriod(const std::span<const double> c, const int degree,
                                        const std::size_t poleCount)
{
    const std::size_t count = c.size() - 1;
    const bool even = degree % 2 == 0;
    const double shift = even ? (1.0 - c[count - 1]) / 2.0 : 0.0;
    const double offset = even ? -0.5 : 0.0;
    const double d = static_cast<double>(count) / static_cast<double>(poleCount);
    PeriodLayout layout;
    layout.Parameters.reserve(count);
    for (std::size_t k = 0; k < count; ++k) {
        layout.Parameters.push_back(c[k] + shift);
    }
    layout.Knots.reserve(poleCount + 1);
    layout.Knots.push_back(0.0);
    for (std::size_t j = 1; j < poleCount; ++j) {
        const double position = static_cast<double>(j) * d + offset;
        const auto i = static_cast<std::size_t>(position);
        const double alpha = position - static_cast<double>(i);
        layout.Knots.push_back((1.0 - alpha) * c[i] + alpha * c[i + 1] + shift);
    }
    layout.Knots.push_back(1.0);
    return layout;
}

// ---------------------------------------------------------------------------
// The cyclic systems

/// The periodic basis functions at t: the degree + 1 that are non-zero,
/// into `basis`, and the distinct pole index of the first, so that basis
/// function a multiplies pole (first + a) mod n. `knots` is the full knot
/// vector of the unclamped equivalent with n + p poles.
[[nodiscard]] std::size_t PeriodicBasis(const std::span<const double> knots, const int degree,
                                        const std::size_t poleCount, const double t,
                                        const std::span<double> basis)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t span = FindSpan(knots, degree, poleCount + p, t);
    BasisFunctions(knots, degree, span, t, basis);
    return span - p;
}

/// The distinct poles of the periodic interpolant: the cyclic collocation
/// system sum_i B_i(t_k) P_i = Q_k over the periodic basis functions B_i.
/// LayOutPeriod puts point k in knot span k, where the non-zero basis
/// functions are B_k, ..., B_{k+p} and the largest value is the central
/// one, B_{k+p/2}; the equation of point k is therefore placed in row
/// k + p/2, which puts that value on the diagonal. The leading block of the
/// bordered matrix is then an ordinary collocation matrix of consecutive
/// points against consecutive basis functions, totally positive and
/// satisfying the Schoenberg-Whitney conditions, so elimination without
/// pivoting is stable; with the equations in point order instead, the block
/// is nearly triangular with the smallest basis values as pivots, and its
/// back substitution grows exponentially.
[[nodiscard]] std::optional<std::vector<Vec3>>
PeriodicInterpolatingPoles(const std::span<const Point3> q, const std::span<const double> t,
                           const std::span<const double> knots, const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = q.size();
    const std::size_t centre = p / 2;
    CyclicBandMatrix matrix(n, p);
    std::vector<Vec3> rhs(n);
    std::array<double, BSplineMaxDegree + 1> basis; // filled by PeriodicBasis
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t first =
            PeriodicBasis(knots, degree, n, t[k], std::span{basis}.first(p + 1));
        const std::size_t row = (k + centre) % n;
        for (std::size_t a = 0; a <= p; ++a) {
            matrix.At(row, (first + a) % n) = basis[a];
        }
        rhs[row] = ToVec(q[k]);
    }
    return matrix.Solve(rhs);
}

/// The distinct poles of the periodic least-squares fit: every pole is
/// free, so the normal equations N^T N P = N^T Q run over all of them, with
/// N_{ki} = B_i(t_k). Two periodic basis functions overlap only when their
/// indices are within p cyclically, so N^T N is a cyclic band matrix; it is
/// positive definite because every knot span holds a parameter.
[[nodiscard]] std::optional<std::vector<Vec3>>
PeriodicLeastSquaresPoles(const std::span<const Point3> q, const std::span<const double> t,
                          const std::span<const double> knots, const int degree,
                          const std::size_t poleCount)
{
    const auto p = static_cast<std::size_t>(degree);
    CyclicBandMatrix normal(poleCount, p);
    std::vector<Vec3> rhs(poleCount);
    std::array<double, BSplineMaxDegree + 1> basis; // filled by PeriodicBasis
    for (std::size_t k = 0; k < q.size(); ++k) {
        const std::size_t first =
            PeriodicBasis(knots, degree, poleCount, t[k], std::span{basis}.first(p + 1));
        for (std::size_t a = 0; a <= p; ++a) {
            const std::size_t i = (first + a) % poleCount;
            rhs[i] += basis[a] * ToVec(q[k]);
            for (std::size_t b = 0; b <= p; ++b) {
                normal.At(i, (first + b) % poleCount) += basis[a] * basis[b];
            }
        }
    }
    return normal.Solve(rhs);
}

/// The solved poles as points, or InvalidInput when the system broke down:
/// the parameters of the points were too close for it, which only rounding
/// can bring about once zero chords are rejected.
[[nodiscard]] Result<std::vector<Point3>> PolesFrom(std::optional<std::vector<Vec3>> solved,
                                                    const std::string_view where)
{
    std::vector<Point3> poles;
    if (solved) {
        poles.reserve(solved->size());
        for (const Vec3& v : *solved) {
            if (!std::isfinite(v.X) || !std::isfinite(v.Y) || !std::isfinite(v.Z)) {
                solved.reset();
                break;
            }
            poles.push_back(ToPoint(v));
        }
    }
    if (!solved) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} +
                        ": the parameters of the points are too close for the periodic "
                        "fitting system");
    }
    return poles;
}

// ---------------------------------------------------------------------------
// The fits

/// The periodic fit with poleCount distinct poles: the interpolant when
/// there are as many poles as points, least squares otherwise. The points
/// and the pole count have been checked.
[[nodiscard]] Result<CurveFit> PeriodicFit(const std::span<const Point3> points, const int degree,
                                           const std::size_t poleCount,
                                           const Parametrization parametrization,
                                           const std::string_view where)
{
    STEREON_TRY(std::vector<double> c, ClosedParameters(points, parametrization, where));
    PeriodLayout layout = LayOutPeriod(c, degree, poleCount);
    std::vector<double> knots = Detail::PeriodicKnotVector(layout.Knots, degree);
    const bool interpolating = poleCount == points.size();
    STEREON_TRY(const std::vector<Point3> poles,
                PolesFrom(interpolating
                              ? PeriodicInterpolatingPoles(points, layout.Parameters, knots, degree)
                              : PeriodicLeastSquaresPoles(points, layout.Parameters, knots, degree,
                                                          poleCount),
                          where));
    BSplineCurve curve{degree, std::move(knots), Detail::WrappedPoles(poles, degree), {}};
    double deviation = 0.0;
    for (std::size_t k = 0; k < points.size(); ++k) {
        deviation = std::max(deviation, Distance(Eval(curve, layout.Parameters[k]), points[k]));
    }
    return CurveFit{std::move(curve), std::move(layout.Parameters), deviation};
}

[[nodiscard]] Result<CurveFit> FitCountPeriodic(const std::span<const Point3> points,
                                                const int degree, const std::size_t poleCount,
                                                const Parametrization parametrization)
{
    constexpr std::string_view where = "FitLeastSquaresPeriodic";
    STEREON_TRY_VOID(CheckClosedPoints(points, degree, where));
    if (poleCount < static_cast<std::size_t>(degree) + 1 || poleCount > points.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "FitLeastSquaresPeriodic: the pole count must be between degree + 1 and "
                    "the number of points");
    }
    return PeriodicFit(points, degree, poleCount, parametrization, where);
}

} // namespace

Result<CurveFit> InterpolatePeriodic(const std::span<const Point3> points, const int degree,
                                     const Parametrization parametrization)
{
    constexpr std::string_view where = "InterpolatePeriodic";
    STEREON_TRY_VOID(CheckClosedPoints(points, degree, where));
    return PeriodicFit(points, degree, points.size(), parametrization, where);
}

Result<CurveFit2d> InterpolatePeriodic(const std::span<const Point2> points, const int degree,
                                       const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return InterpolatePeriodic(lifted, degree, parametrization);
    });
}

Result<CurveFit> FitLeastSquaresPeriodic(const std::span<const Point3> points, const int degree,
                                         const std::size_t poleCount,
                                         const Parametrization parametrization)
{
    return FitCountPeriodic(points, degree, poleCount, parametrization);
}

Result<CurveFit2d> FitLeastSquaresPeriodic(const std::span<const Point2> points, const int degree,
                                           const std::size_t poleCount,
                                           const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return FitCountPeriodic(lifted, degree, poleCount, parametrization);
    });
}

Result<CurveFit> ApproximatePeriodic(const std::span<const Point3> points, const int degree,
                                     const ApproximationBudget budget,
                                     const Parametrization parametrization)
{
    constexpr std::string_view where = "ApproximatePeriodic";
    if (!std::isfinite(budget.Deviation) || !(budget.Deviation > 0.0)) {
        return Fail(ErrorCode::InvalidInput,
                    "ApproximatePeriodic: the budget must be finite and positive");
    }
    STEREON_TRY_VOID(CheckClosedPoints(points, degree, where));
    STEREON_TRY(const std::size_t count,
                Detail::FewestPolesWithin(
                    static_cast<std::size_t>(degree) + 1, points.size(), budget.Deviation,
                    [&](const std::size_t poles) -> Result<double> {
                        STEREON_TRY(CurveFit fit,
                                    PeriodicFit(points, degree, poles, parametrization, where));
                        return fit.MaxDeviation;
                    },
                    where));
    return PeriodicFit(points, degree, count, parametrization, where);
}

Result<CurveFit2d> ApproximatePeriodic(const std::span<const Point2> points, const int degree,
                                       const ApproximationBudget budget,
                                       const Parametrization parametrization)
{
    return InPlane(points, [&](const std::span<const Point3> lifted) {
        return ApproximatePeriodic(lifted, degree, budget, parametrization);
    });
}

} // namespace Stereon

namespace Stereon::Detail {

std::vector<double> PeriodicKnotVector(const std::span<const double> periodKnots, const int degree)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = periodKnots.size() - 1;
    const double period = periodKnots[n] - periodKnots[0];
    std::vector<double> knots;
    knots.reserve(n + 2 * p + 1);
    for (std::size_t r = p; r >= 1; --r) {
        knots.push_back(periodKnots[n - r] - period);
    }
    knots.insert(knots.end(), periodKnots.begin(), periodKnots.end());
    for (std::size_t r = 1; r <= p; ++r) {
        knots.push_back(periodKnots[r] + period);
    }
    return knots;
}

std::vector<Point3> WrappedPoles(const std::span<const Point3> poles, const int degree)
{
    std::vector<Point3> wrapped(poles.begin(), poles.end());
    wrapped.insert(wrapped.end(), poles.begin(),
                   poles.begin() + static_cast<std::ptrdiff_t>(degree));
    return wrapped;
}

} // namespace Stereon::Detail
