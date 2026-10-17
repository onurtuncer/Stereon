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
#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Surface;

import :BSplineSurface;
import :Fitting;

namespace Stereon {

// See src/curve/Fitting.cpp: instantiated here before any use in a function.
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

/// The points of line `fixed` across one direction of the grid: along u
/// (varying k at fixed l) or along v.
[[nodiscard]] std::vector<Point3> GridLine(const std::span<const Point3> grid,
                                           const std::size_t countU, const std::size_t countV,
                                           const bool alongU, const std::size_t fixed)
{
    std::vector<Point3> line;
    const std::size_t count = alongU ? countU : countV;
    for (std::size_t index = 0; index < count; ++index) {
        line.push_back(grid[alongU ? index * countV + fixed : fixed * countV + index]);
    }
    return line;
}

/// The mean over the lines of one direction of their parameters, skipping
/// lines with a zero chord; empty when no line is usable (A9.3).
[[nodiscard]] std::vector<double> MeanParameters(const std::span<const Point3> grid,
                                                 const std::size_t countU, const std::size_t countV,
                                                 const bool alongU,
                                                 const Parametrization parametrization)
{
    const std::size_t count = alongU ? countU : countV;
    const std::size_t lines = alongU ? countV : countU;
    std::vector<double> mean(count, 0.0);
    std::size_t used = 0;
    for (std::size_t fixed = 0; fixed < lines; ++fixed) {
        const std::vector<Point3> line = GridLine(grid, countU, countV, alongU, fixed);
        const std::vector<double> t = Detail::ParametersOf(line, parametrization);
        if (t.empty()) {
            continue;
        }
        for (std::size_t k = 0; k < count; ++k) {
            mean[k] += t[k];
        }
        ++used;
    }
    if (used == 0) {
        return {};
    }
    for (double& t : mean) {
        t /= static_cast<double>(used);
    }
    mean.front() = 0.0;
    mean.back() = 1.0;
    return mean;
}

[[nodiscard]] Result<void> CheckGrid(const std::span<const Point3> grid, const std::size_t countU,
                                     const std::size_t countV, const int degreeU, const int degreeV,
                                     const std::string_view where)
{
    const bool degrees =
        degreeU >= 1 && degreeU <= BSplineMaxDegree && degreeV >= 1 && degreeV <= BSplineMaxDegree;
    if (!degrees || grid.size() != countU * countV ||
        countU < static_cast<std::size_t>(degreeU) + 1 ||
        countV < static_cast<std::size_t>(degreeV) + 1) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} + ": the degrees must be between 1 and 25 and the grid "
                                         "countU x countV with at least degree + 1 points per "
                                         "direction");
    }
    for (const Point3& p : grid) {
        STEREON_TRY_VOID(Detail::CheckFinite(p, "every point", where));
    }
    return {};
}

/// The grid's shared parameters along u and along v.
struct GridParameters {
    std::vector<double> U;
    std::vector<double> V;
};

/// The grid validated and given parameters, or why it cannot be fitted.
[[nodiscard]] Result<GridParameters>
CheckedGridParameters(const std::span<const Point3> grid, const std::size_t countU,
                      const std::size_t countV, const int degreeU, const int degreeV,
                      const Parametrization parametrization, const std::string_view where)
{
    STEREON_TRY_VOID(CheckGrid(grid, countU, countV, degreeU, degreeV, where));
    GridParameters t{MeanParameters(grid, countU, countV, true, parametrization),
                     MeanParameters(grid, countU, countV, false, parametrization)};
    if (t.U.empty() || t.V.empty()) {
        return Fail(ErrorCode::InvalidInput,
                    std::string{where} + ": every line of one direction has coincident points");
    }
    return t;
}

/// How one direction of a tensor-product fit turns a line of points into a
/// line of poles.
using LineFit = std::function<std::vector<Point3>(std::span<const Point3>)>;

/// The poles of a tensor-product fit: `alongU` on every grid line along u
/// gives poleCountU x countV intermediate points, then `alongV` on every
/// line of those along v gives the poles, row by row (A9.4 for
/// interpolation, section 9.4.3 for least squares).
[[nodiscard]] std::vector<Point3> TensorPoles(const std::span<const Point3> grid,
                                              const std::size_t countU, const std::size_t countV,
                                              const std::size_t poleCountU,
                                              const std::size_t poleCountV, const LineFit& alongU,
                                              const LineFit& alongV)
{
    std::vector<Point3> stage(poleCountU * countV);
    for (std::size_t l = 0; l < countV; ++l) {
        const std::vector<Point3> poles = alongU(GridLine(grid, countU, countV, true, l));
        for (std::size_t k = 0; k < poleCountU; ++k) {
            stage[k * countV + l] = poles[k];
        }
    }
    std::vector<Point3> poles;
    poles.reserve(poleCountU * poleCountV);
    for (std::size_t k = 0; k < poleCountU; ++k) {
        const std::vector<Point3> line = alongV(GridLine(stage, poleCountU, countV, false, k));
        poles.insert(poles.end(), line.begin(), line.end());
    }
    return poles;
}

/// The fit with its deviation measured over the grid.
[[nodiscard]] SurfaceFit Measured(BSplineSurface surface, const std::span<const Point3> grid,
                                  const std::size_t countU, const std::size_t countV,
                                  GridParameters t)
{
    double deviation = 0.0;
    for (std::size_t k = 0; k < countU; ++k) {
        for (std::size_t l = 0; l < countV; ++l) {
            deviation = std::max(
                deviation, Distance(Eval(surface, Point2{t.U[k], t.V[l]}), grid[k * countV + l]));
        }
    }
    return SurfaceFit{std::move(surface), std::move(t.U), std::move(t.V), deviation};
}

/// One direction of a least-squares surface fit: knots and a line fit for
/// poleCount poles at these parameters, interpolating when there are as
/// many poles as points.
struct Direction {
    std::vector<double> Knots;
    LineFit Fit;
};

[[nodiscard]] Direction LeastSquaresDirection(const std::span<const double> t, const int degree,
                                              const std::size_t poleCount)
{
    if (poleCount == t.size()) {
        Direction direction{Detail::AveragedKnots(t, degree), {}};
        direction.Fit = [t, degree, knots = direction.Knots](const std::span<const Point3> line) {
            return Detail::InterpolatingPoles(line, t, knots, degree);
        };
        return direction;
    }
    Direction direction{Detail::SpreadKnots(t, degree, poleCount), {}};
    direction.Fit = [t, degree, poleCount,
                     knots = direction.Knots](const std::span<const Point3> line) {
        return Detail::LeastSquaresPoles(line, t, knots, degree, poleCount);
    };
    return direction;
}

} // namespace

Result<SurfaceFit> InterpolateSurface(const std::span<const Point3> grid, const std::size_t countU,
                                      const std::size_t countV, const int degreeU,
                                      const int degreeV, const Parametrization parametrization)
{
    STEREON_TRY(GridParameters t, CheckedGridParameters(grid, countU, countV, degreeU, degreeV,
                                                        parametrization, "InterpolateSurface"));
    std::vector<double> knotsU = Detail::AveragedKnots(t.U, degreeU);
    std::vector<double> knotsV = Detail::AveragedKnots(t.V, degreeV);
    // Along u for every column of the grid, then along v for every row of
    // the result: the tensor product of two curve interpolations.
    std::vector<Point3> poles = TensorPoles(
        grid, countU, countV, countU, countV,
        [&](const std::span<const Point3> line) {
            return Detail::InterpolatingPoles(line, t.U, knotsU, degreeU);
        },
        [&](const std::span<const Point3> line) {
            return Detail::InterpolatingPoles(line, t.V, knotsV, degreeV);
        });
    BSplineSurface surface{degreeU,           degreeV,          std::move(knotsU),
                           std::move(knotsV), std::move(poles), {}};
    return Measured(std::move(surface), grid, countU, countV, std::move(t));
}

Result<SurfaceFit> FitSurfaceLeastSquares(const std::span<const Point3> grid,
                                          const std::size_t countU, const std::size_t countV,
                                          const int degreeU, const int degreeV,
                                          const std::size_t poleCountU,
                                          const std::size_t poleCountV,
                                          const Parametrization parametrization)
{
    STEREON_TRY(GridParameters t, CheckedGridParameters(grid, countU, countV, degreeU, degreeV,
                                                        parametrization, "FitSurfaceLeastSquares"));
    if (poleCountU < static_cast<std::size_t>(degreeU) + 1 || poleCountU > countU ||
        poleCountV < static_cast<std::size_t>(degreeV) + 1 || poleCountV > countV) {
        return Fail(ErrorCode::InvalidInput,
                    "FitSurfaceLeastSquares: each pole count must be between its degree + 1 and "
                    "its grid count");
    }
    Direction u = LeastSquaresDirection(t.U, degreeU, poleCountU);
    Direction v = LeastSquaresDirection(t.V, degreeV, poleCountV);
    std::vector<Point3> poles =
        TensorPoles(grid, countU, countV, poleCountU, poleCountV, u.Fit, v.Fit);
    BSplineSurface surface{degreeU,          degreeV, std::move(u.Knots), std::move(v.Knots),
                           std::move(poles), {}};
    return Measured(std::move(surface), grid, countU, countV, std::move(t));
}

Result<SurfaceFit> ApproximateSurface(const std::span<const Point3> grid, const std::size_t countU,
                                      const std::size_t countV, const int degreeU,
                                      const int degreeV, const ApproximationBudget budget,
                                      const Parametrization parametrization)
{
    constexpr std::string_view where = "ApproximateSurface";
    if (!std::isfinite(budget.Deviation) || !(budget.Deviation > 0.0)) {
        return Fail(ErrorCode::InvalidInput,
                    "ApproximateSurface: the budget must be finite and positive");
    }
    STEREON_TRY_VOID(CheckGrid(grid, countU, countV, degreeU, degreeV, where));
    const auto deviationWith = [&](const std::size_t poleCountU,
                                   const std::size_t poleCountV) -> Result<double> {
        STEREON_TRY(const SurfaceFit fit,
                    FitSurfaceLeastSquares(grid, countU, countV, degreeU, degreeV, poleCountU,
                                           poleCountV, parametrization));
        return fit.MaxDeviation;
    };
    // Along u first, with v interpolated so that only u limits the fit;
    // then along v with the u count found.
    STEREON_TRY(const std::size_t poleCountU,
                Detail::FewestPolesWithin(
                    static_cast<std::size_t>(degreeU) + 1, countU, budget.Deviation,
                    [&](const std::size_t count) { return deviationWith(count, countV); }, where));
    STEREON_TRY(const std::size_t poleCountV,
                Detail::FewestPolesWithin(
                    static_cast<std::size_t>(degreeV) + 1, countV, budget.Deviation,
                    [&](const std::size_t count) { return deviationWith(poleCountU, count); },
                    where));
    return FitSurfaceLeastSquares(grid, countU, countV, degreeU, degreeV, poleCountU, poleCountV,
                                  parametrization);
}

} // namespace Stereon
