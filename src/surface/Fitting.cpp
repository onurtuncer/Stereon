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
#include <cstddef>
#include <span>
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
                                     const std::size_t countV, const int degreeU, const int degreeV)
{
    constexpr std::string_view where = "InterpolateSurface";
    const bool degrees =
        degreeU >= 1 && degreeU <= BSplineMaxDegree && degreeV >= 1 && degreeV <= BSplineMaxDegree;
    if (!degrees || grid.size() != countU * countV ||
        countU < static_cast<std::size_t>(degreeU) + 1 ||
        countV < static_cast<std::size_t>(degreeV) + 1) {
        return Fail(ErrorCode::InvalidInput,
                    "InterpolateSurface: the degrees must be between 1 and 25 and the grid "
                    "countU x countV with at least degree + 1 points per direction");
    }
    for (const Point3& p : grid) {
        STEREON_TRY_VOID(Detail::CheckFinite(p, "every point", where));
    }
    return {};
}

} // namespace

Result<SurfaceFit> InterpolateSurface(const std::span<const Point3> grid, const std::size_t countU,
                                      const std::size_t countV, const int degreeU,
                                      const int degreeV, const Parametrization parametrization)
{
    STEREON_TRY_VOID(CheckGrid(grid, countU, countV, degreeU, degreeV));
    std::vector<double> tu = MeanParameters(grid, countU, countV, true, parametrization);
    std::vector<double> tv = MeanParameters(grid, countU, countV, false, parametrization);
    if (tu.empty() || tv.empty()) {
        return Fail(ErrorCode::InvalidInput,
                    "InterpolateSurface: every line of one direction has coincident points");
    }
    std::vector<double> knotsU = Detail::AveragedKnots(tu, degreeU);
    std::vector<double> knotsV = Detail::AveragedKnots(tv, degreeV);
    // Along u for every column of the grid, then along v for every row of
    // the result: the tensor product of two curve interpolations.
    std::vector<Point3> stage(grid.size());
    for (std::size_t l = 0; l < countV; ++l) {
        const std::vector<Point3> poles = Detail::InterpolatingPoles(
            GridLine(grid, countU, countV, true, l), tu, knotsU, degreeU);
        for (std::size_t k = 0; k < countU; ++k) {
            stage[k * countV + l] = poles[k];
        }
    }
    BSplineSurface surface{degreeU, degreeV, std::move(knotsU), std::move(knotsV), {}, {}};
    for (std::size_t k = 0; k < countU; ++k) {
        const std::vector<Point3> poles = Detail::InterpolatingPoles(
            GridLine(stage, countU, countV, false, k), tv, surface.KnotsV, degreeV);
        surface.Poles.insert(surface.Poles.end(), poles.begin(), poles.end());
    }
    double deviation = 0.0;
    for (std::size_t k = 0; k < countU; ++k) {
        for (std::size_t l = 0; l < countV; ++l) {
            deviation = std::max(
                deviation, Distance(Eval(surface, Point2{tu[k], tv[l]}), grid[k * countV + l]));
        }
    }
    return SurfaceFit{std::move(surface), std::move(tu), std::move(tv), deviation};
}

} // namespace Stereon
