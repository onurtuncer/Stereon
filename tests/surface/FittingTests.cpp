// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Surface fitting. Interpolation must pass through every grid point,
// reproduce a surface of its own space and cope with the collapsed rows of
// a sphere's poles; least squares must satisfy its normal equations over the
// interior poles and keep the corners; approximation must meet its budget
// and spend poles only in the direction that needs them.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::ApproximationBudget;
using Stereon::BSplineSurface;
using Stereon::ErrorCode;
using Stereon::Parametrization;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SurfaceFit;
using Stereon::Vec3;

namespace {

/// A jittered grid of countU x countV points on a wavy surface.
std::vector<Point3> RandomGrid(const std::size_t countU, const std::size_t countV)
{
    std::vector<Point3> grid;
    for (std::size_t k = 0; k < countU; ++k) {
        for (std::size_t l = 0; l < countV; ++l) {
            const double x = static_cast<double>(k) + *rc::gen::inRange(-20, 20) / 100.0;
            const double y = static_cast<double>(l) + *rc::gen::inRange(-20, 20) / 100.0;
            grid.push_back({x, y, std::sin(x) * std::cos(y)});
        }
    }
    return grid;
}

/// A 7 x 9 grid on a sphere from pole to pole: rows of constant latitude,
/// the first and last a single point repeated.
std::vector<Point3> SphereGrid()
{
    std::vector<Point3> grid;
    constexpr int around = 9;
    constexpr int along = 7;
    for (int k = 0; k < along; ++k) {
        const double latitude = -std::numbers::pi / 2 + k * std::numbers::pi / (along - 1);
        for (int l = 0; l < around; ++l) {
            const double longitude = l * std::numbers::pi / (around - 1);
            grid.push_back({std::cos(latitude) * std::cos(longitude),
                            std::cos(latitude) * std::sin(longitude), std::sin(latitude)});
        }
    }
    return grid;
}

double Scale(const std::span<const Point3> points)
{
    double scale = 1.0;
    for (const Point3& p : points) {
        scale = std::max({scale, std::abs(p.X), std::abs(p.Y), std::abs(p.Z)});
    }
    return scale;
}

/// The number of poles along u and along v.
struct PoleCounts {
    std::size_t U;
    std::size_t V;
};

PoleCounts PoleCountsOf(const BSplineSurface& s)
{
    return {s.KnotsU.size() - static_cast<std::size_t>(s.DegreeU) - 1,
            s.KnotsV.size() - static_cast<std::size_t>(s.DegreeV) - 1};
}

} // namespace

TEST_CASE("surface interpolation passes through every grid point", "[surface][fitting][property]")
{
    rc::prop("random wavy grids of every degree", [] {
        const int p = *rc::gen::inRange(1, 4);
        const int q = *rc::gen::inRange(1, 4);
        const auto countU = static_cast<std::size_t>(p + 1 + *rc::gen::inRange(0, 8));
        const auto countV = static_cast<std::size_t>(q + 1 + *rc::gen::inRange(0, 8));
        const std::vector<Point3> grid = RandomGrid(countU, countV);
        const auto fit = Stereon::InterpolateSurface(grid, countU, countV, p, q);
        RC_ASSERT(fit.has_value());
        RC_ASSERT(fit->MaxDeviation <= 1e-10);
        RC_ASSERT(fit->ParametersU.size() == countU);
        RC_ASSERT(fit->ParametersV.size() == countV);
    });
}

TEST_CASE("surface interpolation reproduces a surface of its own space", "[surface][fitting]")
{
    // (u, v, u^2 v - v) has degree (2, 1); sampled at uniform parameters
    // and interpolated at those parameters it comes back exactly.
    const auto surface = [](const double u, const double v) { return Point3{u, v, u * u * v - v}; };
    std::vector<Point3> grid;
    for (int k = 0; k <= 6; ++k) {
        for (int l = 0; l <= 4; ++l) {
            grid.push_back(surface(k / 6.0, l / 4.0));
        }
    }
    const SurfaceFit fit = *Stereon::InterpolateSurface(grid, 7, 5, 2, 1, Parametrization::Uniform);
    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 10; ++j) {
            const Point2 uv{i / 10.0, j / 10.0};
            REQUIRE(Distance(Eval(fit.Surface, uv), surface(uv.X, uv.Y)) <= 1e-12);
        }
    }
    // So does least squares with fewer poles, since the surface lies in the
    // smaller space too and the fit is its unique minimiser.
    const SurfaceFit lsq =
        *Stereon::FitSurfaceLeastSquares(grid, 7, 5, 2, 1, 4, 3, Parametrization::Uniform);
    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 10; ++j) {
            const Point2 uv{i / 10.0, j / 10.0};
            REQUIRE(Distance(Eval(lsq.Surface, uv), surface(uv.X, uv.Y)) <= 1e-12);
        }
    }
}

TEST_CASE("surface interpolation copes with a sphere's collapsed poles", "[surface][fitting]")
{
    // The first and last rows are single points repeated, so their chords
    // vanish and they are left out of the parameter means.
    const std::vector<Point3> grid = SphereGrid();
    const auto fit = Stereon::InterpolateSurface(grid, 7, 9, 3, 3);
    REQUIRE(fit.has_value());
    REQUIRE(fit->MaxDeviation <= 1e-12);
    // Least squares along the collapsed rows fits a constant to constant
    // data, exactly.
    const auto lsq = Stereon::FitSurfaceLeastSquares(grid, 7, 9, 3, 3, 5, 6);
    REQUIRE(lsq.has_value());
    REQUIRE(lsq->MaxDeviation <= 0.1);
    REQUIRE(Distance(Eval(lsq->Surface, Point2{0.0, 0.5}), grid.front()) <= 1e-12);
    REQUIRE(Distance(Eval(lsq->Surface, Point2{1.0, 0.5}), grid.back()) <= 1e-12);
}

TEST_CASE("surface least squares satisfies its normal equations", "[surface][fitting][property]")
{
    rc::prop("the residual is orthogonal to every interior tensor-product basis function", [] {
        const int p = *rc::gen::inRange(1, 4);
        const int q = *rc::gen::inRange(1, 4);
        const auto pu = static_cast<std::size_t>(p);
        const auto qv = static_cast<std::size_t>(q);
        // Grids with at least one more line than the fewest poles, so that
        // both directions are least squares rather than interpolation.
        const std::size_t countU = pu + 2 + *rc::gen::inRange<std::size_t>(0, 10);
        const std::size_t countV = qv + 2 + *rc::gen::inRange<std::size_t>(0, 10);
        const std::size_t nu = pu + 1 + *rc::gen::inRange<std::size_t>(0, countU - pu - 1);
        const std::size_t nv = qv + 1 + *rc::gen::inRange<std::size_t>(0, countV - qv - 1);
        const std::vector<Point3> grid = RandomGrid(countU, countV);
        const auto fit = Stereon::FitSurfaceLeastSquares(grid, countU, countV, p, q, nu, nv);
        RC_ASSERT(fit.has_value());
        const BSplineSurface& s = fit->Surface;
        RC_ASSERT(s.Poles.size() == nu * nv);
        // The corners are the corner points.
        RC_ASSERT(s.Poles.front() == grid.front());
        RC_ASSERT(s.Poles[nv - 1] == grid[countV - 1]);
        RC_ASSERT(s.Poles[(nu - 1) * nv] == grid[(countU - 1) * countV]);
        RC_ASSERT(s.Poles.back() == grid.back());
        // The gradient of the squared residual over the interior grid points
        // with respect to pole (i, j): sum over (k, l) of
        // N_i(u_k) M_j(v_l) (S(u_k, v_l) - Q_kl); zero for interior poles.
        std::vector<Vec3> gradient(nu * nv, Vec3{});
        std::vector<double> bu(pu + 1);
        std::vector<double> bv(qv + 1);
        for (std::size_t k = 1; k + 1 < countU; ++k) {
            for (std::size_t l = 1; l + 1 < countV; ++l) {
                const double u = fit->ParametersU[k];
                const double v = fit->ParametersV[l];
                const std::size_t spanU = Stereon::FindSpan(s.KnotsU, p, nu, u);
                const std::size_t spanV = Stereon::FindSpan(s.KnotsV, q, nv, v);
                Stereon::BasisFunctions(s.KnotsU, p, spanU, u, bu);
                Stereon::BasisFunctions(s.KnotsV, q, spanV, v, bv);
                const Vec3 residual = Eval(s, Point2{u, v}) - grid[k * countV + l];
                for (std::size_t a = 0; a < bu.size(); ++a) {
                    for (std::size_t b = 0; b < bv.size(); ++b) {
                        gradient[(spanU - pu + a) * nv + (spanV - qv + b)] +=
                            bu[a] * bv[b] * residual;
                    }
                }
            }
        }
        const double scale = std::max(Scale(grid), Scale(s.Poles));
        for (std::size_t i = 1; i + 1 < nu; ++i) {
            for (std::size_t j = 1; j + 1 < nv; ++j) {
                RC_ASSERT(Length(gradient[i * nv + j]) <= 1e-9 * scale);
            }
        }
    });
}

TEST_CASE("surface least squares with every pole is the interpolant", "[surface][fitting]")
{
    std::vector<Point3> grid;
    for (int k = 0; k < 6; ++k) {
        for (int l = 0; l < 8; ++l) {
            grid.push_back({0.7 * k, 0.5 * l, std::sin(0.7 * k) * std::cos(0.5 * l)});
        }
    }
    const auto interpolant = Stereon::InterpolateSurface(grid, 6, 8, 3, 2);
    const auto full = Stereon::FitSurfaceLeastSquares(grid, 6, 8, 3, 2, 6, 8);
    REQUIRE(full->Surface == interpolant->Surface);
    // One direction interpolated, the other not: still through every point
    // of the interpolated lines' ends.
    const auto mixed = Stereon::FitSurfaceLeastSquares(grid, 6, 8, 3, 2, 6, 5);
    REQUIRE(mixed.has_value());
    REQUIRE(PoleCountsOf(mixed->Surface).U == 6);
    REQUIRE(PoleCountsOf(mixed->Surface).V == 5);
}

TEST_CASE("surface approximation meets its budget", "[surface][fitting][property]")
{
    rc::prop("the deviation is within the budget", [] {
        const int p = *rc::gen::inRange(1, 4);
        const int q = *rc::gen::inRange(1, 4);
        const auto countU = static_cast<std::size_t>(p + 1 + *rc::gen::inRange(0, 12));
        const auto countV = static_cast<std::size_t>(q + 1 + *rc::gen::inRange(0, 12));
        const std::vector<Point3> grid = RandomGrid(countU, countV);
        const double budget = *rc::gen::inRange(1, 300) / 1000.0;
        const auto fit =
            Stereon::ApproximateSurface(grid, countU, countV, p, q, ApproximationBudget{budget});
        RC_ASSERT(fit.has_value());
        RC_ASSERT(fit->MaxDeviation <= budget);
        const PoleCounts counts = PoleCountsOf(fit->Surface);
        RC_ASSERT(counts.U <= countU && counts.V <= countV);
    });
}

TEST_CASE("surface approximation spends poles only where they are needed", "[surface][fitting]")
{
    // Straight along u, a sine along v: four cubic poles reproduce the
    // straight lines exactly, so u gets the fewest poles whatever v needs.
    std::vector<Point3> grid;
    for (int k = 0; k < 10; ++k) {
        for (int l = 0; l < 20; ++l) {
            grid.push_back({0.5 * k, 0.4 * l, std::sin(0.4 * l)});
        }
    }
    const auto fit = Stereon::ApproximateSurface(grid, 10, 20, 3, 3, ApproximationBudget{1e-6});
    REQUIRE(fit.has_value());
    REQUIRE(fit->MaxDeviation <= 1e-6);
    const PoleCounts counts = PoleCountsOf(fit->Surface);
    REQUIRE(counts.U == 4);
    REQUIRE(counts.V > 4);
}

TEST_CASE("surface fitting validates its input", "[surface][fitting][factory]")
{
    const std::vector<Point3> coincident(12, Point3{});
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::InterpolateSurface(coincident, 3, 3, 1, 1).error().GetCode() ==
            invalid); // 9 != 12
    REQUIRE(Stereon::InterpolateSurface(coincident, 4, 3, 3, 3).error().GetCode() ==
            invalid); // too few
    REQUIRE(Stereon::InterpolateSurface(coincident, 4, 3, 1, 1).error().GetCode() ==
            invalid); // all equal
    REQUIRE(
        Stereon::InterpolateSurface(coincident, 4, 3, 1, 1, Parametrization::Uniform).has_value());
    std::vector<Point3> grid;
    for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 3; ++l) {
            grid.push_back({1.0 * k, 1.0 * l, 0.0});
        }
    }
    REQUIRE(Stereon::FitSurfaceLeastSquares(grid, 4, 3, 1, 1, 2, 2).has_value());
    REQUIRE(Stereon::FitSurfaceLeastSquares(grid, 4, 3, 1, 1, 1, 2).error().GetCode() == invalid);
    REQUIRE(Stereon::FitSurfaceLeastSquares(grid, 4, 3, 1, 1, 5, 2).error().GetCode() == invalid);
    REQUIRE(Stereon::FitSurfaceLeastSquares(grid, 4, 3, 1, 1, 2, 4).error().GetCode() == invalid);
    REQUIRE(
        Stereon::ApproximateSurface(grid, 4, 3, 1, 1, ApproximationBudget{0.0}).error().GetCode() ==
        invalid);
    REQUIRE(
        Stereon::ApproximateSurface(grid, 4, 3, 1, 1, ApproximationBudget{NAN}).error().GetCode() ==
        invalid);
    REQUIRE(Stereon::ApproximateSurface(coincident, 4, 3, 1, 1, ApproximationBudget{1.0})
                .error()
                .GetCode() == invalid);
}
