// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Surface interpolation: through every grid point, reproducing a surface of
// its own space, and coping with the collapsed rows of a sphere's poles.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::ErrorCode;
using Stereon::Parametrization;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::SurfaceFit;

TEST_CASE("surface interpolation passes through every grid point", "[surface][fitting][property]")
{
    rc::prop("random wavy grids of every degree", [] {
        const int p = *rc::gen::inRange(1, 4);
        const int q = *rc::gen::inRange(1, 4);
        const auto countU = static_cast<std::size_t>(p + 1 + *rc::gen::inRange(0, 8));
        const auto countV = static_cast<std::size_t>(q + 1 + *rc::gen::inRange(0, 8));
        std::vector<Point3> grid;
        for (std::size_t k = 0; k < countU; ++k) {
            for (std::size_t l = 0; l < countV; ++l) {
                const double x = static_cast<double>(k) + *rc::gen::inRange(-20, 20) / 100.0;
                const double y = static_cast<double>(l) + *rc::gen::inRange(-20, 20) / 100.0;
                grid.push_back({x, y, std::sin(x) * std::cos(y)});
            }
        }
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
}

TEST_CASE("surface interpolation copes with a sphere's collapsed poles", "[surface][fitting]")
{
    // Rows of constant latitude; the first and last are single points
    // repeated, so their chords vanish and they are left out of the
    // parameter means.
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
    const auto fit = Stereon::InterpolateSurface(grid, along, around, 3, 3);
    REQUIRE(fit.has_value());
    REQUIRE(fit->MaxDeviation <= 1e-12);
}

TEST_CASE("surface interpolation validates its input", "[surface][fitting][factory]")
{
    const std::vector<Point3> grid(12, Point3{});
    const auto invalid = ErrorCode::InvalidInput;
    REQUIRE(Stereon::InterpolateSurface(grid, 3, 3, 1, 1).error().GetCode() == invalid); // 9 != 12
    REQUIRE(Stereon::InterpolateSurface(grid, 4, 3, 3, 3).error().GetCode() == invalid); // too few
    REQUIRE(Stereon::InterpolateSurface(grid, 4, 3, 1, 1).error().GetCode() ==
            invalid); // all equal
    REQUIRE(Stereon::InterpolateSurface(grid, 4, 3, 1, 1, Parametrization::Uniform).has_value());
}
