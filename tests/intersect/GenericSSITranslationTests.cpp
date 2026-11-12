// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// SSI-translated-boundary: far from the origin the generic SSI refused,
// IntersectionFailed near a branch's end on the domain's boundary, where the
// same pair at the origin is answered. The answer, held to a closed form
// written here, must be the same under exact translations (ADR-0053, MR-01)
// and in both orders (MR-02).

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Context;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();
constexpr double tolerance = 1e-6;

/// z = alpha x'^2 + beta y'^2 + lift over [low, low + width]^2, x' = x - 1/2,
/// biquadratic: the Bernstein coefficients of (a + width t)^2, a = low - 1/2,
/// are a^2, a^2 + a width and (a + width)^2. Moved by `shift`.
BSplineSurface Quadric(const double alpha, const double beta, const double lift, const double low,
                       const double width, const Point3& shift)
{
    const double a = low - 0.5;
    const std::array<double, 3> square{a * a, a * a + a * width, (a + width) * (a + width)};
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({shift.X + low + static_cast<double>(i) * width / 2.0,
                             shift.Y + low + static_cast<double>(j) * width / 2.0,
                             shift.Z + lift + alpha * square.at(i) + beta * square.at(j)});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

} // namespace

TEST_CASE("SSI-translated-boundary: far translations answer as at the origin",
          "[intersect][generic][metamorphic][regression]")
{
    // 2x'^2 - 3y'^2 against x'^2 + y'^2 + 1/100 (the bowl over a wider
    // square): x'^2 - 4y'^2 = 1/100, two branches x' = +-sqrt(1/100 + 4y'^2)
    // across the unit square, ending on its sides x = 0 and x = 1.
    for (const Point3 shift :
         {Point3{0.0, 0.0, 0.0}, Point3{16.0, -32.0, 8.0}, Point3{64.0, 0.0, 0.0},
          Point3{0.0, 0.0, 1024.0}, Point3{1024.0, -2048.0, 512.0}, Point3{-4096.0, 256.0, 1024.0},
          Point3{480000.0, -300000.0, 210000.0}}) {
        const Surface a{Quadric(2.0, -3.0, 0.0, 0.0, 1.0, shift)};
        const Surface b{Quadric(1.0, 1.0, 0.01, -0.5, 2.0, shift)};
        for (const bool swapped : {false, true}) {
            INFO("shift (" << shift.X << ", " << shift.Y << ", " << shift.Z << ")"
                           << (swapped ? " swapped" : ""));
            const auto result =
                swapped ? Intersect(b, Box2::Everything(), a, Box2::Everything(), tolerance, ctx)
                        : Intersect(a, Box2::Everything(), b, Box2::Everything(), tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->Branches.size() == 2);
            REQUIRE(IntersectTest::ConnectivityOf(*result, 1e-9).Curves == 2);
            for (const Stereon::IntersectionBranch& branch : result->Branches) {
                REQUIRE(branch.IsCertified());
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Point3 p = Eval(branch.Curve3d, t);
                    const double x = p.X - shift.X - 0.5;
                    const double y = p.Y - shift.Y - 0.5;
                    REQUIRE(std::abs(x * x - 4.0 * y * y - 0.01) <= 8.0 * ctx.GetResolution());
                    REQUIRE(std::abs(p.Z - shift.Z - (2.0 * x * x - 3.0 * y * y)) <=
                            4.0 * ctx.GetResolution());
                }
            }
        }
    }
}
