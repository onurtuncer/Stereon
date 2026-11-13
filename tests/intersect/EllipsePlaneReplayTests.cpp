// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

import Stereon.Intersect;

#include "Checks.hpp"

namespace {
using namespace Stereon;
constexpr double pi = std::numbers::pi;

void CheckCrossings(const Curve& curve, const Surface& plane, const Point3& center,
                    const double radius, const Context& policy)
{
    const auto result = Intersect(curve, plane, policy);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Overlaps.empty());
    // The closed parameter interval may retain both representatives of the
    // periodic seam. Only the two independently known spatial roots are valid.
    REQUIRE(result->Points.size() >= 2);
    REQUIRE(result->Points.size() <= 3);
    bool left = false;
    bool right = false;
    for (const auto& point : result->Points) {
        const Point3 lo = center + Vec3{-radius, 0, 0};
        const Point3 hi = center + Vec3{radius, 0, 0};
        const bool atLo = Distance(point.Position, lo) <= policy.GetResolution();
        const bool atHi = Distance(point.Position, hi) <= policy.GetResolution();
        CHECK((atLo || atHi));
        left = left || atLo;
        right = right || atHi;
        CHECK(point.Decided);
        CHECK(point.Contact == ContactType::Crossing);
        CHECK(Contains(Domain(curve), point.Parameter));
        CHECK(Distance(Eval(curve, point.Parameter), point.Position) <= policy.GetResolution());
        CHECK(Distance(Eval(plane, point.OnSurface), point.Position) <= policy.GetResolution());
    }
    CHECK(left);
    CHECK(right);
}
} // namespace

TEST_CASE("MR macOS ellipse plane shrunk linearization fixture resolves",
          "[intersect][metamorphic][regression][conic-surface]")
{
    // PR-272 macOS run 43273260370, RapidCheck seed 17041817369036985825.
    // Shrunk draws: 0,-16,-40,-40,-13,15,4,2,0,0,-40,-40,-13.
    // Ellipse and horizontal plane share their center; sin(t)=0 gives the
    // two transverse crossings independently of the tiny cos(-pi/2) tilt.
    const double phi = -16.0 * (2.0 * pi) / 64.0;
    for (const double scale : {0.25, 1.0, 4.0}) {
        for (const double offset : {0.0, 8.0}) {
            CAPTURE(scale, offset);
            const Point3 center{scale * (-40.0 / 32.0) + offset, scale * (-40.0 / 32.0),
                                scale * (-13.0 / 32.0)};
            const Ellipse ellipse{center, scale * 1.5, scale * (1.5 * 4.0 / 10.0), Vec3{1, 0, 0},
                                  Vec3{0, std::cos(phi), std::sin(phi)}};
            const Surface plane{Plane{center, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
            const Context policy =
                Context::Default()
                    .WithResolution(Context::DefaultResolution * scale)
                    .WithToleranceCeiling(Context::DefaultToleranceCeiling * scale)
                    .WithScheduler(Scheduler::Serial);
            CheckCrossings(Curve{ellipse}, plane, center, scale * 1.5, policy);
            CheckCrossings(Curve{Reversed(ellipse)}, plane, center, scale * 1.5, policy);
        }
    }
}

TEST_CASE("ellipse plane endpoint proof respects trims and separated supports",
          "[intersect][metamorphic][regression][conic-surface]")
{
    const Point3 center{-1.25, -1.25, -13.0 / 32.0};
    const Ellipse ellipse{center, 1.5, 0.6, Vec3{1, 0, 0}, Vec3{0, 0, -1}};
    const Curve curve{ellipse};
    const Surface plane{Plane{center, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
    const Context policy = Context::Default().WithScheduler(Scheduler::Serial);
    // The seam root is outside this strictly positive angular interval.
    const auto trimmed = Intersect(curve, Interval{0.01, 0.25}, plane, Domain(plane), policy);
    REQUIRE(trimmed.has_value());
    CHECK(trimmed->IsEmpty());
    // Shift beyond the complete ellipse's z extent: no supporting root.
    const Surface apart{Plane{center + Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}};
    const auto separated = Intersect(curve, apart, policy);
    REQUIRE(separated.has_value());
    CHECK(separated->IsEmpty());
    // A bounded interval straddling zero contains exactly the simple root.
    const auto crossing = Intersect(curve, Interval{-0.25, 0.25}, plane, Domain(plane), policy);
    REQUIRE(crossing.has_value());
    REQUIRE(crossing->Points.size() == 1);
    CHECK(crossing->Points.front().Decided);
    CHECK(Distance(crossing->Points.front().Position, center + Vec3{1.5, 0, 0}) <=
          policy.GetResolution());
}
