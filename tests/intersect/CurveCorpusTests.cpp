// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every case of corpus/curve-curve/, in both orders: at the case's known
// closest approach the curves meet exactly when its gap is within the
// resolution. Other meetings of the pair elsewhere are not the case's claim.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"
#include "CurveCorpus.hpp"
#include "CurveWitnessChecks.hpp"

namespace {

/// How far a meeting may be placed from the closest approach, in parameter.
constexpr double Near = 1e-4;

bool MeetsNear(const Stereon::CurveCurveIntersection& r, const double a, const double b,
               const bool swapped)
{
    for (const Stereon::CurveCurvePoint& p : r.Points) {
        const double onA = swapped ? p.ParameterB : p.ParameterA;
        const double onB = swapped ? p.ParameterA : p.ParameterB;
        if (std::abs(onA - a) <= Near && std::abs(onB - b) <= Near) {
            return true;
        }
    }
    return false;
}

} // namespace

TEST_CASE("every curve-curve corpus case meets exactly when its gap is within the resolution",
          "[intersect][curve-curve][corpus]")
{
    const auto files = StereonTest::CurveCorpusFiles();
    REQUIRE(files.size() >= 168);
    for (const auto& file : files) {
        const StereonTest::JsonValue c = StereonTest::ReadJsonFile(file);
        const auto& curves = c["inputs"]["curves"];
        const Stereon::Curve a = StereonTest::CurveFrom(curves[0]);
        const Stereon::Curve b = StereonTest::CurveFrom(curves[1]);
        const Stereon::Interval da = StereonTest::DomainFrom(curves[0]);
        const Stereon::Interval db = StereonTest::DomainFrom(curves[1]);
        const double ta = c["expected"]["parameters"][0].AsNumber();
        const double tb = c["expected"]["parameters"][1].AsNumber();
        const bool meets = c["expected"]["meets"].AsBool();
        for (const bool swapped : {false, true}) {
            INFO(c["id"].AsString() << (swapped ? " swapped" : ""));
            const auto r = swapped ? Intersect(b, db, a, da, IntersectTest::ctx)
                                   : Intersect(a, da, b, db, IntersectTest::ctx);
            INFO(IntersectTest::ErrorText(r));
            REQUIRE(r.has_value());
            CHECK(MeetsNear(*r, ta, tb, swapped) == meets);
            for (const auto& p : r->Points) {
                CHECK(IntersectTest::ValidCurveWitness(p, swapped ? b : a, swapped ? db : da,
                                                       swapped ? a : b, swapped ? da : db,
                                                       IntersectTest::ctx));
            }
        }
    }
}

TEST_CASE("curve corpus witness checks reject corrupt positions and parameters",
          "[intersect][curve-curve][corpus][regression]")
{
    using namespace Stereon;
    const Curve a{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve b{Line{Point3{}, Vec3{0.0, 1.0, 0.0}}};
    const Interval domain{-1.0, 1.0};
    CurveCurvePoint p{0.0, 0.0, Point3{}};
    CHECK(IntersectTest::ValidCurveWitness(p, a, domain, b, domain, IntersectTest::ctx));
    p.Position.X = 1.0;
    CHECK_FALSE(IntersectTest::ValidCurveWitness(p, a, domain, b, domain, IntersectTest::ctx));
    p.Position = Point3{};
    p.ParameterA = 2.0;
    CHECK_FALSE(IntersectTest::ValidCurveWitness(p, a, domain, b, domain, IntersectTest::ctx));
    p.ParameterA = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(IntersectTest::ValidCurveWitness(p, a, domain, b, domain, IntersectTest::ctx));
    p.ParameterA = 0.0;
    p.Position.X = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(IntersectTest::ValidCurveWitness(p, a, domain, b, domain, IntersectTest::ctx));
}
