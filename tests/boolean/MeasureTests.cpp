// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The volumes and areas of the Boolean fixtures' results
// (tests/support/BooleanCases.hpp) against their exact values: the measure
// the OCCT differential tests use (tests/support/ShapeMeasures.hpp) checked
// where the answer is known in closed form, and the Booleans checked by it.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "ShapeMeasures.hpp"

using Stereon::EntityKind;

namespace {

constexpr double pi = std::numbers::pi;

/// A case's exact volume and area, where the result has them.
struct Exact {
    const char* Name;
    double Volume;
    double Area;
};

/// Two boxes of side 10 overlapping in a cube of side 5; a box and a drum of
/// radius 3 through it, 2 below and 3 above; the slab of 3 x 3 cells drilled
/// with holes of radius 2.
const Exact exact[] = {
    {"corner Union", 2000.0 - 125.0, 1200.0 - 150.0},
    {"corner Common", 125.0, 150.0},
    {"corner Cut", 875.0, 600.0},
    {"drum Union", 1000.0 + 45.0 * pi, 600.0 + 30.0 * pi},
    {"drum Common", 90.0 * pi, 78.0 * pi},
    {"drum Cut", 1000.0 - 90.0 * pi, 600.0 + 42.0 * pi},
    {"loose Cut", 1000.0, 600.0},
    {"drilled slab Cut", 9000.0 - 360.0 * pi, 3000.0 + 288.0 * pi},
};

double LargestTolerance(const Stereon::Shape& shape)
{
    double largest = 0.0;
    for (const Stereon::VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        largest = std::max(largest, ToleranceOf(shape, v));
    }
    for (const Stereon::EdgeId e : *AllEntities<EntityKind::Edge>(shape)) {
        largest = std::max(largest, ToleranceOf(shape, e));
    }
    return largest;
}

} // namespace

TEST_CASE("the measure of a Boolean result is its exact volume and area", "[boolean][measure]")
{
    std::size_t checked = 0;
    for (const BooleanCases::Case& c : BooleanCases::Cases()) {
        const auto* known =
            std::ranges::find_if(exact, [&](const Exact& e) { return c.Name == e.Name; });
        if (known == std::end(exact)) {
            continue;
        }
        INFO(c.Name);
        const auto result = Boolean(c.A, c.B, c.Operation, BooleanCases::serial);
        REQUIRE(result.has_value());
        const ShapeMeasures::Measures m = ShapeMeasures::MeasuresOf(result->Shape);
        INFO("volume " << m.Volume << ", area " << m.Area);
        CHECK(m.VolumeEstimate <= 1e-10 * known->Volume);
        CHECK(m.AreaEstimate <= 1e-10 * known->Area);
        CHECK(std::abs(m.Volume - known->Volume) <= 1e-10 * known->Volume);
        CHECK(std::abs(m.Area - known->Area) <= 1e-10 * known->Area);
        ++checked;
    }
    CHECK(checked == std::size(exact));
}

TEST_CASE("a box loosened 0.8 of its tolerance from another is measured within that tolerance",
          "[boolean][measure]")
{
    // The gap closes within the tolerance only, so the union's volume is
    // known to the faces' area times the tolerance: between the two boxes'
    // 1999.2 and the closed 2000. Their common is empty.
    for (const BooleanCases::Case& c : BooleanCases::Cases()) {
        if (c.Name != "loose Union" && c.Name != "loose Common") {
            continue;
        }
        INFO(c.Name);
        const auto result = Boolean(c.A, c.B, c.Operation, BooleanCases::serial);
        REQUIRE(result.has_value());
        const ShapeMeasures::Measures m = ShapeMeasures::MeasuresOf(result->Shape);
        INFO("volume " << m.Volume << ", area " << m.Area);
        if (c.Operation == Stereon::BooleanOperation::Common) {
            CHECK(AllEntities<EntityKind::Face>(result->Shape)->empty());
            CHECK(m.Volume == 0.0);
            continue;
        }
        const double band = m.Area * LargestTolerance(result->Shape);
        CHECK(band <= 1.0 + 1e-9);
        CHECK(m.Volume >= 1999.2 - band);
        CHECK(m.Volume <= 2000.0 + band);
        CHECK(std::abs(m.Area - 1000.0) <= m.EdgeLength * LargestTolerance(result->Shape));
    }
}

TEST_CASE("the measure of the primitives", "[boolean][measure]")
{
    const Stereon::Shape box =
        BooleanCases::BoxFrom(Stereon::Point3{}, Stereon::Point3{1.0, 2.0, 3.0});
    const ShapeMeasures::Measures b = ShapeMeasures::MeasuresOf(box);
    CHECK(std::abs(b.Volume - 6.0) <= 1e-13);
    CHECK(std::abs(b.Area - 22.0) <= 1e-13);
    CHECK(std::abs(b.EdgeLength - 24.0) <= 1e-13);

    const Stereon::Shape cylinder =
        BooleanCases::Must(
            MakeCylinderSolid(Stereon::Point3{1.0, -2.0, 0.5}, Stereon::Vec3{0.0, 1.0, 0.0},
                              Stereon::Vec3{0.0, 0.0, 1.0}, 1.5, 4.0, BooleanCases::serial),
            "a cylinder")
            .Shape;
    const ShapeMeasures::Measures c = ShapeMeasures::MeasuresOf(cylinder);
    CHECK(std::abs(c.Volume - pi * 1.5 * 1.5 * 4.0) <= 1e-12);
    CHECK(std::abs(c.Area - (2.0 * pi * 1.5 * 1.5 + 2.0 * pi * 1.5 * 4.0)) <= 1e-12);
    // Two circles and the seam.
    CHECK(std::abs(c.EdgeLength - (2.0 * 2.0 * pi * 1.5 + 4.0)) <= 1e-12);
}
