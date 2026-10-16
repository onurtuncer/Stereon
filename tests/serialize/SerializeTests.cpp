// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <vector>

import Stereon.Surface;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Cylinder;
using Stereon::Deserialized;
using Stereon::ErrorCode;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Serialized;
using Stereon::ToJson;
using Stereon::Vec3;

namespace {

template <typename T>
T RoundTrip(const T& value)
{
    const std::vector<std::byte> bytes = Serialized(value);
    const auto back = Deserialized<T>(std::span<const std::byte>{bytes});
    REQUIRE(back.has_value());
    return *back;
}

} // namespace

TEST_CASE("geometry structs serialise by reflection and read back exactly", "[serialize]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        REQUIRE(Deserialized<Line>(std::span<const std::byte>{}).error().GetCode() ==
                ErrorCode::NotSupported);
        return;
    }
    const Line line{Point3{1.0, -2.0, 3.5}, Vec3{0.1, 0.2, 0.3}};
    const Circle circle{Point3{0.5, 0.25, -1.0}, 2.75, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Plane plane{Point3{9.0, 8.0, 7.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const Cylinder cylinder{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, false};

    REQUIRE(RoundTrip(line) == line);
    REQUIRE(RoundTrip(circle) == circle);
    REQUIRE(RoundTrip(plane) == plane);
    REQUIRE(RoundTrip(cylinder) == cylinder);

    // Declaration order, 8 bytes per double, 1 per bool.
    REQUIRE(Serialized(line).size() == 6 * 8);
    REQUIRE(Serialized(circle).size() == 10 * 8);
    REQUIRE(Serialized(cylinder).size() == 10 * 8 + 1);

    // Little-endian: 1.0 is 0x3FF0000000000000.
    const std::vector<std::byte> one = Serialized(1.0);
    REQUIRE(one.size() == 8);
    REQUIRE(one[0] == std::byte{0x00});
    REQUIRE(one[6] == std::byte{0xF0});
    REQUIRE(one[7] == std::byte{0x3F});
}

TEST_CASE("malformed input is a parse error, not a crash", "[serialize]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const Cylinder cylinder{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    std::vector<std::byte> bytes = Serialized(cylinder);

    std::vector<std::byte> truncated(bytes.begin(), bytes.end() - 2);
    const auto shortRead = Deserialized<Cylinder>(std::span<const std::byte>{truncated});
    REQUIRE(shortRead.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(shortRead.error().GetMessage()).starts_with("truncated input"));

    std::vector<std::byte> overlong = bytes;
    overlong.push_back(std::byte{0});
    REQUIRE(Deserialized<Cylinder>(std::span<const std::byte>{overlong}).error().GetCode() ==
            ErrorCode::ParseError);

    bytes.back() = std::byte{2}; // the RightHanded flag
    const auto badBool = Deserialized<Cylinder>(std::span<const std::byte>{bytes});
    REQUIRE(std::string(badBool.error().GetMessage()).starts_with("invalid bool byte 2"));
}

TEST_CASE("serialisation round-trips every double bit pattern", "[serialize][property]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    rc::prop("any finite or non-finite coordinates", [] {
        const auto coordinate = [] { return *rc::gen::arbitrary<double>(); };
        const Line line{Point3{coordinate(), coordinate(), coordinate()},
                        Vec3{coordinate(), coordinate(), coordinate()}};
        const std::vector<std::byte> bytes = Serialized(line);
        const auto back = Deserialized<Line>(std::span<const std::byte>{bytes});
        RC_ASSERT(back.has_value());
        RC_ASSERT(Serialized(*back) == bytes); // bit for bit, NaN payloads included
    });
}

TEST_CASE("the JSON view names every member", "[serialize]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        REQUIRE(ToJson(Line{}).empty());
        return;
    }
    const Line line{Point3{1.0, 2.0, 3.0}, Vec3{0.5, 0.0, -1.0}};
    REQUIRE(ToJson(line) == R"({"Origin":{"X":1,"Y":2,"Z":3},"Direction":{"X":0.5,"Y":0,"Z":-1}})");
    const Cylinder cylinder{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, false};
    REQUIRE(ToJson(cylinder).ends_with(
        R"("Radius":1.5,"XAxis":{"X":1,"Y":0,"Z":0},"YAxis":{"X":0,"Y":1,"Z":0},"RightHanded":false})"));
    const Line odd{Point3{std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::infinity(), 0.0},
                   Vec3{}};
    REQUIRE(ToJson(odd).starts_with(R"({"Origin":{"X":"nan","Y":"inf","Z":0})"));
}

TEST_CASE("vectors serialise as a count and their elements", "[serialize]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const BSplineCurve curve{3,
                             {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
                             {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
                              Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}},
                             {1.0, 0.5, 2.0, 0.75, 1.0}};
    REQUIRE(RoundTrip(curve) == curve);
    // The degree (4 bytes), then each vector as 8 bytes of count and its
    // elements: 9 knots, 5 poles of 3 doubles, 5 weights.
    REQUIRE(Serialized(curve).size() == 4 + (8 + 9 * 8) + (8 + 5 * 24) + (8 + 5 * 8));

    const BSplineCurve polynomial{1, {0.0, 0.0, 1.0, 1.0}, {Point3{}, Point3{1.0, 1.0, 1.0}}, {}};
    REQUIRE(RoundTrip(polynomial) == polynomial); // an empty vector reads back empty

    const BSplineSurface surface{1,
                                 1,
                                 {0.0, 0.0, 1.0, 1.0},
                                 {0.0, 0.0, 1.0, 1.0},
                                 {Point3{0.0, 0.0, 0.0}, Point3{0.0, 1.0, 0.0},
                                  Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 1.0}},
                                 {}};
    REQUIRE(RoundTrip(surface) == surface);
    REQUIRE(
        ToJson(polynomial) ==
        R"({"Degree":1,"Knots":[0,0,1,1],"Poles":[{"X":0,"Y":0,"Z":0},{"X":1,"Y":1,"Z":1}],"Weights":[]})");
}

TEST_CASE("a corrupt element count is a parse error, not a huge allocation", "[serialize]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const BSplineCurve curve{1, {0.0, 0.0, 1.0, 1.0}, {Point3{}, Point3{1.0, 1.0, 1.0}}, {}};
    std::vector<std::byte> bytes = Serialized(curve);
    // The knot count follows the 4-byte degree; claim 2^62 knots.
    bytes[4 + 7] = std::byte{0x40};
    const auto back = Deserialized<BSplineCurve>(std::span<const std::byte>{bytes});
    REQUIRE(back.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(back.error().GetMessage()).starts_with("element count"));
}
