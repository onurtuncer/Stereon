// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The analytic handlers on corpus/analytic-ssi/ (RS-0001, the Boolean
// prototype's M3): every plane–plane, plane–cylinder and cylinder–cylinder
// case decided as its nominal configuration, the exact one of the input
// doubles, except that exactly parallel supports within the pair tolerance
// are declared coincident (ADR-0029 section 6, decision D1 of the plan):
// the corpus's one tolerant claim that nominal geometry honours. The pair
// tolerance is the larger of the two (ADR-0001, decision D2). Tier B
// pairs, the quartics, report NotImplemented until M3b.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::ContactType;
using Stereon::ErrorCode;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR) / "analytic-ssi";

Vec3 VectorOf(const JsonValue& v)
{
    return Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Point3 PointOf(const JsonValue& v)
{
    return Point3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

/// The corpus surface as kernel geometry: a plane keeps its normal's
/// direction, a cylinder its axis's; neither keeps the vector's length.
Surface SurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    if (type == "plane") {
        return Surface{IntersectTest::PlaneWithNormal(PointOf(s["point"]), VectorOf(s["normal"]))};
    }
    if (type == "cylinder") {
        return Surface{IntersectTest::CylinderAbout(PointOf(s["point"]), VectorOf(s["axis"]),
                                                    s["radius"].AsNumber())};
    }
    return Surface{}; // Spheres wait for their handlers.
}

/// The configuration a result stands for, in the corpus's words; a tangent
/// ruling of two cylinders is reported without its external/internal side.
std::string ConfigurationOf(const Stereon::SurfaceSurfaceIntersection& r)
{
    if (!r.Overlaps.empty()) {
        return "coincident";
    }
    if (r.Branches.empty()) {
        return "disjoint";
    }
    if (r.Branches.size() == 2) {
        return "two_lines";
    }
    const Stereon::IntersectionBranch& branch = r.Branches.front();
    if (branch.Contact == ContactType::Tangential) {
        return "tangent_line";
    }
    if (branch.Curve3d.Target<Stereon::Circle>() != nullptr) {
        return "circle";
    }
    if (branch.Curve3d.Target<Stereon::Ellipse>() != nullptr) {
        return "ellipse";
    }
    return "line";
}

/// What the kernel must return for the case: the nominal configuration, or
/// the tolerant one where it is a declared coincidence.
std::string ExpectedOf(const JsonValue& expected)
{
    const JsonValue& tolerant = expected["tolerant"];
    if (tolerant.Has("configuration") && tolerant["configuration"].AsString() == "coincident") {
        return "coincident";
    }
    return expected["nominal"]["configuration"].AsString();
}

std::vector<std::filesystem::path> CasesIn(const std::string_view family)
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(corpus / family)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

void CheckCase(const std::filesystem::path& file)
{
    const JsonValue c = StereonTest::ReadJsonFile(file);
    INFO(c["id"].AsString());
    const JsonValue& surfaces = c["inputs"]["surfaces"];
    REQUIRE(surfaces.AsArray().size() == 2);
    const double tolerance =
        std::max(surfaces[0]["tolerance"].AsNumber(), surfaces[1]["tolerance"].AsNumber());
    const auto result =
        Intersect(SurfaceOf(surfaces[0]), Box2::Everything(), SurfaceOf(surfaces[1]),
                  Box2::Everything(), tolerance, IntersectTest::ctx);
    const JsonValue& expected = c["expected"];
    if (expected["nominal"]["tier"].AsString() == "B") {
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::NotImplemented);
        return;
    }
    REQUIRE(result.has_value());
    const std::string want = ExpectedOf(expected);
    const std::string got = ConfigurationOf(*result);
    INFO("expected " << want << ", got " << got);
    REQUIRE(want.starts_with(got));
    if (got == "coincident") {
        const bool nominal = expected["nominal"]["configuration"].AsString() == "coincident";
        REQUIRE((result->Overlaps.front().Evidence == Stereon::CoincidenceEvidence::Exact) ==
                nominal);
    }
}

} // namespace

TEST_CASE("every Tier A plane and cylinder case of the corpus is decided as required",
          "[intersect][corpus]")
{
    std::size_t count = 0;
    for (const std::string_view family :
         {"plane-plane", "plane-cylinder", "cylinder-cylinder", "tier-b"}) {
        for (const std::filesystem::path& file : CasesIn(family)) {
            const JsonValue c = StereonTest::ReadJsonFile(file);
            const JsonValue& surfaces = c["inputs"]["surfaces"];
            const bool planesAndCylinders =
                std::ranges::all_of(surfaces.AsArray(), [](const JsonValue& s) {
                    return s["type"].AsString() == "plane" || s["type"].AsString() == "cylinder";
                });
            if (!planesAndCylinders || c["operation"]["name"].AsString() != "intersect") {
                continue;
            }
            CheckCase(file);
            ++count;
        }
    }
    REQUIRE(count >= 12);
}
