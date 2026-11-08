// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The analytic handlers on corpus/analytic-ssi/ (RS-0001, the Boolean
// prototype's M3): every pair of planes, spheres and cylinders, the
// hand-written seeds and the generated sweeps (tools/generate_corpus.py),
// decided as its nominal configuration, the exact one of the input
// doubles, except that exactly parallel supports within the pair tolerance
// are declared coincident (ADR-0029 section 6, decision D1 of the plan):
// the corpus's one tolerant claim that nominal geometry honours. The pair
// tolerance is the larger of the two (ADR-0001, decision D2). Cylinder
// pairs with skew axes, Tier B, are decided as the configuration the
// corpus names for them (M3b), read off the result's branches: how many
// closed curves they form and at how many nodes they meet.

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
/// direction, a cylinder its axis's; neither keeps the vector's length. A
/// sphere's frame is the coordinate frame.
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
    return Surface{Stereon::Sphere{PointOf(s["centre"]), s["radius"].AsNumber(),
                                   Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
}

/// The kind of curve a branch is, in the corpus's words; a traced branch,
/// a fitted curve, is a "loop", which names any closed conic.
std::string CurveOf(const Stereon::IntersectionBranch& branch)
{
    if (branch.Curve3d.Target<Stereon::Circle>() != nullptr) {
        return "circle";
    }
    if (branch.Curve3d.Target<Stereon::Ellipse>() != nullptr) {
        return "ellipse";
    }
    return branch.Curve3d.Target<Stereon::Line>() != nullptr ? "line" : "loop";
}

/// Whether the configuration found is the one wanted: the same, or the
/// same without the side of a tangency, or a traced loop where the corpus
/// names the closed conic it is (spheres off their poles' axis go to the
/// trace).
bool Matches(const std::string& want, const std::string& got)
{
    if (got == "loop") {
        return want == "circle" || want == "ellipse";
    }
    if (got == "two_loops") {
        return want == "two_circles" || want == "two_ellipses";
    }
    return want.starts_with(got);
}

/// The configuration a result stands for, in the corpus's words; a tangency
/// is reported without its external/internal side, which the expected name
/// adds. Branches cut at a seam count once per curve: a circle is one
/// configuration however many arcs carry it.
std::string ConfigurationOf(const Stereon::SurfaceSurfaceIntersection& r)
{
    if (!r.Overlaps.empty()) {
        return "coincident";
    }
    if (r.Branches.empty()) {
        return r.Points.size() == 1 ? "tangent_point" : r.Points.empty() ? "disjoint" : "points";
    }
    const Stereon::IntersectionBranch& branch = r.Branches.front();
    const std::string curve = CurveOf(branch);
    if (branch.Contact == ContactType::Tangential) {
        return "tangent_" + curve;
    }
    const std::size_t curves =
        curve == "line" ? r.Branches.size() : IntersectTest::ConnectivityOf(r).Curves;
    return curves == 2 ? "two_" + curve + "s" : curve;
}

/// The configuration of two cylinders with skew axes, in the corpus's words:
/// an isolated point is a tangency, loops are counted as the closed curves
/// the branches form, and nodes, where four branch ends meet, mark an
/// internal tangency (one) or two ellipses (two).
std::string SkewConfigurationOf(const Stereon::SurfaceSurfaceIntersection& r)
{
    if (r.IsEmpty()) {
        return "disjoint";
    }
    if (r.Points.size() == 1 && r.Branches.empty()) {
        return "tangent_point";
    }
    const IntersectTest::Connectivity shape = IntersectTest::ConnectivityOf(r);
    if (shape.Nodes == 2) {
        return "two_ellipses";
    }
    if (shape.Nodes == 1) {
        return "internal_tangent";
    }
    if (shape.Curves == 1) {
        return "one_loop";
    }
    return shape.Curves == 2 ? "two_loops" : "unexpected";
}

/// Whether the case is two cylinders whose axes are not parallel.
bool IsSkewCylinders(const JsonValue& surfaces)
{
    if (surfaces[0]["type"].AsString() != "cylinder" ||
        surfaces[1]["type"].AsString() != "cylinder") {
        return false;
    }
    const Vec3 cross = Cross(VectorOf(surfaces[0]["axis"]), VectorOf(surfaces[1]["axis"]));
    return cross != Vec3{};
}

bool HasSphere(const JsonValue& surfaces)
{
    return std::ranges::any_of(surfaces.AsArray(),
                               [](const JsonValue& s) { return s["type"].AsString() == "sphere"; });
}

/// What the kernel must return for the case: the nominal configuration, or
/// the tolerant one where it is a declared coincidence. A pair with a
/// sphere goes to the elementary-pair handler, which also declares a
/// tangency where the profiles meet within the pair tolerance of one
/// (ADR-0029 section 6, ElementaryPairs.cpp): it is held to the tolerant
/// configuration, the one RS-0001's contract proposes, whenever the case
/// gives one.
std::string ExpectedOf(const JsonValue& expected, const JsonValue& surfaces)
{
    const JsonValue& tolerant = expected["tolerant"];
    if (!tolerant.Has("configuration")) {
        return expected["nominal"]["configuration"].AsString();
    }
    const std::string& configuration = tolerant["configuration"].AsString();
    if (configuration == "coincident" || HasSphere(surfaces)) {
        return configuration;
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
    if (expected.Has("kernel")) {
        // The kernel refuses the case, for the reason the case gives.
        REQUIRE_FALSE(result.has_value());
        REQUIRE(EnumName(result.error().GetCode()) == expected["kernel"]["outcome"].AsString());
        return;
    }
    INFO(IntersectTest::ErrorText(result));
    CHECK(result.has_value());
    if (!result) {
        return;
    }
    if (expected["nominal"]["tier"].AsString() == "B") {
        REQUIRE(IsSkewCylinders(surfaces));
        const std::string want = expected["nominal"]["configuration"].AsString();
        const std::string got = SkewConfigurationOf(*result);
        INFO("expected " << want << ", got " << got);
        CHECK(got == want);
        return;
    }
    const std::string want = ExpectedOf(expected, surfaces);
    const std::string got = ConfigurationOf(*result);
    INFO("expected " << want << ", got " << got);
    CHECK(Matches(want, got));
    if (got == "coincident") {
        const bool nominal = expected["nominal"]["configuration"].AsString() == "coincident";
        REQUIRE((result->Overlaps.front().Evidence == Stereon::CoincidenceEvidence::Exact) ==
                nominal);
    }
}

/// Every family of the suite, the sweeps included.
std::vector<std::string> Families()
{
    std::vector<std::string> families;
    for (const auto& entry : std::filesystem::directory_iterator(corpus)) {
        if (entry.is_directory()) {
            families.push_back(entry.path().filename().string());
        }
    }
    std::ranges::sort(families);
    return families;
}

} // namespace

TEST_CASE("every pair case of the corpus is decided as required", "[intersect][corpus]")
{
    std::size_t count = 0;
    for (const std::string& family : Families()) {
        for (const std::filesystem::path& file : CasesIn(family)) {
            const JsonValue c = StereonTest::ReadJsonFile(file);
            const JsonValue& surfaces = c["inputs"]["surfaces"];
            // Tier B pairs other than skew cylinders (a sphere off a
            // cylinder's axis) are decided by the trace, which the corpus
            // names no configuration for yet.
            const bool tierB = c["expected"]["nominal"].Has("tier") &&
                               c["expected"]["nominal"]["tier"].AsString() == "B";
            if (c["operation"]["name"].AsString() != "intersect" ||
                (tierB && !IsSkewCylinders(surfaces))) {
                continue;
            }
            CheckCase(file);
            ++count;
        }
    }
    // The 1024 generated sweeps and the hand-written pairs.
    REQUIRE(count >= 1050);
}
