// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0010's SSI corpus, the corpus of its >= 500-case gate: the surface pairs
// with a closed-form answer that the kernel is run on, in all three suites.
// analytic-ssi holds the plane, sphere and cylinder pairs, curved-ssi the
// pairs with a cone, a torus or an ellipsoid, and nurbs-ssi the pairs of
// B-spline surfaces, which the generic path answers: the path the 500-case
// gate is for (ADR-0010, decisions on acceptance of 2026-10-20) and the
// "near-coincident NURBS" and "pole contact" its Verification names.
//
// The gate runs as the tag [ssi-corpus]: this count, and the tests that hold
// every case of each suite to its closed form (CorpusTests.cpp, and the
// both-ways runs of SurfaceCorpusTests.cpp), run by the Intersection release
// evidence workflow. The count keeps the denominators apart: nurbs-ssi's
// may_fail_closed fixtures may be refused (ADR-0010: a case that exhausts
// the budget fails with the right error), so they are reported apart from
// the fixtures that must be answered.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

import Stereon.Intersect;

using StereonTest::JsonValue;

namespace {

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR);

std::vector<JsonValue> CasesOf(const std::string& suite)
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(corpus / suite)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    std::vector<JsonValue> cases;
    for (const std::filesystem::path& file : files) {
        cases.push_back(StereonTest::ReadJsonFile(file));
    }
    return cases;
}

Stereon::Vec3 VectorOf(const JsonValue& v)
{
    return Stereon::Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

/// Whether CorpusTests.cpp runs an analytic-ssi case against a closed form:
/// a pair, and of Tier B only the skew cylinders, whose configuration the
/// corpus names (a sphere off a cylinder's axis names none yet).
bool HasClosedForm(const JsonValue& c)
{
    if (c["operation"]["name"].AsString() != "intersect") {
        return false;
    }
    const JsonValue& nominal = c["expected"]["nominal"];
    if (!nominal.Has("tier") || nominal["tier"].AsString() != "B") {
        return true;
    }
    const JsonValue& s = c["inputs"]["surfaces"];
    return s[0]["type"].AsString() == "cylinder" && s[1]["type"].AsString() == "cylinder" &&
           Cross(VectorOf(s[0]["axis"]), VectorOf(s[1]["axis"])) != Stereon::Vec3{};
}

} // namespace

TEST_CASE("ADR-0010's SSI corpus holds at least 500 closed-form pairs, B-spline pairs among them",
          "[intersect][corpus][ssi-corpus]")
{
    const auto analytic = std::ranges::count_if(CasesOf("analytic-ssi"), HasClosedForm);
    const auto curved = static_cast<std::ptrdiff_t>(CasesOf("curved-ssi").size());
    const std::vector<JsonValue> nurbs = CasesOf("nurbs-ssi");
    const auto refusable = std::ranges::count_if(
        nurbs, [](const JsonValue& c) { return c["expected"].Has("may_fail_closed"); });
    const auto spline = static_cast<std::ptrdiff_t>(nurbs.size());
    WARN("SSI corpus: " << analytic << " analytic-ssi, " << curved << " curved-ssi, " << spline
                        << " nurbs-ssi (" << spline - refusable << " to be answered, " << refusable
                        << " may fail closed); " << analytic + curved + spline << " in all");
    // The suites' own floors (CorpusTests.cpp, SurfaceCorpusTests.cpp).
    CHECK(analytic >= 1050);
    CHECK(curved >= 600);
    CHECK(spline >= 90);
    CHECK(spline - refusable >= 58);
    REQUIRE(analytic + curved + spline >= 500);
}
