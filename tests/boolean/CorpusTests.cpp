// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Boolean corpus (corpus/boolean/README.md): each case's operands are
// made from their parameters, the operation is run on the serial scheduler,
// and the result's solid and face counts and its validity are compared with
// the case's closed-form expectation. Cases of grids larger than 3 by 3 are
// issue #160's full-size reproducers; they are hidden from the routine CTest
// matrix and run in the "Boolean issue 160 large grids" workflow.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "Fixtures.hpp"

using BooleanTest::CountOf;
using Stereon::BooleanOperation;
using Stereon::EntityKind;
using Stereon::Shape;
using StereonTest::JsonValue;

namespace {

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR) / "boolean";

/// The largest grid the routine matrix runs.
constexpr int RoutineGrid = 3;

/// Every case file of the suite, in path order.
std::vector<std::filesystem::path> CaseFiles()
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(corpus)) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

/// The largest grid among the case's operands.
int LargestGrid(const JsonValue& c)
{
    int largest = 0;
    for (const JsonValue& operand : c["inputs"]["operands"].AsArray()) {
        largest = std::max(largest, static_cast<int>(operand["n"].AsNumber()));
    }
    return largest;
}

Shape OperandOf(const JsonValue& operand)
{
    REQUIRE(operand["type"].AsString() == "drill-grid");
    const JsonValue& offset = operand["offset"];
    return BooleanCases::Drills(
        static_cast<int>(operand["n"].AsNumber()),
        Stereon::Vec3{offset[0].AsNumber(), offset[1].AsNumber(), offset[2].AsNumber()});
}

BooleanOperation OperationOf(const std::string& name)
{
    if (name == "common") {
        return BooleanOperation::Common;
    }
    if (name == "cut") {
        return BooleanOperation::Cut;
    }
    REQUIRE(name == "union");
    return BooleanOperation::Union;
}

/// The case's result: the operand itself for "build", else the Boolean.
Shape ResultOf(const JsonValue& c)
{
    std::map<std::string, Shape> operands;
    for (const JsonValue& operand : c["inputs"]["operands"].AsArray()) {
        operands.emplace(operand["name"].AsString(), OperandOf(operand));
    }
    const JsonValue& operation = c["operation"];
    const JsonValue& names = operation["operands"];
    const Shape& a = operands.at(names[0].AsString());
    if (operation["name"].AsString() == "build") {
        return a;
    }
    const Shape& b = operands.at(names[1].AsString());
    auto result = Boolean(a, b, OperationOf(operation["name"].AsString()), BooleanCases::serial);
    INFO((result ? "" : std::string(result.error().GetMessage())));
    REQUIRE(result.has_value());
    return std::move(result->Shape);
}

void Check(const JsonValue& c)
{
    const Shape result = ResultOf(c);
    const JsonValue& expected = c["expected"];
    CHECK(CountOf<EntityKind::Solid>(result) ==
          static_cast<std::size_t>(expected["solids"].AsNumber()));
    if (expected.Has("faces")) {
        CHECK(CountOf<EntityKind::Face>(result) ==
              static_cast<std::size_t>(expected["faces"].AsNumber()));
    }
    REQUIRE(expected["valid"].AsBool());
    CHECK(ValidateTopology(result).IsValid());
    // At the solid level every shell must close, be oriented outwards and
    // nest; the one finding allowed is the notice that self-intersection
    // was not checked, which every solid-level report holds.
    const Stereon::ValidityReport solid =
        ValidateTopology(result, Stereon::ValidityLevel::Solid, BooleanCases::serial);
    for (const Stereon::ValidityFinding& finding : solid.GetFindings()) {
        INFO(finding.Message);
        CHECK(finding.Rule == Stereon::TopologyRule::SelfIntersectionUnchecked);
    }
}

/// Runs the cases whose largest grid is in [lo, hi]; returns how many.
int RunCases(const int lo, const int hi)
{
    int ran = 0;
    for (const std::filesystem::path& file : CaseFiles()) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        const int grid = LargestGrid(c);
        if (grid < lo || grid > hi) {
            continue;
        }
        INFO(c["id"].AsString());
        Check(c);
        ++ran;
    }
    return ran;
}

} // namespace

TEST_CASE("the Boolean corpus's small cases hold", "[boolean][corpus][issue160]")
{
    CHECK(RunCases(1, RoutineGrid) > 0);
}

TEST_CASE("the Boolean corpus's full-size cases hold",
          "[.][boolean][corpus][regression][issue160][large]")
{
    CHECK(RunCases(RoutineGrid + 1, 1 << 20) > 0);
}
