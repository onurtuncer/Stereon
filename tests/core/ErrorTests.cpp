// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RapidCheck.hpp"
#include <Stereon/Core/Try.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <expected>
#include <memory>
#include <string>
#include <utility>
#include <vector>

import Stereon.Core;

using Stereon::ErrorCode;
using Stereon::Fail;
using Stereon::KernelError;
using Stereon::Result;

namespace {

Result<int> ParsePositive(const int value)
{
    if (value <= 0) {
        return Fail(ErrorCode::InvalidInput, "value must be positive");
    }
    return value;
}

// Counts how far a function got, to show that STEREON_TRY returns early.
struct Trace {
    int Steps = 0;
};

Result<int> SumOfPositives(const int a, const int b, Trace& trace)
{
    STEREON_TRY(const int first, ParsePositive(a));
    ++trace.Steps;
    STEREON_TRY(const int second, ParsePositive(b));
    ++trace.Steps;
    return first + second;
}

Result<void> CheckPositive(const int value)
{
    STEREON_TRY_VOID(ParsePositive(value));
    return {};
}

Result<std::unique_ptr<int>> MakeBoxed(const int value)
{
    STEREON_TRY(const int checked, ParsePositive(value));
    return std::make_unique<int>(checked);
}

Result<int> Unboxed(const int value)
{
    STEREON_TRY(std::unique_ptr<int> boxed, MakeBoxed(value));
    return *boxed;
}

// A different value type, to show that errors propagate across result types.
Result<std::string> Describe(const int value)
{
    STEREON_TRY(const int checked, ParsePositive(value));
    return std::to_string(checked);
}

} // namespace

TEST_CASE("error codes keep their numeric values", "[core][errors]")
{
    // The C ABI exposes these numbers (ADR-0006, ADR-0020). Codes are only
    // ever appended; 0 is left for success.
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::InvalidInput) == 1);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::InvalidResult) == 2);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::ToleranceExceeded) == 3);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::IntersectionFailed) == 4);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::CertificationFailed) == 5);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::DegenerateResult) == 6);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::NotSupported) == 7);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::NotImplemented) == 8);
    STATIC_REQUIRE(static_cast<std::uint32_t>(ErrorCode::ParseError) == 9);
}

TEST_CASE("every error code has a distinct name", "[core][errors]")
{
    REQUIRE(Stereon::ToString(ErrorCode::InvalidInput) == "InvalidInput");
    REQUIRE(Stereon::ToString(ErrorCode::ParseError) == "ParseError");

    std::vector<std::string> names;
    for (std::uint32_t value = 1; value <= 9; ++value) {
        const auto name = Stereon::ToString(static_cast<ErrorCode>(value));
        REQUIRE_FALSE(name.empty());
        REQUIRE(name != "Unknown");
        names.emplace_back(name);
    }
    std::ranges::sort(names);
    REQUIRE(std::ranges::adjacent_find(names) == names.end());

    REQUIRE(Stereon::ToString(static_cast<ErrorCode>(0)) == "Unknown");
}

TEST_CASE("KernelError carries code, message and tolerances", "[core][errors]")
{
    const KernelError plain(ErrorCode::NotImplemented, "boolean on tori");
    REQUIRE(plain.GetCode() == ErrorCode::NotImplemented);
    REQUIRE(plain.GetMessage() == "boolean on tori");
    REQUIRE(plain.GetTolerances().empty());

    const auto grown = KernelError(ErrorCode::ToleranceExceeded, "edge tolerance above ceiling")
                           .WithTolerance(0.02)
                           .WithTolerance(0.01);
    REQUIRE(grown.GetCode() == ErrorCode::ToleranceExceeded);
    REQUIRE(grown.GetTolerances().size() == 2);
    REQUIRE(grown.GetTolerances()[0] == 0.02);
    REQUIRE(grown.GetTolerances()[1] == 0.01);
}

TEST_CASE("Fail builds an unexpected KernelError", "[core][errors]")
{
    const Result<int> result = Fail(ErrorCode::ParseError, "unexpected token");
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(result.error().GetMessage() == "unexpected token");
}

TEST_CASE("STEREON_TRY unwraps values and returns errors early", "[core][errors][try]")
{
    Trace ok;
    const auto sum = SumOfPositives(2, 3, ok);
    REQUIRE(sum == 5);
    REQUIRE(ok.Steps == 2);

    Trace firstFails;
    const auto early = SumOfPositives(-1, 3, firstFails);
    REQUIRE_FALSE(early.has_value());
    REQUIRE(early.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(firstFails.Steps == 0);

    Trace secondFails;
    REQUIRE_FALSE(SumOfPositives(2, -3, secondFails).has_value());
    REQUIRE(secondFails.Steps == 1);
}

TEST_CASE("STEREON_TRY_VOID propagates errors without a value", "[core][errors][try]")
{
    REQUIRE(CheckPositive(4).has_value());

    const auto failed = CheckPositive(0);
    REQUIRE_FALSE(failed.has_value());
    REQUIRE(failed.error().GetMessage() == "value must be positive");
}

TEST_CASE("STEREON_TRY moves move-only values", "[core][errors][try]")
{
    REQUIRE(Unboxed(7) == 7);
    REQUIRE(Unboxed(-7).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("STEREON_TRY converts between result value types", "[core][errors][try]")
{
    REQUIRE(Describe(12) == "12");
    REQUIRE(Describe(-12).error().GetCode() == ErrorCode::InvalidInput);
}

TEST_CASE("errors survive propagation unchanged", "[core][errors][try][property]")
{
    rc::prop("code, message and tolerances pass through STEREON_TRY", [] {
        const auto code = static_cast<ErrorCode>(*rc::gen::inRange<std::uint32_t>(1, 10));
        const auto message = *rc::gen::arbitrary<std::string>();
        // Tolerances are positive lengths in millimetres (ADR-0001).
        const auto tolerances = *rc::gen::container<std::vector<double>>(
            rc::gen::map(rc::gen::inRange(1, 1'000'000), [](const int n) { return n * 1e-6; }));

        const auto source = [&]() -> Result<int> {
            KernelError error(code, message);
            for (const double tolerance : tolerances) {
                error = std::move(error).WithTolerance(tolerance);
            }
            return std::unexpected(std::move(error));
        };
        const auto forward = [&]() -> Result<std::string> {
            STEREON_TRY(const int value, source());
            return std::to_string(value);
        };

        const auto result = forward();
        RC_ASSERT(!result.has_value());
        RC_ASSERT(result.error().GetCode() == code);
        RC_ASSERT(result.error().GetMessage() == message);
        const auto carried = result.error().GetTolerances();
        RC_ASSERT(std::vector<double>(carried.begin(), carried.end()) == tolerances);
    });
}
