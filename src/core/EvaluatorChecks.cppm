// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <string_view>
#include <type_traits>

export module Stereon.Core:EvaluatorChecks;

import :Interval;
import :Scalar;
import :Vector;

/// The preconditions of the geometry evaluators (ADR-0003, ADR-0016):
/// finite parameters and points, parameters inside the domain of a curve or
/// surface that is not periodic in them, and derivative orders that exist.
///
/// The evaluators are templates, and GCC 16.2 drops contracts in templates
/// imported from a module (docs/toolchain/gcc16-feature-probe.md), so these
/// are plain branches that report through EvaluatorViolation(), a function
/// defined in a core implementation unit, which fires in every importer.
/// They are hot-path checks: compiled when STEREON_HOT_CONTRACTS is on (the
/// debug and sanitizer presets) and compiled out in release builds, where the
/// evaluators run in inner loops. Checks apply to double arguments; active
/// scalar types (ADR-0021) are checked through their primal callers.

export namespace Stereon::Detail {

/// True when the evaluator preconditions are checked in this build.
inline constexpr bool EvaluatorChecks = STEREON_HOT_CONTRACTS != 0;

/// Reports a broken evaluator precondition and raises a contract violation;
/// defined in EvaluatorChecks.cpp.
[[noreturn]] void EvaluatorViolation(std::string_view function, std::string_view problem) noexcept;

[[nodiscard]] constexpr bool IsFiniteValue(const double x) noexcept
{
    return x - x == 0.0; // false for infinities and NaN
}

constexpr void Require(const bool holds, const std::string_view function,
                       const std::string_view problem) noexcept
{
    if constexpr (EvaluatorChecks) {
        if (!holds) {
            EvaluatorViolation(function, problem);
        }
    }
}

template <Scalar T>
constexpr void RequireFinite(const T& x, const std::string_view function,
                             const std::string_view problem) noexcept
{
    if constexpr (EvaluatorChecks && std::is_same_v<T, double>) {
        Require(IsFiniteValue(x), function, problem);
    }
}

template <Scalar T>
constexpr void RequireFinite(const BasicPoint2<T>& p, const std::string_view function,
                             const std::string_view problem) noexcept
{
    RequireFinite(p.X, function, problem);
    RequireFinite(p.Y, function, problem);
}

template <Scalar T>
constexpr void RequireFinite(const BasicPoint3<T>& p, const std::string_view function,
                             const std::string_view problem) noexcept
{
    RequireFinite(p.X, function, problem);
    RequireFinite(p.Y, function, problem);
    RequireFinite(p.Z, function, problem);
}

/// t must lie in the domain unless the geometry is periodic in it (period
/// greater than zero), in which case every finite t is a valid parameter.
constexpr void RequireInDomain(const Interval& domain, const double t, const double period,
                               const std::string_view function) noexcept
{
    if constexpr (EvaluatorChecks) {
        Require(period > 0.0 || Contains(domain, t), function, "a parameter outside the domain");
    }
}

} // namespace Stereon::Detail
