// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <concepts>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>

export module Stereon.Robust:Filter;

import Stereon.Core;
import :Sign;

/// Interval filters with exact fallback (ADR-0002). The adaptive predicates
/// carry their own static filters; predicates on curves and surfaces, whose
/// error bounds are not known in closed form, evaluate their quantity in
/// Interval arithmetic instead. An interval that excludes zero decides the
/// sign; one that contains it hands over to an exact or subdivided
/// evaluation.

export namespace Stereon {

/// The sign shared by every number in the interval, or nothing when the
/// interval contains numbers of different signs or is empty. The degenerate
/// interval [0, 0] has sign Zero.
[[nodiscard]] constexpr std::optional<Sign> SignOf(const Interval& i) noexcept
{
    if (i.IsEmpty()) {
        return std::nullopt;
    }
    if (i.GetLower() > 0.0) {
        return Sign::Positive;
    }
    if (i.GetUpper() < 0.0) {
        return Sign::Negative;
    }
    if (i.GetLower() == 0.0 && i.GetUpper() == 0.0) {
        return Sign::Zero;
    }
    return std::nullopt;
}

/// The sign of `estimate` when the interval decides it, otherwise the sign
/// `exact` computes. `exact` is only called when needed.
template <std::invocable Exact>
    requires std::convertible_to<std::invoke_result_t<Exact>, Sign>
[[nodiscard]] constexpr Sign FilteredSign(const Interval& estimate, Exact&& exact)
{
    if (const std::optional<Sign> sign = SignOf(estimate)) {
        return *sign;
    }
    return std::invoke(std::forward<Exact>(exact));
}

} // namespace Stereon
