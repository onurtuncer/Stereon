// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cmath>
#include <numbers>

module Stereon.Core;

import :Frame;
import :Trigonometry;

namespace Stereon {

double Sin(const double x) noexcept
{
    return std::sin(x);
}

double Cos(const double x) noexcept
{
    return std::cos(x);
}

double Atan2(const double y, const double x) noexcept
{
    return std::atan2(y, x);
}

double Sqrt(const double x) noexcept
{
    return std::sqrt(x);
}

double Hypot(const double x, const double y) noexcept
{
    return std::hypot(x, y);
}

double WrappedAngle(const double y, const double x) noexcept
{
    constexpr double period = 2.0 * std::numbers::pi;
    if (x == 0.0 && y == 0.0) {
        return 0.0;
    }
    const double angle = std::atan2(y, x);
    if (angle >= 0.0) {
        return angle;
    }
    const double wrapped = angle + period;
    return wrapped < period ? wrapped : 0.0;
}

} // namespace Stereon
