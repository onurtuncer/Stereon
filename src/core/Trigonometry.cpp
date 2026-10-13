// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cmath>

module Stereon.Core;

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

} // namespace Stereon
