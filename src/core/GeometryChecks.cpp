// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cmath>
#include <expected>
#include <initializer_list>
#include <string>
#include <string_view>
#include <utility>

module Stereon.Core;

import :Error;
import :GeometryChecks;
import :Vector;

namespace Stereon::Detail {

namespace {

/// "<where>: <what> <problem>", concatenated rather than formatted (see
/// EvaluatorChecks.cpp).
[[nodiscard]] std::unexpected<KernelError>
Invalid(const std::string_view where, const std::string_view what, const std::string_view problem)
{
    std::string message{where};
    message += ": ";
    message += what;
    message += ' ';
    message += problem;
    return Fail(ErrorCode::InvalidInput, std::move(message));
}

[[nodiscard]] bool AllFinite(const std::initializer_list<double> values)
{
    for (const double x : values) {
        if (!std::isfinite(x)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool WithinFrameAllowance(const double x)
{
    return std::abs(x) <= FrameAllowance;
}

} // namespace

Result<void> CheckFinite(const double x, const std::string_view what, const std::string_view where)
{
    if (!std::isfinite(x)) {
        return Invalid(where, what, "must be finite");
    }
    return {};
}

Result<void> CheckFinite(const Point2& p, const std::string_view what, const std::string_view where)
{
    if (!AllFinite({p.X, p.Y})) {
        return Invalid(where, what, "must be finite");
    }
    return {};
}

Result<void> CheckFinite(const Point3& p, const std::string_view what, const std::string_view where)
{
    if (!AllFinite({p.X, p.Y, p.Z})) {
        return Invalid(where, what, "must be finite");
    }
    return {};
}

Result<void> CheckFinite(const Vec2& v, const std::string_view what, const std::string_view where)
{
    if (!AllFinite({v.X, v.Y})) {
        return Invalid(where, what, "must be finite");
    }
    return {};
}

Result<void> CheckFinite(const Vec3& v, const std::string_view what, const std::string_view where)
{
    if (!AllFinite({v.X, v.Y, v.Z})) {
        return Invalid(where, what, "must be finite");
    }
    return {};
}

Result<void> CheckPositive(const double x, const std::string_view what,
                           const std::string_view where)
{
    if (!std::isfinite(x) || !(x > 0.0)) {
        return Invalid(where, what, "must be finite and positive");
    }
    return {};
}

Result<void> CheckNonNegative(const double x, const std::string_view what,
                              const std::string_view where)
{
    if (!std::isfinite(x) || x < 0.0) {
        return Invalid(where, what, "must be finite and not negative");
    }
    return {};
}

Result<void> CheckOrthonormal(const Vec2& x, const Vec2& y, const std::string_view where)
{
    if (!AllFinite({x.X, x.Y, y.X, y.Y})) {
        return Invalid(where, "the frame axes", "must be finite");
    }
    if (!WithinFrameAllowance(Dot(x, x) - 1.0) || !WithinFrameAllowance(Dot(y, y) - 1.0) ||
        !WithinFrameAllowance(Dot(x, y))) {
        return Invalid(where, "XAxis and YAxis", "must be orthonormal");
    }
    return {};
}

Result<void> CheckOrthonormal(const Vec3& x, const Vec3& y, const std::string_view where)
{
    if (!AllFinite({x.X, x.Y, x.Z, y.X, y.Y, y.Z})) {
        return Invalid(where, "the frame axes", "must be finite");
    }
    if (!WithinFrameAllowance(Dot(x, x) - 1.0) || !WithinFrameAllowance(Dot(y, y) - 1.0) ||
        !WithinFrameAllowance(Dot(x, y))) {
        return Invalid(where, "XAxis and YAxis", "must be orthonormal");
    }
    return {};
}

} // namespace Stereon::Detail
