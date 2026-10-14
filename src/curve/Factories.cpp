// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Try.hpp>

#include <cstddef>
#include <string_view>

module Stereon.Curve;

/// The checked factories of the built-in curves (ADR-0006, ADR-0016), out of
/// line for the reason given in src/surface/Factories.cpp. Each dimension
/// has its own overload, so that callers need not name it; the checks are
/// written once, for both.

namespace Stereon {

namespace {

template <std::size_t D>
Result<BasicLine<double, D>> MakeLineN(const BasicPointN<double, D>& origin,
                                       const BasicVecN<double, D>& direction)
{
    constexpr std::string_view where = "MakeLine";
    STEREON_TRY_VOID(Detail::CheckFinite(origin, "Origin", where));
    STEREON_TRY_VOID(Detail::CheckFinite(direction, "Direction", where));
    if (!(LengthSquared(direction) > 0.0)) {
        return Fail(ErrorCode::InvalidInput, "MakeLine: Direction must not be zero");
    }
    return BasicLine<double, D>{origin, direction};
}

template <std::size_t D>
Result<BasicCircle<double, D>> MakeCircleN(const BasicPointN<double, D>& center,
                                           const double radius, const BasicVecN<double, D>& xAxis,
                                           const BasicVecN<double, D>& yAxis)
{
    constexpr std::string_view where = "MakeCircle";
    STEREON_TRY_VOID(Detail::CheckFinite(center, "Center", where));
    STEREON_TRY_VOID(Detail::CheckPositive(radius, "Radius", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return BasicCircle<double, D>{center, radius, xAxis, yAxis};
}

template <std::size_t D>
Result<BasicEllipse<double, D>> MakeEllipseN(const BasicPointN<double, D>& center,
                                             const double majorRadius, const double minorRadius,
                                             const BasicVecN<double, D>& xAxis,
                                             const BasicVecN<double, D>& yAxis)
{
    constexpr std::string_view where = "MakeEllipse";
    STEREON_TRY_VOID(Detail::CheckFinite(center, "Center", where));
    STEREON_TRY_VOID(Detail::CheckPositive(minorRadius, "MinorRadius", where));
    if (!Detail::IsFiniteValue(majorRadius) || !(majorRadius >= minorRadius)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeEllipse: MajorRadius must be finite and at least MinorRadius");
    }
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return BasicEllipse<double, D>{center, majorRadius, minorRadius, xAxis, yAxis};
}

} // namespace

Result<Line> MakeLine(const Point3& origin, const Vec3& direction)
{
    return MakeLineN<3>(origin, direction);
}

Result<Line2d> MakeLine(const Point2& origin, const Vec2& direction)
{
    return MakeLineN<2>(origin, direction);
}

Result<Circle> MakeCircle(const Point3& center, const double radius, const Vec3& xAxis,
                          const Vec3& yAxis)
{
    return MakeCircleN<3>(center, radius, xAxis, yAxis);
}

Result<Circle2d> MakeCircle(const Point2& center, const double radius, const Vec2& xAxis,
                            const Vec2& yAxis)
{
    return MakeCircleN<2>(center, radius, xAxis, yAxis);
}

Result<Ellipse> MakeEllipse(const Point3& center, const double majorRadius,
                            const double minorRadius, const Vec3& xAxis, const Vec3& yAxis)
{
    return MakeEllipseN<3>(center, majorRadius, minorRadius, xAxis, yAxis);
}

Result<Ellipse2d> MakeEllipse(const Point2& center, const double majorRadius,
                              const double minorRadius, const Vec2& xAxis, const Vec2& yAxis)
{
    return MakeEllipseN<2>(center, majorRadius, minorRadius, xAxis, yAxis);
}

} // namespace Stereon
