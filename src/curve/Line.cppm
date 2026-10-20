// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>

export module Stereon.Curve:Line;

import Stereon.Robust;
import :Concepts;

/// The straight line, in space and in the parameter plane.

export namespace Stereon {

/// Origin + t * Direction for every real t. Direction need not be a unit
/// vector; its length sets the speed of the parametrisation.
template <Scalar T, std::size_t D>
struct BasicLine {
    BasicPointN<T, D> Origin;
    BasicVecN<T, D> Direction;

    constexpr bool operator==(const BasicLine&) const = default;
};

using Line = BasicLine<double, 3>;
using Line2d = BasicLine<double, 2>;

/// The line after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless Origin and Direction are finite and Direction is not
/// zero.
[[nodiscard]] Result<Line> MakeLine(const Point3& origin, const Vec3& direction);
[[nodiscard]] Result<Line2d> MakeLine(const Point2& origin, const Vec2& direction);

template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicPointN<T, D> Eval(const BasicLine<T, D>& line, const T& t)
{
    Detail::RequireFinite(t, "Eval(Line)", "a non-finite parameter");
    return line.Origin + t * line.Direction;
}

/// The direction for n == 1, zero above.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicVecN<T, D> Deriv(const BasicLine<T, D>& line,
                                              [[maybe_unused]] const T& t, const int n)
{
    Detail::RequireFinite(t, "Deriv(Line)", "a non-finite parameter");
    Detail::Require(n >= 1, "Deriv(Line)", "a derivative order below 1");
    return n == 1 ? line.Direction : BasicVecN<T, D>{};
}

/// The whole real line.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain([[maybe_unused]] const BasicLine<T, D>& line) noexcept
{
    return Interval::Everything();
}

/// Everything: a line is unbounded.
template <std::size_t D>
[[nodiscard]] constexpr BoxN<D> BoundingBox([[maybe_unused]] const BasicLine<double, D>& line)
{
    return BoxN<D>::Everything();
}

/// The parameter of the foot of the perpendicular from p.
template <std::size_t D>
[[nodiscard]] constexpr double Project(const BasicLine<double, D>& line,
                                       const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(Line)", "a non-finite point");
    Detail::Require(LengthSquared(line.Direction) > 0.0, "Project(Line)", "a zero direction");
    return Dot(p - line.Origin, line.Direction) / LengthSquared(line.Direction);
}

/// The same line traversed the other way: Eval(Reversed(l), t) == Eval(l, -t).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicLine<T, D> Reversed(const BasicLine<T, D>& line)
{
    return {line.Origin, -line.Direction};
}

/// The image of the line under x, with the parameter kept:
/// Eval(Transformed(l, x), t) == Transformed(Eval(l, t), x).
template <Scalar T>
[[nodiscard]] constexpr BasicLine<T, 3> Transformed(const BasicLine<T, 3>& line,
                                                    const BasicTransform<T>& x)
{
    return {Transformed(line.Origin, x), Transformed(line.Direction, x)};
}

static_assert(CurveGeometry<Line>);
static_assert(Curve2dGeometry<Line2d>);
static_assert(Detail::FitsInline<Line, CurveBufferSize>);
static_assert(Detail::FitsInline<Line2d, Curve2dBufferSize>);

} // namespace Stereon
