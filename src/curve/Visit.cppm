// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>
#include <string_view>
#include <type_traits>
#include <utility>

export module Stereon.Curve:Visit;

import Stereon.Robust;
import :Circle;
import :Concepts;
import :Ellipse;
import :Erased;
import :Line;

/// The fast path of ADR-0003 for curves. Hot loops call
/// VisitBuiltin(curve, f): when the curve wraps a built-in type, f receives
/// that type and its loop runs with direct calls; any other type reaches f as
/// the Curve itself, whose calls go through the virtual functions. The type
/// test runs once per VisitBuiltin, not once per evaluation.

export namespace Stereon {

/// The built-in curve types, in the order the intersection tables use.
using BuiltinCurveTypes = TypeList<Line, Circle, Ellipse>;
using BuiltinCurve2dTypes = TypeList<Line2d, Circle2d, Ellipse2d>;

template <std::size_t D>
using BuiltinCurveTypesN = std::conditional_t<D == 2, BuiltinCurve2dTypes, BuiltinCurveTypes>;

/// A short name for diagnostics, e.g. in "Line x Circle intersection is not
/// implemented": GeometryNameOf(static_cast<const T*>(nullptr)), which each
/// built-in type overloads and argument-dependent lookup finds; "geometry"
/// for any other type. Overloads rather than explicit specializations of
/// GeometryTypeName, because MSVC rejects specializations of an exported
/// template in a module interface and GCC 16.2 emits no definition for them
/// when another module odr-uses one.
template <typename T>
[[nodiscard]] constexpr std::string_view GeometryNameOf(const T*) noexcept
{
    return "geometry";
}

[[nodiscard]] constexpr std::string_view GeometryNameOf(const Line*) noexcept
{
    return "Line";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Circle*) noexcept
{
    return "Circle";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Ellipse*) noexcept
{
    return "Ellipse";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Line2d*) noexcept
{
    return "Line2d";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Circle2d*) noexcept
{
    return "Circle2d";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Ellipse2d*) noexcept
{
    return "Ellipse2d";
}

template <typename T>
inline constexpr std::string_view GeometryTypeName = GeometryNameOf(static_cast<const T*>(nullptr));

} // namespace Stereon

export namespace Stereon {

/// Calls f(geometry) with the wrapped built-in type, or f(curve) for any
/// other type. Every call of f must return the same type.
template <std::size_t D, typename F>
decltype(auto) VisitBuiltin(const BasicErasedCurve<D>& curve, F&& f)
    STEREON_HOT_PRE(!curve.IsEmpty())
{
    return Detail::VisitErased(BuiltinCurveTypesN<D>{}, curve, std::forward<F>(f));
}

/// The position of the curve's type in BuiltinCurveTypes (or
/// BuiltinCurve2dTypes), or their size when it is not built in.
template <std::size_t D>
[[nodiscard]] std::size_t BuiltinIndexOf(const BasicErasedCurve<D>& curve)
    STEREON_HOT_PRE(!curve.IsEmpty())
{
    return Detail::IndexIn(BuiltinCurveTypesN<D>{}, curve);
}

} // namespace Stereon
