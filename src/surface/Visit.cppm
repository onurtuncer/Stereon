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
#include <utility>

export module Stereon.Surface:Visit;

import Stereon.Curve;
import :Cylinder;
import :Erased;
import :Plane;

/// The fast path of ADR-0003 for surfaces; see Stereon.Curve:Visit.

export namespace Stereon {

/// The built-in surface types, in the order the intersection tables use.
using BuiltinSurfaceTypes = TypeList<Plane, Cylinder>;

/// The diagnostic names of the built-in surfaces; see GeometryNameOf in
/// Stereon.Curve:Visit.
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Plane*) noexcept
{
    return "Plane";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const Cylinder*) noexcept
{
    return "Cylinder";
}

/// Calls f(geometry) with the wrapped built-in type, or f(surface) for any
/// other type. Every call of f must return the same type.
template <typename F>
decltype(auto) VisitBuiltin(const Surface& surface, F&& f) STEREON_HOT_PRE(!surface.IsEmpty())
{
    return Detail::VisitErased(BuiltinSurfaceTypes{}, surface, std::forward<F>(f));
}

/// The position of the surface's type in BuiltinSurfaceTypes, or its size
/// when it is not built in.
[[nodiscard]] inline std::size_t BuiltinIndexOf(const Surface& surface)
    STEREON_HOT_PRE(!surface.IsEmpty())
{
    return Detail::IndexIn(BuiltinSurfaceTypes{}, surface);
}

} // namespace Stereon
