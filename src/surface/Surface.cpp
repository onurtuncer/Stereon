// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <string_view>

module Stereon.Surface;

import :Erased;

namespace Stereon::Detail {

void CheckSurfaceParameter(const SurfaceConcept& model, const Point2& uv,
                           const std::string_view function) noexcept
{
    RequireFinite(uv, function, "a non-finite parameter");
    const Box2 domain = model.GetDomain();
    const Vec2 periods = model.GetPeriods();
    RequireInDomain(domain.X, uv.X, periods.X, function);
    RequireInDomain(domain.Y, uv.Y, periods.Y, function);
}

void CheckDerivativeOrder(const int du, const int dv, const std::string_view function) noexcept
{
    Require(du >= 0 && dv >= 0 && du + dv >= 1, function, "a derivative order that does not exist");
}

void CheckSurfacePoint(const Point3& p, const std::string_view function) noexcept
{
    RequireFinite(p, function, "a non-finite point");
}

} // namespace Stereon::Detail
