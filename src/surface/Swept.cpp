// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <cmath>
#include <utility>

module Stereon.Surface;

import :BSplineSurface;
import :Swept;

namespace Stereon {
namespace {

bool Finite(const Vec3& v)
{
    return std::isfinite(v.X) && std::isfinite(v.Y) && std::isfinite(v.Z);
}

bool FiniteRange(const Interval range)
{
    return std::isfinite(range.GetLower()) && std::isfinite(range.GetUpper()) &&
           range.GetLower() < range.GetUpper() &&
           std::isfinite(range.GetUpper() - range.GetLower());
}

Result<void> CheckProfile(const BSplineCurve& profile)
{
    STEREON_TRY_VOID(
        MakeBSplineCurve(profile.Degree, profile.Knots, profile.Poles, profile.Weights));
    return {};
}

Result<BSplineSurface> CheckedSurface(BSplineSurface surface)
{
    return MakeBSplineSurface(surface.DegreeU, surface.DegreeV, std::move(surface.KnotsU),
                              std::move(surface.KnotsV), std::move(surface.Poles),
                              std::move(surface.Weights));
}

} // namespace

Result<BSplineSurface> MakeExtrusionSurface(const BSplineCurve& profile, const Vec3& direction,
                                            const Interval range)
{
    STEREON_TRY_VOID(CheckProfile(profile));
    if (!Finite(direction) || direction == Vec3{} || !FiniteRange(range)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeExtrusionSurface: expected a finite nonzero direction and finite range");
    }
    const double low = range.GetLower();
    const double high = range.GetUpper();
    BSplineSurface surface{profile.Degree, 1, profile.Knots, {low, low, high, high}, {}, {}};
    for (std::size_t i = 0; i < profile.Poles.size(); ++i) {
        surface.Poles.push_back(profile.Poles[i] + low * direction);
        surface.Poles.push_back(profile.Poles[i] + high * direction);
        if (IsRational(profile)) {
            surface.Weights.push_back(profile.Weights[i]);
            surface.Weights.push_back(profile.Weights[i]);
        }
    }
    return CheckedSurface(std::move(surface));
}

Result<BSplineSurface> MakeRevolutionSurface(const BSplineCurve& profile, const Point3& origin,
                                             const Vec3& axis, const Interval angles)
{
    STEREON_TRY_VOID(CheckProfile(profile));
    const double length = std::hypot(axis.X, axis.Y, axis.Z);
    if (!Finite(origin - Point3{}) || !Finite(axis) || !std::isfinite(length) || length == 0.0 ||
        !FiniteRange(angles) || angles.GetUpper() - angles.GetLower() > CirclePeriod) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeRevolutionSurface: expected a finite axis and at most one turn");
    }
    const Vec3 normal = axis / length;
    const Detail::UnitArc arc = Detail::MakeUnitArc(angles.GetLower(), angles.GetUpper());
    BSplineSurface surface{2, profile.Degree, arc.Knots, profile.Knots, {}, {}};
    for (std::size_t i = 0; i < arc.Poles.size(); ++i) {
        for (std::size_t j = 0; j < profile.Poles.size(); ++j) {
            const Vec3 offset = profile.Poles[j] - origin;
            const Vec3 axial = Dot(offset, normal) * normal;
            const Vec3 radial = offset - axial;
            surface.Poles.push_back(origin + axial + arc.Poles[i].X * radial +
                                    arc.Poles[i].Y * Cross(normal, radial));
            const double weight = IsRational(profile) ? profile.Weights[j] : 1.0;
            surface.Weights.push_back(arc.Weights[i] * weight);
        }
    }
    return CheckedSurface(std::move(surface));
}

} // namespace Stereon
