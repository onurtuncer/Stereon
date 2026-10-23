// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md): every module
// implementation unit of a geometry library includes it, because one that
// calls core's inline functions with contracts and parses none itself
// leaves __tu_has_violation undefined, which the macOS linker rejects.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <numbers>
#include <string_view>

module Stereon.Surface;

/// The checked factories of the built-in surfaces (ADR-0006, ADR-0016). They
/// are defined here, out of line, because inline non-template functions in
/// a module interface that instantiate another module's templates (here
/// Result and the checks of Stereon.Core) corrupt GCC 16.2's module file
/// (docs/toolchain/gcc16-feature-probe.md).

namespace Stereon {

Result<Plane> MakePlane(const Point3& origin, const Vec3& xAxis, const Vec3& yAxis)
{
    constexpr std::string_view where = "MakePlane";
    STEREON_TRY_VOID(Detail::CheckFinite(origin, "Origin", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return Plane{origin, xAxis, yAxis};
}

Result<Cylinder> MakeCylinder(const Point3& origin, const double radius, const Vec3& xAxis,
                              const Vec3& yAxis, const bool rightHanded)
{
    constexpr std::string_view where = "MakeCylinder";
    STEREON_TRY_VOID(Detail::CheckFinite(origin, "Origin", where));
    STEREON_TRY_VOID(Detail::CheckPositive(radius, "Radius", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return Cylinder{origin, radius, xAxis, yAxis, rightHanded};
}

Result<Cone> MakeCone(const Point3& origin, const double radius, const double semiAngle,
                      const Vec3& xAxis, const Vec3& yAxis, const bool rightHanded)
{
    constexpr std::string_view where = "MakeCone";
    STEREON_TRY_VOID(Detail::CheckFinite(origin, "Origin", where));
    if (!Detail::IsFiniteValue(radius) || !(radius >= 0.0)) {
        return Fail(ErrorCode::InvalidInput, "MakeCone: Radius must be finite and not negative");
    }
    const double magnitude = semiAngle < 0.0 ? -semiAngle : semiAngle;
    if (!(magnitude > 0.0 && magnitude < std::numbers::pi / 2)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeCone: SemiAngle must be strictly between 0 and pi/2 in magnitude");
    }
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return Cone{origin, radius, semiAngle, xAxis, yAxis, rightHanded};
}

Result<Sphere> MakeSphere(const Point3& center, const double radius, const Vec3& xAxis,
                          const Vec3& yAxis, const bool rightHanded)
{
    constexpr std::string_view where = "MakeSphere";
    STEREON_TRY_VOID(Detail::CheckFinite(center, "Center", where));
    STEREON_TRY_VOID(Detail::CheckPositive(radius, "Radius", where));
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return Sphere{center, radius, xAxis, yAxis, rightHanded};
}

Result<Torus> MakeTorus(const Point3& center, const double majorRadius, const double minorRadius,
                        const Vec3& xAxis, const Vec3& yAxis, const bool rightHanded)
{
    constexpr std::string_view where = "MakeTorus";
    STEREON_TRY_VOID(Detail::CheckFinite(center, "Center", where));
    STEREON_TRY_VOID(Detail::CheckPositive(minorRadius, "MinorRadius", where));
    if (!Detail::IsFiniteValue(majorRadius) || !(majorRadius >= minorRadius)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeTorus: MajorRadius must be finite and at least MinorRadius "
                    "(ring and horn tori are supported)");
    }
    STEREON_TRY_VOID(Detail::CheckOrthonormal(xAxis, yAxis, where));
    return Torus{center, majorRadius, minorRadius, xAxis, yAxis, rightHanded};
}

} // namespace Stereon
