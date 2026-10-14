// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <string_view>

export module Stereon.Core:GeometryChecks;

import :Error;
import :Vector;

/// Validation of geometry that arrives from outside the kernel (ADR-0006,
/// ADR-0016). The checked factories of the geometry types (MakeEllipse,
/// MakeBSplineCurve, ...) and the readers of stored geometry call these and
/// return their KernelError, InvalidInput, to the caller; a contract would be
/// wrong here, because bad data is not a programming error.
///
/// Every message names the function and the offending value, e.g.
/// "MakeSphere: Radius must be positive". The functions are defined in
/// GeometryChecks.cpp, which keeps the message building out of inline code
/// (docs/toolchain/gcc16-feature-probe.md).

export namespace Stereon {

/// The relative allowance on a frame's orthonormality, dimensionless
/// (ADR-0028: a relative comparison allowance with scale 1, the length of a
/// unit vector). A unit vector computed in double carries a relative error
/// of a few units of 2^-53 per operation; 2^-40, about 1e-12, admits axes
/// that have been normalised, rotated and transformed many times, and
/// rejects any frame that is visibly skewed or scaled.
inline constexpr double FrameAllowance = 0x1p-40;

} // namespace Stereon

export namespace Stereon::Detail {

/// InvalidInput unless x is finite.
[[nodiscard]] Result<void> CheckFinite(double x, std::string_view what, std::string_view where);
[[nodiscard]] Result<void> CheckFinite(const Point2& p, std::string_view what,
                                       std::string_view where);
[[nodiscard]] Result<void> CheckFinite(const Point3& p, std::string_view what,
                                       std::string_view where);
[[nodiscard]] Result<void> CheckFinite(const Vec2& v, std::string_view what,
                                       std::string_view where);
[[nodiscard]] Result<void> CheckFinite(const Vec3& v, std::string_view what,
                                       std::string_view where);

/// InvalidInput unless x is finite and greater than zero.
[[nodiscard]] Result<void> CheckPositive(double x, std::string_view what, std::string_view where);

/// InvalidInput unless x and y are finite unit vectors at right angles, each
/// within FrameAllowance: |x.x - 1|, |y.y - 1| and |x.y| at most the
/// allowance.
[[nodiscard]] Result<void> CheckOrthonormal(const Vec2& x, const Vec2& y, std::string_view where);
[[nodiscard]] Result<void> CheckOrthonormal(const Vec3& x, const Vec3& y, std::string_view where);

} // namespace Stereon::Detail
