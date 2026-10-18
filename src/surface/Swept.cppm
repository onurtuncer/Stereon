// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Surface:Swept;

import Stereon.Curve;
import :BSplineSurface;

export namespace Stereon {

/// Exact tensor-product extrusion: S(u, v) = C(u) + v * direction.
/// Keeps the profile parameter in u and the caller's finite range in v.
/// The direction is a displacement per unit v, not normalised. InvalidInput
/// for an invalid profile, zero/non-finite direction, invalid range or
/// non-finite output poles. Degenerate but finite patches are permitted.
[[nodiscard]] Result<BSplineSurface> MakeExtrusionSurface(const BSplineCurve& profile,
                                                          const Vec3& direction, Interval range);

/// Exact rational revolution about the line through origin along axis.
/// The axis may have any finite nonzero length and is normalised internally.
/// u is the rational circle parameter (equal to angle at the arc knots,
/// not between them); v retains the profile parameter. The angular range
/// must be finite, increasing and at most one turn. Profiles need not be
/// planar. Poles on the axis are allowed; regularity is not guaranteed.
/// InvalidInput for malformed data or non-finite generated poles/weights.
[[nodiscard]] Result<BSplineSurface> MakeRevolutionSurface(const BSplineCurve& profile,
                                                           const Point3& origin, const Vec3& axis,
                                                           Interval angles);

} // namespace Stereon
