// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Core:Frame;

import :Scalar;
import :Trigonometry;
import :Vector;

/// The arithmetic shared by the angle-parametrised curves and surfaces: the
/// axis of a frame of either handedness, the radial direction at an angle
/// and its derivatives, and angles wrapped to one period. It lives in core
/// because a geometry module must not build its functions out of another
/// geometry module's templates (docs/toolchain/gcc16-feature-probe.md).

export namespace Stereon {

/// The third axis of the frame (x, y): Cross(x, y) when rightHanded, its
/// negative otherwise. A surface of revolution stores the flag rather than
/// the axis, which keeps it small and lets Reversed() turn the normal round
/// by flipping it (see BasicCylinder).
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> FrameAxis(const BasicVec3<T>& x, const BasicVec3<T>& y,
                                               const bool rightHanded)
{
    const BasicVec3<T> z = Cross(x, y);
    return rightHanded ? z : -z;
}

/// cos(u + n pi/2) * x + sin(u + n pi/2) * y: the unit radial direction at
/// angle u for n = 0, and its n-th derivative in u, by quarter turns
/// (QuarterTurned), so that no rounded multiple of pi/2 enters an argument.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> RadialDirection(const BasicVec3<T>& x, const BasicVec3<T>& y,
                                                     const T& u, const int n)
{
    const BasicCosSin<T> turned = QuarterTurned(Cos(u), Sin(u), n);
    return turned.Cos * x + turned.Sin * y;
}

/// The polar angle of (x, y) in [0, 2 pi), and 0 for the origin, where
/// every angle is equally valid. atan2 returns (-pi, pi]; a tiny negative
/// angle plus 2 pi rounds to 2 pi itself, which is the same direction as 0
/// and outside the half-open period, so it becomes 0.
[[nodiscard]] double WrappedAngle(double y, double x) noexcept;

} // namespace Stereon
