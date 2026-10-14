// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <concepts>

export module Stereon.Surface:Concepts;

import Stereon.Curve;

/// The geometry concept of ADR-0003 for surfaces. A surface maps a parameter
/// point uv = (u, v) to a point in space; as for curves, the operations are
/// free functions found by argument-dependent lookup and are checked on the
/// double instance of a type.

export namespace Stereon {

/// The tier-1 operations of a surface (ADR-0003):
///
/// - Eval(s, uv): the point at parameter uv;
/// - Deriv(s, uv, du, dv): the mixed partial derivative of order du in u and
///   dv in v, du + dv >= 1;
/// - Domain(s): the parameter box, which may be unbounded;
/// - BoundingBox(s): a box containing the surface over its domain;
/// - Project(s, p): the parameter of a point of the surface closest to p;
/// - Reversed(s): the same surface with the opposite normal;
/// - Transformed(s, x): its image under a rigid motion or a similarity.
/// A surface periodic in u or v provides Periods(s): the period in each
/// parameter, 0 for none. Optional, as PeriodicCurve is for curves.
template <typename S>
concept PeriodicSurface = requires(const S s) {
    { Periods(s) } -> std::same_as<Vec2>;
};

template <typename S>
concept SurfaceGeometry =
    std::copyable<S> && requires(const S s, const Point2 uv, const int du, const int dv,
                                 const Point3 p, const Transform x) {
        { Eval(s, uv) } -> std::same_as<Point3>;
        { Deriv(s, uv, du, dv) } -> std::same_as<Vec3>;
        { Domain(s) } -> std::same_as<Box2>;
        { BoundingBox(s) } -> std::same_as<Box3>;
        { Project(s, p) } -> std::same_as<Point2>;
        { Reversed(s) } -> std::same_as<S>;
        { Transformed(s, x) } -> std::same_as<S>;
    };

} // namespace Stereon
