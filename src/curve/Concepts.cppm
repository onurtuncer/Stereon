// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <concepts>
#include <cstddef>
#include <type_traits>

export module Stereon.Curve:Concepts;

import Stereon.Robust;

/// The geometry concepts of ADR-0003 for curves: what a type must provide to
/// be wrapped in a Curve (in space) or a Curve2d (in a surface's parameter
/// plane). Geometry types are aggregates with no base class; the operations
/// are free functions found by argument-dependent lookup.
///
/// The concepts are checked on the double instance of a type. Eval, Deriv,
/// Reversed and Transformed of the built-in types are also generic over the
/// Scalar concept (ADR-0021); Domain, BoundingBox and Project bound and
/// decide, which stays double, like Interval and Box3.

export namespace Stereon {

/// The point, vector and box types of dimension 2 or 3.
template <Scalar T, std::size_t D>
    requires(D == 2 || D == 3)
using BasicPointN = std::conditional_t<D == 2, BasicPoint2<T>, BasicPoint3<T>>;

template <Scalar T, std::size_t D>
    requires(D == 2 || D == 3)
using BasicVecN = std::conditional_t<D == 2, BasicVec2<T>, BasicVec3<T>>;

template <std::size_t D>
    requires(D == 2 || D == 3)
using BoxN = std::conditional_t<D == 2, Box2, Box3>;

/// The tier-1 operations of a curve of dimension D (ADR-0003):
///
/// - Eval(c, t): the point at parameter t;
/// - Deriv(c, t, n): the n-th derivative with respect to t, n >= 1;
/// - Domain(c): the parameter interval, which may be unbounded;
/// - BoundingBox(c): a box containing the curve over its domain;
/// - Project(c, p): the parameter of a point of the curve closest to p;
/// - Reversed(c): the same point set traversed the other way.
template <typename C, std::size_t D>
concept CurveGeometryN = std::copyable<C> && requires(const C c, const double t, const int n,
                                                      const BasicPointN<double, D> p) {
    { Eval(c, t) } -> std::same_as<BasicPointN<double, D>>;
    { Deriv(c, t, n) } -> std::same_as<BasicVecN<double, D>>;
    { Domain(c) } -> std::same_as<Interval>;
    { BoundingBox(c) } -> std::same_as<BoxN<D>>;
    { Project(c, p) } -> std::same_as<double>;
    { Reversed(c) } -> std::same_as<C>;
};

/// A curve that is periodic in its parameter provides Period(c), the length
/// of one period, so that every finite parameter is valid for it; any other
/// curve must be evaluated inside Domain(c). Optional: the erased wrappers
/// report 0, "not periodic", for a type without it.
template <typename C>
concept PeriodicCurve = requires(const C c) {
    { Period(c) } -> std::same_as<double>;
};

/// A curve in space. It can also be moved: Transformed(c, x) applies a rigid
/// motion or a similarity.
template <typename C>
concept CurveGeometry = CurveGeometryN<C, 3> && requires(const C c, const Transform x) {
    { Transformed(c, x) } -> std::same_as<C>;
};

/// A curve in the (u, v) plane of a surface: a pcurve (ADR-0009). There is
/// no 2D transform type yet, so it has no Transformed.
template <typename C>
concept Curve2dGeometry = CurveGeometryN<C, 2>;

} // namespace Stereon
