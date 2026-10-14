// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Surface:Conversion;

import Stereon.Curve;
import :BSplineSurface;
import :Cone;
import :Cylinder;
import :Erased;
import :Plane;
import :Sphere;
import :Torus;

/// Exact conversion of the analytic surfaces to NURBS (ADR-0003, tier 3),
/// over a finite parameter box (The NURBS Book, chapter 8). The plane
/// becomes a bilinear patch with its own parametrisation; the cylinder,
/// cone, sphere and torus become surfaces of revolution: a rational
/// quadratic circle arc in u swept over a NURBS meridian in v, which
/// agree with the original parameters at every knot (docs/nurbs.rst).

export namespace Stereon {

/// The surface over `range` as a B-spline surface. InvalidInput unless the
/// range is finite and non-empty in both parameters and inside the domain;
/// in a periodic parameter, any range at most one period long.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Plane& plane, Box2 range);
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Cylinder& cylinder, Box2 range);
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Cone& cone, Box2 range);
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Sphere& sphere, Box2 range);
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Torus& torus, Box2 range);
/// A B-spline surface restricted to `range`.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const BSplineSurface& surface, Box2 range);
/// Any built-in surface through the wrapper; NotSupported for other types.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const Surface& surface, Box2 range);

} // namespace Stereon
