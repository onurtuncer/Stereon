// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Build:Primitives;

import Stereon.Topo;

/// Primitive solids (the Boolean prototype's M1, and ADR-0048 §4): a box, a
/// cylinder, a sphere and an ellipsoid, each one solid of one closed shell, built through
/// ShapeBuilder's controlled operations and published by Commit, so the result is certified at the
/// Structural level and checked at the Geometric level (ADR-0001's
/// containment invariant) like any operation's.
///
/// A frame is an origin and two orthonormal axes, x and y, as for the
/// analytic surfaces; z is their cross product. Every edge is Exact (a line,
/// a circle or an ellipse, ADR-0009) with exact pcurves on its faces,
/// including the seams, used twice by one face, except the poles of a
/// sphere or an ellipsoid, which are Degenerate edges. Vertex and edge
/// tolerances are the Context's resolution: native construction is
/// resolution-tight (ADR-0001). Faces are oriented outward.

export namespace Stereon {

/// The box with a corner at `origin`, extending `size.X` along `xAxis`,
/// `size.Y` along `yAxis` and `size.Z` along their cross product: 8 vertices,
/// 12 line edges, 6 planar faces. InvalidInput for a frame that is not
/// orthonormal within FrameAllowance, or a size that is not finite or not
/// larger than the Context's resolution in each direction.
[[nodiscard]] Result<OperationResult> MakeBoxSolid(const Point3& origin, const Vec3& xAxis,
                                                   const Vec3& yAxis, const Vec3& size,
                                                   const Context& context);

/// The axis-aligned box with the corners of `box`, as above.
[[nodiscard]] Result<OperationResult> MakeBoxSolid(const Box3& box, const Context& context);

/// The cylinder whose bottom disk is centred at `baseCentre` in the plane of
/// `xAxis` and `yAxis`, extending `height` along their cross product: 2
/// vertices on the seam at angle 0, a bottom and a top circle, the seam, a
/// lateral cylindrical face and two planar disks. InvalidInput for a frame
/// that is not orthonormal, or a radius or height that is not finite or not
/// larger than the Context's resolution.
[[nodiscard]] Result<OperationResult> MakeCylinderSolid(const Point3& baseCentre, const Vec3& xAxis,
                                                        const Vec3& yAxis, double radius,
                                                        double height, const Context& context);

/// The sphere of `radius` about `center`, its poles along the cross product
/// of `xAxis` and `yAxis`: 2 vertices at the poles, an Exact meridian seam
/// at angle 0 from the south pole to the north, a semicircle, two
/// Degenerate edges at the poles, and one spherical face. InvalidInput for
/// a frame that is not orthonormal, or a radius that is not finite or not
/// larger than the Context's resolution.
[[nodiscard]] Result<OperationResult> MakeSphereSolid(const Point3& center, const Vec3& xAxis,
                                                      const Vec3& yAxis, double radius,
                                                      const Context& context);

/// The ellipsoid with the semi-axes `a`, `b` and `c` along `xAxis`, `yAxis`
/// and their cross product (ADR-0048 §4): the sphere's construction, with
/// the seam a half ellipse. InvalidInput as for the sphere, every semi-axis
/// larger than the resolution. Booleans do not take it yet: their
/// classification and clipping know planes and cylinders.
[[nodiscard]] Result<OperationResult> MakeEllipsoidSolid(const Point3& center, const Vec3& xAxis,
                                                         const Vec3& yAxis, double a, double b,
                                                         double c, const Context& context);

} // namespace Stereon
