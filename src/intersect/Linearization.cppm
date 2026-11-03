// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>

export module Stereon.Intersect:Linearization;

import Stereon.Surface;
import :Dispatch;

/// The linearisation of intersection results for shape sensitivities
/// (ADR-0021, section 6, and section 9's stage 4, moved to v0.2): a result
/// is differentiated through the equations that define it, never through
/// the solver that found it.
///
/// A point where a curve meets a surface solves F(z; p) = C(t; p) -
/// S(u, v; p) = 0 in z = (t, u, v): three equations in three unknowns, so
/// no tracking constraint is needed. With F_z = [C', -S_u, -S_v], a
/// design change p' moves the point by F_z z' = S_p p' - C_p p', and its
/// position by x' = C_p p' + C' t'. F_z is singular exactly where the curve
/// is tangent to the surface: its determinant is C' . (S_u x S_v), and
/// |det F_z| / (|C'| |S_u x S_v|) is the sine of the angle between them.
///
/// Whether a derivative exists is read from what the primal solver proved,
/// not guessed again: a crossing the solver decided is Valid, or NearEvent
/// when its angle is below NearEventSine; a Tangential point or an overlap
/// is Nondifferentiable, the crossing being a branch event there; a point
/// the solver left undecided is Indeterminate. A result with any point that
/// carries no derivative carries none as a whole (ADR-0021, section 5).
///
/// The pairs so far are a line against a plane, a sphere or a cylinder, the
/// parameters bound to coordinates of the line's origin and direction and
/// to the surface's origin or radius. The surface's partials come from its
/// own Deriv, so a new surface needs only the velocities of its defining
/// values; the others answer NotSupported until they bind theirs.

export namespace Stereon {

/// What a design parameter moves on a line: a coordinate of its origin or
/// of its direction, each a Length.
enum class LineComponent : std::uint8_t {
    OriginX,
    OriginY,
    OriginZ,
    DirectionX,
    DirectionY,
    DirectionZ,
};

/// A design parameter bound to one coordinate of a line. Several bindings
/// of one parameter add: a translation moves all three origin coordinates.
struct LineBinding {
    ParameterId Parameter;
    LineComponent Component = LineComponent::OriginX;
};

/// What a design parameter moves on a surface, each a Length: a coordinate
/// of its origin (a sphere's centre), or its radius (a sphere's or a
/// cylinder's). The frame is fixed: turning a surface is a later binding.
enum class SurfaceComponent : std::uint8_t {
    OriginX,
    OriginY,
    OriginZ,
    Radius,
};

struct SurfaceBinding {
    ParameterId Parameter;
    SurfaceComponent Component = SurfaceComponent::OriginX;
};

/// The values a product gives per point of a curve–surface intersection,
/// in this order: the position's velocity x' (three), then t', u' and v'.
inline constexpr std::size_t CurveSurfacePointValues = 6;

/// Below this sine of the angle between curve and surface a derivative is
/// NearEvent: its gain 1 / sin exceeds 2^20, so a relative change in a
/// defining value of 2^-20 moves the point by its own distance from the
/// tangency. An algorithm constant (ADR-0024), not a tolerance.
inline constexpr double NearEventSine = 0x1p-20;

/// The linearisation of one point where a curve meets a surface.
struct CurveSurfaceLinearization {
    /// F = C(t) - S(u, v) at the point's parameters: zero to rounding at a
    /// converged point.
    Vec3 Residual;
    /// The columns of F_z = [C'(t), -S_u, -S_v]: the curve's tangent and
    /// the surface's partials. Named, not an array of Vec3: GCC 16.2 cannot
    /// load an interface holding one (docs/toolchain/gcc16-feature-probe.md,
    /// "recursive lazy load").
    Vec3 Tangent;
    Vec3 Su;
    Vec3 Sv;
    /// |det F_z| / (|C'| |S_u x S_v|): the sine of the angle between the
    /// curve and the surface, 0 at a tangency.
    double Conditioning;
    /// Whether the point has a derivative, from what the solver proved and
    /// the conditioning.
    DerivativeStatus Status;
};

/// The linearisation of a point of Intersect(line, surface), for any
/// surface. InvalidInput for an empty surface or a point whose parameters
/// are not finite.
[[nodiscard]] Result<CurveSurfaceLinearization> Linearize(const Line& line, const Surface& surface,
                                                          const CurveSurfacePoint& point);

/// The forward product of Intersect(line, surface) for the seed of
/// `sensitivity`: for each point of `intersection`, in its order,
/// CurveSurfacePointValues values. `intersection` is the result of
/// intersecting these two geometries. NotSupported (a status) for a
/// surface other than a plane, a sphere or a cylinder. InvalidInput for a
/// binding whose parameter is not of the context's set or is not a Length,
/// a radius binding on a plane, or a seed that was never set.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Line& line, std::span<const LineBinding> lineBindings, const Surface& surface,
             std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);

/// The reverse product: for the adjoint of `sensitivity`, a covector of
/// CurveSurfacePointValues components per point under the Euclidean
/// pairing, one value per parameter of the set, accumulated over the
/// points in their order (ADR-0021, section 7). InvalidInput as
/// IntersectJvp, and for an adjoint of another length.
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const Line& line, std::span<const LineBinding> lineBindings, const Surface& surface,
             std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);

} // namespace Stereon
