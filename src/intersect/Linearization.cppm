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
/// The pairs so far are a line, a circle, an ellipse or a B-spline curve
/// against a plane, a
/// cylinder, a cone, a sphere, a torus, an ellipsoid or a B-spline surface,
/// the parameters bound to coordinates of the line's origin and direction,
/// to a conic's centre and radii (ConicBinding), to the curve's poles and
/// weights (PoleBinding, as EvalJvp takes them),
/// and to the surface's origin, its radii, a cone's semi-angle, an
/// ellipsoid's semi-axes and a B-spline surface's poles and weights. The
/// surface's partials come from its
/// own Deriv, so a new surface needs only the velocities of its defining
/// values; the others answer NotSupported until they bind theirs.
///
/// A point of a branch where two surfaces meet has four unknowns and three
/// equations, S_A = S_B; the fourth is the normal-plane tracking
/// constraint (BranchPointLinearization), so the point moves within the
/// plane through it normal to the branch (the owner's choice of
/// 2026-11-08, docs/sensitivity.rst).

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

/// What a design parameter moves on a circle or an ellipse: a coordinate of
/// its centre; a circle's Radius; an ellipse's MajorRadius, along its XAxis,
/// and MinorRadius, along its YAxis; each a Length. The frame is fixed, as
/// a surface's is.
enum class ConicComponent : std::uint8_t {
    CenterX,
    CenterY,
    CenterZ,
    Radius,
    MajorRadius,
    MinorRadius,
};

/// A design parameter bound to one defining value of a circle or an
/// ellipse. Several bindings of one parameter add.
struct ConicBinding {
    ParameterId Parameter;
    ConicComponent Component = ConicComponent::CenterX;
};

/// What a design parameter moves on a surface: a coordinate of its origin
/// (a sphere's or a torus's centre, a cone's reference point); its Radius
/// (a sphere's, a cylinder's, a cone's reference radius, a torus's major
/// radius); a torus's MinorRadius; an ellipsoid's SemiAxisA, SemiAxisB and
/// SemiAxisC, along its XAxis, YAxis and axis; each a Length; or a cone's
/// SemiAngle, an Angle. A B-spline surface's poles move instead: a
/// coordinate of the binding's Pole (PoleX, PoleY, PoleZ), a Length, or its
/// weight (PoleWeight), Dimensionless, as PoleVelocity moves them. The
/// frame is fixed: turning a surface is a later binding.
enum class SurfaceComponent : std::uint8_t {
    OriginX,
    OriginY,
    OriginZ,
    Radius,
    MinorRadius,
    SemiAngle,
    SemiAxisA,
    SemiAxisB,
    SemiAxisC,
    PoleX,
    PoleY,
    PoleZ,
    PoleWeight,
};

/// A design parameter bound to one defining value of a surface. Pole is
/// the index into a B-spline surface's Poles, (i, j) row by row, for the
/// pole components, and unused otherwise.
struct SurfaceBinding {
    ParameterId Parameter;
    SurfaceComponent Component = SurfaceComponent::OriginX;
    std::size_t Pole = 0;
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
/// The same for a point of Intersect(curve, surface). InvalidInput also for
/// a point outside the curve's domain.
[[nodiscard]] Result<CurveSurfaceLinearization>
Linearize(const BSplineCurve& curve, const Surface& surface, const CurveSurfacePoint& point);
/// The same for a point of Intersect(circle, surface) or Intersect(ellipse,
/// surface), the point's parameter the conic's angle.
[[nodiscard]] Result<CurveSurfaceLinearization>
Linearize(const Circle& circle, const Surface& surface, const CurveSurfacePoint& point);
[[nodiscard]] Result<CurveSurfaceLinearization>
Linearize(const Ellipse& ellipse, const Surface& surface, const CurveSurfacePoint& point);

/// The forward product of Intersect(line, surface) for the seed of
/// `sensitivity`: for each point of `intersection`, in its order,
/// CurveSurfacePointValues values. `intersection` is the result of
/// intersecting these two geometries. NotSupported (a status) for a
/// surface other than a plane, a cylinder, a cone, a sphere, a torus, an
/// ellipsoid or a B-spline surface.
/// InvalidInput for a binding whose parameter is not of the context's set,
/// is not a Length (an Angle for a SemiAngle, Dimensionless for a
/// PoleWeight), or moves a value the surface does not have (a pole it does
/// not have included), or for a seed that was never set.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Line& line, std::span<const LineBinding> lineBindings, const Surface& surface,
             std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);
/// The same for a B-spline curve, its poles and weights moved by
/// `poleBindings` as EvalJvp moves them. InvalidInput also for a binding
/// naming a pole the curve does not have, a weight moved by anything but a
/// dimensionless parameter, or a point outside the curve's domain.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const BSplineCurve& curve, std::span<const PoleBinding> poleBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
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
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const BSplineCurve& curve, std::span<const PoleBinding> poleBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);

/// The products for a circle or an ellipse, its centre and radii moved by
/// `conicBindings`: C(t) = O + r (cos t X + sin t Y) moves by a unit along
/// a centre coordinate and by cos t X + sin t Y with r; an ellipse's
/// O + a cos t X + b sin t Y by cos t X with a and sin t Y with b.
/// InvalidInput as for a line, and also for a radius the conic does not
/// have (MajorRadius on a circle, Radius on an ellipse) or a point whose
/// angle is not finite.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Circle& circle, std::span<const ConicBinding> conicBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Ellipse& ellipse, std::span<const ConicBinding> conicBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const Circle& circle, std::span<const ConicBinding> conicBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const Ellipse& ellipse, std::span<const ConicBinding> conicBindings,
             const Surface& surface, std::span<const SurfaceBinding> surfaceBindings,
             const CurveSurfaceIntersection& intersection, const SensitivityContext& sensitivity);

/// A point of a surface–surface intersection to differentiate: a branch of
/// the result, by its index, at a parameter of its Range (ADR-0009's
/// shared edge parameter).
struct BranchSample {
    std::size_t Branch = 0;
    double Parameter = 0.0;
};

/// The values a product gives per sample of a surface–surface
/// intersection, in this order: the position's velocity x' (three), then
/// u_A', v_A', u_B' and v_B'.
inline constexpr std::size_t BranchPointValues = 7;

/// The linearisation of a branch point. The point solves four equations in
/// z = (u_A, v_A, u_B, v_B): S_A(u_A, v_A) - S_B(u_B, v_B) = 0, and the
/// tracking constraint T . (S_A(u_A, v_A) - X0) = 0 with X0 and T the
/// branch's point and unit tangent at the sample's parameter (the
/// normal-plane constraint): as the geometry moves, the point stays in the
/// plane through X0 normal to the branch. F_z has the rows [S_Au, S_Av,
/// -S_Bu, -S_Bv] and [T . S_Au, T . S_Av, 0, 0]; it is singular where the
/// surfaces are tangent, or where the branch runs inside the plane.
struct BranchPointLinearization {
    /// The point on both surfaces, refined from the branch's curves by
    /// Newton's method on the equations above, and its parameters.
    Point3 Position;
    Point2 OnA;
    Point2 OnB;
    /// S_A - S_B at the refined parameters: zero to rounding.
    Vec3 Residual;
    /// T, and the columns of the surfaces' partials. Named, not an array
    /// (see CurveSurfaceLinearization).
    Vec3 Tangent;
    Vec3 SuA;
    Vec3 SvA;
    Vec3 SuB;
    Vec3 SvB;
    /// |T . (n_A x n_B)| over unit normals: the sine of the angle between
    /// the surfaces when T is the branch's direction, 0 at a tangency.
    double Conditioning;
    DerivativeStatus Status;
};

/// The linearisation of `branch` of Intersect(a, b) at its parameter t.
/// The status: a Tangential or Overlap branch is Nondifferentiable; an
/// uncertified one, or a point at a pole of either surface's parameters,
/// Indeterminate; a crossing NearEvent below NearEventSine, else Valid.
/// InvalidInput for an empty surface, a parameter outside the branch's
/// Range, or a point Newton's method does not bring onto both surfaces
/// within the branch's bound.
[[nodiscard]] Result<BranchPointLinearization>
Linearize(const Surface& a, const Surface& b, const IntersectionBranch& branch, double t);

/// The forward product of Intersect(a, b) for the seed of `sensitivity`:
/// for each sample, in order, BranchPointValues values. The bindings of
/// either surface move it as for a curve–surface intersection; one
/// parameter may move both. A result with an overlap is Nondifferentiable
/// as a whole, and so is any sample on a branch that is not a crossing.
/// NotSupported (a status) unless both surfaces are planes, cylinders,
/// cones, spheres, tori or ellipsoids. InvalidInput as Linearize, for a sample naming
/// a branch the result does not have, a binding as IntersectJvp refuses it
/// above, or a seed that was never set.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Surface& a, std::span<const SurfaceBinding> bindingsA, const Surface& b,
             std::span<const SurfaceBinding> bindingsB,
             const SurfaceSurfaceIntersection& intersection, std::span<const BranchSample> samples,
             const SensitivityContext& sensitivity);

/// The reverse product: for an adjoint of BranchPointValues components per
/// sample, one value per parameter of the set. InvalidInput as IntersectJvp,
/// and for an adjoint of another length.
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const Surface& a, std::span<const SurfaceBinding> bindingsA, const Surface& b,
             std::span<const SurfaceBinding> bindingsB,
             const SurfaceSurfaceIntersection& intersection, std::span<const BranchSample> samples,
             const SensitivityContext& sensitivity);

/// What a design parameter moves on a plane curve: a line's origin and
/// direction coordinates; a circle's, an ellipse's or a hyperbola's centre
/// coordinates, a circle's Radius, an ellipse's or a hyperbola's MajorRadius
/// and MinorRadius; a parabola's vertex coordinates (CenterX, CenterY) and
/// FocalLength; a B-spline's pole coordinates (PoleX, PoleY) and weight
/// (PoleWeight) of the binding's Pole. Each a Length but a weight,
/// Dimensionless; frames fixed.
enum class PlaneCurveComponent : std::uint8_t {
    OriginX,
    OriginY,
    DirectionX,
    DirectionY,
    CenterX,
    CenterY,
    Radius,
    MajorRadius,
    MinorRadius,
    PoleX,
    PoleY,
    PoleWeight,
    FocalLength,
};

/// A design parameter bound to one defining value of a plane curve. Pole
/// is the index into a B-spline's Poles for the pole components, and unused
/// otherwise.
struct PlaneCurveBinding {
    ParameterId Parameter;
    PlaneCurveComponent Component = PlaneCurveComponent::OriginX;
    std::size_t Pole = 0;
};

/// The values a product gives per point where two plane curves meet, in this
/// order: the position's velocity x' (two), then s' and t', the velocities
/// of the parameters on A and on B.
inline constexpr std::size_t PlanarPointValues = 4;

/// The linearisation of one point where two plane curves meet. The point
/// solves F(s, t) = A(s) - B(t) = 0, two equations in two unknowns, with
/// F_z = [A'(s), -B'(t)], singular exactly where the curves are tangent: the
/// sine of the angle between them, |A' x B'| / (|A'| |B'|), is the
/// conditioning. A design change p' moves the parameters by F_z z' =
/// B_p p' - A_p p' and the point by x' = A_p p' + A' s'.
struct PlanarLinearization {
    /// A(s) - B(t) at the point's parameters: zero to rounding.
    Vec2 Residual;
    Vec2 TangentA;
    Vec2 TangentB;
    double Conditioning;
    DerivativeStatus Status;
};

/// The linearisation of a point of Intersect(a, b) for plane curves:
/// Valid or NearEvent for a decided crossing, Nondifferentiable for a
/// tangential one, Indeterminate for one the solver left undecided or whose
/// system is singular. InvalidInput for an empty curve or a point whose
/// parameters are not finite or lie off a curve's domain.
[[nodiscard]] Result<PlanarLinearization> Linearize(const Curve2d& a, const Curve2d& b,
                                                    const PlanarPoint& point);

/// The forward product of Intersect(a, b) for plane curves, each moved by
/// its bindings: for each point of `intersection`, in its order,
/// PlanarPointValues values. An overlap makes the result Nondifferentiable.
/// NotSupported (a status) for a curve of a kind no binding reaches, which
/// no Curve2d of the built-in kinds is. InvalidInput for a binding whose parameter is not
/// of the context's set or of the wrong quantity, or that moves a value the
/// curve does not have (a pole past the last included); for a point off a
/// curve or off the current geometry; or for a seed that was never set.
[[nodiscard]] Result<SensitivityResult>
IntersectJvp(const Curve2d& a, std::span<const PlaneCurveBinding> bindingsA, const Curve2d& b,
             std::span<const PlaneCurveBinding> bindingsB, const PlanarIntersection& intersection,
             const SensitivityContext& sensitivity);

/// The reverse product: for an adjoint of PlanarPointValues components per
/// point, one value per parameter of the set, accumulated over the points in
/// their order. InvalidInput as IntersectJvp, and for an adjoint of another
/// length.
[[nodiscard]] Result<SensitivityResult>
IntersectVjp(const Curve2d& a, std::span<const PlaneCurveBinding> bindingsA, const Curve2d& b,
             std::span<const PlaneCurveBinding> bindingsB, const PlanarIntersection& intersection,
             const SensitivityContext& sensitivity);

/// The products of Intersect(a, b) for two space curves (ADR-0052). Two space curves meet only
/// where their closest approach is within the resolution, and a generic change of either curve
/// undoes the meeting: a meeting is itself an event, of codimension one,
/// with no derivative. A result with a point or an overlap is therefore
/// Nondifferentiable, with that reason; a proved separation is Valid and
/// carries no values (the forward product) or zeros (the reverse one), since
/// it stays a separation under a small enough change. Tracking the closest
/// approach smoothly is a separate product, ClosestApproach, still to come;
/// no binding type is needed until then. InvalidInput for an empty curve,
/// for a seed that was never set (forward), or for an adjoint that is not
/// empty (reverse).
[[nodiscard]] Result<SensitivityResult> IntersectJvp(const Curve& a, const Curve& b,
                                                     const CurveCurveIntersection& intersection,
                                                     const SensitivityContext& sensitivity);
[[nodiscard]] Result<SensitivityResult> IntersectVjp(const Curve& a, const Curve& b,
                                                     const CurveCurveIntersection& intersection,
                                                     const SensitivityContext& sensitivity);

} // namespace Stereon
