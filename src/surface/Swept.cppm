// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <memory>
#include <string_view>
#include <utility>

export module Stereon.Surface:Swept;

import Stereon.Curve;
import :BSplineSurface;
import :Concepts;

/// Swept surfaces: a profile curve moved along a direction (extrusion) or
/// turned about an axis (revolution), as STEP's surface_of_linear_extrusion
/// and surface_of_revolution and OCCT's Geom_SurfaceOfLinearExtrusion and
/// Geom_SurfaceOfRevolution.
///
/// Two forms. ExtrusionSurface and RevolutionSurface wrap any profile in an
/// erased Curve and evaluate the sweep itself: an extrusion is unbounded
/// along its direction, a revolution runs a full turn and is periodic in
/// the angle. MakeExtrusionSurface and MakeRevolutionSurface build the exact
/// tensor-product NURBS form of a B-spline profile over a finite range, the
/// form faces will carry; ToBSpline of a wrapper over a range is that.
///
/// Projection onto a wrapper works through the profile's exact NURBS form
/// where the kernel knows one over a finite domain (a circle, an ellipse, a
/// B-spline): an extrusion projects the point and that form along the
/// direction onto a plane and finds the closest point of the plane curve;
/// a revolution projects onto the NURBS form of its full turn. A line
/// profile, which has no finite domain, extrudes to a plane and is handled
/// as one. For other profiles, those of a user's type or an unbounded
/// conic, Project is approximate: the profile's own projection of the
/// point, with the sweep parameter chosen for it, and documented as such.

export namespace Stereon::Detail {

/// What an extrusion keeps for projection, shared by copies: the profile's
/// NURBS form and that form flattened along the direction, with a projector
/// on the flattened form. Defined in Swept.cpp.
struct ExtrusionProfile;
/// What a revolution keeps for projection: the NURBS form of its full turn
/// and a projector on it. Defined in Swept.cpp.
struct RevolvedPiece;

/// Null when the profile's type has no exact NURBS form over a finite
/// domain that the kernel knows.
[[nodiscard]] std::shared_ptr<const ExtrusionProfile> MakeExtrusionProfile(const Curve& profile,
                                                                           const Vec3& direction);
[[nodiscard]] std::shared_ptr<const RevolvedPiece>
MakeRevolvedPiece(const Curve& profile, const Point3& origin, const Vec3& axis);

/// The finite, non-zero vector v divided by its length, for every such v:
/// scaled first by the power of two that brings its largest component into
/// [1, 2), so that the squared length neither overflows (an axis of 1e200)
/// nor underflows (one of 1e-200). The scaling is exact, so wherever
/// Normalized(v) does not overflow or underflow the result is the same bit
/// for bit, and v and 2^k v give the same unit vector.
[[nodiscard]] Vec3 UnitVector(const Vec3& v);

} // namespace Stereon::Detail

export namespace Stereon {

/// S(u, v) = C(u) + v d: the profile C in u, the whole line of the direction
/// d in v, with d a displacement per unit v, not normalised.
class ExtrusionSurface {
public:
    /// Unchecked: the profile must not be empty and the direction must be
    /// finite and non-zero. MakeExtrusion checks data from outside the
    /// kernel.
    ExtrusionSurface(Curve profile, const Vec3& direction)
        : m_Profile(std::move(profile))
        , m_Direction(direction)
        , m_Projection(Detail::MakeExtrusionProfile(m_Profile, m_Direction))
    {}

    /// The profile, borrowed from this surface: lvalues only (ADR-0046).
    [[nodiscard]] const Curve& GetProfile() const& noexcept { return m_Profile; }
    const Curve& GetProfile() const&& = delete;
    [[nodiscard]] const Vec3& GetDirection() const noexcept { return m_Direction; }
    /// The projection data, or null for a profile without a known NURBS
    /// form over a finite domain.
    [[nodiscard]] const Detail::ExtrusionProfile* GetProjection() const noexcept
    {
        return m_Projection.get();
    }

private:
    Curve m_Profile;
    Vec3 m_Direction;
    std::shared_ptr<const Detail::ExtrusionProfile> m_Projection;
};

/// S(u, v) = O + a(v) + cos u r(v) + sin u (n x r(v)), the profile point
/// C(v) = O + a(v) + r(v) split into its part a along the unit axis n and
/// its radial part r, turned by the angle u about the axis: u in [0, 2 pi],
/// periodic; v the profile's parameter. The turn is right-handed about n.
class RevolutionSurface {
public:
    /// Unchecked: the profile must not be empty, the origin finite and the
    /// axis finite and non-zero; it is normalised (Detail::UnitVector), at
    /// any scale. MakeRevolution checks data from outside the kernel.
    RevolutionSurface(Curve profile, const Point3& origin, const Vec3& axis)
        : m_Profile(std::move(profile))
        , m_Origin(origin)
        , m_Axis(Detail::UnitVector(axis))
        , m_Piece(Detail::MakeRevolvedPiece(m_Profile, m_Origin, m_Axis))
    {}

    /// The tag of the constructor below.
    struct UnitAxis {};

    /// Unchecked, as above, with an axis already of unit length that is
    /// kept bit for bit: normalising a unit vector again can move it by an
    /// ulp. For restoring a surface's state, as the geometry archive does.
    RevolutionSurface(Curve profile, const Point3& origin, const Vec3& unitAxis, UnitAxis)
        : m_Profile(std::move(profile))
        , m_Origin(origin)
        , m_Axis(unitAxis)
        , m_Piece(Detail::MakeRevolvedPiece(m_Profile, m_Origin, m_Axis))
    {}

    /// The profile, borrowed from this surface: lvalues only (ADR-0046).
    [[nodiscard]] const Curve& GetProfile() const& noexcept { return m_Profile; }
    const Curve& GetProfile() const&& = delete;
    [[nodiscard]] const Point3& GetOrigin() const noexcept { return m_Origin; }
    /// Unit length.
    [[nodiscard]] const Vec3& GetAxis() const& noexcept { return m_Axis; }
    const Vec3& GetAxis() const&& = delete;
    /// The NURBS form of the full turn, or null for a profile without a
    /// known NURBS form over a finite domain.
    [[nodiscard]] const Detail::RevolvedPiece* GetPiece() const noexcept { return m_Piece.get(); }

private:
    Curve m_Profile;
    Point3 m_Origin;
    Vec3 m_Axis;
    std::shared_ptr<const Detail::RevolvedPiece> m_Piece;
};

/// The surfaces after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless the profile is not empty, the origin is finite, the
/// axis is finite and non-zero, of any length, and the direction is finite
/// with a squared length that is a finite, non-zero double: an extrusion
/// keeps its direction as given and divides by that squared length.
[[nodiscard]] Result<ExtrusionSurface> MakeExtrusion(const Curve& profile, const Vec3& direction);
[[nodiscard]] Result<RevolutionSurface> MakeRevolution(const Curve& profile, const Point3& origin,
                                                       const Vec3& axis);

[[nodiscard]] Point3 Eval(const ExtrusionSurface& surface, const Point2& uv);
/// S^(k,0) = C^(k)(u), S^(0,1) = d, everything else zero.
[[nodiscard]] Vec3 Deriv(const ExtrusionSurface& surface, const Point2& uv, int du, int dv);
/// The profile's domain in u, the whole real line in v.
[[nodiscard]] Box2 Domain(const ExtrusionSurface& surface);
/// (Period(profile), 0): periodic in u as the profile is, never in v.
[[nodiscard]] Vec2 Periods(const ExtrusionSurface& surface);
/// The profile's box swept along the direction: unbounded in every axis the
/// direction has a component along.
[[nodiscard]] Box3 BoundingBox(const ExtrusionSurface& surface);
/// The parameters of a closest point (see the comment at the top for which
/// profiles this is exact for).
[[nodiscard]] Point2 Project(const ExtrusionSurface& surface, const Point3& p);
/// The profile reversed, the same direction: the normal flips with S_u.
[[nodiscard]] ExtrusionSurface Reversed(const ExtrusionSurface& surface);
/// The profile and the direction transformed.
[[nodiscard]] ExtrusionSurface Transformed(const ExtrusionSurface& surface, const Transform& x);
/// The exact NURBS form over `range` (MakeExtrusionSurface of the profile's
/// form over range.X): range.X must lie in the profile's domain or, for a
/// periodic profile, span at most one period, and range.Y must be finite.
/// NotSupported for a profile without a known NURBS form.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const ExtrusionSurface& surface, Box2 range);

[[nodiscard]] Point3 Eval(const RevolutionSurface& surface, const Point2& uv);
/// S^(k,l): the radial part of C^(l)(v) turned by u and k quarter turns,
/// plus the axial part of C^(l)(v) when k is 0.
[[nodiscard]] Vec3 Deriv(const RevolutionSurface& surface, const Point2& uv, int du, int dv);
/// [0, 2 pi] in u, the profile's domain in v.
[[nodiscard]] Box2 Domain(const RevolutionSurface& surface);
/// (2 pi, Period(profile)).
[[nodiscard]] Vec2 Periods(const RevolutionSurface& surface);
/// The box of the NURBS form of the full turn, or, for a profile without
/// one, the box of the cylinder about the axis that holds the profile's
/// box turned.
[[nodiscard]] Box3 BoundingBox(const RevolutionSurface& surface);
/// The parameters of a closest point (see the comment at the top for which
/// profiles this is exact for).
[[nodiscard]] Point2 Project(const RevolutionSurface& surface, const Point3& p);
/// The same surface turned the other way about its axis, which flips the
/// normal: Eval(Reversed(s), (2 pi - u, v)) is Eval(s, (u, v)).
[[nodiscard]] RevolutionSurface Reversed(const RevolutionSurface& surface);
/// The profile, origin and axis transformed; a reflection also turns the
/// axis round, so that the turn's sense follows the map.
[[nodiscard]] RevolutionSurface Transformed(const RevolutionSurface& surface, const Transform& x);
/// The exact NURBS form over `range` (MakeRevolutionSurface of the
/// profile's form over range.Y through the angles range.X): range.X must be
/// at most one turn and range.Y must lie in the profile's domain or, for a
/// periodic profile, span at most one period. Its u is the rational arc
/// parameter, equal to the angle at the arc knots. NotSupported for a
/// profile without a known NURBS form.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const RevolutionSurface& surface, Box2 range);

[[nodiscard]] constexpr std::string_view GeometryNameOf(const ExtrusionSurface*) noexcept
{
    return "ExtrusionSurface";
}
[[nodiscard]] constexpr std::string_view GeometryNameOf(const RevolutionSurface*) noexcept
{
    return "RevolutionSurface";
}

static_assert(SurfaceGeometry<ExtrusionSurface>);
static_assert(PeriodicSurface<ExtrusionSurface>);
static_assert(SurfaceGeometry<RevolutionSurface>);
static_assert(PeriodicSurface<RevolutionSurface>);

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
