// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdint>
#include <limits>
#include <vector>

export module Stereon.Intersect:Outcome;

import Stereon.Surface;

/// The outcomes of an intersection (ADR-0029, accepted 2026-10-20): what a
/// complete, successful intersection returns. Success means the operation's
/// completeness and accuracy obligations are met over the requested domains
/// (ADR-0010), so a result with no components is a proved separation, never
/// "no seeds found". Anything less is a KernelError (IntersectionFailed,
/// CertificationFailed, NotSupported), never a partial success.
///
/// Components of different kinds may coexist. Exact analytic geometry stays
/// exact: a plane-plane branch is a Line, not a B-spline (ADR-0009).

export namespace Stereon {

/// How two surfaces meet along a branch (ADR-0010, ADR-0029 section 1).
enum class ContactType : std::uint8_t {
    /// Transverse: the surfaces cross.
    Crossing,
    /// The surfaces are tangent along the branch.
    Tangential,
    /// The surfaces coincide; see OverlapRegion. A branch of this kind is
    /// the boundary of an overlap.
    Overlap,
};

/// Whether the fragments on either side of a branch have the same
/// inside/outside state relative to the other solid (ADR-0029 section 3).
/// Established separately from the contact type: a tangency alone implies
/// neither relation.
enum class SideTransition : std::uint8_t {
    /// Not established. Classification must not propagate across it.
    Unresolved,
    /// The same state on both sides: a touch without a crossing.
    SameSide,
    /// Opposite states: a regular transverse crossing implies it.
    OppositeSide,
};

/// How a coincidence was established (ADR-0029 section 6).
enum class CoincidenceEvidence : std::uint8_t {
    /// The supports coincide exactly in the input doubles.
    Exact,
    /// The supports are exactly parallel and within the pair tolerance of
    /// each other: declared, logged, and no surface moved (Distance says how
    /// far apart they are).
    Tolerant,
};

/// An isolated point of contact (a tangent sphere's touch, say): where it
/// is, its parameters on both surfaces, and a bound on its distance from
/// the exact contact; infinite when not certified.
struct IntersectionPoint {
    Point3 Position;
    Point2 OnA;
    Point2 OnB;
    double Bound = std::numeric_limits<double>::infinity();
};

/// One branch of an intersection: its 3D curve, its pcurve on each surface
/// on the shared parameter (ADR-0009), the parameter range of the branch
/// inside the requested domains, how the surfaces meet along it, the side
/// transition across it, and the certified bound on its deviation from both
/// surfaces (ADR-0009 section 4; infinite when not certified).
struct IntersectionBranch {
    Curve Curve3d;
    Curve2d OnA;
    Curve2d OnB;
    Interval Range;
    ContactType Contact = ContactType::Crossing;
    SideTransition Transition = SideTransition::Unresolved;
    double Bound = std::numeric_limits<double>::infinity();

    /// True when the bound is a proved one.
    [[nodiscard]] bool IsCertified() const noexcept
    {
        return Bound < std::numeric_limits<double>::infinity();
    }
};

/// An affine map of the parameter plane, uv' = Offset + u U + v V: how a
/// point of an overlap on surface A corresponds to its point on B.
struct ParameterMap {
    Vec2 U{1.0, 0.0};
    Vec2 V{0.0, 1.0};
    Vec2 Offset{0.0, 0.0};

    [[nodiscard]] Point2 operator()(const Point2& uv) const noexcept
    {
        return Point2{Offset.X + uv.X * U.X + uv.Y * V.X, Offset.Y + uv.X * U.Y + uv.Y * V.Y};
    }
};

/// A region where the surfaces coincide (ADR-0029 section 1): its domain on
/// each surface, the correspondence of their parameters, whether their
/// normals agree, how the coincidence was established, and its boundary
/// branches inside the requested domains (none when the overlap fills
/// them). Coincident supports do not make overlapping faces: the faces'
/// trims decide that.
struct OverlapRegion {
    Box2 OnA;
    Box2 OnB;
    ParameterMap AToB;
    bool SameOrientation = true;
    CoincidenceEvidence Evidence = CoincidenceEvidence::Exact;
    /// The distance between tolerantly coincident supports; 0 when exact.
    double Distance = 0.0;
    /// True when the overlap fills both requested domains.
    bool FullDomain = false;
    std::vector<IntersectionBranch> Boundary;
};

/// A complete intersection of two surfaces over the requested domains.
struct SurfaceSurfaceIntersection {
    std::vector<IntersectionPoint> Points;
    std::vector<IntersectionBranch> Branches;
    std::vector<OverlapRegion> Overlaps;

    /// True when the surfaces are proved separate over the requested
    /// domains: a successful result with no component.
    [[nodiscard]] bool IsEmpty() const noexcept
    {
        return Points.empty() && Branches.empty() && Overlaps.empty();
    }
};

/// A point where two curves meet: its parameter on each and its position.
struct CurveCurvePoint {
    double ParameterA = 0.0;
    double ParameterB = 0.0;
    Point3 Position;
};

/// A stretch where two curves coincide: its parameter interval on each, and
/// whether they run the same way along it.
struct CurveCurveOverlap {
    Interval OnA;
    Interval OnB;
    bool SameDirection = true;
};

/// A complete intersection of two curves: points and coincident stretches.
struct CurveCurveIntersection {
    std::vector<CurveCurvePoint> Points;
    std::vector<CurveCurveOverlap> Overlaps;

    [[nodiscard]] bool IsEmpty() const noexcept { return Points.empty() && Overlaps.empty(); }
};

} // namespace Stereon
