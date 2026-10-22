// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Intersect:Dispatch;

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
/// point of an overlap on surface A corresponds to its point on B. On a
/// periodic surface the image is taken modulo the period. Value-initialised
/// it is the zero map, not the identity: default member initialisers of
/// core's vector template here make GCC 16.2 write a module file its
/// importers cannot read once a std::vector holds the struct
/// (docs/toolchain/gcc16-feature-probe.md, "Bad file data").
struct ParameterMap {
    Vec2 U;
    Vec2 V;
    Vec2 Offset;
};

/// The image of uv under map. A free function defined in Dispatch.cpp
/// rather than a member: an importer of the module could not load the
/// member function's out-of-line declaration (the same compiler bug).
[[nodiscard]] Point2 Mapped(const ParameterMap& map, const Point2& uv) noexcept;

/// A region where the surfaces coincide (ADR-0029 section 1): its domain on
/// each surface, the correspondence of their parameters, whether their
/// normals agree, how the coincidence was established, and its boundary
/// branches inside the requested domains (none when the overlap fills
/// them), as indices into the result's Branches. Coincident supports do
/// not make overlapping faces: the faces' trims decide that.
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
    std::vector<std::size_t> Boundary;
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

/// The type-pair dispatch of ADR-0003 for binary operations: one handler per
/// ordered pair of built-in types, in a table generated at compile time from
/// BuiltinCurveTypes and BuiltinSurfaceTypes. A pair with a closed-form
/// method gets an analytic handler (ADR-0010): plane–plane, plane–cylinder
/// and cylinder–cylinder so far (the Boolean prototype's M3). Every other
/// pair has a stub, which reports NotImplemented naming both types. Types
/// that are not built in report NotSupported: the generic subdivision path
/// needs tier 2 and 3 capabilities that no user type can declare yet.
///
/// The table is constant data, not a registry: nothing is added at run time,
/// so there is no global mutable state (ADR-0007).
///
/// Tier C (ADR-0015): the table is scaffolding for line-by-line review; the
/// analytic handlers are in the implementation units named below.

export namespace Stereon::Detail {

/// A handler: the two geometries, each with the parameter domain the
/// intersection is requested over (ADR-0029 section 1), and for surfaces
/// the pair tolerance within which exactly parallel supports count as
/// coincident (ADR-0029 section 6).
using CurveHandler = Result<CurveCurveIntersection> (*)(const Curve&, const Interval&, const Curve&,
                                                        const Interval&, const Context&);
using SurfaceHandler = Result<SurfaceSurfaceIntersection> (*)(const Surface&, const Box2&,
                                                              const Surface&, const Box2&, double,
                                                              const Context&);

/// The analytic handlers (ADR-0010's table, RS-0001 Tier A): every
/// configuration decision is the sign of a polynomial in the input doubles
/// and the tolerance, settled by Stereon.Robust (RS-0001 G5). Defined in
/// PlanePlane.cpp, PlaneCylinder.cpp and CylinderCylinder.cpp.
Result<SurfaceSurfaceIntersection> PlanePlane(const Surface&, const Box2&, const Surface&,
                                              const Box2&, double, const Context&);
Result<SurfaceSurfaceIntersection> PlaneCylinder(const Surface&, const Box2&, const Surface&,
                                                 const Box2&, double, const Context&);
Result<SurfaceSurfaceIntersection> CylinderPlane(const Surface&, const Box2&, const Surface&,
                                                 const Box2&, double, const Context&);
Result<SurfaceSurfaceIntersection> CylinderCylinder(const Surface&, const Box2&, const Surface&,
                                                    const Box2&, double, const Context&);

template <typename A, typename B>
Result<CurveCurveIntersection> CurveStub(const Curve&, const Interval&, const Curve&,
                                         const Interval&, const Context&)
{
    // The names are read at compile time.
    constexpr std::string_view nameA = GeometryTypeName<A>;
    constexpr std::string_view nameB = GeometryTypeName<B>;
    std::string message(nameA);
    message += " x ";
    message += nameB;
    message += " intersection is not implemented yet (ADR-0010)";
    return Fail(ErrorCode::NotImplemented, std::move(message));
}

template <typename A, typename B>
Result<SurfaceSurfaceIntersection> SurfaceStub(const Surface&, const Box2&, const Surface&,
                                               const Box2&, double, const Context&)
{
    // The names are read at compile time.
    constexpr std::string_view nameA = GeometryTypeName<A>;
    constexpr std::string_view nameB = GeometryTypeName<B>;
    std::string message(nameA);
    message += " x ";
    message += nameB;
    message += " intersection is not implemented yet (ADR-0010)";
    return Fail(ErrorCode::NotImplemented, std::move(message));
}

/// One row of the table: type A against each of Ts.
template <typename Handler, template <typename, typename> typename Stub, typename A, typename... Ts>
consteval std::array<Handler, sizeof...(Ts)> MakePairRow(TypeList<Ts...>)
{
    return {Stub<A, Ts>::Call...};
}

/// table[i][j] handles type i of Ts against type j, generated by pack
/// expansion over the type list.
template <typename Handler, template <typename, typename> typename Stub, typename... Ts>
consteval auto MakePairTable(TypeList<Ts...> types)
{
    return std::array<std::array<Handler, sizeof...(Ts)>, sizeof...(Ts)>{
        MakePairRow<Handler, Stub, Ts>(types)...};
}

/// The handler of the ordered surface pair (A, B): an analytic one where
/// ADR-0010's table has it, the stub otherwise. A chain of if constexpr
/// rather than specialisations: GCC 16.2 and MSVC mishandle explicit
/// specialisations in module interfaces (docs/toolchain/gcc16-feature-probe.md).
template <typename A, typename B>
consteval SurfaceHandler SurfaceHandlerFor()
{
    if constexpr (std::same_as<A, Plane> && std::same_as<B, Plane>) {
        return &PlanePlane;
    } else if constexpr (std::same_as<A, Plane> && std::same_as<B, Cylinder>) {
        return &PlaneCylinder;
    } else if constexpr (std::same_as<A, Cylinder> && std::same_as<B, Plane>) {
        return &CylinderPlane;
    } else if constexpr (std::same_as<A, Cylinder> && std::same_as<B, Cylinder>) {
        return &CylinderCylinder;
    } else {
        return &SurfaceStub<A, B>;
    }
}

/// Class templates wrapping the handlers, so MakePairTable can take them as
/// a template template parameter.
template <typename A, typename B>
struct CurveStubOf {
    static constexpr CurveHandler Call = &CurveStub<A, B>;
};
template <typename A, typename B>
struct SurfaceHandlerOf {
    static constexpr SurfaceHandler Call = SurfaceHandlerFor<A, B>();
};

inline constexpr auto CurvePairTable =
    MakePairTable<CurveHandler, CurveStubOf>(BuiltinCurveTypes{});
inline constexpr auto SurfacePairTable =
    MakePairTable<SurfaceHandler, SurfaceHandlerOf>(BuiltinSurfaceTypes{});

} // namespace Stereon::Detail

export namespace Stereon {

/// The intersection of two curves over the requested parameter domains
/// (ADR-0029): a complete result, or the reason there is none. InvalidInput
/// for an empty curve or a domain that is empty or NaN-bounded; NotSupported
/// for a type that is not built in. Built-in pairs dispatch through the pair
/// table.
[[nodiscard]] Result<CurveCurveIntersection> Intersect(const Curve& a, const Interval& domainA,
                                                       const Curve& b, const Interval& domainB,
                                                       const Context& ctx);

/// Over the curves' whole domains.
[[nodiscard]] Result<CurveCurveIntersection> Intersect(const Curve& a, const Curve& b,
                                                       const Context& ctx);

/// The intersection of two surfaces over the requested parameter domains
/// (ADR-0029); see Intersect(Curve, Interval, Curve, Interval). A face's
/// boxes bound what a face-level operation asks for, and `tolerance` is the
/// pair's: the larger of the two faces' (ADR-0001). Exactly parallel
/// supports within it are declared coincident, with nothing moved
/// (ADR-0029 section 6); every other decision is made on the input doubles
/// as they are. InvalidInput unless the tolerance is finite and not
/// negative.
///
/// A branch's range is clipped to the requested domains with outward
/// rounding on the nominal curve: exactly for lines and for the angular
/// parameter of a cylinder's circles and ellipses, and otherwise only where
/// the branch's parameter box leaves the domain entirely. A face-level
/// operation restricts a branch to its loops itself. A box should include
/// the face's tolerance for the clip to keep every point within it.
[[nodiscard]] Result<SurfaceSurfaceIntersection> Intersect(const Surface& a, const Box2& domainA,
                                                           const Surface& b, const Box2& domainB,
                                                           double tolerance, const Context& ctx);

/// With the Context's resolution as the pair tolerance: the smallest any
/// entity has (ADR-0001).
[[nodiscard]] Result<SurfaceSurfaceIntersection> Intersect(const Surface& a, const Box2& domainA,
                                                           const Surface& b, const Box2& domainB,
                                                           const Context& ctx);

/// Over the surfaces' whole domains, at the Context's resolution.
[[nodiscard]] Result<SurfaceSurfaceIntersection> Intersect(const Surface& a, const Surface& b,
                                                           const Context& ctx);

} // namespace Stereon
