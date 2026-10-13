// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <string_view>
#include <vector>

export module Stereon.Intersect:Dispatch;

import Stereon.Surface;

/// The type-pair dispatch of ADR-0003 for binary operations: one handler per
/// ordered pair of built-in types, in a table generated at compile time from
/// BuiltinCurveTypes and BuiltinSurfaceTypes. A pair with a closed-form
/// method gets an analytic handler (ADR-0010); every pair has a stub until
/// then, which reports NotImplemented naming both types. Types that are not
/// built in report NotSupported: the generic subdivision path needs tier 2
/// and 3 capabilities that no user type can declare yet.
///
/// The table is constant data, not a registry: nothing is added at run time,
/// so there is no global mutable state (ADR-0007).
///
/// Tier C (ADR-0015): this is scaffolding for line-by-line review. It
/// contains no intersection algorithm.

export namespace Stereon {

/// How two surfaces meet along a branch (ADR-0010, ADR-0013).
enum class ContactType : std::uint8_t {
    Crossing,
    Tangential,
    Overlap,
};

/// A point where two curves meet: its parameter on each and its position.
struct CurveCurvePoint {
    double ParameterA = 0.0;
    double ParameterB = 0.0;
    Point3 Position;
};

struct CurveCurveIntersection {
    std::vector<CurveCurvePoint> Points;
};

/// One branch of a surface-surface intersection: its 3D curve, its pcurve on
/// each surface on the shared parameter (ADR-0009), and how the surfaces
/// meet along it. The certificate of ADR-0009 joins it with the first
/// analytic handler.
struct IntersectionBranch {
    Curve Curve3d;
    Curve2d OnA;
    Curve2d OnB;
    ContactType Contact = ContactType::Crossing;
};

struct SurfaceSurfaceIntersection {
    std::vector<IntersectionBranch> Branches;
};

} // namespace Stereon

export namespace Stereon::Detail {

using CurveHandler = Result<CurveCurveIntersection> (*)(const Curve&, const Curve&, const Context&);
using SurfaceHandler = Result<SurfaceSurfaceIntersection> (*)(const Surface&, const Surface&,
                                                              const Context&);

template <typename A, typename B>
Result<CurveCurveIntersection> CurveStub(const Curve&, const Curve&, const Context&)
{
    // constexpr copies: GCC 16.2 emits no definition for the explicit
    // specializations of GeometryTypeName when another module odr-uses them.
    constexpr std::string_view nameA = GeometryTypeName<A>;
    constexpr std::string_view nameB = GeometryTypeName<B>;
    return Fail(
        ErrorCode::NotImplemented,
        std::format("{} x {} intersection is not implemented yet (ADR-0010)", nameA, nameB));
}

template <typename A, typename B>
Result<SurfaceSurfaceIntersection> SurfaceStub(const Surface&, const Surface&, const Context&)
{
    // constexpr copies: GCC 16.2 emits no definition for the explicit
    // specializations of GeometryTypeName when another module odr-uses them.
    constexpr std::string_view nameA = GeometryTypeName<A>;
    constexpr std::string_view nameB = GeometryTypeName<B>;
    return Fail(
        ErrorCode::NotImplemented,
        std::format("{} x {} intersection is not implemented yet (ADR-0010)", nameA, nameB));
}

/// table[i][j] handles type i of Ts against type j.
template <typename Handler, template <typename, typename> typename Stub, typename... Ts>
consteval auto MakePairTable(TypeList<Ts...>)
{
    constexpr std::size_t n = sizeof...(Ts);
    std::array<std::array<Handler, n>, n> table{};
    template for (constexpr std::size_t i : Indices<n>)
    {
        template for (constexpr std::size_t j : Indices<n>)
        {
            table[i][j] = Stub<Ts...[i], Ts...[j]>::Call;
        }
    }
    return table;
}

/// Class templates wrapping the stubs, so MakePairTable can take them as a
/// template template parameter.
template <typename A, typename B>
struct CurveStubOf {
    static constexpr CurveHandler Call = &CurveStub<A, B>;
};
template <typename A, typename B>
struct SurfaceStubOf {
    static constexpr SurfaceHandler Call = &SurfaceStub<A, B>;
};

inline constexpr auto CurvePairTable =
    MakePairTable<CurveHandler, CurveStubOf>(BuiltinCurveTypes{});
inline constexpr auto SurfacePairTable =
    MakePairTable<SurfaceHandler, SurfaceStubOf>(BuiltinSurfaceTypes{});

} // namespace Stereon::Detail

export namespace Stereon {

/// The points where two curves meet. Built-in pairs dispatch through the
/// pair table; other types are NotSupported.
[[nodiscard]] inline Result<CurveCurveIntersection> Intersect(const Curve& a, const Curve& b,
                                                              const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty curve");
    }
    const std::size_t i = BuiltinIndexOf(a);
    const std::size_t j = BuiltinIndexOf(b);
    if (i == BuiltinCurveTypes::Size || j == BuiltinCurveTypes::Size) {
        return Fail(ErrorCode::NotSupported,
                    "curve intersection needs built-in curve types until user types can declare "
                    "the capabilities of tiers 2 and 3 (ADR-0003)");
    }
    return Detail::CurvePairTable[i][j](a, b, ctx);
}

/// The branches along which two surfaces meet; see Intersect(Curve, Curve).
[[nodiscard]] inline Result<SurfaceSurfaceIntersection>
Intersect(const Surface& a, const Surface& b, const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty surface");
    }
    const std::size_t i = BuiltinIndexOf(a);
    const std::size_t j = BuiltinIndexOf(b);
    if (i == BuiltinSurfaceTypes::Size || j == BuiltinSurfaceTypes::Size) {
        return Fail(ErrorCode::NotSupported,
                    "surface intersection needs built-in surface types until user types can "
                    "declare the capabilities of tiers 2 and 3 (ADR-0003)");
    }
    return Detail::SurfacePairTable[i][j](a, b, ctx);
}

} // namespace Stereon
