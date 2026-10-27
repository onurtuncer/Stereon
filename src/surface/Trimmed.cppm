// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <memory>
#include <string_view>
#include <utility>

export module Stereon.Surface:Trimmed;

import Stereon.Curve;
import :BSplineSurface;
import :Concepts;
import :Erased;

/// Rectangular trimmed surfaces: a surface restricted to a box of its
/// parameter plane (STEP's rectangular_trimmed_surface, OCCT's
/// Geom_RectangularTrimmedSurface), such as a quarter of a cylinder or a
/// patch of a plane.
///
/// A TrimmedSurface wraps any surface in an erased Surface and keeps its
/// parametrisation: Eval(trimmed, uv) is Eval(basis, uv) for uv in the box,
/// so a parameter means the same on both. Reversing it runs u backwards,
/// u -> lower + upper - u, which flips the normal, rather than reversing the
/// basis, whose own parametrisation under Reversed is its type's business.
/// A trimmed surface is periodic in neither parameter: its domain is the
/// box, whatever the basis repeats.
///
/// Projection must find the closest point within the box, not the basis's
/// global closest point clamped to it: a point beside a quarter cylinder
/// is closest to an interior point of the quarter that clamping the
/// cylinder's own projection never gives. For the built-in types the
/// trimmed surface keeps the exact NURBS form of its patch (ToBSpline over
/// the box) and projects onto that with the B-spline surface projection,
/// then maps the result back to the basis's parameter. For a user's
/// surface type, which has no known NURBS form, Project compares the
/// basis's global projection, brought into the box, with the box's
/// corners: correct when the closest point on the basis lies in the box or
/// at a corner, approximate otherwise.

export namespace Stereon::Detail {

/// The exact NURBS form of a trimmed patch and its projector, shared by
/// the copies of one trimmed surface. Defined in Trimmed.cpp.
struct TrimmedPatch;

/// The patch of a basis over a box, or null when the basis's type has no
/// exact NURBS form the kernel knows.
[[nodiscard]] std::shared_ptr<const TrimmedPatch> MakePatch(const Surface& basis, Box2 range);

} // namespace Stereon::Detail

export namespace Stereon {

class TrimmedSurface {
public:
    /// The basis over `range`, reversed or not. Unchecked: each interval of
    /// the range must be finite and non-empty, and inside the basis's domain
    /// in that direction or, where the basis is periodic, at most one period
    /// wide. MakeTrimmedSurface checks data from outside the kernel.
    TrimmedSurface(Surface basis, const Box2 range, const bool reversed = false)
        : m_Basis(std::move(basis))
        , m_Range(range)
        , m_Reversed(reversed)
        , m_Patch(Detail::MakePatch(m_Basis, m_Range))
    {}

    [[nodiscard]] const Surface& GetBasis() const noexcept { return m_Basis; }
    [[nodiscard]] Box2 GetRange() const noexcept { return m_Range; }
    [[nodiscard]] bool IsReversed() const noexcept { return m_Reversed; }
    /// The exact NURBS patch, or null for a user's basis type.
    [[nodiscard]] const Detail::TrimmedPatch* GetPatch() const noexcept { return m_Patch.get(); }

    /// The basis's parameters at this surface's parameters uv.
    [[nodiscard]] Point2 BasisParameter(const Point2& uv) const noexcept
    {
        return m_Reversed ? Point2{Reflect(uv.X), uv.Y} : uv;
    }

    /// lower + upper - u in the u range: the parameter of the same point
    /// when the sense is flipped, exact at the ends and kept in the range
    /// where rounding the sum would leave it.
    [[nodiscard]] double Reflect(const double u) const noexcept
    {
        const double lower = m_Range.X.GetLower();
        const double upper = m_Range.X.GetUpper();
        if (u == lower) {
            return upper;
        }
        if (u == upper) {
            return lower;
        }
        return std::clamp(lower + upper - u, lower, upper);
    }

    /// The same patch with the opposite normal. The patch is shared, not
    /// rebuilt.
    [[nodiscard]] TrimmedSurface Reversed() const
    {
        TrimmedSurface copy = *this;
        copy.m_Reversed = !m_Reversed;
        return copy;
    }

private:
    Surface m_Basis;
    Box2 m_Range;
    bool m_Reversed = false;
    std::shared_ptr<const Detail::TrimmedPatch> m_Patch;
};

/// The trimmed surface after validating data from outside the kernel
/// (ADR-0006): InvalidInput unless the basis is not empty and each interval
/// of the range is finite, non-empty, and inside the basis's domain in that
/// direction or, where the basis is periodic, at most one period wide (a
/// patch may cross a seam: a cylinder trimmed to u in [5, 7] is fine).
[[nodiscard]] Result<TrimmedSurface> MakeTrimmedSurface(const Surface& basis, Box2 range);

[[nodiscard]] Point3 Eval(const TrimmedSurface& surface, const Point2& uv);
/// The partial derivative; odd orders in u change sign on a reversed surface.
[[nodiscard]] Vec3 Deriv(const TrimmedSurface& surface, const Point2& uv, int du, int dv);
[[nodiscard]] Box2 Domain(const TrimmedSurface& surface);
/// A box containing the patch: its NURBS form's pole box, or the basis's
/// box for a user's type.
[[nodiscard]] Box3 BoundingBox(const TrimmedSurface& surface);
/// The parameters of a closest point of the patch (see the comment above
/// for a user's basis type). Ties go to the smaller parameters, u first
/// (ADR-0007).
[[nodiscard]] Point2 Project(const TrimmedSurface& surface, const Point3& p);
[[nodiscard]] TrimmedSurface Reversed(const TrimmedSurface& surface);
/// The trimmed surface moved: the basis transformed, and the range rescaled
/// as the basis's Transformed rescales its parameters under a similarity (a
/// plane's both, a cylinder's or cone's v, a sweep's as its profile's), so
/// that the patch is the image of the patch.
[[nodiscard]] TrimmedSurface Transformed(const TrimmedSurface& surface, const Transform& x);
/// The exact NURBS form of the whole patch, with the normal of the trimmed
/// surface. NotSupported for a user's basis type.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const TrimmedSurface& surface);
/// The exact NURBS form of the part of the patch over `range`, which must
/// lie in the trimmed surface's domain, parametrised as the basis's piece
/// and reversed in u when the trimmed surface is. NotSupported for a
/// user's basis type.
[[nodiscard]] Result<BSplineSurface> ToBSpline(const TrimmedSurface& surface, Box2 range);

[[nodiscard]] constexpr std::string_view GeometryNameOf(const TrimmedSurface*) noexcept
{
    return "TrimmedSurface";
}

static_assert(SurfaceGeometry<TrimmedSurface>);

} // namespace Stereon
