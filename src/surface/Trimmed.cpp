// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

module Stereon.Surface;

import :BSplineSurface;
import :Cone;
import :Conversion;
import :Cylinder;
import :Erased;
import :Plane;
import :Swept;
import :Trimmed;
import :Visit;

namespace Stereon::Detail {

struct TrimmedPatch {
    /// The patch, exactly.
    BSplineSurface Surface;
    SurfaceProjector Projector;
    /// True when the patch is parametrised as the basis (a plane or a
    /// B-spline); otherwise a point of the patch is mapped back to the basis
    /// by projecting onto it.
    bool SameParameter = false;
};

} // namespace Stereon::Detail

namespace Stereon {

namespace {

/// True when the basis's NURBS form keeps its parametrisation: a plane or a
/// B-spline surface (src/surface/Conversion.cppm); the rational arcs of the
/// other analytic surfaces agree with it only at their knots.
bool PatchKeepsParameter(const Surface& basis)
{
    return VisitBuiltin(basis, []<typename G>(const G&) {
        return std::is_same_v<G, Plane> || std::is_same_v<G, BSplineSurface>;
    });
}

/// uv brought into the box: by whole periods where the basis has them,
/// then clamped (NearestInRange).
Point2 IntoBox(const Point2& uv, const Box2& range, const Vec2& periods)
{
    return {NearestInRange(uv.X, range.X, periods.X), NearestInRange(uv.Y, range.Y, periods.Y)};
}

/// The basis parameters of the patch's closest point to p, or of the
/// basis's global closest point brought into the box for a user type.
Point2 InteriorCandidate(const TrimmedSurface& surface, const Point3& p)
{
    const Surface& basis = surface.GetBasis();
    const Box2 range = surface.GetRange();
    const Detail::TrimmedPatch* patch = surface.GetPatch();
    if (patch == nullptr) {
        return IntoBox(Project(basis, p), range, Periods(basis));
    }
    const Point2 st = patch->Projector.Project(p);
    if (patch->SameParameter) {
        return IntoBox(st, range, {});
    }
    // The patch's point lies on the basis, an analytic surface that does not
    // cross itself, so projecting it back finds its own parameters, up to
    // the period, which the box fixes.
    return IntoBox(Project(basis, Eval(patch->Surface, st)), range, Periods(basis));
}

/// True when a comes before b: the smaller u, then the smaller v (ADR-0007).
bool Earlier(const Point2& a, const Point2& b)
{
    return a.X < b.X || (a.X == b.X && a.Y < b.Y);
}

Point2 ProjectTrimmed(const TrimmedSurface& surface, const Point3& p)
{
    const Box2& r = surface.GetRange();
    const std::array<Point2, 5> candidates{
        InteriorCandidate(surface, p),          Point2{r.X.GetLower(), r.Y.GetLower()},
        Point2{r.X.GetUpper(), r.Y.GetLower()}, Point2{r.X.GetLower(), r.Y.GetUpper()},
        Point2{r.X.GetUpper(), r.Y.GetUpper()},
    };
    Point2 best = candidates[0];
    double bestDistance = DistanceSquared(Eval(surface.GetBasis(), best), p);
    for (std::size_t k = 1; k < candidates.size(); ++k) {
        const double d = DistanceSquared(Eval(surface.GetBasis(), candidates[k]), p);
        if (d < bestDistance || (d == bestDistance && Earlier(candidates[k], best))) {
            best = candidates[k];
            bestDistance = d;
        }
    }
    return surface.BasisParameter(best); // reflecting is its own inverse
}

/// InvalidInput unless the interval is finite, non-empty and inside the
/// domain or, for a periodic direction, at most one period wide.
Result<void> CheckTrim(const Interval range, const Interval domain, const double period,
                       const char* direction)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        return Fail(ErrorCode::InvalidInput, std::string("MakeTrimmedSurface: the ") + direction +
                                                 " range must be finite and non-empty");
    }
    if (period > 0.0 ? high - low > period : !Contains(domain, range)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string("MakeTrimmedSurface: the ") + direction +
                        (period > 0.0 ? " range spans more than one period"
                                      : " range must lie in the basis's domain"));
    }
    return {};
}

} // namespace

namespace Detail {

std::shared_ptr<const TrimmedPatch> MakePatch(const Surface& basis, const Box2 range)
{
    if (basis.IsEmpty()) {
        return nullptr;
    }
    Result<BSplineSurface> patch = ToBSpline(basis, range);
    if (!patch) {
        return nullptr;
    }
    SurfaceProjector projector(*patch);
    return std::make_shared<const TrimmedPatch>(
        TrimmedPatch{*std::move(patch), std::move(projector), PatchKeepsParameter(basis)});
}

} // namespace Detail

Result<TrimmedSurface> MakeTrimmedSurface(const Surface& basis, const Box2 range)
{
    if (basis.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "MakeTrimmedSurface: the basis surface is empty");
    }
    const Box2 domain = Domain(basis);
    const Vec2 periods = Periods(basis);
    STEREON_TRY_VOID(CheckTrim(range.X, domain.X, periods.X, "u"));
    STEREON_TRY_VOID(CheckTrim(range.Y, domain.Y, periods.Y, "v"));
    return TrimmedSurface(basis, range);
}

Point3 Eval(const TrimmedSurface& surface, const Point2& uv)
{
    Detail::RequireFinite(uv, "Eval(TrimmedSurface)", "a non-finite parameter");
    Detail::RequireInDomain(surface.GetRange().X, uv.X, 0.0, "Eval(TrimmedSurface)");
    Detail::RequireInDomain(surface.GetRange().Y, uv.Y, 0.0, "Eval(TrimmedSurface)");
    return Eval(surface.GetBasis(), surface.BasisParameter(uv));
}

Vec3 Deriv(const TrimmedSurface& surface, const Point2& uv, const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(TrimmedSurface)", "a non-finite parameter");
    Detail::RequireInDomain(surface.GetRange().X, uv.X, 0.0, "Deriv(TrimmedSurface)");
    Detail::RequireInDomain(surface.GetRange().Y, uv.Y, 0.0, "Deriv(TrimmedSurface)");
    const Vec3 d = Deriv(surface.GetBasis(), surface.BasisParameter(uv), du, dv);
    return surface.IsReversed() && du % 2 == 1 ? -d : d;
}

Box2 Domain(const TrimmedSurface& surface)
{
    return surface.GetRange();
}

Box3 BoundingBox(const TrimmedSurface& surface)
{
    const Detail::TrimmedPatch* patch = surface.GetPatch();
    return patch != nullptr ? BoundingBox(patch->Surface) : BoundingBox(surface.GetBasis());
}

Point2 Project(const TrimmedSurface& surface, const Point3& p)
{
    Detail::RequireFinite(p, "Project(TrimmedSurface)", "a non-finite point");
    return ProjectTrimmed(surface, p);
}

TrimmedSurface Reversed(const TrimmedSurface& surface)
{
    return surface.Reversed();
}

namespace {

/// The factors by which a similarity of the scale rescales the basis's
/// parameters, as the basis's own Transformed does: a plane's both, each a
/// distance along an axis; a cylinder's and a cone's v, a distance along
/// the axis; an extrusion's u and a revolution's v as their profiles'; a
/// trimmed basis's as its own basis's; angles, knots and every other
/// builtin's not at all. A user's type is taken to keep its parameters.
[[nodiscard]] Vec2 ParameterScalesOf(const Surface& basis, const double scale)
{
    if (basis.Target<Plane>() != nullptr) {
        return {scale, scale};
    }
    if (basis.Target<Cylinder>() != nullptr || basis.Target<Cone>() != nullptr) {
        return {1.0, scale};
    }
    if (const ExtrusionSurface* extrusion = basis.Target<ExtrusionSurface>()) {
        return {Detail::ParameterScaleOf(extrusion->GetProfile(), scale), 1.0};
    }
    if (const RevolutionSurface* revolution = basis.Target<RevolutionSurface>()) {
        return {1.0, Detail::ParameterScaleOf(revolution->GetProfile(), scale)};
    }
    if (const TrimmedSurface* trimmed = basis.Target<TrimmedSurface>()) {
        return ParameterScalesOf(trimmed->GetBasis(), scale);
    }
    return {1.0, 1.0};
}

[[nodiscard]] Interval Scaled(const Interval& range, const double k)
{
    return {k * range.GetLower(), k * range.GetUpper()};
}

} // namespace

TrimmedSurface Transformed(const TrimmedSurface& surface, const Transform& x)
{
    const Vec2 k = ParameterScalesOf(surface.GetBasis(), Detail::SimilarityScaleOf(x));
    const Box2 range = surface.GetRange();
    return {Transformed(surface.GetBasis(), x), Box2{Scaled(range.X, k.X), Scaled(range.Y, k.Y)},
            surface.IsReversed()};
}

Result<BSplineSurface> ToBSpline(const TrimmedSurface& surface)
{
    const Detail::TrimmedPatch* patch = surface.GetPatch();
    if (patch == nullptr) {
        return Fail(ErrorCode::NotSupported,
                    "ToBSpline(TrimmedSurface): the basis's type has no exact NURBS form the "
                    "kernel knows");
    }
    return surface.IsReversed() ? Reversed(patch->Surface) : patch->Surface;
}

Result<BSplineSurface> ToBSpline(const TrimmedSurface& surface, const Box2 range)
{
    const Box2 domain = surface.GetRange();
    if (!Contains(domain, range) || !(range.X.GetLower() < range.X.GetUpper()) ||
        !(range.Y.GetLower() < range.Y.GetUpper())) {
        return Fail(ErrorCode::InvalidInput,
                    "ToBSpline(TrimmedSurface): the range must be non-empty and lie in the "
                    "trimmed surface's domain");
    }
    const Box2 basisRange = surface.IsReversed()
                                ? Box2{Interval{surface.Reflect(range.X.GetUpper()),
                                                surface.Reflect(range.X.GetLower())},
                                       range.Y}
                                : range;
    STEREON_TRY(BSplineSurface piece, ToBSpline(surface.GetBasis(), basisRange));
    return surface.IsReversed() ? Reversed(piece) : piece;
}

} // namespace Stereon
