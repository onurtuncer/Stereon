// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineCurve;
import :Circle;
import :Conversion;
import :Ellipse;
import :Erased;
import :Hyperbola;
import :Line;
import :Trimmed;
import :Visit;

namespace Stereon::Detail {

struct TrimmedPiece {
    /// The piece in space, exactly; a plane curve lifted to z = 0.
    BSplineCurve Curve;
    CurveProjector Projector;
    /// True when the piece is parametrised as the basis (a line, a parabola
    /// or a B-spline); otherwise a point of the piece is mapped back to the
    /// basis by projecting onto it.
    bool SameParameter = false;
};

} // namespace Stereon::Detail

namespace Stereon {

namespace {

// --- Plane curves lifted to z = 0, parametrised as before --------------------

Point3 Lift(const Point2& p)
{
    return {p.X, p.Y, 0.0};
}

Vec3 Lift(const Vec2& v)
{
    return {v.X, v.Y, 0.0};
}

Line Lift(const Line2d& g)
{
    return {Lift(g.Origin), Lift(g.Direction)};
}

Circle Lift(const Circle2d& g)
{
    return {Lift(g.Center), g.Radius, Lift(g.XAxis), Lift(g.YAxis)};
}

Ellipse Lift(const Ellipse2d& g)
{
    return {Lift(g.Center), g.MajorRadius, g.MinorRadius, Lift(g.XAxis), Lift(g.YAxis)};
}

Hyperbola Lift(const Hyperbola2d& g)
{
    return {Lift(g.Center), g.MajorRadius, g.MinorRadius, Lift(g.XAxis), Lift(g.YAxis)};
}

Parabola Lift(const Parabola2d& g)
{
    return {Lift(g.Vertex), g.FocalLength, Lift(g.XAxis), Lift(g.YAxis)};
}

BSplineCurve Lift(const BSplineCurve2d& g)
{
    BSplineCurve lifted{g.Degree, g.Knots, {}, g.Weights};
    lifted.Poles.reserve(g.Poles.size());
    for (const Point2& p : g.Poles) {
        lifted.Poles.push_back(Lift(p));
    }
    return lifted;
}

/// The basis as a curve in space: itself, a built-in plane curve lifted, or
/// nothing for a user's plane curve type.
std::optional<Curve> InSpace(const Curve& basis)
{
    return basis;
}

std::optional<Curve> InSpace(const Curve2d& basis)
{
    // A trimmed plane curve lifts as its basis does, trimmed the same way.
    if (const auto* trimmed = basis.Target<TrimmedCurve2d>()) {
        std::optional<Curve> inner = InSpace(trimmed->GetBasis());
        if (!inner) {
            return std::nullopt;
        }
        return Curve(TrimmedCurve(*std::move(inner), trimmed->GetRange(), trimmed->IsReversed()));
    }
    return VisitBuiltin(basis, []<typename G>(const G& g) -> std::optional<Curve> {
        if constexpr (std::is_same_v<G, Curve2d>) {
            return std::nullopt;
        } else {
            return Curve(Lift(g));
        }
    });
}

/// A point of the piece in the trimmed curve's dimension.
Point3 InDimension(const Point3& p, const Point3& /*tag*/)
{
    return p;
}

Point2 InDimension(const Point3& p, const Point2& /*tag*/)
{
    return {p.X, p.Y};
}

/// The basis parameter of the piece's closest point to p, or of the
/// basis's global closest point brought into the range for a user type.
template <std::size_t D>
double InteriorCandidate(const BasicTrimmedCurve<D>& curve, const BasicPointN<double, D>& p)
{
    const BasicErasedCurve<D>& basis = curve.GetBasis();
    const Interval range = curve.GetRange();
    const Detail::TrimmedPiece* piece = curve.GetPiece();
    if (piece == nullptr) {
        return NearestInRange(Project(basis, p), range, Period(basis));
    }
    Point3 lifted;
    if constexpr (D == 2) {
        lifted = Lift(p);
    } else {
        lifted = p;
    }
    const double s = piece->Projector.Project(lifted);
    if (piece->SameParameter) {
        return std::clamp(s, range.GetLower(), range.GetUpper());
    }
    // The piece's point lies on the basis, which (a line or a conic) does
    // not cross itself, so projecting it back finds its own parameter.
    const BasicPointN<double, D> onPiece = InDimension(Eval(piece->Curve, s), p);
    return NearestInRange(Project(basis, onPiece), range, Period(basis));
}

template <std::size_t D>
double ProjectTrimmed(const BasicTrimmedCurve<D>& curve, const BasicPointN<double, D>& p)
{
    const Interval range = curve.GetRange();
    const std::array<double, 3> candidates{range.GetLower(), InteriorCandidate(curve, p),
                                           range.GetUpper()};
    double best = candidates[0];
    double bestDistance = DistanceSquared(Eval(curve.GetBasis(), best), p);
    for (std::size_t k = 1; k < candidates.size(); ++k) {
        const double d = DistanceSquared(Eval(curve.GetBasis(), candidates[k]), p);
        if (d < bestDistance) {
            best = candidates[k];
            bestDistance = d;
        }
    }
    return curve.IsReversed() ? curve.Reflect(best) : best;
}

template <std::size_t D>
Result<BasicTrimmedCurve<D>> MakeTrimmed(const BasicErasedCurve<D>& basis, const Interval range)
{
    if (basis.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "MakeTrimmedCurve: the basis curve is empty");
    }
    const double low = range.GetLower();
    const double high = range.GetUpper();
    if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeTrimmedCurve: the range must be finite and non-empty");
    }
    const double period = Period(basis);
    if (period > 0.0 ? high - low > period : !Contains(Domain(basis), range)) {
        return Fail(ErrorCode::InvalidInput,
                    period > 0.0 ? "MakeTrimmedCurve: the range spans more than one period"
                                 : "MakeTrimmedCurve: the range must lie in the basis's domain");
    }
    return BasicTrimmedCurve<D>(basis, range);
}

} // namespace

namespace Detail {

template <std::size_t D>
std::shared_ptr<const TrimmedPiece> MakePieceFor(const BasicErasedCurve<D>& basis,
                                                 const Interval range)
{
    if (basis.IsEmpty()) {
        return nullptr;
    }
    const std::optional<Curve> space = InSpace(basis);
    if (!space) {
        return nullptr;
    }
    Result<BSplineCurve> piece = ToBSpline(*space, range);
    if (!piece) {
        return nullptr;
    }
    CurveProjector projector(*piece);
    return std::make_shared<const TrimmedPiece>(
        TrimmedPiece{*std::move(piece), std::move(projector), Detail::KeepsParameter(basis)});
}

std::shared_ptr<const TrimmedPiece> MakePiece(const Curve& basis, const Interval range)
{
    return MakePieceFor(basis, range);
}

std::shared_ptr<const TrimmedPiece> MakePiece(const Curve2d& basis, const Interval range)
{
    return MakePieceFor(basis, range);
}

} // namespace Detail

Result<TrimmedCurve> MakeTrimmedCurve(const Curve& basis, const Interval range)
{
    return MakeTrimmed(basis, range);
}

Result<TrimmedCurve2d> MakeTrimmedCurve(const Curve2d& basis, const Interval range)
{
    return MakeTrimmed(basis, range);
}

Box3 BoundingBox(const TrimmedCurve& curve)
{
    const Detail::TrimmedPiece* piece = curve.GetPiece();
    return piece != nullptr ? BoundingBox(piece->Curve) : BoundingBox(curve.GetBasis());
}

Box2 BoundingBox(const TrimmedCurve2d& curve)
{
    const Detail::TrimmedPiece* piece = curve.GetPiece();
    if (piece == nullptr) {
        return BoundingBox(curve.GetBasis());
    }
    const Box3 box = BoundingBox(piece->Curve);
    return {box.X, box.Y};
}

double Project(const TrimmedCurve& curve, const Point3& p)
{
    Detail::RequireFinite(p, "Project(TrimmedCurve)", "a non-finite point");
    return ProjectTrimmed(curve, p);
}

double Project(const TrimmedCurve2d& curve, const Point2& p)
{
    Detail::RequireFinite(p, "Project(TrimmedCurve2d)", "a non-finite point");
    return ProjectTrimmed(curve, p);
}

TrimmedCurve Transformed(const TrimmedCurve& curve, const Transform& x)
{
    return {Transformed(curve.GetBasis(), x), curve.GetRange(), curve.IsReversed()};
}

Result<BSplineCurve> ToBSpline(const TrimmedCurve& curve)
{
    const Detail::TrimmedPiece* piece = curve.GetPiece();
    if (piece == nullptr) {
        return Fail(ErrorCode::NotSupported,
                    "ToBSpline(TrimmedCurve): the basis's type has no exact NURBS form the "
                    "kernel knows");
    }
    return curve.IsReversed() ? Reversed(piece->Curve) : piece->Curve;
}

Result<BSplineCurve> ToBSpline(const TrimmedCurve& curve, const Interval range)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    if (!std::isfinite(low) || !std::isfinite(high) || !(low < high) ||
        !Contains(curve.GetRange(), range)) {
        return Fail(ErrorCode::InvalidInput,
                    "ToBSpline(TrimmedCurve): the range must be finite, non-empty and inside the "
                    "trimmed curve's range");
    }
    if (!curve.IsReversed()) {
        // The trimmed curve's parameter is its basis's.
        return ToBSpline(curve.GetBasis(), range);
    }
    // Reversed, the curve at t is the basis at Reflect(t), so the part over
    // [low, high] is the basis's over [Reflect(high), Reflect(low)], run
    // backwards. Reversed(piece) would put it on that reflected domain; its
    // knots are reflected into `range` instead, an affine reparametrisation,
    // with the ends exact and the rest kept in the range where rounding
    // would leave it.
    STEREON_TRY(const BSplineCurve piece,
                ToBSpline(curve.GetBasis(), Interval{curve.Reflect(high), curve.Reflect(low)}));
    BSplineCurve form = Reversed(piece);
    const std::size_t n = piece.Knots.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double knot = piece.Knots[n - 1 - i];
        form.Knots[i] = knot == piece.Knots.front()  ? high
                        : knot == piece.Knots.back() ? low
                                                     : std::clamp(curve.Reflect(knot), low, high);
    }
    return form;
}

} // namespace Stereon
