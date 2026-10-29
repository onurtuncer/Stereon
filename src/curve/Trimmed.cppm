// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>

export module Stereon.Curve:Trimmed;

import Stereon.Robust;
import :BSplineCurve;
import :Concepts;
import :Erased;

/// Trimmed curves: a curve restricted to a parameter range, such as a line
/// segment or a circular arc (STEP's trimmed_curve, OCCT's
/// Geom_TrimmedCurve, IGES entities 100, 104 and 110).
///
/// A TrimmedCurve wraps any curve in an erased Curve (or Curve2d) and keeps
/// its parametrisation: Eval(trimmed, t) is Eval(basis, t) for t in the
/// range, so a parameter means the same on both. Reversing it runs the
/// range backwards, t -> lower + upper - t, rather than reversing the basis,
/// whose own parametrisation under Reversed is its type's business.
///
/// Projection must find the closest point within the range, not the basis's
/// global closest point clamped to it: an arc's closest point to a point
/// near its chord is an interior point the basis's projection never
/// returns. For the built-in types the trimmed curve keeps the exact NURBS
/// form of its piece (ToBSpline over the range; a plane curve lifted to
/// z = 0) and projects onto that with the certified B-spline projection,
/// then maps the result back to the basis's parameter. For a user's curve
/// type, which has no known NURBS form, Project compares the basis's global
/// projection, brought into the range, with the two ends: correct when the
/// closest point on the basis lies in the range or at an end, approximate
/// otherwise.

export namespace Stereon::Detail {

/// The exact NURBS form of a trimmed piece in space and its projector,
/// shared by the copies of one trimmed curve. Defined in Trimmed.cpp.
struct TrimmedPiece;

/// The piece of a basis over a range, or null when the basis's type has no
/// exact NURBS form the kernel knows.
[[nodiscard]] std::shared_ptr<const TrimmedPiece> MakePiece(const Curve& basis, Interval range);
[[nodiscard]] std::shared_ptr<const TrimmedPiece> MakePiece(const Curve2d& basis, Interval range);

/// The scale of a similarity, the factor by which it multiplies lengths:
/// the cube root of its determinant's magnitude.
[[nodiscard]] double SimilarityScaleOf(const Transform& x);

/// The factor by which a similarity of the scale rescales the basis's
/// parameter, as the basis's own Transformed does: a parabola's, a distance
/// along its axis, by the scale; a trimmed basis's as its own basis's;
/// every other builtin's not at all, an angle, a knot, or a line's distance
/// along a direction that is carried unscaled. A user's type is taken to
/// keep its parameter. For Transformed(TrimmedCurve) and the trimmed
/// surfaces swept from a profile.
[[nodiscard]] double ParameterScaleOf(const Curve& basis, double scale);

} // namespace Stereon::Detail

export namespace Stereon {

template <std::size_t D>
class BasicTrimmedCurve {
public:
    /// The basis over `range`, reversed or not. Unchecked: the range must be
    /// finite and non-empty, and inside the basis's domain or, for a
    /// periodic basis, at most one period long. MakeTrimmedCurve checks data
    /// from outside the kernel.
    BasicTrimmedCurve(BasicErasedCurve<D> basis, const Interval range, const bool reversed = false)
        : m_Basis(std::move(basis))
        , m_Range(range)
        , m_Reversed(reversed)
        , m_Piece(Detail::MakePiece(m_Basis, m_Range))
    {}

    /// The basis, borrowed from this trimmed curve: lvalues only (ADR-0046).
    [[nodiscard]] const BasicErasedCurve<D>& GetBasis() const& noexcept { return m_Basis; }
    const BasicErasedCurve<D>& GetBasis() const&& = delete;
    [[nodiscard]] Interval GetRange() const noexcept { return m_Range; }
    [[nodiscard]] bool IsReversed() const noexcept { return m_Reversed; }
    /// The exact NURBS piece, or null for a user's basis type.
    [[nodiscard]] const Detail::TrimmedPiece* GetPiece() const noexcept { return m_Piece.get(); }

    /// The basis's parameter at this curve's parameter t.
    [[nodiscard]] double BasisParameter(const double t) const noexcept
    {
        return m_Reversed ? Reflect(t) : t;
    }

    /// lower + upper - t, the parameter of the same point when the sense is
    /// flipped: exact at the ends, and kept in the range where rounding the
    /// sum would leave it.
    [[nodiscard]] double Reflect(const double t) const noexcept
    {
        const double lower = m_Range.GetLower();
        const double upper = m_Range.GetUpper();
        if (t == lower) {
            return upper;
        }
        if (t == upper) {
            return lower;
        }
        return std::clamp(lower + upper - t, lower, upper);
    }

    /// The same piece run the other way. The piece is shared, not rebuilt.
    [[nodiscard]] BasicTrimmedCurve Reversed() const
    {
        BasicTrimmedCurve copy = *this;
        copy.m_Reversed = !m_Reversed;
        return copy;
    }

private:
    BasicErasedCurve<D> m_Basis;
    Interval m_Range;
    bool m_Reversed = false;
    std::shared_ptr<const Detail::TrimmedPiece> m_Piece;
};

using TrimmedCurve = BasicTrimmedCurve<3>;
using TrimmedCurve2d = BasicTrimmedCurve<2>;

/// The trimmed curve after validating data from outside the kernel
/// (ADR-0006): InvalidInput unless the basis is not empty, the range is
/// finite and non-empty, and it lies in the basis's domain or, for a
/// periodic basis, spans at most one period (an arc may cross the seam:
/// a circle trimmed to [5, 7] is fine).
[[nodiscard]] Result<TrimmedCurve> MakeTrimmedCurve(const Curve& basis, Interval range);
[[nodiscard]] Result<TrimmedCurve2d> MakeTrimmedCurve(const Curve2d& basis, Interval range);

template <std::size_t D>
[[nodiscard]] BasicPointN<double, D> Eval(const BasicTrimmedCurve<D>& curve, const double t)
{
    Detail::RequireFinite(t, "Eval(TrimmedCurve)", "a non-finite parameter");
    Detail::RequireInDomain(curve.GetRange(), t, 0.0, "Eval(TrimmedCurve)");
    return Eval(curve.GetBasis(), curve.BasisParameter(t));
}

/// The n-th derivative; odd orders change sign on a reversed curve.
template <std::size_t D>
[[nodiscard]] BasicVecN<double, D> Deriv(const BasicTrimmedCurve<D>& curve, const double t,
                                         const int n)
{
    Detail::RequireFinite(t, "Deriv(TrimmedCurve)", "a non-finite parameter");
    Detail::RequireInDomain(curve.GetRange(), t, 0.0, "Deriv(TrimmedCurve)");
    const BasicVecN<double, D> d = Deriv(curve.GetBasis(), curve.BasisParameter(t), n);
    return curve.IsReversed() && n % 2 == 1 ? -d : d;
}

template <std::size_t D>
[[nodiscard]] Interval Domain(const BasicTrimmedCurve<D>& curve)
{
    return curve.GetRange();
}

template <std::size_t D>
[[nodiscard]] BasicTrimmedCurve<D> Reversed(const BasicTrimmedCurve<D>& curve)
{
    return curve.Reversed();
}

/// A box containing the piece: its NURBS form's pole box, or the basis's
/// box for a user's type.
[[nodiscard]] Box3 BoundingBox(const TrimmedCurve& curve);
[[nodiscard]] Box2 BoundingBox(const TrimmedCurve2d& curve);

/// The parameter of a closest point of the piece (see the comment above for
/// a user's basis type). Ties go to the smaller parameter (ADR-0007).
[[nodiscard]] double Project(const TrimmedCurve& curve, const Point3& p);
[[nodiscard]] double Project(const TrimmedCurve2d& curve, const Point2& p);

/// The trimmed curve moved: the basis transformed, and the range rescaled
/// as the basis's Transformed rescales its parameter under a similarity
/// (Detail::ParameterScaleOf), so that the piece is the image of the piece.
[[nodiscard]] TrimmedCurve Transformed(const TrimmedCurve& curve, const Transform& x);

/// The exact NURBS form of the piece, run as the trimmed curve runs.
/// NotSupported for a user's basis type.
[[nodiscard]] Result<BSplineCurve> ToBSpline(const TrimmedCurve& curve);

/// The exact NURBS form of the part of the trimmed curve over `range`, a
/// subrange of its own, run as the trimmed curve runs, with domain `range`
/// and the trimmed curve's parameter at every knot (everywhere when
/// Detail::KeepsParameter holds). A reversed curve's form is its basis's
/// piece run backwards with its knots reflected into `range`. A trimmed
/// basis converts the same way, to any depth. InvalidInput unless the range
/// is finite, non-empty and inside the trimmed curve's; NotSupported when
/// the innermost basis has no exact NURBS form the kernel knows.
[[nodiscard]] Result<BSplineCurve> ToBSpline(const TrimmedCurve& curve, Interval range);

} // namespace Stereon
