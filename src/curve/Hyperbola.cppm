// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>

export module Stereon.Curve:Hyperbola;

import Stereon.Robust;
import :Concepts;

/// The hyperbola and the parabola, in space and in the parameter plane: the
/// open conics of STEP and OCCT's Geom package.

export namespace Stereon {

/// Center + MajorRadius cosh t * XAxis + MinorRadius sinh t * YAxis for every
/// real t: the branch of the hyperbola (x/a)^2 - (y/b)^2 = 1 on the side of
/// +XAxis, with asymptotes along a XAxis +- b YAxis. XAxis and YAxis are
/// orthonormal; both radii are positive, in either order. The conventions
/// are STEP's hyperbola and OCCT's Geom_Hyperbola.
template <Scalar T, std::size_t D>
struct BasicHyperbola {
    BasicPointN<T, D> Center;
    T MajorRadius;
    T MinorRadius;
    BasicVecN<T, D> XAxis;
    BasicVecN<T, D> YAxis;

    constexpr bool operator==(const BasicHyperbola&) const = default;
};

using Hyperbola = BasicHyperbola<double, 3>;
using Hyperbola2d = BasicHyperbola<double, 2>;

/// Vertex + t^2 / (4 FocalLength) * XAxis + t * YAxis for every real t: the
/// parabola y^2 = 4 f x in its frame, opening towards +XAxis, with its focus
/// at Vertex + FocalLength * XAxis. XAxis and YAxis are orthonormal;
/// FocalLength is positive. The conventions are STEP's parabola and OCCT's
/// Geom_Parabola.
template <Scalar T, std::size_t D>
struct BasicParabola {
    BasicPointN<T, D> Vertex;
    T FocalLength;
    BasicVecN<T, D> XAxis;
    BasicVecN<T, D> YAxis;

    constexpr bool operator==(const BasicParabola&) const = default;
};

using Parabola = BasicParabola<double, 3>;
using Parabola2d = BasicParabola<double, 2>;

/// The conics after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless every value is finite, the radii or the focal length
/// are positive and the axes are orthonormal within FrameAllowance.
[[nodiscard]] Result<Hyperbola> MakeHyperbola(const Point3& center, double majorRadius,
                                              double minorRadius, const Vec3& xAxis,
                                              const Vec3& yAxis);
[[nodiscard]] Result<Hyperbola2d> MakeHyperbola(const Point2& center, double majorRadius,
                                                double minorRadius, const Vec2& xAxis,
                                                const Vec2& yAxis);
[[nodiscard]] Result<Parabola> MakeParabola(const Point3& vertex, double focalLength,
                                            const Vec3& xAxis, const Vec3& yAxis);
[[nodiscard]] Result<Parabola2d> MakeParabola(const Point2& vertex, double focalLength,
                                              const Vec2& xAxis, const Vec2& yAxis);

// The hyperbola.

template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicPointN<T, D> Eval(const BasicHyperbola<T, D>& hyperbola, const T& t)
{
    Detail::RequireFinite(t, "Eval(Hyperbola)", "a non-finite parameter");
    return hyperbola.Center + (hyperbola.MajorRadius * Cosh(t)) * hyperbola.XAxis +
           (hyperbola.MinorRadius * Sinh(t)) * hyperbola.YAxis;
}

/// cosh and sinh swap under differentiation: even orders give
/// a cosh t X + b sinh t Y, odd orders a sinh t X + b cosh t Y.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicVecN<T, D> Deriv(const BasicHyperbola<T, D>& hyperbola, const T& t,
                                              const int n)
{
    Detail::RequireFinite(t, "Deriv(Hyperbola)", "a non-finite parameter");
    Detail::Require(n >= 1, "Deriv(Hyperbola)", "a derivative order below 1");
    const T c = Cosh(t);
    const T s = Sinh(t);
    const bool odd = n % 2 == 1;
    return (hyperbola.MajorRadius * (odd ? s : c)) * hyperbola.XAxis +
           (hyperbola.MinorRadius * (odd ? c : s)) * hyperbola.YAxis;
}

/// The whole real line.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain([[maybe_unused]] const BasicHyperbola<T, D>& hyperbola)
{
    return Interval::Everything();
}

/// Everything: a hyperbola is unbounded.
template <std::size_t D>
[[nodiscard]] constexpr BoxN<D> BoundingBox([[maybe_unused]] const BasicHyperbola<double, D>& h)
{
    return BoxN<D>::Everything();
}

/// The parameter of a point of the hyperbola closest to p. In the frame,
/// with p's in-plane coordinates (x, y), the stationary points of the
/// squared distance satisfy (a^2 + b^2) sinh t cosh t = a x sinh t + b y
/// cosh t; with s = e^t this is the quartic
/// (a^2 + b^2) s^4 - 2 (a x + b y) s^3 + 2 (a x - b y) s - (a^2 + b^2) = 0,
/// whose positive roots Detail::PolynomialRoots isolates on a range that
/// must contain the closest point (docs/analytic.rst).
[[nodiscard]] double Project(const Hyperbola& hyperbola, const Point3& p);
[[nodiscard]] double Project(const Hyperbola2d& hyperbola, const Point2& p);

/// The same branch traversed the other way: Eval(Reversed(h), -t) ==
/// Eval(h, t).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicHyperbola<T, D> Reversed(const BasicHyperbola<T, D>& hyperbola)
{
    return {hyperbola.Center, hyperbola.MajorRadius, hyperbola.MinorRadius, hyperbola.XAxis,
            -hyperbola.YAxis};
}

/// The image under a similarity; the parameter is kept and the radii scale.
template <Scalar T>
[[nodiscard]] constexpr BasicHyperbola<T, 3> Transformed(const BasicHyperbola<T, 3>& hyperbola,
                                                         const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(hyperbola.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(hyperbola.YAxis, x);
    const T scale = Length(xAxis);
    return {Transformed(hyperbola.Center, x), hyperbola.MajorRadius * scale,
            hyperbola.MinorRadius * scale, xAxis / scale, yAxis / scale};
}

// The parabola.

template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicPointN<T, D> Eval(const BasicParabola<T, D>& parabola, const T& t)
{
    Detail::RequireFinite(t, "Eval(Parabola)", "a non-finite parameter");
    return parabola.Vertex + (t * t / (T{4.0} * parabola.FocalLength)) * parabola.XAxis +
           t * parabola.YAxis;
}

/// t / (2f) X + Y for n = 1, X / (2f) for n = 2, zero above.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicVecN<T, D> Deriv(const BasicParabola<T, D>& parabola, const T& t,
                                              const int n)
{
    Detail::RequireFinite(t, "Deriv(Parabola)", "a non-finite parameter");
    Detail::Require(n >= 1, "Deriv(Parabola)", "a derivative order below 1");
    const T twiceFocal = T{2.0} * parabola.FocalLength;
    if (n == 1) {
        return (t / twiceFocal) * parabola.XAxis + parabola.YAxis;
    }
    return n == 2 ? parabola.XAxis / twiceFocal : BasicVecN<T, D>{};
}

/// The whole real line.
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr Interval Domain([[maybe_unused]] const BasicParabola<T, D>& parabola)
{
    return Interval::Everything();
}

/// Everything: a parabola is unbounded.
template <std::size_t D>
[[nodiscard]] constexpr BoxN<D> BoundingBox([[maybe_unused]] const BasicParabola<double, D>& p)
{
    return BoxN<D>::Everything();
}

/// The parameter of a point of the parabola closest to p: with p's
/// in-plane coordinates (x, y), the real roots of the cubic
/// t^3 / (8 f^2) + t (1 - x / (2f)) - y = 0, isolated on a range that must
/// contain the closest point; the nearest wins, the smaller parameter on a
/// tie (docs/analytic.rst).
[[nodiscard]] double Project(const Parabola& parabola, const Point3& p);
[[nodiscard]] double Project(const Parabola2d& parabola, const Point2& p);

/// The same parabola traversed the other way: Eval(Reversed(c), -t) ==
/// Eval(c, t).
template <Scalar T, std::size_t D>
[[nodiscard]] constexpr BasicParabola<T, D> Reversed(const BasicParabola<T, D>& parabola)
{
    return {parabola.Vertex, parabola.FocalLength, parabola.XAxis, -parabola.YAxis};
}

/// The image under a similarity with scale k: the focal length becomes
/// k f and the parameter is rescaled to k t, because t is the distance
/// along YAxis.
template <Scalar T>
[[nodiscard]] constexpr BasicParabola<T, 3> Transformed(const BasicParabola<T, 3>& parabola,
                                                        const BasicTransform<T>& x)
{
    const BasicVec3<T> xAxis = Transformed(parabola.XAxis, x);
    const BasicVec3<T> yAxis = Transformed(parabola.YAxis, x);
    const T scale = Length(xAxis);
    return {Transformed(parabola.Vertex, x), parabola.FocalLength * scale, xAxis / scale,
            yAxis / scale};
}

static_assert(CurveGeometry<Hyperbola>);
static_assert(Curve2dGeometry<Hyperbola2d>);
static_assert(CurveGeometry<Parabola>);
static_assert(Curve2dGeometry<Parabola2d>);
static_assert(Detail::FitsInline<Hyperbola, GeometryBufferSize>);
static_assert(Detail::FitsInline<Parabola, GeometryBufferSize>);

} // namespace Stereon
