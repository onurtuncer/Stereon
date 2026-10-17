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
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineCurve;
import :Ellipse;
import :Erased;
import :Hyperbola;
import :Interchange;

namespace Stereon {

namespace {

// --- Piecewise polynomials ---------------------------------------------------

/// The binomial coefficients C(n, 0..n), exact in double for any degree a
/// curve has.
std::vector<double> Binomials(const std::size_t n)
{
    std::vector<double> row(n + 1, 1.0);
    for (std::size_t k = 1; k < n; ++k) {
        row[k] = row[k - 1] * static_cast<double>(n - k + 1) / static_cast<double>(k);
    }
    return row;
}

/// The Bernstein coefficients of sum_j a_j s^j on [0, h], degree n:
/// with s = h w, b_k = sum_{j <= k} C(k, j) / C(n, j) a_j h^j.
template <std::size_t D>
std::vector<BasicVecN<double, D>> BernsteinOf(const std::vector<BasicVecN<double, D>>& power,
                                              const double h, const std::size_t n)
{
    std::vector<BasicVecN<double, D>> scaled(n + 1, BasicVecN<double, D>{});
    double hj = 1.0;
    for (std::size_t j = 0; j < power.size(); ++j) {
        scaled[j] = hj * power[j];
        hj *= h;
    }
    const std::vector<double> ofN = Binomials(n);
    std::vector<BasicVecN<double, D>> bernstein(n + 1, BasicVecN<double, D>{});
    for (std::size_t k = 0; k <= n; ++k) {
        const std::vector<double> ofK = Binomials(k);
        for (std::size_t j = 0; j <= k; ++j) {
            bernstein[k] += (ofK[j] / ofN[j]) * scaled[j];
        }
    }
    return bernstein;
}

template <std::size_t D>
Result<BasicPiecewiseConversion<D>>
PiecewiseToBSpline(const std::span<const double> breakpoints,
                   const std::span<const std::vector<BasicVecN<double, D>>> pieces)
{
    constexpr std::string_view where = "PiecewisePolynomialToBSpline";
    if (pieces.empty()) {
        return Fail(ErrorCode::InvalidInput, "PiecewisePolynomialToBSpline: there are no pieces");
    }
    if (breakpoints.size() != pieces.size() + 1) {
        return Fail(ErrorCode::InvalidInput,
                    "PiecewisePolynomialToBSpline: there must be one more breakpoint than pieces");
    }
    for (std::size_t i = 0; i < breakpoints.size(); ++i) {
        STEREON_TRY_VOID(Detail::CheckFinite(breakpoints[i], "a breakpoint", where));
        if (i > 0 && !(breakpoints[i - 1] < breakpoints[i])) {
            return Fail(ErrorCode::InvalidInput,
                        "PiecewisePolynomialToBSpline: the breakpoints must increase strictly");
        }
    }
    std::size_t degree = 1;
    for (const std::vector<BasicVecN<double, D>>& piece : pieces) {
        if (piece.empty()) {
            return Fail(ErrorCode::InvalidInput,
                        "PiecewisePolynomialToBSpline: a piece has no coefficients");
        }
        for (const BasicVecN<double, D>& a : piece) {
            STEREON_TRY_VOID(Detail::CheckFinite(a, "a coefficient", where));
        }
        degree = std::max(degree, piece.size() - 1);
    }

    std::vector<double> knots(degree + 1, breakpoints.front());
    for (std::size_t i = 1; i + 1 < breakpoints.size(); ++i) {
        knots.insert(knots.end(), degree, breakpoints[i]);
    }
    knots.insert(knots.end(), degree + 1, breakpoints.back());

    using PointD = BasicPointN<double, D>;
    std::vector<PointD> poles;
    poles.reserve(pieces.size() * degree + 1);
    double maxGap = 0.0;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        const double h = breakpoints[i + 1] - breakpoints[i];
        const std::vector<BasicVecN<double, D>> b = BernsteinOf<D>(pieces[i], h, degree);
        const PointD start = PointD{} + b[0];
        if (i == 0) {
            poles.push_back(start);
        } else {
            // The shared pole of a join: the mean of where the pieces meet.
            maxGap = std::max(maxGap, Distance(poles.back(), start));
            poles.back() = poles.back() + 0.5 * (start - poles.back());
        }
        for (std::size_t k = 1; k <= degree; ++k) {
            poles.push_back(PointD{} + b[k]);
        }
    }
    STEREON_TRY(auto curve,
                MakeBSplineCurve(static_cast<int>(degree), std::move(knots), std::move(poles), {}));
    return BasicPiecewiseConversion<D>{std::move(curve), maxGap};
}

// --- Implicit conics -----------------------------------------------------------

constexpr double Largest(const ImplicitConic& k)
{
    return std::max(
        {std::abs(k.A), std::abs(k.B), std::abs(k.C), std::abs(k.D), std::abs(k.E), std::abs(k.F)});
}

Result<Curve2d> NotACurve(const std::string_view what)
{
    std::string message = "ConicFromImplicit: the equation describes ";
    message += what;
    message += ", not a conic the kernel can hold";
    return Fail(ErrorCode::InvalidInput, std::move(message));
}

/// The parabola of l1 x'^2 + d1 x' + d2 y' + F = 0 in the eigenframe
/// (e1, e2), the second eigenvalue taken as zero.
Result<Curve2d> ParabolaFrom(const double l1, const Vec2& e1, const Vec2& e2, const double d1,
                             const double d2, const double f)
{
    // The coefficients are normalised to a largest magnitude of 1, so the
    // linear term along the axis is compared with 1 (docs/interchange.rst).
    if (std::abs(d2) <= ConicAllowance) {
        return NotACurve("two parallel lines, one line or no point");
    }
    // l1 (x' - x0)^2 + d2 (y' - y0) = 0.
    const double x0 = -d1 / (2.0 * l1);
    const double y0 = (d1 * d1 / (4.0 * l1) - f) / d2;
    const double k = -l1 / d2;
    const Point2 vertex = Point2{} + x0 * e1 + y0 * e2;
    const Vec2 x = k > 0.0 ? e2 : -e2;
    STEREON_TRY(auto parabola,
                MakeParabola(vertex, std::abs(d2) / (4.0 * std::abs(l1)), x, Perpendicular(x)));
    return Curve2d(parabola);
}

} // namespace

Result<PiecewiseConversion>
PiecewisePolynomialToBSpline(const std::span<const double> breakpoints,
                             const std::span<const std::vector<Vec3>> pieces)
{
    return PiecewiseToBSpline<3>(breakpoints, pieces);
}

Result<PiecewiseConversion2d>
PiecewisePolynomialToBSpline(const std::span<const double> breakpoints,
                             const std::span<const std::vector<Vec2>> pieces)
{
    return PiecewiseToBSpline<2>(breakpoints, pieces);
}

double ImplicitResidual(const ImplicitConic& k, const Point2& p)
{
    const double scale = Largest(k);
    if (!(scale > 0.0)) {
        return 0.0;
    }
    const double x = p.X;
    const double y = p.Y;
    return (k.A * x * x + k.B * x * y + k.C * y * y + k.D * x + k.E * y + k.F) / scale;
}

Result<Curve2d> ConicFromImplicit(const ImplicitConic& conic, const std::optional<Point2> branch)
{
    constexpr std::string_view where = "ConicFromImplicit";
    for (const double c : {conic.A, conic.B, conic.C, conic.D, conic.E, conic.F}) {
        STEREON_TRY_VOID(Detail::CheckFinite(c, "a coefficient", where));
    }
    if (branch) {
        STEREON_TRY_VOID(Detail::CheckFinite(*branch, "the branch point", where));
    }
    const double scale = Largest(conic);
    if (scale == 0.0) {
        return NotACurve("all of the plane");
    }
    const double a = conic.A / scale;
    const double b = conic.B / scale;
    const double c = conic.C / scale;
    const double d = conic.D / scale;
    const double e = conic.E / scale;
    const double f = conic.F / scale;

    // The quadratic part x^T M x, M = [[a, b/2], [b/2, c]], in its
    // eigenframe: the rotation by theta with tan(2 theta) = b / (a - c).
    const double theta = 0.5 * std::atan2(b, a - c);
    const double cs = std::cos(theta);
    const double sn = std::sin(theta);
    Vec2 e1{cs, sn};
    Vec2 e2{-sn, cs};
    double l1 = a * cs * cs + b * sn * cs + c * sn * sn;
    double l2 = a * sn * sn - b * sn * cs + c * cs * cs;
    const double largest = std::max(std::abs(l1), std::abs(l2));
    if (largest == 0.0) {
        return NotACurve(d == 0.0 && e == 0.0 ? "no point" : "a line");
    }
    if (std::abs(l1) < std::abs(l2)) {
        std::swap(l1, l2);
        std::swap(e1, e2);
    }
    // The linear part in the eigenframe.
    const double g1 = d * e1.X + e * e1.Y;
    const double g2 = d * e2.X + e * e2.Y;

    if (std::abs(l2) <= ConicAllowance * largest) {
        return ParabolaFrom(l1, e1, e2, g1, g2, f);
    }

    // A central conic: l1 x'^2 + l2 y'^2 + f' = 0 about the centre.
    const double cx = -g1 / (2.0 * l1);
    const double cy = -g2 / (2.0 * l2);
    const Point2 center = Point2{} + cx * e1 + cy * e2;
    const double shifted = f + 0.5 * (g1 * cx + g2 * cy);
    const double cancelled =
        std::max({std::abs(f), std::abs(0.5 * g1 * cx), std::abs(0.5 * g2 * cy)});
    const bool ellipse = (l1 > 0.0) == (l2 > 0.0);
    if (std::abs(shifted) <= ConicAllowance * cancelled) {
        return NotACurve(ellipse ? "a single point" : "two crossing lines");
    }
    if (ellipse) {
        if ((shifted > 0.0) == (l1 > 0.0)) {
            return NotACurve("no point");
        }
        // |l1| >= |l2|, so the radius along e2 is the larger.
        const double major = std::sqrt(-shifted / l2);
        const double minor = std::sqrt(-shifted / l1);
        STEREON_TRY(auto result, MakeEllipse(center, major, minor, e2, Perpendicular(e2)));
        return Curve2d(result);
    }
    // The transverse axis is the one along which -f' / l is positive.
    const bool alongFirst = -shifted / l1 > 0.0;
    const double transverse = alongFirst ? l1 : l2;
    const double conjugate = alongFirst ? l2 : l1;
    Vec2 x = alongFirst ? e1 : e2;
    if (branch && Dot(*branch - center, x) < 0.0) {
        x = -x;
    }
    STEREON_TRY(auto result, MakeHyperbola(center, std::sqrt(-shifted / transverse),
                                           std::sqrt(shifted / conjugate), x, Perpendicular(x)));
    return Curve2d(result);
}

} // namespace Stereon
