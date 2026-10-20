// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

export module Stereon.Surface:BSplineSurface;

import Stereon.Curve;
import :Concepts;

/// B-spline and NURBS surfaces.

export namespace Stereon {

/// One of a surface's two parameter directions.
enum class ParameterDirection : std::uint8_t {
    U,
    V,
};

/// The tensor-product surface
/// S(u, v) = sum_ij N_i(u) M_j(v) w_ij P_ij / sum_ij N_i(u) M_j(v) w_ij
/// over [KnotsU[p], KnotsU[n+1]] x [KnotsV[q], KnotsV[m+1]], with N_i the
/// degree-p basis of KnotsU and M_j the degree-q basis of KnotsV (see
/// BSplineBasis). Poles and weights are stored row by row: pole (i, j) is
/// Poles[i * PoleCountV(s) + j], i along u. With no weights the surface is
/// polynomial; with weights, all positive, it is a NURBS surface, which can
/// represent the quadrics exactly. The invariants are those
/// MakeBSplineSurface checks.
///
/// As for BasicBSplineCurve, poles and weights may be any Scalar (ADR-0021)
/// while knots and parameters stay double.
template <Scalar T>
struct BasicBSplineSurface {
    int DegreeU = 1;
    int DegreeV = 1;
    std::vector<double> KnotsU;
    std::vector<double> KnotsV;
    std::vector<BasicPoint3<T>> Poles;
    std::vector<T> Weights;

    bool operator==(const BasicBSplineSurface&) const = default;
};

using BSplineSurface = BasicBSplineSurface<double>;

/// The number of poles along u, n + 1 = KnotsU.size() - DegreeU - 1.
template <Scalar T>
[[nodiscard]] constexpr std::size_t PoleCountU(const BasicBSplineSurface<T>& surface)
{
    return surface.KnotsU.size() - static_cast<std::size_t>(surface.DegreeU) - 1;
}

/// The number of poles along v, m + 1 = KnotsV.size() - DegreeV - 1.
template <Scalar T>
[[nodiscard]] constexpr std::size_t PoleCountV(const BasicBSplineSurface<T>& surface)
{
    return surface.KnotsV.size() - static_cast<std::size_t>(surface.DegreeV) - 1;
}

/// True for a NURBS surface, one with weights.
template <Scalar T>
[[nodiscard]] constexpr bool IsRational(const BasicBSplineSurface<T>& surface) noexcept
{
    return !surface.Weights.empty();
}

/// The surface after validating data from outside the kernel (ADR-0006):
/// InvalidInput unless each knot vector suits its degree and pole count
/// (as for MakeBSplineCurve), there are PoleCountU * PoleCountV finite
/// poles, and the weights are absent or one positive finite weight per pole.
[[nodiscard]] Result<BSplineSurface>
MakeBSplineSurface(int degreeU, int degreeV, std::vector<double> knotsU, std::vector<double> knotsV,
                   std::vector<Point3> poles, std::vector<double> weights = {});

/// [KnotsU[p], KnotsU[n+1]] x [KnotsV[q], KnotsV[m+1]].
template <Scalar T>
[[nodiscard]] constexpr Box2 Domain(const BasicBSplineSurface<T>& surface)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    return {Interval{surface.KnotsU[p], surface.KnotsU[PoleCountU(surface)]},
            Interval{surface.KnotsV[q], surface.KnotsV[PoleCountV(surface)]}};
}

} // namespace Stereon

namespace Stereon::Detail {

/// Work arrays for the basis functions of one direction.
inline constexpr std::size_t SurfaceBasisCount = BSplineMaxDegree + 1;

/// The evaluator preconditions (hot-path checks, ADR-0016).
template <Scalar T>
constexpr void CheckSurfaceEvaluation(const BasicBSplineSurface<T>& surface, const Point2& uv,
                                      const std::string_view function)
{
    if constexpr (EvaluatorChecks) {
        const auto p = static_cast<std::size_t>(surface.DegreeU);
        const auto q = static_cast<std::size_t>(surface.DegreeV);
        Require(surface.DegreeU >= 1 && surface.DegreeU <= BSplineMaxDegree &&
                    surface.DegreeV >= 1 && surface.DegreeV <= BSplineMaxDegree,
                function, "an unsupported degree");
        Require(surface.KnotsU.size() > 2 * p + 1 && surface.KnotsV.size() > 2 * q + 1 &&
                    surface.Poles.size() == PoleCountU(surface) * PoleCountV(surface) &&
                    (surface.Weights.empty() || surface.Weights.size() == surface.Poles.size()),
                function, "inconsistent knot, pole or weight counts");
        RequireFinite(uv, function, "a non-finite parameter");
        Require(Contains(Domain(surface), uv), function, "a parameter outside the domain");
    }
}

template <Scalar T>
[[nodiscard]] constexpr T SurfaceWeightOf(const BasicBSplineSurface<T>& surface,
                                          const std::size_t index)
{
    return surface.Weights.empty() ? T{1.0} : surface.Weights[index];
}

/// Small tables live on the stack, larger ones on the heap: the derivative
/// tables of the rational recurrence are (du + 1)(dv + 1) entries, 9 for the
/// second-order partials. Vectors zero themselves on construction, so the
/// capacity is kept small.
inline constexpr std::size_t SurfaceTableCapacity = 16;

/// Basis work arrays for degrees below 8 in both directions, which covers
/// the surfaces met in practice; arrays sized for degree 25 would put 10 KB
/// on the stack per evaluation.
inline constexpr std::size_t SmallSurfaceDegree = 8;

template <typename V>
class SurfaceTable {
public:
    explicit SurfaceTable(const std::size_t size)
    {
        if (size > SurfaceTableCapacity) {
            m_Heap.resize(size);
        }
    }
    [[nodiscard]] V& operator[](const std::size_t i)
    {
        return m_Heap.empty() ? m_Stack[i] : m_Heap[i];
    }

private:
    std::array<V, SurfaceTableCapacity> m_Stack{};
    std::vector<V> m_Heap;
};

[[nodiscard]] constexpr double SurfaceBinomial(const std::size_t n, const std::size_t k)
{
    double c = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        c = c * static_cast<double>(n - k + i) / static_cast<double>(i);
    }
    return c;
}

/// The active window of a surface at a parameter pair: the knot spans and
/// the triangles of basis functions of every degree in each direction,
/// Nu[j * (p + 1) + r] = N_{spanU-j+r,j}(u) for j <= p, r <= j, and Nv
/// likewise (BasisFunctionsOfAllDegrees). The derivatives come from these
/// levels and differenced poles, as for curves, so no basis derivative
/// (algorithm A2.3) is needed. ScaleU[k * (p + 1) + i] is the factor of the
/// (k + 1)-th difference of eq. 3.8, (p - k) / (u_(spanU+i+1) -
/// u_(spanU-p+i+k+1)), and ScaleV likewise: they depend only on the knots,
/// so each is divided once per evaluation, not once per row of poles.
template <std::size_t Capacity>
struct SurfaceWindow {
    std::size_t SpanU = 0;
    std::size_t SpanV = 0;
    std::array<double, Capacity * Capacity> Nu; // filled by BasisFunctionsOfAllDegrees
    std::array<double, Capacity * Capacity> Nv;
    std::array<double, Capacity * Capacity> ScaleU; // filled by DifferenceScales
    std::array<double, Capacity * Capacity> ScaleV;
};

/// The difference factors for orders 1 to `orders` - 1 in a direction of
/// degree n, into scales[k * (n + 1) + i].
template <std::size_t Capacity>
void DifferenceScales(const std::span<const double> knots, const std::size_t span,
                      const std::size_t n, const std::size_t orders,
                      std::array<double, Capacity * Capacity>& scales)
{
    for (std::size_t k = 0; k + 1 < orders; ++k) {
        const std::size_t degree = n - k;
        for (std::size_t i = 0; i < degree; ++i) {
            scales[k * (n + 1) + i] =
                static_cast<double>(degree) / (knots[span + i + 1] - knots[span - n + i + k + 1]);
        }
    }
}

template <std::size_t Capacity, Scalar T>
void FillWindow(const BasicBSplineSurface<T>& surface, const Point2& uv, const std::size_t ku,
                const std::size_t kv, SurfaceWindow<Capacity>& window)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    window.SpanU = FindSpan(surface.KnotsU, surface.DegreeU, PoleCountU(surface), uv.X);
    window.SpanV = FindSpan(surface.KnotsV, surface.DegreeV, PoleCountV(surface), uv.Y);
    BasisFunctionsOfAllDegrees(surface.KnotsU, surface.DegreeU, window.SpanU, uv.X,
                               std::span{window.Nu}.first((p + 1) * (p + 1)));
    BasisFunctionsOfAllDegrees(surface.KnotsV, surface.DegreeV, window.SpanV, uv.Y,
                               std::span{window.Nv}.first((q + 1) * (q + 1)));
    DifferenceScales<Capacity>(surface.KnotsU, window.SpanU, p, ku, window.ScaleU);
    DifferenceScales<Capacity>(surface.KnotsV, window.SpanV, q, kv, window.ScaleV);
}

/// A window of n + 1 homogeneous poles along one direction, differenced in
/// place: after k differences entries 0 to n - k are the poles of the k-th
/// derivative, a B-spline of degree n - k on the same knots (eq. 3.8).
template <std::size_t Capacity, Scalar T>
struct PoleWindow {
    std::array<BasicVec3<T>, Capacity> Points;
    std::array<T, Capacity> Weights; // used by rational surfaces only
};

/// One more difference: the poles of order k to those of order k + 1,
/// Q_i = (n - k) (Q_(i+1) - Q_i) / (t_(span+i+1) - t_(span-n+i+k+1)), over
/// the knots of the direction with degree n, with the factors from
/// DifferenceScales. Every denominator is positive on a non-empty span.
template <bool Rational, std::size_t Capacity, Scalar T>
void DifferenceOnce(const std::array<double, Capacity * Capacity>& scales, const std::size_t n,
                    const std::size_t k, PoleWindow<Capacity, T>& q)
{
    const std::size_t degree = n - k;
    for (std::size_t i = 0; i < degree; ++i) {
        const double scale = scales[k * (n + 1) + i];
        q.Points[i] = T{scale} * (q.Points[i + 1] - q.Points[i]);
        if constexpr (Rational) {
            q.Weights[i] = T{scale} * (q.Weights[i + 1] - q.Weights[i]);
        }
    }
}

/// The combination with one level of a basis triangle of stride n + 1:
/// sum_i N_i q_i over the degree + 1 functions of that level, and the same
/// for the weights.
template <bool Rational, std::size_t Capacity, Scalar T>
void CombineLevel(const std::array<double, Capacity * Capacity>& triangle, const std::size_t n,
                  const std::size_t degree, const PoleWindow<Capacity, T>& q, BasicVec3<T>& sum,
                  T& weight)
{
    sum = {};
    weight = T{};
    for (std::size_t i = 0; i <= degree; ++i) {
        const T basis{triangle[degree * (n + 1) + i]};
        sum += basis * q.Points[i];
        if constexpr (Rational) {
            weight = weight + basis * q.Weights[i];
        }
    }
}

/// The pole of the window with the largest basis product N_r M_s at the
/// point: the evaluation works relative to it. Derivatives do not depend on
/// a translation, but their rounding does: in absolute coordinates a
/// difference of poles near the point carries an error of the order of
/// eps |P|, and where the true derivative is small that error dominates.
/// Next to a collapsed edge, such as the pole of a NURBS sphere, S_u is of
/// the order of the distance d to the edge, so absolute coordinates lose
/// curvature like eps / d^2; relative to the collapsed pole every term
/// shrinks with d as well.
template <std::size_t Capacity, Scalar T>
[[nodiscard]] BasicVec3<T> WindowOrigin(const BasicBSplineSurface<T>& surface,
                                        const SurfaceWindow<Capacity>& window)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    std::size_t r = 0;
    std::size_t c = 0;
    for (std::size_t i = 1; i <= p; ++i) {
        r = window.Nu[p * (p + 1) + i] > window.Nu[p * (p + 1) + r] ? i : r;
    }
    for (std::size_t j = 1; j <= q; ++j) {
        c = window.Nv[q * (q + 1) + j] > window.Nv[q * (q + 1) + c] ? j : c;
    }
    const std::size_t index = (window.SpanU - p + r) * PoleCountV(surface) + (window.SpanV - q + c);
    return ToVec(surface.Poles[index]);
}

/// Stage one of the tensor-product contraction: for each of the p + 1
/// active rows of poles and each v order l < kv, the l-th v derivative of
/// the row's curve, rows[r * kv + l], from the row's poles differenced l
/// times and level q - l of the v triangle; and its weights. The poles are
/// taken relative to `origin` (WindowOrigin). A polynomial surface skips
/// the weights.
template <bool Rational, std::size_t Capacity, std::size_t Rows, Scalar T>
void ContractRows(const BasicBSplineSurface<T>& surface, const SurfaceWindow<Capacity>& window,
                  const BasicVec3<T>& origin, const std::size_t kv,
                  std::array<BasicVec3<T>, Rows>& rows, std::array<T, Rows>& rowWeights)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    PoleWindow<Capacity, T> row;
    for (std::size_t r = 0; r <= p; ++r) {
        const std::size_t first = (window.SpanU - p + r) * countV + (window.SpanV - q);
        for (std::size_t s = 0; s <= q; ++s) {
            if constexpr (Rational) {
                row.Weights[s] = surface.Weights[first + s];
                row.Points[s] = row.Weights[s] * (ToVec(surface.Poles[first + s]) - origin);
            } else {
                row.Points[s] = ToVec(surface.Poles[first + s]) - origin;
            }
        }
        for (std::size_t l = 0; l < kv; ++l) {
            CombineLevel<Rational>(window.Nv, q, q - l, row, rows[r * kv + l],
                                   rowWeights[r * kv + l]);
            if (l + 1 < kv) {
                DifferenceOnce<Rational>(window.ScaleV, q, l, row);
            }
        }
    }
}

/// Stage two: for each v order l, the column rows[. * kv + l] as poles in
/// u, combined with level p - k of the u triangle and differenced between
/// orders: A^(k,l) for k < ku and k + l <= maxTotal, and likewise w^(k,l).
/// A polynomial surface has w = 1 and no weight derivatives.
template <bool Rational, std::size_t Capacity, std::size_t Rows, Scalar T>
void CombineRows(const BasicBSplineSurface<T>& surface, const SurfaceWindow<Capacity>& window,
                 const std::array<BasicVec3<T>, Rows>& rows, const std::array<T, Rows>& rowWeights,
                 const std::size_t ku, const std::size_t kv, const std::size_t maxTotal,
                 SurfaceTable<BasicVec3<T>>& a, SurfaceTable<T>& w)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    PoleWindow<Capacity, T> column;
    for (std::size_t l = 0; l < kv && l <= maxTotal; ++l) {
        for (std::size_t r = 0; r <= p; ++r) {
            column.Points[r] = rows[r * kv + l];
            column.Weights[r] = rowWeights[r * kv + l];
        }
        for (std::size_t k = 0; k < ku && k + l <= maxTotal; ++k) {
            T weight{};
            CombineLevel<Rational>(window.Nu, p, p - k, column, a[k * kv + l], weight);
            w[k * kv + l] = Rational ? weight : T{k == 0 && l == 0 ? 1.0 : 0.0};
            if (k + 1 < ku && k + 1 + l <= maxTotal) {
                DifferenceOnce<Rational>(window.ScaleU, p, k, column);
            }
        }
    }
}

/// The homogeneous derivatives A^(k,l) = sum N_i^(k) M_j^(l) w P and
/// w^(k,l) = sum N_i^(k) M_j^(l) w for k < ku, l < kv and k + l <=
/// maxTotal, into tables indexed k * kv + l: the tensor-product form of
/// A3.6, contracted one direction at a time, which costs about
/// (p + 1)(q + 1) kv + ku kv (p + 1) multiply-adds instead of
/// ku kv (p + 1)(q + 1), with every derivative from differenced poles and
/// the basis triangles rather than from basis derivatives. Rows bounds
/// (p + 1) kv.
template <std::size_t Capacity, std::size_t Rows, Scalar T>
[[nodiscard]] BasicVec3<T>
HomogeneousSurfaceDerivativesWith(const BasicBSplineSurface<T>& surface, const Point2& uv,
                                  const std::size_t ku, const std::size_t kv,
                                  const std::size_t maxTotal, SurfaceTable<BasicVec3<T>>& a,
                                  SurfaceTable<T>& w)
{
    SurfaceWindow<Capacity> window;
    FillWindow(surface, uv, ku, kv, window);
    const BasicVec3<T> origin = WindowOrigin(surface, window);
    std::array<BasicVec3<T>, Rows> rows;
    std::array<T, Rows> rowWeights;
    if (IsRational(surface)) {
        ContractRows<true>(surface, window, origin, kv, rows, rowWeights);
        CombineRows<true>(surface, window, rows, rowWeights, ku, kv, maxTotal, a, w);
    } else {
        ContractRows<false>(surface, window, origin, kv, rows, rowWeights);
        CombineRows<false>(surface, window, rows, rowWeights, ku, kv, maxTotal, a, w);
    }
    return origin;
}

/// Orders up to second in each direction, the common case of Newton steps
/// and curvature, keep the row buffer small.
inline constexpr std::size_t LowSurfaceOrders = 3;

/// The homogeneous derivatives relative to the returned origin: A^(0,0) / w
/// is S - origin, and every derivative is that of S.
template <Scalar T>
[[nodiscard]] BasicVec3<T>
HomogeneousSurfaceDerivatives(const BasicBSplineSurface<T>& surface, const Point2& uv,
                              const std::size_t ku, const std::size_t kv,
                              const std::size_t maxTotal, SurfaceTable<BasicVec3<T>>& a,
                              SurfaceTable<T>& w)
{
    constexpr std::size_t small = SmallSurfaceDegree;
    constexpr std::size_t large = SurfaceBasisCount;
    if (surface.DegreeU < static_cast<int>(small) && surface.DegreeV < static_cast<int>(small)) {
        if (kv <= LowSurfaceOrders) {
            return HomogeneousSurfaceDerivativesWith<small, small * LowSurfaceOrders>(
                surface, uv, ku, kv, maxTotal, a, w);
        } else {
            return HomogeneousSurfaceDerivativesWith<small, small * small>(surface, uv, ku, kv,
                                                                           maxTotal, a, w);
        }
    } else {
        return HomogeneousSurfaceDerivativesWith<large, large * large>(surface, uv, ku, kv,
                                                                       maxTotal, a, w);
    }
}

/// The rational mixed derivative S^(du,dv) from the homogeneous ones
/// (The NURBS Book, eq. 4.20, algorithm A4.4): A = w S, differentiated by
/// Leibniz's rule in both parameters and solved for the highest term,
/// S^(k,l) = (A^(k,l) - sum_{j>=1} C(l,j) w^(0,j) S^(k,l-j)
///            - sum_{i>=1} C(k,i) w^(i,0) S^(k-i,l)
///            - sum_{i,j>=1} C(k,i) C(l,j) w^(i,j) S^(k-i,l-j)) / w.
template <Scalar T>
void RationalSurfaceDerivatives(SurfaceTable<BasicVec3<T>>& a, SurfaceTable<T>& w,
                                const std::size_t ku, const std::size_t kv, const std::size_t du,
                                const std::size_t dv, const std::size_t maxTotal,
                                std::span<BasicVec3<T>> s)
{
    const auto at = [&](const std::size_t k, const std::size_t l) { return k * (dv + 1) + l; };
    for (std::size_t k = 0; k <= du; ++k) {
        for (std::size_t l = 0; l <= dv && k + l <= maxTotal; ++l) {
            BasicVec3<T> v = (k < ku && l < kv) ? a[k * kv + l] : BasicVec3<T>{};
            for (std::size_t i = 0; i <= k && i < ku; ++i) {
                for (std::size_t j = 0; j <= l && j < kv; ++j) {
                    if (i + j > 0) {
                        const double c = SurfaceBinomial(k, i) * SurfaceBinomial(l, j);
                        v -= (T{c} * w[i * kv + j]) * s[at(k - i, l - j)];
                    }
                }
            }
            s[at(k, l)] = v / w[0];
        }
    }
}

} // namespace Stereon::Detail

export namespace Stereon {

template <Scalar T>
[[nodiscard]] BasicPoint3<T> Eval(const BasicBSplineSurface<T>& surface, const Point2& uv)
{
    Detail::CheckSurfaceEvaluation(surface, uv, "Eval(BSplineSurface)");
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    const std::size_t spanU = FindSpan(surface.KnotsU, surface.DegreeU, PoleCountU(surface), uv.X);
    const std::size_t spanV = FindSpan(surface.KnotsV, surface.DegreeV, countV, uv.Y);
    std::array<double, Detail::SurfaceBasisCount> nu; // filled by BasisFunctions
    std::array<double, Detail::SurfaceBasisCount> nv;
    BasisFunctions(surface.KnotsU, surface.DegreeU, spanU, uv.X, std::span{nu}.first(p + 1));
    BasisFunctions(surface.KnotsV, surface.DegreeV, spanV, uv.Y, std::span{nv}.first(q + 1));
    BasicVec3<T> numerator{};
    T denominator{};
    for (std::size_t r = 0; r <= p; ++r) {
        for (std::size_t s = 0; s <= q; ++s) {
            const std::size_t index = (spanU - p + r) * countV + (spanV - q + s);
            const T weighted = T{nu[r] * nv[s]} * Detail::SurfaceWeightOf(surface, index);
            numerator += weighted * ToVec(surface.Poles[index]);
            denominator = denominator + weighted;
        }
    }
    return IsRational(surface) ? ToPoint(numerator / denominator) : ToPoint(numerator);
}

} // namespace Stereon

export namespace Stereon::Detail {

/// Every partial S^(k,l) with k <= du, l <= dv and k + l <= maxTotal into
/// out[k * (dv + 1) + l]: the homogeneous derivatives, then the quotient
/// rule for a rational surface.
template <Scalar T>
void SurfaceDerivativesUpTo(const BasicBSplineSurface<T>& surface, const Point2& uv,
                            const std::size_t du, const std::size_t dv, const std::size_t maxTotal,
                            std::span<BasicVec3<T>> out)
{
    const std::size_t ku = std::min(du, static_cast<std::size_t>(surface.DegreeU)) + 1;
    const std::size_t kv = std::min(dv, static_cast<std::size_t>(surface.DegreeV)) + 1;
    SurfaceTable<BasicVec3<T>> a(ku * kv);
    SurfaceTable<T> w(ku * kv);
    const BasicVec3<T> origin = HomogeneousSurfaceDerivatives(surface, uv, ku, kv, maxTotal, a, w);
    if (IsRational(surface)) {
        RationalSurfaceDerivatives(a, w, ku, kv, du, dv, maxTotal, out);
    } else {
        for (std::size_t k = 0; k <= du; ++k) {
            for (std::size_t l = 0; l <= dv && k + l <= maxTotal; ++l) {
                out[k * (dv + 1) + l] = (k < ku && l < kv) ? a[k * kv + l] : BasicVec3<T>{};
            }
        }
    }
    out[0] += origin; // the point back in absolute coordinates
}

/// The point and the two first partials at uv, into out[0], out[1] =
/// S^(0,1) and out[2] = S^(1,0), the layout of TotalDerivatives for n = 1:
/// what normals, first fundamental forms and a tessellator ask for most. One
/// pass over the (p + 1)(q + 1) active poles with the basis functions and
/// their first derivatives (A2.3) costs three multiply-adds per pole, where
/// the general path above differences the poles once per order. The poles
/// are taken relative to the one with the largest basis product, as
/// WindowOrigin does and for the same rounding reason, and a rational
/// surface applies the quotient rule once (A4.4 at first order).
template <Scalar T>
void FirstDerivatives(const BasicBSplineSurface<T>& surface, const Point2& uv,
                      std::span<BasicVec3<T>> out)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    const std::size_t spanU = FindSpan(surface.KnotsU, surface.DegreeU, PoleCountU(surface), uv.X);
    const std::size_t spanV = FindSpan(surface.KnotsV, surface.DegreeV, countV, uv.Y);
    std::array<double, 2 * SurfaceBasisCount> nu; // N then N', filled by BasisDerivatives
    std::array<double, 2 * SurfaceBasisCount> nv;
    BasisDerivatives(surface.KnotsU, surface.DegreeU, spanU, uv.X, 1,
                     std::span{nu}.first(2 * (p + 1)));
    BasisDerivatives(surface.KnotsV, surface.DegreeV, spanV, uv.Y, 1,
                     std::span{nv}.first(2 * (q + 1)));
    std::size_t r0 = 0;
    std::size_t c0 = 0;
    for (std::size_t r = 1; r <= p; ++r) {
        r0 = nu[r] > nu[r0] ? r : r0;
    }
    for (std::size_t c = 1; c <= q; ++c) {
        c0 = nv[c] > nv[c0] ? c : c0;
    }
    const BasicVec3<T> origin = ToVec(surface.Poles[(spanU - p + r0) * countV + (spanV - q + c0)]);
    const bool rational = IsRational(surface);
    BasicVec3<T> a{};
    BasicVec3<T> au{};
    BasicVec3<T> av{};
    T w{};
    T wu{};
    T wv{};
    for (std::size_t r = 0; r <= p; ++r) {
        const std::size_t first = (spanU - p + r) * countV + (spanV - q);
        for (std::size_t c = 0; c <= q; ++c) {
            const double n = nu[r] * nv[c];
            const double nU = nu[p + 1 + r] * nv[c];
            const double nV = nu[r] * nv[q + 1 + c];
            BasicVec3<T> d = ToVec(surface.Poles[first + c]) - origin;
            if (rational) {
                const T weight = surface.Weights[first + c];
                d = weight * d;
                w = w + T{n} * weight;
                wu = wu + T{nU} * weight;
                wv = wv + T{nV} * weight;
            }
            a += T{n} * d;
            au += T{nU} * d;
            av += T{nV} * d;
        }
    }
    if (rational) {
        const BasicVec3<T> relative = a / w; // S - origin
        out[0] = relative + origin;
        out[2] = (au - wu * relative) / w;
        out[1] = (av - wv * relative) / w;
    } else {
        out[0] = a + origin;
        out[2] = au;
        out[1] = av;
    }
}

} // namespace Stereon::Detail

export namespace Stereon {

/// The point and every partial derivative up to order du in u and dv in v,
/// into out[k * (dv + 1) + l] = S^(k,l)(uv); out[0] is the position vector.
/// Cheaper than separate calls, because all come from one pair of basis
/// tables (A3.6, and A4.4 for a rational surface). `out` holds
/// (du + 1)(dv + 1) entries.
template <Scalar T>
void Derivatives(const BasicBSplineSurface<T>& surface, const Point2& uv, const int du,
                 const int dv, std::span<BasicVec3<T>> out)
{
    Detail::CheckSurfaceEvaluation(surface, uv, "Derivatives(BSplineSurface)");
    Detail::Require(du >= 0 && dv >= 0 &&
                        out.size() >= static_cast<std::size_t>((du + 1) * (dv + 1)),
                    "Derivatives(BSplineSurface)", "a negative order or too small an output");
    const auto ou = static_cast<std::size_t>(du);
    const auto ov = static_cast<std::size_t>(dv);
    Detail::SurfaceDerivativesUpTo(surface, uv, ou, ov, ou + ov, out);
}

/// The point and every partial derivative of total order k + l <= n, into
/// out[k * (n + 1) + l] = S^(k,l)(uv); the entries with k + l > n are not
/// written. What curvature (n = 2) and Newton steps need, without the
/// higher mixed partials that Derivatives(surface, uv, n, n, out) also
/// computes. `out` holds (n + 1)^2 entries.
template <Scalar T>
void TotalDerivatives(const BasicBSplineSurface<T>& surface, const Point2& uv, const int n,
                      std::span<BasicVec3<T>> out)
{
    Detail::CheckSurfaceEvaluation(surface, uv, "TotalDerivatives(BSplineSurface)");
    Detail::Require(n >= 0 && out.size() >= static_cast<std::size_t>((n + 1) * (n + 1)),
                    "TotalDerivatives(BSplineSurface)", "a negative order or too small an output");
    if (n == 1) {
        Detail::FirstDerivatives(surface, uv, out);
        return;
    }
    const auto order = static_cast<std::size_t>(n);
    Detail::SurfaceDerivativesUpTo(surface, uv, order, order, order, out);
}

/// The mixed partial derivative of order du in u and dv in v, du + dv >= 1.
/// A polynomial surface's derivatives beyond its degree in either parameter
/// are zero.
template <Scalar T>
[[nodiscard]] BasicVec3<T> Deriv(const BasicBSplineSurface<T>& surface, const Point2& uv,
                                 const int du, const int dv)
{
    Detail::CheckSurfaceEvaluation(surface, uv, "Deriv(BSplineSurface)");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(BSplineSurface)",
                    "a derivative order that does not exist");
    if (!IsRational(surface) && (du > surface.DegreeU || dv > surface.DegreeV)) {
        return {};
    }
    if (du + dv == 1) {
        std::array<BasicVec3<T>, 3> first;
        Detail::FirstDerivatives(surface, uv, std::span<BasicVec3<T>>{first});
        return first[du == 1 ? 2 : 1];
    }
    const auto size = static_cast<std::size_t>((du + 1) * (dv + 1));
    if (size <= Detail::SurfaceTableCapacity) {
        std::array<BasicVec3<T>, Detail::SurfaceTableCapacity> all;
        Derivatives(surface, uv, du, dv, std::span{all}.first(size));
        return all[size - 1];
    }
    std::vector<BasicVec3<T>> all(size); // high orders: rare
    Derivatives(surface, uv, du, dv, std::span{all});
    return all[size - 1];
}

/// The box of the poles, widened, as for B-spline curves, by a bound on
/// Eval's rounding; it contains the exact surface (convex hull property)
/// and every evaluated point.
[[nodiscard]] Box3 BoundingBox(const BSplineSurface& surface);

/// The parameter of a point of the surface closest to p, over the whole
/// domain. The surface is split into its Bezier patches; patches whose pole
/// box is farther than the best point so far are skipped. On each other
/// patch the closest points of its four boundary curves are found exactly
/// (Project(BSplineCurve)), and interior stationary points are searched by
/// subdividing the two polynomials whose common roots they are: a
/// sub-patch on which either keeps one sign holds none and is discarded,
/// and the small sub-patches that remain are refined by Newton's method.
/// The closest candidate wins (docs/nurbs.rst).
[[nodiscard]] Point2 Project(const BSplineSurface& surface, const Point3& p);

} // namespace Stereon

export namespace Stereon::Detail {

struct SurfaceProjectorData;

} // namespace Stereon::Detail

export namespace Stereon {

/// Projects many points onto one B-spline surface, as CurveProjector does
/// for curves. Construction prepares the Bezier patches, their coordinate
/// polynomials, the point-independent parts of the two stationarity
/// polynomials, and projectors for every patch edge; a query visits patches
/// and sub-patches closest first and stops at the first that cannot hold a
/// closer point. The answers are those of Project(surface, p), which builds
/// a projector for one point. Copies share the prepared data.
class SurfaceProjector {
public:
    explicit SurfaceProjector(const BSplineSurface& surface);

    /// The parameter of a point of the surface closest to p.
    [[nodiscard]] Point2 Project(const Point3& p) const;

private:
    std::shared_ptr<const Detail::SurfaceProjectorData> m_Data;
};

/// The parameters of a closest point of the surface to each point, in the
/// points' order: what Project gives one by one, computed on the Context's
/// scheduler with its thread count (ADR-0007). Each point's answer goes to
/// its own slot, so the result is the same, bit for bit, for the serial
/// scheduler and for every thread count.
[[nodiscard]] std::vector<Point2> ProjectAll(const SurfaceProjector& projector,
                                             std::span<const Point3> points,
                                             const Context& context);

/// The same surface with the opposite normal, by reversing u:
/// Eval(Reversed(s), (a + b - u, v)) == Eval(s, (u, v)) on the u domain
/// [a, b], up to the rounding of the reflected knots.
template <Scalar T>
[[nodiscard]] BasicBSplineSurface<T> Reversed(const BasicBSplineSurface<T>& surface)
{
    const Interval domain = Domain(surface).X;
    const double sum = domain.GetLower() + domain.GetUpper();
    BasicBSplineSurface<T> reversed = surface;
    std::ranges::reverse(reversed.KnotsU);
    for (double& u : reversed.KnotsU) {
        u = sum - u;
    }
    const std::size_t countU = PoleCountU(surface);
    const std::size_t countV = PoleCountV(surface);
    for (std::size_t i = 0; i < countU; ++i) {
        for (std::size_t j = 0; j < countV; ++j) {
            const std::size_t from = (countU - 1 - i) * countV + j;
            reversed.Poles[i * countV + j] = surface.Poles[from];
            if (IsRational(surface)) {
                reversed.Weights[i * countV + j] = surface.Weights[from];
            }
        }
    }
    return reversed;
}

/// The image of the surface under any affine map, by mapping the poles
/// (affine invariance, as for curves).
template <Scalar T>
[[nodiscard]] BasicBSplineSurface<T> Transformed(const BasicBSplineSurface<T>& surface,
                                                 const BasicTransform<T>& x)
{
    BasicBSplineSurface<T> moved = surface;
    for (BasicPoint3<T>& pole : moved.Poles) {
        pole = Transformed(pole, x);
    }
    return moved;
}

/// The surface with the knot value inserted `times` times in one
/// direction (A5.3): the same surface with more poles. InvalidInput unless
/// the value lies in that direction's domain and its multiplicity stays at
/// most the degree.
[[nodiscard]] Result<BSplineSurface> InsertKnot(const BSplineSurface& surface,
                                                ParameterDirection direction, double value,
                                                int times = 1);

/// The surface's Bezier patches, row by row in u then v: one clamped
/// single-span surface per pair of non-empty knot spans of the domain.
[[nodiscard]] std::vector<BSplineSurface> ToBezierPatches(const BSplineSurface& surface);

/// The iso-parametric curve at `value` of the given parameter: for
/// direction U, the curve v -> S(value, v), and for V, u -> S(u, value),
/// with the other direction's degree, knots and parametrisation.
/// InvalidInput unless value lies in that direction's domain.
[[nodiscard]] Result<BSplineCurve> IsoCurve(const BSplineSurface& surface,
                                            ParameterDirection direction, double value);

static_assert(SurfaceGeometry<BSplineSurface>);

} // namespace Stereon
