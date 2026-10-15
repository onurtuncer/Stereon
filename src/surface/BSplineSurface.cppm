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

/// The homogeneous derivatives A^(k,l) = sum N_i^(k) M_j^(l) w P and
/// w^(k,l) = sum N_i^(k) M_j^(l) w for k < ku, l < kv (algorithm A3.6 on the
/// homogeneous form), into tables indexed k * kv + l.
template <std::size_t Capacity, Scalar T>
void HomogeneousSurfaceDerivativesWith(const BasicBSplineSurface<T>& surface, const Point2& uv,
                                       const std::size_t ku, const std::size_t kv,
                                       SurfaceTable<BasicVec3<T>>& a, SurfaceTable<T>& w)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    const std::size_t spanU = FindSpan(surface.KnotsU, surface.DegreeU, PoleCountU(surface), uv.X);
    const std::size_t spanV = FindSpan(surface.KnotsV, surface.DegreeV, countV, uv.Y);
    // Filled by BasisDerivatives for the orders asked.
    std::array<double, Capacity * Capacity> nu;
    std::array<double, Capacity * Capacity> nv;
    BasisDerivatives(surface.KnotsU, surface.DegreeU, spanU, uv.X, static_cast<int>(ku) - 1,
                     std::span{nu}.first(ku * (p + 1)));
    BasisDerivatives(surface.KnotsV, surface.DegreeV, spanV, uv.Y, static_cast<int>(kv) - 1,
                     std::span{nv}.first(kv * (q + 1)));
    for (std::size_t k = 0; k < ku; ++k) {
        for (std::size_t l = 0; l < kv; ++l) {
            BasicVec3<T> sum{};
            T weight{};
            for (std::size_t r = 0; r <= p; ++r) {
                for (std::size_t s = 0; s <= q; ++s) {
                    const std::size_t index = (spanU - p + r) * countV + (spanV - q + s);
                    const T coefficient = T{nu[k * (p + 1) + r] * nv[l * (q + 1) + s]} *
                                          SurfaceWeightOf(surface, index);
                    sum += coefficient * ToVec(surface.Poles[index]);
                    weight = weight + coefficient;
                }
            }
            a[k * kv + l] = sum;
            w[k * kv + l] = weight;
        }
    }
}

template <Scalar T>
void HomogeneousSurfaceDerivatives(const BasicBSplineSurface<T>& surface, const Point2& uv,
                                   const std::size_t ku, const std::size_t kv,
                                   SurfaceTable<BasicVec3<T>>& a, SurfaceTable<T>& w)
{
    if (surface.DegreeU < static_cast<int>(SmallSurfaceDegree) &&
        surface.DegreeV < static_cast<int>(SmallSurfaceDegree)) {
        HomogeneousSurfaceDerivativesWith<SmallSurfaceDegree>(surface, uv, ku, kv, a, w);
    } else {
        HomogeneousSurfaceDerivativesWith<SurfaceBasisCount>(surface, uv, ku, kv, a, w);
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
                                const std::size_t dv, std::span<BasicVec3<T>> s)
{
    const auto at = [&](const std::size_t k, const std::size_t l) { return k * (dv + 1) + l; };
    for (std::size_t k = 0; k <= du; ++k) {
        for (std::size_t l = 0; l <= dv; ++l) {
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
    const std::size_t ku = std::min(ou, static_cast<std::size_t>(surface.DegreeU)) + 1;
    const std::size_t kv = std::min(ov, static_cast<std::size_t>(surface.DegreeV)) + 1;
    Detail::SurfaceTable<BasicVec3<T>> a(ku * kv);
    Detail::SurfaceTable<T> w(ku * kv);
    Detail::HomogeneousSurfaceDerivatives(surface, uv, ku, kv, a, w);
    if (IsRational(surface)) {
        Detail::RationalSurfaceDerivatives(a, w, ku, kv, ou, ov, out);
        return;
    }
    for (std::size_t k = 0; k <= ou; ++k) {
        for (std::size_t l = 0; l <= ov; ++l) {
            out[k * (ov + 1) + l] = (k < ku && l < kv) ? a[k * kv + l] : BasicVec3<T>{};
        }
    }
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
