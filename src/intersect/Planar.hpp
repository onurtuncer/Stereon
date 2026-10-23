// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the plane curve handlers share: a formula written once over a
// scalar field and read as its enclosure, as its rounded value or exactly,
// so that a configuration decision is the sign of a polynomial in the
// input doubles (RS-0001 G5) settled on the heap only when the enclosure
// straddles zero; the placement of a periodic parameter in a requested
// domain; the ordering of a result; and its two sides exchanged.
//
// Included in the purview of the handlers' implementation units, after
// `import :Dispatch;` and `import Stereon.Robust;`, with <algorithm>,
// <cmath>, <cstddef>, <optional>, <utility> and <vector> in their global
// module fragment, and with Analytic.hpp. Everything is in an anonymous
// namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

// --- The scalar fields --------------------------------------------------

/// Outward-rounded intervals: the enclosure of a formula. Each unit uses the
/// operations its formulas need, hence maybe_unused.
struct Rounded {
    using Scalar = Interval;
    [[maybe_unused]] static Interval Of(const double x) { return Interval{x}; }
    [[maybe_unused]] static Interval Diff(const double a, const double b)
    {
        return Interval{a} - Interval{b};
    }
    [[maybe_unused]] static Interval Add(const Interval& a, const Interval& b) { return a + b; }
    [[maybe_unused]] static Interval Sub(const Interval& a, const Interval& b) { return a - b; }
    [[maybe_unused]] static Interval Mul(const Interval& a, const Interval& b) { return a * b; }
};

/// Doubles: the rounded value of a formula, for the geometry it then
/// constructs (never for a decision).
struct Plain {
    using Scalar = double;
    [[maybe_unused]] static double Of(const double x) { return x; }
    [[maybe_unused]] static double Diff(const double a, const double b) { return a - b; }
    [[maybe_unused]] static double Add(const double a, const double b) { return a + b; }
    [[maybe_unused]] static double Sub(const double a, const double b) { return a - b; }
    [[maybe_unused]] static double Mul(const double a, const double b) { return a * b; }
};

/// Expansions on the heap: the exact value of a formula.
struct Exact {
    using Scalar = DynamicExpansion;
    [[maybe_unused]] static DynamicExpansion Of(const double x) { return DynamicExpansion{x}; }
    [[maybe_unused]] static DynamicExpansion Diff(const double a, const double b)
    {
        return DynamicExpansion{TwoDiff(a, b)};
    }
    [[maybe_unused]] static DynamicExpansion Add(const DynamicExpansion& a,
                                                 const DynamicExpansion& b)
    {
        return Sum(a, b);
    }
    [[maybe_unused]] static DynamicExpansion Sub(const DynamicExpansion& a,
                                                 const DynamicExpansion& b)
    {
        return Difference(a, b);
    }
    [[maybe_unused]] static DynamicExpansion Mul(const DynamicExpansion& a,
                                                 const DynamicExpansion& b)
    {
        return Product(a, b);
    }
};

/// The sign of a formula, a callable taking a field and returning its
/// Scalar: from the enclosure when that decides it, exactly otherwise.
template <typename Formula>
[[nodiscard]] Sign SignOfFormula(Formula&& formula)
{
    return FilteredSign(formula(Rounded{}), [&] { return formula(Exact{}).GetSign(); });
}

/// A vector over a field.
template <typename F>
struct Pair2 {
    typename F::Scalar X;
    typename F::Scalar Y;
};

template <typename F>
[[nodiscard]] Pair2<F> Of2(F, const Vec2& v)
{
    return {F::Of(v.X), F::Of(v.Y)};
}

/// p - q exactly.
template <typename F>
[[nodiscard]] Pair2<F> Diff2(F, const Point2& p, const Point2& q)
{
    return {F::Diff(p.X, q.X), F::Diff(p.Y, q.Y)};
}

template <typename F>
[[nodiscard]] typename F::Scalar Dot2(F, const Pair2<F>& a, const Pair2<F>& b)
{
    return F::Add(F::Mul(a.X, b.X), F::Mul(a.Y, b.Y));
}

/// The determinant |a b|: a x b.
template <typename F>
[[nodiscard]] typename F::Scalar Cross2(F, const Pair2<F>& a, const Pair2<F>& b)
{
    return F::Sub(F::Mul(a.X, b.Y), F::Mul(a.Y, b.X));
}

// --- Parameters in their domains ----------------------------------------

/// The parameter of a periodic curve placed in the requested domain: t plus
/// the multiple of the period that lands inside, the smallest when several
/// do, or nothing when none does.
[[nodiscard]] inline std::optional<double> InPeriodicDomain(const double t, const Interval& domain)
{
    if (!domain.IsBounded()) {
        if (!std::isfinite(domain.GetLower())) {
            return std::isfinite(domain.GetUpper()) && t > domain.GetUpper()
                       ? t - CirclePeriod * std::ceil((t - domain.GetUpper()) / CirclePeriod)
                       : t;
        }
        return t < domain.GetLower()
                   ? t + CirclePeriod * std::ceil((domain.GetLower() - t) / CirclePeriod)
                   : t;
    }
    const double first = t + CirclePeriod * std::ceil((domain.GetLower() - t) / CirclePeriod);
    for (double candidate = first; candidate <= domain.GetUpper(); candidate += CirclePeriod) {
        if (candidate >= domain.GetLower()) {
            return candidate;
        }
    }
    return std::nullopt;
}

/// Keeps the point when both parameters lie in their domains, the periodic
/// ones placed there.
inline void KeepPoint(PlanarIntersection& result, PlanarPoint point, const Interval& domainA,
                      const bool periodicA, const Interval& domainB, const bool periodicB)
{
    if (periodicA) {
        const std::optional<double> placed = InPeriodicDomain(point.ParameterA, domainA);
        if (!placed) {
            return;
        }
        point.ParameterA = *placed;
    } else if (!Contains(domainA, point.ParameterA)) {
        return;
    }
    if (periodicB) {
        const std::optional<double> placed = InPeriodicDomain(point.ParameterB, domainB);
        if (!placed) {
            return;
        }
        point.ParameterB = *placed;
    } else if (!Contains(domainB, point.ParameterB)) {
        return;
    }
    result.Points.push_back(point);
}

/// Points in order of the parameter on A, then on B (ADR-0029 G1).
inline void Order(PlanarIntersection& result)
{
    std::ranges::sort(result.Points, [](const PlanarPoint& p, const PlanarPoint& q) {
        return p.ParameterA != q.ParameterA ? p.ParameterA < q.ParameterA
                                            : p.ParameterB < q.ParameterB;
    });
}

/// The same result seen from the other side: A and B exchanged.
[[nodiscard]] inline PlanarIntersection SwappedPlanar(PlanarIntersection result)
{
    for (PlanarPoint& point : result.Points) {
        std::swap(point.ParameterA, point.ParameterB);
    }
    for (PlanarOverlap& overlap : result.Overlaps) {
        std::swap(overlap.OnA, overlap.OnB);
    }
    Order(result);
    return result;
}

/// The angle of p in the frame (center, x, y) of a circle or an ellipse
/// with semi-axes a and b, in [0, 2 pi).
[[nodiscard]] inline double AngleIn(const Point2& center, const Vec2& x, const Vec2& y,
                                    const double a, const double b, const Point2& p)
{
    const Vec2 offset = p - center;
    return WrappedAngle(Dot(offset, y) / b, Dot(offset, x) / a);
}

} // namespace
} // namespace Stereon::Detail
