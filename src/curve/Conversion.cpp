// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md): every module
// implementation unit of a geometry library includes it, because one that
// calls core's inline functions with contracts and parses none itself
// leaves __tu_has_violation undefined, which the macOS linker rejects.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <string_view>
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
import :Visit;

namespace Stereon {

namespace {

/// The largest angle one rational quadratic arc spans: a quarter turn keeps
/// its middle weight cos(pi/4) and its middle pole close to the circle
/// (The NURBS Book, A7.1).
constexpr double MaxArcAngle = std::numbers::pi / 2;

/// The largest parameter length of one hyperbolic arc: the middle weight
/// cosh(1/2) ~ 1.13 stays near 1 and the middle pole near the curve.
constexpr double MaxHyperbolicArc = 1.0;

/// InvalidInput unless range is finite, non-empty, and inside the domain,
/// or for a periodic curve (period > 0) at most one period long.
[[nodiscard]] Result<void> CheckRange(const Interval range, const Interval domain,
                                      const double period, const std::string_view where)
{
    const double low = range.GetLower();
    const double high = range.GetUpper();
    if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        std::string message{where};
        message += ": the range must be finite and non-empty";
        return Fail(ErrorCode::InvalidInput, std::move(message));
    }
    const bool fits =
        period > 0.0 ? high - low <= period : Contains(domain, low) && Contains(domain, high);
    if (!fits) {
        std::string message{where};
        message += period > 0.0 ? ": the range must span at most one period"
                                : ": the range must lie in the domain";
        return Fail(ErrorCode::InvalidInput, std::move(message));
    }
    return {};
}

/// The number of segments of at most `limit` that cover [low, high].
[[nodiscard]] std::size_t SegmentCount(const double low, const double high, const double limit)
{
    return std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((high - low) / limit)));
}

/// Knots for n quadratic segments with the given interior break values:
/// each end three times, each break twice.
[[nodiscard]] std::vector<double> QuadraticKnots(const std::vector<double>& breaks)
{
    std::vector<double> knots(3, breaks.front());
    for (std::size_t k = 1; k + 1 < breaks.size(); ++k) {
        knots.insert(knots.end(), 2, breaks[k]);
    }
    knots.insert(knots.end(), 3, breaks.back());
    return knots;
}

[[nodiscard]] std::vector<double> Breaks(const double low, const double high, const std::size_t n)
{
    std::vector<double> breaks;
    const double step = (high - low) / static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k) {
        breaks.push_back(low + static_cast<double>(k) * step);
    }
    breaks.push_back(high); // exactly
    return breaks;
}

/// The unit arc's poles mapped by origin + x * X + y * Y, scaled per axis.
[[nodiscard]] BSplineCurve MappedArc(const Detail::UnitArc& arc, const Point3& origin,
                                     const Vec3& x, const Vec3& y)
{
    BSplineCurve curve{2, arc.Knots, {}, arc.Weights};
    for (const Point2& p : arc.Poles) {
        curve.Poles.push_back(origin + p.X * x + p.Y * y);
    }
    return curve;
}

/// One hyperbolic arc over [t0, t1] as a rational quadratic Bezier curve:
/// with half-width h and midpoint m, the middle pole is H's offset at m
/// divided by cosh h, with weight cosh h. The arc's point at its parameter
/// 1/2 is then H(m), because cosh t0 + cosh t1 = 2 cosh m cosh h and
/// sinh t0 + sinh t1 = 2 sinh m cosh h: the hyperbolic analogue of the
/// circular arc.
void AppendHyperbolicArc(const Hyperbola& h, const double t0, const double t1, BSplineCurve& curve)
{
    const double half = 0.5 * (t1 - t0);
    const double mid = t0 + half;
    const double weight = Cosh(half);
    if (curve.Poles.empty()) {
        curve.Poles.push_back(Eval(h, t0));
        curve.Weights.push_back(1.0);
    }
    curve.Poles.push_back(h.Center + (h.MajorRadius * Cosh(mid) / weight) * h.XAxis +
                          (h.MinorRadius * Sinh(mid) / weight) * h.YAxis);
    curve.Weights.push_back(weight);
    curve.Poles.push_back(Eval(h, t1));
    curve.Weights.push_back(1.0);
}

} // namespace

Result<BSplineCurve> ToBSpline(const Line& line, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(line), 0.0, "ToBSpline(Line)"));
    const double a = range.GetLower();
    const double b = range.GetUpper();
    return BSplineCurve{1, {a, a, b, b}, {Eval(line, a), Eval(line, b)}, {}};
}

Result<BSplineCurve> ToBSpline(const Circle& circle, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(circle), Period(circle), "ToBSpline(Circle)"));
    const Detail::UnitArc arc = Detail::MakeUnitArc(range.GetLower(), range.GetUpper());
    return MappedArc(arc, circle.Center, circle.Radius * circle.XAxis,
                     circle.Radius * circle.YAxis);
}

Result<BSplineCurve> ToBSpline(const Ellipse& ellipse, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(ellipse), Period(ellipse), "ToBSpline(Ellipse)"));
    // An ellipse is the affine image of the unit circle, and NURBS are
    // affinely invariant: the circle's weights carry over unchanged.
    const Detail::UnitArc arc = Detail::MakeUnitArc(range.GetLower(), range.GetUpper());
    return MappedArc(arc, ellipse.Center, ellipse.MajorRadius * ellipse.XAxis,
                     ellipse.MinorRadius * ellipse.YAxis);
}

Result<BSplineCurve> ToBSpline(const Hyperbola& hyperbola, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(hyperbola), 0.0, "ToBSpline(Hyperbola)"));
    const std::vector<double> breaks =
        Breaks(range.GetLower(), range.GetUpper(),
               SegmentCount(range.GetLower(), range.GetUpper(), MaxHyperbolicArc));
    BSplineCurve curve{2, QuadraticKnots(breaks), {}, {}};
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        AppendHyperbolicArc(hyperbola, breaks[k], breaks[k + 1], curve);
    }
    return curve;
}

Result<BSplineCurve> ToBSpline(const Parabola& parabola, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(parabola), 0.0, "ToBSpline(Parabola)"));
    // A quadratic in t: the Bezier form on [a, b] has the tangent
    // intersection as its middle pole, and the same parametrisation.
    const double a = range.GetLower();
    const double b = range.GetUpper();
    const Point3 start = Eval(parabola, a);
    return BSplineCurve{2,
                        {a, a, a, b, b, b},
                        {start, start + (0.5 * (b - a)) * Deriv(parabola, a, 1), Eval(parabola, b)},
                        {}};
}

Result<BSplineCurve> ToBSpline(const BSplineCurve& curve, const Interval range)
{
    STEREON_TRY_VOID(CheckRange(range, Domain(curve), 0.0, "ToBSpline(BSplineCurve)"));
    BSplineCurve piece = curve;
    if (range.GetLower() > Domain(piece).GetLower()) {
        STEREON_TRY(auto halves, Split(piece, range.GetLower()));
        piece = std::move(halves.second);
    }
    if (range.GetUpper() < Domain(piece).GetUpper()) {
        STEREON_TRY(auto halves, Split(piece, range.GetUpper()));
        piece = std::move(halves.first);
    }
    return piece;
}

Result<BSplineCurve> ToBSpline(const Curve& curve, const Interval range)
{
    return VisitBuiltin(curve, [&]<typename G>(const G& geometry) -> Result<BSplineCurve> {
        if constexpr (std::is_same_v<G, Curve>) {
            return Fail(ErrorCode::NotSupported,
                        "ToBSpline: the curve's type has no exact NURBS form the kernel knows");
        } else {
            return ToBSpline(geometry, range);
        }
    });
}

} // namespace Stereon

namespace Stereon::Detail {

UnitArc MakeUnitArc(const double start, const double end)
{
    const std::vector<double> breaks = Breaks(start, end, SegmentCount(start, end, MaxArcAngle));
    UnitArc arc{QuadraticKnots(breaks), {}, {}};
    arc.Poles.push_back({Cos(start), Sin(start)});
    arc.Weights.push_back(1.0);
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        const double half = 0.5 * (breaks[k + 1] - breaks[k]);
        const double mid = breaks[k] + half;
        const double weight = Cos(half);
        arc.Poles.push_back({Cos(mid) / weight, Sin(mid) / weight});
        arc.Weights.push_back(weight);
        arc.Poles.push_back({Cos(breaks[k + 1]), Sin(breaks[k + 1])});
        arc.Weights.push_back(1.0);
    }
    // A full turn closes exactly: cos and sin of start + 2 pi round
    // differently from those of start.
    if (end - start == CirclePeriod) {
        arc.Poles.back() = arc.Poles.front();
    }
    return arc;
}

} // namespace Stereon::Detail
