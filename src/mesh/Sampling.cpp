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
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Mesh;

import :Sampling;

namespace Stereon {

// Instantiated at namespace scope before any function uses them
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));
static_assert(sizeof(Vec2) == 2 * sizeof(double));

namespace {

/// Segments are split at most this deep below a knot span's first four:
/// 2^-40 of the span, far below anything a display resolves.
constexpr int MaxDepth = 40;

/// Each knot span starts as this many segments, so a span whose midpoint
/// happens to lie on its chord (an S-shaped cubic) is still refined.
constexpr int InitialSegments = 4;

/// The interior points of a segment checked against its chord: k/8.
constexpr int ChordChecks = 8;

[[nodiscard]] Result<void> CheckTolerance(const DisplayTolerance& tolerance,
                                          const std::string_view function)
{
    const bool positive = std::isfinite(tolerance.Chord) && std::isfinite(tolerance.Angle) &&
                          tolerance.Chord > 0.0 && // STEREON-LINT-ALLOW(epsilon): input validation
                          tolerance.Angle > 0.0;   // STEREON-LINT-ALLOW(epsilon): input validation
    if (!positive || tolerance.MaxPoints < 2) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(function) + ": the display tolerance must be finite and positive "
                                            "and allow at least two points");
    }
    return {};
}

/// `range` clipped to `domain` when the domain is not periodic, or the
/// domain itself; InvalidInput unless the result is finite and non-empty.
[[nodiscard]] Result<Interval> SampledRange(const Interval domain, const bool periodic,
                                            const std::optional<Interval> range,
                                            const std::string_view function)
{
    double lower = domain.GetLower();
    double upper = domain.GetUpper();
    if (range) {
        const bool usable = std::isfinite(range->GetLower()) && std::isfinite(range->GetUpper());
        lower = !usable || periodic ? range->GetLower() : std::max(lower, range->GetLower());
        upper = !usable || periodic ? range->GetUpper() : std::min(upper, range->GetUpper());
    }
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(function) + ": the range to sample must be finite and non-empty; "
                                            "unbounded geometry needs one from the caller");
    }
    return Interval{lower, upper};
}

template <std::size_t D>
[[nodiscard]] double DistanceToSegment(const BasicPointN<double, D>& p,
                                       const BasicPointN<double, D>& a,
                                       const BasicPointN<double, D>& b)
{
    const BasicVecN<double, D> v = b - a;
    const BasicVecN<double, D> w = p - a;
    const double length2 = Dot(v, v);
    if (!(length2 > 0.0)) { // STEREON-LINT-ALLOW(epsilon): a degenerate segment, not a tolerance
        return Length(w);
    }
    const double s = std::clamp(Dot(w, v) / length2, 0.0, 1.0);
    return Length(w - s * v);
}

/// The angle between two vectors, zero when either vanishes.
template <std::size_t D>
[[nodiscard]] double Turn(const BasicVecN<double, D>& a, const BasicVecN<double, D>& b)
{
    const double lengths = Length(a) * Length(b);
    if (!(lengths > 0.0)) { // STEREON-LINT-ALLOW(epsilon): a vanishing half, not a tolerance
        return 0.0;
    }
    return std::acos(std::clamp(Dot(a, b) / lengths, -1.0, 1.0));
}

/// One segment of a polyline being refined.
struct Segment {
    double Start = 0.0;
    double End = 0.0;
    int Depth = 0;
};

/// True when the curve between the segment's ends is drawn well enough by
/// the chord from `p0` to `p1`.
template <std::size_t D, typename Eval>
[[nodiscard]] bool IsFlat(const Eval& eval, const Segment& s, const BasicPointN<double, D>& p0,
                          const BasicPointN<double, D>& p1, const DisplayTolerance& tolerance)
{
    const double width = s.End - s.Start;
    for (int k = 1; k < ChordChecks; ++k) {
        const BasicPointN<double, D> q = eval(s.Start + width * k / ChordChecks);
        if (!(DistanceToSegment<D>(q, p0, p1) <= tolerance.Chord)) {
            return false;
        }
    }
    const BasicPointN<double, D> middle = eval(s.Start + width / 2);
    return Turn<D>(middle - p0, p1 - middle) <= tolerance.Angle;
}

/// Appends the refined polyline of one segment, its start already in `out`.
template <std::size_t D, typename Eval>
[[nodiscard]] bool Refine(const Eval& eval, const Segment& first, const DisplayTolerance& tolerance,
                          BasicPolyline<D>& out)
{
    std::vector<Segment> stack{first};
    while (!stack.empty()) {
        const Segment s = stack.back();
        stack.pop_back();
        const BasicPointN<double, D> p0 = out.Points.back();
        const BasicPointN<double, D> p1 = eval(s.End);
        const double middle = s.Start + (s.End - s.Start) / 2;
        const bool splittable = s.Depth < MaxDepth && s.Start < middle && middle < s.End;
        if (splittable && !IsFlat<D>(eval, s, p0, p1, tolerance)) {
            stack.push_back({middle, s.End, s.Depth + 1});
            stack.push_back({s.Start, middle, s.Depth + 1});
            continue;
        }
        if (out.Points.size() >= tolerance.MaxPoints) {
            return false;
        }
        out.Parameters.push_back(s.End);
        out.Points.push_back(p1);
    }
    return true;
}

/// The polyline of `eval` through every breakpoint, ascending and unique,
/// the range's ends first and last.
template <std::size_t D, typename Eval>
[[nodiscard]] Result<BasicPolyline<D>>
SampleBetween(const Eval& eval, const std::vector<double>& breaks,
              const DisplayTolerance& tolerance, const std::string_view function)
{
    BasicPolyline<D> out;
    out.Parameters.push_back(breaks.front());
    out.Points.push_back(eval(breaks.front()));
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double width = breaks[i + 1] - breaks[i];
        for (int k = 0; k < InitialSegments; ++k) {
            const double start = breaks[i] + width * k / InitialSegments;
            const double end = k + 1 == InitialSegments
                                   ? breaks[i + 1]
                                   : breaks[i] + width * (k + 1) / InitialSegments;
            if (!Refine<D>(eval, Segment{start, end, 0}, tolerance, out)) {
                return Fail(ErrorCode::NotSupported,
                            std::string(function) + ": the display tolerance needs more than " +
                                std::to_string(tolerance.MaxPoints) + " points");
            }
        }
    }
    return out;
}

/// The range's ends and the values strictly inside it, sorted and unique.
[[nodiscard]] std::vector<double> WithEnds(std::vector<double> values, const Interval range)
{
    std::erase_if(values,
                  [&](const double x) { return !(range.GetLower() < x && x < range.GetUpper()); });
    values.push_back(range.GetLower());
    values.push_back(range.GetUpper());
    std::ranges::sort(values);
    const auto [first, last] = std::ranges::unique(values);
    values.erase(first, last);
    return values;
}

/// The knots of a B-spline basis, also shifted by one period each way for
/// a periodic basis, so that a trim across the seam keeps them.
template <typename Knots>
void AddKnots(const Knots& knots, const double period, std::vector<double>& out)
{
    for (const double k : knots) {
        out.push_back(k);
        if (period > 0.0) { // STEREON-LINT-ALLOW(epsilon): Period() is 0 when not periodic
            out.push_back(k - period);
            out.push_back(k + period);
        }
    }
}

template <std::size_t D>
void CollectKnots(const BasicErasedCurve<D>& curve, std::vector<double>& out, const int depth)
{
    if (const auto* bspline = curve.template Target<BasicBSplineCurve<double, D>>()) {
        AddKnots(bspline->Knots, Period(curve), out);
        return;
    }
    const auto* trimmed = curve.template Target<BasicTrimmedCurve<D>>();
    if (trimmed == nullptr || depth >= MaxDepth) {
        return;
    }
    std::vector<double> basis;
    CollectKnots(trimmed->GetBasis(), basis, depth + 1);
    for (const double s : basis) {
        const Interval range = trimmed->GetRange();
        if (range.GetLower() <= s && s <= range.GetUpper()) {
            out.push_back(trimmed->IsReversed() ? trimmed->Reflect(s) : s);
        }
    }
}

template <std::size_t D>
[[nodiscard]] Result<BasicPolyline<D>> SampleAnyCurve(const BasicErasedCurve<D>& curve,
                                                      const DisplayTolerance& tolerance,
                                                      const std::optional<Interval> range)
{
    constexpr std::string_view function = "SampleCurve";
    if (curve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "SampleCurve: the curve is empty");
    }
    STEREON_TRY_VOID(CheckTolerance(tolerance, function));
    const bool periodic = Period(curve) > 0.0; // STEREON-LINT-ALLOW(epsilon): 0 means not periodic
    STEREON_TRY(const Interval sampled, SampledRange(Domain(curve), periodic, range, function));
    const auto eval = [&](const double t) { return Eval(curve, t); };
    std::vector<double> knots;
    CollectKnots(curve, knots, 0);
    return SampleBetween<D>(eval, WithEnds(std::move(knots), sampled), tolerance, function);
}

/// The distinct interior knots of a surface in u and v: a B-spline
/// surface's, or a trimmed surface's B-spline basis's.
void SurfaceKnots(const Surface& surface, std::vector<double>& u, std::vector<double>& v,
                  const int depth)
{
    const Vec2 periods = Periods(surface);
    if (const auto* bspline = surface.Target<BSplineSurface>()) {
        AddKnots(bspline->KnotsU, periods.X, u);
        AddKnots(bspline->KnotsV, periods.Y, v);
        return;
    }
    const auto* trimmed = surface.Target<TrimmedSurface>();
    if (trimmed == nullptr || depth >= MaxDepth) {
        return;
    }
    std::vector<double> basisU;
    SurfaceKnots(trimmed->GetBasis(), basisU, v, depth + 1);
    const Interval range = trimmed->GetRange().X;
    for (const double s : basisU) {
        if (range.GetLower() <= s && s <= range.GetUpper()) {
            u.push_back(trimmed->IsReversed() ? trimmed->Reflect(s) : s);
        }
    }
}

/// `count` evenly spaced values over the range, ends included, and the knots.
[[nodiscard]] std::vector<double> IsoValues(const Interval range, const std::size_t count,
                                            std::vector<double> knots)
{
    const double width = range.GetUpper() - range.GetLower();
    for (std::size_t i = 1; i + 1 < count; ++i) {
        knots.push_back(range.GetLower() +
                        width * static_cast<double>(i) / static_cast<double>(count - 1));
    }
    return WithEnds(std::move(knots), range);
}

/// Built by member, not by aggregate initialisation of a type from another
/// module unit (docs/toolchain/gcc16-feature-probe.md).
[[nodiscard]] IsoLine MakeIsoLine(const IsoDirection direction, const double fixed, Polyline line)
{
    IsoLine iso;
    iso.Direction = direction;
    iso.Fixed = fixed;
    iso.Line = std::move(line);
    return iso;
}

} // namespace

std::vector<double> Breakpoints(const Curve& curve, const Interval range)
{
    std::vector<double> knots;
    CollectKnots(curve, knots, 0);
    return WithEnds(std::move(knots), range);
}

std::vector<double> Breakpoints(const Curve2d& curve, const Interval range)
{
    std::vector<double> knots;
    CollectKnots(curve, knots, 0);
    return WithEnds(std::move(knots), range);
}

Result<Polyline> SampleCurve(const Curve& curve, const DisplayTolerance& tolerance,
                             const std::optional<Interval> range)
{
    return SampleAnyCurve(curve, tolerance, range);
}

Result<Polyline2d> SampleCurve(const Curve2d& curve, const DisplayTolerance& tolerance,
                               const std::optional<Interval> range)
{
    return SampleAnyCurve(curve, tolerance, range);
}

Result<std::vector<IsoLine>> SampleSurface(const Surface& surface,
                                           const DisplayTolerance& tolerance,
                                           const std::optional<Box2> range,
                                           const std::size_t linesPerDirection)
{
    constexpr std::string_view function = "SampleSurface";
    if (surface.IsEmpty() || linesPerDirection < 2) {
        return Fail(ErrorCode::InvalidInput,
                    "SampleSurface: needs a surface and at least two lines per direction");
    }
    STEREON_TRY_VOID(CheckTolerance(tolerance, function));
    const Box2 domain = Domain(surface);
    const Vec2 periods = Periods(surface);
    const bool periodicU = periods.X > 0.0; // STEREON-LINT-ALLOW(epsilon): 0 means not periodic
    const bool periodicV = periods.Y > 0.0; // STEREON-LINT-ALLOW(epsilon): 0 means not periodic
    STEREON_TRY(const Interval u,
                SampledRange(domain.X, periodicU,
                             range ? std::optional<Interval>{range->X} : std::nullopt, function));
    STEREON_TRY(const Interval v,
                SampledRange(domain.Y, periodicV,
                             range ? std::optional<Interval>{range->Y} : std::nullopt, function));
    std::vector<double> knotsU;
    std::vector<double> knotsV;
    SurfaceKnots(surface, knotsU, knotsV, 0);
    const std::vector<double> breaksU = WithEnds(knotsU, u);
    const std::vector<double> breaksV = WithEnds(knotsV, v);

    std::vector<IsoLine> lines;
    for (const double fixed : IsoValues(v, linesPerDirection, knotsV)) {
        const auto eval = [&](const double s) { return Eval(surface, Point2{s, fixed}); };
        STEREON_TRY(Polyline line, SampleBetween<3>(eval, breaksU, tolerance, function));
        lines.push_back(MakeIsoLine(IsoDirection::AlongU, fixed, std::move(line)));
    }
    for (const double fixed : IsoValues(u, linesPerDirection, knotsU)) {
        const auto eval = [&](const double s) { return Eval(surface, Point2{fixed, s}); };
        STEREON_TRY(Polyline line, SampleBetween<3>(eval, breaksV, tolerance, function));
        lines.push_back(MakeIsoLine(IsoDirection::AlongV, fixed, std::move(line)));
    }
    return lines;
}

} // namespace Stereon
