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
#include <limits>
#include <optional>
#include <span>
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

/// The largest magnitude among the point's coordinates.
template <std::size_t D>
[[nodiscard]] double MaxAbs(const BasicPointN<double, D>& p)
{
    double most = 0.0;
    for (std::size_t i = 0; i < D; ++i) {
        most = std::max(most, std::abs(p[i]));
    }
    return most;
}

/// A span this many units in the last place of the coordinates is a point
/// for display: 64 roundings.
constexpr double RoundingSpan =
    64.0 * std::numeric_limits<double>::epsilon(); // STEREON-LINT-ALLOW(epsilon): display rounding,
                                                   // not a tolerance

/// One segment of a polyline being refined.
struct Segment {
    double Start = 0.0;
    double End = 0.0;
    int Depth = 0;
};

template <std::size_t D>
[[nodiscard]] bool FinitePoint(const BasicPointN<double, D>& point)
{
    for (std::size_t i = 0; i < D; ++i) {
        if (!std::isfinite(point[i])) {
            return false;
        }
    }
    return true;
}

/// True when the curve between the segment's ends is drawn well enough by
/// the chord from `p0` to `p1`.
template <std::size_t D, typename Eval>
[[nodiscard]] bool IsFlat(const Eval& eval, const Segment& s, const BasicPointN<double, D>& p0,
                          const BasicPointN<double, D>& p1, const DisplayTolerance& tolerance)
{
    if (!FinitePoint<D>(p0) || !FinitePoint<D>(p1)) {
        return false;
    }
    const double width = s.End - s.Start;
    for (int k = 1; k < ChordChecks; ++k) {
        const BasicPointN<double, D> q = eval(s.Start + width * k / ChordChecks);
        if (!FinitePoint<D>(q) || !(DistanceToSegment<D>(q, p0, p1) <= tolerance.Chord)) {
            return false;
        }
    }
    const BasicPointN<double, D> middle = eval(s.Start + width / 2);
    if (!FinitePoint<D>(middle)) {
        return false;
    }
    // A segment within rounding of a point has no direction to turn by. At
    // a sphere's or an ellipsoid's pole the iso-line of constant latitude is
    // a circle of radius cos(pi/2) R, about 1e-16 R: its halves turned
    // sharply, every pole's line was refined until the subdivision gave up,
    // and no whole sphere could be drawn. The scale is the rounding of the
    // points' coordinates, so the angle still bounds every segment with a
    // length of its own, however loose the chord tolerance.
    const BasicVecN<double, D> first = middle - p0;
    const BasicVecN<double, D> second = p1 - middle;
    const double scale = std::max({1.0, MaxAbs<D>(p0), MaxAbs<D>(p1)});
    if (Length(first) + Length(second) <= RoundingSpan * scale) {
        return true;
    }
    return Turn<D>(first, second) <= tolerance.Angle;
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
        if (!IsFlat<D>(eval, s, p0, p1, tolerance)) {
            if (!splittable) {
                return false;
            }
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
                            std::string(function) +
                                ": the display tolerance could not be met within "
                                "the point budget or subdivision limit");
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

namespace {

/// Where the line of constant `fixed` in coordinate `across` (0 for u, 1
/// for v) crosses the loops, in the other coordinate, ascending. A segment
/// counts when its ends lie on opposite sides, one end on the line counted
/// with the side below, so a vertex on the line is crossed once and a
/// segment along it not at all.
[[nodiscard]] std::vector<double> Crossings(std::span<const std::vector<Point2>> loops,
                                            const int across, const double fixed)
{
    const auto fixedOf = [&](const Point2& p) { return across == 0 ? p.X : p.Y; };
    const auto alongOf = [&](const Point2& p) { return across == 0 ? p.Y : p.X; };
    std::vector<double> out;
    for (const std::vector<Point2>& loop : loops) {
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const Point2& p = loop[i];
            const Point2& q = loop[(i + 1) % loop.size()];
            if ((fixedOf(p) <= fixed) == (fixedOf(q) <= fixed)) {
                continue;
            }
            const double s = (fixed - fixedOf(p)) / (fixedOf(q) - fixedOf(p));
            out.push_back(alongOf(p) + s * (alongOf(q) - alongOf(p)));
        }
    }
    std::ranges::sort(out);
    return out;
}

/// The stretches of a line inside the loops, by the even-odd rule: the
/// crossings in pairs, each kept inside `range`.
[[nodiscard]] std::vector<Interval> Inside(const std::vector<double>& crossings,
                                           const Interval range)
{
    std::vector<Interval> out;
    for (std::size_t i = 0; i + 1 < crossings.size(); i += 2) {
        const double lower = std::max(crossings[i], range.GetLower());
        const double upper = std::min(crossings[i + 1], range.GetUpper());
        if (lower < upper) {
            out.push_back(Interval{lower, upper});
        }
    }
    return out;
}

/// The values strictly inside the range.
[[nodiscard]] std::vector<double> Interior(std::vector<double> values, const Interval range)
{
    std::erase_if(values,
                  [&](const double x) { return !(range.GetLower() < x && x < range.GetUpper()); });
    return values;
}

/// InvalidInput unless there is a surface, a bounded range, two lines per
/// direction and loops of three points each.
[[nodiscard]] Result<void> CheckTrimmedInput(const Surface& surface,
                                             std::span<const std::vector<Point2>> loops,
                                             const Box2& range, const std::size_t linesPerDirection)
{
    if (surface.IsEmpty() || linesPerDirection < 2 || range.IsEmpty() || !range.X.IsBounded() ||
        !range.Y.IsBounded()) {
        return Fail(ErrorCode::InvalidInput,
                    "SampleTrimmedSurface: needs a surface, a bounded range and at least two lines "
                    "per direction");
    }
    const bool tooShort =
        std::ranges::any_of(loops, [](const std::vector<Point2>& loop) { return loop.size() < 3; });
    if (tooShort) {
        return Fail(ErrorCode::InvalidInput,
                    "SampleTrimmedSurface: a loop needs at least three points");
    }
    return {};
}

/// The iso-lines across one direction, `across` 1 for lines along u at
/// fixed v and 0 for lines along v at fixed u: each fixed value's stretches
/// inside the loops, sampled and appended to `lines`.
struct IsoFamily {
    int Across;
    Interval Fixed;
    Interval Along;
    std::vector<double> FixedValues;
    const std::vector<double>* AlongKnots;
};

[[nodiscard]] Result<void> AppendTrimmed(const Surface& surface, const DisplayTolerance& tolerance,
                                         std::span<const std::vector<Point2>> loops,
                                         const IsoFamily& family, std::vector<IsoLine>& lines)
{
    const bool alongU = family.Across == 1;
    for (const double fixed : family.FixedValues) {
        for (const Interval stretch :
             Inside(Crossings(loops, family.Across, fixed), family.Along)) {
            const auto eval = [&](const double s) {
                return Eval(surface, alongU ? Point2{s, fixed} : Point2{fixed, s});
            };
            STEREON_TRY(Polyline line, SampleBetween<3>(eval, WithEnds(*family.AlongKnots, stretch),
                                                        tolerance, "SampleTrimmedSurface"));
            lines.push_back(MakeIsoLine(alongU ? IsoDirection::AlongU : IsoDirection::AlongV, fixed,
                                        std::move(line)));
        }
    }
    return {};
}

} // namespace

Result<std::vector<IsoLine>> SampleTrimmedSurface(const Surface& surface,
                                                  const DisplayTolerance& tolerance,
                                                  const std::span<const std::vector<Point2>> loops,
                                                  const Box2 range,
                                                  const std::size_t linesPerDirection)
{
    STEREON_TRY_VOID(CheckTrimmedInput(surface, loops, range, linesPerDirection));
    STEREON_TRY_VOID(CheckTolerance(tolerance, "SampleTrimmedSurface"));
    std::vector<double> knotsU;
    std::vector<double> knotsV;
    SurfaceKnots(surface, knotsU, knotsV, 0);
    // Along u at each interior v, then along v at each interior u.
    const IsoFamily families[] = {
        {1, range.Y, range.X, Interior(IsoValues(range.Y, linesPerDirection, knotsV), range.Y),
         &knotsU},
        {0, range.X, range.Y, Interior(IsoValues(range.X, linesPerDirection, knotsU), range.X),
         &knotsV},
    };
    std::vector<IsoLine> lines;
    for (const IsoFamily& family : families) {
        STEREON_TRY_VOID(AppendTrimmed(surface, tolerance, loops, family, lines));
    }
    return lines;
}

} // namespace Stereon
