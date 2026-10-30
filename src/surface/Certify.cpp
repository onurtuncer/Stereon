// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

module Stereon.Surface;

import :Certify;
import :Cylinder;
import :Plane;

namespace Stereon {

// Instantiated at namespace scope before any function uses them
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));
static_assert(Scalar<IntervalScalar>);

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();
constexpr double TwoPi = 2.0 * std::numbers::pi;

/// [x, x] widened by two units in the last place each way: more than the
/// libm functions' rounding.
[[nodiscard]] Interval Around(const double x)
{
    const double lower = std::nextafter(std::nextafter(x, -Infinity), -Infinity);
    const double upper = std::nextafter(std::nextafter(x, Infinity), Infinity);
    return {lower, upper};
}

/// Whether some point (phase + 2 k) pi, k an integer, lies within x, with
/// pi as an enclosure: the extrema of the trigonometric functions. x is
/// narrower than one period, so a few k suffice.
[[nodiscard]] bool MeetsLattice(const Interval& x, const double phase)
{
    static const Interval pi = Around(std::numbers::pi);
    const double first = std::floor((x.GetLower() / std::numbers::pi - phase) / 2.0) - 1.0;
    const double last = std::ceil((x.GetUpper() / std::numbers::pi - phase) / 2.0) + 1.0;
    for (double k = first; k <= last; k += 1.0) {
        if (Intersects(Interval{phase + 2.0 * k} * pi, x)) {
            return true;
        }
    }
    return false;
}

/// cos over x: the values at the ends, and +-1 where an extremum lies
/// within. A span of a period or more is the whole range.
[[nodiscard]] Interval CosOver(const Interval& x)
{
    if (!(x.GetUpper() - x.GetLower() < TwoPi)) {
        return {-1.0, 1.0};
    }
    Interval r = Hull(Around(std::cos(x.GetLower())), Around(std::cos(x.GetUpper())));
    if (MeetsLattice(x, 0.0)) {
        r = Hull(r, 1.0);
    }
    if (MeetsLattice(x, 1.0)) {
        r = Hull(r, -1.0);
    }
    return Intersection(r, Interval{-1.0, 1.0});
}

/// sin over x, as CosOver with the extrema a quarter turn on.
[[nodiscard]] Interval SinOver(const Interval& x)
{
    if (!(x.GetUpper() - x.GetLower() < TwoPi)) {
        return {-1.0, 1.0};
    }
    Interval r = Hull(Around(std::sin(x.GetLower())), Around(std::sin(x.GetUpper())));
    if (MeetsLattice(x, 0.5)) {
        r = Hull(r, 1.0);
    }
    if (MeetsLattice(x, 1.5)) {
        r = Hull(r, -1.0);
    }
    return Intersection(r, Interval{-1.0, 1.0});
}

} // namespace

IntervalScalar::IntervalScalar()
    : Value(0.0)
{}

IntervalScalar::IntervalScalar(const double x)
    : Value(x)
{}

IntervalScalar::IntervalScalar(const Interval& value)
    : Value(value)
{}

bool operator<(const IntervalScalar& a, const IntervalScalar& b)
{
    return a.Value.GetUpper() < b.Value.GetLower();
}

bool operator>(const IntervalScalar& a, const IntervalScalar& b)
{
    return b < a;
}

IntervalScalar operator+(const IntervalScalar& a, const IntervalScalar& b)
{
    return IntervalScalar{a.Value + b.Value};
}

IntervalScalar operator-(const IntervalScalar& a, const IntervalScalar& b)
{
    return IntervalScalar{a.Value - b.Value};
}

IntervalScalar operator-(const IntervalScalar& a)
{
    return IntervalScalar{-a.Value};
}

IntervalScalar operator*(const IntervalScalar& a, const IntervalScalar& b)
{
    return IntervalScalar{a.Value * b.Value};
}

IntervalScalar operator/(const IntervalScalar& a, const IntervalScalar& b)
{
    return IntervalScalar{a.Value / b.Value};
}

IntervalScalar sin(const IntervalScalar& x)
{
    return IntervalScalar{SinOver(x.Value)};
}

IntervalScalar cos(const IntervalScalar& x)
{
    return IntervalScalar{CosOver(x.Value)};
}

IntervalScalar sqrt(const IntervalScalar& x)
{
    const Interval nonNegative = Intersection(x.Value, Interval{0.0, Infinity});
    return nonNegative.IsEmpty() ? IntervalScalar{} : IntervalScalar{Sqrt(nonNegative)};
}

namespace {

using IPoint2 = BasicPoint2<IntervalScalar>;
using IVec2 = BasicVec2<IntervalScalar>;
using IPoint3 = BasicPoint3<IntervalScalar>;
using IVec3 = BasicVec3<IntervalScalar>;

// --- The geometry over intervals ----------------------------------------

[[nodiscard]] IPoint2 Wide(const Point2& p)
{
    return {IntervalScalar{p.X}, IntervalScalar{p.Y}};
}
[[nodiscard]] IVec2 Wide(const Vec2& v)
{
    return {IntervalScalar{v.X}, IntervalScalar{v.Y}};
}
[[nodiscard]] IPoint3 Wide(const Point3& p)
{
    return {IntervalScalar{p.X}, IntervalScalar{p.Y}, IntervalScalar{p.Z}};
}
[[nodiscard]] IVec3 Wide(const Vec3& v)
{
    return {IntervalScalar{v.X}, IntervalScalar{v.Y}, IntervalScalar{v.Z}};
}

template <std::size_t D>
[[nodiscard]] BasicLine<IntervalScalar, D> Enclosed(const BasicLine<double, D>& l)
{
    return {Wide(l.Origin), Wide(l.Direction)};
}
template <std::size_t D>
[[nodiscard]] BasicCircle<IntervalScalar, D> Enclosed(const BasicCircle<double, D>& c)
{
    return {Wide(c.Center), IntervalScalar{c.Radius}, Wide(c.XAxis), Wide(c.YAxis)};
}
template <std::size_t D>
[[nodiscard]] BasicEllipse<IntervalScalar, D> Enclosed(const BasicEllipse<double, D>& e)
{
    return {Wide(e.Center), IntervalScalar{e.MajorRadius}, IntervalScalar{e.MinorRadius},
            Wide(e.XAxis), Wide(e.YAxis)};
}
[[nodiscard]] BasicPlane<IntervalScalar> Enclosed(const Plane& p)
{
    return {Wide(p.Origin), Wide(p.XAxis), Wide(p.YAxis)};
}
[[nodiscard]] BasicCylinder<IntervalScalar> Enclosed(const Cylinder& c)
{
    return {Wide(c.Origin), IntervalScalar{c.Radius}, Wide(c.XAxis), Wide(c.YAxis), c.RightHanded};
}
template <std::size_t D>
[[nodiscard]] BasicBSplineCurve<IntervalScalar, D> Enclosed(const BasicBSplineCurve<double, D>& c)
{
    BasicBSplineCurve<IntervalScalar, D> wide;
    wide.Degree = c.Degree;
    wide.Knots = c.Knots;
    wide.Poles.reserve(c.Poles.size());
    for (const auto& pole : c.Poles) {
        wide.Poles.push_back(Wide(pole));
    }
    wide.Weights.assign(c.Weights.begin(), c.Weights.end());
    return wide;
}

// --- Traces: what the enclosure needs of each operand -----------------

/// The highest derivative the certificate takes: the Taylor form is of the
/// third order with its remainder from the fourth derivative.
constexpr int TaylorOrder = 4;

/// A curve's position vector and derivatives, [n] the n-th, enclosed at a
/// parameter or over a span.
template <std::size_t D>
using CurveJet = std::array<BasicVecN<IntervalScalar, D>, TaylorOrder + 1>;

/// A curve's jets for a span: at its midpoint, to the order TaylorOrder - 1,
/// and over it, to TaylorOrder.
template <std::size_t D>
struct SpanJets {
    CurveJet<D> AtMidpoint;
    CurveJet<D> Over;
};

/// A curve in D dimensions: its jets for a span, and the parameters where it
/// may be less smooth than they need (a B-spline's knots; none for the
/// closed forms). The jets are evaluated together, since each order over
/// the span is narrowed by the one above it and the midpoint's.
template <std::size_t D>
struct Trace {
    std::function<SpanJets<D>(const Interval&)> Jets;
    std::vector<double> Breaks;
};

/// The 3D curve and the pcurve.
using CurveTrace = Trace<3>;
using PcurveTrace = Trace<2>;

/// The surface's point and partial derivatives at enclosed parameters:
/// Partial[du][dv] for 1 <= du + dv <= the order asked for, at most
/// TaylorOrder, the rest zero.
struct SurfaceJet {
    IPoint3 Point;
    std::array<std::array<IVec3, TaylorOrder + 1>, TaylorOrder + 1> Partial;
};

/// The surface: its jet at enclosed parameters, to an order.
struct SurfaceTrace {
    std::function<SurfaceJet(const IPoint2&, int)> Jet;
};

/// The jet of a closed form at x, a parameter or a span, to the order `top`.
template <std::size_t D, typename G>
[[nodiscard]] CurveJet<D> ClosedFormJet(const G& g, const IntervalScalar& x, const int top)
{
    CurveJet<D> jet{};
    jet[0] = ToVec(Eval(g, x));
    for (int n = 1; n <= top; ++n) {
        jet[static_cast<std::size_t>(n)] = Deriv(g, x, n);
    }
    return jet;
}

/// A circle's and an ellipse's jets from one enclosure of cos x and sin x:
/// their derivatives turn the radii by quarter turns, as Deriv computes
/// them, and each evaluated on its own would enclose the same cos and sin
/// again.
template <std::size_t D>
[[nodiscard]] CurveJet<D> ClosedFormJet(const BasicCircle<IntervalScalar, D>& g,
                                        const IntervalScalar& x, const int top)
{
    const IntervalScalar c = cos(x);
    const IntervalScalar s = sin(x);
    CurveJet<D> jet{};
    for (int n = 0; n <= top; ++n) {
        const BasicCosSin<IntervalScalar> turned = QuarterTurned(c, s, n);
        jet[static_cast<std::size_t>(n)] = g.Radius * (turned.Cos * g.XAxis + turned.Sin * g.YAxis);
    }
    jet[0] = jet[0] + ToVec(g.Center);
    return jet;
}

template <std::size_t D>
[[nodiscard]] CurveJet<D> ClosedFormJet(const BasicEllipse<IntervalScalar, D>& g,
                                        const IntervalScalar& x, const int top)
{
    const IntervalScalar c = cos(x);
    const IntervalScalar s = sin(x);
    CurveJet<D> jet{};
    for (int n = 0; n <= top; ++n) {
        const BasicCosSin<IntervalScalar> turned = QuarterTurned(c, s, n);
        jet[static_cast<std::size_t>(n)] =
            (g.MajorRadius * turned.Cos) * g.XAxis + (g.MinorRadius * turned.Sin) * g.YAxis;
    }
    jet[0] = jet[0] + ToVec(g.Center);
    return jet;
}

template <std::size_t D, typename G>
[[nodiscard]] Trace<D> ClosedFormTrace(const G& g)
{
    return {[g](const Interval& span) {
                const double m = 0.5 * (span.GetLower() + span.GetUpper());
                return SpanJets<D>{ClosedFormJet<D>(g, IntervalScalar{m}, TaylorOrder - 1),
                                   ClosedFormJet<D>(g, IntervalScalar{span}, TaylorOrder)};
            },
            {}};
}

template <typename G>
[[nodiscard]] SurfaceTrace SurfaceOf(const G& g)
{
    return {[g](const IPoint2& uv, const int top) {
        SurfaceJet jet{Eval(g, uv), {}};
        for (int du = 0; du <= top; ++du) {
            for (int dv = du == 0 ? 1 : 0; du + dv <= top; ++dv) {
                jet.Partial[static_cast<std::size_t>(du)][static_cast<std::size_t>(dv)] =
                    Deriv(g, uv, du, dv);
            }
        }
        return jet;
    }};
}

/// A cylinder's jet from one enclosure of cos u and sin u: its u-partials
/// are the radius turned by quarter turns, as Deriv computes them, and the
/// axis its only other partial; each partial evaluated on its own would
/// enclose the same cos and sin again.
[[nodiscard]] SurfaceTrace SurfaceOf(const BasicCylinder<IntervalScalar>& g)
{
    return {[g](const IPoint2& uv, const int top) {
        const IntervalScalar c = cos(uv.X);
        const IntervalScalar s = sin(uv.X);
        const IVec3 axis = Axis(g);
        const auto radial = [&](const int n) {
            const BasicCosSin<IntervalScalar> turned = QuarterTurned(c, s, n);
            return g.Radius * (turned.Cos * g.XAxis + turned.Sin * g.YAxis);
        };
        SurfaceJet jet{g.Origin + radial(0) + uv.Y * axis, {}};
        for (std::size_t du = 1; du <= static_cast<std::size_t>(top); ++du) {
            jet.Partial[du][0] = radial(static_cast<int>(du));
        }
        jet.Partial[0][1] = axis;
        return jet;
    }};
}

// --- B-splines over spans ------------------------------------------------

/// t - m for every t in the span, rounded outward: the differences of the
/// ends need not be exact when the span holds 0.
[[nodiscard]] Interval OffsetOver(const Interval& span, const double m)
{
    return {(Interval{span.GetLower()} - Interval{m}).GetLower(),
            (Interval{span.GetUpper()} - Interval{m}).GetUpper()};
}

/// The derivative curve of a non-rational B-spline (The NURBS Book, eq.
/// 3.8): degree p - 1 on the knots less their first and last, with poles
/// p (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1}), outward-rounded. Nothing for
/// a knot difference that is not positive, which a valid knot vector of
/// degree >= 1 never has inside its domain.
template <std::size_t D>
[[nodiscard]] std::optional<BasicBSplineCurve<IntervalScalar, D>>
DerivativeCurve(const BasicBSplineCurve<IntervalScalar, D>& c)
{
    const auto p = static_cast<std::size_t>(c.Degree);
    if (p == 0 || c.Poles.size() < 2) {
        return std::nullopt;
    }
    BasicBSplineCurve<IntervalScalar, D> d;
    d.Degree = c.Degree - 1;
    d.Knots.assign(c.Knots.begin() + 1, c.Knots.end() - 1);
    d.Poles.reserve(c.Poles.size() - 1);
    const IntervalScalar degree{static_cast<double>(p)};
    for (std::size_t i = 0; i + 1 < c.Poles.size(); ++i) {
        const Interval step = Interval{c.Knots[i + p + 1]} - Interval{c.Knots[i + 1]};
        if (!(step.GetLower() > 0.0)) {
            return std::nullopt;
        }
        const auto difference = ToVec(c.Poles[i + 1]) - ToVec(c.Poles[i]);
        d.Poles.push_back(ToPoint(difference * (degree / IntervalScalar{step})));
    }
    return d;
}

/// The largest magnitude of a pole coordinate, for the rounding margin.
template <std::size_t D>
[[nodiscard]] double PoleMagnitude(const BasicBSplineCurve<IntervalScalar, D>& c)
{
    double magnitude = 0.0;
    for (const auto& pole : c.Poles) {
        for (std::size_t i = 0; i < D; ++i) {
            magnitude = std::max({magnitude, std::abs(pole[i].Value.GetLower()),
                                  std::abs(pole[i].Value.GetUpper())});
        }
    }
    return magnitude;
}

/// The hull of the poles active over the span: every point of the curve
/// over it lies within (the convex hull property, per knot span).
template <std::size_t D>
[[nodiscard]] BasicPointN<IntervalScalar, D> HullOver(const BasicBSplineCurve<IntervalScalar, D>& c,
                                                      const Interval& span)
{
    const Interval domain = Domain(c);
    const double lower = std::clamp(span.GetLower(), domain.GetLower(), domain.GetUpper());
    const double upper = std::clamp(span.GetUpper(), domain.GetLower(), domain.GetUpper());
    const std::size_t first = FindSpan(c.Knots, c.Degree, c.Poles.size(), lower);
    const std::size_t last = FindSpan(c.Knots, c.Degree, c.Poles.size(), upper);
    const auto p = static_cast<std::size_t>(c.Degree);
    BasicPointN<IntervalScalar, D> hull = c.Poles[first - p];
    for (std::size_t i = first - p; i <= last; ++i) {
        for (std::size_t k = 0; k < D; ++k) {
            hull[k] = IntervalScalar{Hull(hull[k].Value, c.Poles[i][k].Value)};
        }
    }
    return hull;
}

/// Eval's basis functions are rounded doubles, not enclosures: the point is
/// widened by the margin BSplineCurve's PoleBoxMargin derives,
/// 8 (p + 1)^2 units of 2^-53 of the largest pole coordinate.
template <std::size_t D>
[[nodiscard]] BasicPointN<IntervalScalar, D>
EnclosedPoint(const BasicBSplineCurve<IntervalScalar, D>& c, const double magnitude, const double t)
{
    const auto terms = static_cast<double>(c.Degree + 1);
    const double margin = 8.0 * terms * terms * 0x1p-53 * magnitude;
    BasicPointN<IntervalScalar, D> point = Eval(c, t);
    for (std::size_t k = 0; k < D; ++k) {
        point[k] = IntervalScalar{Inflated(point[k].Value, margin)};
    }
    return point;
}

//// A non-rational B-spline with its derivative curves up to TaylorOrder, all
/// over intervals, and what it encloses over a parameter span. The hull of
/// the poles active over a span holds for any span but is as wide as the
/// knot spans it covers; within one knot span the spline is a polynomial, so
/// a mean-value form from the next derivative, p(m) + p'(span) (t - m),
/// narrows with the span. Both are enclosures, so their meet is one.
template <std::size_t D>
struct SplineEnclosure {
    /// Derivs[n], the n-th derivative curve, Derivs[0] the spline itself;
    /// none beyond the degree, where the derivative is zero within a knot
    /// span.
    std::vector<BasicBSplineCurve<IntervalScalar, D>> Derivs;
    /// The largest magnitude of each derivative curve's poles.
    std::vector<double> Magnitudes;

    /// True when the span lies within one knot span of the spline.
    [[nodiscard]] bool WithinOneSpan(const Interval& span) const
    {
        const BasicBSplineCurve<IntervalScalar, D>& c = Derivs.front();
        const Interval domain = Domain(c);
        const double lower = std::clamp(span.GetLower(), domain.GetLower(), domain.GetUpper());
        const double upper = std::clamp(span.GetUpper(), domain.GetLower(), domain.GetUpper());
        return FindSpan(c.Knots, c.Degree, c.Poles.size(), lower) ==
               FindSpan(c.Knots, c.Degree, c.Poles.size(), upper);
    }

    /// The n-th derivative at t; zero beyond the degree.
    [[nodiscard]] BasicVecN<IntervalScalar, D> At(const double t, const std::size_t n) const
    {
        if (n >= Derivs.size()) {
            return {};
        }
        // Of degree 0, the derivative is its active pole, exactly; Eval
        // takes degree 1 and above.
        return Derivs[n].Degree == 0 ? ToVec(HullOver(Derivs[n], Interval{t}))
                                     : ToVec(EnclosedPoint(Derivs[n], Magnitudes[n], t));
    }

    /// The jet at t, to the order TaylorOrder - 1.
    [[nodiscard]] CurveJet<D> JetAt(const double t) const
    {
        CurveJet<D> jet{};
        for (std::size_t n = 0; n < static_cast<std::size_t>(TaylorOrder); ++n) {
            jet[n] = At(t, n);
        }
        return jet;
    }

    /// The jets for the span: at its midpoint m, and over it from the top
    /// order down, each order the hull of its active poles, met within one
    /// knot span with the mean-value form from m and the order above.
    [[nodiscard]] SpanJets<D> Jets(const Interval& span) const
    {
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        SpanJets<D> jets{JetAt(m), {}};
        const std::size_t top = std::min(Derivs.size() - 1, static_cast<std::size_t>(TaylorOrder));
        const bool within = WithinOneSpan(span);
        const IntervalScalar offset{OffsetOver(span, m)};
        for (std::size_t n = top + 1; n-- > 0;) {
            jets.Over[n] = ToVec(HullOver(Derivs[n], span));
            if (within && n < top) {
                // n < top <= TaylorOrder, an order the midpoint's jet has.
                jets.Over[n] = Meet(jets.Over[n], jets.AtMidpoint[n] + jets.Over[n + 1] * offset);
            }
        }
        return jets;
    }

    /// The distinct knots: where the spline may be less smooth.
    [[nodiscard]] std::vector<double> Breaks() const
    {
        const auto& knots = Derivs.front().Knots;
        std::vector<double> breaks(knots.begin(), knots.end());
        breaks.erase(std::ranges::unique(breaks).begin(), breaks.end());
        return breaks;
    }

    /// The componentwise meet of two enclosures of the same vector; the
    /// first when they do not meet, which two rigorous enclosures never do.
    [[nodiscard]] static BasicVecN<IntervalScalar, D> Meet(BasicVecN<IntervalScalar, D> a,
                                                           const BasicVecN<IntervalScalar, D>& b)
    {
        for (std::size_t k = 0; k < D; ++k) {
            const Interval met = Intersection(a[k].Value, b[k].Value);
            if (!met.IsEmpty()) {
                a[k] = IntervalScalar{met};
            }
        }
        return a;
    }
};

template <std::size_t D>
[[nodiscard]] std::optional<Trace<D>> SplineTrace(const BasicBSplineCurve<double, D>& spline)
{
    if (IsRational(spline)) {
        return std::nullopt;
    }
    SplineEnclosure<D> e;
    e.Derivs.push_back(Enclosed(spline));
    while (e.Derivs.size() <= static_cast<std::size_t>(TaylorOrder)) {
        std::optional<BasicBSplineCurve<IntervalScalar, D>> next = DerivativeCurve(e.Derivs.back());
        if (!next) {
            break;
        }
        e.Derivs.push_back(std::move(*next));
    }
    if (e.Derivs.size() < 2) {
        return std::nullopt;
    }
    for (const BasicBSplineCurve<IntervalScalar, D>& c : e.Derivs) {
        e.Magnitudes.push_back(PoleMagnitude(c));
    }
    const auto c = std::make_shared<const SplineEnclosure<D>>(std::move(e));
    return Trace<D>{[c](const Interval& span) { return c->Jets(span); }, c->Breaks()};
}

// --- The erased operands --------------------------------------------------

[[nodiscard]] std::optional<CurveTrace> TraceOf(const Curve& curve)
{
    if (const Line* line = curve.Target<Line>()) {
        return ClosedFormTrace<3>(Enclosed(*line));
    }
    if (const Circle* circle = curve.Target<Circle>()) {
        return ClosedFormTrace<3>(Enclosed(*circle));
    }
    if (const Ellipse* ellipse = curve.Target<Ellipse>()) {
        return ClosedFormTrace<3>(Enclosed(*ellipse));
    }
    if (const BSplineCurve* spline = curve.Target<BSplineCurve>()) {
        return SplineTrace(*spline);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<PcurveTrace> TraceOf(const Curve2d& pcurve)
{
    if (const Line2d* line = pcurve.Target<Line2d>()) {
        return ClosedFormTrace<2>(Enclosed(*line));
    }
    if (const Circle2d* circle = pcurve.Target<Circle2d>()) {
        return ClosedFormTrace<2>(Enclosed(*circle));
    }
    if (const Ellipse2d* ellipse = pcurve.Target<Ellipse2d>()) {
        return ClosedFormTrace<2>(Enclosed(*ellipse));
    }
    if (const BSplineCurve2d* spline = pcurve.Target<BSplineCurve2d>()) {
        return SplineTrace(*spline);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<SurfaceTrace> TraceOf(const Surface& surface)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        return SurfaceOf(Enclosed(*plane));
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        return SurfaceOf(Enclosed(*cylinder));
    }
    return std::nullopt;
}

// --- The certificate ---------------------------------------------------

/// The least magnitude of a number in the interval: 0 when it holds 0.
[[nodiscard]] double LeastMagnitude(const Interval& i)
{
    if (Contains(i, 0.0)) {
        return 0.0;
    }
    return std::min(std::abs(i.GetLower()), std::abs(i.GetUpper()));
}

/// The length of a vector of enclosures: an upper bound, rounded up, and a
/// lower bound, rounded down. Outward rounding of the sum of squares can
/// put its lower end a step below zero, so the roots are taken of the
/// doubles rather than through Interval's Sqrt.
[[nodiscard]] Interval Length(const IVec3& v)
{
    Interval squares{0.0};
    Interval leastSquares{0.0};
    for (const IntervalScalar& component : {v.X, v.Y, v.Z}) {
        squares = squares + Squared(component.Value);
        leastSquares = leastSquares + Squared(Interval{LeastMagnitude(component.Value)});
    }
    const double upper = std::nextafter(std::sqrt(squares.GetUpper()), Infinity);
    const double lower =
        std::nextafter(std::sqrt(std::max(0.0, leastSquares.GetLower())), -Infinity);
    return {std::max(0.0, lower), upper};
}

///// A truncated Taylor series in s, [k] the coefficient of s^k.
using Series = std::array<IntervalScalar, TaylorOrder + 1>;

/// 1 / (i! j!) enclosed, for i + j <= TaylorOrder.
[[nodiscard]] const std::array<std::array<IntervalScalar, TaylorOrder + 1>, TaylorOrder + 1>&
ReciprocalFactorials()
{
    static const auto table = [] {
        std::array<std::array<IntervalScalar, TaylorOrder + 1>, TaylorOrder + 1> r{};
        double fi = 1.0;
        for (std::size_t i = 0; i <= static_cast<std::size_t>(TaylorOrder); ++i) {
            fi *= i == 0 ? 1.0 : static_cast<double>(i);
            double fj = 1.0;
            for (std::size_t j = 0; i + j <= static_cast<std::size_t>(TaylorOrder); ++j) {
                fj *= j == 0 ? 1.0 : static_cast<double>(j);
                r[i][j] = IntervalScalar{Interval{1.0} / Interval{fi * fj}};
            }
        }
        return r;
    }();
    return table;
}

/// The coefficients of a curve's Taylor series: its jet's [k] over k!.
template <typename V>
[[nodiscard]] std::array<V, TaylorOrder + 1> Coefficients(std::array<V, TaylorOrder + 1> jet)
{
    const auto& reciprocal = ReciprocalFactorials();
    for (std::size_t k = 2; k <= static_cast<std::size_t>(TaylorOrder); ++k) {
        jet[k] = jet[k] * reciprocal[k][0];
    }
    return jet;
}

/// a b to the order `top`, where a's coefficients below `lowA` and b's below
/// `lowB` are zero.
[[nodiscard]] Series Times(const Series& a, const std::size_t lowA, const Series& b,
                           const std::size_t lowB, const std::size_t top)
{
    Series product{};
    for (std::size_t i = lowA; i + lowB <= top; ++i) {
        for (std::size_t j = lowB; i + j <= top; ++j) {
            product[i + j] = product[i + j] + a[i] * b[j];
        }
    }
    return product;
}

[[nodiscard]] bool IsZero(const IVec3& v)
{
    return v.X.Value == Interval{0.0} && v.Y.Value == Interval{0.0} && v.Z.Value == Interval{0.0};
}

/// The highest powers of du and dv that the partials up to the order `top`
/// need: those of the partials that are not exactly zero.
[[nodiscard]] std::pair<std::size_t, std::size_t> HighestPowers(const SurfaceJet& jet,
                                                                const std::size_t top)
{
    std::size_t highU = 0;
    std::size_t highV = 0;
    for (std::size_t i = 0; i <= top; ++i) {
        for (std::size_t j = 0; i + j <= top; ++j) {
            if (i + j > 0 && !IsZero(jet.Partial[i][j])) {
                highU = std::max(highU, i);
                highV = std::max(highV, j);
            }
        }
    }
    return {highU, highV};
}

/// The Taylor coefficients g_k = g^(k) / k!, from <= k <= top, of g(s) =
/// S(q(s)), from the surface's jet at q(0) and the coefficients q_k of the
/// pcurve's series: S(q0 + d) is the sum over i + j >= 1 of
/// S_ij du^i dv^j / (i! j!), du and dv the series of q less q0, which is the
/// chain rule to every order at once (Faa di Bruno's formula). It holds at
/// each parameter, so the jet and the coefficients enclosed over a span
/// enclose the g_k at every parameter of it. Partials that are exactly zero,
/// most of a cylinder's and a plane's, are skipped, and the powers of du and
/// dv are raised only as far as the partials left need.
[[nodiscard]] std::array<IVec3, TaylorOrder + 1>
Composed(const SurfaceJet& jet, const CurveJet<2>& q, const std::size_t from, const std::size_t top)
{
    const auto [highU, highV] = HighestPowers(jet, top);
    std::array<Series, TaylorOrder + 1> du{};
    std::array<Series, TaylorOrder + 1> dv{};
    du[0][0] = IntervalScalar{1.0};
    dv[0][0] = IntervalScalar{1.0};
    for (std::size_t k = 1; k <= top; ++k) {
        du[1][k] = q[k].X;
        dv[1][k] = q[k].Y;
    }
    for (std::size_t i = 2; i <= highU; ++i) {
        du[i] = Times(du[i - 1], i - 1, du[1], 1, top);
    }
    for (std::size_t j = 2; j <= highV; ++j) {
        dv[j] = Times(dv[j - 1], j - 1, dv[1], 1, top);
    }
    const auto& reciprocal = ReciprocalFactorials();
    std::array<IVec3, TaylorOrder + 1> g{};
    for (std::size_t i = 0; i <= highU; ++i) {
        for (std::size_t j = 0; j <= highV && i + j <= top; ++j) {
            if (i + j == 0 || IsZero(jet.Partial[i][j])) {
                continue;
            }
            // du^i dv^j: one of the factors is 1 unless the partial is mixed.
            const Series w = i == 0 ? dv[j] : j == 0 ? du[i] : Times(du[i], i, dv[j], j, top);
            const IVec3 scaled = jet.Partial[i][j] * reciprocal[i][j];
            for (std::size_t k = std::max(i + j, from); k <= top; ++k) {
                g[k] = g[k] + scaled * w[k];
            }
        }
    }
    return g;
}

/// What a span proves about |c(t) - S(p(t))| over it.
struct SpanBounds {
    /// A bound over the whole span: the smaller of the mean-value form's and
    /// the Taylor form's, the second only where the first misses the target.
    double Upper;
    /// A deviation that certainly occurs, at the midpoint: the least
    /// magnitude of f(m)'s enclosure. A witness that no certificate below it
    /// exists.
    double Witness;
};

/// Bounds on f(t) = c(t) - S(p(t)) over the span from its midpoint m, s = t - m.
///
/// The mean-value form, |f(m) + f'(span) s|, shrinks with the square of
/// the span. Where it misses the target, the Taylor form of the third order
/// with its remainder in integral form, which holds for vectors:
///     sum over k <= 3 of |f_k(m)| |s|^k + sup |f_4| |s|^4,
/// f_k = f^(k) / k!. The coefficients up to the third are evaluated at m,
/// where c and S(p) nearly cancel and their enclosures are tight; only the
/// remainder is enclosed over the span. That matters: f's derivatives are
/// differences of the curve's and the surface's, each varying over the
/// span by much more than their difference, so whichever order is enclosed
/// over the span carries that variation. Enclosing the third, as the form
/// of issue #155 first did, left a remainder of about 0.5 r^4 on the
/// ellipse's fitted pcurve, just over the target at a knot span's width, and
/// every knot span was halved; with the fourth enclosed, a knot span
/// certifies whole. Valid where c and S(p) are four times differentiable
/// over the span: within one polynomial piece of each B-spline, which
/// Bisect's seeding at the knots guarantees.
[[nodiscard]] SpanBounds BoundsOver(const CurveTrace& c, const PcurveTrace& p,
                                    const SurfaceTrace& s, const Interval& span,
                                    const double target)
{
    constexpr auto order = static_cast<std::size_t>(TaylorOrder);
    const double m = 0.5 * (span.GetLower() + span.GetUpper());
    const SpanJets<3> curve = c.Jets(span);
    const SpanJets<2> pcurve = p.Jets(span);
    const CurveJet<3> curveAt = Coefficients(curve.AtMidpoint);
    const CurveJet<2> pcurveAt = Coefficients(pcurve.AtMidpoint);
    const SurfaceJet at = s.Jet(ToPoint(pcurveAt[0]), TaylorOrder - 1);
    const IVec3 f0 = curveAt[0] - ToVec(at.Point);
    const CurveJet<3> curveOver = Coefficients(curve.Over);
    const CurveJet<2> pcurveOver = Coefficients(pcurve.Over);
    const SurfaceJet over = s.Jet(ToPoint(pcurveOver[0]), TaylorOrder);
    const IVec3 slope = curveOver[1] - (over.Partial[1][0] * pcurveOver[1].X +
                                        over.Partial[0][1] * pcurveOver[1].Y);
    const Interval offset = OffsetOver(span, m);
    SpanBounds bounds{Length(f0 + slope * IntervalScalar{offset}).GetUpper(),
                      Length(f0).GetLower()};
    if (bounds.Upper <= target) {
        return bounds;
    }
    const std::array<IVec3, TaylorOrder + 1> surfaceAt = Composed(at, pcurveAt, 1, order - 1);
    const IVec3 remainder = curveOver[order] - Composed(over, pcurveOver, order, order)[order];
    const Interval reach{std::max(std::abs(offset.GetLower()), std::abs(offset.GetUpper()))};
    Interval taylor{Length(f0).GetUpper()};
    Interval power{1.0};
    for (std::size_t k = 1; k <= order; ++k) {
        power = power * reach;
        const IVec3 fk = k < order ? curveAt[k] - surfaceAt[k] : remainder;
        taylor = taylor + Interval{Length(fk).GetUpper()} * power;
    }
    bounds.Upper = std::min(bounds.Upper, taylor.GetUpper());
    return bounds;
}

/// The breaks of either operand strictly inside the range, sorted, without
/// repeats: the seams of the starting spans.
[[nodiscard]] std::vector<double> BreaksWithin(const CurveTrace& c, const PcurveTrace& p,
                                               const Interval& range)
{
    std::vector<double> breaks;
    for (const std::vector<double>* list : {&c.Breaks, &p.Breaks}) {
        for (const double b : *list) {
            if (range.GetLower() < b && b < range.GetUpper()) {
                breaks.push_back(b);
            }
        }
    }
    std::ranges::sort(breaks);
    breaks.erase(std::ranges::unique(breaks).begin(), breaks.end());
    return breaks;
}

[[nodiscard]] std::string Number(const double x)
{
    std::array<char, 32> buffer{};
    const std::to_chars_result r = std::to_chars(buffer.data(), buffer.data() + buffer.size(), x);
    return std::string(buffer.data(), r.ptr);
}

/// Spans bisected until each is within the target, starting from the range
/// cut at every break of the curve and the pcurve, so that each span lies
/// in one polynomial piece of both and the Taylor form holds on it; each
/// span's bound is the smaller of that form's and the mean-value form's,
/// both rigorous. A stack, so the memory stays logarithmic in the number of
/// spans of each starting piece. A span whose midpoint deviates by
/// more than the target for certain ends the search at once: no
/// subdivision can certify what a point already refutes. A span that can
/// no longer be split (its midpoint is one of its ends) or the budget
/// running out ends it too, which is where a target below the rounding of
/// the enclosures leads.
[[nodiscard]] Result<DeviationCertificate> Bisect(const CurveTrace& c, const PcurveTrace& p,
                                                  const SurfaceTrace& s, const Interval& range,
                                                  const double target)
{
    constexpr std::size_t SpanBudget = std::size_t{1} << 18;
    // The starting pieces, pushed last first so that the first is checked
    // first, as when the range was one piece.
    std::vector<double> cuts = BreaksWithin(c, p, range);
    cuts.insert(cuts.begin(), range.GetLower());
    cuts.push_back(range.GetUpper());
    std::vector<Interval> work;
    for (std::size_t k = cuts.size() - 1; k > 0; --k) {
        work.emplace_back(cuts[k - 1], cuts[k]);
    }
    DeviationCertificate certificate;
    for (std::size_t spans = 0; !work.empty(); ++spans) {
        const Interval span = work.back();
        work.pop_back();
        const SpanBounds local = BoundsOver(c, p, s, span, target);
        if (local.Upper <= target) {
            certificate.Bound = std::max(certificate.Bound, local.Upper);
            ++certificate.Spans;
            continue;
        }
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        if (local.Witness > target) {
            return Fail(ErrorCode::CertificationFailed,
                        "CertifyDeviation: at t = " + Number(m) + " the deviation is at least " +
                            Number(local.Witness) + " mm, above the target " + Number(target) +
                            " mm");
        }
        if (spans >= SpanBudget || !(span.GetLower() < m && m < span.GetUpper())) {
            return Fail(ErrorCode::CertificationFailed,
                        "CertifyDeviation: over [" + Number(span.GetLower()) + ", " +
                            Number(span.GetUpper()) + "] the deviation is bounded by " +
                            Number(local.Upper) + " mm, above the target " + Number(target) +
                            " mm, after " + Number(static_cast<double>(spans)) + " spans");
        }
        work.emplace_back(m, span.GetUpper());
        work.emplace_back(span.GetLower(), m);
    }
    return certificate;
}

// --- The closed-form certificate ------------------------------------------

/// A curve in space of the form K0 + Kc cos t + Ks sin t + Kt t, with each
/// coefficient enclosed: what lines, circles and ellipses are, and what a
/// plane composed with any of their plane forms is, and a cylinder composed
/// with a ruling or a unit-speed turn. The deviation of two such curves is
/// the same form in the coefficient differences, so a bound over a whole
/// range needs no subdivision: the analytic identity the generic
/// mean-value forms rediscover span by span (issue #155) is used directly.
struct TrigForm {
    IVec3 Constant;
    IVec3 Cosine;
    IVec3 Sine;
    IVec3 Linear;
};

[[nodiscard]] IVec3 ToVector(const IPoint3& p)
{
    return {p.X, p.Y, p.Z};
}

[[nodiscard]] std::optional<TrigForm> TrigFormOf(const Curve& curve)
{
    if (const Line* line = curve.Target<Line>()) {
        const BasicLine<IntervalScalar, 3> l = Enclosed(*line);
        return TrigForm{ToVector(l.Origin), {}, {}, l.Direction};
    }
    if (const Circle* circle = curve.Target<Circle>()) {
        const BasicCircle<IntervalScalar, 3> c = Enclosed(*circle);
        return TrigForm{ToVector(c.Center), c.Radius * c.XAxis, c.Radius * c.YAxis, {}};
    }
    if (const Ellipse* ellipse = curve.Target<Ellipse>()) {
        const BasicEllipse<IntervalScalar, 3> e = Enclosed(*ellipse);
        return TrigForm{ToVector(e.Center), e.MajorRadius * e.XAxis, e.MinorRadius * e.YAxis, {}};
    }
    return std::nullopt;
}

/// The plane curve as K0 + Kc cos t + Ks sin t + Kt t in the plane.
struct PlaneTrigForm {
    IVec2 Constant;
    IVec2 Cosine;
    IVec2 Sine;
    IVec2 Linear;
};

[[nodiscard]] std::optional<PlaneTrigForm> PlaneTrigFormOf(const Curve2d& pcurve)
{
    if (const Line2d* line = pcurve.Target<Line2d>()) {
        const BasicLine<IntervalScalar, 2> l = Enclosed(*line);
        return PlaneTrigForm{{l.Origin.X, l.Origin.Y}, {}, {}, l.Direction};
    }
    if (const Circle2d* circle = pcurve.Target<Circle2d>()) {
        const BasicCircle<IntervalScalar, 2> c = Enclosed(*circle);
        return PlaneTrigForm{{c.Center.X, c.Center.Y}, c.Radius * c.XAxis, c.Radius * c.YAxis, {}};
    }
    if (const Ellipse2d* ellipse = pcurve.Target<Ellipse2d>()) {
        const BasicEllipse<IntervalScalar, 2> e = Enclosed(*ellipse);
        return PlaneTrigForm{
            {e.Center.X, e.Center.Y}, e.MajorRadius * e.XAxis, e.MinorRadius * e.YAxis, {}};
    }
    return std::nullopt;
}

/// The plane O + u X + v Y along a plane curve: affine, so each coefficient
/// maps through the plane's axes.
[[nodiscard]] TrigForm OnPlane(const Plane& plane, const PlaneTrigForm& p)
{
    const BasicPlane<IntervalScalar> s = Enclosed(plane);
    const auto map = [&s](const IVec2& w) { return w.X * s.XAxis + w.Y * s.YAxis; };
    return {ToVector(s.Origin) + map(p.Constant), map(p.Cosine), map(p.Sine), map(p.Linear)};
}

/// The cylinder O + R (cos u X + sin u Y) + v A along a line (u0 + a t,
/// v0 + b t) of its parameter plane, when that is of the form: a ruling,
/// a = 0, is O + R (cos u0 X + sin u0 Y) + v0 A + b t A; a unit-speed turn,
/// a = +-1 and b = 0, is O + v0 A + R (cos u0 X + sin u0 Y) cos t
/// + a R (cos u0 Y - sin u0 X) sin t. Other lines are helices, not of the
/// form, and are left to the generic certificate.
[[nodiscard]] std::optional<TrigForm> OnCylinder(const Cylinder& cylinder, const Line2d& line)
{
    const double a = line.Direction.X;
    const double b = line.Direction.Y;
    const BasicCylinder<IntervalScalar> s = Enclosed(cylinder);
    const IVec3 axis = Axis(s);
    const IntervalScalar u0{line.Origin.X};
    const IntervalScalar v0{line.Origin.Y};
    const IntervalScalar c = cos(u0);
    const IntervalScalar n = sin(u0);
    const IVec3 radial = s.Radius * (c * s.XAxis + n * s.YAxis);
    if (a == 0.0) { // STEREON-LINT-ALLOW(epsilon): the exact form, not a tolerance
        return TrigForm{ToVector(s.Origin) + radial + v0 * axis, {}, {}, IntervalScalar{b} * axis};
    }
    if ((a == 1.0 || a == -1.0) && b == 0.0) { // STEREON-LINT-ALLOW(epsilon): the exact form
        const IVec3 turned = s.Radius * (c * s.YAxis - n * s.XAxis);
        return TrigForm{ToVector(s.Origin) + v0 * axis, radial, IntervalScalar{a} * turned, {}};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<TrigForm> TrigFormOf(const Surface& surface, const Curve2d& pcurve)
{
    if (const Plane* plane = surface.Target<Plane>()) {
        if (const std::optional<PlaneTrigForm> p = PlaneTrigFormOf(pcurve)) {
            return OnPlane(*plane, *p);
        }
        return std::nullopt;
    }
    if (const Cylinder* cylinder = surface.Target<Cylinder>()) {
        if (const Line2d* line = pcurve.Target<Line2d>()) {
            return OnCylinder(*cylinder, *line);
        }
    }
    return std::nullopt;
}

/// A bound on sup |c(t) - S(p(t))| over the range when both sides are trig
/// forms: |dK0| + |dKc| + |dKs| + |dKt| max |t|, every term rounded up,
/// since |cos t| and |sin t| are at most 1. Nothing when either side is not
/// of the form.
[[nodiscard]] std::optional<double> ClosedFormBound(const Curve& curve, const Surface& surface,
                                                    const Curve2d& pcurve, const Interval& range)
{
    const std::optional<TrigForm> c = TrigFormOf(curve);
    if (!c) {
        return std::nullopt;
    }
    const std::optional<TrigForm> s = TrigFormOf(surface, pcurve);
    if (!s) {
        return std::nullopt;
    }
    const double reach = std::max(std::abs(range.GetLower()), std::abs(range.GetUpper()));
    const Interval total = Interval{Length(c->Constant - s->Constant).GetUpper()} +
                           Interval{Length(c->Cosine - s->Cosine).GetUpper()} +
                           Interval{Length(c->Sine - s->Sine).GetUpper()} +
                           Interval{Length(c->Linear - s->Linear).GetUpper()} * Interval{reach};
    return total.GetUpper();
}

/// A B-spline is evaluated only within its domain; the closed forms
/// anywhere.
[[nodiscard]] bool CoversRange(const Curve& curve, const Interval& range)
{
    const BSplineCurve* spline = curve.Target<BSplineCurve>();
    return spline == nullptr || Contains(Domain(*spline), range);
}

[[nodiscard]] bool CoversRange(const Curve2d& pcurve, const Interval& range)
{
    const BSplineCurve2d* spline = pcurve.Target<BSplineCurve2d>();
    return spline == nullptr || Contains(Domain(*spline), range);
}

} // namespace

Result<DeviationCertificate> CertifyDeviation(const Curve& curve, const Surface& surface,
                                              const Curve2d& pcurve, const Interval& range,
                                              const double target)
{
    if (range.IsEmpty() || !range.IsBounded()) {
        return Fail(ErrorCode::InvalidInput, "CertifyDeviation: the range must be bounded");
    }
    if (!(std::isfinite(target) && target > 0.0)) {
        return Fail(ErrorCode::InvalidInput,
                    "CertifyDeviation: the target must be finite and positive");
    }
    if (!CoversRange(curve, range) || !CoversRange(pcurve, range)) {
        return Fail(ErrorCode::InvalidInput,
                    "CertifyDeviation: the range leaves a B-spline's domain");
    }
    const std::optional<CurveTrace> c = TraceOf(curve);
    const std::optional<PcurveTrace> p = TraceOf(pcurve);
    const std::optional<SurfaceTrace> s = TraceOf(surface);
    if (!c || !p || !s) {
        return Fail(ErrorCode::NotSupported,
                    "CertifyDeviation: interval evaluation covers lines, circles, ellipses and "
                    "non-rational B-splines on planes and cylinders");
    }
    // Analytic curves on a plane, or on a cylinder along a ruling or a
    // unit-speed turn, have a bound over the whole range in one step; one
    // that exceeds the target goes to the bisection, which can also find a
    // witness against it.
    if (const std::optional<double> bound = ClosedFormBound(curve, surface, pcurve, range);
        bound && *bound <= target) {
        return DeviationCertificate{*bound, 1};
    }
    return Bisect(*c, *p, *s, range, target);
}

} // namespace Stereon
