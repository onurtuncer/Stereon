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
    const Interval pi = Around(std::numbers::pi);
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

/// The 3D curve: its point at a parameter, its derivative over a span, its
/// first and second derivatives at a parameter (DAt(t, n), n = 1, 2), its
/// third derivative over a span, and the parameters where it may be less
/// smooth than that (a B-spline's knots; none for the closed forms).
struct CurveTrace {
    std::function<IPoint3(double)> At;
    std::function<IVec3(const Interval&)> D1Over;
    std::function<IVec3(double, int)> DAt;
    std::function<IVec3(const Interval&)> D3Over;
    std::vector<double> Breaks;
};

/// The pcurve: its point at a parameter, and its points and derivatives
/// over a span; its first and second derivatives at a parameter, its second
/// and third over a span, and its breaks, as for the curve.
struct PcurveTrace {
    std::function<IPoint2(double)> At;
    std::function<IPoint2(const Interval&)> Over;
    std::function<IVec2(const Interval&)> D1Over;
    std::function<IVec2(double, int)> DAt;
    std::function<IVec2(const Interval&)> D2Over;
    std::function<IVec2(const Interval&)> D3Over;
    std::vector<double> Breaks;
};

/// The surface: its point and partial derivatives at enclosed parameters.
struct SurfaceTrace {
    std::function<IPoint3(const IPoint2&)> At;
    std::function<IVec3(const IPoint2&, int, int)> Deriv;
};

template <typename G>
[[nodiscard]] CurveTrace ClosedFormCurve(const G& g)
{
    return {[g](const double t) { return Eval(g, IntervalScalar{t}); },
            [g](const Interval& span) { return Deriv(g, IntervalScalar{span}, 1); },
            [g](const double t, const int n) { return Deriv(g, IntervalScalar{t}, n); },
            [g](const Interval& span) { return Deriv(g, IntervalScalar{span}, 3); },
            {}};
}

template <typename G>
[[nodiscard]] PcurveTrace ClosedFormPcurve(const G& g)
{
    return {[g](const double t) { return Eval(g, IntervalScalar{t}); },
            [g](const Interval& span) { return Eval(g, IntervalScalar{span}); },
            [g](const Interval& span) { return Deriv(g, IntervalScalar{span}, 1); },
            [g](const double t, const int n) { return Deriv(g, IntervalScalar{t}, n); },
            [g](const Interval& span) { return Deriv(g, IntervalScalar{span}, 2); },
            [g](const Interval& span) { return Deriv(g, IntervalScalar{span}, 3); },
            {}};
}

template <typename G>
[[nodiscard]] SurfaceTrace SurfaceOf(const G& g)
{
    return {[g](const IPoint2& uv) { return Eval(g, uv); },
            [g](const IPoint2& uv, const int du, const int dv) { return Deriv(g, uv, du, dv); }};
}

// --- B-splines over spans ------------------------------------------------

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

/// A non-rational B-spline with its first two derivative curves, all over
/// intervals, and what it encloses over a parameter span. The hull of the
/// poles active over a span holds for any span but is as wide as the knot
/// spans it covers; within one knot span the spline is a polynomial, so a
/// mean-value form from the next derivative, p(m) + p'(span) (t - m),
/// narrows with the span. Both are enclosures, so their meet is one.
template <std::size_t D>
struct SplineEnclosure {
    BasicBSplineCurve<IntervalScalar, D> C0;
    BasicBSplineCurve<IntervalScalar, D> C1;
    std::optional<BasicBSplineCurve<IntervalScalar, D>> C2;
    /// The third derivative curve; none below degree 3, where the third
    /// derivative is zero within a knot span.
    std::optional<BasicBSplineCurve<IntervalScalar, D>> C3;
    double Magnitude0 = 0.0;
    double Magnitude1 = 0.0;
    double Magnitude2 = 0.0;

    /// True when the span lies within one knot span of C0.
    [[nodiscard]] bool WithinOneSpan(const Interval& span) const
    {
        const Interval domain = Domain(C0);
        const double lower = std::clamp(span.GetLower(), domain.GetLower(), domain.GetUpper());
        const double upper = std::clamp(span.GetUpper(), domain.GetLower(), domain.GetUpper());
        return FindSpan(C0.Knots, C0.Degree, C0.Poles.size(), lower) ==
               FindSpan(C0.Knots, C0.Degree, C0.Poles.size(), upper);
    }

    [[nodiscard]] static IntervalScalar Offset(const Interval& span)
    {
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        return IntervalScalar{Interval{span.GetLower() - m, span.GetUpper() - m}};
    }

    /// The derivative over the span.
    [[nodiscard]] BasicVecN<IntervalScalar, D> D1Over(const Interval& span) const
    {
        BasicVecN<IntervalScalar, D> hull = ToVec(HullOver(C1, span));
        if (!C2 || !WithinOneSpan(span)) {
            return hull;
        }
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        const BasicVecN<IntervalScalar, D> refined =
            ToVec(EnclosedPoint(C1, Magnitude1, m)) + ToVec(HullOver(*C2, span)) * Offset(span);
        return Meet(hull, refined);
    }

    /// The n-th derivative at t, n = 1 or 2; zero where the degree is below n.
    [[nodiscard]] BasicVecN<IntervalScalar, D> DAt(const double t, const int n) const
    {
        if (n == 1) {
            return ToVec(EnclosedPoint(C1, Magnitude1, t));
        }
        return C2 ? ToVec(EnclosedPoint(*C2, Magnitude2, t)) : BasicVecN<IntervalScalar, D>{};
    }

    /// The second derivative over the span, narrowed like D1Over.
    [[nodiscard]] BasicVecN<IntervalScalar, D> D2Over(const Interval& span) const
    {
        if (!C2) {
            return {};
        }
        BasicVecN<IntervalScalar, D> hull = ToVec(HullOver(*C2, span));
        if (!C3 || !WithinOneSpan(span)) {
            return hull;
        }
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        return Meet(hull, ToVec(EnclosedPoint(*C2, Magnitude2, m)) +
                              ToVec(HullOver(*C3, span)) * Offset(span));
    }

    /// The third derivative over the span: the hull of its derivative
    /// curve's active poles, constant within a knot span at degree 3.
    [[nodiscard]] BasicVecN<IntervalScalar, D> D3Over(const Interval& span) const
    {
        return C3 ? ToVec(HullOver(*C3, span)) : BasicVecN<IntervalScalar, D>{};
    }

    /// The distinct knots: where the spline may be less smooth.
    [[nodiscard]] std::vector<double> Breaks() const
    {
        std::vector<double> breaks(C0.Knots.begin(), C0.Knots.end());
        breaks.erase(std::ranges::unique(breaks).begin(), breaks.end());
        return breaks;
    }

    /// The points over the span.
    [[nodiscard]] BasicPointN<IntervalScalar, D> Over(const Interval& span) const
    {
        const double m = 0.5 * (span.GetLower() + span.GetUpper());
        const BasicPointN<IntervalScalar, D> refined =
            EnclosedPoint(C0, Magnitude0, m) + D1Over(span) * Offset(span);
        return ToPoint(Meet(ToVec(HullOver(C0, span)), ToVec(refined)));
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
[[nodiscard]] std::optional<std::shared_ptr<const SplineEnclosure<D>>>
EnclosureOf(const BasicBSplineCurve<double, D>& spline)
{
    if (IsRational(spline)) {
        return std::nullopt;
    }
    SplineEnclosure<D> e{Enclosed(spline), {}, std::nullopt, std::nullopt};
    std::optional<BasicBSplineCurve<IntervalScalar, D>> first = DerivativeCurve(e.C0);
    if (!first) {
        return std::nullopt;
    }
    e.C1 = std::move(*first);
    e.C2 = DerivativeCurve(e.C1);
    if (e.C2) {
        e.C3 = DerivativeCurve(*e.C2);
        e.Magnitude2 = PoleMagnitude(*e.C2);
    }
    e.Magnitude0 = PoleMagnitude(e.C0);
    e.Magnitude1 = PoleMagnitude(e.C1);
    return std::make_shared<const SplineEnclosure<D>>(std::move(e));
}

[[nodiscard]] std::optional<CurveTrace> SplineCurve(const BSplineCurve& spline)
{
    const auto e = EnclosureOf(spline);
    if (!e) {
        return std::nullopt;
    }
    const std::shared_ptr<const SplineEnclosure<3>> c = *e;
    return CurveTrace{[c](const double t) { return EnclosedPoint(c->C0, c->Magnitude0, t); },
                      [c](const Interval& span) { return c->D1Over(span); },
                      [c](const double t, const int n) { return c->DAt(t, n); },
                      [c](const Interval& span) { return c->D3Over(span); }, c->Breaks()};
}

[[nodiscard]] std::optional<PcurveTrace> SplinePcurve(const BSplineCurve2d& spline)
{
    const auto e = EnclosureOf(spline);
    if (!e) {
        return std::nullopt;
    }
    const std::shared_ptr<const SplineEnclosure<2>> c = *e;
    return PcurveTrace{[c](const double t) { return EnclosedPoint(c->C0, c->Magnitude0, t); },
                       [c](const Interval& span) { return c->Over(span); },
                       [c](const Interval& span) { return c->D1Over(span); },
                       [c](const double t, const int n) { return c->DAt(t, n); },
                       [c](const Interval& span) { return c->D2Over(span); },
                       [c](const Interval& span) { return c->D3Over(span); },
                       c->Breaks()};
}

// --- The erased operands --------------------------------------------------

[[nodiscard]] std::optional<CurveTrace> TraceOf(const Curve& curve)
{
    if (const Line* line = curve.Target<Line>()) {
        return ClosedFormCurve(Enclosed(*line));
    }
    if (const Circle* circle = curve.Target<Circle>()) {
        return ClosedFormCurve(Enclosed(*circle));
    }
    if (const Ellipse* ellipse = curve.Target<Ellipse>()) {
        return ClosedFormCurve(Enclosed(*ellipse));
    }
    if (const BSplineCurve* spline = curve.Target<BSplineCurve>()) {
        return SplineCurve(*spline);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<PcurveTrace> TraceOf(const Curve2d& pcurve)
{
    if (const Line2d* line = pcurve.Target<Line2d>()) {
        return ClosedFormPcurve(Enclosed(*line));
    }
    if (const Circle2d* circle = pcurve.Target<Circle2d>()) {
        return ClosedFormPcurve(Enclosed(*circle));
    }
    if (const Ellipse2d* ellipse = pcurve.Target<Ellipse2d>()) {
        return ClosedFormPcurve(Enclosed(*ellipse));
    }
    if (const BSplineCurve2d* spline = pcurve.Target<BSplineCurve2d>()) {
        return SplinePcurve(*spline);
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

/// What a span proves about |c(t) - S(p(t))| over it.
struct SpanBounds {
    /// A bound over the whole span: the mean-value form
    /// f(m) + f'(span) (t - m) at the span's midpoint m.
    double Upper;
    /// A deviation that certainly occurs, at the midpoint: the least
    /// magnitude of f(m)'s enclosure. A witness that no certificate below it
    /// exists.
    double Witness;
};

[[nodiscard]] SpanBounds BoundsOver(const CurveTrace& c, const PcurveTrace& p,
                                    const SurfaceTrace& s, const Interval& span)
{
    const double m = 0.5 * (span.GetLower() + span.GetUpper());
    const IVec3 atMidpoint = c.At(m) - s.At(p.At(m));
    const IPoint2 uv = p.Over(span);
    const IVec2 duv = p.D1Over(span);
    const IVec3 slope = c.D1Over(span) - (s.Deriv(uv, 1, 0) * duv.X + s.Deriv(uv, 0, 1) * duv.Y);
    const IntervalScalar offset{Interval{span.GetLower() - m, span.GetUpper() - m}};
    const IVec3 f = atMidpoint + slope * offset;
    return {Length(f).GetUpper(), Length(atMidpoint).GetLower()};
}

/// A bound over the span from the third-order Taylor form at its midpoint
/// m, with s = t - m and |s| <= r, the span's half width:
///     |f(m)| + |f1(m)| r + |f2(m)| r^2 / 2 + sup |f3| r^3 / 6,
/// fk the k-th derivative of f(t) = c(t) - S(p(t)), the remainder by its
/// integral form, so that it holds for vectors. The
/// point terms are evaluated at m, where c and S(p) nearly cancel and their
/// enclosures are tight; only the cubic term is enclosed over the span, so
/// the bound shrinks with the cube of the span where the mean-value form
/// shrinks with its square, which made the bisection rediscover each knot
/// span in some forty pieces (issue #155). Valid where c and S(p) are three
/// times differentiable over the span: within one polynomial piece of each
/// B-spline, which Bisect's seeding at the knots guarantees.
[[nodiscard]] double TaylorBound(const CurveTrace& c, const PcurveTrace& p, const SurfaceTrace& s,
                                 const Interval& span)
{
    const double m = 0.5 * (span.GetLower() + span.GetUpper());
    const IPoint2 uv = p.At(m);
    const IVec2 d1 = p.DAt(m, 1);
    const IVec2 d2 = p.DAt(m, 2);
    const IVec3 s10 = s.Deriv(uv, 1, 0);
    const IVec3 s01 = s.Deriv(uv, 0, 1);
    const IntervalScalar two{2.0};
    const IntervalScalar three{3.0};
    const IVec3 f0 = c.At(m) - s.At(uv);
    const IVec3 f1 = c.DAt(m, 1) - (s10 * d1.X + s01 * d1.Y);
    const IVec3 f2 =
        c.DAt(m, 2) - (s.Deriv(uv, 2, 0) * (d1.X * d1.X) + s.Deriv(uv, 1, 1) * (two * d1.X * d1.Y) +
                       s.Deriv(uv, 0, 2) * (d1.Y * d1.Y) + s10 * d2.X + s01 * d2.Y);
    // The third derivative of S(p(t)) over the span, by the chain rule.
    const IPoint2 box = p.Over(span);
    const IVec2 e1 = p.D1Over(span);
    const IVec2 e2 = p.D2Over(span);
    const IVec2 e3 = p.D3Over(span);
    const IVec3 g3 = s.Deriv(box, 3, 0) * (e1.X * e1.X * e1.X) +
                     s.Deriv(box, 2, 1) * (three * e1.X * e1.X * e1.Y) +
                     s.Deriv(box, 1, 2) * (three * e1.X * e1.Y * e1.Y) +
                     s.Deriv(box, 0, 3) * (e1.Y * e1.Y * e1.Y) +
                     s.Deriv(box, 2, 0) * (three * e1.X * e2.X) +
                     s.Deriv(box, 1, 1) * (three * (e1.X * e2.Y + e2.X * e1.Y)) +
                     s.Deriv(box, 0, 2) * (three * e1.Y * e2.Y) + s.Deriv(box, 1, 0) * e3.X +
                     s.Deriv(box, 0, 1) * e3.Y;
    const IVec3 f3 = c.D3Over(span) - g3;
    const Interval r{0.5 * (span.GetUpper() - span.GetLower())};
    const Interval total = Interval{Length(f0).GetUpper()} + Interval{Length(f1).GetUpper()} * r +
                           Interval{Length(f2).GetUpper()} * r * r / Interval{2.0} +
                           Interval{Length(f3).GetUpper()} * r * r * r / Interval{6.0};
    return total.GetUpper();
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
        SpanBounds local = BoundsOver(c, p, s, span);
        if (local.Upper > target) {
            local.Upper = std::min(local.Upper, TaylorBound(c, p, s, span));
        }
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
