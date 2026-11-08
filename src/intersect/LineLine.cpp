// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Exact.hpp"

// Two lines meet where they come within the Context's resolution of each
// other, the rule the B-spline pair follows (SplineCurvePair.cpp), so a
// Line and the same segment as a degree-1 B-spline meet alike. Two
// segments' squared distance is convex over the parameter box, so its
// minimum is one point, or a stretch for parallel supports. Every decision,
// whether the minimum lies inside the box and whether it is within the
// resolution, is the sign of a polynomial of degree at most six in the
// input doubles and the resolution, decided exactly; only the placement of
// the meeting is rounded.

namespace Stereon::Detail {
namespace {

// A parameter, an exact quotient with a positive denominator, or an
// infinite end of a domain.
struct LineParameter {
    DynamicExpansion Numerator;
    DynamicExpansion Denominator{1.0};
    int Infinity = 0;
};

template <std::size_t N, std::size_t D>
[[nodiscard]] LineParameter Quotient(const Expansion<N>& n, const Expansion<D>& d)
{
    DynamicExpansion numerator{n};
    DynamicExpansion denominator{d};
    if (denominator.GetSign() == Sign::Negative) {
        return {Negated(numerator), Negated(denominator), 0};
    }
    return {std::move(numerator), std::move(denominator), 0};
}

[[nodiscard]] LineParameter Parameter(const double t)
{
    if (std::isinf(t)) {
        return {DynamicExpansion{}, DynamicExpansion{1.0}, t > 0.0 ? 1 : -1};
    }
    return {DynamicExpansion{t}, DynamicExpansion{1.0}, 0};
}

[[nodiscard]] Sign CompareParameters(const LineParameter& a, const LineParameter& b)
{
    if (a.Infinity != 0 || b.Infinity != 0) {
        return Stereon::SignOf(static_cast<double>(a.Infinity - b.Infinity));
    }
    return Difference(Product(a.Numerator, b.Denominator), Product(b.Numerator, a.Denominator))
        .GetSign();
}

[[nodiscard]] Interval Enclosure(const DynamicExpansion& value)
{
    Interval result{0.0};
    for (const double term : value.GetComponents()) {
        result = result + Interval{term};
    }
    return result;
}

[[nodiscard]] Interval Enclosure(const LineParameter& t)
{
    if (t.Infinity != 0) {
        return Interval{t.Infinity * std::numeric_limits<double>::infinity()};
    }
    return Enclosure(t.Numerator) / Enclosure(t.Denominator);
}

[[nodiscard]] double Estimate(const LineParameter& t)
{
    if (t.Infinity != 0) {
        return t.Infinity * std::numeric_limits<double>::infinity();
    }
    return t.Numerator.GetEstimate() / t.Denominator.GetEstimate();
}

[[nodiscard]] std::array<double, 3> Components(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

[[nodiscard]] std::array<double, 3> Coordinates(const Point3& p)
{
    return {p.X, p.Y, p.Z};
}

[[nodiscard]] ExactVector<1> Exactly(const Point3& p)
{
    return Exactly(Vec3{p.X, p.Y, p.Z});
}

// line.Origin + t line.Direction, exactly.
[[nodiscard]] ExactVector<3> PointAt(const Line& line, const double t)
{
    const auto origin = Coordinates(line.Origin);
    const auto direction = Components(line.Direction);
    ExactVector<3> point;
    for (std::size_t i = 0; i < 3; ++i) {
        point.Hull[i] = Interval{origin[i]} + Interval{direction[i]} * Interval{t};
        point.Value[i] = [o = origin[i], d = direction[i], t] {
            return Sum(Expansion<1>{o}, TwoProduct(d, t));
        };
    }
    return point;
}

// p - q, exactly.
template <std::size_t A, std::size_t B>
[[nodiscard]] ExactVector<A + B> Minus(const ExactVector<A>& p, const ExactVector<B>& q)
{
    ExactVector<A + B> difference;
    for (std::size_t i = 0; i < 3; ++i) {
        difference.Hull[i] = p.Hull[i] - q.Hull[i];
        difference.Value[i] = [pv = p.Value[i], qv = q.Value[i]] { return Difference(pv(), qv()); };
    }
    return difference;
}

// The sign of numerator - bound denominator: how the quotient compares
// with a domain's bound, for a positive denominator.
template <std::size_t N, std::size_t D>
[[nodiscard]] Sign CompareQuotient(const ExactScalar<N>& numerator,
                                   const ExactScalar<D>& denominator, const double bound)
{
    return FilteredSign(numerator.Hull - Interval{bound} * denominator.Hull, [&] {
        return Difference(numerator.Value(), Scaled(denominator.Value(), bound)).GetSign();
    });
}

// Whether numerator / denominator lies in the domain.
template <std::size_t N, std::size_t D>
[[nodiscard]] bool InDomain(const ExactScalar<N>& numerator, const ExactScalar<D>& denominator,
                            const Interval& domain)
{
    const double lower = domain.GetLower();
    const double upper = domain.GetUpper();
    return (std::isinf(lower) ||
            CompareQuotient(numerator, denominator, lower) != Sign::Negative) &&
           (std::isinf(upper) || CompareQuotient(numerator, denominator, upper) != Sign::Positive);
}

[[nodiscard]] std::vector<double> FiniteEnds(const Interval& domain)
{
    std::vector<double> ends;
    for (const double end : {domain.GetLower(), domain.GetUpper()}) {
        if (std::isfinite(end) && std::ranges::find(ends, end) == ends.end()) {
            ends.push_back(end);
        }
    }
    return ends;
}

[[nodiscard]] bool HasFiniteParameter(const Interval& domain)
{
    const double infinity = std::numeric_limits<double>::infinity();
    return domain.GetLower() != infinity && domain.GetUpper() != -infinity;
}

void IncludeBounds(ExponentRange& range, const Interval& da, const Interval& db)
{
    for (const double bound : {da.GetLower(), da.GetUpper(), db.GetLower(), db.GetUpper()}) {
        if (std::isfinite(bound)) {
            range.Include(bound);
        }
    }
}

// Reject unsafe expansion arithmetic before constructing any polynomial.
// Scalable inputs are explicitly unsupported in this handler. The finite
// domain bounds and the resolution enter the decisions, so they are in the
// range too; the decisions have degree at most six.
[[nodiscard]] Result<void> CheckLines(const Line& a, const Line& b, const Interval& da,
                                      const Interval& db, const double resolution)
{
    ExponentRange range;
    for (const auto values : {Coordinates(a.Origin), Components(a.Direction), Coordinates(b.Origin),
                              Components(b.Direction)}) {
        for (const double value : values) {
            if (!std::isfinite(value)) {
                return Fail(ErrorCode::InvalidInput, "LineLine: non-finite line data");
            }
            range.Include(value);
        }
    }
    if (a.Direction == Vec3{} || b.Direction == Vec3{}) {
        return Fail(ErrorCode::InvalidInput, "LineLine: zero direction");
    }
    if (!HasFiniteParameter(da) || !HasFiniteParameter(db)) {
        return Fail(ErrorCode::InvalidInput, "LineLine: domain contains no finite parameter");
    }
    IncludeBounds(range, da, db);
    range.Include(resolution);
    if (CheckDomain(range, 6).Status != DomainStatus::Inside) {
        return Fail(ErrorCode::NotSupported, "LineLine: outside the degree-six exact domain");
    }
    return {};
}

// Bound the returned point's distance from the exact point of A at the
// meeting. L1 is a conservative Euclidean bound; outward interval sums
// include both parameter and coordinate rounding.
[[nodiscard]] double PositionError(const Line& line, const LineParameter& t, const Point3& p)
{
    const auto origin = Coordinates(line.Origin);
    const auto direction = Components(line.Direction);
    const auto position = Coordinates(p);
    const Interval parameter = Enclosure(t);
    Interval error{0.0};
    for (std::size_t i = 0; i < 3; ++i) {
        const Interval delta =
            Interval{origin[i]} + Interval{direction[i]} * parameter - Interval{position[i]};
        error = error + Interval{std::max(std::abs(delta.GetLower()), std::abs(delta.GetUpper()))};
    }
    return error.GetUpper();
}

// The meeting at ta on A and tb on B, placed on A: B's point lies within
// the resolution of A's, which the decision proved.
[[nodiscard]] Result<CurveCurveIntersection> PointResult(const Line& a, const LineParameter& ta,
                                                         const LineParameter& tb,
                                                         const Interval& da, const Interval& db,
                                                         const Context& ctx)
{
    const double t = std::clamp(Estimate(ta), da.GetLower(), da.GetUpper());
    const double s = std::clamp(Estimate(tb), db.GetLower(), db.GetUpper());
    if (!std::isfinite(t) || !std::isfinite(s)) {
        return Fail(ErrorCode::CertificationFailed, "LineLine: unrepresentable contact parameter");
    }
    const Point3 p = Eval(a, t);
    if (!(PositionError(a, ta, p) <= ctx.GetResolution() &&
          PositionError(a, Parameter(t), p) <= ctx.GetResolution())) {
        return Fail(ErrorCode::CertificationFailed,
                    "LineLine: contact rounding exceeds resolution");
    }
    CurveCurveIntersection result;
    result.Points.push_back({t, s, p});
    return result;
}

// Where p lies within the resolution of the line over the domain, the
// parameter of its projection on the line, which then lies in the domain.
template <std::size_t N>
[[nodiscard]] std::optional<LineParameter> Projected(const ExactVector<N>& p, const Line& line,
                                                     const Interval& domain,
                                                     const double resolution)
{
    const auto offset = Minus(p, Exactly(line.Origin));
    const auto direction = Exactly(line.Direction);
    const auto numerator = ExactDot(offset, direction);
    const auto denominator = ExactDot(direction, direction);
    if (!InDomain(numerator, denominator, domain)) {
        return std::nullopt;
    }
    const auto across = ExactCross(offset, direction);
    if (Compare(ExactDot(across, across), Exactly(resolution), denominator) == Sign::Positive) {
        return std::nullopt;
    }
    return Quotient(numerator.Value(), denominator.Value());
}

// Whether two points lie within the resolution of each other.
[[nodiscard]] bool Close(const ExactVector<3>& p, const ExactVector<3>& q, const double resolution)
{
    const auto offset = Minus(p, q);
    return Compare(ExactDot(offset, offset), Exactly(resolution), Exactly(1.0)) != Sign::Positive;
}

struct LineEnd {
    LineParameter A;
    LineParameter B;
};

// The closest approach on the domains' boundary: an end of one segment
// within the resolution of the other's interior, or of the other's end.
// The minimum lies on the boundary when the supports' closest approach
// lies outside the parameter box, so the segments meet exactly when one
// of these does; the nearest is placed.
[[nodiscard]] Result<CurveCurveIntersection>
Boundary(const Line& a, const Line& b, const Interval& da, const Interval& db, const Context& ctx)
{
    const double resolution = ctx.GetResolution();
    std::vector<LineEnd> found;
    const std::vector<double> endsA = FiniteEnds(da);
    const std::vector<double> endsB = FiniteEnds(db);
    for (const double e : endsA) {
        if (auto s = Projected(PointAt(a, e), b, db, resolution)) {
            found.push_back({Parameter(e), std::move(*s)});
        }
        for (const double f : endsB) {
            if (Close(PointAt(a, e), PointAt(b, f), resolution)) {
                found.push_back({Parameter(e), Parameter(f)});
            }
        }
    }
    for (const double f : endsB) {
        if (auto t = Projected(PointAt(b, f), a, da, resolution)) {
            found.push_back({std::move(*t), Parameter(f)});
        }
    }
    if (found.empty()) {
        return CurveCurveIntersection{};
    }
    const auto gap = [&](const LineEnd& end) {
        return Distance(Eval(a, Estimate(end.A)), Eval(b, Estimate(end.B)));
    };
    const auto nearest = std::ranges::min_element(
        found, [&](const LineEnd& p, const LineEnd& q) { return gap(p) < gap(q); });
    return PointResult(a, nearest->A, nearest->B, da, db, ctx);
}

// A source parameter's projection on a parallel target.
[[nodiscard]] LineParameter Project(const Line& source, const double t, const Line& target,
                                    const bool same)
{
    if (std::isinf(t)) {
        return Parameter(same ? t : -t);
    }
    const auto direction = Exactly(target.Direction);
    return Quotient(ExactDot(Minus(PointAt(source, t), Exactly(target.Origin)), direction).Value(),
                    ExactDot(direction, direction).Value());
}

// Parallel supports within the resolution of each other: an overlap over
// the stretch both domains cover, a point where it is one parameter, and
// otherwise their ends' closest approach across the gap.
[[nodiscard]] Result<CurveCurveIntersection>
Parallel(const Line& a, const Line& b, const Interval& da, const Interval& db, const Context& ctx)
{
    const bool same =
        SignOf(ExactDot(Exactly(a.Direction), Exactly(b.Direction))) == Sign::Positive;
    const auto fromA = [&](const double t) -> LineEnd {
        return {Parameter(t), Project(a, t, b, same)};
    };
    const auto fromB = [&](const double s) -> LineEnd {
        return {Project(b, s, a, same), Parameter(s)};
    };
    LineEnd lo = fromA(da.GetLower());
    LineEnd hi = fromA(da.GetUpper());
    LineEnd bl = fromB(same ? db.GetLower() : db.GetUpper());
    LineEnd bh = fromB(same ? db.GetUpper() : db.GetLower());
    if (CompareParameters(lo.A, bl.A) == Sign::Negative) {
        lo = std::move(bl);
    }
    if (CompareParameters(hi.A, bh.A) == Sign::Positive) {
        hi = std::move(bh);
    }
    const Sign order = CompareParameters(lo.A, hi.A);
    if (order == Sign::Positive) {
        return Boundary(a, b, da, db, ctx);
    }
    if (order == Sign::Zero) {
        return PointResult(a, lo.A, lo.B, da, db, ctx);
    }
    // The stretch was decided exactly above. Outward rounding only encloses
    // the endpoint parameters and cannot create a false overlap.
    const Interval onA =
        Intersection(da, Interval{Enclosure(lo.A).GetLower(), Enclosure(hi.A).GetUpper()});
    const Interval onB = Intersection(db, Interval{Enclosure(same ? lo.B : hi.B).GetLower(),
                                                   Enclosure(same ? hi.B : lo.B).GetUpper()});
    CurveCurveIntersection result;
    result.Overlaps.push_back({onA, onB, same});
    return result;
}

// Supports that are not parallel: the closest approach of the supports,
// where it lies in both domains, is the segments' minimum; otherwise the
// minimum lies on the boundary.
[[nodiscard]] Result<CurveCurveIntersection> Nonparallel(const Line& a, const Line& b,
                                                         const ExactVector<4>& normal,
                                                         const Interval& da, const Interval& db,
                                                         const Context& ctx)
{
    const auto w = ExactOffset(b.Origin, a.Origin);
    const auto squared = ExactDot(normal, normal);
    const auto onA = ExactDot(ExactCross(w, Exactly(b.Direction)), normal);
    const auto onB = ExactDot(ExactCross(w, Exactly(a.Direction)), normal);
    if (!InDomain(onA, squared, da) || !InDomain(onB, squared, db)) {
        return Boundary(a, b, da, db, ctx);
    }
    // The supports' distance is |w . n| / |n|.
    if (CompareSquared(ExactDot(w, normal), Exactly(ctx.GetResolution()), squared) ==
        Sign::Positive) {
        return CurveCurveIntersection{};
    }
    return PointResult(a, Quotient(onA.Value(), squared.Value()),
                       Quotient(onB.Value(), squared.Value()), da, db, ctx);
}

} // namespace

Result<CurveCurveIntersection> LineLine(const Curve& a, const Interval& domainA, const Curve& b,
                                        const Interval& domainB, const Context& ctx)
{
    const Line* la = a.Target<Line>();
    const Line* lb = b.Target<Line>();
    if (la == nullptr || lb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineLine: dispatched for other curve types");
    }
    STEREON_TRY_VOID(CheckLines(*la, *lb, domainA, domainB, ctx.GetResolution()));
    const auto normal = ExactCross(Exactly(la->Direction), Exactly(lb->Direction));
    if (!IsExactlyZero(normal)) {
        return Nonparallel(*la, *lb, normal, domainA, domainB, ctx);
    }
    // Parallel supports: their distance is |w x a| / |a|. ADR-0029: only
    // exactly parallel supports overlap; nearly parallel ones meet at a point.
    const auto across = ExactCross(ExactOffset(lb->Origin, la->Origin), Exactly(la->Direction));
    if (Compare(ExactDot(across, across), Exactly(ctx.GetResolution()),
                ExactDot(Exactly(la->Direction), Exactly(la->Direction))) == Sign::Positive) {
        return CurveCurveIntersection{};
    }
    return Parallel(*la, *lb, domainA, domainB, ctx);
}

} // namespace Stereon::Detail
