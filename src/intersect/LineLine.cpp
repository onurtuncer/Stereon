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
#include <optional>
#include <utility>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Exact.hpp"

namespace Stereon::Detail {
namespace {

// Finite parameters are exact quotients, never rounded before clipping.
// Degree at most two in numerator and denominator; comparisons have degree
// at most four. Eight expansion components suffice for each cross product.
struct LineParameter {
    Expansion<8> Numerator;
    Expansion<8> Denominator{1.0};
    int Infinity = 0;
};

template <std::size_t N, std::size_t D>
[[nodiscard]] LineParameter Quotient(const Expansion<N>& n, const Expansion<D>& d)
{
    if (d.GetSign() == Sign::Negative) {
        return {Negated(n), Negated(d), 0};
    }
    return {n, d, 0};
}

[[nodiscard]] LineParameter Parameter(const double t)
{
    if (std::isinf(t)) {
        return {{}, Expansion<8>{1.0}, t > 0.0 ? 1 : -1};
    }
    return {Expansion<8>{t}, Expansion<8>{1.0}, 0};
}

[[nodiscard]] Sign CompareParameters(const LineParameter& a, const LineParameter& b)
{
    if (a.Infinity != 0 || b.Infinity != 0) {
        return Stereon::SignOf(static_cast<double>(a.Infinity - b.Infinity));
    }
    return Difference(Product(a.Numerator, b.Denominator), Product(b.Numerator, a.Denominator))
        .GetSign();
}

[[nodiscard]] bool InDomain(const LineParameter& t, const Interval& domain)
{
    return CompareParameters(t, Parameter(domain.GetLower())) != Sign::Negative &&
           CompareParameters(t, Parameter(domain.GetUpper())) != Sign::Positive;
}

[[nodiscard]] Interval Enclosure(const Expansion<8>& value)
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
// Scalable inputs are explicitly unsupported in this first handler. Including
// finite domain bounds protects the degree-four endpoint comparisons too.
[[nodiscard]] Result<void> CheckLines(const Line& a, const Line& b, const Interval& da,
                                      const Interval& db)
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
    if (CheckDomain(range, 4).Status != DomainStatus::Inside) {
        return Fail(ErrorCode::NotSupported, "LineLine: outside the degree-four exact domain");
    }
    return {};
}

// Bound the returned point's distance from the exact contact on each input.
// L1 is a conservative Euclidean bound; outward interval sums include both
// parameter and coordinate rounding. No model tolerance changes the topology.
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

[[nodiscard]] Result<CurveCurveIntersection>
PointResult(const Line& a, const Line& b, const LineParameter& ta, const LineParameter& tb,
            const Interval& da, const Interval& db, const Context& ctx)
{
    const double t = std::clamp(Estimate(ta), da.GetLower(), da.GetUpper());
    const double s = std::clamp(Estimate(tb), db.GetLower(), db.GetUpper());
    if (!std::isfinite(t) || !std::isfinite(s)) {
        return Fail(ErrorCode::CertificationFailed, "LineLine: unrepresentable contact parameter");
    }
    const Point3 p = Eval(a, t);
    const double bound = std::max(PositionError(a, ta, p), PositionError(b, tb, p));
    // Also cover the caller's evaluation at the rounded parameter on B.
    const double rounded = PositionError(b, Parameter(s), p);
    if (!(bound <= ctx.GetResolution() && rounded <= ctx.GetResolution())) {
        return Fail(ErrorCode::CertificationFailed,
                    "LineLine: contact rounding exceeds resolution");
    }
    CurveCurveIntersection result;
    result.Points.push_back({t, s, p});
    return result;
}

struct LineEnd {
    LineParameter A;
    LineParameter B;
};

// Map a source endpoint onto the target's parameter using one non-zero
// coordinate: (source.origin - target.origin + source.direction * t) / target.direction.
[[nodiscard]] LineParameter MapEnd(const double sourceOrigin, const double sourceDirection,
                                   const double targetOrigin, const double targetDirection,
                                   const double t)
{
    if (std::isinf(t)) {
        const bool same = (sourceDirection > 0.0) == (targetDirection > 0.0);
        return Parameter(same ? t : -t);
    }
    return Quotient(Sum(TwoDiff(sourceOrigin, targetOrigin), TwoProduct(sourceDirection, t)),
                    Expansion<1>{targetDirection});
}

[[nodiscard]] Result<CurveCurveIntersection> OverlapResult(const Line& a, const Line& b,
                                                           const LineEnd& lo, const LineEnd& hi,
                                                           const bool same, const Interval& da,
                                                           const Interval& db, const Context& ctx)
{
    const Sign order = CompareParameters(lo.A, hi.A);
    if (order == Sign::Positive) {
        return CurveCurveIntersection{};
    }
    if (order == Sign::Zero) {
        return PointResult(a, b, lo.A, lo.B, da, db, ctx);
    }
    // Topology was decided exactly above. Outward rounding only encloses
    // the endpoint parameters and cannot create a false overlap.
    const Interval onA =
        Intersection(da, Interval{Enclosure(lo.A).GetLower(), Enclosure(hi.A).GetUpper()});
    const Interval onB = Intersection(db, Interval{Enclosure(same ? lo.B : hi.B).GetLower(),
                                                   Enclosure(same ? hi.B : lo.B).GetUpper()});
    CurveCurveIntersection result;
    result.Overlaps.push_back({onA, onB, same});
    return result;
}

[[nodiscard]] Result<CurveCurveIntersection>
Collinear(const Line& a, const Line& b, const Interval& da, const Interval& db, const Context& ctx)
{
    const auto ad = Components(a.Direction);
    const auto bd = Components(b.Direction);
    const auto ao = Coordinates(a.Origin);
    const auto bo = Coordinates(b.Origin);
    std::size_t axis = 0;
    while (ad[axis] == 0.0) {
        ++axis;
    }
    const bool same = (ad[axis] > 0.0) == (bd[axis] > 0.0);
    const auto fromA = [&](const double t) -> LineEnd {
        return {Parameter(t), MapEnd(ao[axis], ad[axis], bo[axis], bd[axis], t)};
    };
    const auto fromB = [&](const double s) -> LineEnd {
        return {MapEnd(bo[axis], bd[axis], ao[axis], ad[axis], s), Parameter(s)};
    };
    const LineEnd al = fromA(da.GetLower());
    const LineEnd ah = fromA(da.GetUpper());
    const LineEnd bl = fromB(same ? db.GetLower() : db.GetUpper());
    const LineEnd bh = fromB(same ? db.GetUpper() : db.GetLower());
    const LineEnd lo = CompareParameters(al.A, bl.A) == Sign::Negative ? bl : al;
    const LineEnd hi = CompareParameters(ah.A, bh.A) == Sign::Positive ? bh : ah;
    return OverlapResult(a, b, lo, hi, same, da, db, ctx);
}

[[nodiscard]] Result<CurveCurveIntersection> Nonparallel(const Line& a, const Line& b,
                                                         const ExactVector<4>& normal,
                                                         const Interval& da, const Interval& db,
                                                         const Context& ctx)
{
    const auto w = ExactOffset(b.Origin, a.Origin);
    if (SignOf(ExactDot(w, normal)) != Sign::Zero) {
        return CurveCurveIntersection{}; // Skew, however close.
    }
    // Choose the largest exact cross component for well-conditioned division.
    std::size_t axis = 0;
    for (std::size_t i = 1; i < 3; ++i) {
        if (std::abs(normal.Value[i]().GetEstimate()) >
            std::abs(normal.Value[axis]().GetEstimate())) {
            axis = i;
        }
    }
    const auto denominator = normal.Value[axis]();
    const LineParameter ta =
        Quotient(ExactCross(w, Exactly(b.Direction)).Value[axis](), denominator);
    const LineParameter tb =
        Quotient(ExactCross(w, Exactly(a.Direction)).Value[axis](), denominator);
    if (!InDomain(ta, da) || !InDomain(tb, db)) {
        return CurveCurveIntersection{};
    }
    return PointResult(a, b, ta, tb, da, db, ctx);
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
    STEREON_TRY_VOID(CheckLines(*la, *lb, domainA, domainB));
    const auto normal = ExactCross(Exactly(la->Direction), Exactly(lb->Direction));
    if (!IsExactlyZero(normal)) {
        return Nonparallel(*la, *lb, normal, domainA, domainB, ctx);
    }
    if (!IsExactlyZero(ExactCross(ExactOffset(lb->Origin, la->Origin), Exactly(la->Direction)))) {
        return CurveCurveIntersection{}; // Parallel distinct supports.
    }
    return Collinear(*la, *lb, domainA, domainB, ctx);
}

} // namespace Stereon::Detail
