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
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Planar.hpp"

/// A line against a plane, a sphere or a cylinder: v0.2's first
/// curve–surface handlers, ADR-0010's analytic table carried to the
/// line–quadric pairs. Every configuration decision is the sign of a
/// polynomial in the input doubles (RS-0001 G5), written once over the
/// scalar field of Planar.hpp and read as an enclosure first and exactly
/// only when the enclosure straddles zero: whether the line runs along the
/// plane or the axis, whether it then lies in the surface or apart from it,
/// and the discriminant of the quadratic it cuts out of a quadric, which
/// settles how many points there are and whether one is a touch. The
/// crossing of a plane has a rational parameter, so whether it lies in the
/// requested domains is decided exactly too, as the line–line handler
/// decides; a quadric's roots are irrational and are placed by their
/// rounded values, like the plane curve handlers' (PlanarLines.cpp).

namespace Stereon::Detail {

namespace {

/// Whether the rational num / den, den of sign `denSign`, lies in `range`:
/// the sign of num - b den, turned by den's, against each finite bound b,
/// decided exactly.
template <typename Num, typename Den>
[[nodiscard]] bool RationalInRange(const Num& num, const Den& den, const Sign denSign,
                                   const Interval& range)
{
    const auto beyond = [&](const double bound) {
        const Sign sign = SignOfFormula([&](auto f) {
            using F = decltype(f);
            return F::Sub(num(f), F::Mul(F::Of(bound), den(f)));
        });
        return denSign == Sign::Positive ? sign : -sign;
    };
    if (std::isfinite(range.GetLower()) && beyond(range.GetLower()) == Sign::Negative) {
        return false;
    }
    return !(std::isfinite(range.GetUpper()) && beyond(range.GetUpper()) == Sign::Positive);
}

/// Keeps a point of a quadric when its parameter lies in the curve's domain
/// and its surface parameters in the box: the angle u taken modulo 2 pi and
/// placed in the box's range, as a plane curve handler places a circle's,
/// the height or latitude v as it is.
void KeepPoint(CurveSurfaceIntersection& result, CurveSurfacePoint point, const Interval& domain,
               const Box2& box)
{
    if (!Contains(domain, point.Parameter) || !Contains(box.Y, point.OnSurface.Y)) {
        return;
    }
    const std::optional<double> u = InPeriodicDomain(point.OnSurface.X, box.X);
    if (!u) {
        return;
    }
    point.OnSurface.X = *u;
    result.Points.push_back(point);
}

void Order(CurveSurfaceIntersection& result)
{
    std::ranges::sort(result.Points, [](const CurveSurfacePoint& p, const CurveSurfacePoint& q) {
        return p.Parameter < q.Parameter;
    });
}

/// A root placed again from the line re-centred at its rounded point. The
/// coefficients of the line's quadratic carry the rounding of the line's
/// origin's distance from the surface, |w|² eps in C, which places a root
/// off the surface by about eps |w|² / (r sin), hundreds of units in the
/// last place of the point when the origin is far. Re-centred at the
/// point, the same coefficients are formed from a small offset, so the
/// re-centred quadratic's root nearest zero, s, is placed to the point's
/// own rounding, and t + s is the root to that rounding. For a crossing, a
/// step as large as half the distance to the other root is refused: it
/// would not be the same root.
template <typename CoefficientsFrom>
[[nodiscard]] double Polished(const Line& line, const CoefficientsFrom& coefficientsFrom,
                              const double t, const ContactType contact, const double other)
{
    const auto [a, b, c] = coefficientsFrom(Plain{}, Eval(line, t));
    const LineQuadratic q{a, b, c};
    if (!(q.A > 0.0)) {
        return t;
    }
    double s = -q.B / q.A;
    if (contact == ContactType::Crossing) {
        const std::array<double, 2> roots = RootsOf(q, std::max(0.0, q.B * q.B - q.A * q.C));
        s = std::abs(roots[0]) <= std::abs(roots[1]) ? roots[0] : roots[1];
        if (!(std::abs(s) < 0.5 * std::abs(other - t))) {
            return t;
        }
    }
    const double polished = t + s;
    return std::isfinite(polished) ? polished : t;
}

/// The quadratic A t² + 2 B t + C = 0 that a line cuts out of a quadric,
/// its coefficients over a scalar field for a line through a given origin,
/// decided and solved: the sign of the discriminant B² - A C exactly, from
/// the line's own origin, which settles the count and a touch; the roots
/// from the coefficients and the discriminant to a few roundings
/// (QuadraticOf), which keeps two crossings apart however far the origin,
/// polished from the line re-centred at each (Polished), and
/// placed by KeepPoint with the surface parameters `parameters` gives for a
/// point. A root that doubles cannot place, which happens only when the
/// leading coefficient is exactly positive yet rounds to zero, fails rather
/// than dropping the point (ADR-0029: a result with no component is a
/// proved separation).
template <typename CoefficientsFrom, typename Parameters>
[[nodiscard]] Result<void> CutQuadric(CurveSurfaceIntersection& result, const Line& line,
                                      const Interval& domain, const Box2& box,
                                      const CoefficientsFrom& coefficientsFrom,
                                      const Parameters& parameters, const std::string_view handler)
{
    const Sign sign = SignOfFormula([&](auto f) {
        using F = decltype(f);
        const auto [a, b, c] = coefficientsFrom(f, line.Origin);
        return F::Sub(F::Mul(b, b), F::Mul(a, c));
    });
    if (sign == Sign::Negative) {
        return {};
    }
    const auto [q, discriminant] =
        QuadraticOf([&](auto f) { return coefficientsFrom(f, line.Origin); });
    const auto keep = [&](const double rounded, const ContactType contact,
                          const double other) -> Result<void> {
        if (!std::isfinite(rounded)) {
            return Fail(ErrorCode::IntersectionFailed,
                        std::string(handler) + ": a root could not be placed in double precision");
        }
        const double t = Polished(line, coefficientsFrom, rounded, contact, other);
        const Point3 p = Eval(line, t);
        KeepPoint(result, {t, parameters(p), p, contact, true}, domain, box);
        return {};
    };
    if (sign == Sign::Zero) {
        const double t = -q.B / q.A;
        return keep(t, ContactType::Tangential, t);
    }
    const std::array<double, 2> roots = RootsOf(q, discriminant);
    STEREON_TRY_VOID(keep(roots[0], ContactType::Crossing, roots[1]));
    return keep(roots[1], ContactType::Crossing, roots[0]);
}

} // namespace

// --- Line against plane ---------------------------------------------------

Result<CurveSurfaceIntersection> LinePlane(const Curve& a, const Interval& domain, const Surface& b,
                                           const Box2& box, const Context& /*ctx*/)
{
    const Line* line = a.Target<Line>();
    const Plane* plane = b.Target<Plane>();
    if (line == nullptr || plane == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LinePlane: dispatched for other types");
    }
    CurveSurfaceIntersection result;
    // With n = X x Y, the line O + t d meets the plane through P where
    // n . (O + t d - P) = 0: t = n . (P - O) / n . d.
    const auto normal = [&](auto f) {
        return Cross3(f, Of3(f, plane->XAxis), Of3(f, plane->YAxis));
    };
    const auto den = [&](auto f) { return Dot3(f, normal(f), Of3(f, line->Direction)); };
    const auto num = [&](auto f) {
        return Dot3(f, normal(f), Diff3(f, plane->Origin, line->Origin));
    };
    const Sign slope = SignOfFormula(den);
    if (slope != Sign::Zero) {
        if (!RationalInRange(num, den, slope, domain)) {
            return result;
        }
        // The point's parameters on the plane, X . (O + t d - P) and the
        // same with Y, as rationals over the same denominator.
        const auto along = [&](const Vec3& axis) {
            return [&, axis](auto f) {
                using F = decltype(f);
                const Triple<F> x = Of3(f, axis);
                return F::Add(F::Mul(Dot3(f, x, Diff3(f, line->Origin, plane->Origin)), den(f)),
                              F::Mul(num(f), Dot3(f, x, Of3(f, line->Direction))));
            };
        };
        if (!RationalInRange(along(plane->XAxis), den, slope, box.X) ||
            !RationalInRange(along(plane->YAxis), den, slope, box.Y)) {
            return result;
        }
        const Vec3 n = Cross(plane->XAxis, plane->YAxis);
        const double t = Dot(n, plane->Origin - line->Origin) / Dot(n, line->Direction);
        if (!std::isfinite(t)) {
            return Fail(ErrorCode::IntersectionFailed,
                        "LinePlane: the crossing could not be placed in double precision");
        }
        const Point3 p = Eval(*line, t);
        result.Points.push_back({t, Project(*plane, p), p, ContactType::Crossing, true});
        return result;
    }
    // Along the plane: in it when the line's origin is at height zero,
    // decided exactly; apart otherwise, however close.
    if (SignOfFormula(num) != Sign::Zero) {
        return result;
    }
    const Line2d onPlane{Project(*plane, line->Origin), Vec2{Dot(line->Direction, plane->XAxis),
                                                             Dot(line->Direction, plane->YAxis)}};
    const Interval range = Intersection(domain, ClipLine(onPlane, box));
    if (!range.IsEmpty()) {
        result.Overlaps.push_back({range, Curve2d(onPlane)});
    }
    return result;
}

// --- Line against sphere --------------------------------------------------

Result<CurveSurfaceIntersection> LineSphere(const Curve& a, const Interval& domain,
                                            const Surface& b, const Box2& box,
                                            const Context& /*ctx*/)
{
    const Line* line = a.Target<Line>();
    const Sphere* sphere = b.Target<Sphere>();
    if (line == nullptr || sphere == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineSphere: dispatched for other types");
    }
    CurveSurfaceIntersection result;
    // |O + t d - C|² = r²: A = d . d, B = d . (O - C), C = |O - C|² - r²,
    // for the line through `origin` with the line's direction.
    const auto coefficientsFrom = [&](auto f, const Point3& origin) {
        using F = decltype(f);
        const Triple<F> d = Of3(f, line->Direction);
        const Triple<F> w = Diff3(f, origin, sphere->Center);
        const auto r = F::Of(sphere->Radius);
        return std::array{Dot3(f, d, d), Dot3(f, d, w), F::Sub(Dot3(f, w, w), F::Mul(r, r))};
    };
    STEREON_TRY_VOID(CutQuadric(
        result, *line, domain, box, coefficientsFrom,
        [&](const Point3& p) { return Project(*sphere, p); }, "LineSphere"));
    Order(result);
    return result;
}

// --- Line against cylinder ------------------------------------------------

Result<CurveSurfaceIntersection> LineCylinder(const Curve& a, const Interval& domain,
                                              const Surface& b, const Box2& box,
                                              const Context& /*ctx*/)
{
    const Line* line = a.Target<Line>();
    const Cylinder* cylinder = b.Target<Cylinder>();
    if (line == nullptr || cylinder == nullptr) {
        return Fail(ErrorCode::InvalidInput, "LineCylinder: dispatched for other types");
    }
    CurveSurfaceIntersection result;
    // Across the axis n = X x Y (a unit for a valid cylinder, but taken as
    // it is): the parts of d and of w = O - Q across it are d x n and w x n,
    // scaled by |n|, and the line meets the cylinder where |w x n + t d x n|
    // = r |n|: A = |d x n|², B = (d x n) . (w x n), C = |w x n|² - r² |n|².
    // The cross products keep the cancellation out of the doubles that
    // place the roots: (d . d)(n . n) - (d . n)² would lose them for a line
    // nearly along the axis.
    const auto coefficientsFrom = [&](auto f, const Point3& origin) {
        using F = decltype(f);
        const Triple<F> n = Cross3(f, Of3(f, cylinder->XAxis), Of3(f, cylinder->YAxis));
        const Triple<F> dn = Cross3(f, Of3(f, line->Direction), n);
        const Triple<F> wn = Cross3(f, Diff3(f, origin, cylinder->Origin), n);
        const auto r = F::Of(cylinder->Radius);
        return std::array{Dot3(f, dn, dn), Dot3(f, dn, wn),
                          F::Sub(Dot3(f, wn, wn), F::Mul(F::Mul(r, r), Dot3(f, n, n)))};
    };
    const auto coefficients = [&](auto f) { return coefficientsFrom(f, line->Origin); };
    // Along the axis when A = |d x n|² is zero, decided exactly: then on the
    // cylinder when C is zero, the line's distance from the axis being r,
    // and apart otherwise, inside or outside.
    const Sign across = SignOfFormula([&](auto f) { return coefficients(f)[0]; });
    if (across == Sign::Zero) {
        if (SignOfFormula([&](auto f) { return coefficients(f)[2]; }) != Sign::Zero) {
            return result;
        }
        // The ruling: u constant, v = v0 + t (d . axis), the angle placed in
        // the box's range as a point's would be.
        Point2 start = Project(*cylinder, line->Origin);
        const std::optional<double> u = InPeriodicDomain(start.X, box.X);
        if (!u) {
            return result;
        }
        start.X = *u;
        const Line2d onCylinder{start, Vec2{0.0, Dot(line->Direction, Axis(*cylinder))}};
        const Interval range = Intersection(domain, ClipLine(onCylinder, box));
        if (!range.IsEmpty()) {
            result.Overlaps.push_back({range, Curve2d(onCylinder)});
        }
        return result;
    }
    STEREON_TRY_VOID(CutQuadric(
        result, *line, domain, box, coefficientsFrom,
        [&](const Point3& p) { return Project(*cylinder, p); }, "LineCylinder"));
    Order(result);
    return result;
}

} // namespace Stereon::Detail
