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
#include <cstdint>
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

/// A B-spline curve, polynomial or rational, against a plane, a sphere or a
/// cylinder: the surface's implicit equation f = 0 composed with the curve.
/// On each Bezier segment, in homogeneous coordinates (X, W) = (w P, w),
/// W^k f(X / W) is a polynomial g of degree p for the plane and 2p for the
/// quadrics, with the sign of f since every weight is positive. Its
/// Bernstein coefficients are formed in outward-rounded interval
/// arithmetic from the segment's poles, which are themselves enclosed:
/// the segments come from Boehm's knot insertion done in intervals, not
/// from the rounded ToBezierSegments, so that every proof below holds for
/// the curve as given (ADR-0029: a result with no component is a proved
/// separation).
///
/// The segment's parameter is subdivided at halves. A piece is dropped only
/// when g's coefficients are all strictly of one sign, which proves g has no
/// root there (the Bernstein convex hull). A piece whose derivative
/// coefficients are strictly of one sign holds at most one root, a simple
/// one, so a crossing: proved present and bisected to rounding when g's end
/// values have strictly opposite signs, reported as undecided when an end
/// value straddles zero. A piece whose second differences are strictly of
/// one sign holds a single extremum; it is split on around it, and the
/// sub-piece still holding it at MaxDepth is reported as an undecided
/// touch: a tangency, a near miss or two crossings within rounding.
/// Anything else, a contact of higher order, a curve along a quadric, or a
/// search past PieceBudget, fails with IntersectionFailed: an unresolved
/// piece is never dropped (ADR-0010, completeness).
///
/// A curve whose every pole lies exactly in the plane lies in it, an
/// overlap; along a quadric there is no such exact decision yet, and the
/// search fails. The surface parameters of a point are those of its rounded
/// position, kept or dropped against the box as the line handlers do.

namespace Stereon::Detail {

namespace {

/// The deepest subdivision of a segment: below 2^-52 of its width a piece's
/// ends are roundings of the segment's parameter.
constexpr int MaxDepth = std::numeric_limits<double>::digits - 1;

/// The pieces one segment's search may visit. Each isolated root or touch
/// costs a chain of at most 2 MaxDepth pieces, so this leaves room for
/// dozens per segment; a search that needs more is a curve running along
/// the surface, and fails.
constexpr std::size_t PieceBudget = std::size_t{1} << 12;

// --- Bezier segments, enclosed -------------------------------------------

/// A homogeneous pole (w P, w), enclosed.
struct HPoint {
    Interval X;
    Interval Y;
    Interval Z;
    Interval W;
};

[[nodiscard]] HPoint Blend(const HPoint& a, const HPoint& b, const Interval& alpha)
{
    const Interval beta = Interval{1.0} - alpha;
    return {a.X * beta + b.X * alpha, a.Y * beta + b.Y * alpha, a.Z * beta + b.Z * alpha,
            a.W * beta + b.W * alpha};
}

/// A Bezier segment: the stretch [T0, T1] of the curve's parameter, exact
/// knots, and its p + 1 homogeneous poles.
struct Segment {
    double T0;
    double T1;
    std::vector<HPoint> Poles;
};

[[nodiscard]] Result<void> CheckSpline(const BSplineCurve& curve)
{
    constexpr std::string_view where = "Intersect(BSplineCurve, surface)";
    STEREON_TRY_VOID(CheckKnotVector(curve.Degree, curve.Knots, curve.Poles.size(), where));
    STEREON_TRY_VOID(CheckWeights(curve.Weights, curve.Poles.size(), where));
    for (const Point3& pole : curve.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return {};
}

[[nodiscard]] std::vector<HPoint> HomogeneousPoles(const BSplineCurve& curve)
{
    std::vector<HPoint> poles;
    poles.reserve(curve.Poles.size());
    for (std::size_t i = 0; i < curve.Poles.size(); ++i) {
        const Interval w{curve.Weights.empty() ? 1.0 : curve.Weights[i]};
        const Point3& p = curve.Poles[i];
        poles.push_back({w * Interval{p.X}, w * Interval{p.Y}, w * Interval{p.Z}, w});
    }
    return poles;
}

/// Boehm's insertion of u, in the domain [u_p, u_n+1], once: each new pole
/// a combination of two old ones, its coefficient enclosed, so the poles
/// enclose those of the exactly refined curve. At the domain's upper end
/// the last span stands in, the formula holding on its closed span.
[[nodiscard]] Result<void> InsertKnot(std::vector<double>& knots, std::vector<HPoint>& poles,
                                      const std::size_t p, const double u)
{
    const std::size_t n = poles.size() - 1;
    const auto above = static_cast<std::size_t>(std::ranges::upper_bound(knots, u) - knots.begin());
    const std::size_t k = std::min(above - 1, n);
    std::vector<HPoint> refined(poles.begin(),
                                poles.begin() + static_cast<std::ptrdiff_t>(k - p + 1));
    for (std::size_t i = k - p + 1; i <= k; ++i) {
        const Interval width = Interval{knots[i + p]} - Interval{knots[i]};
        if (!(width.GetLower() > 0.0)) {
            return Fail(ErrorCode::NotSupported,
                        "Intersect(BSplineCurve, surface): a knot span too narrow to refine in "
                        "interval arithmetic");
        }
        refined.push_back(
            Blend(poles[i - 1], poles[i], (Interval{u} - Interval{knots[i]}) / width));
    }
    refined.insert(refined.end(), poles.begin() + static_cast<std::ptrdiff_t>(k), poles.end());
    knots.insert(knots.begin() + static_cast<std::ptrdiff_t>(k + 1), u);
    poles = std::move(refined);
    return {};
}

/// The curve's Bezier segments: every distinct knot of the domain, its ends
/// included, raised to multiplicity p, after which the p + 1 poles of each
/// non-empty span are its Bernstein coefficients.
[[nodiscard]] Result<std::vector<Segment>> SegmentsOf(const BSplineCurve& curve)
{
    const auto p = static_cast<std::size_t>(curve.Degree);
    std::vector<double> knots = curve.Knots;
    std::vector<HPoint> poles = HomogeneousPoles(curve);
    std::vector<double> values(knots.begin() + static_cast<std::ptrdiff_t>(p),
                               knots.end() - static_cast<std::ptrdiff_t>(p));
    const auto [first, last] = std::ranges::unique(values);
    values.erase(first, last);
    for (const double u : values) {
        while (std::ranges::count(knots, u) < static_cast<std::ptrdiff_t>(p)) {
            STEREON_TRY_VOID(InsertKnot(knots, poles, p, u));
        }
    }
    const Interval domain = Domain(curve);
    std::vector<Segment> segments;
    for (std::size_t j = p; j < poles.size(); ++j) {
        if (knots[j] < knots[j + 1] && domain.GetLower() <= knots[j] &&
            knots[j + 1] <= domain.GetUpper()) {
            segments.push_back({knots[j],
                                knots[j + 1],
                                {poles.begin() + static_cast<std::ptrdiff_t>(j - p),
                                 poles.begin() + static_cast<std::ptrdiff_t>(j + 1)}});
        }
    }
    return segments;
}

// --- The implicit equations in Bernstein form ----------------------------

/// Bernstein coefficients on [0, 1], enclosed.
using Coefficients = std::vector<Interval>;

/// The binomial coefficient C(n, k), exact for n <= 2 BSplineMaxDegree.
[[nodiscard]] double Binomial(const std::size_t n, const std::size_t k)
{
    double value = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        value = value * static_cast<double>(n + 1 - i) / static_cast<double>(i);
    }
    return value;
}

/// The product of two Bernstein polynomials, of degrees p and q:
/// c_k = sum_{i + j = k} C(p, i) C(q, j) / C(p + q, k) a_i b_j.
[[nodiscard]] Coefficients Product(const Coefficients& a, const Coefficients& b)
{
    const std::size_t p = a.size() - 1;
    const std::size_t q = b.size() - 1;
    Coefficients c(p + q + 1, Interval{0.0});
    for (std::size_t i = 0; i <= p; ++i) {
        for (std::size_t j = 0; j <= q; ++j) {
            const Interval factor =
                Interval{Binomial(p, i) * Binomial(q, j)} / Interval{Binomial(p + q, i + j)};
            c[i + j] = c[i + j] + factor * (a[i] * b[j]);
        }
    }
    return c;
}

/// The coordinates of X - W c over a segment's poles: the curve's offset
/// from c, homogeneous, coordinate by coordinate.
struct Offsets {
    Coefficients X;
    Coefficients Y;
    Coefficients Z;
    Coefficients W;
};

[[nodiscard]] Offsets OffsetsFrom(const std::vector<HPoint>& poles, const Point3& c)
{
    Offsets o;
    for (const HPoint& h : poles) {
        o.X.push_back(h.X - h.W * Interval{c.X});
        o.Y.push_back(h.Y - h.W * Interval{c.Y});
        o.Z.push_back(h.Z - h.W * Interval{c.Z});
        o.W.push_back(h.W);
    }
    return o;
}

/// a x b, enclosed.
struct IVec3 {
    Interval X;
    Interval Y;
    Interval Z;
};

[[nodiscard]] IVec3 CrossOf(const Vec3& a, const Vec3& b)
{
    const IVec3 u{Interval{a.X}, Interval{a.Y}, Interval{a.Z}};
    const IVec3 v{Interval{b.X}, Interval{b.Y}, Interval{b.Z}};
    return {u.Y * v.Z - u.Z * v.Y, u.Z * v.X - u.X * v.Z, u.X * v.Y - u.Y * v.X};
}

/// x² + y² + z² - s w², s enclosing r² (and |n|² for a cylinder).
[[nodiscard]] Coefficients SumOfSquares(const Coefficients& x, const Coefficients& y,
                                        const Coefficients& z, const Coefficients& w,
                                        const Interval& s)
{
    Coefficients g = Product(x, x);
    const Coefficients yy = Product(y, y);
    const Coefficients zz = Product(z, z);
    const Coefficients ww = Product(w, w);
    for (std::size_t k = 0; k < g.size(); ++k) {
        g[k] = g[k] + yy[k] + zz[k] - s * ww[k];
    }
    return g;
}

/// The plane through O with normal n = X x Y: n . (X - W O), degree p.
[[nodiscard]] Coefficients PlaneForm(const std::vector<HPoint>& poles, const Plane& plane)
{
    const Offsets o = OffsetsFrom(poles, plane.Origin);
    const IVec3 n = CrossOf(plane.XAxis, plane.YAxis);
    Coefficients g;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        g.push_back(n.X * o.X[i] + n.Y * o.Y[i] + n.Z * o.Z[i]);
    }
    return g;
}

/// The sphere: |X - W C|² - r² W², degree 2p.
[[nodiscard]] Coefficients SphereForm(const std::vector<HPoint>& poles, const Sphere& sphere)
{
    const Offsets o = OffsetsFrom(poles, sphere.Center);
    const Interval r{sphere.Radius};
    return SumOfSquares(o.X, o.Y, o.Z, o.W, r * r);
}

/// The cylinder about n = X x Y, taken as it is as LineCylinder takes it:
/// |(X - W O) x n|² - r² |n|² W², degree 2p.
[[nodiscard]] Coefficients CylinderForm(const std::vector<HPoint>& poles, const Cylinder& cylinder)
{
    const Offsets o = OffsetsFrom(poles, cylinder.Origin);
    const IVec3 n = CrossOf(cylinder.XAxis, cylinder.YAxis);
    Offsets e;
    for (std::size_t i = 0; i < o.W.size(); ++i) {
        e.X.push_back(o.Y[i] * n.Z - o.Z[i] * n.Y);
        e.Y.push_back(o.Z[i] * n.X - o.X[i] * n.Z);
        e.Z.push_back(o.X[i] * n.Y - o.Y[i] * n.X);
    }
    const Interval r{cylinder.Radius};
    const Interval scale = r * r * (n.X * n.X + n.Y * n.Y + n.Z * n.Z);
    return SumOfSquares(e.X, e.Y, e.Z, o.W, scale);
}

// --- Isolating the roots on a segment -----------------------------------

/// A sign that is proved, or not.
enum class Strict : std::uint8_t { Positive, Negative, Unknown };

[[nodiscard]] Strict StrictSign(const Interval& x)
{
    if (x.GetLower() > 0.0) {
        return Strict::Positive;
    }
    return x.GetUpper() < 0.0 ? Strict::Negative : Strict::Unknown;
}

/// The sign every coefficient strictly has, or Unknown.
[[nodiscard]] Strict SignOfAll(const Coefficients& c)
{
    if (c.empty()) {
        return Strict::Unknown;
    }
    const Strict sign = StrictSign(c.front());
    for (const Interval& x : c) {
        if (StrictSign(x) != sign) {
            return Strict::Unknown;
        }
    }
    return sign;
}

/// The derivative's Bernstein coefficients up to the positive factor of
/// the degree and the piece's width: enough for its sign.
[[nodiscard]] Coefficients Differences(const Coefficients& c)
{
    Coefficients d;
    for (std::size_t k = 0; k + 1 < c.size(); ++k) {
        d.push_back(c[k + 1] - c[k]);
    }
    return d;
}

/// A piece [S0, S1] of the segment's parameter in [0, 1], with g's
/// Bernstein coefficients over it. The ends are dyadic, so halving them is
/// exact.
struct Piece {
    double S0;
    double S1;
    Coefficients G;
};

[[nodiscard]] double MidOf(const Piece& piece)
{
    return 0.5 * piece.S0 + 0.5 * piece.S1;
}

/// de Casteljau at the middle, enclosed.
[[nodiscard]] std::pair<Piece, Piece> Halves(const Piece& piece)
{
    const Interval half{0.5};
    const std::size_t m = piece.G.size();
    Coefficients work = piece.G;
    Coefficients left(m);
    Coefficients right(m);
    left[0] = work[0];
    right[m - 1] = work[m - 1];
    for (std::size_t level = 1; level < m; ++level) {
        for (std::size_t i = 0; i + level < m; ++i) {
            work[i] = (work[i] + work[i + 1]) * half;
        }
        left[level] = work[0];
        right[m - 1 - level] = work[m - 1 - level];
    }
    const double mid = MidOf(piece);
    return {Piece{piece.S0, mid, std::move(left)}, Piece{mid, piece.S1, std::move(right)}};
}

/// A root found on a segment: where it is placed, the stretch [Lower,
/// Upper] proved to hold it, and its contact.
struct LocalRoot {
    double S;
    double Lower;
    double Upper;
    ContactType Contact;
    bool Decided;
};

/// One segment's search: its stretch of the curve's parameter, for
/// messages, its roots in order, and the pieces visited.
struct Search {
    std::string_view Handler;
    double T0;
    double T1;
    std::vector<LocalRoot> Roots;
    std::size_t Pieces = 0;
};

[[nodiscard]] KernelError Unresolved(const Search& search, const Piece& piece,
                                     const std::string_view why)
{
    const double t = search.T0 + (search.T1 - search.T0) * MidOf(piece);
    return Fail(ErrorCode::IntersectionFailed,
                std::string(search.Handler) + ": " + std::string(why) + " near t = " + Number(t) +
                    " (ADR-0010: an unresolved piece fails, it is never dropped)")
        .error();
}

/// The single simple root of a monotone piece whose ends have the strict,
/// opposite signs `start` and its negation, bisected to rounding.
[[nodiscard]] LocalRoot Bisected(Piece piece, const Strict start)
{
    for (;;) {
        const double mid = MidOf(piece);
        if (!(piece.S0 < mid && mid < piece.S1)) {
            break;
        }
        std::pair<Piece, Piece> halves = Halves(piece);
        const Strict atMid = StrictSign(halves.first.G.back());
        if (atMid == Strict::Unknown) {
            return {mid, piece.S0, piece.S1, ContactType::Crossing, true};
        }
        piece = atMid == start ? std::move(halves.second) : std::move(halves.first);
    }
    return {MidOf(piece), piece.S0, piece.S1, ContactType::Crossing, true};
}

/// A piece on which g is strictly monotone: no root, one proved crossing,
/// or, when an end value straddles zero, a crossing at or within rounding
/// of that end, which may lie in the neighbouring piece: undecided, with
/// the whole piece as its stretch so that it merges with the neighbour's.
void Monotone(const Piece& piece, Search& search)
{
    const Strict start = StrictSign(piece.G.front());
    const Strict end = StrictSign(piece.G.back());
    if (start != Strict::Unknown && end != Strict::Unknown) {
        if (start != end) {
            search.Roots.push_back(Bisected(piece, start));
        }
        return;
    }
    double s = MidOf(piece);
    if (start == Strict::Unknown && end != Strict::Unknown) {
        s = piece.S0;
    } else if (start != Strict::Unknown) {
        s = piece.S1;
    }
    search.Roots.push_back({s, piece.S0, piece.S1, ContactType::Crossing, false});
}

/// A piece of a piece proved strictly convex or concave, which its
/// sub-pieces inherit: the slope is monotone, so a sub-piece whose slope
/// has the same strict sign at both ends is monotone. Only the sub-piece
/// holding the single extremum is split on, and once both its end slopes
/// straddle zero, or at MaxDepth, it is an undecided touch: at most two
/// roots, within rounding of the extremum.
/// Convexity is proved once, where the coefficients are wide: halving them
/// keeps an absolute error near the rounding of the first, while g shrinks
/// towards a touch, so it could not be proved again at depth.
[[nodiscard]] Result<void> Convex(const Piece& piece, const int depth, Search& search)
{
    if (++search.Pieces > PieceBudget) {
        return std::unexpected(Unresolved(search, piece, "the subdivision budget is spent"));
    }
    if (SignOfAll(piece.G) != Strict::Unknown) {
        return {};
    }
    const Coefficients slope = Differences(piece.G);
    const Strict start = StrictSign(slope.front());
    const Strict end = StrictSign(slope.back());
    if (start != Strict::Unknown && start == end) {
        Monotone(piece, search);
        return {};
    }
    // Both end slopes within rounding of zero: the slope, monotone, is too
    // throughout, and the piece is the extremum's to rounding.
    const bool flat = start == Strict::Unknown && end == Strict::Unknown;
    if (flat || depth == MaxDepth) {
        search.Roots.push_back({MidOf(piece), piece.S0, piece.S1, ContactType::Tangential, false});
        return {};
    }
    const std::pair<Piece, Piece> halves = Halves(piece);
    STEREON_TRY_VOID(Convex(halves.first, depth + 1, search));
    return Convex(halves.second, depth + 1, search);
}

[[nodiscard]] Result<void> Isolate(const Piece& piece, const int depth, Search& search)
{
    if (++search.Pieces > PieceBudget) {
        return std::unexpected(Unresolved(search, piece, "the subdivision budget is spent"));
    }
    if (SignOfAll(piece.G) != Strict::Unknown) {
        return {}; // Proved apart: g lies in the hull of its coefficients.
    }
    const Coefficients slope = Differences(piece.G);
    if (SignOfAll(slope) != Strict::Unknown) {
        Monotone(piece, search);
        return {};
    }
    if (SignOfAll(Differences(slope)) != Strict::Unknown) {
        return Convex(piece, depth, search);
    }
    if (depth == MaxDepth) {
        return std::unexpected(Unresolved(search, piece,
                                          "neither a crossing nor a touch could be isolated (a "
                                          "contact of higher order, or the curve along the "
                                          "surface)"));
    }
    const std::pair<Piece, Piece> halves = Halves(piece);
    STEREON_TRY_VOID(Isolate(halves.first, depth + 1, search));
    return Isolate(halves.second, depth + 1, search);
}

// --- From segments to the curve -----------------------------------------

/// A root on the curve: its parameter, the stretch proved to hold it, and
/// its contact.
struct Found {
    double T;
    Interval Bracket;
    ContactType Contact;
    bool Decided;
};

/// The segment's parameter at s, its ends exactly.
[[nodiscard]] double ParameterAt(const Search& search, const double s)
{
    if (s == 0.0) {
        return search.T0;
    }
    return s == 1.0 ? search.T1 : search.T0 + (search.T1 - search.T0) * s;
}

/// The curve's parameters over [s0, s1], outward-rounded and kept in the
/// segment.
[[nodiscard]] Interval BracketOf(const Search& search, const double s0, const double s1)
{
    const Interval width = Interval{search.T1} - Interval{search.T0};
    const double lower = (Interval{search.T0} + width * Interval{s0}).GetLower();
    const double upper = (Interval{search.T0} + width * Interval{s1}).GetUpper();
    return Interval{std::max(lower, search.T0), std::min(upper, search.T1)};
}

[[nodiscard]] Result<void> RootsOfSegment(const Segment& segment, Coefficients g,
                                          const std::string_view handler, std::vector<Found>& found)
{
    Search search{handler, segment.T0, segment.T1, {}, 0};
    STEREON_TRY_VOID(Isolate(Piece{0.0, 1.0, std::move(g)}, 0, search));
    for (const LocalRoot& root : search.Roots) {
        found.push_back({ParameterAt(search, root.S), BracketOf(search, root.Lower, root.Upper),
                         root.Contact, root.Decided});
    }
    return {};
}

/// Undecided roots whose stretches meet are one: a root at a piece's or a
/// segment's end is seen from both sides. A touch wins over a crossing.
void Merge(std::vector<Found>& found)
{
    std::ranges::stable_sort(found, [](const Found& a, const Found& b) {
        return a.Bracket.GetLower() < b.Bracket.GetLower();
    });
    std::vector<Found> merged;
    for (const Found& f : found) {
        if (!merged.empty() && (!f.Decided || !merged.back().Decided) &&
            Intersects(f.Bracket, merged.back().Bracket)) {
            Found& last = merged.back();
            last.Bracket = Hull(last.Bracket, f.Bracket);
            if (f.Contact == ContactType::Tangential) {
                last.Contact = ContactType::Tangential;
                last.T = f.T;
            }
            last.Decided = false;
            continue;
        }
        merged.push_back(f);
    }
    found = std::move(merged);
}

/// The root when its stretch reaches the requested domain, placed inside
/// it; a root whose stretch leaves the domain is undecided, never dropped
/// on its rounded parameter alone.
[[nodiscard]] std::optional<Found> InDomain(Found found, const Interval& domain)
{
    if (!Intersects(found.Bracket, domain)) {
        return std::nullopt;
    }
    if (!Contains(domain, found.T)) {
        found.T = std::clamp(found.T, std::max(domain.GetLower(), found.Bracket.GetLower()),
                             std::min(domain.GetUpper(), found.Bracket.GetUpper()));
        found.Decided = false;
    }
    return found;
}

/// A periodic curve is intersected over its own domain: its images a
/// period apart are not enumerated yet.
[[nodiscard]] Result<void> CheckPeriodicDomain(const BSplineCurve& spline, const Interval& domain,
                                               const std::string_view handler)
{
    if (Period(spline) > 0.0 && !Contains(Domain(spline), domain)) {
        return Fail(ErrorCode::NotSupported,
                    std::string(handler) +
                        ": a periodic B-spline is intersected within its own domain only");
    }
    return {};
}

/// Every root of g over the segments that reach the domain, merged and
/// placed in the domain.
template <typename Form>
[[nodiscard]] Result<std::vector<Found>> RootsOf(const BSplineCurve& spline, const Interval& domain,
                                                 Form&& form, const std::string_view handler)
{
    STEREON_TRY_VOID(CheckPeriodicDomain(spline, domain, handler));
    STEREON_TRY(const std::vector<Segment> segments, SegmentsOf(spline));
    std::vector<Found> found;
    for (const Segment& segment : segments) {
        if (Intersects(Interval{segment.T0, segment.T1}, domain)) {
            STEREON_TRY_VOID(RootsOfSegment(segment, form(segment.Poles), handler, found));
        }
    }
    Merge(found);
    std::vector<Found> kept;
    for (const Found& f : found) {
        if (const std::optional<Found> placed = InDomain(f, domain)) {
            kept.push_back(*placed);
        }
    }
    return kept;
}

/// The points of the roots whose surface parameters lie in the box, the
/// angle u of a periodic surface placed there.
template <typename Parameters>
[[nodiscard]] CurveSurfaceIntersection PointsOf(const BSplineCurve& spline,
                                                const std::vector<Found>& found, const Box2& box,
                                                Parameters&& parameters, const bool periodicU)
{
    CurveSurfaceIntersection result;
    for (const Found& f : found) {
        const Point3 p = Eval(spline, f.T);
        Point2 uv = parameters(p);
        const std::optional<double> u =
            periodicU ? InPeriodicDomain(uv.X, box.X)
                      : (Contains(box.X, uv.X) ? std::optional<double>{uv.X} : std::nullopt);
        if (!u || !Contains(box.Y, uv.Y)) {
            continue;
        }
        uv.X = *u;
        result.Points.push_back({f.T, uv, p, f.Contact, f.Decided});
    }
    return result;
}

// --- A curve in a plane --------------------------------------------------

/// Whether every pole lies exactly in the plane, n . (P - O) = 0 decided on
/// the input doubles; the curve is then in the plane, in the hull of its
/// poles with positive weights.
[[nodiscard]] bool InPlane(const BSplineCurve& spline, const Plane& plane)
{
    return std::ranges::all_of(spline.Poles, [&](const Point3& pole) {
        return SignOfFormula([&](auto f) {
                   return Dot3(f, Cross3(f, Of3(f, plane.XAxis), Of3(f, plane.YAxis)),
                               Diff3(f, pole, plane.Origin));
               }) == Sign::Zero;
    });
}

/// The overlap of a curve in the plane: its pcurve is the spline of the
/// projected poles, a plane's parameters being affine. Clipping it to a box
/// its poles leave is not implemented.
[[nodiscard]] Result<CurveSurfaceIntersection>
AlongPlane(const BSplineCurve& spline, const Plane& plane, const Interval& domain, const Box2& box)
{
    BSplineCurve2d pcurve{spline.Degree, spline.Knots, {}, spline.Weights};
    for (const Point3& pole : spline.Poles) {
        const Point2 q = Project(plane, pole);
        if (!Contains(box.X, q.X) || !Contains(box.Y, q.Y)) {
            return Fail(ErrorCode::NotImplemented,
                        "BSplineCurve x Plane: a curve in the plane that leaves the requested box "
                        "is not clipped yet");
        }
        pcurve.Poles.push_back(q);
    }
    CurveSurfaceIntersection result;
    const Interval range = Intersection(domain, Domain(spline));
    if (!range.IsEmpty()) {
        result.Overlaps.push_back({range, Curve2d(std::move(pcurve))});
    }
    return result;
}

} // namespace

Result<CurveSurfaceIntersection> SplinePlane(const Curve& a, const Interval& domain,
                                             const Surface& b, const Box2& box,
                                             const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Plane* plane = b.Target<Plane>();
    if (spline == nullptr || plane == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplinePlane: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    if (InPlane(*spline, *plane)) {
        return AlongPlane(*spline, *plane, domain, box);
    }
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return PlaneForm(poles, *plane); },
                    "BSplineCurve x Plane"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*plane, p); }, false);
}

Result<CurveSurfaceIntersection> SplineSphere(const Curve& a, const Interval& domain,
                                              const Surface& b, const Box2& box,
                                              const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Sphere* sphere = b.Target<Sphere>();
    if (spline == nullptr || sphere == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineSphere: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found,
                RootsOf(
                    *spline, domain,
                    [&](const std::vector<HPoint>& poles) { return SphereForm(poles, *sphere); },
                    "BSplineCurve x Sphere"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*sphere, p); }, true);
}

Result<CurveSurfaceIntersection> SplineCylinder(const Curve& a, const Interval& domain,
                                                const Surface& b, const Box2& box,
                                                const Context& /*ctx*/)
{
    const BSplineCurve* spline = a.Target<BSplineCurve>();
    const Cylinder* cylinder = b.Target<Cylinder>();
    if (spline == nullptr || cylinder == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SplineCylinder: dispatched for other types");
    }
    STEREON_TRY_VOID(CheckSpline(*spline));
    STEREON_TRY(const std::vector<Found> found, RootsOf(
                                                    *spline, domain,
                                                    [&](const std::vector<HPoint>& poles) {
                                                        return CylinderForm(poles, *cylinder);
                                                    },
                                                    "BSplineCurve x Cylinder"));
    return PointsOf(
        *spline, found, box, [&](const Point3& p) { return Project(*cylinder, p); }, true);
}

} // namespace Stereon::Detail
