// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the B-spline and conic handlers share to isolate the roots of a
// polynomial in Bernstein form with proof: a B-spline's Bezier segments,
// enclosed in interval arithmetic by Boehm's knot insertion, and the
// subdivision of a polynomial's enclosed Bernstein coefficients until every
// piece is proved apart, a proved crossing, an isolated touch, or left
// unresolved (Lane and Riesenfeld's subdivision, with the convex-hull
// property as the proof of separation).
//
// Included in the purview of the handlers' implementation units, after
// `import :Dispatch;` and `import Stereon.Robust;`, with <algorithm>,
// <cmath>, <cstddef>, <cstdint>, <limits>, <optional>, <string>,
// <string_view>, <utility> and <vector> in their global module fragment, and
// after Analytic.hpp (Number) and <Stereon/Core/Try.hpp>. Everything is in an
// anonymous namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

/// The deepest subdivision of a segment: below 2^-52 of its width a piece's
/// ends are roundings of the segment's parameter.
constexpr int MaxDepth = std::numeric_limits<double>::digits - 1;

/// The pieces one segment's search may visit. Each isolated root or touch
/// costs a chain of at most 2 MaxDepth pieces, so this leaves room for
/// dozens per segment; a search that needs more is a curve running along
/// the other geometry, and fails.
constexpr std::size_t PieceBudget = std::size_t{1} << 12;

// --- Bezier segments, enclosed -------------------------------------------

/// A homogeneous pole (w P, w), enclosed; Z is zero in the plane.
struct HPoint {
    Interval X;
    Interval Y;
    Interval Z;
    Interval W;
};

[[maybe_unused]] [[nodiscard]] inline HPoint Blend(const HPoint& a, const HPoint& b,
                                                   const Interval& alpha)
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

[[maybe_unused]] [[nodiscard]] inline HPoint Homogeneous(const Point3& p, const Interval& w)
{
    return {w * Interval{p.X}, w * Interval{p.Y}, w * Interval{p.Z}, w};
}

[[maybe_unused]] [[nodiscard]] inline HPoint Homogeneous(const Point2& p, const Interval& w)
{
    return {w * Interval{p.X}, w * Interval{p.Y}, Interval{0.0}, w};
}

template <typename Spline>
[[nodiscard]] std::vector<HPoint> HomogeneousPoles(const Spline& curve)
{
    std::vector<HPoint> poles;
    poles.reserve(curve.Poles.size());
    for (std::size_t i = 0; i < curve.Poles.size(); ++i) {
        const Interval w{curve.Weights.empty() ? 1.0 : curve.Weights[i]};
        poles.push_back(Homogeneous(curve.Poles[i], w));
    }
    return poles;
}

/// Boehm's insertion of u, in the domain [u_p, u_n+1], once: each new pole
/// a combination of two old ones, its coefficient enclosed, so the poles
/// enclose those of the exactly refined curve. At the domain's upper end
/// the last span stands in, the formula holding on its closed span.
[[maybe_unused]] [[nodiscard]] inline Result<void> InsertKnot(std::vector<double>& knots,
                                                              std::vector<HPoint>& poles,
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
                        "Intersect(BSplineCurve): a knot span too narrow to refine in interval "
                        "arithmetic");
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
template <typename Spline>
[[nodiscard]] Result<std::vector<Segment>> SegmentsOf(const Spline& curve)
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

// --- Polynomials in Bernstein form ---------------------------------------

/// Bernstein coefficients on [0, 1], enclosed.
using Coefficients = std::vector<Interval>;

/// The binomial coefficient C(n, k), exact for n <= 2 BSplineMaxDegree.
[[maybe_unused]] [[nodiscard]] inline double Binomial(const std::size_t n, const std::size_t k)
{
    double value = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        value = value * static_cast<double>(n + 1 - i) / static_cast<double>(i);
    }
    return value;
}

/// The product of two Bernstein polynomials, of degrees p and q:
/// c_k = sum_{i + j = k} C(p, i) C(q, j) / C(p + q, k) a_i b_j.
[[maybe_unused]] [[nodiscard]] inline Coefficients Product(const Coefficients& a,
                                                           const Coefficients& b)
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

// --- Isolating the roots on a segment -----------------------------------

/// A sign that is proved, or not.
enum class Strict : std::uint8_t { Positive, Negative, Unknown };

[[maybe_unused]] [[nodiscard]] inline Strict StrictSign(const Interval& x)
{
    if (x.GetLower() > 0.0) {
        return Strict::Positive;
    }
    return x.GetUpper() < 0.0 ? Strict::Negative : Strict::Unknown;
}

/// The sign every coefficient strictly has, or Unknown.
[[maybe_unused]] [[nodiscard]] inline Strict SignOfAll(const Coefficients& c)
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
[[maybe_unused]] [[nodiscard]] inline Coefficients Differences(const Coefficients& c)
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

[[maybe_unused]] [[nodiscard]] inline double MidOf(const Piece& piece)
{
    return 0.5 * piece.S0 + 0.5 * piece.S1;
}

/// de Casteljau at the middle, enclosed.
[[maybe_unused]] [[nodiscard]] inline std::pair<Piece, Piece> Halves(const Piece& piece)
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
/// Upper] proved to hold it, its contact, whether it is proved a single
/// simple crossing (Decided), and whether it was isolated at all: a piece
/// left at MaxDepth that is neither monotone nor convex holds a contact of
/// higher order, or a stretch of the curve along the other geometry, and is
/// not (Isolated false), so that a caller may settle it exactly or fail.
struct LocalRoot {
    double S;
    double Lower;
    double Upper;
    ContactType Contact;
    bool Decided;
    bool Isolated = true;
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

[[maybe_unused]] [[nodiscard]] inline KernelError
Unresolved(const Search& search, const Piece& piece, const std::string_view why)
{
    const double t = search.T0 + (search.T1 - search.T0) * MidOf(piece);
    return Fail(ErrorCode::IntersectionFailed,
                std::string(search.Handler) + ": " + std::string(why) + " near t = " + Number(t) +
                    " (ADR-0010: an unresolved piece fails, it is never dropped)")
        .error();
}

/// The halving of a piece: de Casteljau on its enclosed coefficients, or,
/// for a polynomial known exactly, its coefficients recomputed over each
/// half (see the conic handlers), so that they stay as tight as the
/// rounding of their own values.
struct DeCasteljau {
    [[maybe_unused]] [[nodiscard]] std::pair<Piece, Piece> operator()(const Piece& piece) const
    {
        return Halves(piece);
    }
};

/// The single simple root of a monotone piece whose ends have the strict,
/// opposite signs `start` and its negation, bisected to rounding.
template <typename Split = DeCasteljau>
[[nodiscard]] LocalRoot Bisected(Piece piece, const Strict start, const Split& split = {})
{
    for (;;) {
        const double mid = MidOf(piece);
        if (!(piece.S0 < mid && mid < piece.S1)) {
            break;
        }
        std::pair<Piece, Piece> halves = split(piece);
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
template <typename Split = DeCasteljau>
void Monotone(const Piece& piece, Search& search, const Split& split = {})
{
    const Strict start = StrictSign(piece.G.front());
    const Strict end = StrictSign(piece.G.back());
    if (start != Strict::Unknown && end != Strict::Unknown) {
        if (start != end) {
            search.Roots.push_back(Bisected(piece, start, split));
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
template <typename Split = DeCasteljau>
[[nodiscard]] Result<void> Convex(const Piece& piece, const int depth, Search& search,
                                  const Split& split = {})
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
        Monotone(piece, search, split);
        return {};
    }
    // Both end slopes within rounding of zero: the slope, monotone, is too
    // throughout, and the piece is the extremum's to rounding.
    const bool flat = start == Strict::Unknown && end == Strict::Unknown;
    if (flat || depth == MaxDepth) {
        search.Roots.push_back({MidOf(piece), piece.S0, piece.S1, ContactType::Tangential, false});
        return {};
    }
    const std::pair<Piece, Piece> halves = split(piece);
    STEREON_TRY_VOID(Convex(halves.first, depth + 1, search, split));
    return Convex(halves.second, depth + 1, search, split);
}

/// Every piece apart, monotone, convex, or at MaxDepth unisolated: a
/// contact of higher order or the curve along the other geometry, kept as a
/// root that is not Isolated for the caller to settle or fail on.
template <typename Split = DeCasteljau>
[[nodiscard]] Result<void> Isolate(const Piece& piece, const int depth, Search& search,
                                   const Split& split = {})
{
    if (++search.Pieces > PieceBudget) {
        return std::unexpected(Unresolved(search, piece, "the subdivision budget is spent"));
    }
    if (SignOfAll(piece.G) != Strict::Unknown) {
        return {}; // Proved apart: g lies in the hull of its coefficients.
    }
    const Coefficients slope = Differences(piece.G);
    if (SignOfAll(slope) != Strict::Unknown) {
        Monotone(piece, search, split);
        return {};
    }
    if (SignOfAll(Differences(slope)) != Strict::Unknown) {
        return Convex(piece, depth, search, split);
    }
    if (depth == MaxDepth) {
        search.Roots.push_back(
            {MidOf(piece), piece.S0, piece.S1, ContactType::Tangential, false, false});
        return {};
    }
    const std::pair<Piece, Piece> halves = split(piece);
    STEREON_TRY_VOID(Isolate(halves.first, depth + 1, search, split));
    return Isolate(halves.second, depth + 1, search, split);
}

// --- From segments to the curve -----------------------------------------

/// A root on the curve: its parameter, the stretch proved to hold it, its
/// contact, and whether it is decided and isolated (see LocalRoot).
struct Found {
    double T;
    Interval Bracket;
    ContactType Contact;
    bool Decided;
    bool Isolated = true;
};

/// The segment's parameter at s, its ends exactly.
[[maybe_unused]] [[nodiscard]] inline double ParameterAt(const Search& search, const double s)
{
    if (s == 0.0) {
        return search.T0;
    }
    return s == 1.0 ? search.T1 : search.T0 + (search.T1 - search.T0) * s;
}

/// The curve's parameters over [s0, s1], outward-rounded and kept in the
/// segment.
[[maybe_unused]] [[nodiscard]] inline Interval BracketOf(const Search& search, const double s0,
                                                         const double s1)
{
    const Interval width = Interval{search.T1} - Interval{search.T0};
    const double lower = (Interval{search.T0} + width * Interval{s0}).GetLower();
    const double upper = (Interval{search.T0} + width * Interval{s1}).GetUpper();
    return Interval{std::max(lower, search.T0), std::min(upper, search.T1)};
}

template <typename Split = DeCasteljau>
[[nodiscard]] Result<void> RootsOfSegment(const Segment& segment, Coefficients g,
                                          const std::string_view handler, std::vector<Found>& found,
                                          const Split& split = {})
{
    Search search{handler, segment.T0, segment.T1, {}, 0};
    STEREON_TRY_VOID(Isolate(Piece{0.0, 1.0, std::move(g)}, 0, search, split));
    for (const LocalRoot& root : search.Roots) {
        found.push_back({ParameterAt(search, root.S), BracketOf(search, root.Lower, root.Upper),
                         root.Contact, root.Decided, root.Isolated});
    }
    return {};
}

/// Undecided roots whose stretches meet are one: a root at a piece's or a
/// segment's end is seen from both sides. A touch wins over a crossing.
[[maybe_unused]] inline void Merge(std::vector<Found>& found)
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
            last.Isolated = last.Isolated && f.Isolated;
            continue;
        }
        merged.push_back(f);
    }
    found = std::move(merged);
}

/// The root when its stretch reaches the requested domain, placed inside
/// it; a root whose stretch leaves the domain is undecided, never dropped
/// on its rounded parameter alone.
[[maybe_unused]] [[nodiscard]] inline std::optional<Found> InDomain(Found found,
                                                                    const Interval& domain)
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

/// The invariants of a plane B-spline that the enclosures rely on, checked
/// as MakeBSplineCurve does: a valid knot vector, weights absent or
/// positive, and finite poles.
[[maybe_unused]] [[nodiscard]] inline Result<void> CheckPlaneSpline(const BSplineCurve2d& curve,
                                                                    const std::string_view where)
{
    STEREON_TRY_VOID(CheckKnotVector(curve.Degree, curve.Knots, curve.Poles.size(), where));
    STEREON_TRY_VOID(CheckWeights(curve.Weights, curve.Poles.size(), where));
    for (const Point2& pole : curve.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return {};
}

/// A periodic curve is intersected over its own domain: its images a
/// period apart are not enumerated yet.
template <typename Spline>
[[nodiscard]] Result<void> CheckPeriodicDomain(const Spline& spline, const Interval& domain,
                                               const std::string_view handler)
{
    if (Period(spline) > 0.0 && !Contains(Domain(spline), domain)) {
        return Fail(ErrorCode::NotSupported,
                    std::string(handler) +
                        ": a periodic B-spline is intersected within its own domain only");
    }
    return {};
}

} // namespace
} // namespace Stereon::Detail
