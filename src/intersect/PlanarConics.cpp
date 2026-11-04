// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
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
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Bernstein.hpp"
#include "Planar.hpp"

/// Two conics in the plane, circles, ellipses, hyperbolas and parabolas,
/// against each other (the Boolean prototype's M6), except two circles,
/// which PlanarLines.cpp decides in closed form. One conic, A, is taken in
/// a rational chart: its point is the homogeneous W(u) C + Alpha(u) X +
/// Beta(u) Y for u in [-1, 1], each of W, Alpha and Beta a quadratic in u
/// whose coefficients are products of the input doubles. An ellipse has two
/// charts, u = tan(s/2) about s = 0 and about s = pi; a hyperbola one,
/// u = tanh(t/2), whose ends are its points at infinity; a parabola two,
/// t = L u and t = L / u for a power of two L near 4f, the second's u = 0
/// at infinity. The other conic, B, is taken by its implicit equation in
/// its frame (b² x² ± a² y² = a² b², or y² = 4 f x), so that B's equation
/// along A's chart is a polynomial Q of degree four in u, formed exactly on
/// heap expansions. The open conic is taken in the chart, so that a
/// hyperbola's two branches meet only as B: a point of B's other branch,
/// 2a away from the one parametrised, is dropped by the sign of its
/// abscissa.
///
/// Q's roots at a chart's points at infinity are divided out exactly
/// (synthetic division by u - r, exact for a double r), so that Q is not
/// zero there. Q is then shifted to [0, 1] exactly and its Bernstein
/// coefficients enclosed; Bernstein.hpp's subdivision proves every piece
/// apart, a single simple crossing, an isolated touch, or leaves it
/// unresolved. A root it cannot decide is settled exactly when it lies at
/// a double u*: Q(u*) = 0 decided on expansions, the multiplicity m read
/// from Q's exact Taylor coefficients at u*, and Q / (u - u*)^m proved free
/// of zeros over the root's bracket by interval evaluation, which proves
/// u* the only root there, of multiplicity m: a crossing when m is odd
/// (the curves change sides), a tangency when even, both Decided. Any other
/// undecided root is reported bracketed to rounding, Decided false; a root
/// left unresolved, or one within rounding of a point at infinity, fails
/// with IntersectionFailed. A Q that vanishes identically is a coincidence,
/// whose overlaps are the requested domains' common stretches under the
/// affine correspondence of the two parameters.

namespace Stereon::Detail {

namespace {

// --- Conics --------------------------------------------------------------

enum class Shape : std::uint8_t { Ellipse, Hyperbola, Parabola };

/// A conic in its frame: a circle or an ellipse with semi-axes A and B, a
/// hyperbola with semi-axes A and B, or a parabola with its vertex as
/// Center and its focal length as A.
struct Conic {
    Shape Kind;
    Point2 Center;
    double A;
    double B;
    Vec2 X;
    Vec2 Y;
};

[[nodiscard]] std::optional<Conic> ConicOf(const Curve2d& curve)
{
    if (const Circle2d* c = curve.Target<Circle2d>()) {
        return Conic{Shape::Ellipse, c->Center, c->Radius, c->Radius, c->XAxis, c->YAxis};
    }
    if (const Ellipse2d* e = curve.Target<Ellipse2d>()) {
        return Conic{Shape::Ellipse, e->Center, e->MajorRadius, e->MinorRadius, e->XAxis, e->YAxis};
    }
    if (const Hyperbola2d* h = curve.Target<Hyperbola2d>()) {
        return Conic{Shape::Hyperbola, h->Center, h->MajorRadius,
                     h->MinorRadius,   h->XAxis,  h->YAxis};
    }
    if (const Parabola2d* p = curve.Target<Parabola2d>()) {
        return Conic{Shape::Parabola, p->Vertex, p->FocalLength, 0.0, p->XAxis, p->YAxis};
    }
    return std::nullopt;
}

/// Which conic of a pair is taken in its chart: the one whose implicit
/// equation is the worse, a hyperbola's (two branches) before a
/// parabola's, before an ellipse's.
[[nodiscard]] int ChartRank(const Shape shape)
{
    switch (shape) {
    case Shape::Hyperbola:
        return 0;
    case Shape::Parabola:
        return 1;
    case Shape::Ellipse:
        return 2;
    }
    return 2;
}

/// Whether B, not A, is taken in its chart: the lower ChartRank, and
/// between conics of one rank the first in the order of their data, so
/// that a pair and the pair exchanged are one computation.
[[nodiscard]] bool ChartedSecond(const Conic& a, const Conic& b)
{
    if (ChartRank(a.Kind) != ChartRank(b.Kind)) {
        return ChartRank(b.Kind) < ChartRank(a.Kind);
    }
    const auto data = [](const Conic& c) {
        return std::array{c.A, c.B, c.Center.X, c.Center.Y, c.X.X, c.X.Y, c.Y.X, c.Y.Y};
    };
    return data(b) < data(a);
}

/// The parameter of a point of the conic: the angle of an ellipse, the
/// t = arsinh(y / b) of a hyperbola, the t = y of a parabola, y the
/// point's ordinate in the frame.
[[nodiscard]] double ParameterOn(const Conic& c, const Point2& p)
{
    switch (c.Kind) {
    case Shape::Ellipse:
        return AngleIn(c.Center, c.X, c.Y, c.A, c.B, p);
    case Shape::Hyperbola:
        return std::asinh(Dot(p - c.Center, c.Y) / c.B);
    case Shape::Parabola:
        return Dot(p - c.Center, c.Y);
    }
    return 0.0;
}

// --- Polynomials over a field ---------------------------------------------

/// A polynomial over a field, coefficients in increasing degree.
template <typename F>
using Polynomial = std::vector<typename F::Scalar>;

using ExactPolynomial = Polynomial<Exact>;

template <typename F>
[[nodiscard]] Polynomial<F> Multiply(F, const Polynomial<F>& p, const Polynomial<F>& q)
{
    Polynomial<F> r(p.size() + q.size() - 1, F::Of(0.0));
    for (std::size_t i = 0; i < p.size(); ++i) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            r[i + j] = F::Add(r[i + j], F::Mul(p[i], q[j]));
        }
    }
    return r;
}

template <typename F>
[[nodiscard]] Polynomial<F> Scaled(F, const Polynomial<F>& p, const typename F::Scalar& k)
{
    Polynomial<F> r;
    r.reserve(p.size());
    for (const auto& c : p) {
        r.push_back(F::Mul(c, k));
    }
    return r;
}

template <typename F>
[[nodiscard]] Polynomial<F> Combine(F, const Polynomial<F>& p, const Polynomial<F>& q,
                                    const bool subtract)
{
    Polynomial<F> r(std::max(p.size(), q.size()), F::Of(0.0));
    for (std::size_t i = 0; i < r.size(); ++i) {
        const auto pi = i < p.size() ? p[i] : F::Of(0.0);
        const auto qi = i < q.size() ? q[i] : F::Of(0.0);
        r[i] = subtract ? F::Sub(pi, qi) : F::Add(pi, qi);
    }
    return r;
}

// --- A conic in its chart, against the other's equation ------------------

/// The two charts of an ellipse or a parabola; a hyperbola has the first.
enum class Chart : std::uint8_t { Near, Far };

/// A chart's homogeneous point W C + Alpha X + Beta Y, quadratics in u.
template <typename F>
struct ChartPoint {
    Polynomial<F> W;
    Polynomial<F> Alpha;
    Polynomial<F> Beta;
};

/// An ellipse about s = 0 (h = 1) or s = pi (h = -1): with u = tan(s/2),
/// (1 + u²) (cos s, sin s) = (1 - u², 2u), each times h.
/// A hyperbola with u = tanh(t/2): (1 - u²) (cosh t, sinh t) = (1 + u², 2u).
/// A parabola, times 4f: t = L u gives (4f, L² u², 4f L u), and t = L / u,
/// times u², (4f u², L², 4f L u).
template <typename F>
[[nodiscard]] ChartPoint<F> ChartOf(F, const Conic& c, const Chart chart, const double scale)
{
    const auto zero = F::Of(0.0);
    const auto one = F::Of(1.0);
    const auto a = F::Of(c.A);
    if (c.Kind == Shape::Ellipse) {
        const double h = chart == Chart::Near ? 1.0 : -1.0;
        const auto ha = F::Of(h * c.A);
        return {{one, zero, one}, {ha, zero, F::Sub(zero, ha)}, {zero, F::Of(2.0 * h * c.B), zero}};
    }
    if (c.Kind == Shape::Hyperbola) {
        return {{one, zero, F::Sub(zero, one)}, {a, zero, a}, {zero, F::Of(2.0 * c.B), zero}};
    }
    const auto focal = F::Mul(F::Of(4.0), a);
    const auto ll = F::Of(scale * scale);
    Polynomial<F> beta{zero, F::Mul(focal, F::Of(scale)), zero};
    if (chart == Chart::Near) {
        return {{focal, zero, zero}, {zero, zero, ll}, std::move(beta)};
    }
    return {{zero, zero, focal}, {ll, zero, zero}, std::move(beta)};
}

/// B's implicit equation on the homogeneous coordinates (xi, eta, w) of A's
/// point in B's frame: b² xi² ± a² eta² - a² b² w², or eta² - 4 f xi w.
template <typename F>
[[nodiscard]] Polynomial<F> ImplicitOf(F, const Conic& b, const Polynomial<F>& xi,
                                       const Polynomial<F>& eta, const Polynomial<F>& w)
{
    const Polynomial<F> etaEta = Multiply(F{}, eta, eta);
    if (b.Kind == Shape::Parabola) {
        const auto focal = F::Mul(F::Of(4.0), F::Of(b.A));
        return Combine(F{}, etaEta, Scaled(F{}, Multiply(F{}, xi, w), focal), true);
    }
    const auto aa = F::Mul(F::Of(b.A), F::Of(b.A));
    const auto bb = F::Mul(F::Of(b.B), F::Of(b.B));
    const Polynomial<F> xx = Scaled(F{}, Multiply(F{}, xi, xi), bb);
    const Polynomial<F> yy = Scaled(F{}, etaEta, aa);
    const Polynomial<F> ww = Scaled(F{}, Multiply(F{}, w, w), F::Mul(aa, bb));
    return Combine(F{}, Combine(F{}, xx, yy, b.Kind == Shape::Hyperbola), ww, true);
}

/// Q(u): B's equation along A's chart, of degree four. The coordinates of
/// A's point in B's frame are (H - W c_B) . e for e = X_B, Y_B.
template <typename F>
[[nodiscard]] Polynomial<F> Composite(F, const Conic& a, const ChartPoint<F>& h, const Conic& b)
{
    const Pair2<F> d = Diff2(F{}, a.Center, b.Center);
    const auto coordinate = [&](const Vec2& axis) {
        const Pair2<F> e = Of2(F{}, axis);
        const Polynomial<F> centre = Scaled(F{}, h.W, Dot2(F{}, d, e));
        const Polynomial<F> along = Scaled(F{}, h.Alpha, Dot2(F{}, Of2(F{}, a.X), e));
        const Polynomial<F> across = Scaled(F{}, h.Beta, Dot2(F{}, Of2(F{}, a.Y), e));
        return Combine(F{}, Combine(F{}, centre, along, false), across, false);
    };
    return ImplicitOf(F{}, b, coordinate(b.X), coordinate(b.Y), h.W);
}

// --- Exact polynomials ---------------------------------------------------

/// e, enclosed: the sum of its non-overlapping components, outward-rounded.
[[nodiscard]] Interval Enclose(const DynamicExpansion& e)
{
    Interval sum{0.0};
    for (const double c : e.GetComponents()) {
        sum = sum + Interval{c};
    }
    return sum;
}

[[nodiscard]] DynamicExpansion ValueAt(const ExactPolynomial& p, const double u)
{
    DynamicExpansion value;
    for (std::size_t i = p.size(); i-- > 0;) {
        value = Sum(Scaled(value, u), p[i]);
    }
    return value;
}

/// The quotient of p by u - r, by synthetic division: exact, each step a
/// product by the double r and a sum. The remainder is p(r).
[[nodiscard]] ExactPolynomial Deflated(const ExactPolynomial& p, const double r)
{
    ExactPolynomial q(p.size() - 1);
    DynamicExpansion carry;
    for (std::size_t i = p.size(); i-- > 1;) {
        carry = Sum(Scaled(carry, r), p[i]);
        q[i - 1] = carry;
    }
    return q;
}

[[nodiscard]] bool IsZero(const ExactPolynomial& p)
{
    return std::ranges::all_of(p, [](const DynamicExpansion& c) { return c.IsZero(); });
}

/// The Taylor coefficients c_k of p at the double r, exactly: the
/// remainders of repeated synthetic division by u - r.
[[nodiscard]] std::vector<DynamicExpansion> TaylorAt(ExactPolynomial p, const double r)
{
    std::vector<DynamicExpansion> taylor;
    while (!p.empty()) {
        taylor.push_back(ValueAt(p, r));
        p = Deflated(p, r);
    }
    return taylor;
}

/// The Bernstein coefficients on sigma in [0, 1] of p(u0 + w sigma), for
/// doubles u0 and w, enclosed: the shifted coefficients c_k w^k are exact,
/// from the Taylor coefficients at u0, and only the change of basis
/// b_i = sum_{j <= i} C(i, j) / C(n, j) q_j is outward-rounded. Each piece
/// of the subdivision is formed so, from p itself, rather than by halving
/// the enclosures of its parent: the coefficients then stay as tight as the
/// rounding of their own size, however small the piece, so that a near
/// miss or a contact of high order at a double is resolved rather than
/// lost below the rounding of the whole chart's coefficients.
[[nodiscard]] Coefficients BernsteinOn(const ExactPolynomial& p, const double u0, const double w)
{
    std::vector<DynamicExpansion> q = TaylorAt(p, u0);
    for (std::size_t k = 1; k < q.size(); ++k) {
        for (std::size_t i = 0; i < k; ++i) {
            q[k] = Scaled(q[k], w);
        }
    }
    const std::size_t n = q.size() - 1;
    Coefficients b(n + 1, Interval{0.0});
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            b[i] = b[i] + Interval{Binomial(i, j)} / Interval{Binomial(n, j)} * Enclose(q[j]);
        }
    }
    return b;
}

/// The chart's piece [S0, S1] of s in [0, 1] is u in [2 S0 - 1, 2 S1 - 1],
/// both exact for the dyadic ends of the subdivision.
[[nodiscard]] Piece ExactPiece(const ExactPolynomial& p, const double s0, const double s1)
{
    return {s0, s1, BernsteinOn(p, 2.0 * s0 - 1.0, 2.0 * (s1 - s0))};
}

/// The halving of a piece of the chart, from the exact polynomial.
struct ExactSplit {
    const ExactPolynomial& P;
    [[nodiscard]] std::pair<Piece, Piece> operator()(const Piece& piece) const
    {
        const double mid = MidOf(piece);
        return {ExactPiece(P, piece.S0, mid), ExactPiece(P, mid, piece.S1)};
    }
};

// --- Settling a root exactly ---------------------------------------------

/// The multiplicity of the root u* of p when p / (u - u*)^m has no zero on
/// the bracket, so that u* is the only root there; nothing otherwise. The
/// Taylor coefficients c_k of p at u* are the remainders of repeated
/// synthetic division, exact; the quotient sum_{k >= m} c_k h^(k - m),
/// h = u - u*, is enclosed over the bracket by Horner's rule in intervals.
[[nodiscard]] std::optional<std::size_t> MultiplicityAt(const ExactPolynomial& p, const double root,
                                                        const Interval& bracket)
{
    const std::vector<DynamicExpansion> taylor = TaylorAt(p, root);
    if (!taylor.front().IsZero()) {
        return std::nullopt;
    }
    const auto first = std::ranges::find_if(taylor, [](const auto& c) { return !c.IsZero(); });
    if (first == taylor.end()) {
        return std::nullopt;
    }
    const auto m = static_cast<std::size_t>(first - taylor.begin());
    const Interval h = bracket - Interval{root};
    Interval quotient{0.0};
    for (std::size_t k = taylor.size(); k-- > m;) {
        quotient = quotient * h + Enclose(taylor[k]);
    }
    return Contains(quotient, 0.0) ? std::nullopt : std::optional<std::size_t>{m};
}

/// The doubles where an undecided root may lie exactly: its placement, the
/// ends and the middle of its bracket, and zero.
[[nodiscard]] std::vector<double> Candidates(const Found& found)
{
    const double lower = found.Bracket.GetLower();
    const double upper = found.Bracket.GetUpper();
    std::vector<double> candidates{found.T, lower, upper, 0.5 * lower + 0.5 * upper};
    if (Contains(found.Bracket, 0.0)) {
        candidates.push_back(0.0);
    }
    return candidates;
}

/// An undecided or unisolated root, settled exactly when it lies at a
/// double (see MultiplicityAt): Decided, a crossing for odd multiplicity
/// and a tangency for even. An unisolated root that is not fails.
[[nodiscard]] Result<void> Settle(const ExactPolynomial& p, Found& found,
                                  const std::string_view handler)
{
    if (found.Decided) {
        return {};
    }
    for (const double u : Candidates(found)) {
        if (!Contains(found.Bracket, u)) {
            continue;
        }
        if (const std::optional<std::size_t> m = MultiplicityAt(p, u, found.Bracket)) {
            // u is the bracket's only root: the bracket shrinks to it.
            found = {u, Interval{u}, *m % 2 == 1 ? ContactType::Crossing : ContactType::Tangential,
                     true, true};
            return {};
        }
    }
    if (!found.Isolated) {
        return Fail(ErrorCode::IntersectionFailed,
                    std::string(handler) +
                        ": a contact of higher order that is not at a double could not be "
                        "isolated near chart parameter " +
                        Number(found.T) + " (ADR-0010: it is never dropped)");
    }
    return {};
}

// --- The charts ----------------------------------------------------------

/// A chart's points at infinity: a hyperbola's ends, a far parabola's 0.
[[nodiscard]] std::vector<double> InfinitePoints(const Shape shape, const Chart chart)
{
    if (shape == Shape::Hyperbola) {
        return {-1.0, 1.0};
    }
    if (shape == Shape::Parabola && chart == Chart::Far) {
        return {0.0};
    }
    return {};
}

/// A proved crossing's bracket, whose ends have strictly opposite exact
/// signs, narrowed by bisection on the doubles with Q's exact sign to two
/// adjacent doubles, or to the root itself when it is one.
void Tighten(const ExactPolynomial& p, Found& found)
{
    double a = found.Bracket.GetLower();
    double b = found.Bracket.GetUpper();
    const Sign atA = ValueAt(p, a).GetSign();
    if (atA == Sign::Zero || ValueAt(p, b).GetSign() != -atA) {
        return;
    }
    for (double m = 0.5 * a + 0.5 * b; a < m && m < b; m = 0.5 * a + 0.5 * b) {
        const Sign atM = ValueAt(p, m).GetSign();
        if (atM == Sign::Zero) {
            a = m;
            b = m;
            break;
        }
        (atM == atA ? a : b) = m;
    }
    found.T = 0.5 * a + 0.5 * b;
    found.Bracket = Interval{a, b};
}

/// The roots of Q on the chart u in [-1, 1], s = (u + 1) / 2 being the
/// subdivision's parameter, each with the exact ends u = 2 s - 1 of the
/// pieces that hold it.
[[nodiscard]] Result<std::vector<Found>> IsolateOnChart(const ExactPolynomial& p,
                                                        const std::string_view handler)
{
    Search search{handler, -1.0, 1.0, {}, 0};
    STEREON_TRY_VOID(Isolate(ExactPiece(p, 0.0, 1.0), 0, search, ExactSplit{p}));
    std::vector<Found> found;
    for (const LocalRoot& root : search.Roots) {
        found.push_back({2.0 * root.S - 1.0,
                         Interval{2.0 * root.Lower - 1.0, 2.0 * root.Upper - 1.0}, root.Contact,
                         root.Decided, root.Isolated});
    }
    Merge(found);
    return found;
}

/// The roots of Q on the chart, as chart parameters u with their brackets,
/// each settled or tightened, its roots at infinity divided out first.
[[nodiscard]] Result<std::vector<Found>> RootsOnChart(ExactPolynomial p, const Shape shape,
                                                      const Chart chart,
                                                      const std::string_view handler)
{
    const std::vector<double> infinite = InfinitePoints(shape, chart);
    for (const double r : infinite) {
        while (p.size() > 1 && ValueAt(p, r).IsZero()) {
            p = Deflated(p, r);
        }
    }
    STEREON_TRY(std::vector<Found> found, IsolateOnChart(p, handler));
    for (Found& f : found) {
        if (f.Decided) {
            Tighten(p, f);
        }
        STEREON_TRY_VOID(Settle(p, f, handler));
        for (const double r : infinite) {
            if (Contains(f.Bracket, r)) {
                return Fail(ErrorCode::IntersectionFailed,
                            std::string(handler) +
                                ": a meeting within rounding of a point at infinity, beyond the "
                                "range of the parameter in doubles");
            }
        }
    }
    return found;
}

/// A's parameter at chart parameter u.
[[nodiscard]] double ChartParameter(const Shape shape, const Chart chart, const double scale,
                                    const double u)
{
    if (shape == Shape::Ellipse) {
        return 2.0 * std::atan(u) + (chart == Chart::Far ? std::numbers::pi : 0.0);
    }
    if (shape == Shape::Hyperbola) {
        return 2.0 * std::atanh(u);
    }
    return chart == Chart::Near ? scale * u : scale / u;
}

/// x two doubles further out on each side: a margin for the rounding of
/// the chart's map.
[[nodiscard]] Interval Padded(const double lower, const double upper)
{
    const auto out = [](double x, const double towards) {
        x = std::nextafter(x, towards);
        return std::nextafter(x, towards);
    };
    return Interval{out(lower, -Infinity), out(upper, Infinity)};
}

/// The root in A's parameter, its bracket mapped end by end (each map is
/// monotone over a bracket, which holds no point at infinity).
[[nodiscard]] Found InParameter(Found f, const Shape shape, const Chart chart, const double scale)
{
    const double lower = ChartParameter(shape, chart, scale, f.Bracket.GetLower());
    const double upper = ChartParameter(shape, chart, scale, f.Bracket.GetUpper());
    f.T = ChartParameter(shape, chart, scale, f.T);
    f.Bracket = Padded(std::min(lower, upper), std::max(lower, upper));
    return f;
}

/// Roots seen from two charts at their common ends are one: brackets that
/// meet are merged, an ellipse's taken modulo its period. Two proved roots
/// are the same root proved twice; with an undecided one the merged root is
/// undecided, a touch when either is.
void MergeCharts(std::vector<Found>& found, const bool periodic)
{
    // Brackets of roots at negative angles, the near chart's lower half, turn
    // to meet the far chart's; the parameter is placed in the domain later.
    for (Found& f : found) {
        if (periodic && f.T < 0.0) {
            f.Bracket = f.Bracket + Interval{CirclePeriod};
        }
    }
    std::ranges::stable_sort(found, [](const Found& a, const Found& b) {
        return a.Bracket.GetLower() < b.Bracket.GetLower();
    });
    std::vector<Found> merged;
    for (const Found& f : found) {
        if (merged.empty() || !Intersects(f.Bracket, merged.back().Bracket)) {
            merged.push_back(f);
            continue;
        }
        Found& last = merged.back();
        if (last.Decided && f.Decided) {
            continue;
        }
        last.Bracket = Hull(last.Bracket, f.Bracket);
        last.Decided = false;
        if (f.Contact == ContactType::Tangential) {
            last.Contact = ContactType::Tangential;
            last.T = f.T;
        }
    }
    found = std::move(merged);
}

// --- Two conics ----------------------------------------------------------

/// The two curves with their conics and requested domains.
struct ConicPair {
    const Curve2d& CurveA;
    const Conic& A;
    const Interval& DomainA;
    const Curve2d& CurveB;
    const Conic& B;
    const Interval& DomainB;
    std::string_view Handler;
};

/// The point of a root: A's point at its parameter, B's parameter of it,
/// kept when both lie in the requested domains. A point on a hyperbola B's
/// other branch is not on B: its abscissa is at most -a, 2a from the
/// branch, far beyond rounding; one within a / 2 of the axis would be
/// neither, and fails.
[[nodiscard]] Result<void> KeepRoot(PlanarIntersection& result, const ConicPair& pair,
                                    const Found& f)
{
    const Point2 p = Eval(pair.CurveA, f.T);
    if (pair.B.Kind == Shape::Hyperbola) {
        const double abscissa = Dot(p - pair.B.Center, pair.B.X);
        if (std::abs(abscissa) < 0.5 * pair.B.A) {
            return Fail(ErrorCode::IntersectionFailed,
                        std::string(pair.Handler) + ": a root off both branches of the hyperbola");
        }
        if (abscissa < 0.0) {
            return {};
        }
    }
    KeepPoint(result, {f.T, ParameterOn(pair.B, p), p, f.Contact, f.Decided}, pair.DomainA,
              pair.A.Kind == Shape::Ellipse, pair.DomainB, pair.B.Kind == Shape::Ellipse);
    return {};
}

/// Coincident conics, Q vanishing identically: A lies on B's conic. Two
/// ellipses then coincide, with B's angle sense * s + offset for A's s; two
/// hyperbolas or two parabolas coincide, with B's parameter sense * t +
/// offset, unless A is the other branch of B's hyperbola, a proved
/// separation. Offset and sense are read from A's point and tangent at 0.
[[nodiscard]] Result<PlanarIntersection> Coincident(const ConicPair& pair)
{
    PlanarIntersection result;
    if (pair.A.Kind != pair.B.Kind) {
        return Fail(ErrorCode::IntersectionFailed,
                    std::string(pair.Handler) +
                        ": the conics' equations agree although their kinds differ");
    }
    const Point2 start = Eval(pair.CurveA, 0.0);
    if (pair.B.Kind == Shape::Hyperbola && Dot(start - pair.B.Center, pair.B.X) < 0.0) {
        return result;
    }
    const double offset = ParameterOn(pair.B, start);
    const double sense =
        Dot(Deriv(pair.CurveA, 0.0, 1), Deriv(pair.CurveB, offset, 1)) > 0.0 ? 1.0 : -1.0;
    if (pair.A.Kind == Shape::Ellipse) {
        PeriodicOverlaps(result, sense, offset, pair.DomainA, pair.DomainB);
    } else {
        // t_B = sense t_A + offset, so t_A = -sense offset + sense t_B.
        AffineOverlap(result, -sense * offset, sense, pair.DomainA, pair.DomainB);
    }
    return result;
}

/// The power of two nearest below 4f: the parabola's charts meet at
/// t = ±L, and the chart's coefficients stay exact.
[[nodiscard]] double ChartScale(const Conic& c)
{
    return c.Kind == Shape::Parabola ? std::ldexp(1.0, std::ilogb(4.0 * c.A)) : 1.0;
}

[[nodiscard]] Result<PlanarIntersection> Meet(const ConicPair& pair)
{
    const double scale = ChartScale(pair.A);
    const ExactPolynomial near =
        Composite(Exact{}, pair.A, ChartOf(Exact{}, pair.A, Chart::Near, scale), pair.B);
    if (IsZero(near)) {
        return Coincident(pair);
    }
    std::vector<Found> found;
    for (const Chart chart : {Chart::Near, Chart::Far}) {
        if (chart == Chart::Far && pair.A.Kind == Shape::Hyperbola) {
            continue;
        }
        ExactPolynomial q =
            chart == Chart::Near
                ? near
                : Composite(Exact{}, pair.A, ChartOf(Exact{}, pair.A, chart, scale), pair.B);
        STEREON_TRY(const std::vector<Found> roots,
                    RootsOnChart(std::move(q), pair.A.Kind, chart, pair.Handler));
        for (const Found& f : roots) {
            found.push_back(InParameter(f, pair.A.Kind, chart, scale));
        }
    }
    MergeCharts(found, pair.A.Kind == Shape::Ellipse);
    PlanarIntersection result;
    for (const Found& f : found) {
        STEREON_TRY_VOID(KeepRoot(result, pair, f));
    }
    Order(result);
    return result;
}

} // namespace

Result<PlanarIntersection> ConicConic2d(const Curve2d& a, const Interval& domainA, const Curve2d& b,
                                        const Interval& domainB, const Context& ctx)
{
    const std::optional<Conic> ca = ConicOf(a);
    const std::optional<Conic> cb = ConicOf(b);
    if (!ca || !cb) {
        return Fail(ErrorCode::InvalidInput, "ConicConic2d: dispatched for other curve types");
    }
    if (ChartedSecond(*ca, *cb)) {
        STEREON_TRY(PlanarIntersection result, ConicConic2d(b, domainB, a, domainA, ctx));
        return SwappedPlanar(std::move(result));
    }
    const std::string handler = PairName(a, b);
    return Meet({a, *ca, domainA, b, *cb, domainB, handler});
}

} // namespace Stereon::Detail
