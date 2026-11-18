// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The subdivision machinery the curve-pair handlers share, in the plane and
// in space: an arc, a dyadic piece of one Bezier segment of a B-spline with
// its homogeneous poles enclosed (Bernstein.hpp's SegmentsOf), halved by
// de Casteljau in interval arithmetic, evaluated with its derivatives over
// a parameter interval, and the bookkeeping of a pair's parameter box:
// Krawczyk's test outcome, where a root's bracket lies against its
// segment, and the curve parameter a piece-local bracket stands for.
//
// A plane curve's poles have Z zero, which every step here keeps exactly
// zero, so the plane handlers see the arithmetic they always had.
//
// Included after Bernstein.hpp, in the same purview and with the same
// standard headers, plus <array> and <inplace_vector>. Everything is in an
// anonymous namespace.

#pragma once

namespace Stereon::Detail {
namespace {

/// The most poles an arc or one of its hodographs has: a Bezier segment of
/// a B-spline of the highest degree.
constexpr std::size_t ArcCapacity = static_cast<std::size_t>(BSplineMaxDegree) + 1;

/// A piece [S0, S1] of segment Index's local parameter, dyadic, with its
/// homogeneous poles over it, enclosed.
struct Arc {
    std::size_t Index;
    double S0;
    double S1;
    std::vector<HPoint> Poles;
};

[[maybe_unused]] [[nodiscard]] inline HPoint Mid(const HPoint& a, const HPoint& b)
{
    const Interval half{0.5};
    return {(a.X + b.X) * half, (a.Y + b.Y) * half, (a.Z + b.Z) * half, (a.W + b.W) * half};
}

/// de Casteljau at the middle, enclosed, reduced in the right half's own
/// storage as Bernstein.hpp's Halves is.
[[maybe_unused]] [[nodiscard]] inline std::pair<Arc, Arc> Halves(const Arc& arc)
{
    const std::size_t m = arc.Poles.size();
    std::vector<HPoint> left(m);
    std::vector<HPoint> right = arc.Poles;
    left[0] = right[0];
    for (std::size_t level = 1; level < m; ++level) {
        for (std::size_t i = 0; i + level < m; ++i) {
            right[i] = Mid(right[i], right[i + 1]);
        }
        left[level] = right[0];
    }
    const double mid = 0.5 * arc.S0 + 0.5 * arc.S1;
    return {Arc{arc.Index, arc.S0, mid, std::move(left)},
            Arc{arc.Index, mid, arc.S1, std::move(right)}};
}

/// b + s (c - b): de Casteljau's step in the form whose enclosure shrinks
/// with the poles' spread, for s anywhere, inside [0, 1] or not.
[[maybe_unused]] [[nodiscard]] inline Interval Step(const Interval& b, const Interval& c,
                                                    const Interval& s)
{
    return b + s * (c - b);
}

/// The homogeneous point of the poles' Bezier polynomial at s, enclosed.
/// The reduction's work array is on the stack (std::inplace_vector): the
/// pair searches evaluate an arc and its hodographs at every box they
/// visit, and none of those evaluations may cost a heap allocation.
[[maybe_unused]] [[nodiscard]] inline HPoint EvalH(const std::vector<HPoint>& poles,
                                                   const Interval& s)
{
    STEREON_HOT_ASSERT(poles.size() <= ArcCapacity);
    std::inplace_vector<HPoint, ArcCapacity> work(poles.begin(), poles.end());
    for (std::size_t level = 1; level < work.size(); ++level) {
        for (std::size_t i = 0; i + level < work.size(); ++i) {
            work[i] = {Step(work[i].X, work[i + 1].X, s), Step(work[i].Y, work[i + 1].Y, s),
                       Step(work[i].Z, work[i + 1].Z, s), Step(work[i].W, work[i + 1].W, s)};
        }
    }
    return work.front();
}

/// The hodograph's poles, p (H_i+1 - H_i), of degree p - 1; one zero pole
/// for a constant.
[[maybe_unused]] [[nodiscard]] inline std::vector<HPoint>
Hodograph(const std::vector<HPoint>& poles)
{
    const Interval p{static_cast<double>(poles.size() - 1)};
    std::vector<HPoint> d;
    d.reserve(poles.size());
    for (std::size_t i = 0; i + 1 < poles.size(); ++i) {
        d.push_back({p * (poles[i + 1].X - poles[i].X), p * (poles[i + 1].Y - poles[i].Y),
                     p * (poles[i + 1].Z - poles[i].Z), p * (poles[i + 1].W - poles[i].W)});
    }
    if (d.empty()) {
        d.push_back({Interval{0.0}, Interval{0.0}, Interval{0.0}, Interval{0.0}});
    }
    return d;
}

/// The middle of a pole's enclosure, as a point.
[[maybe_unused]] [[nodiscard]] inline Point3 Nominal3(const HPoint& h)
{
    const double w = Midpoint(h.W);
    return {Midpoint(h.X) / w, Midpoint(h.Y) / w, Midpoint(h.Z) / w};
}

// --- A pair's parameter box ----------------------------------------------

/// The outcome of Krawczyk's test on a box.
enum class Test : std::uint8_t { NoRoot, OneRoot, Unknown };

/// The parameter box of a pair: S on A's piece, T on B's, local to each
/// piece ([0, 1] is the piece itself).
struct ParameterBox {
    Interval S;
    Interval T;
};

[[maybe_unused]] [[nodiscard]] inline bool InsideInterior(const Interval& inner,
                                                          const Interval& outer)
{
    return outer.GetLower() < inner.GetLower() && inner.GetUpper() < outer.GetUpper();
}

/// Krawczyk's outcome from K(X): no root when K(X) misses X, exactly one
/// when it lies inside X's interior, otherwise unknown.
[[maybe_unused]] [[nodiscard]] inline Test Classify(const std::optional<ParameterBox>& k,
                                                    const ParameterBox& x)
{
    if (!k) {
        return Test::Unknown;
    }
    if (!Intersects(k->S, x.S) || !Intersects(k->T, x.T)) {
        return Test::NoRoot;
    }
    return InsideInterior(k->S, x.S) && InsideInterior(k->T, x.T) ? Test::OneRoot : Test::Unknown;
}

/// K(X) for a two-by-two system G, from G at X's middle and its Jacobian J
/// enclosed over X, preconditioned by the inverse of J's midpoint:
/// K(X) = m - Y G(m) + (I - Y J(X)) (X - m). Nothing when J's midpoint is
/// singular.
[[maybe_unused]] [[nodiscard]] inline std::optional<ParameterBox>
KrawczykStep(const ParameterBox& x, const std::array<Interval, 2>& g,
             const std::array<Interval, 4>& j)
{
    const double ms = Midpoint(x.S);
    const double mt = Midpoint(x.T);
    const double c11 = Midpoint(j[0]);
    const double c12 = Midpoint(j[1]);
    const double c21 = Midpoint(j[2]);
    const double c22 = Midpoint(j[3]);
    const double det = c11 * c22 - c12 * c21;
    if (!(std::abs(det) > 0.0) || !std::isfinite(det)) {
        return std::nullopt;
    }
    const std::array<Interval, 4> y{Interval{c22 / det}, Interval{-c12 / det}, Interval{-c21 / det},
                                    Interval{c11 / det}};
    const Interval ds = x.S - Interval{ms};
    const Interval dt = x.T - Interval{mt};
    const auto row = [&](const std::size_t r, const double m, const double unitS,
                         const double unitT) {
        const Interval yg = y[2 * r] * g[0] + y[2 * r + 1] * g[1];
        const Interval rs = Interval{unitS} - (y[2 * r] * j[0] + y[2 * r + 1] * j[2]);
        const Interval rt = Interval{unitT} - (y[2 * r] * j[1] + y[2 * r + 1] * j[3]);
        return Interval{m} - yg + rs * ds + rt * dt;
    };
    return ParameterBox{row(0, ms, 1.0, 0.0), row(1, mt, 0.0, 1.0)};
}

// --- Brackets ------------------------------------------------------------

/// A piece-local box in a segment's local parameter.
[[maybe_unused]] [[nodiscard]] inline Interval InSegment(const Arc& arc, const Interval& local)
{
    return Interval{arc.S0} + (Interval{arc.S1} - Interval{arc.S0}) * local;
}

/// A segment's local parameter in the curve's.
[[maybe_unused]] [[nodiscard]] inline Interval OnCurve(const Segment& segment, const Interval& s)
{
    return Interval{segment.T0} + (Interval{segment.T1} - Interval{segment.T0}) * s;
}

/// Where a root's bracket in a segment's parameter lies against [0, 1]:
/// inside (a root of the curve), beyond (of the polynomial continued, not
/// of the segment), or across an end (undecided).
enum class Within : std::uint8_t { Inside, Beyond, Across };

[[maybe_unused]] [[nodiscard]] inline Within WithinSegment(const Interval& s)
{
    if (Contains(Interval{0.0, 1.0}, s)) {
        return Within::Inside;
    }
    return Intersects(s, Interval{0.0, 1.0}) ? Within::Across : Within::Beyond;
}

/// The segments that reach the requested domain.
[[maybe_unused]] [[nodiscard]] inline std::vector<Segment> Reaching(std::vector<Segment> segments,
                                                                    const Interval& domain)
{
    std::erase_if(segments,
                  [&](const Segment& s) { return !Intersects(Interval{s.T0, s.T1}, domain); });
    return segments;
}

/// A root's parameter in the requested domain: its bracket's middle, or,
/// for a bracket across the domain's end, the end, undecided; nothing for
/// a bracket outside the domain.
[[maybe_unused]] [[nodiscard]] inline std::optional<double>
Placed(const Interval& bracket, const Interval& domain, bool& decided)
{
    if (!Intersects(bracket, domain)) {
        return std::nullopt;
    }
    const double t = Midpoint(bracket);
    if (!Contains(domain, bracket)) {
        decided = false;
    }
    return std::clamp(t, domain.GetLower(), domain.GetUpper());
}

} // namespace
} // namespace Stereon::Detail
