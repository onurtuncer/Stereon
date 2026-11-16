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
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <queue>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

module Stereon.Intersect;

import Stereon.Robust;
import :Dispatch;
import :ClosestApproach;

#include "ApproachSide.hpp"

/// The primal closest approach (ClosestApproach.cppm): two lines in closed
/// form, decided parallel exactly; a line and a conic, or two conics, by
/// branch and bound on phi = |A(s) - B(t)|^2 / 2, then an active-set Newton
/// refinement and a proof that phi is strictly convex over the one cluster
/// of boxes that may hold its minimum.

namespace Stereon {

namespace {

using Detail::ApproachSide;

constexpr double Epsilon = std::numeric_limits<double>::epsilon();
constexpr double Tau = 2.0 * std::numbers::pi;

/// Boxes the search may process before it gives up proving the minimum: a
/// budget, not a tolerance (ADR-0024). A unique minimum needs a few
/// thousand; a set of minimisers exhausts it.
constexpr std::size_t BoxBudget = 200000;
/// Leaves the clustering compares pairwise.
constexpr std::size_t LeafBudget = 4000;
/// A leaf's width as a fraction of its dimension's range: 2^-22.
constexpr double LeafFraction = 0x1p-22;

// --- The parameter ranges ----------------------------------------------------

/// One parameter's searched range: its bounds, whether each is an end of the
/// requested domain (a line's range is cut to where its minimiser can lie,
/// which is no end), and whether it wraps (a whole conic).
struct Range {
    double Lo = 0.0;
    double Hi = 0.0;
    bool LoIsEnd = false;
    bool HiIsEnd = false;
    bool Periodic = false;
    /// The end a range of one point is: the domain's end it sits on.
    ApproachEnd Pin = ApproachEnd::Lower;

    [[nodiscard]] double Width() const { return Hi - Lo; }
};

/// A conic's range: the whole conic for a domain 2 pi wide or unbounded,
/// otherwise its arc, both ends real.
[[nodiscard]] Range ConicRange(const Interval& domain)
{
    if (!domain.IsBounded() || domain.GetUpper() - domain.GetLower() >= Tau) {
        const double lo = std::isfinite(domain.GetLower()) ? domain.GetLower() : 0.0;
        return {lo, lo + Tau, false, false, true};
    }
    return {domain.GetLower(), domain.GetUpper(), true, true, false};
}

/// A line's range against a conic of the given centre and reach: for any
/// point of the conic the nearest point of the whole line lies within
/// [m, M], so the minimiser does; the domain's ends inside that window stay
/// ends.
[[nodiscard]] Range LineRange(const ApproachSide& line, const Interval& domain,
                              const Point3& center, const double reach)
{
    const double dd = Dot(line.D, line.D);
    const double mid = Dot(center - line.Origin, line.D) / dd;
    const double half = reach / std::sqrt(dd);
    const double m = mid - half * (1.0 + 4.0 * Epsilon) - 4.0 * Epsilon * std::abs(mid);
    const double big = mid + half * (1.0 + 4.0 * Epsilon) + 4.0 * Epsilon * std::abs(mid);
    if (domain.GetLower() > big) {
        return {domain.GetLower(), domain.GetLower(), true, true, false};
    }
    if (domain.GetUpper() < m) {
        return {domain.GetUpper(), domain.GetUpper(), true, true, false, ApproachEnd::Upper};
    }
    const bool loEnd = domain.GetLower() >= m;
    const bool hiEnd = domain.GetUpper() <= big;
    return {loEnd ? domain.GetLower() : m, hiEnd ? domain.GetUpper() : big, loEnd, hiEnd, false};
}

// --- phi, its derivatives and its bounds ------------------------------------

/// phi = |d|^2 / 2 with d = A(s) - B(t), its gradient (d . A', -d . B') and
/// its Hessian.
struct Local {
    Vec3 Gap;
    double Phi = 0.0;
    std::array<double, 2> Grad{};
    std::array<double, 3> Hess{}; // ss, st, tt
};

[[nodiscard]] Local LocalAt(const ApproachSide& a, const ApproachSide& b, const double s,
                            const double t)
{
    Local l;
    l.Gap = (a.Origin - b.Origin) + (a.Offset(s) - b.Offset(t));
    const Vec3 a1 = a.D1(s);
    const Vec3 b1 = b.D1(t);
    l.Phi = 0.5 * Dot(l.Gap, l.Gap);
    l.Grad = {Dot(l.Gap, a1), -Dot(l.Gap, b1)};
    l.Hess = {Dot(a1, a1) + Dot(l.Gap, a.D2(s)), -Dot(a1, b1), Dot(b1, b1) - Dot(l.Gap, b.D2(t))};
    return l;
}

/// The smaller eigenvalue of a symmetric 2 x 2 matrix.
[[nodiscard]] double SmallerEigenvalue(const std::array<double, 3>& h)
{
    const double mean = 0.5 * (h[0] + h[2]);
    const double half = 0.5 * (h[0] - h[2]);
    return mean - std::hypot(half, h[1]);
}

/// A bound on the sum of the absolute entries of phi's Hessian, at least its
/// spectral norm, where |d| <= reach.
[[nodiscard]] double SecondBound(const ApproachSide& a, const ApproachSide& b, const double reach)
{
    return a.L1 * a.L1 + reach * a.L2 + b.L1 * b.L1 + reach * b.L2 + 2.0 * a.L1 * b.L1;
}

/// A bound on how fast phi's Hessian changes: its entries' gradients summed,
/// where |d| <= reach.
[[nodiscard]] double ThirdBound(const ApproachSide& a, const ApproachSide& b, const double reach)
{
    const double ss = 3.0 * a.L1 * a.L2 + reach * a.L3 + a.L2 * b.L1;
    const double st = a.L2 * b.L1 + a.L1 * b.L2;
    const double tt = 3.0 * b.L1 * b.L2 + reach * b.L3 + a.L1 * b.L2;
    return ss + 2.0 * st + tt;
}

// --- Branch and bound --------------------------------------------------------

struct Box {
    double S0, S1, T0, T1;
    double Lower; // a lower bound of phi over the box
    double AtCentre;

    [[nodiscard]] double Cs() const { return 0.5 * (S0 + S1); }
    [[nodiscard]] double Ct() const { return 0.5 * (T0 + T1); }
};

struct ByLower {
    bool operator()(const Box& x, const Box& y) const { return x.Lower > y.Lower; }
};

struct Search {
    const ApproachSide& A;
    const ApproachSide& B;
    Range RangeS;
    Range RangeT;
    double Slack = 0.0; // the rounding of phi: boxes within it of the best stay
    double Best = std::numeric_limits<double>::infinity();
    double BestS = 0.0;
    double BestT = 0.0;
    std::vector<Box> Leaves{};
    std::size_t Processed = 0;
    bool Exhausted = false;

    /// phi at the box's centre and a lower bound over it: phi(c) - |g_s| h_s
    /// - |g_t| h_t - M2 (h_s^2 + h_t^2) / 2, M2 at the largest |d| the box
    /// can reach.
    [[nodiscard]] Box Bounded(const double s0, const double s1, const double t0, const double t1)
    {
        Box box{s0, s1, t0, t1, 0.0, 0.0};
        const double hs = 0.5 * (s1 - s0);
        const double ht = 0.5 * (t1 - t0);
        const Local l = LocalAt(A, B, box.Cs(), box.Ct());
        const double reach = std::sqrt(2.0 * l.Phi) + A.L1 * hs + B.L1 * ht;
        box.AtCentre = l.Phi;
        box.Lower = l.Phi - std::abs(l.Grad[0]) * hs - std::abs(l.Grad[1]) * ht -
                    0.5 * SecondBound(A, B, reach) * (hs * hs + ht * ht);
        Offer(box.Cs(), box.Ct(), l.Phi);
        // The corner the gradient descends to: on a box against an end the
        // minimum sits there, and the best value found must reach it, or
        // every box along that end would stay.
        const double s = l.Grad[0] > 0.0 ? s0 : s1;
        const double t = l.Grad[1] > 0.0 ? t0 : t1;
        Offer(s, t, LocalAt(A, B, s, t).Phi);
        ++Processed;
        return box;
    }

    void Offer(const double s, const double t, const double phi)
    {
        if (phi < Best) {
            Best = phi;
            BestS = s;
            BestT = t;
        }
    }

    [[nodiscard]] bool IsLeaf(const Box& box) const
    {
        return box.S1 - box.S0 <= LeafFraction * RangeS.Width() &&
               box.T1 - box.T0 <= LeafFraction * RangeT.Width();
    }

    /// Splits along the dimension the bound is loosest in.
    void Split(const Box& box, std::priority_queue<Box, std::vector<Box>, ByLower>& queue)
    {
        const bool alongS = A.L1 * (box.S1 - box.S0) >= B.L1 * (box.T1 - box.T0);
        const auto keep = [&](const Box& child) {
            if (child.Lower <= Best + Slack) {
                queue.push(child);
            }
        };
        if (alongS) {
            keep(Bounded(box.S0, box.Cs(), box.T0, box.T1));
            keep(Bounded(box.Cs(), box.S1, box.T0, box.T1));
        } else {
            keep(Bounded(box.S0, box.S1, box.T0, box.Ct()));
            keep(Bounded(box.S0, box.S1, box.Ct(), box.T1));
        }
    }

    void Run()
    {
        std::priority_queue<Box, std::vector<Box>, ByLower> queue;
        queue.push(Bounded(RangeS.Lo, RangeS.Hi, RangeT.Lo, RangeT.Hi));
        while (!queue.empty()) {
            const Box box = queue.top();
            queue.pop();
            if (box.Lower > Best + Slack) {
                break;
            }
            if (Processed > BoxBudget || Leaves.size() > LeafBudget) {
                Exhausted = true;
                return;
            }
            if (IsLeaf(box)) {
                Leaves.push_back(box);
            } else {
                Split(box, queue);
            }
        }
        std::erase_if(Leaves, [&](const Box& b) { return b.Lower > Best + Slack; });
    }

    /// Whether a leaf touches the seam of a whole conic's range.
    [[nodiscard]] bool TouchesSeam() const
    {
        return std::ranges::any_of(Leaves, [&](const Box& b) {
            return (RangeS.Periodic && (b.S0 <= RangeS.Lo || b.S1 >= RangeS.Hi)) ||
                   (RangeT.Periodic && (b.T0 <= RangeT.Lo || b.T1 >= RangeT.Hi));
        });
    }
};

/// Disjoint sets of leaves, with path halving.
struct Groups {
    std::vector<std::size_t> Parent;

    explicit Groups(const std::size_t n)
        : Parent(n)
    {
        for (std::size_t i = 0; i < n; ++i) {
            Parent[i] = i;
        }
    }
    [[nodiscard]] std::size_t Root(std::size_t i)
    {
        while (Parent[i] != i) {
            i = Parent[i] = Parent[Parent[i]];
        }
        return i;
    }
};

/// Whether two leaves touch, within a gap.
[[nodiscard]] bool Touch(const Box& x, const Box& y, const double gs, const double gt)
{
    return x.S0 <= y.S1 + gs && y.S0 <= x.S1 + gs && x.T0 <= y.T1 + gt && y.T0 <= x.T1 + gt;
}

/// The hull of two boxes, keeping the smaller bounds.
[[nodiscard]] Box Joined(const Box& h, const Box& x)
{
    return {std::min(h.S0, x.S0), std::max(h.S1, x.S1),       std::min(h.T0, x.T0),
            std::max(h.T1, x.T1), std::min(h.Lower, x.Lower), std::min(h.AtCentre, x.AtCentre)};
}

/// The connected groups of leaves, as the hull of each.
[[nodiscard]] std::vector<Box> Clusters(const std::vector<Box>& leaves, const Range& rs,
                                        const Range& rt)
{
    Groups groups(leaves.size());
    const double gs = 1e-9 * rs.Width();
    const double gt = 1e-9 * rt.Width();
    for (std::size_t i = 0; i < leaves.size(); ++i) {
        for (std::size_t j = i + 1; j < leaves.size(); ++j) {
            if (Touch(leaves[i], leaves[j], gs, gt)) {
                groups.Parent[groups.Root(i)] = groups.Root(j);
            }
        }
    }
    std::vector<Box> hulls;
    std::vector<std::size_t> index(leaves.size(), leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i) {
        const std::size_t r = groups.Root(i);
        if (index[r] == leaves.size()) {
            index[r] = hulls.size();
            hulls.push_back(leaves[i]);
        } else {
            hulls[index[r]] = Joined(hulls[index[r]], leaves[i]);
        }
    }
    return hulls;
}

// --- Refinement ------------------------------------------------------------

/// A point and which bounds it rests on.
struct Refined {
    double S = 0.0;
    double T = 0.0;
    ApproachEnd EndS = ApproachEnd::Interior;
    ApproachEnd EndT = ApproachEnd::Interior;
};

/// Clamped to the range's real ends; which end it rests on.
[[nodiscard]] std::pair<double, ApproachEnd> Clamp(const double x, const Range& r)
{
    if (r.LoIsEnd && r.HiIsEnd && r.Lo == r.Hi) {
        return {r.Lo, r.Pin};
    }
    if (r.LoIsEnd && x <= r.Lo) {
        return {r.Lo, ApproachEnd::Lower};
    }
    if (r.HiIsEnd && x >= r.Hi) {
        return {r.Hi, ApproachEnd::Upper};
    }
    return {x, ApproachEnd::Interior};
}

/// Releases an end the gradient pulls the point away from (KKT).
[[nodiscard]] ApproachEnd Released(const ApproachEnd end, const double grad)
{
    if ((end == ApproachEnd::Lower && grad < 0.0) || (end == ApproachEnd::Upper && grad > 0.0)) {
        return ApproachEnd::Interior;
    }
    return end;
}

/// The Newton step on the free parameters, or the steepest descent where
/// the Hessian is not positive.
[[nodiscard]] std::array<double, 2> Step(const Local& l, const bool freeS, const bool freeT)
{
    const double hs = l.Hess[0];
    const double ht = l.Hess[2];
    if (freeS && freeT) {
        const double det = hs * ht - l.Hess[1] * l.Hess[1];
        if (hs > 0.0 && det > 0.0) {
            return {(-ht * l.Grad[0] + l.Hess[1] * l.Grad[1]) / det,
                    (l.Hess[1] * l.Grad[0] - hs * l.Grad[1]) / det};
        }
        const double scale = 1.0 / std::max(std::abs(hs) + std::abs(ht), 1e-300);
        return {-scale * l.Grad[0], -scale * l.Grad[1]};
    }
    const auto one = [](const double h, const double g) {
        return h > 0.0 ? -g / h : -g / std::max(std::abs(h), 1e-300);
    };
    return {freeS ? one(hs, l.Grad[0]) : 0.0, freeT ? one(ht, l.Grad[1]) : 0.0};
}

/// Active-set Newton with backtracking from (s, t), within the ranges' real
/// ends: on a strictly convex region, the minimiser.
[[nodiscard]] Refined Refine(const ApproachSide& a, const ApproachSide& b, const Range& rs,
                             const Range& rt, double s, double t)
{
    Refined z{s, t};
    for (int iteration = 0; iteration < 200; ++iteration) {
        const Local l = LocalAt(a, b, z.S, z.T);
        z.EndS = Released(z.EndS, l.Grad[0]);
        z.EndT = Released(z.EndT, l.Grad[1]);
        const auto [ds, dt] =
            Step(l, z.EndS == ApproachEnd::Interior, z.EndT == ApproachEnd::Interior);
        double alpha = 1.0;
        Refined next = z;
        for (int halving = 0; halving < 60; ++halving, alpha *= 0.5) {
            std::tie(next.S, next.EndS) = Clamp(z.S + alpha * ds, rs);
            std::tie(next.T, next.EndT) = Clamp(z.T + alpha * dt, rt);
            if (LocalAt(a, b, next.S, next.T).Phi <= l.Phi) {
                break;
            }
        }
        const bool still = std::abs(next.S - z.S) <= 4.0 * Epsilon * (1.0 + std::abs(z.S)) &&
                           std::abs(next.T - z.T) <= 4.0 * Epsilon * (1.0 + std::abs(z.T));
        z = next;
        if (still && iteration > 2) {
            break;
        }
    }
    return z;
}

/// Whether a parameter stays where it rests: free, a range of one point,
/// or an end its gradient pushes outward by more than `moved`, the most the
/// gradient changes over the hull.
[[nodiscard]] bool Held(const ApproachEnd end, const double grad, const bool point,
                        const double moved)
{
    if (end == ApproachEnd::Interior || point) {
        return true;
    }
    const double outward = end == ApproachEnd::Lower ? grad : -grad;
    return outward > moved;
}

/// The smaller curvature of phi over the free parameters; infinite when
/// none is free.
[[nodiscard]] double FreeCurvature(const std::array<double, 3>& hess, const bool freeS,
                                   const bool freeT)
{
    if (freeS && freeT) {
        return SmallerEigenvalue(hess);
    }
    if (freeS || freeT) {
        return freeS ? hess[0] : hess[2];
    }
    return std::numeric_limits<double>::infinity();
}

/// Whether the hull holds one minimiser, around the refined point. A
/// parameter resting on an end whose gradient pushes it outward by more than
/// the gradient can change over the hull keeps every minimiser on that end;
/// over the parameters left free, phi is strictly convex when its Hessian
/// there outweighs the change a third-derivative bound allows out to the
/// hull's farthest corner. A zero-width dimension is no dimension.
[[nodiscard]] bool ProvedUnique(const ApproachSide& a, const ApproachSide& b, const Box& hull,
                                const Refined& z)
{
    const double rs = std::max(std::abs(hull.S0 - z.S), std::abs(hull.S1 - z.S));
    const double rt = std::max(std::abs(hull.T0 - z.T), std::abs(hull.T1 - z.T));
    const double rho = std::hypot(rs, rt);
    const Local l = LocalAt(a, b, z.S, z.T);
    const double reach = std::sqrt(2.0 * l.Phi) + (a.L1 + b.L1) * rho;
    const double rounding = 64.0 * Epsilon * (std::abs(l.Hess[0]) + std::abs(l.Hess[2]));
    const double moved = 1.25 * SecondBound(a, b, reach) * rho + rounding;
    if (!Held(z.EndS, l.Grad[0], hull.S1 == hull.S0, moved) ||
        !Held(z.EndT, l.Grad[1], hull.T1 == hull.T0, moved)) {
        return false;
    }
    const bool freeS = z.EndS == ApproachEnd::Interior && hull.S1 > hull.S0;
    const bool freeT = z.EndT == ApproachEnd::Interior && hull.T1 > hull.T0;
    return FreeCurvature(l.Hess, freeS, freeT) > 1.25 * ThirdBound(a, b, reach) * rho + rounding;
}

// --- Assembling the result -------------------------------------------------

/// A whole conic's parameter in [Lo, Lo + 2 pi).
[[nodiscard]] double Wrapped(const double x, const Range& r)
{
    if (!r.Periodic) {
        return x;
    }
    const double k = std::floor((x - r.Lo) / Tau);
    return x - k * Tau;
}

[[nodiscard]] ClosestApproachResult Assemble(const ApproachSide& a, const ApproachSide& b,
                                             const Range& rs, const Range& rt, const Refined& z,
                                             const ApproachUniqueness uniqueness,
                                             const double lower, const Context& ctx)
{
    ClosestApproachResult r;
    r.ParameterA = Wrapped(z.S, rs);
    r.ParameterB = Wrapped(z.T, rt);
    r.OnA = a.At(r.ParameterA);
    r.OnB = b.At(r.ParameterB);
    r.Gap = r.OnA - r.OnB;
    r.Distance = Length(r.Gap);
    r.LowerBound = std::min(r.Distance, std::max(lower, 0.0));
    r.EndA = z.EndS;
    r.EndB = z.EndT;
    r.Uniqueness = uniqueness;
    r.Meets = r.Distance <= ctx.GetResolution();
    return r;
}

// --- Exact decisions ---------------------------------------------------------

using Exact = DynamicExpansion;

[[nodiscard]] std::array<Exact, 3> ExactOf(const Vec3& v)
{
    return {Exact{v.X}, Exact{v.Y}, Exact{v.Z}};
}

[[nodiscard]] std::array<Exact, 3> ExactDifference(const Point3& p, const Point3& q)
{
    return {Difference(Exact{p.X}, Exact{q.X}), Difference(Exact{p.Y}, Exact{q.Y}),
            Difference(Exact{p.Z}, Exact{q.Z})};
}

[[nodiscard]] std::array<Exact, 3> ExactCross(const std::array<Exact, 3>& u,
                                              const std::array<Exact, 3>& v)
{
    return {Difference(Product(u[1], v[2]), Product(u[2], v[1])),
            Difference(Product(u[2], v[0]), Product(u[0], v[2])),
            Difference(Product(u[0], v[1]), Product(u[1], v[0]))};
}

[[nodiscard]] Exact ExactDot(const std::array<Exact, 3>& u, const std::array<Exact, 3>& v)
{
    return Sum(Sum(Product(u[0], v[0]), Product(u[1], v[1])), Product(u[2], v[2]));
}

[[nodiscard]] bool IsZero(const std::array<Exact, 3>& v)
{
    return std::ranges::all_of(v, [](const Exact& x) { return x.GetSign() == Sign::Zero; });
}

/// The line lies on the circle's axis: its direction is perpendicular to
/// both frame vectors and the circle's centre is on it, exactly.
[[nodiscard]] bool OnAxis(const ApproachSide& line, const ApproachSide& circle)
{
    const auto d = ExactOf(line.D);
    return ExactDot(d, ExactOf(circle.X)).GetSign() == Sign::Zero &&
           ExactDot(d, ExactOf(circle.Y)).GetSign() == Sign::Zero &&
           IsZero(ExactCross(ExactDifference(circle.Origin, line.Origin), d));
}

/// Two circles share their axis: the second's frame vectors and the
/// offset of its centre are perpendicular to the first's normal, written as
/// triple products, and its centre is on the first's axis.
[[nodiscard]] bool Coaxial(const ApproachSide& p, const ApproachSide& q)
{
    const auto n = ExactCross(ExactOf(p.X), ExactOf(p.Y));
    const auto offset = ExactDifference(q.Origin, p.Origin);
    return ExactDot(n, ExactOf(q.X)).GetSign() == Sign::Zero &&
           ExactDot(n, ExactOf(q.Y)).GetSign() == Sign::Zero &&
           ExactDot(offset, ExactOf(p.X)).GetSign() == Sign::Zero &&
           ExactDot(offset, ExactOf(p.Y)).GetSign() == Sign::Zero;
}

// --- Two lines -------------------------------------------------------------

/// A coordinate in the overlap [lo, hi] of two parallel lines' ranges, for a
/// representative minimiser: its middle, or its finite end, or zero.
[[nodiscard]] double Common(const double lo, const double hi)
{
    if (std::isfinite(lo) && std::isfinite(hi)) {
        return 0.5 * (lo + hi);
    }
    if (std::isfinite(lo) || std::isfinite(hi)) {
        return std::isfinite(lo) ? lo : hi;
    }
    return 0.0;
}

/// Two parallel lines: phi depends only on q_A(s) - q_B(t), their
/// coordinates along A's direction. Ranges that overlap in more than a
/// point make a set of minimisers; disjoint ranges a unique pair of ends.
[[nodiscard]] ClosestApproachResult ParallelLines(const ApproachSide& a, const Range& rs,
                                                  const ApproachSide& b, const Range& rt,
                                                  const Context& ctx)
{
    const double qa0 = Dot(a.Origin - b.Origin, a.D);
    const double ka = Dot(a.D, a.D);
    const double kb = Dot(b.D, a.D);
    // Each line's range of q, with the parameter at each end.
    const auto span = [](const double base, const double k, const Range& r) {
        const double x = base + k * r.Lo;
        const double y = base + k * r.Hi;
        return k > 0.0 ? std::array{x, y, r.Lo, r.Hi} : std::array{y, x, r.Hi, r.Lo};
    };
    const auto qa = span(qa0, ka, rs);
    const auto qb = span(0.0, kb, rt);
    const double scale =
        16.0 * Epsilon * (std::abs(qa0) + ka * (std::abs(rs.Lo) + std::abs(rs.Hi)) + 1.0);
    const double gapAB = qb[0] - qa[1]; // B entirely beyond A
    const double gapBA = qa[0] - qb[1]; // A entirely beyond B
    Refined z;
    ApproachUniqueness uniqueness = ApproachUniqueness::Unique;
    if (gapAB > scale) {
        z = {qa[3], qb[2]};
    } else if (gapBA > scale) {
        z = {qa[2], qb[3]};
    } else {
        const bool apart = gapAB > -scale || gapBA > -scale;
        uniqueness = apart ? ApproachUniqueness::Undecided : ApproachUniqueness::NotUnique;
        const double q = Common(std::max(qa[0], qb[0]), std::min(qa[1], qb[1]));
        z = {std::clamp((q - qa0) / ka, rs.Lo, rs.Hi), std::clamp(q / kb, rt.Lo, rt.Hi)};
    }
    z.EndS = Clamp(z.S, rs).second;
    z.EndT = Clamp(z.T, rt).second;
    const double lower = uniqueness == ApproachUniqueness::Undecided
                             ? 0.0
                             : Length(Cross(a.Origin - b.Origin, a.D)) / std::sqrt(ka);
    return Assemble(a, b, rs, rt, z, uniqueness, lower, ctx);
}

/// A line's range is its domain: the convex quadratic needs no window.
[[nodiscard]] Range WholeLineRange(const Interval& domain)
{
    return {domain.GetLower(), domain.GetUpper(), std::isfinite(domain.GetLower()),
            std::isfinite(domain.GetUpper()), false};
}

[[nodiscard]] ClosestApproachResult TwoLines(const ApproachSide& a, const Interval& domainA,
                                             const ApproachSide& b, const Interval& domainB,
                                             const Context& ctx)
{
    const Range rs = WholeLineRange(domainA);
    const Range rt = WholeLineRange(domainB);
    if (IsZero(ExactCross(ExactOf(a.D), ExactOf(b.D)))) {
        return ParallelLines(a, rs, b, rt, ctx);
    }
    // Strictly convex: Newton from any start in range reaches the unique
    // minimiser, in one step when no end binds.
    const auto start = [](const Range& r) {
        return std::clamp(0.0, std::isfinite(r.Lo) ? r.Lo : -1e300,
                          std::isfinite(r.Hi) ? r.Hi : 1e300);
    };
    const Refined z = Refine(a, b, rs, rt, start(rs), start(rt));
    return Assemble(a, b, rs, rt, z, ApproachUniqueness::Unique, 0.0, ctx);
}

// --- A line and a conic, two conics -----------------------------------------

/// The search over ranges whose seams have been moved by the given
/// offsets.
[[nodiscard]] Search Searched(const ApproachSide& a, const Range& rs, const ApproachSide& b,
                              const Range& rt, const double shiftS, const double shiftT)
{
    Range s = rs;
    Range t = rt;
    if (s.Periodic) {
        s.Lo += shiftS;
        s.Hi += shiftS;
    }
    if (t.Periodic) {
        t.Lo += shiftT;
        t.Hi += shiftT;
    }
    Search search{a, b, s, t};
    // The size of the configuration: each conic's reach, each line's
    // searched stretch, the distance between their origins.
    const auto size = [](const ApproachSide& side, const Range& r) {
        return side.IsLine() ? side.L1 * std::max(std::abs(r.Lo), std::abs(r.Hi)) : side.Reach;
    };
    const double extent = size(a, s) + size(b, t) + Length(a.Origin - b.Origin) + 1.0;
    search.Slack = 256.0 * Epsilon * extent * extent;
    search.Run();
    return search;
}

/// A range narrowed to around a cluster's hull, its ends kept where they
/// are the range's.
[[nodiscard]] Range Around(const Range& r, const double lo, const double hi)
{
    const double margin = hi - lo;
    Range narrowed = r;
    narrowed.Lo = std::max(r.Lo, lo - margin);
    narrowed.Hi = std::min(r.Hi, hi + margin);
    narrowed.LoIsEnd = r.LoIsEnd && narrowed.Lo == r.Lo;
    narrowed.HiIsEnd = r.HiIsEnd && narrowed.Hi == r.Hi;
    narrowed.Periodic = false;
    return narrowed;
}

/// Whether the one cluster holds one minimiser: proved on its hull, or,
/// when the hull is too wide for the bounds, by searching again around it
/// with finer leaves, which shrinks it. Every point outside the hull was
/// proved worse, so the minimiser is inside whatever narrower search
/// follows.
[[nodiscard]] bool ProvedAlone(const ApproachSide& a, const ApproachSide& b, const Search& search,
                               Box hull, Refined& z)
{
    Range rs = search.RangeS;
    Range rt = search.RangeT;
    for (int pass = 0; pass < 4; ++pass) {
        if (ProvedUnique(a, b, hull, z)) {
            return true;
        }
        rs = Around(rs, hull.S0, hull.S1);
        rt = Around(rt, hull.T0, hull.T1);
        Search finer{a, b, rs, rt};
        finer.Slack = search.Slack;
        finer.Run();
        if (finer.Exhausted || finer.Leaves.empty()) {
            return false;
        }
        const std::vector<Box> hulls = Clusters(finer.Leaves, rs, rt);
        if (hulls.size() != 1) {
            return false;
        }
        hull = hulls.front();
        const Refined refined = Refine(a, b, search.RangeS, search.RangeT, z.S, z.T);
        z = refined;
    }
    return ProvedUnique(a, b, hull, z);
}

/// The search with the seams of whole conics moved off every leaf.
[[nodiscard]] Search SeamFree(const ApproachSide& a, const Range& rs, const ApproachSide& b,
                              const Range& rt)
{
    constexpr double half = std::numbers::pi;
    for (const auto& [x, y] :
         {std::pair{0.0, 0.0}, std::pair{half, 0.0}, std::pair{0.0, half}, std::pair{half, half}}) {
        Search search = Searched(a, rs, b, rt, x, y);
        if (search.Exhausted || !search.TouchesSeam()) {
            return search;
        }
    }
    Search search = Searched(a, rs, b, rt, 0.0, 0.0);
    search.Exhausted = true; // a cluster around the whole conic: no proof
    return search;
}

[[nodiscard]] ClosestApproachResult Searched(const ApproachSide& a, const Range& rs,
                                             const ApproachSide& b, const Range& rt,
                                             const bool provedSet, const Context& ctx)
{
    const Search search = SeamFree(a, rs, b, rt);
    const double lowest = search.Leaves.empty()
                              ? search.Best
                              : std::ranges::min(search.Leaves, {}, &Box::Lower).Lower;
    const double lower = std::sqrt(2.0 * std::max(lowest, 0.0));
    if (search.Exhausted || provedSet) {
        const Refined z = Refine(a, b, search.RangeS, search.RangeT, search.BestS, search.BestT);
        const ApproachUniqueness u =
            provedSet ? ApproachUniqueness::NotUnique : ApproachUniqueness::Undecided;
        return Assemble(a, b, search.RangeS, search.RangeT, z, u, provedSet ? lower : 0.0, ctx);
    }
    const std::vector<Box> hulls = Clusters(search.Leaves, search.RangeS, search.RangeT);
    const Box& first = hulls.front();
    Refined best = Refine(a, b, search.RangeS, search.RangeT, first.Cs(), first.Ct());
    if (hulls.size() == 1 && ProvedAlone(a, b, search, first, best)) {
        return Assemble(a, b, search.RangeS, search.RangeT, best, ApproachUniqueness::Unique, lower,
                        ctx);
    }
    // Several candidates the search could not separate: the best of them.
    for (const Box& hull : hulls) {
        const Refined z = Refine(a, b, search.RangeS, search.RangeT, hull.Cs(), hull.Ct());
        if (LocalAt(a, b, z.S, z.T).Phi < LocalAt(a, b, best.S, best.T).Phi) {
            best = z;
        }
    }
    return Assemble(a, b, search.RangeS, search.RangeT, best, ApproachUniqueness::Undecided, lower,
                    ctx);
}

/// Whether the minimisers are proved a set: a line on a circle's axis, or
/// two coaxial whole circles.
[[nodiscard]] bool ProvedSet(const ApproachSide& a, const Range& rs, const ApproachSide& b,
                             const Range& rt)
{
    using K = ApproachSide::Kind;
    if (a.Shape == K::Line || b.Shape == K::Line) {
        const ApproachSide& line = a.Shape == K::Line ? a : b;
        const ApproachSide& other = a.Shape == K::Line ? b : a;
        return other.Shape == K::Circle && OnAxis(line, other);
    }
    return a.Shape == K::Circle && b.Shape == K::Circle && rs.Periodic && rt.Periodic &&
           Coaxial(a, b);
}

/// The ranges of a line against a conic, or of two conics, and whether
/// their minimisers are proved a set.
[[nodiscard]] ClosestApproachResult WithConic(const ApproachSide& a, const Interval& domainA,
                                              const ApproachSide& b, const Interval& domainB,
                                              const Context& ctx)
{
    using K = ApproachSide::Kind;
    const Range rs =
        a.Shape == K::Line ? LineRange(a, domainA, b.Origin, b.Reach) : ConicRange(domainA);
    const Range rt =
        b.Shape == K::Line ? LineRange(b, domainB, a.Origin, a.Reach) : ConicRange(domainB);
    return Searched(a, rs, b, rt, ProvedSet(a, rs, b, rt), ctx);
}

[[nodiscard]] Result<void> CheckDomain(const Interval& domain)
{
    if (domain.IsEmpty() || std::isnan(domain.GetLower()) || std::isnan(domain.GetUpper())) {
        return Fail(ErrorCode::InvalidInput, "ClosestApproach: an empty domain");
    }
    return {};
}

} // namespace

Result<ClosestApproachResult> ClosestApproach(const Curve& a, const Interval& domainA,
                                              const Curve& b, const Interval& domainB,
                                              const Context& ctx)
{
    STEREON_TRY(const ApproachSide sideA, Detail::SideOf(a, "ClosestApproach"));
    STEREON_TRY(const ApproachSide sideB, Detail::SideOf(b, "ClosestApproach"));
    STEREON_TRY_VOID(CheckDomain(domainA));
    STEREON_TRY_VOID(CheckDomain(domainB));
    if (sideA.Shape == ApproachSide::Kind::Line && sideB.Shape == ApproachSide::Kind::Line) {
        return TwoLines(sideA, domainA, sideB, domainB, ctx);
    }
    return WithConic(sideA, domainA, sideB, domainB, ctx);
}

Result<ClosestApproachResult> ClosestApproach(const Curve& a, const Curve& b, const Context& ctx)
{
    return ClosestApproach(a, Interval::Everything(), b, Interval::Everything(), ctx);
}

} // namespace Stereon
