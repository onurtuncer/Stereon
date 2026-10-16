// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor: this unit calls core's inline functions that carry
// contracts (Box::operator[], Interval's constructor) and parses no contract
// of its own, which GCC 16.2 cannot link without it
// (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <memory>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Surface;

import :BSplineSurface;
import :Conversion;

/// The B-spline surface operations that are not evaluators. Like the curve
/// operations they work on the homogeneous form (w P, w), in which knot
/// insertion along u acts on every column of poles alike and along v on
/// every row (The NURBS Book, A5.3).

namespace Stereon {

namespace {

// ---------------------------------------------------------------------------
// Homogeneous form

/// A pole in homogeneous coordinates: (w x, w y, w z, w).
using HPoint = std::array<double, 4>;

struct HSurface {
    int DegreeU = 1;
    int DegreeV = 1;
    std::vector<double> KnotsU;
    std::vector<double> KnotsV;
    std::size_t CountU = 0;
    std::size_t CountV = 0;
    std::vector<HPoint> Poles; // row-major, i along u
    bool Rational = false;
};

[[nodiscard]] HPoint Lift(const Point3& p, const double w)
{
    return {w * p.X, w * p.Y, w * p.Z, w};
}

[[nodiscard]] HSurface ToHomogeneous(const BSplineSurface& s)
{
    HSurface h{s.DegreeU,     s.DegreeV,     s.KnotsU, s.KnotsV,
               PoleCountU(s), PoleCountV(s), {},       IsRational(s)};
    h.Poles.reserve(s.Poles.size());
    for (std::size_t i = 0; i < s.Poles.size(); ++i) {
        h.Poles.push_back(Lift(s.Poles[i], Detail::SurfaceWeightOf(s, i)));
    }
    return h;
}

[[nodiscard]] Point3 Euclidean(const HPoint& q, const bool rational)
{
    return rational ? Point3{q[0] / q[3], q[1] / q[3], q[2] / q[3]} : Point3{q[0], q[1], q[2]};
}

[[nodiscard]] BSplineSurface FromHomogeneous(HSurface h)
{
    BSplineSurface s{h.DegreeU, h.DegreeV, std::move(h.KnotsU), std::move(h.KnotsV), {}, {}};
    s.Poles.reserve(h.Poles.size());
    for (const HPoint& q : h.Poles) {
        s.Poles.push_back(Euclidean(q, h.Rational));
        if (h.Rational) {
            s.Weights.push_back(q[3]);
        }
    }
    return s;
}

[[nodiscard]] HPoint Blend(const HPoint& a, const HPoint& b, const double alpha,
                           const bool rational)
{
    return {alpha * a[0] + (1.0 - alpha) * b[0], alpha * a[1] + (1.0 - alpha) * b[1],
            alpha * a[2] + (1.0 - alpha) * b[2],
            rational ? alpha * a[3] + (1.0 - alpha) * b[3] : 1.0};
}

/// Boehm's algorithm (A5.1) on one line of poles: the line after inserting t
/// r times into `knots`, for knot index k (the last with knots[k] <= t) and
/// multiplicity s. The surface applies it to every row or column.
[[nodiscard]] std::vector<HPoint> InsertIntoLine(const std::vector<HPoint>& line,
                                                 const std::vector<double>& knots, const int degree,
                                                 const double t, const std::size_t k,
                                                 const std::size_t s, const std::size_t r,
                                                 const bool rational)
{
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t n = line.size() - 1;
    std::vector<HPoint> q(n + 1 + r);
    for (std::size_t i = 0; i <= k - p; ++i) {
        q[i] = line[i];
    }
    for (std::size_t i = k - s; i <= n; ++i) {
        q[i + r] = line[i];
    }
    std::vector<HPoint> work(line.begin() + static_cast<std::ptrdiff_t>(k - p),
                             line.begin() + static_cast<std::ptrdiff_t>(k - s + 1));
    for (std::size_t j = 1; j <= r; ++j) {
        const std::size_t l = k - p + j;
        for (std::size_t i = 0; i + j + s <= p; ++i) {
            const double alpha = (t - knots[l + i]) / (knots[i + k + 1] - knots[l + i]);
            work[i] = Blend(work[i + 1], work[i], alpha, rational);
        }
        q[l] = work[0];
        q[k + r - j - s] = work[p - j - s];
    }
    const std::size_t l = k - p + r;
    for (std::size_t i = l + 1; i + s < k; ++i) {
        q[i] = work[i - l];
    }
    return q;
}

/// The line of poles i (along u, at v index j fixed) or j (along v).
[[nodiscard]] std::vector<HPoint> LineOf(const HSurface& h, const ParameterDirection direction,
                                         const std::size_t fixed)
{
    std::vector<HPoint> line;
    if (direction == ParameterDirection::U) {
        for (std::size_t i = 0; i < h.CountU; ++i) {
            line.push_back(h.Poles[i * h.CountV + fixed]);
        }
    } else {
        for (std::size_t j = 0; j < h.CountV; ++j) {
            line.push_back(h.Poles[fixed * h.CountV + j]);
        }
    }
    return line;
}

/// Inserts t `times` times into the knots of one direction; times plus the
/// current multiplicity is at most the degree.
void InsertInPlace(HSurface& h, const ParameterDirection direction, const double t,
                   const std::size_t times)
{
    const bool alongU = direction == ParameterDirection::U;
    std::vector<double>& knots = alongU ? h.KnotsU : h.KnotsV;
    const int degree = alongU ? h.DegreeU : h.DegreeV;
    const auto k = static_cast<std::size_t>(std::ranges::upper_bound(knots, t) - knots.begin()) - 1;
    const std::size_t s = KnotMultiplicity(knots, t);
    const std::size_t lines = alongU ? h.CountV : h.CountU;
    const std::size_t countU = alongU ? h.CountU + times : h.CountU;
    const std::size_t countV = alongU ? h.CountV : h.CountV + times;
    std::vector<HPoint> poles(countU * countV);
    for (std::size_t fixed = 0; fixed < lines; ++fixed) {
        const std::vector<HPoint> line =
            InsertIntoLine(LineOf(h, direction, fixed), knots, degree, t, k, s, times, h.Rational);
        for (std::size_t index = 0; index < line.size(); ++index) {
            poles[alongU ? index * countV + fixed : fixed * countV + index] = line[index];
        }
    }
    knots.insert(knots.begin() + static_cast<std::ptrdiff_t>(k + 1), times, t);
    h.Poles = std::move(poles);
    h.CountU = countU;
    h.CountV = countV;
}

/// Raises every knot value of one direction's domain to the degree, so that
/// every pair of spans bounds a Bezier patch.
void RaiseToDegree(HSurface& h, const ParameterDirection direction, const Interval domain)
{
    const bool alongU = direction == ParameterDirection::U;
    std::vector<double> values;
    for (const double u : alongU ? h.KnotsU : h.KnotsV) {
        if (Contains(domain, u) && (values.empty() || values.back() != u)) {
            values.push_back(u);
        }
    }
    const auto degree = static_cast<std::size_t>(alongU ? h.DegreeU : h.DegreeV);
    for (const double u : values) {
        const std::size_t s = KnotMultiplicity(alongU ? h.KnotsU : h.KnotsV, u);
        if (s < degree) {
            InsertInPlace(h, direction, u, degree - s);
        }
    }
}

/// One Bezier patch in homogeneous form: (p + 1)(q + 1) poles, row-major,
/// over [U0, U1] x [V0, V1].
struct Patch {
    Box2 Range;
    std::vector<HPoint> Poles;
};

/// The non-empty spans [knots[i], knots[i+1]] inside the domain, as the
/// index i of their first knot.
[[nodiscard]] std::vector<std::size_t> SpansOf(const std::vector<double>& knots,
                                               const std::size_t degree, const std::size_t count,
                                               const Interval domain)
{
    std::vector<std::size_t> spans;
    for (std::size_t i = degree; i < count; ++i) {
        if (knots[i] < knots[i + 1] && knots[i] >= domain.GetLower() &&
            knots[i + 1] <= domain.GetUpper()) {
            spans.push_back(i);
        }
    }
    return spans;
}

[[nodiscard]] std::vector<Patch> BezierPatchesOf(const BSplineSurface& s)
{
    const Box2 domain = Domain(s);
    HSurface h = ToHomogeneous(s);
    RaiseToDegree(h, ParameterDirection::U, domain.X);
    RaiseToDegree(h, ParameterDirection::V, domain.Y);
    const auto p = static_cast<std::size_t>(h.DegreeU);
    const auto q = static_cast<std::size_t>(h.DegreeV);
    std::vector<Patch> patches;
    for (const std::size_t i : SpansOf(h.KnotsU, p, h.CountU, domain.X)) {
        for (const std::size_t j : SpansOf(h.KnotsV, q, h.CountV, domain.Y)) {
            Patch patch{
                {Interval{h.KnotsU[i], h.KnotsU[i + 1]}, Interval{h.KnotsV[j], h.KnotsV[j + 1]}},
                {}};
            for (std::size_t r = 0; r <= p; ++r) {
                for (std::size_t c = 0; c <= q; ++c) {
                    patch.Poles.push_back(h.Poles[(i - p + r) * h.CountV + (j - q + c)]);
                }
            }
            patches.push_back(std::move(patch));
        }
    }
    return patches;
}

/// Keeps the lines [first, first + count) of the poles in one direction:
/// rows (fixed u index) for U, columns for V.
void KeepLines(HSurface& h, const ParameterDirection direction, const std::size_t first,
               const std::size_t count)
{
    const bool alongU = direction == ParameterDirection::U;
    const std::size_t countU = alongU ? count : h.CountU;
    const std::size_t countV = alongU ? h.CountV : count;
    std::vector<HPoint> poles;
    poles.reserve(countU * countV);
    for (std::size_t i = 0; i < countU; ++i) {
        for (std::size_t j = 0; j < countV; ++j) {
            const std::size_t from = alongU ? (first + i) * h.CountV + j : i * h.CountV + first + j;
            poles.push_back(h.Poles[from]);
        }
    }
    h.Poles = std::move(poles);
    h.CountU = countU;
    h.CountV = countV;
}

/// Cuts one direction down to [low, high], inside its domain, as curve
/// splitting does: each cut value is raised to multiplicity p, after which
/// the poles and knots on its far side are dropped and the value is
/// repeated once more to clamp the new end.
void RestrictDirection(HSurface& h, const ParameterDirection direction, const double low,
                       const double high)
{
    const bool alongU = direction == ParameterDirection::U;
    const auto p = static_cast<std::size_t>(alongU ? h.DegreeU : h.DegreeV);
    std::vector<double>* knots = alongU ? &h.KnotsU : &h.KnotsV;
    const auto raise = [&](const double value) {
        const std::size_t s = KnotMultiplicity(*knots, value);
        if (s < p) {
            InsertInPlace(h, direction, value, p - s);
        }
        return static_cast<std::size_t>(std::ranges::upper_bound(*knots, value) - knots->begin()) -
               1;
    };
    if (low > (*knots)[p]) {
        const std::size_t k = raise(low);
        const std::size_t lines = alongU ? h.CountU : h.CountV;
        KeepLines(h, direction, k - p, lines - (k - p));
        knots->erase(knots->begin(), knots->begin() + static_cast<std::ptrdiff_t>(k - p + 1));
        knots->insert(knots->begin(), low);
    }
    const std::size_t count = alongU ? h.CountU : h.CountV;
    if (high < (*knots)[count]) {
        const std::size_t k = raise(high);
        KeepLines(h, direction, 0, k - p + 1);
        knots->erase(knots->begin() + static_cast<std::ptrdiff_t>(k + 1), knots->end());
        knots->push_back(high);
    }
}

// ---------------------------------------------------------------------------
// Bivariate Bernstein polynomials on [0, 1]^2, for the projection

/// sum_ij C[i (N + 1) + j] B_{i,M}(s) B_{j,N}(t).
struct Bernstein2 {
    std::size_t M = 0;
    std::size_t N = 0;
    std::vector<double> C;

    [[nodiscard]] double& At(const std::size_t i, const std::size_t j)
    {
        return C[i * (N + 1) + j];
    }
    [[nodiscard]] double At(const std::size_t i, const std::size_t j) const
    {
        return C[i * (N + 1) + j];
    }
};

[[nodiscard]] Bernstein2 Zero2(const std::size_t m, const std::size_t n)
{
    return {m, n, std::vector<double>((m + 1) * (n + 1), 0.0)};
}

/// The product in the tensor-product Bernstein basis of degree
/// (Ma + Mb, Na + Nb), the bivariate form of the univariate product
/// formula, applied in each parameter.
/// C(n, 0), ..., C(n, n).
[[nodiscard]] std::vector<double> BinomialRow(const std::size_t n)
{
    std::vector<double> row(n + 1, 1.0);
    for (std::size_t k = 1; k <= n; ++k) {
        row[k] = row[k - 1] * static_cast<double>(n - k + 1) / static_cast<double>(k);
    }
    return row;
}

/// The factors C(m, i) C(n, j) / C(m + n, i + j) of the product formula in
/// one parameter, indexed i * (n + 1) + j.
[[nodiscard]] std::vector<double> ProductFactors(const std::size_t m, const std::size_t n)
{
    const std::vector<double> rm = BinomialRow(m);
    const std::vector<double> rn = BinomialRow(n);
    const std::vector<double> rmn = BinomialRow(m + n);
    std::vector<double> factors((m + 1) * (n + 1));
    for (std::size_t i = 0; i <= m; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            factors[i * (n + 1) + j] = rm[i] * rn[j] / rmn[i + j];
        }
    }
    return factors;
}

[[nodiscard]] Bernstein2 Product(const Bernstein2& a, const Bernstein2& b)
{
    Bernstein2 c = Zero2(a.M + b.M, a.N + b.N);
    const std::vector<double> fu = ProductFactors(a.M, b.M);
    const std::vector<double> fv = ProductFactors(a.N, b.N);
    for (std::size_t i = 0; i <= a.M; ++i) {
        for (std::size_t k = 0; k <= b.M; ++k) {
            const double cu = fu[i * (b.M + 1) + k];
            for (std::size_t j = 0; j <= a.N; ++j) {
                const double aij = cu * a.At(i, j);
                for (std::size_t l = 0; l <= b.N; ++l) {
                    c.At(i + k, j + l) += aij * fv[j * (b.N + 1) + l] * b.At(k, l);
                }
            }
        }
    }
    return c;
}

/// The partial derivative in s (alongU) or t, one degree lower in it.
[[nodiscard]] Bernstein2 Derivative(const Bernstein2& a, const bool alongU)
{
    if (alongU) {
        Bernstein2 d = Zero2(a.M - 1, a.N);
        for (std::size_t i = 0; i < a.M; ++i) {
            for (std::size_t j = 0; j <= a.N; ++j) {
                d.At(i, j) = static_cast<double>(a.M) * (a.At(i + 1, j) - a.At(i, j));
            }
        }
        return d;
    }
    Bernstein2 d = Zero2(a.M, a.N - 1);
    for (std::size_t i = 0; i <= a.M; ++i) {
        for (std::size_t j = 0; j < a.N; ++j) {
            d.At(i, j) = static_cast<double>(a.N) * (a.At(i, j + 1) - a.At(i, j));
        }
    }
    return d;
}

void Accumulate(Bernstein2& into, const Bernstein2& term, const double factor)
{
    for (std::size_t k = 0; k < into.C.size(); ++k) {
        into.C[k] += factor * term.C[k];
    }
}

/// The two halves in one parameter at 1/2, by de Casteljau's algorithm on
/// every row (alongU) or column, each rescaled to [0, 1].
[[nodiscard]] std::pair<Bernstein2, Bernstein2> Halves(const Bernstein2& a, const bool alongU)
{
    std::pair<Bernstein2, Bernstein2> halves{a, a};
    const std::size_t degree = alongU ? a.M : a.N;
    const std::size_t lines = alongU ? a.N + 1 : a.M + 1;
    std::vector<double> work(degree + 1);
    for (std::size_t line = 0; line < lines; ++line) {
        const auto at = [&](Bernstein2& b, const std::size_t k) -> double& {
            return alongU ? b.At(k, line) : b.At(line, k);
        };
        for (std::size_t k = 0; k <= degree; ++k) {
            work[k] = alongU ? a.At(k, line) : a.At(line, k);
        }
        at(halves.first, 0) = work[0];
        at(halves.second, degree) = work[degree];
        for (std::size_t r = 1; r <= degree; ++r) {
            for (std::size_t k = 0; k + r <= degree; ++k) {
                work[k] = 0.5 * (work[k] + work[k + 1]);
            }
            at(halves.first, r) = work[0];
            at(halves.second, degree - r) = work[degree - r];
        }
    }
    return halves;
}

/// True when every coefficient has the same strict sign: the polynomial then
/// has no zero on the patch (convex hull property of the Bernstein form).
[[nodiscard]] bool OneSigned(const Bernstein2& a)
{
    const bool positive = std::ranges::all_of(a.C, [](const double c) { return c > 0.0; });
    return positive || std::ranges::all_of(a.C, [](const double c) { return c < 0.0; });
}

/// A (sub-)patch in the projection's search: its homogeneous coordinate
/// polynomials (x, y, z, w), the two stationarity polynomials F and G, and
/// its parameter range on the surface.
struct SearchPatch {
    std::array<Bernstein2, 4> Coordinates;
    Bernstein2 F;
    Bernstein2 G;
    Box2 Range;
};

/// The point-independent parts of the stationarity polynomials. With the
/// homogeneous patch A / w and H = A_s w - A w_s,
/// F = (A - w p) . H = sum_d A_d H_d - sum_d p_d (w H_d), and G likewise in
/// t: up to the positive factor 2 / w^3 the partial derivatives of the
/// squared distance |A / w - p|^2, so the interior stationary points are
/// their common zeros. The point enters linearly, so a projector computes
/// the constant part and the three linear parts once per patch.
struct StationarityParts {
    Bernstein2 Constant;
    std::array<Bernstein2, 3> Linear;
};

[[nodiscard]] StationarityParts PartsOf(const std::array<Bernstein2, 4>& a, const bool alongU)
{
    const Bernstein2& w = a[3];
    const Bernstein2 dw = Derivative(w, alongU);
    StationarityParts parts;
    for (std::size_t d = 0; d < 3; ++d) {
        Bernstein2 velocity = Product(Derivative(a[d], alongU), w);
        Accumulate(velocity, Product(a[d], dw), -1.0);
        const Bernstein2 term = Product(a[d], velocity);
        if (d == 0) {
            parts.Constant = term;
        } else {
            Accumulate(parts.Constant, term, 1.0);
        }
        parts.Linear[d] = Product(w, velocity);
    }
    return parts;
}

/// The stationarity polynomial for the point p from its parts.
[[nodiscard]] Bernstein2 For(const StationarityParts& parts, const Point3& p)
{
    Bernstein2 result = parts.Constant;
    for (std::size_t d = 0; d < 3; ++d) {
        Accumulate(result, parts.Linear[d], -p[d]);
    }
    return result;
}

/// The homogeneous coordinate polynomials (x, y, z, w) of a Bezier patch.
[[nodiscard]] std::array<Bernstein2, 4> CoordinatesOf(const Patch& patch, const std::size_t p,
                                                      const std::size_t q)
{
    std::array<Bernstein2, 4> a{Zero2(p, q), Zero2(p, q), Zero2(p, q), Zero2(p, q)};
    for (std::size_t i = 0; i <= p; ++i) {
        for (std::size_t j = 0; j <= q; ++j) {
            for (std::size_t d = 0; d < 4; ++d) {
                a[d].At(i, j) = patch.Poles[i * (q + 1) + j][d];
            }
        }
    }
    return a;
}

/// The squared distance from p to the box of the (sub-)patch's poles, a
/// lower bound on its squared distance to the (sub-)patch.
[[nodiscard]] double LowerBoundSquared(const std::array<Bernstein2, 4>& a, const Point3& p)
{
    double total = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        double low = std::numeric_limits<double>::infinity();
        double high = -low;
        for (std::size_t k = 0; k < a[d].C.size(); ++k) {
            low = std::min(low, a[d].C[k] / a[3].C[k]);
            high = std::max(high, a[d].C[k] / a[3].C[k]);
        }
        const double gap = std::max({low - p[d], p[d] - high, 0.0});
        total += gap * gap;
    }
    return total;
}

/// The four quarters of a search patch, in a fixed order.
[[nodiscard]] std::array<SearchPatch, 4> Quarters(const SearchPatch& patch)
{
    const double um = Midpoint(patch.Range.X);
    const double vm = Midpoint(patch.Range.Y);
    const Interval lowU{patch.Range.X.GetLower(), um};
    const Interval highU{um, patch.Range.X.GetUpper()};
    const Interval lowV{patch.Range.Y.GetLower(), vm};
    const Interval highV{vm, patch.Range.Y.GetUpper()};
    std::array<SearchPatch, 4> quarters{
        SearchPatch{{}, {}, {}, Box2{lowU, lowV}}, SearchPatch{{}, {}, {}, Box2{highU, lowV}},
        SearchPatch{{}, {}, {}, Box2{lowU, highV}}, SearchPatch{{}, {}, {}, Box2{highU, highV}}};
    const auto split = [&quarters](const Bernstein2& a, const auto member) {
        const auto [left, right] = Halves(a, true);
        const auto [ll, lh] = Halves(left, false);
        const auto [rl, rh] = Halves(right, false);
        member(quarters[0]) = ll;
        member(quarters[1]) = rl;
        member(quarters[2]) = lh;
        member(quarters[3]) = rh;
    };
    for (std::size_t d = 0; d < 4; ++d) {
        split(patch.Coordinates[d],
              [d](SearchPatch& s) -> Bernstein2& { return s.Coordinates[d]; });
    }
    split(patch.F, [](SearchPatch& s) -> Bernstein2& { return s.F; });
    split(patch.G, [](SearchPatch& s) -> Bernstein2& { return s.G; });
    return quarters;
}

// ---------------------------------------------------------------------------
// Projection

/// The depth of sub-patch subdivision before Newton's method takes over: a
/// leaf is 1/64 of its Bezier patch in each direction. Leaves are reached
/// only where both F and G change sign and the pole box is not farther than
/// the best point, which in practice means near a stationary point.
constexpr int LeafDepth = 6;

/// Newton steps from a leaf's centre. Convergence near a regular
/// stationary point is quadratic, a handful of steps from inside its basin;
/// the bound keeps the work finite from starts outside it.
constexpr int NewtonSteps = 32;

struct Closest {
    Point2 Parameter;
    double DistanceSquared;
};

void Consider(const BSplineSurface& s, const Point3& p, const Point2& uv, Closest& best)
{
    const double d = DistanceSquared(Eval(s, uv), p);
    const bool earlier =
        uv.X < best.Parameter.X || (uv.X == best.Parameter.X && uv.Y < best.Parameter.Y);
    if (d < best.DistanceSquared || (d == best.DistanceSquared && earlier)) {
        best = {uv, d};
    }
}

[[nodiscard]] Point2 Clamped(const Point2& uv, const Box2& box)
{
    return {std::clamp(uv.X, box.X.GetLower(), box.X.GetUpper()),
            std::clamp(uv.Y, box.Y.GetLower(), box.Y.GetUpper())};
}

/// Newton's method on (S - p) . S_u = 0, (S - p) . S_v = 0 from uv, kept
/// inside `box`: the Jacobian is the Hessian of half the squared distance,
/// [[S_u.S_u + r.S_uu, S_u.S_v + r.S_uv], [., S_v.S_v + r.S_vv]]. Each step
/// takes every partial it needs from one Derivatives call.
[[nodiscard]] Point2 Newton(const BSplineSurface& s, const Point3& p, Point2 uv, const Box2& box)
{
    std::array<Vec3, 9> d; // d[k * 3 + l] = S^(k,l)
    for (int step = 0; step < NewtonSteps; ++step) {
        Derivatives(s, uv, 2, 2, std::span<Vec3>{d});
        const Vec3 r = d[0] - ToVec(p);
        const Vec3& su = d[3];
        const Vec3& sv = d[1];
        const double a = Dot(su, su) + Dot(r, d[6]);
        const double b = Dot(su, sv) + Dot(r, d[4]);
        const double c = Dot(sv, sv) + Dot(r, d[2]);
        const double f = Dot(r, su);
        const double g = Dot(r, sv);
        const double det = a * c - b * b;
        if (!(det != 0.0) || !std::isfinite(det)) {
            break;
        }
        const Point2 next =
            Clamped(Point2{uv.X + (b * g - c * f) / det, uv.Y + (b * f - a * g) / det}, box);
        if (next == uv) {
            break;
        }
        uv = next;
    }
    return uv;
}

/// The quarters of a sub-patch in increasing order of their distance bound,
/// ties in their fixed order, so the closest region is searched first and
/// its result prunes the others.
[[nodiscard]] std::array<std::size_t, 4> ClosestFirst(const std::array<SearchPatch, 4>& quarters,
                                                      const Point3& p)
{
    std::array<std::pair<double, std::size_t>, 4> order;
    for (std::size_t k = 0; k < 4; ++k) {
        order[k] = {LowerBoundSquared(quarters[k].Coordinates, p), k};
    }
    std::ranges::sort(order);
    return {order[0].second, order[1].second, order[2].second, order[3].second};
}

void Search(const BSplineSurface& s, const Point3& p, const SearchPatch& patch, const int depth,
            Closest& best)
{
    if (LowerBoundSquared(patch.Coordinates, p) > best.DistanceSquared || OneSigned(patch.F) ||
        OneSigned(patch.G)) {
        return;
    }
    if (depth == LeafDepth) {
        const Point2 centre{Midpoint(patch.Range.X), Midpoint(patch.Range.Y)};
        Consider(s, p, Newton(s, p, centre, patch.Range), best);
        return;
    }
    const std::array<SearchPatch, 4> quarters = Quarters(patch);
    for (const std::size_t k : ClosestFirst(quarters, p)) {
        Search(s, p, quarters[k], depth + 1, best);
    }
}

/// The Bezier curve of one boundary of a patch: a row (fixed u index) or
/// column of its poles, over the patch's range in the other parameter.
[[nodiscard]] BSplineCurve EdgeCurve(const Patch& patch, const std::size_t p, const std::size_t q,
                                     const bool alongU, const std::size_t fixed,
                                     const bool rational)
{
    const std::size_t count = alongU ? p + 1 : q + 1;
    const Interval range = alongU ? patch.Range.X : patch.Range.Y;
    BSplineCurve edge{static_cast<int>(count) - 1, {}, {}, {}};
    edge.Knots.assign(count, range.GetLower());
    edge.Knots.insert(edge.Knots.end(), count, range.GetUpper());
    for (std::size_t k = 0; k < count; ++k) {
        const HPoint& h = patch.Poles[alongU ? k * (q + 1) + fixed : fixed * (q + 1) + k];
        edge.Poles.push_back(Euclidean(h, rational));
        if (rational) {
            edge.Weights.push_back(h[3]);
        }
    }
    return edge;
}

// ---------------------------------------------------------------------------
// Checked construction and the other operations

[[nodiscard]] std::size_t CountFor(const std::vector<double>& knots, const int degree)
{
    const auto needed = static_cast<std::size_t>(degree) + 1;
    return knots.size() > needed ? knots.size() - needed : 0;
}

} // namespace

Result<BSplineSurface> MakeBSplineSurface(const int degreeU, const int degreeV,
                                          std::vector<double> knotsU, std::vector<double> knotsV,
                                          std::vector<Point3> poles, std::vector<double> weights)
{
    constexpr std::string_view where = "MakeBSplineSurface";
    const std::size_t countU = CountFor(knotsU, degreeU);
    const std::size_t countV = CountFor(knotsV, degreeV);
    STEREON_TRY_VOID(Detail::CheckKnotVector(degreeU, knotsU, countU, where));
    STEREON_TRY_VOID(Detail::CheckKnotVector(degreeV, knotsV, countV, where));
    if (poles.size() != countU * countV) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeBSplineSurface: the poles must form a grid of PoleCountU x PoleCountV");
    }
    for (const Point3& pole : poles) {
        STEREON_TRY_VOID(Detail::CheckFinite(pole, "every pole", where));
    }
    STEREON_TRY_VOID(Detail::CheckWeights(weights, poles.size(), where));
    return BSplineSurface{degreeU,           degreeV,          std::move(knotsU),
                          std::move(knotsV), std::move(poles), std::move(weights)};
}

Box3 BoundingBox(const BSplineSurface& surface)
{
    Box3 box;
    double magnitude = 0.0;
    for (const Point3& pole : surface.Poles) {
        box = Hull(box, pole);
        magnitude = std::max({magnitude, std::abs(pole.X), std::abs(pole.Y), std::abs(pole.Z)});
    }
    // As Detail::PoleBoxMargin for curves, with the p + 1 and q + 1 terms
    // of both directions: 8 (p + q + 2)^2 units of 2^-53.
    const auto terms = static_cast<double>(surface.DegreeU + surface.DegreeV + 2);
    return Inflated(box, 8.0 * terms * terms * 0x1p-53 * magnitude);
}

} // namespace Stereon

namespace Stereon::Detail {

/// One Bezier patch prepared for projection: its range, coordinate
/// polynomials, the parts of F and G, and projectors for its four edges
/// (v = V0, v = V1, u = U0, u = U1).
struct PreparedPatch {
    Box2 Range;
    std::array<Bernstein2, 4> Coordinates;
    StationarityParts F;
    StationarityParts G;
    std::vector<CurveProjector> Edges;
};

struct SurfaceProjectorData {
    BSplineSurface Surface;
    std::vector<PreparedPatch> Patches;
};

} // namespace Stereon::Detail

namespace Stereon {

namespace {

[[nodiscard]] Detail::PreparedPatch Prepared(const Patch& patch, const BSplineSurface& surface)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const bool rational = IsRational(surface);
    std::array<Bernstein2, 4> a = CoordinatesOf(patch, p, q);
    Detail::PreparedPatch prepared{patch.Range, a, PartsOf(a, true), PartsOf(a, false), {}};
    prepared.Edges.emplace_back(EdgeCurve(patch, p, q, true, 0, rational));
    prepared.Edges.emplace_back(EdgeCurve(patch, p, q, true, q, rational));
    prepared.Edges.emplace_back(EdgeCurve(patch, p, q, false, 0, rational));
    prepared.Edges.emplace_back(EdgeCurve(patch, p, q, false, p, rational));
    return prepared;
}

/// The closest points of the patch's four boundary curves, which the
/// search for interior stationary points does not cover.
void ConsiderEdges(const BSplineSurface& s, const Point3& point, const Detail::PreparedPatch& patch,
                   Closest& best)
{
    const Box2& r = patch.Range;
    Consider(s, point, Point2{patch.Edges[0].Project(point), r.Y.GetLower()}, best);
    Consider(s, point, Point2{patch.Edges[1].Project(point), r.Y.GetUpper()}, best);
    Consider(s, point, Point2{r.X.GetLower(), patch.Edges[2].Project(point)}, best);
    Consider(s, point, Point2{r.X.GetUpper(), patch.Edges[3].Project(point)}, best);
}

[[nodiscard]] std::vector<std::pair<double, std::size_t>>
PatchesByBound(const std::vector<Detail::PreparedPatch>& patches, const Point3& p)
{
    std::vector<std::pair<double, std::size_t>> order;
    order.reserve(patches.size());
    for (std::size_t i = 0; i < patches.size(); ++i) {
        order.emplace_back(LowerBoundSquared(patches[i].Coordinates, p), i);
    }
    std::ranges::sort(order);
    return order;
}

/// The edges and the interior of one prepared patch, for the point p.
void SearchPatchFor(const BSplineSurface& surface, const Detail::PreparedPatch& patch,
                    const Point3& p, Closest& best)
{
    ConsiderEdges(surface, p, patch, best);
    const SearchPatch search{patch.Coordinates, For(patch.F, p), For(patch.G, p), patch.Range};
    Search(surface, p, search, 0, best);
}

/// The squared distance from p to the box of a raw patch's poles.
[[nodiscard]] double PatchBound(const Patch& patch, const Point3& p)
{
    double total = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        double low = std::numeric_limits<double>::infinity();
        double high = -low;
        for (const HPoint& q : patch.Poles) {
            low = std::min(low, q[d] / q[3]);
            high = std::max(high, q[d] / q[3]);
        }
        const double gap = std::max({low - p[d], p[d] - high, 0.0});
        total += gap * gap;
    }
    return total;
}

} // namespace

SurfaceProjector::SurfaceProjector(const BSplineSurface& surface)
{
    auto data = std::make_shared<Detail::SurfaceProjectorData>();
    for (const Patch& patch : BezierPatchesOf(surface)) {
        data->Patches.push_back(Prepared(patch, surface));
    }
    data->Surface = surface;
    m_Data = std::move(data);
}

Point2 SurfaceProjector::Project(const Point3& p) const
{
    Detail::RequireFinite(p, "Project(BSplineSurface)", "a non-finite point");
    const BSplineSurface& surface = m_Data->Surface;
    const Box2 domain = Domain(surface);
    const Point2 corner{domain.X.GetLower(), domain.Y.GetLower()};
    Closest best{corner, DistanceSquared(Eval(surface, corner), p)};
    for (const auto& [bound, index] : PatchesByBound(m_Data->Patches, p)) {
        if (bound > best.DistanceSquared) {
            break; // every later patch is at least as far
        }
        SearchPatchFor(surface, m_Data->Patches[index], p, best);
    }
    return best.Parameter;
}

Point2 Project(const BSplineSurface& surface, const Point3& p)
{
    // One point: prepare only the patches the best-first search reaches,
    // rather than all of them as a SurfaceProjector does for many points.
    Detail::RequireFinite(p, "Project(BSplineSurface)", "a non-finite point");
    const std::vector<Patch> patches = BezierPatchesOf(surface);
    std::vector<std::pair<double, std::size_t>> order;
    order.reserve(patches.size());
    for (std::size_t i = 0; i < patches.size(); ++i) {
        order.emplace_back(PatchBound(patches[i], p), i);
    }
    std::ranges::sort(order);
    const Box2 domain = Domain(surface);
    const Point2 corner{domain.X.GetLower(), domain.Y.GetLower()};
    Closest best{corner, DistanceSquared(Eval(surface, corner), p)};
    for (const auto& [bound, index] : order) {
        if (bound > best.DistanceSquared) {
            break;
        }
        SearchPatchFor(surface, Prepared(patches[index], surface), p, best);
    }
    return best.Parameter;
}

Result<BSplineSurface> InsertKnot(const BSplineSurface& surface, const ParameterDirection direction,
                                  const double value, const int times)
{
    const bool alongU = direction == ParameterDirection::U;
    const Interval domain = alongU ? Domain(surface).X : Domain(surface).Y;
    if (!std::isfinite(value) || !Contains(domain, value)) {
        return Fail(ErrorCode::InvalidInput, "InsertKnot: the knot must lie in the domain");
    }
    const std::size_t s = KnotMultiplicity(alongU ? surface.KnotsU : surface.KnotsV, value);
    const auto degree = static_cast<std::size_t>(alongU ? surface.DegreeU : surface.DegreeV);
    if (times < 1 || s + static_cast<std::size_t>(times) > degree) {
        return Fail(ErrorCode::InvalidInput,
                    "InsertKnot: times must be positive and keep the multiplicity at most the "
                    "degree");
    }
    HSurface h = ToHomogeneous(surface);
    InsertInPlace(h, direction, value, static_cast<std::size_t>(times));
    return FromHomogeneous(std::move(h));
}

std::vector<BSplineSurface> ToBezierPatches(const BSplineSurface& surface)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    std::vector<BSplineSurface> patches;
    for (const Patch& patch : BezierPatchesOf(surface)) {
        HSurface h{surface.DegreeU, surface.DegreeV, {},          {},
                   p + 1,           q + 1,           patch.Poles, IsRational(surface)};
        h.KnotsU.assign(p + 1, patch.Range.X.GetLower());
        h.KnotsU.insert(h.KnotsU.end(), p + 1, patch.Range.X.GetUpper());
        h.KnotsV.assign(q + 1, patch.Range.Y.GetLower());
        h.KnotsV.insert(h.KnotsV.end(), q + 1, patch.Range.Y.GetUpper());
        patches.push_back(FromHomogeneous(std::move(h)));
    }
    return patches;
}

Result<BSplineSurface> ToBSpline(const BSplineSurface& surface, const Box2 range)
{
    const Box2 domain = Domain(surface);
    const bool finite = std::isfinite(range.X.GetLower()) && std::isfinite(range.X.GetUpper()) &&
                        std::isfinite(range.Y.GetLower()) && std::isfinite(range.Y.GetUpper());
    const bool inside = range.X.GetLower() < range.X.GetUpper() &&
                        range.Y.GetLower() < range.Y.GetUpper() && Contains(domain, range);
    if (!finite || !inside) {
        return Fail(ErrorCode::InvalidInput,
                    "ToBSpline(BSplineSurface): the range must be finite, non-empty and inside "
                    "the domain");
    }
    if (range == domain) {
        return surface; // unchanged, without the rounding of a homogeneous round trip
    }
    HSurface h = ToHomogeneous(surface);
    RestrictDirection(h, ParameterDirection::U, range.X.GetLower(), range.X.GetUpper());
    RestrictDirection(h, ParameterDirection::V, range.Y.GetLower(), range.Y.GetUpper());
    return FromHomogeneous(std::move(h));
}

Result<BSplineCurve> IsoCurve(const BSplineSurface& surface, const ParameterDirection direction,
                              const double value)
{
    const bool fixedU = direction == ParameterDirection::U;
    const Interval domain = fixedU ? Domain(surface).X : Domain(surface).Y;
    if (!std::isfinite(value) || !Contains(domain, value)) {
        return Fail(ErrorCode::InvalidInput, "IsoCurve: the value must lie in the domain");
    }
    const HSurface h = ToHomogeneous(surface);
    const std::vector<double>& knots = fixedU ? h.KnotsU : h.KnotsV;
    const int degree = fixedU ? h.DegreeU : h.DegreeV;
    const std::size_t count = fixedU ? h.CountU : h.CountV;
    const std::size_t span = FindSpan(knots, degree, count, value);
    std::array<double, Detail::SurfaceBasisCount> basis; // filled by BasisFunctions
    const auto order = static_cast<std::size_t>(degree);
    BasisFunctions(knots, degree, span, value, std::span{basis}.first(order + 1));
    // Each pole of the curve combines one line of the surface's poles.
    const std::size_t lines = fixedU ? h.CountV : h.CountU;
    BSplineCurve curve{fixedU ? h.DegreeV : h.DegreeU, fixedU ? h.KnotsV : h.KnotsU, {}, {}};
    for (std::size_t line = 0; line < lines; ++line) {
        HPoint sum{};
        for (std::size_t r = 0; r <= order; ++r) {
            const std::size_t i = span - order + r;
            const HPoint& q = h.Poles[fixedU ? i * h.CountV + line : line * h.CountV + i];
            for (std::size_t d = 0; d < 4; ++d) {
                sum[d] += basis[r] * q[d];
            }
        }
        curve.Poles.push_back(Euclidean(sum, h.Rational));
        if (h.Rational) {
            curve.Weights.push_back(sum[3]);
        }
    }
    return curve;
}

} // namespace Stereon
