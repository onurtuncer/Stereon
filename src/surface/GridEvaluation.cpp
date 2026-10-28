// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/BSplineCurve.cpp.
#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <cstddef>
#include <span>
#include <vector>

module Stereon.Surface;

import :BSplineSurface;
import :GridEvaluation;

/// Eval(surface, uv) finds two spans, computes two basis triangles and sums
/// (p + 1)(q + 1) weighted poles. On a grid the first two depend on u alone
/// or v alone, so they are done once per row and once per column here, and
/// only the sum remains per point: in the order of Eval, with the product
/// nu[r] * nv[s] formed first and then weighted, as Eval does. Along v the
/// sum runs in the lanes of a vector when consecutive columns share a span,
/// each lane performing the scalar operations of its point.

namespace Stereon {

namespace {

using Lanes = Simd::Vec<double>;
constexpr std::size_t LaneCount = Simd::Width<double>;

/// The spans and basis functions of a list of parameters in one direction.
/// Basis[s * Count + k] is the s-th non-zero basis function at parameter k,
/// so that a run of parameters loads as one vector per s.
struct Bases {
    std::size_t Count = 0;
    std::vector<std::size_t> Spans;
    std::vector<double> Basis;
};

[[nodiscard]] Bases BasesOf(const std::span<const double> knots, const int degree,
                            const std::size_t poleCount, const std::span<const double> parameters)
{
    const auto p = static_cast<std::size_t>(degree);
    Bases bases{parameters.size(), {}, std::vector<double>((p + 1) * parameters.size())};
    bases.Spans.reserve(parameters.size());
    std::array<double, Detail::SurfaceBasisCount> values; // filled by BasisFunctions
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        const std::size_t span = FindSpan(knots, degree, poleCount, parameters[k]);
        bases.Spans.push_back(span);
        BasisFunctions(knots, degree, span, parameters[k], std::span{values}.first(p + 1));
        for (std::size_t s = 0; s <= p; ++s) {
            bases.Basis[s * bases.Count + k] = values[s];
        }
    }
    return bases;
}

/// One row of the grid: the u parameter's span and basis functions.
struct Row {
    std::size_t Span;
    std::array<double, Detail::SurfaceBasisCount> Basis;
};

/// The point at column k of the row, by Eval's sum.
[[nodiscard]] Point3 PointAt(const BSplineSurface& surface, const Row& row, const Bases& v,
                             const std::size_t k)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    const std::size_t spanV = v.Spans[k];
    Vec3 numerator{};
    double denominator = 0.0;
    for (std::size_t r = 0; r <= p; ++r) {
        for (std::size_t s = 0; s <= q; ++s) {
            const std::size_t index = (row.Span - p + r) * countV + (spanV - q + s);
            const double weighted =
                (row.Basis[r] * v.Basis[s * v.Count + k]) * Detail::SurfaceWeightOf(surface, index);
            numerator += weighted * ToVec(surface.Poles[index]);
            denominator = denominator + weighted;
        }
    }
    return IsRational(surface) ? ToPoint(numerator / denominator) : ToPoint(numerator);
}

/// The points at columns k, ..., k + LaneCount - 1 of the row, which share
/// the span `spanV`, by Eval's sum in every lane.
void PointsInLanes(const BSplineSurface& surface, const Row& row, const Bases& v,
                   const std::size_t k, const std::size_t spanV, Point3* points)
{
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    const auto q = static_cast<std::size_t>(surface.DegreeV);
    const std::size_t countV = PoleCountV(surface);
    std::array<Lanes, 3> numerator;
    numerator.fill(Simd::Broadcast(0.0));
    Lanes denominator = Simd::Broadcast(0.0);
    for (std::size_t r = 0; r <= p; ++r) {
        const Lanes nu = Simd::Broadcast(row.Basis[r]);
        for (std::size_t s = 0; s <= q; ++s) {
            const std::size_t index = (row.Span - p + r) * countV + (spanV - q + s);
            const Lanes weighted = (nu * Simd::Load(v.Basis.data() + s * v.Count + k)) *
                                   Simd::Broadcast(Detail::SurfaceWeightOf(surface, index));
            const Point3& pole = surface.Poles[index];
            numerator[0] = numerator[0] + weighted * Simd::Broadcast(pole.X);
            numerator[1] = numerator[1] + weighted * Simd::Broadcast(pole.Y);
            numerator[2] = numerator[2] + weighted * Simd::Broadcast(pole.Z);
            denominator = denominator + weighted;
        }
    }
    std::array<std::array<double, LaneCount>, 3> coordinates; // filled by Store
    for (std::size_t d = 0; d < 3; ++d) {
        const Lanes value = IsRational(surface) ? numerator[d] / denominator : numerator[d];
        Simd::Store(value, coordinates[d].data());
    }
    for (std::size_t lane = 0; lane < LaneCount; ++lane) {
        points[lane] = {coordinates[0][lane], coordinates[1][lane], coordinates[2][lane]};
    }
}

/// True when the LaneCount columns from k share a span.
[[nodiscard]] bool SharesSpan(const Bases& v, const std::size_t k)
{
    for (std::size_t lane = 1; lane < LaneCount; ++lane) {
        if (v.Spans[k + lane] != v.Spans[k]) {
            return false;
        }
    }
    return true;
}

void CheckGrid(const BSplineSurface& surface, const std::span<const double> u,
               const std::span<const double> v, const std::span<Point3> points)
{
    Detail::Require(points.size() == u.size() * v.size(), "EvalGrid(BSplineSurface)",
                    "as many points as pairs of parameters");
    const Box2 domain = Domain(surface);
    for (const double s : u) {
        Detail::CheckSurfaceEvaluation(surface, Point2{s, domain.Y.GetLower()},
                                       "EvalGrid(BSplineSurface)");
    }
    for (const double t : v) {
        Detail::CheckSurfaceEvaluation(surface, Point2{domain.X.GetLower(), t},
                                       "EvalGrid(BSplineSurface)");
    }
}

} // namespace

void EvalGrid(const BSplineSurface& surface, const std::span<const double> u,
              const std::span<const double> v, const std::span<Point3> points)
{
    CheckGrid(surface, u, v, points);
    const Bases columns = BasesOf(surface.KnotsV, surface.DegreeV, PoleCountV(surface), v);
    const auto p = static_cast<std::size_t>(surface.DegreeU);
    for (std::size_t i = 0; i < u.size(); ++i) {
        Row row{FindSpan(surface.KnotsU, surface.DegreeU, PoleCountU(surface), u[i]), {}};
        BasisFunctions(surface.KnotsU, surface.DegreeU, row.Span, u[i],
                       std::span{row.Basis}.first(p + 1));
        Point3* out = points.data() + i * v.size();
        std::size_t k = 0;
        for (; k + LaneCount <= v.size(); k += LaneCount) {
            if (SharesSpan(columns, k)) {
                PointsInLanes(surface, row, columns, k, columns.Spans[k], out + k);
                continue;
            }
            for (std::size_t lane = 0; lane < LaneCount; ++lane) {
                out[k + lane] = PointAt(surface, row, columns, k + lane);
            }
        }
        for (; k < v.size(); ++k) {
            out[k] = PointAt(surface, row, columns, k);
        }
    }
}

} // namespace Stereon
