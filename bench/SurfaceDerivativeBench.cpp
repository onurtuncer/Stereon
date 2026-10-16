// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Where the time of a B-spline surface derivative goes, stage by stage:
// the span searches, the basis triangles of every degree (the derivatives
// come from them and differenced poles), the tensor-product contraction
// (every partial up to (2, 2) against those of total order at most 2), and
// the curvature arithmetic on top. The same bicubic 12 x 12 surface and
// 4096 parameter pairs as the OCCT harness (bench/occt/OcctBench.cpp), so
// the stages add up to its B1.6 surface row.

#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

import Stereon.Surface;

namespace {

using Stereon::BSplineSurface;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec3;

constexpr int ParameterCount = 4096;

/// A clamped bicubic surface over a 12 x 12 net of a wavy sheet, as in the
/// OCCT harness.
BSplineSurface SomeSurface(const bool rational)
{
    constexpr int count = 12;
    constexpr int degree = 3;
    std::vector<double> knots(degree + 1, 0.0);
    for (int i = 1; i < count - degree; ++i) {
        knots.push_back(static_cast<double>(i));
    }
    knots.insert(knots.end(), degree + 1, static_cast<double>(count - degree));
    BSplineSurface s{degree, degree, knots, knots, {}, {}};
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            s.Poles.push_back({2.0 * i, 2.0 * j, std::sin(0.7 * i) * std::cos(0.5 * j) * 3.0});
            if (rational) {
                s.Weights.push_back(1.0 + 0.3 * ((i + j) % 3));
            }
        }
    }
    return s;
}

std::vector<Point2> SurfaceParameters(const BSplineSurface& s)
{
    const Stereon::Box2 d = Domain(s);
    std::vector<Point2> uv;
    for (int i = 0; i < ParameterCount; ++i) {
        const double a = (i % 64) / 64.0;
        const double b = (i / 64) / 64.0;
        uv.push_back({d.X.GetLower() + a * (d.X.GetUpper() - d.X.GetLower()),
                      d.Y.GetLower() + b * (d.Y.GetUpper() - d.Y.GetLower())});
    }
    return uv;
}

void Items(benchmark::State& state)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(ParameterCount));
}

/// The two span searches alone.
void Spans(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const std::vector<Point2> uv = SurfaceParameters(s);
    const std::size_t countU = Stereon::PoleCountU(s);
    const std::size_t countV = Stereon::PoleCountV(s);
    for (auto _ : state) {
        std::size_t sum = 0;
        for (const Point2& p : uv) {
            sum += Stereon::FindSpan(s.KnotsU, s.DegreeU, countU, p.X);
            sum += Stereon::FindSpan(s.KnotsV, s.DegreeV, countV, p.Y);
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

/// The span searches and the basis triangles in both directions:
/// everything before the poles are touched.
void Triangles(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const std::vector<Point2> uv = SurfaceParameters(s);
    const std::size_t countU = Stereon::PoleCountU(s);
    const std::size_t countV = Stereon::PoleCountV(s);
    const auto p = static_cast<std::size_t>(s.DegreeU);
    const auto q = static_cast<std::size_t>(s.DegreeV);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& point : uv) {
            std::array<double, 64> nu; // filled by BasisFunctionsOfAllDegrees
            std::array<double, 64> nv;
            const std::size_t spanU = Stereon::FindSpan(s.KnotsU, s.DegreeU, countU, point.X);
            const std::size_t spanV = Stereon::FindSpan(s.KnotsV, s.DegreeV, countV, point.Y);
            Stereon::BasisFunctionsOfAllDegrees(s.KnotsU, s.DegreeU, spanU, point.X,
                                                std::span{nu}.first((p + 1) * (p + 1)));
            Stereon::BasisFunctionsOfAllDegrees(s.KnotsV, s.DegreeV, spanV, point.Y,
                                                std::span{nv}.first((q + 1) * (q + 1)));
            sum += nu[p * (p + 1)] + nv[q * (q + 1)];
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

template <bool Rational>
void D0(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const std::vector<Point2> uv = SurfaceParameters(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            sum += Eval(s, p).Z;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

/// Every partial up to (du, dv).
template <bool Rational, int Du, int Dv>
void Square(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const std::vector<Point2> uv = SurfaceParameters(s);
    std::array<Vec3, (Du + 1) * (Dv + 1)> out;
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            Stereon::Derivatives(s, p, Du, Dv, std::span<Vec3>{out});
            sum += out.back().Z;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

/// Every partial of total order at most 2.
template <bool Rational>
void Total2(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const std::vector<Point2> uv = SurfaceParameters(s);
    std::array<Vec3, 9> out;
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            Stereon::TotalDerivatives(s, p, 2, std::span<Vec3>{out});
            sum += out[6].Z;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

template <bool Rational>
void Curvatures(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const std::vector<Point2> uv = SurfaceParameters(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            const auto k = Stereon::PrincipalCurvatures(s, p);
            sum += k->Min + k->Max + k->MaxDirection.X;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state);
}

} // namespace

BENCHMARK(Spans)->Name("SurfaceD2/Spans");
BENCHMARK(Triangles)->Name("SurfaceD2/SpansAndTriangles");
BENCHMARK(D0<false>)->Name("SurfaceD2/D0");
BENCHMARK(Square<false, 1, 1>)->Name("SurfaceD2/Square11");
BENCHMARK(Square<false, 2, 2>)->Name("SurfaceD2/Square22");
BENCHMARK(Total2<false>)->Name("SurfaceD2/Total2");
BENCHMARK(Curvatures<false>)->Name("SurfaceD2/Curvatures");
BENCHMARK(D0<true>)->Name("SurfaceD2/D0-rational");
BENCHMARK(Square<true, 2, 2>)->Name("SurfaceD2/Square22-rational");
BENCHMARK(Total2<true>)->Name("SurfaceD2/Total2-rational");
BENCHMARK(Curvatures<true>)->Name("SurfaceD2/Curvatures-rational");
