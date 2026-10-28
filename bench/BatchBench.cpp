// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Batch and grid evaluation against the point evaluators they reproduce bit
// for bit (docs/nurbs.rst, *Batch and grid evaluation*): what EvalAll's
// lanes and EvalGrid's shared bases buy over Eval in a loop, on the
// geometry of the OCCT rows B1.1 and B1.2 (bench/occt/OcctBench.cpp), which
// compare the same calls with OCCT's GeomGridEval. Every benchmark writes
// its points to an array allocated once, so the loops differ only in how
// they evaluate.

#include <benchmark/benchmark.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

import Stereon.Surface;

namespace {

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Point2;
using Stereon::Point3;

constexpr std::size_t ParameterCount = 4096;
constexpr std::size_t GridSide = 256;

/// The clamped 24-pole curve of the OCCT rows.
BSplineCurve SomeCurve(const int degree, const bool rational)
{
    constexpr int poles = 24;
    BSplineCurve c{degree, {}, {}, {}};
    c.Knots.assign(static_cast<std::size_t>(degree) + 1, 0.0);
    for (int i = 1; i < poles - degree; ++i) {
        c.Knots.push_back(static_cast<double>(i));
    }
    c.Knots.insert(c.Knots.end(), static_cast<std::size_t>(degree) + 1,
                   static_cast<double>(poles - degree));
    for (int i = 0; i < poles; ++i) {
        c.Poles.push_back({std::cos(0.5 * i) * 10.0, std::sin(0.5 * i) * 10.0, 0.8 * i});
        if (rational) {
            c.Weights.push_back(1.0 + 0.5 * std::sin(1.7 * i) * std::sin(1.7 * i));
        }
    }
    return c;
}

/// The clamped bicubic 12 x 12 surface of the OCCT rows.
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

/// n parameters spread evenly over an interval, ends included.
std::vector<double> Spread(const Stereon::Interval& d, const std::size_t n)
{
    std::vector<double> t(n);
    for (std::size_t k = 0; k < n; ++k) {
        t[k] = d.GetLower() +
               (d.GetUpper() - d.GetLower()) * static_cast<double>(k) / static_cast<double>(n - 1);
    }
    return t;
}

template <bool Rational>
void CurvePointByPoint(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, Rational);
    const std::vector<double> t = Spread(Domain(c), ParameterCount);
    std::vector<Point3> points(t.size());
    for (auto _ : state) {
        for (std::size_t k = 0; k < t.size(); ++k) {
            points[k] = Eval(c, t[k]);
        }
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(t.size()));
}

template <bool Rational>
void CurveBatch(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, Rational);
    const std::vector<double> t = Spread(Domain(c), ParameterCount);
    std::vector<Point3> points(t.size());
    for (auto _ : state) {
        Stereon::EvalAll(c, t, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(t.size()));
}

template <bool Rational>
void SurfacePointByPoint(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const Stereon::Box2 d = Domain(s);
    const std::vector<double> u = Spread(d.X, GridSide);
    const std::vector<double> v = Spread(d.Y, GridSide);
    std::vector<Point3> points(u.size() * v.size());
    for (auto _ : state) {
        for (std::size_t i = 0; i < u.size(); ++i) {
            for (std::size_t j = 0; j < v.size(); ++j) {
                points[i * v.size() + j] = Eval(s, Point2{u[i], v[j]});
            }
        }
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(points.size()));
}

template <bool Rational>
void SurfaceGrid(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const Stereon::Box2 d = Domain(s);
    const std::vector<double> u = Spread(d.X, GridSide);
    const std::vector<double> v = Spread(d.Y, GridSide);
    std::vector<Point3> points(u.size() * v.size());
    for (auto _ : state) {
        Stereon::EvalGrid(s, u, v, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(points.size()));
}

} // namespace

BENCHMARK(CurvePointByPoint<false>)->Name("Batch/CurveD0/cubic/PointByPoint");
BENCHMARK(CurveBatch<false>)->Name("Batch/CurveD0/cubic/EvalAll");
BENCHMARK(CurvePointByPoint<true>)->Name("Batch/CurveD0/cubic-rational/PointByPoint");
BENCHMARK(CurveBatch<true>)->Name("Batch/CurveD0/cubic-rational/EvalAll");
BENCHMARK(SurfacePointByPoint<false>)->Name("Batch/SurfaceGrid/256/bicubic/PointByPoint");
BENCHMARK(SurfaceGrid<false>)->Name("Batch/SurfaceGrid/256/bicubic/EvalGrid");
BENCHMARK(SurfacePointByPoint<true>)->Name("Batch/SurfaceGrid/256/bicubic-rational/PointByPoint");
BENCHMARK(SurfaceGrid<true>)->Name("Batch/SurfaceGrid/256/bicubic-rational/EvalGrid");
