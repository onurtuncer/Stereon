// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The hot paths the C++26 work of docs/performance.rst, *C++26 features
// for speed*, and the inline-storage pass after it, *Inline storage on hot
// paths*, change, measured before and after: per-degree kernels in the
// basis functions, Eval and EvalGrid; the deferred exact values of the
// analytic intersection handlers; the projection's candidate lists; ray
// casts; knot insertion, Bezier extraction and degree elevation.
// Stereon only; the OCCT comparison of the same operations is B1 in
// bench/occt/OcctBench.cpp.

#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

import Stereon.Surface;
import Stereon.Intersect;
import Stereon.Boolean;

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Point3;
using Stereon::Vec3;

namespace {

constexpr std::size_t ParameterCount = 4096;

BSplineCurve SomeCurve(const int degree, const bool rational)
{
    const std::size_t poles = 24;
    const auto p = static_cast<std::size_t>(degree);
    std::vector<double> knots(poles + p + 1);
    for (std::size_t i = 0; i < knots.size(); ++i) {
        knots[i] =
            i <= p
                ? 0.0
                : (i >= poles ? 1.0 : static_cast<double>(i - p) / static_cast<double>(poles - p));
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        const auto x = static_cast<double>(i);
        points.push_back({x, 0.5 * x * x / 24.0 - 3.0, (i % 3 == 0) ? 1.0 : -0.5});
        if (rational) {
            weights.push_back(1.0 + 0.25 * static_cast<double>(i % 4));
        }
    }
    return *Stereon::MakeBSplineCurve(degree, knots, points, weights);
}

BSplineSurface SomeSurface(const bool rational)
{
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0};
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            poles.push_back({1.0 * i, 1.0 * j, 0.1 * ((3 * i + 5 * j) % 7)});
            if (rational) {
                weights.push_back(1.0 + 0.125 * ((i + j) % 4));
            }
        }
    }
    return *Stereon::MakeBSplineSurface(3, 3, knots, knots, poles, weights);
}

std::vector<double> Spread(const Stereon::Interval& d, const std::size_t n)
{
    std::vector<double> t(n);
    for (std::size_t i = 0; i < n; ++i) {
        t[i] = d.GetLower() + (d.GetUpper() - d.GetLower()) * (static_cast<double>(i) + 0.5) /
                                  static_cast<double>(n);
    }
    return t;
}

void Items(benchmark::State& state, const std::size_t perIteration)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(perIteration));
}

void BasisFunctions(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(static_cast<int>(state.range(0)), false);
    const std::vector<double> t = Spread(Domain(c), ParameterCount);
    std::array<double, 32> values{};
    const auto p = static_cast<std::size_t>(c.Degree);
    for (auto _ : state) {
        for (const double s : t) {
            const std::size_t span = Stereon::FindSpan(c.Knots, c.Degree, c.Poles.size(), s);
            Stereon::BasisFunctions(c.Knots, c.Degree, span, s, std::span{values}.first(p + 1));
            benchmark::DoNotOptimize(values.data());
        }
    }
    Items(state, t.size());
}

template <bool Rational>
void CurveEval(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(static_cast<int>(state.range(0)), Rational);
    const std::vector<double> t = Spread(Domain(c), ParameterCount);
    for (auto _ : state) {
        for (const double s : t) {
            benchmark::DoNotOptimize(Eval(c, s));
        }
    }
    Items(state, t.size());
}

template <bool Rational>
void SurfaceEval(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const Stereon::Box2 d = Domain(s);
    const std::vector<double> u = Spread(d.X, 64);
    const std::vector<double> v = Spread(d.Y, 64);
    for (auto _ : state) {
        for (const double a : u) {
            for (const double b : v) {
                benchmark::DoNotOptimize(Eval(s, Stereon::Point2{a, b}));
            }
        }
    }
    Items(state, u.size() * v.size());
}

template <bool Rational>
void SurfaceGrid(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const Stereon::Box2 d = Domain(s);
    const std::vector<double> u = Spread(d.X, 256);
    const std::vector<double> v = Spread(d.Y, 256);
    std::vector<Point3> points(u.size() * v.size());
    for (auto _ : state) {
        Stereon::EvalGrid(s, u, v, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    Items(state, points.size());
}

void CurveProjection(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    std::vector<Point3> queries;
    for (int i = 0; i < 256; ++i) {
        queries.push_back({0.1 * i, 0.03 * i - 2.0, (i % 5) * 0.4 - 1.0});
    }
    const Stereon::CurveProjector projector(c);
    for (auto _ : state) {
        for (const Point3& q : queries) {
            benchmark::DoNotOptimize(projector.Project(q));
        }
    }
    Items(state, queries.size());
}

const Stereon::Surface AboutZ{
    Stereon::Cylinder{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
const Stereon::Surface AboutZShifted{
    Stereon::Cylinder{Point3{3.0, 1.0, 0.0}, 4.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
const Stereon::Surface Tilted{
    Stereon::Plane{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.8, 0.6}}};
const Stereon::Surface Level{
    Stereon::Plane{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};

void Intersection(benchmark::State& state, const Stereon::Surface& a, const Stereon::Surface& b)
{
    const Stereon::Context ctx = Stereon::Context::Default();
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::Intersect(a, b, ctx));
    }
    Items(state, 1);
}

void PlanePlane(benchmark::State& state)
{
    Intersection(state, Tilted, Level);
}

void PlaneCylinder(benchmark::State& state)
{
    Intersection(state, Level, AboutZ);
}

void CylinderCylinderParallel(benchmark::State& state)
{
    Intersection(state, AboutZ, AboutZShifted);
}

const Stereon::Surface Skew{
    Stereon::Cylinder{Point3{2.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.6, 0.8}, true}};

void CylinderCylinderSkew(benchmark::State& state)
{
    Intersection(state, AboutZ, Skew);
}

void RayMeets(benchmark::State& state, const Stereon::Surface& surface)
{
    const Stereon::Context ctx = Stereon::Context::Default();
    const Point3 origin{0.3, 0.2, -7.0};
    const Vec3 d{0.1, 0.2, 1.0};
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::Detail::RayMeets(surface, origin, d, ctx));
    }
    Items(state, 1);
}

void RayMeetsPlane(benchmark::State& state)
{
    RayMeets(state, Level);
}

void RayMeetsCylinder(benchmark::State& state)
{
    RayMeets(state, AboutZ);
}

void CurveKnotInsertion(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, false);
    const std::vector<double> t = Spread(Domain(c), 100);
    for (auto _ : state) {
        for (const double s : t) {
            benchmark::DoNotOptimize(Stereon::InsertKnot(c, s));
        }
    }
    Items(state, t.size());
}

void SurfaceKnotInsertion(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const std::vector<double> u = Spread(Domain(s).X, 16);
    for (auto _ : state) {
        for (const double t : u) {
            benchmark::DoNotOptimize(Stereon::InsertKnot(s, Stereon::ParameterDirection::U, t));
        }
    }
    Items(state, u.size());
}

void BezierSegments(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, false);
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::ToBezierSegments(c));
    }
    Items(state, 1);
}

void DegreeElevation(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, false);
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::ElevatedDegree(c, 5));
    }
    Items(state, 1);
}

} // namespace

BENCHMARK(BasisFunctions)->Name("HotPath/BasisFunctions")->Arg(1)->Arg(2)->Arg(3)->Arg(5)->Arg(7);
BENCHMARK(CurveEval<false>)->Name("HotPath/CurveEval")->Arg(1)->Arg(2)->Arg(3)->Arg(5)->Arg(7);
BENCHMARK(CurveEval<true>)->Name("HotPath/CurveEval/rational")->Arg(3);
BENCHMARK(SurfaceEval<false>)->Name("HotPath/SurfaceEval/bicubic");
BENCHMARK(SurfaceEval<true>)->Name("HotPath/SurfaceEval/bicubic-rational");
BENCHMARK(SurfaceGrid<false>)->Name("HotPath/EvalGrid/256/bicubic");
BENCHMARK(SurfaceGrid<true>)->Name("HotPath/EvalGrid/256/bicubic-rational");
BENCHMARK(CurveProjection)->Name("HotPath/CurveProjection/rational-cubic");
BENCHMARK(PlanePlane)->Name("HotPath/Intersect/PlanePlane");
BENCHMARK(PlaneCylinder)->Name("HotPath/Intersect/PlaneCylinder");
BENCHMARK(CylinderCylinderParallel)->Name("HotPath/Intersect/CylinderCylinder/parallel");
BENCHMARK(CylinderCylinderSkew)->Name("HotPath/Intersect/CylinderCylinder/skew");
BENCHMARK(RayMeetsPlane)->Name("HotPath/RayMeets/plane");
BENCHMARK(RayMeetsCylinder)->Name("HotPath/RayMeets/cylinder");
BENCHMARK(CurveKnotInsertion)->Name("HotPath/CurveKnotInsertion");
BENCHMARK(SurfaceKnotInsertion)->Name("HotPath/SurfaceKnotInsertion/U");
BENCHMARK(BezierSegments)->Name("HotPath/BezierSegments/cubic-24");
BENCHMARK(DegreeElevation)->Name("HotPath/DegreeElevation/3to5");
