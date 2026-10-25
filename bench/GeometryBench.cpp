// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Stereon's side of benchmark B0.2 (docs/benchmarks.md) and the ADR-0003
// verification: curve evaluation through the dispatch paths, against direct
// calls on the concrete type.
//
// - Direct: Eval(line, t) on the concrete type.
// - VisitBuiltin: one type test, then the loop runs on the concrete type.
//   ADR-0003 asks for this within 10% of direct.
// - VisitPerCall: VisitBuiltin around every single evaluation, the worst way
//   to use it.
// - DirectScalar: Eval(line, t) one evaluation at a time, kept from being
//   vectorised: the baseline a virtual call per evaluation is compared with.
// - Erased: Eval(curve, t), one virtual call per evaluation. ADR-0003 asks
//   for this to add at most 5 ns per call over a direct call.
//
// The "/1e8" runs make ADR-0003's 10^8 evaluations in one measurement.
//
// The cubic B-spline of the ADR's verification is included. The OCCT
// comparison, B0.2 proper, is in bench/occt/OcctBench.cpp, which runs the
// same curves along the same paths beside OCCT's. The copy benchmarks
// measure the wrapper at buffer sizes 64, 96 and 128 bytes against the built-in types' sizes, which
// is how ADR-0003 fixes GeometryBufferSize.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <numbers>
#include <type_traits>
#include <vector>

import Stereon.Surface;

namespace {

using Stereon::Circle;
using Stereon::Cone;
using Stereon::Curve;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Line;
using Stereon::Point3;
using Stereon::Torus;
using Stereon::Vec3;

constexpr std::size_t SampleCount = 1 << 16;

/// Iterations of SampleCount evaluations that make at least 10^8.
constexpr auto HundredMillion =
    static_cast<benchmark::IterationCount>((100'000'000 + SampleCount - 1) / SampleCount);

std::vector<double> Parameters()
{
    std::vector<double> t(SampleCount);
    for (std::size_t i = 0; i < SampleCount; ++i) {
        t[i] = static_cast<double>(i) * (2.0 * std::numbers::pi / SampleCount);
    }
    return t;
}

Line SomeLine()
{
    return {Point3{1.0, 2.0, 3.0}, Vec3{0.5, -1.0, 2.0}};
}

Circle SomeCircle()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Ellipse SomeEllipse()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

/// The cubic B-spline of ADR-0003's dispatch benchmark: eight poles over
/// [0, 2 pi], the range Parameters() covers.
Stereon::BSplineCurve SomeBSpline()
{
    constexpr double end = 2.0 * std::numbers::pi;
    return {3,
            {0.0, 0.0, 0.0, 0.0, 0.2 * end, 0.4 * end, 0.6 * end, 0.8 * end, end, end, end, end},
            {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
             Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}, Point3{5.0, 2.0, 1.0},
             Point3{6.0, -1.0, 0.0}, Point3{7.0, 0.0, 1.0}},
            {}};
}

Cylinder SomeCylinder()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Cone SomeCone()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, 0.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

Torus SomeTorus()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, 0.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
}

void SetItems(benchmark::State& state)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(SampleCount));
}

template <typename G>
double SumX(const G& geometry, const std::vector<double>& params)
{
    double total = 0.0;
    for (const double t : params) {
        total += Eval(geometry, t).X;
    }
    return total;
}

template <typename G>
void Direct(benchmark::State& state, const G geometry)
{
    const std::vector<double> params = Parameters();
    for (auto _ : state) {
        benchmark::DoNotOptimize(SumX(geometry, params));
    }
    SetItems(state);
}

template <typename G>
void DirectScalar(benchmark::State& state, const G geometry)
{
    const std::vector<double> params = Parameters();
    for (auto _ : state) {
        for (const double t : params) {
            auto p = Eval(geometry, t);
            benchmark::DoNotOptimize(p);
        }
    }
    SetItems(state);
}

void Visit(benchmark::State& state, const Curve curve)
{
    const std::vector<double> params = Parameters();
    for (auto _ : state) {
        benchmark::DoNotOptimize(
            VisitBuiltin(curve, [&params](const auto& g) { return SumX(g, params); }));
    }
    SetItems(state);
}

void VisitPerCall(benchmark::State& state, const Curve curve)
{
    const std::vector<double> params = Parameters();
    for (auto _ : state) {
        double total = 0.0;
        for (const double t : params) {
            total += VisitBuiltin(curve, [t](const auto& g) { return Eval(g, t).X; });
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state);
}

void Erased(benchmark::State& state, const Curve curve)
{
    const std::vector<double> params = Parameters();
    for (auto _ : state) {
        benchmark::DoNotOptimize(SumX(curve, params));
    }
    SetItems(state);
}

/// Copying a wrapper whose buffer is Size bytes. The model is inline when it
/// fits and shared otherwise, so the counter "shared" shows which.
template <typename G, std::size_t Size>
void CopyWrapper(benchmark::State& state)
{
    using Model = Stereon::Detail::CurveModel<3, G, Size>;
    Stereon::Detail::SmallBuffer<Stereon::Detail::CurveConcept<3>, Size> source;
    if constexpr (std::is_same_v<G, Line>) {
        source.template Emplace<Model>(SomeLine());
    } else {
        source.template Emplace<Model>(SomeCircle());
    }
    for (auto _ : state) {
        auto copy = source;
        benchmark::DoNotOptimize(copy);
    }
    state.counters["sizeof"] = static_cast<double>(sizeof(G));
    state.counters["shared"] = Stereon::Detail::FitsInline<G, Size> ? 0.0 : 1.0;
}

/// Copying a Surface-style wrapper whose buffer is Size bytes. Cone and Torus
/// take 96 bytes, so with the model's table pointer they fit only from 104.
template <typename G, std::size_t Size>
void CopySurfaceWrapper(benchmark::State& state)
{
    using Model = Stereon::Detail::SurfaceModel<G, Size>;
    Stereon::Detail::SmallBuffer<Stereon::Detail::SurfaceConcept, Size> source;
    if constexpr (std::is_same_v<G, Cylinder>) {
        source.template Emplace<Model>(SomeCylinder());
    } else if constexpr (std::is_same_v<G, Cone>) {
        source.template Emplace<Model>(SomeCone());
    } else {
        source.template Emplace<Model>(SomeTorus());
    }
    for (auto _ : state) {
        auto copy = source;
        benchmark::DoNotOptimize(copy);
    }
    state.counters["sizeof"] = static_cast<double>(sizeof(G));
    state.counters["shared"] = Stereon::Detail::FitsInline<G, Size> ? 0.0 : 1.0;
}

} // namespace

BENCHMARK_CAPTURE(Direct, Line, SomeLine());
BENCHMARK_CAPTURE(Visit, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(VisitPerCall, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(Erased, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(Direct, Circle, SomeCircle());
BENCHMARK_CAPTURE(Visit, Circle, Curve(SomeCircle()));
BENCHMARK_CAPTURE(VisitPerCall, Circle, Curve(SomeCircle()));
BENCHMARK_CAPTURE(Erased, Circle, Curve(SomeCircle()));
BENCHMARK_CAPTURE(Direct, BSpline, SomeBSpline());
BENCHMARK_CAPTURE(Visit, BSpline, Curve(SomeBSpline()));
BENCHMARK_CAPTURE(VisitPerCall, BSpline, Curve(SomeBSpline()));
BENCHMARK_CAPTURE(Erased, BSpline, Curve(SomeBSpline()));
BENCHMARK_CAPTURE(Direct, Ellipse, SomeEllipse());
BENCHMARK_CAPTURE(Visit, Ellipse, Curve(SomeEllipse()));
BENCHMARK_CAPTURE(Erased, Ellipse, Curve(SomeEllipse()));
BENCHMARK_CAPTURE(DirectScalar, Line, SomeLine());
BENCHMARK_CAPTURE(DirectScalar, Circle, SomeCircle());
BENCHMARK_CAPTURE(DirectScalar, BSpline, SomeBSpline());

BENCHMARK_CAPTURE(Direct, Line / 1e8, SomeLine())->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Visit, Line / 1e8, Curve(SomeLine()))->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Erased, Line / 1e8, Curve(SomeLine()))->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Direct, Circle / 1e8, SomeCircle())->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Visit, Circle / 1e8, Curve(SomeCircle()))->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Erased, Circle / 1e8, Curve(SomeCircle()))->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Direct, BSpline / 1e8, SomeBSpline())->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Visit, BSpline / 1e8, Curve(SomeBSpline()))->Iterations(HundredMillion);
BENCHMARK_CAPTURE(Erased, BSpline / 1e8, Curve(SomeBSpline()))->Iterations(HundredMillion);

BENCHMARK(CopyWrapper<Line, 64>);
BENCHMARK(CopyWrapper<Line, 96>);
BENCHMARK(CopyWrapper<Line, 128>);
BENCHMARK(CopyWrapper<Circle, 64>);
BENCHMARK(CopyWrapper<Circle, 96>);
BENCHMARK(CopyWrapper<Circle, 128>);
BENCHMARK(CopySurfaceWrapper<Cylinder, 96>);
BENCHMARK(CopySurfaceWrapper<Cylinder, 128>);
BENCHMARK(CopySurfaceWrapper<Cone, 96>);
BENCHMARK(CopySurfaceWrapper<Cone, 128>);
BENCHMARK(CopySurfaceWrapper<Torus, 96>);
BENCHMARK(CopySurfaceWrapper<Torus, 128>);
