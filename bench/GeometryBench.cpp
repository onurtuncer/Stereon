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
// - Erased: Eval(curve, t), one virtual call per evaluation. ADR-0003 asks
//   for this within 2x of direct.
//
// The cubic B-spline of the ADR's verification and the OCCT side
// (Geom_Curve::D0 through occ::handle) join when the B-spline type and the
// OCCT build exist. The copy benchmarks measure the wrapper at buffer sizes
// 64, 96 and 128 bytes against the built-in types' sizes, which is how
// ADR-0003 fixes GeometryBufferSize.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <numbers>
#include <type_traits>
#include <vector>

import Stereon.Surface;

namespace {

using Stereon::Circle;
using Stereon::Curve;
using Stereon::Line;
using Stereon::Point3;
using Stereon::Vec3;

constexpr std::size_t SampleCount = 1 << 16;

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

} // namespace

BENCHMARK_CAPTURE(Direct, Line, SomeLine());
BENCHMARK_CAPTURE(Visit, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(VisitPerCall, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(Erased, Line, Curve(SomeLine()));
BENCHMARK_CAPTURE(Direct, Circle, SomeCircle());
BENCHMARK_CAPTURE(Visit, Circle, Curve(SomeCircle()));
BENCHMARK_CAPTURE(VisitPerCall, Circle, Curve(SomeCircle()));
BENCHMARK_CAPTURE(Erased, Circle, Curve(SomeCircle()));

BENCHMARK(CopyWrapper<Line, 64>);
BENCHMARK(CopyWrapper<Line, 96>);
BENCHMARK(CopyWrapper<Line, 128>);
BENCHMARK(CopyWrapper<Circle, 64>);
BENCHMARK(CopyWrapper<Circle, 96>);
BENCHMARK(CopyWrapper<Circle, 128>);
