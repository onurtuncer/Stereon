// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Phase 1 benchmarks of docs/benchmarks.md, Stereon against OpenCASCADE
// on the same geometry:
//
//   B1.1  B-spline curve and surface evaluation, D0 to D2
//   B1.3  point projection onto B-spline curves and surfaces
//   B1.4  curve interpolation and approximation
//   B1.5  knot insertion and Bezier decomposition
//
// Each Stereon benchmark has an OCCT twin with the same name after the
// slash, so the report pairs them. OCCT gets its best usage where it has
// one: projections reuse an initialised GeomAPI projector and call Perform,
// as an application projecting many points would, and Stereon's reuse a
// CurveProjector or SurfaceProjector; evaluation with derivatives is one
// D2 or D1 call against one Derivatives call. Both libraries are
// release builds of the same compiler; tools/occt/build_occt.sh records the
// OCCT build. Run:
//
//   build/<bench preset with OCCT>/bench/stereon_occt_bench --benchmark_repetitions=5
//
// OCCT is LGPL-2.1 and linked into this executable only (ADR-0008).

#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <benchmark/benchmark.h>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

import Stereon.Surface;

namespace {

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Point2;
using Stereon::Point3;

/// Parameters, query points and fitting points per benchmark iteration.
constexpr int ParameterCount = 4096;
constexpr int QueryCount = 256;
constexpr int FitPointCount = 1000;

gp_Pnt ToOcct(const Point3& p)
{
    return {p.X, p.Y, p.Z};
}

/// A clamped curve of the given degree with 24 poles on a helix-like path,
/// uniform interior knots, and weights when rational.
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

/// A clamped bicubic surface over a 12 x 12 net of a wavy sheet.
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

/// Distinct values and multiplicities, 1-based, as OCCT stores knots.
struct OcctKnots {
    NCollection_Array1<double> Values;
    NCollection_Array1<int> Multiplicities;
};

OcctKnots ToOcctKnots(const std::vector<double>& knots)
{
    std::vector<double> values;
    std::vector<int> counts;
    for (const double u : knots) {
        if (values.empty() || values.back() != u) {
            values.push_back(u);
            counts.push_back(0);
        }
        ++counts.back();
    }
    const int n = static_cast<int>(values.size());
    OcctKnots result{NCollection_Array1<double>(1, n), NCollection_Array1<int>(1, n)};
    for (int i = 0; i < n; ++i) {
        result.Values(i + 1) = values[static_cast<std::size_t>(i)];
        result.Multiplicities(i + 1) = counts[static_cast<std::size_t>(i)];
    }
    return result;
}

occ::handle<Geom_BSplineCurve> ToOcct(const BSplineCurve& c)
{
    const int n = static_cast<int>(c.Poles.size());
    NCollection_Array1<gp_Pnt> poles(1, n);
    NCollection_Array1<double> weights(1, n);
    for (int i = 0; i < n; ++i) {
        poles(i + 1) = ToOcct(c.Poles[static_cast<std::size_t>(i)]);
        weights(i + 1) = c.Weights.empty() ? 1.0 : c.Weights[static_cast<std::size_t>(i)];
    }
    const OcctKnots k = ToOcctKnots(c.Knots);
    return new Geom_BSplineCurve(poles, weights, k.Values, k.Multiplicities, c.Degree);
}

occ::handle<Geom_BSplineSurface> ToOcct(const BSplineSurface& s)
{
    const auto nu = static_cast<int>(Stereon::PoleCountU(s));
    const auto nv = static_cast<int>(Stereon::PoleCountV(s));
    NCollection_Array2<gp_Pnt> poles(1, nu, 1, nv);
    NCollection_Array2<double> weights(1, nu, 1, nv);
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nv; ++j) {
            const auto k = static_cast<std::size_t>(i * nv + j);
            poles(i + 1, j + 1) = ToOcct(s.Poles[k]);
            weights(i + 1, j + 1) = s.Weights.empty() ? 1.0 : s.Weights[k];
        }
    }
    const OcctKnots u = ToOcctKnots(s.KnotsU);
    const OcctKnots v = ToOcctKnots(s.KnotsV);
    return new Geom_BSplineSurface(poles, weights, u.Values, v.Values, u.Multiplicities,
                                   v.Multiplicities, s.DegreeU, s.DegreeV);
}

std::vector<double> CurveParameters(const BSplineCurve& c)
{
    const Stereon::Interval d = Domain(c);
    std::vector<double> t(ParameterCount);
    for (int i = 0; i < ParameterCount; ++i) {
        t[static_cast<std::size_t>(i)] =
            d.GetLower() + (d.GetUpper() - d.GetLower()) * (i + 0.5) / ParameterCount;
    }
    return t;
}

std::vector<Point2> SurfaceParameters(const BSplineSurface& s)
{
    const Stereon::Box2 d = Domain(s);
    std::vector<Point2> uv;
    constexpr int side = 64;
    for (int i = 0; i < side; ++i) {
        for (int j = 0; j < side; ++j) {
            uv.push_back({d.X.GetLower() + (d.X.GetUpper() - d.X.GetLower()) * (i + 0.5) / side,
                          d.Y.GetLower() + (d.Y.GetUpper() - d.Y.GetLower()) * (j + 0.5) / side});
        }
    }
    return uv;
}

/// Query points scattered around the geometry's bounding box.
std::vector<Point3> Queries(const Stereon::Box3& box)
{
    std::vector<Point3> points;
    std::uint32_t state = 12345;
    const auto next = [&state] {
        state = state * 1664525u + 1013904223u; // a fixed LCG: same points every run
        return static_cast<double>(state >> 8) / static_cast<double>(1u << 24);
    };
    for (int i = 0; i < QueryCount; ++i) {
        points.push_back(
            {box.X.GetLower() + next() * (box.X.GetUpper() - box.X.GetLower()),
             box.Y.GetLower() + next() * (box.Y.GetUpper() - box.Y.GetLower()),
             box.Z.GetLower() - 2.0 + next() * (box.Z.GetUpper() - box.Z.GetLower() + 4.0)});
    }
    return points;
}

std::vector<Point3> FitPoints()
{
    std::vector<Point3> points;
    for (int k = 0; k < FitPointCount; ++k) {
        const double s = 0.02 * k;
        points.push_back({s, std::sin(s) * 5.0, std::cos(0.7 * s) * 3.0});
    }
    return points;
}

void Items(benchmark::State& state, const std::size_t perIteration)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(perIteration));
}

// ---------------------------------------------------------------------------
// B1.1: evaluation

template <int Degree, bool Rational>
void CurveD0(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(Degree, Rational);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            sum += Eval(c, u).X;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

template <int Degree, bool Rational>
void OcctCurveD0(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(Degree, Rational);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            sum += occt->Value(u).X();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

template <int Degree, bool Rational>
void CurveD2(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(Degree, Rational);
    const std::vector<double> t = CurveParameters(c);
    std::array<Stereon::Vec3, 3> d{};
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            Derivatives(c, u, 2, std::span<Stereon::Vec3>{d});
            sum += d[0].X + d[1].Y + d[2].Z;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

template <int Degree, bool Rational>
void OcctCurveD2(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(Degree, Rational);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            gp_Pnt p;
            gp_Vec d1;
            gp_Vec d2;
            occt->D2(u, p, d1, d2);
            sum += p.X() + d1.Y() + d2.Z();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

template <bool Rational>
void SurfaceD0(benchmark::State& state)
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
    Items(state, uv.size());
}

template <bool Rational>
void OcctSurfaceD0(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const occ::handle<Geom_BSplineSurface> occt = ToOcct(s);
    const std::vector<Point2> uv = SurfaceParameters(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            sum += occt->Value(p.X, p.Y).Z();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, uv.size());
}

template <bool Rational>
void SurfaceD1(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const std::vector<Point2> uv = SurfaceParameters(s);
    std::array<Stereon::Vec3, 4> d{}; // (0,0), (0,1), (1,0), (1,1)
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            Derivatives(s, p, 1, 1, std::span<Stereon::Vec3>{d});
            sum += d[0].Z + d[2].Z + d[1].Z;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, uv.size());
}

template <bool Rational>
void OcctSurfaceD1(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const occ::handle<Geom_BSplineSurface> occt = ToOcct(s);
    const std::vector<Point2> uv = SurfaceParameters(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            gp_Pnt point;
            gp_Vec du;
            gp_Vec dv;
            occt->D1(p.X, p.Y, point, du, dv);
            sum += point.Z() + du.Z() + dv.Z();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, uv.size());
}

// ---------------------------------------------------------------------------
// B1.3: projection

void CurveProjection(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const std::vector<Point3> q = Queries(BoundingBox(c));
    const Stereon::CurveProjector projector(c);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point3& p : q) {
            sum += projector.Project(p);
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, q.size());
}

void OcctCurveProjection(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    const std::vector<Point3> q = Queries(BoundingBox(c));
    GeomAPI_ProjectPointOnCurve projector;
    projector.Init(occt, occt->FirstParameter(), occt->LastParameter());
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point3& p : q) {
            projector.Perform(ToOcct(p));
            sum += projector.NbPoints() > 0 ? projector.LowerDistanceParameter() : 0.0;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, q.size());
}

void SurfaceProjection(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const std::vector<Point3> q = Queries(BoundingBox(s));
    const Stereon::SurfaceProjector projector(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point3& p : q) {
            sum += projector.Project(p).X;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, q.size());
}

void OcctSurfaceProjection(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const occ::handle<Geom_BSplineSurface> occt = ToOcct(s);
    const std::vector<Point3> q = Queries(BoundingBox(s));
    double u0 = 0.0;
    double u1 = 0.0;
    double v0 = 0.0;
    double v1 = 0.0;
    occt->Bounds(u0, u1, v0, v1);
    GeomAPI_ProjectPointOnSurf projector;
    projector.Init(occt, u0, u1, v0, v1);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point3& p : q) {
            projector.Perform(ToOcct(p));
            double u = 0.0;
            double v = 0.0;
            if (projector.NbPoints() > 0) {
                projector.LowerDistanceParameters(u, v);
            }
            sum += u;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, q.size());
}

// ---------------------------------------------------------------------------
// B1.4: fitting

void Interpolation(benchmark::State& state)
{
    const std::vector<Point3> points = FitPoints();
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::Interpolate(points, 3));
    }
    Items(state, points.size());
}

void OcctInterpolation(benchmark::State& state)
{
    const std::vector<Point3> points = FitPoints();
    occ::handle<NCollection_HArray1<gp_Pnt>> occtPoints =
        new NCollection_HArray1<gp_Pnt>(1, static_cast<int>(points.size()));
    for (std::size_t k = 0; k < points.size(); ++k) {
        occtPoints->SetValue(static_cast<int>(k) + 1, ToOcct(points[k]));
    }
    for (auto _ : state) {
        GeomAPI_Interpolate interpolate(occtPoints, false, 1e-7);
        interpolate.Perform();
        benchmark::DoNotOptimize(interpolate.Curve());
    }
    Items(state, points.size());
}

void Approximation(benchmark::State& state)
{
    const std::vector<Point3> points = FitPoints();
    for (auto _ : state) {
        benchmark::DoNotOptimize(
            Stereon::Approximate(points, 3, Stereon::ApproximationBudget{1e-3}));
    }
    Items(state, points.size());
}

void OcctApproximation(benchmark::State& state)
{
    const std::vector<Point3> points = FitPoints();
    NCollection_Array1<gp_Pnt> occtPoints(1, static_cast<int>(points.size()));
    for (std::size_t k = 0; k < points.size(); ++k) {
        occtPoints(static_cast<int>(k) + 1) = ToOcct(points[k]);
    }
    for (auto _ : state) {
        GeomAPI_PointsToBSpline approximate(occtPoints, 3, 3, GeomAbs_C2, 1e-3);
        benchmark::DoNotOptimize(approximate.Curve());
    }
    Items(state, points.size());
}

// ---------------------------------------------------------------------------
// B1.5: knot operations

void KnotInsertion(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    for (auto _ : state) {
        BSplineCurve refined = c;
        for (int i = 0; i < 100; ++i) {
            refined = *Stereon::InsertKnot(refined, 0.2 + 0.2 * i);
        }
        benchmark::DoNotOptimize(refined);
    }
    Items(state, 100);
}

void OcctKnotInsertion(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    for (auto _ : state) {
        occ::handle<Geom_BSplineCurve> refined = occ::down_cast<Geom_BSplineCurve>(occt->Copy());
        for (int i = 0; i < 100; ++i) {
            refined->InsertKnot(0.2 + 0.2 * i);
        }
        benchmark::DoNotOptimize(refined);
    }
    Items(state, 100);
}

void BezierDecomposition(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::ToBezierSegments(c));
    }
    Items(state, 1);
}

void OcctBezierDecomposition(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    for (auto _ : state) {
        // OCCT builds the arcs on demand: ask for every one, as
        // ToBezierSegments returns them all.
        GeomConvert_BSplineCurveToBezierCurve convert(occt);
        for (int i = 1; i <= convert.NbArcs(); ++i) {
            benchmark::DoNotOptimize(convert.Arc(i));
        }
    }
    Items(state, 1);
}

} // namespace

BENCHMARK(CurveD0<3, false>)->Name("B1.1/CurveD0/cubic/Stereon");
BENCHMARK(OcctCurveD0<3, false>)->Name("B1.1/CurveD0/cubic/OCCT");
BENCHMARK(CurveD0<3, true>)->Name("B1.1/CurveD0/cubic-rational/Stereon");
BENCHMARK(OcctCurveD0<3, true>)->Name("B1.1/CurveD0/cubic-rational/OCCT");
BENCHMARK(CurveD0<7, false>)->Name("B1.1/CurveD0/degree7/Stereon");
BENCHMARK(OcctCurveD0<7, false>)->Name("B1.1/CurveD0/degree7/OCCT");
BENCHMARK(CurveD2<3, false>)->Name("B1.1/CurveD2/cubic/Stereon");
BENCHMARK(OcctCurveD2<3, false>)->Name("B1.1/CurveD2/cubic/OCCT");
BENCHMARK(CurveD2<3, true>)->Name("B1.1/CurveD2/cubic-rational/Stereon");
BENCHMARK(OcctCurveD2<3, true>)->Name("B1.1/CurveD2/cubic-rational/OCCT");
BENCHMARK(SurfaceD0<false>)->Name("B1.1/SurfaceD0/bicubic/Stereon");
BENCHMARK(OcctSurfaceD0<false>)->Name("B1.1/SurfaceD0/bicubic/OCCT");
BENCHMARK(SurfaceD0<true>)->Name("B1.1/SurfaceD0/bicubic-rational/Stereon");
BENCHMARK(OcctSurfaceD0<true>)->Name("B1.1/SurfaceD0/bicubic-rational/OCCT");
BENCHMARK(SurfaceD1<false>)->Name("B1.1/SurfaceD1/bicubic/Stereon");
BENCHMARK(OcctSurfaceD1<false>)->Name("B1.1/SurfaceD1/bicubic/OCCT");
BENCHMARK(CurveProjection)->Name("B1.3/CurveProjection/Stereon");
BENCHMARK(OcctCurveProjection)->Name("B1.3/CurveProjection/OCCT");
BENCHMARK(SurfaceProjection)->Name("B1.3/SurfaceProjection/Stereon");
BENCHMARK(OcctSurfaceProjection)->Name("B1.3/SurfaceProjection/OCCT");
BENCHMARK(Interpolation)->Name("B1.4/Interpolation/1000/Stereon");
BENCHMARK(OcctInterpolation)->Name("B1.4/Interpolation/1000/OCCT");
BENCHMARK(Approximation)->Name("B1.4/Approximation/1000/Stereon");
BENCHMARK(OcctApproximation)->Name("B1.4/Approximation/1000/OCCT");
BENCHMARK(KnotInsertion)->Name("B1.5/KnotInsertion/100/Stereon");
BENCHMARK(OcctKnotInsertion)->Name("B1.5/KnotInsertion/100/OCCT");
BENCHMARK(BezierDecomposition)->Name("B1.5/BezierDecomposition/Stereon");
BENCHMARK(OcctBezierDecomposition)->Name("B1.5/BezierDecomposition/OCCT");
