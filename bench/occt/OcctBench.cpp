// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The curve and surface benchmarks of docs/benchmarks.md, Stereon against
// OpenCASCADE on the same geometry:
//
//   B0.2  curve dispatch, 10^8 evaluations of a line, a circle and a cubic B-spline
//   B1.1  B-spline curve and surface evaluation, D0 to D2, and a batch of curve points
//   B1.2  grid evaluation of a B-spline surface
//   B1.3  point projection onto B-spline curves and surfaces
//   B1.4  curve interpolation and approximation
//   B1.5  knot insertion, degree elevation and Bezier decomposition
//   B1.6  local differential properties and arc length
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

#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <GeomGridEval_BSplineCurve.hxx>
#include <GeomGridEval_BSplineSurface.hxx>
#include <GeomLProp_CLProps.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <benchmark/benchmark.h>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
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
// B0.2: curve dispatch
//
// 10^8 evaluations per measurement of the curves of bench/GeometryBench.cpp,
// a line, a circle and ADR-0003's cubic B-spline, along three paths:
//
//   Direct    the concrete type: Eval(line, t) against EvalD0 through
//             occ::handle<Geom_Line>, where EvalD0 is final, so not virtual
//   Dispatch  VisitBuiltin, one type test and then the concrete loop, against
//             EvalD0 through occ::handle<Geom_Curve>, a virtual call per
//             evaluation: the comparison B0.2 names
//   Virtual   the erased Curve, a virtual call per evaluation, against the
//             same OCCT virtual call: like for like
//
// B0.2's other target, VisitBuiltin within 1.1x of the direct call, compares
// two Stereon rows; the report has a table for it. OCCT's line is
// parametrised by arc length and Stereon's by Origin + t Direction, so OCCT's
// parameters are scaled by |Direction|. The counter sum_x_mm, the x
// coordinates summed over one pass, shows both evaluate the same points.

constexpr std::size_t DispatchSampleCount = 1 << 16;

/// Iterations of DispatchSampleCount evaluations that make at least 10^8.
constexpr auto HundredMillion = static_cast<benchmark::IterationCount>(
    (100'000'000 + DispatchSampleCount - 1) / DispatchSampleCount);

Stereon::Line DispatchLine()
{
    return {Point3{1.0, 2.0, 3.0}, Stereon::Vec3{0.5, -1.0, 2.0}};
}

Stereon::Circle DispatchCircle()
{
    return {Point3{1.0, 2.0, 3.0}, 2.5, Stereon::Vec3{1.0, 0.0, 0.0}, Stereon::Vec3{0.0, 1.0, 0.0}};
}

/// ADR-0003's cubic B-spline: eight poles over [0, 2 pi], the range the
/// parameters cover.
BSplineCurve DispatchBSpline()
{
    constexpr double end = 2.0 * std::numbers::pi;
    return {3,
            {0.0, 0.0, 0.0, 0.0, 0.2 * end, 0.4 * end, 0.6 * end, 0.8 * end, end, end, end, end},
            {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
             Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}, Point3{5.0, 2.0, 1.0},
             Point3{6.0, -1.0, 0.0}, Point3{7.0, 0.0, 1.0}},
            {}};
}

gp_Vec ToOcct(const Stereon::Vec3& v)
{
    return {v.X, v.Y, v.Z};
}

occ::handle<Geom_Line> ToOcct(const Stereon::Line& line)
{
    return new Geom_Line(ToOcct(line.Origin), gp_Dir(ToOcct(line.Direction)));
}

occ::handle<Geom_Circle> ToOcct(const Stereon::Circle& circle)
{
    const gp_Vec x = ToOcct(circle.XAxis);
    const gp_Vec y = ToOcct(circle.YAxis);
    return new Geom_Circle(gp_Ax2(ToOcct(circle.Center), gp_Dir(x.Crossed(y)), gp_Dir(x)),
                           circle.Radius);
}

/// OCCT's parameter per unit of Stereon's: |Direction| for a line, which OCCT
/// parametrises by arc length, 1 for the others.
double OcctParameterScale(const Stereon::Line& line)
{
    return ToOcct(line.Direction).Magnitude();
}

template <typename G>
double OcctParameterScale(const G&)
{
    return 1.0;
}

/// DispatchSampleCount parameters over [0, 2 pi), times scale.
std::vector<double> DispatchParameters(const double scale)
{
    std::vector<double> t(DispatchSampleCount);
    for (std::size_t i = 0; i < DispatchSampleCount; ++i) {
        t[i] = scale * static_cast<double>(i) * (2.0 * std::numbers::pi / DispatchSampleCount);
    }
    return t;
}

template <typename G>
double SumX(const G& geometry, const std::vector<double>& t)
{
    double sum = 0.0;
    for (const double u : t) {
        sum += Eval(geometry, u).X;
    }
    return sum;
}

/// The same sum through an OCCT handle: occ::handle<Geom_Line> calls the
/// final EvalD0 directly, occ::handle<Geom_Curve> through the vtable.
template <typename C>
double OcctSumX(const occ::handle<C>& curve, const std::vector<double>& t)
{
    double sum = 0.0;
    for (const double u : t) {
        sum += curve->EvalD0(u).X();
    }
    return sum;
}

void DispatchItems(benchmark::State& state, const double sum)
{
    Items(state, DispatchSampleCount);
    state.counters["sum_x_mm"] = sum;
}

template <typename G>
void DispatchDirect(benchmark::State& state, const G geometry)
{
    const std::vector<double> t = DispatchParameters(1.0);
    for (auto _ : state) {
        benchmark::DoNotOptimize(SumX(geometry, t));
    }
    DispatchItems(state, SumX(geometry, t));
}

template <typename G>
void DispatchVisit(benchmark::State& state, const G geometry)
{
    const Stereon::Curve curve(geometry);
    const std::vector<double> t = DispatchParameters(1.0);
    const auto pass = [&] {
        return VisitBuiltin(curve, [&t](const auto& g) { return SumX(g, t); });
    };
    for (auto _ : state) {
        benchmark::DoNotOptimize(pass());
    }
    DispatchItems(state, pass());
}

template <typename G>
void DispatchErased(benchmark::State& state, const G geometry)
{
    const Stereon::Curve curve(geometry);
    const std::vector<double> t = DispatchParameters(1.0);
    for (auto _ : state) {
        benchmark::DoNotOptimize(SumX(curve, t));
    }
    DispatchItems(state, SumX(curve, t));
}

/// The OCCT twin of Direct: the concrete handle. The curve is built here, not
/// at registration, so no OCCT object exists before main().
template <typename G>
void OcctDispatchDirect(benchmark::State& state, const G geometry)
{
    const auto curve = ToOcct(geometry);
    const std::vector<double> t = DispatchParameters(OcctParameterScale(geometry));
    for (auto _ : state) {
        benchmark::DoNotOptimize(OcctSumX(curve, t));
    }
    DispatchItems(state, OcctSumX(curve, t));
}

/// The OCCT twin of Dispatch and Virtual: the base-class handle.
template <typename G>
void OcctDispatchVirtual(benchmark::State& state, const G geometry)
{
    const occ::handle<Geom_Curve> curve = ToOcct(geometry);
    const std::vector<double> t = DispatchParameters(OcctParameterScale(geometry));
    for (auto _ : state) {
        benchmark::DoNotOptimize(OcctSumX(curve, t));
    }
    DispatchItems(state, OcctSumX(curve, t));
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
    // The point and the two first partials, as OCCT's D1 below: not the
    // mixed partial that Derivatives(s, p, 1, 1, d) would also give.
    std::array<Stereon::Vec3, 4> d{}; // (0,0), (0,1), (1,0); d[3] is not written
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            TotalDerivatives(s, p, 1, std::span<Stereon::Vec3>{d});
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
// B1.1 batch and B1.2 grid: many points in one call
//
// OCCT's GeomGridEval returns a new array per call, whose gp_Pnt elements
// are zeroed on construction; Stereon's evaluators write into a span, so
// its side allocates the same zeroed array per iteration too.

/// The parameters of CurveParameters as OCCT's 1-based array.
NCollection_Array1<double> ToOcctArray(const std::vector<double>& values)
{
    NCollection_Array1<double> array(1, static_cast<int>(values.size()));
    for (std::size_t k = 0; k < values.size(); ++k) {
        array(static_cast<int>(k) + 1) = values[k];
    }
    return array;
}

/// n parameters of the interval at the midpoints of n equal steps.
std::vector<double> GridParameters(const Stereon::Interval& d, const int n)
{
    std::vector<double> t(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        t[static_cast<std::size_t>(i)] =
            d.GetLower() + (d.GetUpper() - d.GetLower()) * (i + 0.5) / n;
    }
    return t;
}

constexpr int GridSide = 256;

template <bool Rational>
void CurveD0Batch(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, Rational);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        std::vector<Point3> points(t.size());
        Stereon::EvalAll(c, t, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    Items(state, t.size());
}

/// EvalAllFast, with the conversion to span polynomials inside the timed
/// loop: the cost a caller pays who converts for every call.
template <bool Rational>
void CurveD0BatchFast(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, Rational);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        std::vector<Point3> points(t.size());
        Stereon::EvalAllFast(c, t, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    Items(state, t.size());
}

template <bool Rational>
void OcctCurveD0Batch(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, Rational);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    const std::vector<double> t = CurveParameters(c);
    const NCollection_Array1<double> parameters = ToOcctArray(t);
    const GeomGridEval_BSplineCurve evaluator(occt);
    for (auto _ : state) {
        const NCollection_Array1<gp_Pnt> points = evaluator.EvaluateGrid(parameters);
        benchmark::DoNotOptimize(&points.First());
        benchmark::ClobberMemory();
    }
    Items(state, t.size());
}

template <bool Rational>
void SurfaceGrid(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const Stereon::Box2 d = Domain(s);
    const std::vector<double> u = GridParameters(d.X, GridSide);
    const std::vector<double> v = GridParameters(d.Y, GridSide);
    for (auto _ : state) {
        std::vector<Point3> points(u.size() * v.size());
        Stereon::EvalGrid(s, u, v, points);
        benchmark::DoNotOptimize(points.data());
        benchmark::ClobberMemory();
    }
    Items(state, u.size() * v.size());
}

template <bool Rational>
void OcctSurfaceGrid(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(Rational);
    const occ::handle<Geom_BSplineSurface> occt = ToOcct(s);
    const Stereon::Box2 d = Domain(s);
    const NCollection_Array1<double> u = ToOcctArray(GridParameters(d.X, GridSide));
    const NCollection_Array1<double> v = ToOcctArray(GridParameters(d.Y, GridSide));
    const GeomGridEval_BSplineSurface evaluator(occt);
    for (auto _ : state) {
        const NCollection_Array2<gp_Pnt> points = evaluator.EvaluateGrid(u, v);
        benchmark::DoNotOptimize(&points.Value(1, 1));
        benchmark::ClobberMemory();
    }
    Items(state, static_cast<std::size_t>(GridSide) * GridSide);
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
        // A named handle: DoNotOptimize's const& overload is deprecated for
        // types that are not trivially copyable, and -Werror stops the build.
        auto curve = interpolate.Curve();
        benchmark::DoNotOptimize(curve);
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
        auto curve = approximate.Curve();
        benchmark::DoNotOptimize(curve);
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

void DegreeElevation(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    for (auto _ : state) {
        benchmark::DoNotOptimize(*Stereon::ElevatedDegree(c, 5));
    }
    Items(state, 1);
}

void OcctDegreeElevation(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    for (auto _ : state) {
        occ::handle<Geom_BSplineCurve> raised = occ::down_cast<Geom_BSplineCurve>(occt->Copy());
        raised->IncreaseDegree(5);
        benchmark::DoNotOptimize(raised);
    }
    Items(state, 1);
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

// ---------------------------------------------------------------------------
// B1.6: local properties and arc length

/// Normalised tolerance for the property evaluators of both libraries:
/// OCCT treats derivatives shorter than this as zero, Stereon tests exactly.
constexpr double PropertyResolution = 1e-9;

/// Abscissas per iteration of the point-at-length benchmarks.
constexpr int AbscissaCount = 64;

/// The relative accuracy asked of OCCT's arc length, matching Stereon's
/// internal target of about 1e-12 (Detail::ArcLengthRelativeTarget).
constexpr double LengthTolerance = 1e-12;

/// Tangent, principal normal and curvature at every parameter, from one
/// evaluation in both libraries.
void CurveProperties(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const std::vector<double> t = CurveParameters(c);
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            const auto frame = Stereon::FrenetFrame(c, u);
            sum += frame->Curvature + frame->Normal.X + frame->Tangent.X;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

void OcctCurveProperties(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const occ::handle<Geom_BSplineCurve> occt = ToOcct(c);
    const std::vector<double> t = CurveParameters(c);
    GeomLProp_CLProps props(occt, 2, PropertyResolution);
    for (auto _ : state) {
        double sum = 0.0;
        for (const double u : t) {
            props.SetParameter(u);
            sum += props.Curvature();
            gp_Dir normal;
            props.Normal(normal);
            gp_Dir tangent;
            props.Tangent(tangent);
            sum += normal.X() + tangent.X();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, t.size());
}

/// Principal curvatures and directions at every parameter pair.
void SurfaceProperties(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const std::vector<Point2> uv = SurfaceParameters(s);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            const auto k = Stereon::PrincipalCurvatures(s, p);
            sum += k->Min + k->Max + k->MaxDirection.X;
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, uv.size());
}

void OcctSurfaceProperties(benchmark::State& state)
{
    const BSplineSurface s = SomeSurface(false);
    const occ::handle<Geom_BSplineSurface> occt = ToOcct(s);
    const std::vector<Point2> uv = SurfaceParameters(s);
    GeomLProp_SLProps props(occt, 2, PropertyResolution);
    for (auto _ : state) {
        double sum = 0.0;
        for (const Point2& p : uv) {
            props.SetParameters(p.X, p.Y);
            sum += props.MinCurvature() + props.MaxCurvature();
            gp_Dir maxDirection;
            gp_Dir minDirection;
            props.CurvatureDirections(maxDirection, minDirection);
            sum += maxDirection.X();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, uv.size());
}

/// The length of the whole curve.
void CurveLength(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    for (auto _ : state) {
        benchmark::DoNotOptimize(Stereon::ArcLength(c, Domain(c))->Length);
    }
    Items(state, 1);
    // The length itself, to compare the two libraries' answers.
    state.counters["length_mm"] = Stereon::ArcLength(c, Domain(c))->Length;
}

void OcctCurveLength(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const GeomAdaptor_Curve occt(ToOcct(c));
    for (auto _ : state) {
        benchmark::DoNotOptimize(GCPnts_AbscissaPoint::Length(
            occt, occt.FirstParameter(), occt.LastParameter(), LengthTolerance));
    }
    Items(state, 1);
    state.counters["length_mm"] = GCPnts_AbscissaPoint::Length(
        occt, occt.FirstParameter(), occt.LastParameter(), LengthTolerance);
}

/// Parameters at evenly spaced lengths from the start, each from scratch.
void PointAtLength(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const double length = Stereon::ArcLength(c, Domain(c))->Length;
    const double start = Domain(c).GetLower();
    for (auto _ : state) {
        double sum = 0.0;
        for (int i = 1; i <= AbscissaCount; ++i) {
            sum += *Stereon::ParameterAtArcLength(c, start, length * i / (AbscissaCount + 1));
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, AbscissaCount);
}

void OcctPointAtLength(benchmark::State& state)
{
    const BSplineCurve c = SomeCurve(3, true);
    const GeomAdaptor_Curve occt(ToOcct(c));
    const double start = occt.FirstParameter();
    const double length =
        GCPnts_AbscissaPoint::Length(occt, start, occt.LastParameter(), LengthTolerance);
    for (auto _ : state) {
        double sum = 0.0;
        for (int i = 1; i <= AbscissaCount; ++i) {
            const GCPnts_AbscissaPoint point(LengthTolerance, occt,
                                             length * i / (AbscissaCount + 1), start);
            sum += point.Parameter();
        }
        benchmark::DoNotOptimize(sum);
    }
    Items(state, AbscissaCount);
}

} // namespace

// B0.2: each curve along each path, 10^8 evaluations per measurement.
#define STEREON_B0_2(curve, make)                                                                  \
    BENCHMARK_CAPTURE(DispatchDirect, curve, make())                                               \
        ->Name("B0.2/Direct/" #curve "/Stereon")                                                   \
        ->Iterations(HundredMillion);                                                              \
    BENCHMARK_CAPTURE(OcctDispatchDirect, curve, make())                                           \
        ->Name("B0.2/Direct/" #curve "/OCCT")                                                      \
        ->Iterations(HundredMillion);                                                              \
    BENCHMARK_CAPTURE(DispatchVisit, curve, make())                                                \
        ->Name("B0.2/Dispatch/" #curve "/Stereon")                                                 \
        ->Iterations(HundredMillion);                                                              \
    BENCHMARK_CAPTURE(OcctDispatchVirtual, curve, make())                                          \
        ->Name("B0.2/Dispatch/" #curve "/OCCT")                                                    \
        ->Iterations(HundredMillion);                                                              \
    BENCHMARK_CAPTURE(DispatchErased, curve, make())                                               \
        ->Name("B0.2/Virtual/" #curve "/Stereon")                                                  \
        ->Iterations(HundredMillion);                                                              \
    BENCHMARK_CAPTURE(OcctDispatchVirtual, curve, make())                                          \
        ->Name("B0.2/Virtual/" #curve "/OCCT")                                                     \
        ->Iterations(HundredMillion)

STEREON_B0_2(line, DispatchLine);
STEREON_B0_2(circle, DispatchCircle);
STEREON_B0_2(bspline, DispatchBSpline);

#undef STEREON_B0_2

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
BENCHMARK(CurveD0Batch<false>)->Name("B1.1/CurveD0Batch/cubic/Stereon");
BENCHMARK(OcctCurveD0Batch<false>)->Name("B1.1/CurveD0Batch/cubic/OCCT");
BENCHMARK(CurveD0Batch<true>)->Name("B1.1/CurveD0Batch/cubic-rational/Stereon");
BENCHMARK(OcctCurveD0Batch<true>)->Name("B1.1/CurveD0Batch/cubic-rational/OCCT");
BENCHMARK(CurveD0BatchFast<false>)->Name("B1.1/CurveD0BatchFast/cubic/Stereon");
BENCHMARK(OcctCurveD0Batch<false>)->Name("B1.1/CurveD0BatchFast/cubic/OCCT");
BENCHMARK(CurveD0BatchFast<true>)->Name("B1.1/CurveD0BatchFast/cubic-rational/Stereon");
BENCHMARK(OcctCurveD0Batch<true>)->Name("B1.1/CurveD0BatchFast/cubic-rational/OCCT");
BENCHMARK(SurfaceGrid<false>)->Name("B1.2/SurfaceGrid/256/bicubic/Stereon");
BENCHMARK(OcctSurfaceGrid<false>)->Name("B1.2/SurfaceGrid/256/bicubic/OCCT");
BENCHMARK(SurfaceGrid<true>)->Name("B1.2/SurfaceGrid/256/bicubic-rational/Stereon");
BENCHMARK(OcctSurfaceGrid<true>)->Name("B1.2/SurfaceGrid/256/bicubic-rational/OCCT");
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
BENCHMARK(DegreeElevation)->Name("B1.5/DegreeElevation/Stereon");
BENCHMARK(OcctDegreeElevation)->Name("B1.5/DegreeElevation/OCCT");
BENCHMARK(BezierDecomposition)->Name("B1.5/BezierDecomposition/Stereon");
BENCHMARK(OcctBezierDecomposition)->Name("B1.5/BezierDecomposition/OCCT");
BENCHMARK(CurveProperties)->Name("B1.6/CurveProperties/rational-cubic/Stereon");
BENCHMARK(OcctCurveProperties)->Name("B1.6/CurveProperties/rational-cubic/OCCT");
BENCHMARK(SurfaceProperties)->Name("B1.6/SurfaceProperties/bicubic/Stereon");
BENCHMARK(OcctSurfaceProperties)->Name("B1.6/SurfaceProperties/bicubic/OCCT");
BENCHMARK(CurveLength)->Name("B1.6/ArcLength/rational-cubic/Stereon");
BENCHMARK(OcctCurveLength)->Name("B1.6/ArcLength/rational-cubic/OCCT");
BENCHMARK(PointAtLength)->Name("B1.6/PointAtLength/64/Stereon");
BENCHMARK(OcctPointAtLength)->Name("B1.6/PointAtLength/64/OCCT");
