// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Writes, as JSON on standard output, the data behind the figures of the
// manual (docs/fitting.rst, docs/properties.rst): every curve and number in them is what
// Stereon computes, not a sketch. docs/figures/make_fitting_figures.py runs
// this program and renders the SVG files next to it:
//
//   python docs/figures/make_fitting_figures.py build/<preset>/tests/stereon_figure_data
//
// Output goes through std::printf rather than std::format: formatting
// string_views in a program that also links module code defines a libstdc++
// static twice on MinGW (docs/toolchain/gcc16-feature-probe.md).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <numbers>
#include <span>
#include <vector>

import Stereon.Surface;

using namespace Stereon;

namespace {

void Number(const double x)
{
    std::printf("%.10g", x);
}

void Points(const std::span<const Point3> points)
{
    std::printf("[");
    for (std::size_t i = 0; i < points.size(); ++i) {
        std::printf(i == 0 ? "[" : ",[");
        Number(points[i].X);
        std::printf(",");
        Number(points[i].Y);
        std::printf(",");
        Number(points[i].Z);
        std::printf("]");
    }
    std::printf("]");
}

void Numbers(const std::span<const double> values)
{
    std::printf("[");
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            std::printf(",");
        }
        Number(values[i]);
    }
    std::printf("]");
}

/// The curve sampled at `count` evenly spaced parameters of its domain.
std::vector<Point3> Sampled(const BSplineCurve& curve, const int count)
{
    const Interval d = Domain(curve);
    std::vector<Point3> samples;
    for (int i = 0; i <= count; ++i) {
        const double t = d.GetLower() + i * (d.GetUpper() - d.GetLower()) / count;
        samples.push_back(Eval(curve, std::min(t, d.GetUpper())));
    }
    return samples;
}

/// A curve with a sharp turn and uneven spacing: where the three
/// parametrisations differ most.
std::vector<Point3> TurningPoints()
{
    return {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.1, 0.0}, Point3{2.0, 0.0, 0.0},
            Point3{2.2, 0.6, 0.0}, Point3{2.3, 2.4, 0.0}, Point3{3.0, 2.5, 0.0},
            Point3{5.5, 2.6, 0.0}, Point3{6.0, 1.0, 0.0}};
}

/// Measurements of a smooth curve with a deterministic, noise-like error.
std::vector<Point3> NoisyPoints()
{
    std::vector<Point3> points;
    constexpr int count = 60;
    for (int k = 0; k < count; ++k) {
        const double s = 2.0 * std::numbers::pi * k / (count - 1);
        const double noise = 0.04 * std::sin(37.0 * s) + 0.03 * std::cos(91.0 * s);
        points.push_back({s, std::sin(s) + 0.3 * std::sin(3.0 * s) + noise, 0.0});
    }
    return points;
}

void Parametrizations()
{
    const std::vector<Point3> points = TurningPoints();
    std::printf("\"parametrization\":{\"points\":");
    Points(points);
    const char* names[] = {"uniform", "chord", "centripetal"};
    const Parametrization kinds[] = {Parametrization::Uniform, Parametrization::ChordLength,
                                     Parametrization::Centripetal};
    for (int i = 0; i < 3; ++i) {
        const CurveFit fit = *Interpolate(points, 3, kinds[i]);
        std::printf(",\"%s\":{\"curve\":", names[i]);
        Points(Sampled(fit.Curve, 400));
        std::printf(",\"parameters\":");
        Numbers(fit.Parameters);
        std::printf(",\"knots\":");
        Numbers(fit.Curve.Knots);
        std::printf("}");
    }
    std::printf("}");
}

void LeastSquares()
{
    const std::vector<Point3> points = NoisyPoints();
    std::printf("\"leastSquares\":{\"points\":");
    Points(points);
    std::printf(",\"fits\":[");
    const std::size_t counts[] = {6, 10, 18};
    for (std::size_t i = 0; i < 3; ++i) {
        const CurveFit fit = *FitLeastSquares(points, 3, counts[i]);
        std::printf(i == 0 ? "{" : ",{");
        std::printf("\"poleCount\":%zu,\"deviation\":", counts[i]);
        Number(fit.MaxDeviation);
        std::printf(",\"curve\":");
        Points(Sampled(fit.Curve, 400));
        std::printf(",\"poles\":");
        Points(fit.Curve.Poles);
        std::printf("}");
    }
    const CurveFit interpolant = *Interpolate(points, 3);
    std::printf("],\"interpolant\":");
    Points(Sampled(interpolant.Curve, 1200));
    std::printf("}");
}

void ApproximationSearch()
{
    const std::vector<Point3> points = NoisyPoints();
    constexpr double budget = 0.08;
    std::printf("\"approximation\":{\"budget\":");
    Number(budget);
    std::printf(",\"counts\":[");
    std::vector<double> deviations;
    for (std::size_t count = 4; count <= points.size(); ++count) {
        std::printf(count == 4 ? "%zu" : ",%zu", count);
        deviations.push_back(count == points.size()
                                 ? Interpolate(points, 3)->MaxDeviation
                                 : FitLeastSquares(points, 3, count)->MaxDeviation);
    }
    std::printf("],\"deviations\":");
    Numbers(deviations);
    const CurveFit chosen = *Approximate(points, 3, ApproximationBudget{budget});
    std::printf(",\"chosen\":%zu,\"chosenDeviation\":", chosen.Curve.Poles.size());
    Number(chosen.MaxDeviation);
    std::printf("}");
}

/// The collocation matrix N_i(t_k) of a cubic interpolation of 14 points.
void Collocation()
{
    std::vector<Point3> points = NoisyPoints();
    points.resize(14);
    const CurveFit fit = *Interpolate(points, 3);
    const BSplineCurve& c = fit.Curve;
    std::printf("\"collocation\":{\"rows\":[");
    std::vector<double> basis(4);
    for (std::size_t k = 0; k < points.size(); ++k) {
        std::vector<double> row(c.Poles.size(), 0.0);
        const std::size_t span = FindSpan(c.Knots, 3, c.Poles.size(), fit.Parameters[k]);
        BasisFunctions(c.Knots, 3, span, fit.Parameters[k], basis);
        for (std::size_t a = 0; a < 4; ++a) {
            row[span - 3 + a] = basis[a];
        }
        if (k > 0) {
            std::printf(",");
        }
        Numbers(row);
    }
    std::printf("]}");
}

/// A point of a closed cam profile, r = 1 + 0.3 cos 2a + 0.15 sin 3a, with
/// a radial error added.
Point3 CamPoint(const double a, const double noise)
{
    const double r = 1.0 + 0.3 * std::cos(2.0 * a) + 0.15 * std::sin(3.0 * a) + noise;
    return {r * std::cos(a), r * std::sin(a), 0.0};
}

/// Twelve points of the cam at uneven angles, the first at the dent at the
/// top where the profile turns sharply, interpolated as a periodic cubic
/// and as an open cubic through the first point repeated; and fifty noisy
/// samples fitted by periodic least squares and interpolation.
void PeriodicFigure()
{
    std::vector<Point3> points;
    constexpr int few = 12;
    for (int k = 0; k < few; ++k) {
        const double a =
            std::numbers::pi / 2 + 2.0 * std::numbers::pi * (k + 0.35 * std::sin(1.7 * k)) / few;
        points.push_back(CamPoint(a, 0.0));
    }
    const CurveFit closed = *InterpolatePeriodic(points, 3);
    std::vector<Point3> repeated = points;
    repeated.push_back(points.front());
    const CurveFit open = *Interpolate(repeated, 3);
    std::vector<Point3> noisy;
    constexpr int many = 50;
    for (int k = 0; k < many; ++k) {
        const double a = 2.0 * std::numbers::pi * k / many;
        noisy.push_back(CamPoint(a, 0.03 * std::sin(23.0 * a) + 0.02 * std::cos(41.0 * a)));
    }
    constexpr std::size_t poleCount = 8;
    const CurveFit lsq = *FitLeastSquaresPeriodic(noisy, 3, poleCount);
    const CurveFit through = *InterpolatePeriodic(noisy, 3);
    std::printf("\"periodic\":{\"points\":");
    Points(points);
    std::printf(",\"closed\":");
    Points(Sampled(closed.Curve, 400));
    std::printf(",\"open\":");
    Points(Sampled(open.Curve, 400));
    std::printf(",\"noisy\":");
    Points(noisy);
    std::printf(",\"poleCount\":%zu,\"deviation\":", poleCount);
    Number(lsq.MaxDeviation);
    std::printf(",\"leastSquares\":");
    Points(Sampled(lsq.Curve, 400));
    // The distinct poles and the first again, so the polygon closes.
    std::printf(",\"poles\":");
    Points(std::span<const Point3>{lsq.Curve.Poles}.first(poleCount + 1));
    std::printf(",\"interpolant\":");
    Points(Sampled(through.Curve, 1200));
    std::printf("}");
}

/// A grid on a sphere from pole to pole, and the interpolating surface.
void SurfaceGrid()
{
    constexpr int along = 7;
    constexpr int around = 9;
    std::vector<Point3> grid;
    for (int k = 0; k < along; ++k) {
        const double latitude = -std::numbers::pi / 2 + k * std::numbers::pi / (along - 1);
        for (int l = 0; l < around; ++l) {
            const double longitude = l * std::numbers::pi / (around - 1);
            grid.push_back({std::cos(latitude) * std::cos(longitude),
                            std::cos(latitude) * std::sin(longitude), std::sin(latitude)});
        }
    }
    const SurfaceFit fit = *InterpolateSurface(grid, along, around, 3, 3);
    std::printf("\"surface\":{\"grid\":");
    Points(grid);
    std::printf(",\"deviation\":");
    Number(fit.MaxDeviation);
    std::printf(",\"mesh\":[");
    constexpr int n = 24;
    for (int i = 0; i <= n; ++i) {
        std::vector<Point3> row;
        for (int j = 0; j <= n; ++j) {
            row.push_back(
                Eval(fit.Surface, Point2{i / static_cast<double>(n), j / static_cast<double>(n)}));
        }
        if (i > 0) {
            std::printf(",");
        }
        Points(row);
    }
    std::printf("]}");
}

/// The curvature comb of a uniformly parametrised interpolant, and stations
/// at equal steps of its parameter and of its arc length.
void PropertiesFigure()
{
    const BSplineCurve curve = Interpolate(TurningPoints(), 3, Parametrization::Uniform)->Curve;
    const Interval d = Domain(curve);
    const double length = ArcLength(curve, d)->Length;
    // Teeth at equal steps of arc length, the longest 1.5 mm long.
    constexpr int teeth = 240;
    constexpr double longestTooth = 1.5;
    std::vector<double> parameters;
    double kappaMax = 0.0;
    for (int i = 0; i <= teeth; ++i) {
        const double s = std::min(i * length / teeth, length);
        parameters.push_back(*ParameterAtArcLength(curve, d.GetLower(), s));
        kappaMax = std::max(kappaMax, *Curvature(curve, parameters.back()));
    }
    std::vector<Point3> base;
    std::vector<Point3> tips;
    for (const double t : parameters) {
        const Point3 p = Eval(curve, t);
        base.push_back(p);
        // The tooth points away from the centre of curvature; a straight
        // stretch has no normal and no tooth.
        const auto normal = PrincipalNormal(curve, t);
        const double tooth = longestTooth * *Curvature(curve, t) / kappaMax;
        tips.push_back(normal ? p - tooth * *normal : p);
    }
    constexpr int stations = 12;
    std::vector<Point3> byParameter;
    std::vector<Point3> byLength;
    std::vector<double> lengths;
    for (int i = 0; i <= stations; ++i) {
        const double t = d.GetLower() + i * (d.GetUpper() - d.GetLower()) / stations;
        byParameter.push_back(Eval(curve, std::min(t, d.GetUpper())));
        lengths.push_back(
            ArcLength(curve, Interval{d.GetLower(), std::min(t, d.GetUpper())})->Length);
        const double s = std::min(i * length / stations, length);
        byLength.push_back(Eval(curve, *ParameterAtArcLength(curve, d.GetLower(), s)));
    }
    std::printf("\"properties\":{\"points\":");
    Points(TurningPoints());
    std::printf(",\"curve\":");
    Points(base);
    std::printf(",\"comb\":");
    Points(tips);
    std::printf(",\"length\":");
    Number(length);
    std::printf(",\"by_parameter\":");
    Points(byParameter);
    std::printf(",\"lengths_at_parameter\":");
    Numbers(lengths);
    std::printf(",\"by_length\":");
    Points(byLength);
    std::printf("}");
}

} // namespace

int main()
{
    std::printf("{");
    Parametrizations();
    std::printf(",");
    LeastSquares();
    std::printf(",");
    ApproximationSearch();
    std::printf(",");
    Collocation();
    std::printf(",");
    PeriodicFigure();
    std::printf(",");
    SurfaceGrid();
    std::printf(",");
    PropertiesFigure();
    std::printf("}\n");
    return 0;
}
