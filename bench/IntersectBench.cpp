// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The intersection benchmarks of docs/benchmarks.md, Phase 2, Stereon's side;
// the OCCT twins are bench/occt/OcctIntersectBench.cpp, on the cases of
// bench/IntersectCases.hpp:
//
//   B2.1  curve-curve on corpus/planar-curves/ (lines and conics; splines and
//         rational arcs) against Geom2dAPI_InterCurveCurve, and lines through
//         a rational NURBS sphere and cylinder against GeomAPI_IntCS;
//   B2.2  surface-surface on the SSI corpus, nurbs-ssi and curved-ssi by
//         family, against GeomAPI_IntSS;
//   B2.3  the Tier A pairs of planes, spheres and cylinders of analytic-ssi,
//         by pair, against IntAna_QuadQuadGeo.
//
// An iteration intersects every case of the row once; the time per item is
// the time per case. The counters score the row's answers once, outside the
// timed loop, against the closed forms (IntersectCases.hpp). The Benchmarks
// workflow runs this file in stereon_bench for the trend; the OCCT reference
// run compiles it into stereon_occt_bench beside its twins, so that both
// libraries run interleaved in one process.

#include "Json.hpp"
#include <benchmark/benchmark.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "IntersectCases.hpp"

namespace {

namespace B2 = StereonBench::B2;

/// A row's cases and score, made on first use: benchmark functions run many
/// times (each repetition, and the iteration estimate before it).
template <typename Cases>
struct Row {
    Cases Items;
    B2::Score Score;
};

const Row<std::vector<B2::PlanarCase>>& PlanarRow(const bool splines)
{
    static std::map<bool, Row<std::vector<B2::PlanarCase>>> rows;
    if (const auto found = rows.find(splines); found != rows.end()) {
        return found->second;
    }
    Row<std::vector<B2::PlanarCase>> row{B2::PlanarCases(splines), {}};
    for (const B2::PlanarCase& c : row.Items) {
        row.Score.Cases += 1.0;
        const auto result = Intersect(c.A, c.DomainA, c.B, c.DomainB, B2::Ctx());
        if (!result) {
            row.Score.Unresolved += 1.0;
            row.Score.Expected += static_cast<double>(c.Roots.size());
            continue;
        }
        std::vector<Stereon::Point2> points;
        for (const Stereon::PlanarPoint& p : result->Points) {
            points.push_back(p.Position);
        }
        B2::ScoreRoots(c.Roots, c.Within, points, row.Score);
    }
    return rows.emplace(splines, std::move(row)).first->second;
}

void CurveCurve(benchmark::State& state, const bool splines)
{
    const auto& row = PlanarRow(splines);
    for (auto _ : state) {
        for (const B2::PlanarCase& c : row.Items) {
            auto result = Intersect(c.A, c.DomainA, c.B, c.DomainB, B2::Ctx());
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(row.Items.size()));
    B2::Report(state, row.Score);
}

struct LineRow {
    Stereon::Surface Surface;
    std::vector<B2::LineCase> Lines;
    B2::Score Score;
};

const LineRow& LineRowOf(const bool sphere)
{
    static std::map<bool, LineRow> rows;
    if (const auto found = rows.find(sphere); found != rows.end()) {
        return found->second;
    }
    LineRow row{
        Stereon::Surface{sphere ? B2::CorpusSpline("nurbs-ssi/sphere-plane/circle-origin-x1", 0)
                                : B2::CorpusSpline("nurbs-ssi/cylinder-plane/circle-origin-x1", 0)},
        sphere ? B2::LinesThroughSphere() : B2::LinesThroughCylinder(),
        {}};
    for (const B2::LineCase& line : row.Lines) {
        row.Score.Cases += 1.0;
        const auto result = Intersect(Stereon::Curve{line.Line}, row.Surface, B2::Ctx());
        if (!result) {
            row.Score.Unresolved += 1.0;
            row.Score.Expected += static_cast<double>(line.Roots.size());
            continue;
        }
        std::vector<Stereon::Point3> points;
        for (const Stereon::CurveSurfacePoint& p : result->Points) {
            points.push_back(p.Position);
        }
        B2::ScoreRoots(line.Roots, std::vector<double>(line.Roots.size(), B2::RootReach), points,
                       row.Score);
    }
    return rows.emplace(sphere, std::move(row)).first->second;
}

void CurveSurface(benchmark::State& state, const bool sphere)
{
    const LineRow& row = LineRowOf(sphere);
    std::vector<Stereon::Curve> curves;
    for (const B2::LineCase& line : row.Lines) {
        curves.emplace_back(line.Line);
    }
    for (auto _ : state) {
        for (const Stereon::Curve& curve : curves) {
            auto result = Intersect(curve, row.Surface, B2::Ctx());
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(curves.size()));
    B2::Report(state, row.Score);
}

/// The surface-pair rows by name: "nurbs/<family>", "curved/<family>" or
/// "analytic/<pair>".
const Row<std::vector<B2::SsiCase>>& SsiRow(const std::string& name)
{
    static std::map<std::string, Row<std::vector<B2::SsiCase>>> rows;
    if (const auto found = rows.find(name); found != rows.end()) {
        return found->second;
    }
    const std::string kind = name.substr(0, name.find('/'));
    const std::string which = name.substr(name.find('/') + 1);
    std::vector<B2::SsiCase> cases =
        kind == "analytic" ? B2::PairCases(which) : B2::FramedCases(kind + "-ssi", which);
    const B2::Score score = B2::ScoreSsi(cases, B2::StereonAnswers(cases));
    return rows.emplace(name, Row<std::vector<B2::SsiCase>>{std::move(cases), score}).first->second;
}

void SurfaceSurface(benchmark::State& state, const std::string& name)
{
    const auto& row = SsiRow(name);
    for (auto _ : state) {
        for (const B2::SsiCase& c : row.Items) {
            auto result = Intersect(c.A, c.BoxA, c.B, c.BoxB, c.Tolerance, B2::Ctx());
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(row.Items.size()));
    B2::Report(state, row.Score);
}

/// Every B2 row, registered by name at start-up.
void Register(const std::string& name, void (*row)(benchmark::State&, bool), const bool which)
{
    benchmark::RegisterBenchmark(name.c_str(), row, which);
}

void Register(const std::string& name, const std::string& cases)
{
    benchmark::RegisterBenchmark(name.c_str(), SurfaceSurface, cases);
}

[[maybe_unused]] const bool registered = [] {
    Register("B2.1/CurveCurve/conics/Stereon", CurveCurve, false);
    Register("B2.1/CurveCurve/splines/Stereon", CurveCurve, true);
    Register("B2.1/CurveSurface/line-nurbs-sphere/Stereon", CurveSurface, true);
    Register("B2.1/CurveSurface/line-nurbs-cylinder/Stereon", CurveSurface, false);
    for (const std::string& family : B2::NurbsFamilies) {
        Register("B2.2/SSI/nurbs-" + family + "/Stereon", "nurbs/" + family);
    }
    for (const std::string& family : B2::CurvedFamilies) {
        Register("B2.2/SSI/curved-" + family + "/Stereon", "curved/" + family);
    }
    for (const std::string& pair : B2::AnalyticPairs) {
        Register("B2.3/AnalyticSSI/" + pair + "/Stereon", "analytic/" + pair);
    }
    return true;
}();

} // namespace
