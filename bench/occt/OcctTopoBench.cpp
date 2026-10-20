// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The OpenCASCADE side of B0.3 to B0.5 (docs/benchmarks.md): the same
// n x n sheet of planar quads as TopoDS faces sharing their edges and
// vertices, built with BRep_Builder and BRepBuilderAPI; the adjacency
// queries through TopExp::MapShapesAndAncestors (the map must be built
// before the first query, so a query session pays for it: the rows report
// the build and the queries together, and the build alone); memory per
// face as the resident set the build adds, measured for OCCT and for
// Stereon in the same process; and one face replaced through
// BRepTools_ReShape, OCCT's way of editing a shape without rebuilding it.
//
// OCCT is LGPL-2.1 and is used here as a reference only (ADR-0008).

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRep_Builder.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <benchmark/benchmark.h>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "ResidentMemory.hpp"

namespace {

/// The ancestor map of TopExp, under OCCT 8's name for it.
using AncestorMap = NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
                                               TopTools_ShapeMapHasher>;

using StereonBench::ResidentBytes;

std::size_t SideFor(const std::size_t faces)
{
    std::size_t n = 1;
    while ((n + 1) * (n + 1) <= faces) {
        ++n;
    }
    return n;
}

/// The n x n sheet as one shell of planar faces sharing edges and vertices.
TopoDS_Shell MakeOcctSheet(const std::size_t n)
{
    const std::size_t side = n + 1;
    std::vector<TopoDS_Vertex> vertices;
    vertices.reserve(side * side);
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertices.push_back(BRepBuilderAPI_MakeVertex(
                gp_Pnt(static_cast<double>(i), static_cast<double>(j), 0.0)));
        }
    }
    const auto at = [&](const std::size_t i, const std::size_t j) -> const TopoDS_Vertex& {
        return vertices[j * side + i];
    };
    std::vector<TopoDS_Edge> horizontal;
    std::vector<TopoDS_Edge> vertical;
    horizontal.reserve(n * side);
    vertical.reserve(n * side);
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            horizontal.push_back(BRepBuilderAPI_MakeEdge(at(i, j), at(i + 1, j)));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertical.push_back(BRepBuilderAPI_MakeEdge(at(i, j), at(i, j + 1)));
        }
    }
    const gp_Pln plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    BRep_Builder builder;
    TopoDS_Shell shell;
    builder.MakeShell(shell);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            BRepBuilderAPI_MakeWire wire(horizontal[j * n + i], vertical[j * side + i + 1],
                                         TopoDS::Edge(horizontal[(j + 1) * n + i].Reversed()),
                                         TopoDS::Edge(vertical[j * side + i].Reversed()));
            builder.Add(shell, BRepBuilderAPI_MakeFace(plane, wire.Wire(), true).Face());
        }
    }
    return shell;
}

/// The faces of every edge, summed, through the ancestor map.
std::size_t FacesOfEveryEdge(const TopoDS_Shape& shape)
{
    AncestorMap map;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, map);
    std::size_t total = 0;
    for (int k = 1; k <= map.Extent(); ++k) {
        total += static_cast<std::size_t>(map.FindFromIndex(k).Extent());
    }
    return total;
}

std::size_t EdgesOfEveryVertex(const TopoDS_Shape& shape)
{
    AncestorMap map;
    TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, map);
    std::size_t total = 0;
    for (int k = 1; k <= map.Extent(); ++k) {
        total += static_cast<std::size_t>(map.FindFromIndex(k).Extent());
    }
    return total;
}

void OcctSheetBuild(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        const TopoDS_Shell shell = MakeOcctSheet(n);
        benchmark::DoNotOptimize(shell.NbChildren());
    }
    state.counters["faces"] = static_cast<double>(n * n);
}

void OcctFacesOfEdges(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    std::size_t total = 0;
    for (auto _ : state) {
        total = FacesOfEveryEdge(shell);
        benchmark::DoNotOptimize(total);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["incidences"] = static_cast<double>(total);
}

void OcctEdgesOfVertices(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    std::size_t total = 0;
    for (auto _ : state) {
        total = EdgesOfEveryVertex(shell);
        benchmark::DoNotOptimize(total);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["incidences"] = static_cast<double>(total);
}

/// The map alone, the part of a query session Stereon does not pay.
void OcctAncestorMap(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    for (auto _ : state) {
        AncestorMap map;
        TopExp::MapShapesAndAncestors(shell, TopAbs_EDGE, TopAbs_FACE, map);
        benchmark::DoNotOptimize(map.Extent());
    }
    state.counters["faces"] = static_cast<double>(n * n);
}

/// Resident bytes per face the sheet adds: the measure both kernels can
/// give in one process. The counters carry the result; the timed loop is
/// the build.
void OcctResidentPerFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    double perFace = 0.0;
    for (auto _ : state) {
        const std::size_t before = ResidentBytes();
        const TopoDS_Shell shell = MakeOcctSheet(n);
        const std::size_t after = ResidentBytes();
        benchmark::DoNotOptimize(shell.NbChildren());
        perFace = static_cast<double>(after - before) / static_cast<double>(n * n);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["resident_bytes_per_face"] = perFace;
}

/// The Stereon sheet of bench/TopoBench.cpp, here for the resident-set
/// comparison in the same process.
Stereon::Shape MakeStereonSheet(const std::size_t n)
{
    using Stereon::EdgeCurveKind;
    using Stereon::EdgeUseSense;
    using Stereon::Point2;
    using Stereon::Point3;
    using Stereon::Vec2;
    using Stereon::Vec3;
    Stereon::ShapeBuilder b;
    const Stereon::SurfaceRef plane = *b.AddSurface(
        Stereon::Surface(Stereon::Plane{Point3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    const std::size_t side = n + 1;
    std::vector<Stereon::BuilderVertexId> vertices;
    vertices.reserve(side * side);
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertices.push_back(
                *b.AddVertex(Point3{static_cast<double>(i), static_cast<double>(j), 0.0}, 1e-7));
        }
    }
    const auto at = [&](const std::size_t i, const std::size_t j) {
        return vertices[j * side + i];
    };
    std::vector<Stereon::BuilderEdgeId> horizontal;
    std::vector<Stereon::BuilderEdgeId> vertical;
    const auto edge = [&](const Stereon::BuilderVertexId from, const Stereon::BuilderVertexId to,
                          const Point3& origin, const Vec3& direction) {
        return *b.AddEdge({.Kind = EdgeCurveKind::Exact,
                           .Start = from,
                           .End = to,
                           .Curve = *b.AddCurve(Stereon::Curve(Stereon::Line{origin, direction})),
                           .Range = Stereon::Interval(0.0, 1.0),
                           .Tolerance = 1e-7});
    };
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            horizontal.push_back(edge(at(i, j), at(i + 1, j),
                                      Point3{static_cast<double>(i), static_cast<double>(j), 0.0},
                                      Vec3{1, 0, 0}));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertical.push_back(edge(at(i, j), at(i, j + 1),
                                    Point3{static_cast<double>(i), static_cast<double>(j), 0.0},
                                    Vec3{0, 1, 0}));
        }
    }
    const Stereon::BuilderShellId shell = *b.AddShell();
    const auto pcurve = [&](const Point2& from, const Vec2& direction) {
        return *b.AddPCurve(Stereon::Curve2d(Stereon::Line2d{from, direction}));
    };
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const Stereon::BuilderFaceId face = *b.AddFace(plane);
            const Stereon::BuilderLoopId loop = *b.AddLoop(face);
            const auto x = static_cast<double>(i);
            const auto y = static_cast<double>(j);
            static_cast<void>(*b.AppendCoedge(loop, horizontal[j * n + i], EdgeUseSense::Same,
                                              pcurve(Point2{x, y}, Vec2{1, 0})));
            static_cast<void>(*b.AppendCoedge(loop, vertical[j * side + i + 1], EdgeUseSense::Same,
                                              pcurve(Point2{x + 1, y}, Vec2{0, 1})));
            static_cast<void>(*b.AppendCoedge(loop, horizontal[(j + 1) * n + i],
                                              EdgeUseSense::Opposite,
                                              pcurve(Point2{x, y + 1}, Vec2{1, 0})));
            static_cast<void>(*b.AppendCoedge(loop, vertical[j * side + i], EdgeUseSense::Opposite,
                                              pcurve(Point2{x, y}, Vec2{0, 1})));
            static_cast<void>(*b.AddFaceToShell(shell, face));
        }
    }
    return *Stereon::Freeze(std::move(b));
}

void StereonResidentPerFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    double perFace = 0.0;
    for (auto _ : state) {
        const std::size_t before = ResidentBytes();
        const Stereon::Shape shape = MakeStereonSheet(n);
        const std::size_t after = ResidentBytes();
        benchmark::DoNotOptimize(shape.Count<Stereon::EntityKind::Face>());
        perFace = static_cast<double>(after - before) / static_cast<double>(n * n);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["resident_bytes_per_face"] = perFace;
}

/// One face replaced by a copy on another plane, through BRepTools_ReShape:
/// the shell is rebuilt around the untouched faces.
void OcctReplaceFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    // The face in the middle of the sheet.
    TopoDS_Face middle;
    {
        std::size_t k = 0;
        for (TopExp_Explorer faces(shell, TopAbs_FACE); faces.More(); faces.Next(), ++k) {
            if (k == (n / 2) * n + n / 2) {
                middle = TopoDS::Face(faces.Current());
                break;
            }
        }
    }
    const gp_Pln lifted(gp_Pnt(0, 0, 1), gp_Dir(0, 0, 1));
    double residentDelta = 0.0;
    for (auto _ : state) {
        const std::size_t before = ResidentBytes();
        TopoDS_Wire wire;
        for (TopExp_Explorer wires(middle, TopAbs_WIRE); wires.More(); wires.Next()) {
            wire = TopoDS::Wire(wires.Current());
        }
        const TopoDS_Face replacement = BRepBuilderAPI_MakeFace(lifted, wire, true).Face();
        BRepTools_ReShape reshape;
        reshape.Replace(middle, replacement);
        const TopoDS_Shape edited = reshape.Apply(shell);
        benchmark::DoNotOptimize(edited.NbChildren());
        residentDelta = static_cast<double>(ResidentBytes() - before);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["resident_delta_bytes"] = residentDelta;
}

constexpr std::int64_t Quick = 10'000;
constexpr std::int64_t Body = 1'000'000;

} // namespace

BENCHMARK(OcctSheetBuild)
    ->Name("B0.3/SheetBuild/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctFacesOfEdges)
    ->Name("B0.3/FacesOfEdges/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctEdgesOfVertices)
    ->Name("B0.3/EdgesOfVertices/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctAncestorMap)
    ->Name("B0.3/AncestorMap/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctResidentPerFace)
    ->Name("B0.4/ResidentPerFace/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(StereonResidentPerFace)
    ->Name("B0.4/ResidentPerFace/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctReplaceFace)
    ->Name("B0.5/ReplaceSurface/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
