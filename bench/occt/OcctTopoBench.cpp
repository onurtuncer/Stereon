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
// face as the live heap the sheet holds (bench/Memory.hpp); and one face
// replaced or removed through BRepTools_ReShape, OCCT's way of editing a
// shape without rebuilding it. The Stereon twins are bench/TopoBench.cpp,
// compiled into this executable too, so that both libraries' rows run
// interleaved in one process.
//
// OCCT is LGPL-2.1 and is used here as a reference only (ADR-0008).

#include "../Memory.hpp"
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

namespace {

/// The ancestor map of TopExp, under OCCT 8's name for it.
using AncestorMap = NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
                                               TopTools_ShapeMapHasher>;

using StereonBench::HeapBytes;

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

/// The live heap the sheet holds, per face and per edge, with the shell
/// alive: the B0.4 measure, taken the same way for both kernels
/// (TopoBench.cpp's HeapPerFace). The counters carry the result; the timed
/// loop is the build.
void OcctHeapPerFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    double bytes = 0.0;
    for (auto _ : state) {
        const std::size_t before = HeapBytes();
        const TopoDS_Shell shell = MakeOcctSheet(n);
        bytes = static_cast<double>(HeapBytes() - before);
        benchmark::DoNotOptimize(shell.NbChildren());
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["heap_bytes_per_face"] = bytes / static_cast<double>(n * n);
    state.counters["heap_bytes_per_edge"] = bytes / static_cast<double>(2 * n * (n + 1));
}

/// The face in the middle of the sheet, as TopoBench.cpp's MiddleFace.
TopoDS_Face MiddleFace(const TopoDS_Shell& shell, const std::size_t n)
{
    std::size_t k = 0;
    for (TopExp_Explorer faces(shell, TopAbs_FACE); faces.More(); faces.Next(), ++k) {
        if (k == (n / 2) * n + n / 2) {
            return TopoDS::Face(faces.Current());
        }
    }
    return {};
}

/// One face replaced by a copy on another plane, through BRepTools_ReShape:
/// the shell is rebuilt around the untouched faces.
void OcctReplaceFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    const TopoDS_Face middle = MiddleFace(shell, n);
    const gp_Pln lifted(gp_Pnt(0, 0, 1), gp_Dir(0, 0, 1));
    double heapDelta = 0.0;
    for (auto _ : state) {
        const std::size_t before = HeapBytes();
        TopoDS_Wire wire;
        for (TopExp_Explorer wires(middle, TopAbs_WIRE); wires.More(); wires.Next()) {
            wire = TopoDS::Wire(wires.Current());
        }
        const TopoDS_Face replacement = BRepBuilderAPI_MakeFace(lifted, wire, true).Face();
        BRepTools_ReShape reshape;
        reshape.Replace(middle, replacement);
        const TopoDS_Shape edited = reshape.Apply(shell);
        heapDelta = static_cast<double>(HeapBytes() - before);
        benchmark::DoNotOptimize(edited.NbChildren());
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["heap_delta_bytes"] = heapDelta;
}

/// One face removed through BRepTools_ReShape: the twin of TopoBench.cpp's
/// DeleteFace.
void OcctDeleteFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const TopoDS_Shell shell = MakeOcctSheet(n);
    const TopoDS_Face middle = MiddleFace(shell, n);
    double heapDelta = 0.0;
    for (auto _ : state) {
        const std::size_t before = HeapBytes();
        BRepTools_ReShape reshape;
        reshape.Remove(middle);
        const TopoDS_Shape edited = reshape.Apply(shell);
        heapDelta = static_cast<double>(HeapBytes() - before);
        benchmark::DoNotOptimize(edited.NbChildren());
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["heap_delta_bytes"] = heapDelta;
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
BENCHMARK(OcctHeapPerFace)
    ->Name("B0.4/HeapPerFace/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctReplaceFace)
    ->Name("B0.5/ReplaceSurface/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(OcctDeleteFace)
    ->Name("B0.5/DeleteFace/OCCT")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
