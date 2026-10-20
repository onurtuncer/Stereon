// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The topology store on a sheet of N planar quads (B0.3 to B0.5 of
// docs/benchmarks.md, the ADR-0004 and ADR-0005 verification at the
// topology level; the OCCT side is bench/occt/OcctTopoBench.cpp):
//
// - building the sheet through ShapeBuilder and freezing it;
// - the adjacency queries over the whole body: the faces of every edge
//   and the edges of every vertex (B0.3);
// - bytes per face of the frozen tables, and the live heap the sheet holds
//   per face, the measure the OCCT twin takes the same way (B0.4);
// - one local edit on the frozen body, opened in a builder and frozen
//   again: replacing one face's surface, deleting one face, each reporting
//   the pages copied as a fraction of the body (B0.5, ADR-0004 asks < 1%);
// - 100 local edits in a row without Compact(), reporting the memory of
//   the body plus the edited copy relative to the body (ADR-0005 asks
//   < 1.2x), then Compact() of the result.
//
// The sheet has (n + 1)^2 vertices, 2 n (n + 1) edges and 4 n^2 coedges
// for n^2 faces; every interior edge is shared by two faces. The size is
// the benchmark argument: 10^4 faces for a quick run, 10^6 for the body
// of the ADRs (about 1 GB of tables), which the reference run uses.
//
// The OCCT reference run compiles this file into stereon_occt_bench beside
// bench/occt/OcctTopoBench.cpp, so the twins run interleaved in one process.

#include "Memory.hpp"
#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

import Stereon.Topo;

namespace {

using Stereon::BuilderEdgeId;
using Stereon::BuilderFaceId;
using Stereon::BuilderLoopId;
using Stereon::BuilderVertexId;
using Stereon::EdgeCurveKind;
using Stereon::EdgeId;
using Stereon::EdgeUseSense;
using Stereon::EntityKind;
using Stereon::FaceId;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::Vec2;
using Stereon::Vec3;
using Stereon::VertexId;

/// The side of the sheet for a face count: the largest n with n^2 <= faces.
std::size_t SideFor(const std::size_t faces)
{
    std::size_t n = 1;
    while ((n + 1) * (n + 1) <= faces) {
        ++n;
    }
    return n;
}

/// A sheet of n x n unit quads in the plane z = 0, one shell of no solid,
/// Exact line edges with exact pcurves.
Shape MakeSheet(const std::size_t n)
{
    ShapeBuilder b;
    const Stereon::SurfaceRef plane = *b.AddSurface(
        Stereon::Surface(Stereon::Plane{Point3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    const std::size_t side = n + 1;
    std::vector<BuilderVertexId> vertices;
    vertices.reserve(side * side);
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertices.push_back(
                *b.AddVertex(Point3{static_cast<double>(i), static_cast<double>(j), 0.0}, 1e-7));
        }
    }
    // Horizontal edges (i, j) -> (i + 1, j), then vertical (i, j) -> (i, j + 1).
    const auto at = [&](const std::size_t i, const std::size_t j) {
        return vertices[j * side + i];
    };
    std::vector<BuilderEdgeId> horizontal;
    std::vector<BuilderEdgeId> vertical;
    horizontal.reserve(n * side);
    vertical.reserve(n * side);
    const auto line = [&](const Point3& from, const Vec3& direction) {
        return *b.AddCurve(Stereon::Curve(Stereon::Line{from, direction}));
    };
    for (std::size_t j = 0; j < side; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            horizontal.push_back(*b.AddEdge({
                .Kind = EdgeCurveKind::Exact,
                .Start = at(i, j),
                .End = at(i + 1, j),
                .Curve = line(Point3{static_cast<double>(i), static_cast<double>(j), 0.0},
                              Vec3{1, 0, 0}),
                .Range = Interval(0.0, 1.0),
                .Tolerance = 1e-7,
            }));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < side; ++i) {
            vertical.push_back(*b.AddEdge({
                .Kind = EdgeCurveKind::Exact,
                .Start = at(i, j),
                .End = at(i, j + 1),
                .Curve = line(Point3{static_cast<double>(i), static_cast<double>(j), 0.0},
                              Vec3{0, 1, 0}),
                .Range = Interval(0.0, 1.0),
                .Tolerance = 1e-7,
            }));
        }
    }
    const Stereon::BuilderShellId shell = *b.AddShell();
    const auto pcurve = [&](const Point2& from, const Vec2& direction) {
        return *b.AddPCurve(Stereon::Curve2d(Stereon::Line2d{from, direction}));
    };
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const BuilderFaceId face = *b.AddFace(plane);
            const BuilderLoopId loop = *b.AddLoop(face);
            const auto x = static_cast<double>(i);
            const auto y = static_cast<double>(j);
            // Counter-clockwise: bottom along, right along, top against, left against.
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

/// The bytes of every table and pool of a snapshot: the pages they hold.
std::size_t TableBytes(const Shape& shape)
{
    const Stereon::Detail::TopologyTables& t = shape.GetSnapshot()->GetTables();
    return t.Vertices.GetPageBytes() + t.Edges.GetPageBytes() + t.Coedges.GetPageBytes() +
           t.Loops.GetPageBytes() + t.Faces.GetPageBytes() + t.Shells.GetPageBytes() +
           t.Solids.GetPageBytes() + t.Curves.GetPageBytes() + t.PCurves.GetPageBytes() +
           t.Surfaces.GetPageBytes();
}

std::size_t PageCount(const Shape& shape)
{
    const Stereon::Detail::TopologyTables& t = shape.GetSnapshot()->GetTables();
    return t.Vertices.GetPageCount() + t.Edges.GetPageCount() + t.Coedges.GetPageCount() +
           t.Loops.GetPageCount() + t.Faces.GetPageCount() + t.Shells.GetPageCount() +
           t.Solids.GetPageCount() + t.Curves.GetPageCount() + t.PCurves.GetPageCount() +
           t.Surfaces.GetPageCount();
}

/// The pages two snapshots share, over every table and pool.
std::size_t SharedPages(const Shape& a, const Shape& b)
{
    const Stereon::Detail::TopologyTables& x = a.GetSnapshot()->GetTables();
    const Stereon::Detail::TopologyTables& y = b.GetSnapshot()->GetTables();
    return SharedPageCount(x.Vertices, y.Vertices) + SharedPageCount(x.Edges, y.Edges) +
           SharedPageCount(x.Coedges, y.Coedges) + SharedPageCount(x.Loops, y.Loops) +
           SharedPageCount(x.Faces, y.Faces) + SharedPageCount(x.Shells, y.Shells) +
           SharedPageCount(x.Solids, y.Solids) + SharedPageCount(x.Curves, y.Curves) +
           SharedPageCount(x.PCurves, y.PCurves) + SharedPageCount(x.Surfaces, y.Surfaces);
}

/// The faces of every edge, summed: B0.3's first query over the body.
std::size_t FacesOfEveryEdge(const Shape& shape)
{
    std::size_t total = 0;
    for (const EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
        total += (*Stereon::FacesOf(shape, e)).size();
    }
    return total;
}

/// The edges of every vertex, summed: B0.3's second query.
std::size_t EdgesOfEveryVertex(const Shape& shape)
{
    std::size_t total = 0;
    for (const VertexId v : *Stereon::AllEntities<EntityKind::Vertex>(shape)) {
        for ([[maybe_unused]] const EdgeId e : *Stereon::EdgesOf(shape, v)) {
            ++total;
        }
    }
    return total;
}

/// The face in the middle of the sheet, far from its boundary.
FaceId MiddleFace(const Shape& shape, const std::size_t n)
{
    const std::uint32_t index = static_cast<std::uint32_t>((n / 2) * n + n / 2);
    return shape.GetSnapshot()->MakeId<EntityKind::Face>(index);
}

void SheetBuild(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        const Shape shape = MakeSheet(n);
        benchmark::DoNotOptimize(shape.Count<EntityKind::Face>());
    }
    state.counters["faces"] = static_cast<double>(n * n);
}

void FacesOfEdges(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const Shape shape = MakeSheet(n);
    std::size_t total = 0;
    for (auto _ : state) {
        total = FacesOfEveryEdge(shape);
        benchmark::DoNotOptimize(total);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["edges"] = static_cast<double>(shape.Count<EntityKind::Edge>());
    state.counters["incidences"] = static_cast<double>(total);
}

void EdgesOfVertices(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const Shape shape = MakeSheet(n);
    std::size_t total = 0;
    for (auto _ : state) {
        total = EdgesOfEveryVertex(shape);
        benchmark::DoNotOptimize(total);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["vertices"] = static_cast<double>(shape.Count<EntityKind::Vertex>());
    state.counters["incidences"] = static_cast<double>(total);
}

void BytesPerFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const Shape shape = MakeSheet(n);
    for (auto _ : state) {
        benchmark::DoNotOptimize(TableBytes(shape));
    }
    const double faces = static_cast<double>(n * n);
    state.counters["faces"] = faces;
    state.counters["bytes_per_face"] = static_cast<double>(TableBytes(shape)) / faces;
    // Where those bytes are, per face: each table's pages, so that the gap to
    // B0.4's target is attributed rather than estimated.
    const Stereon::Detail::TopologyTables& t = shape.GetSnapshot()->GetTables();
    state.counters["vertex_bytes_per_face"] =
        static_cast<double>(t.Vertices.GetPageBytes()) / faces;
    state.counters["edge_bytes_per_face"] = static_cast<double>(t.Edges.GetPageBytes()) / faces;
    state.counters["coedge_bytes_per_face"] = static_cast<double>(t.Coedges.GetPageBytes()) / faces;
    state.counters["loop_face_bytes_per_face"] =
        static_cast<double>(t.Loops.GetPageBytes() + t.Faces.GetPageBytes()) / faces;
    state.counters["curve_bytes_per_face"] = static_cast<double>(t.Curves.GetPageBytes()) / faces;
    state.counters["pcurve_bytes_per_face"] = static_cast<double>(t.PCurves.GetPageBytes()) / faces;
    state.counters["surface_bytes_per_face"] =
        static_cast<double>(t.Surfaces.GetPageBytes()) / faces;
    state.counters["bytes_per_edge"] = static_cast<double>(TableBytes(shape)) /
                                       static_cast<double>(shape.Count<EntityKind::Edge>());
}

/// The live heap the sheet holds, per face and per edge, with the shape
/// alive: the twin of OcctTopoBench.cpp's HeapPerFace, measured the same
/// way. The builder's working state is gone by then, as the OCCT side's
/// vertex and edge lists are. The timed loop is the build.
void HeapPerFace(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    double bytes = 0.0;
    std::size_t edges = 0;
    for (auto _ : state) {
        const std::size_t before = StereonBench::HeapBytes();
        const Shape shape = MakeSheet(n);
        bytes = static_cast<double>(StereonBench::HeapBytes() - before);
        edges = shape.Count<EntityKind::Edge>();
        benchmark::DoNotOptimize(edges);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["heap_bytes_per_face"] = bytes / static_cast<double>(n * n);
    state.counters["heap_bytes_per_edge"] = bytes / static_cast<double>(edges);
}

/// One edit of the frozen body: opened, edited, frozen again. The pages
/// copied are the pages of the result not shared with the body; the heap
/// delta is what the edited snapshot adds, the measure the OCCT twin gives.
template <typename Edit>
void OneEdit(benchmark::State& state, Edit edit)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const Shape shape = MakeSheet(n);
    std::size_t copied = 0;
    double heapDelta = 0.0;
    for (auto _ : state) {
        const std::size_t before = StereonBench::HeapBytes();
        ShapeBuilder b = *ShapeBuilder::Open(shape);
        edit(b, shape, n);
        const Shape edited = *Stereon::Freeze(std::move(b));
        heapDelta = static_cast<double>(StereonBench::HeapBytes() - before);
        copied = PageCount(edited) - SharedPages(shape, edited);
        benchmark::DoNotOptimize(copied);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["pages"] = static_cast<double>(PageCount(shape));
    state.counters["pages_copied"] = static_cast<double>(copied);
    state.counters["copied_fraction"] =
        static_cast<double>(copied) / static_cast<double>(PageCount(shape));
    state.counters["heap_delta_bytes"] = heapDelta;
}

void ReplaceSurface(ShapeBuilder& b, const Shape& shape, const std::size_t n)
{
    const Stereon::SurfaceRef moved = *b.AddSurface(
        Stereon::Surface(Stereon::Plane{Point3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    static_cast<void>(*b.SetSurface(b.Handle(MiddleFace(shape, n)), moved));
}

void DeleteFace(ShapeBuilder& b, const Shape& shape, const std::size_t n)
{
    static_cast<void>(*b.RemoveFace(b.Handle(MiddleFace(shape, n))));
}

void EditReplaceSurface(benchmark::State& state)
{
    OneEdit(state, ReplaceSurface);
}

void EditDeleteFace(benchmark::State& state)
{
    OneEdit(state, DeleteFace);
}

/// 100 local edits in a row, half surface replacements and half face
/// deletions spread over the sheet, each opened on the previous result
/// and frozen, no Compact(): the memory of the body plus the final copy
/// over the body alone, then Compact() of the result.
void EditSequence(benchmark::State& state)
{
    constexpr std::size_t edits = 100;
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    const Shape body = MakeSheet(n);
    const std::size_t bodyBytes = TableBytes(body);
    double ratio = 0.0;
    for (auto _ : state) {
        Shape current = body;
        for (std::size_t k = 0; k < edits; ++k) {
            ShapeBuilder b = *ShapeBuilder::Open(current);
            // Faces spread along a diagonal, never one twice.
            const std::size_t i = 1 + (k * 7) % (n - 2);
            const std::size_t j = 1 + (k * 11) % (n - 2);
            const FaceId face = current.GetSnapshot()->MakeId<EntityKind::Face>(
                static_cast<std::uint32_t>(j * n + i));
            if (k % 2 == 0) {
                const Stereon::SurfaceRef moved = *b.AddSurface(Stereon::Surface(Stereon::Plane{
                    Point3{0, 0, static_cast<double>(k)}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
                static_cast<void>(*b.SetSurface(b.Handle(face), moved));
            } else {
                static_cast<void>(*b.RemoveFace(b.Handle(face)));
            }
            current = *Stereon::Freeze(std::move(b));
        }
        // What the two snapshots hold together: the body's pages plus the
        // pages the edits copied.
        const std::size_t together = bodyBytes + (PageCount(current) - SharedPages(body, current)) *
                                                     (bodyBytes / PageCount(body));
        ratio = static_cast<double>(together) / static_cast<double>(bodyBytes);
        benchmark::DoNotOptimize(ratio);
    }
    state.counters["faces"] = static_cast<double>(n * n);
    state.counters["memory_ratio"] = ratio;
}

void CompactAfterEdits(benchmark::State& state)
{
    const std::size_t n = SideFor(static_cast<std::size_t>(state.range(0)));
    Shape current = MakeSheet(n);
    for (std::size_t k = 0; k < 50; ++k) {
        ShapeBuilder b = *ShapeBuilder::Open(current);
        const std::size_t i = 1 + (k * 7) % (n - 2);
        const std::size_t j = 1 + (k * 11) % (n - 2);
        static_cast<void>(*b.RemoveFace(b.Handle(current.GetSnapshot()->MakeId<EntityKind::Face>(
            static_cast<std::uint32_t>(j * n + i)))));
        current = *Stereon::Freeze(std::move(b));
    }
    for (auto _ : state) {
        const Stereon::Compacted compacted = *Stereon::Compact(current);
        benchmark::DoNotOptimize(compacted.Shape.Count<EntityKind::Face>());
    }
    state.counters["faces"] = static_cast<double>(n * n);
}

constexpr std::int64_t Quick = 10'000;
constexpr std::int64_t Body = 1'000'000;

} // namespace

BENCHMARK(SheetBuild)
    ->Name("B0.3/SheetBuild/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(FacesOfEdges)
    ->Name("B0.3/FacesOfEdges/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(EdgesOfVertices)
    ->Name("B0.3/EdgesOfVertices/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BytesPerFace)->Name("B0.4/BytesPerFace/Stereon")->Arg(Quick)->Arg(Body);
BENCHMARK(HeapPerFace)
    ->Name("B0.4/HeapPerFace/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(EditReplaceSurface)
    ->Name("B0.5/ReplaceSurface/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(EditDeleteFace)
    ->Name("B0.5/DeleteFace/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(EditSequence)
    ->Name("ADR-0005/EditSequence/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(CompactAfterEdits)
    ->Name("ADR-0004/Compact/Stereon")
    ->Arg(Quick)
    ->Arg(Body)
    ->Unit(benchmark::kMillisecond);
