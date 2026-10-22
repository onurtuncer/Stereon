// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "BodiesStd.hpp"

#include <array>
#include <limits>
#include <random>
import Stereon.Topo;
#include "Bodies.hpp"

using namespace Stereon;
using TopoTest::Must;

namespace {
TopoTest::Body UnownedHoleSheet()
{
    TopoTest::Body body;
    const auto face = body.AddFace(TopoTest::XyPlane());
    const std::array<Point3, 4> outer{{{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 3, 0}}};
    const std::array<Point3, 4> inner{{{1, 1, 0}, {1, 2, 0}, {2, 2, 0}, {2, 1, 0}}};
    TopoTest::AddBoundaryPolygon(body, face, outer);
    TopoTest::AddBoundaryPolygon(body, face, inner);
    return body;
}
} // namespace

TEST_CASE("M5 split and join preserve every nonmanifold use", "[topo][euler][m5]")
{
    for (const std::size_t valence : {1U, 2U, 3U, 7U}) {
        auto body = TopoTest::MakeFan(valence);
        auto& b = body.Builder;
        const auto count = b.Count<EntityKind::Coedge>();
        std::vector<std::vector<BuilderCoedgeId>> cycles;
        for (const auto loop : body.Loops) {
            cycles.push_back(b.CoedgesOf(loop));
        }
        const auto split =
            Must(b.SplitEdgeMakeVertex(body.Edges[0], 0.375, {0, 0, 0.375}, TopoTest::Tight));
        REQUIRE(b.Count<EntityKind::Coedge>() == count + valence);
        REQUIRE(b.Validate().IsValid());
        Must(b.JoinEdgeKillVertex(body.Edges[0], split.Upper));
        REQUIRE(b.Count<EntityKind::Coedge>() == count);
        REQUIRE(b.Validate().IsValid());
        for (std::size_t i = 0; i < cycles.size(); ++i) {
            REQUIRE(b.CoedgesOf(body.Loops[i]) == cycles[i]);
        }
        REQUIRE_FALSE(b.GetTopologyEvents().empty());
    }
}

TEST_CASE("M5 wire splitting has a randomized inverse", "[topo][euler][m5][property]")
{
    std::mt19937 generator(527);
    std::uniform_real_distribution<double> parameter(0.01, 0.99);
    for (int iteration = 0; iteration < 100; ++iteration) {
        auto body = TopoTest::MakeWire(1);
        auto& b = body.Builder;
        const double t = parameter(generator);
        const auto split =
            Must(b.SplitEdgeMakeVertex(body.Edges[0], t, {t, t, 0}, TopoTest::Tight));
        REQUIRE(b.Validate().IsValid());
        Must(b.JoinEdgeKillVertex(body.Edges[0], split.Upper));
        REQUIRE(b.Count<EntityKind::Vertex>() == 2);
        REQUIRE(b.Count<EntityKind::Edge>() == 1);
        REQUIRE(b.Validate().IsValid());
    }
}

TEST_CASE("M5 invalid split parameters and positions leave the builder unchanged",
          "[topo][euler][m5]")
{
    auto body = TopoTest::MakeSquareSheet();
    auto& b = body.Builder;
    const auto epoch = b.GetEpoch();
    for (const double t : {0.0, 1.0, -1.0, std::numeric_limits<double>::infinity()}) {
        REQUIRE_FALSE(b.SplitEdgeMakeVertex(body.Edges[0], t, {t, 0, 0}, TopoTest::Tight));
    }
    REQUIRE_FALSE(b.SplitEdgeMakeVertex(body.Edges[0], 0.5, {0.5, 1, 0}, TopoTest::Tight));
    REQUIRE(b.GetEpoch() == epoch);
    REQUIRE(b.GetTopologyEvents().empty());
    const auto split =
        Must(b.SplitEdgeMakeVertex(body.Edges[0], 0.5, {0.5, 0, 0}, TopoTest::Tight));
    const auto after = b.GetEpoch();
    REQUIRE_FALSE(b.JoinEdgeKillVertex(split.Upper, body.Edges[0]));
    REQUIRE(b.GetEpoch() == after);
    Must(b.JoinEdgeKillVertex(body.Edges[0], split.Upper));
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("M5 ring and hole inverse pairs restore cycles", "[topo][euler][m5]")
{
    auto body = UnownedHoleSheet();
    auto& b = body.Builder;
    const auto outer = b.CoedgesOf(body.Loops[0]);
    const auto inner = b.CoedgesOf(body.Loops[1]);
    const auto face = Must(b.MakeFaceKillRingHole(body.Loops[1]));
    REQUIRE(b.Count<EntityKind::Face>() == 2);
    REQUIRE(b.Validate().IsValid());
    Must(b.KillFaceMakeRingHole(face, body.Faces[0]));
    REQUIRE(b.CoedgesOf(body.Loops[0]) == outer);
    REQUIRE(b.CoedgesOf(body.Loops[1]) == inner);
    EdgeGeometry geometry;
    geometry.Kind = EdgeCurveKind::Boundary;
    geometry.Range = Interval{0, 1};
    geometry.Tolerance = TopoTest::Tight;
    geometry.PCurve = Must(b.AddPCurve(Curve2d(Line2d{{3, 0}, {-2, 2}})));
    const auto bridge = Must(b.MakeEdgeKillRing(outer[0], inner[0], geometry));
    REQUIRE(b.Count<EntityKind::Loop>() == 1);
    REQUIRE(b.Validate().IsValid());
    const auto ring = Must(b.KillEdgeMakeRing(bridge.Edge));
    REQUIRE(b.Count<EntityKind::Loop>() == 2);
    REQUIRE(b.Validate().IsValid());
    REQUIRE(b.CoedgesOf(body.Loops[0]) == outer);
    auto restored = b.CoedgesOf(ring);
    REQUIRE(restored.size() == inner.size());
    REQUIRE(std::is_permutation(restored.begin(), restored.end(), inner.begin(), inner.end()));
}

TEST_CASE("M5 vertex merging uses fixed representatives and sorted inputs", "[topo][euler][m5]")
{
    for (const bool reverse : {false, true}) {
        ShapeBuilder b;
        const auto a = Must(b.AddVertex({0, 0, 0}, 0.001));
        const auto c = Must(b.AddVertex({0.0008, 0, 0}, 0.001));
        const auto d = Must(b.AddVertex({0.0016, 0, 0}, 0.001));
        std::array vertices{a, c, d};
        if (reverse) {
            std::ranges::reverse(vertices);
        }
        const auto merged = Must(b.MergeVertices(vertices, Context::Default()));
        REQUIRE(merged.Entries.size() == 3);
        REQUIRE(merged.Entries[0].Representative == a);
        REQUIRE(merged.Entries[1].Representative == a);
        REQUIRE(merged.Entries[2].Representative == d);
        REQUIRE(b.Count<EntityKind::Vertex>() == 2);
        REQUIRE(b.Validate().IsValid());
    }
}

TEST_CASE("M5 merge policy failures and transactions roll back topology events",
          "[topo][euler][m5]")
{
    ShapeBuilder b;
    const auto a = Must(b.AddVertex({0, 0, 0}, 0.0001));
    const auto c = Must(b.AddVertex({0.0005, 0, 0}, 0.001));
    const std::array vertices{a, c};
    REQUIRE_FALSE(b.MergeVertices(
        vertices, Context::Default().WithToleranceGrowthPolicy(ToleranceGrowth::Fail)));
    REQUIRE(b.Count<EntityKind::Vertex>() == 2);
    REQUIRE(b.GetTopologyEvents().empty());
    const auto result = b.Transaction([&](ShapeBuilder& candidate) -> Result<void> {
        Must(candidate.MergeVertices(vertices, Context::Default()));
        return Fail(ErrorCode::InvalidInput, "injected rollback");
    });
    REQUIRE_FALSE(result);
    REQUIRE(b.Count<EntityKind::Vertex>() == 2);
    REQUIRE(b.GetTopologyEvents().empty());
    Must(b.MergeVertices(vertices, Context::Default()));
    REQUIRE(b.Count<EntityKind::Vertex>() == 1);
    REQUIRE(b.GetTopologyEvents().size() == 1);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("M5 splitting preserves seam sides and closed edge incidence", "[topo][euler][m5]")
{
    for (const bool seam : {false, true}) {
        auto body = TopoTest::MakeCylinder();
        auto& b = body.Builder;
        const auto edge = body.Edges[seam ? 2 : 0];
        const double t = seam ? 2 : std::numbers::pi;
        const Point3 position = seam ? Point3{1.5, 0, 2} : Point3{-1.5, 0, 0};
        const auto split = Must(b.SplitEdgeMakeVertex(edge, t, position, TopoTest::Tight));
        REQUIRE(b.Validate().IsValid());
        Must(b.JoinEdgeKillVertex(edge, split.Upper));
        REQUIRE(b.Count<EntityKind::Edge>() == 3);
        REQUIRE(b.Count<EntityKind::Vertex>() == 2);
        REQUIRE(b.Validate().IsValid());
        const auto published = Must(Commit(std::move(b), Context::Default()));
        REQUIRE_FALSE(published.TopologyEvents.empty());
    }
}

TEST_CASE("M5 merge checks attached endpoint traces and preserves incidence", "[topo][euler][m5]")
{
    auto body = TopoTest::MakeWire(1);
    auto& b = body.Builder;
    // The isolated member joins an existing endpoint without moving its trace.
    const auto near = Must(b.AddVertex({0.001, 0, 0}, 0.001));
    const std::array vertices{body.Vertices[0], near};
    Must(b.MergeVertices(vertices, Context::Default()));
    REQUIRE(b.Count<EntityKind::Vertex>() == 2);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("M5 restricted edge caches and primary sides survive splitting",
          "[topo][euler][m5][cache]")
{
    auto body = TopoTest::MakeHinge();
    Must(body.Builder.CacheEdgeCurve(body.Edges[0]));
    const auto original = std::move(body).Freeze();
    const auto edge = original.Resolve(body.Edges[0]);
    const double bound = CacheBoundOf(original, edge);
    auto b = Must(ShapeBuilder::Open(original));
    const auto lower = b.Handle(edge);
    const auto split = Must(b.SplitEdgeMakeVertex(lower, 0.375, {0, 0.375, 0}, TopoTest::Tight));
    const auto divided = Must(Freeze(std::move(b)));
    const auto low = divided.Resolve(lower);
    const auto high = divided.Resolve(split.Upper);
    REQUIRE(CacheBoundOf(divided, low) == bound);
    REQUIRE(CacheBoundOf(divided, high) == bound);
    REQUIRE(CacheCertificateOf(divided, high) == CacheCertificate::Certified);
    REQUIRE(RangeOf(divided, low) == Interval(0, 0.375));
    REQUIRE(RangeOf(divided, high) == Interval(0.375, 1));
    REQUIRE(RangeOf(original, edge) == Interval(0, 1));
    REQUIRE(ValidateTopology(original).IsValid());
    auto join = Must(ShapeBuilder::Open(divided));
    bool recached = false;
    SECTION("a common certificate is retained") {}
    SECTION("different cache representations are dropped")
    {
        Must(join.CacheEdgeCurve(join.Handle(high)));
        recached = true;
    }
    const auto kept = join.Handle(low);
    Must(join.JoinEdgeKillVertex(kept, join.Handle(high)));
    const auto united = Must(Freeze(std::move(join)));
    REQUIRE(RangeOf(united, united.Resolve(kept)) == Interval(0, 1));
    REQUIRE(CacheCertificateOf(united, united.Resolve(kept)) ==
            (recached ? CacheCertificate::None : CacheCertificate::Certified));
}

TEST_CASE("M5 join refuses a third edge at the join vertex", "[topo][euler][m5]")
{
    auto body = TopoTest::MakeWire(1);
    auto& b = body.Builder;
    const auto split =
        Must(b.SplitEdgeMakeVertex(body.Edges[0], 0.5, {0.5, 0.5, 0}, TopoTest::Tight));
    EdgeSpec branch;
    branch.Kind = EdgeCurveKind::Spatial;
    branch.Start = split.Vertex;
    branch.End = Must(b.AddVertex({0.5, 0.5, 1}, TopoTest::Tight));
    branch.Curve = Must(b.AddCurve(Curve(Line{{0.5, 0.5, 0}, {0, 0, 1}})));
    branch.Range = Interval{0, 1};
    branch.Tolerance = TopoTest::Tight;
    Must(b.AddEdge(branch));
    const auto epoch = b.GetEpoch();
    REQUIRE_FALSE(b.JoinEdgeKillVertex(body.Edges[0], split.Upper));
    REQUIRE(b.GetEpoch() == epoch);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("M5 bridge inverse holds at every pair of loop attachment sites",
          "[topo][euler][m5][property]")
{
    const std::array<Point2, 4> outerEnds{{{3, 0}, {3, 3}, {0, 3}, {0, 0}}};
    const std::array<Point2, 4> innerEnds{{{1, 2}, {2, 2}, {2, 1}, {1, 1}}};
    for (std::size_t a = 0; a < 4; ++a) {
        for (std::size_t c = 0; c < 4; ++c) {
            auto body = UnownedHoleSheet();
            auto& b = body.Builder;
            const auto outer = b.CoedgesOf(body.Loops[0]);
            const auto inner = b.CoedgesOf(body.Loops[1]);
            EdgeGeometry geometry;
            geometry.Kind = EdgeCurveKind::Boundary;
            geometry.Range = Interval{0, 1};
            geometry.Tolerance = TopoTest::Tight;
            geometry.PCurve =
                Must(b.AddPCurve(Curve2d(Line2d{outerEnds[a], innerEnds[c] - outerEnds[a]})));
            const auto bridge = Must(b.MakeEdgeKillRing(outer[a], inner[c], geometry));
            REQUIRE(b.Validate().IsValid());
            const auto ring = Must(b.KillEdgeMakeRing(bridge.Edge));
            REQUIRE(b.Validate().IsValid());
            REQUIRE(b.CoedgesOf(body.Loops[0]) == outer);
            const auto restored = b.CoedgesOf(ring);
            REQUIRE(
                std::is_permutation(restored.begin(), restored.end(), inner.begin(), inner.end()));
            const auto face = Must(b.MakeFaceKillRingHole(ring));
            Must(b.KillFaceMakeRingHole(face, body.Faces[0]));
            REQUIRE(b.Validate().IsValid());
            REQUIRE(b.Count<EntityKind::Face>() == 1);
        }
    }
}

TEST_CASE("M5 endpoint containment failure rolls back earlier merge clusters", "[topo][euler][m5]")
{
    ShapeBuilder b;
    const auto first = Must(b.AddVertex({5, 0, 0}, 0.001));
    const auto firstMember = Must(b.AddVertex({5.0005, 0, 0}, 0.001));
    const auto second = Must(b.AddVertex({0, 0, 0}, 0.001));
    const auto secondMember = Must(b.AddVertex({0.0008, 0, 0}, 0.001));
    EdgeSpec edge;
    edge.Kind = EdgeCurveKind::Spatial;
    edge.Start = secondMember;
    edge.End = Must(b.AddVertex({0.0016, 0, 1}, 0.001));
    edge.Curve = Must(b.AddCurve(Curve(Line{{0.0016, 0, 0}, {0, 0, 1}})));
    edge.Range = Interval{0, 1};
    edge.Tolerance = 0.001;
    Must(b.AddEdge(edge));
    const std::array vertices{first, firstMember, second, secondMember};
    // The trace is within secondMember's tolerance, but outside second's
    // tolerance. The first cluster has already merged in the candidate.
    REQUIRE_FALSE(b.MergeVertices(vertices, Context::Default()));
    REQUIRE(b.Count<EntityKind::Vertex>() == 5);
    REQUIRE(b.GetTopologyEvents().empty());
    REQUIRE(b.Validate().IsValid());
    Must(b.SetTolerance(firstMember, 0.001, "the rolled-back handle remains live"));
}

TEST_CASE("M5 separating a ring refuses to disconnect its shell", "[topo][euler][m5]")
{
    auto body = TopoTest::MakeHoleSheet();
    auto& b = body.Builder;
    const auto epoch = b.GetEpoch();
    const auto loops = b.Count<EntityKind::Loop>();
    REQUIRE_FALSE(b.MakeFaceKillRingHole(body.Loops[1]));
    REQUIRE(b.GetEpoch() == epoch);
    REQUIRE(b.Count<EntityKind::Face>() == 1);
    REQUIRE(b.Count<EntityKind::Loop>() == loops);
    REQUIRE(b.GetTopologyEvents().empty());
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("M5 ring separation preserves an existing path through the shell", "[topo][euler][m5]")
{
    auto body = UnownedHoleSheet();
    const auto shape = std::move(body).Freeze();
    const auto source = shape.Resolve(body.Faces[0]);
    auto b = Must(ShapeBuilder::Open(shape));
    const auto support = Must(b.AddFace(Must(b.AddSurface(SurfaceOf(shape, source)))));
    const auto shell = Must(b.AddShell());
    Must(b.AddFaceToShell(shell, b.Handle(source)));
    Must(b.AddFaceToShell(shell, support));
    // A second face shares both rings. After separation it supplies the
    // edge-connected path between the retained face and the extracted face.
    for (const auto loop : *LoopsOf(shape, source)) {
        const auto added = Must(b.AddLoop(support));
        std::vector<CoedgeId> uses;
        for (const auto use : *CoedgesOf(shape, loop)) {
            uses.push_back(use);
        }
        for (auto it = uses.rbegin(); it != uses.rend(); ++it) {
            const auto pcurve = Must(b.AddPCurve(*PCurveOf(shape, *it)));
            Must(b.AppendCoedge(added, b.Handle(EdgeOf(shape, *it)), EdgeUseSense::Opposite,
                                pcurve));
        }
    }
    REQUIRE(b.Validate().IsValid());
    const auto split = Must(b.MakeFaceKillRingHole(b.Handle(shape.Resolve(body.Loops[1]))));
    REQUIRE(b.Count<EntityKind::Face>() == 3);
    REQUIRE(b.Validate().IsValid());
    Must(b.KillFaceMakeRingHole(split, b.Handle(source)));
    REQUIRE(b.Count<EntityKind::Face>() == 2);
    REQUIRE(b.Validate().IsValid());
}
