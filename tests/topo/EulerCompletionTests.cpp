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
    auto body = TopoTest::MakeHoleSheet();
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
