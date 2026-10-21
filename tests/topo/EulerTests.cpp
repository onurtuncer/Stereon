// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Euler operators of ADR-0027 (the Boolean prototype's M5): a
// tetrahedron built by the classical sequence MVFS, MEV, MEF and torn down
// by the inverses KEF, KEV, KVFS to nothing; each inverse restoring the
// loop it touched up to the surviving coedges; refusals that change
// nothing; transactions that keep their candidate only on success; and
// random MEF/KEF and MEV/KEV pairs that leave a valid body as it was.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <random>
#include <utility>
#include <vector>

import Stereon.Topo;

using Stereon::BuilderCoedgeId;
using Stereon::BuilderEdgeId;
using Stereon::BuilderFaceId;
using Stereon::BuilderLoopId;
using Stereon::BuilderVertexId;
using Stereon::EdgeGeometry;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::Point3;
using Stereon::ShapeBuilder;
using Stereon::Vec3;

namespace {

template <typename T>
T Must(Stereon::Result<T> result)
{
    REQUIRE(result.has_value());
    return *std::move(result);
}

void Must(const Stereon::Result<void>& result)
{
    REQUIRE(result.has_value());
}

/// A loop's coedge slots as a cycle, rotated to start at the smallest, so
/// two cycles compare equal whatever their anchors.
std::vector<std::uint32_t> Cycle(const ShapeBuilder& b, const BuilderLoopId loop)
{
    std::vector<std::uint32_t> slots;
    for (const BuilderCoedgeId c : b.CoedgesOf(loop)) {
        slots.push_back(c.GetIndex());
    }
    std::ranges::rotate(slots, std::ranges::min_element(slots));
    return slots;
}

/// Live entities of every kind, and the epoch.
struct Counts {
    std::uint32_t V, E, C, L, F, S;
    std::uint64_t Epoch;
    bool operator==(const Counts&) const = default;
};

Counts CountsOf(const ShapeBuilder& b)
{
    return {b.Count<EntityKind::Vertex>(),
            b.Count<EntityKind::Edge>(),
            b.Count<EntityKind::Coedge>(),
            b.Count<EntityKind::Loop>(),
            b.Count<EntityKind::Face>(),
            b.Count<EntityKind::Shell>(),
            b.GetEpoch()};
}

/// The tetrahedron (v0, v1, v2, v3) of the classical sequence, with the
/// handles each step returned. Spatial edges on one plane: the topology is
/// under test, not the geometry.
struct Tetrahedron {
    ShapeBuilder B;
    EdgeGeometry Edge;
    Stereon::VertexFaceShell Base;
    Stereon::EdgeVertex E1, E2, E4;
    Stereon::EdgeFace E3, E5, E6;
};

EdgeGeometry SpatialEdge(ShapeBuilder& b)
{
    EdgeGeometry g;
    g.Kind = Stereon::EdgeCurveKind::Spatial;
    g.Curve = Must(b.AddCurve(Stereon::Curve(Stereon::Line{Point3{}, Vec3{1.0, 0.0, 0.0}})));
    g.Range = Stereon::Interval{0.0, 1.0};
    return g;
}

/// MVFS, then two MEVs and an MEF: a triangle and its back, a closed shell
/// of two faces.
void MakePillow(Tetrahedron& t)
{
    ShapeBuilder& b = t.B;
    t.Edge = SpatialEdge(b);
    const Stereon::SurfaceRef plane = Must(b.AddSurface(
        Stereon::Surface(Stereon::Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}})));
    t.Base = Must(b.MakeVertexFaceShell(Point3{}, 0.0, plane, false, Must(b.AddSolid())));
    t.E1 = Must(b.MakeEdgeVertex(t.Base.Loop, t.Base.Vertex, Point3{1, 0, 0}, 0.0, t.Edge));
    t.E2 = Must(b.MakeEdgeVertex(t.E1.Out, Point3{0, 1, 0}, 0.0, t.Edge));
    // From v2, where E2's out-use ends, to v0, where E1's back-use ends.
    t.E3 = Must(b.MakeEdgeFace(t.E2.Out, t.E1.Back, t.Edge));
}

Tetrahedron MakeTetrahedron()
{
    Tetrahedron t;
    MakePillow(t);
    ShapeBuilder& b = t.B;
    // The apex, from v2 in the back face, then its edges to v1 and v0.
    t.E4 = Must(b.MakeEdgeVertex(t.E3.New, Point3{0, 0, 1}, 0.0, t.Edge));
    t.E5 = Must(b.MakeEdgeFace(t.E4.Out, t.E2.Back, t.Edge));
    t.E6 = Must(b.MakeEdgeFace(t.E4.Out, t.E1.Back, t.Edge));
    return t;
}

} // namespace

TEST_CASE("MVFS and KVFS make and remove a vertex, face and shell", "[topo][euler]")
{
    ShapeBuilder b;
    const auto plane = Must(b.AddSurface(
        Stereon::Surface(Stereon::Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}})));
    const auto made = Must(b.MakeVertexFaceShell(Point3{1, 2, 3}, 0.0, plane));
    REQUIRE(b.Count<EntityKind::Vertex>() == 1);
    REQUIRE(b.Count<EntityKind::Face>() == 1);
    REQUIRE(b.Count<EntityKind::Loop>() == 1);
    REQUIRE(b.Count<EntityKind::Shell>() == 1);
    REQUIRE(b.CoedgesOf(made.Loop).empty());
    Must(b.KillVertexFaceShell(made.Vertex, made.Face));
    REQUIRE(b.Count<EntityKind::Vertex>() == 0);
    REQUIRE(b.Count<EntityKind::Face>() == 0);
    REQUIRE(b.Count<EntityKind::Loop>() == 0);
    REQUIRE(b.Count<EntityKind::Shell>() == 0);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("the classical sequence builds a valid tetrahedron", "[topo][euler]")
{
    Tetrahedron t = MakeTetrahedron();
    const ShapeBuilder& b = t.B;
    REQUIRE(b.Count<EntityKind::Vertex>() == 4);
    REQUIRE(b.Count<EntityKind::Edge>() == 6);
    REQUIRE(b.Count<EntityKind::Face>() == 4);
    REQUIRE(b.Count<EntityKind::Coedge>() == 12);
    // Three coedges to each face.
    for (const BuilderLoopId l : {t.Base.Loop, t.E3.Loop, t.E5.Loop, t.E6.Loop}) {
        REQUIRE(b.CoedgesOf(l).size() == 3);
    }
    const Stereon::ValidityReport report = b.Validate();
    INFO(report.ToText());
    REQUIRE(report.IsValid());
    const auto shape = Stereon::Freeze(std::move(t.B));
    REQUIRE(shape.has_value());
}

TEST_CASE("each inverse restores the loop it touched, down to nothing", "[topo][euler]")
{
    Tetrahedron t = MakeTetrahedron();
    ShapeBuilder& b = t.B;
    // The back face's loop after each step of the construction, undone in
    // reverse order.
    Must(b.KillEdgeFace(t.E6.Edge, t.E6.Face));
    REQUIRE(Cycle(b, t.E3.Loop).size() == 4);
    Must(b.KillEdgeFace(t.E5.Edge, t.E5.Face));
    const std::vector<std::uint32_t> withStrut = Cycle(b, t.E3.Loop);
    REQUIRE(withStrut.size() == 5);
    Must(b.KillEdgeVertex(t.E4.Edge, t.E4.Vertex));
    std::vector<std::uint32_t> back{t.E3.New.GetIndex(), t.E2.Back.GetIndex(),
                                    t.E1.Back.GetIndex()};
    std::ranges::rotate(back, std::ranges::min_element(back));
    REQUIRE(Cycle(b, t.E3.Loop) == back);

    Must(b.KillEdgeFace(t.E3.Edge, t.E3.Face));
    std::vector<std::uint32_t> base{t.E1.Out.GetIndex(), t.E2.Out.GetIndex(), t.E2.Back.GetIndex(),
                                    t.E1.Back.GetIndex()};
    std::ranges::rotate(base, std::ranges::min_element(base));
    REQUIRE(Cycle(b, t.Base.Loop) == base);
    Must(b.KillEdgeVertex(t.E2.Edge, t.E2.Vertex));
    Must(b.KillEdgeVertex(t.E1.Edge, t.E1.Vertex));
    REQUIRE(b.CoedgesOf(t.Base.Loop).empty());
    Must(b.KillVertexFaceShell(t.Base.Vertex, t.Base.Face));
    for (const std::uint32_t n :
         {b.Count<EntityKind::Vertex>(), b.Count<EntityKind::Edge>(), b.Count<EntityKind::Coedge>(),
          b.Count<EntityKind::Loop>(), b.Count<EntityKind::Face>(), b.Count<EntityKind::Shell>()}) {
        REQUIRE(n == 0);
    }
}

TEST_CASE("an operator that refuses changes nothing", "[topo][euler]")
{
    Tetrahedron t;
    MakePillow(t);
    ShapeBuilder& b = t.B;
    const Counts before = CountsOf(b);
    const auto refused = [&](const ErrorCode code) {
        REQUIRE(code == ErrorCode::InvalidInput);
        REQUIRE(CountsOf(b) == before);
    };
    // v2 has two edges: E2's is no strut to kill.
    refused(b.KillEdgeVertex(t.E2.Edge, t.E2.Vertex).error().GetCode());
    // Two coedges of different loops.
    refused(b.MakeEdgeFace(t.E1.Out, t.E3.New, t.Edge).error().GetCode());
    // A loop that is not empty, and a vertex that is not isolated.
    refused(b.MakeEdgeVertex(t.Base.Loop, t.Base.Vertex, Point3{}, 0.0, t.Edge).error().GetCode());
    // A face that is not the only one of its shell.
    refused(b.KillVertexFaceShell(t.Base.Vertex, t.Base.Face).error().GetCode());
    // A degenerate edge from MEV, and an edge geometry without its curve.
    EdgeGeometry degenerate = t.Edge;
    degenerate.Kind = Stereon::EdgeCurveKind::Degenerate;
    degenerate.Curve = {};
    degenerate.PCurve = Must(
        b.AddPCurve(Stereon::Curve2d(Stereon::Line2d{Stereon::Point2{}, Stereon::Vec2{1.0, 0.0}})));
    const Counts withPCurve = CountsOf(b);
    REQUIRE(withPCurve == before);
    refused(b.MakeEdgeVertex(t.E1.Out, Point3{}, 0.0, degenerate).error().GetCode());
    EdgeGeometry noCurve = t.Edge;
    noCurve.Curve = {};
    refused(b.MakeEdgeVertex(t.E1.Out, Point3{}, 0.0, noCurve).error().GetCode());

    // A strut: both uses in one face, so KEF has no second face to join.
    const auto strut = Must(b.MakeEdgeVertex(t.E1.Out, Point3{2, 0, 0}, 0.0, t.Edge));
    const Counts withStrut = CountsOf(b);
    REQUIRE(b.KillEdgeFace(strut.Edge, t.Base.Face).error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(CountsOf(b) == withStrut);
}

TEST_CASE("a transaction keeps its candidate only on success", "[topo][euler]")
{
    Tetrahedron t;
    MakePillow(t);
    ShapeBuilder& b = t.B;
    const Counts before = CountsOf(b);
    const std::vector<std::uint32_t> loop = Cycle(b, t.Base.Loop);

    // An edit that fails late, after an MEV and an MEF.
    Stereon::EdgeVertex aborted;
    Stereon::CurveRef abortedCurve;
    const auto failed = b.Transaction([&](ShapeBuilder& candidate) -> Stereon::Result<void> {
        abortedCurve = SpatialEdge(candidate).Curve;
        aborted = Must(candidate.MakeEdgeVertex(t.E1.Out, Point3{2, 0, 0}, 0.0, t.Edge));
        Must(candidate.MakeEdgeFace(aborted.Out, t.E2.Out, t.Edge));
        return Stereon::Fail(ErrorCode::DegenerateResult, "a late failure");
    });
    REQUIRE(failed.error().GetCode() == ErrorCode::DegenerateResult);
    Counts after = CountsOf(b);
    // Nothing changed but the epoch, which advanced: cursors made before
    // the failed edit are stale (ADR-0038).
    REQUIRE(after.Epoch > before.Epoch);
    after.Epoch = before.Epoch;
    REQUIRE(after == before);
    REQUIRE(Cycle(b, t.Base.Loop) == loop);

    // The aborted edit's handles stay invalid, even once this builder
    // allocates again where the candidate had (ADR-0040).
    const auto fresh = Must(b.MakeEdgeVertex(t.E1.Out, Point3{3, 0, 0}, 0.0, t.Edge));
    REQUIRE(fresh.Edge != aborted.Edge);
    REQUIRE(fresh.Vertex != aborted.Vertex);
    REQUIRE(b.KillEdgeVertex(aborted.Edge, aborted.Vertex).error().GetCode() ==
            ErrorCode::InvalidInput);
    Must(b.KillEdgeVertex(fresh.Edge, fresh.Vertex));
    REQUIRE(Cycle(b, t.Base.Loop) == loop);

    // So does its geometry: a curve added later takes another index, and
    // the aborted reference is refused.
    const EdgeGeometry later = SpatialEdge(b);
    REQUIRE(later.Curve.GetIndex() != abortedCurve.GetIndex());
    EdgeGeometry stale = t.Edge;
    stale.Curve = abortedCurve;
    REQUIRE(b.MakeEdgeVertex(t.E1.Out, Point3{4, 0, 0}, 0.0, stale).error().GetCode() ==
            ErrorCode::InvalidInput);

    // An edit that succeeds: its handles are this builder's.
    const auto kept = b.Transaction([&](ShapeBuilder& candidate) {
        return candidate.MakeEdgeVertex(t.E1.Out, Point3{2, 0, 0}, 0.0, t.Edge);
    });
    REQUIRE(kept.has_value());
    REQUIRE(b.Count<EntityKind::Edge>() == before.E + 1);
    Must(b.KillEdgeVertex(kept->Edge, kept->Vertex));
    REQUIRE(Cycle(b, t.Base.Loop) == loop);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("random operator and inverse pairs leave the body as it was", "[topo][euler]")
{
    Tetrahedron t = MakeTetrahedron();
    ShapeBuilder& b = t.B;
    const std::vector<BuilderLoopId> loops{t.Base.Loop, t.E3.Loop, t.E5.Loop, t.E6.Loop};
    std::mt19937 random(20261020);
    for (int trial = 0; trial < 200; ++trial) {
        const BuilderLoopId loop = loops[random() % loops.size()];
        const std::vector<BuilderCoedgeId> coedges = b.CoedgesOf(loop);
        const std::vector<std::uint32_t> before = Cycle(b, loop);
        const BuilderCoedgeId a = coedges[random() % coedges.size()];
        if (random() % 2 == 0) {
            const auto strut = Must(b.MakeEdgeVertex(a, Point3{0, 0, 2}, 0.0, t.Edge));
            REQUIRE(b.Validate().IsValid());
            Must(b.KillEdgeVertex(strut.Edge, strut.Vertex));
        } else {
            BuilderCoedgeId c = coedges[random() % coedges.size()];
            if (c == a) {
                c = coedges[(static_cast<std::size_t>(std::ranges::find(coedges, a) -
                                                      coedges.begin()) +
                             1) %
                            coedges.size()];
            }
            const auto split = Must(b.MakeEdgeFace(a, c, t.Edge));
            REQUIRE(b.Validate().IsValid());
            Must(b.KillEdgeFace(split.Edge, split.Face));
        }
        REQUIRE(Cycle(b, loop) == before);
    }
    REQUIRE(b.Count<EntityKind::Face>() == 4);
    REQUIRE(b.Count<EntityKind::Edge>() == 6);
    REQUIRE(b.Validate().IsValid());
}

TEST_CASE("throwing transactions reserve aborted identities without publishing edits",
          "[topo][euler]")
{
    Tetrahedron t;
    MakePillow(t);
    ShapeBuilder& b = t.B;
    const auto before = CountsOf(b);
    const auto loop = Cycle(b, t.Base.Loop);
    Stereon::EdgeVertex aborted;
    Stereon::CurveRef curve;
    Stereon::PCurveRef pcurve;
    Stereon::SurfaceRef surface;
    const auto edit = [&](ShapeBuilder& candidate) -> Stereon::Result<void> {
        curve = SpatialEdge(candidate).Curve;
        pcurve = Must(candidate.AddPCurve(Stereon::Curve2d(Stereon::Line2d{})));
        surface = Must(candidate.AddSurface(Stereon::Surface(Stereon::Plane{})));
        aborted = Must(candidate.MakeEdgeVertex(t.E1.Out, Point3{2, 0, 0}, 0.0, t.Edge));
        throw std::bad_alloc{};
    };
    REQUIRE_THROWS_AS(b.Transaction(edit), std::bad_alloc);
    auto after = CountsOf(b);
    REQUIRE(after.Epoch > before.Epoch);
    after.Epoch = before.Epoch;
    REQUIRE(after == before);
    REQUIRE(Cycle(b, t.Base.Loop) == loop);

    // Reservations survive a move, and both entity and geometry indices
    // remain unavailable after a thrown failure.
    ShapeBuilder moved(std::move(b));
    const auto fresh = Must(moved.MakeEdgeVertex(t.E1.Out, Point3{3, 0, 0}, 0.0, t.Edge));
    REQUIRE(fresh.Vertex.GetIndex() > aborted.Vertex.GetIndex());
    REQUIRE(fresh.Edge.GetIndex() > aborted.Edge.GetIndex());
    REQUIRE(fresh.Out.GetIndex() > aborted.Back.GetIndex());
    REQUIRE(SpatialEdge(moved).Curve.GetIndex() > curve.GetIndex());
    REQUIRE(Must(moved.AddPCurve(Stereon::Curve2d(Stereon::Line2d{}))).GetIndex() >
            pcurve.GetIndex());
    REQUIRE(Must(moved.AddSurface(Stereon::Surface(Stereon::Plane{}))).GetIndex() >
            surface.GetIndex());
    REQUIRE_FALSE(moved.KillEdgeVertex(aborted.Edge, aborted.Vertex));
    REQUIRE_FALSE(moved.AddFace(surface));
    REQUIRE(moved.Validate().IsValid());
}

TEST_CASE("nested failed transactions never reissue an aborted slot", "[topo][euler][property]")
{
    ShapeBuilder b;
    std::vector<BuilderVertexId> aborted;
    std::mt19937 random(20261021);
    for (int trial = 0; trial < 32; ++trial) {
        const auto epoch = b.GetEpoch();
        const auto failed = b.Transaction([&](ShapeBuilder& outer) -> Stereon::Result<void> {
            const auto inner = [&](ShapeBuilder& candidate) -> Stereon::Result<void> {
                const auto count = 1 + random() % 8;
                for (unsigned i = 0; i < count; ++i) {
                    aborted.push_back(Must(candidate.AddVertex(Point3{}, 0.0)));
                }
                throw std::bad_alloc{};
            };
            REQUIRE_THROWS_AS(outer.Transaction(inner), std::bad_alloc);
            return Stereon::Fail(ErrorCode::InvalidInput, "discard outer too");
        });
        REQUIRE_FALSE(failed);
        REQUIRE(b.GetEpoch() > epoch);
        const auto fresh = Must(b.AddVertex(Point3{}, 0.0));
        for (const auto old : aborted) {
            REQUIRE(fresh.GetIndex() > old.GetIndex());
            REQUIRE_FALSE(b.RemoveVertex(old));
        }
        Must(b.RemoveVertex(fresh));
        REQUIRE(b.Count<EntityKind::Vertex>() == 0);
        REQUIRE(b.Validate().IsValid());
    }
}
