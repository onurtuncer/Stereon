// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Shape-producing operations (ADR-0006, ADR-0011, ADR-0030): Commit() and
// CompactShape() publish an OperationResult whose History accounts for every
// entity, whose tolerance changes are the edit's, and whose warnings say what
// the edit lost; nothing is published when validation fails. And the cached
// 3D curves of Intersection and Trim edges (ADR-0009): built from the
// primary side, uncertified, dropped whenever what they were built from
// changes (ADR-0004).

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CacheCertificate;
using Stereon::Context;
using Stereon::CurveRef;
using Stereon::DiagnosticCode;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::HistoryEntry;
using Stereon::HistoryRelation;
using Stereon::OperationResult;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::ToleranceGrowth;

namespace {

/// The largest distance between the cache and the primary side S(c(t)) over
/// `samples` parameters spread evenly over the edge: the regression check of
/// ADR-0009, never a certificate.
double SampledDeviation(const Shape& shape, const EdgeId e, const Stereon::Surface& surface,
                        const Stereon::Curve2d& pcurve, const int samples)
{
    const Stereon::Curve& cache = Deref(Stereon::CachedCurveOf(shape, e));
    const Interval range = Stereon::RangeOf(shape, e);
    double worst = 0.0;
    for (int k = 0; k <= samples; ++k) {
        const double t = k == samples
                             ? range.GetUpper()
                             : range.GetLower() + Width(range) * static_cast<double>(k) / samples;
        worst = std::max(worst, Distance(Eval(cache, t), Eval(surface, Eval(pcurve, t))));
    }
    return worst;
}

/// The id of slot `index` of kind K in `shape`: for a result reached through`r
/// another builder than the one whose handles the test holds.
template <EntityKind K>
Stereon::Id<K> At(const Shape& shape, const std::uint32_t index)
{
    return shape.GetSnapshot()->template MakeId<K>(index);
}

OperationResult MustCommit(ShapeBuilder&& builder, const Context& context = Context::Default())
{
    return Must(Stereon::Commit(std::move(builder), context));
}

/// The error of a Commit() that must fail with `code`.
std::string CommitError(ShapeBuilder&& builder, const Context& context, const ErrorCode code)
{
    const auto result = Stereon::Commit(std::move(builder), context);
    REQUIRE_FALSE(result.has_value());
    const std::string message(result.error().GetMessage());
    INFO(message);
    REQUIRE(result.error().GetCode() == code);
    return message;
}

/// The one entry of an input entity.
HistoryEntry Fate(const OperationResult& result, const EntityKind kind, const std::uint32_t index)
{
    const std::vector<HistoryEntry> entries = result.History.OfSource(kind, index);
    REQUIRE(entries.size() == 1);
    return entries.front();
}

} // namespace

// --- Cached edge curves (ADR-0009) -------------------------------------------------

TEST_CASE("an Intersection edge's cache follows its primary side and is not certified",
          "[topo][cache]")
{
    Body hinge = MakeHinge();
    const CurveRef cache = Must(hinge.Builder.CacheEdgeCurve(hinge.Edges[0]));
    REQUIRE(cache.IsValid());
    const Shape shape = std::move(hinge).Freeze();
    REQUIRE(Stereon::ValidateTopology(shape).IsValid());
    const EdgeId e = shape.Resolve(hinge.Edges[0]);
    REQUIRE(Stereon::CacheCertificateOf(shape, e) == CacheCertificate::Uncertified);
    REQUIRE(std::isinf(Stereon::CacheBoundOf(shape, e)));
    REQUIRE_FALSE(Stereon::CurveOf(shape, e).has_value()); // a cache is not the edge's own curve
    const Stereon::Curve& curve = Deref(Stereon::CachedCurveOf(shape, e));
    REQUIRE(Domain(curve) == Interval(0.0, 1.0));
    // The primary side is the line (0, t, 0) on face 0's plane.
    const Stereon::Curve2d hingeOnA = Segment2({0, 0}, {0, 1});
    REQUIRE(SampledDeviation(shape, e, XyPlane(), hingeOnA, 1000) <= 1e-12);
    REQUIRE(Stereon::Describe(shape, e).find("cached, uncertified") != std::string::npos);
    // The other kinds keep no cache.
    const EdgeId side = shape.Resolve(hinge.Edges[1]);
    REQUIRE_FALSE(Stereon::CachedCurveOf(shape, side).has_value());
    REQUIRE(Stereon::CacheCertificateOf(shape, side) == CacheCertificate::None);
}

TEST_CASE("a cache of a curved primary side passes through every sample", "[topo][cache]")
{
    // Face 0 on a cylinder of radius 1.5 about the y axis, the hinge's
    // pcurve there u = t, v = 0.25: a quarter of a radian of a circle.
    const Stereon::Surface cylinder(
        Stereon::Cylinder{{0, 0, 0}, 1.5, Vec3{1, 0, 0}, Vec3{0, 0, 1}});
    const Stereon::Curve2d arc = Segment2({0, 0.25}, {1, 0.25});
    Body hinge = MakeHinge(cylinder, arc, [](Body&) {});
    static_cast<void>(Must(hinge.Builder.CacheEdgeCurve(hinge.Edges[0])));
    const Shape shape = std::move(hinge).Freeze();
    const EdgeId e = shape.Resolve(hinge.Edges[0]);
    // At the interpolated samples the cache is the primary side...
    REQUIRE(SampledDeviation(shape, e, cylinder, arc, ShapeBuilder::CacheSamples - 1) <= 1e-12);
    // ...and between them close to it: the regression check, not a bound.
    const double between = SampledDeviation(shape, e, cylinder, arc, 1000);
    REQUIRE(between <= 1e-6);
    // Rebuilding gives the same cache, bit for bit (ADR-0009: deterministic).
    Body again = MakeHinge(cylinder, arc, [](Body&) {});
    static_cast<void>(Must(again.Builder.CacheEdgeCurve(again.Edges[0])));
    const Shape rebuilt = std::move(again).Freeze();
    const Stereon::Curve& first = Deref(Stereon::CachedCurveOf(shape, e));
    const Stereon::Curve& second =
        Deref(Stereon::CachedCurveOf(rebuilt, rebuilt.Resolve(again.Edges[0])));
    for (const double t : {0.0, 0.1234, 0.5, 0.987, 1.0}) {
        REQUIRE(Eval(first, t) == Eval(second, t));
    }
}

TEST_CASE("CacheEdgeCurve refuses edges that keep no cache or cannot build one", "[topo][cache]")
{
    Body box = MakeBox();
    const auto refused = box.Builder.CacheEdgeCurve(box.Edges[0]); // Exact
    REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE(std::string(refused.error().GetMessage()).find("Exact edge keeps no cached curve") !=
            std::string::npos);

    // An Intersection edge with no use yet has no primary side.
    ShapeBuilder b;
    const BuilderVertexId v = Must(b.AddVertex({0, 0, 0}, 0.0));
    const BuilderVertexId w = Must(b.AddVertex({1, 0, 0}, 0.0));
    const auto bare = Must(b.AddEdge(
        {.Kind = EdgeCurveKind::Intersection, .Start = v, .End = w, .Range = Interval(0.0, 1.0)}));
    REQUIRE(b.CacheEdgeCurve(bare).error().GetCode() == ErrorCode::InvalidInput);

    // A primary pcurve outside its surface's bounded domain.
    const Stereon::Surface patch(*Stereon::MakeTrimmedSurface(
        XyPlane(), Stereon::Box2{Interval(0.0, 1.0), Interval(0.0, 1.0)}));
    Body outside = MakeHinge(patch, Segment2({0, 0}, {0, 3}), [](Body&) {});
    REQUIRE(outside.Builder.CacheEdgeCurve(outside.Edges[0]).error().GetCode() ==
            ErrorCode::InvalidInput);
}

TEST_CASE("a cache is dropped whenever a pcurve, the primary or a surface of its edge changes",
          "[topo][cache][edit]")
{
    const auto cached = [](Body& body) {
        static_cast<void>(Must(body.Builder.CacheEdgeCurve(body.Edges[0])));
    };
    SECTION("a new use of the edge brings a new pcurve")
    {
        // The cache is built after face 0 and before face 1 uses the hinge.
        Body hinge = MakeHinge(XyPlane(), Segment2({0, 0}, {0, 1}), cached);
        const Shape shape = std::move(hinge).Freeze();
        const EdgeId e = shape.Resolve(hinge.Edges[0]);
        REQUIRE_FALSE(Stereon::CachedCurveOf(shape, e).has_value());
        REQUIRE(Stereon::CacheCertificateOf(shape, e) == CacheCertificate::None);
    }

    Body hinge = MakeHinge();
    cached(hinge);
    const Shape base = std::move(hinge).Freeze();
    const EdgeId e = base.Resolve(hinge.Edges[0]);
    REQUIRE(Stereon::CachedCurveOf(base, e).has_value());
    const auto requireDropped = [&](ShapeBuilder&& edit) {
        const OperationResult result = MustCommit(std::move(edit));
        REQUIRE_FALSE(Stereon::CachedCurveOf(result.Shape, At<EntityKind::Edge>(result.Shape, 0))
                          .has_value());
        REQUIRE(result.Warnings.size() == 1);
        REQUIRE(result.Warnings[0].Code == DiagnosticCode::CacheDropped);
        REQUIRE(result.Warnings[0].Entities.at(0).find("Edge#0") != std::string::npos);
        // The source keeps its cache: snapshots are immutable.
        REQUIRE(Stereon::CachedCurveOf(base, e).has_value());
    };
    SECTION("the primary passes to the other side")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(base));
        Must(edit.SetPrimaryCoedge(edit.Handle(e), edit.Handle(base.Resolve(hinge.Coedges[4]))));
        requireDropped(std::move(edit));
    }
    SECTION("a face of the edge takes another surface")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(base));
        const Stereon::SurfaceRef moved = Must(edit.AddSurface(
            Stereon::Surface(Stereon::Plane{{0, 0, 2}, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
        Must(edit.SetSurface(edit.Handle(base.Resolve(hinge.Faces[1])), moved));
        requireDropped(std::move(edit));
    }
    SECTION("a use of the edge goes")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(base));
        Must(edit.RemoveFace(edit.Handle(base.Resolve(hinge.Faces[1]))));
        // The hinge, left with one use, is now the boundary of face 0.
        Stereon::Testing::TopologyFaults::SetKind(edit, edit.Handle(e), EdgeCurveKind::Boundary);
        for (const std::size_t side : {std::size_t{4}, std::size_t{5}, std::size_t{6}}) {
            Must(edit.RemoveEdge(edit.Handle(base.Resolve(hinge.Edges[side]))));
        }
        for (const std::size_t vertex : {std::size_t{4}, std::size_t{5}}) {
            Must(edit.RemoveVertex(edit.Handle(base.Resolve(hinge.Vertices[vertex]))));
        }
        requireDropped(std::move(edit));
    }
    SECTION("what the cache was not built from may change")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(base));
        Must(edit.SetTolerance(edit.Handle(e), 1e-6, "a looser hinge"));
        const OperationResult result = MustCommit(std::move(edit));
        REQUIRE(Stereon::CachedCurveOf(result.Shape, At<EntityKind::Edge>(result.Shape, 0))
                    .has_value());
    }
}

TEST_CASE("Compact keeps an edge's cache", "[topo][cache][compact]")
{
    Body hinge = MakeHinge();
    static_cast<void>(Must(hinge.Builder.CacheEdgeCurve(hinge.Edges[0])));
    const Shape shape = std::move(hinge).Freeze();
    const OperationResult compacted = Must(Stereon::CompactShape(shape, Context::Default()));
    const std::vector<HistoryEntry> hingeFate = compacted.History.OfSource(EntityKind::Edge, 0);
    REQUIRE(hingeFate.size() == 1);
    const EdgeId after =
        compacted.Shape.GetSnapshot()->MakeId<EntityKind::Edge>(hingeFate[0].Target);
    const Stereon::Curve& cache = Deref(Stereon::CachedCurveOf(compacted.Shape, after));
    REQUIRE(Stereon::CacheCertificateOf(compacted.Shape, after) == CacheCertificate::Uncertified);
    REQUIRE(Eval(cache, 0.3) ==
            Eval(Deref(Stereon::CachedCurveOf(shape, shape.Resolve(hinge.Edges[0]))), 0.3));
}

// --- OperationResult (ADR-0006) ----------------------------------------------------

TEST_CASE("a new body is all Generated, checked and published together", "[topo][operation]")
{
    Body hinge = MakeHinge();
    const OperationResult result = MustCommit(std::move(hinge.Builder));
    REQUIRE(Stereon::ValidateTopology(result.Shape).IsValid());
    REQUIRE(result.Checked == Stereon::ValidityLevel::Structural);
    REQUIRE(result.ToleranceChanges.empty());
    REQUIRE(result.Warnings.empty());
    REQUIRE(result.History.Count(HistoryRelation::Unchanged) == 0);
    REQUIRE(result.History.Count(EntityKind::Vertex, HistoryRelation::Generated) ==
            result.Shape.Count<EntityKind::Vertex>());
    REQUIRE(result.History.Count(EntityKind::Edge, HistoryRelation::Generated) == 7);
    REQUIRE(result.History.Count(EntityKind::Coedge, HistoryRelation::Generated) == 8);
    REQUIRE(result.History.Count(EntityKind::Face, HistoryRelation::Generated) == 2);
    for (const HistoryEntry& entry : result.History.GetEntries()) {
        REQUIRE(entry.Source == Stereon::NoEntity);
        REQUIRE(entry.Target != Stereon::NoEntity);
    }
}

TEST_CASE("an edit's History: deleted, modified, unchanged and generated", "[topo][operation]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    ShapeBuilder edit = Must(ShapeBuilder::Open(a));
    // Face 0 is the bottom, its loop 0 and coedges 0 to 3.
    Must(edit.RemoveFace(edit.Handle(a.Resolve(box.Faces[0]))));
    const BuilderVertexId added = Must(edit.AddVertex({9, 9, 9}, 1e-7));
    const OperationResult result = MustCommit(std::move(edit));

    REQUIRE(Fate(result, EntityKind::Face, 0).Relation == HistoryRelation::Deleted);
    REQUIRE(Fate(result, EntityKind::Face, 0).Target == Stereon::NoEntity);
    REQUIRE(Fate(result, EntityKind::Loop, 0).Relation == HistoryRelation::Deleted);
    REQUIRE(result.History.Count(EntityKind::Coedge, HistoryRelation::Deleted) == 4);
    // The other faces, and every vertex, continue unchanged at their index.
    for (std::uint32_t f = 1; f < 6; ++f) {
        const HistoryEntry fate = Fate(result, EntityKind::Face, f);
        REQUIRE(fate.Relation == HistoryRelation::Unchanged);
        REQUIRE(fate.Target == f);
    }
    REQUIRE(result.History.Count(EntityKind::Vertex, HistoryRelation::Unchanged) == 8);
    // The shell lost a face: modified. The solid still bounds the same shell.
    REQUIRE(Fate(result, EntityKind::Shell, 0).Relation == HistoryRelation::Modified);
    REQUIRE(Fate(result, EntityKind::Solid, 0).Relation == HistoryRelation::Unchanged);
    // The new vertex is the operation's creation root.
    const std::vector<HistoryEntry> made =
        result.History.OfTarget(EntityKind::Vertex, added.GetIndex());
    REQUIRE(made.size() == 1);
    REQUIRE(made[0].Relation == HistoryRelation::Generated);
    REQUIRE(made[0].Source == Stereon::NoEntity);
    // Every input entity has a disposition; the entries are in their order.
    const auto entries = result.History.GetEntries();
    REQUIRE(std::is_sorted(entries.begin(), entries.end(), [](const auto& x, const auto& y) {
        return std::pair(static_cast<int>(x.Kind), x.Source) <
               std::pair(static_cast<int>(y.Kind), y.Source);
    }));
}

TEST_CASE("a surface replaced: the face is modified, nothing else", "[topo][operation]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    ShapeBuilder edit = Must(ShapeBuilder::Open(a));
    const Stereon::SurfaceRef moved = Must(
        edit.AddSurface(Stereon::Surface(Stereon::Plane{{0, 0, 7}, Vec3{1, 0, 0}, Vec3{0, 1, 0}})));
    Must(edit.SetSurface(edit.Handle(a.Resolve(box.Faces[2])), moved));
    const OperationResult result = MustCommit(std::move(edit));
    REQUIRE(Fate(result, EntityKind::Face, 2).Relation == HistoryRelation::Modified);
    REQUIRE(result.History.Count(HistoryRelation::Modified) == 1);
    REQUIRE(result.History.Count(HistoryRelation::Deleted) == 0);
    REQUIRE(result.History.Count(HistoryRelation::Generated) == 0);
    REQUIRE(result.Warnings.empty()); // the box's Exact edges keep no cache
}

TEST_CASE("tolerance changes are reported, and the Context decides what may grow",
          "[topo][operation][tolerance]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    const auto edgeOf = [&](ShapeBuilder& edit, const std::size_t k) {
        return edit.Handle(a.Resolve(box.Edges[k]));
    };

    SECTION("a growth: reported, warned about, and the edge modified")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(a));
        Must(edit.SetTolerance(edgeOf(edit, 3), 1e-6, "closing a gap"));
        Must(edit.SetTolerance(edgeOf(edit, 3), 1e-5, "and a wider one"));
        const OperationResult result = MustCommit(std::move(edit));
        REQUIRE(result.ToleranceChanges.size() == 1); // one per entity
        const Stereon::ToleranceChange& change = result.ToleranceChanges[0];
        REQUIRE(change.Kind == EntityKind::Edge);
        REQUIRE(change.Index == 3);
        REQUIRE(change.Old == 1e-7);
        REQUIRE(change.New == 1e-5);
        REQUIRE(change.Reason == "closing a gap; and a wider one");
        REQUIRE(result.Warnings.size() == 1);
        REQUIRE(result.Warnings[0].Code == DiagnosticCode::ToleranceGrown);
        REQUIRE(result.Warnings[0].Entities.at(0).find("Edge#3") != std::string::npos);
        REQUIRE(Fate(result, EntityKind::Edge, 3).Relation == HistoryRelation::Modified);
        REQUIRE(Stereon::ToleranceOf(result.Shape, At<EntityKind::Edge>(result.Shape, 3)) == 1e-5);
    }
    SECTION("a reduction is reported without a warning; a change undone is not reported")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(a));
        Must(edit.SetTolerance(edgeOf(edit, 0), 1e-8, "tighter"));
        Must(edit.SetTolerance(edgeOf(edit, 1), 1e-3, "looser"));
        Must(edit.SetTolerance(edgeOf(edit, 1), 1e-7, "back"));
        const OperationResult result = MustCommit(std::move(edit));
        REQUIRE(result.ToleranceChanges.size() == 1);
        REQUIRE(result.ToleranceChanges[0].Index == 0);
        REQUIRE(result.Warnings.empty());
        REQUIRE(Fate(result, EntityKind::Edge, 1).Relation == HistoryRelation::Unchanged);
    }
    SECTION("growth fails when the Context says so")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(a));
        Must(edit.SetTolerance(edgeOf(edit, 3), 1e-6, "closing a gap"));
        const Context strict = Context::Default().WithToleranceGrowthPolicy(ToleranceGrowth::Fail);
        REQUIRE(CommitError(std::move(edit), strict, ErrorCode::ToleranceExceeded)
                    .find("closing a gap") != std::string::npos);
    }
    SECTION("no tolerance above the ceiling")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(a));
        Must(edit.SetTolerance(edit.Handle(a.Resolve(box.Vertices[0])), 0.5, "huge"));
        REQUIRE(CommitError(std::move(edit), Context::Default(), ErrorCode::ToleranceExceeded)
                    .find("ceiling") != std::string::npos);
    }
    SECTION("a tolerance that is not finite and non-negative is refused")
    {
        ShapeBuilder edit = Must(ShapeBuilder::Open(a));
        REQUIRE(edit.SetTolerance(edgeOf(edit, 0), -1.0, "negative").error().GetCode() ==
                ErrorCode::InvalidInput);
        REQUIRE(edit.SetTolerance(edgeOf(edit, 0), std::numeric_limits<double>::infinity(), "")
                    .error()
                    .GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("nothing is published when validation fails", "[topo][operation]")
{
    SECTION("a builder that is not structurally valid")
    {
        Body sheet = MakeSquareSheet();
        sheet.AddLoop(sheet.Faces[0]); // an empty loop: S14
        const std::string message =
            CommitError(std::move(sheet.Builder), Context::Default(), ErrorCode::InvalidInput);
        REQUIRE(message.starts_with("Commit: "));
        REQUIRE(message.find("S14 EmptyLoop") != std::string::npos);
        // The builder was consumed.
        REQUIRE_FALSE(sheet.Builder.GetToken().IsValid()); // NOLINT(bugprone-use-after-move)
    }
    SECTION("a builder that was moved from")
    {
        ShapeBuilder b;
        ShapeBuilder taken(std::move(b));
        REQUIRE(CommitError(std::move(b), Context::Default(), ErrorCode::InvalidInput)
                    .find("moved from") != std::string::npos);
    }
    SECTION("an empty Shape cannot be compacted")
    {
        REQUIRE(Stereon::CompactShape(Shape{}, Context::Default()).error().GetCode() ==
                ErrorCode::InvalidInput);
    }
}

TEST_CASE("CompactShape: every live entity continues unchanged at its new index",
          "[topo][operation][compact]")
{
    Body box = MakeBox();
    const Shape a = std::move(box).Freeze();
    ShapeBuilder edit = Must(ShapeBuilder::Open(a));
    Must(edit.RemoveFace(edit.Handle(a.Resolve(box.Faces[0]))));
    const Shape holed = MustCommit(std::move(edit)).Shape;
    const OperationResult result = Must(Stereon::CompactShape(holed, Context::Default()));
    const Stereon::Compacted reference = Must(Stereon::Compact(holed));
    REQUIRE(result.History.Count(HistoryRelation::Unchanged) == result.History.GetEntries().size());
    REQUIRE(result.History.Count(EntityKind::Face, HistoryRelation::Unchanged) == 5);
    REQUIRE(result.History.Count(EntityKind::Coedge, HistoryRelation::Unchanged) == 20);
    for (const HistoryEntry& entry : result.History.GetEntries()) {
        if (entry.Kind == EntityKind::Face) {
            REQUIRE(entry.Target == reference.Remap.NewIndex<EntityKind::Face>(entry.Source));
        }
    }
    // The tombstones left by the deletion have no entries: they were no
    // live entity of the input.
    REQUIRE(result.History.OfSource(EntityKind::Face, 0).empty());
    REQUIRE(result.Shape.Count<EntityKind::Face>() == 5);
    REQUIRE(Stereon::ValidateTopology(result.Shape).IsValid());
}
