// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the winding number and the seeds stand on, on sheets built by hand:
// a polyline pcurve is sampled at its corners, a sample count the cap
// limits reports the deviation it reaches rather than the target, and a
// region far narrower than a face's extent still gets a seed inside it.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::ctx;
using Catch::Matchers::WithinRel;
using Stereon::BuilderEdgeId;
using Stereon::BuilderFaceId;
using Stereon::BuilderLoopId;
using Stereon::BuilderVertexId;
using Stereon::Curve2d;
using Stereon::EdgeCurveKind;
using Stereon::EdgeUseSense;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ShapeBuilder;

namespace {

template <typename T>
T Must(Stereon::Result<T> result)
{
    if (!result) {
        FAIL("unexpected error: " << std::string(result.error().GetMessage()));
    }
    return std::move(*result);
}

void Must(const Stereon::Result<void>& result)
{
    if (!result) {
        FAIL("unexpected error: " << std::string(result.error().GetMessage()));
    }
}

constexpr double Tight = 1e-5;

/// A sheet face in the plane z = 0 with the identity frame, in a shell of
/// no solid, built one loop at a time.
struct Sheet {
    ShapeBuilder Builder;
    BuilderFaceId Face;

    Sheet()
    {
        const auto plane = Stereon::Surface(
            Stereon::Plane{{0, 0, 0}, Stereon::Vec3{1, 0, 0}, Stereon::Vec3{0, 1, 0}});
        Face = Must(Builder.AddFace(Must(Builder.AddSurface(plane)), false));
        Must(Builder.AddFaceToShell(Must(Builder.AddShell()), Face));
    }

    BuilderVertexId Vertex(const Point2& p)
    {
        return Must(Builder.AddVertex({p.X, p.Y, 0}, Tight));
    }

    /// A Boundary edge from `a` to `b` over `range` along `pcurve`, used once
    /// by `loop` along its parameter.
    BuilderEdgeId Edge(const BuilderLoopId loop, const BuilderVertexId a, const BuilderVertexId b,
                       const Interval& range, const Curve2d& pcurve)
    {
        const BuilderEdgeId e = Must(Builder.AddEdge({.Kind = EdgeCurveKind::Boundary,
                                                      .Start = a,
                                                      .End = b,
                                                      .Range = range,
                                                      .Tolerance = Tight}));
        Must(Builder.AppendCoedge(loop, e, EdgeUseSense::Same, Must(Builder.AddPCurve(pcurve))));
        return e;
    }

    /// A closed polygon of line edges as a new loop.
    void Polygon(std::span<const Point2> corners)
    {
        const BuilderLoopId loop = Must(Builder.AddLoop(Face));
        std::vector<BuilderVertexId> vertices;
        for (const Point2& p : corners) {
            vertices.push_back(Vertex(p));
        }
        for (std::size_t i = 0; i < corners.size(); ++i) {
            const std::size_t j = (i + 1) % corners.size();
            Edge(loop, vertices[i], vertices[j], Interval(0.0, 1.0),
                 Curve2d(Stereon::Line2d{corners[i], corners[j] - corners[i]}));
        }
    }

    Shape Freeze() && { return Must(Stereon::Freeze(std::move(Builder))); }
};

/// The shape's edge of the given index, in creation order.
Stereon::EdgeId EdgeAt(const Shape& shape, const std::size_t index)
{
    const auto edges = AllEntities<Stereon::EntityKind::Edge>(shape);
    REQUIRE(edges.has_value());
    REQUIRE(edges->size() > index);
    return *std::next(edges->begin(), static_cast<std::ptrdiff_t>(index));
}

} // namespace

TEST_CASE("a polyline pcurve is sampled at each of its corners", "[boolean][winding]")
{
    // A triangle whose apex is the corner of a degree-one B-spline: one chord
    // from (0, 0) to (2, 0) would miss the apex (1, 1) by 1 however fine the
    // target, since the polyline curves nowhere between its knots.
    Sheet sheet;
    const BuilderLoopId loop = Must(sheet.Builder.AddLoop(sheet.Face));
    const BuilderVertexId a = sheet.Vertex({0, 0});
    const BuilderVertexId b = sheet.Vertex({2, 0});
    sheet.Edge(loop, a, b, Interval(0.0, 2.0),
               Curve2d(Stereon::BSplineCurve2d{
                   1, {0.0, 0.0, 1.0, 2.0, 2.0}, {Point2{0, 0}, Point2{1, 1}, Point2{2, 0}}, {}}));
    sheet.Edge(loop, b, a, Interval(0.0, 1.0), Curve2d(Stereon::Line2d{{2, 0}, {-2, 0}}));
    const Shape shape = std::move(sheet).Freeze();

    const auto sampling = Stereon::Detail::WindingSamplingOf(shape, EdgeAt(shape, 0), 1e-5);
    REQUIRE(sampling.has_value());
    REQUIRE(sampling->Parameters == std::vector<double>{0.0, 1.0, 2.0});
    REQUIRE(sampling->Deviation == 0.0);
}

TEST_CASE("a capped sampling reports the deviation it reaches, not the target",
          "[boolean][winding]")
{
    // The unit circle at a target of 1e-14 asks for some 4e7 chords; the cap
    // of 2^20 leaves each about 6e-6 long, within h^2 rate / 8 of the circle.
    Sheet sheet;
    const BuilderLoopId loop = Must(sheet.Builder.AddLoop(sheet.Face));
    const BuilderVertexId v = sheet.Vertex({1, 0});
    sheet.Edge(loop, v, v, Interval(0.0, 2.0 * std::numbers::pi),
               Curve2d(Stereon::Circle2d{{0, 0}, 1.0, {1, 0}, {0, 1}}));
    const Shape shape = std::move(sheet).Freeze();

    constexpr double target = 1e-14;
    const auto sampling = Stereon::Detail::WindingSamplingOf(shape, EdgeAt(shape, 0), target);
    REQUIRE(sampling.has_value());
    const auto segments = static_cast<double>(sampling->Parameters.size() - 1);
    REQUIRE(segments == 1 << 20);
    REQUIRE(sampling->Deviation > 100.0 * target);
    // The plane's |Su| + |Sv| = 2 times the circle's |p''| bound 2.
    const double h = 2.0 * std::numbers::pi / segments;
    REQUIRE_THAT(sampling->Deviation, WithinRel(h * h * 4.0 / 8.0, 1e-12));
    // Uncapped, the target holds.
    const auto easy = Stereon::Detail::WindingSamplingOf(shape, EdgeAt(shape, 0), 1e-3);
    REQUIRE(easy.has_value());
    REQUIRE(easy->Deviation <= 1e-3);
}

TEST_CASE("a seed's ray takes no sign the certified filter does not settle", "[boolean][classify]")
{
    // A cylinder of radius 0.1 about z, and rays along x at heights about
    // 0.1. At exactly 0.1 the line is tangent in the input doubles, but
    // |w|^2 - r^2 rounds, so a discriminant formed in doubles may come out
    // negative and be read as a certain miss; the enclosure straddles zero
    // and the ray is ambiguous instead, so another direction is tried.
    const Stereon::Surface drill{
        Stereon::Cylinder{Point3{0, 0, 0}, 0.1, Stereon::Vec3{1, 0, 0}, Stereon::Vec3{0, 1, 0}}};
    const Stereon::Vec3 along{1, 0, 0};
    REQUIRE_FALSE(Stereon::Detail::RayMeets(drill, Point3{-5, 0.1, 0}, along).has_value());
    // Clear of the tangent the filter decides: two crossings, or none.
    const auto through = Stereon::Detail::RayMeets(drill, Point3{-5, 0.05, 0}, along);
    REQUIRE(through.has_value());
    REQUIRE(through->size() == 2);
    const auto past = Stereon::Detail::RayMeets(drill, Point3{-5, 0.2, 0}, along);
    REQUIRE(past.has_value());
    REQUIRE(past->empty());
    // Along the axis, and parallel to a plane: ambiguous.
    REQUIRE_FALSE(
        Stereon::Detail::RayMeets(drill, Point3{0.5, 0, 0}, Stereon::Vec3{0, 0, 1}).has_value());
    const Stereon::Surface floor{
        Stereon::Plane{Point3{0, 0, 0}, Stereon::Vec3{1, 0, 0}, Stereon::Vec3{0, 1, 0}}};
    REQUIRE_FALSE(Stereon::Detail::RayMeets(floor, Point3{0, 0, 1}, along).has_value());
    const auto down = Stereon::Detail::RayMeets(floor, Point3{0, 0, 1}, Stereon::Vec3{0, 0, -1});
    REQUIRE(down.has_value());
    REQUIRE(down->size() == 1);
}

TEST_CASE("a narrow region gets a seed strictly inside it", "[boolean][classify]")
{
    // A square annulus between +-1 and +-0.99: a hundredth wide, against an
    // extent of 2. Offsets fixed at a twentieth and more all land in the hole.
    Sheet sheet;
    const std::vector<Point2> outer{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const std::vector<Point2> inner{{-0.99, -0.99}, {-0.99, 0.99}, {0.99, 0.99}, {0.99, -0.99}};
    sheet.Polygon(outer);
    sheet.Polygon(inner);
    const Shape shape = std::move(sheet).Freeze();
    const auto faces = AllEntities<Stereon::EntityKind::Face>(shape);
    REQUIRE(faces.has_value());
    const Stereon::FaceId face = *faces->begin();

    const auto seed = Stereon::Detail::SeedPointIn(shape, face, ctx);
    REQUIRE(seed.has_value());
    REQUIRE(seed->has_value());
    const Point2 p = **seed;
    INFO(p.X << ", " << p.Y);
    REQUIRE(std::max(std::abs(p.X), std::abs(p.Y)) < 1.0);
    REQUIRE(std::max(std::abs(p.X), std::abs(p.Y)) > 0.99);
    const auto inside = FaceContains(shape, face, p, ctx);
    REQUIRE(inside.has_value());
    REQUIRE(inside->has_value());
    REQUIRE(**inside);
}
