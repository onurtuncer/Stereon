// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Tolerance variants (the Boolean prototype's M9; ADR-0001's verification):
// union, cut and common on boxes and cylinders whose tolerances are loosened
// from 1e-4 to 1e-2 mm (decision D4), the 1e-2 variant raising the ceiling
// in its Context (decision D5). Every result is valid at the Geometric
// level; faces within the pair tolerance of each other, the larger of the
// two (decision D2), are one face, and beyond it are two; every tolerance
// that grows is logged, and refused under ToleranceGrowth::Fail or above the
// ceiling. Then the merge-policy test: vertex chains at 0.8x the pair
// tolerance cluster by ADR-0001's rules 1 to 4 exactly, byte-identically
// across thread counts and repeated runs, and a permuted input changes the
// clusters only as rule 1 predicts.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::CountOf;
using Stereon::BooleanOperation;
using Stereon::Box3;
using Stereon::Context;
using Stereon::DiagnosticCode;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::OperationResult;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ToleranceGrowth;
using Stereon::Vec3;

namespace {

/// The loosened tolerances of ADR-0001's verification.
constexpr double Loosened[] = {1e-4, 1e-3, 1e-2};

/// The Context a variant runs in: the default, its ceiling raised above the
/// loosened tolerance and any growth from it (decision D5).
Context ContextFor(const double tolerance)
{
    return Context::Default().WithToleranceCeiling(
        std::max(Context::DefaultToleranceCeiling, 4.0 * tolerance));
}

/// A box whose vertices and edges have the given tolerance.
Shape LooseBox(const Point3& lo, const Point3& hi, const double tolerance)
{
    const auto box =
        MakeBoxSolid(Box3::FromCorners(lo, hi), ContextFor(tolerance).WithResolution(tolerance));
    REQUIRE(box.has_value());
    return box->Shape;
}

/// A cylinder about z whose vertices and edges have the given tolerance.
Shape LooseDrum(const Point3& base, const double radius, const double height,
                const double tolerance)
{
    const auto drum = MakeCylinderSolid(base, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, radius,
                                        height, ContextFor(tolerance).WithResolution(tolerance));
    REQUIRE(drum.has_value());
    return drum->Shape;
}

OperationResult Run(const Shape& a, const Shape& b, const BooleanOperation operation,
                    const Context& context)
{
    auto result = Boolean(a, b, operation, context);
    const std::string error = result ? std::string() : std::string(result.error().GetMessage());
    INFO(std::string(Stereon::ToString(operation)) << ": " << error);
    REQUIRE(result.has_value());
    REQUIRE(Stereon::ValidateTopology(result->Shape).IsValid());
    REQUIRE(result->Checked == Stereon::ValidityLevel::Geometric);
    return std::move(*result);
}

/// Every tolerance that grew is logged, and warned of with its reason.
void RequireGrowthLogged(const OperationResult& result)
{
    std::size_t grown = 0;
    for (const Stereon::ToleranceChange& change : result.ToleranceChanges) {
        REQUIRE_FALSE(change.Reason.empty());
        grown += change.New > change.Old ? 1 : 0;
    }
    std::size_t warned = 0;
    for (const Stereon::Diagnostic& warning : result.Warnings) {
        warned += warning.Code == DiagnosticCode::ToleranceGrown ? 1 : 0;
    }
    REQUIRE(warned == grown);
}

/// The largest vertex and edge tolerance of a shape.
double LargestTolerance(const Shape& shape)
{
    double largest = 0.0;
    for (const Stereon::VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        largest = std::max(largest, ToleranceOf(shape, v));
    }
    for (const Stereon::EdgeId e : *AllEntities<EntityKind::Edge>(shape)) {
        largest = std::max(largest, ToleranceOf(shape, e));
    }
    return largest;
}

} // namespace

TEST_CASE("loosened boxes and cylinders give the results the exact ones do", "[boolean][tolerance]")
{
    for (const double t : Loosened) {
        CAPTURE(t);
        const Context context = ContextFor(t);
        const Shape a = LooseBox(Point3{}, Point3{10.0, 10.0, 10.0}, t);
        const Shape corner = LooseBox(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}, t);
        REQUIRE(CountOf<EntityKind::Face>(Run(a, corner, BooleanOperation::Union, context).Shape) ==
                12);
        REQUIRE(CountOf<EntityKind::Face>(
                    Run(a, corner, BooleanOperation::Common, context).Shape) == 6);
        REQUIRE(CountOf<EntityKind::Face>(Run(a, corner, BooleanOperation::Cut, context).Shape) ==
                9);

        const Shape drum = LooseDrum(Point3{5.0, 5.0, -2.0}, 3.0, 15.0, t);
        for (const auto& [operation, faces] :
             {std::pair{BooleanOperation::Union, 10}, std::pair{BooleanOperation::Common, 3},
              std::pair{BooleanOperation::Cut, 7}}) {
            const OperationResult result = Run(a, drum, operation, context);
            REQUIRE(CountOf<EntityKind::Face>(result.Shape) == static_cast<std::size_t>(faces));
            RequireGrowthLogged(result);
            REQUIRE(LargestTolerance(result.Shape) <= context.GetToleranceCeiling());
        }
    }
}

TEST_CASE("faces within the pair tolerance are one face, beyond it two", "[boolean][tolerance]")
{
    for (const double t : Loosened) {
        CAPTURE(t);
        const Context context = ContextFor(t);
        // A at the resolution, B loosened to t: the pair tolerance is t.
        const Shape a = LooseBox(Point3{}, Point3{10.0, 10.0, 10.0}, context.GetResolution());
        const Shape near = LooseBox(Point3{10.0 + 0.8 * t, 0.0, 0.0}, Point3{20.0, 10.0, 10.0}, t);
        const OperationResult united = Run(a, near, BooleanOperation::Union, context);
        // As the exact face-to-face union: the shared face gone from both,
        // B's as well as A's though only B's is loose.
        REQUIRE(CountOf<EntityKind::Face>(united.Shape) == 10);
        REQUIRE(CountOf<EntityKind::Solid>(united.Shape) == 1);
        // A's vertices identified with B's, 0.8 t away: grown to B's
        // tolerance at least, and logged.
        REQUIRE_FALSE(united.ToleranceChanges.empty());
        RequireGrowthLogged(united);
        REQUIRE(LargestTolerance(united.Shape) >= t);

        const Shape far = LooseBox(Point3{10.0 + 1.2 * t, 0.0, 0.0}, Point3{20.0, 10.0, 10.0}, t);
        const OperationResult apart = Run(a, far, BooleanOperation::Union, context);
        REQUIRE(CountOf<EntityKind::Face>(apart.Shape) == 12);
        REQUIRE(CountOf<EntityKind::Solid>(apart.Shape) == 2);
        REQUIRE(apart.ToleranceChanges.empty());
    }
}

TEST_CASE("a tolerance that would grow is refused under the Fail policy or above the ceiling",
          "[boolean][tolerance]")
{
    const double t = 1e-3;
    const Context context = ContextFor(t);
    const Shape a = LooseBox(Point3{}, Point3{10.0, 10.0, 10.0}, context.GetResolution());
    const Shape near = LooseBox(Point3{10.0 + 0.8 * t, 0.0, 0.0}, Point3{20.0, 10.0, 10.0}, t);

    const auto refused = Boolean(a, near, BooleanOperation::Union,
                                 context.WithToleranceGrowthPolicy(ToleranceGrowth::Fail));
    REQUIRE_FALSE(refused.has_value());
    REQUIRE(refused.error().GetCode() == ErrorCode::ToleranceExceeded);

    // Inputs loosened to 1e-2 under a ceiling of 5e-3: above it before any
    // growth (decision D5: the 1e-2 variant raises the ceiling).
    const Shape loose = LooseBox(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}, 1e-2);
    const auto above =
        Boolean(a, loose, BooleanOperation::Union, Context::Default().WithToleranceCeiling(5e-3));
    REQUIRE_FALSE(above.has_value());
    REQUIRE(above.error().GetCode() == ErrorCode::ToleranceExceeded);
}

// --- The merge-policy chain (ADR-0001, rules 1 to 4) -------------------------------

namespace {

/// The merge of vertices placed at `xs` along x, each with tolerance t, in
/// the order given (rule 1: candidates by index, which is creation order):
/// each vertex's representative's index, and the bits of the kept vertices'
/// positions and tolerances.
struct Merged {
    std::vector<std::size_t> Representative;
    std::vector<std::uint64_t> Bits;
};

Merged MergeChain(std::span<const double> xs, const double t, const Context& context)
{
    Stereon::ShapeBuilder builder;
    std::vector<Stereon::BuilderVertexId> vertices;
    for (const double x : xs) {
        const auto v = builder.AddVertex(Point3{x, 0.0, 0.0}, t);
        REQUIRE(v.has_value());
        vertices.push_back(*v);
    }
    const auto merged = builder.MergeVertices(vertices, context);
    const std::string error = merged ? std::string() : std::string(merged.error().GetMessage());
    INFO(error);
    REQUIRE(merged.has_value());
    Merged result;
    for (const Stereon::BuilderVertexId v : vertices) {
        std::size_t representative = vertices.size();
        for (const Stereon::VertexMergeEntry& entry : merged->Entries) {
            if (entry.Source == v) {
                representative = static_cast<std::size_t>(
                    std::ranges::find(vertices, entry.Representative) - vertices.begin());
            }
        }
        result.Representative.push_back(
            representative == vertices.size()
                ? static_cast<std::size_t>(std::ranges::find(vertices, v) - vertices.begin())
                : representative);
    }
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (result.Representative[i] == i) {
            const Point3 p = builder.PositionOf(vertices[i]);
            result.Bits.push_back(std::bit_cast<std::uint64_t>(p.X));
        }
    }
    return result;
}

/// The thread counts of ADR-0007's determinism checks, and 0 (the hardware's).
constexpr unsigned ThreadCounts[] = {1U, 2U, 8U, 32U, 0U};

} // namespace

TEST_CASE("vertex chains at 0.8x the pair tolerance cluster by the representative rule",
          "[boolean][tolerance][merge][determinism]")
{
    for (const double t : Loosened) {
        CAPTURE(t);
        const Context context = ContextFor(t);
        const double step = 0.8 * t;
        // A - B - C - D collinear: B joins A; C, 1.6 t from A, opens its own
        // cluster though it is 0.8 t from B (rule 3: within the tolerance of
        // the representative, not of any member); D joins C.
        const std::vector<double> chain{0.0, step, 2.0 * step, 3.0 * step};
        const Merged forward = MergeChain(chain, t, context);
        REQUIRE(forward.Representative == std::vector<std::size_t>{0, 0, 2, 2});
        // Rule 2: the representatives keep their positions, never averaged.
        REQUIRE(forward.Bits ==
                std::vector<std::uint64_t>{std::bit_cast<std::uint64_t>(0.0),
                                           std::bit_cast<std::uint64_t>(2.0 * step)});

        // Byte-identical across thread counts and repeated runs.
        for (const unsigned threads : ThreadCounts) {
            CAPTURE(threads);
            const Context parallel =
                context.WithScheduler(Stereon::Scheduler::Parallel).WithThreadCount(threads);
            for (int run = 0; run < 2; ++run) {
                const Merged again = MergeChain(chain, t, parallel);
                REQUIRE(again.Representative == forward.Representative);
                REQUIRE(again.Bits == forward.Bits);
            }
        }

        // Three in a row, created the other way round: rule 1 orders the
        // candidates by index, so the cluster forms about the other end.
        const std::vector<double> three{0.0, step, 2.0 * step};
        REQUIRE(MergeChain(three, t, context).Representative == std::vector<std::size_t>{0, 0, 2});
        const std::vector<double> reversed{2.0 * step, step, 0.0};
        const Merged back = MergeChain(reversed, t, context);
        REQUIRE(back.Representative == std::vector<std::size_t>{0, 0, 2});
        REQUIRE(back.Bits == std::vector<std::uint64_t>{std::bit_cast<std::uint64_t>(2.0 * step),
                                                        std::bit_cast<std::uint64_t>(0.0)});
    }
}

TEST_CASE("a loosened Boolean is byte-identical across thread counts and repeated runs",
          "[boolean][tolerance][determinism]")
{
    const double t = 1e-3;
    const Context context = ContextFor(t);
    const Shape a = LooseBox(Point3{}, Point3{10.0, 10.0, 10.0}, context.GetResolution());
    const Shape near = LooseBox(Point3{10.0 + 0.8 * t, 0.0, 0.0}, Point3{20.0, 10.0, 10.0}, t);
    const auto fingerprint = [&](const Context& c) {
        const OperationResult result = Run(a, near, BooleanOperation::Union, c);
        std::vector<std::uint64_t> bits;
        for (const Stereon::VertexId v : *AllEntities<EntityKind::Vertex>(result.Shape)) {
            const Point3 p = PositionOf(result.Shape, v);
            for (const double x : {p.X, p.Y, p.Z, ToleranceOf(result.Shape, v)}) {
                bits.push_back(std::bit_cast<std::uint64_t>(x));
            }
        }
        for (const Stereon::EdgeId e : *AllEntities<EntityKind::Edge>(result.Shape)) {
            const Stereon::Interval range = RangeOf(result.Shape, e);
            for (const double x :
                 {range.GetLower(), range.GetUpper(), ToleranceOf(result.Shape, e)}) {
                bits.push_back(std::bit_cast<std::uint64_t>(x));
            }
            bits.push_back(StartOf(result.Shape, e).GetIndex());
            bits.push_back(EndOf(result.Shape, e).GetIndex());
        }
        for (const Stereon::ToleranceChange& change : result.ToleranceChanges) {
            bits.push_back(change.Index);
            bits.push_back(std::bit_cast<std::uint64_t>(change.New));
        }
        return bits;
    };
    const std::vector<std::uint64_t> reference = fingerprint(context);
    for (const unsigned threads : ThreadCounts) {
        CAPTURE(threads);
        const Context parallel =
            context.WithScheduler(Stereon::Scheduler::Parallel).WithThreadCount(threads);
        for (int run = 0; run < 2; ++run) {
            REQUIRE(fingerprint(parallel) == reference);
        }
    }
}
