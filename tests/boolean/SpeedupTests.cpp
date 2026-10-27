// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The many-face fixture of the Boolean prototype's M10 (ADR-0007): a slab
// with a grid of cylindrical holes, cut by one operand that is every
// cylinder of the grid. The result must not depend on the thread count, at
// 1, 2, 8 and 32 threads and over repeated runs; and, hidden from the
// default run since it measures the machine, the cut is timed serially and
// on the pool, stage by stage.

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::CountOf;
using Stereon::BooleanOperation;
using Stereon::Box3;
using Stereon::Context;
using Stereon::EntityKind;
using Stereon::OperationResult;
using Stereon::Point3;
using Stereon::Scheduler;
using Stereon::Shape;
using Stereon::Vec3;

namespace {

const Context serial = Context::Default().WithScheduler(Scheduler::Serial);

/// The slab [0, 10n] x [0, 10n] x [0, 10].
Shape Slab(const int n)
{
    const double side = 10.0 * n;
    const auto slab = MakeBoxSolid(Box3::FromCorners(Point3{}, Point3{side, side, 10.0}), serial);
    REQUIRE(slab.has_value());
    return slab->Shape;
}

/// n x n cylinders of radius 2 at the centres of the slab's 10 x 10 cells,
/// from z = -5 to 15, as one operand: a union of disjoint solids, built up
/// a row at a time.
Shape Drills(const int n)
{
    const auto drill = [](const int i, const int j) {
        const auto made =
            MakeCylinderSolid(Point3{10.0 * i + 5.0, 10.0 * j + 5.0, -5.0}, Vec3{1.0, 0.0, 0.0},
                              Vec3{0.0, 1.0, 0.0}, 2.0, 20.0, serial);
        REQUIRE(made.has_value());
        return made->Shape;
    };
    Shape all = drill(0, 0);
    for (int k = 1; k < n * n; ++k) {
        auto united = Boolean(all, drill(k % n, k / n), BooleanOperation::Union, serial);
        const std::string error = united ? std::string() : std::string(united.error().GetMessage());
        INFO(error);
        REQUIRE(united.has_value());
        all = std::move(united->Shape);
    }
    return all;
}

/// The bits of a result: every vertex's position and tolerance and every
/// edge's range, ends and tolerance, in index order.
std::vector<std::uint64_t> Fingerprint(const Shape& shape)
{
    std::vector<std::uint64_t> bits;
    for (const Stereon::VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        const Point3 p = PositionOf(shape, v);
        for (const double x : {p.X, p.Y, p.Z, ToleranceOf(shape, v)}) {
            bits.push_back(std::bit_cast<std::uint64_t>(x));
        }
    }
    for (const Stereon::EdgeId e : *AllEntities<EntityKind::Edge>(shape)) {
        const Stereon::Interval range = RangeOf(shape, e);
        for (const double x : {range.GetLower(), range.GetUpper(), ToleranceOf(shape, e)}) {
            bits.push_back(std::bit_cast<std::uint64_t>(x));
        }
        bits.push_back(StartOf(shape, e).GetIndex());
        bits.push_back(EndOf(shape, e).GetIndex());
    }
    return bits;
}

OperationResult Cut(const Shape& slab, const Shape& drills, const Context& context)
{
    auto result = Boolean(slab, drills, BooleanOperation::Cut, context);
    const std::string error = result ? std::string() : std::string(result.error().GetMessage());
    INFO(error);
    REQUIRE(result.has_value());
    return std::move(*result);
}

double Seconds(const std::chrono::steady_clock::duration d)
{
    return std::chrono::duration<double>(d).count();
}

} // namespace

TEST_CASE("a slab drilled with a grid of holes is the same at every thread count",
          "[boolean][determinism]")
{
    const int n = 3;
    const Shape slab = Slab(n);
    const Shape drills = Drills(n);
    REQUIRE(CountOf<EntityKind::Solid>(drills) == static_cast<std::size_t>(n * n));
    const OperationResult reference = Cut(slab, drills, serial);
    // The slab's six faces, and one hole's face per drill.
    REQUIRE(CountOf<EntityKind::Face>(reference.Shape) == static_cast<std::size_t>(6 + n * n));
    REQUIRE(CountOf<EntityKind::Solid>(reference.Shape) == 1);
    const std::vector<std::uint64_t> bits = Fingerprint(reference.Shape);
    for (const unsigned threads : {1U, 2U, 8U, 32U, 0U}) {
        CAPTURE(threads);
        const Context parallel =
            Context::Default().WithScheduler(Scheduler::Parallel).WithThreadCount(threads);
        for (int run = 0; run < 2; ++run) {
            const OperationResult again = Cut(slab, drills, parallel);
            REQUIRE(Fingerprint(again.Shape) == bits);
            REQUIRE(again.ToleranceChanges.size() == reference.ToleranceChanges.size());
        }
    }
}

TEST_CASE("the drilled slab's speed-up, stage by stage", "[.][boolean][speedup]")
{
    for (const int n : {7, 10}) {
        using Clock = std::chrono::steady_clock;
        const auto built = Clock::now();
        const Shape slab = Slab(n);
        const Shape drills = Drills(n);
        std::printf("holes %4d: the drills built in %7.3f s\n", n * n,
                    Seconds(Clock::now() - built));
        for (const unsigned threads : {1U, 8U}) {
            const Context context = threads == 1U ? serial
                                                  : Context::Default()
                                                        .WithScheduler(Scheduler::Parallel)
                                                        .WithThreadCount(threads);
            const auto t0 = Clock::now();
            const auto imprinted = ImprintBodies(slab, drills, context);
            REQUIRE(imprinted.has_value());
            const auto t1 = Clock::now();
            const auto classified = Classify(*imprinted, context);
            REQUIRE(classified.has_value());
            const auto t2 = Clock::now();
            const OperationResult cut = Cut(slab, drills, context);
            const auto t3 = Clock::now();
            std::printf("holes %4d threads %2u: imprint %7.3f s, classify %7.3f s, whole cut "
                        "%7.3f s (%zu faces)\n",
                        n * n, threads, Seconds(t1 - t0), Seconds(t2 - t1), Seconds(t3 - t2),
                        CountOf<EntityKind::Face>(cut.Shape));
        }
    }
}
