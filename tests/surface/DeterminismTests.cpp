// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The determinism guarantee of ADR-0007 (G1) on the operations that run in
// parallel: batch projection onto B-spline curves and surfaces gives the
// same bits on the serial scheduler and at 1, 2, 8 and 32 threads, and the
// same bits as projecting one point at a time. The nightly workflow repeats
// the comparison across processes on more points (tests/determinism/).

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

import Stereon.Surface;

#include "DeterminismCases.hpp"

using Stereon::Context;
using Stereon::CurveProjector;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Scheduler;
using Stereon::SurfaceProjector;

namespace {

constexpr std::size_t PointCount = 1500;

Context Serial()
{
    return Context::Default().WithScheduler(Scheduler::Serial);
}

Context Threads(const unsigned count)
{
    return Context::Default().WithScheduler(Scheduler::Parallel).WithThreadCount(count);
}

} // namespace

TEST_CASE("batch projection gives the same bits for every scheduler and thread count",
          "[surface][determinism]")
{
    const std::vector<Point3> points = DeterminismTest::Points(PointCount);
    const CurveProjector curve(DeterminismTest::Curve());
    const SurfaceProjector surface(DeterminismTest::Surface());

    // The serial scheduler is the reference, and it answers as Project does
    // one point at a time.
    const std::vector<std::uint64_t> reference =
        DeterminismTest::Run(curve, surface, points, Serial());
    std::vector<double> oneByOne;
    std::vector<Point2> oneByOneUv;
    for (const Point3& p : points) {
        oneByOne.push_back(curve.Project(p));
        oneByOneUv.push_back(surface.Project(p));
    }
    REQUIRE(DeterminismTest::Bits(oneByOne, oneByOneUv) == reference);

    // Every thread count, twice, and the hardware's count (0).
    for (const unsigned count : {1U, 2U, 8U, 32U, 0U}) {
        INFO(count << " threads");
        REQUIRE(DeterminismTest::Run(curve, surface, points, Threads(count)) == reference);
        REQUIRE(DeterminismTest::Run(curve, surface, points, Threads(count)) == reference);
    }
}

TEST_CASE("batch projection of no points and of fewer points than a chunk",
          "[surface][determinism]")
{
    const CurveProjector curve(DeterminismTest::Curve());
    const SurfaceProjector surface(DeterminismTest::Surface());
    const std::vector<Point3> none;
    REQUIRE(ProjectAll(curve, none, Threads(8)).empty());
    REQUIRE(ProjectAll(surface, none, Serial()).empty());
    const std::vector<Point3> few = DeterminismTest::Points(5);
    REQUIRE(DeterminismTest::Run(curve, surface, few, Threads(32)) ==
            DeterminismTest::Run(curve, surface, few, Serial()));
}
