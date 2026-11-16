// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The skew cylinder handler's choice of the cylinder whose angle
// parametrises the curve does not depend on the operands' order (ADR-0053,
// MR-02): two skew cylinders of the analytic corpus
// (analytic-ssi/cylinder-cylinder/skew-tilted-non-unit), answered in one
// order and refused in the other before ("exceeds the sample budget"), are
// answered in both, and the answer is held to both cylinders' closed-form
// equations.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

/// The distance of p from the cylinder through o about the unit axis a.
double OffCylinder(const Point3& p, const Point3& o, const Vec3& a, const double radius)
{
    const Vec3 w = p - o;
    return std::abs(Length(w - Dot(w, a) * a) - radius);
}

} // namespace

TEST_CASE("two skew cylinders are answered in either order, the parametrisation chosen canonically",
          "[intersect][cylinder][skew][metamorphic]")
{
    const Point3 oa{0.1, -0.2, 0.3};
    const Vec3 axisA{0.0, 0.0, 2.0};
    const Point3 ob{-3.1, 2.2, 1.3};
    const Vec3 axisB{3.0, 4.0, 0.0};
    const Surface a{IntersectTest::CylinderAbout(oa, axisA, 3.0)};
    const Surface b{IntersectTest::CylinderAbout(ob, axisB, 2.0)};
    const double tolerance = 1e-3;
    const Box2 all = Box2::Everything();
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "B, A" : "A, B"));
        const auto r = swapped ? Intersect(b, all, a, all, tolerance, IntersectTest::ctx)
                               : Intersect(a, all, b, all, tolerance, IntersectTest::ctx);
        INFO(IntersectTest::ErrorText(r));
        REQUIRE(r.has_value());
        REQUIRE(r->Overlaps.empty());
        // The corpus's closed form: one loop.
        REQUIRE(IntersectTest::ConnectivityOf(*r, 1e-6).Curves == 1);
        const double reach = 4.0 * IntersectTest::ctx.GetResolution();
        for (const Stereon::IntersectionBranch& branch : r->Branches) {
            REQUIRE(branch.IsCertified());
            for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                const Point3 p = Eval(branch.Curve3d, t);
                REQUIRE(OffCylinder(p, oa, Normalized(axisA), 3.0) <= reach);
                REQUIRE(OffCylinder(p, ob, Normalized(axisB), 2.0) <= reach);
            }
        }
    }
}
