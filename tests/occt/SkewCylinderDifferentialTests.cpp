// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests of the skew cylinder–cylinder intersection (the
// Boolean prototype's M3b) against OpenCASCADE's IntAna_IntQuadQuad, which
// solves the same quartic. OCCT cuts the curve at the branch points and
// hands out the branches v± separately, so curve counts differ by design;
// what must agree is the point set: every sample of OCCT's curves lies
// within a sample spacing of Stereon's branches, and the other way round.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <IntAna_Curve.hxx>
#include <IntAna_IntQuadQuad.hxx>
#include <IntAna_Quadric.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Intersect;

#include "OcctConversion.hpp"

using Stereon::Cylinder;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

using StereonTest::ToOcct;

const Stereon::Context ctx = Stereon::Context::Default();

gp_Cylinder ToOcct(const Cylinder& cylinder)
{
    const Vec3 axis = Axis(cylinder);
    return gp_Cylinder(gp_Ax3(ToOcct(cylinder.Origin), gp_Dir(axis.X, axis.Y, axis.Z),
                              gp_Dir(cylinder.XAxis.X, cylinder.XAxis.Y, cylinder.XAxis.Z)),
                       cylinder.Radius);
}

std::vector<Point3> StereonSamples(const Stereon::SurfaceSurfaceIntersection& result)
{
    std::vector<Point3> points;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        for (int k = 0; k < 720; ++k) {
            points.push_back(Eval(branch.Curve3d, 2.0 * std::numbers::pi * k / 720.0));
        }
    }
    return points;
}

std::vector<Point3> OcctSamples(const IntAna_IntQuadQuad& occt)
{
    std::vector<Point3> points;
    for (int i = 1; i <= occt.NbCurve(); ++i) {
        const IntAna_Curve& curve = occt.Curve(i);
        double theta1 = 0.0;
        double theta2 = 0.0;
        curve.Domain(theta1, theta2);
        for (int k = 0; k <= 200; ++k) {
            const gp_Pnt p = curve.Value(theta1 + (theta2 - theta1) * k / 200.0);
            points.emplace_back(p.X(), p.Y(), p.Z());
        }
    }
    return points;
}

/// The largest distance from a point of `from` to its nearest point of `to`.
double OneSidedDistance(const std::vector<Point3>& from, const std::vector<Point3>& to)
{
    double worst = 0.0;
    for (const Point3& p : from) {
        double nearest = std::numeric_limits<double>::infinity();
        for (const Point3& q : to) {
            nearest = std::min(nearest, Distance(p, q));
        }
        worst = std::max(worst, nearest);
    }
    return worst;
}

const Cylinder aboutZ{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};

Cylinder AboutX(const Point3& origin, const double radius)
{
    return Cylinder{origin, radius, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true};
}

} // namespace

TEST_CASE("skew cylinder-cylinder intersections trace the same curves as OCCT",
          "[occt][intersect][skew]")
{
    const std::vector<Cylinder> others{
        AboutX(Point3{}, 2.0),              // Through, meeting axes: two loops.
        AboutX(Point3{0.0, 4.0, 0.0}, 2.0), // Entering on one side: one loop.
        Cylinder{Point3{0.0, 0.5, 0.0}, 1.5, Normalized(Vec3{0.0, -4.0, 3.0}), Vec3{1.0, 0.0, 0.0},
                 true}, // Oblique, through: two loops.
    };
    for (const Cylinder& other : others) {
        const auto result = Intersect(Surface(aboutZ), Surface(other), ctx);
        REQUIRE(result.has_value());
        REQUIRE_FALSE(result->Branches.empty());
        const IntAna_IntQuadQuad occt(ToOcct(aboutZ), IntAna_Quadric(ToOcct(other)),
                                      Precision::Confusion());
        REQUIRE(occt.IsDone());
        REQUIRE(occt.NbCurve() > 0);
        const std::vector<Point3> ours = StereonSamples(*result);
        const std::vector<Point3> theirs = OcctSamples(occt);
        // Within a sample spacing: 720 samples round a loop of radius 3.
        REQUIRE(OneSidedDistance(theirs, ours) < 0.05);
        REQUIRE(OneSidedDistance(ours, theirs) < 0.1);
    }
}
