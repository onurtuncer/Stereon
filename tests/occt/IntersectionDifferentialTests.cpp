// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Differential tests of the analytic surface intersections (the Boolean
// prototype's M3) against OpenCASCADE's IntAna_QuadQuadGeo, the same closed
// forms by another hand: for each plane–plane, plane–cylinder and parallel
// cylinder–cylinder configuration the two kernels must agree on the kind of
// result and, where there are curves, every point of Stereon's branches
// must lie on one of OCCT's.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Ellipse.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_ResultType.hxx>
#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

import Stereon.Intersect;

#include "OcctConversion.hpp"

using Stereon::Cylinder;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

const Stereon::Context ctx = Stereon::Context::Default();

gp_Dir ToDir(const Vec3& v)
{
    return gp_Dir(v.X, v.Y, v.Z);
}

gp_Pln ToOcct(const Plane& plane)
{
    return gp_Pln(
        gp_Ax3(ToOcct(plane.Origin), ToDir(Cross(plane.XAxis, plane.YAxis)), ToDir(plane.XAxis)));
}

gp_Cylinder ToOcct(const Cylinder& cylinder)
{
    return gp_Cylinder(
        gp_Ax3(ToOcct(cylinder.Origin), ToDir(Axis(cylinder)), ToDir(cylinder.XAxis)),
        cylinder.Radius);
}

/// The distance from p to the nearest of OCCT's result curves.
double DistanceToResult(const IntAna_QuadQuadGeo& occt, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (int i = 1; i <= occt.NbSolutions(); ++i) {
        switch (occt.TypeInter()) {
        case IntAna_Line:
            nearest = std::min(nearest, occt.Line(i).Distance(ToOcct(p)));
            break;
        case IntAna_Circle:
            nearest = std::min(nearest, occt.Circle(i).Distance(ToOcct(p)));
            break;
        case IntAna_Ellipse: {
            GeomAPI_ProjectPointOnCurve projection(ToOcct(p), new Geom_Ellipse(occt.Ellipse(i)));
            nearest = std::min(nearest, projection.LowerDistance());
            break;
        }
        default:
            break;
        }
    }
    return nearest;
}

/// Every sample of every branch lies on one of OCCT's curves.
void RequireOnOcct(const Stereon::SurfaceSurfaceIntersection& result,
                   const IntAna_QuadQuadGeo& occt, const double tolerance)
{
    REQUIRE(occt.IsDone());
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const double lower = std::isinf(branch.Range.GetLower()) ? -7.0 : branch.Range.GetLower();
        const double upper = std::isinf(branch.Range.GetUpper()) ? 7.0 : branch.Range.GetUpper();
        for (int k = 0; k <= 40; ++k) {
            const double t = lower + (upper - lower) * k / 40.0;
            REQUIRE(DistanceToResult(occt, Eval(branch.Curve3d, t)) <= tolerance);
        }
    }
}

const Plane xy{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
const Cylinder aboutZ{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
const Cylinder aboutX{Point3{}, 5.0, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}, true};

Plane Tilted(const Point3& origin, const Vec3& normal)
{
    const Vec3 n = Normalized(normal);
    const Vec3 seed = std::abs(n.X) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
    const Vec3 x = Normalized(Cross(seed, n));
    return Plane{origin, x, Cross(n, x)};
}

} // namespace

TEST_CASE("plane-plane intersections agree with OCCT", "[occt][intersect]")
{
    const Plane tilted = Tilted(Point3{1.0, 2.0, 3.0}, Vec3{1.0, 1.0, 2.0});
    const auto line = Intersect(Surface(xy), Surface(tilted), ctx);
    REQUIRE(line->Branches.size() == 1);
    const IntAna_QuadQuadGeo occtLine(ToOcct(xy), ToOcct(tilted), Precision::Angular(),
                                      Precision::Confusion());
    REQUIRE(occtLine.TypeInter() == IntAna_Line);
    RequireOnOcct(*line, occtLine, 1e-9);

    const Plane lifted{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    REQUIRE(Intersect(Surface(xy), Surface(lifted), ctx)->IsEmpty());
    REQUIRE(
        IntAna_QuadQuadGeo(ToOcct(xy), ToOcct(lifted), Precision::Angular(), Precision::Confusion())
            .TypeInter() == IntAna_Empty);

    const Plane same{Point3{3.0, 4.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{-1.0, 0.0, 0.0}};
    REQUIRE(Intersect(Surface(xy), Surface(same), ctx)->Overlaps.size() == 1);
    REQUIRE(
        IntAna_QuadQuadGeo(ToOcct(xy), ToOcct(same), Precision::Angular(), Precision::Confusion())
            .TypeInter() == IntAna_Same);
}

TEST_CASE("plane-cylinder intersections agree with OCCT", "[occt][intersect]")
{
    struct Case {
        Plane Cut;
        Cylinder Body;
        IntAna_ResultType Type;
        std::size_t Branches;
    };
    const std::vector<Case> cases{
        {Plane{Point3{0.0, 0.0, 3.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutZ,
         IntAna_Circle, 1},
        {Tilted(Point3{0.0, 1.0, 0.0}, Vec3{0.0, 1.0, 1.0}), aboutZ, IntAna_Ellipse, 1},
        {Tilted(Point3{2.0, 0.0, 1.0}, Vec3{1.0, 0.0, 3.0}), aboutZ, IntAna_Ellipse, 1},
        {Plane{Point3{0.0, 0.0, 2.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Line, 2},
        {Plane{Point3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Line, 1},
        {Plane{Point3{0.0, 0.0, 6.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}, aboutX,
         IntAna_Empty, 0},
    };
    for (const Case& c : cases) {
        const auto result = Intersect(Surface(c.Cut), Surface(c.Body), ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.size() == c.Branches);
        const IntAna_QuadQuadGeo occt(ToOcct(c.Cut), ToOcct(c.Body), Precision::Angular(),
                                      Precision::Confusion());
        REQUIRE(occt.TypeInter() == c.Type);
        REQUIRE(static_cast<std::size_t>(std::max(occt.NbSolutions(), 0)) == c.Branches);
        RequireOnOcct(*result, occt, 1e-9);
    }
}

TEST_CASE("parallel cylinder-cylinder intersections agree with OCCT", "[occt][intersect]")
{
    struct Case {
        Cylinder Other;
        IntAna_ResultType Type;
        std::size_t Branches;
    };
    const Cylinder a{Point3{}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true};
    const std::vector<Case> cases{
        {Cylinder{Point3{5.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 2},
        {Cylinder{Point3{6.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 1},
        {Cylinder{Point3{1.0, 0.0, 0.0}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Line, 1},
        {Cylinder{Point3{10.0, 0.0, 0.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true},
         IntAna_Empty, 0},
        {Cylinder{Point3{0.0, 0.0, 9.0}, 3.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}, true},
         IntAna_Same, 0},
    };
    for (const Case& c : cases) {
        const auto result = Intersect(Surface(a), Surface(c.Other), ctx);
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.size() == c.Branches);
        REQUIRE(result->Overlaps.empty() == (c.Type != IntAna_Same));
        const IntAna_QuadQuadGeo occt(ToOcct(a), ToOcct(c.Other), Precision::Confusion());
        REQUIRE(occt.TypeInter() == c.Type);
        RequireOnOcct(*result, occt, 1e-9);
    }
}
