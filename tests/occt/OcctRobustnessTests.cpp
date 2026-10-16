// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Cases where OCCT 8.0.1 gives a wrong answer and Stereon a right one,
// found by the adversarial hunt in RobustnessHunt.cpp. Each test checks
// Stereon against a reference independent of both kernels (exact
// arithmetic, a refined polyline, or a dense sample) and records how far
// OCCT is off. CI pins OCCT 8.0.1 (tools/occt/build_occt.sh), so the OCCT
// figures are reproducible; a later OCCT that fixes a case will fail its
// OCCT check, which should then be retired, not loosened.

#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

import Stereon.Surface;

#include "OcctConversion.hpp"
#include "ReferenceLengths.hpp"

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using StereonTest::PolylineLength;
using StereonTest::PolylineReference;
using StereonTest::ToOcct;

namespace {

/// The tolerance asked of OCCT's arc-length routines in every case.
constexpr double OcctTolerance = 1e-10;

} // namespace

TEST_CASE("the point at a length on a polyline with a stall and a short span",
          "[occt][robustness][curve]")
{
    // A polyline with a repeated point (zero speed on a span) and a span
    // 3.75e-5 wide: GCPnts_AbscissaPoint lands 0.68 % of the length away.
    const BSplineCurve curve{1,
                             {0.0, 0.0, 3.7507542114219575e-05, 1.6168619916534712,
                              1.6168849627863298, 1.6174873544359984, 5.0360710934077577,
                              5.0360710934077577},
                             {Point3{5.3062002650719631, -8.4691886009582724, 1.4037049161974267},
                              Point3{-0.66736793906986591, 3.9355543731033968, -9.2111164209898266},
                              Point3{-8.2174567020354221, -0.65455466621195413, 2.6068698691036669},
                              Point3{-4.6409254184243824, 0.93915723417673114, 4.9105598686212204},
                              Point3{-4.6409254184243824, 0.93915723417673114, 4.9105598686212204},
                              Point3{-6.3772998899613782, 2.92557592571832, -4.724074108776608}},
                             {}};
    const double target = 36.366490109214404;
    const Interval d = Domain(curve);
    const double total = PolylineLength(curve, d.GetUpper());

    const double ours = *Stereon::ParameterAtArcLength(curve, d.GetLower(), target);
    CHECK(std::abs(PolylineLength(curve, ours) - target) <= 1e-10 * total);

    const GeomAdaptor_Curve adaptor(ToOcct(curve));
    const GCPnts_AbscissaPoint occt(OcctTolerance, adaptor, target, d.GetLower());
    REQUIRE(occt.IsDone());
    const double occtError = std::abs(PolylineLength(curve, occt.Parameter()) - target) / total;
    CHECK(occtError > 1e-3); // OCCT 8.0.1: about 6.8e-3
}

TEST_CASE("the length of a rational cubic with weights over six orders of magnitude",
          "[occt][robustness][curve]")
{
    // Weights from 1.1e-3 to 471: the speed varies by orders of magnitude
    // along the curve. GCPnts_AbscissaPoint::Length is off by 1.4 %, and so
    // is the point it finds at a given length.
    const std::vector<Point3> poles{
        Point3{5.8689671605911453, 5.9450512340941248, -6.932323430012735},
        Point3{-1.1655730438767105, 7.5887693818371993, 7.7635696098965639},
        Point3{2.3665850768645207, 8.915181968703088, 1.9412140772297288},
        Point3{1.4902564531983149, 8.077004401390198, -6.2061902053727618},
        Point3{-5.2392474061459415, 6.5054728950109748, 5.0641468150448414},
        Point3{-4.9060971355083254, 2.7275189141375797, 9.2833315158292926},
        Point3{-7.5953010386799562, 7.0540168776767551, 3.7160245713162752},
        Point3{7.9887065358740372, 5.266514728939697, 2.9253933909562591},
        Point3{-1.784317146339216, -8.9808445221584083, 1.4715857410935147}};
    const std::vector<double> weights{
        0.0011163267166932475, 470.6209147013704,   161.78160387556696,
        0.0019569416614698058, 0.89682975549589061, 0.20921122426219615,
        0.0011346594117175202, 378.34193132462798,  28.360278609557295};
    const BSplineCurve curve{
        3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 6.0, 6.0}, poles, weights};
    const Interval d = Domain(curve);
    const double reference = PolylineReference(curve, d.GetLower(), d.GetUpper());

    const auto ours = Stereon::ArcLength(curve, d);
    CHECK(std::abs(ours->Length - reference) <= 1e-9 * reference);
    CHECK(ours->ErrorEstimate <= 1e-9 * reference);

    const GeomAdaptor_Curve adaptor(ToOcct(curve));
    const double occt =
        GCPnts_AbscissaPoint::Length(adaptor, d.GetLower(), d.GetUpper(), OcctTolerance);
    CHECK(std::abs(occt - reference) / reference > 1e-3); // OCCT 8.0.1: about 1.4e-2
}

TEST_CASE("the closest point on the flank of a narrow bump", "[occt][robustness][surface]")
{
    // A bicubic sheet over a 14 x 14 net with one pole raised 29.5 mm: a
    // narrow bump. GeomAPI_ProjectPointOnSurf stops at a local minimum
    // 0.278 mm away; the closest point is 0.230 mm away.
    constexpr int n = 14;
    std::vector<double> knots(4, 0.0);
    for (int i = 1; i < n - 3; ++i) {
        knots.push_back(static_cast<double>(i));
    }
    knots.insert(knots.end(), 4, static_cast<double>(n - 3));
    BSplineSurface surface{3, 3, knots, knots, {}, {}};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double z = i == 9 && j == 10 ? 29.502900254827399 : 0.0;
            surface.Poles.push_back({static_cast<double>(i), static_cast<double>(j), z});
        }
    }
    const Point3 p{9.5837785119550869, 10.532644902983105, 8.6806190559984273};

    const double ours = Distance(Eval(surface, Project(surface, p)), p);
    // No point of a dense grid over the bump's support is closer.
    double sampled = ours + 1.0;
    constexpr int samples = 2000;
    for (int i = 0; i <= samples; ++i) {
        for (int j = 0; j <= samples; ++j) {
            const Point2 uv{4.0 + 6.0 * i / samples, 5.0 + 6.0 * j / samples};
            sampled = std::min(sampled, Distance(Eval(surface, uv), p));
        }
    }
    CHECK(ours <= sampled + 1e-9);
    CHECK(std::abs(ours - 0.23042243723396058) <= 1e-9);

    GeomAPI_ProjectPointOnSurf occt(ToOcct(p), ToOcct(surface));
    REQUIRE(occt.NbPoints() > 0);
    CHECK(occt.LowerDistance() > ours + 0.04); // OCCT 8.0.1: 0.278
}
