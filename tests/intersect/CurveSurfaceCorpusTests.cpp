// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every case of corpus/curve-surface/ (tools/generate_curve_surface_corpus.py):
// a curve of each built-in kind crossing a surface of each built-in kind at a
// point P known in closed form, at 90, 20 and 1 degrees to the tangent plane.
// The result has a decided crossing at P whose curve parameter is the case's,
// and every point it reports lies on the curve, in the requested domain, and
// on the surface at its reported parameters. Other crossings of the pair are
// not the case's claim. A case marked may_fail_closed (a known limit, named
// in the case) may instead fail with IntersectionFailed, never answer
// wrongly.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"
#include "CurveCorpus.hpp"

using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

/// How far the crossing may be from P, and its parameter from the case's,
/// per unit of scale: the inputs are rounded to 1e-12, and a grazing
/// crossing moves by that over the sine of its angle.
constexpr double Reach = 1e-9;

/// How far a reported point may be from the curve and the surface at its
/// parameters, per unit of scale.
constexpr double OnBoth = 1e-9;

std::vector<double> NumbersOf(const JsonValue& v)
{
    std::vector<double> numbers;
    for (const JsonValue& n : v.AsArray()) {
        numbers.push_back(n.AsNumber());
    }
    return numbers;
}

/// The case's surface, as tools/generate_curve_surface_corpus.py writes it.
Surface SurfaceOf(const JsonValue& s)
{
    using StereonTest::PointFrom;
    using StereonTest::VectorFrom;
    const std::string& type = s["type"].AsString();
    if (type == "bspline") {
        std::vector<Point3> poles;
        for (const JsonValue& p : s["poles"].AsArray()) {
            poles.push_back(PointFrom(p));
        }
        auto made = Stereon::MakeBSplineSurface(
            static_cast<int>(s["degree_u"].AsNumber()), static_cast<int>(s["degree_v"].AsNumber()),
            NumbersOf(s["knots_u"]), NumbersOf(s["knots_v"]), std::move(poles),
            s.Has("weights") ? NumbersOf(s["weights"]) : std::vector<double>{});
        REQUIRE(made.has_value());
        return Surface{*std::move(made)};
    }
    const Vec3 x = VectorFrom(s["x_axis"]);
    const Vec3 y = VectorFrom(s["y_axis"]);
    if (type == "plane") {
        return Surface{Stereon::Plane{PointFrom(s["origin"]), x, y}};
    }
    if (type == "sphere") {
        return Surface{Stereon::Sphere{PointFrom(s["center"]), s["radius"].AsNumber(), x, y, true}};
    }
    if (type == "cylinder") {
        return Surface{
            Stereon::Cylinder{PointFrom(s["origin"]), s["radius"].AsNumber(), x, y, true}};
    }
    if (type == "cone") {
        return Surface{Stereon::Cone{PointFrom(s["origin"]), s["radius"].AsNumber(),
                                     s["semi_angle"].AsNumber(), x, y, true}};
    }
    if (type == "torus") {
        return Surface{Stereon::Torus{PointFrom(s["center"]), s["major"].AsNumber(),
                                      s["minor"].AsNumber(), x, y, true}};
    }
    REQUIRE(type == "ellipsoid");
    return Surface{Stereon::Ellipsoid{PointFrom(s["center"]), s["a"].AsNumber(), s["b"].AsNumber(),
                                      s["c"].AsNumber(), x, y, true}};
}

std::vector<std::filesystem::path> CaseFiles()
{
    const std::filesystem::path root = std::filesystem::path(STEREON_CORPUS_DIR) / "curve-surface";
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

} // namespace

TEST_CASE("every curve-surface corpus case crosses at its point",
          "[intersect][curve-surface][corpus]")
{
    const auto files = CaseFiles();
    // Seven curve kinds, eight surface kinds, three angles.
    REQUIRE(files.size() == 168);
    std::size_t refused = 0;
    for (const auto& file : files) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        INFO(c["id"].AsString());
        const JsonValue& inputs = c["inputs"];
        const Stereon::Curve curve = StereonTest::CurveFrom(inputs["curve"]);
        const Stereon::Interval domain = StereonTest::DomainFrom(inputs["curve"]);
        const Surface surface = SurfaceOf(inputs["surface"]);
        const Point3 at = StereonTest::PointFrom(c["expected"]["point"]);
        const double t = c["expected"]["parameter"].AsNumber();
        const double scale = 1.0 + Distance(at, Point3{});

        const auto r = Intersect(curve, domain, surface, Domain(surface), IntersectTest::ctx);
        INFO(IntersectTest::ErrorText(r));
        if (!r && c["expected"].Has("may_fail_closed")) {
            // A known limit (the case's known_limit): refused, never wrong.
            CHECK(r.error().GetCode() == Stereon::ErrorCode::IntersectionFailed);
            ++refused;
            continue;
        }
        REQUIRE(r.has_value());
        REQUIRE(r->Overlaps.empty());
        const auto crossing = std::ranges::find_if(r->Points, [&](const auto& p) {
            return std::abs(p.Parameter - t) <= Reach * scale &&
                   Distance(p.Position, at) <= Reach * scale;
        });
        REQUIRE(crossing != r->Points.end());
        CHECK(crossing->Decided);
        CHECK(crossing->Contact == Stereon::ContactType::Crossing);
        for (const auto& p : r->Points) {
            INFO("at the curve's " << p.Parameter);
            CHECK(Contains(domain, p.Parameter));
            CHECK(Distance(Eval(curve, p.Parameter), p.Position) <= OnBoth * scale);
            CHECK(Distance(Eval(surface, p.OnSurface), p.Position) <= OnBoth * scale);
        }
    }
    // Fifteen cases are marked: the grazing ones on B-spline surfaces and a
    // line grazing a cone.
    CHECK(refused <= 15);
}
