// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every case of corpus/curved-ssi/, in both orders: the pairs with a cone, a
// torus or an ellipsoid (ADR-0010's corpus), whose answers come from closed
// forms in tools/generate_surface_corpus.py. A result has the case's curves
// and points, passes through every witness and stays on both surfaces. A case
// that names may_fail_closed may instead refuse with an error, never answer
// wrongly.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::ContactType;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR) / "curved-ssi";

/// Segments per branch when a witness is looked for.
constexpr std::size_t Samples = 2048;
/// How near a witness must lie to the result, per unit of the case's scale:
/// well above the sagitta of the samples, far below any feature of a case.
constexpr double WitnessReach = 1e-3;

Point3 PointOf(const JsonValue& v)
{
    return Point3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Vec3 VectorOf(const JsonValue& v)
{
    return Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Surface SurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    const Vec3 x = VectorOf(s["x_axis"]);
    const Vec3 y = VectorOf(s["y_axis"]);
    if (type == "plane") {
        return Surface{Stereon::Plane{PointOf(s["origin"]), x, y}};
    }
    if (type == "sphere") {
        return Surface{Stereon::Sphere{PointOf(s["center"]), s["radius"].AsNumber(), x, y, true}};
    }
    if (type == "cylinder") {
        return Surface{Stereon::Cylinder{PointOf(s["origin"]), s["radius"].AsNumber(), x, y, true}};
    }
    if (type == "cone") {
        return Surface{Stereon::Cone{PointOf(s["origin"]), s["radius"].AsNumber(),
                                     s["semi_angle"].AsNumber(), x, y, true}};
    }
    if (type == "torus") {
        return Surface{Stereon::Torus{PointOf(s["center"]), s["major"].AsNumber(),
                                      s["minor"].AsNumber(), x, y, true}};
    }
    return Surface{Stereon::Ellipsoid{PointOf(s["center"]), s["a"].AsNumber(), s["b"].AsNumber(),
                                      s["c"].AsNumber(), x, y, true}};
}

Box2 BoxOf(const JsonValue& s)
{
    if (!s.Has("box")) {
        return Box2::Everything();
    }
    const JsonValue& box = s["box"];
    return Box2{Interval{box[0][0].AsNumber(), box[0][1].AsNumber()},
                Interval{box[1][0].AsNumber(), box[1][1].AsNumber()}};
}

double DistanceToSegment(const Point3& p, const Point3& a, const Point3& b)
{
    const Vec3 along = b - a;
    const double square = Dot(along, along);
    if (square == 0.0) {
        return Distance(p, a);
    }
    const double s = std::clamp(Dot(p - a, along) / square, 0.0, 1.0);
    return Distance(p, a + s * along);
}

/// How far a point is from the nearest component of a result.
double DistanceToResult(const Stereon::SurfaceSurfaceIntersection& r, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Stereon::IntersectionPoint& point : r.Points) {
        nearest = std::min(nearest, Distance(p, point.Position));
    }
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        const std::vector<double> t = IntersectTest::SamplesOf(branch.Range, Samples);
        Point3 previous = Eval(branch.Curve3d, t.front());
        for (std::size_t k = 1; k < t.size(); ++k) {
            const Point3 next = Eval(branch.Curve3d, t[k]);
            nearest = std::min(nearest, DistanceToSegment(p, previous, next));
            previous = next;
        }
    }
    return nearest;
}

void CheckResult(const Stereon::SurfaceSurfaceIntersection& r, const JsonValue& expected,
                 const Surface& a, const Surface& b, const double tolerance, const double scale)
{
    if (expected["overlap"].AsBool()) {
        CHECK_FALSE(r.Overlaps.empty());
        return;
    }
    CHECK(r.Overlaps.empty());
    CHECK(r.Points.size() == static_cast<std::size_t>(expected["points"].AsNumber()));
    CHECK(IntersectTest::ConnectivityOf(r, 2.0 * tolerance).Curves ==
          static_cast<std::size_t>(expected["curves"].AsNumber()));
    const ContactType contact = expected["contact"].AsString() == "tangential"
                                    ? ContactType::Tangential
                                    : ContactType::Crossing;
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        CHECK(branch.Contact == contact);
        CHECK(branch.IsCertified());
        CHECK(IntersectTest::PcurveDeviation(branch, a, b, 100) <= branch.Bound + tolerance);
    }
    for (const JsonValue& witness : expected["witnesses"].AsArray()) {
        INFO("witness " << witness[0].AsNumber() << ' ' << witness[1].AsNumber() << ' '
                        << witness[2].AsNumber());
        CHECK(DistanceToResult(r, PointOf(witness)) <= WitnessReach * scale);
    }
}

bool FailsClosed(const ErrorCode code)
{
    return code == ErrorCode::IntersectionFailed || code == ErrorCode::CertificationFailed ||
           code == ErrorCode::NotSupported || code == ErrorCode::NotImplemented;
}

void CheckCase(const JsonValue& c, const bool swapped)
{
    INFO(c["id"].AsString() << (swapped ? " swapped" : ""));
    const JsonValue& surfaces = c["inputs"]["surfaces"];
    const JsonValue& first = surfaces[swapped ? 1 : 0];
    const JsonValue& second = surfaces[swapped ? 0 : 1];
    const Surface a = SurfaceOf(first);
    const Surface b = SurfaceOf(second);
    const double tolerance = c["inputs"]["tolerance"].AsNumber();
    const auto result = Intersect(a, BoxOf(first), b, BoxOf(second), tolerance, IntersectTest::ctx);
    const JsonValue& expected = c["expected"];
    if (!result) {
        INFO(IntersectTest::ErrorText(result));
        CHECK((expected.Has("may_fail_closed") && FailsClosed(result.error().GetCode())));
        return;
    }
    CheckResult(*result, expected, a, b, tolerance, c["inputs"]["scale"].AsNumber());
}

std::vector<std::filesystem::path> CaseFiles()
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(corpus)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

} // namespace

TEST_CASE("each curved-ssi configuration, tilted, is answered as its closed form says",
          "[intersect][surface-surface][corpus]")
{
    // One placement of each configuration, the one that found the most, in
    // one order: the whole suite takes minutes in a debug build.
    std::size_t count = 0;
    for (const std::filesystem::path& file : CaseFiles()) {
        if (file.stem().string().ends_with("-tilted-x1")) {
            CheckCase(StereonTest::ReadJsonFile(file), false);
            ++count;
        }
    }
    REQUIRE(count >= 100);
}

TEST_CASE("every curved-ssi corpus case is answered as its closed form says",
          "[.][exhaustive][surface-surface][corpus]")
{
    const std::vector<std::filesystem::path> files = CaseFiles();
    REQUIRE(files.size() >= 600);
    for (const std::filesystem::path& file : files) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        CheckCase(c, false);
        CheckCase(c, true);
    }
}
