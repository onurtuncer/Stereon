// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every case of corpus/planar-curves/, in both orders: each unordered pair of
// Line2d, Circle2d, Ellipse2d, Hyperbola2d, Parabola2d and BSplineCurve2d
// (polynomial and rational), crossing, touching with order two and four,
// 3 2^-40 either side of a tangency and coinciding, whose answers come from
// closed forms in tools/generate_planar_corpus.py, never from the kernel. A
// result has the case's points, each within its reach of the expected
// position with the expected contact, on both curves at its parameters, and
// the case's overlaps, whose ends are the expected ends and whose stretches
// lie on both curves.
//
// Two gates. Safety: an answer is never wrong, and a refusal is a fail-closed
// error (IntersectionFailed or NotImplemented), never another one; this one
// must pass. Completeness: only a case that names may_fail_closed, a
// limitation the v0.2 ledger documents, is refused; it is [!mayfail] while
// the kernel refuses cases the ledger does not list (the test's comment).

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::ContactType;
using Stereon::Curve2d;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::PlanarIntersection;
using Stereon::PlanarPoint;
using Stereon::Point2;
using Stereon::Vec2;
using StereonTest::JsonValue;

namespace {

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR);

/// Segments of an overlap's stretch when its points are looked for on the
/// other curve.
constexpr std::size_t Samples = 512;

Point2 PointOf(const JsonValue& v)
{
    return Point2{v[0].AsNumber(), v[1].AsNumber()};
}

Vec2 VectorOf(const JsonValue& v)
{
    return Vec2{v[0].AsNumber(), v[1].AsNumber()};
}

std::vector<double> NumbersOf(const JsonValue& v)
{
    std::vector<double> numbers;
    for (const JsonValue& n : v.AsArray()) {
        numbers.push_back(n.AsNumber());
    }
    return numbers;
}

/// A curve as tools/generate_planar_corpus.py writes it, through the
/// validating factories.
Curve2d CurveOf(const JsonValue& c)
{
    const std::string& type = c["type"].AsString();
    if (type == "line") {
        auto made = Stereon::MakeLine(PointOf(c["origin"]), VectorOf(c["direction"]));
        REQUIRE(made.has_value());
        return Curve2d{*made};
    }
    if (type == "bspline") {
        std::vector<Point2> poles;
        for (const JsonValue& p : c["poles"].AsArray()) {
            poles.push_back(PointOf(p));
        }
        auto made = Stereon::MakeBSplineCurve(
            static_cast<int>(c["degree"].AsNumber()), NumbersOf(c["knots"]), std::move(poles),
            c.Has("weights") ? NumbersOf(c["weights"]) : std::vector<double>{});
        REQUIRE(made.has_value());
        return Curve2d{*std::move(made)};
    }
    const Vec2 x = VectorOf(c["x_axis"]);
    const Vec2 y = VectorOf(c["y_axis"]);
    if (type == "circle") {
        auto made = Stereon::MakeCircle(PointOf(c["center"]), c["radius"].AsNumber(), x, y);
        REQUIRE(made.has_value());
        return Curve2d{*made};
    }
    if (type == "ellipse") {
        auto made = Stereon::MakeEllipse(PointOf(c["center"]), c["major"].AsNumber(),
                                         c["minor"].AsNumber(), x, y);
        REQUIRE(made.has_value());
        return Curve2d{*made};
    }
    if (type == "hyperbola") {
        auto made = Stereon::MakeHyperbola(PointOf(c["center"]), c["major"].AsNumber(),
                                           c["minor"].AsNumber(), x, y);
        REQUIRE(made.has_value());
        return Curve2d{*made};
    }
    auto made = Stereon::MakeParabola(PointOf(c["vertex"]), c["focal"].AsNumber(), x, y);
    REQUIRE(made.has_value());
    return Curve2d{*made};
}

Interval DomainOf(const JsonValue& c, const Curve2d& curve)
{
    if (!c.Has("domain")) {
        return Domain(curve);
    }
    return Interval{c["domain"][0].AsNumber(), c["domain"][1].AsNumber()};
}

double DistanceToSegment(const Point2& p, const Point2& a, const Point2& b)
{
    const Vec2 along = b - a;
    const double square = Dot(along, along);
    if (square == 0.0) {
        return Distance(p, a);
    }
    const double s = std::clamp(Dot(p - a, along) / square, 0.0, 1.0);
    return Distance(p, a + s * along);
}

/// How far a point is from a curve's stretch, by its polyline.
double DistanceToStretch(const Curve2d& curve, const Interval& range, const Point2& p)
{
    const std::vector<double> t = IntersectTest::SamplesOf(range, Samples);
    double nearest = std::numeric_limits<double>::infinity();
    Point2 previous = Eval(curve, t.front());
    for (std::size_t k = 1; k < t.size(); ++k) {
        const Point2 next = Eval(curve, t[k]);
        nearest = std::min(nearest, DistanceToSegment(p, previous, next));
        previous = next;
    }
    return nearest;
}

void CheckOverlaps(const PlanarIntersection& r, const JsonValue& expected, const Curve2d& a,
                   const Curve2d& b)
{
    CHECK(r.Points.empty());
    REQUIRE_FALSE(r.Overlaps.empty());
    double reach = 0.0;
    for (const JsonValue& overlap : expected.AsArray()) {
        reach = std::max(reach, overlap["within"].AsNumber());
    }
    // Every expected end is an end of a stretch found.
    for (const JsonValue& overlap : expected.AsArray()) {
        for (const char* end : {"from", "to"}) {
            const Point2 p = PointOf(overlap[end]);
            INFO("overlap end " << p.X << ' ' << p.Y);
            CHECK(std::ranges::any_of(r.Overlaps, [&](const Stereon::PlanarOverlap& o) {
                return std::min(Distance(Eval(a, o.OnA.GetLower()), p),
                                Distance(Eval(a, o.OnA.GetUpper()), p)) <= reach;
            }));
        }
    }
    // Every stretch found lies on both curves. The sagitta of the polyline
    // is far below the reach at the corpus's curvatures.
    for (const Stereon::PlanarOverlap& o : r.Overlaps) {
        for (const double s : {0.0, 0.125, 0.5, 0.875, 1.0}) {
            const double ta = o.OnA.GetLower() + s * Width(o.OnA);
            const double tb = o.OnB.GetLower() + s * Width(o.OnB);
            CHECK(DistanceToStretch(b, o.OnB, Eval(a, ta)) <= reach);
            CHECK(DistanceToStretch(a, o.OnA, Eval(b, tb)) <= reach);
        }
    }
}

void CheckPoints(const PlanarIntersection& r, const JsonValue& expected, const Curve2d& a,
                 const Curve2d& b)
{
    CHECK(r.Overlaps.empty());
    const auto& points = expected["points"].AsArray();
    CHECK(r.Points.size() == points.size());
    const bool decided = expected.Has("decided") && expected["decided"].AsBool();
    std::vector<bool> used(r.Points.size(), false);
    double reach = 0.0;
    for (const JsonValue& point : points) {
        const Point2 position = PointOf(point["position"]);
        const double within = point["within"].AsNumber();
        reach = std::max(reach, within);
        const ContactType contact = point["contact"].AsString() == "tangential"
                                        ? ContactType::Tangential
                                        : ContactType::Crossing;
        INFO("expected " << point["contact"].AsString() << " point " << position.X << ' '
                         << position.Y << " within " << within);
        std::size_t best = r.Points.size();
        for (std::size_t k = 0; k < r.Points.size(); ++k) {
            if (!used[k] && Distance(r.Points[k].Position, position) <= within &&
                (best == r.Points.size() ||
                 Distance(r.Points[k].Position, position) <
                     Distance(r.Points[best].Position, position))) {
                best = k;
            }
        }
        CHECK(best < r.Points.size());
        if (best == r.Points.size()) {
            continue;
        }
        used[best] = true;
        CHECK(r.Points[best].Contact == contact);
        if (decided) {
            CHECK(r.Points[best].Decided);
        }
    }
    for (const PlanarPoint& p : r.Points) {
        INFO("found point " << p.Position.X << ' ' << p.Position.Y << " at " << p.ParameterA << ", "
                            << p.ParameterB);
        CHECK(Distance(Eval(a, p.ParameterA), p.Position) <= reach);
        CHECK(Distance(Eval(b, p.ParameterB), p.Position) <= reach);
    }
}

bool FailsClosed(const ErrorCode code)
{
    return code == ErrorCode::IntersectionFailed || code == ErrorCode::NotImplemented;
}

/// How many cases were answered, refused as their may_fail_closed allows,
/// and refused without it, for the log.
struct Tally {
    std::size_t Answered = 0;
    std::size_t Refused = 0;
    std::size_t Undocumented = 0;
};

enum class Gate : std::uint8_t { Safety, Completeness };

void CheckCase(const JsonValue& c, const bool swapped, const Gate gate, Tally& tally)
{
    INFO(c["id"].AsString() << (swapped ? " swapped" : ""));
    const JsonValue& curves = c["inputs"]["curves"];
    const JsonValue& first = curves[swapped ? 1 : 0];
    const JsonValue& second = curves[swapped ? 0 : 1];
    const Curve2d a = CurveOf(first);
    const Curve2d b = CurveOf(second);
    const auto result = Intersect(a, DomainOf(first, a), b, DomainOf(second, b), IntersectTest::ctx);
    const JsonValue& expected = c["expected"];
    if (!result) {
        INFO(IntersectTest::ErrorText(result));
        const bool documented = expected.Has("may_fail_closed");
        ++(documented ? tally.Refused : tally.Undocumented);
        if (gate == Gate::Safety) {
            CHECK(FailsClosed(result.error().GetCode()));
        } else {
            CHECK((documented && FailsClosed(result.error().GetCode())));
        }
        return;
    }
    ++tally.Answered;
    if (gate == Gate::Completeness) {
        return;
    }
    if (!expected["overlaps"].AsArray().empty()) {
        CheckOverlaps(*result, expected["overlaps"], a, b);
    } else {
        CheckPoints(*result, expected, a, b);
    }
}

std::vector<std::filesystem::path> CaseFiles()
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry :
         std::filesystem::recursive_directory_iterator(corpus / "planar-curves")) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

/// One placement of each configuration: the rotated one at scale 1, or the
/// far one where a rotation's rounding would change the answer.
std::vector<std::filesystem::path> OnePlacementEach()
{
    std::vector<std::filesystem::path> picked;
    for (const std::filesystem::path& file : CaseFiles()) {
        const std::string stem = file.stem().string();
        bool pick = stem.ends_with("-rotated-x1");
        if (stem.ends_with("-far-x1")) {
            std::filesystem::path rotated = file;
            rotated.replace_filename(stem.substr(0, stem.size() - 6) + "rotated-x1.json");
            pick = !std::filesystem::exists(rotated);
        }
        if (pick) {
            picked.push_back(file);
        }
    }
    return picked;
}

void RunAll(const std::vector<std::filesystem::path>& files, const bool bothWays, const Gate gate)
{
    Tally tally;
    for (const std::filesystem::path& file : files) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        CheckCase(c, false, gate, tally);
        if (bothWays) {
            CheckCase(c, true, gate, tally);
        }
    }
    WARN(tally.Answered << " answered, " << tally.Refused << " refused as documented, "
                        << tally.Undocumented << " refused without a documented limitation");
}

} // namespace

TEST_CASE("no planar-curves configuration, rotated or else far, is answered wrongly",
          "[intersect][planar][corpus]")
{
    const std::vector<std::filesystem::path> files = OnePlacementEach();
    REQUIRE(files.size() >= 140);
    RunAll(files, false, Gate::Safety);
}

// [!mayfail] until the kernel answers what the ledger does not list as a
// limitation (2026-11-12, main 227d79a; nothing is answered wrongly):
// every contact of order four between a B-spline and a line or a conic is
// refused, the subdivision budget spent at the touch (line-spline,
// {circle,ellipse,hyperbola,parabola}-{spline,rational}/osculating), although
// the ledger says touches are reported bracketed; and two B-splines at the far
// placement are refused as "neither apart nor a single crossing": the near
// pairs (rational-rational/near-far-*, spline-rational/near-far-*) and one
// transversal crossing (spline-rational/crossing-far-x1). Remove the tag once
// they are answered or the ledger lists them and their cases name
// may_fail_closed.
TEST_CASE("only the planar-curves configurations the ledger lists are refused",
          "[intersect][planar][corpus][!mayfail]")
{
    RunAll(OnePlacementEach(), false, Gate::Completeness);
}

TEST_CASE("no planar-curves corpus case is answered wrongly, both ways",
          "[.][exhaustive][intersect][planar][corpus]")
{
    const std::vector<std::filesystem::path> files = CaseFiles();
    REQUIRE(files.size() >= 800);
    RunAll(files, true, Gate::Safety);
}

// [!mayfail] for the refusals listed above, at every placement and scale.
TEST_CASE("only the planar-curves corpus cases the ledger lists are refused, both ways",
          "[.][exhaustive][intersect][planar][corpus][!mayfail]")
{
    RunAll(CaseFiles(), true, Gate::Completeness);
}
