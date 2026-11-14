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
#include <numbers>
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
                (best == r.Points.size() || Distance(r.Points[k].Position, position) <
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
    const auto result =
        Intersect(a, DomainOf(first, a), b, DomainOf(second, b), IntersectTest::ctx);
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

// --- MR-04: equivalent reparameterization ---------------------------------

/// One operand in an equivalent chart, built exactly from the case's data,
/// not by the kernel: a line's origin moved two directions on and its
/// direction doubled, t = 2 t' + 2; a circle's frame turned a quarter,
/// t = t' + pi/2; an ellipse's turned half, t = t' + pi; a hyperbola's or a
/// parabola's y axis negated, t = -t'; a B-spline reversed (knots
/// reflected, poles and weights reversed), t = a + b - t'. Every step is
/// exact on dyadic data: a negation, a doubling or a reflection of dyadic
/// knots.
struct Reparam {
    Curve2d Curve;
    Interval Domain;
    double (*Back)(double, double) = nullptr; ///< t from t' and the spline's a + b
    double Sum = 0.0;
};

double BackLine(const double t, double /*unused*/)
{
    return 2.0 * t + 2.0;
}

double BackCircle(const double t, double /*unused*/)
{
    return t + std::numbers::pi / 2.0;
}

double BackEllipse(const double t, double /*unused*/)
{
    return t + std::numbers::pi;
}

double BackNegated(const double t, double /*unused*/)
{
    return -t;
}

double BackSpline(const double t, const double sum)
{
    return sum - t;
}

/// The requested domain's ends, when the case bounds it.
struct Bounds {
    bool Bounded = false;
    double Lo = 0.0;
    double Hi = 0.0;
};

Bounds BoundsOf(const JsonValue& c)
{
    if (!c.Has("domain")) {
        return {};
    }
    return {true, c["domain"][0].AsNumber(), c["domain"][1].AsNumber()};
}

Reparam ReparamLine(const JsonValue& c, const Bounds& d)
{
    const Point2 o = PointOf(c["origin"]);
    const Vec2 v = VectorOf(c["direction"]);
    Reparam r;
    r.Curve = Curve2d{*Stereon::MakeLine(o + 2.0 * v, 2.0 * v)};
    r.Back = &BackLine;
    r.Domain =
        d.Bounded ? Interval{(d.Lo - 2.0) / 2.0, (d.Hi - 2.0) / 2.0} : Interval::Everything();
    return r;
}

Reparam ReparamSpline(const JsonValue& c, const Bounds& d)
{
    const std::vector<double> knots = NumbersOf(c["knots"]);
    const double sum = knots.front() + knots.back();
    std::vector<double> reflected;
    for (auto it = knots.rbegin(); it != knots.rend(); ++it) {
        reflected.push_back(sum - *it);
    }
    std::vector<Point2> poles;
    for (const JsonValue& p : c["poles"].AsArray()) {
        poles.insert(poles.begin(), PointOf(p));
    }
    std::vector<double> weights =
        c.Has("weights") ? NumbersOf(c["weights"]) : std::vector<double>{};
    std::ranges::reverse(weights);
    auto made =
        Stereon::MakeBSplineCurve(static_cast<int>(c["degree"].AsNumber()), std::move(reflected),
                                  std::move(poles), std::move(weights));
    REQUIRE(made.has_value());
    Reparam r;
    r.Curve = Curve2d{*std::move(made)};
    r.Back = &BackSpline;
    r.Sum = sum;
    r.Domain = d.Bounded ? Interval{sum - d.Hi, sum - d.Lo} : Domain(r.Curve);
    return r;
}

/// A closed conic turned by `turn` (a quarter for a circle, half for an
/// ellipse): its chart shifted by that angle.
Reparam ReparamTurned(Curve2d curve, double (*back)(double, double), const double turn,
                      const Bounds& d)
{
    Reparam r;
    r.Curve = std::move(curve);
    r.Back = back;
    r.Domain = d.Bounded ? Interval{d.Lo - turn, d.Hi - turn} : Domain(r.Curve);
    return r;
}

Reparam ReparamConic(const JsonValue& c, const Bounds& d)
{
    const std::string& type = c["type"].AsString();
    const Vec2 x = VectorOf(c["x_axis"]);
    const Vec2 y = VectorOf(c["y_axis"]);
    if (type == "circle") {
        return ReparamTurned(
            Curve2d{*Stereon::MakeCircle(PointOf(c["center"]), c["radius"].AsNumber(), y, -x)},
            &BackCircle, std::numbers::pi / 2.0, d);
    }
    if (type == "ellipse") {
        return ReparamTurned(
            Curve2d{*Stereon::MakeEllipse(PointOf(c["center"]), c["major"].AsNumber(),
                                          c["minor"].AsNumber(), -x, -y)},
            &BackEllipse, std::numbers::pi, d);
    }
    Reparam r;
    if (type == "hyperbola") {
        r.Curve = Curve2d{*Stereon::MakeHyperbola(PointOf(c["center"]), c["major"].AsNumber(),
                                                  c["minor"].AsNumber(), x, -y)};
    } else {
        r.Curve =
            Curve2d{*Stereon::MakeParabola(PointOf(c["vertex"]), c["focal"].AsNumber(), x, -y)};
    }
    r.Back = &BackNegated;
    r.Domain = d.Bounded ? Interval{-d.Hi, -d.Lo} : Interval::Everything();
    return r;
}

Reparam ReparamOf(const JsonValue& c)
{
    const std::string& type = c["type"].AsString();
    const Bounds d = BoundsOf(c);
    if (type == "line") {
        return ReparamLine(c, d);
    }
    if (type == "bspline") {
        return ReparamSpline(c, d);
    }
    return ReparamConic(c, d);
}

/// The largest reach of the case's expected points and overlaps.
double ReachOf(const JsonValue& expected)
{
    double reach = 0.0;
    for (const JsonValue& p : expected["points"].AsArray()) {
        reach = std::max(reach, p["within"].AsNumber());
    }
    for (const JsonValue& o : expected["overlaps"].AsArray()) {
        reach = std::max(reach, o["within"].AsNumber());
    }
    return reach;
}

/// Every parameter on the moved operand maps back, by the closed-form
/// change of parameter, to the original curve at the same place.
void CheckMappedBack(const PlanarIntersection& result, const Reparam& moved,
                     const Curve2d& original, const bool moveFirst, const double reach)
{
    for (const PlanarPoint& p : result.Points) {
        const double t = moveFirst ? p.ParameterA : p.ParameterB;
        CHECK(Distance(Eval(original, moved.Back(t, moved.Sum)), p.Position) <= reach);
    }
    for (const Stereon::PlanarOverlap& o : result.Overlaps) {
        const Interval& range = moveFirst ? o.OnA : o.OnB;
        for (const double t : {range.GetLower(), range.GetUpper()}) {
            CHECK(Distance(Eval(original, moved.Back(t, moved.Sum)), Eval(moved.Curve, t)) <=
                  reach);
        }
    }
}

/// A refusal never establishes a relation. Safety checks its error class;
/// the separate completeness gate requires a documented fixture limitation.
void CheckRefusal(const JsonValue& expected, const ErrorCode code, const Gate gate, Tally& tally)
{
    const bool documented = expected.Has("may_fail_closed");
    CHECK(FailsClosed(code));
    if (gate == Gate::Completeness) {
        CHECK(documented);
    }
    ++(documented ? tally.Refused : tally.Undocumented);
}

/// MR-04 for one case, one order, one operand moved to its equivalent
/// chart: the result in the new chart is the case's closed-form answer,
/// and every parameter on the moved operand maps back, by the closed-form
/// change of parameter, to a point of the original curve at the same
/// place.
void CheckReparam(const JsonValue& c, const bool swapped, const bool moveFirst, const Gate gate,
                  Tally& tally)
{
    INFO(c["id"].AsString() << (swapped ? " swapped" : "") << (moveFirst ? " A'" : " B'")
                            << " reparametrized");
    const JsonValue& curves = c["inputs"]["curves"];
    const JsonValue& first = curves[swapped ? 1 : 0];
    const JsonValue& second = curves[swapped ? 0 : 1];
    const Curve2d originalA = CurveOf(first);
    const Curve2d originalB = CurveOf(second);
    const Reparam moved = ReparamOf(moveFirst ? first : second);
    const Curve2d& a = moveFirst ? moved.Curve : originalA;
    const Curve2d& b = moveFirst ? originalB : moved.Curve;
    const Interval domainA = moveFirst ? moved.Domain : DomainOf(first, originalA);
    const Interval domainB = moveFirst ? DomainOf(second, originalB) : moved.Domain;
    const auto result = Intersect(a, domainA, b, domainB, IntersectTest::ctx);
    const JsonValue& expected = c["expected"];
    if (!result) {
        INFO(IntersectTest::ErrorText(result));
        CheckRefusal(expected, result.error().GetCode(), gate, tally);
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
    CheckMappedBack(*result, moved, moveFirst ? originalA : originalB, moveFirst,
                    ReachOf(expected));
}

void RunReparam(const std::vector<std::filesystem::path>& files, const bool bothWays,
                const Gate gate = Gate::Safety)
{
    Tally tally;
    for (const std::filesystem::path& file : files) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        for (const bool swapped : {false, true}) {
            if (swapped && !bothWays) {
                continue;
            }
            CheckReparam(c, swapped, true, gate, tally);
            CheckReparam(c, swapped, false, gate, tally);
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

TEST_CASE("MR-04: each planar-curves configuration answers alike in an equivalent chart",
          "[intersect][planar][corpus][metamorphic]")
{
    // One placement of each configuration, one order, each operand in turn
    // in its equivalent chart (ReparamOf).
    const std::vector<std::filesystem::path> files = OnePlacementEach();
    REQUIRE(files.size() >= 140);
    RunReparam(files, false);
}

TEST_CASE("MR-04: every planar-curves corpus case answers alike in an equivalent chart, both ways",
          "[.][exhaustive][intersect][planar][corpus][metamorphic]")
{
    const std::vector<std::filesystem::path> files = CaseFiles();
    REQUIRE(files.size() >= 800);
    RunReparam(files, true);
}

// Incomplete evidence is recorded as MR implemented, never verified. This
// test exposes the documented completeness debt separately from safety.
TEST_CASE("MR-04: only documented planar configurations refuse an equivalent chart",
          "[intersect][planar][corpus][metamorphic][!mayfail]")
{
    RunReparam(OnePlacementEach(), false, Gate::Completeness);
}

TEST_CASE("MR-04: only documented planar corpus cases refuse equivalent charts, both ways",
          "[.][exhaustive][intersect][planar][corpus][metamorphic][!mayfail]")
{
    RunReparam(CaseFiles(), true, Gate::Completeness);
}
