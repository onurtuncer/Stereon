// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The cases of the Phase 2 benchmarks (docs/benchmarks.md, B2.1 to B2.3) and
// the closed-form answers both libraries are scored against, shared by the
// Stereon rows (bench/IntersectBench.cpp) and their OCCT twins
// (bench/occt/OcctIntersectBench.cpp), so that both solve the same problems
// and are judged by the same rules:
//
//   B2.1  corpus/planar-curves/ (every pair without an overlap), and lines
//         through the exact rational NURBS sphere and cylinder of
//         corpus/nurbs-ssi/, whose roots are a quadratic's;
//   B2.2  corpus/nurbs-ssi/ and corpus/curved-ssi/ by family, every case but
//         coincident surfaces (GeomAPI_IntSS has no overlap outcome);
//   B2.3  corpus/analytic-ssi/'s Tier A pairs of planes, spheres and
//         cylinders, by pair.
//
// Each library's answer is reduced to points and sampled curves (Found) and
// scored once, outside the timed loop: roots or branches missed and spurious
// against the closed form, cases left unresolved (a refusal, OCCT's IsDone()
// false or an exception), configurations that are not exact, and the largest
// distance of a sampled point from either surface (max_deviation, mm). The
// counters are named without the _mm suffix that marks a result both
// libraries must share (tools/check_release_benchmarks.py): they differ by
// design.
//
// Include after Json.hpp, <algorithm>, <cmath>, <cstddef>, <filesystem>,
// <limits>, <map>, <optional>, <stdexcept>, <string>, <utility>, <vector>
// and `import Stereon.Intersect;`.

#pragma once

namespace StereonBench::B2 {

using Stereon::Box2;
using Stereon::Curve2d;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::JsonValue;

/// Samples per curve of a result, and the reach of a witness or a sample,
/// and of joined curve ends, per unit of a case's scale: the reach of
/// tests/intersect/SurfaceCorpusTests.cpp and tests/occt/CorpusDifferentialTests.cpp.
inline constexpr int CurveSamples = 64;
inline constexpr double WitnessReach = 1e-3;
inline constexpr double JoinReach = 1e-4;
/// An unbounded line is sampled this far either side of the case's anchor,
/// per unit of scale; a parameter past Unbounded counts as unbounded (OCCT
/// trims the lines of an analytic intersection at about 1.5e9).
inline constexpr double LineWindow = 10.0;
inline constexpr double Unbounded = 1e7;
/// How near a root of B2.1's lines a found point must be.
inline constexpr double RootReach = 1e-7;

inline const Stereon::Context& Ctx()
{
    static const Stereon::Context context = Stereon::Context::Default();
    return context;
}

template <typename R>
auto Must(R made, const char* what)
{
    if (!made) {
        throw std::runtime_error(what);
    }
    return *std::move(made);
}

inline std::vector<JsonValue> CasesOf(const std::string& suite)
{
    const std::filesystem::path root = std::filesystem::path(STEREON_CORPUS_DIR) / suite;
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    std::vector<JsonValue> cases;
    for (const std::filesystem::path& file : files) {
        cases.push_back(StereonTest::ReadJsonFile(file));
    }
    return cases;
}

/// The family of a case ID: "planar-curves/line-circle/x" -> "line-circle".
inline std::string FamilyOf(const std::string& id)
{
    const std::size_t first = id.find('/');
    return id.substr(first + 1, id.rfind('/') - first - 1);
}

inline Point2 Point2Of(const JsonValue& v)
{
    return Point2{v[0].AsNumber(), v[1].AsNumber()};
}
inline Vec2 Vec2Of(const JsonValue& v)
{
    return Vec2{v[0].AsNumber(), v[1].AsNumber()};
}
inline Point3 Point3Of(const JsonValue& v)
{
    return Point3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}
inline Vec3 Vec3Of(const JsonValue& v)
{
    return Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}
inline std::vector<double> NumbersOf(const JsonValue& v)
{
    std::vector<double> numbers;
    for (const JsonValue& n : v.AsArray()) {
        numbers.push_back(n.AsNumber());
    }
    return numbers;
}

/// What a row scored, summed over its cases.
struct Score {
    double Cases = 0.0;
    double Expected = 0.0;   ///< roots, or connected curves and points, the closed form has
    double Missed = 0.0;     ///< of those, not found (witnesses not on the result included)
    double Spurious = 0.0;   ///< found, not in the closed form
    double Unresolved = 0.0; ///< cases with no answer
    double Inexact =
        0.0; ///< B2.2, B2.3: answered cases whose configuration is not the closed form's
    double MaxDeviation = 0.0;
};

// --- B2.1, plane curves -------------------------------------------------------

struct PlanarCase {
    Curve2d A;
    Interval DomainA;
    Curve2d B;
    Interval DomainB;
    std::vector<Point2> Roots;
    std::vector<double> Within;
};

inline Curve2d PlanarCurveOf(const JsonValue& c)
{
    const std::string& type = c["type"].AsString();
    if (type == "line") {
        return Curve2d{
            Must(Stereon::MakeLine(Point2Of(c["origin"]), Vec2Of(c["direction"])), "line")};
    }
    if (type == "bspline") {
        std::vector<Point2> poles;
        for (const JsonValue& p : c["poles"].AsArray()) {
            poles.push_back(Point2Of(p));
        }
        return Curve2d{Must(Stereon::MakeBSplineCurve(static_cast<int>(c["degree"].AsNumber()),
                                                      NumbersOf(c["knots"]), std::move(poles),
                                                      c.Has("weights") ? NumbersOf(c["weights"])
                                                                       : std::vector<double>{}),
                            "spline")};
    }
    const Vec2 x = Vec2Of(c["x_axis"]);
    const Vec2 y = Vec2Of(c["y_axis"]);
    if (type == "circle") {
        return Curve2d{Must(
            Stereon::MakeCircle(Point2Of(c["center"]), c["radius"].AsNumber(), x, y), "circle")};
    }
    if (type == "ellipse") {
        return Curve2d{Must(Stereon::MakeEllipse(Point2Of(c["center"]), c["major"].AsNumber(),
                                                 c["minor"].AsNumber(), x, y),
                            "ellipse")};
    }
    if (type == "hyperbola") {
        return Curve2d{Must(Stereon::MakeHyperbola(Point2Of(c["center"]), c["major"].AsNumber(),
                                                   c["minor"].AsNumber(), x, y),
                            "hyperbola")};
    }
    return Curve2d{Must(Stereon::MakeParabola(Point2Of(c["vertex"]), c["focal"].AsNumber(), x, y),
                        "parabola")};
}

inline Interval PlanarDomainOf(const JsonValue& c, const Curve2d& curve)
{
    return c.Has("domain") ? Interval{c["domain"][0].AsNumber(), c["domain"][1].AsNumber()}
                           : Domain(curve);
}

/// The pairs of corpus/planar-curves/ without an overlap: of the families with
/// a spline or a rational arc, or of the others (lines and conics).
inline std::vector<PlanarCase> PlanarCases(const bool splines)
{
    std::vector<PlanarCase> cases;
    for (const JsonValue& c : CasesOf("planar-curves")) {
        const std::string family = FamilyOf(c["id"].AsString());
        const bool spline = family.find("spline") != std::string::npos ||
                            family.find("rational") != std::string::npos;
        const JsonValue& expected = c["expected"];
        if (spline != splines ||
            (expected.Has("overlaps") && !expected["overlaps"].AsArray().empty())) {
            continue;
        }
        const JsonValue& curves = c["inputs"]["curves"];
        PlanarCase made{PlanarCurveOf(curves[0]), {}, PlanarCurveOf(curves[1]), {}, {}, {}};
        made.DomainA = PlanarDomainOf(curves[0], made.A);
        made.DomainB = PlanarDomainOf(curves[1], made.B);
        for (const JsonValue& point : expected["points"].AsArray()) {
            made.Roots.push_back(Point2Of(point["position"]));
            made.Within.push_back(point["within"].AsNumber());
        }
        cases.push_back(std::move(made));
    }
    return cases;
}

/// Points closer than `reach` to one kept before them are one point: OCCT
/// reports a tangency or a root on a seam more than once.
template <typename P>
std::vector<P> Distinct(const std::vector<P>& points, const double reach)
{
    std::vector<P> kept;
    for (const P& p : points) {
        if (std::ranges::none_of(kept, [&](const P& q) { return Distance(p, q) <= reach; })) {
            kept.push_back(p);
        }
    }
    return kept;
}

/// Roots matched one to one, each to an unused found point within its reach.
template <typename P>
void ScoreRoots(const std::vector<P>& roots, const std::vector<double>& within,
                const std::vector<P>& found, Score& score)
{
    const double smallest = within.empty() ? RootReach : *std::ranges::min_element(within);
    const std::vector<P> distinct = Distinct(found, smallest);
    std::vector<bool> used(distinct.size(), false);
    score.Expected += static_cast<double>(roots.size());
    for (std::size_t r = 0; r < roots.size(); ++r) {
        std::size_t best = distinct.size();
        for (std::size_t k = 0; k < distinct.size(); ++k) {
            if (!used[k] && Distance(distinct[k], roots[r]) <= within[r] &&
                (best == distinct.size() ||
                 Distance(distinct[k], roots[r]) < Distance(distinct[best], roots[r]))) {
                best = k;
            }
        }
        if (best == distinct.size()) {
            score.Missed += 1.0;
        } else {
            used[best] = true;
        }
    }
    score.Spurious += static_cast<double>(std::ranges::count(used, false));
}

// --- B2.1, lines through NURBS surfaces ---------------------------------------

struct LineCase {
    Stereon::Line Line;
    std::vector<Point3> Roots;
};

/// A B-spline surface of a corpus/nurbs-ssi/ case.
inline Stereon::BSplineSurface CorpusSpline(const std::string& id, const std::size_t index)
{
    const JsonValue c =
        StereonTest::ReadJsonFile(std::filesystem::path(STEREON_CORPUS_DIR) / (id + ".json"));
    const JsonValue& s = c["inputs"]["surfaces"][index];
    std::vector<Point3> poles;
    for (const JsonValue& p : s["poles"].AsArray()) {
        poles.push_back(Point3Of(p));
    }
    return Must(Stereon::MakeBSplineSurface(
                    static_cast<int>(s["degree_u"].AsNumber()),
                    static_cast<int>(s["degree_v"].AsNumber()), NumbersOf(s["knots_u"]),
                    NumbersOf(s["knots_v"]), std::move(poles),
                    s.Has("weights") ? NumbersOf(s["weights"]) : std::vector<double>{}),
                "corpus spline");
}

inline constexpr int LineCount = 64;

/// A unit vector across `d`.
inline Vec3 Across(const Vec3& d)
{
    return Normalized(Cross(d, std::abs(d.X) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0}));
}

/// The fractional part of k times the golden ratio's conjugate: a spread of
/// offsets in [0, 1).
inline double Spread(const int k)
{
    const double golden = 0.6180339887498949;
    return k * golden - std::floor(k * golden);
}

/// Lines through or past the unit sphere about the origin, the sphere of
/// nurbs-ssi/sphere-plane: directions on a Fibonacci spiral, and three in
/// four at a distance below 0.9 from the centre (two crossings, at
/// 3 -+ sqrt(1 - s^2) along the unit direction from the line's origin), the
/// fourth beyond 1.1 (none). Nothing is near a tangency.
inline std::vector<LineCase> LinesThroughSphere()
{
    std::vector<LineCase> lines;
    for (int k = 0; k < LineCount; ++k) {
        const double z = 1.0 - (2.0 * k + 1.0) / LineCount;
        const double phi = 2.399963229728653 * k;
        const double rho = std::sqrt(1.0 - z * z);
        const Vec3 d{rho * std::cos(phi), rho * std::sin(phi), z};
        const bool hits = k % 4 != 3;
        const double s = hits ? 0.9 * Spread(k) : 1.1 + 0.3 * Spread(k);
        const Point3 closest = Point3{} + s * Across(d);
        LineCase line{Stereon::Line{closest - 3.0 * d, d}, {}};
        if (hits) {
            const double half = std::sqrt(1.0 - s * s);
            line.Roots = {closest - half * d, closest + half * d};
        }
        lines.push_back(std::move(line));
    }
    return lines;
}

/// Lines through or past the cylinder of radius 1 about the z axis, |z| <= 3,
/// of nurbs-ssi/cylinder-plane: at most 37 degrees from the horizontal and
/// within |z| <= 1 at their closest to the axis, so both crossings lie
/// within |z| <= 1.75, far inside the patch.
inline std::vector<LineCase> LinesThroughCylinder()
{
    std::vector<LineCase> lines;
    for (int k = 0; k < LineCount; ++k) {
        const double phi = 2.399963229728653 * k;
        const double rise = 0.6 * (2.0 * Spread(k + 7) - 1.0);
        const double level = 2.0 * Spread(k + 3) - 1.0;
        const double flat = std::sqrt(1.0 - rise * rise);
        const Vec3 h{std::cos(phi), std::sin(phi), 0.0};
        const Vec3 n{-std::sin(phi), std::cos(phi), 0.0};
        const Vec3 d = flat * h + Vec3{0.0, 0.0, rise};
        const bool hits = k % 4 != 3;
        const double s = hits ? 0.9 * Spread(k) : 1.1 + 0.3 * Spread(k);
        const Point3 closest = Point3{0.0, 0.0, level} + s * n;
        LineCase line{Stereon::Line{closest - 3.0 * d, d}, {}};
        if (hits) {
            const double t = std::sqrt(1.0 - s * s) / flat;
            line.Roots = {closest - t * d, closest + t * d};
        }
        lines.push_back(std::move(line));
    }
    return lines;
}

// --- B2.2 and B2.3, surface pairs ---------------------------------------------

/// A closed form's shape: connected curves and isolated points, or an
/// overlap. Counted is false where curves meet at nodes, which ends alone
/// cannot count: the result must only have a curve.
struct Shape {
    std::size_t Curves = 0;
    std::size_t Points = 0;
    bool Overlap = false;
    bool Counted = true;
};

struct SsiCase {
    std::string Id;
    Surface A;
    Box2 BoxA;
    Surface B;
    Box2 BoxB;
    double Tolerance = 0.0;
    double Scale = 1.0;
    Point3 Anchor;
    std::vector<Shape> Readings; ///< the closed-form answers (analytic-ssi: nominal, tolerant)
    std::vector<Point3> Witnesses;
};

/// A surface of curved-ssi or nurbs-ssi (corpus/curved-ssi/README.md, Inputs),
/// its box, and its reference point.
inline std::pair<Surface, Point3> FramedSurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    if (type == "bspline") {
        std::vector<Point3> poles;
        for (const JsonValue& p : s["poles"].AsArray()) {
            poles.push_back(Point3Of(p));
        }
        const Point3 first = poles.front();
        return {
            Surface{Must(Stereon::MakeBSplineSurface(
                             static_cast<int>(s["degree_u"].AsNumber()),
                             static_cast<int>(s["degree_v"].AsNumber()), NumbersOf(s["knots_u"]),
                             NumbersOf(s["knots_v"]), std::move(poles),
                             s.Has("weights") ? NumbersOf(s["weights"]) : std::vector<double>{}),
                         "spline")},
            first};
    }
    const Vec3 x = Vec3Of(s["x_axis"]);
    const Vec3 y = Vec3Of(s["y_axis"]);
    const Point3 at = Point3Of(s[s.Has("center") ? "center" : "origin"]);
    if (type == "plane") {
        return {Surface{Stereon::Plane{at, x, y}}, at};
    }
    if (type == "sphere") {
        return {Surface{Stereon::Sphere{at, s["radius"].AsNumber(), x, y, true}}, at};
    }
    if (type == "cylinder") {
        return {Surface{Stereon::Cylinder{at, s["radius"].AsNumber(), x, y, true}}, at};
    }
    if (type == "cone") {
        return {Surface{Stereon::Cone{at, s["radius"].AsNumber(), s["semi_angle"].AsNumber(), x, y,
                                      true}},
                at};
    }
    if (type == "torus") {
        return {
            Surface{Stereon::Torus{at, s["major"].AsNumber(), s["minor"].AsNumber(), x, y, true}},
            at};
    }
    return {Surface{Stereon::Ellipsoid{at, s["a"].AsNumber(), s["b"].AsNumber(), s["c"].AsNumber(),
                                       x, y, true}},
            at};
}

inline Box2 BoxOf(const JsonValue& s)
{
    if (!s.Has("box")) {
        return Box2::Everything();
    }
    const JsonValue& box = s["box"];
    return Box2{Interval{box[0][0].AsNumber(), box[0][1].AsNumber()},
                Interval{box[1][0].AsNumber(), box[1][1].AsNumber()}};
}

/// The cases of one family of curved-ssi or nurbs-ssi, coincident surfaces left out.
inline std::vector<SsiCase> FramedCases(const std::string& suite, const std::string& family)
{
    std::vector<SsiCase> cases;
    for (const JsonValue& c : CasesOf(suite + "/" + family)) {
        const JsonValue& expected = c["expected"];
        if (expected["overlap"].AsBool()) {
            continue;
        }
        const JsonValue& surfaces = c["inputs"]["surfaces"];
        auto [a, anchor] = FramedSurfaceOf(surfaces[0]);
        Surface b = FramedSurfaceOf(surfaces[1]).first;
        SsiCase made{c["id"].AsString(),
                     std::move(a),
                     BoxOf(surfaces[0]),
                     std::move(b),
                     BoxOf(surfaces[1]),
                     c["inputs"]["tolerance"].AsNumber(),
                     c["inputs"]["scale"].AsNumber(),
                     anchor,
                     {},
                     {}};
        made.Readings.push_back(Shape{static_cast<std::size_t>(expected["curves"].AsNumber()),
                                      static_cast<std::size_t>(expected["points"].AsNumber())});
        for (const JsonValue& w : expected["witnesses"].AsArray()) {
            made.Witnesses.push_back(Point3Of(w));
        }
        cases.push_back(std::move(made));
    }
    return cases;
}

/// The shape of an analytic-ssi configuration (corpus/analytic-ssi/README.md).
inline std::optional<Shape> ShapeOf(const std::string& configuration)
{
    if (configuration == "disjoint") {
        return Shape{};
    }
    if (configuration == "coincident") {
        return Shape{0, 0, true};
    }
    if (configuration.starts_with("tangent_point")) {
        return Shape{0, 1};
    }
    if (configuration == "two_lines" || configuration == "two_circles" ||
        configuration == "two_loops") {
        return Shape{2, 0};
    }
    if (configuration == "internal_tangent" || configuration == "two_ellipses") {
        return Shape{1, 0, false, false};
    }
    if (configuration == "line" || configuration == "circle" || configuration == "ellipse" ||
        configuration == "one_loop" || configuration.starts_with("tangent_")) {
        return Shape{1, 0};
    }
    return std::nullopt;
}

/// An analytic-ssi surface: a plane keeps its normal's direction, a
/// cylinder its axis's; a sphere's frame is the coordinate frame.
inline Surface PairSurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    const auto frameAbout = [](const Vec3& direction) {
        const Vec3 n = Normalized(direction);
        const Vec3 x = Across(n);
        return std::pair{x, Cross(n, x)};
    };
    if (type == "plane") {
        const auto [x, y] = frameAbout(Vec3Of(s["normal"]));
        return Surface{Stereon::Plane{Point3Of(s["point"]), x, y}};
    }
    if (type == "cylinder") {
        const auto [x, y] = frameAbout(Vec3Of(s["axis"]));
        return Surface{Stereon::Cylinder{Point3Of(s["point"]), s["radius"].AsNumber(), x, y, true}};
    }
    return Surface{Stereon::Sphere{Point3Of(s["centre"]), s["radius"].AsNumber(),
                                   Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
}

inline Point3 PairAnchorOf(const JsonValue& s)
{
    return Point3Of(s[s.Has("centre") ? "centre" : "point"]);
}

/// The Tier A pair cases of analytic-ssi whose two surfaces are of the given
/// types, in alphabetical order ("cylinder-sphere", "plane-plane", ...), the
/// cases the kernel must refuse (expected.kernel) left out.
inline std::vector<SsiCase> PairCases(const std::string& pair)
{
    std::vector<SsiCase> cases;
    for (const JsonValue& c : CasesOf("analytic-ssi")) {
        const JsonValue& expected = c["expected"];
        if (c["operation"]["name"].AsString() != "intersect" || expected.Has("kernel") ||
            expected["nominal"]["tier"].AsString() != "A") {
            continue;
        }
        const JsonValue& s = c["inputs"]["surfaces"];
        std::string first = s[0]["type"].AsString();
        std::string second = s[1]["type"].AsString();
        if (second < first) {
            std::swap(first, second);
        }
        if (first + "-" + second != pair) {
            continue;
        }
        const std::string nominal = expected["nominal"]["configuration"].AsString();
        const std::string tolerant = expected["tolerant"].Has("configuration")
                                         ? expected["tolerant"]["configuration"].AsString()
                                         : nominal;
        SsiCase made{c["id"].AsString(),
                     PairSurfaceOf(s[0]),
                     Box2::Everything(),
                     PairSurfaceOf(s[1]),
                     Box2::Everything(),
                     std::max(s[0]["tolerance"].AsNumber(), s[1]["tolerance"].AsNumber()),
                     1.0,
                     PairAnchorOf(s[0]),
                     {},
                     {}};
        for (const std::string& configuration : {nominal, tolerant}) {
            if (const std::optional<Shape> shape = ShapeOf(configuration)) {
                made.Readings.push_back(*shape);
            }
        }
        cases.push_back(std::move(made));
    }
    return cases;
}

/// A library's answer to a surface pair: sampled curves, points, an overlap,
/// or no answer.
struct Found {
    bool Resolved = false;
    bool Overlap = false;
    std::vector<std::vector<Point3>> Curves;
    std::vector<Point3> Points;
};

/// The parameters at which a branch is sampled: its range, or a window about
/// the anchor along an unbounded line, or [-10, 10] of an unbounded conic.
inline std::vector<double> BranchParameters(const Stereon::IntersectionBranch& branch,
                                            const Point3& anchor, const double scale)
{
    double lower = branch.Range.GetLower();
    double upper = branch.Range.GetUpper();
    if (const auto* line = branch.Curve3d.Target<Stereon::Line>(); line != nullptr) {
        const double speed = Length(line->Direction);
        const double centre = Dot(anchor - line->Origin, line->Direction) / (speed * speed);
        const double half = LineWindow * scale / speed;
        lower = std::max(lower, centre - half);
        upper = std::min(upper, centre + half);
    } else {
        lower = std::isinf(lower) ? -10.0 : lower;
        upper = std::isinf(upper) ? 10.0 : upper;
    }
    std::vector<double> t;
    for (int k = 0; k <= CurveSamples; ++k) {
        t.push_back(lower + (upper - lower) * k / CurveSamples);
    }
    return t;
}

inline Found FoundOf(const Stereon::Result<Stereon::SurfaceSurfaceIntersection>& result,
                     const Point3& anchor, const double scale)
{
    Found found;
    if (!result) {
        return found;
    }
    found.Resolved = true;
    found.Overlap = !result->Overlaps.empty();
    for (const Stereon::IntersectionPoint& p : result->Points) {
        found.Points.push_back(p.Position);
    }
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        std::vector<Point3> samples;
        for (const double t : BranchParameters(branch, anchor, scale)) {
            samples.push_back(Eval(branch.Curve3d, t));
        }
        found.Curves.push_back(std::move(samples));
    }
    return found;
}

/// The connected curves the sampled curves form, ends within `join` joined.
inline std::size_t ConnectedCurves(const std::vector<std::vector<Point3>>& curves,
                                   const double join)
{
    std::vector<Point3> clusters;
    std::vector<std::size_t> parent;
    const auto clusterOf = [&](const Point3& p) {
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (Distance(clusters[i], p) <= join) {
                return i;
            }
        }
        clusters.push_back(p);
        parent.push_back(parent.size());
        return clusters.size() - 1;
    };
    const auto find = [&](std::size_t i) {
        while (parent[i] != i) {
            i = parent[i] = parent[parent[i]];
        }
        return i;
    };
    for (const std::vector<Point3>& curve : curves) {
        const std::size_t a = clusterOf(curve.front());
        const std::size_t b = clusterOf(curve.back());
        parent[find(a)] = find(b);
    }
    std::size_t roots = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i) {
        if (find(i) == i) {
            ++roots;
        }
    }
    return roots;
}

inline double DistanceToSegment(const Point3& p, const Point3& a, const Point3& b)
{
    const Vec3 along = b - a;
    const double square = Dot(along, along);
    if (square == 0.0) {
        return Distance(p, a);
    }
    const double s = std::clamp(Dot(p - a, along) / square, 0.0, 1.0);
    return Distance(p, a + s * along);
}

inline double DistanceToFound(const Found& found, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Point3& q : found.Points) {
        nearest = std::min(nearest, Distance(p, q));
    }
    for (const std::vector<Point3>& curve : found.Curves) {
        for (std::size_t k = 1; k < curve.size(); ++k) {
            nearest = std::min(nearest, DistanceToSegment(p, curve[k - 1], curve[k]));
        }
    }
    return nearest;
}

inline double DistanceToSurface(const Surface& s, const Point3& p)
{
    return Distance(p, Eval(s, Project(s, p)));
}

/// The largest distance of a sampled point from either surface.
inline double DeviationOf(const SsiCase& c, const Found& found)
{
    double worst = 0.0;
    const auto measure = [&](const Point3& p) {
        worst = std::max({worst, DistanceToSurface(c.A, p), DistanceToSurface(c.B, p)});
    };
    std::ranges::for_each(found.Points, measure);
    for (const std::vector<Point3>& curve : found.Curves) {
        std::ranges::for_each(curve, measure);
    }
    return worst;
}

/// Whether the answer has the shape: the overlap, or the curves and points.
inline bool HasShape(const Found& found, const std::size_t curves, const Shape& want)
{
    if (want.Overlap || found.Overlap) {
        return want.Overlap == found.Overlap;
    }
    const bool counted = want.Counted ? curves == want.Curves : curves > 0;
    return counted && found.Points.size() == want.Points;
}

inline void ScoreCase(const SsiCase& c, const Found& found, Score& score)
{
    score.Cases += 1.0;
    const Shape& want = c.Readings.front();
    const double expected = static_cast<double>(want.Curves + want.Points);
    score.Expected += expected;
    if (!found.Resolved) {
        score.Unresolved += 1.0;
        return;
    }
    const std::size_t curves = ConnectedCurves(found.Curves, JoinReach * c.Scale);
    if (std::ranges::none_of(c.Readings,
                             [&](const Shape& s) { return HasShape(found, curves, s); })) {
        score.Inexact += 1.0;
        const double got = static_cast<double>(curves + found.Points.size());
        score.Missed += std::max(0.0, expected - got);
        score.Spurious += std::max(0.0, got - expected);
    }
    for (const Point3& w : c.Witnesses) {
        score.Missed += DistanceToFound(found, w) > WitnessReach * c.Scale ? 1.0 : 0.0;
    }
    score.MaxDeviation = std::max(score.MaxDeviation, DeviationOf(c, found));
}

/// Stereon's answers to a row's cases, once.
inline std::vector<Found> StereonAnswers(const std::vector<SsiCase>& cases)
{
    std::vector<Found> answers;
    for (const SsiCase& c : cases) {
        answers.push_back(
            FoundOf(Intersect(c.A, c.BoxA, c.B, c.BoxB, c.Tolerance, Ctx()), c.Anchor, c.Scale));
    }
    return answers;
}

inline Score ScoreSsi(const std::vector<SsiCase>& cases, const std::vector<Found>& answers)
{
    Score score;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        ScoreCase(cases[i], answers[i], score);
    }
    return score;
}

/// The row's counters (docs/benchmarks.md: branches found; missed and
/// spurious; deviation; exactness).
template <typename State>
void Report(State& state, const Score& score)
{
    state.counters["cases"] = score.Cases;
    state.counters["expected"] = score.Expected;
    state.counters["missed"] = score.Missed;
    state.counters["spurious"] = score.Spurious;
    state.counters["unresolved"] = score.Unresolved;
    state.counters["inexact"] = score.Inexact;
    state.counters["max_deviation"] = score.MaxDeviation;
}

/// The families of the B2.2 rows, and the pairs of the B2.3 rows.
inline const std::vector<std::string> NurbsFamilies{"bump-plane",     "cylinder-cylinder",
                                                    "cylinder-plane", "plane-plane",
                                                    "saddle-plane",   "sphere-plane"};
inline const std::vector<std::string> CurvedFamilies{
    "cone-cone",          "cone-ellipsoid",  "cone-torus",          "cylinder-cone",
    "cylinder-ellipsoid", "cylinder-torus",  "ellipsoid-ellipsoid", "plane-cone",
    "plane-ellipsoid",    "plane-torus",     "sphere-cone",         "sphere-ellipsoid",
    "sphere-torus",       "torus-ellipsoid", "torus-torus"};
inline const std::vector<std::string> AnalyticPairs{"cylinder-cylinder", "cylinder-plane",
                                                    "cylinder-sphere",   "plane-plane",
                                                    "plane-sphere",      "sphere-sphere"};

} // namespace StereonBench::B2
