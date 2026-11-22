// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The SSI corpus against OpenCASCADE's general surface intersection (ADR-0010,
// Verification: completeness checked against ground truth and against OCCT):
// every pair case of corpus/analytic-ssi/, corpus/curved-ssi/ and
// corpus/nurbs-ssi/ is intersected by GeomAPI_IntSS's engine, GeomInt_IntSS
// with GeomAPI_IntSS's arguments (Approx on, no pcurve approximation), which
// also returns the isolated points GeomAPI_IntSS does not expose, at the
// case's tolerance.
//
// The closed form is the oracle, not OCCT. Each case is decided in two steps:
//
// 1. OCCT against the closed form: the number of connected curves its lines
//    form (ends within JoinReach joined), the number of isolated points,
//    every witness of the case on OCCT's result, and every sample of OCCT's
//    curves on both surfaces (WitnessReach). Where OCCT disagrees, the closed
//    form wins: the case is an OCCT discrepancy, recorded with the reason,
//    never a kernel failure. analytic-ssi gives a configuration instead of
//    counts, and both its exact (nominal) and its tolerant reading count as
//    a closed-form answer.
// 2. Where OCCT agrees with the closed form, Stereon against OCCT, both ways:
//    every sample of Stereon's branches and every Stereon point lies on
//    OCCT's result, and every sample of OCCT's on Stereon's. A difference
//    there is a kernel failure: the closed form and OCCT agree on it. A
//    refusal (IntersectionFailed and the like) is allowed where the case
//    says may_fail_closed, as in tests/intersect/SurfaceCorpusTests.cpp.
//
// Not compared: coincident surfaces, which GeomAPI_IntSS answers with no
// curve (it has no overlap outcome), analytic-ssi's Tier B pairs other than
// skew cylinders, which have no closed-form configuration, and arrangements.
// Each test writes its tally and every discrepancy as Markdown to the file
// named by STEREON_OCCT_CORPUS_REPORT, when it is set (the OCCT workflow puts
// it in the job summary).
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008).

#include "Json.hpp"
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomEval_EllipsoidSurface.hxx>
#include <GeomInt_IntSS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Standard_Failure.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "../intersect/Checks.hpp"
#include "OcctConversion.hpp"

using Stereon::Box2;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

using StereonTest::ToOcct;

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR);

/// How near a witness or a sample must lie to a result or a surface, per unit
/// of the case's scale: tests/intersect/SurfaceCorpusTests.cpp's reach, well
/// above the sagitta of the samples and OCCT's approximation of a walking
/// line, far below any feature of a case.
constexpr double WitnessReach = 1e-3;
/// Ends of OCCT's curves nearer than this, per unit of scale, are one point:
/// OCCT cuts a walking line where its patches or its approximation end.
constexpr double JoinReach = 1e-4;
/// Samples per OCCT curve, and segments per Stereon branch when a point is
/// looked for on it.
constexpr int CurveSamples = 64;
constexpr std::size_t BranchSamples = 64;
constexpr std::size_t BranchSegments = 2048;
/// An unbounded line is sampled this far either side of the case's anchor,
/// per unit of scale; an unbounded conic (a parabola, a hyperbola) over this
/// stretch of its parameter about its vertex.
constexpr double LineWindow = 10.0;
constexpr double ConicWindow = 3.0;
/// A parameter beyond this is unbounded: OCCT trims the lines of an
/// analytic intersection at about 1.5e9, far past the size box (ADR-0001).
constexpr double Unbounded = 1e7;

Point3 PointOf(const JsonValue& v)
{
    return Point3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Vec3 VectorOf(const JsonValue& v)
{
    return Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Point3 PointOf(const gp_Pnt& p)
{
    return Point3{p.X(), p.Y(), p.Z()};
}

gp_Dir ToDir(const Vec3& v)
{
    return gp_Dir(v.X, v.Y, v.Z);
}

std::vector<double> NumbersOf(const JsonValue& v)
{
    std::vector<double> numbers;
    for (const JsonValue& n : v.AsArray()) {
        numbers.push_back(n.AsNumber());
    }
    return numbers;
}

// --- The cases' surfaces, in both kernels -------------------------------------

/// The two kernels' versions of one surface of a case, over the case's box.
struct CaseSurface {
    Surface Stereon;
    Box2 Box = Box2::Everything();
    occ::handle<Geom_Surface> Occt;
    /// A point of the surface near the case: where unbounded curves are sampled.
    Point3 Anchor;
};

Stereon::BSplineSurface SplineOf(const JsonValue& s)
{
    std::vector<Point3> poles;
    for (const JsonValue& p : s["poles"].AsArray()) {
        poles.push_back(PointOf(p));
    }
    auto made = Stereon::MakeBSplineSurface(
        static_cast<int>(s["degree_u"].AsNumber()), static_cast<int>(s["degree_v"].AsNumber()),
        NumbersOf(s["knots_u"]), NumbersOf(s["knots_v"]), std::move(poles),
        s.Has("weights") ? NumbersOf(s["weights"]) : std::vector<double>{});
    REQUIRE(made.has_value());
    return *std::move(made);
}

/// A surface of curved-ssi or nurbs-ssi (corpus/curved-ssi/README.md, Inputs).
/// Stereon's frames and parametrisations of the analytic surfaces are
/// OCCT's (STEP's), so a box carries over as it is.
CaseSurface FramedSurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    if (type == "bspline") {
        const Stereon::BSplineSurface spline = SplineOf(s);
        return {Surface{spline}, Box2::Everything(), ToOcct(spline), spline.Poles.front()};
    }
    const Vec3 x = VectorOf(s["x_axis"]);
    const Vec3 y = VectorOf(s["y_axis"]);
    const Point3 at = PointOf(s[s.Has("center") ? "center" : "origin"]);
    const gp_Ax3 frame(ToOcct(at), ToDir(Cross(x, y)), ToDir(x));
    CaseSurface made{Surface{Stereon::Plane{at, x, y}}, Box2::Everything(), nullptr, at};
    if (type == "plane") {
        made.Occt = new Geom_Plane(frame);
    } else if (type == "sphere") {
        const double r = s["radius"].AsNumber();
        made.Stereon = Surface{Stereon::Sphere{at, r, x, y, true}};
        made.Occt = new Geom_SphericalSurface(frame, r);
    } else if (type == "cylinder") {
        const double r = s["radius"].AsNumber();
        made.Stereon = Surface{Stereon::Cylinder{at, r, x, y, true}};
        made.Occt = new Geom_CylindricalSurface(frame, r);
    } else if (type == "cone") {
        const double r = s["radius"].AsNumber();
        const double angle = s["semi_angle"].AsNumber();
        made.Stereon = Surface{Stereon::Cone{at, r, angle, x, y, true}};
        made.Occt = new Geom_ConicalSurface(frame, angle, r);
    } else if (type == "torus") {
        const double major = s["major"].AsNumber();
        const double minor = s["minor"].AsNumber();
        made.Stereon = Surface{Stereon::Torus{at, major, minor, x, y, true}};
        made.Occt = new Geom_ToroidalSurface(frame, major, minor);
    } else {
        REQUIRE(type == "ellipsoid");
        const double a = s["a"].AsNumber();
        const double b = s["b"].AsNumber();
        const double c = s["c"].AsNumber();
        made.Stereon = Surface{Stereon::Ellipsoid{at, a, b, c, x, y, true}};
        made.Occt = new GeomEval_EllipsoidSurface(frame, a, b, c);
    }
    if (s.Has("box")) {
        const JsonValue& box = s["box"];
        made.Box = Box2{Interval{box[0][0].AsNumber(), box[0][1].AsNumber()},
                        Interval{box[1][0].AsNumber(), box[1][1].AsNumber()}};
        made.Occt = new Geom_RectangularTrimmedSurface(made.Occt, made.Box.X.GetLower(),
                                                       made.Box.X.GetUpper(), made.Box.Y.GetLower(),
                                                       made.Box.Y.GetUpper());
    }
    return made;
}

/// A surface of analytic-ssi (corpus/analytic-ssi/README.md, Inputs): a
/// plane keeps its normal's direction, a cylinder its axis's; a sphere's
/// frame is the coordinate frame. OCCT picks its own X direction from the
/// normal or the axis, which changes the parametrisation, not the surface.
CaseSurface PairSurfaceOf(const JsonValue& s)
{
    const std::string& type = s["type"].AsString();
    if (type == "plane") {
        const Point3 at = PointOf(s["point"]);
        const Vec3 n = VectorOf(s["normal"]);
        return {Surface{IntersectTest::PlaneWithNormal(at, n)}, Box2::Everything(),
                new Geom_Plane(gp_Pln(ToOcct(at), ToDir(n))), at};
    }
    if (type == "cylinder") {
        const Point3 at = PointOf(s["point"]);
        const Vec3 axis = VectorOf(s["axis"]);
        const double r = s["radius"].AsNumber();
        return {Surface{IntersectTest::CylinderAbout(at, axis, r)}, Box2::Everything(),
                new Geom_CylindricalSurface(gp_Ax3(ToOcct(at), ToDir(axis)), r), at};
    }
    const Point3 centre = PointOf(s["centre"]);
    const double r = s["radius"].AsNumber();
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    return {Surface{Stereon::Sphere{centre, r, ex, ey, true}}, Box2::Everything(),
            new Geom_SphericalSurface(gp_Ax3(ToOcct(centre), ToDir(Cross(ex, ey)), ToDir(ex)), r),
            centre};
}

// --- OCCT's result ------------------------------------------------------------

struct OcctCurve {
    occ::handle<Geom_Curve> Curve;
    double First = 0.0;
    double Last = 0.0;
};

struct OcctResult {
    bool Done = false;
    std::string Failure;
    std::vector<OcctCurve> Curves;
    std::vector<Point3> Points;
};

OcctResult IntersectWithOcct(const occ::handle<Geom_Surface>& a, const occ::handle<Geom_Surface>& b,
                             const double tolerance)
{
    OcctResult result;
    try {
        // GeomAPI_IntSS(a, b, tolerance) is this call (GeomAPI_IntSS.lxx).
        const GeomInt_IntSS occt(a, b, tolerance, true);
        if (!occt.IsDone()) {
            result.Failure = "IsDone() is false";
            return result;
        }
        for (int i = 1; i <= occt.NbLines(); ++i) {
            const occ::handle<Geom_Curve>& c = occt.Line(i);
            result.Curves.push_back({c, c->FirstParameter(), c->LastParameter()});
        }
        for (int i = 1; i <= occt.NbPoints(); ++i) {
            result.Points.push_back(PointOf(occt.Point(i)));
        }
        result.Done = true;
    } catch (const std::exception& e) {
        result.Failure = std::string("exception: ") + e.what();
    }
    return result;
}

/// The line under a curve, or null.
const Geom_Line* LineOf(const occ::handle<Geom_Curve>& curve)
{
    const Geom_Curve* basis = curve.get();
    if (const auto* trimmed = dynamic_cast<const Geom_TrimmedCurve*>(basis)) {
        basis = trimmed->BasisCurve().get();
    }
    return dynamic_cast<const Geom_Line*>(basis);
}

bool IsUnbounded(const OcctCurve& c)
{
    return std::abs(c.First) >= Unbounded || std::abs(c.Last) >= Unbounded;
}

/// The stretch of a curve's parameter that is sampled: the whole range of a
/// bounded curve, a window about the anchor on an unbounded line, a window
/// about the vertex on an unbounded conic.
std::pair<double, double> SampledRange(const OcctCurve& c, const Point3& anchor, const double scale)
{
    if (!IsUnbounded(c)) {
        return {c.First, c.Last};
    }
    double centre = 0.0;
    double half = ConicWindow;
    if (const Geom_Line* line = LineOf(c.Curve)) {
        const gp_Lin lin = line->Lin();
        centre = gp_Vec(lin.Location(), ToOcct(anchor)).Dot(gp_Vec(lin.Direction()));
        half = LineWindow * scale;
    }
    return {std::max(c.First, centre - half), std::min(c.Last, centre + half)};
}

std::vector<Point3> SamplesOf(const OcctCurve& c, const Point3& anchor, const double scale)
{
    const auto [first, last] = SampledRange(c, anchor, scale);
    std::vector<Point3> samples;
    for (int k = 0; k <= CurveSamples; ++k) {
        samples.push_back(PointOf(c.Curve->Value(first + (last - first) * k / CurveSamples)));
    }
    return samples;
}

/// The distance from p to a curve of OCCT's: an unbounded line taken whole;
/// otherwise the nearest foot of a perpendicular over the sampled range or
/// the nearer end, or the nearest of dense samples if the projection fails.
double DistanceTo(const OcctCurve& c, const Point3& p, const Point3& anchor, const double scale)
{
    const gp_Pnt q = ToOcct(p);
    if (const Geom_Line* line = LineOf(c.Curve); line != nullptr && IsUnbounded(c)) {
        return line->Lin().Distance(q);
    }
    const auto [first, last] = SampledRange(c, anchor, scale);
    double nearest = std::min(c.Curve->Value(first).Distance(q), c.Curve->Value(last).Distance(q));
    try {
        GeomAPI_ProjectPointOnCurve projection(q, c.Curve, first, last);
        if (projection.NbPoints() > 0) {
            nearest = std::min(nearest, projection.LowerDistance());
        }
    } catch (const std::exception&) {
        for (int k = 0; k <= 20000; ++k) {
            nearest =
                std::min(nearest, c.Curve->Value(first + (last - first) * k / 20000.0).Distance(q));
        }
    }
    return nearest;
}

double DistanceTo(const OcctResult& r, const Point3& p, const Point3& anchor, const double scale)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Point3& point : r.Points) {
        nearest = std::min(nearest, Distance(p, point));
    }
    for (const OcctCurve& c : r.Curves) {
        nearest = std::min(nearest, DistanceTo(c, p, anchor, scale));
    }
    return nearest;
}

/// The connected curves OCCT's lines form, ends within `join` joined; an
/// unbounded curve's ends at infinity join nothing.
std::size_t ConnectedCurves(const OcctResult& r, const double join)
{
    std::vector<Point3> clusters;
    std::vector<std::size_t> ends;
    const auto clusterOf = [&](const Point3& p) {
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (Distance(clusters[i], p) <= join) {
                return i;
            }
        }
        clusters.push_back(p);
        return clusters.size() - 1;
    };
    for (const OcctCurve& c : r.Curves) {
        if (IsUnbounded(c)) {
            clusters.push_back(Point3{});
            ends.push_back(clusters.size() - 1);
            ends.push_back(clusters.size() - 1);
            continue;
        }
        ends.push_back(clusterOf(PointOf(c.Curve->Value(c.First))));
        ends.push_back(clusterOf(PointOf(c.Curve->Value(c.Last))));
    }
    std::vector<std::size_t> parent(clusters.size());
    for (std::size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
    const auto find = [&](std::size_t i) {
        while (parent[i] != i) {
            i = parent[i] = parent[parent[i]];
        }
        return i;
    };
    for (std::size_t e = 0; e + 1 < ends.size(); e += 2) {
        parent[find(ends[e])] = find(ends[e + 1]);
    }
    std::size_t roots = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i) {
        if (find(i) == i) {
            ++roots;
        }
    }
    return roots;
}

double DistanceToSurface(const Surface& s, const Point3& p)
{
    return Distance(p, Eval(s, Project(s, p)));
}

// --- Stereon's result ---------------------------------------------------------

/// A branch's parameters: its range, or a window about the anchor along an
/// unbounded line (lines are the only unbounded branches of the corpus but
/// for a parabola or hyperbola of a cone, sampled over [-10, 10]).
std::vector<double> ParametersOf(const Stereon::IntersectionBranch& branch, const Point3& anchor,
                                 const double scale, const std::size_t count)
{
    const auto* line = branch.Curve3d.Target<Stereon::Line>();
    if (line == nullptr || branch.Range.IsBounded()) {
        return IntersectTest::SamplesOf(branch.Range, count);
    }
    const double speed = Length(line->Direction);
    const double centre = Dot(anchor - line->Origin, line->Direction) / (speed * speed);
    const double half = LineWindow * scale / speed;
    const double lower = std::isinf(branch.Range.GetLower())
                             ? centre - half
                             : std::max(branch.Range.GetLower(), centre - half);
    const double upper = std::isinf(branch.Range.GetUpper())
                             ? centre + half
                             : std::min(branch.Range.GetUpper(), centre + half);
    return IntersectTest::SamplesOf(Interval{lower, upper}, count);
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

/// The distance from p to Stereon's result: an unbounded line taken whole,
/// any other branch as a polyline of BranchSegments segments.
double DistanceTo(const Stereon::SurfaceSurfaceIntersection& r, const Point3& p,
                  const Point3& anchor, const double scale)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const Stereon::IntersectionPoint& point : r.Points) {
        nearest = std::min(nearest, Distance(p, point.Position));
    }
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        if (const auto* line = branch.Curve3d.Target<Stereon::Line>();
            line != nullptr && !branch.Range.IsBounded()) {
            const Vec3 u = Normalized(line->Direction);
            nearest = std::min(nearest, Length(Cross(p - line->Origin, u)));
            continue;
        }
        const std::vector<double> t = ParametersOf(branch, anchor, scale, BranchSegments);
        Point3 previous = Eval(branch.Curve3d, t.front());
        for (std::size_t k = 1; k < t.size(); ++k) {
            const Point3 next = Eval(branch.Curve3d, t[k]);
            nearest = std::min(nearest, DistanceToSegment(p, previous, next));
            previous = next;
        }
    }
    return nearest;
}

std::string Format(const Point3& p)
{
    std::ostringstream text;
    text.precision(10);
    text << '(' << p.X << ", " << p.Y << ", " << p.Z << ')';
    return text.str();
}

std::string Format(const double d)
{
    std::ostringstream text;
    text.precision(3);
    text << d;
    return text.str();
}

/// Both ways within `reach`: every sample of Stereon's result on OCCT's, every
/// sample of OCCT's on Stereon's. The first point that is not, or "".
std::string Difference(const Stereon::SurfaceSurfaceIntersection& ours, const OcctResult& theirs,
                       const Point3& anchor, const double scale, const double reach)
{
    for (const Stereon::IntersectionPoint& point : ours.Points) {
        const double d = DistanceTo(theirs, point.Position, anchor, scale);
        if (d > reach) {
            return "Stereon's point " + Format(point.Position) + " is " + Format(d) +
                   " from OCCT's result";
        }
    }
    for (const Stereon::IntersectionBranch& branch : ours.Branches) {
        for (const double t : ParametersOf(branch, anchor, scale, BranchSamples)) {
            const Point3 p = Eval(branch.Curve3d, t);
            const double d = DistanceTo(theirs, p, anchor, scale);
            if (d > reach) {
                return "Stereon's branch point " + Format(p) + " is " + Format(d) +
                       " from OCCT's result";
            }
        }
    }
    for (const Point3& point : theirs.Points) {
        const double d = DistanceTo(ours, point, anchor, scale);
        if (d > reach) {
            return "OCCT's point " + Format(point) + " is " + Format(d) + " from Stereon's result";
        }
    }
    for (const OcctCurve& c : theirs.Curves) {
        for (const Point3& p : SamplesOf(c, anchor, scale)) {
            const double d = DistanceTo(ours, p, anchor, scale);
            if (d > reach) {
                return "OCCT's curve point " + Format(p) + " is " + Format(d) +
                       " from Stereon's result";
            }
        }
    }
    return "";
}

bool FailsClosed(const ErrorCode code)
{
    return code == ErrorCode::IntersectionFailed || code == ErrorCode::CertificationFailed ||
           code == ErrorCode::NotSupported || code == ErrorCode::NotImplemented;
}

// --- Judging OCCT against the closed form -------------------------------------

/// The closed form's shape: connected curves and isolated points, or an
/// overlap. `Counted` is false where the curves meet at nodes (two loops
/// touching, two ellipses crossing), which ends alone cannot count: there
/// the result must only have a curve.
struct Shape {
    std::size_t Curves = 0;
    std::size_t Points = 0;
    bool Overlap = false;
    bool Counted = true;
};

/// Why OCCT's result is not the shape, or "".
std::string Mismatch(const OcctResult& r, const Shape& want, const CaseSurface& a,
                     const CaseSurface& b, const double scale)
{
    if (want.Overlap) {
        return r.Curves.empty() && r.Points.empty() ? ""
                                                    : "curves or points on coincident surfaces";
    }
    const double reach = WitnessReach * scale;
    std::string why;
    const std::size_t curves = ConnectedCurves(r, JoinReach * scale);
    if (want.Counted ? curves != want.Curves : curves == 0) {
        why += std::to_string(curves) + " connected curves (" + std::to_string(r.Curves.size()) +
               " lines) for " + std::to_string(want.Curves) + "; ";
    }
    if (r.Points.size() != want.Points) {
        why +=
            std::to_string(r.Points.size()) + " points for " + std::to_string(want.Points) + "; ";
    }
    double worst = 0.0;
    Point3 at;
    for (const OcctCurve& c : r.Curves) {
        for (const Point3& p : SamplesOf(c, a.Anchor, scale)) {
            const double d =
                std::max(DistanceToSurface(a.Stereon, p), DistanceToSurface(b.Stereon, p));
            if (d > worst) {
                worst = d;
                at = p;
            }
        }
    }
    for (const Point3& p : r.Points) {
        const double d = std::max(DistanceToSurface(a.Stereon, p), DistanceToSurface(b.Stereon, p));
        if (d > worst) {
            worst = d;
            at = p;
        }
    }
    if (worst > reach) {
        why +=
            "a point of its result, " + Format(at) + ", is " + Format(worst) + " from a surface; ";
    }
    return why;
}

/// The witnesses of a case that OCCT's result misses, as a reason, or "".
std::string MissedWitnesses(const OcctResult& r, const JsonValue& witnesses, const Point3& anchor,
                            const double scale)
{
    std::string why;
    for (const JsonValue& witness : witnesses.AsArray()) {
        const Point3 p = PointOf(witness);
        const double d = DistanceTo(r, p, anchor, scale);
        if (d > WitnessReach * scale) {
            why += "witness " + Format(p) + " is " + Format(d) + " from it; ";
        }
    }
    return why;
}

// --- The tally ----------------------------------------------------------------

/// How one case came out.
enum class Outcome : std::uint8_t {
    Agree,           ///< OCCT agrees with the closed form, Stereon with OCCT
    StereonRefused,  ///< OCCT agrees with the closed form; Stereon fails closed, as allowed
    OtherReading,    ///< OCCT agrees with the closed-form reading Stereon is not held to
    OcctDiscrepancy, ///< OCCT disagrees with the closed form
    OcctFailed,      ///< OCCT reports no result (IsDone false, or an exception)
    NotCompared,     ///< coincident surfaces, or a Stereon refusal the case requires
    KernelFailure,   ///< the closed form and OCCT agree, Stereon does not
};

constexpr const char* NameOf(const Outcome o)
{
    switch (o) {
    case Outcome::Agree:
        return "agree";
    case Outcome::StereonRefused:
        return "Stereon fails closed";
    case Outcome::OtherReading:
        return "OCCT takes the other reading";
    case Outcome::OcctDiscrepancy:
        return "OCCT discrepancy";
    case Outcome::OcctFailed:
        return "OCCT failed";
    case Outcome::NotCompared:
        return "not compared";
    case Outcome::KernelFailure:
        return "kernel failure";
    }
    return "?";
}

constexpr Outcome Outcomes[] = {
    Outcome::Agree,      Outcome::StereonRefused, Outcome::OtherReading, Outcome::OcctDiscrepancy,
    Outcome::OcctFailed, Outcome::NotCompared,    Outcome::KernelFailure};

struct Tally {
    std::map<std::string, std::map<Outcome, std::size_t>> ByFamily;
    std::vector<std::pair<std::string, std::string>> Notes; // (case, outcome: reason)

    void Add(const std::string& id, const Outcome outcome, const std::string& reason = "")
    {
        const std::string family = id.substr(0, id.rfind('/'));
        ++ByFamily[family][outcome];
        if (outcome != Outcome::Agree && outcome != Outcome::NotCompared) {
            Notes.emplace_back(id, std::string(NameOf(outcome)) + ": " + reason);
        }
    }

    [[nodiscard]] std::size_t Count(const Outcome outcome) const
    {
        std::size_t n = 0;
        for (const auto& [family, counts] : ByFamily) {
            const auto found = counts.find(outcome);
            n += found == counts.end() ? 0 : found->second;
        }
        return n;
    }

    [[nodiscard]] std::size_t Total() const
    {
        std::size_t n = 0;
        for (const Outcome o : Outcomes) {
            n += Count(o);
        }
        return n;
    }

    /// The tally as Markdown, appended to STEREON_OCCT_CORPUS_REPORT if set.
    void Write(const std::string& title) const
    {
        const char* path = std::getenv("STEREON_OCCT_CORPUS_REPORT");
        if (path == nullptr || *path == '\0') {
            return;
        }
        std::ofstream out(path, std::ios::app);
        out << "### " << title << "\n\n| Family |";
        for (const Outcome o : Outcomes) {
            out << ' ' << NameOf(o) << " |";
        }
        out << "\n| --- |";
        for ([[maybe_unused]] const Outcome o : Outcomes) {
            out << " ---: |";
        }
        out << '\n';
        for (const auto& [family, counts] : ByFamily) {
            out << "| " << family << " |";
            for (const Outcome o : Outcomes) {
                const auto found = counts.find(o);
                out << ' ' << (found == counts.end() ? 0 : found->second) << " |";
            }
            out << '\n';
        }
        out << "| **all** |";
        for (const Outcome o : Outcomes) {
            out << " **" << Count(o) << "** |";
        }
        out << "\n\n";
        if (!Notes.empty()) {
            out << "<details><summary>" << Notes.size()
                << " cases other than agreement</summary>\n\n| Case | Outcome |\n| --- | --- |\n";
            for (const auto& [id, note] : Notes) {
                std::string cell = note;
                std::ranges::replace(cell, '|', '/');
                out << "| `" << id << "` | " << cell << " |\n";
            }
            out << "\n</details>\n\n";
        }
    }
};

std::vector<std::filesystem::path> CaseFiles(const std::string& suite)
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(corpus / suite)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

// --- curved-ssi and nurbs-ssi -------------------------------------------------

void CompareFramedCase(const JsonValue& c, Tally& tally)
{
    const std::string& id = c["id"].AsString();
    INFO(id);
    const JsonValue& surfaces = c["inputs"]["surfaces"];
    const CaseSurface a = FramedSurfaceOf(surfaces[0]);
    const CaseSurface b = FramedSurfaceOf(surfaces[1]);
    const double tolerance = c["inputs"]["tolerance"].AsNumber();
    const double scale = c["inputs"]["scale"].AsNumber();
    const JsonValue& expected = c["expected"];

    const OcctResult occt = IntersectWithOcct(a.Occt, b.Occt, tolerance);
    if (!occt.Done) {
        tally.Add(id, Outcome::OcctFailed, occt.Failure);
        return;
    }
    const Shape want{static_cast<std::size_t>(expected["curves"].AsNumber()),
                     static_cast<std::size_t>(expected["points"].AsNumber()),
                     expected["overlap"].AsBool()};
    if (want.Overlap) {
        const std::string why = Mismatch(occt, want, a, b, scale);
        tally.Add(id, why.empty() ? Outcome::NotCompared : Outcome::OcctDiscrepancy, why);
        return;
    }
    const std::string why = Mismatch(occt, want, a, b, scale) +
                            MissedWitnesses(occt, expected["witnesses"], a.Anchor, scale);
    if (!why.empty()) {
        tally.Add(id, Outcome::OcctDiscrepancy, why);
        return;
    }

    const auto ours = Intersect(a.Stereon, a.Box, b.Stereon, b.Box, tolerance, IntersectTest::ctx);
    if (!ours) {
        INFO(IntersectTest::ErrorText(ours));
        CHECK(FailsClosed(ours.error().GetCode()));
        CHECK(expected.Has("may_fail_closed"));
        tally.Add(id, Outcome::StereonRefused, std::string(EnumName(ours.error().GetCode())));
        return;
    }
    const std::string difference = Difference(*ours, occt, a.Anchor, scale, WitnessReach * scale);
    INFO(difference);
    CHECK(difference.empty());
    tally.Add(id, difference.empty() ? Outcome::Agree : Outcome::KernelFailure, difference);
}

Tally CompareFramedSuite(const std::string& suite)
{
    Tally tally;
    for (const std::filesystem::path& file : CaseFiles(suite)) {
        CompareFramedCase(StereonTest::ReadJsonFile(file), tally);
    }
    tally.Write(suite + " against GeomAPI_IntSS");
    return tally;
}

// --- analytic-ssi -------------------------------------------------------------

/// The shape of an analytic-ssi configuration (corpus/analytic-ssi/README.md,
/// Configurations), or none for a Tier B quartic, which names no shape.
std::optional<Shape> ShapeOf(const std::string& configuration)
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

bool HasSphere(const JsonValue& surfaces)
{
    return std::ranges::any_of(surfaces.AsArray(),
                               [](const JsonValue& s) { return s["type"].AsString() == "sphere"; });
}

/// The configuration tests/intersect/CorpusTests.cpp holds the kernel to: the
/// nominal one, or the tolerant one where it is a declared coincidence or
/// the pair has a sphere.
std::string RequiredOf(const JsonValue& expected, const JsonValue& surfaces)
{
    const JsonValue& tolerant = expected["tolerant"];
    if (!tolerant.Has("configuration")) {
        return expected["nominal"]["configuration"].AsString();
    }
    const std::string& configuration = tolerant["configuration"].AsString();
    if (configuration == "coincident" || HasSphere(surfaces)) {
        return configuration;
    }
    return expected["nominal"]["configuration"].AsString();
}

void ComparePairCase(const JsonValue& c, Tally& tally)
{
    const std::string& id = c["id"].AsString();
    INFO(id);
    const JsonValue& surfaces = c["inputs"]["surfaces"];
    REQUIRE(surfaces.AsArray().size() == 2);
    const CaseSurface a = PairSurfaceOf(surfaces[0]);
    const CaseSurface b = PairSurfaceOf(surfaces[1]);
    const double tolerance =
        std::max(surfaces[0]["tolerance"].AsNumber(), surfaces[1]["tolerance"].AsNumber());
    const JsonValue& expected = c["expected"];
    const std::string nominal = expected["nominal"]["configuration"].AsString();
    const std::string tolerant = expected["tolerant"].Has("configuration")
                                     ? expected["tolerant"]["configuration"].AsString()
                                     : nominal;
    const std::optional<Shape> nominalShape = ShapeOf(nominal);
    if (!nominalShape) {
        tally.Add(id, Outcome::NotCompared);
        return;
    }
    const std::optional<Shape> tolerantShape = ShapeOf(tolerant);

    const OcctResult occt = IntersectWithOcct(a.Occt, b.Occt, tolerance);
    if (!occt.Done) {
        tally.Add(id, Outcome::OcctFailed, occt.Failure);
        return;
    }
    const std::string whyNominal = Mismatch(occt, *nominalShape, a, b, 1.0);
    const std::string whyTolerant =
        tolerantShape ? Mismatch(occt, *tolerantShape, a, b, 1.0) : whyNominal;
    std::vector<std::string> readings;
    if (whyNominal.empty()) {
        readings.push_back(nominal);
    }
    if (whyTolerant.empty()) {
        readings.push_back(tolerant);
    }
    if (readings.empty()) {
        std::string why = "nominal " + nominal + ": " + whyNominal;
        if (tolerant != nominal) {
            why += "tolerant " + tolerant + ": " + whyTolerant;
        }
        tally.Add(id, Outcome::OcctDiscrepancy, why);
        return;
    }
    if (expected.Has("kernel")) {
        tally.Add(id, Outcome::NotCompared); // Stereon refuses the case by contract
        return;
    }
    const std::string required = RequiredOf(expected, surfaces);
    if (required == "coincident") {
        tally.Add(id, Outcome::NotCompared);
        return;
    }
    if (std::ranges::find(readings, required) == readings.end()) {
        tally.Add(id, Outcome::OtherReading,
                  "OCCT answers " + readings.front() + ", Stereon is held to " + required);
        return;
    }
    const auto ours = Intersect(a.Stereon, a.Box, b.Stereon, b.Box, tolerance, IntersectTest::ctx);
    INFO(IntersectTest::ErrorText(ours));
    CHECK(ours.has_value());
    if (!ours) {
        tally.Add(id, Outcome::KernelFailure, IntersectTest::ErrorText(ours));
        return;
    }
    const std::string difference = Difference(*ours, occt, a.Anchor, 1.0, WitnessReach);
    INFO(difference);
    CHECK(difference.empty());
    tally.Add(id, difference.empty() ? Outcome::Agree : Outcome::KernelFailure, difference);
}

} // namespace

TEST_CASE("every analytic-ssi pair case agrees with its closed form and GeomAPI_IntSS",
          "[occt][intersect][corpus]")
{
    Tally tally;
    for (const std::filesystem::path& file : CaseFiles("analytic-ssi")) {
        const JsonValue c = StereonTest::ReadJsonFile(file);
        if (c["operation"]["name"].AsString() == "intersect") {
            ComparePairCase(c, tally);
        }
    }
    tally.Write("analytic-ssi against GeomAPI_IntSS");
    // The 1024 generated sweeps and the hand-written pairs.
    REQUIRE(tally.Total() >= 1050);
    CHECK(tally.Count(Outcome::KernelFailure) == 0);
}

TEST_CASE("every curved-ssi case agrees with its closed form and GeomAPI_IntSS",
          "[occt][intersect][corpus]")
{
    const Tally tally = CompareFramedSuite("curved-ssi");
    REQUIRE(tally.Total() >= 600);
    CHECK(tally.Count(Outcome::KernelFailure) == 0);
}

TEST_CASE("every nurbs-ssi case agrees with its closed form and GeomAPI_IntSS",
          "[occt][intersect][corpus]")
{
    const Tally tally = CompareFramedSuite("nurbs-ssi");
    REQUIRE(tally.Total() >= 90);
    CHECK(tally.Count(Outcome::KernelFailure) == 0);
}
