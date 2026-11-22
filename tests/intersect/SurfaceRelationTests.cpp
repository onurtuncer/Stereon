// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0053 relations MR-01 to MR-04 for every ordered pair of built-in
// surfaces: Plane, Cylinder, Cone, Sphere, Torus, Ellipsoid and
// BSplineSurface, 28 unordered pairs and 49 ordered ones.
//
// Each unordered pair has one transverse configuration whose answer is
// known in closed form, never from the kernel:
//
// - the fifteen pairs with a cone, a torus or an ellipsoid: a case of
//   corpus/curved-ssi/, at its origin placement and scale 1;
// - two B-spline surfaces: corpus/nurbs-ssi/cylinder-cylinder/two-loops;
// - the six pairs of planes, spheres and cylinders, and a B-spline plane
//   against each analytic surface: written here, the closed form beside
//   each.
//
// The claim for a configuration: no overlap; its count of isolated points;
// its count of connected curves; every branch crossing and certified within
// the resolution, its pcurves on both surfaces; every sample of every
// branch, and every point, on both analytic surfaces by their closed-form
// distances; and every witness (a point known to lie on the intersection)
// near the result. The relations, each run in both operand orders:
//
// - MR-02: the operands exchanged (every relation below also runs in both
//   orders, so each ordered pair has each relation);
// - MR-01: two cyclic permutations of the axes and a quarter turn about z,
//   exact, two of them with the translation (8, -4, 2);
// - MR-03: scales 1/4 and 4, exact, with the pair tolerance, the
//   resolution and the tolerance ceiling scaled;
// - MR-04: each operand in turn in an equivalent chart: an analytic
//   surface's frame turned half about its own axis (X and Y negated, so
//   u' = u - pi, or (u, v)' = -(u, v) on a plane, its box mapped alike); a
//   B-spline's u reversed (knots u' = a + b - u, rows of poles and weights
//   reversed). The pcurves must lie on the recharted surface.

#include "Json.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::Cone;
using Stereon::ContactType;
using Stereon::Context;
using Stereon::Cylinder;
using Stereon::Ellipsoid;
using Stereon::IntersectionBranch;
using Stereon::IntersectionPoint;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::SurfaceSurfaceIntersection;
using Stereon::ToleranceCause;
using Stereon::Torus;
using Stereon::Vec3;
using StereonTest::JsonValue;

namespace {

constexpr double pi = std::numbers::pi;

const std::filesystem::path corpus = std::filesystem::path(STEREON_CORPUS_DIR);

enum class Type : std::uint8_t { Plane, Sphere, Cylinder, Cone, Torus, Ellipsoid, Spline };

/// One operand: a surface's defining data, its requested box (none for the
/// whole chart). Lengths: R (radius, major or a), R2 (minor or b), R3 (c).
struct Spec {
    Type Kind = Type::Plane;
    Point3 Origin;
    Vec3 X{1, 0, 0};
    Vec3 Y{0, 1, 0};
    double R = 0.0;
    double R2 = 0.0;
    double R3 = 0.0;
    double Angle = 0.0;
    std::optional<Box2> Box;
    int DegreeU = 0;
    int DegreeV = 0;
    std::vector<double> KnotsU;
    std::vector<double> KnotsV;
    std::vector<Point3> Poles; // row by row, a row per u
    std::vector<double> Weights;
};

struct Expected {
    std::size_t Curves = 0;
    std::size_t Points = 0;
    std::vector<Point3> Witnesses;
};

struct Fixture {
    std::string Name;
    Spec A;
    Spec B;
    double Tolerance = 1e-5;
    Expected Want;
};

// --- The corpus -----------------------------------------------------------

Point3 PointOf(const JsonValue& v)
{
    return Point3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

Vec3 VectorOf(const JsonValue& v)
{
    return Vec3{v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

std::vector<double> NumbersOf(const JsonValue& v)
{
    std::vector<double> numbers;
    for (const JsonValue& n : v.AsArray()) {
        numbers.push_back(n.AsNumber());
    }
    return numbers;
}

Spec SpecOf(const JsonValue& s)
{
    Spec spec;
    const std::string& type = s["type"].AsString();
    if (s.Has("box")) {
        const JsonValue& box = s["box"];
        spec.Box = Box2{Interval{box[0][0].AsNumber(), box[0][1].AsNumber()},
                        Interval{box[1][0].AsNumber(), box[1][1].AsNumber()}};
    }
    if (type == "bspline") {
        spec.Kind = Type::Spline;
        spec.DegreeU = static_cast<int>(s["degree_u"].AsNumber());
        spec.DegreeV = static_cast<int>(s["degree_v"].AsNumber());
        spec.KnotsU = NumbersOf(s["knots_u"]);
        spec.KnotsV = NumbersOf(s["knots_v"]);
        for (const JsonValue& p : s["poles"].AsArray()) {
            spec.Poles.push_back(PointOf(p));
        }
        if (s.Has("weights")) {
            spec.Weights = NumbersOf(s["weights"]);
        }
        return spec;
    }
    spec.X = VectorOf(s["x_axis"]);
    spec.Y = VectorOf(s["y_axis"]);
    if (type == "plane") {
        spec.Kind = Type::Plane;
        spec.Origin = PointOf(s["origin"]);
    } else if (type == "sphere") {
        spec.Kind = Type::Sphere;
        spec.Origin = PointOf(s["center"]);
        spec.R = s["radius"].AsNumber();
    } else if (type == "cylinder") {
        spec.Kind = Type::Cylinder;
        spec.Origin = PointOf(s["origin"]);
        spec.R = s["radius"].AsNumber();
    } else if (type == "cone") {
        spec.Kind = Type::Cone;
        spec.Origin = PointOf(s["origin"]);
        spec.R = s["radius"].AsNumber();
        spec.Angle = s["semi_angle"].AsNumber();
    } else if (type == "torus") {
        spec.Kind = Type::Torus;
        spec.Origin = PointOf(s["center"]);
        spec.R = s["major"].AsNumber();
        spec.R2 = s["minor"].AsNumber();
    } else {
        spec.Kind = Type::Ellipsoid;
        spec.Origin = PointOf(s["center"]);
        spec.R = s["a"].AsNumber();
        spec.R2 = s["b"].AsNumber();
        spec.R3 = s["c"].AsNumber();
    }
    return spec;
}

/// A corpus case at its origin placement and scale 1, with its expected
/// answer, which tools/generate_surface_corpus.py or
/// tools/generate_nurbs_corpus.py derives in closed form.
Fixture FromCorpus(const std::string& suite, const std::string& family, const std::string& config)
{
    const JsonValue c =
        StereonTest::ReadJsonFile(corpus / suite / family / (config + "-origin-x1.json"));
    const JsonValue& surfaces = c["inputs"]["surfaces"];
    const JsonValue& expected = c["expected"];
    // Only answered, transverse configurations are used here.
    REQUIRE_FALSE(expected["overlap"].AsBool());
    REQUIRE_FALSE(expected.Has("may_fail_closed"));
    REQUIRE(expected["contact"].AsString() == "crossing");
    REQUIRE(c["inputs"]["scale"].AsNumber() == 1.0);
    Fixture f;
    f.Name = suite + "/" + family + "/" + config;
    f.A = SpecOf(surfaces[0]);
    f.B = SpecOf(surfaces[1]);
    f.Tolerance = c["inputs"]["tolerance"].AsNumber();
    f.Want.Curves = static_cast<std::size_t>(expected["curves"].AsNumber());
    f.Want.Points = static_cast<std::size_t>(expected["points"].AsNumber());
    for (const JsonValue& w : expected["witnesses"].AsArray()) {
        f.Want.Witnesses.push_back(PointOf(w));
    }
    return f;
}

// --- Configurations written here ------------------------------------------

const Vec3 ex{1, 0, 0};
const Vec3 ey{0, 1, 0};
const Vec3 ez{0, 0, 1};

Spec PlaneSpec(const Point3& origin, const Vec3& x, const Vec3& y, std::optional<Box2> box)
{
    Spec s;
    s.Kind = Type::Plane;
    s.Origin = origin;
    s.X = x;
    s.Y = y;
    s.Box = box;
    return s;
}

Spec RoundSpec(const Type kind, const Point3& origin, const Vec3& x, const Vec3& y, const double r,
               std::optional<Box2> box = std::nullopt)
{
    Spec s;
    s.Kind = kind;
    s.Origin = origin;
    s.X = x;
    s.Y = y;
    s.R = r;
    s.Box = box;
    return s;
}

/// A cylinder's or a cone's box: a whole turn over the heights given.
Box2 TurnBox(const double low, const double high)
{
    return Box2{Interval{0.0, 2.0 * pi}, Interval{low, high}};
}

/// The bilinear B-spline that is exactly the plane z = 1/2 over [-6, 6]^2.
Spec FlatSpline()
{
    Spec s;
    s.Kind = Type::Spline;
    s.DegreeU = 1;
    s.DegreeV = 1;
    s.KnotsU = {0, 0, 1, 1};
    s.KnotsV = {0, 0, 1, 1};
    for (const double x : {-6.0, 6.0}) {
        for (const double y : {-6.0, 6.0}) {
            s.Poles.push_back(Point3{x, y, 0.5});
        }
    }
    return s;
}

std::vector<Fixture> Written()
{
    const double r3 = std::sqrt(0.75);
    const Box2 square{Interval{-4.0, 4.0}, Interval{-4.0, 4.0}};
    std::vector<Fixture> out;
    // z = 0 against x = 1/2 over y in [-3, 3]: the segment x = 1/2, z = 0.
    out.push_back(
        {"plane-plane",
         PlaneSpec(Point3{}, ex, ey, square),
         PlaneSpec(Point3{0.5, 0, 0}, ey, ez, Box2{Interval{-3.0, 3.0}, Interval{-4.0, 4.0}}),
         1e-5,
         {1, 0, {Point3{0.5, 0, 0}, Point3{0.5, 2, 0}, Point3{0.5, -2.5, 0}}}});
    // z = 1/2 cuts the unit sphere in the circle of radius sqrt(3/4).
    out.push_back({"plane-sphere",
                   PlaneSpec(Point3{0, 0, 0.5}, ex, ey, square),
                   RoundSpec(Type::Sphere, Point3{}, ex, ey, 1.0),
                   1e-5,
                   {1, 0, {Point3{r3, 0, 0.5}, Point3{0, -r3, 0.5}, Point3{-r3, 0, 0.5}}}});
    // The plane through (0, 0, 1/2) with axes x and (0, 3/5, 4/5), normal
    // (0, -4/5, 3/5), cuts the unit cylinder about z in the ellipse
    // z = 1/2 + 4y/3.
    out.push_back({"plane-cylinder",
                   PlaneSpec(Point3{0, 0, 0.5}, ex, Vec3{0, 0.6, 0.8}, square),
                   RoundSpec(Type::Cylinder, Point3{}, ex, ey, 1.0, TurnBox(-4.0, 4.0)),
                   1e-5,
                   {1,
                    0,
                    {Point3{1, 0, 0.5}, Point3{-1, 0, 0.5}, Point3{0, 1, 0.5 + 4.0 / 3.0},
                     Point3{0, -1, 0.5 - 4.0 / 3.0}}}});
    // Unit spheres a unit apart along x meet in the circle x = 1/2, radius
    // sqrt(3/4).
    out.push_back({"sphere-sphere",
                   RoundSpec(Type::Sphere, Point3{}, ex, ey, 1.0),
                   RoundSpec(Type::Sphere, Point3{1, 0, 0}, ex, ey, 1.0),
                   1e-5,
                   {1, 0, {Point3{0.5, r3, 0}, Point3{0.5, 0, -r3}, Point3{0.5, -r3, 0}}}});
    // The unit cylinder about z and the unit sphere about (3/2, 0, 0):
    // z^2 = 3x - 9/4 on the cylinder, one loop over x in [3/4, 1], away from
    // the sphere's poles.
    out.push_back(
        {"sphere-cylinder",
         RoundSpec(Type::Sphere, Point3{1.5, 0, 0}, ex, ey, 1.0),
         RoundSpec(Type::Cylinder, Point3{}, ex, ey, 1.0, TurnBox(-4.0, 4.0)),
         1e-5,
         {1, 0, {Point3{1, 0, r3}, Point3{1, 0, -r3}, Point3{0.75, std::sqrt(0.4375), 0}}}});
    // The unit cylinder about z and the cylinder of radius 1/2 about x: two
    // loops, x = +-sqrt(1 - y^2) with y^2 + z^2 = 1/4.
    out.push_back({"cylinder-cylinder",
                   RoundSpec(Type::Cylinder, Point3{}, ex, ey, 1.0, TurnBox(-4.0, 4.0)),
                   RoundSpec(Type::Cylinder, Point3{}, ey, ez, 0.5, TurnBox(-4.0, 4.0)),
                   1e-5,
                   {2, 0, {Point3{1, 0, 0.5}, Point3{-1, 0, -0.5}, Point3{r3, 0.5, 0}}}});
    // A B-spline plane z = 1/2 against each analytic surface: the plane
    // section's closed form (as in SplineAnalyticTests, now also turned).
    const double h = 0.5;
    out.push_back({"bspline-plane",
                   FlatSpline(),
                   PlaneSpec(Point3{0.5, 0, 0}, ey, ez, std::nullopt),
                   1e-5,
                   {1, 0, {Point3{0.5, 0, h}, Point3{0.5, 5, h}}}});
    out.push_back({"bspline-cylinder",
                   FlatSpline(),
                   RoundSpec(Type::Cylinder, Point3{}, ex, ey, 2.0),
                   1e-5,
                   {1, 0, {Point3{2, 0, h}, Point3{0, -2, h}}}});
    {
        const double a = 0.5;
        const double rho = 1.0 + h * std::tan(a);
        Spec cone = RoundSpec(Type::Cone, Point3{}, ex, ey, 1.0);
        cone.Angle = a;
        out.push_back({"bspline-cone",
                       FlatSpline(),
                       cone,
                       1e-5,
                       {1, 0, {Point3{rho, 0, h}, Point3{0, rho, h}}}});
    }
    {
        const double rho = std::sqrt(4.0 - h * h);
        out.push_back({"bspline-sphere",
                       FlatSpline(),
                       RoundSpec(Type::Sphere, Point3{}, ex, ey, 2.0),
                       1e-5,
                       {1, 0, {Point3{rho, 0, h}, Point3{-rho, 0, h}}}});
    }
    {
        const double w = std::sqrt(1.0 - h * h);
        Spec torus = RoundSpec(Type::Torus, Point3{}, ex, ey, 3.0);
        torus.R2 = 1.0;
        out.push_back({"bspline-torus",
                       FlatSpline(),
                       torus,
                       1e-5,
                       {2, 0, {Point3{3.0 - w, 0, h}, Point3{0, 3.0 + w, h}}}});
    }
    {
        const double k = std::sqrt(1.0 - (h / 1.5) * (h / 1.5));
        Spec ellipsoid = RoundSpec(Type::Ellipsoid, Point3{}, ex, ey, 3.0);
        ellipsoid.R2 = 2.0;
        ellipsoid.R3 = 1.5;
        out.push_back({"bspline-ellipsoid",
                       FlatSpline(),
                       ellipsoid,
                       1e-5,
                       {1, 0, {Point3{3.0 * k, 0, h}, Point3{0, 2.0 * k, h}}}});
    }
    return out;
}

/// Every unordered pair of built-in surfaces, once.
std::vector<Fixture> Fixtures()
{
    std::vector<Fixture> out = Written();
    const std::array<std::array<const char*, 2>, 15> curved{{
        {"cone-cone", "coaxial-unequal"},
        {"cone-ellipsoid", "both-nappes"},
        {"cone-torus", "coaxial-two"},
        {"cylinder-cone", "cone-axis-inside"},
        {"cylinder-ellipsoid", "between-sides"},
        {"cylinder-torus", "over-the-tube"},
        {"ellipsoid-ellipsoid", "shifted-x-2"},
        {"plane-cone", "tilted-20"},
        {"plane-ellipsoid", "oblique-1p0"},
        {"plane-torus", "two-ovals"},
        {"sphere-cone", "off-axis-loop"},
        {"sphere-ellipsoid", "around-long-axis"},
        {"sphere-torus", "across-the-tube"},
        {"torus-ellipsoid", "through-the-tube"},
        {"torus-torus", "shifted-1"},
    }};
    for (const auto& [family, config] : curved) {
        out.push_back(FromCorpus("curved-ssi", family, config));
    }
    out.push_back(FromCorpus("nurbs-ssi", "cylinder-cylinder", "two-loops"));
    REQUIRE(out.size() == 28);
    return out;
}

// --- Transformations -------------------------------------------------------

/// p -> Offset + Scale Q p, with Q an exact rotation given by its columns.
struct Motion {
    std::array<Vec3, 3> Columns{ex, ey, ez};
    Vec3 Offset{};
    double Scale = 1.0;

    [[nodiscard]] Vec3 Turn(const Vec3& v) const
    {
        return v.X * Columns[0] + v.Y * Columns[1] + v.Z * Columns[2];
    }
    [[nodiscard]] Point3 Of(const Point3& p) const
    {
        return Point3{} + Offset + Scale * Turn(p - Point3{});
    }
    [[nodiscard]] Context Policy() const
    {
        return Context::Default()
            .WithResolution(Context::DefaultResolution * Scale)
            .WithToleranceCeiling(Context::DefaultToleranceCeiling * Scale);
    }
};

Interval Scaled(const Interval& i, const double s)
{
    return Interval{s * i.GetLower(), s * i.GetUpper()};
}

Spec Moved(Spec s, const Motion& m)
{
    if (s.Kind == Type::Spline) {
        for (Point3& p : s.Poles) {
            p = m.Of(p);
        }
        return s;
    }
    s.Origin = m.Of(s.Origin);
    s.X = m.Turn(s.X);
    s.Y = m.Turn(s.Y);
    s.R *= m.Scale;
    s.R2 *= m.Scale;
    s.R3 *= m.Scale;
    if (s.Box) {
        // A plane's parameters are lengths; a cylinder's and a cone's v.
        if (s.Kind == Type::Plane) {
            s.Box = Box2{Scaled(s.Box->X, m.Scale), Scaled(s.Box->Y, m.Scale)};
        } else if (s.Kind == Type::Cylinder || s.Kind == Type::Cone) {
            s.Box = Box2{s.Box->X, Scaled(s.Box->Y, m.Scale)};
        }
    }
    return s;
}

/// The same surface in an equivalent chart (see the file's comment).
Spec Recharted(Spec s)
{
    if (s.Kind == Type::Spline) {
        const double a = s.KnotsU.front();
        const double b = s.KnotsU.back();
        std::vector<double> knots;
        for (auto k = s.KnotsU.rbegin(); k != s.KnotsU.rend(); ++k) {
            knots.push_back(a + b - *k);
        }
        s.KnotsU = std::move(knots);
        const std::size_t columns = s.KnotsV.size() - static_cast<std::size_t>(s.DegreeV) - 1;
        const std::size_t rows = s.Poles.size() / columns;
        std::vector<Point3> poles;
        std::vector<double> weights;
        for (std::size_t i = rows; i-- > 0;) {
            for (std::size_t j = 0; j < columns; ++j) {
                poles.push_back(s.Poles[i * columns + j]);
                if (!s.Weights.empty()) {
                    weights.push_back(s.Weights[i * columns + j]);
                }
            }
        }
        s.Poles = std::move(poles);
        s.Weights = std::move(weights);
        if (s.Box) {
            s.Box =
                Box2{Interval{a + b - s.Box->X.GetUpper(), a + b - s.Box->X.GetLower()}, s.Box->Y};
        }
        return s;
    }
    s.X = -s.X;
    s.Y = -s.Y;
    if (s.Box) {
        if (s.Kind == Type::Plane) {
            s.Box = Box2{Interval{-s.Box->X.GetUpper(), -s.Box->X.GetLower()},
                         Interval{-s.Box->Y.GetUpper(), -s.Box->Y.GetLower()}};
        } else {
            s.Box = Box2{Interval{s.Box->X.GetLower() - pi, s.Box->X.GetUpper() - pi}, s.Box->Y};
        }
    }
    return s;
}

Surface SurfaceOf(const Spec& s)
{
    switch (s.Kind) {
    case Type::Plane:
        return Surface{Plane{s.Origin, s.X, s.Y}};
    case Type::Sphere:
        return Surface{Sphere{s.Origin, s.R, s.X, s.Y, true}};
    case Type::Cylinder:
        return Surface{Cylinder{s.Origin, s.R, s.X, s.Y, true}};
    case Type::Cone:
        return Surface{Cone{s.Origin, s.R, s.Angle, s.X, s.Y, true}};
    case Type::Torus:
        return Surface{Torus{s.Origin, s.R, s.R2, s.X, s.Y, true}};
    case Type::Ellipsoid:
        return Surface{Ellipsoid{s.Origin, s.R, s.R2, s.R3, s.X, s.Y, true}};
    case Type::Spline:
        break;
    }
    auto made =
        Stereon::MakeBSplineSurface(s.DegreeU, s.DegreeV, s.KnotsU, s.KnotsV, s.Poles, s.Weights);
    REQUIRE(made.has_value());
    return Surface{*std::move(made)};
}

/// The closed-form distance of a point from an analytic surface (for the
/// ellipsoid, its algebraic value times the largest semi-axis over two, an
/// estimate from above near it); none for a B-spline.
std::optional<double> Off(const Spec& s, const Point3& p)
{
    const Vec3 w = p - s.Origin;
    const Vec3 z = Cross(s.X, s.Y);
    const double x = Dot(w, s.X);
    const double y = Dot(w, s.Y);
    const double h = Dot(w, z);
    const double rho = std::hypot(x, y);
    switch (s.Kind) {
    case Type::Plane:
        return std::abs(h);
    case Type::Sphere:
        return std::abs(Length(w) - s.R);
    case Type::Cylinder:
        return std::abs(rho - s.R);
    case Type::Cone:
        return std::abs(rho - std::abs(s.R + h * std::tan(s.Angle))) * std::cos(s.Angle);
    case Type::Torus:
        return std::abs(std::hypot(rho - s.R, h) - s.R2);
    case Type::Ellipsoid: {
        const double f = x * x / (s.R * s.R) + y * y / (s.R2 * s.R2) + h * h / (s.R3 * s.R3) - 1.0;
        return std::abs(f) * std::max({s.R, s.R2, s.R3}) / 2.0;
    }
    case Type::Spline:
        break;
    }
    return std::nullopt;
}

// --- The claim ---------------------------------------------------------------

/// Segments per branch when a witness is looked for.
constexpr std::size_t Samples = 2048;
/// How near a witness must lie to the result, per unit of scale: well
/// above the sagitta of the samples and the witnesses' ten decimals, far
/// below any feature of a configuration.
constexpr double WitnessReach = 1e-3;

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

double DistanceToResult(const SurfaceSurfaceIntersection& r, const Point3& p)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const IntersectionPoint& point : r.Points) {
        nearest = std::min(nearest, Distance(p, point.Position));
    }
    for (const IntersectionBranch& branch : r.Branches) {
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

void CheckOn(const Spec& s, const Point3& p, const double reach)
{
    if (const std::optional<double> off = Off(s, p)) {
        CHECK(*off <= reach);
    }
}

enum class Chart : std::uint8_t { Given, OfA, OfB };

/// The claim for one configuration under one motion, order and chart.
void Claim(const Fixture& f, const Motion& m, const bool swapped, const Chart chart = Chart::Given)
{
    Spec a = Moved(f.A, m);
    Spec b = Moved(f.B, m);
    if (chart == Chart::OfA) {
        a = Recharted(a);
    } else if (chart == Chart::OfB) {
        b = Recharted(b);
    }
    const Spec& first = swapped ? b : a;
    const Spec& second = swapped ? a : b;
    INFO(f.Name << (swapped ? " swapped" : "") << " scale=" << m.Scale << " offset=(" << m.Offset.X
                << ", " << m.Offset.Y << ", " << m.Offset.Z << ") turn=(" << m.Columns[0].X << ", "
                << m.Columns[0].Y << ", " << m.Columns[0].Z
                << ") chart=" << static_cast<int>(chart));
    const Surface sa = SurfaceOf(first);
    const Surface sb = SurfaceOf(second);
    const Context policy = m.Policy();
    const double resolution = policy.GetResolution();
    const double tolerance = f.Tolerance * m.Scale;
    const auto r = Intersect(sa, first.Box.value_or(Box2::Everything()), sb,
                             second.Box.value_or(Box2::Everything()), tolerance, policy);
    // A refusal of a configuration that is answered at the origin is a
    // relation failure, never a pass.
    INFO(IntersectTest::ErrorText(r));
    REQUIRE(r.has_value());
    CHECK(r->Overlaps.empty());
    CHECK(r->Points.size() == f.Want.Points);
    CHECK(IntersectTest::ConnectivityOf(*r, 2.0 * tolerance).Curves == f.Want.Curves);
    const double pointReach = 4.0 * (tolerance + resolution);
    for (const IntersectionPoint& point : r->Points) {
        CheckOn(first, point.Position, pointReach);
        CheckOn(second, point.Position, pointReach);
    }
    for (const IntersectionBranch& branch : r->Branches) {
        CHECK(branch.Contact == ContactType::Crossing);
        CHECK(branch.Range.IsBounded());
        CHECK(branch.IsCertified());
        CHECK(branch.BoundCause == ToleranceCause::CertificateWidening);
        CHECK(branch.Bound <= resolution);
        CHECK(IntersectTest::PcurveDeviation(branch, sa, sb, 100) <= branch.Bound + tolerance);
        const double reach = 4.0 * (branch.Bound + tolerance);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            const Point3 p = Eval(branch.Curve3d, t);
            CheckOn(first, p, reach);
            CheckOn(second, p, reach);
        }
    }
    for (const Point3& w : f.Want.Witnesses) {
        const Point3 placed = m.Of(w);
        INFO("witness (" << placed.X << ", " << placed.Y << ", " << placed.Z << ")");
        CHECK(DistanceToResult(*r, placed) <= WitnessReach * m.Scale);
    }
}

const std::array<Vec3, 3> Cyclic{ey, ez, ex};
const std::array<Vec3, 3> Cyclic2{ez, ex, ey};
const std::array<Vec3, 3> QuarterZ{ey, Vec3{-1, 0, 0}, ez};
const Vec3 Far{8, -4, 2};

} // namespace

TEST_CASE("MR-02: surface-surface answers hold with the operands exchanged",
          "[intersect][surface-surface][metamorphic]")
{
    for (const Fixture& f : Fixtures()) {
        Claim(f, {}, false);
        Claim(f, {}, true);
    }
}

TEST_CASE("MR-01: surface-surface answers move with exact rigid motions",
          "[intersect][surface-surface][metamorphic]")
{
    const std::array<Motion, 3> motions{Motion{Cyclic, Far, 1.0}, Motion{Cyclic2, Vec3{}, 1.0},
                                        Motion{QuarterZ, Far, 1.0}};
    for (const Fixture& f : Fixtures()) {
        for (const Motion& m : motions) {
            Claim(f, m, false);
            Claim(f, m, true);
        }
    }
}

TEST_CASE("MR-03: surface-surface answers scale with the geometry and its tolerances",
          "[intersect][surface-surface][metamorphic]")
{
    for (const Fixture& f : Fixtures()) {
        for (const double scale : {0.25, 4.0}) {
            Motion m;
            m.Scale = scale;
            Claim(f, m, false);
            Claim(f, m, true);
        }
    }
}

TEST_CASE("MR-04: surface-surface answers keep their place in equivalent surface charts",
          "[intersect][surface-surface][metamorphic]")
{
    for (const Fixture& f : Fixtures()) {
        for (const Chart chart : {Chart::OfA, Chart::OfB}) {
            Claim(f, {}, false, chart);
            Claim(f, {}, true, chart);
        }
    }
}
