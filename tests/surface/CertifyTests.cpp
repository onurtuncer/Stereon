// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Certificates (ADR-0009 section 4, the Boolean prototype's M4):
// IntervalScalar is a Scalar whose every result encloses the exact one, its
// sine and cosine enclose their extrema, the geometry templates evaluate
// over it, and CertifyDeviation bounds a curve's deviation from a surface
// along a pcurve over a whole range, refusing what it cannot bound.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

import Stereon.Surface;

using Catch::Matchers::WithinAbs;
using Stereon::Circle;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::IntervalScalar;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;

const Surface xy{Plane{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
const Surface aboutZ{Cylinder{Point3{}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};

bool Encloses(const IntervalScalar& s, const double x)
{
    return Contains(s.Value, x);
}

} // namespace

TEST_CASE("IntervalScalar is a Scalar whose results enclose the exact ones", "[surface][certify]")
{
    STATIC_REQUIRE(Stereon::Scalar<IntervalScalar>);
    REQUIRE(IntervalScalar{} == IntervalScalar{0.0});
    const IntervalScalar a{0.1};
    const IntervalScalar b{0.2};
    REQUIRE(Encloses(a + b, 0.1 + 0.2));
    REQUIRE(Encloses(a + b, 0.3)); // The exact sum, which the double 0.1 + 0.2 misses.
    REQUIRE(Encloses(a * b, 0.02));
    REQUIRE(Encloses(a / b, 0.5));
    REQUIRE(Encloses(-a, -0.1));
    REQUIRE(Encloses(sqrt(IntervalScalar{2.0}), std::sqrt(2.0)));
    // The certain ordering: undecided when the intervals meet.
    REQUIRE(a < b);
    REQUIRE_FALSE(IntervalScalar{Interval{0.0, 1.0}} < IntervalScalar{Interval{0.5, 2.0}});
    REQUIRE_FALSE(IntervalScalar{Interval{0.5, 2.0}} < IntervalScalar{Interval{0.0, 1.0}});
}

TEST_CASE("the trigonometric enclosures cover the extrema inside a span", "[surface][certify]")
{
    const auto over = [](const double lower, const double upper) {
        return IntervalScalar{Interval{lower, upper}};
    };
    // A span through pi/2 reaches 1; one through 3 pi/2 reaches -1.
    REQUIRE(sin(over(1.0, 2.0)).Value.GetUpper() >= 1.0);
    REQUIRE(sin(over(4.0, 5.0)).Value.GetLower() <= -1.0);
    REQUIRE(cos(over(3.0, 3.2)).Value.GetLower() <= -1.0);
    REQUIRE(cos(over(-0.1, 0.1)).Value.GetUpper() >= 1.0);
    // A span short of an extremum stays short of it.
    REQUIRE(cos(over(0.1, 0.2)).Value.GetUpper() < 1.0);
    REQUIRE(sin(over(0.1, 0.2)).Value.GetLower() > 0.0);
    // Negative angles, as a branch starting before a cylinder's seam has.
    REQUIRE(sin(over(-2.0, -1.0)).Value.GetLower() <= -1.0);
    REQUIRE(cos(over(-3.2, -3.0)).Value.GetLower() <= -1.0);
    REQUIRE(cos(over(-0.2, -0.1)).Value.GetUpper() < 1.0);
    // A whole period is everything; a point is tight.
    REQUIRE(sin(over(0.0, 7.0)).Value == Interval{-1.0, 1.0});
    REQUIRE(Width(sin(IntervalScalar{0.7}).Value) < 1e-15);
    REQUIRE(Encloses(sin(IntervalScalar{0.7}), std::sin(0.7)));
    // Every sample of a span lies within its enclosure.
    for (int k = 0; k <= 100; ++k) {
        const double t = 0.3 + 2.5 * k / 100.0;
        REQUIRE(Encloses(sin(over(0.3, 2.8)), std::sin(t)));
        REQUIRE(Encloses(cos(over(0.3, 2.8)), std::cos(t)));
    }
}

TEST_CASE("a circle evaluated over a parameter interval encloses its points", "[surface][certify]")
{
    const Circle circle{Point3{1.0, 2.0, 3.0}, 2.5, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const Stereon::BasicCircle<IntervalScalar, 3> wide{
        {IntervalScalar{1.0}, IntervalScalar{2.0}, IntervalScalar{3.0}},
        IntervalScalar{2.5},
        {IntervalScalar{0.0}, IntervalScalar{1.0}, IntervalScalar{0.0}},
        {IntervalScalar{0.0}, IntervalScalar{0.0}, IntervalScalar{1.0}}};
    const Interval span{0.4, 1.1};
    const auto p = Eval(wide, IntervalScalar{span});
    const auto d = Deriv(wide, IntervalScalar{span}, 1);
    for (int k = 0; k <= 50; ++k) {
        const double t = span.GetLower() + Width(span) * k / 50.0;
        const Point3 q = Eval(circle, t);
        const Vec3 v = Deriv(circle, t, 1);
        REQUIRE((Encloses(p.X, q.X) && Encloses(p.Y, q.Y) && Encloses(p.Z, q.Z)));
        REQUIRE((Encloses(d.X, v.X) && Encloses(d.Y, v.Y) && Encloses(d.Z, v.Z)));
    }
}

TEST_CASE("a line on a plane and a circle on a cylinder are certified to tiny bounds",
          "[surface][certify]")
{
    const Curve line{Line{Point3{1.0, 2.0, 0.0}, Normalized(Vec3{3.0, 4.0, 0.0})}};
    const Curve2d onPlane{Line2d{Point2{1.0, 2.0}, Vec2{0.6, 0.8}}};
    const auto flat = CertifyDeviation(line, xy, onPlane, Interval{-10.0, 10.0}, 1e-6);
    REQUIRE(flat.has_value());
    REQUIRE(flat->Bound <= 1e-13);

    const Curve circle{
        Circle{Point3{0.0, 0.0, 3.0}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Curve2d onCylinder{Line2d{Point2{0.0, 3.0}, Vec2{1.0, 0.0}}};
    const auto round = CertifyDeviation(circle, aboutZ, onCylinder, Interval{0.0, 2.0 * pi}, 1e-5);
    REQUIRE(round.has_value());
    // A circle along a unit-speed turn is of the closed form, so the whole
    // turn is proved at once and the bound is the rounding of the
    // coefficients, where the bisection needed 4096 spans to get within the
    // target (issue #155).
    REQUIRE(round->Bound <= 1e-13);
    REQUIRE(round->Spans == 1);
}

namespace {

/// The largest |c(t) - S(p(t))| at n + 1 parameters of the range: what
/// every proved bound must be at least.
double SampledDeviation(const Curve& c, const Surface& s, const Curve2d& p, const Interval& range,
                        const int n = 2000)
{
    double worst = 0.0;
    for (int k = 0; k <= n; ++k) {
        const double t = range.GetLower() + Width(range) * k / n;
        worst = std::max(worst, Distance(Eval(c, t), Eval(s, Eval(p, t))));
    }
    return worst;
}

} // namespace

TEST_CASE("the closed-form certificate never bounds below the deviation", "[surface][certify]")
{
    // Pcurves moved off their curves by known amounts: the one-step bound
    // must cover the sampled deviation at every amount, and the result must
    // agree with what the bisection proves or refutes.
    const Curve circle{
        Circle{Point3{0.0, 0.0, 3.0}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Curve ellipse{Stereon::Ellipse{Point3{1.0, -2.0, 0.0}, 4.0, 1.5, Vec3{0.6, 0.8, 0.0},
                                         Vec3{-0.8, 0.6, 0.0}}};
    for (const double off : {0.0, 1e-9, 1e-7, 3e-6, 1e-4}) {
        INFO("offset " << off);
        const Interval turn{-1.0, 2.0 * pi - 1.0};
        const std::vector<std::pair<Curve, std::pair<Surface, Curve2d>>> cases{
            {circle, {aboutZ, Curve2d{Line2d{Point2{0.0, 3.0 + off}, Vec2{1.0, 0.0}}}}},
            {circle, {aboutZ, Curve2d{Line2d{Point2{off, 3.0}, Vec2{1.0, 0.0}}}}},
            {ellipse,
             {xy, Curve2d{Stereon::Ellipse2d{Point2{1.0 + off, -2.0}, 4.0, 1.5, Vec2{0.6, 0.8},
                                             Vec2{-0.8, 0.6}}}}},
            {ellipse,
             {xy, Curve2d{Stereon::Ellipse2d{Point2{1.0, -2.0}, 4.0 + off, 1.5, Vec2{0.6, 0.8},
                                             Vec2{-0.8, 0.6}}}}},
        };
        for (const auto& [curve, onto] : cases) {
            const auto& [surface, pcurve] = onto;
            const double sampled = SampledDeviation(curve, surface, pcurve, turn);
            const auto loose = CertifyDeviation(curve, surface, pcurve, turn, 1.0);
            REQUIRE(loose.has_value());
            REQUIRE(loose->Bound >= sampled);
            const auto tight = CertifyDeviation(curve, surface, pcurve, turn, 1e-5);
            REQUIRE(tight.has_value() == (sampled <= 1e-5));
            if (tight) {
                REQUIRE(tight->Bound >= sampled);
                REQUIRE(tight->Bound <= 1e-5);
            }
        }
    }
}

TEST_CASE("a helix on a cylinder is left to the bisection", "[surface][certify]")
{
    // u = t / 2 is a line of the parameter plane but not a unit-speed turn:
    // its image is not of the closed form, so the bisection proves it.
    const Curve circle{
        Circle{Point3{0.0, 0.0, 3.0}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Curve2d slow{Line2d{Point2{0.0, 3.0}, Vec2{0.5, 0.0}}};
    const auto result = CertifyDeviation(circle, aboutZ, slow, Interval{0.0, 1.0}, 1e-5);
    REQUIRE_FALSE(result.has_value()); // the curves part: the circle runs twice as fast
    REQUIRE(result.error().GetCode() == ErrorCode::CertificationFailed);
}

TEST_CASE("a pcurve off its curve fails certification, or is bounded above its error",
          "[surface][certify]")
{
    // The circle sits at height 3; the pcurve claims 3.001.
    const Curve circle{
        Circle{Point3{0.0, 0.0, 3.0}, 5.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}}};
    const Curve2d wrong{Line2d{Point2{0.0, 3.001}, Vec2{1.0, 0.0}}};
    const auto failed = CertifyDeviation(circle, aboutZ, wrong, Interval{0.0, 2.0 * pi}, 1e-5);
    REQUIRE_FALSE(failed.has_value());
    REQUIRE(failed.error().GetCode() == ErrorCode::CertificationFailed);
    const auto loose = CertifyDeviation(circle, aboutZ, wrong, Interval{0.0, 2.0 * pi}, 2e-3);
    REQUIRE(loose.has_value());
    REQUIRE(loose->Bound >= 1e-3);
    REQUIRE(loose->Bound <= 2e-3);
}

TEST_CASE("B-spline curves and pcurves are certified through their poles", "[surface][certify]")
{
    // A cubic through samples of a line, as a pcurve and as a 3D curve.
    std::vector<Point2> samples2d;
    std::vector<Point3> samples3d;
    for (int k = 0; k <= 8; ++k) {
        const double t = k / 8.0;
        samples2d.emplace_back(1.0 + 3.0 * t, 2.0 + 4.0 * t);
        samples3d.emplace_back(1.0 + 3.0 * t, 2.0 + 4.0 * t, 0.0);
    }
    const auto fit2d = Interpolate(samples2d, 3, Stereon::Parametrization::Uniform);
    const auto fit3d = Interpolate(samples3d, 3, Stereon::Parametrization::Uniform);
    REQUIRE((fit2d.has_value() && fit3d.has_value()));
    const Curve line{Line{Point3{1.0, 2.0, 0.0}, Vec3{3.0, 4.0, 0.0}}};
    const Curve2d onPlane{Line2d{Point2{1.0, 2.0}, Vec2{3.0, 4.0}}};
    const Interval unit{0.0, 1.0};
    const auto splinePcurve = CertifyDeviation(line, xy, Curve2d{fit2d->Curve}, unit, 1e-6);
    REQUIRE(splinePcurve.has_value());
    REQUIRE(splinePcurve->Bound <= 1e-6);
    const auto splineCurve = CertifyDeviation(Curve{fit3d->Curve}, xy, onPlane, unit, 1e-6);
    REQUIRE(splineCurve.has_value());
    REQUIRE(splineCurve->Bound <= 1e-6);
    REQUIRE(splineCurve->Bound < 1e-9);
}

TEST_CASE("certification covers parallels on cones, spheres and tori", "[surface][certify]")
{
    const Vec3 x{1.0, 0.0, 0.0};
    const Vec3 y{0.0, 1.0, 0.0};
    const Interval turn{0.0, 2.0 * std::numbers::pi};
    const Curve2d parallel{Line2d{Point2{0.0, 0.5}, Vec2{1.0, 0.0}}};
    const auto check = [&](const Surface& surface, const Circle& circle) {
        const auto certificate = CertifyDeviation(Curve{circle}, surface, parallel, turn, 1e-9);
        INFO((certificate ? std::string() : std::string(certificate.error().GetMessage())));
        REQUIRE(certificate.has_value());
        REQUIRE(certificate->Bound <= 1e-9);
    };
    // The sphere's latitude 0.5, the torus's tube angle 0.5, the cone's
    // generator distance 0.5.
    check(Surface{Stereon::Sphere{Point3{}, 2.0, x, y, true}},
          Circle{Point3{0.0, 0.0, 2.0 * std::sin(0.5)}, 2.0 * std::cos(0.5), x, y});
    check(Surface{Stereon::Torus{Point3{}, 3.0, 1.0, x, y, true}},
          Circle{Point3{0.0, 0.0, std::sin(0.5)}, 3.0 + std::cos(0.5), x, y});
    const double a = std::numbers::pi / 6.0;
    check(Surface{Stereon::Cone{Point3{}, 1.0, a, x, y, true}},
          Circle{Point3{0.0, 0.0, 0.5 * std::cos(a)}, 1.0 + 0.5 * std::sin(a), x, y});
}

TEST_CASE("certification refuses what it cannot bound", "[surface][certify]")
{
    const Curve line{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve2d onPlane{Line2d{Point2{}, Vec2{1.0, 0.0}}};
    // A parabola as the pcurve has no interval trace.
    const Curve2d parabola{Stereon::Parabola2d{Point2{}, 1.0, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}};
    REQUIRE(CertifyDeviation(line, xy, parabola, Interval{0.0, 1.0}, 1e-6).error().GetCode() ==
            ErrorCode::NotSupported);
    REQUIRE(CertifyDeviation(line, xy, onPlane, Interval::Everything(), 1e-6).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(CertifyDeviation(line, xy, onPlane, Interval{0.0, 1.0}, 0.0).error().GetCode() ==
            ErrorCode::InvalidInput);
}

TEST_CASE("certification handles the edges of its input", "[surface][certify]")
{
    const Curve line{Line{Point3{1.0, 2.0, 0.0}, Vec3{3.0, 4.0, 0.0}}};
    const Curve2d onPlane{Line2d{Point2{1.0, 2.0}, Vec2{3.0, 4.0}}};
    // A range of one parameter: one span, the point's own enclosure.
    const auto point = CertifyDeviation(line, xy, onPlane, Interval{0.5, 0.5}, 1e-6);
    REQUIRE(point.has_value());
    REQUIRE(point->Spans == 1);
    REQUIRE(point->Bound <= 1e-14);
    // A range too narrow to split below a target under the rounding floor:
    // the failure names the span and the bound it reached.
    const auto floor = CertifyDeviation(line, xy, onPlane, Interval{0.5, 0.5 + 1e-15}, 1e-300);
    REQUIRE_FALSE(floor.has_value());
    REQUIRE(floor.error().GetCode() == ErrorCode::CertificationFailed);
    REQUIRE(std::string(floor.error().GetMessage()).find("after") != std::string::npos);
    // A pcurve off its curve fails on the first span, by a witness point.
    const Curve2d wrong{Line2d{Point2{1.0, 2.5}, Vec2{3.0, 4.0}}};
    const auto witness = CertifyDeviation(line, xy, wrong, Interval{0.0, 1.0}, 1e-6);
    REQUIRE_FALSE(witness.has_value());
    REQUIRE(std::string(witness.error().GetMessage()).find("at least") != std::string::npos);
}

TEST_CASE("B-splines of low degree are certified, rational ones refused, and the range checked",
          "[surface][certify]")
{
    std::vector<Point2> samples;
    for (int k = 0; k <= 6; ++k) {
        samples.emplace_back(1.0 + 3.0 * k / 6.0, 2.0 + 4.0 * k / 6.0);
    }
    const Curve line{Line{Point3{1.0, 2.0, 0.0}, Vec3{3.0, 4.0, 0.0}}};
    for (const int degree : {1, 2}) {
        const auto fit = Interpolate(samples, degree, Stereon::Parametrization::Uniform);
        REQUIRE(fit.has_value());
        const auto certificate =
            CertifyDeviation(line, xy, Curve2d{fit->Curve}, Interval{0.0, 1.0}, 1e-6);
        REQUIRE(certificate.has_value());
        REQUIRE(certificate->Bound <= 1e-6);
        // Outside the spline's domain: refused, never evaluated.
        REQUIRE(CertifyDeviation(line, xy, Curve2d{fit->Curve}, Interval{0.0, 1.5}, 1e-6)
                    .error()
                    .GetCode() == ErrorCode::InvalidInput);
    }
    const auto rational = Stereon::MakeBSplineCurve(
        1, {0.0, 0.0, 1.0, 1.0}, std::vector<Point2>{Point2{1.0, 2.0}, Point2{4.0, 6.0}},
        std::vector<double>{1.0, 2.0});
    REQUIRE(rational.has_value());
    REQUIRE(CertifyDeviation(line, xy, Curve2d{*rational}, Interval{0.0, 1.0}, 1e-6)
                .error()
                .GetCode() == ErrorCode::NotSupported);
}

namespace {

/// A biquadratic bump over the unit square with an interior knot in u, so
/// that a pcurve crosses a knot line; rational when asked.
Stereon::BSplineSurface Bump(const bool rational)
{
    std::vector<Stereon::Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            poles.push_back({i / 3.0, j / 2.0, (i == 1 || i == 2) && j == 1 ? 1.0 : 0.0});
            weights.push_back(j == 1 ? 2.0 : 1.0);
        }
    }
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knotsU, knotsV, poles,
                                                  rational ? weights : std::vector<double>{});
    REQUIRE(made.has_value());
    return *made;
}

/// The polyline of `count` chords of the surface along v = 1/2, and its
/// parameters, on a shared parameter u.
std::pair<Stereon::BSplineCurve, Stereon::BSplineCurve2d>
ChordsAlong(const Stereon::Surface& surface, const int count, const double offset)
{
    std::vector<Stereon::Point3> points;
    std::vector<Stereon::Point2> parameters;
    std::vector<double> knots{0.0};
    for (int k = 0; k <= count; ++k) {
        const double u = static_cast<double>(k) / count;
        const Stereon::Point3 p = Eval(surface, Stereon::Point2{u, 0.5});
        points.push_back({p.X, p.Y, p.Z + offset});
        parameters.push_back({u, 0.5});
        knots.push_back(u);
    }
    knots.push_back(1.0);
    const auto curve = Stereon::MakeBSplineCurve(1, knots, points);
    const auto pcurve = Stereon::MakeBSplineCurve(1, knots, parameters);
    REQUIRE(curve.has_value());
    REQUIRE(pcurve.has_value());
    return {*curve, *pcurve};
}

} // namespace

TEST_CASE("chords of a B-spline surface certify across a knot line", "[surface][certify]")
{
    for (const bool rational : {false, true}) {
        const Stereon::Surface surface{Bump(rational)};
        const auto [curve, pcurve] = ChordsAlong(surface, 1024, 0.0);
        const auto c =
            CertifyDeviation(Curve{curve}, surface, Curve2d{pcurve}, Interval{0.0, 1.0}, 1e-5);
        REQUIRE(c.has_value());
        REQUIRE(c->Bound <= 1e-5);
    }
}

TEST_CASE("a certificate on a B-spline surface refuses a curve off it", "[surface][certify]")
{
    for (const bool rational : {false, true}) {
        const Stereon::Surface surface{Bump(rational)};
        const auto [curve, pcurve] = ChordsAlong(surface, 256, 1e-3);
        const auto c =
            CertifyDeviation(Curve{curve}, surface, Curve2d{pcurve}, Interval{0.0, 1.0}, 1e-5);
        REQUIRE_FALSE(c.has_value());
        REQUIRE(c.error().GetCode() == ErrorCode::CertificationFailed);
    }
}

TEST_CASE("repeated knots do not turn unavailable derivatives into zero",
          "[surface][certify][regression]")
{
    const std::vector<double> ku{0, 0, 0, .5, .5, 1, 1, 1};
    const std::vector<double> kv{0, 0, 1, 1};
    const std::vector<double> z{0, 0, .25, .5, 1};
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < z.size(); ++i) {
        poles.push_back({static_cast<double>(i) * .25, 0, z[i]});
        poles.push_back({static_cast<double>(i) * .25, 1, z[i]});
    }
    const auto surface = Stereon::MakeBSplineSurface(2, 1, ku, kv, poles);
    REQUIRE(surface.has_value());
    const auto curve = Stereon::MakeBSplineCurve(
        1, {0, 0, 1, 1}, std::vector<Point3>{{.125, .5, 0}, {.375, .5, .125}});
    const auto pc =
        Stereon::MakeBSplineCurve(1, {0, 0, 1, 1}, std::vector<Point2>{{.125, .5}, {.375, .5}});
    REQUIRE(curve.has_value());
    REQUIRE(pc.has_value());
    const auto certificate =
        CertifyDeviation(Curve{*curve}, Surface{*surface}, Curve2d{*pc}, Interval{0, 1}, 1e-6);
    REQUIRE_FALSE(certificate.has_value());
}

namespace {

/// A surface whose x (alongU) or y (across) is exactly linear in one
/// parameter over knots a short span `h` apart, so that the partial along
/// it is one number, `slope`, at every parameter: the oracle is exact.
struct SteepPatch {
    Stereon::Surface Patch;
    bool AlongU;
    double Slope;
    int PartialDegrees; // the sum of the partial spline's two degrees
    double EndU;
    double EndV;
};

/// The poles are in proportion to the Greville abscissae, small dyadic
/// numbers for these knots, so every pole difference and knot difference is
/// exact in doubles: 3 u / h on one span, u / h with a knot at h.
[[nodiscard]] SteepPatch Steep(const bool interiorKnot, const bool alongU, const double h)
{
    const std::vector<double> shortKnots =
        interiorKnot ? std::vector<double>{0, 0, 0, h, 2 * h, 2 * h, 2 * h}
                     : std::vector<double>{0, 0, 0, 0, h, h, h, h};
    const std::vector<double> steps =
        interiorKnot ? std::vector<double>{0, 0.5, 1.5, 2} : std::vector<double>{0, 1, 2, 3};
    const std::vector<double> unitKnots{0, 0, 0, 0, 1, 1, 1, 1};
    const int shortDegree = interiorKnot ? 2 : 3;
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            poles.push_back(alongU ? Point3{steps[i], static_cast<double>(j), 0.0}
                                   : Point3{static_cast<double>(i), steps[j], 0.0});
        }
    }
    auto surface = alongU
                       ? Stereon::MakeBSplineSurface(shortDegree, 3, shortKnots, unitKnots, poles)
                       : Stereon::MakeBSplineSurface(3, shortDegree, unitKnots, shortKnots, poles);
    REQUIRE(surface.has_value());
    const double end = shortKnots.back();
    return {Stereon::Surface{*surface},
            alongU,
            (interiorKnot ? 1.0 : 3.0) / h,
            shortDegree - 1 + 3,
            alongU ? end : 1.0,
            alongU ? 1.0 : end};
}

} // namespace

// The margin of a partial derivative's enclosure was sized by the largest
// pole of the surface, not of the partial. Over a knot span h the partial's
// poles are about 1 / h times the surface's, so on a short span the margin
// was a small fraction of the rounding it is there to cover, and the
// enclosure held the exact value only by the outward steps of the interval
// sum, which no bound relates to the basis functions' rounding.
//
// The property is the one the margin is defined by: Detail::PoleBoxMargin
// bounds "how far Eval's rounding can move a computed point", relative to
// the largest pole coordinate of the spline evaluated, and SplineEnclosure
// keeps that magnitude per derivative curve. CertifyDeviation cannot show
// it: its bound meets this enclosure with the hull of the poles and adds
// the roundings of everything after, so the test reads the enclosure itself
// through Detail::HomogeneousJetAt.
//
// A2.2 reaches a basis function of degree p through p recurrence steps of
// at most six rounded operations each, all on non-negative numbers, so its
// relative error is within 6 p units of 2^-53 to first order, and that of a
// product of two within 6 (p + q). The partial here is `slope` exactly, and
// its poles all equal it, so its computed sum can be that far from it,
// relatively: the enclosure has to reach that far both ways.
TEST_CASE("a partial's enclosure covers the worst rounding of its own poles",
          "[surface][certify][regression]")
{
    // The span the defect was found on, 2^-20, among shorter and longer ones.
    for (const double h : {0x1p-1, 0x1p-10, 0x1p-20, 0x1p-30, 0x1p-44}) {
        for (const bool interiorKnot : {false, true}) {
            for (const bool alongU : {true, false}) {
                const SteepPatch patch = Steep(interiorKnot, alongU, h);
                const double reach = 6.0 * patch.PartialDegrees * 0x1p-53 * patch.Slope;
                CAPTURE(h, interiorKnot, alongU, patch.Slope, reach);
                for (int a = 1; a < 8; ++a) {
                    for (int b = 1; b < 8; ++b) {
                        const double u = patch.EndU * a / 8.0;
                        const double v = patch.EndV * b / 8.0;
                        const auto jet = Stereon::Detail::HomogeneousJetAt(patch.Patch, u, v);
                        REQUIRE(jet.has_value());
                        const Interval partial = alongU ? jet->AlongU[0] : jet->AlongV[1];
                        CAPTURE(u, v, partial.GetLower() - patch.Slope,
                                partial.GetUpper() - patch.Slope);
                        // Sound: the exact value is inside.
                        REQUIRE(Contains(partial, patch.Slope));
                        // Proved so: with room for the worst rounding.
                        REQUIRE(partial.GetLower() <= patch.Slope - reach);
                        REQUIRE(partial.GetUpper() >= patch.Slope + reach);
                        // And still a rounding margin, not a loose box.
                        REQUIRE(Width(partial) <= 0x1p-40 * patch.Slope);
                    }
                }
            }
        }
    }
}

TEST_CASE("a steep patch's own value is enclosed within its poles' margin", "[surface][certify]")
{
    // x = 3 u / h exactly, a dyadic number at these parameters: the value's
    // enclosure holds it, within the margin of the surface's own poles.
    const double h = 0x1p-20;
    const SteepPatch patch = Steep(false, true, h);
    for (int a = 0; a <= 8; ++a) {
        const double u = h * a / 8.0;
        const auto jet = Stereon::Detail::HomogeneousJetAt(patch.Patch, u, 0.5);
        REQUIRE(jet.has_value());
        const double x = 3.0 * a / 8.0;
        REQUIRE(Contains(jet->Value[0], x));
        REQUIRE(Contains(jet->Value[3], 1.0));
        REQUIRE(Width(jet->Value[0]) <= 0x1p-40);
    }
}

// --- An exact oracle for the homogeneous jet --------------------------------
//
// The splines below have integer poles and weights, and knots 0, 1 and 2
// only, so every knot difference is 1 or 2 and every basis function at a
// parameter that is a multiple of 1/4 is a dyadic number of a few bits:
// Cox-de Boor's recurrence and the sums over the poles are then exact in
// doubles (at most 8 fractional bits a direction at degree 4, on poles
// below 2^18). Scaling the knots by a power of two scales the partials by
// its inverse, exactly, so the same oracle answers for knot spans of any
// size.

namespace {

using Pole4 = std::array<double, 4>;

/// A tensor-product spline of homogeneous poles (w x, w y, w z, w) on
/// normalised knots.
struct Net {
    int DegreeU = 0;
    int DegreeV = 0;
    std::vector<double> KnotsU;
    std::vector<double> KnotsV;
    std::size_t CountU = 0;
    std::size_t CountV = 0;
    std::vector<Pole4> Poles;
};

/// Clamped knots of the degree on [0, 1], or on [0, 2] with a knot at 1 of
/// the given multiplicity.
[[nodiscard]] std::vector<double> NormalisedKnots(const int degree, const int interior)
{
    std::vector<double> knots(static_cast<std::size_t>(degree + 1), 0.0);
    knots.insert(knots.end(), static_cast<std::size_t>(interior), 1.0);
    knots.insert(knots.end(), static_cast<std::size_t>(degree + 1), interior == 0 ? 1.0 : 2.0);
    return knots;
}

/// The span FindSpan takes: knots[span] <= s < knots[span + 1], the last
/// one that is not empty at the end of the domain.
[[nodiscard]] std::size_t OracleSpan(const std::vector<double>& knots, const std::size_t count,
                                     const double s)
{
    std::size_t span = count - 1;
    while (knots[span] > s || knots[span] == knots[span + 1]) {
        --span;
    }
    return span;
}

/// N_{i,d}(s) on the span, by Cox-de Boor's recurrence from the definition;
/// a term over an empty support is zero.
[[nodiscard]] double OracleBasis(const std::vector<double>& k, const std::size_t i, const int d,
                                 const std::size_t span, const double s)
{
    if (d == 0) {
        return i == span ? 1.0 : 0.0;
    }
    const auto step = static_cast<std::size_t>(d);
    double value = 0.0;
    if (k[i + step] > k[i]) {
        value += (s - k[i]) / (k[i + step] - k[i]) * OracleBasis(k, i, d - 1, span, s);
    }
    if (k[i + step + 1] > k[i + 1]) {
        value += (k[i + step + 1] - s) / (k[i + step + 1] - k[i + 1]) *
                 OracleBasis(k, i + 1, d - 1, span, s);
    }
    return value;
}

/// The partial along one direction (The NURBS Book, eq. 3.8, along each
/// line of poles).
[[nodiscard]] Net PartialAlong(const Net& s, const bool alongU)
{
    const auto d = static_cast<std::size_t>(alongU ? s.DegreeU : s.DegreeV);
    const std::vector<double>& knots = alongU ? s.KnotsU : s.KnotsV;
    Net t = s;
    (alongU ? t.DegreeU : t.DegreeV) -= 1;
    (alongU ? t.CountU : t.CountV) -= 1;
    (alongU ? t.KnotsU : t.KnotsV).assign(knots.begin() + 1, knots.end() - 1);
    t.Poles.clear();
    for (std::size_t i = 0; i < t.CountU; ++i) {
        for (std::size_t j = 0; j < t.CountV; ++j) {
            const std::size_t k = alongU ? i : j;
            const double scale = static_cast<double>(d) / (knots[k + d + 1] - knots[k + 1]);
            const Pole4& low = s.Poles[i * s.CountV + j];
            const Pole4& high = s.Poles[(alongU ? i + 1 : i) * s.CountV + (alongU ? j : j + 1)];
            Pole4 pole{};
            for (std::size_t c = 0; c < pole.size(); ++c) {
                pole[c] = (high[c] - low[c]) * scale;
            }
            t.Poles.push_back(pole);
        }
    }
    return t;
}

/// What a spline is at one parameter, per homogeneous component: its value,
/// the sum of the magnitudes of the terms that make it (the most a relative
/// error of the basis functions can weigh on), and the largest magnitude
/// among the poles that take part.
struct ExactSum {
    Pole4 Value{};
    Pole4 Mass{};
    Pole4 Largest{};
};

[[nodiscard]] ExactSum SumAt(const Net& s, const double u, const double v)
{
    const std::size_t spanU = OracleSpan(s.KnotsU, s.CountU, u);
    const std::size_t spanV = OracleSpan(s.KnotsV, s.CountV, v);
    const auto p = static_cast<std::size_t>(s.DegreeU);
    const auto q = static_cast<std::size_t>(s.DegreeV);
    ExactSum sum;
    for (std::size_t i = spanU - p; i <= spanU; ++i) {
        for (std::size_t j = spanV - q; j <= spanV; ++j) {
            const double basis = OracleBasis(s.KnotsU, i, s.DegreeU, spanU, u) *
                                 OracleBasis(s.KnotsV, j, s.DegreeV, spanV, v);
            const Pole4& pole = s.Poles[i * s.CountV + j];
            for (std::size_t c = 0; c < pole.size(); ++c) {
                sum.Value[c] += basis * pole[c];
                sum.Mass[c] += basis * std::abs(pole[c]);
                sum.Largest[c] = std::max(sum.Largest[c], std::abs(pole[c]));
            }
        }
    }
    return sum;
}

/// A surface on the normalised knots scaled by spanU and spanV, powers of
/// two, with the oracle's net of its homogeneous poles.
struct OraclePatch {
    Stereon::Surface Patch;
    Net Homogeneous;
    double SpanU;
    double SpanV;
};

/// Integer poles from a fixed sequence: x and y of either sign below 2^8, z
/// about 2^12 so that one component's poles dwarf the others', and, when
/// rational, weights from 1 to 4.
[[nodiscard]] OraclePatch MakeOraclePatch(const int degreeU, const int degreeV, const int interiorU,
                                          const int interiorV, const bool rational,
                                          const double spanU, const double spanV)
{
    Net net{degreeU,
            degreeV,
            NormalisedKnots(degreeU, interiorU),
            NormalisedKnots(degreeV, interiorV),
            0,
            0,
            {}};
    net.CountU = net.KnotsU.size() - static_cast<std::size_t>(degreeU) - 1;
    net.CountV = net.KnotsV.size() - static_cast<std::size_t>(degreeV) - 1;
    std::uint32_t state = 0x9E3779B9U + static_cast<std::uint32_t>(degreeU * 64 + degreeV * 16 +
                                                                   interiorU * 4 + interiorV);
    const auto next = [&state](const std::uint32_t below) {
        state = state * 1664525U + 1013904223U;
        return static_cast<double>((state >> 8U) % below);
    };
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < net.CountU * net.CountV; ++i) {
        const Point3 pole{next(513) - 256.0, next(513) - 256.0, 4096.0 + next(65)};
        const double w = rational ? 1.0 + next(4) : 1.0;
        poles.push_back(pole);
        if (rational) {
            weights.push_back(w);
        }
        net.Poles.push_back({w * pole.X, w * pole.Y, w * pole.Z, w});
    }
    std::vector<double> knotsU = net.KnotsU;
    std::vector<double> knotsV = net.KnotsV;
    for (double& k : knotsU) {
        k *= spanU;
    }
    for (double& k : knotsV) {
        k *= spanV;
    }
    auto surface =
        Stereon::MakeBSplineSurface(degreeU, degreeV, std::move(knotsU), std::move(knotsV),
                                    std::move(poles), std::move(weights));
    REQUIRE(surface.has_value());
    return {Stereon::Surface{*surface}, std::move(net), spanU, spanV};
}

/// What a sweep found: how many enclosed components it compared, how many
/// did not cover the exact value with the rounding the basis functions may
/// have, how many were wider than a rounding margin, and the first of each.
struct Tally {
    int Compared = 0;
    int Uncovered = 0;
    int Wide = 0;
    std::string FirstUncovered;
    std::string FirstWide;
};

/// Compare the components of one enclosed spline with the oracle's. `scale`
/// takes the oracle's normalised parameters to the surface's: 1 for the
/// value, the inverse of the knot span for a partial.
void Compare(Tally& tally, const std::array<Interval, 4>& enclosed, const Net& net, const double s,
             const double t, const double scale, const std::string& where)
{
    const ExactSum exact = SumAt(net, s, t);
    // Six rounded operations a recurrence step, a step a degree, in each
    // direction: the relative error of a product of two basis functions.
    const double rounding = 6.0 * (net.DegreeU + net.DegreeV) * 0x1p-53;
    // Narrow on the spline's largest pole of any component, which is all a
    // margin has to be. Poles that are exactly zero are enclosed within a
    // subnormal.
    const double widest = 0x1p-40 * *std::ranges::max_element(exact.Largest) * scale + 0x1p-900;
    for (std::size_t c = 0; c < enclosed.size(); ++c) {
        const double value = exact.Value[c] * scale;
        const double reach = rounding * exact.Mass[c] * scale;
        const auto describe = [&] {
            return where + ", component " + std::to_string(c) + ": [" +
                   std::to_string((enclosed[c].GetLower() - value) / 0x1p-53) + ", " +
                   std::to_string((enclosed[c].GetUpper() - value) / 0x1p-53) +
                   "] units of 2^-53 about the exact value, to reach " +
                   std::to_string(reach / 0x1p-53) + " and stay within " +
                   std::to_string(widest / 0x1p-53);
        };
        ++tally.Compared;
        if (!(enclosed[c].GetLower() <= value - reach && enclosed[c].GetUpper() >= value + reach)) {
            if (tally.Uncovered++ == 0) {
                tally.FirstUncovered = describe();
            }
        }
        if (!(Width(enclosed[c]) <= widest)) {
            if (tally.Wide++ == 0) {
                tally.FirstWide = describe();
            }
        }
    }
}

} // namespace

// The whole of what the mean-value forms start from: the value and both
// first partials, all four homogeneous components, of polynomial and
// rational surfaces of degrees 1 to 4 each way (so partials of degree 0 to
// 3), Bezier and with a knot of multiplicity 1 or the degree, over knot
// spans from 1 down to 2^-45 in either direction or both, at span ends,
// span interiors and on the knot. Each enclosure has to hold the exact
// value with room either way for the worst rounding of the basis functions
// on the terms that make it, and stay as narrow as a rounding margin on its
// own poles.
//
// The second partials are not here: the certificate never evaluates them at
// a point, only takes the hull of their poles, which are interval
// enclosures and need no margin.
TEST_CASE("a B-spline surface's jet covers its rounding on every span, polynomial or rational",
          "[surface][certify][regression]")
{
    struct Spans {
        double U;
        double V;
    };
    const std::array<Spans, 6> spans{{{1.0, 1.0},
                                      {0x1p-12, 1.0},
                                      {1.0, 0x1p-30},
                                      {0x1p-30, 0x1p-12},
                                      {0x1p-45, 0x1p-45},
                                      {0x1p+10, 0x1p-20}}};
    Tally tally;
    for (int degreeU = 1; degreeU <= 4; ++degreeU) {
        for (int degreeV = 1; degreeV <= 4; ++degreeV) {
            for (const int interiorU : {0, 1, degreeU}) {
                for (const int interiorV : {0, 1, degreeV}) {
                    for (const bool rational : {false, true}) {
                        for (const Spans& span : spans) {
                            const OraclePatch patch = MakeOraclePatch(
                                degreeU, degreeV, interiorU, interiorV, rational, span.U, span.V);
                            const Net alongU = PartialAlong(patch.Homogeneous, true);
                            const Net alongV = PartialAlong(patch.Homogeneous, false);
                            const double endU = interiorU == 0 ? 1.0 : 2.0;
                            const double endV = interiorV == 0 ? 1.0 : 2.0;
                            const std::string which =
                                "degrees " + std::to_string(degreeU) + " by " +
                                std::to_string(degreeV) + ", interior knots " +
                                std::to_string(interiorU) + " and " + std::to_string(interiorV) +
                                (rational ? ", rational" : ", polynomial") + ", spans 2^" +
                                std::to_string(std::ilogb(span.U)) + " by 2^" +
                                std::to_string(std::ilogb(span.V));
                            for (int a = 0; a <= 4; ++a) {
                                for (int b = 0; b <= 4; ++b) {
                                    const double s = endU * a / 4.0;
                                    const double t = endV * b / 4.0;
                                    const auto jet = Stereon::Detail::HomogeneousJetAt(
                                        patch.Patch, s * span.U, t * span.V);
                                    if (!jet.has_value()) {
                                        FAIL(which << ": no jet");
                                    }
                                    const std::string at = which + ", at (" + std::to_string(s) +
                                                           ", " + std::to_string(t) + ")";
                                    Compare(tally, jet->Value, patch.Homogeneous, s, t, 1.0,
                                            at + ", the value");
                                    Compare(tally, jet->AlongU, alongU, s, t, 1.0 / span.U,
                                            at + ", the partial along u");
                                    Compare(tally, jet->AlongV, alongV, s, t, 1.0 / span.V,
                                            at + ", the partial along v");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // 4 x 4 degrees, 3 x 3 knot vectors, 2 kinds, 6 spans, 5 x 5 parameters;
    // 3 splines of 4 components each.
    REQUIRE(tally.Compared == 16 * 9 * 2 * 6 * 25 * 12);
    CAPTURE(tally.Uncovered, tally.FirstUncovered, tally.Wide, tally.FirstWide);
    CHECK(tally.Uncovered == 0);
    CHECK(tally.Wide == 0);
}

// Poles of both signs that cancel: x = (2 u / h - 1)^2 over one span has
// the partial 4 (2 u / h - 1) / h, zero at the middle of the span as the
// difference of two poles 4 / h in magnitude. The rounding the enclosure
// has to cover there is that of those poles, not of the value it encloses.
TEST_CASE("a partial that vanishes by cancellation is enclosed to its poles' rounding",
          "[surface][certify][regression]")
{
    for (const double h : {1.0, 0x1p-10, 0x1p-20, 0x1p-40}) {
        const std::vector<double> shortKnots{0, 0, 0, h, h, h};
        const std::vector<double> unitKnots{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
        std::vector<Point3> poles;
        for (const double x : {1.0, -1.0, 1.0}) {
            for (int j = 0; j <= 5; ++j) {
                poles.push_back(Point3{x, static_cast<double>(j), 0.0});
            }
        }
        const auto surface = Stereon::MakeBSplineSurface(2, 5, shortKnots, unitKnots, poles);
        REQUIRE(surface.has_value());
        const auto jet = Stereon::Detail::HomogeneousJetAt(Stereon::Surface{*surface}, h / 2, 0.5);
        REQUIRE(jet.has_value());
        // The partial is of degree 1 by 5, on the poles -4 / h and 4 / h.
        const double reach = 6.0 * 6.0 * 0x1p-53 * (4.0 / h);
        CAPTURE(h, reach, jet->AlongU[0].GetLower(), jet->AlongU[0].GetUpper());
        REQUIRE(jet->AlongU[0].GetLower() <= -reach);
        REQUIRE(jet->AlongU[0].GetUpper() >= reach);
        REQUIRE(Width(jet->AlongU[0]) <= 0x1p-40 * (4.0 / h));
    }
}
