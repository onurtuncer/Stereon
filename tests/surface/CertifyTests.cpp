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

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
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
    REQUIRE(round->Bound <= 1e-5);
    REQUIRE(round->Spans > 1); // The bound is the enclosures', so near the target, not tiny. // The
                               // spans end up finer than the target needs.
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

TEST_CASE("certification refuses what it cannot bound", "[surface][certify]")
{
    const Curve line{Line{Point3{}, Vec3{1.0, 0.0, 0.0}}};
    const Curve2d onPlane{Line2d{Point2{}, Vec2{1.0, 0.0}}};
    const Surface sphere{
        Stereon::Sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
    REQUIRE(CertifyDeviation(line, sphere, onPlane, Interval{0.0, 1.0}, 1e-6).error().GetCode() ==
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
