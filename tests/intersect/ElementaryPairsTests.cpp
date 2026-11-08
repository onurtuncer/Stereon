// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The closed forms of the analytic table's remaining pairs (ADR-0010,
// ADR-0029): coaxial surfaces of revolution meeting in exact circles,
// touching along them or at a point, and coinciding; a plane through a
// cone's apex; a plane containing a torus's axis. The configurations a
// fillet makes, a torus tangent to a cylinder and a plane, are the ones
// the trace alone could not resolve.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Catch::Matchers::WithinAbs;
using IntersectTest::ctx;
using Stereon::Box2;
using Stereon::Circle;
using Stereon::Cone;
using Stereon::ContactType;
using Stereon::Interval;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::SideTransition;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
const Vec3 X{1.0, 0.0, 0.0};
const Vec3 Y{0.0, 1.0, 0.0};
const Vec3 Z{0.0, 0.0, 1.0};

const Surface torus{Torus{Point3{}, 3.0, 1.0, X, Y, true}};

Surface Horizontal(const double z)
{
    return Surface{Plane{Point3{0.0, 0.0, z}, X, Y}};
}

Surface CylinderZ(const double r)
{
    return Surface{IntersectTest::CylinderAbout(Point3{}, Z, r)};
}

/// The circles of a result: every branch an exact circle, certified, its
/// pcurves agreeing with it; their radii and heights in order.
struct CircleAt {
    double Radius;
    double Z;
};

std::vector<CircleAt> CirclesOf(const Stereon::SurfaceSurfaceIntersection& result, const Surface& a,
                                const Surface& b, const ContactType contact)
{
    std::vector<CircleAt> circles;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        const Circle* c = branch.Curve3d.Target<Circle>();
        REQUIRE(c != nullptr);
        REQUIRE(branch.Contact == contact);
        REQUIRE(branch.Transition == (contact == ContactType::Tangential
                                          ? SideTransition::SameSide
                                          : SideTransition::OppositeSide));
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 100) <= 1e-12);
        const bool seen = std::ranges::any_of(circles, [&](const CircleAt& k) {
            return std::abs(k.Radius - c->Radius) < 1e-12 && std::abs(k.Z - c->Center.Z) < 1e-12;
        });
        if (!seen) {
            circles.push_back({c->Radius, c->Center.Z});
        }
    }
    return circles;
}

/// The total angle the branches of one circle cover: a whole turn when the
/// boxes restrict nothing, in pieces cut at the charts' seams.
double TurnOf(const Stereon::SurfaceSurfaceIntersection& result)
{
    double total = 0.0;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        total += branch.Range.GetUpper() - branch.Range.GetLower();
    }
    return total;
}

} // namespace

TEST_CASE("a fillet torus touches the cylinder round it along one circle",
          "[intersect][elementary]")
{
    for (const double r : {4.0, 2.0}) {
        const Surface cylinder = CylinderZ(r);
        const auto result = Intersect(torus, cylinder, ctx);
        INFO("r = " << r << ": " << IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        const auto circles = CirclesOf(*result, torus, cylinder, ContactType::Tangential);
        REQUIRE(circles.size() == 1);
        REQUIRE_THAT(circles[0].Radius, WithinAbs(r, 1e-14));
        REQUIRE_THAT(circles[0].Z, WithinAbs(0.0, 1e-14));
        REQUIRE_THAT(TurnOf(*result), WithinAbs(2.0 * pi, 1e-12));
        // The cylinder first: the same circle.
        REQUIRE(Intersect(cylinder, torus, ctx)->Branches.front().Contact ==
                ContactType::Tangential);
    }
}

TEST_CASE("a cylinder through the tube cuts two circles", "[intersect][elementary]")
{
    const Surface cylinder = CylinderZ(3.5);
    const auto result = Intersect(cylinder, torus, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    const auto circles = CirclesOf(*result, cylinder, torus, ContactType::Crossing);
    REQUIRE(circles.size() == 2);
    REQUIRE_THAT(circles[0].Z, WithinAbs(-std::sqrt(0.75), 1e-14));
    REQUIRE_THAT(circles[1].Z, WithinAbs(std::sqrt(0.75), 1e-14));
}

TEST_CASE("a plane across the axis touches the torus's top or cuts two circles",
          "[intersect][elementary]")
{
    const auto top = Intersect(Horizontal(1.0), torus, ctx);
    REQUIRE(top.has_value());
    const auto touch = CirclesOf(*top, Horizontal(1.0), torus, ContactType::Tangential);
    REQUIRE(touch.size() == 1);
    REQUIRE_THAT(touch[0].Radius, WithinAbs(3.0, 1e-14));
    REQUIRE(top->Branches.front().OnA.Target<Stereon::Circle2d>() != nullptr);

    const auto cut = Intersect(torus, Horizontal(0.5), ctx);
    REQUIRE(cut.has_value());
    const auto circles = CirclesOf(*cut, torus, Horizontal(0.5), ContactType::Crossing);
    REQUIRE(circles.size() == 2);
    REQUIRE(Intersect(Horizontal(1.5), torus, ctx)->IsEmpty());
}

TEST_CASE("a sphere seated in a cone touches it along a circle", "[intersect][elementary]")
{
    // The apex at the origin, 30 degrees: a sphere of radius 2 centred 4 up
    // the axis touches every generator at distance 4 sin 30 = 2.
    const Surface cone{Cone{Point3{}, 0.0, pi / 6.0, X, Y, true}};
    const Surface sphere{Sphere{Point3{0.0, 0.0, 4.0}, 2.0, X, Y, true}};
    const auto result = Intersect(cone, sphere, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    const auto circles = CirclesOf(*result, cone, sphere, ContactType::Tangential);
    REQUIRE(circles.size() == 1);
    REQUIRE_THAT(circles[0].Radius, WithinAbs(std::sqrt(3.0), 1e-12));
    REQUIRE_THAT(circles[0].Z, WithinAbs(3.0, 1e-12));
}

TEST_CASE("a sphere through a cone's apex meets it there and in a circle",
          "[intersect][elementary]")
{
    const Surface cone{Cone{Point3{}, 0.0, pi / 6.0, X, Y, true}};
    const Surface sphere{Sphere{Point3{0.0, 0.0, 2.0}, 2.0, X, Y, true}};
    const auto result = Intersect(sphere, cone, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Points.size() == 1);
    REQUIRE(Distance(result->Points.front().Position, Point3{}) <= 1e-12);
    REQUIRE(CirclesOf(*result, sphere, cone, ContactType::Crossing).size() == 1);
}

TEST_CASE("spheres meet in a circle, touch at a point, or coincide", "[intersect][elementary]")
{
    const Surface a{Sphere{Point3{}, 2.0, X, Y, true}};
    const Surface crossing{Sphere{Point3{0.0, 0.0, 3.0}, 2.0, X, Y, true}};
    const auto circle = Intersect(a, crossing, ctx);
    INFO(IntersectTest::ErrorText(circle));
    REQUIRE(circle.has_value());
    const auto circles = CirclesOf(*circle, a, crossing, ContactType::Crossing);
    REQUIRE(circles.size() == 1);
    REQUIRE_THAT(circles[0].Radius, WithinAbs(std::sqrt(4.0 - 2.25), 1e-14));
    REQUIRE_THAT(circles[0].Z, WithinAbs(1.5, 1e-14));

    const Surface touching{Sphere{Point3{0.0, 0.0, 4.0}, 2.0, X, Y, true}};
    const auto point = Intersect(a, touching, ctx);
    REQUIRE(point.has_value());
    REQUIRE(point->Branches.empty());
    REQUIRE(point->Points.size() == 1);
    REQUIRE(Distance(point->Points.front().Position, Point3{0.0, 0.0, 2.0}) <= 1e-14);

    REQUIRE(Intersect(a, Surface{Sphere{Point3{0.0, 0.0, 5.0}, 2.0, X, Y, true}}, ctx)->IsEmpty());
    REQUIRE(Intersect(a, Surface{Sphere{Point3{0.1, 0.0, 0.0}, 0.5, X, Y, true}}, ctx)->IsEmpty());

    const auto same = Intersect(a, Surface{Sphere{Point3{}, 2.0, X, Y, true}}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Overlaps.size() == 1);
    REQUIRE(same->Overlaps.front().FullDomain);
    REQUIRE(same->Overlaps.front().Evidence == Stereon::CoincidenceEvidence::Exact);
    REQUIRE(same->Overlaps.front().Distance == 0.0);
}

TEST_CASE("a tangent circle declared within the tolerance is certified within it",
          "[intersect][elementary]")
{
    // Regression (corpus analytic-ssi/sweep-sphere-cylinder-coaxial): a
    // sphere centred on a cylinder's axis, the radii 5e-4 apart, under a
    // pair tolerance of 1e-3, is declared tangent along a circle (ADR-0029
    // section 6); that circle lies up to 5e-4 from one surface, and its
    // certification against the resolution failed the whole intersection.
    for (const double radius : {2.0 - 5e-4, 2.0 + 5e-4}) {
        const auto touch = Intersect(Surface{Sphere{Point3{}, 2.0, X, Y, true}}, Box2::Everything(),
                                     CylinderZ(radius), Box2::Everything(), 1e-3, ctx);
        INFO(IntersectTest::ErrorText(touch));
        REQUIRE(touch.has_value());
        // One circle, in arcs cut at the surfaces' seams.
        REQUIRE(IntersectTest::ConnectivityOf(*touch).Curves == 1);
        for (const Stereon::IntersectionBranch& arc : touch->Branches) {
            REQUIRE(arc.Contact == ContactType::Tangential);
            REQUIRE(arc.IsCertified());
            // The bound is what the circle is from the surfaces: about the
            // gap, within the tolerance.
            REQUIRE(arc.Bound >= 4e-4);
            REQUIRE(arc.Bound <= 1e-3);
        }
    }
}

TEST_CASE("coaxial surfaces of revolution within the tolerance are declared, not exact",
          "[intersect][elementary]")
{
    // Regression (corpus analytic-ssi/sphere-sphere/concentric-near-equal):
    // concentric spheres whose radii differ by less than the pair tolerance
    // were declared coincident with Exact evidence and no distance. ADR-0029
    // section 6: a coincidence declared within the tolerance says so, and
    // how far apart the surfaces are, as for parallel planes and cylinders.
    const Surface a{Sphere{Point3{}, 2.0, X, Y, true}};
    const Surface b{Sphere{Point3{}, 2.0005, X, Y, true}};
    const auto declared = Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-3, ctx);
    REQUIRE(declared.has_value());
    REQUIRE(declared->Overlaps.size() == 1);
    const Stereon::OverlapRegion& overlap = declared->Overlaps.front();
    REQUIRE(overlap.Evidence == Stereon::CoincidenceEvidence::Tolerant);
    REQUIRE_THAT(overlap.Distance, WithinAbs(0.0005, 1e-15));
    // Beyond the tolerance they are two spheres apart, one inside the other.
    REQUIRE(Intersect(a, Box2::Everything(), b, Box2::Everything(), 1e-4, ctx)->IsEmpty());
    // Coaxial cylinders 4e-4 apart in radius, within 1e-3, which the
    // cylinder pair's own handler declares: the same evidence.
    const auto wide = Intersect(CylinderZ(2.0), Box2::Everything(), CylinderZ(2.0004),
                                Box2::Everything(), 1e-3, ctx);
    REQUIRE(wide.has_value());
    REQUIRE(wide->Overlaps.size() == 1);
    REQUIRE(wide->Overlaps.front().Evidence == Stereon::CoincidenceEvidence::Tolerant);
    REQUIRE_THAT(wide->Overlaps.front().Distance, WithinAbs(0.0004, 1e-15));
}

TEST_CASE("spheres whose poles are off their centres' line go to the trace",
          "[intersect][elementary]")
{
    const Surface a{Sphere{Point3{}, 2.0, X, Y, true}};
    const Surface b{Sphere{Point3{3.0, 0.0, 0.0}, 2.0, X, Y, true}};
    const auto result = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_FALSE(result->Branches.empty());
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 200) <= ctx.GetResolution());
    }
    REQUIRE(IntersectTest::ConnectivityOf(*result).Curves == 1);
}

TEST_CASE("a plane touches a sphere at a point", "[intersect][elementary]")
{
    const Surface sphere{Sphere{Point3{}, 2.0, X, Y, true}};
    const Vec3 toward = 2.0 * Normalized(Vec3{1.0, 2.0, 2.0});
    const Surface oblique{IntersectTest::PlaneWithNormal(Point3{} + toward, Vec3{1.0, 2.0, 2.0})};
    const auto result = Intersect(oblique, sphere, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.empty());
    REQUIRE(result->Points.size() == 1);
}

TEST_CASE("a plane through a cone's apex cuts two generators, touches one, or only the apex",
          "[intersect][elementary]")
{
    const Surface cone{Cone{Point3{}, 0.0, pi / 6.0, X, Y, true}};
    // Containing the axis: two generators.
    const Surface vertical{Plane{Point3{}, X, Z}};
    const auto two = Intersect(vertical, cone, ctx);
    INFO(IntersectTest::ErrorText(two));
    REQUIRE(two.has_value());
    REQUIRE(two->Branches.size() == 2);
    for (const Stereon::IntersectionBranch& branch : two->Branches) {
        REQUIRE(branch.Curve3d.Target<Stereon::Line>() != nullptr);
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(IntersectTest::PcurveDeviation(branch, vertical, cone, 50) <= 1e-12);
    }
    // Tilted by 60 degrees from the axis: tangent along one generator.
    const Surface tangent{IntersectTest::PlaneWithNormal(
        Point3{}, Vec3{std::cos(pi / 6.0), 0.0, -std::sin(pi / 6.0)})};
    const auto one = Intersect(tangent, cone, ctx);
    INFO(IntersectTest::ErrorText(one));
    REQUIRE(one.has_value());
    REQUIRE(one->Branches.size() == 1);
    REQUIRE(one->Branches.front().Contact == ContactType::Tangential);
    // Across the axis: the apex alone.
    const auto apex = Intersect(Horizontal(0.0), cone, ctx);
    REQUIRE(apex.has_value());
    REQUIRE(apex->Branches.empty());
    REQUIRE(apex->Points.size() == 1);
}

TEST_CASE("a plane containing the torus's axis cuts two meridian circles",
          "[intersect][elementary]")
{
    const Surface meridian{Plane{Point3{}, X, Z}};
    const auto result = Intersect(meridian, torus, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 2);
    for (const Stereon::IntersectionBranch& branch : result->Branches) {
        const Circle* c = branch.Curve3d.Target<Circle>();
        REQUIRE(c != nullptr);
        REQUIRE(c->Radius == 1.0);
        REQUIRE_THAT(std::abs(c->Center.X), WithinAbs(3.0, 1e-14));
        REQUIRE(branch.IsCertified());
        REQUIRE(IntersectTest::PcurveDeviation(branch, meridian, torus, 100) <= 1e-12);
    }
}

TEST_CASE("coaxial cones cross in circles and coincide when alike", "[intersect][elementary]")
{
    const Surface a{Cone{Point3{}, 0.0, pi / 6.0, X, Y, true}};
    const Surface b{Cone{Point3{0.0, 0.0, 1.0}, 0.0, pi / 4.0, X, Y, true}};
    const auto result = Intersect(a, b, ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    // rho = z tan 30 = |z - 1| tan 45: at z below and above the second apex.
    const auto circles = CirclesOf(*result, a, b, ContactType::Crossing);
    REQUIRE(circles.size() == 2);

    const auto same = Intersect(a, Surface{Cone{Point3{}, 0.0, pi / 6.0, X, Y, true}}, ctx);
    REQUIRE(same.has_value());
    REQUIRE(same->Overlaps.size() == 1);
}

TEST_CASE("a sphere off the axis touches a cylinder, a cone or a torus at a point",
          "[intersect][elementary]")
{
    // Against the wall of a cylinder of radius 3, outside and inside.
    const Surface cylinder = CylinderZ(3.0);
    for (const auto& [x, r] : {std::pair{4.0, 1.0}, std::pair{2.0, 1.0}}) {
        const Surface ball{Sphere{Point3{x, 0.0, 1.0}, r, X, Y, true}};
        const auto result = Intersect(ball, cylinder, ctx);
        INFO("x = " << x << ": " << IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        REQUIRE(result->Branches.empty());
        REQUIRE(result->Points.size() == 1);
        REQUIRE(Distance(result->Points.front().Position, Point3{3.0, 0.0, 1.0}) <= 1e-14);
    }
    // Clear of the wall: proved separate.
    REQUIRE(Intersect(Surface{Sphere{Point3{4.5, 0.0, 1.0}, 1.0, X, Y, true}}, cylinder, ctx)
                ->IsEmpty());
    // On the tube of the torus, from outside it, beside the axis.
    const Surface onTube{Sphere{Point3{0.0, 5.5, 0.0}, 1.5, X, Y, true}};
    const auto tube = Intersect(torus, onTube, ctx);
    INFO(IntersectTest::ErrorText(tube));
    REQUIRE(tube.has_value());
    REQUIRE(tube->Points.size() == 1);
    REQUIRE(Distance(tube->Points.front().Position, Point3{0.0, 4.0, 0.0}) <= 1e-14);
    // Against a generator of a cone, the apex at the origin.
    const Surface cone{Cone{Point3{}, 0.0, pi / 4.0, X, Y, true}};
    const Point3 foot{2.0, 0.0, 2.0};
    const Surface ball{Sphere{foot + 0.5 * Normalized(Vec3{1.0, 0.0, -1.0}), 0.5, X, Y, true}};
    const auto touch = Intersect(cone, ball, ctx);
    INFO(IntersectTest::ErrorText(touch));
    REQUIRE(touch.has_value());
    REQUIRE(touch->Points.size() == 1);
    REQUIRE(Distance(touch->Points.front().Position, foot) <= 1e-14);
}

TEST_CASE("boxes clip the circles of coaxial surfaces", "[intersect][elementary]")
{
    const Surface cylinder = CylinderZ(3.5);
    const Box2 half{Interval{0.0, pi}, Interval::Everything()};
    const auto result = Intersect(cylinder, half, torus, Box2::Everything(), ctx);
    INFO(IntersectTest::ErrorText(result));
    REQUIRE(result.has_value());
    REQUIRE_THAT(TurnOf(*result), WithinAbs(2.0 * pi, 1e-12)); // Half of each of two circles.
    const Box2 upper{Interval{0.0, 2.0 * pi}, Interval{0.0, 5.0}};
    const auto above = Intersect(cylinder, upper, torus, Box2::Everything(), ctx);
    REQUIRE(above.has_value());
    REQUIRE(CirclesOf(*above, cylinder, torus, ContactType::Crossing).size() == 1);
}
