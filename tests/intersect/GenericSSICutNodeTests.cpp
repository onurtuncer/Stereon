// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Transverse nodes on a cut of the requested domain, or on an end of a
// surface's domain, in the generic SSI (ADR-0010): only the half-branches
// inside the searched domain are answered, each from the node on the cut to
// its own end. Each answer is held to a closed form written here: the node,
// the ends, the curve every sample lies on, every pcurve inside the request.
// The same answer, moved accordingly, under exchange of the surfaces,
// reversal and transposition of parameters (the request moved with them), a
// far translation and scaling by powers of two (ADR-0053, MR-01 to MR-04).
// A branch running along the cut, and a node on cuts of both surfaces, fail
// closed.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Context;
using Stereon::ErrorCode;
using Stereon::Interval;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();
constexpr double tolerance = 1e-6;
constexpr double Pi = std::numbers::pi;

/// z = alpha x'^2 + beta y'^2, x' = x - 1/2, over x in [x0, x1] and y in
/// [y0, y1], biquadratic, x = u and y = v: the Bernstein coefficients of
/// (a + w t)^2 over t in [0, 1] are a^2, a (a + w) and (a + w)^2.
BSplineSurface Quadric(const double alpha, const double beta, const Box2& over)
{
    const auto squares = [](const Interval& range) {
        const double a = range.GetLower() - 0.5;
        const double b = range.GetUpper() - 0.5;
        return std::array<double, 3>{a * a, a * b, b * b};
    };
    const std::array<double, 3> sx = squares(over.X);
    const std::array<double, 3> sy = squares(over.Y);
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const double x = over.X.GetLower() + 0.5 * static_cast<double>(i) * Width(over.X);
            const double y = over.Y.GetLower() + 0.5 * static_cast<double>(j) * Width(over.Y);
            poles.push_back({x, y, alpha * sx.at(i) + beta * sy.at(j)});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

const Box2 Unit{Interval{0.0, 1.0}, Interval{0.0, 1.0}};

/// The crossing quadric 2 x'^2 - 3 y'^2 over the unit square.
BSplineSurface Crossing()
{
    return Quadric(2.0, -3.0, Unit);
}

/// The paraboloid x'^2 + y'^2, wider than the crossing quadric.
BSplineSurface Bowl()
{
    return Quadric(1.0, 1.0, Box2{Interval{-0.5, 1.5}, Interval{-0.5, 1.5}});
}

/// z = (x - 1/2)(y - 1/2), bilinear: its branches run along x = 1/2 and
/// y = 1/2.
BSplineSurface Saddle()
{
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.25}, {0.0, 1.0, -0.25}, {1.0, 0.0, -0.25}, {1.0, 1.0, 0.25}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = 0 over [-1, 2]^2, bilinear.
BSplineSurface Flat()
{
    const std::vector<Point3> poles{
        {-1.0, -1.0, 0.0}, {-1.0, 2.0, 0.0}, {2.0, -1.0, 0.0}, {2.0, 2.0, 0.0}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// The rational quadratic arc of the unit circle from angle a0 to a1.
struct CircleArc {
    std::array<std::array<double, 2>, 3> Poles;
    std::array<double, 3> Weights;
};

/// The arc's ends, its tangents' meeting point, and the cosine of half the
/// span as the middle weight.
CircleArc ArcOf(const double a0, const double a1)
{
    const double half = 0.5 * (a1 - a0);
    const double middle = 0.5 * (a0 + a1);
    const double c = std::cos(half);
    return {{{{std::cos(a0), std::sin(a0)},
              {std::cos(middle) / c, std::sin(middle) / c},
              {std::cos(a1), std::sin(a1)}}},
            {1.0, c, 1.0}};
}

/// The unit cylinder about the z axis over an arc, z in [-1, 1]; or, with
/// aboutX, about the x axis, the arc in (y, z), x in [-1, 1].
BSplineSurface Cylinder(const CircleArc& arc, const bool aboutX)
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < 3; ++i) {
        for (const double along : {-1.0, 1.0}) {
            const auto [c, s] = arc.Poles.at(i);
            poles.push_back(aboutX ? Point3{along, c, s} : Point3{c, s, along});
            weights.push_back(arc.Weights.at(i));
        }
    }
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 1, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

// --- Transformations of the input -----------------------------------------

/// The surface with every pole moved by f: a NURBS surface follows exactly.
BSplineSurface Mapped(BSplineSurface s, const std::function<Point3(const Point3&)>& f)
{
    for (Point3& p : s.Poles) {
        p = f(p);
    }
    return s;
}

/// The number of poles along v.
std::size_t CountV(const BSplineSurface& s)
{
    return s.KnotsV.size() - static_cast<std::size_t>(s.DegreeV) - 1;
}

/// The same surface, u running backwards over [0, 1]: poles reversed along u.
BSplineSurface ReversedU(const BSplineSurface& s)
{
    const std::size_t nv = CountV(s);
    const std::size_t nu = s.Poles.size() / nv;
    std::vector<Point3> poles(s.Poles.size());
    std::vector<double> weights(s.Weights.size());
    for (std::size_t i = 0; i < nu; ++i) {
        for (std::size_t j = 0; j < nv; ++j) {
            poles[(nu - 1 - i) * nv + j] = s.Poles[i * nv + j];
            if (!weights.empty()) {
                weights[(nu - 1 - i) * nv + j] = s.Weights[i * nv + j];
            }
        }
    }
    std::vector<double> knots;
    for (auto k = s.KnotsU.rbegin(); k != s.KnotsU.rend(); ++k) {
        knots.push_back(1.0 - *k);
    }
    const auto made =
        Stereon::MakeBSplineSurface(s.DegreeU, s.DegreeV, knots, s.KnotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// The same surface, u and v exchanged.
BSplineSurface Transposed(const BSplineSurface& s)
{
    const std::size_t nv = CountV(s);
    const std::size_t nu = s.Poles.size() / nv;
    std::vector<Point3> poles(s.Poles.size());
    std::vector<double> weights(s.Weights.size());
    for (std::size_t i = 0; i < nu; ++i) {
        for (std::size_t j = 0; j < nv; ++j) {
            poles[j * nu + i] = s.Poles[i * nv + j];
            if (!weights.empty()) {
                weights[j * nu + i] = s.Weights[i * nv + j];
            }
        }
    }
    const auto made =
        Stereon::MakeBSplineSurface(s.DegreeV, s.DegreeU, s.KnotsV, s.KnotsU, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// A translation and a scale applied to space, and its inverse.
struct Motion {
    double Scale = 1.0;
    Point3 Shift{0.0, 0.0, 0.0};
    /// The point moved.
    [[nodiscard]] Point3 Forward(const Point3& p) const
    {
        return {Scale * p.X + Shift.X, Scale * p.Y + Shift.Y, Scale * p.Z + Shift.Z};
    }
    /// The moved point taken back.
    [[nodiscard]] Point3 Back(const Point3& p) const
    {
        return {(p.X - Shift.X) / Scale, (p.Y - Shift.Y) / Scale, (p.Z - Shift.Z) / Scale};
    }
};

/// The context with its resolution and ceiling scaled by a power of two,
/// the ceiling moved first when both grow.
Context Scaled(const double scale)
{
    const double resolution = ctx.GetResolution() * scale;
    const double ceiling = ctx.GetToleranceCeiling() * scale;
    return scale >= 1.0 ? ctx.WithToleranceCeiling(ceiling).WithResolution(resolution)
                        : ctx.WithResolution(resolution).WithToleranceCeiling(ceiling);
}

// --- The closed-form oracle -------------------------------------------------

/// What a case must answer: the node, the ends of the branches in the
/// request, and how far a point is from the closed-form curve.
struct Expected {
    Point3 Node;
    std::vector<Point3> Ends;
    std::function<double(const Point3&)> OffCurve;
};

/// 2 x'^2 - 3 y'^2 = x'^2 + y'^2: x' = +-2 y', z = 5 y'^2, its ends on the
/// unit square's sides x = 0 and x = 1 at y = 1/2 -+ 1/4. `keep` picks the
/// ends in the request.
Expected QuadricCrossing(const std::function<bool(const Point3&)>& keep)
{
    constexpr double z = 5.0 / 16.0;
    Expected e{{0.5, 0.5, 0.0}, {}, [](const Point3& p) {
                   const double x = p.X - 0.5;
                   const double y = p.Y - 0.5;
                   return std::max(std::abs(p.Z - 5.0 * y * y),
                                   std::abs(std::abs(x) - 2.0 * std::abs(y)));
               }};
    for (const Point3 end :
         {Point3{1.0, 0.75, z}, Point3{0.0, 0.25, z}, Point3{1.0, 0.25, z}, Point3{0.0, 0.75, z}}) {
        if (keep(end)) {
            e.Ends.push_back(end);
        }
    }
    return e;
}

/// x^2 + y^2 = 1 = y^2 + z^2 through (0, 1, 0): x = +-z, each half leaving
/// the upright cylinder's arc where x = +-sqrt(1/2).
Expected CylinderCrossing(const std::function<bool(const Point3&)>& keep)
{
    const double s = std::sqrt(0.5);
    Expected e{{0.0, 1.0, 0.0}, {}, [](const Point3& p) {
                   return std::max(std::abs(p.X * p.X + p.Y * p.Y - 1.0),
                                   std::abs(p.Y * p.Y + p.Z * p.Z - 1.0));
               }};
    for (const Point3 end :
         {Point3{s, s, s}, Point3{s, s, -s}, Point3{-s, s, s}, Point3{-s, s, -s}}) {
        if (keep(end)) {
            e.Ends.push_back(end);
        }
    }
    return e;
}

/// The result is the closed form, moved: one certified branch per end, from
/// the node, every sample on the curve, every pcurve sample in its request.
void CheckCut(const Stereon::SurfaceSurfaceIntersection& r, const Surface& a, const Box2& boxA,
              const Surface& b, const Box2& boxB, const Expected& e, const Motion& motion,
              const double resolution)
{
    REQUIRE(r.Overlaps.empty());
    REQUIRE(r.Points.empty());
    REQUIRE(r.Branches.size() == e.Ends.size());
    REQUIRE(IntersectTest::ConnectivityOf(r, 1e-9 * motion.Scale).Curves == 1);
    const double reach = 4.0 * resolution / motion.Scale;
    const double largest =
        std::max({std::abs(motion.Shift.X), std::abs(motion.Shift.Y), std::abs(motion.Shift.Z)}) +
        2.0 * motion.Scale;
    const double rounding = 4.0 * std::numeric_limits<double>::epsilon() * largest;
    std::vector<bool> matched(e.Ends.size(), false);
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= resolution);
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 64) <=
                branch.Bound + 1e-12 * motion.Scale + rounding);
        Point3 first = motion.Back(Eval(branch.Curve3d, branch.Range.GetLower()));
        Point3 last = motion.Back(Eval(branch.Curve3d, branch.Range.GetUpper()));
        if (Distance(last, e.Node) < Distance(first, e.Node)) {
            std::swap(first, last);
        }
        INFO("branch from (" << first.X << ", " << first.Y << ", " << first.Z << ") to (" << last.X
                             << ", " << last.Y << ", " << last.Z << ")");
        REQUIRE(Distance(first, e.Node) <= reach);
        bool found = false;
        for (std::size_t k = 0; k < e.Ends.size(); ++k) {
            if (!matched[k] && Distance(last, e.Ends[k]) <= reach) {
                matched[k] = found = true;
                break;
            }
        }
        REQUIRE(found);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            REQUIRE(e.OffCurve(motion.Back(Eval(branch.Curve3d, t))) <= reach);
            REQUIRE(Contains(boxA, Eval(branch.OnA, t)));
            REQUIRE(Contains(boxB, Eval(branch.OnB, t)));
        }
    }
}

/// Intersect the moved surfaces over their requests with the tolerance and
/// resolution scaled alike, in both orders, and check each.
void CheckMoved(const BSplineSurface& sa, const Box2& boxA, const BSplineSurface& sb,
                const Box2& boxB, const Expected& e, const Motion& motion)
{
    const auto move = [&](const Point3& p) { return motion.Forward(p); };
    const Surface a{Mapped(sa, move)};
    const Surface b{Mapped(sb, move)};
    const Context scaled = Scaled(motion.Scale);
    for (const bool swapped : {false, true}) {
        INFO("scale " << motion.Scale << " shift (" << motion.Shift.X << ", " << motion.Shift.Y
                      << ", " << motion.Shift.Z << ")" << (swapped ? " swapped" : ""));
        const auto result = swapped ? Intersect(b, boxB, a, boxA, tolerance * motion.Scale, scaled)
                                    : Intersect(a, boxA, b, boxB, tolerance * motion.Scale, scaled);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        if (swapped) {
            CheckCut(*result, b, boxB, a, boxA, e, motion, scaled.GetResolution());
        } else {
            CheckCut(*result, a, boxA, b, boxB, e, motion, scaled.GetResolution());
        }
    }
}

const Box2 Everything = Box2::Everything();
const Box2 LeftHalf{Interval{0.0, 0.5}, Interval{0.0, 1.0}};
const Box2 RightHalf{Interval{0.5, 1.0}, Interval{0.0, 1.0}};
const Box2 LowerHalf{Interval{0.0, 1.0}, Interval{0.0, 0.5}};
const Box2 LowerLeft{Interval{0.0, 0.5}, Interval{0.0, 0.5}};

} // namespace

/// A requested cut through the node of two biquadratics, along u, along v and
/// at a corner of both.
TEST_CASE("a node on a cut of the requested domain answers the half-branches inside it",
          "[intersect][generic][node][cut]")
{
    // Before nodes on cuts were settled, these failed closed.
    CheckMoved(Crossing(), LeftHalf, Bowl(), Everything,
               QuadricCrossing([](const Point3& p) { return p.X < 0.5; }), {});
    CheckMoved(Crossing(), RightHalf, Bowl(), Everything,
               QuadricCrossing([](const Point3& p) { return p.X > 0.5; }), {});
    CheckMoved(Crossing(), LowerHalf, Bowl(), Everything,
               QuadricCrossing([](const Point3& p) { return p.Y < 0.5; }), {});
    CheckMoved(Crossing(), LowerLeft, Bowl(), Everything,
               QuadricCrossing([](const Point3& p) { return p.X < 0.5 && p.Y < 0.5; }), {});
}

/// The surface's own domain ends at the node: an end that does not close.
TEST_CASE("a node on an end of a surface's domain answers the half-branches on it",
          "[intersect][generic][node][cut]")
{
    const BSplineSurface half = Quadric(2.0, -3.0, Box2{Interval{0.0, 0.5}, Interval{0.0, 1.0}});
    CheckMoved(half, Everything, Bowl(), Everything,
               QuadricCrossing([](const Point3& p) { return p.X < 0.5; }), {});
}

/// A rational node on a cut: two exact cylinders, the upright one cut at the
/// node's angle.
TEST_CASE("a rational node on a cut answers the half-branches inside it",
          "[intersect][generic][node][cut][rational]")
{
    const BSplineSurface upright = Cylinder(ArcOf(Pi / 4.0, 3.0 * Pi / 4.0), false);
    const BSplineSurface lying = Cylinder(ArcOf(-Pi / 3.0, Pi / 3.0), true);
    // u = 1/2 is the angle 90 degrees, (0, 1, z); u < 1/2 holds x > 0.
    CheckMoved(upright, LeftHalf, lying, Everything,
               CylinderCrossing([](const Point3& p) { return p.X > 0.0; }), {});
}

/// MR-01 to MR-04 on the cut fixtures, the request moved with the
/// parameters.
TEST_CASE("cut node answers move with the surfaces", "[intersect][generic][node][cut][metamorphic]")
{
    const auto left = [](const Point3& p) { return p.X < 0.5; };
    SECTION("MR-01: translations, exact in binary")
    {
        for (const Point3 shift : {Point3{16.0, -32.0, 8.0}, Point3{1024.0, -2048.0, 512.0},
                                   Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(Crossing(), LeftHalf, Bowl(), Everything, QuadricCrossing(left),
                       {1.0, shift});
        }
    }
    SECTION("MR-03: scaling by powers of two, with the tolerance and resolution")
    {
        for (const int k : {-6, 3, 10}) {
            CheckMoved(Crossing(), LeftHalf, Bowl(), Everything, QuadricCrossing(left),
                       {std::ldexp(1.0, k), Point3{}});
        }
    }
    SECTION("MR-04: reversed and transposed parameters, the request with them")
    {
        CheckMoved(ReversedU(Crossing()), RightHalf, Bowl(), Everything, QuadricCrossing(left), {});
        CheckMoved(Transposed(Crossing()), LowerHalf, Transposed(Bowl()), Everything,
                   QuadricCrossing(left), {});
    }
    // MR-02, the exchange of the surfaces, is in every CheckMoved.
}

/// Sub-cases the proof does not settle, failing closed.
TEST_CASE("nodes on cuts the proof does not settle fail closed",
          "[intersect][generic][node][cut][special]")
{
    const auto failsClosed = [](const Surface& a, const Box2& boxA, const Surface& b,
                                const Box2& boxB) {
        for (const bool swapped : {false, true}) {
            const auto result = swapped ? Intersect(b, boxB, a, boxA, tolerance, ctx)
                                        : Intersect(a, boxA, b, boxB, tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    };
    SECTION("a branch along the cut: the saddle's line x = 1/2 on the cut u = 1/2")
    {
        failsClosed(Surface{Saddle()}, LeftHalf, Surface{Flat()}, Everything);
    }
    SECTION("a node on cuts of both surfaces at once")
    {
        const Box2 bowlLeft{Interval{-0.5, 0.5}, Interval{-0.5, 1.5}};
        failsClosed(Surface{Crossing()}, LeftHalf, Surface{Bowl()}, bowlLeft);
    }
}
