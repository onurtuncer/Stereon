// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Isolated touches in the generic SSI (ADR-0010): polynomial and rational
// B-spline surfaces meeting at one point where they share a tangent plane and
// part on every side. Each answer is held to a closed form written here: one
// point, no branch, its parameters on both surfaces, and a bound within the
// resolution for an exact touch. The same answer, moved accordingly, under
// exchange of the surfaces, reversal and transposition of parameters, a far
// translation and scaling by powers of two (ADR-0053, MR-01 to MR-04). Near
// events: a gap within the tolerance with the surfaces apart is a declared
// touch (ADR-0050); a small loop of the zero set wider than the resolution is
// never reported as a point; four tolerances apart is proved empty. A touch
// on a knot line is not settled and fails closed.

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
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();
/// A pair tolerance ten resolutions wide, as for the tangencies.
constexpr double tolerance = 1e-4;
constexpr double Pi = std::numbers::pi;

/// z = alpha x'^2 + beta y'^2 over the unit square, x' = x - 1/2, biquadratic:
/// the Bernstein coefficients of x'^2 over [0, 1] are 1/4, -1/4 and 1/4.
BSplineSurface Bowl(const double alpha, const double beta, const double lift = 0.0)
{
    constexpr std::array<double, 3> square{0.25, -0.25, 0.25};
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({static_cast<double>(i) / 2.0, static_cast<double>(j) / 2.0,
                             lift + alpha * square.at(i) + beta * square.at(j)});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// A bilinear plane through `origin`, spanned by x and y over the unit
/// square, lx and ly long.
BSplineSurface PlanePatch(const Point3& origin, const Point3& x, const Point3& y, const double lx,
                          const double ly)
{
    std::vector<Point3> poles;
    for (const double i : {0.0, 1.0}) {
        for (const double j : {0.0, 1.0}) {
            poles.push_back({origin.X + i * lx * x.X + j * ly * y.X,
                             origin.Y + i * lx * x.Y + j * ly * y.Y,
                             origin.Z + i * lx * x.Z + j * ly * y.Z});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = height over [-1, 2]^2: wider than the bowls.
BSplineSurface Wide(const double height)
{
    return PlanePatch({-1.0, -1.0, height}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 3.0, 3.0);
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

/// A patch of the unit sphere about `centre`: longitudes of the arc `turn`
/// about the z axis, latitudes from -40 to 40 degrees, exact, the products of
/// the two arcs' poles and weights.
BSplineSurface SpherePatch(const Point3& centre, const CircleArc& turn)
{
    const CircleArc latitude = ArcOf(-2.0 * Pi / 9.0, 2.0 * Pi / 9.0);
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const auto [cx, cy] = turn.Poles.at(i);
            const auto [rho, z] = latitude.Poles.at(j);
            poles.push_back({centre.X + rho * cx, centre.Y + rho * cy, centre.Z + z});
            weights.push_back(turn.Weights.at(i) * latitude.Weights.at(j));
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// The unit sphere about the origin near (1, 0, 0): longitudes -40 to 40.
BSplineSurface EastCap()
{
    return SpherePatch({0.0, 0.0, 0.0}, ArcOf(-2.0 * Pi / 9.0, 2.0 * Pi / 9.0));
}

/// The unit sphere about (2, 0, 0) near (1, 0, 0): longitudes 140 to 220.
BSplineSurface WestCap()
{
    return SpherePatch({2.0, 0.0, 0.0}, ArcOf(7.0 * Pi / 9.0, 11.0 * Pi / 9.0));
}

/// The plane x = 1 over y, z in [-1.5, 1.5].
BSplineSurface Wall()
{
    return PlanePatch({1.0, -1.5, -1.5}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, 3.0, 3.0);
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

/// The same surface, u running backwards: poles reversed along u.
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

/// The result is the closed form, moved: one point at the touch, no branch,
/// its bound within `bound`, its parameters on both surfaces at it.
void CheckTouch(const Stereon::SurfaceSurfaceIntersection& r, const Surface& a, const Surface& b,
                const Point3& touch, const Motion& motion, const double resolution,
                const double bound)
{
    REQUIRE(r.Overlaps.empty());
    REQUIRE(r.Branches.empty());
    REQUIRE(r.Points.size() == 1);
    const Stereon::IntersectionPoint& p = r.Points.front();
    REQUIRE(p.Bound <= bound);
    const double largest =
        std::max({std::abs(motion.Shift.X), std::abs(motion.Shift.Y), std::abs(motion.Shift.Z)}) +
        4.0 * motion.Scale;
    const double rounding = 8.0 * std::numeric_limits<double>::epsilon() * largest;
    INFO("point (" << p.Position.X << ", " << p.Position.Y << ", " << p.Position.Z << ")");
    REQUIRE(Distance(motion.Back(p.Position), touch) <= 4.0 * resolution / motion.Scale);
    REQUIRE(Distance(Eval(a, p.OnA), p.Position) <= p.Bound + rounding);
    REQUIRE(Distance(Eval(b, p.OnB), p.Position) <= p.Bound + rounding);
}

/// Intersect the moved surfaces with the tolerance and resolution scaled
/// alike, in both orders, and check each against the moved touch.
void CheckMoved(const BSplineSurface& sa, const BSplineSurface& sb, const Point3& touch,
                const Motion& motion)
{
    const auto move = [&](const Point3& p) { return motion.Forward(p); };
    const Surface a{Mapped(sa, move)};
    const Surface b{Mapped(sb, move)};
    const Context scaled = Scaled(motion.Scale);
    for (const bool swapped : {false, true}) {
        INFO("scale " << motion.Scale << " shift (" << motion.Shift.X << ", " << motion.Shift.Y
                      << ", " << motion.Shift.Z << ")" << (swapped ? " swapped" : ""));
        const auto result = swapped ? Intersect(b, Box2::Everything(), a, Box2::Everything(),
                                                tolerance * motion.Scale, scaled)
                                    : Intersect(a, Box2::Everything(), b, Box2::Everything(),
                                                tolerance * motion.Scale, scaled);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE(result.has_value());
        CheckTouch(*result, swapped ? b : a, swapped ? a : b, touch, motion, scaled.GetResolution(),
                   scaled.GetResolution());
    }
}

const Point3 BowlBottom{0.5, 0.5, 0.0};
const Point3 East{1.0, 0.0, 0.0};

} // namespace

/// A polynomial touch at a bowl's bottom, and two bowls back to back.
TEST_CASE("a polynomial bowl touching a plane, or another bowl, is one point",
          "[intersect][generic][touch]")
{
    // Before isolated touches were settled, these failed closed.
    CheckMoved(Bowl(1.0, 2.0), Wide(0.0), BowlBottom, {});
    CheckMoved(Bowl(1.0, 1.0), Bowl(-2.0, -1.0), BowlBottom, {});
}

/// Rational touches: an exact sphere on a plane, and two exact spheres.
TEST_CASE("an exact sphere touching a plane, or another sphere, is one point",
          "[intersect][generic][touch][rational]")
{
    CheckMoved(EastCap(), Wall(), East, {});
    CheckMoved(EastCap(), WestCap(), East, {});
}

/// MR-01 to MR-04 on the touch fixtures.
TEST_CASE("touch answers move with the surfaces", "[intersect][generic][touch][metamorphic]")
{
    SECTION("MR-01: translations, exact in binary")
    {
        for (const Point3 shift : {Point3{16.0, -32.0, 8.0}, Point3{1024.0, -2048.0, 512.0},
                                   Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(Bowl(1.0, 2.0), Wide(0.0), BowlBottom, {1.0, shift});
            CheckMoved(EastCap(), WestCap(), East, {1.0, shift});
        }
    }
    SECTION("MR-03: scaling by powers of two, with the tolerance and resolution")
    {
        for (const int k : {-6, 3, 10}) {
            CheckMoved(Bowl(1.0, 1.0), Bowl(-2.0, -1.0), BowlBottom,
                       {std::ldexp(1.0, k), Point3{}});
            CheckMoved(EastCap(), Wall(), East, {std::ldexp(1.0, k), Point3{}});
        }
    }
    SECTION("MR-04: reversed and transposed parameters leave the point in place")
    {
        CheckMoved(ReversedU(Bowl(1.0, 2.0)), Transposed(Wide(0.0)), BowlBottom, {});
        CheckMoved(Transposed(ReversedU(EastCap())), ReversedU(WestCap()), East, {});
    }
    // MR-02, the exchange of the surfaces, is in every CheckMoved.
}

/// Near events: gaps within and beyond the tolerance, on either side.
TEST_CASE("a gap within the tolerance apart is a declared touch; a wider loop is never a point",
          "[intersect][generic][touch][near-event]")
{
    const Surface bowl{Bowl(1.0, 2.0)};
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "swapped" : "in order"));
        const auto run = [&](const double height) {
            const Surface flat{Wide(height)};
            return swapped ? Intersect(flat, Box2::Everything(), bowl, Box2::Everything(),
                                       tolerance, ctx)
                           : Intersect(bowl, Box2::Everything(), flat, Box2::Everything(),
                                       tolerance, ctx);
        };
        SECTION("half a tolerance below the bottom: apart, a declared touch")
        {
            const auto result = run(-0.5 * tolerance);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->Branches.empty());
            REQUIRE(result->Points.size() == 1);
            const Stereon::IntersectionPoint& p = result->Points.front();
            REQUIRE(p.Bound <= tolerance);
            REQUIRE(p.Bound >= 0.5 * tolerance - 4.0 * ctx.GetResolution());
            REQUIRE(std::abs(p.Position.X - 0.5) <= 4.0 * ctx.GetResolution());
            REQUIRE(std::abs(p.Position.Y - 0.5) <= 4.0 * ctx.GetResolution());
        }
        SECTION("four tolerances below: proved apart")
        {
            const auto result = run(-4.0 * tolerance);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->IsEmpty());
        }
        SECTION("half a tolerance above: a loop of radius about 0.007, never a point")
        {
            const double height = 0.5 * tolerance;
            const auto result = run(height);
            INFO(IntersectTest::ErrorText(result));
            if (result) {
                REQUIRE(result->Points.empty());
                REQUIRE_FALSE(result->Branches.empty());
                for (const Stereon::IntersectionBranch& branch : result->Branches) {
                    for (const double t : IntersectTest::SamplesOf(branch.Range, 32)) {
                        const Point3 q = Eval(branch.Curve3d, t);
                        const double x = q.X - 0.5;
                        const double y = q.Y - 0.5;
                        REQUIRE(std::abs(x * x + 2.0 * y * y - height) <=
                                4.0 * ctx.GetResolution());
                    }
                }
            } else {
                REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
            }
        }
    }
}

/// A sub-case the proof does not settle, failing closed.
TEST_CASE("a touch on a knot line of the bowl fails closed", "[intersect][generic][touch][special]")
{
    // z = x'^2 + y'^2 in two Bezier pieces, x in [0, 1/2] and [1/2, 1]: the
    // bottom is on the knot line, in no single patch.
    std::vector<Point3> poles;
    const std::array<double, 5> xs{0.0, 0.25, 0.5, 0.75, 1.0};
    const std::array<double, 5> squares{0.25, 0.0, 0.0, 0.0, 0.25};
    const std::array<double, 3> ys{0.0, 0.5, 1.0};
    const std::array<double, 3> ysq{0.25, -0.25, 0.25};
    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({xs.at(i), ys.at(j), squares.at(i) + ysq.at(j)});
        }
    }
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knotsU, knotsV, poles, {});
    REQUIRE(made.has_value());
    const Surface pieced{*made};
    const Surface flat{Wide(0.0)};
    for (const bool swapped : {false, true}) {
        const auto result =
            swapped
                ? Intersect(flat, Box2::Everything(), pieced, Box2::Everything(), tolerance, ctx)
                : Intersect(pieced, Box2::Everything(), flat, Box2::Everything(), tolerance, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}
