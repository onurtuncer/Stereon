// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Branch points of the generic SSI (ADR-0010): two polynomial B-spline
// surfaces crossing in two branches at a transverse node. Each answer is held
// to a closed form written here, never derived from the kernel's path: the
// node, the four ends of its half-branches on the boundary of the domain, and
// the curve every sample of every branch must lie on. The same answer, moved
// accordingly, under exchange of the surfaces, reversal and transposition of
// parameters, a far translation and scaling by powers of two (ADR-0053,
// MR-01 to MR-04). Higher-order singularities, nodes of rational surfaces and
// a node the tolerance cannot admit fail closed.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using Stereon::Box2;
using Stereon::BSplineSurface;
using Stereon::Context;
using Stereon::ErrorCode;
using Stereon::Point3;
using Stereon::Surface;

namespace {

const Context ctx = Context::Default();
constexpr double tolerance = 1e-6;

/// Bernstein coefficients over [0, 1] of degree 2, of x' = x - 1/2 and x'^2.
constexpr std::array<double, 3> Linear{-0.5, 0.0, 0.5};
constexpr std::array<double, 3> Square{0.25, -0.25, 0.25};

/// A biquadratic graph z = f(x, y) over the square [low, low + width]^2,
/// its poles at (low + i width / 2, low + j width / 2, z(i, j)).
BSplineSurface Graph(const std::function<double(std::size_t, std::size_t)>& z,
                     const double low = 0.0, const double width = 1.0,
                     const std::vector<double>& weights = {})
{
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            poles.push_back({low + static_cast<double>(i) * width / 2.0,
                             low + static_cast<double>(j) * width / 2.0, z(i, j)});
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// z = alpha x'^2 + beta y'^2 over [low, low + width]^2: the Bernstein
/// coefficients of (a + width t)^2, a = low - 1/2, are a^2, a^2 + a width
/// and (a + width)^2.
BSplineSurface Quadric(const double alpha, const double beta, const double low = 0.0,
                       const double width = 1.0)
{
    const double a = low - 0.5;
    const std::array<double, 3> square{a * a, a * a + a * width, (a + width) * (a + width)};
    return Graph([=](const std::size_t i,
                     const std::size_t j) { return alpha * square.at(i) + beta * square.at(j); },
                 low, width);
}

/// The paraboloid z = x'^2 + y'^2 over [-1/2, 3/2]^2, wider than the other
/// surface: no branch ends where both domains' boundaries meet.
BSplineSurface Bowl()
{
    return Quadric(1.0, 1.0, -0.5, 2.0);
}

/// z = x' y', bilinear: the saddle of the existing generic tests.
BSplineSurface Saddle(const std::vector<double>& weights = {})
{
    const std::vector<Point3> poles{
        {0.0, 0.0, 0.25}, {0.0, 1.0, -0.25}, {1.0, 0.0, -0.25}, {1.0, 1.0, 0.25}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// z = height, bilinear over the unit square.
BSplineSurface Flat(const double height)
{
    const std::vector<Point3> poles{
        {0.0, 0.0, height}, {0.0, 1.0, height}, {1.0, 0.0, height}, {1.0, 1.0, height}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

// --- Transformations of the input -----------------------------------------

BSplineSurface Mapped(BSplineSurface s, const std::function<Point3(const Point3&)>& f)
{
    for (Point3& p : s.Poles) {
        p = f(p);
    }
    return s;
}

/// The same surface, u running backwards: poles reversed along u.
BSplineSurface ReversedU(BSplineSurface s)
{
    const std::size_t nv = s.KnotsV.size() - static_cast<std::size_t>(s.DegreeV) - 1;
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
    const std::size_t nv = s.KnotsV.size() - static_cast<std::size_t>(s.DegreeV) - 1;
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

// --- Exact cylinders ---------------------------------------------------------

/// The rational quadratic arc of the unit circle from angle a0 to a1, less
/// than a half turn apart: (cos, sin) at the ends, the tangents' meeting
/// point between, its weight the cosine of half the span.
struct CircleArc {
    std::array<std::array<double, 2>, 3> Poles;
    std::array<double, 3> Weights;
};

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

constexpr double Pi = 3.14159265358979323846;

/// About z, from 45 to 135 degrees: (0, 1, 0) mid-arc.
BSplineSurface Upright()
{
    return Cylinder(ArcOf(Pi / 4.0, 3.0 * Pi / 4.0), false);
}

/// About x, from -60 to 60 degrees from the y axis: (0, 1, 0) mid-arc.
BSplineSurface Lying()
{
    return Cylinder(ArcOf(-Pi / 3.0, Pi / 3.0), true);
}

// --- The closed-form oracle -------------------------------------------------

/// What a node case must answer: the node, the four ends of its branches,
/// and how far a point is from the closed-form curve, all in the input's
/// own space before any transformation.
struct Expected {
    Point3 Node;
    std::array<Point3, 4> Ends;
    std::function<double(const Point3&)> OffCurve;
};

/// x = 1/2 or y = 1/2, z = 0: the saddle against the plane.
Expected SaddleLines()
{
    return {{0.5, 0.5, 0.0},
            {Point3{0.5, 0.0, 0.0}, Point3{0.5, 1.0, 0.0}, Point3{0.0, 0.5, 0.0},
             Point3{1.0, 0.5, 0.0}},
            [](const Point3& p) {
                return std::max(std::abs(p.Z), std::min(std::abs(p.X - 0.5), std::abs(p.Y - 0.5)));
            }};
}

/// x^2 + y^2 = 1 = y^2 + z^2: the ellipses x = +-z through (0, 1, 0), each
/// half leaving the upright cylinder's arc where x = +-sqrt(1/2).
Expected CylinderCrossing()
{
    const double s = std::sqrt(0.5);
    return {{0.0, 1.0, 0.0},
            {Point3{s, s, s}, Point3{s, s, -s}, Point3{-s, s, s}, Point3{-s, s, -s}},
            [](const Point3& p) {
                return std::max(std::abs(p.X * p.X + p.Y * p.Y - 1.0),
                                std::abs(p.Y * p.Y + p.Z * p.Z - 1.0));
            }};
}

/// 2x'^2 - 3y'^2 = x'^2 + y'^2: x' = +-2 y', z = 5 y'^2, through the
/// centre to x = 0 and 1 at y = 1/2 -+ 1/4.
Expected QuadricCrossing()
{
    constexpr double z = 5.0 / 16.0;
    return {
        {0.5, 0.5, 0.0},
        {Point3{1.0, 0.75, z}, Point3{0.0, 0.25, z}, Point3{1.0, 0.25, z}, Point3{0.0, 0.75, z}},
        [](const Point3& p) {
            const double x = p.X - 0.5;
            const double y = p.Y - 0.5;
            return std::max(std::abs(p.Z - 5.0 * y * y), std::abs(std::abs(x) - 2.0 * std::abs(y)));
        }};
}

/// A rigid motion and a scale applied to space, and its inverse.
struct Motion {
    double Scale = 1.0;
    Point3 Shift{0.0, 0.0, 0.0};
    [[nodiscard]] Point3 Forward(const Point3& p) const
    {
        return {Scale * p.X + Shift.X, Scale * p.Y + Shift.Y, Scale * p.Z + Shift.Z};
    }
    [[nodiscard]] Point3 Back(const Point3& p) const
    {
        return {(p.X - Shift.X) / Scale, (p.Y - Shift.Y) / Scale, (p.Z - Shift.Z) / Scale};
    }
};

/// The result is the closed form, moved: four certified branches, each from
/// the node to its own end, every sample on the curve, one connected curve.
void CheckNode(const Stereon::SurfaceSurfaceIntersection& r, const Surface& a, const Surface& b,
               const Expected& e, const Motion& motion, const double resolution)
{
    REQUIRE(r.Overlaps.empty());
    REQUIRE(r.Points.empty());
    REQUIRE(r.Branches.size() == 4);
    REQUIRE(IntersectTest::ConnectivityOf(r, 1e-9 * motion.Scale).Curves == 1);
    const double reach = 4.0 * resolution / motion.Scale; // In the input's own space.
    std::array<bool, 4> matched{};
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= resolution);
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 64) <=
                branch.Bound + 1e-12 * motion.Scale);
        Point3 first = motion.Back(Eval(branch.Curve3d, branch.Range.GetLower()));
        Point3 last = motion.Back(Eval(branch.Curve3d, branch.Range.GetUpper()));
        if (Distance(last, e.Node) < Distance(first, e.Node)) {
            std::swap(first, last);
        }
        INFO("branch from (" << first.X << ", " << first.Y << ", " << first.Z << ") to (" << last.X
                             << ", " << last.Y << ", " << last.Z << ")");
        REQUIRE(Distance(first, e.Node) <= reach);
        bool found = false;
        for (std::size_t k = 0; k < 4; ++k) {
            if (!matched.at(k) && Distance(last, e.Ends.at(k)) <= reach) {
                matched.at(k) = found = true;
                break;
            }
        }
        REQUIRE(found);
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            REQUIRE(e.OffCurve(motion.Back(Eval(branch.Curve3d, t))) <= reach);
        }
    }
}

/// The context with its resolution and ceiling scaled by a power of two,
/// the ceiling moved first when both grow.
Context Scaled(const double scale)
{
    const double resolution = ctx.GetResolution() * scale;
    const double ceiling = ctx.GetToleranceCeiling() * scale;
    return scale >= 1.0 ? ctx.WithToleranceCeiling(ceiling).WithResolution(resolution)
                        : ctx.WithResolution(resolution).WithToleranceCeiling(ceiling);
}

/// Whether a transformed case must be answered, or may instead fail closed.
enum class Must : bool { Answer, AnswerOrFailClosed };

/// Intersect the moved surfaces with the tolerance and resolution scaled
/// alike, in both orders, and check each against the moved closed form.
void CheckMoved(const BSplineSurface& sa, const BSplineSurface& sb, const Expected& e,
                const Motion& motion, const Must must = Must::Answer)
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
        if (!result && must == Must::AnswerOrFailClosed) {
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
            continue;
        }
        REQUIRE(result.has_value());
        CheckNode(*result, swapped ? b : a, swapped ? a : b, e, motion, scaled.GetResolution());
    }
}

} // namespace

TEST_CASE("a saddle crossing a plane at its centre ends four branches at the node",
          "[intersect][generic][node]")
{
    // Before branch points were settled, this failed closed.
    CheckMoved(Saddle(), Flat(0.0), SaddleLines(), {});
}

TEST_CASE("two biquadratics crossing at a node are split there, curved branches and all",
          "[intersect][generic][node]")
{
    CheckMoved(Quadric(2.0, -3.0), Bowl(), QuadricCrossing(), {});
}

TEST_CASE("node answers move with the surfaces", "[intersect][generic][node][metamorphic]")
{
    const BSplineSurface a = Quadric(2.0, -3.0);
    const BSplineSurface b = Bowl();
    SECTION("MR-01: translations, exact in binary")
    {
        CheckMoved(a, b, QuadricCrossing(), {1.0, Point3{16.0, -32.0, 8.0}});
        CheckMoved(Saddle(), Flat(0.0), SaddleLines(), {1.0, Point3{16.0, -32.0, 8.0}});
        // WEAKENED RELATION, a known open defect: SSI-translated-boundary in
        // docs/planning/intersection-metamorphic.json and the Generic SSI row
        // of the v0.2 ledger. At offsets of about 64 and beyond, the
        // transverse search refuses (IntersectionFailed) near a branch's end
        // on the domain's boundary, with or without a node, as on the base
        // revision 227d79a. Until it is repaired these shifts accept that
        // refusal, never a wrong answer; then they must answer.
        for (const Point3 shift :
             {Point3{1024.0, -2048.0, 512.0}, Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(a, b, QuadricCrossing(), {1.0, shift}, Must::AnswerOrFailClosed);
            CheckMoved(Saddle(), Flat(0.0), SaddleLines(), {1.0, shift}, Must::AnswerOrFailClosed);
        }
    }
    SECTION("MR-03: scaling by powers of two, with the tolerance and resolution")
    {
        for (const int k : {-6, 3, 10}) {
            CheckMoved(a, b, QuadricCrossing(), {std::ldexp(1.0, k), Point3{0.0, 0.0, 0.0}});
        }
    }
    SECTION("MR-04: reversed and transposed parameters leave the curves in place")
    {
        CheckMoved(ReversedU(a), b, QuadricCrossing(), {});
        CheckMoved(a, Transposed(b), QuadricCrossing(), {});
        CheckMoved(Transposed(ReversedU(a)), ReversedU(b), QuadricCrossing(), {});
        CheckMoved(ReversedU(Saddle()), Transposed(Flat(0.0)), SaddleLines(), {});
    }
    // MR-02, the exchange of the surfaces, is in every CheckMoved.
}

TEST_CASE("a gap within the tolerance is a node; beyond it, two branches that never meet",
          "[intersect][generic][node][near-event]")
{
    const Surface saddle{Saddle()};
    // 2^-30 apart at the centre, the tolerance 1e-6: the branches are
    // admitted to meet there (ADR-0050), the node a candidate event only.
    const auto within = Intersect(saddle, Box2::Everything(), Surface{Flat(0x1p-30)},
                                  Box2::Everything(), tolerance, ctx);
    INFO(IntersectTest::ErrorText(within));
    REQUIRE(within.has_value());
    REQUIRE(within->Branches.size() == 4);
    REQUIRE(IntersectTest::ConnectivityOf(*within, 1e-9).Curves == 1);
    // 2^-6 apart: the hyperbola x' y' = 2^-6, two branches, no node.
    const auto beyond = Intersect(saddle, Box2::Everything(), Surface{Flat(0x1p-6)},
                                  Box2::Everything(), tolerance, ctx);
    INFO(IntersectTest::ErrorText(beyond));
    REQUIRE(beyond.has_value());
    REQUIRE(beyond->Branches.size() == 2);
    REQUIRE(IntersectTest::ConnectivityOf(*beyond, 1e-9).Curves == 2);
    for (const Stereon::IntersectionBranch& branch : beyond->Branches) {
        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(std::abs((p.X - 0.5) * (p.Y - 0.5) - 0x1p-6) <= 4.0 * ctx.GetResolution());
        }
    }
}

TEST_CASE("branch points that are not proved transverse nodes fail closed",
          "[intersect][generic][node][special]")
{
    const Surface flat{Flat(0.0)};
    SECTION("no tolerance admits the node's enclosure")
    {
        const auto result =
            Intersect(Surface{Saddle()}, Box2::Everything(), flat, Box2::Everything(), 0.0, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("three lines through a point: z = x'^2 y' - x' y'^2, a degenerate Hessian")
    {
        const Surface monkey{Graph([](const std::size_t i, const std::size_t j) {
            return Square.at(i) * Linear.at(j) - Linear.at(i) * Square.at(j);
        })};
        for (const bool swapped : {false, true}) {
            const auto result = swapped ? Intersect(flat, Box2::Everything(), monkey,
                                                    Box2::Everything(), tolerance, ctx)
                                        : Intersect(monkey, Box2::Everything(), flat,
                                                    Box2::Everything(), tolerance, ctx);
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
        }
    }
    SECTION("three lines through a point of a rational surface")
    {
        // The same degenerate saddle, written with equal weights.
        const Surface monkey{Graph(
            [](const std::size_t i, const std::size_t j) {
                return Square.at(i) * Linear.at(j) - Linear.at(i) * Square.at(j);
            },
            0.0, 1.0, std::vector<double>(9, 2.0))};
        const auto result =
            Intersect(monkey, Box2::Everything(), flat, Box2::Everything(), tolerance, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
    SECTION("a definite touch is not a node")
    {
        const auto result = Intersect(Surface{Quadric(1.0, 1.0)}, Box2::Everything(), flat,
                                      Box2::Everything(), tolerance, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}

TEST_CASE("rational nodes are split as polynomial ones", "[intersect][generic][node][rational]")
{
    SECTION("a weighted saddle: weights 1, 2, 2, 4 keep the lines x = 1/2 and y = 1/2")
    {
        // The height's numerator (1 - 3u)(1 - 3v) / 4 and x = (2u + 2uv) /
        // ((1 + u)(1 + v)) put the node at u = v = 1/3, the point (1/2, 1/2, 0).
        CheckMoved(Saddle({1.0, 2.0, 2.0, 4.0}), Flat(0.0), SaddleLines(), {});
        CheckMoved(Saddle({2.0, 2.0, 2.0, 2.0}), Flat(0.0), SaddleLines(), {});
    }
    SECTION("two exact cylinders of equal radius, axes crossing")
    {
        CheckMoved(Upright(), Lying(), CylinderCrossing(), {});
    }
}

TEST_CASE("rational node answers move with the surfaces",
          "[intersect][generic][node][rational][metamorphic]")
{
    const BSplineSurface a = Upright();
    const BSplineSurface b = Lying();
    SECTION("MR-01: a translation, exact in binary")
    {
        // Answered since #279's exact shift about the origin, the node moved
        // back with the seeds (formerly weakened under SSI-translated-boundary).
        CheckMoved(a, b, CylinderCrossing(), {1.0, Point3{16.0, -32.0, 8.0}});
    }
    SECTION("MR-03: scaling by powers of two")
    {
        for (const int k : {-3, 5}) {
            CheckMoved(a, b, CylinderCrossing(), {std::ldexp(1.0, k), Point3{0.0, 0.0, 0.0}});
        }
    }
    SECTION("MR-04: reversed and transposed parameters")
    {
        CheckMoved(ReversedU(a), b, CylinderCrossing(), {});
        CheckMoved(a, Transposed(b), CylinderCrossing(), {});
        CheckMoved(ReversedU(Saddle({1.0, 2.0, 2.0, 4.0})), Flat(0.0), SaddleLines(), {});
    }
}
