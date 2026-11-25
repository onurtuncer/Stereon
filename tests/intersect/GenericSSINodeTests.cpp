// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Branch points of the generic SSI (ADR-0010): polynomial or rational B-spline
// surfaces crossing in two branches at a transverse node. Each answer is held
// to a closed form written here, never derived from the kernel's path: the
// node, the four ends of its half-branches on the boundary of the domain, and
// the curve every sample of every branch must lie on. The same answer, moved
// accordingly, under exchange of the surfaces, reversal and transposition of
// parameters, a far translation and scaling by powers of two (ADR-0053,
// MR-01 to MR-04). A node on a knot line of either surface, on a crossing of
// two, or on the seam of a closed chart is answered the same way, the Bezier
// patches around it settled together; a crease there within the tolerance is
// not told apart from a smooth surface. Higher-order singularities, a crease
// beyond the tolerance, a seam not closed in the doubles and a node the
// tolerance cannot admit fail closed.

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

constexpr double Pi = std::numbers::pi;

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

// --- Surfaces with knot lines and seams through the node --------------------

/// The Bernstein coefficients of x and of x'^2, x' = x - 1/2, over the
/// spans between the breaks: three per span, the ends shared.
struct Pieces {
    std::vector<double> Points;
    std::vector<double> Squares;
};

Pieces PiecesOver(const std::vector<double>& breaks)
{
    Pieces p{{breaks.front()}, {(breaks.front() - 0.5) * (breaks.front() - 0.5)}};
    for (std::size_t k = 1; k < breaks.size(); ++k) {
        const double a = breaks[k - 1] - 0.5;
        const double b = breaks[k] - 0.5;
        p.Points.insert(p.Points.end(), {0.5 * (breaks[k - 1] + breaks[k]), breaks[k]});
        p.Squares.insert(p.Squares.end(), {a * b, b * b});
    }
    return p;
}

/// Quadratic knots with a double knot at each inner value.
std::vector<double> DoubleKnots(const std::vector<double>& values)
{
    std::vector<double> knots{values.front()};
    for (const double v : values) {
        knots.insert(knots.end(), {v, v});
    }
    knots.push_back(values.back());
    return knots;
}

/// z = alpha x'^2 + beta y'^2 + kink |x - breaksX[1]| over [0, 1]^2, in
/// Bezier pieces between the breaks along x and y, joined at double knots at
/// the given parameters: the same paraboloid as Quadric's when the kink is
/// zero, however its parameters run. The kink is a crease along the first
/// inner break in x, the slope jumping by twice it.
BSplineSurface Pieced(const std::array<double, 2>& form, const std::vector<double>& breaksX,
                      const std::vector<double>& knotsX, const std::vector<double>& breaksY,
                      const std::vector<double>& knotsY, const double kink = 0.0)
{
    const Pieces x = PiecesOver(breaksX);
    const Pieces y = PiecesOver(breaksY);
    std::vector<Point3> poles;
    for (std::size_t i = 0; i < x.Points.size(); ++i) {
        for (std::size_t j = 0; j < y.Points.size(); ++j) {
            const double crease =
                breaksX.size() > 2 ? kink * std::abs(x.Points[i] - breaksX[1]) : 0.0;
            poles.push_back({x.Points[i], y.Points[j],
                             form[0] * x.Squares[i] + form[1] * y.Squares[j] + crease});
        }
    }
    const auto made =
        Stereon::MakeBSplineSurface(2, 2, DoubleKnots(knotsX), DoubleKnots(knotsY), poles, {});
    REQUIRE(made.has_value());
    return *made;
}

const std::vector<double> OneSpan{0.0, 1.0};
/// A break through the node's x or y, at a parameter that is not its
/// middle: the two sides run at different speeds.
const std::vector<double> BreakAtNode{0.0, 0.5, 1.0};
const std::vector<double> KnotOffMiddle{0.0, 0.375, 1.0};

/// The crossing quadric of QuadricCrossing with a knot line along x = 1/2,
/// through the node; with `both`, another along y = 1/2, crossing it there.
BSplineSurface KnottedQuadric(const bool both, const double kink = 0.0)
{
    return Pieced({2.0, -3.0}, BreakAtNode, KnotOffMiddle, both ? BreakAtNode : OneSpan,
                  both ? KnotOffMiddle : OneSpan, kink);
}

/// The unit cylinder over several arcs joined at double knots, as Cylinder
/// builds one.
BSplineSurface CylinderOf(const std::vector<CircleArc>& arcs, const bool aboutX)
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    std::vector<double> values{0.0};
    for (std::size_t a = 0; a < arcs.size(); ++a) {
        for (std::size_t i = a == 0 ? 0 : 1; i < 3; ++i) {
            for (const double along : {-1.0, 1.0}) {
                const auto [c, s] = arcs[a].Poles.at(i);
                poles.push_back(aboutX ? Point3{along, c, s} : Point3{c, s, along});
                weights.push_back(arcs[a].Weights.at(i));
            }
        }
        values.push_back(static_cast<double>(a + 1) / static_cast<double>(arcs.size()));
    }
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made =
        Stereon::MakeBSplineSurface(2, 1, DoubleKnots(values), knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// Upright and Lying, each in two arcs meeting at (0, 1, 0)'s parameter.
BSplineSurface UprightInTwo()
{
    return CylinderOf({ArcOf(Pi / 4.0, Pi / 2.0), ArcOf(Pi / 2.0, 3.0 * Pi / 4.0)}, false);
}

BSplineSurface LyingInTwo()
{
    return CylinderOf({ArcOf(-Pi / 3.0, 0.0), ArcOf(0.0, Pi / 3.0)}, true);
}

/// A full turn of the unit cylinder, the nine-pole rational circle with its
/// middle weights rounded as a file would hold them: about z from (1, 0),
/// or about x from (y, z) = (1, 0), over [-3, 3] along the axis. Its chart
/// closes on the seam at angle zero, the first and last poles the same.
BSplineSurface FullTurn(const bool aboutX)
{
    const std::array<std::array<double, 2>, 9> circle{{{1.0, 0.0},
                                                       {1.0, 1.0},
                                                       {0.0, 1.0},
                                                       {-1.0, 1.0},
                                                       {-1.0, 0.0},
                                                       {-1.0, -1.0},
                                                       {0.0, -1.0},
                                                       {1.0, -1.0},
                                                       {1.0, 0.0}}};
    // sqrt(1/2) to fifteen decimals, as tools/generate_nurbs_corpus.py writes it.
    const double rounded = std::round(0.5 * std::numbers::sqrt2 * 1e15) / 1e15;
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < circle.size(); ++i) {
        for (const double along : {-3.0, 3.0}) {
            const auto [c, s] = circle.at(i);
            poles.push_back(aboutX ? Point3{along, c, s} : Point3{c, s, along});
            weights.push_back(i % 2 == 0 ? 1.0 : rounded);
        }
    }
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 1, DoubleKnots({0.0, 0.25, 0.5, 0.75, 1.0}),
                                                  knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// The length of the ellipse x = z on x^2 + y^2 = 1, (cos, sin, cos): the
/// integral of sqrt(1 + sin^2) over a turn, by Simpson's rule.
double EllipseLength()
{
    constexpr std::size_t n = 4096;
    const double step = 2.0 * Pi / static_cast<double>(n);
    double sum = 0.0;
    for (std::size_t k = 0; k <= n; ++k) {
        const double s = std::sin(static_cast<double>(k) * step);
        const double weight = k == 0 || k == n ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0);
        sum += weight * std::sqrt(1.0 + s * s);
    }
    return sum * step / 3.0;
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
    // The samples are evaluated in doubles at the placed coordinates, each
    // point rounded to their spacing: 9.1e-13 near x = -4096, where the
    // measured gap of a branch certified to 2.7e-13 is 1.37e-12, a gap of
    // one and a half spacings. Four spacings of the largest coordinate,
    // beyond the 1e-12 kept for the evaluation itself.
    const double largest =
        std::max({std::abs(motion.Shift.X), std::abs(motion.Shift.Y), std::abs(motion.Shift.Z)}) +
        2.0 * motion.Scale;
    const double rounding = 4.0 * std::numeric_limits<double>::epsilon() * largest;
    std::array<bool, 4> matched{};
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

/// Intersect the moved surfaces with the tolerance and resolution scaled
/// alike, in both orders, and check each against the moved closed form.
void CheckMoved(const BSplineSurface& sa, const BSplineSurface& sb, const Expected& e,
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
        // #279 repairs SSI-translated-boundary; nodes must move with the arcs.
        for (const Point3 shift :
             {Point3{1024.0, -2048.0, 512.0}, Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(a, b, QuadricCrossing(), {1.0, shift});
            CheckMoved(Saddle(), Flat(0.0), SaddleLines(), {1.0, shift});
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
    // Weights 2^(i+j) reparameterize the bilinear surface by
    // x = 2u/(1+u), y = 2v/(1+v), preserving z = (x-1/2)(y-1/2).
    // The critical parameters move from (1/2,1/2) to (1/3,1/3), while
    // the geometric oracle stays the same and W has nonzero partials.
    for (const bool rational : {false, true}) {
        const Surface saddle{rational ? Saddle({1.0, 2.0, 2.0, 4.0}) : Saddle()};
        for (const bool swapped : {false, true}) {
            INFO("rational " << rational << " swapped " << swapped);
            for (const double height : {0x1p-30, 0x1p-6}) {
                INFO("height " << height);
                const Surface flat{Flat(height)};
                const auto result = swapped ? Intersect(flat, Box2::Everything(), saddle,
                                                        Box2::Everything(), tolerance, ctx)
                                            : Intersect(saddle, Box2::Everything(), flat,
                                                        Box2::Everything(), tolerance, ctx);
                INFO(IntersectTest::ErrorText(result));
                REQUIRE(result.has_value());
                REQUIRE(result->Points.empty());
                REQUIRE(result->Overlaps.empty());
                if (height < tolerance) {
                    // ADR-0050 admits the near crossing as a candidate node.
                    // CheckNode checks its location, all four distinct ends,
                    // connectivity, certificates and both pcurves.
                    CheckNode(*result, swapped ? flat : saddle, swapped ? saddle : flat,
                              SaddleLines(), {}, ctx.GetResolution());
                } else {
                    // x' y' = 2^-6: two disconnected hyperbola branches.
                    REQUIRE(result->Branches.size() == 2);
                    REQUIRE(IntersectTest::ConnectivityOf(*result, 1e-9).Curves == 2);
                    for (const Stereon::IntersectionBranch& branch : result->Branches) {
                        REQUIRE(branch.IsCertified());
                        REQUIRE(branch.Bound <= ctx.GetResolution());
                        REQUIRE(IntersectTest::PcurveDeviation(branch, swapped ? flat : saddle,
                                                               swapped ? saddle : flat,
                                                               64) <= branch.Bound + 1e-12);
                        for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                            const Point3 p = Eval(branch.Curve3d, t);
                            REQUIRE(std::abs((p.X - 0.5) * (p.Y - 0.5) - height) <=
                                    4.0 * ctx.GetResolution());
                            REQUIRE(std::abs(p.Z - height) <= 4.0 * ctx.GetResolution());
                        }
                    }
                }
            }
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
        // Equal weights exercise the rational representation of a polynomial.
        // Weights 2^(i+j) also exercise a genuinely varying denominator:
        // W=(1+u)^2(1+v)^2, x=2u/(1+u), y=2v/(1+v).
        // Thus z=x'^2 y'-x' y'^2 still has three crossing lines and a
        // zero Hessian at (x,y)=(1/2,1/2), now at u=v=1/3.
        for (const std::vector<double>& weights :
             {std::vector<double>(9, 2.0),
              std::vector<double>{1.0, 2.0, 4.0, 2.0, 4.0, 8.0, 4.0, 8.0, 16.0}}) {
            INFO("last weight " << weights.back());
            const Surface monkey{Graph(
                [](const std::size_t i, const std::size_t j) {
                    return Square.at(i) * Linear.at(j) - Linear.at(i) * Square.at(j);
                },
                0.0, 1.0, weights)};
            for (const bool swapped : {false, true}) {
                INFO("swapped " << swapped);
                const auto result = swapped ? Intersect(flat, Box2::Everything(), monkey,
                                                        Box2::Everything(), tolerance, ctx)
                                            : Intersect(monkey, Box2::Everything(), flat,
                                                        Box2::Everything(), tolerance, ctx);
                INFO(IntersectTest::ErrorText(result));
                REQUIRE_FALSE(result.has_value());
                REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
            }
        }
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

// --- Nodes on knot lines and seams -----------------------------------------

namespace {

/// Two full turns of unit cylinders, about z and about x: the ellipses
/// x = z and x = -z, crossing at (0, 1, 0) and (0, -1, 0). The result is
/// the whole of both and nothing twice (its length is theirs), every sample
/// on both cylinders, in one connected curve whose branches end, four each,
/// at the two nodes.
void CheckFullTurns(const Stereon::SurfaceSurfaceIntersection& r, const Surface& a,
                    const Surface& b, const Motion& motion, const double resolution)
{
    REQUIRE(r.Overlaps.empty());
    REQUIRE(r.Points.empty());
    const IntersectTest::Connectivity shape = IntersectTest::ConnectivityOf(r, 1e-9 * motion.Scale);
    REQUIRE(shape.Curves == 1);
    REQUIRE(shape.Nodes == 2);
    const double reach = 4.0 * resolution / motion.Scale;
    double length = 0.0;
    std::array<std::size_t, 2> endsAt{};
    for (const Stereon::IntersectionBranch& branch : r.Branches) {
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= resolution);
        REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 64) <=
                branch.Bound + 1e-12 * motion.Scale);
        const std::vector<double> t = IntersectTest::SamplesOf(branch.Range, 1024);
        Point3 previous = motion.Back(Eval(branch.Curve3d, t.front()));
        for (std::size_t k = 1; k < t.size(); ++k) {
            const Point3 p = motion.Back(Eval(branch.Curve3d, t[k]));
            REQUIRE(std::abs(p.X * p.X + p.Y * p.Y - 1.0) <= reach);
            REQUIRE(std::abs(p.Y * p.Y + p.Z * p.Z - 1.0) <= reach);
            length += Distance(previous, p);
            previous = p;
        }
        for (const double end : {t.front(), t.back()}) {
            const Point3 p = motion.Back(Eval(branch.Curve3d, end));
            endsAt.at(0) += Distance(p, Point3{0.0, 1.0, 0.0}) <= reach ? 1U : 0U;
            endsAt.at(1) += Distance(p, Point3{0.0, -1.0, 0.0}) <= reach ? 1U : 0U;
        }
    }
    REQUIRE(endsAt.at(0) == 4);
    REQUIRE(endsAt.at(1) == 4);
    REQUIRE(std::abs(length - 2.0 * EllipseLength()) <= 1e-4);
}

void CheckFullTurnsMoved(const BSplineSurface& sa, const BSplineSurface& sb, const Motion& motion)
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
        CheckFullTurns(*result, swapped ? b : a, swapped ? a : b, motion, scaled.GetResolution());
    }
}

void RequireFailsClosed(const Surface& a, const Surface& b, const double within)
{
    for (const bool swapped : {false, true}) {
        INFO("swapped " << swapped);
        const auto result =
            swapped ? Intersect(b, Box2::Everything(), a, Box2::Everything(), within, ctx)
                    : Intersect(a, Box2::Everything(), b, Box2::Everything(), within, ctx);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}

} // namespace

TEST_CASE("a node on a knot line is split there as one inside a patch is",
          "[intersect][generic][node][knot]")
{
    // Before the patches around such a node were settled together, each of
    // these failed closed: no single patch holds the node inside it.
    SECTION("polynomial: a knot line of one surface through the node")
    {
        CheckMoved(KnottedQuadric(false), Bowl(), QuadricCrossing(), {});
    }
    SECTION("polynomial: two knot lines crossing at the node")
    {
        CheckMoved(KnottedQuadric(true), Bowl(), QuadricCrossing(), {});
    }
    SECTION("polynomial: knot lines of both surfaces through the node")
    {
        // The bowl over [-1/2, 3/2]^2 in four pieces meeting above the node.
        const std::vector<double> wide{-0.5, 0.5, 1.5};
        const std::vector<double> knots{0.0, 0.625, 1.0};
        CheckMoved(KnottedQuadric(true), Pieced({1.0, 1.0}, wide, knots, wide, knots),
                   QuadricCrossing(), {});
    }
    SECTION("rational: exact cylinders in two arcs each")
    {
        CheckMoved(UprightInTwo(), Lying(), CylinderCrossing(), {});
        CheckMoved(Upright(), LyingInTwo(), CylinderCrossing(), {});
        CheckMoved(UprightInTwo(), LyingInTwo(), CylinderCrossing(), {});
    }
}

TEST_CASE("nodes on a knot line and on a seam of two full cylinders are both split",
          "[intersect][generic][node][knot][seam]")
{
    // (0, 1, 0) is on a knot line of the cylinder about z and on the seam of
    // the one about x; (0, -1, 0) on a knot line of each. This is
    // corpus/nurbs-ssi/cylinder-cylinder/node-origin-x1.
    CheckFullTurnsMoved(FullTurn(false), FullTurn(true), {});
}

TEST_CASE("answers at nodes on knot lines and seams move with the surfaces",
          "[intersect][generic][node][knot][seam][metamorphic]")
{
    const BSplineSurface a = KnottedQuadric(true);
    const BSplineSurface b = Bowl();
    SECTION("MR-01: translations, exact in binary")
    {
        for (const Point3 shift : {Point3{16.0, -32.0, 8.0}, Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(a, b, QuadricCrossing(), {1.0, shift});
            CheckMoved(UprightInTwo(), LyingInTwo(), CylinderCrossing(), {1.0, shift});
        }
        CheckFullTurnsMoved(FullTurn(false), FullTurn(true), {1.0, Point3{16.0, -32.0, 8.0}});
    }
    SECTION("MR-03: scaling by powers of two, with the tolerance and resolution")
    {
        for (const int k : {-6, 10}) {
            const Motion scaled{std::ldexp(1.0, k), Point3{0.0, 0.0, 0.0}};
            CheckMoved(a, b, QuadricCrossing(), scaled);
            CheckMoved(UprightInTwo(), LyingInTwo(), CylinderCrossing(), scaled);
        }
        CheckFullTurnsMoved(FullTurn(false), FullTurn(true),
                            {std::ldexp(1.0, 5), Point3{0.0, 0.0, 0.0}});
    }
    SECTION("MR-04: reversed and transposed parameters leave the curves in place")
    {
        CheckMoved(ReversedU(a), b, QuadricCrossing(), {});
        CheckMoved(Transposed(a), Transposed(b), QuadricCrossing(), {});
        CheckMoved(Transposed(ReversedU(UprightInTwo())), ReversedU(LyingInTwo()),
                   CylinderCrossing(), {});
        CheckFullTurnsMoved(ReversedU(FullTurn(false)), Transposed(FullTurn(true)), {});
    }
    // MR-02, the exchange of the surfaces, is in every CheckMoved.
}

TEST_CASE("near a node on a knot line, only what the tolerance admits is a node",
          "[intersect][generic][node][knot][near-event]")
{
    const Surface bowl{Bowl()};
    SECTION("a crease along the knot line, within the tolerance: a node")
    {
        // The slope jumps by 2^-39 across x = 1/2: the branches bend within
        // 2^-40 of the smooth surface's, far inside the resolution.
        CheckMoved(KnottedQuadric(false, 0x1p-40), Bowl(), QuadricCrossing(), {});
    }
    SECTION("a crease beyond the tolerance: two branches touching on it, not a node")
    {
        // x'^2 + |x'| / 128 = 4 y'^2: a parabola on each side of the crease,
        // tangent to it at the centre. Not a transverse node; fails closed.
        RequireFailsClosed(Surface{KnottedQuadric(false, 0x1p-7)}, bowl, tolerance);
    }
    SECTION("a gap within the tolerance is a node; beyond it, two branches")
    {
        for (const double height : {0x1p-30, 0x1p-6}) {
            INFO("height " << height);
            const auto lift = [&](const Point3& p) { return Point3{p.X, p.Y, p.Z + height}; };
            const Surface lifted{Mapped(Bowl(), lift)};
            const Surface knotted{KnottedQuadric(true)};
            const auto result =
                Intersect(knotted, Box2::Everything(), lifted, Box2::Everything(), tolerance, ctx);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->Points.empty());
            const std::size_t curves = height < tolerance ? 1 : 2;
            REQUIRE(result->Branches.size() == (curves == 1 ? 4U : 2U));
            REQUIRE(IntersectTest::ConnectivityOf(*result, 1e-9).Curves == curves);
            for (const Stereon::IntersectionBranch& branch : result->Branches) {
                REQUIRE(branch.IsCertified());
                REQUIRE(branch.Bound <= ctx.GetResolution());
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    // x'^2 - 4 y'^2 = height, on the lifted bowl.
                    const Point3 p = Eval(branch.Curve3d, t);
                    const double x = p.X - 0.5;
                    const double y = p.Y - 0.5;
                    REQUIRE(std::abs(p.Z - (x * x + y * y + height)) <= 4.0 * ctx.GetResolution());
                    REQUIRE(std::abs(x * x - 4.0 * y * y - height) <= 4.0 * ctx.GetResolution());
                }
            }
        }
    }
    SECTION("a node beside the knot line is split, or refused where neither proof reaches")
    {
        // The knot line 2^k beside the node. From 2^-12 on the node is proved
        // inside its own patch; to 2^-20 it is within the box given up
        // around the line. Between (2^-14 and 2^-16 were refused when this
        // was written), the box would need the surfaces further apart than
        // the tolerance and the patch has no room: a refusal is admitted
        // there, a wrong answer nowhere.
        for (const int k : {-6, -16, -20, -40}) {
            INFO("knot line 2^" << k << " beside the node");
            const BSplineSurface beside = Pieced({2.0, -3.0}, {0.0, 0.5 + std::ldexp(1.0, k), 1.0},
                                                 KnotOffMiddle, OneSpan, OneSpan);
            if (k != -16) {
                CheckMoved(beside, Bowl(), QuadricCrossing(), {});
                continue;
            }
            const Surface a{beside};
            const auto result =
                Intersect(a, Box2::Everything(), bowl, Box2::Everything(), tolerance, ctx);
            if (result) {
                CheckNode(*result, a, bowl, QuadricCrossing(), {}, ctx.GetResolution());
            } else {
                REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
            }
        }
    }
}

TEST_CASE("what is not a proved transverse node on a knot line or seam fails closed",
          "[intersect][generic][node][knot][seam][special]")
{
    SECTION("no tolerance admits the node")
    {
        RequireFailsClosed(Surface{KnottedQuadric(true)}, Surface{Bowl()}, 0.0);
    }
    SECTION("a tangency along a curve through the knot line: a degenerate form")
    {
        // x'^2 - 3 y'^2 = x'^2 + y'^2: the line y' = 0, doubled.
        RequireFailsClosed(
            Surface{Pieced({1.0, -3.0}, BreakAtNode, KnotOffMiddle, OneSpan, OneSpan)},
            Surface{Bowl()}, tolerance);
    }
    SECTION("a definite touch on the knot line")
    {
        RequireFailsClosed(
            Surface{Pieced({3.0, 2.0}, BreakAtNode, KnotOffMiddle, OneSpan, OneSpan)},
            Surface{Bowl()}, tolerance);
    }
    SECTION("a seam whose two ends differ in the doubles is not a closed chart")
    {
        BSplineSurface open = FullTurn(true);
        for (std::size_t k = open.Poles.size() - 2; k < open.Poles.size(); ++k) {
            open.Poles[k].Z = std::numeric_limits<double>::denorm_min();
        }
        RequireFailsClosed(Surface{FullTurn(false)}, Surface{open}, tolerance);
    }
}
