// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Tangencies along a curve in the generic SSI (ADR-0010): polynomial and
// rational B-spline surfaces touching along a line, a parabola and a circle.
// Each answer is held to a closed form written here, never derived from the
// kernel's path: one tangential branch, its two ends, the curve every sample
// lies on, and a certificate within the pair tolerance. The same answer,
// moved accordingly, under exchange of the surfaces, reversal and
// transposition of parameters, a far translation and scaling by powers of two
// (ADR-0053, MR-01 to MR-04). Near events: a gap within the tolerance either
// way is declared a tangency (ADR-0029 section 6, ADR-0050), never two
// crossings or nothing; beyond it, two crossings or a proved separation.
// What is not proved fails closed: a valley that ends where both surfaces'
// domains end together.

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
/// A pair tolerance ten resolutions wide, so that a tangency declared within
/// it can be told from one certified within the resolution.
constexpr double tolerance = 1e-4;
constexpr double Pi = std::numbers::pi;

// --- Polynomial graphs over the unit square ----------------------------------

double Binomial(const std::size_t n, const std::size_t k)
{
    double c = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        c = c * static_cast<double>(n + 1 - i) / static_cast<double>(i);
    }
    return c;
}

/// The Bernstein coefficients of degree n over t in [0, 1] of (t - 1/2)^p,
/// p <= n: its power coefficients a_r = C(p, r) (-1/2)^(p - r), and t^r has
/// the coefficients C(i, r) / C(n, r).
std::vector<double> BernsteinOfPower(const std::size_t n, const std::size_t p)
{
    std::vector<double> b(n + 1, 0.0);
    for (std::size_t r = 0; r <= p; ++r) {
        const double a = Binomial(p, r) * std::pow(-0.5, static_cast<double>(p - r));
        for (std::size_t i = r; i <= n; ++i) {
            b[i] += a * Binomial(i, r) / Binomial(n, r);
        }
    }
    return b;
}

/// The graph z = sum c[p][q] x'^p y'^q, x' = x - 1/2 and y' = y - 1/2, over
/// the unit square, of degrees (n, m): x = u and y = v, the poles evenly
/// spaced. Weights along v only, when given, reparameterize y and keep the
/// graph: z does not depend on the pole's row j.
BSplineSurface PolynomialGraph(const std::size_t n, const std::size_t m,
                               const std::vector<std::vector<double>>& c,
                               const std::vector<double>& weightsV = {})
{
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= m; ++j) {
            double z = 0.0;
            for (std::size_t p = 0; p < c.size(); ++p) {
                for (std::size_t q = 0; q < c[p].size(); ++q) {
                    z += c[p][q] * BernsteinOfPower(n, p)[i] * BernsteinOfPower(m, q)[j];
                }
            }
            poles.push_back({static_cast<double>(i) / static_cast<double>(n),
                             static_cast<double>(j) / static_cast<double>(m), z});
            if (!weightsV.empty()) {
                weights.push_back(weightsV.at(j));
            }
        }
    }
    const auto knots = [](const std::size_t d) {
        std::vector<double> k(d + 1, 0.0);
        k.insert(k.end(), d + 1, 1.0);
        return k;
    };
    const auto made = Stereon::MakeBSplineSurface(static_cast<int>(n), static_cast<int>(m),
                                                  knots(n), knots(m), poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

/// z = 4 x'^2, linear in y: tangent to z = 0 along x = 1/2.
BSplineSurface Trough(const std::vector<double>& weightsV = {})
{
    return PolynomialGraph(2, 1, {{0.0}, {0.0}, {4.0}}, weightsV);
}

/// z = (y' - x'^2)^2 = y'^2 - 2 x'^2 y' + x'^4: tangent to z = 0 along the
/// parabola y = 1/2 + (x - 1/2)^2.
BSplineSurface ParabolicValley()
{
    return PolynomialGraph(4, 2, {{0.0, 0.0, 1.0}, {}, {0.0, -2.0}, {}, {1.0}});
}

/// z = height over [-1, 3]^2, bilinear: wider than every other surface
/// here, so that a tangency's ends lie on the other surface's boundary.
BSplineSurface Wide(const double height)
{
    const std::vector<Point3> poles{
        {-1.0, -1.0, height}, {-1.0, 3.0, height}, {3.0, -1.0, height}, {3.0, 3.0, height}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

/// z = height over exactly the unit square: the trough's own domain.
BSplineSurface Flush(const double height)
{
    const std::vector<Point3> poles{
        {0.0, 0.0, height}, {0.0, 1.0, height}, {1.0, 0.0, height}, {1.0, 1.0, height}};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(1, 1, knots, knots, poles, {});
    REQUIRE(made.has_value());
    return *made;
}

// --- Rational surfaces -------------------------------------------------------

/// The rational quadratic arc of the unit circle from angle a0 to a1, less
/// than a half turn apart.
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

/// The unit cylinder about the x axis, x in [0, 1], over the arc of (y, z)
/// from -150 to -30 degrees: it lies on z = -1 along y = 0.
BSplineSurface LyingCylinder()
{
    const CircleArc arc = ArcOf(-5.0 * Pi / 6.0, -Pi / 6.0);
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < 3; ++i) {
        for (const double x : {0.0, 1.0}) {
            const auto [y, z] = arc.Poles.at(i);
            poles.push_back({x, y, z});
            weights.push_back(arc.Weights.at(i));
        }
    }
    const std::vector<double> knotsU{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const std::vector<double> knotsV{0.0, 0.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 1, knotsU, knotsV, poles, weights);
    REQUIRE(made.has_value());
    return *made;
}

constexpr double TorusR = 2.0;
constexpr double Torusr = 0.5;

/// A patch of the torus about the z axis, major radius 2 and minor 1/2: a
/// quarter turn about the axis in u, from 0 to 90 degrees, and the minor
/// circle from 30 to 150 degrees in v, over its top. It lies on z = 1/2
/// along the circle x^2 + y^2 = 4. The poles are the products of the two
/// arcs' and the weights the products of theirs: exact.
BSplineSurface TorusPatch()
{
    const CircleArc turn = ArcOf(0.0, Pi / 2.0);
    const CircleArc minor = ArcOf(Pi / 6.0, 5.0 * Pi / 6.0);
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const auto [cx, cy] = turn.Poles.at(i);
            const auto [rho, z] = minor.Poles.at(j);
            const double radius = TorusR + Torusr * rho;
            poles.push_back({radius * cx, radius * cy, Torusr * z});
            weights.push_back(turn.Weights.at(i) * minor.Weights.at(j));
        }
    }
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const auto made = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, weights);
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
    [[nodiscard]] Point3 Forward(const Point3& p) const
    {
        return {Scale * p.X + Shift.X, Scale * p.Y + Shift.Y, Scale * p.Z + Shift.Z};
    }
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

/// What a tangency case must answer, in the input's own space: the ends of
/// its one branch, and how far a point is from the closed-form curve.
struct Expected {
    std::array<Point3, 2> Ends;
    std::function<double(const Point3&)> OffCurve;
};

Expected TroughLine(const double z = 0.0)
{
    return {{Point3{0.5, 0.0, z}, Point3{0.5, 1.0, z}},
            [z](const Point3& p) { return std::max(std::abs(p.X - 0.5), std::abs(p.Z - z)); }};
}

Expected Parabola()
{
    return {{Point3{0.0, 0.75, 0.0}, Point3{1.0, 0.75, 0.0}}, [](const Point3& p) {
                const double x = p.X - 0.5;
                return std::max(std::abs(p.Y - 0.5 - x * x), std::abs(p.Z));
            }};
}

Expected CylinderLine()
{
    return {{Point3{0.0, 0.0, -1.0}, Point3{1.0, 0.0, -1.0}},
            [](const Point3& p) { return std::max(std::abs(p.Y), std::abs(p.Z + 1.0)); }};
}

Expected TorusCircle()
{
    return {{Point3{TorusR, 0.0, Torusr}, Point3{0.0, TorusR, Torusr}}, [](const Point3& p) {
                return std::max(std::abs(std::hypot(p.X, p.Y) - TorusR), std::abs(p.Z - Torusr));
            }};
}

/// The result is the closed form, moved: one tangential branch between the
/// two ends, certified within the pair tolerance on both surfaces, every
/// sample on the curve within the tolerance's reach across the valley. A gap
/// g <= tolerance between the surfaces lets the valley's point lie off the
/// exact curve by sqrt(2 g / h_aa) across it; here g is zero, and the
/// samples are held to a few resolutions.
void CheckTangency(const Stereon::SurfaceSurfaceIntersection& r, const Surface& a, const Surface& b,
                   const Expected& e, const Motion& motion, const double resolution,
                   const double pairTolerance)
{
    REQUIRE(r.Overlaps.empty());
    REQUIRE(r.Points.empty());
    REQUIRE(r.Branches.size() == 1);
    const Stereon::IntersectionBranch& branch = r.Branches.front();
    REQUIRE(branch.Contact == Stereon::ContactType::Tangential);
    REQUIRE(branch.Transition == Stereon::SideTransition::Unresolved);
    REQUIRE(branch.IsCertified());
    REQUIRE(branch.Bound <= std::max(resolution, pairTolerance));
    if (branch.Bound > resolution) {
        REQUIRE(branch.BoundCause == Stereon::ToleranceCause::TolerantCoincidence);
    } else {
        REQUIRE(branch.BoundCause == Stereon::ToleranceCause::CertificateWidening);
    }
    const double largest =
        std::max({std::abs(motion.Shift.X), std::abs(motion.Shift.Y), std::abs(motion.Shift.Z)}) +
        4.0 * motion.Scale;
    const double rounding = 4.0 * std::numeric_limits<double>::epsilon() * largest;
    REQUIRE(IntersectTest::PcurveDeviation(branch, a, b, 64) <=
            branch.Bound + 1e-12 * motion.Scale + rounding);
    const double reach = 4.0 * resolution / motion.Scale; // In the input's own space.
    Point3 first = motion.Back(Eval(branch.Curve3d, branch.Range.GetLower()));
    Point3 last = motion.Back(Eval(branch.Curve3d, branch.Range.GetUpper()));
    if (Distance(first, e.Ends[1]) < Distance(first, e.Ends[0])) {
        std::swap(first, last);
    }
    INFO("branch from (" << first.X << ", " << first.Y << ", " << first.Z << ") to (" << last.X
                         << ", " << last.Y << ", " << last.Z << ")");
    REQUIRE(Distance(first, e.Ends[0]) <= reach);
    REQUIRE(Distance(last, e.Ends[1]) <= reach);
    for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
        REQUIRE(e.OffCurve(motion.Back(Eval(branch.Curve3d, t))) <= reach);
    }
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
        CheckTangency(*result, swapped ? b : a, swapped ? a : b, e, motion, scaled.GetResolution(),
                      tolerance * motion.Scale);
    }
}

} // namespace

TEST_CASE("a trough lying on a plane is one tangential branch along its bottom line",
          "[intersect][generic][tangent]")
{
    // Before tangencies along a curve were traced, this failed closed.
    CheckMoved(Trough(), Wide(0.0), TroughLine(), {});
}

TEST_CASE("a tangency along a parabola is traced across the patch", "[intersect][generic][tangent]")
{
    CheckMoved(ParabolicValley(), Wide(0.0), Parabola(), {});
}

TEST_CASE("rational surfaces tangent along a line and along a circle",
          "[intersect][generic][tangent][rational]")
{
    SECTION("the trough with its y reparameterized by weights along v")
    {
        CheckMoved(Trough({1.0, 3.0}), Wide(0.0), TroughLine(), {});
    }
    SECTION("an exact cylinder lying on a plane")
    {
        CheckMoved(LyingCylinder(), Wide(-1.0), CylinderLine(), {});
    }
    SECTION("an exact torus patch lying on a plane, along a quarter circle")
    {
        CheckMoved(TorusPatch(), Wide(Torusr), TorusCircle(), {});
    }
}

TEST_CASE("tangency answers move with the surfaces", "[intersect][generic][tangent][metamorphic]")
{
    SECTION("MR-01: translations, exact in binary")
    {
        for (const Point3 shift : {Point3{16.0, -32.0, 8.0}, Point3{1024.0, -2048.0, 512.0},
                                   Point3{-4096.0, 256.0, 1024.0}}) {
            CheckMoved(Trough(), Wide(0.0), TroughLine(), {1.0, shift});
            CheckMoved(TorusPatch(), Wide(Torusr), TorusCircle(), {1.0, shift});
        }
    }
    SECTION("MR-03: scaling by powers of two, with the tolerance and resolution")
    {
        for (const int k : {-6, 3, 10}) {
            CheckMoved(ParabolicValley(), Wide(0.0), Parabola(), {std::ldexp(1.0, k), Point3{}});
            CheckMoved(LyingCylinder(), Wide(-1.0), CylinderLine(), {std::ldexp(1.0, k), Point3{}});
        }
    }
    SECTION("MR-04: reversed and transposed parameters leave the curves in place")
    {
        CheckMoved(ReversedU(Trough()), Wide(0.0), TroughLine(), {});
        CheckMoved(Transposed(Trough()), Wide(0.0), TroughLine(), {});
        CheckMoved(Transposed(ReversedU(ParabolicValley())), Transposed(Wide(0.0)), Parabola(), {});
        CheckMoved(ReversedU(Transposed(TorusPatch())), ReversedU(Wide(Torusr)), TorusCircle(), {});
    }
    // MR-02, the exchange of the surfaces, is in every CheckMoved.
}

TEST_CASE("a gap within the tolerance is declared a tangency; beyond it, crossings or nothing",
          "[intersect][generic][tangent][near-event]")
{
    const Surface trough{Trough()};
    for (const bool swapped : {false, true}) {
        INFO((swapped ? "swapped" : "in order"));
        const auto run = [&](const double height) {
            const Surface flat{Wide(height)};
            return swapped ? Intersect(flat, Box2::Everything(), trough, Box2::Everything(),
                                       tolerance, ctx)
                           : Intersect(trough, Box2::Everything(), flat, Box2::Everything(),
                                       tolerance, ctx);
        };
        SECTION("within the tolerance, above and below: a declared tangency")
        {
            for (const double height : {0.5 * tolerance, -0.5 * tolerance}) {
                INFO("height " << height);
                const auto result = run(height);
                INFO(IntersectTest::ErrorText(result));
                REQUIRE(result.has_value());
                REQUIRE(result->Branches.size() == 1);
                const Stereon::IntersectionBranch& branch = result->Branches.front();
                REQUIRE(branch.Contact == Stereon::ContactType::Tangential);
                REQUIRE(branch.IsCertified());
                REQUIRE(branch.Bound <= tolerance);
                // The valley of either surface's height over the other is
                // x = 1/2, where both normals are the z axis; its point is
                // on one of the two, at z = 0 or at the height.
                for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                    const Point3 p = Eval(branch.Curve3d, t);
                    REQUIRE(std::abs(p.X - 0.5) <= 4.0 * ctx.GetResolution());
                    REQUIRE(std::min(std::abs(p.Z), std::abs(p.Z - height)) <=
                            4.0 * ctx.GetResolution());
                }
                IntersectTest::CheckToleranceEvidence(*result, tolerance);
            }
        }
        SECTION("four tolerances below: proved apart")
        {
            const auto result = run(-4.0 * tolerance);
            INFO(IntersectTest::ErrorText(result));
            REQUIRE(result.has_value());
            REQUIRE(result->IsEmpty());
        }
        SECTION("four tolerances above: two crossing lines, or a refusal, never a tangency")
        {
            const double height = 4.0 * tolerance;
            const auto result = run(height);
            INFO(IntersectTest::ErrorText(result));
            if (result) {
                REQUIRE(result->Branches.size() == 2);
                const double x = 0.5 * std::sqrt(height);
                for (const Stereon::IntersectionBranch& branch : result->Branches) {
                    REQUIRE(branch.Contact == Stereon::ContactType::Crossing);
                    for (const double t : IntersectTest::SamplesOf(branch.Range, 64)) {
                        const Point3 p = Eval(branch.Curve3d, t);
                        REQUIRE(std::abs(std::abs(p.X - 0.5) - x) <= 4.0 * ctx.GetResolution());
                        REQUIRE(std::abs(p.Z - height) <= 4.0 * ctx.GetResolution());
                    }
                }
            } else {
                REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
            }
        }
    }
}

TEST_CASE("a tangency whose ends lie where both domains end fails closed",
          "[intersect][generic][tangent][special]")
{
    // The valley's ends are on the boundary of both surfaces at once: the
    // plane's height over the trough is not defined past its ends, nor the
    // trough's over the plane, so neither is proved a graph over the other
    // there.
    const Surface trough{Trough()};
    const Surface flat{Flush(0.0)};
    for (const bool swapped : {false, true}) {
        const auto result = swapped ? Intersect(flat, trough, ctx) : Intersect(trough, flat, ctx);
        INFO(IntersectTest::ErrorText(result));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().GetCode() == ErrorCode::IntersectionFailed);
    }
}
