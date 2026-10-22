// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"

/// Cylinder–cylinder with skew or intersecting axes (ADR-0010's revision of
/// 2026-10-20, the Boolean prototype's M3b). On cylinder A at angle u, B's
/// equation is a quadratic in the height v, alpha v² + 2 beta(u) v +
/// gamma(u) = 0, whose discriminant D(u) = beta² - alpha gamma is a
/// trigonometric polynomial of degree 2 in u, a quartic P(t) in t =
/// tan(u/2). The curve is the two branches v±(u) where D >= 0, meeting at
/// the branch points where D vanishes.
///
/// Decisions are exact (RS-0001 G5): the axes' distance against the sum
/// and difference of the radii, degree-10 polynomials in the frames
/// through DynamicExpansion, settle whether the cylinders are apart, touch
/// at a point, or are tangent along the curve (the singular case, which
/// is not supported yet); with multiple roots thereby excluded, a Sturm
/// sequence of P counts the branch points and bisection with exact signs
/// finds them. Each loop of the curve is parametrised smoothly, u = m -
/// h cos s over a branch interval, so the fold points are regular, then
/// sampled and interpolated by B-splines: a pcurve on each cylinder and a
/// cached 3D curve on the shared parameter, with the samples doubled until
/// CertifyDeviation certifies both pcurves within the resolution.

namespace Stereon::Detail {

namespace {

constexpr double TwoPi = 2.0 * std::numbers::pi;

// --- Exact vectors --------------------------------------------------------

using DVec = std::array<DynamicExpansion, 3>;

[[nodiscard]] DVec Exact(const Vec3& v)
{
    return {DynamicExpansion{v.X}, DynamicExpansion{v.Y}, DynamicExpansion{v.Z}};
}

[[nodiscard]] DVec ExactOffset(const Point3& p, const Point3& q)
{
    return {TwoDiff(p.X, q.X), TwoDiff(p.Y, q.Y), TwoDiff(p.Z, q.Z)};
}

[[nodiscard]] DynamicExpansion Dot(const DVec& a, const DVec& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[nodiscard]] DVec Cross(const DVec& a, const DVec& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[nodiscard]] DVec Negated(const DVec& a)
{
    return {-a[0], -a[1], -a[2]};
}

/// The cylinder's axis as the frame gives it, ±X x Y, exactly.
[[nodiscard]] DVec ExactAxis(const Cylinder& c)
{
    const DVec cross = Cross(Exact(c.XAxis), Exact(c.YAxis));
    return c.RightHanded ? cross : Negated(cross);
}

// --- Trigonometric polynomials ----------------------------------------

/// A polynomial of degree at most 2 in c = cos u and s = sin u: the
/// coefficients of 1, c, s, c², cs, s².
struct TrigPoly {
    std::array<DynamicExpansion, 6> K;
};

[[nodiscard]] TrigPoly Linear(DynamicExpansion k0, DynamicExpansion kc, DynamicExpansion ks)
{
    return {{std::move(k0), std::move(kc), std::move(ks), {}, {}, {}}};
}

[[nodiscard]] TrigPoly Sum(const TrigPoly& a, const TrigPoly& b)
{
    TrigPoly r;
    for (std::size_t i = 0; i < 6; ++i) {
        r.K[i] = a.K[i] + b.K[i];
    }
    return r;
}

[[nodiscard]] TrigPoly Scaled(const TrigPoly& a, const DynamicExpansion& f)
{
    TrigPoly r;
    for (std::size_t i = 0; i < 6; ++i) {
        r.K[i] = a.K[i] * f;
    }
    return r;
}

/// The product of two polynomials of degree at most 1.
[[nodiscard]] TrigPoly Product(const TrigPoly& a, const TrigPoly& b)
{
    return {{a.K[0] * b.K[0], a.K[0] * b.K[1] + a.K[1] * b.K[0], a.K[0] * b.K[2] + a.K[2] * b.K[0],
             a.K[1] * b.K[1], a.K[1] * b.K[2] + a.K[2] * b.K[1], a.K[2] * b.K[2]}};
}

/// The exact sign at u = pi: c = -1, s = 0.
[[nodiscard]] Sign SignAtPi(const TrigPoly& p)
{
    return (p.K[0] - p.K[1] + p.K[3]).GetSign();
}

/// P(t) = p(c(t), s(t)) (1 + t²)² with c = (1 - t²)/(1 + t²) and s =
/// 2t/(1 + t²): each monomial's image is a polynomial in t with small
/// integer coefficients.
[[nodiscard]] ExactPolynomial Quartic(const TrigPoly& p)
{
    const auto poly = [](std::initializer_list<double> coefficients) {
        std::vector<DynamicExpansion> c;
        for (const double x : coefficients) {
            c.emplace_back(x);
        }
        return ExactPolynomial(std::move(c));
    };
    const std::array<ExactPolynomial, 6> images{
        poly({1.0, 0.0, 2.0, 0.0, 1.0}),  // 1: (1 + t²)²
        poly({1.0, 0.0, 0.0, 0.0, -1.0}), // c: (1 - t²)(1 + t²)
        poly({0.0, 2.0, 0.0, 2.0}),       // s: 2t(1 + t²)
        poly({1.0, 0.0, -2.0, 0.0, 1.0}), // c²: (1 - t²)²
        poly({0.0, 2.0, 0.0, -2.0}),      // cs: 2t(1 - t²)
        poly({0.0, 0.0, 4.0}),            // s²: 4t²
    };
    ExactPolynomial sum;
    for (std::size_t i = 0; i < 6; ++i) {
        sum = Sum(sum, Scaled(images[i], p.K[i]));
    }
    return sum;
}

// --- The pair ---------------------------------------------------------------

/// Both cylinders with their axes as the frames give them, in doubles for
/// the geometry and exactly for the decisions.
struct Pair {
    const Surface& ErasedA;
    const Surface& ErasedB;
    const Cylinder& A;
    const Cylinder& B;
    Vec3 AxisA;
    Vec3 AxisB;
    /// |a x b|² = |a|²|b|² - (a.b)², the quadratic's leading coefficient.
    double Alpha;
};

[[nodiscard]] Vec3 FrameAxis(const Cylinder& c)
{
    const Vec3 cross = Cross(c.XAxis, c.YAxis);
    return c.RightHanded ? cross : -cross;
}

[[nodiscard]] Pair PairOf(const Surface& erasedA, const Surface& erasedB, const Cylinder& a,
                          const Cylinder& b)
{
    const Vec3 axisA = FrameAxis(a);
    const Vec3 axisB = FrameAxis(b);
    return {erasedA, erasedB, a, b, axisA, axisB, LengthSquared(Cross(axisA, axisB))};
}

/// beta(u), gamma(u) and D(u) in doubles, for the points once the
/// configuration is decided: w = cA - cB + rA (cos u XA + sin u YA).
struct Quadratic {
    double Beta;
    double Gamma;
    double Discriminant;
};

[[nodiscard]] Quadratic QuadraticAt(const Pair& p, const double u)
{
    const Vec3 w = (p.A.Origin - p.B.Origin) + RadialTurned(p.A, u, 0);
    const double ab = Dot(p.AxisA, p.AxisB);
    const double bb = LengthSquared(p.AxisB);
    const double wa = Dot(w, p.AxisA);
    const double wb = Dot(w, p.AxisB);
    const double beta = wa * bb - wb * ab;
    const double gamma = LengthSquared(w) * bb - wb * wb - p.B.Radius * p.B.Radius * bb;
    return {beta, gamma, beta * beta - p.Alpha * gamma};
}

/// The exact discriminant D(u) = beta(u)² - alpha gamma(u).
[[nodiscard]] TrigPoly ExactDiscriminant(const Pair& p)
{
    const DVec a = ExactAxis(p.A);
    const DVec b = ExactAxis(p.B);
    const DVec w0 = ExactOffset(p.A.Origin, p.B.Origin);
    const DVec x = Exact(p.A.XAxis);
    const DVec y = Exact(p.A.YAxis);
    const double rA = p.A.Radius;
    const DynamicExpansion ab = Dot(a, b);
    const DynamicExpansion bb = Dot(b, b);
    const DynamicExpansion alpha = Dot(a, a) * bb - ab * ab;
    // The components of w, each linear in (c, s), and w.a, w.b.
    std::array<TrigPoly, 3> w;
    for (std::size_t i = 0; i < 3; ++i) {
        w[i] = Linear(w0[i], x[i] * rA, y[i] * rA);
    }
    const TrigPoly wa = Linear(Dot(w0, a), Dot(x, a) * rA, Dot(y, a) * rA);
    const TrigPoly wb = Linear(Dot(w0, b), Dot(x, b) * rA, Dot(y, b) * rA);
    const TrigPoly beta = Sum(Scaled(wa, bb), Scaled(wb, -ab));
    TrigPoly ww = Product(w[0], w[0]);
    ww = Sum(ww, Product(w[1], w[1]));
    ww = Sum(ww, Product(w[2], w[2]));
    const DynamicExpansion rB2 = TwoProduct(p.B.Radius, p.B.Radius);
    const TrigPoly gamma = Sum(Sum(Scaled(ww, bb), Scaled(Product(wb, wb), DynamicExpansion{-1.0})),
                               Linear(-(rB2 * bb), {}, {}));
    return Sum(Product(beta, beta), Scaled(gamma, -alpha));
}

// --- The configuration -----------------------------------------------------

enum class Configuration : std::uint8_t {
    /// Further apart than the radii's sum: proved separate.
    Apart,
    /// The axes' distance is exactly the radii's sum: one contact point.
    TouchingOutside,
    /// The axes' distance is exactly the radii's difference: the curve is
    /// singular (the surfaces are tangent along it).
    TouchingInside,
    /// Transverse curves.
    Crossing,
};

/// The axes' distance against the radii, exactly: with n = a x b and w =
/// cA - cB, the distance is |w.n|/|n|, so (w.n)² is compared with
/// (rA ± rB)² |n|², where |n|² = alpha.
[[nodiscard]] Configuration ConfigurationOf(const Pair& p)
{
    const DVec a = ExactAxis(p.A);
    const DVec b = ExactAxis(p.B);
    const DVec n = Cross(a, b);
    const DynamicExpansion alpha = Dot(n, n);
    const DynamicExpansion wn = Dot(ExactOffset(p.A.Origin, p.B.Origin), n);
    const DynamicExpansion squared = wn * wn;
    const DynamicExpansion sum = TwoSum(p.A.Radius, p.B.Radius);
    const DynamicExpansion difference = TwoDiff(p.A.Radius, p.B.Radius);
    const Sign outside = (squared - sum * sum * alpha).GetSign();
    if (outside == Sign::Positive) {
        return Configuration::Apart;
    }
    if (outside == Sign::Zero) {
        return Configuration::TouchingOutside;
    }
    const Sign inside = (squared - difference * difference * alpha).GetSign();
    return inside == Sign::Zero ? Configuration::TouchingInside : Configuration::Crossing;
}

/// The one contact point of cylinders touching from outside: on the common
/// perpendicular of the axes, rA from A's. Its bound is the gap between the
/// surfaces at its parameters, by interval evaluation.
[[nodiscard]] IntersectionPoint ContactPoint(const Pair& p)
{
    const Vec3 w = p.A.Origin - p.B.Origin;
    const double aa = LengthSquared(p.AxisA);
    const double bb = LengthSquared(p.AxisB);
    const double ab = Dot(p.AxisA, p.AxisB);
    const double wa = Dot(w, p.AxisA);
    const double wb = Dot(w, p.AxisB);
    const double tA = (ab * wb - bb * wa) / p.Alpha;
    const double tB = (aa * wb - ab * wa) / p.Alpha;
    const Point3 footA = p.A.Origin + tA * p.AxisA;
    const Point3 footB = p.B.Origin + tB * p.AxisB;
    IntersectionPoint point;
    point.Position = footA + p.A.Radius * Normalized(footB - footA);
    point.OnA = Project(p.A, point.Position);
    point.OnB = Project(p.B, point.Position);
    using I = IntervalScalar;
    const BasicCylinder<I> wideA{{I{p.A.Origin.X}, I{p.A.Origin.Y}, I{p.A.Origin.Z}},
                                 I{p.A.Radius},
                                 {I{p.A.XAxis.X}, I{p.A.XAxis.Y}, I{p.A.XAxis.Z}},
                                 {I{p.A.YAxis.X}, I{p.A.YAxis.Y}, I{p.A.YAxis.Z}},
                                 p.A.RightHanded};
    const BasicCylinder<I> wideB{{I{p.B.Origin.X}, I{p.B.Origin.Y}, I{p.B.Origin.Z}},
                                 I{p.B.Radius},
                                 {I{p.B.XAxis.X}, I{p.B.XAxis.Y}, I{p.B.XAxis.Z}},
                                 {I{p.B.YAxis.X}, I{p.B.YAxis.Y}, I{p.B.YAxis.Z}},
                                 p.B.RightHanded};
    const BasicVec3<I> gap = Eval(wideA, BasicPoint2<I>{I{point.OnA.X}, I{point.OnA.Y}}) -
                             Eval(wideB, BasicPoint2<I>{I{point.OnB.X}, I{point.OnB.Y}});
    const Interval squares = Squared(gap.X.Value) + Squared(gap.Y.Value) + Squared(gap.Z.Value);
    point.Bound = std::nextafter(std::sqrt(squares.GetUpper()), Infinity);
    return point;
}

// --- The branch points --------------------------------------------------

/// The real roots of P in (lo, hi], known to hold exactly one, bisected
/// with exact signs to a double.
[[nodiscard]] double Refine(const ExactPolynomial& p, double lo, double hi)
{
    const Sign atLow = p.SignAt(lo);
    if (atLow == Sign::Zero) {
        return lo;
    }
    for (int i = 0; i < 2000; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (!(lo < mid && mid < hi)) {
            break;
        }
        const Sign atMid = p.SignAt(mid);
        if (atMid == Sign::Zero) {
            return mid;
        }
        (atMid == atLow ? lo : hi) = mid;
    }
    return hi;
}

/// The roots of P, each isolated by Sturm counts and refined: P has only
/// simple roots once tangency is excluded.
[[nodiscard]] Result<std::vector<double>> RealRoots(const ExactPolynomial& p)
{
    std::vector<double> roots;
    const SturmSequence sturm(p);
    const int total = sturm.CountRoots(-Infinity, Infinity);
    if (total == 0) {
        return roots;
    }
    double bound = 1.0;
    while (sturm.CountRoots(-bound, bound) < total && bound < 1e300) {
        bound *= 2.0;
    }
    std::vector<std::pair<double, double>> work{{-bound, bound}};
    while (!work.empty()) {
        const auto [lo, hi] = work.back();
        work.pop_back();
        const int count = sturm.CountRoots(lo, hi);
        if (count == 0) {
            continue;
        }
        if (count == 1) {
            roots.push_back(Refine(p, lo, hi));
            continue;
        }
        const double mid = 0.5 * (lo + hi);
        if (!(lo < mid && mid < hi)) {
            return Fail(ErrorCode::IntersectionFailed,
                        "SkewCylinders: two branch points closer than doubles can tell apart");
        }
        work.emplace_back(mid, hi);
        work.emplace_back(lo, mid);
    }
    return roots;
}

/// The branch points on A as angles in [0, 2 pi), sorted: the roots of D,
/// through t = tan(u/2), and u = pi itself when D vanishes there exactly.
[[nodiscard]] Result<std::vector<double>> BranchPoints(const TrigPoly& d,
                                                       const ExactPolynomial& quartic)
{
    STEREON_TRY(const std::vector<double> ts, RealRoots(quartic));
    std::vector<double> us;
    for (const double t : ts) {
        const double u = 2.0 * std::atan(t);
        us.push_back(u < 0.0 ? u + TwoPi : u);
    }
    if (SignAtPi(d) == Sign::Zero) {
        us.push_back(std::numbers::pi);
    }
    std::ranges::sort(us);
    if (us.size() % 2 != 0) {
        return Fail(ErrorCode::IntersectionFailed,
                    "SkewCylinders: an odd number of branch points, which simple roots of the "
                    "discriminant cannot give");
    }
    return us;
}

/// The exact sign of D at an angle: through the quartic at tan(u/2), or
/// at pi itself.
[[nodiscard]] Sign SignOfDiscriminant(const TrigPoly& d, const ExactPolynomial& quartic,
                                      const double u)
{
    const double reduced = u >= std::numbers::pi ? u - TwoPi : u;
    if (reduced == std::numbers::pi || reduced == -std::numbers::pi) {
        return SignAtPi(d);
    }
    return quartic.SignAt(std::tan(0.5 * reduced));
}

// --- The loops -----------------------------------------------------------------

/// One closed loop of the curve: over a branch interval [U1, U2] of A's
/// angle, parametrised u = m - h cos s, or over A's whole turn on one of
/// the branches v± (Side).
struct Loop {
    double U1 = 0.0;
    double U2 = 0.0;
    bool Whole = false;
    double Side = 1.0;
};

/// The loops from the branch points: the intervals where D > 0 (an interval
/// may run past 2 pi), or both branches over the whole turn when D never
/// vanishes and is positive.
[[nodiscard]] std::vector<Loop> LoopsOf(const TrigPoly& d, const ExactPolynomial& quartic,
                                        const std::vector<double>& branchPoints)
{
    std::vector<Loop> loops;
    if (branchPoints.empty()) {
        if (quartic.SignAt(0.0) == Sign::Positive) {
            loops.push_back({0.0, TwoPi, true, 1.0});
            loops.push_back({0.0, TwoPi, true, -1.0});
        }
        return loops;
    }
    for (std::size_t i = 0; i < branchPoints.size(); ++i) {
        const double u1 = branchPoints[i];
        const double u2 =
            i + 1 < branchPoints.size() ? branchPoints[i + 1] : branchPoints[0] + TwoPi;
        if (SignOfDiscriminant(d, quartic, 0.5 * (u1 + u2)) == Sign::Positive) {
            loops.push_back({u1, u2, false, 1.0});
        }
    }
    return loops;
}

/// The loop's point at its parameter s in [0, 2 pi): the angle, the height
/// from the quadratic's root on the right side, and the 3D point.
struct Sample {
    Point2 OnA;
    Point2 OnB;
    Point3 Position;
};

[[nodiscard]] Sample SampleAt(const Pair& p, const Loop& loop, const double s)
{
    const double m = 0.5 * (loop.U1 + loop.U2);
    const double h = 0.5 * (loop.U2 - loop.U1);
    const double u = loop.Whole ? s : m - h * Cos(s);
    const double side = loop.Whole ? loop.Side : (Sin(s) < 0.0 ? -1.0 : 1.0);
    const Quadratic q = QuadraticAt(p, u);
    const double v = (-q.Beta + side * Sqrt(std::max(q.Discriminant, 0.0))) / p.Alpha;
    const Point3 position = Eval(p.A, Point2{u, v});
    return {Point2{u, v}, Project(p.B, position), position};
}

/// The samples of a loop at s = 2 pi k / n, k = 0..n, the last the first
/// again; B's angle unwrapped so that it runs continuously, past the seam
/// where the loop crosses it, and the last sample shows how often the loop
/// winds around B.
struct Samples {
    std::vector<Point2> OnA;
    std::vector<Point2> OnB;
    std::vector<Point3> Positions;
};

[[nodiscard]] Samples SampleLoop(const Pair& p, const Loop& loop, const std::size_t n)
{
    Samples samples;
    double previousB = 0.0;
    for (std::size_t k = 0; k <= n; ++k) {
        const double s = TwoPi * static_cast<double>(k % n) / static_cast<double>(n);
        Sample sample = SampleAt(p, loop, s);
        if (k > 0) {
            sample.OnB.X += TwoPi * std::round((previousB - sample.OnB.X) / TwoPi);
        }
        if (k == n && loop.Whole) {
            sample.OnA.X += TwoPi; // The whole turn, closed by the seam again.
        }
        previousB = sample.OnB.X;
        samples.OnA.push_back(sample.OnA);
        samples.OnB.push_back(sample.OnB);
        samples.Positions.push_back(sample.Position);
    }
    return samples;
}

/// A pcurve through samples on the parameter s in [0, 2 pi]: periodic when
/// the loop closes in the parameter plane, open from the first sample to
/// its image a period on when the loop winds around the surface.
[[nodiscard]] Result<Curve2d> FitPcurve(const std::vector<Point2>& samples)
{
    const bool winds = std::abs(samples.back().X - samples.front().X) > std::numbers::pi;
    const Interval range{0.0, TwoPi};
    if (winds) {
        STEREON_TRY(const CurveFit2d fit, Interpolate(samples, 3, Parametrization::Uniform));
        STEREON_TRY(BSplineCurve2d curve, Reparametrized(fit.Curve, range));
        return Curve2d(std::move(curve));
    }
    const std::span<const Point2> once(samples.data(), samples.size() - 1);
    STEREON_TRY(const CurveFit2d fit, InterpolatePeriodic(once, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d curve, Reparametrized(fit.Curve, range));
    return Curve2d(std::move(curve));
}

[[nodiscard]] Result<Curve> FitCache(const std::vector<Point3>& samples)
{
    const std::span<const Point3> once(samples.data(), samples.size() - 1);
    STEREON_TRY(const CurveFit fit, InterpolatePeriodic(once, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve curve, Reparametrized(fit.Curve, Interval{0.0, TwoPi}));
    return Curve(std::move(curve));
}

/// Whether a parameter range may meet a box's periodic axis: compared one
/// period either way too, since a pcurve may run past the seam.
[[nodiscard]] bool MeetsPeriodic(const Interval& range, const Interval& box)
{
    for (const double shift : {-TwoPi, 0.0, TwoPi}) {
        if (Intersects(Interval{range.GetLower() + shift, range.GetUpper() + shift}, box)) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] bool MayMeet(const Curve2d& pcurve, const Box2& box)
{
    const Box2 hull = BoundingBox(pcurve);
    return MeetsPeriodic(hull.X, box.X) && Intersects(hull.Y, box.Y);
}

/// Whether every sample lies within the model's size box: loops of nearly
/// parallel axes run out of it (RS-0001 Tier B).
[[nodiscard]] bool WithinSizeBox(const Samples& samples)
{
    return std::ranges::all_of(samples.Positions, [](const Point3& position) {
        return Contains(Context::GetSizeBox(), position);
    });
}

/// The samples fitted: both pcurves and the cache, certified against both
/// cylinders within the target. The certification's failure is reported;
/// nothing when the loop lies outside the requested domains.
[[nodiscard]] Result<std::optional<IntersectionBranch>> Fitted(const Pair& p,
                                                               const Samples& samples,
                                                               const Box2& boxA, const Box2& boxB,
                                                               const double target)
{
    STEREON_TRY(Curve2d onA, FitPcurve(samples.OnA));
    STEREON_TRY(Curve2d onB, FitPcurve(samples.OnB));
    if (!MayMeet(onA, boxA) || !MayMeet(onB, boxB)) {
        return std::optional<IntersectionBranch>{};
    }
    STEREON_TRY(Curve cache, FitCache(samples.Positions));
    const Interval range{0.0, TwoPi};
    STEREON_TRY(const DeviationCertificate certA,
                CertifyDeviation(cache, p.ErasedA, onA, range, target));
    STEREON_TRY(const DeviationCertificate certB,
                CertifyDeviation(cache, p.ErasedB, onB, range, target));
    IntersectionBranch branch = Crossing(std::move(cache), std::move(onA), std::move(onB), range);
    branch.Bound = std::max(certA.Bound, certB.Bound);
    return std::optional<IntersectionBranch>{std::move(branch)};
}

/// A loop fitted and certified, the samples doubled until both pcurves are
/// within the resolution; nothing when the loop lies outside the requested
/// domains.
[[nodiscard]] Result<std::optional<IntersectionBranch>>
BranchOf(const Pair& p, const Loop& loop, const Box2& boxA, const Box2& boxB, const Context& ctx)
{
    constexpr std::size_t MostSamples = 4096;
    for (std::size_t n = 64; n <= MostSamples; n *= 2) {
        const Samples samples = SampleLoop(p, loop, n);
        if (!WithinSizeBox(samples)) {
            return Fail(ErrorCode::NotSupported,
                        "SkewCylinders: the intersection leaves the size box; the axes are "
                        "nearly parallel (RS-0001 Tier B)");
        }
        Result<std::optional<IntersectionBranch>> fitted =
            Fitted(p, samples, boxA, boxB, ctx.GetResolution());
        if (fitted || fitted.error().GetCode() != ErrorCode::CertificationFailed ||
            n == MostSamples) {
            return fitted;
        }
    }
    return Fail(ErrorCode::CertificationFailed, "SkewCylinders: the loop could not be certified");
}

} // namespace

Result<SurfaceSurfaceIntersection> SkewCylinders(const Surface& a, const Box2& boxA,
                                                 const Surface& b, const Box2& boxB,
                                                 const Context& ctx)
{
    const Cylinder* ca = a.Target<Cylinder>();
    const Cylinder* cb = b.Target<Cylinder>();
    if (ca == nullptr || cb == nullptr) {
        return Fail(ErrorCode::InvalidInput, "SkewCylinders: dispatched for other surface types");
    }
    const Pair p = PairOf(a, b, *ca, *cb);
    SurfaceSurfaceIntersection result;
    switch (ConfigurationOf(p)) {
    case Configuration::Apart:
        return result; // Proved separate.
    case Configuration::TouchingOutside:
        result.Points.push_back(ContactPoint(p));
        return result;
    case Configuration::TouchingInside:
        return Fail(ErrorCode::NotSupported,
                    "SkewCylinders: the axes' distance equals the radii's difference, so the "
                    "cylinders are tangent along their intersection; the singular curve is not "
                    "supported yet");
    case Configuration::Crossing:
        break;
    }
    const TrigPoly d = ExactDiscriminant(p);
    const ExactPolynomial quartic = Quartic(d);
    STEREON_TRY(const std::vector<double> branchPoints, BranchPoints(d, quartic));
    for (const Loop& loop : LoopsOf(d, quartic, branchPoints)) {
        STEREON_TRY(std::optional<IntersectionBranch> branch, BranchOf(p, loop, boxA, boxB, ctx));
        if (branch) {
            result.Branches.push_back(std::move(*branch));
        }
    }
    return result;
}

} // namespace Stereon::Detail
