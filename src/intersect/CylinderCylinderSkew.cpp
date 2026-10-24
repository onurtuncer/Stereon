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
#include <expected>
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
#include "Exact.hpp"

/// Cylinder–cylinder with skew or intersecting axes (ADR-0010's decisions on
/// acceptance, RS-0001 Tier B, the Boolean prototype's M3b). Along A's
/// ruling at angle u, B's equation is a quadratic in A's height v, so the
/// curve is two sheets v±(u) over the angles where its discriminant is not
/// negative. With σ the angle from the axes' common normal and δ the
/// distance between the axes, the discriminant is r_B² - (r_A cos σ - δ)²:
/// the sheets exist where |r_A cos σ - δ| <= r_B, and the configuration
/// follows from δ against r_A + r_B and |r_A - r_B|, decided exactly on the
/// input doubles.
///
/// - δ > r_A + r_B: proved separate; δ = r_A + r_B: one point of contact.
/// - |r_A - r_B| < δ < r_A + r_B: one closed loop over an arc of σ about
///   the common normal.
/// - δ < |r_A - r_B|: two loops, over two arcs of σ when A is the larger
///   cylinder (B passes through A) and over the whole turn, one per sheet,
///   when A is the smaller (A passes through B).
/// - δ = |r_A - r_B|: internally tangent at a node, a double root of the
///   discriminant, where the loops meet; δ = 0 with equal radii: two
///   ellipses through two nodes.
///
/// A branch is parametrised so that it is smooth: by A's angle where the
/// discriminant stays positive, and otherwise so that the signed square
/// root of the discriminant is smooth across the arc's ends, the folds
/// where the curve's tangent runs along A's ruling. Branches end at nodes
/// and at the chart seams of both cylinders. A branch's 3D curve and both
/// pcurves are cubic B-splines through samples of the closed form, refined
/// until CertifyDeviation proves them within the Context's resolution of
/// both surfaces (decision D3). Where the requested boxes are bounded, a
/// branch is first clipped to them by interval evaluation of the closed
/// form, so that a loop reaching far along nearly parallel axes costs only
/// what lies inside.

namespace Stereon::Detail {

namespace {

constexpr double Pi = std::numbers::pi;

/// Samples per branch: doubled from the first count until the fit
/// certifies. The budget bounds the work, not the geometry (ADR-0010).
constexpr std::size_t InitialSamples = 32;
constexpr std::size_t SampleBudget = std::size_t{1} << 14;
/// How often a track's range is bisected at most when clipping it to the
/// requested boxes, and how many enclosures the clip evaluates at most:
/// the bisection follows the boxes' boundaries, so the depth costs a few
/// evaluations per boundary, and the budget keeps a loose enclosure from
/// turning it into a sweep. Then how many samples scan a piece for B's
/// seam.
constexpr int ClipDepth = 40;
constexpr std::size_t ClipBudget = std::size_t{1} << 16;
constexpr std::size_t SeamScan = 256;
/// A seam that comes within this of a fold's extreme angle only touches
/// the loop there: one cut, not the two that acos would amplify rounding
/// into.
constexpr double FoldTouch = 1e-9;

[[nodiscard]] ExactVector<4> ExactAxisOf(const Cylinder& c)
{
    const ExactVector<4> cross = ExactCross(Exactly(c.XAxis), Exactly(c.YAxis));
    return c.RightHanded ? cross : ExactNegated(cross);
}

// --- Exact decisions ----------------------------------------------------

/// The axes' common normal n = a x b and s = (O_B - O_A) . n: the distance
/// between the axes is |s| / |n|, and the sign of s says on which side of
/// A's axis along n B's axis lies.
struct ExactNormal {
    ExactVector<64> Normal;
    ExactScalar<768> Offset;
};

[[nodiscard]] ExactNormal ExactNormalOf(const Cylinder& a, const Cylinder& b)
{
    const ExactVector<64> normal = ExactCross(ExactAxisOf(a), ExactAxisOf(b));
    return {normal, ExactDot(ExactOffset(b.Origin, a.Origin), normal)};
}

/// The sign of s² - k² |n|²: the distance between the axes against k, a
/// polynomial of degree ten in the input doubles. The interval decides it
/// unless the pair is within rounding of the threshold; the exact stage
/// then runs on the heap, where its capacity as an Expansion would exceed
/// a million components.
[[nodiscard]] Sign CompareDistance(const ExactNormal& e, const ExactScalar<2>& k)
{
    const Interval normSquared =
        Squared(e.Normal.Hull[0]) + Squared(e.Normal.Hull[1]) + Squared(e.Normal.Hull[2]);
    const Interval estimate = Squared(e.Offset.Hull) - Squared(k.Hull) * normSquared;
    return FilteredSign(estimate, [&] {
        DynamicExpansion norm;
        for (const auto& component : e.Normal.Value) {
            const DynamicExpansion c{component};
            norm = Sum(norm, Product(c, c));
        }
        const DynamicExpansion s{e.Offset.Value};
        const DynamicExpansion factor{k.Value};
        return Difference(Product(s, s), Product(Product(factor, factor), norm)).GetSign();
    });
}

// --- The nominal pair ---------------------------------------------------

/// What the closed form needs of the pair, in doubles: the unit axes, the
/// unit common normal from A's axis towards B's (either way when they
/// meet), its angle on A, the distance between the axes, the cosine of the
/// angle between the axes with its sine squared as |â x b̂|², and that sine
/// times r_A, the scale of the discriminant's root.
struct Pair {
    const Surface& ErasedA;
    const Surface& ErasedB;
    const Cylinder& A;
    const Cylinder& B;
    Vec3 AxisA;
    Vec3 AxisB;
    Vec3 Normal;
    Vec3 Offset; ///< O_A - O_B
    double Phi;
    double Distance;
    double Cosine;
    double SineSquared;
    double RootScale;
};

[[nodiscard]] Pair PairOf(const Surface& erasedA, const Surface& erasedB, const Cylinder& a,
                          const Cylinder& b, const Sign side)
{
    const Vec3 axisA = Normalized(Axis(a));
    const Vec3 axisB = Normalized(Axis(b));
    const Vec3 cross = Cross(axisA, axisB);
    const Vec3 normal = side == Sign::Negative ? -Normalized(cross) : Normalized(cross);
    const double distance = side == Sign::Zero ? 0.0 : std::abs(Dot(b.Origin - a.Origin, normal));
    const double sineSquared = LengthSquared(cross);
    return {erasedA,
            erasedB,
            a,
            b,
            axisA,
            axisB,
            normal,
            a.Origin - b.Origin,
            WrappedAngle(Dot(normal, a.YAxis), Dot(normal, a.XAxis)),
            distance,
            Dot(axisA, axisB),
            sineSquared,
            std::sqrt(sineSquared) * a.Radius};
}

template <Scalar T>
[[nodiscard]] BasicVec3<T> Wide(const Vec3& v)
{
    return {T{v.X}, T{v.Y}, T{v.Z}};
}

template <Scalar T>
[[nodiscard]] BasicPoint3<T> Wide(const Point3& p)
{
    return {T{p.X}, T{p.Y}, T{p.Z}};
}

template <Scalar T>
[[nodiscard]] BasicCylinder<T> Wide(const Cylinder& c)
{
    return {Wide<T>(c.Origin), T{c.Radius}, Wide<T>(c.XAxis), Wide<T>(c.YAxis), c.RightHanded};
}

template <Scalar T>
[[nodiscard]] T SqrtOf(const T& x)
{
    using std::sqrt;
    return sqrt(x);
}

// --- The configuration --------------------------------------------------

/// How the cylinders meet, with the cosines of σ where the discriminant
/// vanishes, (δ ∓ r_B) / r_A: exactly ∓1 where the exact decision puts a
/// node, and kept off ±1 where it puts an arc's end strictly inside the
/// turn or the discriminant strictly positive, so that the square roots
/// below never see a negative.
enum class Meeting : std::uint8_t {
    Separate,
    Point,
    Loop,
    TwoLoops,
    TwoTurns,
    Lobes,
    TouchingTurns,
    Ellipses,
};

struct Configuration {
    Meeting Kind;
    double CosLow;
    double CosHigh;
};

[[nodiscard]] Meeting MeetingOf(const Pair& p, const Sign outer, const Sign inner)
{
    if (outer == Sign::Positive) {
        return Meeting::Separate;
    }
    if (outer == Sign::Zero) {
        return Meeting::Point;
    }
    if (inner == Sign::Positive) {
        return Meeting::Loop;
    }
    const bool larger = p.A.Radius > p.B.Radius;
    if (inner == Sign::Negative) {
        return larger ? Meeting::TwoLoops : Meeting::TwoTurns;
    }
    if (p.A.Radius == p.B.Radius) {
        return Meeting::Ellipses;
    }
    return larger ? Meeting::Lobes : Meeting::TouchingTurns;
}

/// x inside the open interval (-1, 1), above 1, or below -1, by the
/// smallest step: the decision's strict inequality, imposed on the rounded
/// quotient.
[[nodiscard]] double Within(const double x)
{
    return std::clamp(x, std::nextafter(-1.0, 0.0), std::nextafter(1.0, 0.0));
}
[[nodiscard]] double Above(const double x)
{
    return std::max(x, std::nextafter(1.0, 2.0));
}
[[nodiscard]] double Below(const double x)
{
    return std::min(x, std::nextafter(-1.0, -2.0));
}

[[nodiscard]] Configuration Configure(const Pair& p, const Sign outer, const Sign inner)
{
    const Meeting kind = MeetingOf(p, outer, inner);
    const double low = (p.Distance - p.B.Radius) / p.A.Radius;
    const double high = (p.Distance + p.B.Radius) / p.A.Radius;
    switch (kind) {
    case Meeting::Loop:
        return {kind, Within(low), Above(high)};
    case Meeting::TwoLoops:
        return {kind, Within(low), Within(high)};
    case Meeting::TwoTurns:
        return {kind, Below(low), Above(high)};
    case Meeting::Lobes:
        return {kind, Within(low), 1.0};
    case Meeting::TouchingTurns:
        return {kind, -1.0, Above(high)};
    case Meeting::Ellipses:
        return {kind, -1.0, 1.0};
    default:
        return {kind, low, high};
    }
}

void LogMeeting(const Configuration& c, const Context& ctx)
{
    switch (c.Kind) {
    case Meeting::Point:
        ctx.Log(LogLevel::Info, "Intersect: cylinders with skew axes touch at one point: the axes "
                                "are the sum of the radii apart (ADR-0029)");
        break;
    case Meeting::Lobes:
    case Meeting::TouchingTurns:
        ctx.Log(LogLevel::Info, "Intersect: cylinders internally tangent: the axes are the "
                                "difference of the radii apart, and the branches meet at the "
                                "point of tangency (ADR-0029)");
        break;
    case Meeting::Ellipses:
        ctx.Log(LogLevel::Info, "Intersect: cylinders of one radius with meeting axes: two "
                                "ellipses, as branches between their two crossings (ADR-0029)");
        break;
    default:
        break;
    }
}

// --- The closed form ----------------------------------------------------

/// B's equation along A's ruling at angle u, α v² + 2 β v + γ = 0 in A's
/// height v, with α the sine squared of the axes' angle: Q is the ruling's
/// point at height 0 relative to B's origin, β = Q.â - cos Q.b̂, and
/// γ = |Q x b̂|² - r_B², the squared distance of that point from B's axis
/// against B's radius.
template <Scalar T>
struct Quadratic {
    T Beta;
    T Gamma;
};

template <Scalar T>
[[nodiscard]] Quadratic<T> QuadraticAt(const Pair& p, const T& u)
{
    const BasicVec3<T> q =
        Wide<T>(p.Offset) +
        T{p.A.Radius} * RadialDirection(Wide<T>(p.A.XAxis), Wide<T>(p.A.YAxis), u, 0);
    const BasicVec3<T> axisB = Wide<T>(p.AxisB);
    const BasicVec3<T> across = Cross(q, axisB);
    const T radius{p.B.Radius};
    return {Dot(q, Wide<T>(p.AxisA)) - T{p.Cosine} * Dot(q, axisB),
            Dot(across, across) - radius * radius};
}

/// The root (-β + R) / α of the quadratic, R the signed square root of its
/// discriminant, by the form that does not cancel: γ / (-β - R) when β and
/// R have the same sign. Nearly parallel axes make α tiny and one root
/// huge; this keeps the other accurate.
[[nodiscard]] double HeightFrom(const Quadratic<double>& q, const double root, const double alpha)
{
    if (q.Beta * root > 0.0) {
        return q.Gamma / (-q.Beta - root);
    }
    return (-q.Beta + root) / alpha;
}

/// Both forms enclose the root; their meet is the tighter where α is small.
[[nodiscard]] IntervalScalar HeightFrom(const Quadratic<IntervalScalar>& q,
                                        const IntervalScalar& root, const IntervalScalar& alpha)
{
    Interval height = ((-q.Beta + root) / alpha).Value;
    const Interval denominator = (-q.Beta - root).Value;
    if (!Contains(denominator, 0.0)) {
        const Interval met = Intersection(height, (q.Gamma / IntervalScalar{denominator}).Value);
        if (!met.IsEmpty()) {
            height = met;
        }
    }
    return IntervalScalar{height};
}

/// The factors of the discriminant over r_A², (CosHigh - cos σ)(cos σ -
/// CosLow), in forms exact at a node: 2 sin²(σ/2) for CosHigh = 1 and
/// 2 cos²(σ/2) for CosLow = -1.
template <Scalar T>
[[nodiscard]] T FactorHigh(const Configuration& c, const T& sigma)
{
    if (c.CosHigh == 1.0) {
        const T s = Sin(sigma * T{0.5});
        return T{2.0} * s * s;
    }
    return T{c.CosHigh} - Cos(sigma);
}

template <Scalar T>
[[nodiscard]] T FactorLow(const Configuration& c, const T& sigma)
{
    if (c.CosLow == -1.0) {
        const T k = Cos(sigma * T{0.5});
        return T{2.0} * k * k;
    }
    return Cos(sigma) - T{c.CosLow};
}

/// The sign of x as a scalar, 0 for zero; for an interval, the hull of the
/// signs it contains.
[[nodiscard]] double SignScalar(const double x)
{
    if (x > 0.0) {
        return 1.0;
    }
    return x < 0.0 ? -1.0 : 0.0;
}
[[nodiscard]] IntervalScalar SignScalar(const IntervalScalar& x)
{
    const double lower = x.Value.GetLower() > 0.0 ? 1.0 : -1.0;
    const double upper = x.Value.GetUpper() < 0.0 ? -1.0 : 1.0;
    return IntervalScalar{Interval{lower, upper}};
}

/// A branch candidate: how its parameter t runs over A's angle, and the
/// sign of the discriminant's root along it.
struct Track {
    enum class Shape : std::uint8_t {
        /// t is A's angle u; the root keeps the sign Sheet.
        Turn,
        /// σ = Mirror (Center - Half cos t) over an arc from fold to fold,
        /// closed over a period of t; the root has the sign of sin t.
        /// Center is 0 for the arc about the common normal.
        FoldLoop,
        /// σ = Mirror Fold (1 - t²), t in [-1, 1]: from the node through the
        /// fold at t = 0 and back; the root has the sign of t.
        FoldLobe,
    };
    Shape Kind;
    Interval Range;
    bool Closed = false;
    double Sheet = 1.0;
    double Mirror = 1.0;
    double Center = 0.0;
    double Half = 0.0;
    double Fold = 0.0;
    /// acos of the configuration's cosines: the arc's ends in σ.
    double AngleLow = 0.0;
    double AngleHigh = 0.0;
};

template <Scalar T>
[[nodiscard]] T SigmaOf(const Pair& p, const Track& k, const T& t)
{
    switch (k.Kind) {
    case Track::Shape::Turn:
        return t - T{p.Phi};
    case Track::Shape::FoldLoop:
        return T{k.Mirror} * (T{k.Center} - T{k.Half} * Cos(t));
    case Track::Shape::FoldLobe:
        return T{k.Mirror * k.Fold} * (T{1.0} - t * t);
    }
    return T{};
}

/// The signed square root of the discriminant over r_A² along the track.
/// A factor vanishing at an arc's end |σ| = ψ is 2 sin((ψ + |σ|)/2)
/// sin((ψ - |σ|)/2), and over a fold-to-fold arc (ψ - |σ|)/2 is Half
/// sin²(t/2) at one end and Half cos²(t/2) at the other; the sines that
/// vanish are taken out with the sign, so that the root is smooth in t.
template <Scalar T>
[[nodiscard]] T RootOf(const Configuration& c, const Track& k, const T& t, const T& sigma)
{
    switch (k.Kind) {
    case Track::Shape::Turn:
        return T{k.Sheet} * SqrtOf(FactorHigh(c, sigma) * FactorLow(c, sigma));
    case Track::Shape::FoldLoop: {
        const T s = Sin(t * T{0.5});
        const T co = Cos(t * T{0.5});
        const T ends = Sin(T{k.Half} * s * s) * Sin(T{k.Half} * co * co);
        const T magnitude = T{k.Mirror} * sigma;
        const T rest = k.Center == 0.0 ? T{2.0} * (T{c.CosHigh} - Cos(sigma))
                                       : T{4.0} * Sin((magnitude + T{k.AngleHigh}) * T{0.5}) *
                                             Sin((magnitude + T{k.AngleLow}) * T{0.5});
        return SignScalar(Sin(t)) * SqrtOf(rest * ends);
    }
    case Track::Shape::FoldLobe: {
        const T magnitude = T{k.Mirror} * sigma;
        const T rest = Sin((magnitude + T{k.Fold}) * T{0.5}) * Sin(T{0.5 * k.Fold} * t * t);
        return SignScalar(t) * T{2.0} * Sin(magnitude * T{0.5}) * SqrtOf(rest);
    }
    }
    return T{};
}

/// A point of the track: its parameters on A and the point in space.
template <Scalar T>
struct Place {
    T U;
    T Height;
    BasicPoint3<T> Point;
};

template <Scalar T>
[[nodiscard]] Place<T> PlaceOf(const Pair& p, const Configuration& c, const Track& k, const T& t)
{
    const T sigma = SigmaOf(p, k, t);
    const T u = k.Kind == Track::Shape::Turn ? t : T{p.Phi} + sigma;
    const T root = T{p.RootScale} * RootOf(c, k, t, sigma);
    const T height = HeightFrom(QuadraticAt(p, u), root, T{p.SineSquared});
    return {u, height, Eval(Wide<T>(p.A), BasicPoint2<T>{u, height})};
}

/// B's parameters of a point.
[[nodiscard]] Point2 OnB(const Pair& p, const Point3& point)
{
    const Vec3 offset = point - p.B.Origin;
    return {WrappedAngle(Dot(offset, p.B.YAxis), Dot(offset, p.B.XAxis)), Dot(offset, p.AxisB)};
}

// --- The tracks of a configuration --------------------------------------

[[nodiscard]] std::vector<Track> TracksOf(const Configuration& c, const Pair& p)
{
    using Shape = Track::Shape;
    const double angleLow = std::acos(c.CosLow);
    const double angleHigh = std::acos(std::clamp(c.CosHigh, -1.0, 1.0));
    const Interval turn{0.0, CirclePeriod};
    switch (c.Kind) {
    case Meeting::Loop:
        return {{.Kind = Shape::FoldLoop, .Range = turn, .Closed = true, .Half = angleLow}};
    case Meeting::TwoLoops: {
        const Track right{.Kind = Shape::FoldLoop,
                          .Range = turn,
                          .Closed = true,
                          .Center = (angleHigh + angleLow) / 2.0,
                          .Half = (angleLow - angleHigh) / 2.0,
                          .AngleLow = angleLow,
                          .AngleHigh = angleHigh};
        Track left = right;
        left.Mirror = -1.0;
        return {right, left};
    }
    case Meeting::Lobes: {
        const Track right{.Kind = Shape::FoldLobe, .Range = Interval{-1.0, 1.0}, .Fold = angleLow};
        Track left = right;
        left.Mirror = -1.0;
        return {right, left};
    }
    case Meeting::TwoTurns:
        return {{.Kind = Shape::Turn, .Range = turn, .Closed = true, .Sheet = 1.0},
                {.Kind = Shape::Turn, .Range = turn, .Closed = true, .Sheet = -1.0}};
    case Meeting::TouchingTurns: {
        // The node is at σ = pi: both sheets run from it round to it.
        const double node = WrappedAngle(Sin(p.Phi + Pi), Cos(p.Phi + Pi));
        const Interval range{node, node + CirclePeriod};
        return {{.Kind = Shape::Turn, .Range = range, .Sheet = 1.0},
                {.Kind = Shape::Turn, .Range = range, .Sheet = -1.0}};
    }
    case Meeting::Ellipses: {
        // Nodes at σ = 0 and pi: each sheet is an ellipse's half between them.
        const Interval first{p.Phi, p.Phi + Pi};
        const Interval second{p.Phi + Pi, p.Phi + CirclePeriod};
        return {{.Kind = Shape::Turn, .Range = first, .Sheet = 1.0},
                {.Kind = Shape::Turn, .Range = second, .Sheet = 1.0},
                {.Kind = Shape::Turn, .Range = first, .Sheet = -1.0},
                {.Kind = Shape::Turn, .Range = second, .Sheet = -1.0}};
    }
    default:
        return {};
    }
}

// --- Clipping to the requested boxes ------------------------------------

/// Whether a box restricts a cylinder's domain at all: a turn or more in u
/// and everything in v restricts nothing.
[[nodiscard]] bool Restricts(const Box2& box)
{
    const bool wholeTurn = !box.X.IsBounded() || Width(box.X) >= CirclePeriod;
    return !wholeTurn || box.Y != Interval::Everything();
}

/// The shifts by whole periods that could bring `angle` onto `range`, both
/// narrower than a period.
template <typename Test>
[[nodiscard]] bool SomeShift(const Interval& angle, const Interval& range, Test&& test)
{
    const double first = std::ceil((range.GetLower() - angle.GetUpper()) / CirclePeriod);
    const double last = std::floor((range.GetUpper() - angle.GetLower()) / CirclePeriod);
    // Two intervals narrower than a turn need at most three candidate shifts.
    // Count iterations with an integer even at the limit of double precision.
    for (int step = 0; step < 3; ++step) {
        const double k = first + static_cast<double>(step);
        if (k > last) {
            break;
        }
        const Interval shifted =
            Interval{angle.GetLower(), angle.GetUpper()} + Interval{k} * Interval{CirclePeriod};
        if (test(shifted)) {
            return true;
        }
    }
    return false;
}

/// Whether some angle of `angle` lies in `range` modulo the period.
[[nodiscard]] bool MeetsAngles(const Interval& angle, const Interval& range)
{
    if (!range.IsBounded() || Width(range) >= CirclePeriod || !angle.IsBounded() ||
        Width(angle) >= CirclePeriod) {
        return true;
    }
    return SomeShift(angle, range, [&](const Interval& s) { return Intersects(s, range); });
}

/// Whether every angle of `angle` lies in `range` modulo the period.
[[nodiscard]] bool AnglesWithin(const Interval& angle, const Interval& range)
{
    if (!range.IsBounded() || Width(range) >= CirclePeriod) {
        return true;
    }
    if (!angle.IsBounded() || Width(angle) >= CirclePeriod) {
        return false;
    }
    return SomeShift(angle, range, [&](const Interval& s) { return Contains(range, s); });
}

enum class Fit : std::uint8_t { Outside, Partly, Inside };

/// Refuse narrow charts whose shift can no longer be counted exactly. A full
/// turn needs no reduction. Never silently discard angular constraints at
/// magnitudes where the double chart cannot resolve individual turns.
[[nodiscard]] Result<void> CheckAngularChart(const Interval& range)
{
    constexpr double MaxChartAngle = 0x1p52;
    if (range.IsBounded() && Width(range) < CirclePeriod &&
        std::max(std::abs(range.GetLower()), std::abs(range.GetUpper())) >= MaxChartAngle) {
        return Fail(ErrorCode::InvalidInput,
                    "CylinderCylinder: angular chart is too large to resolve whole-period shifts");
    }
    return {};
}

/// Combine two oriented half-plane tests for a narrow sector (intersection)
/// or a wide one (union), retaining uncertainty for subdivision.
[[nodiscard]] Fit HalfPlaneFit(const Interval& left, const Interval& right, const bool wide)
{
    const bool leftOut = left.GetUpper() < 0.0;
    const bool rightOut = right.GetUpper() < 0.0;
    const bool leftIn = left.GetLower() >= 0.0;
    const bool rightIn = right.GetLower() >= 0.0;
    if (wide) {
        return leftOut && rightOut ? Fit::Outside : (leftIn || rightIn ? Fit::Inside : Fit::Partly);
    }
    return leftOut || rightOut ? Fit::Outside : (leftIn && rightIn ? Fit::Inside : Fit::Partly);
}

/// Membership of B's angular sector without atan2 or a branch cut. For an
/// arc at most pi wide both oriented boundary half-planes must hold; for
/// a wider arc either suffices. Interval signs prove inclusion or exclusion,
/// otherwise the clip subdivides. A zero-width arc must also point forwards.
[[nodiscard]] Fit SectorFit(const IntervalScalar& x, const IntervalScalar& y, const Interval& range)
{
    if (!range.IsBounded() || Width(range) >= CirclePeriod) {
        return Fit::Inside;
    }
    using I = IntervalScalar;
    const I start{range.GetLower()};
    const I end{range.GetUpper()};
    const Interval left = (Cos(start) * y - Sin(start) * x).Value;
    const Interval right = (Sin(end) * x - Cos(end) * y).Value;
    const Fit sides = HalfPlaneFit(left, right, Width(range) > Pi);
    if (sides == Fit::Outside) {
        return sides;
    }
    if (Width(range) == 0.0) {
        const Interval forward = (Cos(start) * x + Sin(start) * y).Value;
        return forward.GetUpper() < 0.0 ? Fit::Outside : Fit::Partly;
    }
    return sides;
}

/// How the track over the span lies against both requested boxes. Outside
/// is proved; Inside says every enclosed point satisfies every restriction.
[[nodiscard]] Fit FitOf(const Pair& p, const Configuration& c, const Track& k, const Interval& span,
                        const Box2& boxA, const Box2& boxB)
{
    using I = IntervalScalar;
    const Place<I> place = PlaceOf(p, c, k, I{span});
    const BasicVec3<I> offsetB = place.Point - Wide<I>(p.B.Origin);
    const Interval heightB = Dot(offsetB, Wide<I>(p.AxisB)).Value;
    const Fit angleB =
        SectorFit(Dot(offsetB, Wide<I>(p.B.XAxis)), Dot(offsetB, Wide<I>(p.B.YAxis)), boxB.X);
    if (!MeetsAngles(place.U.Value, boxA.X) || !Intersects(place.Height.Value, boxA.Y) ||
        !Intersects(heightB, boxB.Y) || angleB == Fit::Outside) {
        return Fit::Outside;
    }
    const bool inside = AnglesWithin(place.U.Value, boxA.X) &&
                        Contains(boxA.Y, place.Height.Value) && Contains(boxB.Y, heightB) &&
                        angleB == Fit::Inside;
    return inside ? Fit::Inside : Fit::Partly;
}

/// The clip in progress: what it keeps, and how many enclosures it has
/// evaluated against its budget.
struct Clip {
    std::vector<Interval> Kept;
    std::size_t Evaluations = 0;
};

/// Keeps the spans of t over which the track can lie inside the boxes, by
/// bisection to ClipDepth within the budget; adjacent kept spans merge.
void ClipSpan(const Pair& p, const Configuration& c, const Track& k, const Interval& span,
              const Box2& boxA, const Box2& boxB, const int depth, Clip& clip)
{
    ++clip.Evaluations;
    const Fit fit = FitOf(p, c, k, span, boxA, boxB);
    if (fit == Fit::Outside) {
        return;
    }
    const double m = Midpoint(span);
    if (fit == Fit::Partly && depth < ClipDepth && clip.Evaluations < ClipBudget &&
        span.GetLower() < m && m < span.GetUpper()) {
        ClipSpan(p, c, k, Interval{span.GetLower(), m}, boxA, boxB, depth + 1, clip);
        ClipSpan(p, c, k, Interval{m, span.GetUpper()}, boxA, boxB, depth + 1, clip);
        return;
    }
    if (!clip.Kept.empty() && clip.Kept.back().GetUpper() == span.GetLower()) {
        clip.Kept.back() = Interval{clip.Kept.back().GetLower(), span.GetUpper()};
    } else {
        clip.Kept.push_back(span);
    }
}

/// The spans of the track's range to fit: the whole range when the boxes
/// restrict nothing, otherwise what the clip keeps, with a closed track's
/// first and last spans joined across the end of its period.
[[nodiscard]] std::vector<Interval> ClipTrack(const Pair& p, const Configuration& c, const Track& k,
                                              const Box2& boxA, const Box2& boxB)
{
    if (!Restricts(boxA) && !Restricts(boxB)) {
        return {k.Range};
    }
    Clip clip;
    ClipSpan(p, c, k, k.Range, boxA, boxB, 0, clip);
    std::vector<Interval>& kept = clip.Kept;
    const bool wraps = k.Closed && kept.size() >= 2 &&
                       kept.front().GetLower() == k.Range.GetLower() &&
                       kept.back().GetUpper() == k.Range.GetUpper();
    if (wraps) {
        const Interval joined{kept.back().GetLower(), kept.front().GetUpper() + Width(k.Range)};
        kept.back() = joined;
        kept.erase(kept.begin());
    }
    return kept;
}

// --- Cuts at the seams --------------------------------------------------

/// The parameters of a fold-to-fold loop at which σ takes the value
/// `target`, over a period of t on either side of the range: where
/// cos t = x, or at the fold alone where the loop only touches the value.
[[nodiscard]] std::vector<double> LoopParametersAt(const Track& k, const double target)
{
    const double x = (k.Center - k.Mirror * target) / k.Half;
    if (!(std::abs(x) <= 1.0)) {
        return {};
    }
    const bool touches = std::abs(x) >= 1.0 - FoldTouch;
    const double t0 = touches ? (x > 0.0 ? 0.0 : Pi) : std::acos(x);
    std::vector<double> roots;
    for (const double shift : {-CirclePeriod, 0.0, CirclePeriod}) {
        roots.push_back(t0 + shift);
        if (!touches) {
            roots.push_back(CirclePeriod - t0 + shift);
        }
    }
    return roots;
}

/// The parameters of a lobe at which σ takes the value `target`: where
/// 1 - t² = target / (Mirror Fold), or the fold alone where it is touched.
[[nodiscard]] std::vector<double> LobeParametersAt(const Track& k, const double target)
{
    const double squared = 1.0 - target / (k.Mirror * k.Fold);
    if (!(squared >= 0.0 && squared <= 1.0)) {
        return {};
    }
    if (squared <= FoldTouch) {
        return {0.0};
    }
    return {-std::sqrt(squared), std::sqrt(squared)};
}

/// The parameters at which σ(t) takes the value `target`; a turn's are its
/// caller's, since there t is the angle itself.
[[nodiscard]] std::vector<double> ParametersAt(const Track& k, const double target)
{
    switch (k.Kind) {
    case Track::Shape::FoldLoop:
        return LoopParametersAt(k, target);
    case Track::Shape::FoldLobe:
        return LobeParametersAt(k, target);
    default:
        return {};
    }
}

/// Adds a cut strictly inside the span.
void AddCut(const Interval& span, const double t, std::vector<double>& cuts)
{
    if (span.GetLower() < t && t < span.GetUpper()) {
        cuts.push_back(t);
    }
}

/// Where A's angle u = Phi + σ crosses a multiple of the period inside the
/// span: for a turn, at the multiples themselves.
void SeamCutsOnA(const Pair& p, const Track& k, const Interval& span, std::vector<double>& cuts)
{
    const double first = std::ceil(span.GetLower() / CirclePeriod - 1.0);
    const double last = std::floor(span.GetUpper() / CirclePeriod + 1.0);
    for (double n = first; n <= last; n += 1.0) {
        if (k.Kind == Track::Shape::Turn) {
            AddCut(span, n * CirclePeriod, cuts);
            continue;
        }
        for (const double t : ParametersAt(k, n * CirclePeriod - p.Phi)) {
            AddCut(span, t, cuts);
        }
    }
}

/// (P(t) - O_B) . Y_B: zero with (P(t) - O_B) . X_B positive where B's
/// angle crosses its seam.
[[nodiscard]] double SeamSideOnB(const Pair& p, const Configuration& c, const Track& k,
                                 const double t)
{
    return Dot(PlaceOf(p, c, k, t).Point - p.B.Origin, p.B.YAxis);
}

/// The crossing between lo and hi, where the seam side changes sign, to
/// rounding.
[[nodiscard]] double BisectSeam(const Pair& p, const Configuration& c, const Track& k, double lo,
                                double hi)
{
    const bool lowerNegative = SeamSideOnB(p, c, k, lo) < 0.0;
    for (int i = 0; i < 64 && lo < hi; ++i) {
        const double m = 0.5 * (lo + hi);
        if (m <= lo || m >= hi) {
            break;
        }
        ((SeamSideOnB(p, c, k, m) < 0.0) == lowerNegative ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
}

/// Where B's angle crosses its seam inside the span: scanned at SeamScan
/// samples, a jump of more than half a turn between neighbours marks a
/// crossing, which is then bisected.
void SeamCutsOnB(const Pair& p, const Configuration& c, const Track& k, const Interval& span,
                 std::vector<double>& cuts)
{
    const double width = Width(span);
    const auto at = [&](const std::size_t i) {
        return span.GetLower() + width * static_cast<double>(i) / static_cast<double>(SeamScan);
    };
    double previous = OnB(p, PlaceOf(p, c, k, span.GetLower()).Point).X;
    for (std::size_t i = 1; i <= SeamScan; ++i) {
        const double t = at(i);
        const double angle = OnB(p, PlaceOf(p, c, k, t).Point).X;
        if (std::abs(angle - previous) > Pi) {
            const double before = at(i - 1);
            AddCut(span, BisectSeam(p, c, k, before, t), cuts);
        }
        previous = angle;
    }
}

/// The span's ends and the cuts inside it, in order; a cut within rounding
/// of the one before joins it, so that no sliver is fitted on its own.
[[nodiscard]] std::vector<double> CutsOf(const Pair& p, const Configuration& c, const Track& k,
                                         const Interval& span)
{
    std::vector<double> cuts{span.GetLower(), span.GetUpper()};
    SeamCutsOnA(p, k, span, cuts);
    SeamCutsOnB(p, c, k, span, cuts);
    std::ranges::sort(cuts);
    const double sliver = 1e-9 * Width(span); // STEREON-LINT-ALLOW(epsilon): a parameter gap
    std::vector<double> distinct{cuts.front()};
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        if (cuts[i] - distinct.back() > sliver) {
            distinct.push_back(cuts[i]);
        }
    }
    return distinct;
}

// --- Sampling, fitting and certifying -----------------------------------

struct Samples {
    std::vector<Point3> Points;
    std::vector<Point2> OnA;
    std::vector<Point2> OnB;
};

/// A's angles into [0, 2 pi]: the piece crosses no seam inside, so the
/// shift of a whole number of periods that brings its middle sample into
/// the chart brings every sample, and the ends, which may lie on the seam
/// to rounding, are clamped onto it.
void ChartOnA(std::vector<Point2>& onA)
{
    const double shift = CirclePeriod * std::floor(onA[onA.size() / 2].X / CirclePeriod);
    for (Point2& uv : onA) {
        uv.X = std::clamp(uv.X - shift, 0.0, CirclePeriod);
    }
}

/// B's angles, each wrapped into [0, 2 pi) on its own, into one chart: the
/// piece crosses no seam inside, so only an end sample at the seam can have
/// wrapped to the other side of it; it joins its neighbour's side.
void ChartOnB(std::vector<Point2>& onB)
{
    const auto join = [](Point2& end, const Point2& neighbour) {
        if (end.X - neighbour.X > Pi) {
            end.X -= CirclePeriod;
        } else if (neighbour.X - end.X > Pi) {
            end.X += CirclePeriod;
        }
    };
    join(onB.front(), onB[1]);
    join(onB.back(), onB[onB.size() - 2]);
    for (Point2& uv : onB) {
        uv.X = std::clamp(uv.X, 0.0, CirclePeriod);
    }
}

/// The closed form at count + 1 parameters evenly spaced over the range.
[[nodiscard]] Samples SampleTrack(const Pair& p, const Configuration& c, const Track& k,
                                  const Interval& range, const std::size_t count)
{
    Samples s;
    s.Points.reserve(count + 1);
    s.OnA.reserve(count + 1);
    s.OnB.reserve(count + 1);
    const double width = Width(range);
    for (std::size_t i = 0; i <= count; ++i) {
        const double t =
            range.GetLower() + width * static_cast<double>(i) / static_cast<double>(count);
        const Place<double> place = PlaceOf(p, c, k, t);
        s.Points.push_back(place.Point);
        s.OnA.push_back(Point2{place.U, place.Height});
        s.OnB.push_back(OnB(p, place.Point));
    }
    ChartOnA(s.OnA);
    ChartOnB(s.OnB);
    return s;
}

/// Cubic B-splines through the samples, on the shared parameter over the
/// range: the 3D curve and both pcurves, as a transverse branch.
[[nodiscard]] Result<IntersectionBranch> FitBranch(const Samples& s, const Interval& range)
{
    STEREON_TRY(const CurveFit fit, Interpolate(s.Points, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve curve, Reparametrized(fit.Curve, range));
    STEREON_TRY(const CurveFit2d fitA, Interpolate(s.OnA, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d onA, Reparametrized(fitA.Curve, range));
    STEREON_TRY(const CurveFit2d fitB, Interpolate(s.OnB, 3, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d onB, Reparametrized(fitB.Curve, range));
    return Crossing(Curve(std::move(curve)), Curve2d(std::move(onA)), Curve2d(std::move(onB)),
                    range);
}

/// The largest deviation at the midpoints between the samples: of the
/// fitted curve from the closed form, and of each surface through its
/// fitted pcurve from the fitted curve. A guide for the sample count only;
/// the certificate is the proof.
[[nodiscard]] double MidpointDeviation(const Pair& p, const Configuration& c, const Track& k,
                                       const IntersectionBranch& branch, const std::size_t count)
{
    const double width = Width(branch.Range);
    double worst = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        const double t = branch.Range.GetLower() +
                         width * (static_cast<double>(i) + 0.5) / static_cast<double>(count);
        const Point3 on = Eval(branch.Curve3d, t);
        worst = std::max({worst, Distance(on, PlaceOf(p, c, k, t).Point),
                          Distance(on, Eval(p.ErasedA, Eval(branch.OnA, t))),
                          Distance(on, Eval(p.ErasedB, Eval(branch.OnB, t)))});
    }
    return worst;
}

/// A piece that leaves the size box (ADR-0001) is outside the kernel's
/// domain: at that scale the rounding of a B-spline's poles alone exceeds
/// the resolution, so no fit could certify. Refused at once rather than
/// after the sample budget.
[[nodiscard]] Result<void> CheckSizeBox(const Samples& s)
{
    const Box3 sizeBox = Context::GetSizeBox();
    for (const Point3& point : s.Points) {
        if (!Contains(sizeBox, point)) {
            return Fail(ErrorCode::CertificationFailed,
                        "CylinderCylinder: a branch of the skew cylinders leaves the size box at " +
                            Number(point.X) + ", " + Number(point.Y) + ", " + Number(point.Z) +
                            " mm, where no deviation certifies within the resolution (ADR-0001)");
        }
    }
    return {};
}

/// The branch over one piece of a track: the samples are doubled until the
/// midpoints agree with the closed form and the certificate holds, within
/// the sample budget.
[[nodiscard]] Result<void> AddPiece(SurfaceSurfaceIntersection& result, const Pair& p,
                                    const Configuration& c, const Track& k, const Interval& range,
                                    const Context& ctx)
{
    const double target = ctx.GetResolution();
    for (std::size_t count = InitialSamples;; count *= 2) {
        const Samples samples = SampleTrack(p, c, k, range, count);
        STEREON_TRY_VOID(CheckSizeBox(samples));
        STEREON_TRY(IntersectionBranch branch, FitBranch(samples, range));
        if (MidpointDeviation(p, c, k, branch, count) <= target / 4.0) {
            const Result<void> certified = Certify(branch, p.ErasedA, p.ErasedB, ctx);
            if (certified) {
                result.Branches.push_back(std::move(branch));
                return {};
            }
            if (certified.error().GetCode() != ErrorCode::CertificationFailed) {
                return std::unexpected(certified.error());
            }
        }
        if (count >= SampleBudget) {
            return Fail(ErrorCode::CertificationFailed,
                        "CylinderCylinder: a branch of the skew cylinders exceeds the sample "
                        "budget before it certifies within the resolution " +
                            Number(target) + " mm");
        }
    }
}

/// Every piece of the track inside the boxes: the clipped spans, cut at
/// the seams of both cylinders.
[[nodiscard]] Result<void> AddTrack(SurfaceSurfaceIntersection& result, const Pair& p,
                                    const Configuration& c, const Track& k, const Box2& boxA,
                                    const Box2& boxB, const Context& ctx)
{
    for (const Interval& span : ClipTrack(p, c, k, boxA, boxB)) {
        const std::vector<double> cuts = CutsOf(p, c, k, span);
        for (std::size_t i = 1; i < cuts.size(); ++i) {
            STEREON_TRY_VOID(AddPiece(result, p, c, k, Interval{cuts[i - 1], cuts[i]}, ctx));
        }
    }
    return {};
}

// --- The point of external tangency -------------------------------------

/// A proved bound on the point's distance from the cylinder, by interval
/// evaluation of |(P - O) x axis| / |axis| - r with the axis the frame's
/// cross product.
[[nodiscard]] double PointBound(const Point3& point, const Cylinder& c)
{
    using I = IntervalScalar;
    const BasicVec3<I> axis = FrameAxis(Wide<I>(c.XAxis), Wide<I>(c.YAxis), c.RightHanded);
    const BasicVec3<I> cross = Cross(Wide<I>(point) - Wide<I>(c.Origin), axis);
    const Interval residual = (SqrtOf(Dot(cross, cross) / Dot(axis, axis)) - I{c.Radius}).Value;
    return std::max(std::abs(residual.GetLower()), std::abs(residual.GetUpper()));
}

/// The one point where externally tangent cylinders touch: on the common
/// perpendicular, r_A from A's axis along the normal, at the height of A's
/// foot, (w . â - cos w . b̂) / sin² with w = O_B - O_A.
[[nodiscard]] Result<void> AddPoint(SurfaceSurfaceIntersection& result, const Pair& p,
                                    const Box2& boxA, const Box2& boxB)
{
    const Vec3 w = -p.Offset;
    const double height = (Dot(w, p.AxisA) - p.Cosine * Dot(w, p.AxisB)) / p.SineSquared;
    const Point2 onA{p.Phi, height};
    const Point3 point = Eval(p.A, onA);
    const Point2 onB = OnB(p, point);
    if (!MeetsAngles(Interval{onA.X}, boxA.X) || !Contains(boxA.Y, onA.Y) ||
        !MeetsAngles(Interval{onB.X}, boxB.X) || !Contains(boxB.Y, onB.Y)) {
        return {};
    }
    result.Points.push_back(
        {point, onA, onB, std::max(PointBound(point, p.A), PointBound(point, p.B))});
    return {};
}

} // namespace

Result<SurfaceSurfaceIntersection> CylinderCylinderSkew(const Surface& a, const Box2& boxA,
                                                        const Surface& b, const Box2& boxB,
                                                        const double /*tolerance*/,
                                                        const Context& ctx)
{
    const Cylinder* ca = a.Target<Cylinder>();
    const Cylinder* cb = b.Target<Cylinder>();
    if (ca == nullptr || cb == nullptr) {
        return Fail(ErrorCode::InvalidInput,
                    "CylinderCylinderSkew: dispatched for other surface types");
    }
    STEREON_TRY_VOID(CheckAngularChart(boxA.X));
    STEREON_TRY_VOID(CheckAngularChart(boxB.X));
    const ExactNormal exact = ExactNormalOf(*ca, *cb);
    const Sign outer = CompareDistance(exact, ExactSum(ca->Radius, cb->Radius));
    const Sign inner = CompareDistance(exact, ExactDifference(ca->Radius, cb->Radius));
    const Pair p = PairOf(a, b, *ca, *cb, SignOf(exact.Offset));
    const Configuration c = Configure(p, outer, inner);
    LogMeeting(c, ctx);
    SurfaceSurfaceIntersection result;
    if (c.Kind == Meeting::Point) {
        STEREON_TRY_VOID(AddPoint(result, p, boxA, boxB));
        return result;
    }
    for (const Track& k : TracksOf(c, p)) {
        STEREON_TRY_VOID(AddTrack(result, p, c, k, boxA, boxB, ctx));
    }
    return result;
}

} // namespace Stereon::Detail
