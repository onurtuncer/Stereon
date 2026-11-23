// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Issue #195: the chart seams of two skew cylinders at, beside and between
// the folds of their intersection. The cylinders stay where they are and
// only their frames turn about their axes, so the curve is one fixed curve
// in space, known in closed form in the world's coordinates, and each seam
// is a line whose points on the other cylinder are counted independently
// of the handler: by the sign of a product that does not cancel, and by a
// dense sampling of the reference curve that brackets every change of the
// seam's side. Against that the answer is held to its whole contract:
//
// - every crossing of a seam is a cut, where exactly two branches end;
// - no branch crosses a seam inside by more than its bound;
// - the branches chain into the reference's closed curves with no end left
//   over and none shared by a third branch, and their length is the
//   reference's, so that no arc is published twice or left out;
// - both pcurves stay in their charts, without a jump;
// - a refusal is IntersectionFailed or CertificationFailed, and only where
//   two crossings lie within a few resolutions of each other: the pieces
//   between them could not be fitted, and nothing else excuses one.
//
// All of it in both orders of the operands and for both handednesses.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Intersect;

#include "Checks.hpp"

using IntersectTest::ctx;
using IntersectTest::DistanceFrom;
using Stereon::ContactType;
using Stereon::Cylinder;
using Stereon::ErrorCode;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Surface;
using Stereon::Vec3;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double radiusA = 3.0;

/// A is the cylinder of radius 3 about z through the origin; B runs along
/// x through (0, Offset, 0). A point of A at the world's angle u about z,
/// (3 cos u, y, z) with y = 3 sin u, is on B where
/// z² = (y - Low)(High - y), Low and High the nearest and farthest B gets
/// from the plane y = 0: a product that keeps its sign beside a fold,
/// where B's equation itself would cancel.
struct World {
    double Offset;
    double RadiusB;

    [[nodiscard]] double Low() const { return Offset - RadiusB; }
    [[nodiscard]] double High() const { return Offset + RadiusB; }
};

/// A piece of the reference curve: A's angle from From to To on the sheet
/// z >= 0 or z <= 0, on the side x >= 0 or mirrored to x <= 0.
struct Leg {
    double From;
    double To;
    double Sheet;
    double Side;
};

[[nodiscard]] Point3 OnLeg(const World& w, const Leg& leg, const double fraction)
{
    const double u = fraction == 1.0 ? leg.To : leg.From + fraction * (leg.To - leg.From);
    const double y = radiusA * std::sin(u);
    const double squared = std::max((y - w.Low()) * (w.High() - y), 0.0);
    return {leg.Side * radiusA * std::cos(u), y, leg.Sheet * std::sqrt(squared)};
}

/// The reference's closed curves, each a cycle of legs: one loop over the
/// top of A where B reaches past it, two loops where B passes through,
/// and two lobes from the node at the top where B touches A from inside.
[[nodiscard]] std::vector<std::vector<Leg>> LoopsOf(const World& w)
{
    const double low = std::asin(w.Low() / radiusA);
    if (w.High() >= radiusA) {
        const double top = pi / 2.0;
        if (w.High() > radiusA) {
            return {{{low, top, 1.0, 1.0},
                     {top, low, 1.0, -1.0},
                     {low, top, -1.0, -1.0},
                     {top, low, -1.0, 1.0}}};
        }
        return {{{top, low, 1.0, 1.0}, {low, top, -1.0, 1.0}},
                {{top, low, 1.0, -1.0}, {low, top, -1.0, -1.0}}};
    }
    const double high = std::asin(w.High() / radiusA);
    return {{{low, high, 1.0, 1.0}, {high, low, -1.0, 1.0}},
            {{low, high, 1.0, -1.0}, {high, low, -1.0, -1.0}}};
}

/// A dense sample of the reference: its points in order around each closed
/// curve, with the leg and the fraction of each, and its length. Each leg
/// is sampled at both of its ends, so that two samples next to each other
/// on different legs are one point, a fold or the top, and every stretch
/// between two other neighbours lies on one leg.
struct Sample {
    Point3 Point;
    std::size_t Leg;
    double Fraction;
};

struct Reference {
    World Pair;
    std::vector<std::vector<Leg>> Loops;
    std::vector<std::vector<Sample>> Samples;
    double Length = 0.0;
};

[[nodiscard]] Reference ReferenceOf(const World& w)
{
    constexpr std::size_t perLeg = std::size_t{1} << 15;
    Reference r{w, LoopsOf(w), {}, 0.0};
    for (const std::vector<Leg>& loop : r.Loops) {
        std::vector<Sample> samples;
        for (std::size_t leg = 0; leg < loop.size(); ++leg) {
            for (std::size_t i = 0; i <= perLeg; ++i) {
                const double fraction = static_cast<double>(i) / static_cast<double>(perLeg);
                samples.push_back({OnLeg(w, loop[leg], fraction), leg, fraction});
            }
        }
        for (std::size_t i = 0; i < samples.size(); ++i) {
            r.Length += Distance(samples[i].Point, samples[(i + 1) % samples.size()].Point);
        }
        r.Samples.push_back(std::move(samples));
    }
    return r;
}

/// Where a point stands against a cylinder's chart seam, the ruling at
/// angle zero: Ahead on the seam's side of the axis, and Side the signed
/// distance from the plane of the axis and the seam.
struct SeamSide {
    bool Ahead;
    double Side;
};

[[nodiscard]] SeamSide SeamSideOf(const Cylinder& c, const Point3& p)
{
    const Vec3 offset = p - c.Origin;
    return {Dot(offset, c.XAxis) > 0.0, Dot(offset, c.YAxis)};
}

/// The seam crossing between two samples of one leg, by bisection of the
/// seam's side on the reference itself.
[[nodiscard]] Point3 Bracketed(const Reference& r, const Leg& leg, const Cylinder& c, double lo,
                               double hi)
{
    const bool negative = SeamSideOf(c, OnLeg(r.Pair, leg, lo)).Side < 0.0;
    for (int i = 0; i < 100; ++i) {
        const double middle = 0.5 * (lo + hi);
        if (middle <= lo || middle >= hi) {
            break;
        }
        ((SeamSideOf(c, OnLeg(r.Pair, leg, middle)).Side < 0.0) == negative ? lo : hi) = middle;
    }
    return OnLeg(r.Pair, leg, 0.5 * (lo + hi));
}

/// A sample against a seam: on one side of it, exactly on its plane, or
/// behind the axis, where the plane is not the seam.
enum class Standing { Negative, Positive, OnPlane, Behind };

[[nodiscard]] Standing StandingOf(const Cylinder& c, const Point3& p)
{
    const SeamSide at = SeamSideOf(c, p);
    if (!at.Ahead) {
        return Standing::Behind;
    }
    if (at.Side == 0.0) {
        return Standing::OnPlane;
    }
    return at.Side > 0.0 ? Standing::Positive : Standing::Negative;
}

[[nodiscard]] bool OnASide(const Standing s)
{
    return s == Standing::Negative || s == Standing::Positive;
}

/// The crossings of a cylinder's seam along one closed curve of the
/// reference: every change of the seam's side between samples ahead of the
/// axis, walked once around from a sample that is on a side. Between two
/// neighbours of one leg the change is bracketed; otherwise the samples
/// between the two sides are on the seam's plane, or one point of two
/// legs, and the crossing is the sample before the change. A touch, with
/// the same side before and after, is no crossing.
void AddSeamCrossings(const Reference& r, const std::size_t loop, const Cylinder& c,
                      std::vector<Point3>& crossings)
{
    const std::vector<Sample>& samples = r.Samples[loop];
    const std::size_t count = samples.size();
    std::size_t start = 0;
    while (start < count && !OnASide(StandingOf(c, samples[start].Point))) {
        ++start;
    }
    if (start == count) {
        return;
    }
    Standing last = StandingOf(c, samples[start].Point);
    std::size_t lastStep = start;
    for (std::size_t step = start + 1; step <= start + count; ++step) {
        const Sample& here = samples[step % count];
        const Standing now = StandingOf(c, here.Point);
        if (now == Standing::OnPlane) {
            continue;
        }
        if (OnASide(now) && OnASide(last) && now != last) {
            const Sample& before = samples[lastStep % count];
            const bool oneLeg = lastStep + 1 == step && before.Leg == here.Leg;
            crossings.push_back(
                oneLeg ? Bracketed(r, r.Loops[loop][here.Leg], c, before.Fraction, here.Fraction)
                       : samples[(step - 1) % count].Point);
        }
        last = now;
        lastStep = step;
    }
}

[[nodiscard]] std::vector<Point3> SeamCrossingsOf(const Reference& r, const Cylinder& c)
{
    std::vector<Point3> crossings;
    for (std::size_t loop = 0; loop < r.Samples.size(); ++loop) {
        AddSeamCrossings(r, loop, c, crossings);
    }
    return crossings;
}

/// The number of points of a seam on the other cylinder, by the sign of
/// the product alone: the seam of A at the world's height y of its ruling,
/// or the seam of B, the line along x through (y, z) of its ruling, which
/// meets A where x² = (3 - y)(3 + y).
[[nodiscard]] std::size_t CountOf(const double product)
{
    if (product > 0.0) {
        return 2;
    }
    return product == 0.0 ? 1 : 0;
}

[[nodiscard]] std::size_t SeamPointsOfA(const World& w, const Cylinder& a)
{
    const double y = radiusA * a.XAxis.Y;
    return CountOf((y - w.Low()) * (w.High() - y));
}

[[nodiscard]] std::size_t SeamPointsOfB(const World& w, const Cylinder& b)
{
    const double y = w.Offset + w.RadiusB * b.XAxis.Y;
    return CountOf((radiusA - y) * (radiusA + y));
}

/// A, its seam at the world's angle whose cosine and sine are given.
[[nodiscard]] Cylinder FrameOfA(const double cosine, const double sine, const bool rightHanded)
{
    return {Point3{}, radiusA, Vec3{cosine, sine, 0.0}, Vec3{-sine, cosine, 0.0}, rightHanded};
}

/// B, its seam at the angle about x, from y towards z, whose cosine and
/// sine are given.
[[nodiscard]] Cylinder FrameOfB(const World& w, const double cosine, const double sine,
                                const bool rightHanded)
{
    return {Point3{0.0, w.Offset, 0.0}, w.RadiusB, Vec3{0.0, cosine, sine},
            Vec3{0.0, -sine, cosine}, rightHanded};
}

/// The smallest distance between two of the points.
[[nodiscard]] double ClosestPair(const std::vector<Point3>& points)
{
    double closest = 1e300;
    for (std::size_t i = 0; i < points.size(); ++i) {
        for (std::size_t j = i + 1; j < points.size(); ++j) {
            closest = std::min(closest, Distance(points[i], points[j]));
        }
    }
    return closest;
}

constexpr std::size_t perBranch = 256;

/// A branch against one cylinder's seam: its pcurve stays in the chart
/// without a jump, and the curve does not pass from one side of the seam
/// to the other, ahead of the axis, by more than the reach on both sides.
void RequireOneChart(const Stereon::IntersectionBranch& branch, const Cylinder& c,
                     const Stereon::Curve2d& pcurve, const double reach)
{
    double side = 0.0;
    double previous = Eval(pcurve, branch.Range.GetLower()).X;
    for (const double t : IntersectTest::SamplesOf(branch.Range, perBranch)) {
        const double angle = Eval(pcurve, t).X;
        REQUIRE(angle >= -1e-9);
        REQUIRE(angle <= Stereon::CirclePeriod + 1e-9);
        REQUIRE(std::abs(angle - previous) < pi);
        previous = angle;
        const SeamSide at = SeamSideOf(c, Eval(branch.Curve3d, t));
        if (!at.Ahead) {
            side = 0.0;
            continue;
        }
        if (std::abs(at.Side) <= reach) {
            continue;
        }
        const double sign = at.Side > 0.0 ? 1.0 : -1.0;
        REQUIRE((side == 0.0 || side == sign));
        side = sign;
    }
}

/// The ends of the branches chain into closed curves: two ends at every
/// point where branches meet, four at each of the expected nodes.
void RequireClosedChains(const Stereon::SurfaceSurfaceIntersection& result, const std::size_t nodes)
{
    std::vector<Point3> clusters;
    std::vector<std::size_t> degree;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
            const Point3 p = Eval(branch.Curve3d, t);
            std::size_t cluster = 0;
            while (cluster < clusters.size() && Distance(clusters[cluster], p) > 1e-7) {
                ++cluster;
            }
            if (cluster == clusters.size()) {
                clusters.push_back(p);
                degree.push_back(0);
            }
            ++degree[cluster];
        }
    }
    std::size_t four = 0;
    for (const std::size_t d : degree) {
        REQUIRE((d == 2 || d == 4));
        four += d == 4 ? 1U : 0U;
    }
    REQUIRE(four == nodes);
}

/// The whole contract of an answer against the reference.
void RequireContract(const Stereon::SurfaceSurfaceIntersection& result, const Reference& r,
                     const Cylinder& a, const Cylinder& b, const bool swapped,
                     const std::vector<Point3>& crossings, const std::size_t nodes)
{
    const double resolution = ctx.GetResolution();
    const Surface erasedA{a};
    const Surface erasedB{b};
    REQUIRE(result.Points.empty());
    REQUIRE(result.Overlaps.empty());
    REQUIRE_FALSE(result.Branches.empty());
    double length = 0.0;
    std::vector<Point3> ends;
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        REQUIRE(branch.Contact == ContactType::Crossing);
        REQUIRE(branch.IsCertified());
        REQUIRE(branch.Bound <= resolution);
        const Stereon::Curve2d& onA = swapped ? branch.OnB : branch.OnA;
        const Stereon::Curve2d& onB = swapped ? branch.OnA : branch.OnB;
        Point3 previous = Eval(branch.Curve3d, branch.Range.GetLower());
        for (const double t : IntersectTest::SamplesOf(branch.Range, perBranch)) {
            const Point3 p = Eval(branch.Curve3d, t);
            REQUIRE(DistanceFrom(a, p) <= branch.Bound + 1e-12);
            REQUIRE(DistanceFrom(b, p) <= branch.Bound + 1e-12);
            REQUIRE(Distance(p, Eval(erasedA, Eval(onA, t))) <= branch.Bound);
            REQUIRE(Distance(p, Eval(erasedB, Eval(onB, t))) <= branch.Bound);
            length += Distance(previous, p);
            previous = p;
        }
        RequireOneChart(branch, a, onA, 2.0 * resolution);
        RequireOneChart(branch, b, onB, 2.0 * resolution);
        ends.push_back(Eval(branch.Curve3d, branch.Range.GetLower()));
        ends.push_back(Eval(branch.Curve3d, branch.Range.GetUpper()));
    }
    RequireClosedChains(result, nodes);
    REQUIRE(std::abs(length - r.Length) <= 1e-2);
    // Crossings a few resolutions apart are each a cut of their own, where
    // one branch ends and the next begins.
    if (ClosestPair(crossings) > 4.0 * resolution) {
        for (const Point3& crossing : crossings) {
            const auto near = std::ranges::count_if(ends, [&](const Point3& end) {
                return Distance(end, crossing) <= 1.5 * resolution;
            });
            REQUIRE(near == 2);
        }
    }
}

/// A refusal says which of the handler's limits it met: the budget of the
/// seam isolation, two seam points in one unresolved cell, two cuts within
/// the resolution of each other, a seam within rounding of a fold, or a
/// piece that would not certify within the sample budget.
[[nodiscard]] bool IsNamedRefusal(const std::string_view message)
{
    constexpr std::array<std::string_view, 5> limits{
        "exhausted its interval budget", "cannot be told apart within one unresolved cell",
        "within the fit resolution of each other", "within rounding of a fold",
        "exceeds the sample budget"};
    return std::ranges::any_of(limits, [&](const std::string_view limit) {
        return message.find(limit) != std::string_view::npos;
    });
}

/// One pair of frames in one order: the independent count of each seam's
/// crossings, then the answer's contract, or a refusal, which names its
/// limit and which only crossings within a few resolutions of each other
/// excuse. Whether it was answered.
bool Check(const Reference& r, const Cylinder& a, const Cylinder& b, const bool swapped,
           const std::size_t nodes = 0)
{
    const std::vector<Point3> onSeamA = SeamCrossingsOf(r, a);
    const std::vector<Point3> onSeamB = SeamCrossingsOf(r, b);
    // The dense reference and the sign of the product agree on each seam: a
    // tangent seam touches and does not cross.
    REQUIRE(onSeamA.size() == (SeamPointsOfA(r.Pair, a) == 2 ? 2U : 0U));
    REQUIRE(onSeamB.size() == (SeamPointsOfB(r.Pair, b) == 2 ? 2U : 0U));
    std::vector<Point3> crossings = onSeamA;
    crossings.insert(crossings.end(), onSeamB.begin(), onSeamB.end());

    const Surface erasedA{a};
    const Surface erasedB{b};
    const auto result =
        swapped ? Intersect(erasedB, erasedA, ctx) : Intersect(erasedA, erasedB, ctx);
    INFO(IntersectTest::ErrorText(result));
    if (!result) {
        REQUIRE(IsNamedRefusal(result.error().GetMessage()));
        REQUIRE((result.error().GetCode() == ErrorCode::IntersectionFailed ||
                 result.error().GetCode() == ErrorCode::CertificationFailed));
        REQUIRE(ClosestPair(crossings) <= 4.0 * ctx.GetResolution());
        return false;
    }
    RequireContract(*result, r, a, b, swapped, crossings, nodes);
    return true;
}

/// The offsets from a fold, in radians, on either side of it and down to
/// nothing: by powers of two, close together where the two crossings come
/// within the resolution of each other and where the offset falls under
/// the rounding of the fold's own angle. Fewer on the positive side, where
/// a seam inside the fold is refused once its crossings are that close,
/// which costs both parametrisations' budgets.
[[nodiscard]] std::vector<double> OffsetsToZero()
{
    std::vector<double> offsets{0.0};
    for (const int power : {4, 16, 28, 34, 36, 38, 40, 46, 53, 62}) {
        if (power != 38 && power != 46 && power != 53) {
            offsets.push_back(std::ldexp(1.0, -power));
        }
        offsets.push_back(-std::ldexp(1.0, -power));
    }
    return offsets;
}

/// B reaching past A, in one loop whose folds are A's rulings in the plane
/// y = 0; and B passing through A, in two loops with a fold each there.
constexpr World oneLoop{2.0, 2.0};
constexpr World twoLoops{1.0, 1.0};

} // namespace

TEST_CASE("the parametrising cylinder's seam at and beside a fold is cut where it crosses",
          "[intersect][cylinder][skew][regression][issue195]")
{
    // A's seam, the ruling at the world's angle e or pi - e, is tangent to
    // B at e = 0, the fold itself; crosses the curve twice just inside the
    // fold for e > 0, the crossings 2 sqrt(y (High - y)) apart with
    // y = 3 sin e; and misses it for e < 0. B's seam meets nothing where B
    // reaches past A, on the far side of its axis, and crosses each loop at
    // an ordinary place where B passes through. A has the larger radius,
    // so it is the cylinder whose angle parametrises the curve in either
    // order. The handedness and the order of the operands alternate along
    // the offsets, both of each at either fold.
    for (const World& w : {oneLoop, twoLoops}) {
        const Reference r = ReferenceOf(w);
        const bool through = w.High() < radiusA;
        std::size_t answered = 0;
        std::size_t run = 0;
        std::size_t index = 0;
        for (const double e : OffsetsToZero()) {
            ++index;
            for (const double mirror : {1.0, -1.0}) {
                const bool rightHanded = (index % 4 < 2) == (mirror > 0.0);
                const Cylinder a = FrameOfA(mirror * std::cos(e), std::sin(e), rightHanded);
                const Cylinder b = through ? FrameOfB(w, std::cos(1.0), std::sin(1.0), !rightHanded)
                                           : FrameOfB(w, 1.0, 0.0, !rightHanded);
                REQUIRE(SeamPointsOfA(w, a) == (e > 0.0 ? 2U : (e == 0.0 ? 1U : 0U)));
                REQUIRE(SeamPointsOfB(w, b) == (through ? 2U : 0U));
                const bool swapped = (index % 3 == 0) == (mirror > 0.0);
                INFO("B radius " << w.RadiusB << ", offset " << e << ", mirror " << mirror
                                 << ", right-handed " << rightHanded << ", swapped " << swapped);
                ++run;
                const bool ok = Check(r, a, b, swapped);
                answered += ok ? 1U : 0U;
                // Continuity at zero: the tangent seam and every seam that
                // misses the curve are answered, however near the fold.
                REQUIRE((ok || e > 0.0));
            }
        }
        // The refusals are the few offsets whose crossings are within
        // resolutions of each other.
        REQUIRE(answered * 4 >= run * 3);
    }
}

TEST_CASE("the other cylinder's seam at and beside the folds is cut where it crosses",
          "[intersect][cylinder][skew][regression][issue195]")
{
    // B's seam at the angle pi + e about its axis is the line along x
    // through (Offset - r_B cos e, -r_B sin e): at e = 0 the x axis, which
    // passes through both folds, (±3, 0, 0), and beside it a line crossing
    // the curve 2 r_B e from each fold, on one sheet or the other. A's
    // seam crosses the curve at an ordinary place.
    for (const World& w : {oneLoop, twoLoops}) {
        const Reference r = ReferenceOf(w);
        std::size_t index = 0;
        for (const double e : OffsetsToZero()) {
            const bool rightHanded = ++index % 2 == 0;
            const Cylinder a = FrameOfA(std::cos(0.3), std::sin(0.3), !rightHanded);
            const Cylinder b = FrameOfB(w, -std::cos(e), -std::sin(e), rightHanded);
            REQUIRE(SeamPointsOfA(w, a) == 2);
            REQUIRE(SeamPointsOfB(w, b) == 2);
            for (const bool swapped : {false, true}) {
                INFO("B radius " << w.RadiusB << ", offset " << e << ", right-handed "
                                 << rightHanded << ", swapped " << swapped);
                // The crossings are far apart: nothing excuses a refusal.
                REQUIRE(Check(r, a, b, swapped));
            }
        }
    }
}

TEST_CASE("two crossings of a seam closer than the former scan's step are both cut",
          "[intersect][cylinder][skew][regression][issue195]")
{
    // B's seam is the line along x at y = 3 - d, which meets A at
    // x = ±sqrt(d (6 - d)), both on one sheet beside the top of the loop:
    // 0.02 apart at d = 1e-4 and 5e-5 at d = 1e-10, where the 256 samples
    // of the former scan, 0.07 apart along the loop, stepped over both. At
    // d = 1e-13 they are 1.5e-6 apart, under the resolution: no piece
    // between them could be fitted, so the pair is refused, or answered by
    // the parametrisation in which that arc lies within the resolution of
    // a branch's end, which stands for both cuts. At d = 0 the seam is
    // tangent to A at the top of the loop, exactly, and for d < 0 it passes
    // outside: the curve touches the seam or comes within d of it, crosses
    // nothing, and is answered.
    const Reference r = ReferenceOf(oneLoop);
    // The seams that touch or nearly touch run the isolation to its budget
    // in the first parametrisation, so those are asked in one order each.
    std::size_t index = 0;
    for (const double d : {1e-4, 1e-10, 1e-13, 0.0, -1e-15}) {
        const double cosine = (1.0 - d) / 2.0;
        const double sheet = ++index % 2 == 0 ? 1.0 : -1.0;
        const bool rightHanded = index % 4 < 2;
        const double sine = sheet * std::sqrt((1.0 - cosine) * (1.0 + cosine));
        const Cylinder a = FrameOfA(std::cos(0.3), std::sin(0.3), rightHanded);
        const Cylinder b = FrameOfB(oneLoop, cosine, sine, rightHanded);
        REQUIRE(SeamPointsOfB(oneLoop, b) == (d > 0.0 ? 2U : (d == 0.0 ? 1U : 0U)));
        const bool costly = d >= 0.0 && d < 1e-6;
        for (const bool swapped : {false, true}) {
            if (costly && swapped != (index % 2 == 0)) {
                continue;
            }
            INFO("depth " << d << ", sheet " << sheet << ", right-handed " << rightHanded
                          << ", swapped " << swapped);
            const bool ok = Check(r, a, b, swapped);
            REQUIRE((ok || d == 1e-13));
        }
    }
}

TEST_CASE("both seams at one fold are cut or refused, never joined across a seam",
          "[intersect][cylinder][skew][regression][issue195]")
{
    // A's seam beside the fold at (3, 0, 0) while B's passes through or
    // beside the same fold: cuts of both seams within the resolution of
    // each other and of the track's own end.
    const Reference r = ReferenceOf(oneLoop);
    std::size_t answered = 0;
    std::size_t index = 0;
    for (const double ea : {0.0, 0x1p-30, -0x1p-30, 0x1p-50, -0x1p-50, 0x1p-8}) {
        for (const double eb : {0.0, 0x1p-30, -0x1p-30, 0x1p-50, 0x1p-8}) {
            for (const bool rightHanded : {true, false}) {
                const Cylinder a = FrameOfA(std::cos(ea), std::sin(ea), rightHanded);
                const Cylinder b = FrameOfB(oneLoop, -std::cos(eb), -std::sin(eb), !rightHanded);
                const bool swapped = ++index % 2 == 0;
                INFO("A offset " << ea << ", B offset " << eb << ", right-handed " << rightHanded
                                 << ", swapped " << swapped);
                answered += Check(r, a, b, swapped) ? 1U : 0U;
            }
        }
    }
    REQUIRE(answered > 0);
}

TEST_CASE("the seams of internally tangent cylinders beside a lobe's fold are cut where they cross",
          "[intersect][cylinder][skew][regression][issue195]")
{
    // B of radius 1 inside A, touching it at the node (0, 3, 0): two lobes
    // whose folds are A's rulings at y = 1. A's seam beside the fold of
    // one lobe, inside it or outside; B's seam along x through the folds
    // of both, at y = 1, or beside them.
    const World lobes{2.0, 1.0};
    const Reference r = ReferenceOf(lobes);
    const double fold = std::asin(1.0 / 3.0);
    for (const double e : {1e-2, -1e-2, 1e-5, -1e-5, 1e-8, -1e-8, 1e-11, -1e-11}) {
        for (const bool rightHanded : {true, false}) {
            const Cylinder a = FrameOfA(std::cos(fold + e), std::sin(fold + e), rightHanded);
            const Cylinder b = FrameOfB(lobes, -std::cos(e), -std::sin(e), !rightHanded);
            REQUIRE(SeamPointsOfA(lobes, a) == (e > 0.0 ? 2U : 0U));
            REQUIRE(SeamPointsOfB(lobes, b) == 2);
            for (const bool swapped : {false, true}) {
                INFO("offset " << e << ", right-handed " << rightHanded << ", swapped " << swapped);
                const bool ok = Check(r, a, b, swapped, 1);
                REQUIRE((ok || e > 0.0));
            }
        }
    }
}
