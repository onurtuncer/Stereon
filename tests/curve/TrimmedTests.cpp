// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Trimmed curves: a basis restricted to a parameter range. A trimmed curve
// must evaluate as its basis on the range, project to the closest point of
// the piece rather than of the whole basis, and convert exactly.

#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

import Stereon.Curve;

#include "Formatting.hpp"

using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Ellipse;
using Stereon::ErrorCode;
using Stereon::Hyperbola;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::Parabola;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Transform;
using Stereon::TrimmedCurve;
using Stereon::TrimmedCurve2d;
using Stereon::Vec2;
using Stereon::Vec3;

namespace {

constexpr double Tolerance = 1e-11;
constexpr double Pi = std::numbers::pi;

Circle UnitCircle()
{
    return {{0.0, 0.0, 0.0}, 1.0, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
}

TrimmedCurve Trim(const Curve& basis, const Interval range)
{
    return *Stereon::MakeTrimmedCurve(basis, range);
}

/// The parameter of the closest of `samples` points of the curve to p.
template <typename C, typename P>
double Brute(const C& curve, const P& p)
{
    const Interval range = Domain(curve);
    constexpr int Samples = 4000;
    double best = range.GetLower();
    double bestDistance = DistanceSquared(Eval(curve, best), p);
    for (int k = 1; k <= Samples; ++k) {
        const double t = range.GetLower() + Width(range) * k / Samples;
        const double d = DistanceSquared(Eval(curve, t), p);
        if (d < bestDistance) {
            best = t;
            bestDistance = d;
        }
    }
    return best;
}

/// True when Project's answer is at least as close as the sampled minimum.
template <typename C, typename P>
bool ProjectsClosest(const C& curve, const P& p)
{
    const double t = Project(curve, p);
    const Interval range = Domain(curve);
    if (t < range.GetLower() || t > range.GetUpper()) {
        return false;
    }
    const double found = Distance(Eval(curve, t), p);
    const double sampled = Distance(Eval(curve, Brute(curve, p)), p);
    return found <= sampled + 1e-9;
}

} // namespace

TEST_CASE("a trimmed curve evaluates as its basis on its range", "[curve][trimmed]")
{
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{0.5, 2.0});
    REQUIRE(Domain(arc) == Interval{0.5, 2.0});
    for (const double t : {0.5, 1.0, 1.7, 2.0}) {
        REQUIRE(Distance(Eval(arc, t), Eval(UnitCircle(), t)) <= Tolerance);
        REQUIRE(Length(Deriv(arc, t, 1) - Deriv(UnitCircle(), t, 1)) <= Tolerance);
        REQUIRE(Length(Deriv(arc, t, 2) - Deriv(UnitCircle(), t, 2)) <= Tolerance);
    }
    STATIC_REQUIRE(Stereon::CurveGeometry<TrimmedCurve>);
    STATIC_REQUIRE(Stereon::Curve2dGeometry<TrimmedCurve2d>);
}

TEST_CASE("a trimmed curve projects onto its piece, not its basis", "[curve][trimmed]")
{
    // The quarter arc from angle 0 to pi/2. The centre's neighbourhood near
    // the chord projects to the interior; points beyond the ends to an end.
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{0.0, Pi / 2.0});
    REQUIRE(std::abs(Project(arc, Point3{0.3, 0.3, 0.0}) - Pi / 4.0) <= 1e-9);
    REQUIRE(std::abs(Project(arc, Point3{2.0, 2.0, 1.0}) - Pi / 4.0) <= 1e-9);
    REQUIRE(Project(arc, Point3{-1.0, -0.2, 0.0}) == Pi / 2.0);
    REQUIRE(Project(arc, Point3{0.5, -3.0, 0.0}) == 0.0);

    // A segment of a line: the foot of the perpendicular, or an end.
    const TrimmedCurve segment =
        Trim(Curve(Line{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}), Interval{-1.0, 2.0});
    REQUIRE(std::abs(Project(segment, Point3{0.75, 4.0, -1.0}) - 0.75) <= Tolerance);
    REQUIRE(Project(segment, Point3{9.0, 0.0, 0.0}) == 2.0);
    REQUIRE(Project(segment, Point3{-9.0, 1.0, 0.0}) == -1.0);

    // The lower half of an ellipse: the global closest point of (0, 0.1) is
    // the top vertex, outside the range, and clamping it would give an end
    // at distance 3; the bottom vertex, at distance 1.1, is the answer
    // (docs/trimmed.rst).
    const TrimmedCurve lower =
        Trim(Curve(Ellipse{{0.0, 0.0, 0.0}, 3.0, 1.0, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}),
             Interval{Pi, 2.0 * Pi});
    REQUIRE(std::abs(Project(lower, Point3{0.0, 0.1, 0.0}) - 1.5 * Pi) <= 1e-9);
}

TEST_CASE("trimmed curves project to the closest point of the piece", "[curve][trimmed][property]")
{
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 1.0, -2.0}, 0.9);
    const Vec3 x = Transformed(Vec3{1.0, 0.0, 0.0}, tilt);
    const Vec3 y = Transformed(Vec3{0.0, 1.0, 0.0}, tilt);
    const Point3 c{1.0, -2.0, 0.5};
    const std::vector<TrimmedCurve> curves{
        Trim(Curve(Circle{c, 2.0, x, y}), Interval{5.0, 7.0}), // across the seam
        Trim(Curve(Ellipse{c, 3.0, 1.0, x, y}), Interval{0.4, 3.9}),
        Trim(Curve(Hyperbola{c, 2.0, 1.0, x, y}), Interval{-1.2, 0.8}),
        Trim(Curve(Parabola{c, 0.5, x, y}), Interval{-2.0, 1.0}),
        Trim(Curve(*Stereon::MakeBSplineCurve(
                 3, {0, 0, 0, 0, 1, 2, 2, 2, 2},
                 {{0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 1, 0}, {4, 0, 2}}, {})),
             Interval{0.3, 1.6}),
    };
    rc::prop("Project agrees with dense sampling of the piece",
             [&](const double px, const double py, const double pz) {
                 const Point3 p{std::fmod(px, 6.0), std::fmod(py, 6.0), std::fmod(pz, 6.0)};
                 for (const TrimmedCurve& curve : curves) {
                     RC_ASSERT(ProjectsClosest(curve, p));
                     RC_ASSERT(ProjectsClosest(Reversed(curve), p));
                 }
             });
}

TEST_CASE("a circle trimmed across its seam is one arc", "[curve][trimmed]")
{
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{5.0, 7.0});
    REQUIRE(Distance(Eval(arc, 7.0), Point3{std::cos(7.0), std::sin(7.0), 0.0}) <= Tolerance);
    // The point at angle 0 lies inside the arc, at parameter 2 pi.
    REQUIRE(std::abs(Project(arc, Point3{3.0, 0.0, 0.0}) - 2.0 * Pi) <= 1e-9);
    // The arc's length is its angle.
    REQUIRE(std::abs(Stereon::ArcLength(arc, Domain(arc))->Length - 2.0) <= 1e-10);
}

TEST_CASE("a reversed trimmed curve runs its range backwards", "[curve][trimmed]")
{
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{0.5, 2.0});
    const TrimmedCurve back = Reversed(arc);
    REQUIRE(back.IsReversed());
    REQUIRE(Domain(back) == Domain(arc));
    REQUIRE(back.GetPiece() == arc.GetPiece());
    for (const double t : {0.5, 0.9, 2.0}) {
        REQUIRE(Distance(Eval(back, t), Eval(arc, 2.5 - t)) <= Tolerance);
        REQUIRE(Length(Deriv(back, t, 1) + Deriv(arc, 2.5 - t, 1)) <= Tolerance);
        REQUIRE(Length(Deriv(back, t, 2) - Deriv(arc, 2.5 - t, 2)) <= Tolerance);
    }
    const Point3 p{0.2, 0.4, 0.0};
    REQUIRE(std::abs(Project(back, p) - (2.5 - Project(arc, p))) <= Tolerance);
    REQUIRE(!Reversed(back).IsReversed());
}

TEST_CASE("a trimmed curve converts, moves and bounds as its piece", "[curve][trimmed]")
{
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{0.0, Pi / 2.0});

    // The NURBS form is the quarter arc, run as the trimmed curve runs.
    const BSplineCurve nurbs = *Stereon::ToBSpline(arc);
    REQUIRE(Distance(Eval(nurbs, Domain(nurbs).GetLower()), Point3{1.0, 0.0, 0.0}) <= Tolerance);
    REQUIRE(Distance(Eval(nurbs, Domain(nurbs).GetUpper()), Point3{0.0, 1.0, 0.0}) <= Tolerance);
    const BSplineCurve back = *Stereon::ToBSpline(Reversed(arc));
    REQUIRE(Distance(Eval(back, Domain(back).GetLower()), Point3{0.0, 1.0, 0.0}) <= Tolerance);

    // The box holds the piece and not the rest of the circle.
    const Stereon::Box3 box = BoundingBox(arc);
    REQUIRE(box.X.GetLower() >= -Tolerance);
    REQUIRE(box.Y.GetLower() >= -Tolerance);
    REQUIRE(box.X.GetUpper() >= 1.0 - Tolerance);

    const Transform shift = Stereon::MakeTranslation(Vec3{0.0, 0.0, 3.0});
    const TrimmedCurve moved = Transformed(Reversed(arc), shift);
    REQUIRE(moved.IsReversed());
    REQUIRE(Distance(Eval(moved, 0.0), Point3{0.0, 1.0, 3.0}) <= Tolerance);
}

TEST_CASE("arc length of a trimmed B-spline splits at its knots", "[curve][trimmed]")
{
    const BSplineCurve spline =
        *Stereon::MakeBSplineCurve(1, {0, 0, 1, 2, 2}, {{0, 0, 0}, {3, 0, 0}, {3, 4, 0}}, {});
    const TrimmedCurve piece = Trim(Curve(spline), Interval{0.5, 1.5});
    REQUIRE(std::abs(Stereon::ArcLength(piece, Domain(piece))->Length - 3.5) <= 1e-12);
    REQUIRE(std::abs(Stereon::ArcLength(Reversed(piece), Domain(piece))->Length - 3.5) <= 1e-12);

    // Reversed over its whole domain, a B-spline is evaluated at its ends
    // exactly: the reflection of an end is the other end, not a sum that
    // rounds outside the basis's domain.
    const TrimmedCurve whole = Reversed(Trim(Curve(spline), Interval{0.1, 2.0}));
    REQUIRE(Eval(whole, 0.1) == Point3{3.0, 4.0, 0.0});
    REQUIRE(Eval(whole, 2.0) == Eval(spline, 0.1));
}

TEST_CASE("a trimmed curve is wrapped as any curve", "[curve][trimmed][erased]")
{
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{0.0, Pi});
    const Curve wrapped(arc);
    REQUIRE(Domain(wrapped) == Interval{0.0, Pi});
    REQUIRE(std::abs(Project(wrapped, Point3{0.1, -0.5, 0.0}) - 0.0) <= Tolerance);

    // Trimming it again: the wrapped trimmed curve is not a built-in type,
    // but its piece of a circle has an exact form, so the inner trim has one
    // too (audit 2026-10-20, F2) and projects onto it.
    const TrimmedCurve inner = Trim(wrapped, Interval{0.5, 2.5});
    REQUIRE(inner.GetPiece() != nullptr);
    const BSplineCurve nurbs = *Stereon::ToBSpline(inner);
    REQUIRE(Domain(nurbs) == Interval{0.5, 2.5});
    REQUIRE(Distance(Eval(nurbs, 0.5), Eval(UnitCircle(), 0.5)) <= Tolerance);
    REQUIRE(Distance(Eval(nurbs, 2.5), Eval(UnitCircle(), 2.5)) <= Tolerance);
    REQUIRE(std::abs(Project(inner, Point3{0.0, 3.0, 0.0}) - Pi / 2.0) <= 1e-9);
    REQUIRE(Project(inner, Point3{1.0, -0.1, 0.0}) == 0.5);
}

namespace {

/// The NURBS form of `curve` over `range` through the erased wrapper, as a
/// sweep or a trimmed surface asks for it.
BSplineCurve ErasedForm(const TrimmedCurve& curve, const Interval range)
{
    const auto form = Stereon::ToBSpline(Curve(curve), range);
    REQUIRE(form.has_value());
    return *form;
}

/// The form agrees with the trimmed curve at every knot in its domain, and
/// everywhere when the basis keeps its parameter.
void RequireFormOf(const TrimmedCurve& curve, const Interval range)
{
    const BSplineCurve form = ErasedForm(curve, range);
    REQUIRE(Domain(form) == range);
    for (const double k : form.Knots) {
        REQUIRE(Distance(Eval(form, k), Eval(curve, k)) <= Tolerance);
    }
    if (Stereon::Detail::KeepsParameter(Curve(curve))) {
        for (int i = 0; i <= 16; ++i) {
            const double t = range.GetLower() + Width(range) * i / 16.0;
            REQUIRE(Distance(Eval(form, t), Eval(curve, t)) <= Tolerance);
        }
    }
    // Run as the curve runs: the tangents agree in direction at the ends.
    for (const double t : {range.GetLower(), range.GetUpper()}) {
        REQUIRE(Dot(Deriv(form, t, 1), Deriv(curve, t, 1)) > 0.0);
    }
}

} // namespace

TEST_CASE("a trimmed curve converts through the erased wrapper over any subrange",
          "[curve][trimmed][conversion][audit]")
{
    // Audit 2026-10-20, F2: ToBSpline(const Curve&, Interval) knew only the
    // built-in types, so a sweep or trimmed surface of a trimmed profile
    // found no exact form and fell back to approximation. Every basis kind,
    // forwards and reversed, over its whole range and over subranges.
    const Transform tilt = Stereon::MakeRotation(Vec3{1.0, 1.0, -2.0}, 0.9);
    const Vec3 x = Transformed(Vec3{1.0, 0.0, 0.0}, tilt);
    const Vec3 y = Transformed(Vec3{0.0, 1.0, 0.0}, tilt);
    const Point3 c{1.0, -2.0, 0.5};
    const std::vector<TrimmedCurve> curves{
        Trim(Curve(Line{c, Vec3{0.5, -1.0, 2.0}}), Interval{-1.0, 3.0}),
        Trim(Curve(Circle{c, 2.0, x, y}), Interval{5.0, 7.0}), // across the seam
        Trim(Curve(Ellipse{c, 3.0, 1.0, x, y}), Interval{0.4, 3.9}),
        Trim(Curve(Hyperbola{c, 2.0, 1.0, x, y}), Interval{-1.2, 0.8}),
        Trim(Curve(Parabola{c, 0.5, x, y}), Interval{-2.0, 1.0}),
        Trim(Curve(*Stereon::MakeBSplineCurve(
                 3, {0, 0, 0, 0, 1, 2, 2, 2, 2},
                 {{0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 1, 0}, {4, 0, 2}}, {1, 2, 1, 0.5, 1})),
             Interval{0.3, 1.6}),
    };
    for (const TrimmedCurve& forwards : curves) {
        for (const TrimmedCurve& curve : {forwards, Reversed(forwards)}) {
            const Interval whole = Domain(curve);
            const double a = whole.GetLower();
            const double w = Width(whole);
            INFO("range " << a << " + " << w << (curve.IsReversed() ? ", reversed" : ""));
            // The whole range: the same form as the concrete conversion.
            const BSplineCurve concrete = *Stereon::ToBSpline(curve);
            const BSplineCurve erased = ErasedForm(curve, whole);
            REQUIRE(Domain(erased) == Domain(concrete));
            for (const double k : concrete.Knots) {
                REQUIRE(Distance(Eval(erased, k), Eval(concrete, k)) <= Tolerance);
            }
            RequireFormOf(curve, whole);
            // Subranges at either end and inside.
            RequireFormOf(curve, Interval{a, a + 0.4 * w});
            RequireFormOf(curve, Interval{a + 0.6 * w, a + w});
            RequireFormOf(curve, Interval{a + 0.15 * w, a + 0.7 * w});
        }
    }
}

TEST_CASE("the reversed piece of a line is the line, parameter for parameter",
          "[curve][trimmed][conversion][audit]")
{
    // C(t) = (t, 0, 0) on [0, 4], reversed: the trimmed curve at t is the
    // line at 4 - t. Its form over [1, 2] must have domain [1, 2] and run
    // from (3, 0, 0) to (2, 0, 0) on that parameter, not on [2, 3].
    const TrimmedCurve back =
        Reversed(Trim(Curve(Line{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}), Interval{0.0, 4.0}));
    const BSplineCurve form = ErasedForm(back, Interval{1.0, 2.0});
    REQUIRE(Domain(form) == Interval{1.0, 2.0});
    REQUIRE(Eval(form, 1.0) == Point3{3.0, 0.0, 0.0});
    REQUIRE(Eval(form, 2.0) == Point3{2.0, 0.0, 0.0});
    REQUIRE(Distance(Eval(form, 1.25), Point3{2.75, 0.0, 0.0}) <= Tolerance);
    REQUIRE(Stereon::Detail::KeepsParameter(Curve(back)));
}

TEST_CASE("nested trims convert to any depth, in space and in the plane",
          "[curve][trimmed][conversion][audit]")
{
    const Curve circle(UnitCircle());
    Curve nested = circle;
    TrimmedCurve innermost = Trim(circle, Interval{0.0, 6.0});
    for (int depth = 0; depth < 4; ++depth) {
        innermost = Trim(nested, Interval{0.5 + 0.25 * depth, 5.5 - 0.25 * depth});
        if (depth % 2 == 1) {
            innermost = Reversed(innermost);
        }
        nested = Curve(innermost);
        REQUIRE(innermost.GetPiece() != nullptr);
        RequireFormOf(innermost, Domain(innermost));
        RequireFormOf(innermost, Interval{2.0, 3.0});
    }
    REQUIRE_FALSE(Stereon::Detail::KeepsParameter(nested));
    REQUIRE(ProjectsClosest(innermost, Point3{0.3, -0.9, 0.2}));

    // A line trimmed twice keeps its parameter through both trims.
    const Curve line(Line{{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}});
    const TrimmedCurve twice =
        Reversed(Trim(Curve(Trim(line, Interval{-3.0, 3.0})), Interval{-1.0, 2.0}));
    REQUIRE(Stereon::Detail::KeepsParameter(Curve(twice)));
    RequireFormOf(twice, Interval{0.0, 1.5});

    // A plane curve trimmed twice has its piece too, lifted to z = 0.
    const Curve2d arc(Circle2d{{1.0, 1.0}, 2.0, {1.0, 0.0}, {0.0, 1.0}});
    const TrimmedCurve2d inner = *Stereon::MakeTrimmedCurve(
        Curve2d(*Stereon::MakeTrimmedCurve(arc, Interval{0.0, Pi})), Interval{0.5, 2.5});
    REQUIRE(inner.GetPiece() != nullptr);
    REQUIRE(ProjectsClosest(inner, Point2{1.2, 1.5}));
    REQUIRE(ProjectsClosest(Reversed(inner), Point2{-3.0, 0.0}));
    REQUIRE(Stereon::BoundingBox(inner).Y.GetLower() >= 1.0 + 2.0 * std::sin(0.5) - 1e-9);
}

TEST_CASE("the erased conversion refuses empty curves and ranges outside the trim",
          "[curve][trimmed][conversion][audit]")
{
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    REQUIRE(code(Stereon::ToBSpline(Curve{}, Interval{0.0, 1.0})) == ErrorCode::InvalidInput);
    REQUIRE_FALSE(Stereon::Detail::KeepsParameter(Curve{}));
    REQUIRE_FALSE(Stereon::Detail::KeepsParameter(Curve2d{}));
    const TrimmedCurve arc = Trim(Curve(UnitCircle()), Interval{1.0, 2.0});
    for (const TrimmedCurve& curve : {arc, Reversed(arc)}) {
        const Curve wrapped(curve);
        REQUIRE(code(Stereon::ToBSpline(wrapped, Interval{0.5, 1.5})) == ErrorCode::InvalidInput);
        REQUIRE(code(Stereon::ToBSpline(wrapped, Interval{1.5, 2.5})) == ErrorCode::InvalidInput);
        REQUIRE(code(Stereon::ToBSpline(wrapped, Interval{1.5, 1.5})) == ErrorCode::InvalidInput);
        REQUIRE(code(Stereon::ToBSpline(wrapped,
                                        Interval{1.0, std::numeric_limits<double>::infinity()})) ==
                ErrorCode::InvalidInput);
    }
}

TEST_CASE("trimmed plane curves project in the plane", "[curve][trimmed][2d]")
{
    const Circle2d circle{{1.0, 1.0}, 2.0, {1.0, 0.0}, {0.0, 1.0}};
    const TrimmedCurve2d arc = *Stereon::MakeTrimmedCurve(Curve2d(circle), Interval{0.0, Pi});
    REQUIRE(std::abs(Project(arc, Point2{1.2, 1.5}) - std::atan2(0.5, 0.2)) <= 1e-9);
    REQUIRE(Project(arc, Point2{4.0, -1.0}) == 0.0);
    REQUIRE(Distance(Eval(arc, Pi), Point2{-1.0, 1.0}) <= Tolerance);

    const TrimmedCurve2d segment =
        *Stereon::MakeTrimmedCurve(Curve2d(Line2d{{0.0, 0.0}, {0.0, 2.0}}), Interval{0.0, 1.0});
    REQUIRE(std::abs(Project(segment, Point2{5.0, 1.0}) - 0.5) <= Tolerance);
    REQUIRE(Project(Reversed(segment), Point2{5.0, 9.0}) == 0.0);
    const Stereon::Box2 box = BoundingBox(segment);
    REQUIRE(box.Y.GetUpper() <= 2.0 + Tolerance);
    REQUIRE(box.Y.GetLower() >= -Tolerance);
}

TEST_CASE("MakeTrimmedCurve validates the range", "[curve][trimmed][validation]")
{
    const auto code = [](const auto& result) { return result.error().GetCode(); };
    const Curve circle(UnitCircle());
    const Curve line(Line{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    const Curve spline(*Stereon::MakeBSplineCurve(1, {0, 0, 1, 1}, {{0, 0, 0}, {1, 0, 0}}, {}));
    constexpr double Inf = std::numeric_limits<double>::infinity();

    REQUIRE(code(Stereon::MakeTrimmedCurve(Curve{}, Interval{0.0, 1.0})) ==
            ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeTrimmedCurve(line, Interval{1.0, 1.0})) == ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeTrimmedCurve(line, Interval{0.0, Inf})) == ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeTrimmedCurve(circle, Interval{0.0, 7.0})) == ErrorCode::InvalidInput);
    REQUIRE(code(Stereon::MakeTrimmedCurve(spline, Interval{0.5, 1.5})) == ErrorCode::InvalidInput);

    REQUIRE(Stereon::MakeTrimmedCurve(circle, Interval{-3.0, 3.0}));
    REQUIRE(Stereon::MakeTrimmedCurve(circle, Interval{0.0, 2.0 * Pi}));
    REQUIRE(Stereon::MakeTrimmedCurve(line, Interval{-1e6, 1e6}));
    REQUIRE(Stereon::MakeTrimmedCurve(spline, Interval{0.0, 1.0}));
}
