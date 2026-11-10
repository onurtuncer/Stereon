// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Every built-in 3D curve as enclosed Bezier pieces over a bounded domain,
// with the map from a piece's linear parameter back to the curve's, and the
// clipping of an unbounded domain to a ball round a bounded partner. Shared
// by the curve-curve pair (SplineCurvePair.cpp) and curves against B-spline
// surfaces (CurveSplineSurface.cpp). Moved from SplineCurvePair.cpp.
//
// Included in the purview of the handlers' implementation units, after
// Bernstein.hpp, BezierArc.hpp and ConicArc.hpp, with <algorithm>, <cmath>,
// <initializer_list>, <numbers>, <string> and <vector> in their global
// module fragment. Everything is in an anonymous namespace, so each unit
// gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

// --- Every curve as enclosed Bezier pieces --------------------------------

/// How a piece's linear parameter between its ends maps to the curve's.
enum class ParameterMap : std::uint8_t { Same, Angle, Hyperbolic };

/// A curve over a bounded domain: its enclosed pieces, each [T0, T1] of
/// the curve's own parameter at its ends, and the map between them.
struct PieceCurve {
    std::vector<Segment> Segments;
    Interval Domain;
    ParameterMap Map = ParameterMap::Same;
    bool Closed = false;
    std::string Name;
};

/// The curve's parameter at u, linear on a piece between its ends.
[[maybe_unused]] [[nodiscard]] inline double CurveParameter(const PieceCurve& c, const double u)
{
    if (c.Map == ParameterMap::Same) {
        return u;
    }
    const auto piece =
        std::ranges::find_if(c.Segments, [&](const Segment& s) { return s.T0 <= u && u <= s.T1; });
    if (piece == c.Segments.end()) {
        return u;
    }
    const double s = (u - piece->T0) / (piece->T1 - piece->T0);
    return c.Map == ParameterMap::Angle ? AngleAt(piece->T0, piece->T1, s)
                                        : HyperbolicAt(piece->T0, piece->T1, s);
}

/// A one-parameter domain is searched on a piece around it.
[[maybe_unused]] [[nodiscard]] inline Interval BuildRange(const Interval& domain)
{
    return domain.GetLower() == domain.GetUpper()
               ? Interval{domain.GetLower() - 0.5, domain.GetUpper() + 0.5}
               : domain;
}

/// Cuts the range into `count` pieces and makes each with `piece`.
template <typename Piece>
[[maybe_unused]] [[nodiscard]] inline std::vector<Segment> Cut(const Interval& range,
                                                               const double count, Piece&& piece)
{
    const auto n = static_cast<std::size_t>(count);
    const double width = Width(range);
    std::vector<Segment> segments;
    for (std::size_t k = 0; k < n; ++k) {
        const double t0 = range.GetLower() + width * (static_cast<double>(k) / count);
        const double t1 = k + 1 == n
                              ? range.GetUpper()
                              : range.GetLower() + width * (static_cast<double>(k + 1) / count);
        Segment segment = piece(t0, t1);
        segment.T0 = t0;
        segment.T1 = t1;
        segments.push_back(std::move(segment));
    }
    return segments;
}

[[maybe_unused]] [[nodiscard]] inline HPoint PolePoint(const Point3& origin, const Interval& along,
                                                       const Vec3& axis, const Interval& across,
                                                       const Vec3& other)
{
    return {Interval{origin.X} + along * Interval{axis.X} + across * Interval{other.X},
            Interval{origin.Y} + along * Interval{axis.Y} + across * Interval{other.Y},
            Interval{origin.Z} + along * Interval{axis.Z} + across * Interval{other.Z},
            Interval{1.0}};
}

[[maybe_unused]] [[nodiscard]] inline PieceCurve Pieces(const Line& line, const Interval& domain)
{
    const Vec3 none{};
    const auto end = [&](const double t) {
        return PolePoint(line.Origin, Interval{t}, line.Direction, Interval{0.0}, none);
    };
    const Interval range = BuildRange(domain);
    return {{Segment{range.GetLower(),
                     range.GetUpper(),
                     {end(range.GetLower()), end(range.GetUpper())}}},
            domain,
            ParameterMap::Same,
            false,
            "Line"};
}

/// A parabola's arc is the quadratic Bezier whose middle pole is the blossom
/// of its two end parameters; it keeps the parabola's parameter.
[[maybe_unused]] [[nodiscard]] inline PieceCurve Pieces(const Parabola& parabola,
                                                        const Interval& domain)
{
    const Interval scale = Interval{1.0} / (Interval{4.0} * Interval{parabola.FocalLength});
    const auto pole = [&](const double r, const double s) {
        const Interval a{r};
        const Interval b{s};
        return PolePoint(parabola.Vertex, a * b * scale, parabola.XAxis, (a + b) * Interval{0.5},
                         parabola.YAxis);
    };
    const Interval range = BuildRange(domain);
    const double lower = range.GetLower();
    const double upper = range.GetUpper();
    return {{Segment{lower, upper, {pole(lower, lower), pole(lower, upper), pole(upper, upper)}}},
            domain,
            ParameterMap::Same,
            false,
            "Parabola"};
}

/// Whether the domain is a whole turn, whose ends are no ends: the curve
/// runs on through the seam.
[[maybe_unused]] [[nodiscard]] inline bool WholeTurn(const Interval& domain)
{
    return Width(domain) >= 2.0 * std::numbers::pi * (1.0 - 0x1p-40);
}

template <typename Conic>
[[maybe_unused]] [[nodiscard]] inline PieceCurve
ConicPieces(const Conic& conic, const ConicData& data, const Interval& domain, std::string name)
{
    static_cast<void>(conic);
    const Interval range = BuildRange(domain);
    const double count = std::max(1.0, std::ceil(Width(range) / (0.5 * std::numbers::pi)));
    return {Cut(range, count,
                [&](const double t0, const double t1) { return ArcSegment(data, t0, t1); }),
            domain, ParameterMap::Angle, WholeTurn(domain), std::move(name)};
}

/// The widest stretch of a hyperbola's parameter on one piece: its arc's
/// weight, cosh of half of it, stays moderate.
constexpr double HyperbolaPiece = 1.0;

[[maybe_unused]] [[nodiscard]] inline PieceCurve Pieces(const Hyperbola& h, const Interval& domain)
{
    const ConicData data{h.Center, Times(h.MajorRadius, h.XAxis), Times(h.MinorRadius, h.YAxis)};
    const Interval range = BuildRange(domain);
    const double count = std::max(1.0, std::ceil(Width(range) / HyperbolaPiece));
    return {Cut(range, count,
                [&](const double t0, const double t1) { return HyperbolicSegment(data, t0, t1); }),
            domain, ParameterMap::Hyperbolic, false, "Hyperbola"};
}

[[maybe_unused]] [[nodiscard]] inline bool AllFinite(const std::initializer_list<double> values)
{
    return std::ranges::all_of(values, [](const double v) { return std::isfinite(v); });
}

[[maybe_unused]] [[nodiscard]] inline bool FiniteFrame(const Point3& p, const Vec3& x,
                                                       const Vec3& y)
{
    return AllFinite({p.X, p.Y, p.Z, x.X, x.Y, x.Z, y.X, y.Y, y.Z});
}

[[maybe_unused]] [[nodiscard]] inline Result<void>
CheckSpline(const BSplineCurve& curve, const Interval& domain, const std::string_view where)
{
    STEREON_TRY_VOID(CheckKnotVector(curve.Degree, curve.Knots, curve.Poles.size(), where));
    STEREON_TRY_VOID(CheckWeights(curve.Weights, curve.Poles.size(), where));
    for (const Point3& pole : curve.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return CheckPeriodicDomain(curve, domain, where);
}

[[maybe_unused]] [[nodiscard]] inline bool Valid(const Line& l)
{
    return FiniteFrame(l.Origin, l.Direction, Vec3{}) && l.Direction != Vec3{};
}

[[maybe_unused]] [[nodiscard]] inline bool Valid(const Circle& o)
{
    return FiniteFrame(o.Center, o.XAxis, o.YAxis) && AllFinite({o.Radius}) && o.Radius > 0.0;
}

[[maybe_unused]] [[nodiscard]] inline bool Valid(const Ellipse& e)
{
    return FiniteFrame(e.Center, e.XAxis, e.YAxis) && AllFinite({e.MajorRadius, e.MinorRadius}) &&
           e.MinorRadius > 0.0;
}

[[maybe_unused]] [[nodiscard]] inline bool Valid(const Hyperbola& h)
{
    return FiniteFrame(h.Center, h.XAxis, h.YAxis) && AllFinite({h.MajorRadius, h.MinorRadius}) &&
           h.MinorRadius > 0.0;
}

[[maybe_unused]] [[nodiscard]] inline bool Valid(const Parabola& p)
{
    return FiniteFrame(p.Vertex, p.XAxis, p.YAxis) && AllFinite({p.FocalLength}) &&
           p.FocalLength > 0.0;
}

/// Whether the curve, when it is a T, has valid data.
template <typename T>
[[maybe_unused]] [[nodiscard]] inline bool ValidIf(const Curve& c)
{
    const T* curve = c.Target<T>();
    return curve == nullptr || Valid(*curve);
}

[[maybe_unused]] [[nodiscard]] inline bool ValidCurve(const Curve& c)
{
    return ValidIf<Line>(c) && ValidIf<Circle>(c) && ValidIf<Ellipse>(c) && ValidIf<Hyperbola>(c) &&
           ValidIf<Parabola>(c);
}

/// The curve's enclosed pieces over a bounded domain.
[[maybe_unused]] [[nodiscard]] inline Result<PieceCurve> PiecesOf(const Curve& c,
                                                                  const Interval& domain)
{
    if (const BSplineCurve* s = c.Target<BSplineCurve>()) {
        STEREON_TRY(std::vector<Segment> segments, SegmentsOf(*s));
        const bool closed = Period(*s) > 0.0 && Contains(domain, Domain(*s));
        return PieceCurve{Reaching(std::move(segments), domain), domain, ParameterMap::Same, closed,
                          "BSplineCurve"};
    }
    if (const Line* l = c.Target<Line>()) {
        return Pieces(*l, domain);
    }
    if (const Circle* o = c.Target<Circle>()) {
        const ConicData data{o->Center, Times(o->Radius, o->XAxis), Times(o->Radius, o->YAxis)};
        return ConicPieces(*o, data, domain, "Circle");
    }
    if (const Ellipse* e = c.Target<Ellipse>()) {
        const ConicData data{e->Center, Times(e->MajorRadius, e->XAxis),
                             Times(e->MinorRadius, e->YAxis)};
        return ConicPieces(*e, data, domain, "Ellipse");
    }
    if (const Hyperbola* h = c.Target<Hyperbola>()) {
        return Pieces(*h, domain);
    }
    if (const Parabola* p = c.Target<Parabola>()) {
        return Pieces(*p, domain);
    }
    return Fail(ErrorCode::NotImplemented, "Intersect: a curve type with no enclosed pieces");
}

// --- Unbounded domains ----------------------------------------------------

/// A ball holding a bounded curve, with room to spare.
struct Ball {
    Point3 Center;
    double Radius = 0.0;
};

/// The ball of the pieces' poles (they hold the curve, the weights being
/// positive), widened past the resolution so that a curve clipped to it
/// ends more than the resolution from the bounded curve.
[[maybe_unused]] [[nodiscard]] inline Ball BallOf(const std::vector<Segment>& segments,
                                                  const double resolution)
{
    Box3 box = Box3::Empty();
    for (const Segment& s : segments) {
        for (const HPoint& h : s.Poles) {
            const Point3 p = Nominal3(h);
            box.X = Hull(box.X, Interval{p.X});
            box.Y = Hull(box.Y, Interval{p.Y});
            box.Z = Hull(box.Z, Interval{p.Z});
        }
    }
    const Point3 center{Midpoint(box.X), Midpoint(box.Y), Midpoint(box.Z)};
    double radius = 0.0;
    for (const Segment& s : segments) {
        for (const HPoint& h : s.Poles) {
            radius = std::max(radius, Distance(Nominal3(h), center));
        }
    }
    const double scale = std::max({1.0, radius, Length(center - Point3{})});
    return {center, radius + 2.0 * resolution + 0x1p-30 * scale};
}

/// The parameters of an unbounded curve that can lie in the ball.
[[maybe_unused]] [[nodiscard]] inline Result<Interval> Reach(const Curve& c, const Ball& ball)
{
    if (const Line* l = c.Target<Line>()) {
        const double along =
            Dot(ball.Center - l->Origin, l->Direction) / Dot(l->Direction, l->Direction);
        const double half = ball.Radius / Length(l->Direction);
        const double slack = 0x1p-30 * (std::abs(along) + half);
        return Interval{along - half - slack, along + half + slack};
    }
    if (const Parabola* p = c.Target<Parabola>()) {
        // Off the axis by |t|, the axes being orthonormal.
        const double reach = (Distance(p->Vertex, ball.Center) + ball.Radius) * (1.0 + 0x1p-30);
        return Interval{-reach, reach};
    }
    if (const Hyperbola* h = c.Target<Hyperbola>()) {
        // Off the major axis by b |sinh t|, the axes being orthonormal.
        const double reach =
            std::asinh((Distance(h->Center, ball.Center) + ball.Radius) / h->MinorRadius) *
            (1.0 + 0x1p-30);
        return Interval{-reach, reach};
    }
    return Fail(ErrorCode::InvalidInput,
                "Intersect: a curve of bounded parameter given an unbounded domain");
}

/// The requested domain on the curve, checked: nothing past the curve's
/// own domain, and at most one turn of a circle or an ellipse.
[[maybe_unused]] [[nodiscard]] inline Result<Interval> RequestedOf(const Curve& c,
                                                                   const Interval& requested)
{
    if (!ValidCurve(c)) {
        return Fail(ErrorCode::InvalidInput, "Intersect: invalid curve data");
    }
    if (const BSplineCurve* s = c.Target<BSplineCurve>()) {
        STEREON_TRY_VOID(CheckSpline(*s, requested, "BSplineCurve x curve"));
        return Intersection(requested, Domain(*s));
    }
    if (c.Target<Circle>() != nullptr || c.Target<Ellipse>() != nullptr) {
        if (!requested.IsBounded()) {
            return Fail(ErrorCode::InvalidInput, "Intersect: a conic's domain must be bounded");
        }
        if (Width(requested) > 2.0 * std::numbers::pi * (1.0 + 0x1p-40)) {
            return Fail(ErrorCode::NotSupported,
                        "Intersect: a conic's domain of more than one turn");
        }
        return requested;
    }
    return Intersection(requested, Domain(c));
}

} // namespace
} // namespace Stereon::Detail
