// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "BoxClip.hpp"
#include "Planar.hpp"

/// A plane against an ellipsoid (ADR-0048 §3.3), in closed form. In the
/// ellipsoid's frame (e1, e2, e3), e3 its axis as the evaluator rounds it,
/// the map s = (x / A, y / B, z / C) takes the ellipsoid onto the unit
/// sphere and the plane n . (x - O) = 0 onto the plane m . s = d, with
/// m = (A n . e1, B n . e2, C n . e3) and d = n . (O - Ce). The plane meets
/// the unit sphere in a circle, touches it at a point or misses it as d²
/// is below, at or above |m|²: a polynomial in the input doubles, decided
/// exactly (SignOfFormula). The circle maps back onto an ellipse, the
/// affine image of a circle, whose conjugate semi-diameters are turned to
/// its principal axes; that ellipse is the branch's 3D curve, exact, and
/// its pcurve on the plane is the exact Ellipse2d.
///
/// Its pcurve on the ellipsoid is not a conic in (u, v): it is the cubic
/// interpolant through the section's longitudes and latitudes, sampled more
/// densely until CertifyDeviation proves it within the Context's
/// resolution, as the trace fits its pcurves (ADR-0009, D3). The section is
/// cut where it crosses the chart's seam, so each piece's pcurve stays in
/// the base turn [0, 2 pi]. A plane through a pole of the chart, where the
/// longitude has no value, is decided exactly too and fails closed.

namespace Stereon::Detail {

namespace {

constexpr double Pi = std::numbers::pi;

/// Samples of a pcurve piece at first and at most: a work budget, not a
/// geometric tolerance (ADR-0010).
constexpr std::size_t FirstSamples = 32;
constexpr std::size_t SampleBudget = std::size_t{1} << 14;

/// The two surfaces, erased and as themselves.
struct Sides {
    const Surface& A;
    const Surface& B;
    const Plane& P;
    const Ellipsoid& E;
};

/// The plane's normal X x Y over a field.
template <typename F>
[[nodiscard]] Triple<F> NormalOf(F f, const Plane& p)
{
    return Cross3(f, Of3(f, p.XAxis), Of3(f, p.YAxis));
}

/// The ellipsoid's centre's height over the plane, times |n|: n . (Ce - O).
template <typename F>
[[nodiscard]] typename F::Scalar HeightOf(F f, const Sides& s)
{
    return Dot3(f, NormalOf(f, s.P), Diff3(f, s.E.Center, s.P.Origin));
}

/// The normal's components along the ellipsoid's frame, the axis as
/// X x Y: its handedness flips the sign of the third, which no decision
/// below depends on.
template <typename F>
[[nodiscard]] std::array<typename F::Scalar, 3> AlongFrame(F f, const Sides& s)
{
    const Triple<F> n = NormalOf(f, s.P);
    const Triple<F> x = Of3(f, s.E.XAxis);
    const Triple<F> y = Of3(f, s.E.YAxis);
    return {Dot3(f, n, x), Dot3(f, n, y), Dot3(f, n, Cross3(f, x, y))};
}

/// |m|² - d²: positive for a section, zero for a touch, negative apart.
template <typename F>
[[nodiscard]] typename F::Scalar Room(F f, const Sides& s)
{
    using S = typename F::Scalar;
    const std::array<S, 3> k = AlongFrame(f, s);
    const std::array<S, 3> semi{F::Of(s.E.A), F::Of(s.E.B), F::Of(s.E.C)};
    S reach = F::Mul(F::Mul(semi[0], k[0]), F::Mul(semi[0], k[0]));
    for (std::size_t i = 1; i < 3; ++i) {
        reach = F::Add(reach, F::Mul(F::Mul(semi[i], k[i]), F::Mul(semi[i], k[i])));
    }
    const S h = HeightOf(f, s);
    return F::Sub(reach, F::Mul(h, h));
}

/// Whether the plane passes through a pole, Ce +- C e3: n . (Ce - O) +-
/// C n . e3 = 0, decided exactly.
[[nodiscard]] bool ThroughPole(const Sides& s)
{
    for (const double side : {1.0, -1.0}) {
        const Sign sign = SignOfFormula([&](auto f) {
            using F = decltype(f);
            const auto k = AlongFrame(f, s);
            return F::Add(HeightOf(f, s), F::Mul(F::Of(side * s.E.C), k[2]));
        });
        if (sign == Sign::Zero) {
            return true;
        }
    }
    return false;
}

/// The ellipsoid's frame in doubles, the axis as the evaluator takes it.
struct Frame3 {
    Vec3 X;
    Vec3 Y;
    Vec3 Z;
    Vec3 Semi;
};

[[nodiscard]] Frame3 FrameOf(const Ellipsoid& e)
{
    return {e.XAxis, e.YAxis, Axis(e), Vec3{e.A, e.B, e.C}};
}

[[nodiscard]] Vec3 Local(const Frame3& f, const Vec3& w)
{
    return {Dot(w, f.X), Dot(w, f.Y), Dot(w, f.Z)};
}

[[nodiscard]] Vec3 World(const Frame3& f, const Vec3& l)
{
    return l.X * f.X + l.Y * f.Y + l.Z * f.Z;
}

[[nodiscard]] Vec3 Scaled(const Vec3& a, const Vec3& b)
{
    return {a.X * b.X, a.Y * b.Y, a.Z * b.Z};
}

/// A unit vector across the unit vector m.
[[nodiscard]] Vec3 Across(const Vec3& m)
{
    const Vec3 helper = std::abs(m.X) <= std::abs(m.Y) && std::abs(m.X) <= std::abs(m.Z)
                            ? Vec3{1.0, 0.0, 0.0}
                        : std::abs(m.Y) <= std::abs(m.Z) ? Vec3{0.0, 1.0, 0.0}
                                                         : Vec3{0.0, 0.0, 1.0};
    return Normalized(Cross(m, helper));
}

/// The section's ellipse: the image of the circle in the sphere's
/// coordinates, its conjugate semi-diameters turned to principal axes.
[[nodiscard]] Ellipse SectionOf(const Sides& s)
{
    const Frame3 f = FrameOf(s.E);
    const Vec3 n = Cross(s.P.XAxis, s.P.YAxis);
    const Vec3 m = Scaled(Local(f, n), f.Semi);
    const double mm = Dot(m, m);
    const double d = Dot(n, s.P.Origin - s.E.Center);
    // |m|² - d² to a few roundings: it cancels near a touch.
    const double room = std::max(0.0, ValueOfFormula([&](auto g) { return Room(g, s); }));
    const double rho = std::sqrt(room / mm);
    const Vec3 unit = m / std::sqrt(mm);
    const Vec3 e1 = Across(unit);
    const Vec3 e2 = Cross(unit, e1);
    const Vec3 p = rho * Scaled(e1, f.Semi);
    const Vec3 q = rho * Scaled(e2, f.Semi);
    const double t0 = 0.5 * std::atan2(2.0 * Dot(p, q), Dot(p, p) - Dot(q, q));
    Vec3 major = std::cos(t0) * p + std::sin(t0) * q;
    Vec3 minor = -std::sin(t0) * p + std::cos(t0) * q;
    if (Length(minor) > Length(major)) {
        std::swap(major, minor);
        minor = -minor;
    }
    const Point3 center = s.E.Center + World(f, Scaled((d / mm) * m, f.Semi));
    const Vec3 x = Normalized(World(f, major));
    const Vec3 along = World(f, minor);
    const Vec3 y = Normalized(along - Dot(along, x) * x);
    return Ellipse{center, Length(major), Length(minor), x, y};
}

[[nodiscard]] Vec2 InPlane(const Plane& plane, const Vec3& v)
{
    return {Dot(v, plane.XAxis), Dot(v, plane.YAxis)};
}

/// The section's longitude and latitude at the ellipse's angle t, with
/// their derivatives: s = (x / A, y / B, z / C) lies on the unit sphere,
/// u = atan2(s2, s1), v = atan2(s3, |(s1, s2)|).
struct UvJet {
    Point2 At;
    Vec2 Rate;
};

[[nodiscard]] UvJet UvAt(const Frame3& f, const Ellipsoid& e, const Ellipse& section,
                         const double t)
{
    const Vec3 inverse{1.0 / f.Semi.X, 1.0 / f.Semi.Y, 1.0 / f.Semi.Z};
    const Vec3 s = Scaled(Local(f, Eval(section, t) - e.Center), inverse);
    const Vec3 ds = Scaled(Local(f, Deriv(section, t, 1)), inverse);
    const double radial2 = s.X * s.X + s.Y * s.Y;
    const double radial = std::sqrt(radial2);
    const double length2 = radial2 + s.Z * s.Z;
    const double du = (s.X * ds.Y - s.Y * ds.X) / radial2;
    const double dv = (ds.Z * radial2 - s.Z * (s.X * ds.X + s.Y * ds.Y)) / (length2 * radial);
    return {Point2{WrappedAngle(s.Y, s.X), std::atan2(s.Z, radial)}, Vec2{du, dv}};
}

/// Where the ellipse crosses the chart's seam, the half-plane of local
/// y = 0 and x > 0, in [0, 2 pi), in order. A touch of the seam is not a
/// crossing: the pcurve there reaches u = 0 or 2 pi and turns back.
[[nodiscard]] std::vector<double> SeamCrossings(const Frame3& f, const Ellipsoid& e,
                                                const Ellipse& section)
{
    const double c = Dot(section.Center - e.Center, f.Y);
    const double a = section.MajorRadius * Dot(section.XAxis, f.Y);
    const double b = section.MinorRadius * Dot(section.YAxis, f.Y);
    const double r = std::hypot(a, b);
    std::vector<double> crossings;
    if (!(r > std::abs(c))) {
        return crossings;
    }
    const double phi = std::atan2(b, a);
    const double alpha = std::acos(-c / r);
    for (const double t : {phi + alpha, phi - alpha}) {
        const double placed = t - CirclePeriod * std::floor(t / CirclePeriod);
        if (Dot(Eval(section, placed) - e.Center, f.X) > 0.0) {
            crossings.push_back(placed);
        }
    }
    std::ranges::sort(crossings);
    return crossings;
}

/// The pieces of the ellipse's angle between seam crossings: the whole
/// turn when it crosses none, from a crossing round to it when one.
[[nodiscard]] std::vector<Interval> PiecesOf(const std::vector<double>& crossings)
{
    if (crossings.empty()) {
        return {Interval{0.0, CirclePeriod}};
    }
    std::vector<Interval> pieces;
    for (std::size_t i = 0; i < crossings.size(); ++i) {
        const double start = crossings[i];
        const double end =
            i + 1 < crossings.size() ? crossings[i + 1] : crossings.front() + CirclePeriod;
        pieces.push_back(Interval{start, end});
    }
    return pieces;
}

/// The longitudes along a piece, unwrapped from its middle, which lies off
/// the seam, and clamped into the base turn: at a piece's ends on the seam
/// they reach 0 or 2 pi to rounding.
[[nodiscard]] std::vector<UvJet> SamplesOf(const Frame3& f, const Ellipsoid& e,
                                           const Ellipse& section, const Interval& piece,
                                           const std::size_t n)
{
    std::vector<UvJet> samples(n + 1);
    const double width = piece.GetUpper() - piece.GetLower();
    for (std::size_t k = 0; k <= n; ++k) {
        const double t =
            k == n ? piece.GetUpper()
                   : piece.GetLower() + width * (static_cast<double>(k) / static_cast<double>(n));
        samples[k] = UvAt(f, e, section, t);
    }
    const auto unwrap = [&](const std::size_t from, const std::size_t to) {
        const double near = samples[from].At.X;
        double& u = samples[to].At.X;
        u -= CirclePeriod * std::round((u - near) / CirclePeriod);
    };
    const std::size_t middle = n / 2;
    for (std::size_t k = middle + 1; k <= n; ++k) {
        unwrap(k - 1, k);
    }
    for (std::size_t k = middle; k-- > 0;) {
        unwrap(k + 1, k);
    }
    for (UvJet& s : samples) {
        s.At.X = std::clamp(s.At.X, 0.0, CirclePeriod);
    }
    return samples;
}

/// The cubic interpolant through n + 1 samples of the piece, on the
/// ellipse's own angle.
[[nodiscard]] Result<Curve2d> InterpolantOf(const std::vector<UvJet>& samples,
                                            const Interval& piece)
{
    std::vector<Point2> points;
    points.reserve(samples.size());
    for (const UvJet& s : samples) {
        points.push_back(s.At);
    }
    // The fit runs on [0, 1], `width` times faster than the angle.
    const double width = piece.GetUpper() - piece.GetLower();
    const EndDerivatives2d ends{width * samples.front().Rate, width * samples.back().Rate};
    STEREON_TRY(const CurveFit2d fit, Interpolate(points, 3, ends, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d pcurve, Reparametrized(fit.Curve, piece));
    return Curve2d(std::move(pcurve));
}

/// The pcurve on the ellipsoid over a piece: samples doubled until the
/// interpolant is proved within half the Context's resolution of the
/// ellipsoid along the section. Past the budget the section passes within
/// rounding of a pole, where its longitude turns by a half turn in no
/// distance: IntersectionFailed, as the trace answers a branch through a
/// pole, and as a plane exactly through one is answered here. A plane
/// whose frame is rounded can miss the pole it was made through by a unit
/// in the last place, which no exact test on the rounded plane can tell
/// from a plane through it.
[[nodiscard]] Result<Curve2d> FittedPcurve(const Sides& s, const Ellipse& section,
                                           const Interval& piece, const Context& ctx)
{
    const Frame3 f = FrameOf(s.E);
    const double target = 0.5 * ctx.GetResolution();
    for (std::size_t n = FirstSamples; n <= SampleBudget; n *= 2) {
        STEREON_TRY(Curve2d pcurve, InterpolantOf(SamplesOf(f, s.E, section, piece, n), piece));
        const Result<DeviationCertificate> proved =
            CertifyDeviation(Curve(section), s.B, pcurve, piece, target);
        if (proved) {
            return pcurve;
        }
        if (proved.error().GetCode() != ErrorCode::CertificationFailed) {
            return std::unexpected(proved.error());
        }
    }
    return Fail(ErrorCode::IntersectionFailed,
                "PlaneEllipsoid: the section passes within rounding of a pole of the ellipsoid's "
                "chart, and its pcurve could not be certified within the sample budget "
                "(ADR-0009)");
}

/// Whether the ellipsoid's requested box holds its whole chart.
[[nodiscard]] bool WholeChart(const Box2& box)
{
    return Contains(box.X, Interval{0.0, CirclePeriod}) &&
           Contains(box.Y, Interval{-Pi / 2.0, Pi / 2.0});
}

/// The stretches of a piece inside both requested boxes.
[[nodiscard]] Result<std::vector<Interval>>
StretchesOf(const Curve2d& onPlane, const Box2& planeBox, const Curve2d& onEllipsoid,
            const Interval& piece, const Box2& boxA, const Box2& boxB, const Context& ctx)
{
    std::vector<Interval> stretches{piece};
    if (!Contains(boxA, planeBox)) {
        STEREON_TRY(stretches, ClipPcurveToBox(onPlane, piece, boxA, ctx, "PlaneEllipsoid"));
    }
    if (WholeChart(boxB)) {
        return stretches;
    }
    std::vector<Interval> kept;
    for (const Interval& stretch : stretches) {
        STEREON_TRY(const std::vector<Interval> inside,
                    ClipPcurveToBox(onEllipsoid, stretch, boxB, ctx, "PlaneEllipsoid"));
        kept.insert(kept.end(), inside.begin(), inside.end());
    }
    return kept;
}

/// The point where the plane touches the ellipsoid: the foot of the
/// centre's perpendicular in the sphere's coordinates, mapped back.
[[nodiscard]] SurfaceSurfaceIntersection Touch(const Sides& s, const Box2& boxA, const Box2& boxB)
{
    const Frame3 f = FrameOf(s.E);
    const Vec3 n = Cross(s.P.XAxis, s.P.YAxis);
    const Vec3 m = Scaled(Local(f, n), f.Semi);
    const double d = Dot(n, s.P.Origin - s.E.Center);
    const Point3 at = s.E.Center + World(f, Scaled((d / Dot(m, m)) * m, f.Semi));
    const Point2 onA = Project(s.P, at);
    const Point2 onB = Project(s.E, at);
    SurfaceSurfaceIntersection result;
    const std::optional<double> u = InPeriodicDomain(onB.X, boxB.X);
    if (Contains(boxA.X, onA.X) && Contains(boxA.Y, onA.Y) && u && Contains(boxB.Y, onB.Y)) {
        result.Points.push_back({at, onA, Point2{*u, onB.Y}});
    }
    return result;
}

/// The section: one branch per piece between seam crossings and stretch
/// inside both boxes.
[[nodiscard]] Result<SurfaceSurfaceIntersection> Section(const Sides& s, const Box2& boxA,
                                                         const Box2& boxB, const Context& ctx)
{
    if (!WholeChart(boxB) && !Contains(Interval{0.0, CirclePeriod}, boxB.X)) {
        return Fail(ErrorCode::NotImplemented,
                    "PlaneEllipsoid: a parameter box on the ellipsoid outside its base turn "
                    "[0, 2 pi] in u");
    }
    const Ellipse section = SectionOf(s);
    const Ellipse2d onPlane{Project(s.P, section.Center), section.MajorRadius, section.MinorRadius,
                            InPlane(s.P, section.XAxis), InPlane(s.P, section.YAxis)};
    const Box2 planeBox{Around(onPlane.Center.X, section.MajorRadius),
                        Around(onPlane.Center.Y, section.MajorRadius)};
    SurfaceSurfaceIntersection result;
    if (!Intersects(planeBox, boxA)) {
        return result;
    }
    const Frame3 f = FrameOf(s.E);
    for (const Interval& piece : PiecesOf(SeamCrossings(f, s.E, section))) {
        STEREON_TRY(const Curve2d onEllipsoid, FittedPcurve(s, section, piece, ctx));
        STEREON_TRY(const std::vector<Interval> stretches,
                    StretchesOf(Curve2d(onPlane), planeBox, onEllipsoid, piece, boxA, boxB, ctx));
        for (const Interval& stretch : stretches) {
            STEREON_TRY_VOID(Keep(result,
                                  Crossing(Curve(section), Curve2d(onPlane), onEllipsoid, stretch),
                                  s.A, s.B, ctx));
        }
    }
    return result;
}

} // namespace

Result<SurfaceSurfaceIntersection> PlaneEllipsoid(const Surface& a, const Box2& boxA,
                                                  const Surface& b, const Box2& boxB,
                                                  const double /*tolerance*/, const Context& ctx)
{
    const Plane* plane = a.Target<Plane>();
    const Ellipsoid* ellipsoid = b.Target<Ellipsoid>();
    if (plane == nullptr || ellipsoid == nullptr) {
        return Fail(ErrorCode::InvalidInput, "PlaneEllipsoid: dispatched for other surface types");
    }
    const Sides s{a, b, *plane, *ellipsoid};
    const Sign room = SignOfFormula([&](auto f) { return Room(f, s); });
    if (room == Sign::Negative) {
        return SurfaceSurfaceIntersection{}; // Proved separate.
    }
    if (room == Sign::Zero) {
        return Touch(s, boxA, boxB);
    }
    if (ThroughPole(s)) {
        return Fail(ErrorCode::IntersectionFailed,
                    "PlaneEllipsoid: the section passes through a pole of the ellipsoid's chart, "
                    "where its pcurve has no longitude (ADR-0009)");
    }
    return Section(s, boxA, boxB, ctx);
}

Result<SurfaceSurfaceIntersection> EllipsoidPlane(const Surface& a, const Box2& boxA,
                                                  const Surface& b, const Box2& boxB,
                                                  const double tolerance, const Context& ctx)
{
    STEREON_TRY(SurfaceSurfaceIntersection result,
                PlaneEllipsoid(b, boxB, a, boxA, tolerance, ctx));
    return Swapped(std::move(result));
}

} // namespace Stereon::Detail
