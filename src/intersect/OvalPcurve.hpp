// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The pcurve of a planar conic lying on a sphere or an ellipsoid, in the
// surface's longitude and latitude chart: not a conic in (u, v), but the
// cubic interpolant through the conic's longitudes and latitudes, sampled
// more densely until CertifyDeviation proves it within half the Context's
// resolution of the surface (ADR-0009, D3), each piece between the
// crossings of the chart's seam. A conic through a pole of the chart,
// where the longitude has no rate, fails closed.
//
// Included in the purview of the handlers' implementation units, after
// `import :Dispatch;`, with <algorithm>, <cmath>, <cstddef>, <utility> and
// <vector> in their global module fragment. Everything is in an anonymous
// namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

/// A planar conic as C + P cos t + Q sin t: a circle or an ellipse, P and Q
/// conjugate semi-diameters.
struct Oval {
    Point3 Center;
    Vec3 P;
    Vec3 Q;
};

/// A circle as an oval: its radius along each axis.
[[nodiscard]] inline Oval OvalOf(const Circle& c)
{
    return {c.Center, c.Radius * c.XAxis, c.Radius * c.YAxis};
}

/// An ellipse as an oval: its semi-axes along its axes.
[[nodiscard]] inline Oval OvalOf(const Ellipse& e)
{
    return {e.Center, e.MajorRadius * e.XAxis, e.MinorRadius * e.YAxis};
}

/// The chart's frame: the surface's X and Y axes, its axis as the evaluator
/// takes it, its centre and its semi-axes (the radius three times for a
/// sphere).
struct OvalChart {
    Point3 Center;
    Vec3 X;
    Vec3 Y;
    Vec3 Z;
    Vec3 Semi;
};

/// An ellipsoid's chart.
[[nodiscard]] inline OvalChart ChartOf(const Ellipsoid& e)
{
    return {e.Center, e.XAxis, e.YAxis, Axis(e), Vec3{e.A, e.B, e.C}};
}

/// A sphere's chart, its radius as each semi-axis.
[[nodiscard]] inline OvalChart ChartOf(const Sphere& s)
{
    return {s.Center, s.XAxis, s.YAxis, Axis(s), Vec3{s.Radius, s.Radius, s.Radius}};
}

/// A direction in the chart's frame, divided by the semi-axes.
[[nodiscard]] inline Vec3 OvalLocal(const OvalChart& f, const Vec3& w)
{
    return {Dot(w, f.X) / f.Semi.X, Dot(w, f.Y) / f.Semi.Y, Dot(w, f.Z) / f.Semi.Z};
}

/// The longitude and latitude at the conic's angle t, with their rates:
/// s = (x / A, y / B, z / C) lies on the unit sphere, u = atan2(s2, s1),
/// v = atan2(s3, |(s1, s2)|).
struct OvalUv {
    Point2 At;
    Vec2 Rate;
};

/// The longitude and latitude at t, and their rates.
[[nodiscard]] inline OvalUv OvalUvAt(const OvalChart& f, const Oval& o, const double t)
{
    const double c = std::cos(t);
    const double sn = std::sin(t);
    const Vec3 s = OvalLocal(f, (o.Center - f.Center) + c * o.P + sn * o.Q);
    const Vec3 ds = OvalLocal(f, c * o.Q - sn * o.P);
    const double radial2 = s.X * s.X + s.Y * s.Y;
    const double radial = std::sqrt(radial2);
    const double length2 = radial2 + s.Z * s.Z;
    const double du = (s.X * ds.Y - s.Y * ds.X) / radial2;
    const double dv = (ds.Z * radial2 - s.Z * (s.X * ds.X + s.Y * ds.Y)) / (length2 * radial);
    return {Point2{WrappedAngle(s.Y, s.X), std::atan2(s.Z, radial)}, Vec2{du, dv}};
}

/// Where the conic crosses the chart's seam, the half-plane of local y = 0
/// and x > 0, in [0, 2 pi), in order. A touch of the seam is not a
/// crossing: the pcurve there reaches u = 0 or 2 pi and turns back.
[[nodiscard]] inline std::vector<double> OvalSeamCrossings(const OvalChart& f, const Oval& o)
{
    const double c = Dot(o.Center - f.Center, f.Y);
    const double a = Dot(o.P, f.Y);
    const double b = Dot(o.Q, f.Y);
    const double r = std::hypot(a, b);
    std::vector<double> crossings;
    if (!(r > std::abs(c))) {
        return crossings;
    }
    const double phi = std::atan2(b, a);
    const double alpha = std::acos(-c / r);
    for (const double t : {phi + alpha, phi - alpha}) {
        const double placed = t - CirclePeriod * std::floor(t / CirclePeriod);
        const Vec3 at = (o.Center - f.Center) + std::cos(placed) * o.P + std::sin(placed) * o.Q;
        if (Dot(at, f.X) > 0.0) {
            crossings.push_back(placed);
        }
    }
    std::ranges::sort(crossings);
    return crossings;
}

/// The longitudes along a piece, unwrapped from its middle, which lies off
/// the seam, and clamped into the base turn: at a piece's ends on the seam
/// they reach 0 or 2 pi to rounding.
[[nodiscard]] inline std::vector<OvalUv> OvalSamples(const OvalChart& f, const Oval& o,
                                                     const Interval& piece, const std::size_t n)
{
    std::vector<OvalUv> samples(n + 1);
    const double width = piece.GetUpper() - piece.GetLower();
    for (std::size_t k = 0; k <= n; ++k) {
        const double t =
            k == n ? piece.GetUpper()
                   : piece.GetLower() + width * (static_cast<double>(k) / static_cast<double>(n));
        samples[k] = OvalUvAt(f, o, t);
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
    for (OvalUv& s : samples) {
        s.At.X = std::clamp(s.At.X, 0.0, CirclePeriod);
    }
    return samples;
}

/// The pieces of the conic's angle between seam crossings: the whole turn
/// when it crosses none, from a crossing round to the next otherwise.
[[nodiscard]] inline std::vector<Interval> OvalPieces(const std::vector<double>& crossings)
{
    if (crossings.empty()) {
        return {Interval{0.0, CirclePeriod}};
    }
    std::vector<Interval> pieces;
    for (std::size_t i = 0; i < crossings.size(); ++i) {
        const double end =
            i + 1 < crossings.size() ? crossings[i + 1] : crossings.front() + CirclePeriod;
        pieces.emplace_back(crossings[i], end);
    }
    return pieces;
}

/// The cubic interpolant through the samples of a piece, on the conic's own
/// angle; nothing at a sample on a pole to rounding, where the longitude's
/// rate is 0 / 0.
[[nodiscard]] inline Result<std::optional<Curve2d>>
OvalInterpolant(const std::vector<OvalUv>& samples, const Interval& piece)
{
    for (const OvalUv& s : samples) {
        if (!std::isfinite(s.Rate.X) || !std::isfinite(s.Rate.Y)) {
            return std::optional<Curve2d>{};
        }
    }
    std::vector<Point2> points;
    points.reserve(samples.size());
    for (const OvalUv& s : samples) {
        points.push_back(s.At);
    }
    // The fit runs on [0, 1], `width` times faster than the angle.
    const double width = piece.GetUpper() - piece.GetLower();
    const EndDerivatives2d ends{width * samples.front().Rate, width * samples.back().Rate};
    STEREON_TRY(const CurveFit2d fit, Interpolate(points, 3, ends, Parametrization::Uniform));
    STEREON_TRY(BSplineCurve2d pcurve, Reparametrized(fit.Curve, piece));
    return std::optional<Curve2d>{Curve2d(std::move(pcurve))};
}

/// The pcurve on the surface over a piece: samples doubled from 32 to 2^14
/// until the interpolant is proved within half the Context's resolution of
/// the surface along the conic. Past the budget, or at a sample on a pole,
/// the conic passes within rounding of a pole of the chart, where its
/// longitude turns by a half turn in no distance: IntersectionFailed.
[[nodiscard]] inline Result<Curve2d> FittedOvalPcurve(const OvalChart& f, const Oval& o,
                                                      const Curve& conic, const Surface& surface,
                                                      const Interval& piece, const Context& ctx,
                                                      const std::string_view handler)
{
    const double target = 0.5 * ctx.GetResolution();
    for (std::size_t n = 32; n <= (std::size_t{1} << 14); n *= 2) {
        STEREON_TRY(const std::optional<Curve2d> pcurve,
                    OvalInterpolant(OvalSamples(f, o, piece, n), piece));
        if (!pcurve) {
            break;
        }
        const Result<DeviationCertificate> proved =
            CertifyDeviation(conic, surface, *pcurve, piece, target);
        if (proved) {
            return *pcurve;
        }
        if (proved.error().GetCode() != ErrorCode::CertificationFailed) {
            return std::unexpected(proved.error());
        }
    }
    return Fail(ErrorCode::IntersectionFailed,
                std::string(handler) +
                    ": the curve passes within rounding of a pole of the surface's chart, where "
                    "its pcurve has no longitude, and could not be certified (ADR-0009)");
}

} // namespace
} // namespace Stereon::Detail
