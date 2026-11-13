// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Conics as enclosed rational quadratic Bezier arcs, shared by the
// curve-surface and curve-curve searches. Every quantity is enclosed from
// the doubles given, so an arc encloses the conic as given, not a rounded
// NURBS form of it.
//
// Included in the purview of the handlers' implementation units, after
// Bernstein.hpp, with <algorithm> and <cmath> in their global module
// fragment. Everything is in an anonymous namespace, so each unit gets its
// own copy.

#pragma once

namespace Stereon::Detail {
namespace {

/// A vector enclosed.
struct IVec3 {
    Interval X;
    Interval Y;
    Interval Z;
};

/// A circle or an ellipse as C + A cos t + B sin t, the vectors enclosed:
/// the radius times the axis, as given.
struct ConicData {
    Point3 Center;
    IVec3 A;
    IVec3 B;
};

[[maybe_unused]] [[nodiscard]] inline IVec3 Times(const double r, const Vec3& v)
{
    const Interval s{r};
    return {s * Interval{v.X}, s * Interval{v.Y}, s * Interval{v.Z}};
}

/// The homogeneous pole w C + c A + s B of weight w, enclosed.
[[maybe_unused]] [[nodiscard]] inline HPoint ConicPole(const ConicData& conic, const Interval& w,
                                                       const Interval& c, const Interval& s)
{
    return {w * Interval{conic.Center.X} + c * conic.A.X + s * conic.B.X,
            w * Interval{conic.Center.Y} + c * conic.A.Y + s * conic.B.Y,
            w * Interval{conic.Center.Z} + c * conic.A.Z + s * conic.B.Z, w};
}

/// The arc [t0, t1], at most a quarter turn, as a rational quadratic
/// Bezier segment on s in [0, 1]: the end poles at t0 and t1, weight 1; the
/// middle one where the end tangents meet, C + (A cos m + B sin m) / cos h
/// with h the half-angle and m the middle, weight cos h, stored as its
/// homogeneous (cos h C + A cos m + B sin m, cos h). Every quantity is
/// enclosed from the doubles t0 and t1, so the segment encloses the arc of
/// the conic as given, not a rounded NURBS form of it.
[[maybe_unused]] [[nodiscard]] inline Segment ArcSegment(const ConicData& conic, const double t0,
                                                         const double t1)
{
    const IntervalScalar a{Interval{t0}};
    const IntervalScalar b{Interval{t1}};
    const IntervalScalar half{(Interval{t1} - Interval{t0}) / Interval{2.0}};
    const IntervalScalar middle{Interval{t0} + half.Value};
    // sin(0) is exactly zero; a libm rounding enclosure would lose an
    // exact endpoint contact before the polynomial search sees it.
    const Interval one{1.0};
    return {0.0,
            1.0,
            {ConicPole(conic, one, cos(a).Value, t0 == 0.0 ? Interval{0.0} : sin(a).Value),
             ConicPole(conic, cos(half).Value, cos(middle).Value, sin(middle).Value),
             ConicPole(conic, one, cos(b).Value, t1 == 0.0 ? Interval{0.0} : sin(b).Value)}};
}

/// The conic's angle at the arc's Bezier parameter s: stereographically,
/// tan((t - m) / 2) = tan(h / 2) (2 s - 1), exact at the arc's ends.
[[maybe_unused]] [[nodiscard]] inline double AngleAt(const double t0, const double t1,
                                                     const double s)
{
    if (s <= 0.0) {
        return t0;
    }
    if (s >= 1.0) {
        return t1;
    }
    const double half = 0.5 * (t1 - t0);
    return std::clamp(t0 + half + 2.0 * std::atan(std::tan(0.5 * half) * (2.0 * s - 1.0)), t0, t1);
}

/// The hyperbola's arc [t0, t1] as a rational quadratic Bezier segment on
/// s in [0, 1]: as ArcSegment with cosh and sinh for cos and sin, the middle
/// weight cosh h for the half-width h.
[[maybe_unused]] [[nodiscard]] inline Segment HyperbolicSegment(const ConicData& conic,
                                                                const double t0, const double t1)
{
    const IntervalScalar a{Interval{t0}};
    const IntervalScalar b{Interval{t1}};
    const IntervalScalar half{(Interval{t1} - Interval{t0}) / Interval{2.0}};
    const IntervalScalar middle{Interval{t0} + half.Value};
    const Interval one{1.0};
    return {0.0,
            1.0,
            {ConicPole(conic, one, cosh(a).Value, sinh(a).Value),
             ConicPole(conic, cosh(half).Value, cosh(middle).Value, sinh(middle).Value),
             ConicPole(conic, one, cosh(b).Value, sinh(b).Value)}};
}

/// The hyperbola's parameter at the arc's Bezier parameter s:
/// tanh((t - m) / 2) = tanh(h / 2) (2 s - 1), exact at the arc's ends.
[[maybe_unused]] [[nodiscard]] inline double HyperbolicAt(const double t0, const double t1,
                                                          const double s)
{
    if (s <= 0.0) {
        return t0;
    }
    if (s >= 1.0) {
        return t1;
    }
    const double half = 0.5 * (t1 - t0);
    return std::clamp(t0 + half + 2.0 * std::atanh(std::tanh(0.5 * half) * (2.0 * s - 1.0)), t0,
                      t1);
}

} // namespace
} // namespace Stereon::Detail
