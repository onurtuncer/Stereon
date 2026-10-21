// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

export module Stereon.Mesh:Sampling;

import Stereon.Surface;

/// Display sampling (ADR-0035, section 2): curves to polylines and surfaces
/// to iso-lines, for the debug viewer and any other frontend. A display
/// tolerance steers it, a purpose of its own (ADR-0028, ADR-0034 section 5)
/// that never reaches modelling code: nothing here decides topology or
/// geometry, and no modelling tolerance is read.
///
/// A polyline is refined until, on every segment, seven evenly spaced
/// interior points of the curve lie within the chord tolerance of the
/// segment and the turn between its two halves is within the angle
/// tolerance. Between those points the deviation is not bounded: it is a
/// display, not a certificate. Every distinct knot of a B-spline (and of a
/// trimmed B-spline's basis) inside the range is a vertex of the polyline,
/// so a kink at a knot is drawn as one.
///
/// Unbounded geometry (a line, a plane, a cylinder's height) is sampled over
/// the range the caller gives, which the domain clips; without one it is
/// refused, since the scene, not the geometry, knows how far to draw.

export namespace Stereon {

/// How finely to sample for display, in the model's units.
struct DisplayTolerance {
    /// The largest distance between the curve and a segment, millimetres.
    double Chord = 0.01;
    /// The largest turn between consecutive half-segments, radians.
    double Angle = 0.1;
    /// The most points one polyline may have; sampling that needs more is
    /// refused rather than truncated.
    std::size_t MaxPoints = 65536;
};

/// A curve sampled for display: its points and the parameter of each.
template <std::size_t D>
struct BasicPolyline {
    std::vector<double> Parameters;
    std::vector<BasicPointN<double, D>> Points;
};

using Polyline = BasicPolyline<3>;
using Polyline2d = BasicPolyline<2>;

/// Which parameter an iso-line runs along.
enum class IsoDirection : std::uint8_t {
    /// u varies; v is fixed.
    AlongU,
    /// v varies; u is fixed.
    AlongV,
};

/// One iso-line of a sampled surface: its direction, the fixed parameter's
/// value, and the polyline, whose parameters are the varying one's.
struct IsoLine {
    IsoDirection Direction = IsoDirection::AlongU;
    double Fixed = 0.0;
    Polyline Line;
};

/// The parameters a sampled curve must pass through inside `range`: its
/// ends and every distinct knot of a B-spline or of a trimmed curve's
/// B-spline basis, ascending.
[[nodiscard]] std::vector<double> Breakpoints(const Curve& curve, Interval range);
[[nodiscard]] std::vector<double> Breakpoints(const Curve2d& curve, Interval range);

/// The curve over `range` (the domain if none is given; the range is
/// clipped to a non-periodic domain). InvalidInput for a tolerance that is
/// not finite and positive, fewer than two points allowed, or a range that
/// is not finite and non-empty; NotSupported when the tolerance needs more
/// than MaxPoints points or cannot be met before subdivision stops (depth
/// limit, parameter precision, or non-finite evaluations).
[[nodiscard]] Result<Polyline> SampleCurve(const Curve& curve, const DisplayTolerance& tolerance,
                                           std::optional<Interval> range = std::nullopt);
[[nodiscard]] Result<Polyline2d> SampleCurve(const Curve2d& curve,
                                             const DisplayTolerance& tolerance,
                                             std::optional<Interval> range = std::nullopt);

/// Iso-lines of the surface over `range` (the domain if none is given,
/// clipped as for curves): `linesPerDirection` evenly spaced lines each way,
/// ends included, and one more at every distinct interior knot of a
/// B-spline surface or a trimmed one's B-spline basis. The AlongU lines come
/// first, by ascending v. Errors as SampleCurve, and InvalidInput for fewer
/// than two lines per direction.
[[nodiscard]] Result<std::vector<IsoLine>> SampleSurface(const Surface& surface,
                                                         const DisplayTolerance& tolerance,
                                                         std::optional<Box2> range = std::nullopt,
                                                         std::size_t linesPerDirection = 9);

} // namespace Stereon
