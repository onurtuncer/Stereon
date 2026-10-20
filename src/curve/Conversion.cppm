// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <vector>

export module Stereon.Curve:Conversion;

import Stereon.Robust;
import :BSplineCurve;
import :Circle;
import :Ellipse;
import :Erased;
import :Hyperbola;
import :Line;

/// Exact conversion of the analytic curves to NURBS (ADR-0003, tier 3): every
/// built-in curve, restricted to a finite parameter range, is a B-spline or
/// a rational B-spline exactly, not by approximation (The NURBS Book,
/// chapter 7). The result has the same point set and orientation, and its
/// domain is the range itself. Lines and parabolas keep their
/// parametrisation everywhere; the conics converted to rational quadratic
/// arcs agree with the original parameter at every knot, but not between
/// knots, because no rational polynomial is proportional to the angle along
/// a circle (docs/nurbs.rst).

export namespace Stereon {

/// The curve over `range` as a B-spline. InvalidInput unless the range is
/// finite, non-empty and inside the domain; for the periodic circle and
/// ellipse, any finite range at most one period long.
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Line& line, Interval range);
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Circle& circle, Interval range);
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Ellipse& ellipse, Interval range);
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Hyperbola& hyperbola, Interval range);
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Parabola& parabola, Interval range);
/// A B-spline restricted to `range`: the piece between the range's ends.
[[nodiscard]] Result<BSplineCurve> ToBSpline(const BSplineCurve& curve, Interval range);
/// Any built-in curve through the wrapper, and a trimmed curve of one,
/// nested to any depth (ToBSpline(const TrimmedCurve&, Interval));
/// NotSupported for other types, which have no exact representation the
/// kernel knows of, and InvalidInput for an empty wrapper.
[[nodiscard]] Result<BSplineCurve> ToBSpline(const Curve& curve, Interval range);

} // namespace Stereon

export namespace Stereon::Detail {

/// True when ToBSpline keeps the curve's parametrisation everywhere, not
/// only at its knots: a line, a parabola, a B-spline, or a trimmed curve of
/// one, nested to any depth. A conic's rational arcs agree with it only at
/// their knots. False for an empty wrapper.
[[nodiscard]] bool KeepsParameter(const Curve& curve);
[[nodiscard]] bool KeepsParameter(const Curve2d& curve);

/// The unit circle arc from angle `start` to `end` (end > start, at most
/// 2 pi apart) as a rational quadratic B-spline in the plane: split into
/// segments of at most a quarter turn, each the Bezier arc with poles
/// (cos a, sin a), (cos m, sin m) / cos(d/2), (cos b, sin b) and middle
/// weight cos(d/2), for half-angle d/2 and mid-angle m (algorithm A7.1).
/// Its knots are the segment end angles, so it passes through the circle's
/// point of angle t at every knot t. An arc of exactly one turn ends at its
/// first pole, so that it closes exactly.
struct UnitArc {
    std::vector<double> Knots;
    std::vector<Point2> Poles;
    std::vector<double> Weights;
};

[[nodiscard]] UnitArc MakeUnitArc(double start, double end);

} // namespace Stereon::Detail
