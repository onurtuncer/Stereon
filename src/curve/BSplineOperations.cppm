// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <vector>

export module Stereon.Curve:BSplineOperations;

import Stereon.Robust;
import :BSplineBasis;
import :BSplineCurve;
import :Concepts;

/// The B-spline curve operations beyond knot insertion (docs/nurbs.rst,
/// *Degree elevation, knot removal and the other rearrangements*): the
/// continuity and closure queries, clamping and the periodic form, degree
/// elevation, knot removal and concatenation. Each is a free function on
/// BSplineCurve and BSplineCurve2d, defined in BSplineOperations.cpp on the
/// homogeneous form of the curve, so that a rational curve is handled as
/// the polynomial curve (w P, w).

export namespace Stereon {

// ---------------------------------------------------------------------------
// Queries

/// The order of continuity the knot vector guarantees across the domain:
/// p - s for the largest multiplicity s of a knot strictly inside the
/// domain, and the degree when there is none (one polynomial span, which is
/// smoother still). The curve may be smoother than its knots say, when a
/// knot is removable exactly (RemoveKnot); this query does not look.
template <Scalar T, std::size_t D>
[[nodiscard]] int Continuity(const BasicBSplineCurve<T, D>& curve)
{
    const Interval domain = Domain(curve);
    std::size_t largest = 0;
    std::size_t run = 0;
    for (std::size_t i = 0; i < curve.Knots.size(); ++i) {
        run = (i > 0 && curve.Knots[i] == curve.Knots[i - 1]) ? run + 1 : 1;
        if (curve.Knots[i] > domain.GetLower() && curve.Knots[i] < domain.GetUpper()) {
            largest = run > largest ? run : largest;
        }
    }
    return curve.Degree - static_cast<int>(largest);
}

/// The order of continuity across the parameter t: p - s when t is a knot
/// of multiplicity s strictly inside the domain, the degree otherwise (no
/// knot, or a domain end, where there is no other side).
template <Scalar T, std::size_t D>
[[nodiscard]] int Continuity(const BasicBSplineCurve<T, D>& curve, const double t)
{
    const Interval domain = Domain(curve);
    if (!(t > domain.GetLower() && t < domain.GetUpper())) {
        return curve.Degree;
    }
    return curve.Degree - static_cast<int>(KnotMultiplicity(curve.Knots, t));
}

/// Whether the curve returns to its start: the distance between its points
/// at the two ends of the domain is at most `tolerance`. A periodic form
/// (IsPeriodic) is closed up to rounding.
[[nodiscard]] bool IsClosed(const BSplineCurve& curve, double tolerance);
[[nodiscard]] bool IsClosed(const BSplineCurve2d& curve, double tolerance);

/// Whether the curve is in the periodic form of a closed curve, the one the
/// periodic fits produce and AsPeriodic builds: with n + p poles the last p
/// repeat the first p (and their weights), and the knots repeat with the
/// width T of the domain, u_{i+n} = u_i + T up to rounding, so that the
/// curve continues across either end of its domain with the continuity the
/// seam knot's multiplicity gives. Its period is the width of its domain.
/// The rounding allowed a knot is that of a double of its size, and it must
/// stay below 1/256 of the shortest knot interval of the domain: knots so
/// far from the origin that it does not (about 1e13 for unit intervals)
/// cannot resolve the period, and such a curve is not recognised whatever
/// its values are.
[[nodiscard]] bool IsPeriodic(const BSplineCurve& curve);
[[nodiscard]] bool IsPeriodic(const BSplineCurve2d& curve);

// ---------------------------------------------------------------------------
// Clamping and the periodic form

/// The same curve on the same domain with both end knots at multiplicity
/// p + 1, so that it starts at its first pole and ends at its last: the
/// domain ends are inserted up to multiplicity p (Boehm's algorithm) and
/// the knots and poles beyond them, which the curve on the domain never
/// used, are dropped. A clamped curve is returned unchanged.
[[nodiscard]] BSplineCurve Clamped(const BSplineCurve& curve);
[[nodiscard]] BSplineCurve2d Clamped(const BSplineCurve2d& curve);

/// The periodic form of a closed curve, with the period the width of its
/// domain: the curve is clamped, its two end poles, which must be within
/// `tolerance` of each other, are merged at their midpoint, and the knots
/// and poles are written out as IsPeriodic describes, with the seam knot at
/// multiplicity p, so that the curve on the domain is unchanged except for
/// the merge, which moves it by at most half the distance between the end
/// poles. A curve already in the periodic form is returned unchanged.
/// InvalidInput unless the tolerance is finite and non-negative, the end
/// poles are within it and, for a rational curve, the end weights are
/// equal (scale the weights of one side first otherwise). NotSupported when
/// no valid curve holds the periodic form: a knot copy shifted by the
/// period overflows, or the extended knot vector is not of finite width.
[[nodiscard]] Result<BSplineCurve> AsPeriodic(const BSplineCurve& curve, double tolerance);
[[nodiscard]] Result<BSplineCurve2d> AsPeriodic(const BSplineCurve2d& curve, double tolerance);

/// A curve in the periodic form opened at t: the clamped curve over
/// [t, t + Period(curve)] that traces the curve once from t round to t
/// again, with the curve's parametrisation, cut from the form written out
/// over one more period on each side. t may lie in any period. The piece
/// closes exactly, its last pole and weight its first, and is C^0 there;
/// AsPeriodic turns it back into the form.
/// InvalidInput unless the curve is in the periodic form and t is finite.
[[nodiscard]] Result<BSplineCurve> Opened(const BSplineCurve& curve, double t);
[[nodiscard]] Result<BSplineCurve2d> Opened(const BSplineCurve2d& curve, double t);

// ---------------------------------------------------------------------------
// Degree elevation and knot removal

/// The same curve as a B-spline of the higher `degree` (A5.9): the clamped
/// form is split into its Bezier segments, each is raised by repeated
/// elevation Q_i = (i / (q + 1)) P_{i-1} + (1 - i / (q + 1)) P_i, and the
/// segments are joined back with every interior knot at its old
/// multiplicity plus the degree increase. That join is exact in real
/// arithmetic, but it inverts knot insertion, which closely spaced knots
/// make ill-conditioned, so each removal is kept only while the bound on
/// the curve's movement (the one RemoveKnot reports) stays within
/// 2^-40 (1 + max |P|); otherwise the knot keeps its remaining copies.
/// The result therefore has as few poles as the arithmetic allows, is
/// within that bound of the curve everywhere, and is clamped; a curve
/// already of that degree is returned unchanged. InvalidInput unless
/// curve.Degree <= degree <= BSplineMaxDegree.
[[nodiscard]] Result<BSplineCurve> ElevatedDegree(const BSplineCurve& curve, int degree);
[[nodiscard]] Result<BSplineCurve2d> ElevatedDegree(const BSplineCurve2d& curve, int degree);

/// What RemoveKnot did.
template <std::size_t D>
struct BasicKnotRemoval {
    /// The curve with the knot removed `Removed` times.
    BasicBSplineCurve<double, D> Curve;
    /// How many times the knot was removed, 0 to the `times` asked for.
    int Removed = 0;
    /// A bound on how far any point of Curve lies from the point of the
    /// input at the same parameter: the sum, over the removals made, of
    /// the distance by which the poles failed to agree (A5.8), scaled for a
    /// rational curve by (1 + the largest pole norm) / the smallest weight.
    /// At most the tolerance given; 0 when nothing was removed.
    double Deviation = 0.0;
};

using KnotRemoval = BasicKnotRemoval<3>;
using KnotRemoval2d = BasicKnotRemoval<2>;

/// The curve with the knot t removed up to `times` times, while the curve
/// moves by at most `tolerance` (A5.8): each removal computes the poles
/// that reproduce the curve without the knot from both sides, measures how
/// far the two computations disagree, and is kept only while the bound
/// Deviation stays within the tolerance; the first removal that would
/// exceed it stops the sequence, so Removed may be less than `times` and
/// the curve is then returned as far as it got. A tolerance of 0 removes
/// only what the rounding of the poles allows. InvalidInput unless t is a
/// knot strictly inside the domain, 1 <= times <= its multiplicity and the
/// tolerance is finite and non-negative.
[[nodiscard]] Result<KnotRemoval> RemoveKnot(const BSplineCurve& curve, double t, int times,
                                             double tolerance);
[[nodiscard]] Result<KnotRemoval2d> RemoveKnot(const BSplineCurve2d& curve, double t, int times,
                                               double tolerance);

// ---------------------------------------------------------------------------
// Concatenation

/// One curve tracing `first` and then `second`: both are clamped, the one
/// of lower degree is elevated to the other's, the second's domain is
/// shifted to start where the first's ends and its weights are scaled so
/// that the two poles at the junction, which must be within `tolerance`,
/// carry the same weight; those two poles are merged at their midpoint,
/// which moves each curve by at most half their distance. The junction
/// knot has multiplicity p, continuity C^0; RemoveKnot there, with a
/// tolerance, finds whatever continuity the two curves actually have. The
/// parametrisation of each piece is kept. InvalidInput unless the tolerance
/// is finite and non-negative and the end of `first` and the start of
/// `second` are within it. NotSupported when no valid curve holds the join:
/// shifting the second curve's knots onto the first's rounds one of its
/// intervals away (a large knot origin against a small span), or the scaled
/// weights leave [MinWeight, MaxWeight] even after both pieces are rescaled
/// by a common power of two, which changes no point.
[[nodiscard]] Result<BSplineCurve> Concatenated(const BSplineCurve& first,
                                                const BSplineCurve& second, double tolerance);
[[nodiscard]] Result<BSplineCurve2d> Concatenated(const BSplineCurve2d& first,
                                                  const BSplineCurve2d& second, double tolerance);

} // namespace Stereon
