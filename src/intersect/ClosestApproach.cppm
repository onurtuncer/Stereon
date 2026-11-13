// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>

export module Stereon.Intersect:ClosestApproach;

import Stereon.Surface;
import :Dispatch;
import :Linearization;

/// The closest approach of two space curves (ADR-0052, option 1, and its
/// addendum of 2026-11-12): the pair of parameters (s, t) whose points A(s)
/// and B(t) are nearest over the requested domains, with its own forward
/// and reverse products. It is not Intersect: Intersect's products for two
/// space curves stay those of ADR-0052's decision, and ClosestApproach
/// describes the approach whether or not the curves meet.
///
/// The minimum is proved global, or the result says it is not. Two lines
/// are a convex quadratic, decided exactly parallel or not. A line and a
/// conic, or two conics, are searched by branch and bound on |A(s) -
/// B(t)|^2 with second-order lower bounds; the boxes that may hold the
/// minimum form clusters, and a single cluster on which the squared
/// distance is proved strictly convex (its Hessian at the minimiser
/// outweighs a bound on its third derivatives over the cluster) holds the
/// unique minimiser. Exactly parallel lines, a line along a circle's axis
/// and two coaxial full circles are proved to have a set of minimisers.
///
/// At an interior minimiser the parameters solve G(s, t) = (d . A', -d . B')
/// = 0 with d = A(s) - B(t); G_z is the Hessian H = [[A'.A' + d.A'',
/// -A'.B'], [-A'.B', B'.B' - d.B'']], and a design change moves them by
/// H z' = -G_p p'. A parameter at an end of its domain is held there and
/// the other solves its one equation (the owner's decision 2 of the
/// addendum).

export namespace Stereon {

/// Where a parameter of the closest approach lies in its domain.
enum class ApproachEnd : std::uint8_t {
    Interior,
    Lower,
    Upper,
};

/// What the primal proved about the minimiser: unique; a set of
/// minimisers (exactly parallel lines, a line on a circle's axis, coaxial
/// full circles), proved; or undecided, two candidates the search could not
/// separate, or one it could not prove alone (the reported pair is then the
/// best found, its distance an upper bound and LowerBound a lower one).
enum class ApproachUniqueness : std::uint8_t {
    Unique,
    NotUnique,
    Undecided,
};

/// The closest approach of A and B: both closest points, their parameters,
/// the gap d = A(s) - B(t) and its length (the owner's decision 1: no
/// single position). Meets says the distance is within the Context's
/// resolution: the curves meet there, and the products answer
/// Nondifferentiable (decision 4).
struct ClosestApproachResult {
    double ParameterA = 0.0;
    double ParameterB = 0.0;
    Point3 OnA;
    Point3 OnB;
    Vec3 Gap;
    double Distance = 0.0;
    /// A certified lower bound on the minimum distance over the domains.
    double LowerBound = 0.0;
    ApproachEnd EndA = ApproachEnd::Interior;
    ApproachEnd EndB = ApproachEnd::Interior;
    ApproachUniqueness Uniqueness = ApproachUniqueness::Unique;
    bool Meets = false;
};

/// The closest approach of a and b over the given domains. Lines, circles
/// and ellipses (decision 5); a line's domain may be unbounded, a conic's is
/// an interval of its angle, the whole conic when 2 pi wide or wider.
/// NotSupported for any other curve (a B-spline until a subdivision proof
/// exists); InvalidInput for an empty curve, a zero direction, a domain
/// that is empty or not finite where it must be.
[[nodiscard]] Result<ClosestApproachResult> ClosestApproach(const Curve& a, const Interval& domainA,
                                                            const Curve& b, const Interval& domainB,
                                                            const Context& ctx);

/// Over each curve's whole domain.
[[nodiscard]] Result<ClosestApproachResult> ClosestApproach(const Curve& a, const Curve& b,
                                                            const Context& ctx);

/// The design parameters bound to one curve of a closest approach: a
/// line's (LineBinding) or a circle's or an ellipse's (ConicBinding),
/// whichever the curve is; a binding of the other kind is InvalidInput.
struct CurveBindings {
    std::span<const LineBinding> Line;
    std::span<const ConicBinding> Conic;
};

/// The values the products give, in this order: A(s)' (three), B(t)'
/// (three), s', t', the gap's d' (three) and the distance's |d|'.
inline constexpr std::size_t ClosestApproachValues = 12;

/// The forward product of a closest approach. Nondifferentiable when the
/// minimisers are a set, the curves meet, or a parameter sits at an end of
/// its domain with the one-sided gradient within rounding of zero (where
/// the minimiser leaves the end: decision 2); Indeterminate when the primal
/// left the minimiser undecided (decision 3) or its Hessian is singular;
/// NearEvent when the normalised Hessian's smaller eigenvalue, or an end's
/// normalised gradient, is below NearEventSine (decision 6); Valid
/// otherwise. NotSupported for curves other than lines, circles and
/// ellipses. InvalidInput for a binding of another curve kind, of a value
/// the curve does not have or of the wrong quantity; for a result that does
/// not lie on the current curves; or for a seed that was never set.
[[nodiscard]] Result<SensitivityResult>
ClosestApproachJvp(const Curve& a, const CurveBindings& bindingsA, const Curve& b,
                   const CurveBindings& bindingsB, const ClosestApproachResult& approach,
                   const SensitivityContext& sensitivity);

/// The reverse product: for an adjoint of ClosestApproachValues
/// components, one value per parameter of the set. InvalidInput as the
/// forward product, and for an adjoint of another length.
[[nodiscard]] Result<SensitivityResult>
ClosestApproachVjp(const Curve& a, const CurveBindings& bindingsA, const Curve& b,
                   const CurveBindings& bindingsB, const ClosestApproachResult& approach,
                   const SensitivityContext& sensitivity);

} // namespace Stereon
