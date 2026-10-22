// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>

export module Stereon.Surface:Certify;

import Stereon.Curve;
import :Erased;

/// Certificates (ADR-0009, section 4): a bound on a curve's deviation from
/// a surface that holds over a whole parameter range, not at samples. The
/// bound comes from interval evaluation: IntervalScalar is a Scalar
/// (ADR-0021) whose arithmetic is core's outward-rounded Interval, so the
/// geometry templates, given it, evaluate over parameter intervals and
/// return rigorous enclosures. CertifyDeviation bisects the range until the
/// enclosure of |c(t) - S(p(t))| on every span is within the target.

export namespace Stereon {

/// A real number known to lie within an interval, with the Scalar
/// interface: every operation returns an interval containing every result
/// the operands' numbers could give (ADR-0002). The default value is the
/// number 0, so that T{} in generic code is zero. The ordering is the
/// certain one, a < b when every number of a is below every number of b,
/// which a geometric decision on enclosures must read as "undecided" when
/// neither holds.
///
/// Every operation is defined out of line, in Certify.cpp: inline
/// functions here that call core's Interval arithmetic make GCC 16.2 write
/// a module file its importers cannot read
/// (docs/toolchain/gcc16-feature-probe.md, "Bad file data").
struct IntervalScalar {
    Interval Value;

    IntervalScalar();
    // Implicit, so that constants written T{2} in generic code are enclosures
    // of the constant.
    IntervalScalar(double x); // NOLINT(google-explicit-constructor)
    explicit IntervalScalar(const Interval& value);

    bool operator==(const IntervalScalar&) const = default;

    friend bool operator<(const IntervalScalar& a, const IntervalScalar& b);
    friend bool operator>(const IntervalScalar& a, const IntervalScalar& b);
    friend IntervalScalar operator+(const IntervalScalar& a, const IntervalScalar& b);
    friend IntervalScalar operator-(const IntervalScalar& a, const IntervalScalar& b);
    friend IntervalScalar operator-(const IntervalScalar& a);
    friend IntervalScalar operator*(const IntervalScalar& a, const IntervalScalar& b);
    friend IntervalScalar operator/(const IntervalScalar& a, const IntervalScalar& b);

    // Found by ADL from generic code that writes `using std::sin; sin(x);`
    // (Stereon.Core:Trigonometry): enclosures of the functions over the
    // interval.
    friend IntervalScalar sin(const IntervalScalar& x);
    friend IntervalScalar cos(const IntervalScalar& x);
    /// Over the non-negative part of the interval.
    friend IntervalScalar sqrt(const IntervalScalar& x);
};

/// What CertifyDeviation proves: the bound, and the spans it took. Its own
/// type, because a declaration returning Result<double> here makes GCC
/// 16.2 write a module file its importers cannot read ("Bad file data";
/// the Curve module declares that specialisation too).
struct DeviationCertificate {
    /// A bound on sup |c(t) - S(p(t))| over the range, in millimetres.
    double Bound = 0.0;
    /// How many spans the range was bisected into to prove it.
    std::size_t Spans = 0;

    bool operator==(const DeviationCertificate&) const = default;
};

/// A proved bound on the deviation of a 3D curve from a surface along a
/// pcurve, sup over `range` of |c(t) - S(p(t))| (ADR-0009, section 4): the
/// range is bisected until the enclosure on every span is within `target`,
/// and the largest enclosure is the bound. The enclosures are mean-value
/// forms, f(m) + f'(span) (t - m), from the geometry templates evaluated
/// over IntervalScalar; a B-spline's derivative over a span is the hull of
/// its derivative curve's poles.
///
/// Covers lines, circles, ellipses and non-rational B-splines as the curve,
/// lines, circles, ellipses and non-rational B-splines as the pcurve, and
/// planes and cylinders as the surface; NotSupported otherwise. InvalidInput
/// for an empty or unbounded range or a target that is not finite and
/// positive. CertificationFailed when the span budget runs out before every
/// span is within the target, with the message naming the span and the
/// bound it reached; nothing is sampled, so the failure is never a guess.
[[nodiscard]] Result<DeviationCertificate> CertifyDeviation(const Curve& curve,
                                                            const Surface& surface,
                                                            const Curve2d& pcurve,
                                                            const Interval& range, double target);

} // namespace Stereon
