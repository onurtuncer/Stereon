// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Robust:Predicates;

import Stereon.Core;
import :Domain;
import :Expansion;
import :Sign;

/// The exact geometric predicates of ADR-0002: Orient2D, Orient3D, InCircle
/// and InSphere, after Shewchuk (1997). Each is the sign of a determinant,
/// computed adaptively:
///
/// - Stage A evaluates the determinant in floating point and compares it with
///   an error bound, a constant times the "permanent", the same formula with
///   every subtraction replaced by an addition of magnitudes. Almost every
///   input is decided here, at the cost of the plain formula plus a compare.
/// - Stage B computes the determinant of the rounded coordinate differences
///   exactly, as an expansion, and tests the estimate against a tighter
///   bound. When the differences had no roundoff (coordinates of similar
///   magnitude) the result is exact and returned at once. Orient2D and
///   Orient3D add Shewchuk's stage C, a first-order correction in the
///   roundoff of the differences.
/// - The exact stage evaluates the determinant of the raw coordinates with
///   expansion arithmetic, without any translation, and reads the sign off the
///   result. It is reached only for degenerate or nearly degenerate input.
///
/// The constants are Shewchuk's, for round-to-nearest-even doubles (unit
/// roundoff 2^-53) and no fused multiply-add. Shewchuk's analysis also
/// assumes that nothing underflows or overflows; stage A therefore carries an
/// absolute term, UnderflowGuard, and the adaptive and exact stages first
/// check the numerical domain of :Domain and rescale the coordinates by an
/// exact power of two when needed. Coordinates must be finite, and their
/// binary exponents must span less than about 2045 / d orders for a
/// predicate of degree d; otherwise the exact stage reports a contract
/// violation. The manual's "Exact predicates" chapter derives the
/// determinants, the bounds and the domain.

namespace Stereon::Detail {

/// The unit roundoff of Shewchuk's analysis, 2^-53: half of
/// std::numeric_limits<double>::epsilon().
inline constexpr double Epsilon = 0x1p-53;

inline constexpr double ResultBound = (3.0 + 8.0 * Epsilon) * Epsilon;
inline constexpr double Orient2DBoundA = (3.0 + 16.0 * Epsilon) * Epsilon;
inline constexpr double Orient2DBoundB = (2.0 + 12.0 * Epsilon) * Epsilon;
inline constexpr double Orient2DBoundC = (9.0 + 64.0 * Epsilon) * Epsilon * Epsilon;
inline constexpr double Orient3DBoundA = (7.0 + 56.0 * Epsilon) * Epsilon;
inline constexpr double Orient3DBoundB = (3.0 + 28.0 * Epsilon) * Epsilon;
inline constexpr double Orient3DBoundC = (26.0 + 288.0 * Epsilon) * Epsilon * Epsilon;
inline constexpr double InCircleBoundA = (10.0 + 96.0 * Epsilon) * Epsilon;
inline constexpr double InCircleBoundB = (4.0 + 48.0 * Epsilon) * Epsilon;
inline constexpr double InSphereBoundA = (16.0 + 224.0 * Epsilon) * Epsilon;
inline constexpr double InSphereBoundB = (5.0 + 72.0 * Epsilon) * Epsilon;

// --- Floating-point pieces (stage A) ----------------------------------------

/// The 2x2 minor px*qy - py*qx.
[[nodiscard]] constexpr double Minor(const Vec2& p, const Vec2& q) noexcept
{
    return p.X * q.Y - p.Y * q.X;
}

/// The minor with magnitudes added: |px*qy| + |py*qx|.
[[nodiscard]] constexpr double MinorPermanent(const Vec2& p, const Vec2& q) noexcept
{
    return Magnitude(p.X * q.Y) + Magnitude(p.Y * q.X);
}

[[nodiscard]] constexpr Vec2 Xy(const Vec3& v) noexcept
{
    return {v.X, v.Y};
}

/// The 3x3 determinant of p, q, r as rows, expanded along the z column.
[[nodiscard]] constexpr double Det3(const Vec3& p, const Vec3& q, const Vec3& r) noexcept
{
    return p.Z * Minor(Xy(q), Xy(r)) + q.Z * Minor(Xy(r), Xy(p)) + r.Z * Minor(Xy(p), Xy(q));
}

[[nodiscard]] constexpr double Det3Permanent(const Vec3& p, const Vec3& q, const Vec3& r) noexcept
{
    return MinorPermanent(Xy(q), Xy(r)) * Magnitude(p.Z) +
           MinorPermanent(Xy(r), Xy(p)) * Magnitude(q.Z) +
           MinorPermanent(Xy(p), Xy(q)) * Magnitude(r.Z);
}

/// The squared norm: the lifted coordinate of the circle and sphere tests.
[[nodiscard]] constexpr double Lift(const Vec2& v) noexcept
{
    return v.X * v.X + v.Y * v.Y;
}
[[nodiscard]] constexpr double Lift(const Vec3& v) noexcept
{
    return v.X * v.X + v.Y * v.Y + v.Z * v.Z;
}

// --- Exact pieces (stages B and exact) --------------------------------------

[[nodiscard]] constexpr Expansion<4> ExactMinor(const double px, const double py, const double qx,
                                                const double qy) noexcept
{
    return TwoProduct(px, qy) - TwoProduct(py, qx);
}
[[nodiscard]] constexpr Expansion<4> ExactMinor(const Vec2& p, const Vec2& q) noexcept
{
    return ExactMinor(p.X, p.Y, q.X, q.Y);
}
[[nodiscard]] constexpr Expansion<4> ExactMinor(const Point2& p, const Point2& q) noexcept
{
    return ExactMinor(p.X, p.Y, q.X, q.Y);
}
[[nodiscard]] constexpr Expansion<4> ExactMinor(const Point3& p, const Point3& q) noexcept
{
    return ExactMinor(p.X, p.Y, q.X, q.Y);
}

/// The 3x3 determinant of translated points, exactly.
[[nodiscard]] constexpr Expansion<24> ExactDet3(const Vec3& p, const Vec3& q,
                                                const Vec3& r) noexcept
{
    return ExactMinor(Xy(q), Xy(r)) * p.Z + ExactMinor(Xy(r), Xy(p)) * q.Z +
           ExactMinor(Xy(p), Xy(q)) * r.Z;
}

/// The 3x3 determinant of raw points as rows, exactly.
[[nodiscard]] constexpr Expansion<24> ExactDet3(const Point3& p, const Point3& q,
                                                const Point3& r) noexcept
{
    return ExactMinor(q, r) * p.Z + ExactMinor(r, p) * q.Z + ExactMinor(p, q) * r.Z;
}

/// The orientation determinant of raw 2D points: the 3x3 determinant of
/// rows (x, y, 1), which is the sum of the three cyclic minors.
[[nodiscard]] constexpr Expansion<12> ExactOrient2(const Point2& a, const Point2& b,
                                                   const Point2& c) noexcept
{
    return ExactMinor(a, b) + ExactMinor(b, c) + ExactMinor(c, a);
}

/// The orientation determinant of raw 3D points: the 4x4 determinant of
/// rows (x, y, z, 1), expanded along the column of ones.
[[nodiscard]] constexpr Expansion<96> ExactDet4(const Point3& p, const Point3& q, const Point3& r,
                                                const Point3& s) noexcept
{
    return -ExactDet3(q, r, s) + ExactDet3(p, r, s) - ExactDet3(p, q, s) + ExactDet3(p, q, r);
}

[[nodiscard]] constexpr Expansion<4> ExactLift(const double x, const double y) noexcept
{
    return TwoProduct(x, x) + TwoProduct(y, y);
}
[[nodiscard]] constexpr Expansion<6> ExactLift(const double x, const double y,
                                               const double z) noexcept
{
    return ExactLift(x, y) + TwoProduct(z, z);
}
[[nodiscard]] constexpr Expansion<4> ExactLift(const Vec2& v) noexcept
{
    return ExactLift(v.X, v.Y);
}
[[nodiscard]] constexpr Expansion<4> ExactLift(const Point2& p) noexcept
{
    return ExactLift(p.X, p.Y);
}
[[nodiscard]] constexpr Expansion<6> ExactLift(const Vec3& v) noexcept
{
    return ExactLift(v.X, v.Y, v.Z);
}
[[nodiscard]] constexpr Expansion<6> ExactLift(const Point3& p) noexcept
{
    return ExactLift(p.X, p.Y, p.Z);
}

// --- Roundoff of the coordinate differences ---------------------------------

/// The roundoff of diff = fl(a - b), per coordinate.
[[nodiscard]] constexpr Vec2 Tail(const Point2& a, const Point2& b, const Vec2& diff) noexcept
{
    return {TwoDiffTail(a.X, b.X, diff.X), TwoDiffTail(a.Y, b.Y, diff.Y)};
}
[[nodiscard]] constexpr Vec3 Tail(const Point3& a, const Point3& b, const Vec3& diff) noexcept
{
    return {TwoDiffTail(a.X, b.X, diff.X), TwoDiffTail(a.Y, b.Y, diff.Y),
            TwoDiffTail(a.Z, b.Z, diff.Z)};
}
[[nodiscard]] constexpr bool IsZero(const Vec2& v) noexcept
{
    return v.X == 0.0 && v.Y == 0.0;
}
[[nodiscard]] constexpr bool IsZero(const Vec3& v) noexcept
{
    return v.X == 0.0 && v.Y == 0.0 && v.Z == 0.0;
}

} // namespace Stereon::Detail

export namespace Stereon {

// --- Exact stages -----------------------------------------------------------

/// Orient2D without filters: the exact sign, at the cost of about a dozen
/// TwoProducts. The adaptive predicates below fall back to these.
[[nodiscard]] constexpr Sign Orient2DExact(const Point2& a, const Point2& b,
                                           const Point2& c) noexcept
{
    const Detail::DomainCheck domain = Detail::CheckDomain({a, b, c}, 2);
    if (domain.Status == Detail::DomainStatus::Scalable) {
        using Detail::ScaledByPowerOfTwo;
        return Orient2DExact(ScaledByPowerOfTwo(a, domain.Shift),
                             ScaledByPowerOfTwo(b, domain.Shift),
                             ScaledByPowerOfTwo(c, domain.Shift));
    }
    if (domain.Status != Detail::DomainStatus::Inside) {
        Detail::PredicateDomainViolation("Orient2D", domain.Status);
    }
    return Detail::ExactOrient2(a, b, c).GetSign();
}

[[nodiscard]] constexpr Sign Orient3DExact(const Point3& a, const Point3& b, const Point3& c,
                                           const Point3& d) noexcept
{
    const Detail::DomainCheck domain = Detail::CheckDomain({a, b, c, d}, 3);
    if (domain.Status == Detail::DomainStatus::Scalable) {
        using Detail::ScaledByPowerOfTwo;
        const int s = domain.Shift;
        return Orient3DExact(ScaledByPowerOfTwo(a, s), ScaledByPowerOfTwo(b, s),
                             ScaledByPowerOfTwo(c, s), ScaledByPowerOfTwo(d, s));
    }
    if (domain.Status != Detail::DomainStatus::Inside) {
        Detail::PredicateDomainViolation("Orient3D", domain.Status);
    }
    return Detail::ExactDet4(a, b, c, d).GetSign();
}

/// The 4x4 determinant of rows (x, y, x^2 + y^2, 1), expanded along the
/// lifted column.
[[nodiscard]] constexpr Sign InCircleExact(const Point2& a, const Point2& b, const Point2& c,
                                           const Point2& d) noexcept
{
    const Detail::DomainCheck domain = Detail::CheckDomain({a, b, c, d}, 4);
    if (domain.Status == Detail::DomainStatus::Scalable) {
        using Detail::ScaledByPowerOfTwo;
        const int s = domain.Shift;
        return InCircleExact(ScaledByPowerOfTwo(a, s), ScaledByPowerOfTwo(b, s),
                             ScaledByPowerOfTwo(c, s), ScaledByPowerOfTwo(d, s));
    }
    if (domain.Status != Detail::DomainStatus::Inside) {
        Detail::PredicateDomainViolation("InCircle", domain.Status);
    }
    using Detail::ExactLift;
    using Detail::ExactOrient2;
    const auto det = ExactLift(a) * ExactOrient2(b, c, d) - ExactLift(b) * ExactOrient2(a, c, d) +
                     ExactLift(c) * ExactOrient2(a, b, d) - ExactLift(d) * ExactOrient2(a, b, c);
    return det.GetSign();
}

/// The 5x5 determinant of rows (x, y, z, x^2 + y^2 + z^2, 1), expanded
/// along the lifted column. The intermediate expansions reach a capacity of
/// 5 760 components, the same as Shewchuk's insphereexact.
[[nodiscard]] constexpr Sign InSphereExact(const Point3& a, const Point3& b, const Point3& c,
                                           const Point3& d, const Point3& e) noexcept
{
    const Detail::DomainCheck domain = Detail::CheckDomain({a, b, c, d, e}, 5);
    if (domain.Status == Detail::DomainStatus::Scalable) {
        using Detail::ScaledByPowerOfTwo;
        const int s = domain.Shift;
        return InSphereExact(ScaledByPowerOfTwo(a, s), ScaledByPowerOfTwo(b, s),
                             ScaledByPowerOfTwo(c, s), ScaledByPowerOfTwo(d, s),
                             ScaledByPowerOfTwo(e, s));
    }
    if (domain.Status != Detail::DomainStatus::Inside) {
        Detail::PredicateDomainViolation("InSphere", domain.Status);
    }
    using Detail::ExactDet4;
    using Detail::ExactLift;
    const auto ta = ExactLift(a) * ExactDet4(b, c, d, e);
    const auto tb = ExactLift(b) * ExactDet4(a, c, d, e);
    const auto tc = ExactLift(c) * ExactDet4(a, b, d, e);
    const auto td = ExactLift(d) * ExactDet4(a, b, c, e);
    const auto te = ExactLift(e) * ExactDet4(a, b, c, d);
    return (((-ta + tb) - tc) + (td - te)).GetSign();
}

} // namespace Stereon

namespace Stereon::Detail {

// --- Adaptive stages (B, C) -------------------------------------------------

[[nodiscard]] constexpr Sign Orient2DAdaptive(const Point2& a, const Point2& b, const Point2& c,
                                              const double detSum) noexcept
{
    // Stages B and C are exact only inside the domain; outside it the exact
    // stage rescales or reports.
    if (CheckDomain({a, b, c}, 2).Status != DomainStatus::Inside) {
        return Orient2DExact(a, b, c);
    }
    const Vec2 ac = a - c;
    const Vec2 bc = b - c;
    const Expansion<4> det = TwoProduct(ac.X, bc.Y) - TwoProduct(ac.Y, bc.X);
    double estimate = det.GetEstimate();
    double errorBound = Orient2DBoundB * detSum;
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }

    const Vec2 acTail = Tail(a, c, ac);
    const Vec2 bcTail = Tail(b, c, bc);
    if (IsZero(acTail) && IsZero(bcTail)) {
        return det.GetSign(); // the differences were exact, so det is
    }

    errorBound = Orient2DBoundC * detSum + ResultBound * Magnitude(estimate);
    estimate += (ac.X * bcTail.Y + bc.Y * acTail.X) - (ac.Y * bcTail.X + bc.X * acTail.Y);
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }
    return Orient2DExact(a, b, c);
}

[[nodiscard]] constexpr Sign Orient3DAdaptive(const Point3& a, const Point3& b, const Point3& c,
                                              const Point3& d, const double permanent) noexcept
{
    if (CheckDomain({a, b, c, d}, 3).Status != DomainStatus::Inside) {
        return Orient3DExact(a, b, c, d);
    }
    const Vec3 ad = a - d;
    const Vec3 bd = b - d;
    const Vec3 cd = c - d;
    const Expansion<24> det = ExactDet3(ad, bd, cd);
    double estimate = det.GetEstimate();
    double errorBound = Orient3DBoundB * permanent;
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }

    const Vec3 adTail = Tail(a, d, ad);
    const Vec3 bdTail = Tail(b, d, bd);
    const Vec3 cdTail = Tail(c, d, cd);
    if (IsZero(adTail) && IsZero(bdTail) && IsZero(cdTail)) {
        return det.GetSign();
    }

    errorBound = Orient3DBoundC * permanent + ResultBound * Magnitude(estimate);
    estimate +=
        (ad.Z * ((bd.X * cdTail.Y + cd.Y * bdTail.X) - (bd.Y * cdTail.X + cd.X * bdTail.Y)) +
         adTail.Z * (bd.X * cd.Y - bd.Y * cd.X)) +
        (bd.Z * ((cd.X * adTail.Y + ad.Y * cdTail.X) - (cd.Y * adTail.X + ad.X * cdTail.Y)) +
         bdTail.Z * (cd.X * ad.Y - cd.Y * ad.X)) +
        (cd.Z * ((ad.X * bdTail.Y + bd.Y * adTail.X) - (ad.Y * bdTail.X + bd.X * adTail.Y)) +
         cdTail.Z * (ad.X * bd.Y - ad.Y * bd.X));
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }
    return Orient3DExact(a, b, c, d);
}

[[nodiscard]] constexpr Sign InCircleAdaptive(const Point2& a, const Point2& b, const Point2& c,
                                              const Point2& d, const double permanent) noexcept
{
    if (CheckDomain({a, b, c, d}, 4).Status != DomainStatus::Inside) {
        return InCircleExact(a, b, c, d);
    }
    const Vec2 ad = a - d;
    const Vec2 bd = b - d;
    const Vec2 cd = c - d;
    const auto det = ExactLift(ad) * ExactMinor(bd, cd) + ExactLift(bd) * ExactMinor(cd, ad) +
                     ExactLift(cd) * ExactMinor(ad, bd);
    const double estimate = det.GetEstimate();
    const double errorBound = InCircleBoundB * permanent;
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }
    if (IsZero(Tail(a, d, ad)) && IsZero(Tail(b, d, bd)) && IsZero(Tail(c, d, cd))) {
        return det.GetSign();
    }
    return InCircleExact(a, b, c, d);
}

[[nodiscard]] constexpr Sign InSphereAdaptive(const Point3& a, const Point3& b, const Point3& c,
                                              const Point3& d, const Point3& e,
                                              const double permanent) noexcept
{
    if (CheckDomain({a, b, c, d, e}, 5).Status != DomainStatus::Inside) {
        return InSphereExact(a, b, c, d, e);
    }
    const Vec3 ae = a - e;
    const Vec3 be = b - e;
    const Vec3 ce = c - e;
    const Vec3 de = d - e;
    const auto det =
        (ExactLift(de) * ExactDet3(ae, be, ce) - ExactLift(ce) * ExactDet3(de, ae, be)) +
        (ExactLift(be) * ExactDet3(ce, de, ae) - ExactLift(ae) * ExactDet3(be, ce, de));
    const double estimate = det.GetEstimate();
    const double errorBound = InSphereBoundB * permanent;
    if (Magnitude(estimate) >= errorBound) {
        return SignOf(estimate);
    }
    if (IsZero(Tail(a, e, ae)) && IsZero(Tail(b, e, be)) && IsZero(Tail(c, e, ce)) &&
        IsZero(Tail(d, e, de))) {
        return det.GetSign();
    }
    return InSphereExact(a, b, c, d, e);
}

} // namespace Stereon::Detail

export namespace Stereon {

// --- Adaptive predicates ----------------------------------------------------

/// Positive when a, b, c turn counterclockwise, negative when clockwise, zero
/// when they are collinear: the sign of the 2x2 determinant |a - c, b - c|.
[[nodiscard]] constexpr Sign Orient2D(const Point2& a, const Point2& b, const Point2& c) noexcept
{
    const double detLeft = (a.X - c.X) * (b.Y - c.Y);
    const double detRight = (a.Y - c.Y) * (b.X - c.X);
    const double det = detLeft - detRight;

    // Shewchuk's orient2d returns early when the two terms cannot cancel
    // (opposite signs, or one of them zero). The bound test below gives the
    // same answer in those cases, since |det| is then the whole permanent,
    // and it does so without a branch on the sign of random data.
    const double detSum = Detail::Magnitude(detLeft) + Detail::Magnitude(detRight);
    const double errorBound = Detail::Orient2DBoundA * detSum + Detail::UnderflowGuard;
    // Strict, unlike Shewchuk's >=: an overflowed product gives inf >= inf.
    if (Detail::Magnitude(det) > errorBound) {
        return SignOf(det);
    }
    return Detail::Orient2DAdaptive(a, b, c, detSum);
}

/// Positive when d lies below the plane through a, b, c, with a, b, c seen
/// counterclockwise from above; negative when above; zero when the four
/// points are coplanar. The sign of the 3x3 determinant |a - d, b - d, c - d|.
[[nodiscard]] constexpr Sign Orient3D(const Point3& a, const Point3& b, const Point3& c,
                                      const Point3& d) noexcept
{
    const Vec3 ad = a - d;
    const Vec3 bd = b - d;
    const Vec3 cd = c - d;
    const double det = Detail::Det3(ad, bd, cd);
    const double permanent = Detail::Det3Permanent(ad, bd, cd);
    const double errorBound = Detail::Orient3DBoundA * permanent + Detail::UnderflowGuard;
    if (Detail::Magnitude(det) > errorBound) {
        return SignOf(det);
    }
    return Detail::Orient3DAdaptive(a, b, c, d, permanent);
}

/// Positive when d lies inside the circle through a, b, c, negative when
/// outside, zero when the four points are cocircular. a, b, c must be in
/// counterclockwise order (Orient2D(a, b, c) positive), or the sign flips.
[[nodiscard]] constexpr Sign InCircle(const Point2& a, const Point2& b, const Point2& c,
                                      const Point2& d) noexcept
{
    const Vec2 ad = a - d;
    const Vec2 bd = b - d;
    const Vec2 cd = c - d;
    const double aLift = Detail::Lift(ad);
    const double bLift = Detail::Lift(bd);
    const double cLift = Detail::Lift(cd);
    const double det = aLift * Detail::Minor(bd, cd) + bLift * Detail::Minor(cd, ad) +
                       cLift * Detail::Minor(ad, bd);
    const double permanent = Detail::MinorPermanent(bd, cd) * aLift +
                             Detail::MinorPermanent(cd, ad) * bLift +
                             Detail::MinorPermanent(ad, bd) * cLift;
    const double errorBound = Detail::InCircleBoundA * permanent + Detail::UnderflowGuard;
    if (Detail::Magnitude(det) > errorBound) {
        return SignOf(det);
    }
    return Detail::InCircleAdaptive(a, b, c, d, permanent);
}

/// Positive when e lies inside the sphere through a, b, c, d, negative when
/// outside, zero when the five points are cospherical. a, b, c, d must be
/// positively oriented (Orient3D(a, b, c, d) positive), or the sign flips.
[[nodiscard]] constexpr Sign InSphere(const Point3& a, const Point3& b, const Point3& c,
                                      const Point3& d, const Point3& e) noexcept
{
    const Vec3 ae = a - e;
    const Vec3 be = b - e;
    const Vec3 ce = c - e;
    const Vec3 de = d - e;
    const double aLift = Detail::Lift(ae);
    const double bLift = Detail::Lift(be);
    const double cLift = Detail::Lift(ce);
    const double dLift = Detail::Lift(de);
    const double det = (dLift * Detail::Det3(ae, be, ce) - cLift * Detail::Det3(de, ae, be)) +
                       (bLift * Detail::Det3(ce, de, ae) - aLift * Detail::Det3(be, ce, de));
    const double permanent =
        Detail::Det3Permanent(be, ce, de) * aLift + Detail::Det3Permanent(ce, de, ae) * bLift +
        Detail::Det3Permanent(de, ae, be) * cLift + Detail::Det3Permanent(ae, be, ce) * dLift;
    const double errorBound = Detail::InSphereBoundA * permanent + Detail::UnderflowGuard;
    if (Detail::Magnitude(det) > errorBound) {
        return SignOf(det);
    }
    return Detail::InSphereAdaptive(a, b, c, d, e, permanent);
}

} // namespace Stereon
