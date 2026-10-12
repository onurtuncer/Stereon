// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// An independent referee for the robust predicates (ADR-0015): the same
/// determinants evaluated in 128-bit integer arithmetic on integer-valued
/// coordinates, where every operation is exact. Test inputs are integers, or
/// dyadic rationals scaled to integers by a power of two, which leaves the
/// floating-point behaviour of the predicate under test unchanged.
///
/// Magnitude limits, so that nothing overflows the 127 bits available:
/// Orient2D coordinates below 2^60, Orient3D below 2^30, InCircle below
/// 2^28, InSphere below 2^20. Compilers without __int128 (MSVC) get no
/// oracle; the tests that need it are skipped there.

#pragma once

#include <cmath>
#include <cstdint>

#if defined(__SIZEOF_INT128__)
#define STEREON_HAS_EXACT_ORACLE 1
#else
#define STEREON_HAS_EXACT_ORACLE 0
#endif

#if STEREON_HAS_EXACT_ORACLE

namespace Stereon::Tests {

// __extension__ keeps -Wpedantic quiet about the non-standard type.
__extension__ typedef __int128 Wide; // NOLINT(modernize-use-using)

/// The integer a double holds; the double must be integer-valued.
[[nodiscard]] inline Wide ToWide(const double x)
{
    return static_cast<Wide>(x);
}

[[nodiscard]] inline bool IsIntegerValued(const double x)
{
    return std::isfinite(x) && x == std::trunc(x);
}

[[nodiscard]] inline Stereon::Sign SignOfWide(const Wide v)
{
    if (v > 0) {
        return Stereon::Sign::Positive;
    }
    return v < 0 ? Stereon::Sign::Negative : Stereon::Sign::Zero;
}

struct WidePoint2 {
    Wide X;
    Wide Y;
};

struct WidePoint3 {
    Wide X;
    Wide Y;
    Wide Z;
};

[[nodiscard]] inline WidePoint2 Widen(const Stereon::Point2& p)
{
    return {ToWide(p.X), ToWide(p.Y)};
}

[[nodiscard]] inline WidePoint3 Widen(const Stereon::Point3& p)
{
    return {ToWide(p.X), ToWide(p.Y), ToWide(p.Z)};
}

/// px*qy - py*qx.
[[nodiscard]] inline Wide Minor(const WidePoint2& p, const WidePoint2& q)
{
    return p.X * q.Y - p.Y * q.X;
}

[[nodiscard]] inline WidePoint2 Xy(const WidePoint3& p)
{
    return {p.X, p.Y};
}

/// The determinant of rows (x, y, 1).
[[nodiscard]] inline Wide Orient2(const WidePoint2& a, const WidePoint2& b, const WidePoint2& c)
{
    return Minor(a, b) + Minor(b, c) + Minor(c, a);
}

/// The determinant of rows (x, y, z).
[[nodiscard]] inline Wide Det3(const WidePoint3& p, const WidePoint3& q, const WidePoint3& r)
{
    return p.Z * Minor(Xy(q), Xy(r)) + q.Z * Minor(Xy(r), Xy(p)) + r.Z * Minor(Xy(p), Xy(q));
}

/// The determinant of rows (x, y, z, 1).
[[nodiscard]] inline Wide Det4(const WidePoint3& p, const WidePoint3& q, const WidePoint3& r,
                               const WidePoint3& s)
{
    return -Det3(q, r, s) + Det3(p, r, s) - Det3(p, q, s) + Det3(p, q, r);
}

[[nodiscard]] inline Wide Lift(const WidePoint2& p)
{
    return p.X * p.X + p.Y * p.Y;
}

[[nodiscard]] inline Wide Lift(const WidePoint3& p)
{
    return p.X * p.X + p.Y * p.Y + p.Z * p.Z;
}

[[nodiscard]] inline Stereon::Sign ExactOrient2D(const Stereon::Point2& a, const Stereon::Point2& b,
                                                 const Stereon::Point2& c)
{
    return SignOfWide(Orient2(Widen(a), Widen(b), Widen(c)));
}

[[nodiscard]] inline Stereon::Sign ExactOrient3D(const Stereon::Point3& a, const Stereon::Point3& b,
                                                 const Stereon::Point3& c, const Stereon::Point3& d)
{
    return SignOfWide(Det4(Widen(a), Widen(b), Widen(c), Widen(d)));
}

/// The determinant of rows (x, y, x^2 + y^2, 1), expanded along the lifted
/// column.
[[nodiscard]] inline Stereon::Sign ExactInCircle(const Stereon::Point2& a, const Stereon::Point2& b,
                                                 const Stereon::Point2& c, const Stereon::Point2& d)
{
    const WidePoint2 wa = Widen(a);
    const WidePoint2 wb = Widen(b);
    const WidePoint2 wc = Widen(c);
    const WidePoint2 wd = Widen(d);
    const Wide det = Lift(wa) * Orient2(wb, wc, wd) - Lift(wb) * Orient2(wa, wc, wd) +
                     Lift(wc) * Orient2(wa, wb, wd) - Lift(wd) * Orient2(wa, wb, wc);
    return SignOfWide(det);
}

/// The determinant of rows (x, y, z, x^2 + y^2 + z^2, 1), expanded along the
/// lifted column.
[[nodiscard]] inline Stereon::Sign ExactInSphere(const Stereon::Point3& a, const Stereon::Point3& b,
                                                 const Stereon::Point3& c, const Stereon::Point3& d,
                                                 const Stereon::Point3& e)
{
    const WidePoint3 wa = Widen(a);
    const WidePoint3 wb = Widen(b);
    const WidePoint3 wc = Widen(c);
    const WidePoint3 wd = Widen(d);
    const WidePoint3 we = Widen(e);
    const Wide det = -Lift(wa) * Det4(wb, wc, wd, we) + Lift(wb) * Det4(wa, wc, wd, we) -
                     Lift(wc) * Det4(wa, wb, wd, we) + Lift(wd) * Det4(wa, wb, wc, we) -
                     Lift(we) * Det4(wa, wb, wc, wd);
    return SignOfWide(det);
}

} // namespace Stereon::Tests

#endif // STEREON_HAS_EXACT_ORACLE
