// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <algorithm>
#include <bit>
#include <climits>
#include <cstdint>
#include <initializer_list>

export module Stereon.Robust:Domain;

import Stereon.Core;

/// The numerical domain of the exact predicates, and the exact rescaling that
/// brings inputs into it.
///
/// Shewchuk's error bounds and the error-free transformations assume that no
/// operation underflows or overflows. Inside the predicates every value is a
/// polynomial of degree at most d in the coordinates (d = 2 for Orient2D, 3
/// for Orient3D, 4 for InCircle, 5 for InSphere), and every intermediate
/// value, rounded or exact, is a multiple of 2^(d * L), where 2^L is the
/// lowest set bit among the non-zero coordinates: differences, products,
/// roundoff terms and the halves of Split() all keep that granularity. So if
///
///     d * L >= -1022                         (no value is subnormal), and
///     d * (E + 1) + 10 <= 1023               (no value overflows),
///
/// where 2^E bounds the largest coordinate and the 10 bits cover the sums of
/// the at most 2^10 terms of a determinant, every non-zero intermediate is a
/// normal double and the arithmetic is as exact as Shewchuk proves.
///
/// The determinants are homogeneous of degree d in the coordinates, so
/// multiplying every coordinate by 2^s multiplies the determinant by
/// 2^(d * s) and keeps its sign. Inputs outside the domain are scaled by the
/// power of two that brings them in, when one exists; it does when the
/// coordinates' binary exponents span less than about 2045 / d orders, far
/// beyond anything ADR-0001's size box and resolution produce. Inputs that no
/// scaling can bring in are a contract violation, reported from Robust.cpp
/// so that it fires in every importer.

export namespace Stereon::Detail {

inline constexpr int MinNormalExponent = -1022;
inline constexpr int MaxFiniteExponent = 1023;
/// log2 of the most terms any determinant sums, with room for the error
/// bound constants.
inline constexpr int TermHeadroom = 10;

/// The exponent of the lowest set bit of a non-zero finite double.
[[nodiscard]] constexpr int LowestBitExponent(const double x) noexcept
{
    const auto bits = std::bit_cast<std::uint64_t>(x) & 0x7FFF'FFFF'FFFF'FFFFu;
    const auto biased = static_cast<int>(bits >> 52);
    const std::uint64_t fraction = bits & 0x000F'FFFF'FFFF'FFFFu;
    if (biased == 0) {
        return -1074 + std::countr_zero(fraction); // subnormal
    }
    return biased - 1075 + std::countr_zero(fraction | (std::uint64_t{1} << 52));
}

/// The exponent of the highest set bit of a non-zero finite double.
[[nodiscard]] constexpr int HighestBitExponent(const double x) noexcept
{
    const auto bits = std::bit_cast<std::uint64_t>(x) & 0x7FFF'FFFF'FFFF'FFFFu;
    const auto biased = static_cast<int>(bits >> 52);
    if (biased == 0) {
        return -1075 + std::bit_width(bits); // subnormal
    }
    return biased - 1023;
}

/// The lowest and highest bit exponents over the non-zero coordinates.
struct ExponentRange {
    int Lowest = INT_MAX;
    int Highest = INT_MIN;

    constexpr void Include(const double x) noexcept
    {
        if (x != 0.0) {
            Lowest = std::min(Lowest, LowestBitExponent(x));
            Highest = std::max(Highest, HighestBitExponent(x));
        }
    }
    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return Lowest == INT_MAX; }
};

/// Why a set of coordinates cannot be evaluated exactly.
enum class DomainStatus : std::uint8_t {
    Inside,    ///< no scaling needed
    Scalable,  ///< inside after multiplying by 2^Shift
    NotFinite, ///< a coordinate is infinite or NaN
    TooWide,   ///< no power of two brings both ends into range
};

struct DomainCheck {
    DomainStatus Status = DomainStatus::Inside;
    int Shift = 0;
};

/// The ceiling of a / b for b > 0.
[[nodiscard]] constexpr int CeilDiv(const int a, const int b) noexcept
{
    return a >= 0 ? (a + b - 1) / b : -((-a) / b);
}

/// Classifies coordinates for a predicate of the given degree.
[[nodiscard]] constexpr DomainCheck CheckDomain(const ExponentRange& range,
                                                const int degree) noexcept
{
    if (range.IsEmpty()) {
        return {}; // every coordinate is zero
    }
    // The shift s must satisfy degree * (Lowest + s) >= MinNormalExponent and
    // degree * (Highest + s + 2) + TermHeadroom <= MaxFiniteExponent.
    const int lowest = CeilDiv(MinNormalExponent, degree) - range.Lowest;
    // Differences of coordinates below 2^(Highest + 1) stay below 2^(Highest + 2).
    const int highest = (MaxFiniteExponent - TermHeadroom) / degree - 2 - range.Highest;
    if (lowest > highest) {
        return {DomainStatus::TooWide, 0};
    }
    if (lowest <= 0 && 0 <= highest) {
        return {};
    }
    return {DomainStatus::Scalable, lowest > 0 ? lowest : highest};
}

/// x * 2^shift, exactly: the domain check guarantees that neither the result
/// nor any intermediate step is subnormal or infinite.
[[nodiscard]] constexpr double ScaledByPowerOfTwo(double x, int shift) noexcept
{
    constexpr int Step = 1000;
    while (shift > Step) {
        x *= 0x1p1000;
        shift -= Step;
    }
    while (shift < -Step) {
        x *= 0x1p-1000;
        shift += Step;
    }
    const auto biased = static_cast<std::uint64_t>(shift + 1023);
    return x * std::bit_cast<double>(biased << 52);
}

[[nodiscard]] constexpr Point2 ScaledByPowerOfTwo(const Point2& p, const int shift) noexcept
{
    return {ScaledByPowerOfTwo(p.X, shift), ScaledByPowerOfTwo(p.Y, shift)};
}

[[nodiscard]] constexpr Point3 ScaledByPowerOfTwo(const Point3& p, const int shift) noexcept
{
    return {ScaledByPowerOfTwo(p.X, shift), ScaledByPowerOfTwo(p.Y, shift),
            ScaledByPowerOfTwo(p.Z, shift)};
}

[[nodiscard]] constexpr bool IsFinite(const double x) noexcept
{
    return x - x == 0.0; // false for infinities and NaN
}

[[nodiscard]] constexpr DomainCheck CheckDomain(std::initializer_list<Point2> points,
                                                const int degree) noexcept
{
    ExponentRange range;
    for (const Point2& p : points) {
        if (!IsFinite(p.X) || !IsFinite(p.Y)) {
            return {DomainStatus::NotFinite, 0};
        }
        range.Include(p.X);
        range.Include(p.Y);
    }
    return CheckDomain(range, degree);
}

[[nodiscard]] constexpr DomainCheck CheckDomain(std::initializer_list<Point3> points,
                                                const int degree) noexcept
{
    ExponentRange range;
    for (const Point3& p : points) {
        if (!IsFinite(p.X) || !IsFinite(p.Y) || !IsFinite(p.Z)) {
            return {DomainStatus::NotFinite, 0};
        }
        range.Include(p.X);
        range.Include(p.Y);
        range.Include(p.Z);
    }
    return CheckDomain(range, degree);
}

/// Reports coordinates outside every scaled domain and raises a contract
/// violation; defined in Robust.cpp, a module implementation unit.
[[noreturn]] void PredicateDomainViolation(const char* predicate, DomainStatus status) noexcept;

/// The absolute part of the stage-A error bounds: the total roundoff of the
/// few operations of stage A that could have produced subnormal results,
/// each at most 2^-1075, is far below 2^-1022. With it, stage A is sound for
/// any finite input, and inputs whose determinant is tiny go on to the
/// adaptive stages, which check the domain.
inline constexpr double UnderflowGuard = 0x1p-1022;

} // namespace Stereon::Detail
