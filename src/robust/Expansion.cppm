// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

export module Stereon.Robust:Expansion;

import :Sign;

/// Expansion arithmetic (ADR-0002), after Shewchuk, "Adaptive Precision
/// Floating-Point Arithmetic and Fast Robust Geometric Predicates" (1997).
///
/// An expansion is a sum of doubles that represents a real number exactly:
/// the components do not overlap (the least significant bit of each is more
/// significant than the most significant bit of the one before it), they are
/// stored in increasing magnitude, and none of them is zero. The sign of an
/// expansion is therefore the sign of its last component, and the exact
/// value of a determinant comes out of ordinary double arithmetic, with the
/// rounding errors kept instead of dropped.
///
/// Three error-free transformations supply the components: TwoSum (Knuth),
/// FastTwoSum (Dekker) and TwoProduct (Dekker, by splitting each factor into
/// two 26-bit halves). Sum, Scaled and Product combine expansions, each
/// returning an expansion whose capacity is the worst case of the operation,
/// so that a predicate's intermediate sizes follow from its formula at
/// compile time. Every operation is constexpr and has no state.
///
/// The arithmetic assumes round-to-nearest-even, no fused multiply-add (the
/// robust library is built with -ffp-contract=off, and so is every unit that
/// instantiates these templates) and no overflow. Underflow is harmless: the
/// transformations stay exact for subnormal results.
///
/// The manual's "Exact predicates" chapter derives the algorithms.

export namespace Stereon {

/// 2^27 + 1: Split() multiplies by it to cut a 53-bit significand into two
/// 26-bit halves whose product is exact.
inline constexpr double Splitter = 0x1p27 + 1.0;

/// A double as the sum of two non-overlapping halves of at most 26 bits.
struct Halves {
    double Hi;
    double Lo;
};

[[nodiscard]] constexpr Halves Split(const double a) noexcept
{
    const double c = Splitter * a;
    const double big = c - a;
    const double hi = c - big;
    return {hi, a - hi};
}

} // namespace Stereon

namespace Stereon::Detail {

/// |x|, by clearing the sign bit. std::fabs is not a constant expression on
/// every target (MinGW's math.h supplies an inline C fabs at -O3), and GCC
/// 16.2 rejects a contract check in the same constant evaluation as an
/// `if consteval`, so this is the one spelling that is constexpr everywhere
/// and one instruction at run time.
[[nodiscard]] constexpr double Magnitude(const double x) noexcept
{
    return std::bit_cast<double>(std::bit_cast<std::uint64_t>(x) & 0x7FFFFFFFFFFFFFFFu);
}

/// A rounded result with its exact roundoff: Value + Error is the true sum
/// or product.
struct Terms {
    double Value;
    double Error;
};

/// Requires |a| >= |b| (Dekker). Three operations.
[[nodiscard]] constexpr Terms FastTwoSumTerms(const double a, const double b) noexcept
{
    const double x = a + b;
    const double bVirtual = x - a;
    return {x, b - bVirtual};
}

/// Any a and b (Knuth). Six operations.
[[nodiscard]] constexpr Terms TwoSumTerms(const double a, const double b) noexcept
{
    const double x = a + b;
    const double bVirtual = x - a;
    const double aVirtual = x - bVirtual;
    const double bRound = b - bVirtual;
    const double aRound = a - aVirtual;
    return {x, aRound + bRound};
}

[[nodiscard]] constexpr Terms TwoDiffTerms(const double a, const double b) noexcept
{
    const double x = a - b;
    const double bVirtual = a - x;
    const double aVirtual = x + bVirtual;
    const double bRound = bVirtual - b;
    const double aRound = a - aVirtual;
    return {x, aRound + bRound};
}

/// a * b with b already split, for scaling a whole expansion by b.
[[nodiscard]] constexpr Terms TwoProductTerms(const double a, const double b,
                                              const Halves bHalves) noexcept
{
    const double x = a * b;
    const Halves aHalves = Split(a);
    const double error1 = x - (aHalves.Hi * bHalves.Hi);
    const double error2 = error1 - (aHalves.Lo * bHalves.Hi);
    const double error3 = error2 - (aHalves.Hi * bHalves.Lo);
    return {x, (aHalves.Lo * bHalves.Lo) - error3};
}

} // namespace Stereon::Detail

export namespace Stereon {

/// A real number as a sum of at most N doubles; see the module comment.
/// N is a capacity: the length is whatever the arithmetic needed.
template <std::size_t N>
    requires(N > 0)
class Expansion {
public:
    static constexpr std::size_t Capacity = N;

    /// Zero.
    constexpr Expansion() = default;

    /// The expansion of one double; zero has no components.
    constexpr explicit Expansion(const double value) noexcept
    {
        if (value != 0.0) {
            m_Components[0] = value;
            m_Length = 1;
        }
    }

    /// A smaller expansion as a larger one.
    template <std::size_t M>
        requires(M < N)
    constexpr Expansion(const Expansion<M>& e) noexcept // NOLINT(google-explicit-constructor)
        : m_Length(e.GetLength())
    {
        for (std::size_t i = 0; i < m_Length; ++i) {
            m_Components[i] = e[i];
        }
    }

    /// A larger expansion as a smaller one, where the formula bounds the
    /// length more tightly than the capacity arithmetic does.
    template <std::size_t M>
        requires(M > N)
    [[nodiscard]] static constexpr Expansion Narrowed(const Expansion<M>& e)
        STEREON_HOT_PRE(e.GetLength() <= N)
    {
        Expansion h;
        h.m_Length = e.GetLength();
        for (std::size_t i = 0; i < h.m_Length; ++i) {
            h.m_Components[i] = e[i];
        }
        return h;
    }

    [[nodiscard]] constexpr std::size_t GetLength() const noexcept { return m_Length; }
    [[nodiscard]] constexpr bool IsZero() const noexcept { return m_Length == 0; }

    /// The components, smallest magnitude first, none of them zero.
    [[nodiscard]] constexpr std::span<const double> GetComponents() const noexcept
    {
        return std::span<const double>{m_Components.data(), m_Length};
    }
    [[nodiscard]] constexpr double operator[](const std::size_t i) const
        STEREON_HOT_PRE(i < GetLength())
    {
        return m_Components[i];
    }

    /// The exact sign of the value: the sign of the last component.
    [[nodiscard]] constexpr Sign GetSign() const noexcept
    {
        return m_Length == 0 ? Sign::Zero : SignOf(m_Components[m_Length - 1]);
    }

    /// The last component, or zero. It is within a factor of two of the
    /// value: the components below it sum to less than its least significant
    /// bit, which is why its sign is the sign of the value.
    [[nodiscard]] constexpr double GetMostSignificant() const noexcept
    {
        return m_Length == 0 ? 0.0 : m_Components[m_Length - 1];
    }

    /// The components summed in floating point, smallest first (Shewchuk's
    /// estimate): the value with a relative error of a few roundings.
    [[nodiscard]] constexpr double GetEstimate() const noexcept
    {
        double sum = 0.0;
        for (std::size_t i = 0; i < m_Length; ++i) {
            sum += m_Components[i];
        }
        return sum;
    }

    /// Appends a component; a zero is dropped. The arithmetic below builds
    /// expansions with this, in increasing magnitude.
    constexpr void PushBackNonZero(const double c) STEREON_HOT_PRE(c == 0.0 || GetLength() < N)
    {
        if (c != 0.0) {
            m_Components[m_Length++] = c;
        }
    }

private:
    std::array<double, N> m_Components{};
    std::size_t m_Length = 0;
};

// --- Error-free transformations ---------------------------------------------

/// a + b exactly, as [roundoff, sum].
[[nodiscard]] constexpr Expansion<2> TwoSum(const double a, const double b) noexcept
{
    const Detail::Terms t = Detail::TwoSumTerms(a, b);
    Expansion<2> h;
    h.PushBackNonZero(t.Error);
    h.PushBackNonZero(t.Value);
    return h;
}

/// a + b exactly when |a| >= |b|, three operations instead of six.
[[nodiscard]] constexpr Expansion<2> FastTwoSum(const double a, const double b)
    STEREON_HOT_PRE(Detail::Magnitude(a) >= Detail::Magnitude(b))
{
    const Detail::Terms t = Detail::FastTwoSumTerms(a, b);
    Expansion<2> h;
    h.PushBackNonZero(t.Error);
    h.PushBackNonZero(t.Value);
    return h;
}

/// a - b exactly.
[[nodiscard]] constexpr Expansion<2> TwoDiff(const double a, const double b) noexcept
{
    const Detail::Terms t = Detail::TwoDiffTerms(a, b);
    Expansion<2> h;
    h.PushBackNonZero(t.Error);
    h.PushBackNonZero(t.Value);
    return h;
}

/// The roundoff of x = fl(a - b): a - b == x + TwoDiffTail(a, b, x).
[[nodiscard]] constexpr double TwoDiffTail(const double a, const double b, const double x) noexcept
{
    const double bVirtual = a - x;
    const double aVirtual = x + bVirtual;
    const double bRound = bVirtual - b;
    const double aRound = a - aVirtual;
    return aRound + bRound;
}

/// a * b exactly, as [roundoff, product].
[[nodiscard]] constexpr Expansion<2> TwoProduct(const double a, const double b) noexcept
{
    const Detail::Terms t = Detail::TwoProductTerms(a, b, Split(b));
    Expansion<2> h;
    h.PushBackNonZero(t.Error);
    h.PushBackNonZero(t.Value);
    return h;
}

// --- Expansion arithmetic ---------------------------------------------------

/// e + f (Shewchuk's fast_expansion_sum_zeroelim): the components of both
/// are merged by magnitude and summed with a running TwoSum, so that every
/// roundoff becomes a component of the result.
template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<A + B> Sum(const Expansion<A>& e, const Expansion<B>& f) noexcept
{
    if (e.IsZero()) {
        return Expansion<A + B>{f};
    }
    if (f.IsZero()) {
        return Expansion<A + B>{e};
    }
    const std::span<const double> ec = e.GetComponents();
    const std::span<const double> fc = f.GetComponents();
    std::size_t ei = 0;
    std::size_t fi = 0;
    // The running sum starts as the smallest component; the first merge may
    // use the cheaper FastTwoSum because the next component is at least as
    // large.
    double q = Detail::Magnitude(fc[0]) > Detail::Magnitude(ec[0]) ? ec[ei++] : fc[fi++];
    bool first = true;

    Expansion<A + B> h;
    const auto absorb = [&](const double next) {
        const Detail::Terms t =
            first ? Detail::FastTwoSumTerms(next, q) : Detail::TwoSumTerms(q, next);
        first = false;
        h.PushBackNonZero(t.Error);
        q = t.Value;
    };
    while (ei < ec.size() && fi < fc.size()) {
        absorb(Detail::Magnitude(fc[fi]) > Detail::Magnitude(ec[ei]) ? ec[ei++] : fc[fi++]);
    }
    while (ei < ec.size()) {
        absorb(ec[ei++]);
    }
    while (fi < fc.size()) {
        absorb(fc[fi++]);
    }
    h.PushBackNonZero(q);
    return h;
}

template <std::size_t A>
[[nodiscard]] constexpr Expansion<A> Negated(const Expansion<A>& e) noexcept
{
    Expansion<A> h;
    for (const double c : e.GetComponents()) {
        h.PushBackNonZero(-c);
    }
    return h;
}

/// e - f.
template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<A + B> Difference(const Expansion<A>& e,
                                                    const Expansion<B>& f) noexcept
{
    return Sum(e, Negated(f));
}

/// e * b (Shewchuk's scale_expansion_zeroelim): each component times b as
/// an exact TwoProduct, the pieces summed in order of magnitude.
template <std::size_t A>
[[nodiscard]] constexpr Expansion<2 * A> Scaled(const Expansion<A>& e, const double b) noexcept
{
    Expansion<2 * A> h;
    if (e.IsZero() || b == 0.0) {
        return h;
    }
    const Halves bHalves = Split(b);
    const std::span<const double> ec = e.GetComponents();
    Detail::Terms head = Detail::TwoProductTerms(ec[0], b, bHalves);
    h.PushBackNonZero(head.Error);
    double q = head.Value;
    for (std::size_t i = 1; i < ec.size(); ++i) {
        const Detail::Terms product = Detail::TwoProductTerms(ec[i], b, bHalves);
        const Detail::Terms sum = Detail::TwoSumTerms(q, product.Error);
        h.PushBackNonZero(sum.Error);
        head = Detail::FastTwoSumTerms(product.Value, sum.Value);
        h.PushBackNonZero(head.Error);
        q = head.Value;
    }
    h.PushBackNonZero(q);
    return h;
}

/// e * f: the longer expansion scaled by each component of the shorter one,
/// the pieces summed. The capacity is the worst case 2 * A * B; the length
/// is what the operands needed.
template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<2 * A * B> Product(const Expansion<A>& e,
                                                     const Expansion<B>& f) noexcept
{
    Expansion<2 * A * B> h;
    if (e.IsZero() || f.IsZero()) {
        return h;
    }
    if (f.GetLength() <= e.GetLength()) {
        for (const double c : f.GetComponents()) {
            h = Expansion<2 * A * B>::Narrowed(Sum(h, Scaled(e, c)));
        }
    } else {
        for (const double c : e.GetComponents()) {
            h = Expansion<2 * A * B>::Narrowed(Sum(h, Scaled(f, c)));
        }
    }
    return h;
}

template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<A + B> operator+(const Expansion<A>& e,
                                                   const Expansion<B>& f) noexcept
{
    return Sum(e, f);
}
template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<A + B> operator-(const Expansion<A>& e,
                                                   const Expansion<B>& f) noexcept
{
    return Difference(e, f);
}
template <std::size_t A>
[[nodiscard]] constexpr Expansion<A> operator-(const Expansion<A>& e) noexcept
{
    return Negated(e);
}
template <std::size_t A>
[[nodiscard]] constexpr Expansion<2 * A> operator*(const Expansion<A>& e, const double b) noexcept
{
    return Scaled(e, b);
}
template <std::size_t A>
[[nodiscard]] constexpr Expansion<2 * A> operator*(const double b, const Expansion<A>& e) noexcept
{
    return Scaled(e, b);
}
template <std::size_t A, std::size_t B>
[[nodiscard]] constexpr Expansion<2 * A * B> operator*(const Expansion<A>& e,
                                                       const Expansion<B>& f) noexcept
{
    return Product(e, f);
}

} // namespace Stereon
