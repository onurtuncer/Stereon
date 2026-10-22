// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <span>
#include <vector>

export module Stereon.Robust:Dynamic;

import :Expansion;
import :Sign;

/// Exact arithmetic without a capacity (ADR-0002): DynamicExpansion is
/// Expansion<N> on the heap, for polynomials in the input doubles whose
/// degree no compile-time bound should limit, and ExactPolynomial is a
/// polynomial in one variable with such coefficients, with the Sturm
/// sequence that counts its real roots in an interval exactly. The
/// analytic intersections decide the configurations of degree 10 and more
/// with them (the distance between skew axes against the radii, the
/// branch points of two cylinders).

export namespace Stereon {

/// A real number as a sum of doubles of any length: the components do not
/// overlap, increase in magnitude and are never zero, as in Expansion<N>,
/// so the sign is the last component's. Every operation is exact; the
/// same round-to-nearest, no-FMA, no-underflow assumptions apply.
class DynamicExpansion {
public:
    /// Zero.
    DynamicExpansion() = default;
    /// The expansion of one double.
    explicit DynamicExpansion(double value);
    /// A fixed-capacity expansion.
    template <std::size_t N>
    DynamicExpansion(const Expansion<N>& e) // NOLINT(google-explicit-constructor)
        : m_Components(e.GetComponents().begin(), e.GetComponents().end())
    {}

    [[nodiscard]] std::size_t GetLength() const noexcept { return m_Components.size(); }
    [[nodiscard]] bool IsZero() const noexcept { return m_Components.empty(); }
    [[nodiscard]] std::span<const double> GetComponents() const noexcept { return m_Components; }
    /// The exact sign: the sign of the last component.
    [[nodiscard]] Sign GetSign() const noexcept;
    /// The components summed in floating point, smallest first.
    [[nodiscard]] double GetEstimate() const noexcept;

    /// Appends a component; a zero is dropped.
    void PushBackNonZero(double c);

    friend DynamicExpansion Sum(const DynamicExpansion& e, const DynamicExpansion& f);
    friend DynamicExpansion Negated(const DynamicExpansion& e);
    friend DynamicExpansion Difference(const DynamicExpansion& e, const DynamicExpansion& f);
    friend DynamicExpansion Scaled(const DynamicExpansion& e, double b);
    friend DynamicExpansion Product(const DynamicExpansion& e, const DynamicExpansion& f);

private:
    std::vector<double> m_Components;
};

[[nodiscard]] DynamicExpansion operator+(const DynamicExpansion& e, const DynamicExpansion& f);
[[nodiscard]] DynamicExpansion operator-(const DynamicExpansion& e, const DynamicExpansion& f);
[[nodiscard]] DynamicExpansion operator-(const DynamicExpansion& e);
[[nodiscard]] DynamicExpansion operator*(const DynamicExpansion& e, const DynamicExpansion& f);
[[nodiscard]] DynamicExpansion operator*(const DynamicExpansion& e, double b);

/// A polynomial in one variable with exact coefficients, lowest degree
/// first; the leading coefficient is never zero, and the zero polynomial
/// has no coefficients and degree -1.
class ExactPolynomial {
public:
    ExactPolynomial() = default;
    /// From coefficients, lowest degree first; leading zeros are dropped.
    explicit ExactPolynomial(std::vector<DynamicExpansion> coefficients);

    [[nodiscard]] int GetDegree() const noexcept
    {
        return static_cast<int>(m_Coefficients.size()) - 1;
    }
    [[nodiscard]] bool IsZero() const noexcept { return m_Coefficients.empty(); }
    [[nodiscard]] const DynamicExpansion& operator[](std::size_t k) const
    {
        return m_Coefficients[k];
    }
    [[nodiscard]] std::span<const DynamicExpansion> GetCoefficients() const noexcept
    {
        return m_Coefficients;
    }

    /// The exact sign of the value at x, by Horner's rule.
    [[nodiscard]] Sign SignAt(double x) const;
    /// The sign as x grows without bound: the leading coefficient's.
    [[nodiscard]] Sign SignAtInfinity() const noexcept;
    /// The sign as x falls without bound.
    [[nodiscard]] Sign SignAtMinusInfinity() const noexcept;

    friend ExactPolynomial Sum(const ExactPolynomial& p, const ExactPolynomial& q);
    friend ExactPolynomial Difference(const ExactPolynomial& p, const ExactPolynomial& q);
    friend ExactPolynomial Product(const ExactPolynomial& p, const ExactPolynomial& q);
    friend ExactPolynomial Scaled(const ExactPolynomial& p, const DynamicExpansion& c);
    friend ExactPolynomial Derivative(const ExactPolynomial& p);

private:
    std::vector<DynamicExpansion> m_Coefficients;
};

/// The Sturm sequence of a polynomial, by pseudo-remainders, so that every
/// member keeps exact coefficients: the number of distinct real roots in
/// (a, b] is the difference of the sign variations at a and at b. Exact
/// for a polynomial with simple roots; one with a multiple root counts it
/// once, which the callers exclude first.
class SturmSequence {
public:
    explicit SturmSequence(ExactPolynomial p);

    /// The sign variations of the sequence at x, or at +-infinity for an
    /// infinite x.
    [[nodiscard]] int VariationsAt(double x) const;
    /// The number of distinct real roots in (a, b], a < b, either infinite.
    [[nodiscard]] int CountRoots(double a, double b) const;
    [[nodiscard]] const ExactPolynomial& GetPolynomial() const noexcept
    {
        return m_Members.front();
    }

private:
    std::vector<ExactPolynomial> m_Members;
};

} // namespace Stereon
