// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

module Stereon.Robust;

import :Dynamic;
import :Expansion;

namespace Stereon {

// --- DynamicExpansion ---------------------------------------------------

DynamicExpansion::DynamicExpansion(const double value)
{
    PushBackNonZero(value);
}

Sign DynamicExpansion::GetSign() const noexcept
{
    return m_Components.empty() ? Sign::Zero : SignOf(m_Components.back());
}

double DynamicExpansion::GetEstimate() const noexcept
{
    double sum = 0.0;
    for (const double c : m_Components) {
        sum += c;
    }
    return sum;
}

void DynamicExpansion::PushBackNonZero(const double c)
{
    if (c != 0.0) {
        m_Components.push_back(c);
    }
}

/// e + f (Shewchuk's fast_expansion_sum_zeroelim), as Expansion<N>'s Sum.
DynamicExpansion Sum(const DynamicExpansion& e, const DynamicExpansion& f)
{
    if (e.IsZero()) {
        return f;
    }
    if (f.IsZero()) {
        return e;
    }
    const std::span<const double> ec = e.GetComponents();
    const std::span<const double> fc = f.GetComponents();
    std::size_t ei = 0;
    std::size_t fi = 0;
    double q = Detail::Magnitude(fc[0]) > Detail::Magnitude(ec[0]) ? ec[ei++] : fc[fi++];
    bool first = true;
    DynamicExpansion h;
    h.m_Components.reserve(ec.size() + fc.size());
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

DynamicExpansion Negated(const DynamicExpansion& e)
{
    DynamicExpansion h;
    h.m_Components.reserve(e.GetLength());
    for (const double c : e.GetComponents()) {
        h.m_Components.push_back(-c);
    }
    return h;
}

DynamicExpansion Difference(const DynamicExpansion& e, const DynamicExpansion& f)
{
    return Sum(e, Negated(f));
}

/// e * b (Shewchuk's scale_expansion_zeroelim), as Expansion<N>'s Scaled.
DynamicExpansion Scaled(const DynamicExpansion& e, const double b)
{
    DynamicExpansion h;
    if (e.IsZero() || b == 0.0) {
        return h;
    }
    const Halves bHalves = Split(b);
    const std::span<const double> ec = e.GetComponents();
    h.m_Components.reserve(2 * ec.size());
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

/// e * f: the longer scaled by each component of the shorter, summed.
DynamicExpansion Product(const DynamicExpansion& e, const DynamicExpansion& f)
{
    DynamicExpansion h;
    if (e.IsZero() || f.IsZero()) {
        return h;
    }
    const DynamicExpansion& longer = f.GetLength() <= e.GetLength() ? e : f;
    const DynamicExpansion& shorter = f.GetLength() <= e.GetLength() ? f : e;
    for (const double c : shorter.GetComponents()) {
        h = Sum(h, Scaled(longer, c));
    }
    return h;
}

DynamicExpansion operator+(const DynamicExpansion& e, const DynamicExpansion& f)
{
    return Sum(e, f);
}
DynamicExpansion operator-(const DynamicExpansion& e, const DynamicExpansion& f)
{
    return Difference(e, f);
}
DynamicExpansion operator-(const DynamicExpansion& e)
{
    return Negated(e);
}
DynamicExpansion operator*(const DynamicExpansion& e, const DynamicExpansion& f)
{
    return Product(e, f);
}
DynamicExpansion operator*(const DynamicExpansion& e, const double b)
{
    return Scaled(e, b);
}

// --- ExactPolynomial ----------------------------------------------------

ExactPolynomial::ExactPolynomial(std::vector<DynamicExpansion> coefficients)
    : m_Coefficients(std::move(coefficients))
{
    while (!m_Coefficients.empty() && m_Coefficients.back().IsZero()) {
        m_Coefficients.pop_back();
    }
}

Sign ExactPolynomial::SignAt(const double x) const
{
    DynamicExpansion value;
    for (std::size_t k = m_Coefficients.size(); k-- > 0;) {
        value = Sum(Scaled(value, x), m_Coefficients[k]);
    }
    return value.GetSign();
}

Sign ExactPolynomial::SignAtInfinity() const noexcept
{
    return IsZero() ? Sign::Zero : m_Coefficients.back().GetSign();
}

Sign ExactPolynomial::SignAtMinusInfinity() const noexcept
{
    const Sign leading = SignAtInfinity();
    return GetDegree() % 2 == 0 ? leading : -leading;
}

ExactPolynomial Sum(const ExactPolynomial& p, const ExactPolynomial& q)
{
    const std::size_t n = std::max(p.m_Coefficients.size(), q.m_Coefficients.size());
    std::vector<DynamicExpansion> c(n);
    for (std::size_t k = 0; k < n; ++k) {
        if (k < p.m_Coefficients.size()) {
            c[k] = Sum(c[k], p.m_Coefficients[k]);
        }
        if (k < q.m_Coefficients.size()) {
            c[k] = Sum(c[k], q.m_Coefficients[k]);
        }
    }
    return ExactPolynomial(std::move(c));
}

ExactPolynomial Difference(const ExactPolynomial& p, const ExactPolynomial& q)
{
    return Sum(p, Scaled(q, DynamicExpansion{-1.0}));
}

ExactPolynomial Product(const ExactPolynomial& p, const ExactPolynomial& q)
{
    if (p.IsZero() || q.IsZero()) {
        return {};
    }
    std::vector<DynamicExpansion> c(p.m_Coefficients.size() + q.m_Coefficients.size() - 1);
    for (std::size_t i = 0; i < p.m_Coefficients.size(); ++i) {
        for (std::size_t j = 0; j < q.m_Coefficients.size(); ++j) {
            c[i + j] = Sum(c[i + j], Product(p.m_Coefficients[i], q.m_Coefficients[j]));
        }
    }
    return ExactPolynomial(std::move(c));
}

ExactPolynomial Scaled(const ExactPolynomial& p, const DynamicExpansion& c)
{
    std::vector<DynamicExpansion> scaled;
    scaled.reserve(p.m_Coefficients.size());
    for (const DynamicExpansion& coefficient : p.m_Coefficients) {
        scaled.push_back(Product(coefficient, c));
    }
    return ExactPolynomial(std::move(scaled));
}

ExactPolynomial Derivative(const ExactPolynomial& p)
{
    if (p.GetDegree() < 1) {
        return {};
    }
    std::vector<DynamicExpansion> d;
    d.reserve(p.m_Coefficients.size() - 1);
    for (std::size_t k = 1; k < p.m_Coefficients.size(); ++k) {
        d.push_back(Scaled(p.m_Coefficients[k], static_cast<double>(k)));
    }
    return ExactPolynomial(std::move(d));
}

// --- SturmSequence ------------------------------------------------------

namespace {

/// The pseudo-remainder of a by b with the sign of the true remainder: a is
/// multiplied by b's leading coefficient once per reduction step, and by it
/// once more when that leaves an odd power of a negative coefficient.
[[nodiscard]] ExactPolynomial SignedPseudoRemainder(ExactPolynomial a, const ExactPolynomial& b)
{
    const DynamicExpansion& lead = b[static_cast<std::size_t>(b.GetDegree())];
    int scalings = 0;
    while (!a.IsZero() && a.GetDegree() >= b.GetDegree()) {
        const auto shift = static_cast<std::size_t>(a.GetDegree() - b.GetDegree());
        std::vector<DynamicExpansion> shifted(shift);
        for (const DynamicExpansion& c : b.GetCoefficients()) {
            shifted.push_back(Product(c, a[static_cast<std::size_t>(a.GetDegree())]));
        }
        a = Difference(Scaled(a, lead), ExactPolynomial(std::move(shifted)));
        ++scalings;
    }
    if (scalings % 2 == 1 && lead.GetSign() == Sign::Negative) {
        a = Scaled(a, lead);
    }
    return a;
}

} // namespace

SturmSequence::SturmSequence(ExactPolynomial p)
{
    m_Members.push_back(std::move(p));
    if (m_Members.front().GetDegree() < 1) {
        return;
    }
    m_Members.push_back(Derivative(m_Members.front()));
    while (m_Members.back().GetDegree() > 0) {
        ExactPolynomial next =
            Scaled(SignedPseudoRemainder(m_Members[m_Members.size() - 2], m_Members.back()),
                   DynamicExpansion{-1.0});
        if (next.IsZero()) {
            break;
        }
        m_Members.push_back(std::move(next));
    }
}

int SturmSequence::VariationsAt(const double x) const
{
    int variations = 0;
    Sign previous = Sign::Zero;
    for (const ExactPolynomial& member : m_Members) {
        const Sign sign = std::isinf(x)
                              ? (x > 0.0 ? member.SignAtInfinity() : member.SignAtMinusInfinity())
                              : member.SignAt(x);
        if (sign == Sign::Zero) {
            continue;
        }
        if (previous != Sign::Zero && sign != previous) {
            ++variations;
        }
        previous = sign;
    }
    return variations;
}

int SturmSequence::CountRoots(const double a, const double b) const
{
    return VariationsAt(a) - VariationsAt(b);
}

} // namespace Stereon
