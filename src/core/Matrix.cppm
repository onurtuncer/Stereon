// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <format>

export module Stereon.Core:Matrix;

import :Scalar;
import :Vector;

/// Small square matrices and affine transforms, scalar-generic like the
/// vectors (ADR-0021). Mat2, Mat3, Mat4 and Transform are the double
/// instances. Everything is constexpr and closed-form (cofactor expansion),
/// so results do not depend on a pivoting strategy (ADR-0007).
export namespace Stereon {

/// N×N matrix stored row-major: BasicMat<T, 2>{a, b, c, d} is [[a, b], [c, d]].
/// Default-constructed it is the zero matrix.
template <Scalar T, std::size_t N>
    requires(N >= 1 && N <= 4)
struct BasicMat {
    std::array<T, N * N> Elements{};

    [[nodiscard]] static constexpr BasicMat Identity()
    {
        BasicMat m;
        for (std::size_t i = 0; i < N; ++i) {
            m(i, i) = T{1};
        }
        return m;
    }

    [[nodiscard]] static constexpr BasicMat Diagonal(const std::array<T, N>& diagonal)
    {
        BasicMat m;
        for (std::size_t i = 0; i < N; ++i) {
            m(i, i) = diagonal[i];
        }
        return m;
    }

    constexpr bool operator==(const BasicMat&) const = default;

    constexpr T& operator()(const std::size_t row, const std::size_t col)
        STEREON_HOT_PRE(row < N && col < N)
    {
        return Elements[row * N + col];
    }
    constexpr const T& operator()(const std::size_t row, const std::size_t col) const
        STEREON_HOT_PRE(row < N && col < N)
    {
        return Elements[row * N + col];
    }

    constexpr BasicMat& operator+=(const BasicMat& o)
    {
        for (std::size_t i = 0; i < N * N; ++i) {
            Elements[i] = Elements[i] + o.Elements[i];
        }
        return *this;
    }
    constexpr BasicMat& operator-=(const BasicMat& o)
    {
        for (std::size_t i = 0; i < N * N; ++i) {
            Elements[i] = Elements[i] - o.Elements[i];
        }
        return *this;
    }
    constexpr BasicMat& operator*=(const T& s)
    {
        for (T& e : Elements) {
            e = e * s;
        }
        return *this;
    }

    friend constexpr BasicMat operator+(BasicMat a, const BasicMat& b) { return a += b; }
    friend constexpr BasicMat operator-(BasicMat a, const BasicMat& b) { return a -= b; }
    friend constexpr BasicMat operator-(const BasicMat& a) { return a * T{-1}; }
    friend constexpr BasicMat operator*(BasicMat a, const T& s) { return a *= s; }
    friend constexpr BasicMat operator*(const T& s, BasicMat a) { return a *= s; }

    friend constexpr BasicMat operator*(const BasicMat& a, const BasicMat& b)
    {
        BasicMat r;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                T sum{};
                for (std::size_t k = 0; k < N; ++k) {
                    sum = sum + a(i, k) * b(k, j);
                }
                r(i, j) = sum;
            }
        }
        return r;
    }
};

template <Scalar T>
using BasicMat2 = BasicMat<T, 2>;
template <Scalar T>
using BasicMat3 = BasicMat<T, 3>;
template <Scalar T>
using BasicMat4 = BasicMat<T, 4>;

using Mat2 = BasicMat<double, 2>;
using Mat3 = BasicMat<double, 3>;
using Mat4 = BasicMat<double, 4>;

template <Scalar T>
[[nodiscard]] constexpr BasicVec2<T> operator*(const BasicMat<T, 2>& m, const BasicVec2<T>& v)
{
    return {m(0, 0) * v.X + m(0, 1) * v.Y, m(1, 0) * v.X + m(1, 1) * v.Y};
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> operator*(const BasicMat<T, 3>& m, const BasicVec3<T>& v)
{
    return {m(0, 0) * v.X + m(0, 1) * v.Y + m(0, 2) * v.Z,
            m(1, 0) * v.X + m(1, 1) * v.Y + m(1, 2) * v.Z,
            m(2, 0) * v.X + m(2, 1) * v.Y + m(2, 2) * v.Z};
}

template <Scalar T, std::size_t N>
[[nodiscard]] constexpr BasicMat<T, N> Transposed(const BasicMat<T, N>& m)
{
    BasicMat<T, N> r;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            r(j, i) = m(i, j);
        }
    }
    return r;
}

template <Scalar T, std::size_t N>
[[nodiscard]] constexpr T Trace(const BasicMat<T, N>& m)
{
    T sum{};
    for (std::size_t i = 0; i < N; ++i) {
        sum = sum + m(i, i);
    }
    return sum;
}

} // namespace Stereon

namespace Stereon::Detail {

/// m without row `row` and column `col`.
template <Scalar T, std::size_t N>
    requires(N > 1)
[[nodiscard]] constexpr BasicMat<T, N - 1> MinorOf(const BasicMat<T, N>& m, const std::size_t row,
                                                   const std::size_t col)
{
    BasicMat<T, N - 1> r;
    std::size_t ri = 0;
    for (std::size_t i = 0; i < N; ++i) {
        if (i == row) {
            continue;
        }
        std::size_t rj = 0;
        for (std::size_t j = 0; j < N; ++j) {
            if (j == col) {
                continue;
            }
            r(ri, rj) = m(i, j);
            ++rj;
        }
        ++ri;
    }
    return r;
}

} // namespace Stereon::Detail

export namespace Stereon {

template <Scalar T, std::size_t N>
[[nodiscard]] constexpr T Determinant(const BasicMat<T, N>& m)
{
    if constexpr (N == 1) {
        return m(0, 0);
    } else if constexpr (N == 2) {
        return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
    } else {
        // Expansion along the first row.
        T det{};
        for (std::size_t j = 0; j < N; ++j) {
            const T term = m(0, j) * Determinant(Detail::MinorOf(m, 0, j));
            det = (j % 2 == 0) ? det + term : det - term;
        }
        return det;
    }
}

/// The signed minor that Inverse() is built from; exposed because robust code
/// sometimes wants the adjugate without the division.
template <Scalar T, std::size_t N>
[[nodiscard]] constexpr T Cofactor(const BasicMat<T, N>& m, const std::size_t row,
                                   const std::size_t col) STEREON_HOT_PRE(row < N && col < N)
{
    if constexpr (N == 1) {
        return T{1};
    } else {
        const T minor = Determinant(Detail::MinorOf(m, row, col));
        return ((row + col) % 2 == 0) ? minor : -minor;
    }
}

/// Inverse by the adjugate. A singular matrix is a bug in the caller: code
/// that may meet one checks Determinant() first and reports a KernelError.
template <Scalar T, std::size_t N>
[[nodiscard]] constexpr BasicMat<T, N> Inverse(const BasicMat<T, N>& m)
    STEREON_HOT_PRE(Determinant(m) != T{})
{
    const T det = Determinant(m);
    BasicMat<T, N> r;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            r(i, j) = Cofactor(m, j, i) / det;
        }
    }
    return r;
}

/// An affine map of 3D space: x -> Linear * x + Translation. Default-constructed
/// it is the identity.
template <Scalar T>
struct BasicTransform {
    BasicMat<T, 3> Linear = BasicMat<T, 3>::Identity();
    BasicVec3<T> Translation{};

    [[nodiscard]] static constexpr BasicTransform Identity() { return {}; }

    constexpr bool operator==(const BasicTransform&) const = default;

    /// Composition: (a * b) applied to p is a applied to (b applied to p).
    friend constexpr BasicTransform operator*(const BasicTransform& a, const BasicTransform& b)
    {
        return {a.Linear * b.Linear, a.Linear * b.Translation + a.Translation};
    }
};

using Transform = BasicTransform<double>;

template <Scalar T>
[[nodiscard]] constexpr BasicTransform<T> Inverse(const BasicTransform<T>& t)
    STEREON_HOT_PRE(Determinant(t.Linear) != T{})
{
    const BasicMat<T, 3> inverse = Inverse(t.Linear);
    return {inverse, -(inverse * t.Translation)};
}

/// Points move with the whole map, vectors only with its linear part.
template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Transformed(const BasicPoint3<T>& p,
                                                   const BasicTransform<T>& t)
{
    return ToPoint(t.Linear * ToVec(p) + t.Translation);
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Transformed(const BasicVec3<T>& v, const BasicTransform<T>& t)
{
    return t.Linear * v;
}

template <Scalar T>
[[nodiscard]] constexpr BasicTransform<T> MakeTranslation(const BasicVec3<T>& offset)
{
    return {BasicMat<T, 3>::Identity(), offset};
}

template <Scalar T>
[[nodiscard]] constexpr BasicTransform<T> MakeScaling(const T& factor)
{
    return {BasicMat<T, 3>::Diagonal({factor, factor, factor}), {}};
}

template <Scalar T>
[[nodiscard]] constexpr BasicTransform<T> MakeScaling(const BasicVec3<T>& factors)
{
    return {BasicMat<T, 3>::Diagonal({factors.X, factors.Y, factors.Z}), {}};
}

/// Rotation by `angle` radians about `axis` through the origin, right-handed
/// (Rodrigues' formula). The axis need not be unit length, only non-zero.
template <Scalar T>
[[nodiscard]] constexpr BasicTransform<T> MakeRotation(const BasicVec3<T>& axis, const T& angle)
    STEREON_HOT_PRE(T{} < LengthSquared(axis))
{
    using std::cos;
    using std::sin;
    const BasicVec3<T> k = Normalized(axis);
    const T c = cos(angle);
    const T s = sin(angle);
    const T t = T{1} - c;
    BasicMat<T, 3> r;
    r(0, 0) = c + k.X * k.X * t;
    r(0, 1) = k.X * k.Y * t - k.Z * s;
    r(0, 2) = k.X * k.Z * t + k.Y * s;
    r(1, 0) = k.Y * k.X * t + k.Z * s;
    r(1, 1) = c + k.Y * k.Y * t;
    r(1, 2) = k.Y * k.Z * t - k.X * s;
    r(2, 0) = k.Z * k.X * t - k.Y * s;
    r(2, 1) = k.Z * k.Y * t + k.X * s;
    r(2, 2) = c + k.Z * k.Z * t;
    return {r, {}};
}

/// The homogeneous 4×4 form, for graphics pipelines and other consumers that
/// expect one matrix.
template <Scalar T>
[[nodiscard]] constexpr BasicMat<T, 4> ToMat4(const BasicTransform<T>& t)
{
    BasicMat<T, 4> m = BasicMat<T, 4>::Identity();
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            m(i, j) = t.Linear(i, j);
        }
        m(i, 3) = t.Translation[i];
    }
    return m;
}

} // namespace Stereon

// std::format support: a matrix as nested rows, "[[a, b], [c, d]]"; a transform
// as "{Linear: [[...]], Translation: (...)}".

template <Stereon::Scalar T, std::size_t N>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicMat<T, N>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicMat<T, N>& m, FormatContext& ctx) const
    {
        auto out = ctx.out();
        *out++ = '[';
        for (std::size_t i = 0; i < N; ++i) {
            *out++ = '[';
            for (std::size_t j = 0; j < N; ++j) {
                if (j > 0) {
                    out = std::format_to(out, ", ");
                }
                out = std::format_to(out, "{}", m(i, j));
            }
            *out++ = ']';
            if (i + 1 < N) {
                out = std::format_to(out, ", ");
            }
        }
        *out++ = ']';
        return out;
    }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicTransform<T>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicTransform<T>& t, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "{{Linear: {}, Translation: {}}}", t.Linear,
                              t.Translation);
    }
};
