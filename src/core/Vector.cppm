// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cmath>
#include <cstddef>
#include <format>

export module Stereon.Core:Vector;

import :Scalar;

/// Fixed-size vectors and points. Vectors are displacements, points are
/// positions: Point - Point is a Vec, Point + Vec is a Point, and Point + Point
/// does not compile. The Basic* templates are scalar-generic (ADR-0021); Vec2,
/// Vec3, Point2 and Point3 are their double instances and what most code uses.
///
/// Lengths are in millimetres (ADR-0001). Comparisons are exact; tolerance
/// decisions belong to the code that owns the tolerance.
export namespace Stereon {

template <Scalar T>
struct BasicVec2 {
    T X{};
    T Y{};

    // Member, not hidden friend: a defaulted friend operator== exported from a
    // module crashes GCC 16.2 (docs/toolchain/gcc16-feature-probe.md).
    constexpr bool operator==(const BasicVec2&) const = default;

    constexpr T& operator[](const std::size_t i) STEREON_HOT_PRE(i < 2) { return i == 0 ? X : Y; }
    constexpr const T& operator[](const std::size_t i) const STEREON_HOT_PRE(i < 2)
    {
        return i == 0 ? X : Y;
    }

    constexpr BasicVec2& operator+=(const BasicVec2& o)
    {
        X = X + o.X;
        Y = Y + o.Y;
        return *this;
    }
    constexpr BasicVec2& operator-=(const BasicVec2& o)
    {
        X = X - o.X;
        Y = Y - o.Y;
        return *this;
    }
    constexpr BasicVec2& operator*=(const T& s)
    {
        X = X * s;
        Y = Y * s;
        return *this;
    }
    constexpr BasicVec2& operator/=(const T& s)
    {
        X = X / s;
        Y = Y / s;
        return *this;
    }

    friend constexpr BasicVec2 operator+(BasicVec2 a, const BasicVec2& b) { return a += b; }
    friend constexpr BasicVec2 operator-(BasicVec2 a, const BasicVec2& b) { return a -= b; }
    friend constexpr BasicVec2 operator-(const BasicVec2& a) { return {-a.X, -a.Y}; }
    friend constexpr BasicVec2 operator*(BasicVec2 a, const T& s) { return a *= s; }
    friend constexpr BasicVec2 operator*(const T& s, BasicVec2 a) { return a *= s; }
    friend constexpr BasicVec2 operator/(BasicVec2 a, const T& s) { return a /= s; }
};

template <Scalar T>
struct BasicVec3 {
    T X{};
    T Y{};
    T Z{};

    constexpr bool operator==(const BasicVec3&) const = default;

    constexpr T& operator[](const std::size_t i) STEREON_HOT_PRE(i < 3)
    {
        return i == 0 ? X : (i == 1 ? Y : Z);
    }
    constexpr const T& operator[](const std::size_t i) const STEREON_HOT_PRE(i < 3)
    {
        return i == 0 ? X : (i == 1 ? Y : Z);
    }

    constexpr BasicVec3& operator+=(const BasicVec3& o)
    {
        X = X + o.X;
        Y = Y + o.Y;
        Z = Z + o.Z;
        return *this;
    }
    constexpr BasicVec3& operator-=(const BasicVec3& o)
    {
        X = X - o.X;
        Y = Y - o.Y;
        Z = Z - o.Z;
        return *this;
    }
    constexpr BasicVec3& operator*=(const T& s)
    {
        X = X * s;
        Y = Y * s;
        Z = Z * s;
        return *this;
    }
    constexpr BasicVec3& operator/=(const T& s)
    {
        X = X / s;
        Y = Y / s;
        Z = Z / s;
        return *this;
    }

    friend constexpr BasicVec3 operator+(BasicVec3 a, const BasicVec3& b) { return a += b; }
    friend constexpr BasicVec3 operator-(BasicVec3 a, const BasicVec3& b) { return a -= b; }
    friend constexpr BasicVec3 operator-(const BasicVec3& a) { return {-a.X, -a.Y, -a.Z}; }
    friend constexpr BasicVec3 operator*(BasicVec3 a, const T& s) { return a *= s; }
    friend constexpr BasicVec3 operator*(const T& s, BasicVec3 a) { return a *= s; }
    friend constexpr BasicVec3 operator/(BasicVec3 a, const T& s) { return a /= s; }
};

template <Scalar T>
struct BasicPoint2 {
    T X{};
    T Y{};

    constexpr bool operator==(const BasicPoint2&) const = default;

    constexpr T& operator[](const std::size_t i) STEREON_HOT_PRE(i < 2) { return i == 0 ? X : Y; }
    constexpr const T& operator[](const std::size_t i) const STEREON_HOT_PRE(i < 2)
    {
        return i == 0 ? X : Y;
    }

    constexpr BasicPoint2& operator+=(const BasicVec2<T>& v)
    {
        X = X + v.X;
        Y = Y + v.Y;
        return *this;
    }
    constexpr BasicPoint2& operator-=(const BasicVec2<T>& v)
    {
        X = X - v.X;
        Y = Y - v.Y;
        return *this;
    }

    friend constexpr BasicPoint2 operator+(BasicPoint2 p, const BasicVec2<T>& v) { return p += v; }
    friend constexpr BasicPoint2 operator+(const BasicVec2<T>& v, BasicPoint2 p) { return p += v; }
    friend constexpr BasicPoint2 operator-(BasicPoint2 p, const BasicVec2<T>& v) { return p -= v; }
    friend constexpr BasicVec2<T> operator-(const BasicPoint2& p, const BasicPoint2& q)
    {
        return {p.X - q.X, p.Y - q.Y};
    }
};

template <Scalar T>
struct BasicPoint3 {
    T X{};
    T Y{};
    T Z{};

    constexpr bool operator==(const BasicPoint3&) const = default;

    constexpr T& operator[](const std::size_t i) STEREON_HOT_PRE(i < 3)
    {
        return i == 0 ? X : (i == 1 ? Y : Z);
    }
    constexpr const T& operator[](const std::size_t i) const STEREON_HOT_PRE(i < 3)
    {
        return i == 0 ? X : (i == 1 ? Y : Z);
    }

    constexpr BasicPoint3& operator+=(const BasicVec3<T>& v)
    {
        X = X + v.X;
        Y = Y + v.Y;
        Z = Z + v.Z;
        return *this;
    }
    constexpr BasicPoint3& operator-=(const BasicVec3<T>& v)
    {
        X = X - v.X;
        Y = Y - v.Y;
        Z = Z - v.Z;
        return *this;
    }

    friend constexpr BasicPoint3 operator+(BasicPoint3 p, const BasicVec3<T>& v) { return p += v; }
    friend constexpr BasicPoint3 operator+(const BasicVec3<T>& v, BasicPoint3 p) { return p += v; }
    friend constexpr BasicPoint3 operator-(BasicPoint3 p, const BasicVec3<T>& v) { return p -= v; }
    friend constexpr BasicVec3<T> operator-(const BasicPoint3& p, const BasicPoint3& q)
    {
        return {p.X - q.X, p.Y - q.Y, p.Z - q.Z};
    }
};

using Vec2 = BasicVec2<double>;
using Vec3 = BasicVec3<double>;
using Point2 = BasicPoint2<double>;
using Point3 = BasicPoint3<double>;

// Vector algebra.

template <Scalar T>
[[nodiscard]] constexpr T Dot(const BasicVec2<T>& a, const BasicVec2<T>& b)
{
    return a.X * b.X + a.Y * b.Y;
}

template <Scalar T>
[[nodiscard]] constexpr T Dot(const BasicVec3<T>& a, const BasicVec3<T>& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

/// The 2D cross product a.X * b.Y - a.Y * b.X: positive when b lies
/// counter-clockwise from a.
template <Scalar T>
[[nodiscard]] constexpr T Cross(const BasicVec2<T>& a, const BasicVec2<T>& b)
{
    return a.X * b.Y - a.Y * b.X;
}

/// Right-handed: Cross(X, Y) == Z.
template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Cross(const BasicVec3<T>& a, const BasicVec3<T>& b)
{
    return {a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X};
}

/// a rotated a quarter turn counter-clockwise.
template <Scalar T>
[[nodiscard]] constexpr BasicVec2<T> Perpendicular(const BasicVec2<T>& a)
{
    return {-a.Y, a.X};
}

template <Scalar T>
[[nodiscard]] constexpr T LengthSquared(const BasicVec2<T>& a)
{
    return Dot(a, a);
}

template <Scalar T>
[[nodiscard]] constexpr T LengthSquared(const BasicVec3<T>& a)
{
    return Dot(a, a);
}

template <Scalar T>
[[nodiscard]] constexpr T Length(const BasicVec2<T>& a)
{
    using std::sqrt;
    return sqrt(LengthSquared(a));
}

template <Scalar T>
[[nodiscard]] constexpr T Length(const BasicVec3<T>& a)
{
    using std::sqrt;
    return sqrt(LengthSquared(a));
}

/// The zero vector has no direction; passing it is a bug in the caller, so it
/// is a contract, not a KernelError (ADR-0006).
template <Scalar T>
[[nodiscard]] constexpr BasicVec2<T> Normalized(const BasicVec2<T>& a)
    STEREON_HOT_PRE(T{} < LengthSquared(a))
{
    return a / Length(a);
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Normalized(const BasicVec3<T>& a)
    STEREON_HOT_PRE(T{} < LengthSquared(a))
{
    return a / Length(a);
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec2<T> Lerp(const BasicVec2<T>& a, const BasicVec2<T>& b, const T& t)
{
    return a + (b - a) * t;
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> Lerp(const BasicVec3<T>& a, const BasicVec3<T>& b, const T& t)
{
    return a + (b - a) * t;
}

// Points.

/// The position vector of p, from the origin. Explicit because a point and a
/// displacement are different things; the conversion is never implicit.
template <Scalar T>
[[nodiscard]] constexpr BasicVec2<T> ToVec(const BasicPoint2<T>& p)
{
    return {p.X, p.Y};
}

template <Scalar T>
[[nodiscard]] constexpr BasicVec3<T> ToVec(const BasicPoint3<T>& p)
{
    return {p.X, p.Y, p.Z};
}

/// The point at the tip of v when it starts at the origin.
template <Scalar T>
[[nodiscard]] constexpr BasicPoint2<T> ToPoint(const BasicVec2<T>& v)
{
    return {v.X, v.Y};
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> ToPoint(const BasicVec3<T>& v)
{
    return {v.X, v.Y, v.Z};
}

template <Scalar T>
[[nodiscard]] constexpr T DistanceSquared(const BasicPoint2<T>& p, const BasicPoint2<T>& q)
{
    return LengthSquared(p - q);
}

template <Scalar T>
[[nodiscard]] constexpr T DistanceSquared(const BasicPoint3<T>& p, const BasicPoint3<T>& q)
{
    return LengthSquared(p - q);
}

template <Scalar T>
[[nodiscard]] constexpr T Distance(const BasicPoint2<T>& p, const BasicPoint2<T>& q)
{
    return Length(p - q);
}

template <Scalar T>
[[nodiscard]] constexpr T Distance(const BasicPoint3<T>& p, const BasicPoint3<T>& q)
{
    return Length(p - q);
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint2<T> Lerp(const BasicPoint2<T>& p, const BasicPoint2<T>& q,
                                            const T& t)
{
    return p + (q - p) * t;
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Lerp(const BasicPoint3<T>& p, const BasicPoint3<T>& q,
                                            const T& t)
{
    return p + (q - p) * t;
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint2<T> Midpoint(const BasicPoint2<T>& p, const BasicPoint2<T>& q)
{
    return Lerp(p, q, T{0.5});
}

template <Scalar T>
[[nodiscard]] constexpr BasicPoint3<T> Midpoint(const BasicPoint3<T>& p, const BasicPoint3<T>& q)
{
    return Lerp(p, q, T{0.5});
}

} // namespace Stereon

// std::format support: "(x, y)" and "(x, y, z)", with no format specification
// of their own. Only for scalars that are themselves formattable.

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicVec2<T>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicVec2<T>& v, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "({}, {})", v.X, v.Y);
    }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicVec3<T>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicVec3<T>& v, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "({}, {}, {})", v.X, v.Y, v.Z);
    }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicPoint2<T>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicPoint2<T>& p, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "({}, {})", p.X, p.Y);
    }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct std::formatter<Stereon::BasicPoint3<T>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return Stereon::Detail::EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const Stereon::BasicPoint3<T>& p, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "({}, {}, {})", p.X, p.Y, p.Z);
    }
};
