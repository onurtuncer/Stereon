// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <format>

export module Stereon.Core:Box;

import :Interval;
import :Scalar;
import :Vector;

/// Axis-aligned bounding boxes as products of closed intervals, one per axis.
/// Containment, overlap, union and intersection are the per-axis Interval
/// operations, so a box is empty as soon as one axis is, and Box3{} is the
/// empty box, the identity of Hull. Like Interval, boxes are double only.

namespace Stereon::Detail {

template <typename B, typename P>
[[nodiscard]] constexpr B BoxFromCorners(const P& a, const P& b)
{
    B box;
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        box[i] = Interval{std::min(a[i], b[i]), std::max(a[i], b[i])};
    }
    return box;
}

} // namespace Stereon::Detail

export namespace Stereon {

struct Box2 {
    using Point = Point2;
    using Vec = Vec2;
    static constexpr std::size_t Dimension = 2;

    Interval X;
    Interval Y;

    [[nodiscard]] static constexpr Box2 Empty() { return {}; }
    [[nodiscard]] static constexpr Box2 Everything()
    {
        return {Interval::Everything(), Interval::Everything()};
    }
    /// The box spanned by two opposite corners, in either order.
    [[nodiscard]] static constexpr Box2 FromCorners(const Point2& a, const Point2& b)
    {
        return Detail::BoxFromCorners<Box2>(a, b);
    }

    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return X.IsEmpty() || Y.IsEmpty(); }

    constexpr bool operator==(const Box2&) const = default;

    constexpr Interval& operator[](const std::size_t axis) STEREON_HOT_PRE(axis < 2)
    {
        return axis == 0 ? X : Y;
    }
    constexpr const Interval& operator[](const std::size_t axis) const STEREON_HOT_PRE(axis < 2)
    {
        return axis == 0 ? X : Y;
    }
};

struct Box3 {
    using Point = Point3;
    using Vec = Vec3;
    static constexpr std::size_t Dimension = 3;

    Interval X;
    Interval Y;
    Interval Z;

    [[nodiscard]] static constexpr Box3 Empty() { return {}; }
    [[nodiscard]] static constexpr Box3 Everything()
    {
        return {Interval::Everything(), Interval::Everything(), Interval::Everything()};
    }
    [[nodiscard]] static constexpr Box3 FromCorners(const Point3& a, const Point3& b)
    {
        return Detail::BoxFromCorners<Box3>(a, b);
    }

    [[nodiscard]] constexpr bool IsEmpty() const noexcept
    {
        return X.IsEmpty() || Y.IsEmpty() || Z.IsEmpty();
    }

    constexpr bool operator==(const Box3&) const = default;

    constexpr Interval& operator[](const std::size_t axis) STEREON_HOT_PRE(axis < 3)
    {
        return axis == 0 ? X : (axis == 1 ? Y : Z);
    }
    constexpr const Interval& operator[](const std::size_t axis) const STEREON_HOT_PRE(axis < 3)
    {
        return axis == 0 ? X : (axis == 1 ? Y : Z);
    }
};

template <typename B>
concept AnyBox = std::same_as<B, Box2> || std::same_as<B, Box3>;

template <AnyBox B>
[[nodiscard]] constexpr bool Contains(const B& box, const typename B::Point& p) noexcept
{
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        if (!Contains(box[i], p[i])) {
            return false;
        }
    }
    return true;
}

/// Every box contains the empty box.
template <AnyBox B>
[[nodiscard]] constexpr bool Contains(const B& outer, const B& inner) noexcept
{
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        if (!Contains(outer[i], inner[i])) {
            return false;
        }
    }
    return true;
}

/// Boxes that share only a face, edge or corner intersect.
template <AnyBox B>
[[nodiscard]] constexpr bool Intersects(const B& a, const B& b) noexcept
{
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        if (!Intersects(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

template <AnyBox B>
[[nodiscard]] constexpr B Hull(const B& box, const typename B::Point& p)
{
    B r = box;
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        r[i] = Hull(r[i], p[i]);
    }
    return r;
}

template <AnyBox B>
[[nodiscard]] constexpr B Hull(const B& a, const B& b)
{
    B r;
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        r[i] = Hull(a[i], b[i]);
    }
    return r;
}

/// The empty box, B{}, when the boxes do not intersect.
template <AnyBox B>
[[nodiscard]] constexpr B Intersection(const B& a, const B& b)
{
    B r;
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        r[i] = Intersection(a[i], b[i]);
        if (r[i].IsEmpty()) {
            return {};
        }
    }
    return r;
}

template <AnyBox B>
[[nodiscard]] constexpr typename B::Point MinCorner(const B& box) STEREON_HOT_PRE(!box.IsEmpty())
{
    typename B::Point p{};
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        p[i] = box[i].GetLower();
    }
    return p;
}

template <AnyBox B>
[[nodiscard]] constexpr typename B::Point MaxCorner(const B& box) STEREON_HOT_PRE(!box.IsEmpty())
{
    typename B::Point p{};
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        p[i] = box[i].GetUpper();
    }
    return p;
}

template <AnyBox B>
[[nodiscard]] constexpr typename B::Point Center(const B& box) STEREON_HOT_PRE(!box.IsEmpty())
{
    typename B::Point p{};
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        p[i] = Midpoint(box[i]);
    }
    return p;
}

/// The vector from the minimum to the maximum corner.
template <AnyBox B>
[[nodiscard]] constexpr typename B::Vec Diagonal(const B& box) STEREON_HOT_PRE(!box.IsEmpty())
{
    typename B::Vec v{};
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        v[i] = Width(box[i]);
    }
    return v;
}

/// Every face moved outward by margin; see Inflated(Interval, double).
template <AnyBox B>
[[nodiscard]] constexpr B Inflated(const B& box, const double margin)
{
    B r;
    for (std::size_t i = 0; i < B::Dimension; ++i) {
        r[i] = Inflated(box[i], margin);
        if (r[i].IsEmpty()) {
            return {};
        }
    }
    return r;
}

} // namespace Stereon

// std::format support: "[min corner, max corner]", or "[]" when empty.

namespace Stereon::Detail {

template <AnyBox B>
struct BoxFormatter {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return EmptySpecParser::Parse(ctx);
    }
    template <typename FormatContext>
    auto format(const B& box, FormatContext& ctx) const
    {
        if (box.IsEmpty()) {
            return std::format_to(ctx.out(), "[]");
        }
        return std::format_to(ctx.out(), "[{}, {}]", MinCorner(box), MaxCorner(box));
    }
};

} // namespace Stereon::Detail

template <>
struct std::formatter<Stereon::Box2> : Stereon::Detail::BoxFormatter<Stereon::Box2> {};

template <>
struct std::formatter<Stereon::Box3> : Stereon::Detail::BoxFormatter<Stereon::Box3> {};
