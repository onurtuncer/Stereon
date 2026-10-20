// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/Conversion.cpp.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <format>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Io;

import :Debug;

namespace Stereon::Detail {

namespace {

constexpr std::uint32_t DebugLayoutVersion = 1;

/// The smallest encoding of one item of each section, to refuse a count no
/// payload of this length could hold before reserving for it.
constexpr std::size_t PointBytes = 3 * 8;
constexpr std::size_t ColourBytes = 4;

// ---------------------------------------------------------------------------
// Writing

class Writer {
public:
    void U32(const std::uint32_t value) { Put(value, 4); }
    void U64(const std::uint64_t value) { Put(value, 8); }
    void I64(const std::int64_t value) { Put(static_cast<std::uint64_t>(value), 8); }
    void F64(const double value) { Put(std::bit_cast<std::uint64_t>(value), 8); }
    void Text(const std::string_view text)
    {
        U32(static_cast<std::uint32_t>(text.size()));
        for (const char c : text) {
            m_Bytes.push_back(static_cast<std::byte>(c));
        }
    }
    void Point(const Point3& p)
    {
        F64(p.X);
        F64(p.Y);
        F64(p.Z);
    }
    void Colour(const DebugColour& c)
    {
        for (const std::uint8_t b : {c.R, c.G, c.B, c.A}) {
            m_Bytes.push_back(static_cast<std::byte>(b));
        }
    }
    [[nodiscard]] std::vector<std::byte> Take() { return std::move(m_Bytes); }

private:
    void Put(const std::uint64_t value, const std::size_t count)
    {
        for (std::size_t i = 0; i < count; ++i) {
            m_Bytes.push_back(static_cast<std::byte>((value >> (8 * i)) & 0xFFU));
        }
    }

    std::vector<std::byte> m_Bytes;
};

[[nodiscard]] bool Finite(const Point3& p)
{
    return std::isfinite(p.X) && std::isfinite(p.Y) && std::isfinite(p.Z);
}

[[nodiscard]] bool AllFinite(const std::vector<Point3>& points)
{
    for (const Point3& p : points) {
        if (!Finite(p)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidBox(const Box3& box)
{
    for (const Interval& i : {box.X, box.Y, box.Z}) {
        if (!std::isfinite(i.GetLower()) || !std::isfinite(i.GetUpper()) ||
            !(i.GetLower() <= i.GetUpper())) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidTriangles(const DebugTriangles& t)
{
    if (!AllFinite(t.Vertices) || t.Indices.size() % 3 != 0) {
        return false;
    }
    for (const std::uint32_t index : t.Indices) {
        if (index >= t.Vertices.size()) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool ValidPrimitives(const DebugChunk& chunk)
{
    for (const DebugPoint& p : chunk.Points) {
        if (!Finite(p.Position)) {
            return false;
        }
    }
    for (const DebugPolyline& p : chunk.Polylines) {
        if (!AllFinite(p.Points)) {
            return false;
        }
    }
    for (const DebugBox& b : chunk.Boxes) {
        if (!ValidBox(b.Box)) {
            return false;
        }
    }
    for (const DebugTriangles& t : chunk.Triangles) {
        if (!ValidTriangles(t)) {
            return false;
        }
    }
    for (const DebugLabel& l : chunk.Labels) {
        if (!Finite(l.Position)) {
            return false;
        }
    }
    return true;
}

void WritePrimitives(Writer& w, const DebugChunk& chunk)
{
    w.U64(chunk.Points.size());
    for (const DebugPoint& p : chunk.Points) {
        w.Point(p.Position);
        w.Colour(p.Colour);
    }
    w.U64(chunk.Polylines.size());
    for (const DebugPolyline& p : chunk.Polylines) {
        w.U64(p.Points.size());
        for (const Point3& q : p.Points) {
            w.Point(q);
        }
        w.Colour(p.Colour);
    }
    w.U64(chunk.Boxes.size());
    for (const DebugBox& b : chunk.Boxes) {
        w.Point(Point3{b.Box.X.GetLower(), b.Box.Y.GetLower(), b.Box.Z.GetLower()});
        w.Point(Point3{b.Box.X.GetUpper(), b.Box.Y.GetUpper(), b.Box.Z.GetUpper()});
        w.Colour(b.Colour);
    }
    w.U64(chunk.Triangles.size());
    for (const DebugTriangles& t : chunk.Triangles) {
        w.U64(t.Vertices.size());
        for (const Point3& q : t.Vertices) {
            w.Point(q);
        }
        w.U64(t.Indices.size());
        for (const std::uint32_t index : t.Indices) {
            w.U32(index);
        }
        w.Colour(t.Colour);
    }
    w.U64(chunk.Labels.size());
    for (const DebugLabel& l : chunk.Labels) {
        w.Point(l.Position);
        w.Text(l.Text);
        w.Colour(l.Colour);
    }
}

// ---------------------------------------------------------------------------
// Reading

[[nodiscard]] std::unexpected<KernelError> Malformed(const std::string_view what,
                                                     const std::size_t offset)
{
    return Fail(ErrorCode::ParseError,
                std::format("ReadGeometry: the \"DBUG\" chunk is malformed at payload offset "
                            "{}: {}",
                            offset, what));
}

class Reader {
public:
    explicit Reader(const std::span<const std::byte> bytes)
        : m_Bytes(bytes)
    {}

    [[nodiscard]] bool AtEnd() const { return m_Offset == m_Bytes.size(); }
    [[nodiscard]] std::size_t Offset() const { return m_Offset; }

    [[nodiscard]] Result<std::uint64_t> Unsigned(const std::size_t count)
    {
        if (count > Left()) {
            return Malformed("truncated", m_Offset);
        }
        std::uint64_t value = 0;
        for (std::size_t i = 0; i < count; ++i) {
            value |= std::uint64_t{std::to_integer<std::uint8_t>(m_Bytes[m_Offset + i])} << (8 * i);
        }
        m_Offset += count;
        return value;
    }
    [[nodiscard]] Result<double> F64()
    {
        STEREON_TRY(const std::uint64_t bits, Unsigned(8));
        return std::bit_cast<double>(bits);
    }
    /// A count of items at least `itemBytes` long each; ParseError for more
    /// than what is left could hold.
    [[nodiscard]] Result<std::size_t> Count(const std::size_t itemBytes)
    {
        STEREON_TRY(const std::uint64_t count, Unsigned(8));
        if (count > Left() / itemBytes) {
            return Malformed("a count larger than the rest of the chunk can hold", m_Offset);
        }
        return static_cast<std::size_t>(count);
    }
    [[nodiscard]] Result<std::string> Text()
    {
        STEREON_TRY(const std::uint64_t length, Unsigned(4));
        if (length > Left()) {
            return Malformed("truncated text", m_Offset);
        }
        std::string text(static_cast<std::size_t>(length), '\0');
        for (char& c : text) {
            c = static_cast<char>(std::to_integer<unsigned char>(m_Bytes[m_Offset++]));
        }
        return text;
    }
    [[nodiscard]] Result<Point3> Point()
    {
        Point3 p;
        STEREON_TRY(p.X, F64());
        STEREON_TRY(p.Y, F64());
        STEREON_TRY(p.Z, F64());
        if (!Finite(p)) {
            return Malformed("a coordinate that is not finite", m_Offset);
        }
        return p;
    }
    [[nodiscard]] Result<DebugColour> Colour()
    {
        STEREON_TRY(const std::uint64_t rgba, Unsigned(4));
        DebugColour c;
        c.R = static_cast<std::uint8_t>(rgba & 0xFFU);
        c.G = static_cast<std::uint8_t>((rgba >> 8) & 0xFFU);
        c.B = static_cast<std::uint8_t>((rgba >> 16) & 0xFFU);
        c.A = static_cast<std::uint8_t>((rgba >> 24) & 0xFFU);
        return c;
    }
    [[nodiscard]] Result<std::vector<Point3>> Points()
    {
        STEREON_TRY(const std::size_t count, Count(PointBytes));
        std::vector<Point3> points;
        points.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            STEREON_TRY(const Point3 p, Point());
            points.push_back(p);
        }
        return points;
    }

private:
    [[nodiscard]] std::size_t Left() const { return m_Bytes.size() - m_Offset; }

    std::span<const std::byte> m_Bytes;
    std::size_t m_Offset = 0;
};

[[nodiscard]] Result<void> ReadPoints(Reader& r, DebugChunk& chunk)
{
    STEREON_TRY(const std::size_t count, r.Count(PointBytes + ColourBytes));
    for (std::size_t i = 0; i < count; ++i) {
        DebugPoint p;
        STEREON_TRY(p.Position, r.Point());
        STEREON_TRY(p.Colour, r.Colour());
        chunk.Points.push_back(p);
    }
    return {};
}

[[nodiscard]] Result<void> ReadPolylines(Reader& r, DebugChunk& chunk)
{
    STEREON_TRY(const std::size_t count, r.Count(8 + ColourBytes));
    for (std::size_t i = 0; i < count; ++i) {
        DebugPolyline p;
        STEREON_TRY(p.Points, r.Points());
        STEREON_TRY(p.Colour, r.Colour());
        chunk.Polylines.push_back(std::move(p));
    }
    return {};
}

[[nodiscard]] Result<void> ReadBoxes(Reader& r, DebugChunk& chunk)
{
    STEREON_TRY(const std::size_t count, r.Count(2 * PointBytes + ColourBytes));
    for (std::size_t i = 0; i < count; ++i) {
        STEREON_TRY(const Point3 lower, r.Point());
        STEREON_TRY(const Point3 upper, r.Point());
        if (!(lower.X <= upper.X && lower.Y <= upper.Y && lower.Z <= upper.Z)) {
            return Malformed("a box whose lower corner is above its upper one", r.Offset());
        }
        DebugBox b;
        b.Box = Box3{Interval{lower.X, upper.X}, Interval{lower.Y, upper.Y},
                     Interval{lower.Z, upper.Z}};
        STEREON_TRY(b.Colour, r.Colour());
        chunk.Boxes.push_back(b);
    }
    return {};
}

[[nodiscard]] Result<void> ReadTriangles(Reader& r, DebugChunk& chunk)
{
    STEREON_TRY(const std::size_t count, r.Count(8 + 8 + ColourBytes));
    for (std::size_t i = 0; i < count; ++i) {
        DebugTriangles t;
        STEREON_TRY(t.Vertices, r.Points());
        STEREON_TRY(const std::size_t indices, r.Count(4));
        t.Indices.reserve(indices);
        for (std::size_t k = 0; k < indices; ++k) {
            STEREON_TRY(const std::uint64_t index, r.Unsigned(4));
            t.Indices.push_back(static_cast<std::uint32_t>(index));
        }
        if (!ValidTriangles(t)) {
            return Malformed("a triangle index out of range, or an index count that is not a "
                             "multiple of three",
                             r.Offset());
        }
        STEREON_TRY(t.Colour, r.Colour());
        chunk.Triangles.push_back(std::move(t));
    }
    return {};
}

[[nodiscard]] Result<void> ReadLabels(Reader& r, DebugChunk& chunk)
{
    STEREON_TRY(const std::size_t count, r.Count(PointBytes + 4 + ColourBytes));
    for (std::size_t i = 0; i < count; ++i) {
        DebugLabel l;
        STEREON_TRY(l.Position, r.Point());
        STEREON_TRY(l.Text, r.Text());
        STEREON_TRY(l.Colour, r.Colour());
        chunk.Labels.push_back(std::move(l));
    }
    return {};
}

// ---------------------------------------------------------------------------
// JSON

[[nodiscard]] std::string Quoted(const std::string_view text)
{
    std::string out = "\"";
    for (const char c : text) {
        const auto u = static_cast<unsigned char>(c);
        if (c == '"' || c == '\\') {
            out += '\\';
            out += c;
        } else if (u < 0x20) {
            out += std::format("\\u{:04x}", static_cast<unsigned>(u));
        } else {
            out += c;
        }
    }
    return out + "\"";
}

[[nodiscard]] std::string PointJson(const Point3& p)
{
    return std::format("[{},{},{}]", p.X, p.Y, p.Z);
}

[[nodiscard]] std::string PointsJson(const std::vector<Point3>& points)
{
    std::string out = "[";
    for (std::size_t i = 0; i < points.size(); ++i) {
        out += (i > 0 ? "," : "") + PointJson(points[i]);
    }
    return out + "]";
}

[[nodiscard]] std::string ColourJson(const DebugColour& c)
{
    if (!c.IsSet()) {
        return "null";
    }
    return std::format("[{},{},{},{}]", static_cast<unsigned>(c.R), static_cast<unsigned>(c.G),
                       static_cast<unsigned>(c.B), static_cast<unsigned>(c.A));
}

template <typename T, typename F>
[[nodiscard]] std::string ArrayJson(const std::vector<T>& items, F item)
{
    std::string out = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        out += (i > 0 ? "," : "") + item(items[i]);
    }
    return out + "]";
}

[[nodiscard]] std::string IndicesJson(const std::vector<std::uint32_t>& indices)
{
    std::string out = "[";
    for (std::size_t i = 0; i < indices.size(); ++i) {
        out += (i > 0 ? "," : "") + std::to_string(indices[i]);
    }
    return out + "]";
}

} // namespace

Result<std::vector<std::byte>> DebugPayload(const DebugChunk& chunk)
{
    if (!ValidPrimitives(chunk)) {
        return Fail(ErrorCode::InvalidInput,
                    "WriteGeometry: a debug primitive has a coordinate that is not finite, a box "
                    "with its corners swapped or a triangle index out of range");
    }
    Writer w;
    w.U32(DebugLayoutVersion);
    w.Text(chunk.Name);
    w.Text(chunk.File);
    w.U32(chunk.Line);
    w.I64(chunk.Time);
    WritePrimitives(w, chunk);
    return w.Take();
}

Result<DebugChunk> ReadDebugPayload(const std::span<const std::byte> payload)
{
    Reader r(payload);
    STEREON_TRY(const std::uint64_t version, r.Unsigned(4));
    if (version != DebugLayoutVersion) {
        return Malformed(
            std::format("layout version {}, this reader knows {}", version, DebugLayoutVersion), 0);
    }
    DebugChunk chunk;
    STEREON_TRY(chunk.Name, r.Text());
    STEREON_TRY(chunk.File, r.Text());
    STEREON_TRY(const std::uint64_t line, r.Unsigned(4));
    chunk.Line = static_cast<std::uint32_t>(line);
    STEREON_TRY(const std::uint64_t time, r.Unsigned(8));
    chunk.Time = static_cast<std::int64_t>(time);
    STEREON_TRY_VOID(ReadPoints(r, chunk));
    STEREON_TRY_VOID(ReadPolylines(r, chunk));
    STEREON_TRY_VOID(ReadBoxes(r, chunk));
    STEREON_TRY_VOID(ReadTriangles(r, chunk));
    STEREON_TRY_VOID(ReadLabels(r, chunk));
    if (!r.AtEnd()) {
        return Malformed("bytes after the last section", r.Offset());
    }
    return chunk;
}

std::string DebugJson(const DebugChunk& chunk)
{
    const auto point = [](const DebugPoint& p) {
        return std::format("{{\"at\":{},\"colour\":{}}}", PointJson(p.Position),
                           ColourJson(p.Colour));
    };
    const auto polyline = [](const DebugPolyline& p) {
        return std::format("{{\"points\":{},\"colour\":{}}}", PointsJson(p.Points),
                           ColourJson(p.Colour));
    };
    const auto box = [](const DebugBox& b) {
        const Point3 lower{b.Box.X.GetLower(), b.Box.Y.GetLower(), b.Box.Z.GetLower()};
        const Point3 upper{b.Box.X.GetUpper(), b.Box.Y.GetUpper(), b.Box.Z.GetUpper()};
        return std::format("{{\"lower\":{},\"upper\":{},\"colour\":{}}}", PointJson(lower),
                           PointJson(upper), ColourJson(b.Colour));
    };
    const auto triangles = [](const DebugTriangles& t) {
        return std::format("{{\"vertices\":{},\"indices\":{},\"colour\":{}}}",
                           PointsJson(t.Vertices), IndicesJson(t.Indices), ColourJson(t.Colour));
    };
    const auto label = [](const DebugLabel& l) {
        return std::format("{{\"at\":{},\"text\":{},\"colour\":{}}}", PointJson(l.Position),
                           Quoted(l.Text), ColourJson(l.Colour));
    };
    return std::format("{{\"name\":{},\"file\":{},\"line\":{},\"time\":{},\"points\":{},"
                       "\"polylines\":{},\"boxes\":{},\"triangles\":{},\"labels\":{}}}",
                       Quoted(chunk.Name), Quoted(chunk.File), chunk.Line, chunk.Time,
                       ArrayJson(chunk.Points, point), ArrayJson(chunk.Polylines, polyline),
                       ArrayJson(chunk.Boxes, box), ArrayJson(chunk.Triangles, triangles),
                       ArrayJson(chunk.Labels, label));
}

} // namespace Stereon::Detail
