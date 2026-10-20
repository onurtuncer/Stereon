// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

export module Stereon.Io:Debug;

import Stereon.Core;

/// The optional "DBUG" chunk of a dump (ADR-0035, section 1): who dumped and
/// where, and the debug primitives that are not kernel objects. Optional,
/// so a reader that does not know it skips it with a warning, and adding it
/// left the format version as it was.
///
/// Payload layout, every value little-endian: a layout version (u32, 1);
/// the name, the source file (each a u32 length and UTF-8 bytes), the line
/// (u32) and the time (i64, milliseconds since the Unix epoch); then five
/// sections, each a count (u64) and its items:
///
/// - points: x, y, z (f64 each), colour;
/// - polylines: a point count (u64), the points, colour;
/// - boxes: lower x, y, z, then upper x, y, z (f64 each), colour;
/// - triangle soups: a vertex count (u64), the vertices, an index count
///   (u64, a multiple of three), the indices (u32 each, below the vertex
///   count), colour;
/// - labels: the position, the text (u32 length and UTF-8 bytes), colour;
///
/// a colour being four bytes, red, green, blue and alpha, where alpha 0
/// means "the viewer's default". Coordinates are finite and a box's lower
/// corner is not above its upper one, both checked on write and on read.

export namespace Stereon {

/// A primitive's colour; alpha 0 (the default) leaves it to the viewer.
struct DebugColour {
    std::uint8_t R = 0;
    std::uint8_t G = 0;
    std::uint8_t B = 0;
    std::uint8_t A = 0;

    [[nodiscard]] constexpr bool IsSet() const noexcept { return A != 0; }
    constexpr bool operator==(const DebugColour&) const = default;
};

struct DebugPoint {
    Point3 Position;
    DebugColour Colour;

    bool operator==(const DebugPoint&) const = default;
};

struct DebugPolyline {
    std::vector<Point3> Points;
    DebugColour Colour;

    bool operator==(const DebugPolyline&) const = default;
};

struct DebugBox {
    Box3 Box;
    DebugColour Colour;

    bool operator==(const DebugBox&) const = default;
};

/// Triangles as a vertex list and three indices per triangle.
struct DebugTriangles {
    std::vector<Point3> Vertices;
    std::vector<std::uint32_t> Indices;
    DebugColour Colour;

    bool operator==(const DebugTriangles&) const = default;
};

struct DebugLabel {
    Point3 Position;
    std::string Text;
    DebugColour Colour;

    bool operator==(const DebugLabel&) const = default;
};

/// The content of a "DBUG" chunk.
struct DebugChunk {
    /// The dump's name, as its file is named.
    std::string Name;
    /// Where the dump was made; empty and 0 when not known.
    std::string File;
    std::uint32_t Line = 0;
    /// When, in milliseconds since the Unix epoch.
    std::int64_t Time = 0;

    std::vector<DebugPoint> Points;
    std::vector<DebugPolyline> Polylines;
    std::vector<DebugBox> Boxes;
    std::vector<DebugTriangles> Triangles;
    std::vector<DebugLabel> Labels;

    bool operator==(const DebugChunk&) const = default;
};

} // namespace Stereon

export namespace Stereon::Detail {

/// The chunk's payload; InvalidInput for a coordinate that is not finite,
/// a box with its corners swapped, or a triangle index out of range.
[[nodiscard]] Result<std::vector<std::byte>> DebugPayload(const DebugChunk& chunk);

/// The chunk from its payload; ParseError for anything malformed.
[[nodiscard]] Result<DebugChunk> ReadDebugPayload(std::span<const std::byte> payload);

/// The chunk as JSON, for the stn-dump view.
[[nodiscard]] std::string DebugJson(const DebugChunk& chunk);

} // namespace Stereon::Detail
