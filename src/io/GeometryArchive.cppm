// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

export module Stereon.Io:GeometryArchive;

import Stereon.Surface;

/// The geometry subset of the native .stn format (ADR-0012): curves, plane
/// curves and surfaces of the built-in types, with the model's numerical
/// policy (ADR-0028), saved and read back bit for bit.
///
/// Layout, every value little-endian:
///
/// - header, 32 bytes: the magic "STN\x1A"; the format version (u32); the
///   writing kernel's version (three u16) and two reserved zero bytes; the
///   length unit (u32, 1 = millimetre); the model resolution in that unit
///   (f64); four reserved zero bytes;
/// - chunks, each a tag of four ASCII letters, flags (u32, bit 0 =
///   required), the payload length (u64), the CRC-32 of the payload (u32)
///   and the payload.
///
/// The chunks, all required:
///
/// - "POLI", the numerical policy: angular resolution, tolerance ceiling
///   and frame allowance (f64 each);
/// - "CRVS", "PCRV" and "SRFS": curves in space, curves in the parameter
///   plane and surfaces. Each payload is a record count (u64) and the
///   records: the type's name (u16 length and its bytes, e.g. "Circle"), the
///   value's length (u64) and the value as Serialized writes it.
///
/// Format version 2 adds trimmed curves to "CRVS" and "PCRV", named
/// "TrimmedCurve" and "TrimmedCurve2d". The value of one is its basis as a
/// record of the same chunk (name, length and value, so a basis may itself
/// be trimmed, at most MaxTrimmedNesting deep), then the trim (Lower and
/// Upper f64, Reversed one byte) as Serialized writes it.
///
/// Format version 3 adds trimmed surfaces and sweeps to "SRFS". The value
/// of a "TrimmedSurface" is its basis as a record of the same chunk, then
/// the trim (LowerU, UpperU, LowerV, UpperV f64, Reversed one byte); the
/// value of an "ExtrusionSurface" is its profile as a record of the "CRVS"
/// kind, then the direction (three f64), and of a "RevolutionSurface" its
/// profile likewise, then the origin and the unit axis (three f64 each),
/// which the reader keeps as stored and refuses when not of unit length. A
/// basis or profile may itself be such a composite record, all nested at
/// most MaxTrimmedNesting deep together. Each earlier version is the next
/// without its new records, so a version-1 or version-2 file reads
/// unchanged.
///
/// A reader refuses a newer format version, a damaged chunk, an unknown
/// required chunk, an unknown type in a geometry chunk and any value its
/// type's checked factory rejects, each with a ParseError that says where;
/// it skips an unknown optional chunk with a warning. It never substitutes
/// its own defaults for the stored policy.

export namespace Stereon {

/// The length unit of a model. ADR-0001 fixes millimetres; the field exists
/// so that a file states its unit rather than implying it.
enum class LengthUnit : std::uint32_t {
    Millimetre = 1,
};

/// The numerical policy that gives a model's geometry its meaning
/// (ADR-0028): stored with it and read back as stored.
struct ModelPolicy {
    LengthUnit Unit = LengthUnit::Millimetre;
    /// Linear resolution, in the length unit.
    double Resolution = Context::DefaultResolution;
    /// Angular resolution, radians.
    double AngularResolution = Context::DefaultAngularResolution;
    /// The largest local tolerance an operation may grow, in the length unit.
    double ToleranceCeiling = Context::DefaultToleranceCeiling;
    /// How far the axes of an analytic frame may be from orthonormal when
    /// the geometry was accepted (FrameAllowance).
    double FrameAllowance = Stereon::FrameAllowance;

    bool operator==(const ModelPolicy&) const = default;
};

/// The policy of a Context, for saving geometry built under it.
[[nodiscard]] ModelPolicy PolicyOf(const Context& context);

/// Geometry and the policy it was built under.
struct GeometryArchive {
    ModelPolicy Policy;
    std::vector<Curve> Curves;
    std::vector<Curve2d> Curves2d;
    std::vector<Surface> Surfaces;
};

/// An archive read back, with what the reader skipped.
struct LoadedGeometry {
    GeometryArchive Archive;
    /// The format version the file was written in.
    std::uint32_t FormatVersion = 0;
    /// The version of the kernel that wrote it.
    Version WrittenBy{0, 0, 0};
    /// One line per optional chunk the reader did not know and skipped.
    std::vector<std::string> Warnings;
};

/// The format version this kernel writes and the newest it reads.
inline constexpr std::uint32_t GeometryFormatVersion = 3;

/// How deep the composite records (trimmed curves, trimmed surfaces and
/// sweeps) may be nested in a file: a trimmed curve whose basis is a
/// trimmed curve is depth 2, as is an extrusion of a trimmed curve. Deeper
/// nesting is refused on read, so that a hostile file cannot exhaust the
/// stack.
inline constexpr int MaxTrimmedNesting = 8;

/// The archive as bytes. NotSupported for geometry of a type the format has
/// no record for (the built-in types, trimmed curves and surfaces of them,
/// and extrusions and revolutions of any stored profile have one) and
/// without static reflection; InvalidInput for a policy value that is not
/// finite and positive.
[[nodiscard]] Result<std::vector<std::byte>> WriteGeometry(const GeometryArchive& archive);

/// The archive in bytes; ParseError for anything malformed (see above).
[[nodiscard]] Result<LoadedGeometry> ReadGeometry(std::span<const std::byte> bytes);

/// WriteGeometry to a file, replacing it; InvalidInput if it cannot be
/// written.
[[nodiscard]] Result<void> SaveGeometry(const std::filesystem::path& path,
                                        const GeometryArchive& archive);

/// ReadGeometry from a file; InvalidInput if it cannot be read.
[[nodiscard]] Result<LoadedGeometry> LoadGeometry(const std::filesystem::path& path);

/// A JSON view of an archive's bytes, for debugging and diffs (the
/// stn-dump view of ADR-0012): the header, the policy and every record with
/// its type and members. ParseError as ReadGeometry.
[[nodiscard]] Result<std::string> GeometryToJson(std::span<const std::byte> bytes);

} // namespace Stereon

export namespace Stereon::Detail {

/// The CRC-32 (IEEE 802.3, reflected, polynomial 0xEDB88320) of the bytes,
/// as each chunk stores for its payload.
[[nodiscard]] std::uint32_t Crc32(std::span<const std::byte> bytes) noexcept;

} // namespace Stereon::Detail
