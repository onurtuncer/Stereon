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

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Io;

import :GeometryArchive;
import :Records;

namespace Stereon {

namespace {

// ---------------------------------------------------------------------------
// Layout constants (format version 1)

constexpr std::array<std::byte, 4> Magic{std::byte{'S'}, std::byte{'T'}, std::byte{'N'},
                                         std::byte{0x1A}};
constexpr std::size_t HeaderSize = 32;
constexpr std::uint32_t RequiredFlag = 1;
constexpr std::size_t ChunkHeaderSize = 4 + 4 + 8 + 4;

constexpr std::string_view PolicyTag = "POLI";
constexpr std::string_view CurvesTag = "CRVS";
constexpr std::string_view PlaneCurvesTag = "PCRV";
constexpr std::string_view SurfacesTag = "SRFS";

/// The record names of trimmed curves (format version 2).
template <std::size_t D>
constexpr std::string_view TrimmedName = D == 3 ? "TrimmedCurve" : "TrimmedCurve2d";

// ---------------------------------------------------------------------------
// Writing

void PutBytes(std::vector<std::byte>& out, const std::uint64_t value, const std::size_t count)
{
    for (std::size_t i = 0; i < count; ++i) {
        out.push_back(static_cast<std::byte>((value >> (8 * i)) & 0xFFU));
    }
}

void PutU16(std::vector<std::byte>& out, const std::uint16_t value)
{
    PutBytes(out, value, 2);
}

void PutU32(std::vector<std::byte>& out, const std::uint32_t value)
{
    PutBytes(out, value, 4);
}

void PutU64(std::vector<std::byte>& out, const std::uint64_t value)
{
    PutBytes(out, value, 8);
}

void PutF64(std::vector<std::byte>& out, const double value)
{
    PutU64(out, std::bit_cast<std::uint64_t>(value));
}

void PutTag(std::vector<std::byte>& out, const std::string_view tag)
{
    for (const char c : tag) {
        out.push_back(static_cast<std::byte>(c));
    }
}

void PutChunk(std::vector<std::byte>& out, const std::string_view tag, const bool required,
              const std::vector<std::byte>& payload)
{
    PutTag(out, tag);
    PutU32(out, required ? RequiredFlag : 0U);
    PutU64(out, payload.size());
    PutU32(out, Detail::Crc32(payload));
    out.insert(out.end(), payload.begin(), payload.end());
}

[[nodiscard]] bool PositiveFinite(const double x)
{
    return std::isfinite(x) && x > 0.0;
}

/// Every value finite and positive, in millimetres, and the resolution at
/// most the tolerance ceiling, as Context requires of its settings.
[[nodiscard]] bool ValidPolicy(const ModelPolicy& p)
{
    return p.Unit == LengthUnit::Millimetre && PositiveFinite(p.Resolution) &&
           PositiveFinite(p.AngularResolution) && PositiveFinite(p.ToleranceCeiling) &&
           PositiveFinite(p.FrameAllowance) && p.Resolution <= p.ToleranceCeiling;
}

/// The factories that validate geometry, on write and on read, check frames
/// against the kernel's FrameAllowance and no other; a policy that states
/// another would claim a validation that never happened.
[[nodiscard]] bool AppliedFrameAllowance(const ModelPolicy& p)
{
    return p.FrameAllowance == FrameAllowance;
}

void PutHeader(std::vector<std::byte>& out, const ModelPolicy& policy)
{
    out.insert(out.end(), Magic.begin(), Magic.end());
    PutU32(out, GeometryFormatVersion);
    PutU16(out, static_cast<std::uint16_t>(LibraryVersion.Major));
    PutU16(out, static_cast<std::uint16_t>(LibraryVersion.Minor));
    PutU16(out, static_cast<std::uint16_t>(LibraryVersion.Patch));
    PutU16(out, 0);
    PutU32(out, static_cast<std::uint32_t>(policy.Unit));
    PutF64(out, policy.Resolution);
    PutU32(out, 0);
}

[[nodiscard]] std::vector<std::byte> PolicyPayload(const ModelPolicy& policy)
{
    std::vector<std::byte> payload;
    PutF64(payload, policy.AngularResolution);
    PutF64(payload, policy.ToleranceCeiling);
    PutF64(payload, policy.FrameAllowance);
    return payload;
}

/// One record: the type's name, the value's length and the value.
void PutNamed(std::vector<std::byte>& out, const std::string_view name,
              const std::vector<std::byte>& value)
{
    PutU16(out, static_cast<std::uint16_t>(name.size()));
    PutTag(out, name);
    PutU64(out, value.size());
    out.insert(out.end(), value.begin(), value.end());
}

template <typename T>
void PutRecord(std::vector<std::byte>& out, const T& value)
{
    PutNamed(out, GeometryTypeName<T>, Detail::RecordBytes(value));
}

template <typename Wrapper>
[[nodiscard]] bool PutItem(std::vector<std::byte>& out, const Wrapper& item, int depth);

/// A composite record: an inner item's record (a trimmed curve's or
/// surface's basis, a sweep's profile) followed by `tail` as Serialized
/// writes it. False when the inner item cannot be stored or the nesting is
/// already as deep as the format allows.
template <typename Inner, typename Tail>
[[nodiscard]] bool PutComposite(std::vector<std::byte>& out, const std::string_view name,
                                const Inner& inner, const Tail& tail, const int depth)
{
    if (depth >= MaxTrimmedNesting) {
        return false;
    }
    std::vector<std::byte> value;
    if (!PutItem(value, inner, depth + 1)) {
        return false;
    }
    const std::vector<std::byte> bytes = Detail::RecordBytes(tail);
    value.insert(value.end(), bytes.begin(), bytes.end());
    PutNamed(out, name, value);
    return true;
}

/// The record of a trimmed surface or a sweep (format version 3), if the
/// item is one whose basis or profile can be stored.
[[nodiscard]] bool PutWrapped(std::vector<std::byte>& out, const Surface& item, const int depth)
{
    if (const auto* trimmed = item.Target<TrimmedSurface>()) {
        const Box2 range = trimmed->GetRange();
        return PutComposite(out, GeometryTypeName<TrimmedSurface>, trimmed->GetBasis(),
                            Detail::SurfaceTrimRecord{range.X.GetLower(), range.X.GetUpper(),
                                                      range.Y.GetLower(), range.Y.GetUpper(),
                                                      trimmed->IsReversed()},
                            depth);
    }
    if (const auto* extrusion = item.Target<ExtrusionSurface>()) {
        return PutComposite(out, GeometryTypeName<ExtrusionSurface>, extrusion->GetProfile(),
                            Detail::ExtrusionRecord{extrusion->GetDirection()}, depth);
    }
    if (const auto* revolution = item.Target<RevolutionSurface>()) {
        return PutComposite(
            out, GeometryTypeName<RevolutionSurface>, revolution->GetProfile(),
            Detail::RevolutionRecord{revolution->GetOrigin(), revolution->GetAxis()}, depth);
    }
    return false;
}

/// The record of a trimmed curve (format version 2), if the item is one
/// whose basis can be stored.
template <std::size_t D>
[[nodiscard]] bool PutWrapped(std::vector<std::byte>& out, const BasicErasedCurve<D>& item,
                              const int depth)
{
    const auto* trimmed = item.template Target<BasicTrimmedCurve<D>>();
    if (trimmed == nullptr) {
        return false;
    }
    const Interval range = trimmed->GetRange();
    return PutComposite(
        out, TrimmedName<D>, trimmed->GetBasis(),
        Detail::TrimRecord{range.GetLower(), range.GetUpper(), trimmed->IsReversed()}, depth);
}

/// The item's record: a built-in type's, or a composite one (a trimmed
/// curve, a trimmed surface, a sweep). False for anything else, which the
/// format has no record for.
template <typename Wrapper>
bool PutItem(std::vector<std::byte>& out, const Wrapper& item, const int depth)
{
    const bool builtin = VisitBuiltin(item, [&]<typename G>(const G& geometry) {
        if constexpr (std::is_same_v<G, Wrapper>) {
            return false;
        } else {
            PutRecord(out, geometry);
            return true;
        }
    });
    return builtin || PutWrapped(out, item, depth);
}

/// The payload of a geometry chunk: the record count and the records.
/// NotSupported for a wrapped type that is not built in.
template <typename Wrapper>
[[nodiscard]] Result<std::vector<std::byte>> GeometryPayload(const std::vector<Wrapper>& items,
                                                             const std::string_view what)
{
    std::vector<std::byte> payload;
    PutU64(payload, items.size());
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (items[i].IsEmpty()) {
            return Fail(ErrorCode::InvalidInput,
                        std::format("WriteGeometry: {} {} is empty", what, i));
        }
        if (!PutItem(payload, items[i], 0)) {
            return Fail(ErrorCode::NotSupported,
                        std::format("WriteGeometry: {} {} is of a type the archive has no "
                                    "record for; only the built-in types, trimmed curves and "
                                    "surfaces of them and sweeps of stored profiles, nested "
                                    "at most {} deep, are stored",
                                    what, i, MaxTrimmedNesting));
        }
    }
    return payload;
}

// ---------------------------------------------------------------------------
// Reading

/// Reads bytes in order; every read checks what is left.
class Cursor {
public:
    explicit Cursor(const std::span<const std::byte> bytes)
        : m_Bytes(bytes)
    {}

    [[nodiscard]] std::size_t Offset() const { return m_Offset; }
    [[nodiscard]] std::size_t Left() const { return m_Bytes.size() - m_Offset; }
    [[nodiscard]] bool AtEnd() const { return m_Offset == m_Bytes.size(); }

    [[nodiscard]] Result<std::span<const std::byte>> Take(const std::size_t count,
                                                          const std::string_view what)
    {
        if (count > Left()) {
            return Fail(ErrorCode::ParseError,
                        std::format("ReadGeometry: truncated {}: {} bytes needed at offset {}, "
                                    "{} left",
                                    what, count, m_Offset, Left()));
        }
        const std::span<const std::byte> taken = m_Bytes.subspan(m_Offset, count);
        m_Offset += count;
        return taken;
    }

    [[nodiscard]] Result<std::uint64_t> Unsigned(const std::size_t count,
                                                 const std::string_view what)
    {
        STEREON_TRY(const std::span<const std::byte> taken, Take(count, what));
        std::uint64_t value = 0;
        for (std::size_t i = 0; i < count; ++i) {
            value |= std::uint64_t{std::to_integer<std::uint8_t>(taken[i])} << (8 * i);
        }
        return value;
    }

    [[nodiscard]] Result<double> Double(const std::string_view what)
    {
        STEREON_TRY(const std::uint64_t bits, Unsigned(8, what));
        return std::bit_cast<double>(bits);
    }

private:
    std::span<const std::byte> m_Bytes;
    std::size_t m_Offset = 0;
};

[[nodiscard]] std::string TagOf(const std::span<const std::byte> bytes)
{
    std::string tag;
    for (const std::byte b : bytes) {
        const auto c = std::to_integer<unsigned char>(b);
        tag += (c >= 0x20 && c < 0x7F) ? static_cast<char>(c) : '?';
    }
    return tag;
}

/// The header's fields.
struct Header {
    std::uint32_t FormatVersion = 0;
    Version WrittenBy{0, 0, 0};
    std::uint32_t Unit = 0;
    double Resolution = 0.0;
};

[[nodiscard]] Result<Header> ReadHeader(Cursor& cursor)
{
    STEREON_TRY(const std::span<const std::byte> magic, cursor.Take(Magic.size(), "header"));
    if (!std::ranges::equal(magic, Magic)) {
        return Fail(ErrorCode::ParseError, "ReadGeometry: not a Stereon archive (bad magic)");
    }
    Header h;
    STEREON_TRY(const std::uint64_t version, cursor.Unsigned(4, "header"));
    h.FormatVersion = static_cast<std::uint32_t>(version);
    if (h.FormatVersion == 0 || h.FormatVersion > GeometryFormatVersion) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: format version {} is not one this reader knows "
                                "(it reads 1 to {})",
                                h.FormatVersion, GeometryFormatVersion));
    }
    STEREON_TRY(const std::uint64_t major, cursor.Unsigned(2, "header"));
    STEREON_TRY(const std::uint64_t minor, cursor.Unsigned(2, "header"));
    STEREON_TRY(const std::uint64_t patch, cursor.Unsigned(2, "header"));
    h.WrittenBy = {static_cast<int>(major), static_cast<int>(minor), static_cast<int>(patch)};
    STEREON_TRY_VOID(cursor.Take(2, "header"));
    STEREON_TRY(const std::uint64_t unit, cursor.Unsigned(4, "header"));
    h.Unit = static_cast<std::uint32_t>(unit);
    STEREON_TRY(h.Resolution, cursor.Double("header"));
    STEREON_TRY_VOID(cursor.Take(4, "header"));
    return h;
}

/// One chunk as stored: its tag, whether it is required, and its payload,
/// whose checksum has been verified.
struct Chunk {
    std::string Tag;
    bool Required = false;
    std::span<const std::byte> Payload;
    std::size_t Offset = 0;
};

[[nodiscard]] Result<Chunk> ReadChunk(Cursor& cursor)
{
    Chunk chunk;
    chunk.Offset = cursor.Offset();
    STEREON_TRY(const std::span<const std::byte> tag, cursor.Take(4, "chunk header"));
    chunk.Tag = TagOf(tag);
    STEREON_TRY(const std::uint64_t flags, cursor.Unsigned(4, "chunk header"));
    chunk.Required = (flags & RequiredFlag) != 0;
    STEREON_TRY(const std::uint64_t length, cursor.Unsigned(8, "chunk header"));
    STEREON_TRY(const std::uint64_t crc, cursor.Unsigned(4, "chunk header"));
    if (length > cursor.Left()) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: chunk \"{}\" at offset {} claims {} bytes, {} left",
                                chunk.Tag, chunk.Offset, length, cursor.Left()));
    }
    STEREON_TRY(chunk.Payload, cursor.Take(static_cast<std::size_t>(length), "chunk"));
    if (Detail::Crc32(chunk.Payload) != crc) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: chunk \"{}\" at offset {} is damaged (checksum "
                                "mismatch)",
                                chunk.Tag, chunk.Offset));
    }
    return chunk;
}

/// The records of a geometry chunk: name and serialized value each.
struct RawRecord {
    std::string Name;
    std::span<const std::byte> Value;
};

/// One record at the cursor.
[[nodiscard]] Result<RawRecord> ReadRecord(Cursor& cursor)
{
    STEREON_TRY(const std::uint64_t nameLength, cursor.Unsigned(2, "record"));
    STEREON_TRY(const std::span<const std::byte> name,
                cursor.Take(static_cast<std::size_t>(nameLength), "record"));
    STEREON_TRY(const std::uint64_t valueLength, cursor.Unsigned(8, "record"));
    STEREON_TRY(const std::span<const std::byte> value,
                cursor.Take(static_cast<std::size_t>(valueLength), "record"));
    return RawRecord{TagOf(name), value};
}

/// A composite record split into its inner record and its tail: a trimmed
/// curve's or surface's basis and trim, a sweep's profile and sweep.
struct RawComposite {
    RawRecord Inner;
    std::span<const std::byte> Tail;
};

[[nodiscard]] Result<RawComposite> SplitComposite(const RawRecord& record)
{
    Cursor cursor(record.Value);
    STEREON_TRY(RawRecord inner, ReadRecord(cursor));
    STEREON_TRY(const std::span<const std::byte> tail, cursor.Take(cursor.Left(), "tail"));
    return RawComposite{std::move(inner), tail};
}

[[nodiscard]] Result<std::vector<RawRecord>> ReadRecords(const Chunk& chunk)
{
    Cursor cursor(chunk.Payload);
    STEREON_TRY(const std::uint64_t count, cursor.Unsigned(8, "record count"));
    // Every record takes at least its 10 bytes of name length and value length.
    if (count > cursor.Left() / 10) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: chunk \"{}\" claims {} records in {} bytes",
                                chunk.Tag, count, cursor.Left()));
    }
    std::vector<RawRecord> records;
    for (std::uint64_t i = 0; i < count; ++i) {
        STEREON_TRY(RawRecord record, ReadRecord(cursor));
        records.push_back(std::move(record));
    }
    if (!cursor.AtEnd()) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: chunk \"{}\" has {} bytes after its records",
                                chunk.Tag, cursor.Left()));
    }
    return records;
}

// ---------------------------------------------------------------------------
// Validation: every value read goes through its type's checked factory
// (ADR-0016), as any data from outside the kernel does.

template <std::size_t D>
[[nodiscard]] auto Validated(const BasicLine<double, D>& g)
{
    return MakeLine(g.Origin, g.Direction);
}
template <std::size_t D>
[[nodiscard]] auto Validated(const BasicCircle<double, D>& g)
{
    return MakeCircle(g.Center, g.Radius, g.XAxis, g.YAxis);
}
template <std::size_t D>
[[nodiscard]] auto Validated(const BasicEllipse<double, D>& g)
{
    return MakeEllipse(g.Center, g.MajorRadius, g.MinorRadius, g.XAxis, g.YAxis);
}
template <std::size_t D>
[[nodiscard]] auto Validated(const BasicHyperbola<double, D>& g)
{
    return MakeHyperbola(g.Center, g.MajorRadius, g.MinorRadius, g.XAxis, g.YAxis);
}
template <std::size_t D>
[[nodiscard]] auto Validated(const BasicParabola<double, D>& g)
{
    return MakeParabola(g.Vertex, g.FocalLength, g.XAxis, g.YAxis);
}
template <std::size_t D>
[[nodiscard]] auto Validated(const BasicBSplineCurve<double, D>& g)
{
    return MakeBSplineCurve(g.Degree, g.Knots, g.Poles, g.Weights);
}
[[nodiscard]] auto Validated(const Plane& g)
{
    return MakePlane(g.Origin, g.XAxis, g.YAxis);
}
[[nodiscard]] auto Validated(const Cylinder& g)
{
    return MakeCylinder(g.Origin, g.Radius, g.XAxis, g.YAxis, g.RightHanded);
}
[[nodiscard]] auto Validated(const Cone& g)
{
    return MakeCone(g.Origin, g.Radius, g.SemiAngle, g.XAxis, g.YAxis, g.RightHanded);
}
[[nodiscard]] auto Validated(const Sphere& g)
{
    return MakeSphere(g.Center, g.Radius, g.XAxis, g.YAxis, g.RightHanded);
}
[[nodiscard]] auto Validated(const Torus& g)
{
    return MakeTorus(g.Center, g.MajorRadius, g.MinorRadius, g.XAxis, g.YAxis, g.RightHanded);
}
[[nodiscard]] auto Validated(const BSplineSurface& g)
{
    return MakeBSplineSurface(g.DegreeU, g.DegreeV, g.KnotsU, g.KnotsV, g.Poles, g.Weights);
}

/// One record decoded as T and validated, in a Wrapper.
template <typename T, typename Wrapper>
[[nodiscard]] Result<Wrapper> Decode(const RawRecord& record, const std::string& where)
{
    const Result<T> value = Detail::RecordValue<T>(record.Value);
    if (!value) {
        return Fail(ErrorCode::ParseError, std::format("ReadGeometry: {} ({}): {}", where,
                                                       record.Name, value.error().GetMessage()));
    }
    const auto valid = Validated(*value);
    if (!valid) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} ({}) is not valid geometry: {}", where,
                                record.Name, valid.error().GetMessage()));
    }
    return Wrapper(*valid);
}

template <typename Wrapper, typename... Ts>
[[nodiscard]] Result<Wrapper> DecodeAny(TypeList<Ts...> types, const RawRecord& record,
                                        const std::string& where, int depth, std::uint32_t version);

/// The format versions that introduced the composite records.
constexpr std::uint32_t TrimmedCurvesVersion = 2;
constexpr std::uint32_t SurfaceRecordsVersion = 3;

/// ParseError unless the file's format version has the record's type: no
/// writer of version 1 stored a trimmed curve, so a version-1 file with one
/// was not written as it claims.
[[nodiscard]] Result<void> RequireVersion(const RawRecord& record, const std::string& where,
                                          const std::uint32_t since, const std::uint32_t version)
{
    if (version < since) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} is a \"{}\" record, which format version {} "
                                "does not have (it came in version {})",
                                where, record.Name, version, since));
    }
    return {};
}

/// A composite record's parts, unless it nests deeper than the format
/// allows or does not split: ParseError naming the record either way.
[[nodiscard]] Result<RawComposite> Parts(const RawRecord& record, const std::string& where,
                                         const int depth)
{
    if (depth >= MaxTrimmedNesting) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} nests records more than {} deep", where,
                                MaxTrimmedNesting));
    }
    Result<RawComposite> parts = SplitComposite(record);
    if (!parts) {
        return Fail(ErrorCode::ParseError, std::format("ReadGeometry: {} ({}): {}", where,
                                                       record.Name, parts.error().GetMessage()));
    }
    return parts;
}

/// The composite record's tail decoded as T, or a ParseError naming the
/// record.
template <typename T>
[[nodiscard]] Result<T> TailValue(const RawComposite& parts, const RawRecord& record,
                                  const std::string& where)
{
    Result<T> value = Detail::RecordValue<T>(parts.Tail);
    if (!value) {
        return Fail(ErrorCode::ParseError, std::format("ReadGeometry: {} ({}): {}", where,
                                                       record.Name, value.error().GetMessage()));
    }
    return value;
}

/// Checked before an Interval is made: its constructor requires
/// lower <= upper, and NaN satisfies nothing.
[[nodiscard]] bool FiniteNonEmpty(const double lower, const double upper)
{
    return std::isfinite(lower) && std::isfinite(upper) && lower < upper;
}

template <typename Wrapper>
[[nodiscard]] Result<Wrapper> RangeError(const RawRecord& record, const std::string& where)
{
    return Fail(ErrorCode::ParseError,
                std::format("ReadGeometry: {} ({}) has a range that is not finite and non-empty",
                            where, record.Name));
}

/// The checked factory's result as a Wrapper, or its refusal as a
/// ParseError naming the record.
template <typename Wrapper, typename T>
[[nodiscard]] Result<Wrapper> AsRead(const Result<T>& valid, const RawRecord& record,
                                     const std::string& where)
{
    if (!valid) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} ({}) is not valid geometry: {}", where,
                                record.Name, valid.error().GetMessage()));
    }
    return Wrapper(*valid);
}

/// A trimmed curve's record: its basis decoded as any record of the chunk,
/// then its trim, validated by MakeTrimmedCurve.
template <std::size_t D, typename Types>
[[nodiscard]] Result<BasicErasedCurve<D>> DecodeTrimmed(Types types, const RawRecord& record,
                                                        const std::string& where, const int depth,
                                                        const std::uint32_t version)
{
    STEREON_TRY(const RawComposite parts, Parts(record, where, depth));
    STEREON_TRY(const BasicErasedCurve<D> basis,
                DecodeAny<BasicErasedCurve<D>>(types, parts.Inner, where + ", its basis", depth + 1,
                                               version));
    STEREON_TRY(const Detail::TrimRecord trim, TailValue<Detail::TrimRecord>(parts, record, where));
    if (!FiniteNonEmpty(trim.Lower, trim.Upper)) {
        return RangeError<BasicErasedCurve<D>>(record, where);
    }
    return AsRead<BasicErasedCurve<D>>(MakeTrimmedCurve(basis, Interval{trim.Lower, trim.Upper})
                                           .transform([&](const BasicTrimmedCurve<D>& t) {
                                               return trim.Reversed ? t.Reversed() : t;
                                           }),
                                       record, where);
}

/// A trimmed surface's record (format version 3): its basis decoded as any
/// record of the surface chunk, then its trim, validated by
/// MakeTrimmedSurface.
template <typename Types>
[[nodiscard]] Result<Surface> DecodeTrimmedSurface(Types types, const RawRecord& record,
                                                   const std::string& where, const int depth,
                                                   const std::uint32_t version)
{
    STEREON_TRY(const RawComposite parts, Parts(record, where, depth));
    STEREON_TRY(const Surface basis,
                DecodeAny<Surface>(types, parts.Inner, where + ", its basis", depth + 1, version));
    STEREON_TRY(const Detail::SurfaceTrimRecord trim,
                TailValue<Detail::SurfaceTrimRecord>(parts, record, where));
    if (!FiniteNonEmpty(trim.LowerU, trim.UpperU) || !FiniteNonEmpty(trim.LowerV, trim.UpperV)) {
        return RangeError<Surface>(record, where);
    }
    const Box2 box{Interval{trim.LowerU, trim.UpperU}, Interval{trim.LowerV, trim.UpperV}};
    return AsRead<Surface>(MakeTrimmedSurface(basis, box).transform([&](const TrimmedSurface& t) {
        return trim.Reversed ? t.Reversed() : t;
    }),
                           record, where);
}

/// A sweep's record (format version 3): its profile decoded as any record
/// of the curve chunk, then its sweep, validated by `make`.
template <typename Record, typename Make>
[[nodiscard]] Result<Surface> DecodeSwept(const RawRecord& record, const std::string& where,
                                          const int depth, const std::uint32_t version, Make make)
{
    STEREON_TRY(const RawComposite parts, Parts(record, where, depth));
    STEREON_TRY(const Curve profile, DecodeAny<Curve>(BuiltinCurveTypes{}, parts.Inner,
                                                      where + ", its profile", depth + 1, version));
    STEREON_TRY(const Record sweep, TailValue<Record>(parts, record, where));
    return AsRead<Surface>(make(profile, sweep), record, where);
}

/// A revolution restored from its record: validated by MakeRevolution, then
/// built with the stored axis as it is. The writer stores the unit axis the
/// surface keeps, and normalising it again could move it by an ulp and
/// break the bit-for-bit round trip; an axis that is not of unit length is
/// not a surface's state and is refused.
[[nodiscard]] Result<RevolutionSurface> RestoredRevolution(const Curve& profile,
                                                           const Detail::RevolutionRecord& r)
{
    STEREON_TRY_VOID(MakeRevolution(profile, r.Origin, r.Axis));
    if (std::abs(Dot(r.Axis, r.Axis) - 1.0) > FrameAllowance) {
        return Fail(ErrorCode::InvalidInput, "the axis is not of unit length");
    }
    return RevolutionSurface(profile, r.Origin, r.Axis, RevolutionSurface::UnitAxis{});
}

/// The record decoded as the type of the list it names, or as a composite
/// record of the chunk's kind: a trimmed curve in a curve chunk, a trimmed
/// surface or a sweep in the surface chunk.
template <typename Wrapper, typename... Ts>
Result<Wrapper> DecodeAny(TypeList<Ts...> types, const RawRecord& record, const std::string& where,
                          const int depth, const std::uint32_t version)
{
    if constexpr (std::is_same_v<Wrapper, Curve> || std::is_same_v<Wrapper, Curve2d>) {
        constexpr std::size_t D = std::is_same_v<Wrapper, Curve> ? 3 : 2;
        if (record.Name == TrimmedName<D>) {
            STEREON_TRY_VOID(RequireVersion(record, where, TrimmedCurvesVersion, version));
            return DecodeTrimmed<D>(types, record, where, depth, version);
        }
    }
    if constexpr (std::is_same_v<Wrapper, Surface>) {
        const bool composite = record.Name == GeometryTypeName<TrimmedSurface> ||
                               record.Name == GeometryTypeName<ExtrusionSurface> ||
                               record.Name == GeometryTypeName<RevolutionSurface>;
        if (composite) {
            STEREON_TRY_VOID(RequireVersion(record, where, SurfaceRecordsVersion, version));
        }
        if (record.Name == GeometryTypeName<TrimmedSurface>) {
            return DecodeTrimmedSurface(types, record, where, depth, version);
        }
        if (record.Name == GeometryTypeName<ExtrusionSurface>) {
            return DecodeSwept<Detail::ExtrusionRecord>(
                record, where, depth, version,
                [](const Curve& profile, const Detail::ExtrusionRecord& r) {
                    return MakeExtrusion(profile, r.Direction);
                });
        }
        if (record.Name == GeometryTypeName<RevolutionSurface>) {
            return DecodeSwept<Detail::RevolutionRecord>(record, where, depth, version,
                                                         RestoredRevolution);
        }
    }
    std::optional<Result<Wrapper>> decoded;
    static_cast<void>(((record.Name == GeometryTypeName<Ts>
                            ? (decoded.emplace(Decode<Ts, Wrapper>(record, where)), true)
                            : false) ||
                       ...));
    if (!decoded) {
        return Fail(
            ErrorCode::ParseError,
            std::format("ReadGeometry: {} is of an unknown type \"{}\"", where, record.Name));
    }
    return *std::move(decoded);
}

template <typename Wrapper, typename Types>
[[nodiscard]] Result<std::vector<Wrapper>> DecodeChunk(const Chunk& chunk, Types types,
                                                       const std::uint32_t version)
{
    STEREON_TRY(const std::vector<RawRecord> records, ReadRecords(chunk));
    std::vector<Wrapper> items;
    items.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
        const std::string where = std::format("record {} of chunk \"{}\"", i, chunk.Tag);
        STEREON_TRY(Wrapper item, DecodeAny<Wrapper>(types, records[i], where, 0, version));
        items.push_back(std::move(item));
    }
    return items;
}

[[nodiscard]] Result<void> ReadPolicy(const Chunk& chunk, ModelPolicy& policy)
{
    Cursor cursor(chunk.Payload);
    STEREON_TRY(policy.AngularResolution, cursor.Double("policy"));
    STEREON_TRY(policy.ToleranceCeiling, cursor.Double("policy"));
    STEREON_TRY(policy.FrameAllowance, cursor.Double("policy"));
    if (!cursor.AtEnd()) {
        return Fail(ErrorCode::ParseError, "ReadGeometry: the policy chunk is too long");
    }
    return {};
}

/// Which of the chunks every archive has one of were seen. Each must appear
/// exactly once: a missing one would read as an empty collection and a
/// repeated one would replace what the first decoded, losing geometry with
/// every checksum intact.
struct ChunkState {
    bool Policy = false;
    bool Curves = false;
    bool PlaneCurves = false;
    bool Surfaces = false;
};

/// Marks the chunk's kind seen; ParseError if it was already.
[[nodiscard]] Result<void> FirstOfKind(const Chunk& chunk, bool& seen)
{
    if (seen) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: a second \"{}\" chunk at offset {}; an archive has "
                                "exactly one",
                                chunk.Tag, chunk.Offset));
    }
    seen = true;
    return {};
}

[[nodiscard]] Result<void> ApplyChunk(const Chunk& chunk, LoadedGeometry& loaded, ChunkState& state)
{
    GeometryArchive& archive = loaded.Archive;
    const std::uint32_t version = loaded.FormatVersion;
    if (chunk.Tag == PolicyTag) {
        STEREON_TRY_VOID(FirstOfKind(chunk, state.Policy));
        return ReadPolicy(chunk, archive.Policy);
    }
    if (chunk.Tag == CurvesTag) {
        STEREON_TRY_VOID(FirstOfKind(chunk, state.Curves));
        STEREON_TRY(archive.Curves, DecodeChunk<Curve>(chunk, BuiltinCurveTypes{}, version));
        return {};
    }
    if (chunk.Tag == PlaneCurvesTag) {
        STEREON_TRY_VOID(FirstOfKind(chunk, state.PlaneCurves));
        STEREON_TRY(archive.Curves2d, DecodeChunk<Curve2d>(chunk, BuiltinCurve2dTypes{}, version));
        return {};
    }
    if (chunk.Tag == SurfacesTag) {
        STEREON_TRY_VOID(FirstOfKind(chunk, state.Surfaces));
        STEREON_TRY(archive.Surfaces, DecodeChunk<Surface>(chunk, BuiltinSurfaceTypes{}, version));
        return {};
    }
    if (chunk.Required) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: unknown required chunk \"{}\" at offset {}; a "
                                "newer writer stored data this reader cannot interpret",
                                chunk.Tag, chunk.Offset));
    }
    loaded.Warnings.push_back(std::format("skipped unknown optional chunk \"{}\" ({} bytes)",
                                          chunk.Tag, chunk.Payload.size()));
    return {};
}

[[nodiscard]] Result<void> CheckPolicy(const ModelPolicy& policy)
{
    if (!ValidPolicy(policy)) {
        return Fail(ErrorCode::ParseError,
                    "ReadGeometry: the stored policy has a unit this reader does not know, a "
                    "value that is not finite and positive, or a resolution above its tolerance "
                    "ceiling");
    }
    if (!AppliedFrameAllowance(policy)) {
        // ADR-0028: never claim to restore a policy the reader did not apply.
        return Fail(ErrorCode::NotSupported,
                    std::format("ReadGeometry: the stored frame allowance {} is not the one this "
                                "reader validates frames against, FrameAllowance ({})",
                                policy.FrameAllowance, FrameAllowance));
    }
    return {};
}

[[nodiscard]] std::string Hex(const std::uint32_t x)
{
    return std::format("{:#010x}", x);
}

} // namespace

ModelPolicy PolicyOf(const Context& context)
{
    return {LengthUnit::Millimetre, context.GetResolution(), context.GetAngularResolution(),
            context.GetToleranceCeiling(), FrameAllowance};
}

Result<std::vector<std::byte>> WriteGeometry(const GeometryArchive& archive)
{
    if constexpr (!SerializationAvailable) {
        return Fail(ErrorCode::NotSupported, "WriteGeometry: needs static reflection (ADR-0016)");
    }
    if (!ValidPolicy(archive.Policy)) {
        return Fail(ErrorCode::InvalidInput,
                    "WriteGeometry: every policy value must be finite and positive, in "
                    "millimetres and radians, and the resolution at most the tolerance ceiling");
    }
    if (!AppliedFrameAllowance(archive.Policy)) {
        return Fail(ErrorCode::InvalidInput,
                    std::format("WriteGeometry: the frame allowance {} is not the one this kernel "
                                "validates frames against, FrameAllowance ({})",
                                archive.Policy.FrameAllowance, FrameAllowance));
    }
    STEREON_TRY(const std::vector<std::byte> curves,
                GeometryPayload(archive.Curves, std::string_view{"curve"}));
    STEREON_TRY(const std::vector<std::byte> planeCurves,
                GeometryPayload(archive.Curves2d, std::string_view{"plane curve"}));
    STEREON_TRY(const std::vector<std::byte> surfaces,
                GeometryPayload(archive.Surfaces, std::string_view{"surface"}));
    std::vector<std::byte> out;
    PutHeader(out, archive.Policy);
    PutChunk(out, PolicyTag, true, PolicyPayload(archive.Policy));
    PutChunk(out, CurvesTag, true, curves);
    PutChunk(out, PlaneCurvesTag, true, planeCurves);
    PutChunk(out, SurfacesTag, true, surfaces);
    return out;
}

Result<LoadedGeometry> ReadGeometry(const std::span<const std::byte> bytes)
{
    if constexpr (!SerializationAvailable) {
        return Fail(ErrorCode::NotSupported, "ReadGeometry: needs static reflection (ADR-0016)");
    }
    Cursor cursor(bytes);
    STEREON_TRY(const Header header, ReadHeader(cursor));
    LoadedGeometry loaded;
    loaded.FormatVersion = header.FormatVersion;
    loaded.WrittenBy = header.WrittenBy;
    loaded.Archive.Policy.Unit = static_cast<LengthUnit>(header.Unit);
    loaded.Archive.Policy.Resolution = header.Resolution;
    ChunkState state;
    while (!cursor.AtEnd()) {
        STEREON_TRY(const Chunk chunk, ReadChunk(cursor));
        STEREON_TRY_VOID(ApplyChunk(chunk, loaded, state));
    }
    if (!state.Policy) {
        // ADR-0028: never substitute the reader's defaults for missing policy.
        return Fail(ErrorCode::ParseError, "ReadGeometry: the archive has no policy chunk");
    }
    // Every writer stores all three geometry chunks, empty or not, so one
    // that is missing was cut off or removed, not left out.
    for (const auto& [seen, tag] :
         {std::pair{state.Curves, CurvesTag}, std::pair{state.PlaneCurves, PlaneCurvesTag},
          std::pair{state.Surfaces, SurfacesTag}}) {
        if (!seen) {
            return Fail(ErrorCode::ParseError,
                        std::format("ReadGeometry: the archive has no \"{}\" chunk; it was cut "
                                    "short or the chunk removed",
                                    tag));
        }
    }
    STEREON_TRY_VOID(CheckPolicy(loaded.Archive.Policy));
    return loaded;
}

Result<void> SaveGeometry(const std::filesystem::path& path, const GeometryArchive& archive)
{
    STEREON_TRY(const std::vector<std::byte> bytes, WriteGeometry(archive));
    // Written beside the target and renamed over it only once complete and
    // closed: a failed save leaves an existing file as it was, and an error
    // the stream reports only when it flushes (a full disk) is reported
    // rather than lost in its destructor.
    std::filesystem::path partial = path;
    partial += ".partial";
    std::error_code ignored;
    {
        std::ofstream file(partial, std::ios::binary | std::ios::trunc);
        if (!file.is_open()) {
            // Nothing was created, so nothing is removed.
            return Fail(ErrorCode::InvalidInput,
                        std::format("SaveGeometry: cannot write \"{}\"", path.string()));
        }
        const bool written = Detail::WriteAll(file, bytes);
        file.close();
        if (!written || file.fail()) {
            std::filesystem::remove(partial, ignored);
            return Fail(ErrorCode::InvalidInput,
                        std::format("SaveGeometry: cannot write \"{}\"", path.string()));
        }
    }
    std::error_code error;
    std::filesystem::rename(partial, path, error);
    if (error) {
        std::filesystem::remove(partial, ignored);
        return Fail(ErrorCode::InvalidInput, std::format("SaveGeometry: cannot replace \"{}\": {}",
                                                         path.string(), error.message()));
    }
    return {};
}

Result<LoadedGeometry> LoadGeometry(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return Fail(ErrorCode::InvalidInput,
                    std::format("LoadGeometry: cannot read \"{}\"", path.string()));
    }
    const std::vector<char> chars{std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>()};
    std::vector<std::byte> bytes(chars.size());
    std::ranges::transform(chars, bytes.begin(), [](const char c) {
        return static_cast<std::byte>(static_cast<unsigned char>(c));
    });
    return ReadGeometry(bytes);
}

namespace {

template <typename... Ts>
[[nodiscard]] std::string RecordJsonOf(TypeList<Ts...> types, const RawRecord& record, int depth);

/// The JSON of a composite record's value: the inner record's, decoded
/// with `innerTypes`, under `innerKey`, and the tail's under `tailKey`.
template <typename Tail, typename InnerTypes>
[[nodiscard]] std::string CompositeJson(const RawRecord& record, InnerTypes innerTypes,
                                        const std::string_view innerKey,
                                        const std::string_view tailKey, const int depth)
{
    const Result<RawComposite> parts = SplitComposite(record);
    if (!parts) {
        return "null";
    }
    const std::string tail = Detail::RecordValue<Tail>(parts->Tail)
                                 .transform([](const Tail& t) { return Detail::RecordJson(t); })
                                 .value_or(std::string{"null"});
    return std::format("{{\"{}\":{},\"{}\":{}}}", innerKey,
                       RecordJsonOf(innerTypes, parts->Inner, depth + 1), tailKey, tail);
}

/// The JSON of a composite record's value when its name is one of the
/// chunk's composite kinds, else empty.
template <typename... Ts>
[[nodiscard]] std::string CompositeJsonOf(TypeList<Ts...> types, const RawRecord& record,
                                          const int depth)
{
    if (depth >= MaxTrimmedNesting) {
        return {};
    }
    if constexpr (std::is_same_v<TypeList<Ts...>, BuiltinCurveTypes>) {
        if (record.Name == TrimmedName<3>) {
            return CompositeJson<Detail::TrimRecord>(record, types, "basis", "trim", depth);
        }
    } else if constexpr (std::is_same_v<TypeList<Ts...>, BuiltinCurve2dTypes>) {
        if (record.Name == TrimmedName<2>) {
            return CompositeJson<Detail::TrimRecord>(record, types, "basis", "trim", depth);
        }
    } else if constexpr (std::is_same_v<TypeList<Ts...>, BuiltinSurfaceTypes>) {
        if (record.Name == GeometryTypeName<TrimmedSurface>) {
            return CompositeJson<Detail::SurfaceTrimRecord>(record, types, "basis", "trim", depth);
        }
        if (record.Name == GeometryTypeName<ExtrusionSurface>) {
            return CompositeJson<Detail::ExtrusionRecord>(record, BuiltinCurveTypes{}, "profile",
                                                          "sweep", depth);
        }
        if (record.Name == GeometryTypeName<RevolutionSurface>) {
            return CompositeJson<Detail::RevolutionRecord>(record, BuiltinCurveTypes{}, "profile",
                                                           "sweep", depth);
        }
    }
    return {};
}

/// The JSON of one record, decoded but not validated: the view is for
/// looking at a file. A composite record's value is its parts.
template <typename... Ts>
std::string RecordJsonOf(TypeList<Ts...> types, const RawRecord& record, const int depth)
{
    std::string value = CompositeJsonOf(types, record, depth);
    if (value.empty()) {
        value = "null";
        static_cast<void>(
            ((record.Name == GeometryTypeName<Ts>
                  ? (value = Detail::RecordValue<Ts>(record.Value)
                                 .transform([](const Ts& g) { return Detail::RecordJson(g); })
                                 .value_or(std::string{"null"}),
                     true)
                  : false) ||
             ...));
    }
    return std::format("{{\"type\":\"{}\",\"value\":{}}}", record.Name, value);
}

/// The JSON of one geometry chunk's records, decoded but not validated:
/// the view is for looking at a file, including one that does not load.
template <typename Types>
[[nodiscard]] Result<std::string> RecordsJson(Types types, const Chunk& chunk)
{
    STEREON_TRY(const std::vector<RawRecord> records, ReadRecords(chunk));
    std::string out = "[";
    for (std::size_t i = 0; i < records.size(); ++i) {
        out += (i > 0 ? "," : "") + RecordJsonOf(types, records[i], 0);
    }
    return out + "]";
}

[[nodiscard]] Result<std::string> ChunkJson(const Chunk& chunk)
{
    std::string body;
    if (chunk.Tag == CurvesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinCurveTypes{}, chunk));
    } else if (chunk.Tag == PlaneCurvesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinCurve2dTypes{}, chunk));
    } else if (chunk.Tag == SurfacesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinSurfaceTypes{}, chunk));
    } else {
        body = std::format("{{\"bytes\":{}}}", chunk.Payload.size());
    }
    return std::format("{{\"tag\":\"{}\",\"required\":{},\"crc\":\"{}\",\"content\":{}}}",
                       chunk.Tag, chunk.Required ? "true" : "false",
                       Hex(Detail::Crc32(chunk.Payload)), body);
}

} // namespace

Result<std::string> GeometryToJson(const std::span<const std::byte> bytes)
{
    STEREON_TRY(const LoadedGeometry loaded, ReadGeometry(bytes));
    const ModelPolicy& p = loaded.Archive.Policy;
    std::string out =
        std::format("{{\"format\":{},\"writtenBy\":\"{}.{}.{}\",\"policy\":{{\"unit\":\"mm\","
                    "\"resolution\":{},\"angularResolution\":{},\"toleranceCeiling\":{},"
                    "\"frameAllowance\":{}}},\"chunks\":[",
                    loaded.FormatVersion, loaded.WrittenBy.Major, loaded.WrittenBy.Minor,
                    loaded.WrittenBy.Patch, p.Resolution, p.AngularResolution, p.ToleranceCeiling,
                    p.FrameAllowance);
    Cursor cursor(bytes);
    STEREON_TRY_VOID(ReadHeader(cursor));
    for (bool first = true; !cursor.AtEnd(); first = false) {
        STEREON_TRY(const Chunk chunk, ReadChunk(cursor));
        STEREON_TRY(const std::string json, ChunkJson(chunk));
        out += (first ? "" : ",") + json;
    }
    return out + "]}";
}

} // namespace Stereon

namespace Stereon::Detail {

namespace {

[[nodiscard]] constexpr std::array<std::uint32_t, 256> Crc32Table() noexcept
{
    constexpr std::uint32_t polynomial = 0xEDB88320U;
    std::array<std::uint32_t, 256> table{};
    for (std::uint32_t i = 0; i < 256; ++i) {
        std::uint32_t c = i;
        for (int bit = 0; bit < 8; ++bit) {
            c = (c & 1U) != 0 ? polynomial ^ (c >> 1) : c >> 1;
        }
        table[i] = c;
    }
    return table;
}

constexpr std::array<std::uint32_t, 256> Crc32Lookup = Crc32Table();

} // namespace

bool WriteAll(std::ostream& out, const std::span<const std::byte> bytes)
{
    out.write(reinterpret_cast<const char*>(bytes.data()), // NOLINT: byte I/O
              static_cast<std::streamsize>(bytes.size()));
    out.flush();
    return !out.fail();
}

std::uint32_t Crc32(const std::span<const std::byte> bytes) noexcept
{
    std::uint32_t crc = 0xFFFFFFFFU;
    for (const std::byte b : bytes) {
        crc = Crc32Lookup[(crc ^ std::to_integer<std::uint32_t>(b)) & 0xFFU] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFFU;
}

} // namespace Stereon::Detail
