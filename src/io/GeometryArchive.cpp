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
#include <span>
#include <string>
#include <string_view>
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

[[nodiscard]] bool ValidPolicy(const ModelPolicy& p)
{
    return p.Unit == LengthUnit::Millimetre && PositiveFinite(p.Resolution) &&
           PositiveFinite(p.AngularResolution) && PositiveFinite(p.ToleranceCeiling) &&
           PositiveFinite(p.FrameAllowance);
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

/// A trimmed curve's record, if the item is one whose basis can be stored.
[[nodiscard]] bool PutTrimmed(std::vector<std::byte>& /*out*/, const Surface& /*item*/,
                              int /*depth*/)
{
    return false;
}

template <std::size_t D>
[[nodiscard]] bool PutTrimmed(std::vector<std::byte>& out, const BasicErasedCurve<D>& item,
                              const int depth)
{
    const auto* trimmed = item.template Target<BasicTrimmedCurve<D>>();
    if (trimmed == nullptr || depth >= MaxTrimmedNesting) {
        return false;
    }
    std::vector<std::byte> value;
    if (!PutItem(value, trimmed->GetBasis(), depth + 1)) {
        return false;
    }
    const Interval range = trimmed->GetRange();
    const std::vector<std::byte> trim = Detail::RecordBytes(
        Detail::TrimRecord{range.GetLower(), range.GetUpper(), trimmed->IsReversed()});
    value.insert(value.end(), trim.begin(), trim.end());
    PutNamed(out, TrimmedName<D>, value);
    return true;
}

/// The item's record: a built-in type's, or a trimmed curve's. False for
/// anything else, which the format has no record for.
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
    return builtin || PutTrimmed(out, item, depth);
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
                                    "record for; only the built-in types, and trimmed curves "
                                    "of them nested at most {} deep, are stored",
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

/// A trimmed curve's record split into its basis record and its trim bytes.
struct RawTrimmed {
    RawRecord Basis;
    std::span<const std::byte> Trim;
};

[[nodiscard]] Result<RawTrimmed> SplitTrimmed(const RawRecord& record)
{
    Cursor cursor(record.Value);
    STEREON_TRY(RawRecord basis, ReadRecord(cursor));
    STEREON_TRY(const std::span<const std::byte> trim, cursor.Take(cursor.Left(), "trim"));
    return RawTrimmed{std::move(basis), trim};
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
                                        const std::string& where, int depth);

/// A trimmed curve's record: its basis decoded as any record of the chunk,
/// then its trim, validated by MakeTrimmedCurve.
template <std::size_t D, typename Types>
[[nodiscard]] Result<BasicErasedCurve<D>> DecodeTrimmed(Types types, const RawRecord& record,
                                                        const std::string& where, const int depth)
{
    if (depth >= MaxTrimmedNesting) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} nests trimmed curves more than {} deep", where,
                                MaxTrimmedNesting));
    }
    const Result<RawTrimmed> split = SplitTrimmed(record);
    if (!split) {
        return Fail(ErrorCode::ParseError, std::format("ReadGeometry: {} ({}): {}", where,
                                                       record.Name, split.error().GetMessage()));
    }
    STEREON_TRY(
        const BasicErasedCurve<D> basis,
        DecodeAny<BasicErasedCurve<D>>(types, split->Basis, where + ", its basis", depth + 1));
    const Result<Detail::TrimRecord> trim = Detail::RecordValue<Detail::TrimRecord>(split->Trim);
    if (!trim) {
        return Fail(ErrorCode::ParseError, std::format("ReadGeometry: {} ({}): {}", where,
                                                       record.Name, trim.error().GetMessage()));
    }
    // Checked before the Interval is made: its constructor requires
    // lower <= upper, and NaN satisfies nothing.
    if (!std::isfinite(trim->Lower) || !std::isfinite(trim->Upper) ||
        !(trim->Lower < trim->Upper)) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} ({}) has a range that is not finite and "
                                "non-empty",
                                where, record.Name));
    }
    const auto valid = MakeTrimmedCurve(basis, Interval{trim->Lower, trim->Upper});
    if (!valid) {
        return Fail(ErrorCode::ParseError,
                    std::format("ReadGeometry: {} ({}) is not valid geometry: {}", where,
                                record.Name, valid.error().GetMessage()));
    }
    return BasicErasedCurve<D>(trim->Reversed ? valid->Reversed() : *valid);
}

/// The record decoded as the type of the list it names, or as a trimmed
/// curve in a curve chunk.
template <typename Wrapper, typename... Ts>
Result<Wrapper> DecodeAny(TypeList<Ts...> types, const RawRecord& record, const std::string& where,
                          const int depth)
{
    if constexpr (std::is_same_v<Wrapper, Curve> || std::is_same_v<Wrapper, Curve2d>) {
        constexpr std::size_t D = std::is_same_v<Wrapper, Curve> ? 3 : 2;
        if (record.Name == TrimmedName<D>) {
            return DecodeTrimmed<D>(types, record, where, depth);
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
[[nodiscard]] Result<std::vector<Wrapper>> DecodeChunk(const Chunk& chunk, Types types)
{
    STEREON_TRY(const std::vector<RawRecord> records, ReadRecords(chunk));
    std::vector<Wrapper> items;
    items.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
        const std::string where = std::format("record {} of chunk \"{}\"", i, chunk.Tag);
        STEREON_TRY(Wrapper item, DecodeAny<Wrapper>(types, records[i], where, 0));
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

/// What the chunks of a file set, and which required ones were seen.
struct ChunkState {
    bool Policy = false;
};

[[nodiscard]] Result<void> ApplyChunk(const Chunk& chunk, LoadedGeometry& loaded, ChunkState& state)
{
    GeometryArchive& archive = loaded.Archive;
    if (chunk.Tag == PolicyTag) {
        state.Policy = true;
        return ReadPolicy(chunk, archive.Policy);
    }
    if (chunk.Tag == CurvesTag) {
        STEREON_TRY(archive.Curves, DecodeChunk<Curve>(chunk, BuiltinCurveTypes{}));
        return {};
    }
    if (chunk.Tag == PlaneCurvesTag) {
        STEREON_TRY(archive.Curves2d, DecodeChunk<Curve2d>(chunk, BuiltinCurve2dTypes{}));
        return {};
    }
    if (chunk.Tag == SurfacesTag) {
        STEREON_TRY(archive.Surfaces, DecodeChunk<Surface>(chunk, BuiltinSurfaceTypes{}));
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
                    "ReadGeometry: the stored policy has a unit this reader does not know or a "
                    "value that is not finite and positive");
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
                    "millimetres and radians");
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
    STEREON_TRY_VOID(CheckPolicy(loaded.Archive.Policy));
    return loaded;
}

Result<void> SaveGeometry(const std::filesystem::path& path, const GeometryArchive& archive)
{
    STEREON_TRY(const std::vector<std::byte> bytes, WriteGeometry(archive));
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    file.write(reinterpret_cast<const char*>(bytes.data()), // NOLINT: byte I/O
               static_cast<std::streamsize>(bytes.size()));
    if (!file) {
        return Fail(ErrorCode::InvalidInput,
                    std::format("SaveGeometry: cannot write \"{}\"", path.string()));
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

/// The JSON of one record, decoded but not validated: the view is for
/// looking at a file. A trimmed curve's value is its basis and its trim;
/// `trimmedName` is empty for a chunk that holds none.
template <typename... Ts>
[[nodiscard]] std::string RecordJsonOf(TypeList<Ts...> types, const RawRecord& record,
                                       const std::string_view trimmedName, const int depth)
{
    std::string value = "null";
    if (!trimmedName.empty() && record.Name == trimmedName && depth < MaxTrimmedNesting) {
        const Result<RawTrimmed> split = SplitTrimmed(record);
        if (split) {
            const std::string trim =
                Detail::RecordValue<Detail::TrimRecord>(split->Trim)
                    .transform([](const Detail::TrimRecord& t) { return Detail::RecordJson(t); })
                    .value_or(std::string{"null"});
            value = std::format("{{\"basis\":{},\"trim\":{}}}",
                                RecordJsonOf(types, split->Basis, trimmedName, depth + 1), trim);
        }
    }
    static_cast<void>(
        ((record.Name == GeometryTypeName<Ts>
              ? (value = Detail::RecordValue<Ts>(record.Value)
                             .transform([](const Ts& g) { return Detail::RecordJson(g); })
                             .value_or(std::string{"null"}),
                 true)
              : false) ||
         ...));
    return std::format("{{\"type\":\"{}\",\"value\":{}}}", record.Name, value);
}

/// The JSON of one geometry chunk's records, decoded but not validated:
/// the view is for looking at a file, including one that does not load.
template <typename Types>
[[nodiscard]] Result<std::string> RecordsJson(Types types, const Chunk& chunk,
                                              const std::string_view trimmedName)
{
    STEREON_TRY(const std::vector<RawRecord> records, ReadRecords(chunk));
    std::string out = "[";
    for (std::size_t i = 0; i < records.size(); ++i) {
        out += (i > 0 ? "," : "") + RecordJsonOf(types, records[i], trimmedName, 0);
    }
    return out + "]";
}

[[nodiscard]] Result<std::string> ChunkJson(const Chunk& chunk)
{
    std::string body;
    if (chunk.Tag == CurvesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinCurveTypes{}, chunk, TrimmedName<3>));
    } else if (chunk.Tag == PlaneCurvesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinCurve2dTypes{}, chunk, TrimmedName<2>));
    } else if (chunk.Tag == SurfacesTag) {
        STEREON_TRY(body, RecordsJson(BuiltinSurfaceTypes{}, chunk, std::string_view{}));
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

std::uint32_t Crc32(const std::span<const std::byte> bytes) noexcept
{
    std::uint32_t crc = 0xFFFFFFFFU;
    for (const std::byte b : bytes) {
        crc = Crc32Lookup[(crc ^ std::to_integer<std::uint32_t>(b)) & 0xFFU] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFFU;
}

} // namespace Stereon::Detail
