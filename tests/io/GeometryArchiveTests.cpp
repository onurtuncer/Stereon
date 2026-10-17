// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The geometry archive (ADR-0012, ADR-0028): every built-in type saved and
// read back bit for bit, the numerical policy preserved rather than
// defaulted, damaged and hostile files refused with a ParseError that says
// where, trimmed curves (format version 2) stored with their bases, and the
// committed files of every format version, which every later version must
// still read and reproduce byte for byte.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <vector>

import Stereon.Io;

using Stereon::BSplineCurve;
using Stereon::BSplineCurve2d;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::Cone;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Cylinder;
using Stereon::Ellipse;
using Stereon::Ellipse2d;
using Stereon::ErrorCode;
using Stereon::GeometryArchive;
using Stereon::Hyperbola;
using Stereon::Hyperbola2d;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::LoadedGeometry;
using Stereon::Parabola;
using Stereon::Parabola2d;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::ReadGeometry;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::TrimmedCurve;
using Stereon::TrimmedCurve2d;
using Stereon::Vec2;
using Stereon::Vec3;
using Stereon::WriteGeometry;

namespace {

using Bytes = std::vector<std::byte>;

constexpr Vec3 X{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
constexpr Vec3 Y{2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0};
constexpr Point3 O{1.0, -2.0, 0.5};

/// One of every built-in type, curves, plane curves and surfaces, rational
/// and polynomial, with a non-default policy: the content of the committed
/// version-1 file.
GeometryArchive EveryType()
{
    GeometryArchive a;
    a.Policy.Resolution = 2e-5;
    a.Policy.AngularResolution = 3e-11;
    a.Policy.ToleranceCeiling = 5e-3;
    a.Curves = {
        Curve(Line{O, X}),
        Curve(Circle{O, 2.5, X, Y}),
        Curve(Ellipse{O, 3.0, 1.25, X, Y}),
        Curve(Hyperbola{O, 2.0, 0.75, X, Y}),
        Curve(Parabola{O, 0.5, X, Y}),
        Curve(BSplineCurve{3,
                           {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
                           {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
                            Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}},
                           {1.0, 0.5, 2.0, 0.75, 1.0}})};
    a.Curves2d = {Curve2d(Line2d{Point2{1.0, 2.0}, Vec2{0.6, 0.8}}),
                  Curve2d(Circle2d{Point2{0.5, 0.5}, 0.25, Vec2{1.0, 0.0}, Vec2{0.0, -1.0}}),
                  Curve2d(Ellipse2d{Point2{}, 2.0, 1.0, Vec2{0.0, 1.0}, Vec2{-1.0, 0.0}}),
                  Curve2d(Hyperbola2d{Point2{}, 1.0, 0.5, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}),
                  Curve2d(Parabola2d{Point2{}, 0.25, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}}),
                  Curve2d(BSplineCurve2d{2,
                                         {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                                         {Point2{0.0, 0.0}, Point2{0.5, 1.0}, Point2{1.0, 0.0}},
                                         {}})};
    a.Surfaces = {Surface(Plane{O, X, Y}),
                  Surface(Cylinder{O, 1.5, X, Y, false}),
                  Surface(Cone{O, 1.0, std::numbers::pi / 6.0, X, Y}),
                  Surface(Sphere{O, 2.5, X, Y}),
                  Surface(Torus{O, 3.0, 0.75, X, Y}),
                  Surface(BSplineSurface{1,
                                         2,
                                         {0.0, 0.0, 1.0, 1.0},
                                         {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                                         {Point3{0.0, 0.0, 0.0}, Point3{0.0, 1.0, 0.5},
                                          Point3{0.0, 2.0, 0.0}, Point3{1.0, 0.0, 0.0},
                                          Point3{1.0, 1.0, 1.0}, Point3{1.0, 2.0, 0.0}},
                                         {1.0, 0.5, 1.0, 1.0, 2.0, 1.0}})};
    return a;
}

/// EveryType with trimmed curves added: an arc across a circle's seam, a
/// reversed piece of a B-spline, a trimmed curve of a trimmed curve and a
/// plane arc. The content of the committed version-2 file.
GeometryArchive WithTrimmed()
{
    GeometryArchive a = EveryType();
    const TrimmedCurve arc =
        *Stereon::MakeTrimmedCurve(Curve(Circle{O, 2.5, X, Y}), Stereon::Interval{5.0, 7.0});
    const TrimmedCurve piece =
        Stereon::MakeTrimmedCurve(a.Curves[5], Stereon::Interval{0.25, 0.75})->Reversed();
    const TrimmedCurve inner =
        *Stereon::MakeTrimmedCurve(Curve(piece), Stereon::Interval{0.3, 0.6});
    a.Curves.push_back(Curve(arc));
    a.Curves.push_back(Curve(piece));
    a.Curves.push_back(Curve(inner));
    a.Curves2d.push_back(
        Curve2d(*Stereon::MakeTrimmedCurve(a.Curves2d[1], Stereon::Interval{0.5, 2.0})));
    return a;
}

Bytes Written(const GeometryArchive& a)
{
    auto bytes = WriteGeometry(a);
    REQUIRE(bytes.has_value());
    return *bytes;
}

LoadedGeometry Read(const Bytes& bytes)
{
    auto loaded = ReadGeometry(std::span<const std::byte>{bytes});
    REQUIRE(loaded.has_value());
    return *std::move(loaded);
}

std::uint32_t U32At(const Bytes& b, const std::size_t at)
{
    std::uint32_t x = 0;
    for (std::size_t i = 0; i < 4; ++i) {
        x |= std::uint32_t{std::to_integer<std::uint8_t>(b[at + i])} << (8 * i);
    }
    return x;
}

void PutU32At(Bytes& b, const std::size_t at, const std::uint32_t x)
{
    for (std::size_t i = 0; i < 4; ++i) {
        b[at + i] = static_cast<std::byte>((x >> (8 * i)) & 0xFFU);
    }
}

/// The offset of the first chunk with the given tag; the header is 32
/// bytes and each chunk header 20.
std::size_t ChunkAt(const Bytes& b, const std::string& tag)
{
    std::size_t at = 32;
    while (at + 20 <= b.size()) {
        const std::string here{static_cast<char>(b[at]), static_cast<char>(b[at + 1]),
                               static_cast<char>(b[at + 2]), static_cast<char>(b[at + 3])};
        if (here == tag) {
            return at;
        }
        at += 20 + U32At(b, at + 8); // payload lengths here fit in 32 bits
    }
    FAIL("no chunk " << tag);
    return 0;
}

/// Recomputes a chunk's checksum after its payload was edited.
void Reseal(Bytes& b, const std::size_t chunk)
{
    const std::size_t length = U32At(b, chunk + 8);
    PutU32At(b, chunk + 16,
             Stereon::Detail::Crc32(std::span<const std::byte>{b}.subspan(chunk + 20, length)));
}

/// An extra chunk with the given tag and flags, appended.
Bytes WithChunk(Bytes b, const std::string& tag, const std::uint32_t flags)
{
    const Bytes payload{std::byte{1}, std::byte{2}, std::byte{3}};
    for (const char c : tag) {
        b.push_back(static_cast<std::byte>(c));
    }
    const std::size_t at = b.size();
    b.resize(b.size() + 16);
    PutU32At(b, at, flags);
    PutU32At(b, at + 4, static_cast<std::uint32_t>(payload.size()));
    PutU32At(b, at + 8, 0);
    PutU32At(b, at + 12, Stereon::Detail::Crc32(payload));
    b.insert(b.end(), payload.begin(), payload.end());
    return b;
}

std::filesystem::path GoldenFile(const int version)
{
    return std::filesystem::path(__FILE__).parent_path() / "data" /
           ("geometry-v" + std::to_string(version) + ".stn");
}

Bytes FileBytes(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.good());
    const std::vector<char> chars{std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>()};
    Bytes bytes;
    for (const char c : chars) {
        bytes.push_back(static_cast<std::byte>(static_cast<unsigned char>(c)));
    }
    return bytes;
}

} // namespace

TEST_CASE("every built-in type is saved and read back bit for bit", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        REQUIRE(WriteGeometry(EveryType()).error().GetCode() == ErrorCode::NotSupported);
        return;
    }
    const GeometryArchive original = EveryType();
    const Bytes bytes = Written(original);
    const LoadedGeometry loaded = Read(bytes);
    REQUIRE(loaded.FormatVersion == Stereon::GeometryFormatVersion);
    REQUIRE(loaded.WrittenBy == Stereon::LibraryVersion);
    REQUIRE(loaded.Warnings.empty());
    REQUIRE(loaded.Archive.Curves.size() == original.Curves.size());
    REQUIRE(loaded.Archive.Curves2d.size() == original.Curves2d.size());
    REQUIRE(loaded.Archive.Surfaces.size() == original.Surfaces.size());
    // Writing what was read gives the same bytes: every value, bit for bit.
    REQUIRE(Written(loaded.Archive) == bytes);
    // And the geometry evaluates as before.
    for (std::size_t i = 0; i < original.Curves.size(); ++i) {
        REQUIRE(Eval(loaded.Archive.Curves[i], 0.25) == Eval(original.Curves[i], 0.25));
    }
    for (std::size_t i = 0; i < original.Surfaces.size(); ++i) {
        const Point2 uv{0.25, 0.25};
        REQUIRE(Eval(loaded.Archive.Surfaces[i], uv) == Eval(original.Surfaces[i], uv));
    }
}

TEST_CASE("the numerical policy is read back as stored, never defaulted", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const GeometryArchive original = EveryType();
    const LoadedGeometry loaded = Read(Written(original));
    REQUIRE(loaded.Archive.Policy == original.Policy);
    REQUIRE(loaded.Archive.Policy.Resolution == 2e-5); // not Context::DefaultResolution

    // A policy that is not finite and positive is refused, on write and on read.
    GeometryArchive bad = original;
    bad.Policy.ToleranceCeiling = -1.0;
    REQUIRE(WriteGeometry(bad).error().GetCode() == ErrorCode::InvalidInput);
    Bytes bytes = Written(original);
    const std::size_t policy = ChunkAt(bytes, "POLI");
    bytes[policy + 20 + 7] = std::byte{0xFF}; // angular resolution's sign and exponent: NaN
    bytes[policy + 20 + 6] = std::byte{0xFF};
    Reseal(bytes, policy);
    REQUIRE(ReadGeometry(bytes).error().GetCode() == ErrorCode::ParseError);

    // Without a policy chunk there is no policy to read, and no default to
    // substitute (ADR-0028).
    Bytes noPolicy = Written(original);
    noPolicy[ChunkAt(noPolicy, "POLI")] = std::byte{'p'}; // now an unknown chunk...
    const std::size_t renamed = ChunkAt(noPolicy, "pOLI");
    PutU32At(noPolicy, renamed + 4, 0); // ...and optional, so only its absence matters
    const auto missing = ReadGeometry(noPolicy);
    REQUIRE(missing.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(missing.error().GetMessage()).find("no policy") != std::string::npos);
}

TEST_CASE("damaged archives are refused with a parse error", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const Bytes good = Written(EveryType());
    const auto codeOf = [](const Bytes& b) { return ReadGeometry(b).error().GetCode(); };

    Bytes badMagic = good;
    badMagic[0] = std::byte{'X'};
    REQUIRE(codeOf(badMagic) == ErrorCode::ParseError);

    Bytes newer = good;
    PutU32At(newer, 4, Stereon::GeometryFormatVersion + 1);
    REQUIRE(codeOf(newer) == ErrorCode::ParseError);

    for (const std::size_t cut :
         {std::size_t{0}, std::size_t{10}, std::size_t{40}, good.size() / 2, good.size() - 1}) {
        const Bytes truncated(good.begin(), good.begin() + static_cast<std::ptrdiff_t>(cut));
        REQUIRE(codeOf(truncated) == ErrorCode::ParseError);
    }

    Bytes flipped = good;
    flipped[ChunkAt(flipped, "SRFS") + 40] ^= std::byte{0x10};
    const auto damaged = ReadGeometry(flipped);
    REQUIRE(damaged.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(damaged.error().GetMessage()).find("checksum") != std::string::npos);
}

TEST_CASE("invalid geometry in a file is caught by the checked factories", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    GeometryArchive one;
    one.Curves = {Curve(Circle{O, 2.5, X, Y})};
    Bytes bytes = Written(one);
    // The circle record: count (8), name length (2), "Circle" (6), value
    // length (8), then Center (24) and Radius (8). Make the radius -2.5.
    const std::size_t curves = ChunkAt(bytes, "CRVS");
    bytes[curves + 20 + 8 + 2 + 6 + 8 + 24 + 7] ^= std::byte{0x80};
    Reseal(bytes, curves);
    const auto read = ReadGeometry(bytes);
    REQUIRE(read.error().GetCode() == ErrorCode::ParseError);
    const std::string message{read.error().GetMessage()};
    REQUIRE(message.find("record 0 of chunk \"CRVS\" (Circle) is not valid") != std::string::npos);

    // An unknown type name in a geometry chunk: the geometry cannot be read,
    // so the file is refused.
    Bytes renamed = Written(one);
    renamed[ChunkAt(renamed, "CRVS") + 20 + 8 + 2] = std::byte{'Q'};
    Reseal(renamed, ChunkAt(renamed, "CRVS"));
    REQUIRE(
        std::string(ReadGeometry(renamed).error().GetMessage()).find("unknown type \"Qircle\"") !=
        std::string::npos);
}

TEST_CASE("unknown chunks: required ones are refused, optional ones skipped", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const Bytes good = Written(EveryType());
    const auto required = ReadGeometry(WithChunk(good, "HIST", 1));
    REQUIRE(required.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(required.error().GetMessage()).find("unknown required chunk \"HIST\"") !=
            std::string::npos);

    const LoadedGeometry optional = Read(WithChunk(good, "NOTE", 0));
    REQUIRE(optional.Warnings.size() == 1);
    REQUIRE(optional.Warnings[0].find("NOTE") != std::string::npos);
    REQUIRE(Written(optional.Archive) == good); // the geometry is untouched
}

namespace {

/// A curve type the kernel does not know: the archive has no record for it.
struct UserCurve {
    Line Base;
};

Point3 Eval(const UserCurve& c, const double t)
{
    return Eval(c.Base, t);
}
Vec3 Deriv(const UserCurve& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const UserCurve& c)
{
    return Domain(c.Base);
}
Stereon::Box3 BoundingBox(const UserCurve& c)
{
    return BoundingBox(c.Base);
}
double Project(const UserCurve& c, const Point3& p)
{
    return Project(c.Base, p);
}
UserCurve Reversed(const UserCurve& c)
{
    return {Reversed(c.Base)};
}
UserCurve Transformed(const UserCurve& c, const Stereon::Transform& x)
{
    return {Transformed(c.Base, x)};
}

} // namespace

TEST_CASE("a user's curve type cannot be stored", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    GeometryArchive a;
    a.Curves = {Curve(Line{O, X}), Curve(UserCurve{Line{O, Y}})};
    const auto written = WriteGeometry(a);
    REQUIRE(written.error().GetCode() == ErrorCode::NotSupported);
    REQUIRE(std::string(written.error().GetMessage()).find("curve 1") != std::string::npos);
}

TEST_CASE("files round-trip, and the JSON view shows them", "[io][archive]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const auto path = std::filesystem::temp_directory_path() / "stereon-archive-test.stn";
    REQUIRE(Stereon::SaveGeometry(path, EveryType()).has_value());
    const auto loaded = Stereon::LoadGeometry(path);
    REQUIRE(loaded.has_value());
    REQUIRE(Written(loaded->Archive) == Written(EveryType()));
    std::filesystem::remove(path);
    REQUIRE(Stereon::LoadGeometry(path).error().GetCode() == ErrorCode::InvalidInput);

    const auto json = Stereon::GeometryToJson(Written(EveryType()));
    REQUIRE(json.has_value());
    REQUIRE(json->starts_with(R"({"format":2,)"));
    REQUIRE(json->find(R"("resolution":2e-05)") != std::string::npos);
    REQUIRE(
        json->find(R"({"type":"Circle","value":{"Center":{"X":1,"Y":-2,"Z":0.5},"Radius":2.5,)") !=
        std::string::npos);
    REQUIRE(json->find(R"("tag":"SRFS","required":true)") != std::string::npos);
}

TEST_CASE("the committed files of every format version still read and rewrite byte for byte",
          "[io][archive][compatibility]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // ADR-0012: files written by each release stay readable by every later
    // one. Each version's file is frozen once committed; only the current
    // version's is regenerated, when that version is new, never to make
    // this test pass, by running the hidden "[.regenerate]" case.
    const struct {
        int Version;
        GeometryArchive Content;
    } files[] = {{1, EveryType()}, {2, WithTrimmed()}};
    for (const auto& [version, content] : files) {
        const Bytes golden = FileBytes(GoldenFile(version));
        const LoadedGeometry loaded = Read(golden);
        REQUIRE(loaded.FormatVersion == static_cast<std::uint32_t>(version));
        REQUIRE(loaded.Archive.Policy == content.Policy);
        REQUIRE(loaded.Archive.Curves.size() == content.Curves.size());
        REQUIRE(loaded.Archive.Curves2d.size() == content.Curves2d.size());
        // The format and kernel versions in the header are the writer's, so
        // compare the rest: everything after the header's version fields.
        // Version 1 is version 2 without trimmed curves, so rewriting it as
        // version 2 changes nothing else.
        const Bytes rewritten = Written(loaded.Archive);
        REQUIRE(rewritten.size() == golden.size());
        REQUIRE(std::equal(rewritten.begin() + 16, rewritten.end(), golden.begin() + 16));
        REQUIRE(Written(content).size() == golden.size());
    }
}

TEST_CASE("regenerate the committed file of the current format version", "[.regenerate]")
{
    static_assert(Stereon::GeometryFormatVersion == 2,
                  "a new format version gets a new committed file; earlier ones stay frozen");
    REQUIRE(Stereon::SaveGeometry(GoldenFile(2), WithTrimmed()).has_value());
}

TEST_CASE("trimmed curves are saved with their bases and read back bit for bit",
          "[io][archive][trimmed]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const GeometryArchive original = WithTrimmed();
    const Bytes bytes = Written(original);
    const LoadedGeometry loaded = Read(bytes);
    REQUIRE(Written(loaded.Archive) == bytes);
    for (std::size_t i = 6; i < original.Curves.size(); ++i) {
        const auto* before = original.Curves[i].Target<TrimmedCurve>();
        const auto* after = loaded.Archive.Curves[i].Target<TrimmedCurve>();
        REQUIRE(before != nullptr);
        REQUIRE(after != nullptr);
        REQUIRE(after->GetRange() == before->GetRange());
        REQUIRE(after->IsReversed() == before->IsReversed());
        const Stereon::Interval range = before->GetRange();
        for (const double f : {0.0, 0.3, 1.0}) {
            const double t = range.GetLower() + f * (range.GetUpper() - range.GetLower());
            REQUIRE(Eval(*after, t) == Eval(*before, t));
        }
    }
    // The nested one: a trimmed curve whose basis is a trimmed curve.
    const auto* nested = loaded.Archive.Curves[8].Target<TrimmedCurve>();
    REQUIRE(nested->GetBasis().Target<TrimmedCurve>() != nullptr);
    REQUIRE(nested->GetBasis().Target<TrimmedCurve>()->IsReversed());
    const auto* plane = loaded.Archive.Curves2d.back().Target<TrimmedCurve2d>();
    REQUIRE(plane != nullptr);
    REQUIRE(Eval(*plane, 1.0) == Eval(original.Curves2d[1], 1.0));

    // The JSON view shows the basis and the trim.
    const auto json = Stereon::GeometryToJson(bytes);
    REQUIRE(json.has_value());
    REQUIRE(json->starts_with(R"({"format":2,)"));
    REQUIRE(json->find(R"({"type":"TrimmedCurve","value":{"basis":{"type":"Circle","value":)") !=
            std::string::npos);
    REQUIRE(json->find(R"("trim":{"Lower":5,"Upper":7,"Reversed":false})") != std::string::npos);
    REQUIRE(json->find(R"({"type":"TrimmedCurve2d","value":{"basis":{"type":"Circle2d")") !=
            std::string::npos);
}

namespace {

/// A line trimmed `depth` times, each trim inside the one before.
Curve NestedTrims(const int depth)
{
    Curve curve(Line{O, X});
    for (int k = 0; k < depth; ++k) {
        const double half = 10.0 - k;
        curve = Curve(*Stereon::MakeTrimmedCurve(curve, Stereon::Interval{-half, half}));
    }
    return curve;
}

/// The archive with its one curve record wrapped in one more trimmed-curve
/// record with the same trim, written by hand: the writer refuses to nest
/// so deep, a hostile file need not.
Bytes WrappedOnceMore(Bytes b)
{
    const std::size_t chunk = ChunkAt(b, "CRVS");
    const std::size_t length = U32At(b, chunk + 8);
    const Bytes record(b.begin() + static_cast<std::ptrdiff_t>(chunk + 20 + 8),
                       b.begin() + static_cast<std::ptrdiff_t>(chunk + 20 + length));
    // The record ends with its trim: two doubles and a bool.
    Bytes value = record;
    value.insert(value.end(), record.end() - 17, record.end());
    const std::string name = "TrimmedCurve";
    Bytes payload(8, std::byte{0});
    payload[0] = std::byte{1};
    payload.push_back(static_cast<std::byte>(name.size()));
    payload.push_back(std::byte{0});
    for (const char c : name) {
        payload.push_back(static_cast<std::byte>(c));
    }
    for (std::size_t i = 0; i < 8; ++i) {
        payload.push_back(static_cast<std::byte>((value.size() >> (8 * i)) & 0xFFU));
    }
    payload.insert(payload.end(), value.begin(), value.end());
    Bytes out(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(chunk + 20));
    PutU32At(out, chunk + 8, static_cast<std::uint32_t>(payload.size()));
    out.insert(out.end(), payload.begin(), payload.end());
    out.insert(out.end(), b.begin() + static_cast<std::ptrdiff_t>(chunk + 20 + length), b.end());
    Reseal(out, chunk);
    return out;
}

} // namespace

TEST_CASE("trimmed-curve records are validated and their nesting bounded", "[io][archive][trimmed]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // The deepest nesting the format allows is written and read back...
    GeometryArchive deep;
    deep.Curves = {NestedTrims(Stereon::MaxTrimmedNesting)};
    const Bytes deepest = Written(deep);
    REQUIRE(Written(Read(deepest).Archive) == deepest);
    // ...one more is refused on write, and on read from a file made by hand.
    GeometryArchive deeper;
    deeper.Curves = {NestedTrims(Stereon::MaxTrimmedNesting + 1)};
    REQUIRE(WriteGeometry(deeper).error().GetCode() == ErrorCode::NotSupported);
    const auto tooDeep = ReadGeometry(WrappedOnceMore(deepest));
    REQUIRE(tooDeep.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(tooDeep.error().GetMessage()).find("more than 8 deep") !=
            std::string::npos);
    // The hand-made wrapping is a valid record at a legal depth.
    GeometryArchive shallow;
    shallow.Curves = {NestedTrims(1)};
    REQUIRE(Read(WrappedOnceMore(Written(shallow))).Archive.Curves.size() == 1);

    // A trim edited in the file: the trimmed curve record is the circle arc
    // [5, 7]; its trim is the record's last 17 bytes.
    GeometryArchive arc;
    arc.Curves = {Curve(
        *Stereon::MakeTrimmedCurve(Curve(Circle{O, 2.5, X, Y}), Stereon::Interval{5.0, 7.0}))};
    const Bytes good = Written(arc);
    const std::size_t curves = ChunkAt(good, "CRVS");
    const std::size_t end = curves + 20 + U32At(good, curves + 8);
    const auto withUpper = [&](const double upper) {
        Bytes b = good;
        const auto bits = std::bit_cast<std::uint64_t>(upper);
        for (std::size_t i = 0; i < 8; ++i) {
            b[end - 9 + i] = static_cast<std::byte>((bits >> (8 * i)) & 0xFFU);
        }
        Reseal(b, curves);
        return std::string(ReadGeometry(b).error().GetMessage());
    };
    REQUIRE(withUpper(std::numeric_limits<double>::quiet_NaN()).find("not finite and non-empty") !=
            std::string::npos);
    REQUIRE(withUpper(4.0).find("not finite and non-empty") != std::string::npos);
    // A range longer than a period is caught by MakeTrimmedCurve.
    REQUIRE(withUpper(12.0).find("(TrimmedCurve) is not valid geometry") != std::string::npos);
    // A sense byte that is not 0 or 1.
    Bytes badSense = good;
    badSense[end - 1] = std::byte{2};
    Reseal(badSense, curves);
    REQUIRE(ReadGeometry(badSense).error().GetCode() == ErrorCode::ParseError);
    // A trimmed record in the plane-curve chunk is of an unknown type there.
    GeometryArchive plane;
    plane.Curves2d = {Curve2d(*Stereon::MakeTrimmedCurve(
        Curve2d(Line2d{Point2{0.0, 0.0}, Vec2{1.0, 0.0}}), Stereon::Interval{0.0, 1.0}))};
    Bytes renamed = Written(plane);
    const std::size_t planeCurves = ChunkAt(renamed, "PCRV");
    // "TrimmedCurve2d" -> "TrimmedCurve2x": the name, after count and length.
    renamed[planeCurves + 20 + 8 + 2 + 13] = std::byte{'x'};
    Reseal(renamed, planeCurves);
    REQUIRE(std::string(ReadGeometry(renamed).error().GetMessage()).find("unknown type") !=
            std::string::npos);
}
