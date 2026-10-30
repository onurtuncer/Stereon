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
#include <atomic>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <mutex>
#include <numbers>
#include <ostream>
#include <span>
#include <sstream>
#include <streambuf>
#include <string>
#include <thread>
#include <vector>

import Stereon.Io;

using Stereon::Box2;
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
using Stereon::ExtrusionSurface;
using Stereon::GeometryArchive;
using Stereon::Hyperbola;
using Stereon::Hyperbola2d;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Line2d;
using Stereon::LoadedGeometry;
using Stereon::Parabola;
using Stereon::Parabola2d;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::ReadGeometry;
using Stereon::RevolutionSurface;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Torus;
using Stereon::TrimmedCurve;
using Stereon::TrimmedCurve2d;
using Stereon::TrimmedSurface;
using Stereon::Vec2;
using Stereon::Vec3;
using Stereon::WriteGeometry;

namespace {

/// The message of a failed result, copied: the result is named by the
/// parameter, so the message is read from a named owner (ADR-0046).
template <typename R>
std::string MessageOf(const R& result)
{
    return std::string(result.error().GetMessage());
}

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

/// WithTrimmed with the surface records of format version 3 added: a
/// cylinder patch across its seam, a reversed corner of the B-spline
/// surface, an extrusion of the circle arc (a trimmed profile), a
/// revolution of the B-spline curve about a tilted axis, and a band of that
/// revolution. The content of the committed version-3 file.
GeometryArchive WithSurfaces()
{
    GeometryArchive a = WithTrimmed();
    const TrimmedSurface patch =
        *Stereon::MakeTrimmedSurface(a.Surfaces[1], Box2{Interval{5.0, 7.0}, Interval{-1.0, 2.0}});
    const TrimmedSurface corner =
        Stereon::MakeTrimmedSurface(a.Surfaces[5], Box2{Interval{0.25, 0.75}, Interval{0.0, 0.5}})
            ->Reversed();
    const ExtrusionSurface extrusion = *Stereon::MakeExtrusion(a.Curves[6], Vec3{0.5, -1.0, 2.0});
    const RevolutionSurface revolution =
        *Stereon::MakeRevolution(a.Curves[5], O, Vec3{2.0, -1.0, 3.0});
    const TrimmedSurface band = *Stereon::MakeTrimmedSurface(
        Surface(revolution), Box2{Interval{1.0, 4.0}, Interval{0.2, 0.8}});
    a.Surfaces.push_back(Surface(patch));
    a.Surfaces.push_back(Surface(corner));
    a.Surfaces.push_back(Surface(extrusion));
    a.Surfaces.push_back(Surface(revolution));
    a.Surfaces.push_back(Surface(band));
    return a;
}

/// A parameter inside the surface's domain, away from its edges.
Point2 SomewhereOn(const Surface& s)
{
    const Box2 d = Domain(s);
    const auto inside = [](const Interval& i) {
        return std::isfinite(i.GetLower()) && std::isfinite(i.GetUpper())
                   ? i.GetLower() + 0.3 * (i.GetUpper() - i.GetLower())
                   : 0.3;
    };
    return {inside(d.X), inside(d.Y)};
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
    REQUIRE(MessageOf(ReadGeometry(renamed)).find("unknown type \"Qircle\"") != std::string::npos);
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
    REQUIRE(json->starts_with(R"({"format":3,)"));
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
    } files[] = {{1, EveryType()}, {2, WithTrimmed()}, {3, WithSurfaces()}};
    for (const auto& [version, content] : files) {
        const Bytes golden = FileBytes(GoldenFile(version));
        const LoadedGeometry loaded = Read(golden);
        REQUIRE(loaded.FormatVersion == static_cast<std::uint32_t>(version));
        REQUIRE(loaded.Archive.Policy == content.Policy);
        REQUIRE(loaded.Archive.Curves.size() == content.Curves.size());
        REQUIRE(loaded.Archive.Curves2d.size() == content.Curves2d.size());
        REQUIRE(loaded.Archive.Surfaces.size() == content.Surfaces.size());
        // The format and kernel versions in the header are the writer's, so
        // compare the rest: everything after the header's version fields.
        // Each version is the next without its new records, so rewriting an
        // earlier one as the current version changes nothing else.
        const Bytes rewritten = Written(loaded.Archive);
        REQUIRE(rewritten.size() == golden.size());
        REQUIRE(std::equal(rewritten.begin() + 16, rewritten.end(), golden.begin() + 16));
        REQUIRE(Written(content).size() == golden.size());
    }
}

TEST_CASE("regenerate the committed file of the current format version", "[.regenerate]")
{
    static_assert(Stereon::GeometryFormatVersion == 3,
                  "a new format version gets a new committed file; earlier ones stay frozen");
    REQUIRE(Stereon::SaveGeometry(GoldenFile(3), WithSurfaces()).has_value());
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
    REQUIRE(json->starts_with(R"({"format":3,)"));
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

/// The archive with the one record of chunk `tag` wrapped in one more
/// composite record `name` with the same tail (the record's last `tail`
/// bytes), written by hand: the writer refuses to nest so deep, a hostile
/// file need not.
Bytes WrappedOnceMore(Bytes b, const std::string& tag, const std::string& name,
                      const std::size_t tail)
{
    const std::size_t chunk = ChunkAt(b, tag);
    const std::size_t length = U32At(b, chunk + 8);
    const Bytes record(b.begin() + static_cast<std::ptrdiff_t>(chunk + 20 + 8),
                       b.begin() + static_cast<std::ptrdiff_t>(chunk + 20 + length));
    Bytes value = record;
    value.insert(value.end(), record.end() - static_cast<std::ptrdiff_t>(tail), record.end());
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
    const auto tooDeep = ReadGeometry(WrappedOnceMore(deepest, "CRVS", "TrimmedCurve", 17));
    REQUIRE(tooDeep.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(tooDeep.error().GetMessage()).find("more than 8 deep") !=
            std::string::npos);
    // The hand-made wrapping is a valid record at a legal depth.
    GeometryArchive shallow;
    shallow.Curves = {NestedTrims(1)};
    REQUIRE(
        Read(WrappedOnceMore(Written(shallow), "CRVS", "TrimmedCurve", 17)).Archive.Curves.size() ==
        1);

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
        return MessageOf(ReadGeometry(b));
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
    REQUIRE(MessageOf(ReadGeometry(renamed)).find("unknown type") != std::string::npos);
}

TEST_CASE("swept NURBS surfaces survive the geometry archive", "[io][swept]")
{
    const BSplineCurve profile{2, {0, 0, 0, 1, 1, 1}, {{2, 0, 0}, {3, 1, 1}, {2, 0, 2}}, {1, 2, 1}};
    const auto extruded =
        Stereon::MakeExtrusionSurface(profile, Vec3{0, 1, 2}, Stereon::Interval{-1, 2});
    const auto revolved = Stereon::MakeRevolutionSurface(
        profile, Point3{}, Vec3{0, 0, 1}, Stereon::Interval{0, Stereon::CirclePeriod});
    REQUIRE(extruded.has_value());
    REQUIRE(revolved.has_value());
    GeometryArchive archive;
    archive.Surfaces = {Surface{*extruded}, Surface{*revolved}};
    const auto bytes = WriteGeometry(archive);
    REQUIRE(bytes.has_value());
    const auto loaded = ReadGeometry(*bytes);
    REQUIRE(loaded.has_value());
    const auto writtenAgain = WriteGeometry(loaded->Archive);
    REQUIRE(writtenAgain.has_value());
    REQUIRE(*writtenAgain == *bytes);
}

TEST_CASE("trimmed surfaces and sweeps are saved with their bases and profiles and read back "
          "bit for bit",
          "[io][archive][surfaces]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    const GeometryArchive original = WithSurfaces();
    const Bytes bytes = Written(original);
    const LoadedGeometry loaded = Read(bytes);
    REQUIRE(loaded.FormatVersion == 3);
    REQUIRE(Written(loaded.Archive) == bytes);
    const std::vector<Surface>& after = loaded.Archive.Surfaces;
    REQUIRE(after.size() == original.Surfaces.size());
    for (std::size_t i = 6; i < after.size(); ++i) {
        REQUIRE(Domain(after[i]) == Domain(original.Surfaces[i]));
        const Point2 uv = SomewhereOn(original.Surfaces[i]);
        REQUIRE(Eval(after[i], uv) == Eval(original.Surfaces[i], uv));
        REQUIRE(Deriv(after[i], uv, 1, 1) == Deriv(original.Surfaces[i], uv, 1, 1));
    }
    const auto* patch = after[6].Target<TrimmedSurface>();
    REQUIRE(patch != nullptr);
    REQUIRE(patch->GetRange() == Box2{Interval{5.0, 7.0}, Interval{-1.0, 2.0}});
    REQUIRE_FALSE(patch->IsReversed());
    REQUIRE(patch->GetBasis().Target<Cylinder>() != nullptr);
    const auto* corner = after[7].Target<TrimmedSurface>();
    REQUIRE(corner != nullptr);
    REQUIRE(corner->IsReversed());
    const auto* extrusion = after[8].Target<ExtrusionSurface>();
    REQUIRE(extrusion != nullptr);
    REQUIRE(extrusion->GetDirection() == Vec3{0.5, -1.0, 2.0});
    REQUIRE(extrusion->GetProfile().Target<TrimmedCurve>() != nullptr);
    const auto* revolution = after[9].Target<RevolutionSurface>();
    REQUIRE(revolution != nullptr);
    REQUIRE(revolution->GetOrigin() == O);
    REQUIRE(revolution->GetAxis() == original.Surfaces[9].Target<RevolutionSurface>()->GetAxis());
    REQUIRE(revolution->GetProfile().Target<BSplineCurve>() != nullptr);
    const auto* band = after[10].Target<TrimmedSurface>();
    REQUIRE(band != nullptr);
    REQUIRE(band->GetBasis().Target<RevolutionSurface>() != nullptr);

    // The JSON view shows the parts.
    const auto json = Stereon::GeometryToJson(bytes);
    REQUIRE(json.has_value());
    REQUIRE(
        json->find(R"({"type":"TrimmedSurface","value":{"basis":{"type":"Cylinder","value":)") !=
        std::string::npos);
    REQUIRE(
        json->find(R"("trim":{"LowerU":5,"UpperU":7,"LowerV":-1,"UpperV":2,"Reversed":false})") !=
        std::string::npos);
    REQUIRE(json->find(R"({"type":"ExtrusionSurface","value":{"profile":{"type":"TrimmedCurve")") !=
            std::string::npos);
    REQUIRE(json->find(R"("sweep":{"Direction":{"X":0.5,"Y":-1,"Z":2}})") != std::string::npos);
    REQUIRE(
        json->find(R"({"type":"RevolutionSurface","value":{"profile":{"type":"BSplineCurve")") !=
        std::string::npos);
    REQUIRE(json->find(R"("sweep":{"Origin":{"X":1,"Y":-2,"Z":0.5},"Axis":{"X":)") !=
            std::string::npos);
    REQUIRE(
        json->find(R"({"type":"TrimmedSurface","value":{"basis":{"type":"RevolutionSurface")") !=
        std::string::npos);
}

namespace {

/// A plane trimmed `depth` times, each box inside the one before.
Surface NestedSurfaceTrims(const int depth)
{
    Surface surface(Plane{O, X, Y});
    for (int k = 0; k < depth; ++k) {
        const double half = 10.0 - k;
        surface = Surface(*Stereon::MakeTrimmedSurface(
            surface, Box2{Interval{-half, half}, Interval{-half, half}}));
    }
    return surface;
}

/// The error message of reading `b` with `bytes` written over its surface
/// chunk's last `tail` bytes, resealed.
std::string ReadWithTail(Bytes b, const Bytes& bytes)
{
    const std::size_t surfaces = ChunkAt(b, "SRFS");
    const std::size_t end = surfaces + 20 + U32At(b, surfaces + 8);
    std::ranges::copy(bytes, b.begin() + static_cast<std::ptrdiff_t>(end - bytes.size()));
    Reseal(b, surfaces);
    return MessageOf(ReadGeometry(b));
}

Bytes DoubleBytes(const double x)
{
    Bytes out(8);
    const auto bits = std::bit_cast<std::uint64_t>(x);
    for (std::size_t i = 0; i < 8; ++i) {
        out[i] = static_cast<std::byte>((bits >> (8 * i)) & 0xFFU);
    }
    return out;
}

} // namespace

TEST_CASE("trimmed-surface and sweep records are validated and their nesting bounded",
          "[io][archive][surfaces]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // The deepest nesting the format allows is written and read back...
    GeometryArchive deep;
    deep.Surfaces = {NestedSurfaceTrims(Stereon::MaxTrimmedNesting)};
    const Bytes deepest = Written(deep);
    REQUIRE(Written(Read(deepest).Archive) == deepest);
    // ...one more is refused on write, and on read from a file made by hand.
    GeometryArchive deeper;
    deeper.Surfaces = {NestedSurfaceTrims(Stereon::MaxTrimmedNesting + 1)};
    REQUIRE(WriteGeometry(deeper).error().GetCode() == ErrorCode::NotSupported);
    const auto tooDeep = ReadGeometry(WrappedOnceMore(deepest, "SRFS", "TrimmedSurface", 33));
    REQUIRE(tooDeep.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(tooDeep.error().GetMessage()).find("more than 8 deep") !=
            std::string::npos);
    // Trims and sweeps count together: an extrusion of a curve trimmed as
    // deep as the format allows is one level too many.
    const Vec3 up{0.0, 0.0, 1.0};
    GeometryArchive sweptDeep;
    sweptDeep.Surfaces = {
        Surface(*Stereon::MakeExtrusion(NestedTrims(Stereon::MaxTrimmedNesting), up))};
    REQUIRE(WriteGeometry(sweptDeep).error().GetCode() == ErrorCode::NotSupported);
    GeometryArchive sweptEnough;
    sweptEnough.Surfaces = {
        Surface(*Stereon::MakeExtrusion(NestedTrims(Stereon::MaxTrimmedNesting - 1), up))};
    const Bytes swept = Written(sweptEnough);
    REQUIRE(Written(Read(swept).Archive) == swept);

    // A trim edited in the file: the trimmed-surface record is the cylinder
    // patch u in [5, 7], v in [-1, 2]; its trim is the record's last 33
    // bytes, four doubles and the sense byte.
    GeometryArchive patch;
    patch.Surfaces = {Surface(*Stereon::MakeTrimmedSurface(
        Surface(Cylinder{O, 1.5, X, Y, false}), Box2{Interval{5.0, 7.0}, Interval{-1.0, 2.0}}))};
    const Bytes good = Written(patch);
    const auto withUpperU = [&](const double upper) {
        Bytes tail = DoubleBytes(upper);
        const Bytes rest = DoubleBytes(-1.0);
        tail.insert(tail.end(), rest.begin(), rest.end());
        const Bytes upperV = DoubleBytes(2.0);
        tail.insert(tail.end(), upperV.begin(), upperV.end());
        tail.push_back(std::byte{0});
        return ReadWithTail(good, tail);
    };
    REQUIRE(withUpperU(std::numeric_limits<double>::quiet_NaN()).find("not finite and non-empty") !=
            std::string::npos);
    REQUIRE(withUpperU(4.0).find("not finite and non-empty") != std::string::npos);
    // A box wider than a period is caught by MakeTrimmedSurface.
    REQUIRE(withUpperU(12.0).find("(TrimmedSurface) is not valid geometry") != std::string::npos);
    // A sense byte that is not 0 or 1.
    REQUIRE(ReadWithTail(good, Bytes{std::byte{2}}).find("TrimmedSurface") != std::string::npos);

    // A sweep edited in the file: the extrusion's direction is its record's
    // last 24 bytes, the revolution's origin and axis its last 48.
    GeometryArchive sweeps;
    sweeps.Surfaces = {Surface(*Stereon::MakeExtrusion(Curve(Line{O, X}), up)),
                       Surface(*Stereon::MakeRevolution(Curve(Line{O, X}), O, up))};
    const Bytes goodSweeps = Written(sweeps);
    const Bytes zeros(24, std::byte{0});
    const std::string zeroAxis = ReadWithTail(goodSweeps, zeros);
    CAPTURE(zeroAxis);
    REQUIRE(zeroAxis.find("(RevolutionSurface) is not valid geometry") != std::string::npos);
    GeometryArchive extrusion;
    extrusion.Surfaces = {sweeps.Surfaces[0]};
    const std::string zeroDirection = ReadWithTail(Written(extrusion), zeros);
    CAPTURE(zeroDirection);
    REQUIRE(zeroDirection.find("(ExtrusionSurface) is not valid geometry") != std::string::npos);
    // Origin (NaN, 0, 0), axis (0, 0, 1).
    Bytes nanOrigin = DoubleBytes(std::numeric_limits<double>::quiet_NaN());
    const Bytes sixteen(16, std::byte{0});
    const Bytes one = DoubleBytes(1.0);
    nanOrigin.insert(nanOrigin.end(), sixteen.begin(), sixteen.end());
    nanOrigin.insert(nanOrigin.end(), sixteen.begin(), sixteen.end());
    nanOrigin.insert(nanOrigin.end(), one.begin(), one.end());
    const std::string badOrigin = ReadWithTail(goodSweeps, nanOrigin);
    CAPTURE(badOrigin);
    REQUIRE(badOrigin.find("(RevolutionSurface) is not valid geometry") != std::string::npos);
    // An axis that is not of unit length is not a surface's state: origin O,
    // axis (0, 0, 2).
    Bytes longAxis;
    for (const double x : {O.X, O.Y, O.Z, 0.0, 0.0, 2.0}) {
        const Bytes word = DoubleBytes(x);
        longAxis.insert(longAxis.end(), word.begin(), word.end());
    }
    const std::string notUnit = ReadWithTail(goodSweeps, longAxis);
    CAPTURE(notUnit);
    REQUIRE(notUnit.find("the axis is not of unit length") != std::string::npos);

    // A composite record under an unknown name is of an unknown type.
    Bytes renamed = good;
    const std::size_t surfaces = ChunkAt(renamed, "SRFS");
    // "TrimmedSurface" -> "TrimmedSurfacx": the name, after count and length.
    renamed[surfaces + 20 + 8 + 2 + 13] = std::byte{'x'};
    Reseal(renamed, surfaces);
    REQUIRE(MessageOf(ReadGeometry(renamed)).find("unknown type") != std::string::npos);
}

TEST_CASE("a revolution's unit axis survives the round trip whatever its direction",
          "[io][archive][surfaces]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // The writer stores the unit axis the surface keeps. Normalising it
    // again on read could move it by an ulp, so the reader keeps it as
    // stored: the round trip is bit for bit in every direction.
    const Curve profile(Line{O, X});
    int count = 0;
    for (int i = -4; i <= 4; ++i) {
        for (int j = -4; j <= 4; ++j) {
            for (int k = -4; k <= 4; ++k) {
                const Vec3 axis{0.37 * i + 0.011 * j, 0.59 * j - 0.023 * k, 0.83 * k + 0.007 * i};
                if (!(Dot(axis, axis) > 0.0)) {
                    continue;
                }
                GeometryArchive a;
                a.Surfaces = {Surface(*Stereon::MakeRevolution(profile, O, axis))};
                const Bytes bytes = Written(a);
                REQUIRE(Written(Read(bytes).Archive) == bytes);
                ++count;
            }
        }
    }
    REQUIRE(count == 9 * 9 * 9 - 1);
}

// ---------------------------------------------------------------------------
// Audit 2026-10-20: archive integrity (F1), policy (F8), saving (F7) and
// sweeps of trimmed profiles read back (F2)

namespace {

/// One chunk of an archive: its tag and its bytes [Start, End), header and
/// payload.
struct ChunkSpan {
    std::string Tag;
    std::size_t Start;
    std::size_t End;
};

/// The chunks of an archive in order: the header is 32 bytes and each chunk
/// header 20.
std::vector<ChunkSpan> ChunksOf(const Bytes& b)
{
    std::vector<ChunkSpan> chunks;
    for (std::size_t at = 32; at + 20 <= b.size();) {
        const std::string tag{static_cast<char>(b[at]), static_cast<char>(b[at + 1]),
                              static_cast<char>(b[at + 2]), static_cast<char>(b[at + 3])};
        const std::size_t end = at + 20 + U32At(b, at + 8);
        chunks.push_back({tag, at, end});
        at = end;
    }
    return chunks;
}

Bytes Slice(const Bytes& b, const std::size_t from, const std::size_t to)
{
    return {b.begin() + static_cast<std::ptrdiff_t>(from),
            b.begin() + static_cast<std::ptrdiff_t>(to)};
}

/// The bytes with a required chunk appended, `payload` its content.
Bytes Appended(Bytes b, const std::string& tag, const Bytes& payload)
{
    for (const char c : tag) {
        b.push_back(static_cast<std::byte>(c));
    }
    const std::size_t at = b.size();
    b.resize(b.size() + 16);
    PutU32At(b, at, 1);
    PutU32At(b, at + 4, static_cast<std::uint32_t>(payload.size()));
    PutU32At(b, at + 8, 0);
    PutU32At(b, at + 12, Stereon::Detail::Crc32(payload));
    b.insert(b.end(), payload.begin(), payload.end());
    return b;
}

/// The error of reading `b`, which must fail with `code`.
std::string ReadError(const Bytes& b, const ErrorCode code = ErrorCode::ParseError)
{
    const auto read = ReadGeometry(b);
    REQUIRE_FALSE(read.has_value());
    const std::string message(read.error().GetMessage());
    INFO(message);
    REQUIRE(read.error().GetCode() == code);
    return message;
}

} // namespace

TEST_CASE("an archive missing a chunk, cut at a chunk boundary or repeating one is refused",
          "[io][archive][audit]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // F1: the reader tracked only the policy chunk, so a file cut after a
    // complete chunk read as one with no curves or surfaces, every checksum
    // intact, and a repeated chunk replaced what the first had decoded.
    for (const int version : {1, 2, 3}) {
        INFO("the committed version-" << version << " file");
        const Bytes golden = FileBytes(GoldenFile(version));
        const std::vector<ChunkSpan> chunks = ChunksOf(golden);
        REQUIRE(chunks.size() == 4);
        REQUIRE(chunks.back().End == golden.size());
        if (version == 3) {
            // The boundaries the audit found.
            REQUIRE(chunks[0].End == 76);
            REQUIRE(chunks[1].End == 1670);
            REQUIRE(chunks[2].End == 2327);
        }
        // Cut after each complete chunk but the last: the next is missing.
        for (std::size_t k = 0; k + 1 < chunks.size(); ++k) {
            const std::string message = ReadError(Slice(golden, 0, chunks[k].End));
            REQUIRE(message.find("no \"" + chunks[k + 1].Tag + "\" chunk") != std::string::npos);
        }
        for (const ChunkSpan& gone : chunks) {
            // Each chunk removed whole.
            Bytes without = Slice(golden, 0, gone.Start);
            without.insert(without.end(), golden.begin() + static_cast<std::ptrdiff_t>(gone.End),
                           golden.end());
            const std::string missing = ReadError(without);
            REQUIRE(missing.find(gone.Tag == "POLI"
                                     ? std::string("no policy chunk")
                                     : "no \"" + gone.Tag + "\" chunk") != std::string::npos);
            // Each chunk repeated as it is.
            Bytes twice = golden;
            twice.insert(twice.end(), golden.begin() + static_cast<std::ptrdiff_t>(gone.Start),
                         golden.begin() + static_cast<std::ptrdiff_t>(gone.End));
            REQUIRE(ReadError(twice).find("a second \"" + gone.Tag + "\" chunk") !=
                    std::string::npos);
        }
        // A valid, empty curve chunk after the real one would have erased
        // every curve.
        const Bytes emptyCurves = Appended(golden, "CRVS", Bytes(8, std::byte{0}));
        REQUIRE(ReadError(emptyCurves).find("a second \"CRVS\" chunk") != std::string::npos);
    }

    // Cut anywhere at all, an archive is refused.
    const Bytes good = Written(EveryType());
    for (std::size_t cut = 0; cut < good.size(); ++cut) {
        const auto read = ReadGeometry(Slice(good, 0, cut));
        REQUIRE_FALSE(read.has_value());
        REQUIRE(read.error().GetCode() == ErrorCode::ParseError);
    }
    // Chunks in another order are still one of each, and read.
    const std::vector<ChunkSpan> chunks = ChunksOf(good);
    Bytes shuffled = Slice(good, 0, 32);
    for (const std::size_t k : {std::size_t{3}, std::size_t{1}, std::size_t{0}, std::size_t{2}}) {
        shuffled.insert(shuffled.end(), good.begin() + static_cast<std::ptrdiff_t>(chunks[k].Start),
                        good.begin() + static_cast<std::ptrdiff_t>(chunks[k].End));
    }
    REQUIRE(Written(Read(shuffled).Archive) == good);
}

TEST_CASE("a record newer than the file's format version is refused", "[io][archive][audit]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // F1: the decoder took any record name it knew, whatever version the
    // file claimed; a version-1 file with a trimmed curve was not written by
    // a version-1 writer.
    const auto relabelled = [](const GeometryArchive& a, const std::uint32_t version) {
        Bytes b = Written(a);
        PutU32At(b, 4, version);
        return b;
    };
    const Curve arc(*Stereon::MakeTrimmedCurve(Curve(Circle{O, 2.5, X, Y}), Interval{5.0, 7.0}));
    GeometryArchive trimmed;
    trimmed.Curves = {arc};
    REQUIRE(ReadError(relabelled(trimmed, 1))
                .find("\"TrimmedCurve\" record, which format version 1 does not have") !=
            std::string::npos);
    REQUIRE(Read(relabelled(trimmed, 2)).FormatVersion == 2);
    GeometryArchive plane;
    plane.Curves2d = {Curve2d(*Stereon::MakeTrimmedCurve(
        Curve2d(Line2d{Point2{0.0, 0.0}, Vec2{1.0, 0.0}}), Interval{0.0, 1.0}))};
    REQUIRE(ReadError(relabelled(plane, 1)).find("\"TrimmedCurve2d\" record") != std::string::npos);

    const std::vector<Surface> composites{
        Surface(*Stereon::MakeTrimmedSurface(Surface(Plane{O, X, Y}),
                                             Box2{Interval{0.0, 1.0}, Interval{0.0, 1.0}})),
        Surface(*Stereon::MakeExtrusion(Curve(Circle{O, 2.5, X, Y}), Vec3{0.0, 0.0, 1.0})),
        Surface(*Stereon::MakeRevolution(Curve(Line{O, X}), Point3{}, Vec3{0.0, 0.0, 1.0})),
    };
    for (const Surface& s : composites) {
        GeometryArchive a;
        a.Surfaces = {s};
        for (const std::uint32_t version : {1U, 2U}) {
            REQUIRE(ReadError(relabelled(a, version)).find("which format version") !=
                    std::string::npos);
        }
        REQUIRE(Read(relabelled(a, 3)).Archive.Surfaces.size() == 1);
    }
    // Nested: a version-2 file may hold a trimmed curve, but not as the
    // profile of a sweep.
    GeometryArchive swept;
    swept.Surfaces = {Surface(*Stereon::MakeExtrusion(arc, Vec3{0.0, 0.0, 1.0}))};
    REQUIRE(ReadError(relabelled(swept, 2)).find("\"ExtrusionSurface\" record") !=
            std::string::npos);
}

TEST_CASE("the stored policy is one the reader applies", "[io][archive][audit]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // F8: the frame allowance was stored and read back, but the factories
    // that validate the geometry checked the compiled FrameAllowance, so a
    // file claimed a validation that never happened; and a resolution above
    // the tolerance ceiling, which no Context can have, was accepted.
    GeometryArchive above = EveryType();
    above.Policy.Resolution = 1e-2;
    above.Policy.ToleranceCeiling = 1e-3;
    REQUIRE(WriteGeometry(above).error().GetCode() == ErrorCode::InvalidInput);
    for (const double allowance : {1e-16, 1e-6}) {
        GeometryArchive other = EveryType();
        other.Policy.FrameAllowance = allowance;
        REQUIRE(WriteGeometry(other).error().GetCode() == ErrorCode::InvalidInput);
    }
    // Any Context's policy can be written.
    GeometryArchive fromContext = EveryType();
    fromContext.Policy = Stereon::PolicyOf(
        Stereon::Context::Default().WithToleranceCeiling(1e-1).WithResolution(1e-3));
    REQUIRE(Read(Written(fromContext)).Archive.Policy == fromContext.Policy);

    // The same refused on read, from a file edited after writing. The policy
    // chunk holds the angular resolution, the tolerance ceiling and the frame
    // allowance, in that order.
    const Bytes good = Written(EveryType());
    const auto withPolicyValue = [&](const std::size_t index, const double value) {
        Bytes b = good;
        const std::size_t policy = ChunkAt(b, "POLI");
        const Bytes word = DoubleBytes(value);
        std::ranges::copy(word, b.begin() + static_cast<std::ptrdiff_t>(policy + 20 + 8 * index));
        Reseal(b, policy);
        return b;
    };
    // A ceiling of 1e-6 below the stored resolution of 2e-5.
    REQUIRE(ReadError(withPolicyValue(1, 1e-6)).find("resolution above") != std::string::npos);
    for (const double allowance : {1e-16, 2e-13, 1e-6}) {
        INFO("frame allowance " << allowance);
        REQUIRE(ReadError(withPolicyValue(2, allowance), ErrorCode::NotSupported)
                    .find("frame allowance") != std::string::npos);
    }
    REQUIRE(Read(withPolicyValue(2, Stereon::FrameAllowance)).Archive.Policy == EveryType().Policy);
}

namespace {

/// A stream buffer that takes every byte and fails when flushed: a disk
/// that fills up after the write call returned.
class FailsOnFlush : public std::streambuf {
protected:
    std::streamsize xsputn(const char* /*bytes*/, const std::streamsize count) override
    {
        return count;
    }
    int_type overflow(const int_type c) override { return traits_type::not_eof(c); }
    int sync() override { return -1; }
};

/// A stream buffer that takes nothing.
class FailsOnWrite : public std::streambuf {
protected:
    std::streamsize xsputn(const char* /*bytes*/, std::streamsize /*count*/) override { return 0; }
    int_type overflow(int_type /*c*/) override { return traits_type::eof(); }
};

} // namespace

TEST_CASE("a save that fails is reported and leaves the file it would replace",
          "[io][archive][audit]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // F7: SaveGeometry checked the stream after write() and returned, so an
    // archive still in the stream's buffer could fail to reach the disk
    // without the caller hearing of it.
    const Bytes bytes = Written(EveryType());
    FailsOnFlush flushing;
    std::ostream full(&flushing);
    REQUIRE_FALSE(Stereon::Detail::WriteAll(full, bytes));
    FailsOnWrite writing;
    std::ostream none(&writing);
    REQUIRE_FALSE(Stereon::Detail::WriteAll(none, bytes));
    std::ostringstream fine;
    REQUIRE(Stereon::Detail::WriteAll(fine, bytes));
    REQUIRE(fine.str().size() == bytes.size());

    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "stereon-save-test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const fs::path path = dir / "model.stn";
    GeometryArchive other;
    other.Curves = {Curve(Line{O, X})};
    const auto entries = [&] {
        return std::distance(fs::directory_iterator(dir), fs::directory_iterator());
    };
    // Replacing a file: the new content, and nothing left beside it.
    REQUIRE(Stereon::SaveGeometry(path, EveryType()).has_value());
    REQUIRE(Stereon::SaveGeometry(path, other).has_value());
    REQUIRE(FileBytes(path) == Written(other));
    REQUIRE(entries() == 1);
    // A target that cannot be replaced, a directory: reported, nothing left.
    fs::create_directory(dir / "taken.stn");
    REQUIRE(Stereon::SaveGeometry(dir / "taken.stn", other).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE(fs::is_directory(dir / "taken.stn"));
    REQUIRE(entries() == 2);
    // A directory that does not exist: nothing written.
    REQUIRE(Stereon::SaveGeometry(dir / "missing" / "model.stn", other).error().GetCode() ==
            ErrorCode::InvalidInput);
    REQUIRE_FALSE(fs::exists(dir / "missing"));
    // Saves racing for one target, from two threads: each writes a file of
    // its own and renames it over the target once whole, so the target is
    // always one archive or the other, never one's bytes in the other's
    // file, and nothing is left beside it.
    constexpr int RaceRounds = 25;
    const GeometryArchive every = EveryType();
    const std::vector<std::byte> first = Written(every);
    const std::vector<std::byte> second = Written(other);
    std::atomic<int> failures = 0;
    std::mutex lock;
    std::string firstError;
    {
        const auto race = [&](const GeometryArchive& archive) {
            for (int round = 0; round < RaceRounds; ++round) {
                const auto saved = Stereon::SaveGeometry(path, archive);
                if (!saved) {
                    ++failures;
                    const std::scoped_lock guard(lock);
                    if (firstError.empty()) {
                        firstError = std::string(saved.error().GetMessage());
                    }
                }
            }
        };
        const std::jthread a(race, std::cref(every));
        const std::jthread b(race, std::cref(other));
    }
    // Two renames racing for one target can fail on Windows, where a file
    // being replaced cannot be replaced again at that instant: a failed
    // save is reported, never a torn file, and the other save goes through.
    INFO("first failure: " << firstError << " (" << failures.load() << " of " << 2 * RaceRounds
                           << ")");
    REQUIRE(failures < 2 * RaceRounds);
    const std::vector<std::byte> published = FileBytes(path);
    REQUIRE((published == first || published == second));
    REQUIRE(entries() == 2);
    fs::remove_all(dir);
}

TEST_CASE("sweeps of trimmed profiles read from a file keep their exact forms",
          "[io][archive][audit]")
{
    if constexpr (!Stereon::SerializationAvailable) {
        return;
    }
    // F2: the committed version-3 file extrudes a circle arc, a trimmed
    // curve, which had no exact form through the erased wrapper.
    const LoadedGeometry loaded = Read(FileBytes(GoldenFile(3)));
    const auto* extrusion = loaded.Archive.Surfaces[8].Target<ExtrusionSurface>();
    REQUIRE(extrusion != nullptr);
    REQUIRE(extrusion->GetProfile().Target<TrimmedCurve>() != nullptr);
    REQUIRE(extrusion->GetProjection() != nullptr);
    const Box2 box{Domain(extrusion->GetProfile()), Interval{-1.0, 1.0}};
    REQUIRE(ToBSpline(*extrusion, box).has_value());
    for (const double u : {5.0, 5.7, 6.4, 7.0}) {
        const Point3 p = Eval(*extrusion, Point2{u, 0.5});
        REQUIRE(Distance(Eval(*extrusion, Project(*extrusion, p)), p) <= 1e-9);
    }
    // The trimmed curve of a trimmed curve has its piece too.
    const auto* nested = loaded.Archive.Curves[8].Target<TrimmedCurve>();
    REQUIRE(nested != nullptr);
    REQUIRE(nested->GetPiece() != nullptr);
}

TEST_CASE("horn tori survive geometry archive validation", "[io][archive][regression]")
{
    GeometryArchive original;
    original.Surfaces = {Surface(Torus{O, 1.0, 1.0, X, Y})};
    if constexpr (!Stereon::SerializationAvailable) {
        REQUIRE(WriteGeometry(original).error().GetCode() == ErrorCode::NotSupported);
        return;
    }
    const Bytes bytes = Written(original);
    const auto loaded = Read(bytes);
    REQUIRE(loaded.Archive.Surfaces.size() == 1);
    REQUIRE(Written(loaded.Archive) == bytes);
    const auto* horn = loaded.Archive.Surfaces[0].Target<Torus>();
    REQUIRE(horn != nullptr);
    REQUIRE(horn->MajorRadius == horn->MinorRadius);
}
