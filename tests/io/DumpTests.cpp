// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Dumps for the debug viewer (ADR-0035 verification): the "DBUG" chunk
// round-trips bit for bit, a reader that does not know it loads the
// geometry with one warning, malformed primitives are refused on write and
// on read, STEREON_DUMP_DIR switches dumping on and off, STEREON_DUMP
// records its call site, and a shape is stored with its topology overlays.

#include "../topo/BodiesStd.hpp"
#include <Stereon/Io/Dump.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Io;

#include "../topo/Bodies.hpp"

using Stereon::Circle;
using Stereon::Curve;
using Stereon::DebugBox;
using Stereon::DebugChunk;
using Stereon::DebugColour;
using Stereon::DebugLabel;
using Stereon::DebugPoint;
using Stereon::DebugPolyline;
using Stereon::DebugTriangles;
using Stereon::EdgeCurveKind;
using Stereon::ErrorCode;
using Stereon::GeometryArchive;
using Stereon::Interval;
using Stereon::Point3;
using Stereon::Vec3;

namespace {

using Bytes = std::vector<std::byte>;

constexpr Vec3 X{1.0, 0.0, 0.0};
constexpr Vec3 Y{0.0, 1.0, 0.0};

DebugColour Colour(const std::uint8_t r, const std::uint8_t g, const std::uint8_t b)
{
    DebugColour c;
    c.R = r;
    c.G = g;
    c.B = b;
    c.A = 200;
    return c;
}

/// One primitive of every kind, with colours set and unset.
DebugChunk EveryPrimitive()
{
    DebugChunk d;
    d.Name = "probe";
    d.File = "tests/io/DumpTests.cpp";
    d.Line = 42;
    d.Time = 1'700'000'000'123;
    DebugPoint point;
    point.Position = Point3{1.0, 2.0, 3.0};
    point.Colour = Colour(255, 0, 0);
    d.Points.push_back(point);
    DebugPolyline polyline;
    polyline.Points = {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.5, 0.0}, Point3{2.0, 0.0, -1.0}};
    d.Polylines.push_back(polyline);
    DebugBox box;
    box.Box = Stereon::Box3{Interval{-1.0, 1.0}, Interval{0.0, 0.0}, Interval{2.0, 3.5}};
    box.Colour = Colour(0, 128, 255);
    d.Boxes.push_back(box);
    DebugTriangles triangles;
    triangles.Vertices = {Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}, Point3{0.0, 1.0, 0.0},
                          Point3{1.0, 1.0, 0.0}};
    triangles.Indices = {0, 1, 2, 2, 1, 3};
    d.Triangles.push_back(triangles);
    DebugLabel label;
    label.Position = Point3{0.5, 0.5, 0.0};
    label.Text = "S05 \"quoted\" \\ and a\ttab";
    d.Labels.push_back(label);
    return d;
}

GeometryArchive WithDebug()
{
    GeometryArchive archive;
    archive.Curves.push_back(Curve(Circle{Point3{}, 2.0, X, Y}));
    archive.Debug = EveryPrimitive();
    return archive;
}

Bytes Written(const GeometryArchive& archive)
{
    return TopoTest::Must(Stereon::WriteGeometry(archive));
}

/// The offset of the first chunk with this tag.
std::size_t ChunkOffset(const Bytes& bytes, const std::string_view tag)
{
    for (std::size_t i = 0; i + tag.size() <= bytes.size(); ++i) {
        bool match = true;
        for (std::size_t k = 0; k < tag.size(); ++k) {
            match = match && bytes[i + k] == static_cast<std::byte>(tag[k]);
        }
        if (match) {
            return i;
        }
    }
    FAIL("no chunk " << std::string(tag));
    return 0;
}

/// Changes one payload byte of the chunk at `offset` and fixes its CRC, so
/// that the payload, not the checksum, is what the reader refuses.
void DamagePayload(Bytes& bytes, const std::size_t chunk, const std::size_t at,
                   const std::byte value)
{
    constexpr std::size_t lengthAt = 8;
    constexpr std::size_t crcAt = 16;
    constexpr std::size_t payloadAt = 20;
    std::uint64_t length = 0;
    for (std::size_t i = 0; i < 8; ++i) {
        length |= std::uint64_t{std::to_integer<std::uint8_t>(bytes[chunk + lengthAt + i])}
                  << (8 * i);
    }
    bytes[chunk + payloadAt + at] = value;
    const std::uint32_t crc = Stereon::Detail::Crc32(
        std::span<const std::byte>(bytes).subspan(chunk + payloadAt, length));
    for (std::size_t i = 0; i < 4; ++i) {
        bytes[chunk + crcAt + i] = static_cast<std::byte>((crc >> (8 * i)) & 0xFFU);
    }
}

/// Sets STEREON_DUMP_DIR, or unsets it for an empty value; false if the
/// environment refused. Not a REQUIRE: the guard's destructor calls it.
[[nodiscard]] bool SetDumpDirectory(const std::string& value)
{
#ifdef _WIN32
    return _putenv_s("STEREON_DUMP_DIR", value.c_str()) == 0;
#else
    return value.empty() ? unsetenv("STEREON_DUMP_DIR") == 0
                         : setenv("STEREON_DUMP_DIR", value.c_str(), 1) == 0;
#endif
}

/// Points STEREON_DUMP_DIR at a fresh directory for one test, and removes
/// both afterwards.
class DumpDirectoryGuard {
public:
    DumpDirectoryGuard()
        : m_Path(std::filesystem::temp_directory_path() /
                 ("stereon-dump-test-" + std::to_string(std::random_device{}())))
    {
        REQUIRE(SetDumpDirectory(m_Path.string()));
    }
    DumpDirectoryGuard(const DumpDirectoryGuard&) = delete;
    DumpDirectoryGuard& operator=(const DumpDirectoryGuard&) = delete;
    ~DumpDirectoryGuard()
    {
        static_cast<void>(SetDumpDirectory(""));
        std::error_code ignored;
        std::filesystem::remove_all(m_Path, ignored);
    }

    [[nodiscard]] const std::filesystem::path& GetPath() const { return m_Path; }

private:
    std::filesystem::path m_Path;
};

std::size_t FilesIn(const std::filesystem::path& directory)
{
    std::size_t count = 0;
    for ([[maybe_unused]] const auto& entry : std::filesystem::directory_iterator(directory)) {
        ++count;
    }
    return count;
}

} // namespace

TEST_CASE("a debug chunk round-trips bit for bit", "[io][dump]")
{
    const GeometryArchive archive = WithDebug();
    const Bytes bytes = Written(archive);
    const auto loaded = Stereon::ReadGeometry(bytes);
    REQUIRE(loaded);
    REQUIRE(loaded->Warnings.empty());
    REQUIRE(loaded->Archive.Debug.has_value());
    REQUIRE(*loaded->Archive.Debug == *archive.Debug);
    REQUIRE(loaded->FormatVersion == Stereon::GeometryFormatVersion);
    REQUIRE(Written(loaded->Archive) == bytes);

    // Without a debug chunk the archive is written exactly as before.
    GeometryArchive plain = archive;
    plain.Debug.reset();
    const auto reread = Stereon::ReadGeometry(Written(plain));
    REQUIRE(reread);
    REQUIRE_FALSE(reread->Archive.Debug.has_value());
}

TEST_CASE("a reader that does not know the debug chunk loads the geometry with one warning",
          "[io][dump]")
{
    Bytes bytes = Written(WithDebug());
    const std::size_t chunk = ChunkOffset(bytes, "DBUG");
    bytes[chunk + 3] = std::byte{'X'}; // "DBUX": an optional chunk of an unknown kind
    const auto loaded = Stereon::ReadGeometry(bytes);
    REQUIRE(loaded);
    REQUIRE(loaded->Warnings.size() == 1);
    REQUIRE(loaded->Warnings[0].find("DBUX") != std::string::npos);
    REQUIRE_FALSE(loaded->Archive.Debug.has_value());
    REQUIRE(loaded->Archive.Curves.size() == 1);
}

TEST_CASE("malformed debug primitives are refused on write", "[io][dump]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<GeometryArchive> bad(4, WithDebug());
    bad[0].Debug->Points[0].Position.X = nan;
    bad[1].Debug->Polylines[0].Points[1].Z = std::numeric_limits<double>::infinity();
    bad[2].Debug->Triangles[0].Indices.back() = 4; // four vertices: 0 to 3
    bad[3].Debug->Triangles[0].Indices.pop_back(); // not a multiple of three
    for (const GeometryArchive& archive : bad) {
        const auto written = Stereon::WriteGeometry(archive);
        REQUIRE_FALSE(written);
        REQUIRE(written.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("a damaged debug chunk is refused on read", "[io][dump]")
{
    const Bytes good = Written(WithDebug());
    const std::size_t chunk = ChunkOffset(good, "DBUG");

    Bytes newerLayout = good;
    DamagePayload(newerLayout, chunk, 0, std::byte{2}); // layout version 2
    const auto refused = Stereon::ReadGeometry(newerLayout);
    REQUIRE_FALSE(refused);
    REQUIRE(refused.error().GetCode() == ErrorCode::ParseError);
    REQUIRE(std::string(refused.error().GetMessage()).find("DBUG") != std::string::npos);

    // The name's length (after the u32 layout version) made huge: the
    // reader must not trust it.
    Bytes longName = good;
    DamagePayload(longName, chunk, 7, std::byte{0x7F});
    REQUIRE_FALSE(Stereon::ReadGeometry(longName));

    // A second debug chunk would replace the first.
    Bytes twice = good;
    twice.insert(twice.end(), good.begin() + static_cast<std::ptrdiff_t>(chunk), good.end());
    REQUIRE_FALSE(Stereon::ReadGeometry(twice));
}

TEST_CASE("the JSON view shows the debug chunk", "[io][dump]")
{
    const auto json = Stereon::GeometryToJson(Written(WithDebug()));
    REQUIRE(json);
    for (const std::string_view expected :
         {R"("tag":"DBUG","required":false)", R"("name":"probe")", R"("line":42)",
          R"("indices":[0,1,2,2,1,3])", R"("text":"S05 \"quoted\" \\ and a\u0009tab")",
          R"("colour":[0,128,255,200])", R"("colour":null)"}) {
        INFO(expected);
        REQUIRE(json->find(expected) != std::string::npos);
    }
}

TEST_CASE("dump names are checked", "[io][dump]")
{
    for (const std::string_view good : {"a", "after-insertion", "edge_4.v2", "A0"}) {
        REQUIRE(Stereon::IsDumpName(good));
    }
    for (const std::string_view bad : {"", ".hidden", "a/b", "a\\b", "a b", "..", "naïve"}) {
        INFO(bad);
        REQUIRE_FALSE(Stereon::IsDumpName(bad));
    }
}

TEST_CASE("without STEREON_DUMP_DIR a dump writes nothing", "[io][dump]")
{
    REQUIRE(SetDumpDirectory(""));
    REQUIRE_FALSE(Stereon::DumpDirectory().has_value());
    const Curve circle(Circle{Point3{}, 1.0, X, Y});
    const auto written = Stereon::DumpAt("circle", circle);
    REQUIRE(written);
    REQUIRE_FALSE(written->has_value());
    REQUIRE_FALSE(Stereon::Dump("circle", circle));
    // A bad name is refused even when nothing would be written.
    REQUIRE_FALSE(Stereon::DumpAt("a/b", circle));
}

TEST_CASE("a dump writes <name>.stn with where and when it was made", "[io][dump]")
{
    const DumpDirectoryGuard guard;
    const Curve circle(Circle{Point3{}, 1.0, X, Y});
    const auto written = Stereon::DumpAt("circle", circle);
    REQUIRE(written);
    REQUIRE(written->has_value());
    REQUIRE(**written == guard.GetPath() / "circle.stn");

    const auto loaded = Stereon::LoadGeometry(**written);
    REQUIRE(loaded);
    REQUIRE(loaded->Archive.Curves.size() == 1);
    const DebugChunk& debug = loaded->Archive.Debug.value();
    REQUIRE(debug.Name == "circle");
    REQUIRE(debug.File.find("DumpTests.cpp") != std::string::npos);
    REQUIRE(debug.Line > 0);
    REQUIRE(debug.Time > 0);

    // A second dump of the name replaces the first and leaves no temporary
    // file behind.
    REQUIRE(Stereon::Dump("circle", Curve(Circle{Point3{}, 3.0, X, Y})));
    REQUIRE(FilesIn(guard.GetPath()) == 1);

    REQUIRE_FALSE(Stereon::DumpAt("a b", circle));
    REQUIRE_FALSE(Stereon::Dump(nullptr, circle));
}

TEST_CASE("STEREON_DUMP records the line it is called on", "[io][dump]")
{
    const DumpDirectoryGuard guard;
    const Curve circle(Circle{Point3{}, 1.0, X, Y});
    const std::uint32_t line = __LINE__ + 1;
    STEREON_DUMP("macro", circle);
    const auto loaded = Stereon::LoadGeometry(guard.GetPath() / "macro.stn");
    REQUIRE(loaded);
    REQUIRE(loaded->Archive.Debug->Line == line);
    REQUIRE(loaded->Archive.Debug->File.find("DumpTests.cpp") != std::string::npos);

    STEREON_DUMP("primitives", EveryPrimitive());
    const auto primitives = Stereon::LoadGeometry(guard.GetPath() / "primitives.stn");
    REQUIRE(primitives);
    REQUIRE(primitives->Archive.Debug->Labels == EveryPrimitive().Labels);
    REQUIRE(primitives->Archive.Debug->Name == "primitives");
}

TEST_CASE("a shape is dumped with its geometry and topology overlays", "[io][dump]")
{
    const Stereon::Shape box = TopoTest::MakeBox().Freeze();
    const auto archive = Stereon::ShapeArchive(box);
    REQUIRE(archive);
    REQUIRE(archive->Surfaces.size() == 6);
    REQUIRE(archive->Curves.size() == 12);
    REQUIRE(archive->Curves2d.size() == 24);
    const DebugChunk& debug = archive->Debug.value();
    REQUIRE(debug.Points.size() == 8);
    REQUIRE(debug.Boxes.size() == 8);
    REQUIRE(debug.Labels.size() == 8 + 12);
    REQUIRE(debug.Polylines.size() == 12);
    for (const DebugPolyline& edge : debug.Polylines) {
        REQUIRE(edge.Colour == Stereon::EdgeCurveColour(EdgeCurveKind::Exact));
    }
    for (const DebugBox& tolerance : debug.Boxes) {
        REQUIRE(tolerance.Box.X.GetLower() < tolerance.Box.X.GetUpper());
    }

    // It is an ordinary archive: written, read back and viewed as JSON.
    const auto loaded = Stereon::ReadGeometry(Written(*archive));
    REQUIRE(loaded);
    REQUIRE(loaded->Archive.Debug == archive->Debug);
    REQUIRE(Stereon::GeometryToJson(Written(*archive)));

    // The cylinder's edges are of two kinds, each in its colour.
    const auto cylinder = Stereon::ShapeArchive(TopoTest::MakeCylinder().Freeze());
    REQUIRE(cylinder);
    std::vector<DebugColour> colours;
    for (const DebugPolyline& edge : cylinder->Debug->Polylines) {
        colours.push_back(edge.Colour);
    }
    REQUIRE(std::ranges::count(colours, Stereon::EdgeCurveColour(EdgeCurveKind::Exact)) == 2);
    REQUIRE(std::ranges::count(colours, Stereon::EdgeCurveColour(EdgeCurveKind::Boundary)) == 1);

    // Faces are drawn where they are: each surface trimmed to the box its
    // boundary covers, so the cylinder's lateral face ends at its two circles.
    for (const Stereon::Surface& face : cylinder->Surfaces) {
        REQUIRE(face.Target<Stereon::TrimmedSurface>() != nullptr);
    }
    const Stereon::Box2 lateral = Stereon::Domain(cylinder->Surfaces.front());
    REQUIRE(lateral.Y.GetLower() == 0.0);
    REQUIRE(lateral.Y.GetUpper() == 4.0);

    REQUIRE_FALSE(Stereon::ShapeArchive(Stereon::Shape{}));

    const DumpDirectoryGuard guard;
    REQUIRE(Stereon::Dump("box", box));
    REQUIRE(std::filesystem::exists(guard.GetPath() / "box.stn"));
}
