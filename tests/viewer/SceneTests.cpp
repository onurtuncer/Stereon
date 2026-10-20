// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The debug viewer's scene (ADR-0035 verification): dump files turned into
// the expected display primitives, overlays by layer, unbounded geometry
// placed by the scene's bounded content, and a dump directory followed as
// files are added, changed and removed. No window is opened.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <random>
#include <string>
#include <system_error>
#include <vector>

import Stereon.Viewer.Scene;

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Circle;
using Stereon::Circle2d;
using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::DebugBox;
using Stereon::DebugChunk;
using Stereon::DebugColour;
using Stereon::DebugLabel;
using Stereon::DebugTriangles;
using Stereon::GeometryArchive;
using Stereon::Interval;
using Stereon::Line;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sphere;
using Stereon::Surface;
using Stereon::Vec2;
using Stereon::Vec3;
using Stereon::Viewer::Layer;
using Stereon::Viewer::Primitives;
using Stereon::Viewer::Rgba;
using Stereon::Viewer::SceneOptions;

namespace {

constexpr Vec3 X{1.0, 0.0, 0.0};
constexpr Vec3 Y{0.0, 1.0, 0.0};
constexpr Rgba Red{1.0F, 0.0F, 0.0F, 1.0F};

template <typename T>
std::size_t CountLayer(const std::vector<T>& items, const Layer layer)
{
    return static_cast<std::size_t>(
        std::ranges::count_if(items, [&](const T& item) { return item.Kind == layer; }));
}

GeometryArchive Holding(const Curve& curve)
{
    GeometryArchive archive;
    archive.Curves.push_back(curve);
    return archive;
}

BSplineCurve Cubic()
{
    return BSplineCurve{3,
                        {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0},
                        {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.0}, Point3{2.0, -1.0, 1.0},
                         Point3{3.0, 1.0, 2.0}, Point3{4.0, 0.0, 0.0}},
                        {}};
}

/// A fresh directory, removed afterwards.
class TemporaryDirectory {
public:
    TemporaryDirectory()
        : m_Path(std::filesystem::temp_directory_path() /
                 ("stereon-scene-test-" + std::to_string(std::random_device{}())))
    {
        std::filesystem::create_directories(m_Path);
    }
    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
    ~TemporaryDirectory()
    {
        std::error_code ignored;
        std::filesystem::remove_all(m_Path, ignored);
    }
    [[nodiscard]] const std::filesystem::path& GetPath() const { return m_Path; }

private:
    std::filesystem::path m_Path;
};

/// Saves an archive whose debug chunk names it, as a dump would.
void SaveNamed(const std::filesystem::path& path, GeometryArchive archive, const std::string& name)
{
    archive.Debug.emplace();
    archive.Debug->Name = name;
    archive.Debug->File = "tests/viewer/SceneTests.cpp";
    archive.Debug->Line = 7;
    REQUIRE(Stereon::SaveGeometry(path, archive));
}

/// Moves a file's modification time forward, so that a rewrite within the
/// file system's time resolution still reads as a change.
void Touch(const std::filesystem::path& path)
{
    const auto when = std::filesystem::last_write_time(path) + std::chrono::seconds(2);
    std::filesystem::last_write_time(path, when);
}

} // namespace

TEST_CASE("a curve is drawn with its ticks and its direction", "[viewer][scene]")
{
    SceneOptions options;
    options.Ticks = 7;
    const Primitives p =
        Stereon::Viewer::BuildPrimitives(Holding(Curve(Circle{Point3{}, 2.0, X, Y})), options, Red);
    REQUIRE(CountLayer(p.Segments, Layer::Geometry) > 10);
    REQUIRE(CountLayer(p.Markers, Layer::ParameterTicks) == 7);
    REQUIRE(CountLayer(p.Segments, Layer::Directions) == 3);
    REQUIRE(CountLayer(p.Segments, Layer::ControlPolygon) == 0);
    REQUIRE(CountLayer(p.Markers, Layer::Knots) == 0);
    REQUIRE_FALSE(p.HasUnbounded);
    REQUIRE(p.Warnings.empty());
    for (const auto& s : p.Segments) {
        if (s.Kind == Layer::Geometry) {
            REQUIRE(s.Colour == Red);
        }
    }
}

TEST_CASE("a B-spline curve shows its control polygon and knots", "[viewer][scene]")
{
    const Primitives p =
        Stereon::Viewer::BuildPrimitives(Holding(Curve(Cubic())), SceneOptions{}, Red);
    REQUIRE(CountLayer(p.Segments, Layer::ControlPolygon) == 4); // five poles
    REQUIRE(CountLayer(p.Markers, Layer::Knots) == 3);           // 0, 0.5 and 1
}

TEST_CASE("a B-spline surface shows its iso-lines and control net", "[viewer][scene]")
{
    GeometryArchive archive;
    archive.Surfaces.push_back(Surface(
        BSplineSurface{1,
                       1,
                       {0.0, 0.0, 0.5, 1.0, 1.0},
                       {0.0, 0.0, 1.0, 1.0},
                       {Point3{0.0, 0.0, 0.0}, Point3{0.0, 1.0, 0.0}, Point3{1.0, 0.0, 1.0},
                        Point3{1.0, 1.0, 1.0}, Point3{2.0, 0.0, 0.0}, Point3{2.0, 1.0, 0.0}},
                       {}}));
    const Primitives p = Stereon::Viewer::BuildPrimitives(archive, SceneOptions{}, Red);
    // Three rows of two poles: three segments along v, four along u.
    REQUIRE(CountLayer(p.Segments, Layer::ControlPolygon) == 3 + 4);
    REQUIRE(CountLayer(p.Segments, Layer::Geometry) > 0);
}

TEST_CASE("unbounded geometry waits for a placement", "[viewer][scene]")
{
    const GeometryArchive line = Holding(Curve(Line{Point3{0.0, 5.0, 0.0}, X}));
    const Primitives unplaced = Stereon::Viewer::BuildPrimitives(line, SceneOptions{}, Red);
    REQUIRE(unplaced.HasUnbounded);
    REQUIRE(unplaced.Segments.empty());

    Stereon::Viewer::Placement placement;
    placement.Centre = Point3{3.0, 0.0, 0.0};
    placement.Extent = 2.0;
    const Primitives placed =
        Stereon::Viewer::BuildPrimitives(line, SceneOptions{}, Red, placement);
    REQUIRE(placed.HasUnbounded);
    REQUIRE_FALSE(placed.Segments.empty());
    const Stereon::Box3 bounds = Stereon::Viewer::BoundsOf(placed);
    // Around the point of the line nearest the centre, two units each way.
    REQUIRE(bounds.X.GetLower() >= 1.0 - 1e-9);
    REQUIRE(bounds.X.GetUpper() <= 5.0 + 1e-9);

    GeometryArchive sphere;
    sphere.Surfaces.push_back(Surface(Sphere{Point3{}, 1.0, X, Y}));
    REQUIRE_FALSE(Stereon::Viewer::BuildPrimitives(sphere, SceneOptions{}, Red).HasUnbounded);
}

TEST_CASE("plane curves go to the 2D panel", "[viewer][scene]")
{
    GeometryArchive archive;
    archive.Curves2d.push_back(Curve2d(Circle2d{Point2{0.5, 0.5}, 0.25, Vec2{1, 0}, Vec2{0, 1}}));
    const Primitives p = Stereon::Viewer::BuildPrimitives(archive, SceneOptions{}, Red);
    REQUIRE(p.Segments.empty());
    REQUIRE(CountLayer(p.PlaneSegments, Layer::Geometry) > 10);
    REQUIRE(CountLayer(p.PlaneSegments, Layer::Directions) == 2);
}

TEST_CASE("the debug chunk's primitives are drawn in their own colours", "[viewer][scene]")
{
    GeometryArchive archive;
    DebugChunk& debug = archive.Debug.emplace();
    DebugBox box;
    box.Box = Stereon::Box3{Interval{0.0, 1.0}, Interval{0.0, 2.0}, Interval{0.0, 3.0}};
    box.Colour = DebugColour{0, 255, 0, 255};
    debug.Boxes.push_back(box);
    DebugTriangles triangles;
    triangles.Vertices = {Point3{}, Point3{1, 0, 0}, Point3{0, 1, 0}, Point3{1, 1, 0}};
    triangles.Indices = {0, 1, 2, 2, 1, 3};
    debug.Triangles.push_back(triangles);
    DebugLabel label;
    label.Text = "Vertex#3";
    debug.Labels.push_back(label);

    const Primitives p = Stereon::Viewer::BuildPrimitives(archive, SceneOptions{}, Red);
    REQUIRE(CountLayer(p.Segments, Layer::Debug) == 12);
    REQUIRE(p.Segments.front().Colour == Rgba{0.0F, 1.0F, 0.0F, 1.0F});
    REQUIRE(p.Triangles.size() == 2);
    REQUIRE(p.Triangles.front().Colour == Red); // no colour of its own
    REQUIRE(p.Labels.size() == 1);
    REQUIRE(p.Labels.front().Text == "Vertex#3");

    const Stereon::Box3 bounds = Stereon::Viewer::BoundsOf(p);
    REQUIRE(bounds.Z.GetUpper() == 3.0);
    REQUIRE(Stereon::Viewer::BoundsOf(Primitives{}).IsEmpty());
}

TEST_CASE("object colours are stable and distinct", "[viewer][scene]")
{
    REQUIRE(Stereon::Viewer::ColourFor("edge") == Stereon::Viewer::ColourFor("edge"));
    REQUIRE_FALSE(Stereon::Viewer::ColourFor("edge") == Stereon::Viewer::ColourFor("face"));
}

TEST_CASE("a scene follows its dump directory", "[viewer][scene]")
{
    const TemporaryDirectory directory;
    const std::filesystem::path circle = directory.GetPath() / "circle.stn";
    SaveNamed(circle, Holding(Curve(Circle{Point3{}, 2.0, X, Y})), "b-circle");
    SaveNamed(directory.GetPath() / "cubic.stn", Holding(Curve(Cubic())), "a-cubic");
    // A dump's temporary file and anything not .stn are not objects.
    std::ofstream(directory.GetPath() / "half.stn.partial") << "partial";
    std::ofstream(directory.GetPath() / "notes.txt") << "notes";

    Stereon::Viewer::Scene scene(directory.GetPath());
    REQUIRE(scene.Refresh());
    REQUIRE(scene.GetObjects().size() == 2);
    REQUIRE(scene.GetObjects()[0].Name == "a-cubic");
    REQUIRE(scene.GetObjects()[1].Name == "b-circle");
    REQUIRE(scene.GetObjects()[1].Line == 7);
    REQUIRE_FALSE(scene.Refresh()); // nothing changed

    scene.SetVisible(1, false);
    SaveNamed(circle, Holding(Curve(Circle{Point3{}, 6.0, X, Y})), "b-circle");
    Touch(circle);
    REQUIRE(scene.Refresh());
    REQUIRE_FALSE(scene.GetObjects()[1].Visible); // visibility survives the reload
    REQUIRE(Stereon::Viewer::BoundsOf(scene.GetObjects()[1].Content).X.GetUpper() > 5.0);
    // The hidden circle no longer counts towards the scene's bounds.
    REQUIRE(scene.GetBounds().X.GetUpper() < 4.5);

    std::filesystem::remove(circle);
    REQUIRE(scene.Refresh());
    REQUIRE(scene.GetObjects().size() == 1);
}

TEST_CASE("a damaged file is listed with its error", "[viewer][scene]")
{
    const TemporaryDirectory directory;
    std::ofstream(directory.GetPath() / "broken.stn", std::ios::binary) << "not an archive";
    Stereon::Viewer::Scene scene(directory.GetPath());
    REQUIRE(scene.Refresh());
    REQUIRE(scene.GetObjects().size() == 1);
    REQUIRE(scene.GetObjects()[0].Name == "broken");
    REQUIRE_FALSE(scene.GetObjects()[0].Error.empty());
    REQUIRE(scene.GetObjects()[0].Content.Segments.empty());
}

TEST_CASE("unbounded geometry is drawn where the bounded content is", "[viewer][scene]")
{
    const TemporaryDirectory directory;
    SaveNamed(directory.GetPath() / "circle.stn",
              Holding(Curve(Circle{Point3{10.0, 0.0, 0.0}, 1.0, X, Y})), "circle");
    SaveNamed(directory.GetPath() / "line.stn", Holding(Curve(Line{Point3{}, X})), "line");
    Stereon::Viewer::Scene scene(directory.GetPath());
    REQUIRE(scene.Refresh());
    const auto& line = scene.GetObjects()[1];
    REQUIRE(line.Name == "line");
    REQUIRE_FALSE(line.Content.Segments.empty());
    const Stereon::Box3 bounds = Stereon::Viewer::BoundsOf(line.Content);
    // Around x = 10, about as far as the circle reaches.
    REQUIRE(bounds.X.GetLower() > 7.0);
    REQUIRE(bounds.X.GetUpper() < 13.0);
}
