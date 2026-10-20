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

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

module Stereon.Viewer.Scene;

namespace Stereon::Viewer {

static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

/// The union of two boxes, either of which may be empty.
[[nodiscard]] Box3 Union(const Box3& a, const Box3& b)
{
    if (a.IsEmpty()) {
        return b;
    }
    if (b.IsEmpty()) {
        return a;
    }
    const auto join = [](const Interval& x, const Interval& y) {
        return Interval{std::min(x.GetLower(), y.GetLower()), std::max(x.GetUpper(), y.GetUpper())};
    };
    return Box3{join(a.X, b.X), join(a.Y, b.Y), join(a.Z, b.Z)};
}

/// Where the scene's unbounded geometry goes: around the middle of the
/// bounded content, as far as it reaches.
[[nodiscard]] Placement PlacementFor(const Box3& bounded, const double defaultExtent)
{
    Placement placement;
    placement.Extent = defaultExtent;
    if (bounded.IsEmpty()) {
        return placement;
    }
    const Point3 lo{bounded.X.GetLower(), bounded.Y.GetLower(), bounded.Z.GetLower()};
    const Point3 hi{bounded.X.GetUpper(), bounded.Y.GetUpper(), bounded.Z.GetUpper()};
    placement.Centre = lo + 0.5 * (hi - lo);
    const double reach = 0.5 * Length(hi - lo);
    if (reach > 0.0) {
        placement.Extent = reach;
    }
    return placement;
}

[[nodiscard]] std::filesystem::file_time_type ModifiedTime(const std::filesystem::path& path)
{
    std::error_code ignored;
    return std::filesystem::last_write_time(path, ignored);
}

/// The `.stn` files of the directory: dumps only, never the temporary
/// files a dump is written to before it is renamed.
[[nodiscard]] std::vector<std::filesystem::path> DumpFiles(const std::filesystem::path& directory)
{
    std::vector<std::filesystem::path> files;
    std::error_code error;
    std::filesystem::directory_iterator it(directory, error);
    for (; !error && it != std::filesystem::directory_iterator(); it.increment(error)) {
        if (it->path().extension() == ".stn" && it->is_regular_file(error)) {
            files.push_back(it->path());
        }
    }
    return files;
}

} // namespace

SceneObject LoadObject(const std::filesystem::path& path, const SceneOptions& options)
{
    SceneObject object;
    object.Source = path;
    object.Name = path.stem().string();
    object.Modified = ModifiedTime(path);
    Result<LoadedGeometry> loaded = LoadGeometry(path);
    if (!loaded) {
        object.Error = std::string(loaded.error().GetMessage());
        object.Colour = ColourFor(object.Name);
        return object;
    }
    object.Archive = std::move(loaded->Archive);
    if (object.Archive.Debug) {
        const DebugChunk& debug = *object.Archive.Debug;
        object.Name = debug.Name.empty() ? object.Name : debug.Name;
        object.File = debug.File;
        object.Line = debug.Line;
        object.Time = debug.Time;
    }
    object.Colour = ColourFor(object.Name);
    object.Content = BuildPrimitives(object.Archive, options, object.Colour);
    object.Bounded = BoundsOf(object.Content);
    return object;
}

Scene::Scene(std::filesystem::path directory, SceneOptions options)
    : m_Directory(std::move(directory))
    , m_Options(options)
{}

bool Scene::Refresh()
{
    std::map<std::filesystem::path, SceneObject> previous;
    for (SceneObject& object : m_Objects) {
        previous.emplace(object.Source, std::move(object));
    }
    bool changed = false;
    std::vector<SceneObject> objects;
    for (const std::filesystem::path& path : DumpFiles(m_Directory)) {
        const auto old = previous.find(path);
        if (old != previous.end() && old->second.Modified == ModifiedTime(path)) {
            objects.push_back(std::move(old->second));
            previous.erase(old);
            continue;
        }
        SceneObject object = LoadObject(path, m_Options);
        if (old != previous.end()) {
            object.Visible = old->second.Visible;
            previous.erase(old);
        }
        objects.push_back(std::move(object));
        changed = true;
    }
    changed = changed || !previous.empty(); // files removed
    std::ranges::sort(objects, [](const SceneObject& a, const SceneObject& b) {
        return a.Name != b.Name ? a.Name < b.Name : a.Source < b.Source;
    });
    m_Objects = std::move(objects);
    if (changed) {
        Rebuild();
    }
    return changed;
}

void Scene::Rebuild()
{
    Box3 bounded = Box3::Empty();
    for (const SceneObject& object : m_Objects) {
        bounded = Union(bounded, object.Bounded);
    }
    const Placement placement = PlacementFor(bounded, m_Options.DefaultExtent);
    for (SceneObject& object : m_Objects) {
        if (object.Error.empty() && object.Content.HasUnbounded) {
            object.Content = BuildPrimitives(object.Archive, m_Options, object.Colour, placement);
        }
    }
}

void Scene::SetVisible(const std::size_t index, const bool visible)
{
    if (index < m_Objects.size()) {
        m_Objects[index].Visible = visible;
    }
}

Box3 Scene::GetBounds() const
{
    Box3 bounds = Box3::Empty();
    for (const SceneObject& object : m_Objects) {
        if (object.Visible) {
            bounds = Union(bounds, BoundsOf(object.Content));
        }
    }
    return bounds;
}

} // namespace Stereon::Viewer
