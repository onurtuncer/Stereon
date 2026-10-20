// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

export module Stereon.Viewer.Scene;

export import Stereon.Io;

/// The debug viewer's scene (ADR-0035, section 3): dump files turned into
/// display primitives, headless, so that it is tested without a window and
/// any frontend can draw it. Geometry is sampled through Stereon.Mesh's
/// display sampling; the debug chunk's primitives are copied; overlays
/// (control polygons and nets, knots, parameter ticks, direction arrows)
/// are added by layer, which a frontend shows or hides.
///
/// Unbounded geometry (a line, a plane, a cylinder) is drawn around the
/// point nearest the middle of the scene's bounded content, as far as that
/// content reaches, or over a default extent when there is none.

export namespace Stereon::Viewer {

/// A colour, components in [0, 1].
struct Rgba {
    float R = 0.0F;
    float G = 0.0F;
    float B = 0.0F;
    float A = 1.0F;

    bool operator==(const Rgba&) const = default;
};

/// What a primitive shows, so that a frontend can toggle each kind.
enum class Layer : std::uint8_t {
    /// Curves, surfaces' iso-lines, plane curves.
    Geometry,
    /// B-spline control polygons and nets.
    ControlPolygon,
    /// Points at the distinct knots of B-splines.
    Knots,
    /// Points at evenly spaced parameters of curves.
    ParameterTicks,
    /// Arrows along curves in their parameter's direction.
    Directions,
    /// The dump's own primitives (its "DBUG" chunk).
    Debug,
};

inline constexpr std::size_t LayerCount = 6;

struct Segment {
    Point3 A;
    Point3 B;
    Rgba Colour;
    Layer Kind = Layer::Geometry;
};

struct Marker {
    Point3 Position;
    Rgba Colour;
    Layer Kind = Layer::Geometry;
};

struct Triangle {
    Point3 A;
    Point3 B;
    Point3 C;
    Rgba Colour;
};

struct Label {
    Point3 Position;
    std::string Text;
    Rgba Colour;
};

/// A segment of the 2D panel: plane curves and pcurves in parameter space.
struct PlaneSegment {
    Point2 A;
    Point2 B;
    Rgba Colour;
    Layer Kind = Layer::Geometry;
};

/// Everything one object draws.
struct Primitives {
    std::vector<Segment> Segments;
    std::vector<Marker> Markers;
    std::vector<Triangle> Triangles;
    std::vector<Label> Labels;
    std::vector<PlaneSegment> PlaneSegments;
    /// What could not be drawn, one line each (a curve sampling refused).
    std::vector<std::string> Warnings;
    /// True when some geometry is unbounded and was drawn around a
    /// placement (or skipped, without one).
    bool HasUnbounded = false;
};

/// How to draw.
struct SceneOptions {
    DisplayTolerance Tolerance;
    /// Iso-lines per direction of a surface.
    std::size_t IsoLines = 9;
    /// Ticks per curve on the ParameterTicks layer.
    std::size_t Ticks = 10;
    /// The extent of unbounded geometry when the scene has nothing bounded.
    double DefaultExtent = 10.0;
};

/// Where unbounded geometry is drawn: around `Centre`, `Extent` each way.
struct Placement {
    Point3 Centre;
    double Extent = 10.0;
};

/// The primitives of an archive. Unbounded geometry is drawn by
/// `placement`, or skipped without one (HasUnbounded tells the caller to
/// ask again once the scene's bounds are known). `colour` is the object's;
/// debug primitives with a colour of their own keep it.
[[nodiscard]] Primitives BuildPrimitives(const GeometryArchive& archive,
                                         const SceneOptions& options, Rgba colour,
                                         std::optional<Placement> placement = std::nullopt);

/// The bounding box of the primitives' 3D points; empty if there are none.
[[nodiscard]] Box3 BoundsOf(const Primitives& primitives);

/// A distinct colour for an object, from its name.
[[nodiscard]] Rgba ColourFor(std::string_view name);

/// One dump file in the scene.
struct SceneObject {
    std::filesystem::path Source;
    /// The dump's name and where it was made, from its debug chunk; the
    /// file's stem and nothing for a file without one.
    std::string Name;
    std::string File;
    std::uint32_t Line = 0;
    std::int64_t Time = 0;
    bool Visible = true;
    Rgba Colour;
    Primitives Content;
    /// The file failed to load: why. Content is empty then.
    std::string Error;
    std::filesystem::file_time_type Modified;
    GeometryArchive Archive;
    /// The bounds of what the object draws without a placement: its
    /// bounded content, which places the scene's unbounded geometry.
    Box3 Bounded;
};

/// The objects of a dump directory, kept in step with it by Refresh.
class Scene {
public:
    explicit Scene(std::filesystem::path directory, SceneOptions options = {});

    /// Loads new and changed `.stn` files of the directory (by modification
    /// time), drops removed ones, and rebuilds the primitives. True when
    /// anything changed. A file that fails to load stays listed with its
    /// error. Visibility survives a reload.
    bool Refresh();

    /// The objects, ordered by name.
    [[nodiscard]] std::span<const SceneObject> GetObjects() const noexcept { return m_Objects; }
    void SetVisible(std::size_t index, bool visible);

    /// The bounds of the visible objects' primitives; empty if none.
    [[nodiscard]] Box3 GetBounds() const;
    [[nodiscard]] const std::filesystem::path& GetDirectory() const noexcept { return m_Directory; }
    [[nodiscard]] const SceneOptions& GetOptions() const noexcept { return m_Options; }

private:
    void Rebuild();

    std::filesystem::path m_Directory;
    SceneOptions m_Options;
    std::vector<SceneObject> m_Objects;
};

/// One object loaded from a file, its primitives built without placement.
[[nodiscard]] SceneObject LoadObject(const std::filesystem::path& path,
                                     const SceneOptions& options);

} // namespace Stereon::Viewer
