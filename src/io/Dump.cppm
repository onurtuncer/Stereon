// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <filesystem>
#include <optional>
#include <string_view>

export module Stereon.Io:Dump;

import Stereon.Intersect;
import Stereon.Mesh;
import :Debug;
import :GeometryArchive;

/// Dumps for the debug viewer (ADR-0035, section 1): one ordinary .stn file
/// per call, `<STEREON_DUMP_DIR>/<name>.stn`, holding the geometry chunks
/// and a "DBUG" chunk with the dump's name, where and when it was made, and
/// its debug primitives. With STEREON_DUMP_DIR unset or empty a dump writes
/// nothing and costs one environment lookup. A dump is written to a
/// temporary file and renamed (SaveGeometry), so a reader never sees half a
/// file, and a later dump of the same name replaces the earlier one.
///
/// `STEREON_DUMP(name, object)` in <Stereon/Io/Dump.hpp> records the call's
/// source location; the `Dump(name, object)` overloads, non-template and
/// out of line, can be called from a GDB or LLDB breakpoint.
///
/// Tests, tools and applications dump; kernel algorithms below io cannot
/// (ADR-0035: that would go through a sink on the Context, decided when
/// first needed).

export namespace Stereon {

/// The environment variable that names the dump directory.
inline constexpr std::string_view DumpDirectoryVariable = "STEREON_DUMP_DIR";

/// The dump directory, or nothing when STEREON_DUMP_DIR is unset or empty.
[[nodiscard]] std::optional<std::filesystem::path> DumpDirectory();

/// True for a name a dump can have: non-empty, of `[A-Za-z0-9_.-]` only, and
/// not starting with a dot.
[[nodiscard]] bool IsDumpName(std::string_view name) noexcept;

/// Writes the archive as `<name>.stn` in the dump directory, creating the
/// directory if needed. Its DBUG chunk (added if it has none) gets the name,
/// `where` and the current time. The path written, or nothing when dumping
/// is off; InvalidInput for a name IsDumpName refuses or a file that cannot
/// be written, and WriteGeometry's errors for what it cannot store.
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, GeometryArchive archive,
       SourceLocation where = SourceLocation::Current());

[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, const Curve& curve, SourceLocation where = SourceLocation::Current());
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, const Curve2d& curve,
       SourceLocation where = SourceLocation::Current());
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, const Surface& surface,
       SourceLocation where = SourceLocation::Current());
/// A shape, as ShapeArchive stores it.
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, const Shape& shape, SourceLocation where = SourceLocation::Current());
/// A shape with a validity report's findings drawn on it, as
/// ValidityArchive stores it.
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, const Shape& shape, const ValidityReport& report,
       SourceLocation where = SourceLocation::Current());
/// Debug primitives alone.
[[nodiscard]] Result<std::optional<std::filesystem::path>>
DumpAt(std::string_view name, DebugChunk primitives,
       SourceLocation where = SourceLocation::Current());

/// What a dump of the shape stores, until the native format has topology
/// chunks (ADR-0012): every face as its surface's iso-lines inside its
/// loops, polylines in FaceColour (SampleTrimmedSurface), or, where its
/// loops cannot be sampled, its surface trimmed to the parameter box its
/// boundary's pcurves cover; the curve of every edge trimmed to the edge's
/// range, the pcurve of every coedge, and as debug
/// primitives the topology overlays of ADR-0035 section 4: each vertex as a
/// point with its tolerance as a box, each edge as a polyline coloured by
/// its curve kind (EdgeCurveColour), drawn from its 3D curve or, for a kind
/// without one, from its primary coedge's pcurve on the face's surface, and
/// a Describe() label at each vertex and at the middle of each edge.
/// InvalidInput for an empty shape.
[[nodiscard]] Result<GeometryArchive> ShapeArchive(const Shape& shape);

/// The colour ShapeArchive draws an edge of the kind in.
[[nodiscard]] DebugColour EdgeCurveColour(EdgeCurveKind kind) noexcept;

/// The colour ShapeArchive draws a face's iso-lines in.
[[nodiscard]] DebugColour FaceColour() noexcept;

/// The shape as ShapeArchive stores it, with the report's findings drawn
/// on it (ADR-0035 section 4): every vertex, edge and face a finding names
/// in FindingColour, a coedge as its edge, drawn before the shape's own
/// primitives so that they show on top, and a label at the first of them
/// with the rule, the entities' short names and the message, "ShellOpen
/// (Shell#0, Edge#8): the shell uses the edge 1 times, not twice: it is not
/// closed". The labels are the findings' alone: the shape's Describe()
/// labels would bury them. Entities are found by their Describe() text, as
/// the validator names them; one no longer in the shape is skipped, and the
/// finding's label is still written, at the origin when it names nothing
/// drawn. SelfIntersectionUnchecked, which every Solid-level report holds,
/// is left out; a report with nothing else gives ShapeArchive's archive.
/// Errors as ShapeArchive.
[[nodiscard]] Result<GeometryArchive> ValidityArchive(const Shape& shape,
                                                      const ValidityReport& report);

/// The colour ValidityArchive draws a finding's entities in.
[[nodiscard]] DebugColour FindingColour() noexcept;

/// A surface–surface intersection as a dump: both surfaces, each trimmed to
/// its requested box when the box is bounded (an unbounded one is drawn
/// around the scene), every branch's pcurves on A and on B for the
/// Parameter space window, and as debug primitives every branch's 3D curve
/// over its range in its contact's colour (ContactColour), every isolated
/// point, and a label for each, "branch 0: Crossing, OppositeSide, bound
/// 2.1e-07" or "point 0: bound inf". An overlap region is labelled at the
/// first point of its first boundary branch, or at A's origin. Errors from
/// sampling a branch.
[[nodiscard]] Result<GeometryArchive> IntersectionArchive(const Surface& a, const Box2& boxA,
                                                          const Surface& b, const Box2& boxB,
                                                          const SurfaceSurfaceIntersection& result);

/// A curve–surface intersection as a dump: the curve over its domain when
/// bounded, the surface over its box when bounded, every point as a debug
/// point in its contact's colour with a label, "t = 0.25: Crossing", and
/// every overlap as the curve's polyline over its range in the overlap's
/// colour, with its pcurve for the Parameter space window. Errors from
/// sampling an overlap.
[[nodiscard]] Result<GeometryArchive> IntersectionArchive(const Curve& curve,
                                                          const Interval& domain,
                                                          const Surface& surface, const Box2& box,
                                                          const CurveSurfaceIntersection& result);

/// The colour IntersectionArchive draws a branch or point of the contact
/// in: red for a crossing, orange for a touch, green for an overlap.
[[nodiscard]] DebugColour ContactColour(ContactType contact) noexcept;

/// For a debugger: DumpAt without a source location, true when a file was
/// written. A null name, dumping switched off and any error give false.
bool Dump(const char* name, const Curve& curve);
bool Dump(const char* name, const Curve2d& curve);
bool Dump(const char* name, const Surface& surface);
bool Dump(const char* name, const Shape& shape);

} // namespace Stereon
