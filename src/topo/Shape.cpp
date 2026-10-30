// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <format>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

module Stereon.Topo;

namespace Stereon {

namespace Detail {

void ShapeHandleViolation(const EntityKind kind, const std::uint32_t index,
                          const std::uint32_t generation, const StoreId store,
                          const StoreId expectedStore, const std::uint32_t size, const bool alive)
{
    HandleCheckFailed({
        .Kind = kind,
        .Index = index,
        .Generation = generation,
        .Store = store,
        .ExpectedStore = expectedStore,
        .Size = size,
        .SlotGeneration = 0,
        .SlotAlive = alive,
    });
}

namespace {

/// The registry behind RegisterSnapshot(), built on first use.
struct SnapshotRegistry {
    std::mutex Mutex;
    std::unordered_map<StoreId, const TopologyTables*> Tables;
    /// The text stereon_describe_handle() returns, valid until its next call.
    std::string Text;
};

SnapshotRegistry& Registry()
{
    static SnapshotRegistry registry;
    return registry;
}

} // namespace

void RegisterSnapshot(const StoreId store, const TopologyTables* tables) noexcept
{
    SnapshotRegistry& r = Registry();
    const std::lock_guard lock(r.Mutex);
    r.Tables[store] = tables;
}

void UnregisterSnapshot(const StoreId store) noexcept
{
    SnapshotRegistry& r = Registry();
    const std::lock_guard lock(r.Mutex);
    r.Tables.erase(store);
}

/// The registered description of slot `index` of kind `kind` in `store`, or
/// null when no live snapshot has that store.
[[nodiscard]] const char* DescribeRegistered(const int kind, const std::uint32_t index,
                                             const StoreId store)
{
    if (kind < 0 || static_cast<std::size_t>(kind) >= EntityKindCount) {
        return nullptr;
    }
    SnapshotRegistry& r = Registry();
    const std::lock_guard lock(r.Mutex);
    const auto found = r.Tables.find(store);
    if (found == r.Tables.end()) {
        return nullptr;
    }
    r.Text = DescribeEntity(*found->second, static_cast<EntityKind>(kind), index);
    return r.Text.c_str();
}

/// As DescribeRegistered, into the caller's buffer, NUL-terminated and cut
/// to capacity - 1 characters; the full length, 0 when no live snapshot has
/// the store or the kind is not one. Built and copied under the registry's
/// lock, so nothing of the registry leaves it (ADR-0046, section 7).
[[nodiscard]] std::size_t DescribeRegisteredInto(const int kind, const std::uint32_t index,
                                                 const StoreId store, char* const buffer,
                                                 const std::size_t capacity)
{
    if (kind < 0 || static_cast<std::size_t>(kind) >= EntityKindCount) {
        return 0;
    }
    SnapshotRegistry& r = Registry();
    const std::lock_guard lock(r.Mutex);
    const auto found = r.Tables.find(store);
    if (found == r.Tables.end()) {
        return 0;
    }
    const std::string text = DescribeEntity(*found->second, static_cast<EntityKind>(kind), index);
    if (buffer != nullptr && capacity > 0) {
        const std::size_t copied = std::min(text.size(), capacity - 1);
        std::memcpy(buffer, text.data(), copied);
        buffer[copied] = '\0';
    }
    return text.size();
}

} // namespace Detail

std::optional<Curve> CurveOf(const Shape& shape, const EdgeId e)
{
    const CurveRef curve = Detail::RecordOf(shape, e).Curve;
    if (!curve.IsValid()) {
        return std::nullopt;
    }
    return shape.GetSnapshot()->GetTables().Curves[curve.GetIndex()];
}

std::optional<Curve> CachedCurveOf(const Shape& shape, const EdgeId e)
{
    const CurveRef curve = Detail::RecordOf(shape, e).Cached;
    if (!curve.IsValid()) {
        return std::nullopt;
    }
    return shape.GetSnapshot()->GetTables().Curves[curve.GetIndex()];
}

std::optional<Curve2d> PCurveOf(const Shape& shape, const CoedgeId c)
{
    const PCurveRef pcurve = Detail::RecordOf(shape, c).PCurve;
    if (!pcurve.IsValid()) {
        return std::nullopt;
    }
    return shape.GetSnapshot()->GetTables().PCurves[pcurve.GetIndex()];
}

Surface SurfaceOf(const Shape& shape, const FaceId f)
{
    return shape.GetSnapshot()->GetTables().Surfaces[Detail::RecordOf(shape, f).Surface.GetIndex()];
}

Result<LoopRange> LoopsOf(const Shape& shape, const FaceId face)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, face, "LoopsOf"));
    return LoopRange(shape.GetSnapshot(), Detail::RecordOf(shape, face).FirstLoop.GetIndex());
}

Result<ShellFaceRange> FacesOf(const Shape& shape, const ShellId shell)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, shell, "FacesOf"));
    return ShellFaceRange(shape.GetSnapshot(), Detail::RecordOf(shape, shell).FirstFace.GetIndex());
}

Result<SolidShellRange> ShellsOf(const Shape& shape, const SolidId solid)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, solid, "ShellsOf"));
    return SolidShellRange(shape.GetSnapshot(),
                           Detail::RecordOf(shape, solid).FirstShell.GetIndex());
}

Result<CoedgeRange> CoedgesOf(const Shape& shape, const LoopId loop)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, loop, "CoedgesOf"));
    return CoedgeRange(shape.GetSnapshot(), Detail::RecordOf(shape, loop).Anchor.GetIndex());
}

Result<RadialCoedgeRange> RadialCoedgesOf(const Shape& shape, const EdgeId edge)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, edge, "RadialCoedgesOf"));
    return RadialCoedgeRange(shape.GetSnapshot(),
                             Detail::RecordOf(shape, edge).RadialAnchor.GetIndex());
}

Result<VertexEdgeRange> EdgesOf(const Shape& shape, const VertexId vertex)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, vertex, "EdgesOf"));
    return VertexEdgeRange(shape.GetSnapshot(), vertex.GetIndex(),
                           Detail::RecordOf(shape, vertex).FirstEdge.GetIndex());
}

Result<std::vector<FaceId>> FacesOf(const Shape& shape, const EdgeId edge)
{
    STEREON_TRY_VOID(Detail::CheckQuery(shape, edge, "FacesOf"));
    const TopologySnapshot& snapshot = *shape.GetSnapshot();
    const Detail::TopologyTables& t = snapshot.GetTables();
    std::vector<FaceId> faces;
    const std::uint32_t anchor = t.Edges[edge.GetIndex()].RadialAnchor.GetIndex();
    if (anchor != CoedgeId::Invalid) {
        std::uint32_t c = anchor;
        do {
            faces.push_back(snapshot.MakeId<EntityKind::Face>(
                t.Loops[t.Coedges[c].Loop.GetIndex()].Face.GetIndex()));
            c = t.Coedges[c].Radial.GetIndex();
        } while (c != anchor);
    }
    std::ranges::sort(faces);
    const auto [first, last] = std::ranges::unique(faces);
    faces.erase(first, last);
    return faces;
}

std::string RuleId(const TopologyRule rule)
{
    const int number = static_cast<int>(rule);
    return number > GeometricRuleBase ? std::format("G{:02}", number - GeometricRuleBase)
                                      : std::format("S{:02}", number);
}

std::string_view RuleName(const TopologyRule rule) noexcept
{
    // In rule order: S01 first.
    static constexpr std::array<std::string_view, 17> Names{
        "DanglingLink",    "LoopLinksBroken",   "LoopOpen",          "CoedgeOutsideLoop",
        "RadialRingOpen",  "CoedgeOutsideRing", "EdgeUseCount",      "PrimaryCoedge",
        "LoopOutsideFace", "FaceOutsideShell",  "ShellOutsideSolid", "MissingGeometry",
        "EdgeEndpoints",   "EmptyLoop",         "ShellDisconnected", "VertexIncidence",
        "EdgeCache",
    };
    static constexpr std::array<std::string_view, 2> GeometricNames{
        "VertexContainment",
        "EdgeContainment",
    };
    const auto number = static_cast<std::size_t>(rule);
    constexpr auto geometric = static_cast<std::size_t>(GeometricRuleBase);
    if (number > geometric) {
        const std::size_t index = number - geometric - 1;
        return index < GeometricNames.size() ? GeometricNames[index]
                                             : std::string_view{"UnknownRule"};
    }
    return number - 1 < Names.size() ? Names[number - 1] : std::string_view{"UnknownRule"};
}

std::string ValidityReport::ToText() const
{
    std::string text;
    for (const ValidityFinding& finding : m_Findings) {
        text += RuleId(finding.Rule);
        text += ' ';
        text += RuleName(finding.Rule);
        text += ':';
        for (const std::string& entity : finding.Entities) {
            text += ' ';
            text += entity;
        }
        text += ": ";
        text += finding.Message;
        text += '\n';
    }
    if (m_Count > m_Findings.size()) {
        text += std::format("... and {} more\n", m_Count - m_Findings.size());
    }
    return text;
}

} // namespace Stereon

/// For a debugger (ADR-0004, safeguard 4): the Describe() text of slot
/// `index` of the entity kind `kind` (EntityKind's order) in the live
/// snapshot with store `store`, or null when there is none or the build
/// checks no handles. The GDB and LLDB pretty-printers call it to show the
/// entity behind a handle; the text stays valid until the next call. That
/// text is a process-wide buffer, so this is the one named legacy exception
/// of ADR-0046, section 7: for serialised debugger invocation only, with no
/// concurrent consumers. Every other consumer uses
/// stereon_describe_handle_into().
extern "C" const char* stereon_describe_handle(const int kind, const unsigned int index,
                                               const unsigned int store)
{
    if constexpr (STEREON_CHECKED_HANDLES) {
        return Stereon::Detail::DescribeRegistered(kind, index, store);
    }
    return nullptr;
}

/// The description stereon_describe_handle() gives, copied into the caller's
/// `buffer` of `capacity` bytes, NUL-terminated and cut to fit; returns the
/// full length, 0 when there is no such entity or the build checks no
/// handles. Safe to call concurrently: the copy is made under the registry's
/// lock and nothing of the registry escapes it (ADR-0046, section 7).
extern "C" std::size_t stereon_describe_handle_into(const int kind, const unsigned int index,
                                                    const unsigned int store, char* const buffer,
                                                    const std::size_t capacity)
{
    if constexpr (STEREON_CHECKED_HANDLES) {
        return Stereon::Detail::DescribeRegisteredInto(kind, index, store, buffer, capacity);
    }
    return 0;
}
