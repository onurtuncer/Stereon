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

#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <string>
#include <string_view>

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

} // namespace Detail

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

std::string RuleId(const TopologyRule rule)
{
    return std::format("S{:02}", static_cast<int>(rule));
}

std::string_view RuleName(const TopologyRule rule) noexcept
{
    // In rule order: S01 first.
    static constexpr std::array<std::string_view, 15> Names{
        "DanglingLink",    "LoopLinksBroken",   "LoopOpen",          "CoedgeOutsideLoop",
        "RadialRingOpen",  "CoedgeOutsideRing", "EdgeUseCount",      "PrimaryCoedge",
        "LoopOutsideFace", "FaceOutsideShell",  "ShellOutsideSolid", "MissingGeometry",
        "EdgeEndpoints",   "EmptyLoop",         "ShellDisconnected",
    };
    const auto index = static_cast<std::size_t>(rule) - 1;
    return index < Names.size() ? Names[index] : std::string_view{"UnknownRule"};
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
