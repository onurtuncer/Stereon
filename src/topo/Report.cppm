// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Topo:Report;

/// The findings of the validity checker (ADR-0018): rule IDs, the entities
/// involved, described, and what was wrong. A report is a list, not a
/// boolean, so that a failure says where to look.

export namespace Stereon {

/// The rules of the validity checker (ADR-0018). The number is the stable
/// rule ID: the Structural rules (ADR-0038's Freeze invariants) from 1, S01
/// being DanglingLink; the Geometric rules from 101, G01 being
/// VertexContainment. Solid (L) rules join as that level arrives.
enum class TopologyRule : std::uint8_t {
    /// A link names a slot out of range or a dead one.
    DanglingLink = 1,
    /// Next and Prev of a loop are not mutual inverses.
    LoopLinksBroken,
    /// A loop walk does not return to its anchor within the coedge count,
    /// or reaches a coedge of another loop.
    LoopOpen,
    /// A live coedge is not reached exactly once by the walk of its loop.
    CoedgeOutsideLoop,
    /// A radial walk does not return to its anchor within the coedge count,
    /// or reaches a use of another edge.
    RadialRingOpen,
    /// A live coedge is not reached exactly once by the radial walk of its
    /// edge.
    CoedgeOutsideRing,
    /// The number of coedges does not fit the edge's curve kind (ADR-0009).
    EdgeUseCount,
    /// A primary coedge is missing where the kind needs one, or is not in
    /// the edge's ring.
    PrimaryCoedge,
    /// A loop is not in its face's list exactly once, or the list names a
    /// loop of another face.
    LoopOutsideFace,
    /// A face is not in its shell's list exactly once, or the list names a
    /// face of another shell.
    FaceOutsideShell,
    /// A shell is not in its solid's list exactly once, or the list names a
    /// shell of another solid.
    ShellOutsideSolid,
    /// Geometry the kind requires is missing or names no pool entry.
    MissingGeometry,
    /// A degenerate edge with two different end vertices.
    EdgeEndpoints,
    /// A loop without coedges.
    EmptyLoop,
    /// The faces of a shell are not connected through shared edges.
    ShellDisconnected,
    /// The vertex-incidence index does not match the edges: a vertex's list
    /// names an edge not at it, twice or out of order, or misses one.
    VertexIncidence,
    /// An edge's cached curve and its certificate disagree with each other
    /// or with the edge (ADR-0009): a cache on a kind that keeps none, a
    /// cache that names no non-empty pool entry, a status without a cache,
    /// or a bound that the status does not allow.
    EdgeCache,

    /// G01 (ADR-0001's containment invariant): an end of an edge lies
    /// outside the tolerance sphere of its vertex.
    VertexContainment = 101,
    /// G02 (ADR-0001's containment invariant): the edge's 3D curve or one of
    /// its pcurves mapped through its face's surface leaves the edge's
    /// tolerance tube, or does not cover the edge's range.
    EdgeContainment,
};

/// "S01", "S02", ...
/// The stable ID of a rule: S01 to S17, G01 and G02.
[[nodiscard]] std::string RuleId(TopologyRule rule);

/// The Geometric rules are numbered from one past this.
inline constexpr int GeometricRuleBase = 100;
/// "DanglingLink", ...
[[nodiscard]] std::string_view RuleName(TopologyRule rule) noexcept;

/// One violated rule: the entities involved, described by Describe(), and
/// what was found.
struct ValidityFinding {
    TopologyRule Rule{};
    std::vector<std::string> Entities;
    std::string Message;
};

/// Everything a check found. A corrupt store can break every rule at every
/// entity; the report keeps the first MaxFindings findings and counts the
/// rest, so that a check of a damaged million-face body stays bounded.
class ValidityReport {
public:
    static constexpr std::size_t MaxFindings = 64;

    void Add(ValidityFinding finding)
    {
        ++m_Count;
        if (m_Findings.size() < MaxFindings) {
            m_Findings.push_back(std::move(finding));
        }
    }

    /// No findings at all.
    [[nodiscard]] bool IsValid() const noexcept { return m_Count == 0; }
    [[nodiscard]] std::span<const ValidityFinding> GetFindings() const noexcept
    {
        return m_Findings;
    }
    /// Findings found, including those not kept.
    [[nodiscard]] std::size_t GetCount() const noexcept { return m_Count; }
    /// True when a finding of `rule` was kept.
    [[nodiscard]] bool Has(TopologyRule rule) const noexcept
    {
        for (const ValidityFinding& finding : m_Findings) {
            if (finding.Rule == rule) {
                return true;
            }
        }
        return false;
    }

    /// One line per kept finding: "S03 LoopOpen: Loop#2 [...]: message".
    [[nodiscard]] std::string ToText() const;

private:
    std::vector<ValidityFinding> m_Findings;
    std::size_t m_Count = 0;
};

} // namespace Stereon
