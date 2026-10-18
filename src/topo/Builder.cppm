// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdint>
#include <string>

export module Stereon.Topo:Builder;

import Stereon.Intersect;
import :Records;
import :Report;
import :Shape;

/// ShapeBuilder: the mutable, single-threaded side of ADR-0004, and the
/// controlled construction path of ADR-0027. Its tables are private; every
/// connectivity change goes through an operation below that keeps loop
/// links, radial rings and parent lists consistent together, so no caller
/// writes a link field. Freeze() then runs the independent structural
/// validator and, only if it finds nothing, seals the tables into a
/// certified, immutable Shape without renumbering anything.
///
/// Operations check their handles: a handle that is not a live entity of
/// this builder, or a request the topology cannot take (a coedge of an
/// edge in a loop of another builder, a face added to a second shell), is
/// reported as InvalidInput and changes nothing. Every successful change
/// advances the builder's topology epoch (ADR-0038).

export namespace Stereon {

using BuilderVertexId = BuilderId<EntityKind::Vertex>;
using BuilderEdgeId = BuilderId<EntityKind::Edge>;
using BuilderCoedgeId = BuilderId<EntityKind::Coedge>;
using BuilderLoopId = BuilderId<EntityKind::Loop>;
using BuilderFaceId = BuilderId<EntityKind::Face>;
using BuilderShellId = BuilderId<EntityKind::Shell>;
using BuilderSolidId = BuilderId<EntityKind::Solid>;

/// A new edge. Exact and Spatial edges carry their own 3D curve; every
/// kind has the shared parameter range [t0, t1] of its curve and pcurves
/// (ADR-0009). Start and End are the vertices at t0 and t1; a closed edge
/// has the same vertex at both.
struct EdgeSpec {
    EdgeCurveKind Kind = EdgeCurveKind::Spatial;
    BuilderVertexId Start{};
    BuilderVertexId End{};
    CurveRef Curve{};
    Interval Range{};
    double Tolerance = 0.0;
};

/// The stable identity of one builder (ADR-0038): it moves with the
/// builder, and a moved-from builder has none.
class BuilderToken {
public:
    constexpr BuilderToken() = default;
    constexpr explicit BuilderToken(const StoreId value) noexcept
        : m_Value(value)
    {}
    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Value != 0; }
    [[nodiscard]] constexpr StoreId GetValue() const noexcept { return m_Value; }
    constexpr bool operator==(const BuilderToken&) const = default;

private:
    StoreId m_Value = 0;
};

namespace Testing {
class TopologyFaults;
} // namespace Testing

class ShapeBuilder {
public:
    ShapeBuilder();
    ShapeBuilder(const ShapeBuilder&) = delete;
    ShapeBuilder& operator=(const ShapeBuilder&) = delete;
    /// The token, the epoch and the tables move; the source is left with no
    /// token, and every operation on it fails.
    ShapeBuilder(ShapeBuilder&& other) noexcept;
    ShapeBuilder& operator=(ShapeBuilder&& other) noexcept;
    ~ShapeBuilder();

    // --- Geometry pools ----------------------------------------------------
    // Values are shared by reference: several edges may use one curve, and
    // several faces one surface, each with its own bounds (ADR-0025).

    [[nodiscard]] CurveRef AddCurve(Curve curve);
    [[nodiscard]] PCurveRef AddPCurve(Curve2d pcurve);
    [[nodiscard]] SurfaceRef AddSurface(Surface surface);

    // --- Entities ----------------------------------------------------------

    [[nodiscard]] Result<BuilderVertexId> AddVertex(const Point3& position, double tolerance);

    /// An edge with no uses yet; coedges join its radial ring as loops take
    /// it. Its curve, where the kind has one, must already be in the pool.
    [[nodiscard]] Result<BuilderEdgeId> AddEdge(const EdgeSpec& spec);

    /// A face on a surface of the pool, in no shell and with no loops.
    [[nodiscard]] Result<BuilderFaceId> AddFace(SurfaceRef surface, bool reversed = false);

    /// An empty loop at the end of the face's loop list. A loop must have a
    /// coedge before Freeze().
    [[nodiscard]] Result<BuilderLoopId> AddLoop(BuilderFaceId face);

    /// A use of `edge` at the end of `loop`'s cycle, running along the edge
    /// or against it, with its pcurve on the face's surface. The coedge
    /// joins the edge's radial ring; the first use of an edge of a kind that
    /// evaluates through a surface (Intersection, Trim, Boundary,
    /// Degenerate) becomes the edge's primary coedge (ADR-0009), which
    /// SetPrimaryCoedge() can change explicitly.
    [[nodiscard]] Result<BuilderCoedgeId> AppendCoedge(BuilderLoopId loop, BuilderEdgeId edge,
                                                       EdgeUseSense sense, PCurveRef pcurve);

    /// Makes `coedge`, a use of `edge`, the edge's primary coedge.
    [[nodiscard]] Result<void> SetPrimaryCoedge(BuilderEdgeId edge, BuilderCoedgeId coedge);

    [[nodiscard]] Result<BuilderSolidId> AddSolid();
    /// A shell bounding `solid`, at the end of its shell list.
    [[nodiscard]] Result<BuilderShellId> AddShell(BuilderSolidId solid);
    /// A shell of no solid: a sheet or an open shell.
    [[nodiscard]] Result<BuilderShellId> AddShell();
    /// Adds a face of no shell at the end of the shell's face list.
    [[nodiscard]] Result<void> AddFaceToShell(BuilderShellId shell, BuilderFaceId face);

    // --- Inspection --------------------------------------------------------

    /// The structural check Freeze() will run, without freezing.
    [[nodiscard]] ValidityReport Validate() const;

    /// Live entities of kind K.
    template <EntityKind K>
    [[nodiscard]] std::uint32_t Count() const noexcept
    {
        return m_Tables.template LiveCount<K>();
    }

    [[nodiscard]] BuilderToken GetToken() const noexcept { return m_Token; }
    /// Advanced by every change to connectivity, membership or liveness.
    [[nodiscard]] std::uint64_t GetEpoch() const noexcept { return m_Epoch; }

private:
    friend Result<Shape> Freeze(ShapeBuilder&& builder);
    friend class Testing::TopologyFaults;

    template <EntityKind K>
    [[nodiscard]] bool IsLive(BuilderId<K> h) const noexcept;
    template <EntityKind K>
    [[nodiscard]] Result<void> Require(BuilderId<K> h, const char* operation) const;
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> Allocate(Detail::Record<K> record);
    [[nodiscard]] Result<void> RequireToken(const char* operation) const;
    void Advance() noexcept;

    Detail::TopologyTables m_Tables;
    BuilderToken m_Token;
    std::uint64_t m_Epoch = 0;
};

/// Seals the builder into a Shape: runs the structural validator over every
/// table (ADR-0038's Freeze invariants) and, if it finds nothing, moves the
/// tables into a certified snapshot, keeping every index. Otherwise reports
/// InvalidInput with the findings. Either way the builder is consumed.
[[nodiscard]] Result<Shape> Freeze(ShapeBuilder&& builder);

/// The structural check of a frozen shape, run again by the independent
/// validator: what Freeze() certified, re-derived (ADR-0036).
[[nodiscard]] ValidityReport ValidateTopology(const Shape& shape);

} // namespace Stereon

export namespace Stereon::Detail {

/// The version recorded in every ValidationCertificate.
inline constexpr std::uint32_t StructuralValidatorVersion = 1;

/// The independent structural validator (ADR-0018, ADR-0038): its own
/// bounded walks over the raw tables, checking every link before following
/// it, so it terminates and reads nothing out of range on any input.
[[nodiscard]] ValidityReport ValidateStructure(const TopologyTables& tables);

/// Describe() over raw tables, safe on corrupt ones: it checks every index
/// and bounds every walk, for the validator's findings.
[[nodiscard]] std::string DescribeEntity(const TopologyTables& tables, EntityKind kind,
                                         std::uint32_t index);

} // namespace Stereon::Detail

export namespace Stereon::Testing {

/// Seeded faults for the validator's tests (ADR-0027, section 8): direct
/// writes to link fields that the builder's operations never make. Not for
/// production code, which cannot otherwise reach the tables.
class TopologyFaults {
public:
    static void SetNext(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId next);
    static void SetPrev(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId prev);
    static void SetRadial(ShapeBuilder& b, BuilderCoedgeId c, BuilderCoedgeId radial);
    static void SetCoedgeLoop(ShapeBuilder& b, BuilderCoedgeId c, BuilderLoopId loop);
    static void SetPrimary(ShapeBuilder& b, BuilderEdgeId e, BuilderCoedgeId primary);
    static void SetKind(ShapeBuilder& b, BuilderEdgeId e, EdgeCurveKind kind);
    static void SetEnd(ShapeBuilder& b, BuilderEdgeId e, BuilderVertexId end);
    static void SetPCurve(ShapeBuilder& b, BuilderCoedgeId c, PCurveRef pcurve);
    static void SetFaceShell(ShapeBuilder& b, BuilderFaceId f, BuilderShellId shell);
    static void SetLoopFace(ShapeBuilder& b, BuilderLoopId l, BuilderFaceId face);
    static void SetShellSolid(ShapeBuilder& b, BuilderShellId s, BuilderSolidId solid);
    /// Marks the slot dead without unlinking it.
    static void KillCoedge(ShapeBuilder& b, BuilderCoedgeId c);
    /// Points a link at a slot past the end of its table.
    static void SetNextOutOfRange(ShapeBuilder& b, BuilderCoedgeId c);
};

} // namespace Stereon::Testing
