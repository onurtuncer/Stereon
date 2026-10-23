// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

export module Stereon.Topo:Builder;

import Stereon.Intersect;
import :Operation;
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
///
/// A builder opened on a frozen Shape (Open) starts from the shape's
/// tables and shares every page with it until it writes to one. Deletion
/// (the Remove operations) leaves tombstones: a removed slot is dead, its
/// generation bumped so that every handle to it is stale, and never reused
/// by this builder, so an index live in both the shape and its edit names
/// the same entity. Compact() is the separate, explicit renumbering.

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

/// The geometry of an edge an Euler operator makes: an EdgeSpec without its
/// vertices, which the operator decides, and the pcurve of both new uses.
/// MEV and MEF make the edge's two uses on one surface (a strut in one
/// loop, or the two sides of a face split), so one pcurve serves both; it
/// may be empty only for a Spatial edge.
struct EdgeGeometry {
    EdgeCurveKind Kind = EdgeCurveKind::Spatial;
    CurveRef Curve{};
    Interval Range{};
    double Tolerance = 0.0;
    PCurveRef PCurve{};
};

/// What MVFS makes: an isolated vertex, and a face with one empty loop in a
/// new shell. The loop is the vertex's loop: MEV grows its first edge from
/// that vertex.
struct VertexFaceShell {
    BuilderVertexId Vertex{};
    BuilderFaceId Face{};
    BuilderLoopId Loop{};
    BuilderShellId Shell{};
};

/// What MEV makes: an edge from an existing vertex to a new one, and its
/// two uses in one loop, out along the edge and back against it.
struct EdgeVertex {
    BuilderEdgeId Edge{};
    BuilderVertexId Vertex{};
    BuilderCoedgeId Out{};
    BuilderCoedgeId Back{};
};

/// What MEF makes: an edge splitting a loop, and the new face and loop on
/// its other side. `Kept` is the edge's use in the original loop, `New` its
/// use in the new one.
struct EdgeFace {
    BuilderEdgeId Edge{};
    BuilderFaceId Face{};
    BuilderLoopId Loop{};
    BuilderCoedgeId Kept{};
    BuilderCoedgeId New{};
};

/// SEMV keeps the input edge on the lower range and creates Upper and Vertex.
struct EdgeSplit {
    BuilderEdgeId Upper{};
    BuilderVertexId Vertex{};
};

/// MEKR joins two loops with two uses of a new bridge edge.
struct EdgeRing {
    BuilderEdgeId Edge{};
    BuilderCoedgeId Out{};
    BuilderCoedgeId Back{};
};

/// One input vertex and its deterministic representative after MergeVertices.
struct VertexMergeEntry {
    BuilderVertexId Source{};
    BuilderVertexId Representative{};
};
struct VertexMergeResult {
    std::vector<VertexMergeEntry> Entries;
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

class Compacted;
[[nodiscard]] Result<Compacted> Compact(const Shape& shape);
[[nodiscard]] Result<OperationResult> Commit(ShapeBuilder&& builder, const Context& context);

class ShapeBuilder {
public:
    ShapeBuilder();
    /// A builder on the shape's tables, sharing every page with it until it
    /// writes to one (ADR-0005). Its handles to the shape's entities come
    /// from Handle(); new entities take slots after the shape's last one.
    /// InvalidInput for an empty Shape.
    [[nodiscard]] static Result<ShapeBuilder> Open(const Shape& shape);
    ShapeBuilder(const ShapeBuilder&) = delete;
    ShapeBuilder& operator=(const ShapeBuilder&) = delete;
    /// The token, the epoch and the tables move; the source is left with no
    /// token, and every operation on it fails.
    ShapeBuilder(ShapeBuilder&& other) noexcept;
    ShapeBuilder& operator=(ShapeBuilder&& other) noexcept;
    ~ShapeBuilder();

    // --- Geometry pools ----------------------------------------------------
    // Values are shared by reference: several edges may use one curve, and
    // several faces one surface, each with its own bounds (ADR-0025). An
    // empty wrapper supports nothing and is refused with InvalidInput, as is
    // any addition to a builder that was moved from or frozen.

    [[nodiscard]] Result<CurveRef> AddCurve(Curve curve);
    [[nodiscard]] Result<PCurveRef> AddPCurve(Curve2d pcurve);
    [[nodiscard]] Result<SurfaceRef> AddSurface(Surface surface);

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

    // --- Editing -----------------------------------------------------------
    // Each operation unlinks what it removes from every list and ring it is
    // in and tombstones the slot; nothing else changes.

    /// Replaces the face's surface by another of the pool. The pcurves of its
    /// coedges are the caller's to keep consistent. Drops the cached curve of
    /// every edge the face uses, which was S(c(t)) on the old surface.
    [[nodiscard]] Result<void> SetSurface(BuilderFaceId face, SurfaceRef surface);

    /// Sets an entity's tolerance, the explicit operation through which
    /// ADR-0001 lets a local tolerance change, and logs the change with its
    /// reason. Commit() reports it in the OperationResult, warns when it grew,
    /// and fails with ToleranceExceeded for a tolerance above the Context's
    /// ceiling or, when the Context forbids growth, for any that grew.
    /// InvalidInput unless the tolerance is finite and not negative.
    [[nodiscard]] Result<void> SetTolerance(BuilderVertexId vertex, double tolerance,
                                            std::string_view reason);
    [[nodiscard]] Result<void> SetTolerance(BuilderEdgeId edge, double tolerance,
                                            std::string_view reason);

    /// Builds the cached 3D curve of an Intersection or Trim edge (ADR-0009,
    /// section 5) from its primary coedge: S(c(t)) sampled at CacheSamples
    /// parameters spread evenly over the edge's range and interpolated by a
    /// cubic B-spline on the edge's own parameter, so the cache passes
    /// through the primary side at every sample. Deterministic: the same
    /// pcurve and surface give the same cache. The cache's deviation from
    /// S(c(t)) is then certified within the edge's tolerance by
    /// CertifyDeviation (ADR-0009, section 4) where interval evaluation
    /// covers the types, and the record holds the proven bound; otherwise
    /// the cache is CacheCertificate::Uncertified. Any later change to a
    /// pcurve, the primary coedge or a surface of the edge drops it.
    /// InvalidInput for an edge of another kind, one without a primary
    /// coedge, or a range the pcurve or surface does not cover.
    [[nodiscard]] Result<CurveRef> CacheEdgeCurve(BuilderEdgeId edge);

    /// The samples CacheEdgeCurve() interpolates.
    static constexpr std::size_t CacheSamples = 33;

    /// Removes one use of an edge from its loop and its edge's radial ring.
    /// Refused when the coedge is its edge's primary and the edge has other
    /// uses: the caller reassigns the primary first (ADR-0027, section 4).
    [[nodiscard]] Result<void> RemoveCoedge(BuilderCoedgeId coedge);
    /// Removes a loop with its coedges from its face. A primary coedge among
    /// them passes the primary to the next use in the ring, since every use
    /// of this loop goes; an edge left with no use loses its primary.
    [[nodiscard]] Result<void> RemoveLoop(BuilderLoopId loop);
    /// Removes a face with its loops and coedges from its shell, as
    /// RemoveLoop does for each loop. Its edges and vertices stay.
    [[nodiscard]] Result<void> RemoveFace(BuilderFaceId face);
    /// Removes a shell from its solid; its faces stay, in no shell.
    [[nodiscard]] Result<void> RemoveShell(BuilderShellId shell);
    /// Removes a solid; its shells stay, of no solid.
    [[nodiscard]] Result<void> RemoveSolid(BuilderSolidId solid);
    /// Removes an edge with no uses from its vertices' lists; refused while
    /// a coedge uses it.
    [[nodiscard]] Result<void> RemoveEdge(BuilderEdgeId edge);
    /// Removes a vertex no edge starts or ends at; refused otherwise.
    [[nodiscard]] Result<void> RemoveVertex(BuilderVertexId vertex);

    // --- Euler operators (ADR-0027) -----------------------------------------
    // The classical set, under Stereon's contracts. Each operator checks all
    // of its preconditions before it writes anything, so a refusal
    // (InvalidInput) leaves the builder exactly as it was. Each states its
    // inverse, which restores the connectivity up to the correspondence of
    // the entities that survive: deleted slots are tombstones and never
    // reused, so the restored entities have new handles. Faces keep their
    // shells, and a new face joins the shell of the face it came from. In
    // checked builds each operator then checks the loops, radial rings and
    // vertex lists it touched (ADR-0027, section 9), a contract.

    /// MVFS: an isolated vertex at `position`, and a face on `surface` with
    /// one empty loop, the vertex's, in a new shell bounding `solid` or, if
    /// `solid` is invalid, of no solid. Inverse: KillVertexFaceShell.
    [[nodiscard]] Result<VertexFaceShell> MakeVertexFaceShell(const Point3& position,
                                                              double tolerance, SurfaceRef surface,
                                                              bool reversed = false,
                                                              BuilderSolidId solid = {});

    /// KVFS: removes an isolated vertex, and a face whose only loop is empty
    /// and which is the only face of its shell, with the loop and the shell.
    /// The solid, if any, stays. Inverse of MakeVertexFaceShell.
    [[nodiscard]] Result<void> KillVertexFaceShell(BuilderVertexId vertex, BuilderFaceId face);

    /// MEV: a new vertex at `position` and an edge to it from the vertex at
    /// which `at` ends, as a strut of `at`'s loop: the edge's use along it
    /// (Out) and its use against it (Back) follow `at` in the cycle. The
    /// edge runs from the existing vertex to the new one. Refused for a
    /// Degenerate edge, whose ends are one vertex. Inverse: KillEdgeVertex.
    [[nodiscard]] Result<EdgeVertex> MakeEdgeVertex(BuilderCoedgeId at, const Point3& position,
                                                    double tolerance, const EdgeGeometry& geometry);

    /// MEV on an empty loop: the first edge, from `from`, which must be an
    /// isolated vertex (the loop's vertex, as MVFS made it), to a new one.
    [[nodiscard]] Result<EdgeVertex> MakeEdgeVertex(BuilderLoopId loop, BuilderVertexId from,
                                                    const Point3& position, double tolerance,
                                                    const EdgeGeometry& geometry);

    /// KEV: removes a strut edge and its far vertex: the edge's only uses
    /// are adjacent in one loop and meet at `vertex`, which no other edge
    /// starts or ends at. A loop left with no coedge stays, empty, the loop
    /// of the edge's other vertex. Inverse of MakeEdgeVertex.
    [[nodiscard]] Result<void> KillEdgeVertex(BuilderEdgeId edge, BuilderVertexId vertex);

    /// MEF: an edge from the vertex at which `a` ends to the vertex at which
    /// `b` ends, both coedges of one loop, splitting the loop in two. The
    /// original loop keeps `a`, followed by the edge's use along it (Kept),
    /// and the coedges after `b`; a new loop of a new face, on the same
    /// surface and in the same shell, takes the edge's use against it (New)
    /// and the coedges after `a` up to `b`. With `a` equal to `b` the edge
    /// is closed and its own new loop. The face's other loops stay with it:
    /// moving a hole to the new face is a separate operation. Inverse:
    /// KillEdgeFace.
    [[nodiscard]] Result<EdgeFace> MakeEdgeFace(BuilderCoedgeId a, BuilderCoedgeId b,
                                                const EdgeGeometry& geometry);

    /// KEF: removes an edge whose only two uses lie in two faces on the same
    /// surface, the same way round and in the same shell, and removes
    /// `face`, one of them, joining its loop into the other face's: the
    /// coedges of `face`'s loop follow those of the kept loop in the cycle.
    /// `face` must have that one loop. Inverse of MakeEdgeFace.
    [[nodiscard]] Result<void> KillEdgeFace(BuilderEdgeId edge, BuilderFaceId face);

    /// KEMR: remove a bridge with exactly two opposite uses in one loop,
    /// leaving two nonempty cycles on the same face. The anchor's side stays
    /// in the original loop. Refuses adjacent uses (use KEV for a strut).
    /// Returns the new loop; inverse: MakeEdgeKillRing.
    [[nodiscard]] Result<BuilderLoopId> KillEdgeMakeRing(BuilderEdgeId edge);
    /// MEKR: bridge the ends of a and b in distinct nonempty loops of one
    /// face. Keeps a's loop, deletes b's loop. Inverse: KillEdgeMakeRing.
    [[nodiscard]] Result<EdgeRing> MakeEdgeKillRing(BuilderCoedgeId a, BuilderCoedgeId b,
                                                    const EdgeGeometry& geometry);
    /// KFMRH: move dead's sole nonempty loop into kept and remove dead.
    /// Faces must share surface, orientation and shell. No geometric nesting
    /// is inferred. Inverse: MakeFaceKillRingHole on the moved loop.
    [[nodiscard]] Result<void> KillFaceMakeRingHole(BuilderFaceId dead, BuilderFaceId kept);
    /// MFKRH: move a nonempty loop of a multi-loop face into a new face on
    /// the same surface and shell. Refuses to disconnect a shell; faces
    /// outside a shell are supported. Other loops stay. Inverse: KFMRH.
    [[nodiscard]] Result<BuilderFaceId> MakeFaceKillRingHole(BuilderLoopId loop);
    /// SEMV: split every radial use at t, strictly inside the bounded range.
    /// The caller supplies the split position and vertex tolerance; every
    /// incident trace at t must lie inside it. Shared curves/pcurves keep
    /// their parameterisation and are restricted by the two edge ranges.
    /// Cached certificates remain valid on each restricted range. Supports
    /// wires, seams and arbitrary radial valence, but not Degenerate edges.
    /// Inverse: JoinEdgeKillVertex(lower, upper).
    [[nodiscard]] Result<EdgeSplit> SplitEdgeMakeVertex(BuilderEdgeId edge, double t,
                                                        const Point3& position, double tolerance);
    /// JEKV: join contiguous lower/upper ranges of the same geometry at
    /// their sole degree-two vertex. Every use must be paired adjacently in
    /// the same loop, sense and pcurve; primary uses must correspond.
    /// Keeps lower, removes upper and the common vertex. Different cached
    /// representations are dropped; identical certificates are preserved.
    [[nodiscard]] Result<void> JoinEdgeKillVertex(BuilderEdgeId lower, BuilderEdgeId upper);
    /// ADR-0001: sort unique candidates by index, cluster against each
    /// representative's original position/tolerance (never another member),
    /// and retain its position. Tolerance growth obeys ctx and is logged.
    /// Redirect all incident edges, including closed and non-manifold ones;
    /// reject a merge if the resulting tolerance cannot contain their ends.
    /// Atomic: failure leaves connectivity, geometry and events unchanged.
    [[nodiscard]] Result<VertexMergeResult> MergeVertices(std::span<const BuilderVertexId> vertices,
                                                          const Context& ctx);
    /// Transaction-local low-level events, not final snapshot references.
    /// The enclosing operation supplies semantic history (ADR-0030).
    [[nodiscard]] std::span<const TopologyEvent> GetTopologyEvents() const noexcept
    {
        return m_Events;
    }

    // --- Transactions (ADR-0027, section 9) ---------------------------------

    /// Runs a compound edit on a candidate copy of this builder, which shares
    /// every page with it until it writes to one (ADR-0005), and keeps the
    /// copy only if the edit succeeds: on failure the candidate is discarded
    /// and this builder is exactly as it was. The edit takes the candidate,
    /// a ShapeBuilder&, and returns a Result; handles it obtains are valid in
    /// this builder after a successful transaction, and handles of this
    /// builder are valid in the candidate. After a failed transaction the
    /// handles the edit obtained stay invalid: its slot indices are reserved,
    /// never reused (ADR-0040), and the epoch advances (ADR-0038). This also
    /// holds when the edit throws. Failure cleanup does not allocate; any
    /// missing tombstones are materialised before the next allocation.
    template <typename Edit>
    [[nodiscard]] auto Transaction(Edit&& edit) -> std::invoke_result_t<Edit&, ShapeBuilder&>
    {
        ShapeBuilder candidate(*this, CandidateKey{});
        TransactionScope scope{*this, candidate};
        auto result = edit(candidate);
        if (result.has_value()) {
            *this = std::move(candidate);
            scope.Committed = true;
        }
        return result;
    }

    // --- Inspection --------------------------------------------------------

    /// The builder's handle to an entity of the shape it was opened on: the
    /// same slot, at its current generation. Invalid when the handle names
    /// no live entity of this builder, in checked builds also when it comes
    /// from another shape.
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> Handle(const Id<K> h) const noexcept
    {
        if (!h.IsValid() || !m_Tables.template IsLive<K>(h.GetIndex())) {
            return {};
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            if (h.GetStore() != m_SourceStore) {
                return {};
            }
        }
        return BuilderId<K>::Make(h.GetIndex(),
                                  m_Tables.template Table<K>()[h.GetIndex()].Slot.GetGeneration(),
                                  m_Token.GetValue());
    }

    /// The loop's coedges in cycle order from its anchor; none for an empty
    /// loop or a handle that is not a live loop of this builder.
    [[nodiscard]] std::vector<BuilderCoedgeId> CoedgesOf(BuilderLoopId loop) const;

    /// Connectivity read back during an edit, for the operations that
    /// compose the operators (the Booleans' imprinting, M6): a coedge's
    /// edge, loop and sense, an edge's vertices, a vertex's position, and a
    /// face's loops in order. An invalid handle, Same, the origin or an
    /// empty list for a handle that is not live in this builder. Geometry
    /// beyond positions is read from the frozen Shape the edit started from.
    [[nodiscard]] BuilderEdgeId EdgeOf(BuilderCoedgeId coedge) const noexcept;
    [[nodiscard]] BuilderLoopId LoopOf(BuilderCoedgeId coedge) const noexcept;
    [[nodiscard]] EdgeUseSense SenseOf(BuilderCoedgeId coedge) const noexcept;
    [[nodiscard]] BuilderVertexId StartOf(BuilderEdgeId edge) const noexcept;
    [[nodiscard]] BuilderVertexId EndOf(BuilderEdgeId edge) const noexcept;
    [[nodiscard]] Point3 PositionOf(BuilderVertexId vertex) const noexcept;
    [[nodiscard]] std::vector<BuilderLoopId> LoopsOf(BuilderFaceId face) const;

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

    /// A builder on given tables, for Open() and Compact(): `source` is the
    /// store the tables came from, which Handle() checks.
    ShapeBuilder(Detail::TopologyTables tables, StoreId source, Detail::FreezeKey);

private:
    friend Result<Shape> Freeze(ShapeBuilder&& builder);
    friend Result<OperationResult> Commit(ShapeBuilder&& builder, const Context& context);
    friend class Testing::TopologyFaults;

    /// What CacheEdgeCurve proves about a cache (ADR-0009, section 4).
    struct CacheBound {
        CacheCertificate Status;
        double Bound;
    };
    [[nodiscard]] static CacheBound CertifyCache(const Curve& cache, const Surface& surface,
                                                 const Curve2d& pcurve, const Interval& range,
                                                 double tolerance);

    /// A transaction's candidate: a copy with the same token, so that the
    /// handles of either are valid in the other, which only Transaction()
    /// makes and only one of which survives it.
    struct CandidateKey {};
    ShapeBuilder(const ShapeBuilder& other, CandidateKey);
    /// Reserves the discarded candidate's indices and advances the epoch,
    /// without allocating, including during exception unwinding.
    void Abandon(const ShapeBuilder& candidate) noexcept;
    struct TransactionScope {
        ShapeBuilder& Original;
        const ShapeBuilder& Candidate;
        bool Committed = false;
        ~TransactionScope() noexcept
        {
            if (!Committed) {
                Original.Abandon(Candidate);
            }
        }
    };

    // The Euler operators' shared steps (Euler.cpp).
    [[nodiscard]] Result<void> CheckEdgeGeometry(const EdgeGeometry& geometry,
                                                 const char* operation) const;
    [[nodiscard]] std::uint32_t EndVertexOf(std::uint32_t coedge) const noexcept;
    [[nodiscard]] std::uint32_t SecondOfTwoUses(std::uint32_t edge) const noexcept;
    [[nodiscard]] bool IsFreeEnd(std::uint32_t edge, std::uint32_t vertex) const noexcept;
    [[nodiscard]] std::uint32_t FaceOfUse(std::uint32_t coedge) const noexcept;
    [[nodiscard]] Result<void> CheckJoinable(std::uint32_t dead, std::uint32_t alive,
                                             const char* operation) const;
    [[nodiscard]] std::vector<std::uint32_t> CycleFrom(std::uint32_t coedge) const;
    void RelinkCycle(std::uint32_t loop, const std::vector<std::uint32_t>& cycle);
    [[nodiscard]] BuilderCoedgeId NewUse(std::uint32_t edge, std::uint32_t loop, EdgeUseSense sense,
                                         PCurveRef pcurve);
    [[nodiscard]] Result<EdgeVertex> GrowStrut(std::uint32_t loop, std::uint32_t after,
                                               std::uint32_t from, const Point3& position,
                                               double tolerance, const EdgeGeometry& geometry);
    void KillEdgeAndUses(std::uint32_t edge);
    void CheckNeighbourhood(const std::vector<std::uint32_t>& loops,
                            const std::vector<std::uint32_t>& edges) const;

    [[nodiscard]] std::vector<std::uint32_t> UsesOf(std::uint32_t edge) const;
    [[nodiscard]] bool CanSeparateLoop(std::uint32_t loop) const;
    void MoveLoop(std::uint32_t loop, std::uint32_t face);
    void RepointEdge(std::uint32_t edge, std::uint32_t start, std::uint32_t end);
    [[nodiscard]] Result<void> CheckPointAt(std::uint32_t edge, double t, const Point3& position,
                                            double tolerance) const;
    [[nodiscard]] Result<void> CheckEdgeJoin(std::uint32_t lower, std::uint32_t upper) const;
    void SplitUses(std::uint32_t lower, std::uint32_t upper);
    [[nodiscard]] Result<VertexMergeResult>
    MergeVerticesImpl(std::span<const BuilderVertexId> vertices, const Context& ctx);
    [[nodiscard]] Result<void> MergeVertexCluster(const std::vector<BuilderVertexId>& cluster,
                                                  const Context& ctx);
    [[nodiscard]] Result<void> CheckUsePointAt(std::uint32_t c, double t, const Point3& position,
                                               double tolerance) const;
    [[nodiscard]] Result<void> CheckJoinVertex(std::uint32_t lower, std::uint32_t upper) const;
    [[nodiscard]] Result<void> CheckJoinUses(std::uint32_t lower, std::uint32_t upper) const;
    [[nodiscard]] Result<void> CheckMergeEndpoints(const std::vector<BuilderVertexId>& cluster,
                                                   const Point3& position, double tolerance) const;
    void CheckVertexNeighbourhood(std::uint32_t vertex) const;
    void RecordEvent(EntityKind kind, TopologyEventKind event, std::uint32_t source,
                     std::uint32_t target);

    template <EntityKind K>
    [[nodiscard]] bool IsLive(BuilderId<K> h) const noexcept;
    /// The builder handle of a live slot, at its current generation.
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> HandleOfSlot(std::uint32_t index) const noexcept;
    template <EntityKind K>
    [[nodiscard]] Result<void> Require(BuilderId<K> h, const char* operation) const;
    template <EntityKind K>
    [[nodiscard]] BuilderId<K> Allocate(Detail::Record<K> record);
    template <EntityKind K>
    void Release(std::uint32_t index);
    [[nodiscard]] Result<void> RequireToken(const char* operation) const;
    void Advance() noexcept;

    void LinkAtVertex(std::uint32_t edge, std::uint32_t vertex, bool atStart);
    void UnlinkAtVertex(std::uint32_t edge, std::uint32_t vertex);
    /// Drops the edge's cached curve, if it has one: a pcurve, the primary
    /// coedge or a surface it was built from has changed (ADR-0004).
    void DropCache(std::uint32_t edge);
    template <EntityKind K>
    [[nodiscard]] Result<void> SetToleranceOf(BuilderId<K> h, double tolerance,
                                              std::string_view reason, const char* operation);
    void DetachCoedge(std::uint32_t coedge, bool passPrimaryOn);
    void RemoveLoopAt(std::uint32_t loop);
    template <EntityKind P, EntityKind C>
    void UnlinkChild(std::uint32_t child, Id<P> Detail::Record<C>::* parent,
                     Id<C> Detail::Record<C>::* next, Id<C> Detail::Record<P>::* first,
                     Id<C> Detail::Record<P>::* last);

    Detail::TopologyTables m_Tables;
    BuilderToken m_Token;
    // Exclusive high-water marks, including indices issued by aborted
    // nested edits. Moves leave the source's marks intact for its edit scope.
    std::array<std::size_t, EntityKindCount> m_ReservedEntities{};
    std::array<std::size_t, 3> m_ReservedGeometry{};
    std::uint64_t m_Epoch = 0;
    StoreId m_SourceStore = 0;
    /// The snapshot Open() started from, which Commit()'s History compares
    /// the result with; null for a new builder.
    std::shared_ptr<const TopologySnapshot> m_Source;
    /// Every SetTolerance() call, in order, for Commit() to report.
    std::vector<ToleranceChange> m_ToleranceChanges;
    std::vector<TopologyEvent> m_Events;
};

/// The old-to-new slot map of one Compact() (ADR-0004), which feeds the
/// history of ADR-0011: for every kind, the new index of each old slot,
/// Invalid for a dead one, and the same for the three geometry pools, whose
/// unreferenced entries are dropped.
class RemapTable {
public:
    template <EntityKind K>
    [[nodiscard]] std::uint32_t NewIndex(const std::uint32_t old) const noexcept
    {
        const std::vector<std::uint32_t>& map = m_Entities[static_cast<std::size_t>(K)];
        return old < map.size() ? map[old] : Id<K>::Invalid;
    }
    [[nodiscard]] std::uint32_t NewCurve(const std::uint32_t old) const noexcept
    {
        return old < m_Curves.size() ? m_Curves[old] : CurveRef::Invalid;
    }
    [[nodiscard]] std::uint32_t NewPCurve(const std::uint32_t old) const noexcept
    {
        return old < m_PCurves.size() ? m_PCurves[old] : PCurveRef::Invalid;
    }
    [[nodiscard]] std::uint32_t NewSurface(const std::uint32_t old) const noexcept
    {
        return old < m_Surfaces.size() ? m_Surfaces[old] : SurfaceRef::Invalid;
    }
    /// Slots of kind K in the source, live or dead.
    template <EntityKind K>
    [[nodiscard]] std::size_t OldCount() const noexcept
    {
        return m_Entities[static_cast<std::size_t>(K)].size();
    }

private:
    friend Result<Compacted> Compact(const Shape& shape);

    std::array<std::vector<std::uint32_t>, EntityKindCount> m_Entities;
    std::vector<std::uint32_t> m_Curves;
    std::vector<std::uint32_t> m_PCurves;
    std::vector<std::uint32_t> m_Surfaces;
};

/// What Compact() gives: the compacted shape, certified like any other,
/// and the map from the source's slots to its own.
class Compacted {
public:
    Stereon::Shape Shape;
    RemapTable Remap;
    /// The store of the shape that was compacted, which Translate() checks
    /// in checked builds.
    StoreId SourceStore = 0;

    /// The compacted shape's handle to an entity of the source; invalid for a
    /// handle that named no live entity of the source.
    template <EntityKind K>
    [[nodiscard]] Id<K> Translate(const Id<K> old) const noexcept
    {
        if (!old.IsValid() || Shape.IsEmpty()) {
            return {};
        }
        if constexpr (STEREON_CHECKED_HANDLES) {
            if (old.GetStore() != SourceStore) {
                return {};
            }
        }
        const std::uint32_t index = Remap.template NewIndex<K>(old.GetIndex());
        return index == Id<K>::Invalid ? Id<K>{} : Shape.GetSnapshot()->template MakeId<K>(index);
    }
};

/// The explicit renumbering of ADR-0004: the shape's live entities with
/// consecutive indices in their old order, dead slots and unreferenced
/// geometry dropped, frozen and certified again, with the remap table. The
/// source shape is left as it is. InvalidInput for an empty Shape.
[[nodiscard]] Result<Compacted> Compact(const Shape& shape);

/// Seals the builder into a Shape: runs the structural validator over every
/// table (ADR-0038's Freeze invariants) and, if it finds nothing, moves the
/// tables into a certified snapshot, keeping every index. Otherwise reports
/// InvalidInput with the findings. Either way the builder is consumed.
[[nodiscard]] Result<Shape> Freeze(ShapeBuilder&& builder);

/// The structural check of a frozen shape, run again by the independent
/// validator: what Freeze() certified, re-derived (ADR-0036).
[[nodiscard]] ValidityReport ValidateTopology(const Shape& shape);

/// The checks of a level and of every level below it (ADR-0018). Geometric
/// adds ADR-0001's containment invariant, sampled: each vertex's tolerance
/// sphere holds its edges' ends (G01), and each edge's tolerance tube holds
/// its 3D curve and its pcurves mapped through their surfaces (G02), and a
/// Degenerate edge's traces stay at its vertex (G02). It runs only on a
/// structurally valid shape. Solid has no rules yet and checks as Geometric.
/// An empty Shape is valid at every level.
[[nodiscard]] ValidityReport ValidateTopology(const Shape& shape, ValidityLevel level);

/// Publishes the builder's tables as a shape-producing operation does
/// (ADR-0006, ADR-0026): Freeze() and its structural certification, then the
/// OperationResult. Its History compares the result with the shape the
/// builder was opened on (all Generated for a new builder); its tolerance
/// changes are the builder's SetTolerance() calls, one per entity (the first
/// old value, the last new one, the reasons in order); its warnings are the
/// tolerances that grew and the cached edge curves the edit dropped. Fails,
/// publishing nothing, with InvalidInput when the tables are not
/// structurally valid or, in checked builds, when the shape the builder was
/// opened on no longer is, or when a live tolerance of the result is below
/// the Context's resolution; ToleranceExceeded when a live tolerance of the
/// result, logged or not, passes the Context's ceiling, or one grew under
/// ToleranceGrowth::Fail; and InvalidResult
/// when the History does not account for every entity (a kernel bug) or the
/// result fails the Geometric checks (ADR-0001's containment invariant),
/// which run when the Context's validity level asks for them and always in
/// checked builds. The builder is consumed either way.
[[nodiscard]] Result<OperationResult> Commit(ShapeBuilder&& builder, const Context& context);

/// Compact() as a shape-producing operation: the compacted shape with a
/// History in which every live entity continues Unchanged at its new index
/// (ADR-0011: the remap is composed into the History). It publishes under the
/// Context's policy as Commit() does: ToleranceExceeded for a tolerance above
/// the ceiling, InvalidInput for one below the resolution, and InvalidInput
/// when the Geometric level runs (the Context asks for it, or a checked
/// build) and the shape fails it; Checked says which level ran. In checked
/// builds the input's structure is validated again first, InvalidInput if
/// it fails.
[[nodiscard]] Result<OperationResult> CompactShape(const Shape& shape, const Context& context);

} // namespace Stereon

export namespace Stereon::Detail {

/// The independent structural validator (ADR-0018, ADR-0038): its own
/// bounded walks over the raw tables, checking every link before following
/// it, so it terminates and reads nothing out of range on any input.
[[nodiscard]] ValidityReport ValidateStructure(const TopologyTables& tables);

/// The Geometric checks of ADR-0018 on structurally valid tables: ADR-0001's
/// containment invariant, sampled (G01, G02).
[[nodiscard]] ValidityReport ValidateGeometry(const TopologyTables& tables);

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
    static void SetFirstEdge(ShapeBuilder& b, BuilderVertexId v, BuilderEdgeId first);
    static void SetNextAtStart(ShapeBuilder& b, BuilderEdgeId e, BuilderEdgeId next);
    static void SetNextAtEnd(ShapeBuilder& b, BuilderEdgeId e, BuilderEdgeId next);
    /// Marks the slot dead without unlinking it.
    static void KillCoedge(ShapeBuilder& b, BuilderCoedgeId c);
    /// Points a link at a slot past the end of its table.
    static void SetNextOutOfRange(ShapeBuilder& b, BuilderCoedgeId c);
    /// Empties a geometry pool's entry, which AddCurve, AddPCurve and
    /// AddSurface refuse to store.
    static void EmptyCurve(ShapeBuilder& b, CurveRef curve);
    static void EmptyPCurve(ShapeBuilder& b, PCurveRef pcurve);
    static void EmptySurface(ShapeBuilder& b, SurfaceRef surface);
    /// Sets an edge's cache fields as they are, consistent or not.
    static void SetCache(ShapeBuilder& b, BuilderEdgeId e, CurveRef cached, CacheCertificate status,
                         double bound);
};

} // namespace Stereon::Testing
