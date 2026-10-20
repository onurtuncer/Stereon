// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

module Stereon.Topo;

/// The structural validator of ADR-0018 and ADR-0038. It shares nothing
/// with the ranges it certifies: its own walks over the raw tables, each
/// link checked for bounds and liveness before it is followed, every walk
/// bounded by the size of the table it walks, so that it terminates and
/// reads nothing out of range on any input, corrupt or hostile. A finding
/// stops the walk that made it; the other walks go on, so one fault does
/// not hide another.

namespace Stereon::Detail {

namespace {

constexpr std::uint32_t None = 0xFFFFFFFFu;

[[nodiscard]] bool NeedsPrimary(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Intersection || kind == EdgeCurveKind::Trim ||
           kind == EdgeCurveKind::Boundary || kind == EdgeCurveKind::Degenerate;
}

[[nodiscard]] bool HasOwnCurve(const EdgeCurveKind kind) noexcept
{
    return kind == EdgeCurveKind::Exact || kind == EdgeCurveKind::Spatial;
}

/// The coedge counts ADR-0009 allows for each kind; a two-use Boundary edge
/// must also be a seam, which the caller checks.
[[nodiscard]] bool UseCountFits(const EdgeCurveKind kind, const std::size_t uses) noexcept
{
    switch (kind) {
    case EdgeCurveKind::Intersection:
        return uses >= 2;
    case EdgeCurveKind::Exact:
    case EdgeCurveKind::Degenerate:
        return uses >= 1;
    case EdgeCurveKind::Trim:
        return uses == 2;
    case EdgeCurveKind::Boundary:
        return uses == 1 || uses == 2;
    case EdgeCurveKind::Spatial:
        return true;
    }
    return false;
}

/// A certified cache's bound: a finite length, never negative.
[[nodiscard]] bool IsProvenBound(const double bound) noexcept
{
    return std::isfinite(bound) && bound >= 0.0; // STEREON-LINT-ALLOW(epsilon): a sign test
}

class StructureValidator {
public:
    explicit StructureValidator(const TopologyTables& tables)
        : m_T(tables)
        , m_VertexOk(tables.Vertices.GetSize(), false)
        , m_EdgeOk(tables.Edges.GetSize(), false)
        , m_CoedgeOk(tables.Coedges.GetSize(), false)
        , m_LoopOk(tables.Loops.GetSize(), false)
        , m_RingOk(tables.Edges.GetSize(), false)
        , m_InLoop(tables.Coedges.GetSize(), 0)
        , m_InRing(tables.Coedges.GetSize(), 0)
    {}

    [[nodiscard]] ValidityReport Run()
    {
        CheckLinks();
        CheckLoops();
        CheckRings();
        CheckLoopLists();
        CheckFaceLists();
        CheckShellLists();
        CheckGeometry();
        CheckEdgeCaches();
        CheckShellConnectivity();
        CheckIncidence();
        return std::move(m_Report);
    }

private:
    // --- Reporting -------------------------------------------------------

    void Report(const TopologyRule rule, std::vector<std::string> entities, std::string message)
    {
        m_Report.Add(
            {.Rule = rule, .Entities = std::move(entities), .Message = std::move(message)});
    }

    [[nodiscard]] std::string Name(const EntityKind kind, const std::uint32_t index) const
    {
        return DescribeEntity(m_T, kind, index);
    }

    // --- S01: every link names a live slot ----------------------------------

    /// True when `link` is live, or invalid and allowed to be; reports S01
    /// otherwise.
    template <EntityKind K>
    bool Link(const EntityKind owner, const std::uint32_t index, const Id<K> link,
              const bool required, const char* field)
    {
        if (!link.IsValid() ? !required : m_T.IsLive<K>(link.GetIndex())) {
            return true;
        }
        Report(TopologyRule::DanglingLink, {Name(owner, index)},
               std::string(field) +
                   (link.IsValid() ? " names a slot that is out of range or dead" : " is missing"));
        return false;
    }

    void CheckLinks()
    {
        for (std::uint32_t i = 0; i < m_T.Vertices.GetSize(); ++i) {
            m_VertexOk[i] =
                m_T.Vertices[i].Slot.Alive &&
                Link(EntityKind::Vertex, i, m_T.Vertices[i].FirstEdge, false, "FirstEdge");
        }
        for (std::uint32_t i = 0; i < m_T.Edges.GetSize(); ++i) {
            m_EdgeOk[i] = m_T.Edges[i].Slot.Alive && CheckEdgeLinks(i);
        }
        for (std::uint32_t i = 0; i < m_T.Coedges.GetSize(); ++i) {
            m_CoedgeOk[i] = m_T.Coedges[i].Slot.Alive && CheckCoedgeLinks(i);
        }
        for (std::uint32_t i = 0; i < m_T.Loops.GetSize(); ++i) {
            m_LoopOk[i] = m_T.Loops[i].Slot.Alive && CheckLoopLinks(i);
        }
        CheckParentLinks();
    }

    bool CheckEdgeLinks(const std::uint32_t i)
    {
        const EdgeRecord& e = m_T.Edges[i];
        constexpr EntityKind k = EntityKind::Edge;
        bool ok = Link(k, i, e.Start, true, "Start");
        ok = Link(k, i, e.End, true, "End") && ok;
        ok = Link(k, i, e.RadialAnchor, false, "RadialAnchor") && ok;
        ok = Link(k, i, e.NextAtStart, false, "NextAtStart") && ok;
        ok = Link(k, i, e.NextAtEnd, false, "NextAtEnd") && ok;
        return Link(k, i, e.Primary, false, "Primary") && ok;
    }

    bool CheckCoedgeLinks(const std::uint32_t i)
    {
        const CoedgeRecord& c = m_T.Coedges[i];
        constexpr EntityKind k = EntityKind::Coedge;
        bool ok = Link(k, i, c.Next, true, "Next");
        ok = Link(k, i, c.Prev, true, "Prev") && ok;
        ok = Link(k, i, c.Radial, true, "Radial") && ok;
        ok = Link(k, i, c.Edge, true, "Edge") && ok;
        return Link(k, i, c.Loop, true, "Loop") && ok;
    }

    bool CheckLoopLinks(const std::uint32_t i)
    {
        const LoopRecord& l = m_T.Loops[i];
        constexpr EntityKind k = EntityKind::Loop;
        bool ok = Link(k, i, l.Face, true, "Face");
        ok = Link(k, i, l.Anchor, false, "Anchor") && ok;
        return Link(k, i, l.NextInFace, false, "NextInFace") && ok;
    }

    void CheckParentLinks()
    {
        for (std::uint32_t i = 0; i < m_T.Faces.GetSize(); ++i) {
            const FaceRecord& f = m_T.Faces[i];
            if (f.Slot.Alive) {
                Link(EntityKind::Face, i, f.Shell, false, "Shell");
                Link(EntityKind::Face, i, f.FirstLoop, false, "FirstLoop");
                Link(EntityKind::Face, i, f.NextInShell, false, "NextInShell");
            }
        }
        for (std::uint32_t i = 0; i < m_T.Shells.GetSize(); ++i) {
            const ShellRecord& s = m_T.Shells[i];
            if (s.Slot.Alive) {
                Link(EntityKind::Shell, i, s.Solid, false, "Solid");
                Link(EntityKind::Shell, i, s.FirstFace, false, "FirstFace");
                Link(EntityKind::Shell, i, s.NextInSolid, false, "NextInSolid");
            }
        }
        for (std::uint32_t i = 0; i < m_T.Solids.GetSize(); ++i) {
            if (m_T.Solids[i].Slot.Alive) {
                Link(EntityKind::Solid, i, m_T.Solids[i].FirstShell, false, "FirstShell");
            }
        }
    }

    // --- S02, S03, S04, S14: loop cycles --------------------------------------

    void CheckLoops()
    {
        for (std::uint32_t l = 0; l < m_T.Loops.GetSize(); ++l) {
            if (m_LoopOk[l]) {
                WalkLoop(l);
            }
        }
        for (std::uint32_t c = 0; c < m_T.Coedges.GetSize(); ++c) {
            if (m_CoedgeOk[c]) {
                CheckReciprocity(c);
                CheckInLoop(c);
            }
        }
    }

    void WalkLoop(const std::uint32_t l)
    {
        const std::uint32_t anchor = m_T.Loops[l].Anchor.GetIndex();
        if (anchor == None) {
            Report(TopologyRule::EmptyLoop, {Name(EntityKind::Loop, l)}, "the loop has no coedge");
            return;
        }
        std::uint32_t c = anchor;
        for (std::size_t steps = 0; steps <= m_T.Coedges.GetSize(); ++steps) {
            if (!m_CoedgeOk[c]) {
                return; // reported by S01
            }
            if (!VisitInLoop(l, c)) {
                return;
            }
            c = m_T.Coedges[c].Next.GetIndex();
            if (c == anchor) {
                return;
            }
        }
        Report(TopologyRule::LoopOpen, {Name(EntityKind::Loop, l)},
               "the walk from the anchor does not return to it");
    }

    /// One step of a loop walk at coedge c; false when the walk must stop.
    bool VisitInLoop(const std::uint32_t l, const std::uint32_t c)
    {
        if (m_T.Coedges[c].Loop.GetIndex() != l) {
            Report(TopologyRule::LoopOpen, {Name(EntityKind::Loop, l), Name(EntityKind::Coedge, c)},
                   "the walk reaches a coedge of another loop");
            return false;
        }
        if (m_InLoop[c] != 0) {
            Report(TopologyRule::LoopOpen, {Name(EntityKind::Loop, l), Name(EntityKind::Coedge, c)},
                   "the walk revisits a coedge before returning to the anchor");
            return false;
        }
        m_InLoop[c] = 1;
        return true;
    }

    void CheckReciprocity(const std::uint32_t c)
    {
        const CoedgeRecord& r = m_T.Coedges[c];
        const std::uint32_t next = r.Next.GetIndex();
        const std::uint32_t prev = r.Prev.GetIndex();
        if (m_T.Coedges[next].Prev.GetIndex() != c) {
            Report(TopologyRule::LoopLinksBroken,
                   {Name(EntityKind::Coedge, c), Name(EntityKind::Coedge, next)},
                   "Prev of Next is not the coedge");
        }
        if (m_T.Coedges[prev].Next.GetIndex() != c) {
            Report(TopologyRule::LoopLinksBroken,
                   {Name(EntityKind::Coedge, c), Name(EntityKind::Coedge, prev)},
                   "Next of Prev is not the coedge");
        }
    }

    void CheckInLoop(const std::uint32_t c)
    {
        if (m_InLoop[c] == 0 && m_LoopOk[m_T.Coedges[c].Loop.GetIndex()]) {
            Report(TopologyRule::CoedgeOutsideLoop, {Name(EntityKind::Coedge, c)},
                   "the walk of the coedge's loop does not reach it");
        }
    }

    // --- S05, S06, S07, S08: radial rings ---------------------------------------

    void CheckRings()
    {
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (m_EdgeOk[e]) {
                const std::size_t uses = WalkRing(e);
                if (m_RingOk[e]) {
                    CheckUseCount(e, uses);
                    CheckPrimary(e, uses);
                }
            }
        }
        for (std::uint32_t c = 0; c < m_T.Coedges.GetSize(); ++c) {
            if (m_CoedgeOk[c] && m_InRing[c] == 0 && m_RingOk[m_T.Coedges[c].Edge.GetIndex()]) {
                Report(TopologyRule::CoedgeOutsideRing, {Name(EntityKind::Coedge, c)},
                       "the radial walk of the coedge's edge does not reach it");
            }
        }
    }

    /// The uses of edge e along its ring; sets m_RingOk[e] when the ring
    /// closes over uses of e only.
    std::size_t WalkRing(const std::uint32_t e)
    {
        const std::uint32_t anchor = m_T.Edges[e].RadialAnchor.GetIndex();
        std::size_t uses = 0;
        if (anchor == None) {
            m_RingOk[e] = true;
            return 0;
        }
        for (std::uint32_t c = anchor; uses <= m_T.Coedges.GetSize(); ++uses) {
            if (!m_CoedgeOk[c] || !VisitInRing(e, c)) {
                return uses;
            }
            c = m_T.Coedges[c].Radial.GetIndex();
            if (c == anchor) {
                m_RingOk[e] = true;
                return uses + 1;
            }
        }
        Report(TopologyRule::RadialRingOpen, {Name(EntityKind::Edge, e)},
               "the radial walk does not return to its anchor");
        return uses;
    }

    bool VisitInRing(const std::uint32_t e, const std::uint32_t c)
    {
        if (m_T.Coedges[c].Edge.GetIndex() != e) {
            Report(TopologyRule::RadialRingOpen,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, c)},
                   "the radial walk reaches a use of another edge");
            return false;
        }
        if (m_InRing[c] != 0) {
            Report(TopologyRule::RadialRingOpen,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, c)},
                   "the radial walk revisits a coedge before returning to the anchor");
            return false;
        }
        m_InRing[c] = 1;
        return true;
    }

    /// The face of a coedge whose links are known good, or None.
    [[nodiscard]] std::uint32_t FaceOf(const std::uint32_t c) const
    {
        const std::uint32_t l = m_T.Coedges[c].Loop.GetIndex();
        return m_LoopOk[l] ? m_T.Loops[l].Face.GetIndex() : None;
    }

    void CheckUseCount(const std::uint32_t e, const std::size_t uses)
    {
        const EdgeCurveKind kind = m_T.Edges[e].Kind;
        if (!UseCountFits(kind, uses)) {
            Report(TopologyRule::EdgeUseCount, {Name(EntityKind::Edge, e)},
                   std::string("a ") + std::string(ToString(kind)) + " edge with " +
                       std::to_string(uses) + " coedges");
            return;
        }
        if (kind == EdgeCurveKind::Boundary && uses == 2) {
            const std::uint32_t first = m_T.Edges[e].RadialAnchor.GetIndex();
            const std::uint32_t second = m_T.Coedges[first].Radial.GetIndex();
            if (FaceOf(first) != FaceOf(second)) {
                Report(TopologyRule::EdgeUseCount, {Name(EntityKind::Edge, e)},
                       "a Boundary edge with two coedges on different faces is not a seam");
            }
        }
    }

    void CheckPrimary(const std::uint32_t e, const std::size_t uses)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        const std::uint32_t primary = edge.Primary.GetIndex();
        if (primary == None) {
            if (NeedsPrimary(edge.Kind) && uses > 0) {
                Report(TopologyRule::PrimaryCoedge, {Name(EntityKind::Edge, e)},
                       std::string("a ") + std::string(ToString(edge.Kind)) +
                           " edge without its primary coedge");
            }
            return;
        }
        if (m_T.Coedges[primary].Edge.GetIndex() != e) {
            Report(TopologyRule::PrimaryCoedge,
                   {Name(EntityKind::Edge, e), Name(EntityKind::Coedge, primary)},
                   "the primary coedge is not a use of the edge");
        }
    }

    // --- S09, S10, S11: parent lists -----------------------------------------------

    /// Checks that every live child of kind C with a parent is in that
    /// parent's list exactly once, and that every list names children of
    /// its own parent and ends at the recorded last member.
    template <EntityKind P, EntityKind C, typename Members, typename Next, typename ParentOf>
    void CheckMembership(const TopologyRule rule, Members members, Next next, ParentOf parentOf)
    {
        const auto& parents = m_T.Table<P>();
        const auto& children = m_T.Table<C>();
        std::vector<std::uint8_t> seen(children.GetSize(), 0);
        for (std::uint32_t p = 0; p < parents.GetSize(); ++p) {
            if (parents[p].Slot.Alive) {
                const auto [first, last] = members(parents[p]);
                WalkList<P, C>(rule, p, first, last, seen, next, parentOf);
            }
        }
        for (std::uint32_t c = 0; c < children.GetSize(); ++c) {
            const std::uint32_t parent = parentOf(children[c]);
            if (children[c].Slot.Alive && seen[c] == 0 && m_T.IsLive<P>(parent)) {
                Report(rule, {Name(C, c), Name(P, parent)}, "missing from its parent's list");
            }
        }
    }

    template <EntityKind P, EntityKind C, typename Next, typename ParentOf>
    void WalkList(const TopologyRule rule, const std::uint32_t p, std::uint32_t child,
                  const std::uint32_t last, std::vector<std::uint8_t>& seen, Next next,
                  ParentOf parentOf)
    {
        const auto& children = m_T.Table<C>();
        std::uint32_t previous = None;
        while (m_T.IsLive<C>(child)) {
            if (parentOf(children[child]) != p || seen[child] != 0) {
                Report(rule, {Name(P, p), Name(C, child)},
                       seen[child] != 0 ? "the list names a member twice"
                                        : "the list names a member of another parent");
                return;
            }
            seen[child] = 1;
            previous = child;
            child = next(children[child]);
        }
        if (child == None && previous != last) {
            Report(rule, {Name(P, p)}, "the recorded last member is not the end of the list");
        }
    }

    void CheckLoopLists()
    {
        CheckMembership<EntityKind::Face, EntityKind::Loop>(
            TopologyRule::LoopOutsideFace,
            [](const FaceRecord& f) {
                return std::pair{f.FirstLoop.GetIndex(), f.LastLoop.GetIndex()};
            },
            [](const LoopRecord& l) { return l.NextInFace.GetIndex(); },
            [](const LoopRecord& l) { return l.Face.GetIndex(); });
    }

    void CheckFaceLists()
    {
        CheckMembership<EntityKind::Shell, EntityKind::Face>(
            TopologyRule::FaceOutsideShell,
            [](const ShellRecord& s) {
                return std::pair{s.FirstFace.GetIndex(), s.LastFace.GetIndex()};
            },
            [](const FaceRecord& f) { return f.NextInShell.GetIndex(); },
            [](const FaceRecord& f) { return f.Shell.GetIndex(); });
    }

    void CheckShellLists()
    {
        CheckMembership<EntityKind::Solid, EntityKind::Shell>(
            TopologyRule::ShellOutsideSolid,
            [](const SolidRecord& s) {
                return std::pair{s.FirstShell.GetIndex(), s.LastShell.GetIndex()};
            },
            [](const ShellRecord& s) { return s.NextInSolid.GetIndex(); },
            [](const ShellRecord& s) { return s.Solid.GetIndex(); });
    }

    // --- S12, S13: geometry associations ------------------------------------------

    /// An entry of a geometry pool is a support only when it holds geometry:
    /// an empty wrapper in the pool is as missing as an index past its end.
    void CheckGeometry()
    {
        for (std::uint32_t f = 0; f < m_T.Faces.GetSize(); ++f) {
            const SurfaceRef s = m_T.Faces[f].Surface;
            if (!m_T.Faces[f].Slot.Alive) {
                continue;
            }
            if (!(s.IsValid() && s.GetIndex() < m_T.Surfaces.GetSize())) {
                Report(TopologyRule::MissingGeometry, {Name(EntityKind::Face, f)},
                       "the face has no surface in the pool");
            } else if (m_T.Surfaces[s.GetIndex()].IsEmpty()) {
                Report(TopologyRule::MissingGeometry, {Name(EntityKind::Face, f)},
                       "the face's surface in the pool is empty");
            }
        }
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (m_EdgeOk[e]) {
                CheckEdgeGeometry(e);
            }
        }
        for (std::uint32_t c = 0; c < m_T.Coedges.GetSize(); ++c) {
            if (m_CoedgeOk[c]) {
                CheckPCurve(c);
            }
        }
    }

    void CheckEdgeGeometry(const std::uint32_t e)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        const bool hasCurve = edge.Curve.IsValid() && edge.Curve.GetIndex() < m_T.Curves.GetSize();
        if (HasOwnCurve(edge.Kind) != hasCurve || (edge.Curve.IsValid() && !hasCurve)) {
            Report(TopologyRule::MissingGeometry, {Name(EntityKind::Edge, e)},
                   HasOwnCurve(edge.Kind) ? "the edge has no 3D curve in the pool"
                                          : "the edge kind has no 3D curve of its own");
        } else if (hasCurve && m_T.Curves[edge.Curve.GetIndex()].IsEmpty()) {
            Report(TopologyRule::MissingGeometry, {Name(EntityKind::Edge, e)},
                   "the edge's 3D curve in the pool is empty");
        }
        if (edge.Range.IsEmpty()) {
            Report(TopologyRule::MissingGeometry, {Name(EntityKind::Edge, e)},
                   "the edge has an empty parameter range");
        }
        if (edge.Kind == EdgeCurveKind::Degenerate && edge.Start != edge.End) {
            Report(TopologyRule::EdgeEndpoints, {Name(EntityKind::Edge, e)},
                   "a degenerate edge with two different vertices");
        }
    }

    void CheckPCurve(const std::uint32_t c)
    {
        const PCurveRef p = m_T.Coedges[c].PCurve;
        const EdgeCurveKind kind = m_T.Edges[m_T.Coedges[c].Edge.GetIndex()].Kind;
        const bool inPool = p.IsValid() && p.GetIndex() < m_T.PCurves.GetSize();
        if ((p.IsValid() && !inPool) || (!p.IsValid() && kind != EdgeCurveKind::Spatial)) {
            Report(TopologyRule::MissingGeometry, {Name(EntityKind::Coedge, c)},
                   "the coedge has no pcurve in the pool");
        } else if (inPool && m_T.PCurves[p.GetIndex()].IsEmpty()) {
            Report(TopologyRule::MissingGeometry, {Name(EntityKind::Coedge, c)},
                   "the coedge's pcurve in the pool is empty");
        }
    }

    // --- S17: cached edge curves and their certificates ---------------------------

    void CheckEdgeCaches()
    {
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (m_EdgeOk[e]) {
                CheckEdgeCache(e);
            }
        }
    }

    /// The edge's cache fields agree with one another and with its kind;
    /// whether the cache is close to S(c(t)) is a Geometric matter.
    void CheckEdgeCache(const std::uint32_t e)
    {
        const EdgeRecord& edge = m_T.Edges[e];
        const bool cached = edge.Cached.IsValid();
        const bool certified = edge.CacheStatus == CacheCertificate::Certified;
        const char* problem = nullptr;
        if (cached && !HasCachedCurve(edge.Kind)) {
            problem = "the edge kind keeps no cached curve";
        } else if (cached && (edge.Cached.GetIndex() >= m_T.Curves.GetSize() ||
                              m_T.Curves[edge.Cached.GetIndex()].IsEmpty())) {
            problem = "the cached curve names no non-empty pool entry";
        } else if (cached != (edge.CacheStatus != CacheCertificate::None)) {
            problem = "the cache status does not say whether there is a cache";
        } else if (certified ? !IsProvenBound(edge.CacheBound)
                             : edge.CacheBound != std::numeric_limits<double>::infinity()) {
            problem = "the cache bound does not fit its status";
        }
        if (problem != nullptr) {
            Report(TopologyRule::EdgeCache, {Name(EntityKind::Edge, e)}, problem);
        }
    }

    // --- S15: shells are edge-connected --------------------------------------------

    [[nodiscard]] std::uint32_t Find(std::vector<std::uint32_t>& parent, std::uint32_t f) const
    {
        while (parent[f] != f) {
            parent[f] = parent[parent[f]];
            f = parent[f];
        }
        return f;
    }

    void CheckShellConnectivity()
    {
        std::vector<std::uint32_t> parent(m_T.Faces.GetSize());
        for (std::uint32_t f = 0; f < parent.size(); ++f) {
            parent[f] = f;
        }
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (m_EdgeOk[e] && m_RingOk[e]) {
                JoinFacesOf(e, parent);
            }
        }
        for (std::uint32_t s = 0; s < m_T.Shells.GetSize(); ++s) {
            if (m_T.Shells[s].Slot.Alive) {
                CheckShellComponent(s, parent);
            }
        }
    }

    /// Joins the faces around edge e that lie in one shell, for every shell
    /// among them: each face to the first face of its shell on the ring. A
    /// non-manifold edge may be shared by several shells, and which of them
    /// holds the anchor must not matter.
    void JoinFacesOf(const std::uint32_t e, std::vector<std::uint32_t>& parent) const
    {
        const std::uint32_t anchor = m_T.Edges[e].RadialAnchor.GetIndex();
        if (anchor == None) {
            return;
        }
        std::vector<std::uint32_t> firsts; // the first face met in each shell
        std::uint32_t c = anchor;
        do {
            const std::uint32_t face = FaceOf(c);
            if (face != None) {
                JoinToShellFirst(face, firsts, parent);
            }
            c = m_T.Coedges[c].Radial.GetIndex();
        } while (c != anchor);
    }

    void JoinToShellFirst(const std::uint32_t face, std::vector<std::uint32_t>& firsts,
                          std::vector<std::uint32_t>& parent) const
    {
        for (const std::uint32_t first : firsts) {
            if (m_T.Faces[first].Shell == m_T.Faces[face].Shell) {
                parent[Find(parent, face)] = Find(parent, first);
                return;
            }
        }
        firsts.push_back(face);
    }

    void CheckShellComponent(const std::uint32_t s, std::vector<std::uint32_t>& parent)
    {
        const std::uint32_t first = m_T.Shells[s].FirstFace.GetIndex();
        if (!m_T.IsLive<EntityKind::Face>(first)) {
            return;
        }
        const std::uint32_t root = Find(parent, first);
        std::uint32_t f = m_T.Faces[first].NextInShell.GetIndex();
        for (std::size_t steps = 0; m_T.IsLive<EntityKind::Face>(f) && steps < m_T.Faces.GetSize();
             ++steps) {
            if (Find(parent, f) != root) {
                Report(TopologyRule::ShellDisconnected,
                       {Name(EntityKind::Shell, s), Name(EntityKind::Face, f)},
                       "the face shares no chain of edges with the shell's first face");
                return;
            }
            f = m_T.Faces[f].NextInShell.GetIndex();
        }
    }

    // --- S16: the vertex-incidence index ---------------------------------------

    /// Walks every live vertex's list, then checks that every live edge was
    /// reached from each of its ends once (a closed edge once in all).
    void CheckIncidence()
    {
        std::vector<std::uint8_t> fromStart(m_T.Edges.GetSize(), 0);
        std::vector<std::uint8_t> fromEnd(m_T.Edges.GetSize(), 0);
        for (std::uint32_t v = 0; v < m_T.Vertices.GetSize(); ++v) {
            if (m_VertexOk[v]) {
                WalkIncidence(v, fromStart, fromEnd);
            }
        }
        for (std::uint32_t e = 0; e < m_T.Edges.GetSize(); ++e) {
            if (!m_EdgeOk[e]) {
                continue;
            }
            const EdgeRecord& edge = m_T.Edges[e];
            const bool closed = edge.Start == edge.End;
            if (fromStart[e] != 1 && m_VertexOk[edge.Start.GetIndex()]) {
                Report(TopologyRule::VertexIncidence,
                       {Name(EntityKind::Edge, e), Name(EntityKind::Vertex, edge.Start.GetIndex())},
                       "the edge is not in its start vertex's list exactly once");
            }
            if (!closed && fromEnd[e] != 1 && m_VertexOk[edge.End.GetIndex()]) {
                Report(TopologyRule::VertexIncidence,
                       {Name(EntityKind::Edge, e), Name(EntityKind::Vertex, edge.End.GetIndex())},
                       "the edge is not in its end vertex's list exactly once");
            }
        }
    }

    /// One vertex's list: every edge live, at the vertex, after the one
    /// before it in index, and reached once.
    void WalkIncidence(const std::uint32_t v, std::vector<std::uint8_t>& fromStart,
                       std::vector<std::uint8_t>& fromEnd)
    {
        std::uint32_t previous = None;
        std::uint32_t e = m_T.Vertices[v].FirstEdge.GetIndex();
        for (std::size_t steps = 0; e != None && steps <= m_T.Edges.GetSize(); ++steps) {
            if (!m_EdgeOk[e]) {
                return; // reported by S01
            }
            const EdgeRecord& edge = m_T.Edges[e];
            const bool atStart = edge.Start.GetIndex() == v;
            if (!atStart && edge.End.GetIndex() != v) {
                Report(TopologyRule::VertexIncidence,
                       {Name(EntityKind::Vertex, v), Name(EntityKind::Edge, e)},
                       "the list names an edge that does not start or end at the vertex");
                return;
            }
            if (previous != None && e <= previous) {
                Report(TopologyRule::VertexIncidence,
                       {Name(EntityKind::Vertex, v), Name(EntityKind::Edge, e)},
                       "the list is not in ascending index order");
                return;
            }
            std::uint8_t& reached = atStart ? fromStart[e] : fromEnd[e];
            if (reached != 0) {
                Report(TopologyRule::VertexIncidence,
                       {Name(EntityKind::Vertex, v), Name(EntityKind::Edge, e)},
                       "the list names an edge twice");
                return;
            }
            reached = 1;
            previous = e;
            e = (atStart ? edge.NextAtStart : edge.NextAtEnd).GetIndex();
        }
        if (e != None) {
            Report(TopologyRule::VertexIncidence, {Name(EntityKind::Vertex, v)},
                   "the list does not end within the edge count");
        }
    }

    const TopologyTables& m_T;
    ValidityReport m_Report;
    std::vector<bool> m_VertexOk;
    std::vector<bool> m_EdgeOk;
    std::vector<bool> m_CoedgeOk;
    std::vector<bool> m_LoopOk;
    std::vector<bool> m_RingOk;
    std::vector<std::uint8_t> m_InLoop;
    std::vector<std::uint8_t> m_InRing;
};

} // namespace

ValidityReport ValidateStructure(const TopologyTables& tables)
{
    return StructureValidator(tables).Run();
}

} // namespace Stereon::Detail
