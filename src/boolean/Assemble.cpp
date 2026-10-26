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
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

module Stereon.Boolean;

import :Assemble;

namespace Stereon {

namespace {

constexpr std::uint8_t BodyA = 0;
constexpr std::uint8_t BodyB = 1;

/// An entity of an imprinted body, by body and slot.
using Key = std::pair<std::uint8_t, std::uint32_t>;

// --- Where the imprinted entities came from -------------------------------------

/// Where an entity of an imprinted body came from in its operand: its own
/// slot when the imprint kept it (Open() keeps every index), the slot it was
/// split from, through any chain of splits, when the imprint made it by a
/// split, and nothing when the imprint made it from scratch.
class Origins {
public:
    explicit Origins(const OperationResult& imprinted)
    {
        for (const HistoryEntry& e : imprinted.History.GetEntries()) {
            if (e.Relation == HistoryRelation::Unchanged ||
                e.Relation == HistoryRelation::Modified) {
                m_Kept.insert({static_cast<int>(e.Kind), e.Target});
            }
        }
        for (const TopologyEvent& e : imprinted.TopologyEvents) {
            if (e.Event == TopologyEventKind::Split) {
                m_Parent[{static_cast<int>(e.Kind), e.Target}] = e.Source;
            }
        }
    }

    [[nodiscard]] std::optional<std::uint32_t> Of(const EntityKind kind, std::uint32_t slot) const
    {
        // A chain of splits is at most as long as the events that made it.
        for (std::size_t step = 0; step <= m_Parent.size(); ++step) {
            if (m_Kept.contains({static_cast<int>(kind), slot})) {
                return slot;
            }
            const auto parent = m_Parent.find({static_cast<int>(kind), slot});
            if (parent == m_Parent.end()) {
                return std::nullopt;
            }
            slot = parent->second;
        }
        return std::nullopt;
    }

private:
    std::set<std::pair<int, std::uint32_t>> m_Kept;
    std::map<std::pair<int, std::uint32_t>, std::uint32_t> m_Parent;
};

// --- What the result records of its sources ------------------------------------

/// An entity of an operand: the input's ordinal and its slot.
struct Source {
    std::uint8_t Input = 0;
    std::uint32_t Index = NoEntity;

    auto operator<=>(const Source&) const = default;
};

/// The sources of one entity of the result: the operand entities it
/// continues, and for an intersection edge the faces of each operand that
/// met there.
struct Provenance {
    std::set<Source> Sources;
    std::set<Source> Generators;
};

/// The provenance of every entity of the result, by kind and result slot.
using Ledger = std::map<std::pair<int, std::uint32_t>, Provenance>;

// --- The assembly in progress ------------------------------------------------------

/// A piece of the intersection, an edge in each body: the body whose edge
/// the result takes, and both edges.
struct SharedPair {
    std::uint8_t Representative = BodyA;
    std::array<EdgeId, 2> Edges;
};

struct Work {
    const ImprintedBodies& Bodies;
    const Context& Ctx;
    std::array<Origins, 2> Origin;
    std::array<bool, 2> Reversed{};
    std::array<std::set<std::uint32_t>, 2> Kept;
    ShapeBuilder Builder;
    /// B's vertices at the ends of the intersection's pieces, each the
    /// vertex of A it is identified with.
    std::map<Key, VertexId> Alias;
    /// A's vertices with the B vertices identified with them.
    std::map<std::uint32_t, std::vector<VertexId>> Aliased;
    std::map<Key, std::size_t> PairOf;
    std::vector<SharedPair> Pairs;
    std::map<Key, BuilderVertexId> Vertices;
    std::map<Key, BuilderEdgeId> Edges;
    /// The result's edges that keep a cached curve, to build once their
    /// uses exist.
    std::vector<BuilderEdgeId> Cached;
    /// The result's faces, the body and face each copies, and the result's
    /// edges each uses.
    std::vector<BuilderFaceId> Faces;
    std::vector<std::pair<std::uint8_t, FaceId>> Copied;
    std::vector<std::vector<std::uint32_t>> FaceEdges;
    Ledger Records;

    [[nodiscard]] const Shape& ShapeOf(const std::uint8_t body) const
    {
        return body == BodyA ? Bodies.A.Shape : Bodies.B.Shape;
    }
};

void Record(Work& w, const EntityKind kind, const std::uint32_t target, const std::uint8_t body,
            const std::optional<std::uint32_t> origin)
{
    Provenance& p = w.Records[{static_cast<int>(kind), target}];
    if (origin) {
        p.Sources.insert({body, *origin});
    }
}

// --- The pieces of the intersection ---------------------------------------------

/// Whether a kept face of the body uses the edge.
[[nodiscard]] Result<bool> KeptUses(const Work& w, const std::uint8_t body, const EdgeId edge)
{
    STEREON_TRY(const std::vector<FaceId> faces, FacesOf(w.ShapeOf(body), edge));
    return std::ranges::any_of(faces,
                               [&](const FaceId f) { return w.Kept[body].contains(f.GetIndex()); });
}

/// Identifies B's vertex with A's: within the pair tolerance, the larger
/// of the two (decision D2), of A's, whose position the result keeps.
[[nodiscard]] Result<void> Identify(Work& w, const VertexId a, const VertexId b)
{
    const Shape& sa = w.Bodies.A.Shape;
    const Shape& sb = w.Bodies.B.Shape;
    const double tolerance = std::max(ToleranceOf(sa, a), ToleranceOf(sb, b));
    if (Distance(PositionOf(sa, a), PositionOf(sb, b)) > tolerance) {
        return Fail(ErrorCode::IntersectionFailed,
                    "Boolean: a piece of the intersection ends at " + std::string(Describe(sa, a)) +
                        " in one body and at " + std::string(Describe(sb, b)) +
                        " in the other, further apart than their tolerances");
    }
    if (w.Alias.try_emplace({BodyB, b.GetIndex()}, a).second) {
        w.Aliased[a.GetIndex()].push_back(b);
    }
    return {};
}

/// The pieces of the intersection: each edge pair, the body whose edge the
/// result takes (A's where a kept face of A uses it), and the ends of B's
/// identified with A's. Both edges carry the branch's curve over one range,
/// so they run the same way.
[[nodiscard]] Result<void> PairPieces(Work& w)
{
    for (const SharedEdge& shared : w.Bodies.Shared) {
        STEREON_TRY(const bool usedByA, KeptUses(w, BodyA, shared.A));
        SharedPair pair{usedByA ? BodyA : BodyB, {shared.A, shared.B}};
        w.PairOf[{BodyA, shared.A.GetIndex()}] = w.Pairs.size();
        w.PairOf[{BodyB, shared.B.GetIndex()}] = w.Pairs.size();
        w.Pairs.push_back(pair);
        const Shape& sa = w.Bodies.A.Shape;
        const Shape& sb = w.Bodies.B.Shape;
        STEREON_TRY_VOID(Identify(w, StartOf(sa, shared.A), StartOf(sb, shared.B)));
        STEREON_TRY_VOID(Identify(w, EndOf(sa, shared.A), EndOf(sb, shared.B)));
    }
    return {};
}

// --- Vertices and edges ------------------------------------------------------------

/// The result's vertex for a body's vertex: B's identified with A's is A's,
/// made once, at A's position with the largest tolerance of the vertices
/// identified, the growth logged (ADR-0001, rule 4).
[[nodiscard]] Result<BuilderVertexId> VertexFor(Work& w, std::uint8_t body, VertexId vertex)
{
    if (const auto alias = w.Alias.find({body, vertex.GetIndex()}); alias != w.Alias.end()) {
        body = BodyA;
        vertex = alias->second;
    }
    const Key key{body, vertex.GetIndex()};
    if (const auto made = w.Vertices.find(key); made != w.Vertices.end()) {
        return made->second;
    }
    const Shape& shape = w.ShapeOf(body);
    const double own = ToleranceOf(shape, vertex);
    STEREON_TRY(const BuilderVertexId made, w.Builder.AddVertex(PositionOf(shape, vertex), own));
    w.Vertices[key] = made;
    Record(w, EntityKind::Vertex, made.GetIndex(), body,
           w.Origin[body].Of(EntityKind::Vertex, vertex.GetIndex()));
    double tolerance = own;
    for (const VertexId b :
         body == BodyA ? w.Aliased[vertex.GetIndex()] : std::vector<VertexId>{}) {
        tolerance = std::max(tolerance, ToleranceOf(w.Bodies.B.Shape, b));
        Record(w, EntityKind::Vertex, made.GetIndex(), BodyB,
               w.Origin[BodyB].Of(EntityKind::Vertex, b.GetIndex()));
    }
    if (tolerance > own) {
        STEREON_TRY_VOID(w.Builder.SetTolerance(
            made, tolerance, "Boolean: identified with the other body's vertex (ADR-0001)"));
    }
    return made;
}

/// The faces of each operand an intersection edge's faces came from: the
/// generators of the edge (ADR-0030, section 6).
[[nodiscard]] Result<void> RecordGenerators(Work& w, const std::uint32_t target,
                                            const SharedPair& pair)
{
    Provenance& p = w.Records[{static_cast<int>(EntityKind::Edge), target}];
    for (const std::uint8_t body : {BodyA, BodyB}) {
        STEREON_TRY(const std::vector<FaceId> faces, FacesOf(w.ShapeOf(body), pair.Edges[body]));
        for (const FaceId f : faces) {
            if (const auto origin = w.Origin[body].Of(EntityKind::Face, f.GetIndex())) {
                p.Generators.insert({body, *origin});
            }
        }
    }
    return {};
}

/// The result's edge for an edge of the representative body: its kind,
/// range and curve, between the result's vertices for its ends, at the
/// larger tolerance of the pair's two edges when it is a piece of the
/// intersection.
[[nodiscard]] Result<BuilderEdgeId> MakeEdge(Work& w, const std::uint8_t body, const EdgeId edge,
                                             const std::optional<std::size_t> pair)
{
    const Shape& shape = w.ShapeOf(body);
    EdgeSpec spec;
    spec.Kind = KindOf(shape, edge);
    STEREON_TRY(spec.Start, VertexFor(w, body, StartOf(shape, edge)));
    STEREON_TRY(spec.End, VertexFor(w, body, EndOf(shape, edge)));
    spec.Range = RangeOf(shape, edge);
    spec.Tolerance = ToleranceOf(shape, edge);
    if (const std::optional<Curve> curve = CurveOf(shape, edge)) {
        STEREON_TRY(spec.Curve, w.Builder.AddCurve(*curve));
    }
    STEREON_TRY(const BuilderEdgeId made, w.Builder.AddEdge(spec));
    if (spec.Kind == EdgeCurveKind::Trim || spec.Kind == EdgeCurveKind::Intersection) {
        w.Cached.push_back(made);
    }
    if (!pair) {
        return made;
    }
    const SharedPair& shared = w.Pairs[*pair];
    const std::uint8_t other = body == BodyA ? BodyB : BodyA;
    const double tolerance =
        std::max(spec.Tolerance, ToleranceOf(w.ShapeOf(other), shared.Edges[other]));
    if (tolerance > spec.Tolerance) {
        STEREON_TRY_VOID(w.Builder.SetTolerance(
            made, tolerance, "Boolean: one piece of the intersection with the other body's"));
    }
    STEREON_TRY_VOID(RecordGenerators(w, made.GetIndex(), shared));
    Record(w, EntityKind::Edge, made.GetIndex(), other,
           w.Origin[other].Of(EntityKind::Edge, shared.Edges[other].GetIndex()));
    return made;
}

/// The result's edge for a body's edge: a piece of the intersection is the
/// representative's, made once for both bodies.
[[nodiscard]] Result<BuilderEdgeId> EdgeFor(Work& w, std::uint8_t body, EdgeId edge)
{
    std::optional<std::size_t> pair;
    if (const auto shared = w.PairOf.find({body, edge.GetIndex()}); shared != w.PairOf.end()) {
        pair = shared->second;
        body = w.Pairs[*pair].Representative;
        edge = w.Pairs[*pair].Edges[body];
    }
    const Key key{body, edge.GetIndex()};
    if (const auto made = w.Edges.find(key); made != w.Edges.end()) {
        return made->second;
    }
    STEREON_TRY(const BuilderEdgeId made, MakeEdge(w, body, edge, pair));
    w.Edges[key] = made;
    Record(w, EntityKind::Edge, made.GetIndex(), body,
           w.Origin[body].Of(EntityKind::Edge, edge.GetIndex()));
    return made;
}

// --- Faces ---------------------------------------------------------------------------

[[nodiscard]] EdgeUseSense Flipped(const EdgeUseSense sense)
{
    return sense == EdgeUseSense::Same ? EdgeUseSense::Opposite : EdgeUseSense::Same;
}

/// Copies a loop of a body's face into the result's face: in order, or for
/// a reversed body against it, each use turned round.
[[nodiscard]] Result<void> CopyLoop(Work& w, const std::uint8_t body, const LoopId loop,
                                    const BuilderFaceId face, std::vector<std::uint32_t>& edges)
{
    const Shape& shape = w.ShapeOf(body);
    STEREON_TRY(const BuilderLoopId made, w.Builder.AddLoop(face));
    Record(w, EntityKind::Loop, made.GetIndex(), body,
           w.Origin[body].Of(EntityKind::Loop, loop.GetIndex()));
    STEREON_TRY(const CoedgeRange range, CoedgesOf(shape, loop));
    std::vector<CoedgeId> coedges(range.begin(), range.end());
    if (w.Reversed[body]) {
        std::ranges::reverse(coedges);
    }
    for (const CoedgeId coedge : coedges) {
        const std::optional<Curve2d> pcurve = PCurveOf(shape, coedge);
        if (!pcurve) {
            return Fail(ErrorCode::NotSupported,
                        "Boolean: a use of " + std::string(Describe(shape, EdgeOf(shape, coedge))) +
                            " has no pcurve");
        }
        STEREON_TRY(const BuilderEdgeId edge, EdgeFor(w, body, EdgeOf(shape, coedge)));
        STEREON_TRY(const PCurveRef curve, w.Builder.AddPCurve(*pcurve));
        const EdgeUseSense sense = StoredEdgeUseSense(shape, coedge);
        STEREON_TRY(
            const BuilderCoedgeId use,
            w.Builder.AppendCoedge(made, edge, w.Reversed[body] ? Flipped(sense) : sense, curve));
        Record(w, EntityKind::Coedge, use.GetIndex(), body,
               w.Origin[body].Of(EntityKind::Coedge, coedge.GetIndex()));
        edges.push_back(edge.GetIndex());
    }
    return {};
}

/// Copies a kept face of a body into the result, its orientation turned
/// round for a reversed body.
[[nodiscard]] Result<void> CopyFace(Work& w, const std::uint8_t body, const FaceId face)
{
    const Shape& shape = w.ShapeOf(body);
    STEREON_TRY(const SurfaceRef surface, w.Builder.AddSurface(SurfaceOf(shape, face)));
    STEREON_TRY(const BuilderFaceId made,
                w.Builder.AddFace(surface, IsReversed(shape, face) != w.Reversed[body]));
    Record(w, EntityKind::Face, made.GetIndex(), body,
           w.Origin[body].Of(EntityKind::Face, face.GetIndex()));
    std::vector<std::uint32_t> edges;
    STEREON_TRY(const LoopRange loops, LoopsOf(shape, face));
    for (const LoopId loop : loops) {
        STEREON_TRY_VOID(CopyLoop(w, body, loop, made, edges));
    }
    w.Faces.push_back(made);
    w.Copied.emplace_back(body, face);
    w.FaceEdges.push_back(std::move(edges));
    return {};
}

[[nodiscard]] Result<void> CopyKeptFaces(Work& w)
{
    for (const std::uint8_t body : {BodyA, BodyB}) {
        STEREON_TRY(const auto faces, AllEntities<EntityKind::Face>(w.ShapeOf(body)));
        for (const FaceId face : faces) {
            if (w.Kept[body].contains(face.GetIndex())) {
                STEREON_TRY_VOID(CopyFace(w, body, face));
            }
        }
    }
    for (const BuilderEdgeId edge : w.Cached) {
        STEREON_TRY(const CurveRef cache, w.Builder.CacheEdgeCurve(edge));
        static_cast<void>(cache);
    }
    return {};
}

// --- Shells and solids ---------------------------------------------------------------

/// The result's faces grouped into connected components across the edges
/// they share, each component's faces in order.
[[nodiscard]] std::vector<std::vector<std::size_t>> Components(const Work& w)
{
    std::vector<std::size_t> parent(w.Faces.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t i) {
        while (parent[i] != i) {
            i = parent[i] = parent[parent[i]];
        }
        return i;
    };
    std::map<std::uint32_t, std::size_t> firstUser; // result edge -> a face using it
    for (std::size_t f = 0; f < w.Faces.size(); ++f) {
        for (const std::uint32_t edge : w.FaceEdges[f]) {
            const auto [at, fresh] = firstUser.try_emplace(edge, f);
            if (!fresh) {
                parent[find(f)] = find(at->second);
            }
        }
    }
    std::map<std::size_t, std::size_t> slot; // root -> component
    std::vector<std::vector<std::size_t>> components;
    for (std::size_t f = 0; f < w.Faces.size(); ++f) {
        const auto [at, fresh] = slot.try_emplace(find(f), components.size());
        if (fresh) {
            components.emplace_back();
        }
        components[at->second].push_back(f);
    }
    return components;
}

/// Records the shell's and its solid's sources: the shells and solids of
/// the operands its faces came from.
void RecordShell(Work& w, const std::vector<std::size_t>& faces, const BuilderShellId shell,
                 const BuilderSolidId solid)
{
    for (const std::size_t f : faces) {
        const auto [body, face] = w.Copied[f];
        const Shape& shape = w.ShapeOf(body);
        const ShellId from = ShellOf(shape, face);
        if (!from.IsValid()) {
            continue;
        }
        Record(w, EntityKind::Shell, shell.GetIndex(), body,
               w.Origin[body].Of(EntityKind::Shell, from.GetIndex()));
        if (const SolidId whole = SolidOf(shape, from); whole.IsValid()) {
            Record(w, EntityKind::Solid, solid.GetIndex(), body,
                   w.Origin[body].Of(EntityKind::Solid, whole.GetIndex()));
        }
    }
}

[[nodiscard]] Result<void> AddShellOf(Work& w, const std::vector<std::size_t>& faces,
                                      const BuilderSolidId solid)
{
    STEREON_TRY(const BuilderShellId shell, w.Builder.AddShell(solid));
    for (const std::size_t f : faces) {
        STEREON_TRY_VOID(w.Builder.AddFaceToShell(shell, w.Faces[f]));
    }
    // Every shell and solid of the result has a record, a creation root
    // when no operand shell or solid continues in it.
    w.Records[{static_cast<int>(EntityKind::Shell), shell.GetIndex()}];
    w.Records[{static_cast<int>(EntityKind::Solid), solid.GetIndex()}];
    RecordShell(w, faces, shell, solid);
    return {};
}

/// The shells: a component made only of a reversed body's faces is a void
/// of the one solid the others bound; every other component bounds a solid
/// of its own.
[[nodiscard]] Result<void> MakeShells(Work& w)
{
    std::vector<std::vector<std::size_t>> lumps;
    std::vector<std::vector<std::size_t>> voids;
    for (std::vector<std::size_t>& component : Components(w)) {
        const bool hollow = std::ranges::all_of(
            component, [&](const std::size_t f) { return w.Reversed[w.Copied[f].first]; });
        (hollow ? voids : lumps).push_back(std::move(component));
    }
    if (!voids.empty() && lumps.size() != 1) {
        return Fail(ErrorCode::NotSupported,
                    "Boolean: a void among several solids, which the prototype does not place");
    }
    std::vector<BuilderSolidId> solids;
    for (const std::vector<std::size_t>& lump : lumps) {
        STEREON_TRY(const BuilderSolidId solid, w.Builder.AddSolid());
        STEREON_TRY_VOID(AddShellOf(w, lump, solid));
        solids.push_back(solid);
    }
    for (const std::vector<std::size_t>& hollow : voids) {
        STEREON_TRY_VOID(AddShellOf(w, hollow, solids.front()));
    }
    return {};
}

// --- The History -------------------------------------------------------------------

[[nodiscard]] HistoryEntry EntryOf(const EntityKind kind, const HistoryRelation relation,
                                   const Source source, const std::uint32_t target)
{
    HistoryEntry e;
    e.Kind = kind;
    e.Relation = relation;
    e.Input = source.Input;
    e.Source = source.Index;
    e.Target = target;
    return e;
}

/// An operand entity by kind, input and slot.
using SourceKey = std::tuple<int, std::uint8_t, std::uint32_t>;

/// How many result entities each operand entity continues in, and the
/// operand entities that are pure splits: continued in several result
/// entities, each continuing it alone.
struct Fanout {
    std::map<SourceKey, std::size_t> Count;
    std::set<SourceKey> Splits;
};

[[nodiscard]] Fanout FanoutOf(const Ledger& records)
{
    Fanout fanout;
    std::set<SourceKey> shared; // continued in a result entity with other sources
    for (const auto& [target, p] : records) {
        for (const Source s : p.Sources) {
            const SourceKey key{target.first, s.Input, s.Index};
            ++fanout.Count[key];
            if (p.Sources.size() >= 2) {
                shared.insert(key);
            }
        }
    }
    for (const auto& [key, count] : fanout.Count) {
        if (count >= 2 && !shared.contains(key)) {
            fanout.Splits.insert(key);
        }
    }
    return fanout;
}

/// The groups the History's split and merge events take, numbered as they
/// are met: a split's by its source, a merge's by its target.
struct Groups {
    std::map<SourceKey, std::uint32_t> Splits;
    std::uint32_t Next = 0;
};

/// The entries of one result entity: Generated by its generators, or as a
/// creation root when it has no source; Modified from each source, in its
/// split's group when the source is a pure split, in a merge group when
/// several sources continue in this entity alone.
void EntriesOf(const std::pair<int, std::uint32_t>& target, const Provenance& p,
               const Fanout& fanout, std::vector<HistoryEntry>& entries, Groups& groups)
{
    const auto kind = static_cast<EntityKind>(target.first);
    for (const Source g : p.Generators) {
        HistoryEntry e = EntryOf(kind, HistoryRelation::Generated, g, target.second);
        e.SourceKind = EntityKind::Face;
        entries.push_back(e);
    }
    if (p.Sources.empty() && p.Generators.empty()) {
        entries.push_back(EntryOf(kind, HistoryRelation::Generated, Source{}, target.second));
    }
    const bool merge = p.Sources.size() >= 2 && std::ranges::all_of(p.Sources, [&](const Source s) {
                           return fanout.Count.at({target.first, s.Input, s.Index}) == 1;
                       });
    const std::uint32_t mergeGroup = merge ? groups.Next++ : NoGroup;
    for (const Source s : p.Sources) {
        const SourceKey key{target.first, s.Input, s.Index};
        HistoryEntry e = EntryOf(kind, HistoryRelation::Modified, s, target.second);
        e.Group = mergeGroup;
        if (fanout.Splits.contains(key)) {
            const auto [at, fresh] = groups.Splits.try_emplace(key, groups.Next);
            groups.Next += fresh ? 1 : 0;
            e.Group = at->second;
        }
        entries.push_back(e);
    }
}

/// Deleted entries for every live entity of an operand no result entity
/// continues.
template <EntityKind K>
[[nodiscard]] Result<void> AddDeleted(const Shape& operand, const std::uint8_t input,
                                      const Fanout& fanout, std::vector<HistoryEntry>& entries)
{
    STEREON_TRY(const auto all, AllEntities<K>(operand));
    for (const auto id : all) {
        if (!fanout.Count.contains({static_cast<int>(K), input, id.GetIndex()})) {
            entries.push_back(
                EntryOf(K, HistoryRelation::Deleted, Source{input, id.GetIndex()}, NoEntity));
        }
    }
    return {};
}

[[nodiscard]] Result<History> HistoryOf(const Work& w, const std::array<Shape, 2>& operands)
{
    const Fanout fanout = FanoutOf(w.Records);
    std::vector<HistoryEntry> entries;
    Groups groups;
    for (const auto& [target, p] : w.Records) {
        EntriesOf(target, p, fanout, entries, groups);
    }
    for (const std::uint8_t input : {BodyA, BodyB}) {
        const Shape& operand = operands[input];
        STEREON_TRY_VOID(AddDeleted<EntityKind::Vertex>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Edge>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Coedge>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Loop>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Face>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Shell>(operand, input, fanout, entries));
        STEREON_TRY_VOID(AddDeleted<EntityKind::Solid>(operand, input, fanout, entries));
    }
    return History(std::move(entries));
}

// --- The keep rules ------------------------------------------------------------------

/// The side of the other solid a body's regions must lie on to be kept: for
/// a union outside, for a common inside, and for a cut A outside B and B
/// inside A.
[[nodiscard]] Side KeptSide(const BooleanOperation operation, const std::uint8_t body)
{
    switch (operation) {
    case BooleanOperation::Union:
        return Side::Outside;
    case BooleanOperation::Common:
        return Side::Inside;
    case BooleanOperation::Cut:
        return body == BodyA ? Side::Outside : Side::Inside;
    }
    return Side::Unresolved;
}

[[nodiscard]] std::set<std::uint32_t> KeptFaces(const std::vector<Region>& regions, const Side side)
{
    std::set<std::uint32_t> kept;
    for (const Region& region : regions) {
        if (region.Side == side) {
            for (const FaceId face : region.Faces) {
                kept.insert(face.GetIndex());
            }
        }
    }
    return kept;
}

} // namespace

Result<OperationResult> Boolean(const Shape& a, const Shape& b, const BooleanOperation operation,
                                const Context& ctx)
{
    STEREON_TRY(const ImprintedBodies bodies, ImprintBodies(a, b, ctx));
    STEREON_TRY(const Classification classes, Classify(bodies, ctx));
    Work w{bodies, ctx, {Origins(bodies.A), Origins(bodies.B)}};
    w.Reversed = {false, operation == BooleanOperation::Cut};
    w.Kept = {KeptFaces(classes.A, KeptSide(operation, BodyA)),
              KeptFaces(classes.B, KeptSide(operation, BodyB))};
    STEREON_TRY_VOID(PairPieces(w));
    STEREON_TRY_VOID(CopyKeptFaces(w));
    STEREON_TRY_VOID(MakeShells(w));
    const std::array<Shape, 2> operands{a, b};
    STEREON_TRY(History history, HistoryOf(w, operands));
    return Commit(std::move(w.Builder), ctx, std::move(history), operands);
}

} // namespace Stereon
