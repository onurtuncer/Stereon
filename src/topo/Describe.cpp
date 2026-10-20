// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <string>
#include <string_view>

module Stereon.Topo;

namespace Stereon {

namespace {

using Detail::TopologyTables;

/// The diagnostic name of a curve or surface's type: "Line", "Plane",
/// "geometry" for a type that is not built in, or "empty geometry" for an
/// empty wrapper, which the validator reports and Describe must still name.
template <typename Wrapper>
[[nodiscard]] std::string_view TypeNameOf(const Wrapper& wrapper)
{
    if (wrapper.IsEmpty()) {
        return "empty geometry";
    }
    return VisitBuiltin(
        wrapper, []<typename G>(const G&) -> std::string_view { return GeometryTypeName<G>; });
}

/// The bound of every walk below: one more step than there are coedges, so
/// a walk on a corrupt ring stops.
[[nodiscard]] std::size_t WalkBound(const TopologyTables& t) noexcept
{
    return t.Coedges.GetSize() + 1;
}

/// The face of a coedge, or the invalid index when a link on the way is
/// broken.
[[nodiscard]] std::uint32_t FaceOfCoedge(const TopologyTables& t, const std::uint32_t coedge)
{
    if (!t.IsLive<EntityKind::Coedge>(coedge)) {
        return FaceId::Invalid;
    }
    const std::uint32_t loop = t.Coedges[coedge].Loop.GetIndex();
    if (!t.IsLive<EntityKind::Loop>(loop)) {
        return FaceId::Invalid;
    }
    return t.Loops[loop].Face.GetIndex();
}

[[nodiscard]] std::string IndexText(const std::uint32_t index)
{
    return index == Id<EntityKind::Vertex>::Invalid ? std::string("-") : std::to_string(index);
}

/// The edges in the vertex's incidence list, bounded and checked, so that
/// a corrupt list still describes.
[[nodiscard]] std::size_t DegreeOf(const TopologyTables& t, const std::uint32_t v)
{
    std::size_t degree = 0;
    std::uint32_t e = t.Vertices[v].FirstEdge.GetIndex();
    while (t.IsLive<EntityKind::Edge>(e) && degree < t.Edges.GetSize()) {
        ++degree;
        const Detail::EdgeRecord& edge = t.Edges[e];
        e = (edge.Start.GetIndex() == v ? edge.NextAtStart : edge.NextAtEnd).GetIndex();
    }
    return degree;
}

[[nodiscard]] std::string DescribeVertex(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::VertexRecord& v = t.Vertices[i];
    return std::format("Vertex#{} [({}, {}, {}), tol {}, {} edges]", i, v.Position.X, v.Position.Y,
                       v.Position.Z, v.Tolerance, DegreeOf(t, i));
}

[[nodiscard]] std::string DescribeEdge(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::EdgeRecord& e = t.Edges[i];
    std::string text = std::format("Edge#{} [{}", i, ToString(e.Kind));
    if (e.Curve.IsValid() && e.Curve.GetIndex() < t.Curves.GetSize()) {
        text += ' ';
        text += TypeNameOf(t.Curves[e.Curve.GetIndex()]);
    }
    if (e.Cached.IsValid()) {
        text += ", ";
        text += ToString(e.CacheStatus);
    }
    text += std::format(", tol {}, faces ", e.Tolerance);
    const std::uint32_t anchor = e.RadialAnchor.GetIndex();
    if (anchor == CoedgeId::Invalid) {
        text += "none";
    } else {
        std::uint32_t c = anchor;
        std::size_t steps = 0;
        do {
            if (steps > 0) {
                text += '|';
            }
            text += IndexText(FaceOfCoedge(t, c));
            if (!t.IsLive<EntityKind::Coedge>(c)) {
                break;
            }
            c = t.Coedges[c].Radial.GetIndex();
        } while (c != anchor && ++steps < WalkBound(t));
    }
    text += ']';
    return text;
}

[[nodiscard]] std::string DescribeCoedge(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::CoedgeRecord& c = t.Coedges[i];
    return std::format("Coedge#{} [Edge#{}, Loop#{}, {}]", i, IndexText(c.Edge.GetIndex()),
                       IndexText(c.Loop.GetIndex()),
                       c.Sense == EdgeUseSense::Same ? "same sense" : "opposite sense");
}

[[nodiscard]] std::string DescribeLoop(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::LoopRecord& l = t.Loops[i];
    const std::uint32_t anchor = l.Anchor.GetIndex();
    std::size_t count = 0;
    bool closed = anchor == CoedgeId::Invalid;
    for (std::uint32_t c = anchor; t.IsLive<EntityKind::Coedge>(c) && count < WalkBound(t);) {
        ++count;
        c = t.Coedges[c].Next.GetIndex();
        if (c == anchor) {
            closed = true;
            break;
        }
    }
    return std::format("Loop#{} [Face#{}, {} coedge{}{}]", i, IndexText(l.Face.GetIndex()), count,
                       count == 1 ? "" : "s", closed ? "" : ", open");
}

/// Members of a singly linked list, bounded by the size of the table.
template <typename Table, typename Step>
[[nodiscard]] std::size_t CountList(const Table& table, std::uint32_t first, Step step)
{
    std::size_t count = 0;
    while (first < table.GetSize() && table[first].Slot.Alive && count <= table.GetSize()) {
        ++count;
        first = step(table[first]);
    }
    return count;
}

[[nodiscard]] std::string DescribeFace(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::FaceRecord& f = t.Faces[i];
    std::string_view surface = "no surface";
    if (f.Surface.IsValid() && f.Surface.GetIndex() < t.Surfaces.GetSize()) {
        surface = TypeNameOf(t.Surfaces[f.Surface.GetIndex()]);
    }
    const std::size_t loops = CountList(t.Loops, f.FirstLoop.GetIndex(),
                                        [](const auto& l) { return l.NextInFace.GetIndex(); });
    return std::format("Face#{} [{}, {}, {} loop{}{}]", i, surface,
                       f.Shell.IsValid() ? "Shell#" + std::to_string(f.Shell.GetIndex())
                                         : std::string("no shell"),
                       loops, loops == 1 ? "" : "s", f.Reversed ? ", reversed" : "");
}

[[nodiscard]] std::string DescribeShell(const TopologyTables& t, const std::uint32_t i)
{
    const Detail::ShellRecord& s = t.Shells[i];
    const std::size_t faces = CountList(t.Faces, s.FirstFace.GetIndex(),
                                        [](const auto& f) { return f.NextInShell.GetIndex(); });
    return std::format("Shell#{} [{}, {} face{}]", i,
                       s.Solid.IsValid() ? "Solid#" + std::to_string(s.Solid.GetIndex())
                                         : std::string("no solid"),
                       faces, faces == 1 ? "" : "s");
}

[[nodiscard]] std::string DescribeSolid(const TopologyTables& t, const std::uint32_t i)
{
    const std::size_t shells = CountList(t.Shells, t.Solids[i].FirstShell.GetIndex(),
                                         [](const auto& s) { return s.NextInSolid.GetIndex(); });
    return std::format("Solid#{} [{} shell{}]", i, shells, shells == 1 ? "" : "s");
}

/// The description of slot `index` of the table of kind K, or a note that
/// it is not a live slot.
template <EntityKind K>
[[nodiscard]] std::string DescribeKind(const TopologyTables& t, const std::uint32_t index)
{
    if (!t.IsLive<K>(index)) {
        return std::format("{}#{} [no live slot]", ToString(K), IndexText(index));
    }
    if constexpr (K == EntityKind::Vertex) {
        return DescribeVertex(t, index);
    } else if constexpr (K == EntityKind::Edge) {
        return DescribeEdge(t, index);
    } else if constexpr (K == EntityKind::Coedge) {
        return DescribeCoedge(t, index);
    } else if constexpr (K == EntityKind::Loop) {
        return DescribeLoop(t, index);
    } else if constexpr (K == EntityKind::Face) {
        return DescribeFace(t, index);
    } else if constexpr (K == EntityKind::Shell) {
        return DescribeShell(t, index);
    } else {
        return DescribeSolid(t, index);
    }
}

template <EntityKind K>
[[nodiscard]] std::string DescribeIn(const Shape& shape, const Id<K> h)
{
    if (shape.IsEmpty() || !shape.GetSnapshot()->Contains(h)) {
        return std::format("{}#{} [not in this shape]", ToString(K), IndexText(h.GetIndex()));
    }
    return DescribeKind<K>(shape.GetSnapshot()->GetTables(), h.GetIndex());
}

} // namespace

namespace Detail {

std::string DescribeEntity(const TopologyTables& t, const EntityKind kind,
                           const std::uint32_t index)
{
    switch (kind) {
    case EntityKind::Vertex:
        return DescribeKind<EntityKind::Vertex>(t, index);
    case EntityKind::Edge:
        return DescribeKind<EntityKind::Edge>(t, index);
    case EntityKind::Coedge:
        return DescribeKind<EntityKind::Coedge>(t, index);
    case EntityKind::Loop:
        return DescribeKind<EntityKind::Loop>(t, index);
    case EntityKind::Face:
        return DescribeKind<EntityKind::Face>(t, index);
    case EntityKind::Shell:
        return DescribeKind<EntityKind::Shell>(t, index);
    case EntityKind::Solid:
        return DescribeKind<EntityKind::Solid>(t, index);
    }
    return {};
}

} // namespace Detail

std::string Describe(const Shape& shape, const VertexId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const EdgeId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const CoedgeId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const LoopId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const FaceId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const ShellId h)
{
    return DescribeIn(shape, h);
}
std::string Describe(const Shape& shape, const SolidId h)
{
    return DescribeIn(shape, h);
}

} // namespace Stereon
