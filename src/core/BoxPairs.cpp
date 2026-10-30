// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

module Stereon.Core;

import :Box;
import :BoxPairs;
import :Interval;

namespace Stereon {

namespace {

/// Boxes per leaf: testing a few boxes costs less than another level.
constexpr std::size_t LeafSize = 4;

/// Hits of one query that are more than one in this many of B's places are
/// put in order by marking rather than by sorting.
constexpr std::size_t DenseShare = 16;

/// Deeper than any median-split hierarchy of a list that fits in memory.
constexpr std::size_t StackDepth = 128;

[[nodiscard]] bool IsBounded(const Box3& box)
{
    return box.X.IsBounded() && box.Y.IsBounded() && box.Z.IsBounded();
}

/// The middle of a bounded interval, halved before the sum so that it
/// cannot overflow.
[[nodiscard]] double CentreOf(const Interval& i)
{
    return 0.5 * i.GetLower() + 0.5 * i.GetUpper();
}

/// A node of the hierarchy, in preorder: an inner node's first child is the
/// next node and its second is at Second; a leaf holds Count entries from
/// First.
struct Node {
    Box3 Bounds;
    std::size_t First = 0;
    std::size_t Count = 0;
    std::size_t Second = 0;
};

/// The hierarchy of the bounded, non-empty boxes of a list.
class Hierarchy {
public:
    Hierarchy(const std::span<const Box3> boxes, std::vector<std::size_t> entries)
        : m_Boxes(boxes)
        , m_Entries(std::move(entries))
    {
        if (!m_Entries.empty()) {
            m_Nodes.reserve(2 * (m_Entries.size() / LeafSize + 1));
            static_cast<void>(Build(0, m_Entries.size()));
        }
    }

    /// Calls hit(j) for every entry j whose box meets `box`, in no order.
    template <typename F>
    void Query(const Box3& box, F&& hit) const
    {
        if (m_Nodes.empty()) {
            return;
        }
        std::array<std::size_t, StackDepth> stack; // filled as the walk goes
        std::size_t top = 0;
        stack[top++] = 0;
        while (top > 0) {
            const Node& node = m_Nodes[stack[--top]];
            if (!Intersects(node.Bounds, box)) {
                continue;
            }
            if (node.Count > 0) {
                for (std::size_t k = node.First; k < node.First + node.Count; ++k) {
                    if (Intersects(m_Boxes[m_Entries[k]], box)) {
                        hit(m_Entries[k]);
                    }
                }
                continue;
            }
            stack[top++] = node.Second;
            stack[top++] = static_cast<std::size_t>(&node - m_Nodes.data()) + 1;
        }
    }

private:
    /// The node of entries [begin, end), and below it its children; returns
    /// its index.
    std::size_t Build(const std::size_t begin, const std::size_t end)
    {
        const std::size_t index = m_Nodes.size();
        m_Nodes.emplace_back();
        Box3 bounds = Box3::Empty();
        Box3 centres = Box3::Empty();
        for (std::size_t k = begin; k < end; ++k) {
            const Box3& box = m_Boxes[m_Entries[k]];
            bounds = Hull(bounds, box);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                centres[axis] = Hull(centres[axis], CentreOf(box[axis]));
            }
        }
        std::size_t axis = 0;
        for (std::size_t other = 1; other < 3; ++other) {
            if (Width(centres[other]) > Width(centres[axis])) {
                axis = other;
            }
        }
        if (end - begin <= LeafSize || !(Width(centres[axis]) > 0.0)) {
            m_Nodes[index] = Node{bounds, begin, end - begin, 0};
            return index;
        }
        // The median by centre, ties by place: a strict order, so the halves
        // are the same sets on every standard library.
        const std::size_t middle = begin + (end - begin) / 2;
        const auto offset = [](const std::size_t k) { return static_cast<std::ptrdiff_t>(k); };
        std::nth_element(m_Entries.begin() + offset(begin), m_Entries.begin() + offset(middle),
                         m_Entries.begin() + offset(end),
                         [this, axis](const std::size_t i, const std::size_t j) {
                             const double ci = CentreOf(m_Boxes[i][axis]);
                             const double cj = CentreOf(m_Boxes[j][axis]);
                             return ci < cj || (ci == cj && i < j);
                         });
        static_cast<void>(Build(begin, middle));
        const std::size_t second = Build(middle, end);
        m_Nodes[index] = Node{bounds, 0, 0, second};
        return index;
    }

    std::span<const Box3> m_Boxes;
    std::vector<std::size_t> m_Entries;
    std::vector<Node> m_Nodes;
};

/// Appends (i, j) for every hit j in increasing j. Few hits are sorted; many
/// are marked and read off in one pass over B's places, which costs less
/// than sorting once the hits are more than a small share of them, as when
/// every box meets every other.
void AppendInOrder(const std::size_t i, std::vector<std::size_t>& hits,
                   std::vector<unsigned char>& marks, std::vector<BoxPair>& pairs)
{
    if (hits.size() * DenseShare < marks.size()) {
        std::ranges::sort(hits);
        for (const std::size_t j : hits) {
            pairs.push_back({i, j});
        }
        return;
    }
    for (const std::size_t j : hits) {
        marks[j] = 1;
    }
    for (std::size_t j = 0; j < marks.size(); ++j) {
        if (marks[j] != 0) {
            marks[j] = 0;
            pairs.push_back({i, j});
        }
    }
}

} // namespace

std::vector<BoxPair> OverlappingPairs(const std::span<const Box3> a, const std::span<const Box3> b)
{
    std::vector<std::size_t> bounded;
    std::vector<std::size_t> unbounded;
    for (std::size_t j = 0; j < b.size(); ++j) {
        if (!b[j].IsEmpty()) {
            (IsBounded(b[j]) ? bounded : unbounded).push_back(j);
        }
    }
    const Hierarchy hierarchy(b, std::move(bounded));
    std::vector<BoxPair> pairs;
    std::vector<std::size_t> hits;
    std::vector<unsigned char> marks(b.size(), 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].IsEmpty()) {
            continue;
        }
        hits.clear();
        hierarchy.Query(a[i], [&hits](const std::size_t j) { hits.push_back(j); });
        for (const std::size_t j : unbounded) {
            if (Intersects(a[i], b[j])) {
                hits.push_back(j);
            }
        }
        AppendInOrder(i, hits, marks, pairs);
    }
    return pairs;
}

} // namespace Stereon
