// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Property tests over generated bodies (ADR-0015 asks for them in topology
// code): every generated valid body freezes and reads back with the counts,
// links and multiplicities its construction implies; every generated
// single fault is refused by the rule that covers it.

#include "BodiesStd.hpp"
#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::CoedgeId;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::FaceId;
using Stereon::LoopId;
using Stereon::Shape;
using Stereon::TopologyRule;
using Faults = Stereon::Testing::TopologyFaults;

namespace {

/// The structural facts any frozen shape must satisfy, re-derived through
/// the public ranges and accessors, independently of the validator.
void RequireConsistent(const Shape& shape)
{
    std::map<CoedgeId, int> inLoops;
    std::map<CoedgeId, int> inRings;
    for (const LoopId l : *Stereon::AllEntities<EntityKind::Loop>(shape)) {
        for (const CoedgeId c : *Stereon::CoedgesOf(shape, l)) {
            RC_ASSERT(Stereon::LoopOf(shape, c) == l);
            ++inLoops[c];
        }
        const FaceId f = Stereon::FaceOf(shape, l);
        bool listed = false;
        for (const LoopId m : *Stereon::LoopsOf(shape, f)) {
            listed = listed || m == l;
        }
        RC_ASSERT(listed);
    }
    for (const EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
        for (const CoedgeId c : *Stereon::RadialCoedgesOf(shape, e)) {
            RC_ASSERT(Stereon::EdgeOf(shape, c) == e);
            ++inRings[c];
        }
    }
    const auto coedges = *Stereon::AllEntities<EntityKind::Coedge>(shape);
    RC_ASSERT(inLoops.size() == coedges.size());
    RC_ASSERT(inRings.size() == coedges.size());
    for (const CoedgeId c : coedges) {
        RC_ASSERT(inLoops[c] == 1);
        RC_ASSERT(inRings[c] == 1);
    }
    RC_ASSERT(Stereon::ValidateTopology(shape).IsValid());
}

} // namespace

TEST_CASE("generated prisms freeze with the counts of a closed polyhedron", "[topo][property]")
{
    rc::prop("prism over an n-gon", [] {
        const auto n = static_cast<std::size_t>(*rc::gen::inRange(3, 60));
        Body body = MakePrism(n, *rc::gen::inRange(1, 50) / 10.0);
        const Shape shape = std::move(body).Freeze();
        RC_ASSERT(shape.Count<EntityKind::Vertex>() == 2 * n);
        RC_ASSERT(shape.Count<EntityKind::Edge>() == 3 * n);
        RC_ASSERT(shape.Count<EntityKind::Face>() == n + 2);
        RC_ASSERT(shape.Count<EntityKind::Coedge>() == 6 * n);
        // V - E + F = 2 for a closed polyhedron of genus 0.
        RC_ASSERT(2 * n - 3 * n + (n + 2) == 2);
        for (const EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
            const auto ring = Collect(*Stereon::RadialCoedgesOf(shape, e));
            RC_ASSERT(ring.size() == 2u);
            RC_ASSERT(Stereon::StoredEdgeUseSense(shape, ring[0]) !=
                      Stereon::StoredEdgeUseSense(shape, ring[1]));
        }
        // Every face's loop is a closed chain of vertices.
        for (const LoopId l : *Stereon::AllEntities<EntityKind::Loop>(shape)) {
            const auto loop = Collect(*Stereon::CoedgesOf(shape, l));
            for (std::size_t i = 0; i < loop.size(); ++i) {
                RC_ASSERT(EndOfUse(shape, loop[i]) ==
                          StartOfUse(shape, loop[(i + 1) % loop.size()]));
            }
        }
        RequireConsistent(shape);
    });
}

TEST_CASE("generated fans keep every use of the shared edge", "[topo][property]")
{
    rc::prop("fan of k faces", [] {
        const auto k = static_cast<std::size_t>(*rc::gen::inRange(1, 24));
        Body body = MakeFan(k);
        const Shape shape = std::move(body).Freeze();
        const auto ring = Collect(*Stereon::RadialCoedgesOf(shape, shape.Resolve(body.Edges[0])));
        RC_ASSERT(ring.size() == k);
        std::set<FaceId> faces;
        for (const CoedgeId c : ring) {
            faces.insert(Stereon::FaceOf(shape, Stereon::LoopOf(shape, c)));
        }
        RC_ASSERT(faces.size() == k);
        RC_ASSERT(Collect(*Stereon::FacesOf(shape, shape.Resolve(body.Shells[0]))).size() == k);
        RequireConsistent(shape);
    });
}

TEST_CASE("generated wires have edges and no uses", "[topo][property]")
{
    rc::prop("wire of n edges", [] {
        const auto n = static_cast<std::size_t>(*rc::gen::inRange(1, 100));
        const Shape shape = MakeWire(n).Freeze();
        RC_ASSERT(shape.Count<EntityKind::Edge>() == n);
        RC_ASSERT(shape.Count<EntityKind::Coedge>() == 0u);
        for (const EdgeId e : *Stereon::AllEntities<EntityKind::Edge>(shape)) {
            RC_ASSERT(Stereon::RadialCoedgesOf(shape, e)->empty());
        }
        RequireConsistent(shape);
    });
}

TEST_CASE("every generated single fault is refused by its rule", "[topo][property]")
{
    // Each fault below breaks exactly the invariant its rule states, at a
    // random place in a random prism.
    rc::prop("one fault in a prism", [] {
        const auto n = static_cast<std::size_t>(*rc::gen::inRange(3, 30));
        Body body = MakePrism(n);
        const std::size_t coedges = body.Coedges.size();
        const auto pick = [](const std::size_t size) {
            return static_cast<std::size_t>(*rc::gen::inRange<std::size_t>(0, size));
        };
        const std::size_t c = pick(coedges);
        const int fault = *rc::gen::inRange(0, 7);
        TopologyRule expected{};
        switch (fault) {
        case 0:
            // Prev naming the coedge itself, never its predecessor in a loop
            // of three or more.
            Faults::SetPrev(body.Builder, body.Coedges[c], body.Coedges[c]);
            expected = TopologyRule::LoopLinksBroken;
            break;
        case 1:
            Faults::SetNextOutOfRange(body.Builder, body.Coedges[c]);
            expected = TopologyRule::DanglingLink;
            break;
        case 2: {
            // The ring stepping to a coedge of another edge.
            // An identical prism, frozen, says which coedges share an edge.
            std::size_t other = pick(coedges);
            const Shape probe = MakePrism(n).Freeze();
            const auto uses = Collect(*Stereon::AllEntities<EntityKind::Coedge>(probe));
            while (Stereon::EdgeOf(probe, uses[other]) == Stereon::EdgeOf(probe, uses[c])) {
                other = (other + 1) % coedges;
            }
            Faults::SetRadial(body.Builder, body.Coedges[c], body.Coedges[other]);
            expected = TopologyRule::RadialRingOpen;
            break;
        }
        case 3: {
            const std::size_t loop = pick(body.Loops.size());
            const Shape probe = MakePrism(n).Freeze();
            const auto uses = Collect(*Stereon::AllEntities<EntityKind::Coedge>(probe));
            const std::size_t own = Stereon::LoopOf(probe, uses[c]).GetIndex();
            Faults::SetCoedgeLoop(body.Builder, body.Coedges[c],
                                  body.Loops[loop == own ? (loop + 1) % body.Loops.size() : loop]);
            expected = TopologyRule::CoedgeOutsideLoop;
            break;
        }
        case 4:
            Faults::SetPCurve(body.Builder, body.Coedges[c], Stereon::PCurveRef{});
            expected = TopologyRule::MissingGeometry;
            break;
        case 5:
            body.AddLoop(body.Faces[pick(body.Faces.size())]);
            expected = TopologyRule::EmptyLoop;
            break;
        default: {
            const BuilderShellId other = Must(body.Builder.AddShell());
            Faults::SetFaceShell(body.Builder, body.Faces[pick(body.Faces.size())], other);
            expected = TopologyRule::FaceOutsideShell;
            break;
        }
        }
        const Stereon::ValidityReport report = body.Builder.Validate();
        RC_ASSERT(report.Has(expected));
        RC_ASSERT_FALSE(Stereon::Freeze(std::move(body.Builder)).has_value());
    });
}

TEST_CASE("random corruption never hangs or reads out of range", "[topo][property]")
{
    // Arbitrary link writes: the validator must terminate and, when it
    // approves, the approved store must read back consistently.
    rc::prop("random faults in a box", [] {
        Body body = MakeBox();
        const auto writes = *rc::gen::inRange(0, 12);
        const auto coedge = [&] {
            return body.Coedges[static_cast<std::size_t>(*rc::gen::inRange<std::size_t>(0, 24))];
        };
        for (int w = 0; w < writes; ++w) {
            switch (*rc::gen::inRange(0, 4)) {
            case 0:
                Faults::SetNext(body.Builder, coedge(), coedge());
                break;
            case 1:
                Faults::SetPrev(body.Builder, coedge(), coedge());
                break;
            case 2:
                Faults::SetRadial(body.Builder, coedge(), coedge());
                break;
            default:
                Faults::SetCoedgeLoop(
                    body.Builder, coedge(),
                    body.Loops[static_cast<std::size_t>(*rc::gen::inRange<std::size_t>(0, 6))]);
                break;
            }
        }
        const bool valid = body.Builder.Validate().IsValid();
        auto frozen = Stereon::Freeze(std::move(body.Builder));
        RC_ASSERT(frozen.has_value() == valid);
        if (frozen) {
            RequireConsistent(*frozen);
        }
    });
}
