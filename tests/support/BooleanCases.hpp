// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Boolean fixtures of ADR-0007's checks (the Boolean prototype's M10):
// union, cut and common of two boxes overlapping at a corner, of a box and
// a cylinder through it, of a box and a loosened box 0.8 of its tolerance
// from its face (M9), and a slab drilled with a grid of holes; each result
// as its bit patterns, which the determinism checks compare, and as a G2
// record, which tools/compare_g2.py compares across platforms within the
// vertices' tolerances. Shared by tests/determinism/DeterminismRun.cpp, the
// Boolean tests and bench/BooleanBench.cpp.
//
// Include after `import Stereon.Boolean;`, with <bit>, <cstdint>, <cstdio>,
// <string>, <utility> and <vector> included before it (GCC 16.2 rejects a
// standard header first included after an import of a module that uses
// it). Failures throw std::runtime_error: the callers are a test, a
// benchmark and a command-line check, none of which can go on without the
// fixture.

#pragma once

#include <stdexcept>

namespace BooleanCases {

template <typename T>
T Must(Stereon::Result<T> result, const char* what)
{
    if (!result) {
        throw std::runtime_error(std::string(what) + ": " +
                                 std::string(result.error().GetMessage()));
    }
    return std::move(*result);
}

/// Fixtures are built on the serial scheduler: their making is not under test.
inline const Stereon::Context serial =
    Stereon::Context::Default().WithScheduler(Stereon::Scheduler::Serial);

inline Stereon::Shape BoxFrom(const Stereon::Point3& lo, const Stereon::Point3& hi,
                              const Stereon::Context& context = serial)
{
    return Must(MakeBoxSolid(Stereon::Box3::FromCorners(lo, hi), context), "a box").Shape;
}

/// The slab [0, 10n] x [0, 10n] x [0, 10].
inline Stereon::Shape Slab(const int n)
{
    return BoxFrom(Stereon::Point3{}, Stereon::Point3{10.0 * n, 10.0 * n, 10.0});
}

/// n x n cylinders of radius 2 at the centres of the slab's 10 x 10 cells,
/// from z = -5 to 15, as one operand of n^2 solids: united pairwise, a
/// balanced tree of unions, so that making it costs n^2 log n, not n^4.
inline Stereon::Shape Drills(const int n)
{
    if (n < 1) {
        throw std::runtime_error("Drills: a grid needs at least one hole");
    }
    std::vector<Stereon::Shape> level;
    for (int k = 0; k < n * n; ++k) {
        level.push_back(
            Must(MakeCylinderSolid(
                     Stereon::Point3{10.0 * (k % n) + 5.0, 10.0 * (k / n) + 5.0, -5.0},
                     Stereon::Vec3{1.0, 0.0, 0.0}, Stereon::Vec3{0.0, 1.0, 0.0}, 2.0, 20.0, serial),
                 "a drill")
                .Shape);
    }
    while (level.size() > 1) {
        std::vector<Stereon::Shape> next;
        for (std::size_t i = 0; i + 1 < level.size(); i += 2) {
            next.push_back(
                Must(Boolean(level[i], level[i + 1], Stereon::BooleanOperation::Union, serial),
                     "drills united")
                    .Shape);
        }
        if (level.size() % 2 == 1) {
            next.push_back(std::move(level.back()));
        }
        level = std::move(next);
    }
    // Checked again where the move happens: GCC at -O3 cannot see that the
    // loop above leaves one shape, and warns of a null shared pointer.
    if (level.size() != 1) {
        throw std::runtime_error("Drills: the unions did not leave one operand");
    }
    Stereon::Shape all = std::move(level[0]);
    return all;
}

/// One Boolean of the fixtures: its name, its operands and its operation.
struct Case {
    std::string Name;
    Stereon::Shape A;
    Stereon::Shape B;
    Stereon::BooleanOperation Operation;
};

/// The fixtures, each operation of each pair.
inline std::vector<Case> Cases()
{
    using Stereon::Point3;
    const Stereon::Shape box = BoxFrom(Point3{}, Point3{10.0, 10.0, 10.0});
    const Stereon::Shape corner = BoxFrom(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0});
    const Stereon::Shape drum =
        Must(MakeCylinderSolid(Point3{5.0, 5.0, -2.0}, Stereon::Vec3{1.0, 0.0, 0.0},
                               Stereon::Vec3{0.0, 1.0, 0.0}, 3.0, 15.0, serial),
             "a drum")
            .Shape;
    const double t = 1e-3;
    const Stereon::Shape near = BoxFrom(Point3{10.0 + 0.8 * t, 0.0, 0.0}, Point3{20.0, 10.0, 10.0},
                                        serial.WithResolution(t));
    std::vector<Case> cases;
    for (const Stereon::BooleanOperation op :
         {Stereon::BooleanOperation::Union, Stereon::BooleanOperation::Common,
          Stereon::BooleanOperation::Cut}) {
        const std::string name(Stereon::ToString(op));
        cases.push_back({"corner " + name, box, corner, op});
        cases.push_back({"drum " + name, box, drum, op});
        cases.push_back({"loose " + name, box, near, op});
    }
    cases.push_back({"drilled slab Cut", Slab(3), Drills(3), Stereon::BooleanOperation::Cut});
    return cases;
}

/// The bits of a result: every vertex's position and tolerance and every
/// edge's range, ends and tolerance, in index order, and the tolerance
/// changes. Comparing doubles with == would let -0.0 and 0.0 pass as equal.
inline std::vector<std::uint64_t> Bits(const Stereon::OperationResult& result)
{
    const Stereon::Shape& shape = result.Shape;
    std::vector<std::uint64_t> bits;
    for (const Stereon::VertexId v : *AllEntities<Stereon::EntityKind::Vertex>(shape)) {
        const Stereon::Point3 p = PositionOf(shape, v);
        for (const double x : {p.X, p.Y, p.Z, ToleranceOf(shape, v)}) {
            bits.push_back(std::bit_cast<std::uint64_t>(x));
        }
    }
    for (const Stereon::EdgeId e : *AllEntities<Stereon::EntityKind::Edge>(shape)) {
        const Stereon::Interval range = RangeOf(shape, e);
        for (const double x : {range.GetLower(), range.GetUpper(), ToleranceOf(shape, e)}) {
            bits.push_back(std::bit_cast<std::uint64_t>(x));
        }
        bits.push_back(StartOf(shape, e).GetIndex());
        bits.push_back(EndOf(shape, e).GetIndex());
    }
    for (const Stereon::ToleranceChange& change : result.ToleranceChanges) {
        bits.push_back(change.Index);
        bits.push_back(std::bit_cast<std::uint64_t>(change.New));
    }
    return bits;
}

/// Every case run under `context`, its bits appended in case order.
inline std::vector<std::uint64_t> Run(const std::vector<Case>& cases,
                                      const Stereon::Context& context)
{
    std::vector<std::uint64_t> bits;
    for (const Case& c : cases) {
        const std::vector<std::uint64_t> own =
            Bits(Must(Boolean(c.A, c.B, c.Operation, context), c.Name.c_str()));
        bits.push_back(own.size()); // where one case ends and the next begins
        bits.insert(bits.end(), own.begin(), own.end());
    }
    return bits;
}

template <Stereon::EntityKind K>
std::size_t CountOf(const Stereon::Shape& shape)
{
    return AllEntities<K>(shape)->size();
}

/// The G2 record of one result, one JSON object (ADR-0007, G2): its entity
/// counts, every vertex's position and tolerance in index order, with the
/// doubles in shortest round-trip form and as their bits, which the
/// comparison reads to count the results that are bit-identical.
inline std::string Record(const std::string& name, const Stereon::OperationResult& result)
{
    using Stereon::EntityKind;
    const Stereon::Shape& shape = result.Shape;
    std::string out = "{\"name\": \"" + name + "\", \"counts\": {";
    out += "\"vertex\": " + std::to_string(CountOf<EntityKind::Vertex>(shape));
    out += ", \"edge\": " + std::to_string(CountOf<EntityKind::Edge>(shape));
    out += ", \"coedge\": " + std::to_string(CountOf<EntityKind::Coedge>(shape));
    out += ", \"loop\": " + std::to_string(CountOf<EntityKind::Loop>(shape));
    out += ", \"face\": " + std::to_string(CountOf<EntityKind::Face>(shape));
    out += ", \"shell\": " + std::to_string(CountOf<EntityKind::Shell>(shape));
    out += ", \"solid\": " + std::to_string(CountOf<EntityKind::Solid>(shape));
    out += "}, \"vertices\": [";
    // Four doubles at 17 significant digits and their exponents: room for
    // all, so that -O3's -Wformat-truncation has nothing to prove.
    char buffer[256];
    bool first = true;
    for (const Stereon::VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        const Stereon::Point3 p = PositionOf(shape, v);
        std::snprintf(buffer, sizeof buffer, "%s[%.17g, %.17g, %.17g, %.17g]", first ? "" : ", ",
                      p.X, p.Y, p.Z, ToleranceOf(shape, v));
        out += buffer;
        first = false;
    }
    std::uint64_t hash = 0xCBF29CE484222325ULL;
    for (const std::uint64_t word : Bits(result)) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((word >> (8 * byte)) & 0xFFU)) * 0x100000001B3ULL;
        }
    }
    std::snprintf(buffer, sizeof buffer, "], \"bits\": \"%016llx\"}",
                  static_cast<unsigned long long>(hash));
    out += buffer;
    return out;
}

/// Every case's G2 record under `context`, as a JSON array.
inline std::string Report(const std::vector<Case>& cases, const Stereon::Context& context)
{
    std::string out = "[\n";
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Case& c = cases[i];
        out += "  " + Record(c.Name, Must(Boolean(c.A, c.B, c.Operation, context), c.Name.c_str()));
        out += i + 1 < cases.size() ? ",\n" : "\n";
    }
    return out + "]\n";
}

} // namespace BooleanCases
