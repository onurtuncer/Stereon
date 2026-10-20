// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Violates one contract of Stereon.Core from an importing unit, chosen by the
// argument; CTest expects the violation to be reported. Values come from argc
// so that no check can be settled at compile time. See tests/CMakeLists.txt.
//
// The probe installs its own violation handler and exits normally, because
// the default handler terminates, which CTest reports as a crash on Windows
// (exception 0xC0000409) before it looks at the output.

#include <array>
#include <contracts>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <utility>

import Stereon.Topo;

void handle_contract_violation(const std::contracts::contract_violation& violation)
{
    std::printf("contract violation: %s at %s:%u\n", violation.comment(),
                violation.location().file_name(), violation.location().line());
    std::exit(0);
}

namespace {

/// A curve on the bounded, non-periodic domain [0, 1].
struct Segment {
    Stereon::Line Base;
};
Stereon::Point3 Eval(const Segment& c, const double t)
{
    return Eval(c.Base, t);
}
Stereon::Vec3 Deriv(const Segment& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const Segment&)
{
    return Stereon::Interval{0.0, 1.0};
}
Stereon::Box3 BoundingBox(const Segment& c)
{
    return BoundingBox(c.Base);
}
double Project(const Segment& c, const Stereon::Point3& p)
{
    return Project(c.Base, p);
}
Segment Reversed(const Segment& c)
{
    return {Reversed(c.Base)};
}
Segment Transformed(const Segment& c, const Stereon::Transform& x)
{
    return {Transformed(c.Base, x)};
}

// One function per case; each violates one contract or check. value is
// 2.0, read from argc so that no check can be settled at compile time.

int ProbeInline(const double value)
{
    // Interval's constructor, a non-template inline function, requires
    // lower <= upper.
    const Stereon::Interval reversed{value + 1.0, 1.0};
    return reversed.IsEmpty() ? 1 : 0;
}

int ProbeHandle(double)
{
    // A checked-handle failure is reported from the module's own unit, so
    // it fires at every optimisation level.
    Stereon::GenerationTable<Stereon::EntityKind::Vertex> table;
    const auto handle = table.Allocate();
    table.Release(handle);
    table.Check(handle);
    return table.IsAlive(handle) ? 0 : 1;
}

int ProbePredicate(const double value)
{
    // Collinear points whose coordinates span from 2^-1000 to 2^500: no power
    // of two keeps every intermediate of the exact stage normal.
    const double tiny = 0x1p-1000 * (value - 1.0);
    const double huge = 0x1p500 * (value - 1.0);
    const Stereon::Sign sign = Stereon::Orient2D(Stereon::Point2{tiny, tiny},
                                                 Stereon::Point2{huge, huge}, Stereon::Point2{});
    return sign == Stereon::Sign::Zero ? 0 : 1;
}

int ProbeArena(const double value)
{
    // Copying an arena inside its own edit callback would let the callback's
    // reference write into a page the copy shares.
    Stereon::PagedArena<int, 4> arena;
    arena.PushBack(static_cast<int>(value));
    int shared = 0;
    arena.Modify(0, [&](int& entry) {
        const Stereon::PagedArena<int, 4> snapshot = arena;
        entry = 7;
        shared = snapshot[0];
    });
    return shared == 7 ? 0 : 1;
}

int ProbeArenaPage(const double value)
{
    // Taking a page inside the edit callback: the page is already unique, so
    // the holder would see the callback's write (audit 2026-10-20, F4).
    Stereon::PagedArena<int, 4> arena;
    arena.PushBack(static_cast<int>(value));
    std::shared_ptr<const Stereon::PagedArena<int, 4>::Page> held;
    arena.Modify(0, [&](int& entry) {
        held = arena.GetPage(0);
        entry = 7;
    });
    return (*held)[0] == 7 ? 0 : 1;
}

int ProbeEvaluatorNan(const double value)
{
    // A NaN parameter, through the erased wrapper.
    const Stereon::Curve curve(Stereon::Line{Stereon::Point3{}, Stereon::Vec3{1.0, 0.0, 0.0}});
    const double nan = std::numeric_limits<double>::quiet_NaN() * value;
    return Eval(curve, nan).X == 0.0 ? 0 : 1;
}

int ProbeEvaluatorDomain(const double value)
{
    // A parameter (2.0) outside the domain of a curve that is not periodic.
    const Stereon::Curve curve(
        Segment{Stereon::Line{Stereon::Point3{}, Stereon::Vec3{1.0, 0.0, 0.0}}});
    return Eval(curve, value).X == 0.0 ? 0 : 1;
}

int ProbeTemplate(const double value)
{
    // Normalized, a template, requires a non-zero vector.
    const Stereon::Vec3 unit = Normalized(Stereon::Vec3{value - 2.0, 0.0, 0.0});
    return unit.X == unit.X ? 0 : 1;
}

// The topology accessors (PositionOf, CurveOf, ...) treat a handle that does
// not name a live entity of the shape as a contract violation in checked
// builds (ADR-0006: the O(1) tier); the report comes from the topology
// module's own unit, so it fires at every optimisation level.

#if STEREON_CHECKED_HANDLES
Stereon::Shape Vertices(const double x, const int count, const int removed)
{
    Stereon::ShapeBuilder b;
    for (int i = 0; i < count; ++i) {
        const auto v = *b.AddVertex(Stereon::Point3{x + i, 0.0, 0.0}, 1e-7);
        if (i == removed) {
            static_cast<void>(*b.RemoveVertex(v));
        }
    }
    return *Stereon::Freeze(std::move(b));
}

int ProbeShapeForeign(const double value)
{
    // A vertex of another shape.
    const Stereon::Shape a = Vertices(value, 1, -1);
    const Stereon::Shape b = Vertices(value, 1, -1);
    const Stereon::VertexId foreign = b.GetSnapshot()->MakeId<Stereon::EntityKind::Vertex>(0);
    return PositionOf(a, foreign).X == value ? 0 : 1;
}

int ProbeShapeRange(const double value)
{
    // An index past the end of the shape's vertex table (5: two vertices).
    const Stereon::Shape a = Vertices(value, 2, -1);
    const auto past = Stereon::VertexId::Make(static_cast<std::uint32_t>(value) + 3U,
                                              a.GetSnapshot()->GetStore(), 0);
    return PositionOf(a, past).X == value ? 0 : 1;
}

int ProbeShapeDead(const double value)
{
    // The slot of a vertex the builder removed: a tombstone in the shape.
    const Stereon::Shape a = Vertices(value, 2, 1);
    const Stereon::VertexId dead = a.GetSnapshot()->MakeId<Stereon::EntityKind::Vertex>(1);
    return PositionOf(a, dead).X == value ? 0 : 1;
}
#else
int ProbeShapeForeign(double)
{
    return 2;
}
int ProbeShapeRange(double)
{
    return 2;
}
int ProbeShapeDead(double)
{
    return 2;
}
#endif

struct Probe {
    const char* Name;
    int (*Run)(double);
};

constexpr std::array<Probe, 11> Probes{{
    {"shape-foreign", &ProbeShapeForeign},
    {"shape-range", &ProbeShapeRange},
    {"shape-dead", &ProbeShapeDead},
    {"inline", &ProbeInline},
    {"handle", &ProbeHandle},
    {"predicate", &ProbePredicate},
    {"arena", &ProbeArena},
    {"arena-page", &ProbeArenaPage},
    {"evaluator-nan", &ProbeEvaluatorNan},
    {"evaluator-domain", &ProbeEvaluatorDomain},
    {"template", &ProbeTemplate},
}};

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 2) {
        return 2;
    }
    for (const Probe& probe : Probes) {
        if (std::strcmp(argv[1], probe.Name) == 0) {
            return probe.Run(static_cast<double>(argc));
        }
    }
    return 2;
}
