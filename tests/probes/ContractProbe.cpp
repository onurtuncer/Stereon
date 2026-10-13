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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

import Stereon.Surface;

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

struct Probe {
    const char* Name;
    int (*Run)(double);
};

constexpr std::array<Probe, 7> Probes{{
    {"inline", &ProbeInline},
    {"handle", &ProbeHandle},
    {"predicate", &ProbePredicate},
    {"arena", &ProbeArena},
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
