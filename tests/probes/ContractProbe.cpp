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

#include <contracts>
#include <cstdio>
#include <cstdlib>
#include <cstring>

import Stereon.Robust;

void handle_contract_violation(const std::contracts::contract_violation& violation)
{
    std::printf("contract violation: %s at %s:%u\n", violation.comment(),
                violation.location().file_name(), violation.location().line());
    std::exit(0);
}

int main(const int argc, char** argv)
{
    if (argc != 2) {
        return 2;
    }
    const auto value = static_cast<double>(argc); // 2.0
    if (std::strcmp(argv[1], "inline") == 0) {
        // Interval's constructor, a non-template inline function, requires
        // lower <= upper.
        const Stereon::Interval reversed{value + 1.0, 1.0};
        return reversed.IsEmpty() ? 1 : 0;
    }
    if (std::strcmp(argv[1], "handle") == 0) {
        // A checked-handle failure is reported from the module's own unit,
        // so it fires at every optimisation level.
        Stereon::GenerationTable<Stereon::EntityKind::Vertex> table;
        const auto handle = table.Allocate();
        table.Release(handle);
        table.Check(handle);
        return table.IsAlive(handle) ? 0 : 1;
    }
    if (std::strcmp(argv[1], "predicate") == 0) {
        // Collinear points whose coordinates span from 2^-1000 to 2^500: no
        // power of two keeps every intermediate of the exact stage normal.
        const double tiny = 0x1p-1000 * (value - 1.0);
        const double huge = 0x1p500 * (value - 1.0);
        const Stereon::Sign sign = Stereon::Orient2D(
            Stereon::Point2{tiny, tiny}, Stereon::Point2{huge, huge}, Stereon::Point2{});
        return sign == Stereon::Sign::Zero ? 0 : 1;
    }
    if (std::strcmp(argv[1], "arena") == 0) {
        // Copying an arena inside its own edit callback would let the
        // callback's reference write into a page the copy shares.
        Stereon::PagedArena<int, 4> arena;
        arena.PushBack(argc);
        int shared = 0;
        arena.Modify(0, [&](int& entry) {
            const Stereon::PagedArena<int, 4> snapshot = arena;
            entry = 7;
            shared = snapshot[0];
        });
        return shared == 7 ? 0 : 1;
    }
    if (std::strcmp(argv[1], "template") == 0) {
        // Normalized, a template, requires a non-zero vector.
        const Stereon::Vec3 unit = Normalized(Stereon::Vec3{value - 2.0, 0.0, 0.0});
        return unit.X == unit.X ? 0 : 1;
    }
    return 2;
}
