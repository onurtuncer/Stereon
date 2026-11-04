// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The native file (docs/examples.rst): a curve and a surface saved to a
// .stn archive with the model's numerical policy, read back, and evaluated
// to the same bits. Values are written as their IEEE 754 bit patterns, so a
// round trip loses nothing.

#include <bit>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <print>

import Stereon.Io;

using namespace Stereon;

int main()
{
    const Context ctx = Context::Default();

    GeometryArchive archive;
    archive.Policy = PolicyOf(ctx);
    archive.Curves.push_back(
        Curve{*MakeCircle(Point3{1.0, 2.0, 3.0}, 4.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0})});
    archive.Surfaces.push_back(
        Surface{Cylinder{Point3{}, 2.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}}});

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "stereon-example.stn";
    if (const Result<void> saved = SaveGeometry(path, archive); !saved) {
        std::println(stderr, "save: {}", saved.error().GetMessage());
        return 1;
    }
    const Result<LoadedGeometry> loaded = LoadGeometry(path);
    std::filesystem::remove(path);
    if (!loaded) {
        std::println(stderr, "load: {}", loaded.error().GetMessage());
        return 1;
    }
    std::println("format version {}, written by {}.{}.{}, resolution {} mm", loaded->FormatVersion,
                 loaded->WrittenBy.Major, loaded->WrittenBy.Minor, loaded->WrittenBy.Patch,
                 loaded->Archive.Policy.Resolution);

    // Bit for bit: the same doubles in, the same doubles out.
    const auto bits = [](const double x) { return std::bit_cast<std::uint64_t>(x); };
    const Point3 before = Eval(archive.Curves.front(), 0.7);
    const Point3 after = Eval(loaded->Archive.Curves.front(), 0.7);
    const Point3 sBefore = Eval(archive.Surfaces.front(), Point2{0.3, 1.1});
    const Point3 sAfter = Eval(loaded->Archive.Surfaces.front(), Point2{0.3, 1.1});
    const bool same = bits(before.X) == bits(after.X) && bits(before.Y) == bits(after.Y) &&
                      bits(before.Z) == bits(after.Z) && bits(sBefore.X) == bits(sAfter.X) &&
                      bits(sBefore.Y) == bits(sAfter.Y) && bits(sBefore.Z) == bits(sAfter.Z);
    std::println("curve and surface evaluate to the same bits after the round trip: {}", same);
    return same ? 0 : 1;
}
