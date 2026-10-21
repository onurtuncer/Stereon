// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A public consumer: linking the kernel must not leak profiling definitions.
#if defined(STEREON_ENABLE_PROFILING) || defined(TRACY_ENABLE) || defined(TRACY_ON_DEMAND)
#error "Private profiling definitions leaked into the consumer"
#endif

#include <charconv>
#include <iostream>
#include <string_view>
#include <system_error>

import Stereon.Build;

namespace {

bool BuildPair()
{
    const auto context = Stereon::Context::Default();
    const Stereon::Point3 origin{};
    const Stereon::Vec3 x{1.0, 0.0, 0.0};
    const Stereon::Vec3 y{0.0, 1.0, 0.0};
    const auto box = Stereon::MakeBoxSolid(origin, x, y, {2.0, 3.0, 4.0}, context);
    const auto cylinder = Stereon::MakeCylinderSolid(origin, x, y, 1.0, 2.0, context);
    const auto invalid = Stereon::MakeCylinderSolid(origin, x, y, -1.0, 2.0, context);
    if (!box || !cylinder || invalid ||
        invalid.error().GetCode() != Stereon::ErrorCode::InvalidInput) {
        return false;
    }
    return box->Shape.Count<Stereon::EntityKind::Face>() == 6 &&
           cylinder->Shape.Count<Stereon::EntityKind::Face>() == 3 &&
           box->Checked == Stereon::ValidityLevel::Geometric &&
           cylinder->Checked == Stereon::ValidityLevel::Geometric;
}

} // namespace

int main(const int argc, const char* const* argv)
{
    unsigned iterations = 1;
    if (argc == 2) {
        const std::string_view text(argv[1]);
        const auto [end, error] =
            std::from_chars(text.data(), text.data() + text.size(), iterations);
        if (error != std::errc{} || end != text.data() + text.size() || iterations == 0) {
            std::cerr << "iterations must be a positive integer\n";
            return 2;
        }
    } else if (argc != 1) {
        std::cerr << "usage: stereon_profiling_workload [iterations]\n";
        return 2;
    }
    for (unsigned i = 0; i < iterations; ++i) {
        if (!BuildPair()) {
            std::cerr << "primitive workload failed\n";
            return 1;
        }
    }
    std::cout << iterations << " box/cylinder pairs validated; invalid radii rejected\n";
}
