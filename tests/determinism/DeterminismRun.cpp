// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The determinism check of ADR-0007 (Verification) for the operations that
// run in parallel: batch projection onto a B-spline curve and surface
// (tests/support/DeterminismCases.hpp), serialised as the results' bit
// patterns.
//
//   stereon_determinism serial <out> [points]   one run on the serial scheduler
//   stereon_determinism <threads> <out> [points] one run on that many threads
//   stereon_determinism --check [points]         serial and 1, 2, 8, 32 threads
//                                                in this process, compared
//
// The nightly workflow runs each configuration in its own process, twice,
// and compares the files byte for byte; CTest runs --check on every build.

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string>
#include <vector>

import Stereon.Surface;

#include "DeterminismCases.hpp"

namespace {

using Stereon::Context;
using Stereon::Scheduler;

constexpr std::size_t DefaultPoints = 20'000;

/// FNV-1a over the bytes, for the log: equal digests for equal runs.
std::uint64_t Digest(const std::vector<std::uint64_t>& bits)
{
    std::uint64_t hash = 0xCBF29CE484222325ULL;
    for (const std::uint64_t word : bits) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((word >> (8 * byte)) & 0xFFU)) * 0x100000001B3ULL;
        }
    }
    return hash;
}

/// The Context of a configuration name: "serial" or a thread count.
bool ContextFor(const char* name, Context& context)
{
    if (std::strcmp(name, "serial") == 0) {
        context = Context::Default().WithScheduler(Scheduler::Serial);
        return true;
    }
    char* end = nullptr;
    const unsigned long threads = std::strtoul(name, &end, 10);
    if (end == name || *end != '\0' || threads > 1024) {
        return false;
    }
    context = Context::Default()
                  .WithScheduler(Scheduler::Parallel)
                  .WithThreadCount(static_cast<unsigned>(threads));
    return true;
}

std::size_t PointsFrom(const int argc, char** argv, const int index)
{
    return argc > index ? static_cast<std::size_t>(std::strtoull(argv[index], nullptr, 10))
                        : DefaultPoints;
}

/// Little-endian bytes of the words, whatever the host's order, so that
/// files from two machines compare too (ADR-0007, G2).
bool Write(const char* path, const std::vector<std::uint64_t>& bits)
{
    std::FILE* file = std::fopen(path, "wb");
    if (file == nullptr) {
        return false;
    }
    bool ok = true;
    for (const std::uint64_t word : bits) {
        unsigned char bytes[8];
        for (int byte = 0; byte < 8; ++byte) {
            bytes[byte] = static_cast<unsigned char>((word >> (8 * byte)) & 0xFFU);
        }
        ok = ok && std::fwrite(bytes, 1, sizeof bytes, file) == sizeof bytes;
    }
    return std::fclose(file) == 0 && ok;
}

int Check(const std::size_t count)
{
    const std::vector<Stereon::Point3> points = DeterminismTest::Points(count);
    const Stereon::CurveProjector curve(DeterminismTest::Curve());
    const Stereon::SurfaceProjector surface(DeterminismTest::Surface());
    Context context;
    static_cast<void>(ContextFor("serial", context));
    const std::vector<std::uint64_t> reference =
        DeterminismTest::Run(curve, surface, points, context);
    std::printf("serial: %016llx\n", static_cast<unsigned long long>(Digest(reference)));
    int failures = 0;
    for (const char* name : {"1", "2", "8", "32"}) {
        static_cast<void>(ContextFor(name, context));
        const std::vector<std::uint64_t> bits =
            DeterminismTest::Run(curve, surface, points, context);
        const bool same = bits == reference;
        failures += same ? 0 : 1;
        std::printf("%s threads: %016llx %s\n", name, static_cast<unsigned long long>(Digest(bits)),
                    same ? "identical" : "DIFFERS from serial");
    }
    std::printf("%zu points projected onto a curve and a surface: %s\n", count,
                failures == 0 ? "every configuration identical" : "configurations differ");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc >= 2 && std::strcmp(argv[1], "--check") == 0) {
        return Check(PointsFrom(argc, argv, 2));
    }
    Context context;
    if (argc < 3 || !ContextFor(argv[1], context)) {
        std::fputs("usage: stereon_determinism serial|<threads> <out> [points]\n"
                   "       stereon_determinism --check [points]\n",
                   stderr);
        return 2;
    }
    const std::size_t count = PointsFrom(argc, argv, 3);
    const std::vector<Stereon::Point3> points = DeterminismTest::Points(count);
    const Stereon::CurveProjector curve(DeterminismTest::Curve());
    const Stereon::SurfaceProjector surface(DeterminismTest::Surface());
    const std::vector<std::uint64_t> bits = DeterminismTest::Run(curve, surface, points, context);
    std::printf("%s: %zu points, %016llx\n", argv[1], count,
                static_cast<unsigned long long>(Digest(bits)));
    return Write(argv[2], bits) ? 0 : 1;
}
