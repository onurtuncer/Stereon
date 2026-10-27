// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The determinism check of ADR-0007 (Verification) for the operations that
// run in parallel: batch projection onto a B-spline curve and surface
// (tests/support/DeterminismCases.hpp) and the Booleans' union, cut and
// common (tests/support/BooleanCases.hpp), serialised as the results' bit
// patterns.
//
//   stereon_determinism serial <out> [points]   one run on the serial scheduler
//   stereon_determinism <threads> <out> [points] one run on that many threads
//   stereon_determinism --check [points]         serial and 1, 2, 8, 32 threads
//                                                in this process, compared
//   stereon_determinism --report <out>           the Booleans' G2 records, as
//                                                JSON, serially
//
// The nightly workflow runs each configuration in its own process, twice,
// and compares the files byte for byte; CTest runs --check on every build.
// The G2 workflow writes the report on Linux and on Windows and compares
// the two within the vertices' tolerances (tools/compare_g2.py).

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <span>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
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

/// Everything parallel, run under the configuration: the projections' bits,
/// then the Booleans'.
std::vector<std::uint64_t> RunAll(std::span<const Stereon::Point3> points,
                                  const std::vector<BooleanCases::Case>& cases,
                                  const Context& context)
{
    const Stereon::CurveProjector curve(DeterminismTest::Curve());
    const Stereon::SurfaceProjector surface(DeterminismTest::Surface());
    std::vector<std::uint64_t> bits = DeterminismTest::Run(curve, surface, points, context);
    const std::vector<std::uint64_t> booleans = BooleanCases::Run(cases, context);
    bits.insert(bits.end(), booleans.begin(), booleans.end());
    return bits;
}

int Check(const std::size_t count)
{
    const std::vector<Stereon::Point3> points = DeterminismTest::Points(count);
    const std::vector<BooleanCases::Case> cases = BooleanCases::Cases();
    Context context;
    static_cast<void>(ContextFor("serial", context));
    const std::vector<std::uint64_t> reference = RunAll(points, cases, context);
    std::printf("serial: %016llx\n", static_cast<unsigned long long>(Digest(reference)));
    int failures = 0;
    for (const char* name : {"1", "2", "8", "32"}) {
        static_cast<void>(ContextFor(name, context));
        const std::vector<std::uint64_t> bits = RunAll(points, cases, context);
        const bool same = bits == reference;
        failures += same ? 0 : 1;
        std::printf("%s threads: %016llx %s\n", name, static_cast<unsigned long long>(Digest(bits)),
                    same ? "identical" : "DIFFERS from serial");
    }
    std::printf("%zu points projected onto a curve and a surface, %zu Booleans: %s\n", count,
                cases.size(),
                failures == 0 ? "every configuration identical" : "configurations differ");
    return failures == 0 ? 0 : 1;
}

int Report(const char* path)
{
    const std::string report = BooleanCases::Report(
        BooleanCases::Cases(), Context::Default().WithScheduler(Scheduler::Serial));
    std::FILE* file = std::fopen(path, "wb");
    if (file == nullptr) {
        return 1;
    }
    const bool ok = std::fwrite(report.data(), 1, report.size(), file) == report.size();
    std::printf("G2 records of the Booleans written to %s\n", path);
    return std::fclose(file) == 0 && ok ? 0 : 1;
}

int Main(const int argc, char** argv)
{
    if (argc >= 2 && std::strcmp(argv[1], "--check") == 0) {
        return Check(PointsFrom(argc, argv, 2));
    }
    if (argc >= 3 && std::strcmp(argv[1], "--report") == 0) {
        return Report(argv[2]);
    }
    Context context;
    if (argc < 3 || !ContextFor(argv[1], context)) {
        std::fputs("usage: stereon_determinism serial|<threads> <out> [points]\n"
                   "       stereon_determinism --check [points]\n"
                   "       stereon_determinism --report <out>\n",
                   stderr);
        return 2;
    }
    const std::size_t count = PointsFrom(argc, argv, 3);
    const std::vector<Stereon::Point3> points = DeterminismTest::Points(count);
    const std::vector<std::uint64_t> bits = RunAll(points, BooleanCases::Cases(), context);
    std::printf("%s: %zu points and the Booleans, %016llx\n", argv[1], count,
                static_cast<unsigned long long>(Digest(bits)));
    return Write(argv[2], bits) ? 0 : 1;
}

} // namespace

int main(const int argc, char** argv)
{
    try {
        return Main(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "stereon_determinism: %s\n", e.what());
        return 1;
    }
}
