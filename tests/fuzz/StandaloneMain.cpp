// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A driver for the fuzz targets of tests/fuzz/ without libFuzzer, which GCC
// does not have. Each target exports LLVMFuzzerTestOneInput, the libFuzzer
// entry point, so a Clang build can link it with -fsanitize=fuzzer instead
// (STEREON_LIBFUZZER). This driver runs the target on every file of its
// seed corpus, then on mutations of them: bits flipped, bytes overwritten
// with boundary values, slices dropped, duplicated or inserted. The
// mutations come from a seeded generator, so a run is reproducible.
//
//     stereon_fuzz_<name> [--runs N] [--seed S] [--max-size B] <file or directory>...
//
// A crash, a sanitizer report or an abort by the target is a finding, as
// under libFuzzer. So is a contract violation (ADR-0006: data from outside
// the kernel must yield a KernelError, never a violation): the driver
// counts them, which is how the observe-mode run of ADR-0016 sees them, and
// exits 1 when there was one. The input that caused one is written to
// fuzz-finding-<n>.bin in the working directory.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#if defined(__cpp_contracts) || __has_include(<contracts>)
#include <contracts>
#define STEREON_FUZZ_HANDLER 1
#endif

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size);

namespace {

using Bytes = std::vector<std::uint8_t>;

std::size_t s_Violations = 0;

/// xorshift64*: small, fast and the same on every platform.
class Generator {
public:
    explicit Generator(const std::uint64_t seed)
        : m_State(seed == 0 ? 0x9E3779B97F4A7C15ULL : seed)
    {}
    std::uint64_t Next()
    {
        m_State ^= m_State >> 12;
        m_State ^= m_State << 25;
        m_State ^= m_State >> 27;
        return m_State * 0x2545F4914F6CDD1DULL;
    }
    /// A value in [0, bound), bound > 0.
    std::size_t Below(const std::size_t bound) { return static_cast<std::size_t>(Next() % bound); }

private:
    std::uint64_t m_State;
};

constexpr std::array<std::uint8_t, 8> Boundaries{0x00, 0x01, 0x7F, 0x80, 0xFE, 0xFF, 0x20, 0x40};

void Mutate(Bytes& input, Generator& g, const std::size_t maxSize)
{
    const std::size_t steps = 1 + g.Below(4);
    for (std::size_t step = 0; step < steps; ++step) {
        const std::size_t size = input.size();
        switch (g.Below(6)) {
        case 0: // flip one bit
            if (size > 0) {
                const std::size_t at = g.Below(size);
                input[at] = static_cast<std::uint8_t>(input[at] ^ (1U << g.Below(8)));
            }
            break;
        case 1: // a boundary byte
            if (size > 0) {
                input[g.Below(size)] = Boundaries[g.Below(Boundaries.size())];
            }
            break;
        case 2: // drop a slice
            if (size > 0) {
                const std::size_t at = g.Below(size);
                const std::size_t count = 1 + g.Below(std::min<std::size_t>(size - at, 64));
                input.erase(input.begin() + static_cast<std::ptrdiff_t>(at),
                            input.begin() + static_cast<std::ptrdiff_t>(at + count));
            }
            break;
        case 3: // duplicate a slice
            if (size > 0 && size < maxSize) {
                const std::size_t at = g.Below(size);
                const std::size_t count = 1 + g.Below(std::min<std::size_t>(size - at, 64));
                const Bytes slice(input.begin() + static_cast<std::ptrdiff_t>(at),
                                  input.begin() + static_cast<std::ptrdiff_t>(at + count));
                input.insert(input.begin() + static_cast<std::ptrdiff_t>(g.Below(size + 1)),
                             slice.begin(), slice.end());
            }
            break;
        case 4: // insert random bytes
            if (size < maxSize) {
                const std::size_t count = 1 + g.Below(16);
                for (std::size_t k = 0; k < count; ++k) {
                    input.insert(input.begin() +
                                     static_cast<std::ptrdiff_t>(g.Below(input.size() + 1)),
                                 static_cast<std::uint8_t>(g.Next()));
                }
            }
            break;
        default: // overwrite eight bytes, as a double or a count would be
            if (size >= 8) {
                const std::size_t at = g.Below(size - 7);
                const std::uint64_t value = g.Next();
                for (std::size_t k = 0; k < 8; ++k) {
                    input[at + k] = static_cast<std::uint8_t>(value >> (8 * k));
                }
            }
            break;
        }
    }
    if (input.size() > maxSize) {
        input.resize(maxSize);
    }
}

Bytes ReadFile(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    return Bytes{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

std::vector<Bytes> ReadCorpus(const std::vector<std::filesystem::path>& paths)
{
    std::vector<Bytes> corpus;
    for (const std::filesystem::path& path : paths) {
        if (std::filesystem::is_directory(path)) {
            std::vector<std::filesystem::path> files;
            for (const auto& entry : std::filesystem::directory_iterator(path)) {
                if (entry.is_regular_file()) {
                    files.push_back(entry.path());
                }
            }
            std::ranges::sort(files); // the same order on every platform
            for (const std::filesystem::path& file : files) {
                corpus.push_back(ReadFile(file));
            }
        } else if (std::filesystem::is_regular_file(path)) {
            corpus.push_back(ReadFile(path));
        }
    }
    return corpus;
}

/// Runs one input; true when it caused a contract violation, which is then
/// written out for reproduction.
bool Run(const Bytes& input)
{
    const std::size_t before = s_Violations;
    static_cast<void>(LLVMFuzzerTestOneInput(input.data(), input.size()));
    if (s_Violations == before) {
        return false;
    }
    const std::string name = "fuzz-finding-" + std::to_string(s_Violations) + ".bin";
    std::ofstream out(name, std::ios::binary);
    out.write(reinterpret_cast<const char*>(input.data()), // NOLINT: byte I/O
              static_cast<std::streamsize>(input.size()));
    std::fprintf(stderr, "fuzz: contract violation on an input of %zu bytes, written to %s\n",
                 input.size(), name.c_str());
    return true;
}

struct Options {
    std::size_t Runs = 1000;
    std::uint64_t Seed = 1;
    std::size_t MaxSize = 1 << 16;
    std::vector<std::filesystem::path> Paths;
};

bool Parse(const int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        const bool hasValue = i + 1 < argc;
        if (arg == "--runs" && hasValue) {
            options.Runs = std::strtoull(argv[++i], nullptr, 10);
        } else if (arg == "--seed" && hasValue) {
            options.Seed = std::strtoull(argv[++i], nullptr, 10);
        } else if (arg == "--max-size" && hasValue) {
            options.MaxSize = std::strtoull(argv[++i], nullptr, 10);
        } else if (arg.starts_with("--")) {
            return false;
        } else {
            options.Paths.emplace_back(arg);
        }
    }
    return options.MaxSize > 0;
}

} // namespace

#if defined(STEREON_FUZZ_HANDLER)
// Counts and reports; the evaluation semantic decides what follows. Under
// enforce the program terminates after the handler returns, under observe
// it continues and the count makes the run fail at the end.
void handle_contract_violation(const std::contracts::contract_violation& violation)
{
    ++s_Violations;
    std::fprintf(stderr, "contract violation: %s at %s:%u\n", violation.comment(),
                 violation.location().file_name(), violation.location().line());
}
#endif

int main(const int argc, char** argv)
{
    Options options;
    if (!Parse(argc, argv, options)) {
        std::fprintf(stderr,
                     "usage: %s [--runs N] [--seed S] [--max-size B] <file or directory>...\n",
                     argv[0]);
        return 2;
    }
    const std::vector<Bytes> corpus = ReadCorpus(options.Paths);
    std::size_t findings = 0;
    for (const Bytes& input : corpus) {
        findings += static_cast<std::size_t>(Run(input));
    }
    Generator g(options.Seed);
    for (std::size_t run = 0; run < options.Runs; ++run) {
        Bytes input = corpus.empty() ? Bytes(g.Below(256)) : corpus[g.Below(corpus.size())];
        Mutate(input, g, options.MaxSize);
        findings += static_cast<std::size_t>(Run(input));
    }
    std::printf("fuzz: %zu corpus inputs and %zu mutations, %zu contract violations\n",
                corpus.size(), options.Runs, findings);
    return findings == 0 ? 0 : 1;
}
