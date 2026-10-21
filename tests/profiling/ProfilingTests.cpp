// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Profiling.hpp"

#include <array>
#include <thread>

namespace {

constexpr int Evaluate(const bool valid, [[maybe_unused]] int& evaluations)
{
    STEREON_PROFILE_SCOPE("test.evaluate");
    STEREON_PROFILE_COUNTER("test.evaluations", ++evaluations);
    if (!valid) {
        return -1;
    }
    {
        STEREON_PROFILE_SCOPE("test.nested");
        return 42;
    }
}

constexpr bool ConstantEvaluation()
{
    int evaluations = 0;
    return Evaluate(true, evaluations) == 42 && Evaluate(false, evaluations) == -1 &&
           evaluations == 0;
}

static_assert(ConstantEvaluation());

bool RuntimeEvaluation()
{
    int evaluations = 0;
    if (Evaluate(true, evaluations) != 42 || Evaluate(false, evaluations) != -1) {
        return false;
    }
#if STEREON_ENABLE_PROFILING
    return evaluations == 2;
#else
    // Disabled instrumentation must not even require these names to exist.
    STEREON_PROFILE_SCOPE(UndefinedProfileName());
    STEREON_PROFILE_COUNTER(UndefinedCounterName(), UndefinedCounterValue());
    return evaluations == 0;
#endif
}

} // namespace

int main()
{
    if (!RuntimeEvaluation()) {
        return 1;
    }
    // Separate workers exercise scope nesting and early returns with no
    // collector attached. The same executable can be run during a capture.
    std::array<bool, 4> results{};
    {
        std::array<std::jthread, 4> workers;
        for (std::size_t i = 0; i < workers.size(); ++i) {
            workers[i] = std::jthread([&, i] { results[i] = RuntimeEvaluation(); });
        }
    }
    for (const bool result : results) {
        if (!result) {
            return 1;
        }
    }
    return 0;
}
