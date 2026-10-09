// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Contract macros (ADR-0016). Modules cannot export macros, so this is a
/// header; include it in the global module fragment of a module unit.
///
/// Two kinds of contract:
///
/// - Boundary contracts (STEREON_PRE, STEREON_POST, STEREON_ASSERT) are cheap
///   checks on public entry points. They are compiled in every build and
///   evaluated with the build's semantic (STEREON_CONTRACT_SEMANTIC).
/// - Hot-path contracts (STEREON_HOT_PRE, STEREON_HOT_POST, STEREON_HOT_ASSERT)
///   guard inner loops and evaluators. They are compiled only when
///   STEREON_HOT_CONTRACTS is 1, which release builds turn off.
///
/// The global semantic is chosen per translation unit, not per contract, so
/// this split is how release builds "ignore hot-path checks and enforce
/// boundary checks". Compilers without P2900 contracts (Clang and MSVC in
/// non-blocking CI) get empty expansions, so code using the macros still builds.
///
/// Usage:
///
///     double Eval(const Line& line, double t)
///         STEREON_PRE(IsFinite(t))
///         STEREON_HOT_PRE(InDomain(line, t))
///         STEREON_POST(r : IsFinite(r));

#pragma once

#if defined(__cpp_contracts) && __cpp_contracts >= 202502L
#define STEREON_HAS_CONTRACTS 1
#else
#define STEREON_HAS_CONTRACTS 0
#endif

// Set by CMake on stereon::core; default to the safe choice elsewhere.
#if !defined(STEREON_HOT_CONTRACTS)
#define STEREON_HOT_CONTRACTS 1
#endif

#if STEREON_HAS_CONTRACTS
#define STEREON_PRE(...) pre(__VA_ARGS__)
#define STEREON_POST(...) post(__VA_ARGS__)
#define STEREON_ASSERT(...) contract_assert(__VA_ARGS__)
#else
#define STEREON_PRE(...)
#define STEREON_POST(...)
#define STEREON_ASSERT(...) static_cast<void>(sizeof(static_cast<bool>(__VA_ARGS__)))
#endif

#if STEREON_HAS_CONTRACTS && STEREON_HOT_CONTRACTS
#define STEREON_HOT_PRE(...) pre(__VA_ARGS__)
#define STEREON_HOT_POST(...) post(__VA_ARGS__)
#define STEREON_HOT_ASSERT(...) contract_assert(__VA_ARGS__)
#else
#define STEREON_HOT_PRE(...)
#define STEREON_HOT_POST(...)
#define STEREON_HOT_ASSERT(...) static_cast<void>(sizeof(static_cast<bool>(__VA_ARGS__)))
#endif
