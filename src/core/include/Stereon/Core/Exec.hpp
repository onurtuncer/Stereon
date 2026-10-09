// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Senders and schedulers (ADR-0007): Stereon::Exec.
///
/// Stereon code spells these Stereon::Exec::then, never std::execution:: or
/// stdexec:: directly, so the switch to the standard library is a change in
/// this file only. The names are the standard ones and keep their snake_case
/// spelling (ADR-0017 exempts standard-library names). Only what Stereon uses
/// is exported; add names here as algorithms need them.
///
/// This is a header, not a module partition, because of compiler bugs
/// (docs/toolchain/gcc16-feature-probe.md):
///
/// - GCC 16.2 crashes on stdexec in any translation unit compiled with
///   -fmodules, which CMake adds to every source it scans for modules.
/// - MSVC 19.50 crashes on stdexec reached through a module.
///
/// So a source file that includes this header must be excluded from module
/// scanning with stereon_exec_sources() (cmake/StereonOptions.cmake), and it
/// cannot import Stereon modules. Until the compilers are fixed, parallel code
/// lives in such files behind plain function declarations.

#pragma once

#include <version>

// Use std::execution once the standard library ships it; until then NVIDIA
// stdexec, which implements the same proposal (P2300) under namespace stdexec.
#if defined(__cpp_lib_senders)
#include <execution>
#define STEREON_EXEC_USES_STD 1
#define STEREON_EXEC_NS std::execution
#define STEREON_EXEC_SYNC_NS std::this_thread
#else
#include <stdexec/execution.hpp>
#define STEREON_EXEC_USES_STD 0
#define STEREON_EXEC_NS stdexec
#define STEREON_EXEC_SYNC_NS stdexec
#endif

namespace Stereon::Exec {

/// True when std::execution is used; false while the stdexec polyfill is.
inline constexpr bool UsesStandardLibrary = STEREON_EXEC_USES_STD == 1;

// Concepts
using STEREON_EXEC_NS::operation_state;
using STEREON_EXEC_NS::receiver;
using STEREON_EXEC_NS::scheduler;
using STEREON_EXEC_NS::sender;
using STEREON_EXEC_NS::sender_in;

// Sender factories
using STEREON_EXEC_NS::just;
using STEREON_EXEC_NS::just_error;
using STEREON_EXEC_NS::just_stopped;
using STEREON_EXEC_NS::read_env;
using STEREON_EXEC_NS::schedule;

// Sender adaptors
using STEREON_EXEC_NS::bulk;
using STEREON_EXEC_NS::continues_on;
using STEREON_EXEC_NS::into_variant;
using STEREON_EXEC_NS::let_error;
using STEREON_EXEC_NS::let_stopped;
using STEREON_EXEC_NS::let_value;
using STEREON_EXEC_NS::on;
using STEREON_EXEC_NS::starts_on;
using STEREON_EXEC_NS::stopped_as_optional;
using STEREON_EXEC_NS::then;
using STEREON_EXEC_NS::upon_error;
using STEREON_EXEC_NS::upon_stopped;
using STEREON_EXEC_NS::when_all;

// Execution policies for bulk
using STEREON_EXEC_NS::par;
using STEREON_EXEC_NS::par_unseq;
using STEREON_EXEC_NS::seq;

// Queries and low-level operations
using STEREON_EXEC_NS::connect;
using STEREON_EXEC_NS::get_env;
using STEREON_EXEC_NS::get_scheduler;
using STEREON_EXEC_NS::start;

// Execution contexts. run_loop is the serial, single-threaded context that
// ADR-0007 requires for debugging; the parallel scheduler is the default
// multi-threaded one.
using STEREON_EXEC_NS::get_parallel_scheduler;
using STEREON_EXEC_NS::parallel_scheduler;
using STEREON_EXEC_NS::run_loop;

// Consumers
using STEREON_EXEC_SYNC_NS::sync_wait;

} // namespace Stereon::Exec

#undef STEREON_EXEC_NS
#undef STEREON_EXEC_SYNC_NS
