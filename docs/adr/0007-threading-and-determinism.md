# ADR-0007: Threading and determinism

- **Status:** Proposed
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (scope of the guarantee; fixed-order reductions; seeded randomness)
- **Revised:** 2026-10-25 (review with the owner: the scheduler in the `Context` as GCC 16 allows it, a pool owned through the `Context`, process-wide state that never reaches a result, the G1 envelope, fixed choice sequences, a shared `Context`, and a verification CI can run; see *Review of 2026-10-25*)
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Booleans, tessellation and STEP import parallelise naturally (per face pair, per face, per entity). CAD users also need reproducibility: the same input must give the same output, or regression testing and bug reports become meaningless. Many kernels lose determinism through thread-dependent merge order, global caches and non-associative floating-point reductions.

## Options

1. **Single-threaded kernel; callers parallelise across operations.** Simple; wastes cores on large single operations.
2. **oneTBB task graph.** Mature and fast; an extra dependency outside the standard.
3. **C++26 `std::execution` (senders/receivers).** Standard, composable, supports custom schedulers; not yet shipped by standard libraries, so a polyfill (NVIDIA stdexec) is needed.

## Decision

Adopt option 3, with stdexec as the polyfill until libstdc++ ships `std::execution`.

- **No global mutable state that reaches a result.** Caches, schedulers, tolerances and logging live in a `Context` passed explicitly to every operation. Process-wide state is allowed only for identity and diagnostics, and only when it never reaches geometry, tolerances, entity ordering, serialised output, errors or a decision: today the counter that issues store identities (`Detail::NextStoreId`), the registry through which debuggers resolve handles in checked builds, and the optional profiling backend of ADR-0043. A store identity may name a snapshot; it may never order, key or appear in anything an operation returns or writes.
- **Operation-level parallelism.** Operations run their parallel work through the kernel's parallel loop (`Detail::ParallelFor`) on the `Context`'s scheduler. The `Context` carries the scheduler as a kind (`Serial` or `Parallel`) and a thread count rather than as a stdexec object, because GCC 16.2 cannot compile stdexec in module units (ADR-0016): sender graphs live in plain source files behind a header. The serial scheduler is always available for debugging and is the determinism tests' reference.
- **One pool per `Context`.** The worker threads of a parallel `Context` belong to a pool the `Context` owns through a shared handle, created on first use with its thread count: copies made with `With...` share it, and a different thread count makes a different pool. No operation creates threads of its own. The pool is a scheduler resource, not state that reaches a result. (Done in the boolean prototype's M10: `Context::GetWorkerPool()`, and `ParallelFor` takes the pool; it built one per call before.)
- **A `Context` may be shared.** A `Context` is an immutable value apart from its pool and its sinks; several threads may run operations with one `Context` at once. Its sinks are thread-safe or serialise themselves (ADR-0023), and its pool accepts work from any thread. This answers the question ADR-0020 leaves open for free-threaded Python.
- **Determinism guarantee (G1).** For the same executable on the same hardware (the same FMA availability and math library), the same input, entity indices included, and the same `Context` apart from its scheduler, thread count and sinks, results are bit-identical for any scheduler and any thread count. Entity indices are part of the input because `Compact()` renumbers them and fixed-order merges are keyed by them. Cross-compiler and cross-hardware identity (G2) is a goal, not a v1.0 promise; the test suite compares G2 results within tolerance, not bitwise. Floating-point contraction is off where results are decided (`-ffp-contract=off`, a public option of the robust library and so of every target that links it), so FMA availability does not change those bits.
- **Fixed order, always.** Parallel stages write into per-task buffers; merges happen in a fixed order keyed by input entity index, never by completion order. The parallel loop cuts its range into chunks of a fixed size whatever the thread count. Floating-point reductions are performed in that fixed order; compensated summation may be added for accuracy, but it does not replace the fixed order, because reordering a compensated sum still changes the bits.
- **Deterministic choices.** Any choice that would otherwise be random (ray directions in ADR-0013, symbolic perturbation, hash seeds) is a pure function of the input and the attempt number: either a fixed sequence (Classify's sixteen ray directions, irrational fractions of a turn) or a counter-based generator seeded by input entity indices and the attempt number. There is no global or thread-local random state.
- **No `thread_local` state** that can affect results.
- `std::simd` (xsimd fallback) is used inside evaluators for data parallelism.

## Consequences

- Determinism enables exact regression comparison and reliable bug reproduction.
- Fixed-order merges cost some parallel efficiency; accepted.
- Every parallel algorithm needs a determinism test.
- Sender code cannot appear in module units until GCC compiles stdexec there; the parallel loop's header is the boundary, and new parallel patterns extend it rather than calling stdexec from a module.

## Verification

G1, every pull request and nightly: the determinism cases (`tests/support/DeterminismCases.hpp`) run with the serial scheduler and 1, 2, 8 and 32 threads, twice each, in one process and across processes, and compare their outputs byte for byte; the TSan job runs them too. Each parallel operation adds its case. For the Boolean prototype (its M10) the cases are the many-face fixture's union, cut and common, compared through their serialised result: entity counts, positions, tolerances and the History, in entity order.

Speed-up: a Boolean on a body of more than 1 000 faces (the prototype's box with a grid of cylindrical holes) runs at least 3× faster on the 4 cores of the hosted runner of the OCCT reference workflow than with the serial scheduler, measured as that workflow measures (median of interleaved repetitions, CPU model recorded). A run whose spread is large is repeated, not quoted.

G2: a job runs the same fixtures with the GCC 16 builds on Linux and on Windows (MSYS2 UCRT64) and compares the volumes, areas and entity counts of the results within 1e-9 relative, and every vertex position within the larger of the two vertex tolerances; it reports, without blocking, how many results are bit-identical.

## Review of 2026-10-25

Reviewed with the owner against the code: the parallel loop (`Stereon/Core/Parallel.hpp`, `src/core/Parallel.cpp`), the `Context`, the determinism tests and the process-wide state in `src/`.

- **Held up:** no `thread_local` state; modelling state in the `Context`; a fixed chunking whatever the thread count, with a serial path; a G1 test on every pull request and nightly, across processes and under TSan; contraction off where results are decided.
- **Revised to match what GCC 16 allows:** the `Context` carries a scheduler kind and a thread count, not a stdexec scheduler, and sender graphs stay out of module units (ADR-0016).
- **Decided:** the worker threads belong to a pool owned through the `Context` (today's pool per call is replaced in M10); process-wide state is allowed for identity and diagnostics when it never reaches a result, and the present cases are named; the G1 envelope excludes the scheduler, the thread count and the sinks, and includes entity indices; fixed choice sequences are deterministic choices; one `Context` may be shared between threads.
- **Verification made runnable:** the hosted runners have 4 vCPUs, so the speed-up target is at least 3× on 4 cores, where the earlier text asked for 5× on 8 cores that CI cannot measure; the G2 comparison names what it compares, since shapes have no `.stn` form before v0.3.
- **Status:** stays Proposed. It is accepted once the Boolean prototype's M10 passes the verification above: thread-count identity on Booleans, the speed-up and the G2 job.
