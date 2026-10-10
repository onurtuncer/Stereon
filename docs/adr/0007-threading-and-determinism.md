# ADR-0007: Threading and determinism

- **Status:** Proposed
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (scope of the guarantee; fixed-order reductions; seeded randomness)
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

- **No global mutable state.** All caches, schedulers, tolerances and logging live in a `Context` passed explicitly to every operation.
- **Operation-level parallelism.** Operations express their work as sender graphs on a scheduler supplied by the `Context`; a serial scheduler is always available for debugging.
- **Determinism guarantee (G1).** For the same executable, the same hardware (same FMA availability and math library) and the same `Context` settings, results are bit-identical for any thread count and any scheduler. Cross-compiler and cross-hardware identity (G2) is a goal, not a v1.0 promise; the test suite compares G2 results within tolerance, not bitwise.
- **Fixed order, always.** Parallel stages write into per-task buffers; merges happen in a fixed order keyed by input entity index, never by completion order. Floating-point reductions are performed in that fixed order; compensated summation may be added for accuracy, but it does not replace the fixed order, because reordering a compensated sum still changes the bits.
- **Seeded randomness.** Any pseudo-random choice (ray directions in ADR-0013, symbolic perturbation, hash seeds) comes from a counter-based generator seeded by the input entity indices and the attempt number. There is no global or thread-local random state, so the choice is a pure function of the input.
- **No `thread_local` state** that can affect results.
- `std::simd` (xsimd fallback) is used inside evaluators for data parallelism.

## Consequences

- Determinism enables exact regression comparison and reliable bug reproduction.
- Fixed-order merges cost some parallel efficiency; accepted.
- Every parallel algorithm needs a determinism test.

## Verification

CI runs the corpus with 1, 2, 8 and 32 threads, and with the serial scheduler, on one build and compares serialised outputs byte for byte (G1). Pass: identical output; boolean speed-up ≥ 5× on 8 cores for bodies with more than 1 000 faces. A second job compares the GCC Linux and GCC Windows outputs within tolerance (G2) and reports, without blocking, how many are bit-identical.
